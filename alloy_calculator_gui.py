#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
合金热力学计算工具（UEM-Pycalphad）v3.0
==========================================

单窗口多标签页架构 + LLM 智能辅助

核心重构：
- 单窗口多标签页（6个计算类型各自独立）
- LLM 智能辅助（相预测、成分解析、参数建议、结果解读、模型推荐）
- 统一使用 CalculationResult 数据契约
- 线程安全的 Figure 数据复制
"""
import os

os.environ['SYMENGINE_LLVM'] = '0'

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import matplotlib

matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import threading
from datetime import datetime
import numpy as np
from pathlib import Path
import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import re
import traceback

from thermocal_core import (
	DatabaseManager,
	PseudoBinaryCalculator,
	ThermodynamicPropertiesCalculator,
	LiquidusSurfaceCalculator,
	SolubilityCalculator,
	LiquidusSolidusCalculator,
	CalculationResult,
	CalculationMode, TernaryCalculator
)

# LLM 模块（可选）
try:
	from llm_backend import (LLMConfig, LLMService, PROVIDER_PRESETS,
	                         fetch_ollama_models, fetch_openai_models)
	from llm_prompts import parse_phase_prediction
	LLM_AVAILABLE = True
except ImportError:
	LLM_AVAILABLE = False

# AI 代理加速（可选）
try:
	from ai_surrogate import AdaptiveSurfaceCalculator, SurrogateCache
	AI_SURROGATE_AVAILABLE = True
except ImportError:
	AI_SURROGATE_AVAILABLE = False


# =============================================================================
# 智能相过滤函数
# =============================================================================
def filter_phases_strict(db, components):
	"""智能相过滤函数：严格筛选可用相

	分两步筛选：
	1. 每个子格位必须包含至少一个选定组分或VA（基本CALPHAD条件）
	2. 所有非VA物种都在组分集合内的相直接通过；
	   否则仅保留标准溶液相（FCC_A1, BCC_A2, LIQUID 等），
	   排除含有无关元素的化合物/金属间化合物相。
	"""
	if db is None:
		return []

	allowed_species = set(c.upper() for c in components)
	allowed_species.add('VA')
	allowed_species.add('*')

	# 标准溶液相名称集合（这些相通常在多元数据库中定义了所有元素）
	_solution_names = {
		'LIQUID', 'FCC_A1', 'BCC_A2', 'HCP_A3', 'HCP_ZN',
		'DIAMOND_A4', 'CBCC_A12', 'CUB_A13', 'BCT_A5', 'DHCP'
	}

	candidate_phases = []

	for phase_name, phase_obj in db.phases.items():
		constituents = phase_obj.constituents
		if not constituents:
			continue

		# 第一步：每个子格位必须包含至少一个选定组分/VA
		sublattice_ok = True
		all_species_allowed = True
		has_any_component = False

		for sublattice in constituents:
			sublattice_has_match = False
			for species in sublattice:
				species_name = species.name.upper() if hasattr(species, 'name') else str(species).upper()
				if species_name in allowed_species:
					sublattice_has_match = True
					if species_name not in ('VA', '*'):
						has_any_component = True
				else:
					all_species_allowed = False
			if not sublattice_has_match:
				sublattice_ok = False
				break

		if not sublattice_ok or not has_any_component:
			continue

		# 第二步：如果所有物种都在允许集合内，直接通过
		if all_species_allowed:
			candidate_phases.append(str(phase_name))
			continue

		# 如果含有无关物种，仅保留标准溶液相
		name_upper = str(phase_name).upper()
		is_solution = (name_upper in _solution_names or
		               name_upper.startswith('LIQUID'))
		if is_solution:
			candidate_phases.append(str(phase_name))

	return sorted(candidate_phases)


# 相分类常量（复用 SolubilityCalculator 的经验）
_STANDARD_SOLUTION_PHASES = {
	'LIQUID', 'FCC_A1', 'BCC_A2', 'HCP_A3', 'HCP_ZN',
	'DIAMOND_A4', 'CBCC_A12', 'CUB_A13', 'BCT_A5', 'DHCP'
}
_ORDERED_SUFFIXES = ['_L12', '_L10', '_D0', '_B2', '_ORD', '_INT_']


def classify_phases(db, components):
	"""智能相分类：将候选相按类型分组

	先用 filter_phases_strict 筛选候选相，再按类型分为：
	- solution:  溶液相（FCC_A1, BCC_A2, LIQUID 等）
	- compound:  化合物/金属间化合物（相物种 ⊆ 组分集合）
	- ordered:   排序相（含 _L12, _L10, _D0 等后缀）
	- other:     其他
	"""
	candidates = filter_phases_strict(db, components)
	comp_set = set(c.upper() for c in components if c.upper() != 'VA')

	groups = {'solution': [], 'compound': [], 'ordered': [], 'other': []}

	for phase_name in candidates:
		name_upper = phase_name.upper()

		# 1. 检查排序相
		is_ordered = any(suf in name_upper for suf in _ORDERED_SUFFIXES)
		if is_ordered:
			groups['ordered'].append(phase_name)
			continue

		# 2. 检查溶液相
		is_solution = (name_upper in _STANDARD_SOLUTION_PHASES or
		               name_upper.startswith('LIQUID'))
		if is_solution:
			groups['solution'].append(phase_name)
			continue

		# 3. 检查化合物（相物种 ⊆ 组分集合）
		try:
			phase_obj = db.phases[phase_name]
			phase_species = set()
			for sublattice in phase_obj.constituents:
				for species in sublattice:
					sp_name = species.name if hasattr(species, 'name') else str(species)
					if sp_name.upper() not in ('VA', '*', 'VACANCY'):
						phase_species.add(sp_name.upper())
			if phase_species and phase_species.issubset(comp_set):
				groups['compound'].append(phase_name)
				continue
		except Exception:
			pass

		# 4. 其他
		groups['other'].append(phase_name)

	return groups


# =============================================================================
# 相选择弹窗类（增强版：分组显示 + AI 推荐预选）
# =============================================================================

# 分组显示名称与颜色
_GROUP_LABELS = {
	'solution': ('溶液相', '#2E7D32'),
	'compound': ('化合物', '#1565C0'),
	'ordered':  ('排序相', '#E65100'),
	'other':    ('其他',   '#616161'),
}

# 默认选中规则
_GROUP_DEFAULTS = {
	'solution': True,
	'compound': True,
	'ordered':  False,
	'other':    False,
}


class PhaseSelectorDialog(tk.Toplevel):
	"""相选择弹窗：分组显示候选相，支持 AI 推荐预选"""

	def __init__(self, parent, phases, title="相选择确认",
	             default_selected=True, recommended_phases=None,
	             phase_groups=None):
		super().__init__(parent)
		self.title(title)
		self.result = None
		self.phases = phases
		self.phase_vars = {}
		self.recommended_phases = recommended_phases
		self.phase_groups = phase_groups  # {'solution': [...], 'compound': [...], ...}

		self.transient(parent)
		self.grab_set()
		self.resizable(True, True)

		num_phases = len(phases)
		height = min(700, max(450, 150 + num_phases * 22))
		self.geometry(f"520x{height}")

		self._create_widgets(default_selected)

		self.update_idletasks()
		x = parent.winfo_x() + (parent.winfo_width() - self.winfo_width()) // 2
		y = parent.winfo_y() + (parent.winfo_height() - self.winfo_height()) // 2
		self.geometry(f"+{x}+{y}")

		self.wait_window(self)

	def _create_widgets(self, default_selected):
		"""创建弹窗界面"""
		main_frame = ttk.Frame(self, padding=15)
		main_frame.pack(fill=tk.BOTH, expand=True)

		ttk.Label(main_frame, text="请确认参与计算的相",
		          font=('', 14, 'bold')).pack(pady=(0, 5))

		if self.recommended_phases:
			ttk.Label(main_frame,
			          text=f"AI 推荐了 {len(self.recommended_phases)} 个相（已预选），"
			               f"共 {len(self.phases)} 个候选相。",
			          font=('', 10), foreground='#0066cc', justify=tk.LEFT
			          ).pack(anchor=tk.W, pady=(0, 5))

		# 说明文字
		if self.phase_groups:
			n_sol = len(self.phase_groups.get('solution', []))
			n_cpd = len(self.phase_groups.get('compound', []))
			n_ord = len(self.phase_groups.get('ordered', []))
			desc = f"已智能分类 {len(self.phases)} 个候选相：" \
			       f"{n_sol}溶液相 + {n_cpd}化合物 + {n_ord}排序相"
			ttk.Label(main_frame, text=desc,
			          font=('', 10), foreground='gray').pack(anchor=tk.W, pady=(0, 5))
		else:
			ttk.Label(main_frame,
			          text=f"共 {len(self.phases)} 个候选相，请勾选参与平衡计算的相：",
			          font=('', 10), foreground='gray').pack(anchor=tk.W, pady=(0, 5))

		# 快捷按钮
		btn_frame = ttk.Frame(main_frame)
		btn_frame.pack(fill=tk.X, pady=(0, 8))

		ttk.Button(btn_frame, text="全选", width=8,
		           command=self._select_all).pack(side=tk.LEFT, padx=2)
		ttk.Button(btn_frame, text="全不选", width=8,
		           command=self._deselect_all).pack(side=tk.LEFT, padx=2)

		if self.phase_groups:
			ttk.Button(btn_frame, text="仅溶液相", width=10,
			           command=self._select_solution_only).pack(side=tk.LEFT, padx=2)
			ttk.Button(btn_frame, text="溶液+化合物", width=12,
			           command=self._select_sol_compound).pack(side=tk.LEFT, padx=2)
		else:
			ttk.Button(btn_frame, text="仅液相", width=8,
			           command=self._select_liquid_only).pack(side=tk.LEFT, padx=2)

		ttk.Button(btn_frame, text="反选", width=6,
		           command=self._invert_selection).pack(side=tk.LEFT, padx=2)

		if self.recommended_phases:
			ttk.Button(btn_frame, text="AI推荐", width=8,
			           command=self._select_recommended).pack(side=tk.LEFT, padx=2)

		# 滚动区域
		list_frame = ttk.Frame(main_frame)
		list_frame.pack(fill=tk.BOTH, expand=True)

		canvas = tk.Canvas(list_frame, highlightthickness=0)
		scrollbar = ttk.Scrollbar(list_frame, orient="vertical", command=canvas.yview)
		scrollable_frame = ttk.Frame(canvas)

		scrollable_frame.bind("<Configure>",
		                      lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

		canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
		canvas.configure(yscrollcommand=scrollbar.set)

		canvas.pack(side="left", fill="both", expand=True)
		scrollbar.pack(side="right", fill="y")

		def _on_mousewheel(event):
			canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
		canvas.bind("<MouseWheel>", _on_mousewheel)
		scrollable_frame.bind("<MouseWheel>", _on_mousewheel)

		# 填充相列表
		if self.phase_groups:
			self._create_grouped_list(scrollable_frame, default_selected)
		else:
			self._create_flat_list(scrollable_frame, default_selected)

		# 底部
		ttk.Separator(main_frame, orient='horizontal').pack(fill=tk.X, pady=8)

		bottom_frame = ttk.Frame(main_frame)
		bottom_frame.pack(fill=tk.X)

		self.count_label = ttk.Label(bottom_frame, text="", font=('', 10))
		self.count_label.pack(side=tk.LEFT)
		self._update_count()

		ttk.Button(bottom_frame, text="取消", width=10,
		           command=self._on_cancel).pack(side=tk.RIGHT, padx=5)
		ttk.Button(bottom_frame, text="确定计算", width=12,
		           command=self._on_confirm).pack(side=tk.RIGHT, padx=5)

		for var in self.phase_vars.values():
			var.trace_add('write', lambda *args: self._update_count())

		self.bind('<Return>', lambda e: self._on_confirm())
		self.bind('<Escape>', lambda e: self._on_cancel())

	def _create_grouped_list(self, parent, default_selected):
		"""按分组创建相列表"""
		row_offset = 0
		for group_key in ['solution', 'compound', 'ordered', 'other']:
			group_phases = self.phase_groups.get(group_key, [])
			if not group_phases:
				continue

			label_text, color = _GROUP_LABELS.get(group_key, ('其他', '#616161'))
			group_default = _GROUP_DEFAULTS.get(group_key, False)

			# 分组标题
			header = tk.Frame(parent)
			header.grid(row=row_offset, column=0, columnspan=3, sticky=tk.W,
			            padx=2, pady=(8, 2))
			tk.Label(header, text=f"■ {label_text} ({len(group_phases)})",
			         font=('', 10, 'bold'), fg=color).pack(side=tk.LEFT)

			# 分组全选/全不选
			def _toggle_group(phases=group_phases, select=True):
				for p in phases:
					if p in self.phase_vars:
						self.phase_vars[p].set(select)

			tk.Button(header, text="✓", width=2, relief='flat', fg='green',
			          command=lambda p=group_phases: _toggle_group(p, True)).pack(side=tk.LEFT, padx=2)
			tk.Button(header, text="✗", width=2, relief='flat', fg='red',
			          command=lambda p=group_phases: _toggle_group(p, False)).pack(side=tk.LEFT, padx=1)

			row_offset += 1

			# 相列表（3列）
			num_cols = 3
			for i, phase in enumerate(sorted(group_phases)):
				if self.recommended_phases:
					selected = phase in self.recommended_phases
				else:
					selected = group_default

				var = tk.BooleanVar(value=selected)
				self.phase_vars[phase] = var

				row = row_offset + i // num_cols
				col = i % num_cols

				cb = ttk.Checkbutton(parent, text=phase, variable=var)
				cb.grid(row=row, column=col, sticky=tk.W, padx=8, pady=1)

			row_offset += (len(group_phases) + num_cols - 1) // num_cols

	def _create_flat_list(self, parent, default_selected):
		"""平铺显示（兼容无分组调用）"""
		num_cols = 3
		for i, phase in enumerate(self.phases):
			if self.recommended_phases:
				selected = phase in self.recommended_phases
			else:
				selected = default_selected

			var = tk.BooleanVar(value=selected)
			self.phase_vars[phase] = var

			row = i // num_cols
			col = i % num_cols

			cb = ttk.Checkbutton(parent, text=phase, variable=var)
			cb.grid(row=row, column=col, sticky=tk.W, padx=5, pady=2)

	def _update_count(self):
		count = sum(1 for var in self.phase_vars.values() if var.get())
		self.count_label.config(text=f"已选择: {count}/{len(self.phases)} 个相")

	def _select_all(self):
		for var in self.phase_vars.values():
			var.set(True)

	def _deselect_all(self):
		for var in self.phase_vars.values():
			var.set(False)

	def _select_liquid_only(self):
		for phase, var in self.phase_vars.items():
			var.set('LIQUID' in phase.upper())

	def _select_solution_only(self):
		"""仅选中溶液相"""
		sol_set = set(self.phase_groups.get('solution', []))
		for phase, var in self.phase_vars.items():
			var.set(phase in sol_set)

	def _select_sol_compound(self):
		"""选中溶液相 + 化合物"""
		active = set(self.phase_groups.get('solution', []) +
		             self.phase_groups.get('compound', []))
		for phase, var in self.phase_vars.items():
			var.set(phase in active)

	def _invert_selection(self):
		for var in self.phase_vars.values():
			var.set(not var.get())

	def _select_recommended(self):
		if self.recommended_phases:
			for phase, var in self.phase_vars.items():
				var.set(phase in self.recommended_phases)

	def _on_confirm(self):
		selected = [phase for phase, var in self.phase_vars.items() if var.get()]
		if not selected:
			messagebox.showwarning("警告", "请至少选择一个相！", parent=self)
			return
		self.result = selected
		self.destroy()

	def _on_cancel(self):
		self.result = None
		self.destroy()


# =============================================================================
# 计算标签页抽象基类
# =============================================================================
class CalculatorTab:
	"""计算标签页基类，提供共享基础设施"""

	def __init__(self, parent_notebook, app, title):
		self.app = app
		self.frame = ttk.Frame(parent_notebook)
		parent_notebook.add(self.frame, text=title)

		self.calc_result = None
		self.results_data = {}
		self._results_lock = threading.Lock()

		self._build_layout()

	def _build_layout(self):
		"""构建标签页内部布局"""
		paned = ttk.PanedWindow(self.frame, orient=tk.VERTICAL)
		paned.pack(fill=tk.BOTH, expand=True)

		# 上部：参数条 + 图表区
		top_frame = ttk.Frame(paned)
		paned.add(top_frame, weight=3)

		# 参数条
		params_frame = ttk.Frame(top_frame)
		params_frame.pack(fill=tk.X, padx=5, pady=5)
		self._create_params_bar(params_frame)

		# 图表区
		chart_frame = ttk.Frame(top_frame)
		chart_frame.pack(fill=tk.BOTH, expand=True)
		self._create_chart_area(chart_frame)

		# 下部：数据 + 日志子标签
		bottom_frame = ttk.Frame(paned)
		paned.add(bottom_frame, weight=1)
		self._create_sub_notebook(bottom_frame)

	def _create_params_bar(self, parent):
		"""创建参数条（子类实现）"""
		pass

	def _create_chart_area(self, parent):
		"""创建图表区"""
		self.fig = Figure(figsize=(10, 7))
		self.ax = self.fig.add_subplot(111)
		self.canvas = FigureCanvasTkAgg(self.fig, parent)
		self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
		self.toolbar = NavigationToolbar2Tk(self.canvas, parent)
		self.toolbar.update()

	def _create_sub_notebook(self, parent):
		"""创建数据+日志子标签"""
		sub_nb = ttk.Notebook(parent)
		sub_nb.pack(fill=tk.BOTH, expand=True)

		# 数据页
		data_frame = ttk.Frame(sub_nb)
		sub_nb.add(data_frame, text="数据")

		data_toolbar = ttk.Frame(data_frame)
		data_toolbar.pack(fill=tk.X, pady=2, padx=5)
		ttk.Button(data_toolbar, text="导出数据", command=self.export_data).pack(side=tk.LEFT, padx=2)
		ttk.Button(data_toolbar, text="导出图形", command=self.export_plot).pack(side=tk.LEFT, padx=2)
		ttk.Button(data_toolbar, text="清空", command=lambda: self.data_text.delete('1.0', tk.END)).pack(side=tk.LEFT, padx=2)

		self.data_text = scrolledtext.ScrolledText(data_frame, wrap=tk.NONE,
		                                           font=('Courier New', 10), height=8)
		self.data_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=2)

		# 日志页
		log_frame = ttk.Frame(sub_nb)
		sub_nb.add(log_frame, text="日志")

		log_toolbar = ttk.Frame(log_frame)
		log_toolbar.pack(fill=tk.X, pady=2, padx=5)
		ttk.Button(log_toolbar, text="清空日志",
		           command=lambda: self.log_text.delete('1.0', tk.END)).pack(side=tk.LEFT, padx=2)

		self.log_text = scrolledtext.ScrolledText(log_frame, wrap=tk.WORD,
		                                          font=('Courier New', 9), height=8)
		self.log_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=2)

	def log(self, message):
		"""线程安全日志"""
		timestamp = datetime.now().strftime("%H:%M:%S")
		full_message = f"[{timestamp}] {message}\n"

		def _update():
			try:
				if self.frame.winfo_exists():
					self.log_text.insert(tk.END, full_message)
					self.log_text.see(tk.END)
			except Exception:
				pass

		try:
			self.frame.after(0, _update)
		except Exception:
			pass

	def _get_model_spec(self, model_key, phases=None):
		"""获取模型规格"""
		if model_key == 'RKM':
			return None
		elif model_key == 'UEM1':
			try:
				from pycalphad.models.model_uem import ModelUEM
				from pycalphad import Model
				uem1_liq = getattr(self, 'uem1_liquid_only', None)
				if uem1_liq and uem1_liq.get() and phases:
					return {ph: ModelUEM if 'LIQUID' in ph.upper() else Model
					        for ph in phases}
				else:
					return ModelUEM
			except ImportError:
				self.log("警告: UEM1模型不可用")
				return None
		else:
			return self.app.available_models.get(model_key)

	def _parse_base_alloy(self, alloy_str):
		"""解析基础合金字符串"""
		alloy_str_upper = alloy_str.strip().upper().replace(" ", "")
		if not alloy_str_upper:
			raise ValueError("基础合金字符串不能为空")

		# 尝试 LLM 解析
		if hasattr(self.app, 'llm_service') and self.app.llm_service.is_available():
			try:
				result = self.app.llm_service.parse_composition(
					alloy_str, self.app.available_comps)
				if result:
					self.log(f"AI 解析成分: {alloy_str} → {result}")
					return result
			except Exception:
				pass

		# 标准解析
		matches = re.findall(r'([A-Z]{1,2})(\d*\.?\d*)', alloy_str_upper)

		if not matches:
			raise ValueError(f"无法解析基础合金: '{alloy_str}'")

		base_comps_dict = {}
		available_comps_upper = [c.upper() for c in self.app.available_comps]

		for comp, ratio_str in matches:
			if comp not in available_comps_upper:
				if comp == 'VA':
					continue
				raise ValueError(f"组分 '{comp}' 不在数据库中！")

			ratio = 1.0 if ratio_str == "" else float(ratio_str)
			if ratio <= 0:
				raise ValueError(f"组分 '{comp}' 的比例必须大于 0")

			base_comps_dict[comp] = base_comps_dict.get(comp, 0) + ratio

		if not base_comps_dict:
			raise ValueError("基础合金未包含任何有效组分")

		return base_comps_dict

	def _rebuild_chart(self, projection=None):
		"""重建图表区域"""
		parent = self.canvas.get_tk_widget().master
		self.toolbar.destroy()
		self.canvas.get_tk_widget().destroy()

		self.fig = Figure(figsize=(10, 7))
		if projection:
			self.ax = self.fig.add_subplot(111, projection=projection)
		else:
			self.ax = self.fig.add_subplot(111)

		self.canvas = FigureCanvasTkAgg(self.fig, parent)
		self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
		self.toolbar = NavigationToolbar2Tk(self.canvas, parent)
		self.toolbar.update()

	def _copy_binplot_to_ax(self, result, ax, label='', scan_comp=''):
		"""从后端 Figure 复制 binplot 数据到指定 axes"""
		source_fig = result.figure
		if source_fig is None or len(source_fig.axes) == 0:
			return

		source_ax = source_fig.axes[0]

		try:
			# 1. 复制填充区域 (Collections) - 背景层
			for collection in source_ax.collections:
				try:
					paths = collection.get_paths()
					facecolors = collection.get_facecolors()
					edgecolors = collection.get_edgecolors()
					alpha = collection.get_alpha()

					for i, path in enumerate(paths):
						vertices = path.vertices
						if len(vertices) < 3:
							continue
						fc = facecolors[i % len(facecolors)] if len(facecolors) > 0 else 'none'
						ec = edgecolors[i % len(edgecolors)] if len(edgecolors) > 0 else 'none'
						patch = mpatches.Polygon(vertices, facecolor=fc, edgecolor=ec,
						                         alpha=alpha, zorder=1)
						ax.add_patch(patch)
				except Exception:
					pass

			# 2. 复制线条 (Line2D) - 中间层
			for line in source_ax.get_lines():
				try:
					xdata = line.get_xdata()
					ydata = line.get_ydata()
					if len(xdata) == 0:
						continue
					ax.plot(xdata, ydata, color=line.get_color(),
					        linestyle=line.get_linestyle(),
					        linewidth=line.get_linewidth(),
					        marker=line.get_marker(),
					        label=line.get_label(), zorder=10)
				except Exception:
					pass

			# 3. 复制文本 (Text) - 顶层
			for text in source_ax.texts:
				try:
					content = text.get_text()
					if not content or content.strip() == '':
						continue
					pos = text.get_position()
					bbox_prop = None
					bbox = text.get_bbox_patch()
					if bbox:
						bbox_prop = dict(boxstyle=bbox.get_boxstyle(),
						                 facecolor=bbox.get_facecolor(),
						                 edgecolor=bbox.get_edgecolor(),
						                 alpha=bbox.get_alpha())
					ax.text(pos[0], pos[1], content,
					        fontsize=text.get_fontsize(),
					        color=text.get_color(),
					        ha=text.get_ha(), va=text.get_va(),
					        rotation=text.get_rotation(),
					        bbox=bbox_prop, zorder=100)
				except Exception:
					pass

			# 4. 复制坐标轴设置
			try:
				ax.set_xlim(source_ax.get_xlim())
				ax.set_ylim(source_ax.get_ylim())
				ax.set_xlabel(source_ax.get_xlabel() or f'X({scan_comp})')
				ax.set_ylabel(source_ax.get_ylabel() or 'Temperature (K)')
				ax.set_title(source_ax.get_title() or f'真二元相图 - {label}')
			except Exception:
				pass
		except Exception as e:
			self.log(f"binplot 复制出错: {e}")

	def _plot_phase_map_to_ax(self, result, ax, label='', scan_comp=''):
		"""绘制 phase_map 到指定 axes（平滑版 + 自动标注，支持整数编码）"""
		x_grid = result.x_axis
		t_grid = result.t_grid
		phase_map = result.phase_map
		phase_legend = getattr(result, 'phase_legend', None)

		if x_grid is None or t_grid is None or phase_map is None:
			self.log(f"警告: {label} 的 phase_map 数据不完整")
			return

		try:
			# 根据是否有 phase_legend 区分整数编码 vs 旧字符串格式
			if phase_legend is not None:
				# 整数编码格式：phase_map 为 int32，phase_legend 为 {int: str}
				unique_ids = sorted(set(phase_map.flatten()) - {-1, -2, -3})
				unique_phases = [phase_legend.get(pid, f"Phase-{pid}") for pid in unique_ids]
				# 重映射到连续 0..N-1
				remap = {pid: idx for idx, pid in enumerate(unique_ids)}
				num_map = np.zeros(phase_map.shape, dtype=float)
				for r in range(phase_map.shape[0]):
					for c in range(phase_map.shape[1]):
						num_map[r, c] = remap.get(int(phase_map[r, c]), -1)
			else:
				# 旧字符串格式（向后兼容）
				unique_phases_raw = np.unique(phase_map)
				unique_phases = [p for p in unique_phases_raw if p != 'Error']
				if len(unique_phases) == 0:
					ax.text(0.5, 0.5, "计算失败：所有点均返回Error",
					        ha='center', va='center', transform=ax.transAxes,
					        fontsize=14, color='red',
					        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
					return
				phase_to_num = {phase: i for i, phase in enumerate(unique_phases)}
				num_map = np.zeros(phase_map.shape, dtype=float)
				for r in range(phase_map.shape[0]):
					for c in range(phase_map.shape[1]):
						num_map[r, c] = phase_to_num.get(phase_map[r, c], -1)

			if len(unique_phases) == 0:
				ax.text(0.5, 0.5, "计算失败：所有点均返回Error",
				        ha='center', va='center', transform=ax.transAxes,
				        fontsize=14, color='red',
				        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
				return

			# 平滑处理
			try:
				from scipy.ndimage import zoom
				zoom_factor = 5
				num_map_smooth = zoom(num_map, zoom_factor, order=0)
				x_smooth = np.linspace(x_grid.min(), x_grid.max(), num_map_smooth.shape[1])
				t_smooth = np.linspace(t_grid.min(), t_grid.max(), num_map_smooth.shape[0])
				X, T = np.meshgrid(x_smooth, t_smooth)
				Z = num_map_smooth
			except ImportError:
				X, T = np.meshgrid(x_grid, t_grid)
				Z = num_map
				x_smooth = x_grid
				t_smooth = t_grid

			if len(unique_phases) == 1:
				ax.pcolormesh(X, T, Z, cmap='Blues', alpha=0.5, shading='auto')
				ax.text(0.5, 0.5, unique_phases[0], ha='center', va='center',
				        transform=ax.transAxes, fontsize=16, fontweight='bold',
				        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
			else:
				levels = np.arange(len(unique_phases) + 1) - 0.5
				n_colors = len(unique_phases)
				if n_colors <= 10:
					cmap = plt.cm.get_cmap('tab10', n_colors)
				elif n_colors <= 20:
					cmap = plt.cm.get_cmap('tab20', n_colors)
				else:
					cmap = plt.cm.get_cmap('gist_ncar', n_colors)

				ax.contourf(X, T, Z, levels=levels, cmap=cmap, alpha=0.9)
				ax.contour(X, T, Z, levels=levels, colors='k', linewidths=0.5, alpha=0.5)

				import matplotlib.patheffects as PathEffects

				# 收集所有标签位置，防止重叠
				label_positions = []  # [(x, y, w, h), ...]
				x_span = x_smooth[-1] - x_smooth[0]
				t_span = t_smooth[-1] - t_smooth[0]
				lbl_w = x_span * 0.10
				lbl_h = t_span * 0.04

				for i, phase_name in enumerate(unique_phases):
					mask = (Z == i)
					if not np.any(mask):
						continue
					y_indices, x_indices = np.where(mask)

					# 质心
					center_idx_y = int(np.median(y_indices))
					center_idx_x = int(np.median(x_indices))
					pos_x = x_smooth[center_idx_x]
					pos_y = t_smooth[center_idx_y]

					# 避免重叠：检查已有标签位置
					overlaps = True
					best_x, best_y = pos_x, pos_y
					candidates = [(pos_x, pos_y)]
					# 尝试四分位
					for q in [0.3, 0.7, 0.2, 0.8]:
						qx = x_smooth[int(np.quantile(x_indices, q))]
						qy = t_smooth[int(np.quantile(y_indices, q))]
						candidates.append((qx, qy))

					for cx, cy in candidates:
						no_overlap = True
						for ox, oy, ow, oh in label_positions:
							if (abs(cx - ox) < (lbl_w + ow) / 2 and
							    abs(cy - oy) < (lbl_h + oh) / 2):
								no_overlap = False
								break
						if no_overlap:
							best_x, best_y = cx, cy
							overlaps = False
							break

					if overlaps:
						best_x, best_y = pos_x, pos_y

					display_name = phase_name
					if len(display_name) > 15 and '+' in display_name:
						display_name = display_name.replace(' + ', '\n+ ')

					# 根据区域大小调整字体
					n_pts = np.sum(mask)
					total_pts = Z.size
					ratio = n_pts / total_pts
					fontsize = 8 if ratio < 0.05 else (9 if ratio < 0.15 else 10)

					txt = ax.text(best_x, best_y, display_name,
					              fontsize=fontsize, ha='center', va='center',
					              color='black', fontweight='bold', zorder=100,
					              bbox=dict(boxstyle='round,pad=0.2',
					                        facecolor='white', alpha=0.75,
					                        edgecolor='none'))
					txt.set_path_effects([PathEffects.withStroke(linewidth=2.5,
					                                            foreground='white', alpha=0.9)])

					label_positions.append((best_x, best_y, lbl_w, lbl_h))

			ax.set_xlabel(f'X({scan_comp})', fontsize=12)
			ax.set_ylabel('Temperature (K)', fontsize=12)
			ax.set_title(f'伪二元相区分布图 ({label})', fontsize=14, fontweight='bold')

		except Exception as e:
			self.log(f"phase_map 绘制失败: {e}")
			self.log(traceback.format_exc())

	def export_plot(self):
		"""导出图表"""
		if not hasattr(self, 'fig') or not self.fig.axes:
			messagebox.showwarning("警告", "没有可导出的图表！")
			return

		file_path = filedialog.asksaveasfilename(
			title="保存图表", defaultextension=".png",
			filetypes=[("PNG", "*.png"), ("SVG", "*.svg"), ("PDF", "*.pdf"), ("All", "*.*")]
		)
		if file_path:
			try:
				self.fig.savefig(file_path, dpi=300, bbox_inches='tight')
				self.log(f"图表已导出到: {file_path}")
			except Exception as e:
				messagebox.showerror("错误", f"导出图表失败:\n{e}")

	def export_data(self):
		"""导出数据（子类可覆盖）"""
		messagebox.showinfo("提示", "当前标签页暂无可导出的数据")

	def clear_results(self):
		"""清除结果"""
		self._rebuild_chart()
		self.log_text.delete('1.0', tk.END)
		self.data_text.delete('1.0', tk.END)
		self.calc_result = None
		self.results_data.clear()
		self.canvas.draw()
		self.log("结果已清除")

	def start_calculation(self):
		"""开始计算（子类实现）"""
		pass


# =============================================================================
# 液相线/固相线标签页
# =============================================================================
class LiquidusSolidusTab(CalculatorTab):
	"""液相线/固相线计算标签页 - 支持多模型对比"""

	def __init__(self, parent_notebook, app):
		self.model_vars = {}
		self.uem1_liquid_only = tk.BooleanVar(value=False)
		super().__init__(parent_notebook, app, "液相线/固相线")

	def _get_selected_models(self):
		"""获取当前标签页选中的模型列表"""
		selected = []
		for key, var in self.model_vars.items():
			if var.get():
				if self.app.available_models[key] is None and key != 'RKM':
					pass
				else:
					selected.append(key)
		if not selected:
			from tkinter import messagebox
			messagebox.showwarning("警告", "请至少选择一个计算模型！")
		return selected

	def _create_params_bar(self, parent):
		"""创建参数条"""
		frame = ttk.LabelFrame(parent, text="液相线/固相线参数", padding=5)
		frame.pack(fill=tk.X)
		frame.columnconfigure(1, weight=1)
		frame.columnconfigure(3, weight=1)
		frame.columnconfigure(5, weight=1)

		# 第一行
		ttk.Label(frame, text="基础合金:").grid(row=0, column=0, sticky=tk.W, padx=3)
		self.comps_entry = ttk.Entry(frame, width=15)
		self.comps_entry.grid(row=0, column=1, sticky=tk.EW, padx=3)
		self.comps_entry.insert(0, "AL")

		ttk.Label(frame, text="变化组分:").grid(row=0, column=2, sticky=tk.W, padx=3)
		self.scan_comp_combobox = ttk.Combobox(frame, width=8, state='normal')
		self.scan_comp_combobox.grid(row=0, column=3, sticky=tk.W, padx=3)

		ttk.Label(frame, text="范围:").grid(row=0, column=4, sticky=tk.W, padx=3)
		range_f = ttk.Frame(frame)
		range_f.grid(row=0, column=5, sticky=tk.W, padx=3)
		self.scan_start = ttk.Entry(range_f, width=5)
		self.scan_start.pack(side=tk.LEFT, padx=1)
		self.scan_start.insert(0, "0.0")
		ttk.Label(range_f, text="~").pack(side=tk.LEFT)
		self.scan_end = ttk.Entry(range_f, width=5)
		self.scan_end.pack(side=tk.LEFT, padx=1)
		self.scan_end.insert(0, "1.0")
		ttk.Label(range_f, text="共").pack(side=tk.LEFT, padx=2)
		self.scan_points = ttk.Entry(range_f, width=4)
		self.scan_points.pack(side=tk.LEFT, padx=1)
		self.scan_points.insert(0, "51")
		ttk.Label(range_f, text="点").pack(side=tk.LEFT)

		# 第二行
		ttk.Label(frame, text="温度(K):").grid(row=1, column=0, sticky=tk.W, padx=3)
		temp_f = ttk.Frame(frame)
		temp_f.grid(row=1, column=1, columnspan=3, sticky=tk.W, padx=3)
		self.temp_min = ttk.Entry(temp_f, width=6)
		self.temp_min.pack(side=tk.LEFT, padx=1)
		self.temp_min.insert(0, "873")
		ttk.Label(temp_f, text="~").pack(side=tk.LEFT)
		self.temp_max = ttk.Entry(temp_f, width=6)
		self.temp_max.pack(side=tk.LEFT, padx=1)
		self.temp_max.insert(0, "2073")
		ttk.Label(temp_f, text="步长").pack(side=tk.LEFT, padx=3)
		self.temp_step = ttk.Entry(temp_f, width=4)
		self.temp_step.pack(side=tk.LEFT, padx=1)
		self.temp_step.insert(0, "10")
		ttk.Label(temp_f, text="K").pack(side=tk.LEFT)

		ttk.Button(frame, text="开始计算",
		           command=self.start_calculation).grid(row=1, column=4, columnspan=2,
		                                               sticky=tk.EW, padx=5, pady=3)

		# 第三行 — 模型选择（多选复选框，支持多模型对比）
		ttk.Label(frame, text="模型:").grid(row=2, column=0, sticky=tk.W, padx=3)
		model_f = ttk.Frame(frame)
		model_f.grid(row=2, column=1, columnspan=4, sticky=tk.W, padx=3, pady=2)
		for key in ('RKM', 'Muggianu', 'Toop', 'UEM1'):
			var = tk.BooleanVar(value=(key == 'RKM'))
			self.model_vars[key] = var
			ttk.Checkbutton(model_f, text=self.app.model_labels[key],
			                variable=var).pack(side=tk.LEFT, padx=4)
		self.cb_uem1_liq = ttk.Checkbutton(model_f, text="UEM Liq Only",
		                                   variable=self.uem1_liquid_only,
		                                   state='disabled')
		self.cb_uem1_liq.pack(side=tk.LEFT, padx=4)
		self.model_vars['UEM1'].trace_add('write', self._on_uem1_toggle)

	def _on_uem1_toggle(self, *args):
		if self.model_vars['UEM1'].get():
			self.cb_uem1_liq.config(state='normal')
		else:
			self.cb_uem1_liq.config(state='disabled')
			self.uem1_liquid_only.set(False)

	def update_comps(self):
		"""更新组分下拉列表"""
		self.scan_comp_combobox['values'] = self.app.available_comps
		if self.app.available_comps:
			self.scan_comp_combobox.current(0)

	def start_calculation(self):
		"""液相线/固相线计算"""
		if not self.app.dbe:
			messagebox.showwarning("警告", "请先加载数据库！")
			return

		try:
			base_str = self.comps_entry.get().strip()
			base_dict = self._parse_base_alloy(base_str)
			vary_comp = self.scan_comp_combobox.get().strip().upper()

			if not vary_comp:
				messagebox.showwarning("警告", "请选择变化组分！")
				return

			selected_comps = self.app.get_selected_components()
			if not selected_comps:
				messagebox.showwarning("警告", "请至少选择一个组分！")
				return

			for comp in base_dict.keys():
				if comp not in selected_comps:
					selected_comps.append(comp)
			if vary_comp not in selected_comps:
				selected_comps.append(vary_comp)

			all_comps = sorted(list(set(selected_comps + ['VA'])))

			x_scan = np.linspace(float(self.scan_start.get()),
			                     float(self.scan_end.get()),
			                     int(self.scan_points.get()))

			selected_keys = self._get_selected_models()
			if not selected_keys:
				return

			phase_groups = classify_phases(self.app.dbe, selected_comps)
			candidate_phases = sum(phase_groups.values(), [])
			if not candidate_phases:
				messagebox.showerror("错误", "没有找到适用于当前组分的相！")
				return

			# AI 推荐相
			recommended = None
			if hasattr(self.app, 'llm_service') and self.app.llm_service.is_available():
				recommended = self.app.llm_service.predict_phases(
					selected_comps, candidate_phases, 'liquidus_solidus')

			dialog = PhaseSelectorDialog(self.app.root, candidate_phases,
			                             title="液相线计算 - 相选择",
			                             recommended_phases=recommended,
			                             phase_groups=phase_groups)
			if dialog.result is None:
				self.log("用户取消了计算")
				return

			final_phases = dialog.result
			self.log(f"用户确认使用 {len(final_phases)} 个相")

			inputs = {
				'base_comps_dict': base_dict,
				'varying_comp': vary_comp,
				'study_comps': all_comps,
				'selected_comps': selected_comps,
				'x_scan_range': x_scan,
				'temp_min': float(self.temp_min.get()),
				'temp_max': float(self.temp_max.get()),
				'temp_step': float(self.temp_step.get()),
				'scan_comp': vary_comp,
			}

			threading.Thread(
				target=self._run_liquidus_thread,
				args=(inputs, selected_keys, final_phases),
				daemon=True
			).start()

		except Exception as e:
			messagebox.showerror("输入错误", f"参数解析失败:\n{e}")

	def _run_liquidus_thread(self, inputs, models, phases):
		"""液固相线计算线程"""
		self.results_data['x_scan'] = inputs['x_scan_range']
		self.results_data['scan_comp'] = inputs['scan_comp']

		T_min = inputs['temp_min']
		T_max = inputs['temp_max']
		T_step = inputs['temp_step']

		compositions = []
		base_comp = inputs['base_comps_dict']
		varying_comp = inputs['varying_comp']

		for x_val in inputs['x_scan_range']:
			comp = dict(base_comp)
			remaining = 1.0 - x_val
			if base_comp:
				total_base = sum(base_comp.values())
				for elem, frac in base_comp.items():
					comp[elem] = (frac / total_base) * remaining
			comp[varying_comp] = x_val
			compositions.append(comp)

		self.log(f"成分点数: {len(compositions)}, 温度范围: {T_min}-{T_max}K")
		active_components = inputs.get('study_comps')
		if not active_components:
			active_components = sorted(list(set(inputs['selected_comps'] + ['VA'])))

		for model_key in models:
			self.log(f"开始计算 {self.app.model_labels[model_key]} 模型...")
			try:
				model_spec = self._get_model_spec(model_key, phases)
				calculator = LiquidusSolidusCalculator(
					database=self.app.dbe, logger=self.log)

				result = calculator.calculate_batch(
					compositions=compositions,
					T_range=(T_min, T_max),
					components=active_components,
					phases=phases,
					model_spec=model_spec,
					pdens=self.app.pdens,
					T_step=T_step,
					refine=True,
					use_bisection=True,
					bisection_tol=1.0,
					verbose=True
				)

				if not result.success:
					self.log(f"✗ {self.app.model_labels[model_key]} 计算失败: {result.message}")
					continue

				with self._results_lock:
					self.results_data[model_key] = {
						'type': 'liquidus_solidus',
						'result': result,
						'label': self.app.model_labels[model_key]
					}

				self.log(f"✓ {self.app.model_labels[model_key]} 计算完成")
				self.log(f"  {result.message}")

			except Exception as e:
				self.log(f"✗ {self.app.model_labels[model_key]} 计算失败: {e}")
				self.log(traceback.format_exc())

		if self.frame.winfo_exists():
			self.frame.after_idle(self._plot_results)

		# LLM 结果解读
		self._interpret_results()

	def _interpret_results(self):
		"""LLM 解读计算结果"""
		if not hasattr(self.app, 'llm_service') or not self.app.llm_service.is_available():
			return
		try:
			for key, data in self.results_data.items():
				if key in ('x_scan', 'scan_comp'):
					continue
				result = data.get('result')
				if result and result.success:
					context = {'calc_type': '液相线/固相线', 'model': data.get('label', key)}
					interpretation = self.app.llm_service.interpret_result(
						result.to_dict() if hasattr(result, 'to_dict') else {}, context)
					if interpretation:
						self.app.llm_panel.show_message(f"AI 分析 ({data.get('label', key)}):\n{interpretation}")
					break
		except Exception:
			pass

	def _plot_results(self):
		"""绘制结果"""
		while len(self.fig.axes) > 1:
			self.fig.delaxes(self.fig.axes[-1])
		self.ax.clear()

		x = self.results_data.get('x_scan')
		scan_comp = self.results_data.get('scan_comp', 'C')

		if x is None:
			return

		has_special_plot = False

		for model_key in self.model_vars:
			if model_key not in self.results_data:
				continue

			data = self.results_data[model_key]
			result = data.get('result')
			if result is None:
				continue

			label = data.get('label', model_key)
			mode = result.mode

			if mode == CalculationMode.BINPLOT.value:
				has_special_plot = True
				self._copy_binplot_to_ax(result, self.ax, label, scan_comp)
			elif mode == CalculationMode.PHASE_MAP.value:
				has_special_plot = True
				self._plot_phase_map_to_ax(result, self.ax, label, scan_comp)
			elif mode == CalculationMode.LINE.value:
				if not has_special_plot:
					self._plot_line_mode(result, label, x)

		if not has_special_plot:
			self.ax.set_xlabel(f'X({scan_comp})', fontsize=12)
			self.ax.set_ylabel('Temperature (K)', fontsize=12)
			self.ax.set_title('液相线/固相线对比', fontsize=14, fontweight='bold')

		handles, labels = self.ax.get_legend_handles_labels()
		if handles:
			filtered = [(h, l) for h, l in zip(handles, labels) if not l.startswith('_')]
			if filtered:
				self.ax.legend(*zip(*filtered), loc='best')

		self.ax.grid(True, alpha=0.3)
		self.canvas.draw()

	def _plot_line_mode(self, result, label, x):
		"""绘制 line 模式"""
		liquidus_data = result.liquidus_data
		if liquidus_data is None:
			return

		liquidus = liquidus_data.get('liquidus_K')
		solidus = liquidus_data.get('solidus_K')

		if isinstance(liquidus, (int, float)) and isinstance(solidus, (int, float)):
			self.log(f"  {label}: 液相线={liquidus:.1f}K, 固相线={solidus:.1f}K")
			return

		batch_results = liquidus_data.get('batch_results')
		if batch_results:
			liquidus_arr = np.array([r.get('liquidus_K', np.nan) for r in batch_results])
			solidus_arr = np.array([r.get('solidus_K', np.nan) for r in batch_results])

			if len(liquidus_arr) == len(x):
				self.ax.plot(x, liquidus_arr, '-', label=f"{label} - Liquidus", linewidth=2)
				self.ax.plot(x, solidus_arr, '--', label=f"{label} - Solidus", linewidth=2)

	def export_data(self):
		"""导出液相线数据"""
		if not self.results_data:
			messagebox.showinfo("提示", "暂无数据可导出")
			return

		file_path = filedialog.asksaveasfilename(
			title="保存数据", defaultextension=".csv",
			filetypes=[("CSV", "*.csv"), ("All", "*.*")])
		if not file_path:
			return

		try:
			with open(file_path, 'w', encoding='utf-8') as f:
				f.write(f"# 液相线/固相线计算结果\n")
				f.write(f"# {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

				for key, data in self.results_data.items():
					if key in ('x_scan', 'scan_comp'):
						continue
					if 'result' not in data:
						continue

					result = data['result']
					label = data.get('label', key)

					if result.mode == CalculationMode.BINPLOT.value:
						f.write(f"\n# {label} (Binplot)\n")
						f.write("Line_ID,Label,X,Temperature_K\n")
						if result.figure and result.figure.axes:
							ax = result.figure.axes[0]
							for i, line in enumerate(ax.get_lines()):
								ll = line.get_label()
								if not ll or ll.startswith('_'):
									ll = f"Boundary_{i}"
								for x, y in zip(line.get_xdata(), line.get_ydata()):
									f.write(f"{i},{ll},{x:.6f},{y:.6f}\n")

					elif result.mode == CalculationMode.LINE.value:
						f.write(f"\n# {label} (Line)\n")
						ld = result.liquidus_data
						if ld and ld.get('batch_results'):
							f.write("X,Liquidus_K,Solidus_K\n")
							x_scan = self.results_data.get('x_scan', [])
							for i, br in enumerate(ld['batch_results']):
								xi = x_scan[i] if i < len(x_scan) else i
								f.write(f"{xi:.6f},{br.get('liquidus_K', 'NaN'):.6f},"
								        f"{br.get('solidus_K', 'NaN'):.6f}\n")

			self.log(f"数据已导出到: {file_path}")
		except Exception as e:
			messagebox.showerror("错误", f"导出失败:\n{e}")


# =============================================================================
# 热力学性质标签页
# =============================================================================
class ThermodynamicPropsTab(CalculatorTab):
	"""热力学性质计算标签页 - 支持多模型对比"""

	def __init__(self, parent_notebook, app):
		self.model_vars = {}
		self.uem1_liquid_only = tk.BooleanVar(value=False)
		super().__init__(parent_notebook, app, "热力学性质")

	def _get_selected_models(self):
		"""获取当前标签页选中的模型列表"""
		selected = []
		for key, var in self.model_vars.items():
			if var.get():
				if self.app.available_models[key] is None and key != 'RKM':
					pass
				else:
					selected.append(key)
		if not selected:
			from tkinter import messagebox
			messagebox.showwarning("警告", "请至少选择一个计算模型！")
		return selected

	def _create_chart_area(self, parent):
		"""创建多图表区（Gibbs、混合焓、活度）"""
		self.chart_notebook = ttk.Notebook(parent)
		self.chart_notebook.pack(fill=tk.BOTH, expand=True)

		# Gibbs
		f1 = ttk.Frame(self.chart_notebook)
		self.chart_notebook.add(f1, text="Gibbs自由能")
		self.fig_gibbs = Figure(figsize=(10, 6))
		self.ax_gibbs = self.fig_gibbs.add_subplot(111)
		self.canvas_gibbs = FigureCanvasTkAgg(self.fig_gibbs, f1)
		self.canvas_gibbs.get_tk_widget().pack(fill=tk.BOTH, expand=True)
		NavigationToolbar2Tk(self.canvas_gibbs, f1)

		# 混合焓
		f2 = ttk.Frame(self.chart_notebook)
		self.chart_notebook.add(f2, text="混合焓")
		self.fig_enthalpy = Figure(figsize=(10, 6))
		self.ax_enthalpy = self.fig_enthalpy.add_subplot(111)
		self.canvas_enthalpy = FigureCanvasTkAgg(self.fig_enthalpy, f2)
		self.canvas_enthalpy.get_tk_widget().pack(fill=tk.BOTH, expand=True)
		NavigationToolbar2Tk(self.canvas_enthalpy, f2)

		# 活度
		f3 = ttk.Frame(self.chart_notebook)
		self.chart_notebook.add(f3, text="活度")
		self.fig_activity = Figure(figsize=(10, 6))
		self.ax_activity = self.fig_activity.add_subplot(111)
		self.canvas_activity = FigureCanvasTkAgg(self.fig_activity, f3)
		self.canvas_activity.get_tk_widget().pack(fill=tk.BOTH, expand=True)
		NavigationToolbar2Tk(self.canvas_activity, f3)

		# 保持 base class 兼容
		self.fig = self.fig_gibbs
		self.ax = self.ax_gibbs
		self.canvas = self.canvas_gibbs

	def _create_params_bar(self, parent):
		"""创建参数条"""
		frame = ttk.LabelFrame(parent, text="热力学性质参数", padding=5)
		frame.pack(fill=tk.X)

		row0 = ttk.Frame(frame)
		row0.pack(fill=tk.X, pady=2)

		ttk.Label(row0, text="基础合金:").pack(side=tk.LEFT, padx=3)
		self.comps_entry = ttk.Entry(row0, width=12)
		self.comps_entry.pack(side=tk.LEFT, padx=3)
		self.comps_entry.insert(0, "AL")

		ttk.Label(row0, text="变化组分:").pack(side=tk.LEFT, padx=3)
		self.scan_comp_combobox = ttk.Combobox(row0, width=8)
		self.scan_comp_combobox.pack(side=tk.LEFT, padx=3)

		ttk.Label(row0, text="范围:").pack(side=tk.LEFT, padx=3)
		self.scan_start = ttk.Entry(row0, width=5)
		self.scan_start.pack(side=tk.LEFT, padx=1)
		self.scan_start.insert(0, "0.0")
		ttk.Label(row0, text="~").pack(side=tk.LEFT)
		self.scan_end = ttk.Entry(row0, width=5)
		self.scan_end.pack(side=tk.LEFT, padx=1)
		self.scan_end.insert(0, "1.0")
		ttk.Label(row0, text="共").pack(side=tk.LEFT, padx=2)
		self.scan_points = ttk.Entry(row0, width=4)
		self.scan_points.pack(side=tk.LEFT, padx=1)
		self.scan_points.insert(0, "51")
		ttk.Label(row0, text="点").pack(side=tk.LEFT)

		ttk.Label(row0, text="温度(K):").pack(side=tk.LEFT, padx=(10, 3))
		self.temperature = ttk.Entry(row0, width=6)
		self.temperature.pack(side=tk.LEFT, padx=1)
		self.temperature.insert(0, "1273")

		# 第二行：模型选择 + 开始计算
		row1 = ttk.Frame(frame)
		row1.pack(fill=tk.X, pady=2)
		ttk.Label(row1, text="模型:").pack(side=tk.LEFT, padx=3)
		for key in ('RKM', 'Muggianu', 'Toop', 'UEM1'):
			var = tk.BooleanVar(value=(key == 'RKM'))
			self.model_vars[key] = var
			ttk.Checkbutton(row1, text=self.app.model_labels[key],
			                variable=var).pack(side=tk.LEFT, padx=4)
		self.cb_uem1_liq = ttk.Checkbutton(row1, text="UEM Liq Only",
		                                   variable=self.uem1_liquid_only,
		                                   state='disabled')
		self.cb_uem1_liq.pack(side=tk.LEFT, padx=4)
		self.model_vars['UEM1'].trace_add('write', self._on_uem1_toggle)
		ttk.Button(row1, text="开始计算",
		           command=self.start_calculation).pack(side=tk.RIGHT, padx=5)

	def _on_uem1_toggle(self, *args):
		if self.model_vars['UEM1'].get():
			self.cb_uem1_liq.config(state='normal')
		else:
			self.cb_uem1_liq.config(state='disabled')
			self.uem1_liquid_only.set(False)

	def update_comps(self):
		self.scan_comp_combobox['values'] = self.app.available_comps
		if self.app.available_comps:
			self.scan_comp_combobox.current(0)

	def start_calculation(self):
		"""热力学性质计算"""
		if not self.app.dbe:
			messagebox.showwarning("警告", "请先加载数据库！")
			return

		try:
			base_dict = self._parse_base_alloy(self.comps_entry.get().strip())
			vary_comp = self.scan_comp_combobox.get().strip().upper()

			selected_comps = self.app.get_selected_components()
			if not selected_comps:
				messagebox.showwarning("警告", "请至少选择一个组分！")
				return

			for comp in base_dict.keys():
				if comp not in selected_comps:
					selected_comps.append(comp)
			if vary_comp and vary_comp not in selected_comps:
				selected_comps.append(vary_comp)

			x_scan = np.linspace(float(self.scan_start.get()),
			                     float(self.scan_end.get()),
			                     int(self.scan_points.get()))

			selected_keys = self._get_selected_models()
			if not selected_keys:
				return

			phase_groups = classify_phases(self.app.dbe, selected_comps)
			candidate_phases = sum(phase_groups.values(), [])
			if not candidate_phases:
				messagebox.showerror("错误", "没有找到适用于当前组分的相！")
				return

			# AI 推荐相
			recommended = None
			if hasattr(self.app, 'llm_service') and self.app.llm_service.is_available():
				recommended = self.app.llm_service.predict_phases(
					selected_comps, candidate_phases, 'thermodynamic_properties')

			dialog = PhaseSelectorDialog(self.app.root, candidate_phases,
			                             title="热力学性质 - 相选择",
			                             recommended_phases=recommended,
			                             phase_groups=phase_groups)
			if dialog.result is None:
				self.log("用户取消了计算")
				return

			final_phases = dialog.result
			temp = float(self.temperature.get())

			inputs = {
				'base_comps_dict': base_dict,
				'varying_comp': vary_comp,
				'x_scan_range': x_scan,
				'temperature': temp,
				'selected_comps': selected_comps,
			}

			threading.Thread(
				target=self._run_props_thread,
				args=(inputs, selected_keys, final_phases),
				daemon=True
			).start()

		except Exception as e:
			messagebox.showerror("输入错误", f"参数解析失败:\n{e}")

	def _run_props_thread(self, inputs, models, phases):
		"""热力学性质计算线程"""
		for model_key in models:
			self.log(f"计算 {self.app.model_labels[model_key]} 热力学性质...")
			try:
				model_spec = self._get_model_spec(model_key, phases)
				calculator = ThermodynamicPropertiesCalculator(
					database=self.app.dbe, logger=self.log)

				result = calculator.calculate(
					base_composition=inputs['base_comps_dict'],
					varying_comp=inputs['varying_comp'],
					x_range=(inputs['x_scan_range'][0], inputs['x_scan_range'][-1]),
					x_points=len(inputs['x_scan_range']),
					temperature=inputs['temperature'],
					components=self.app.available_comps,
					phases=phases,
					model_spec=model_spec,
					pdens=self.app.pdens,
					verbose=True
				)

				if not result.success:
					self.log(f"✗ {self.app.model_labels[model_key]} 失败: {result.message}")
					continue

				props_data = result.properties_data
				with self._results_lock:
					self.results_data[f"{model_key}_props"] = {
						'gibbs': props_data.get('gibbs_energy'),
						'enthalpy_mix': props_data.get('enthalpy_mix'),
						'activity': props_data.get('activity'),
						'x_scan': props_data.get('x_values', inputs['x_scan_range']),
						'label': self.app.model_labels[model_key]
					}

				gibbs = props_data.get('gibbs_energy')
				if gibbs is not None:
					valid_gibbs = ~np.isnan(gibbs)
					if np.any(valid_gibbs):
						self.log(f"✓ {self.app.model_labels[model_key]} 完成:")
						self.log(f"  Gibbs能: {np.nanmin(gibbs):.2e} ~ {np.nanmax(gibbs):.2e} J/mol")

			except Exception as e:
				self.log(f"✗ {self.app.model_labels[model_key]} 失败: {e}")
				self.log(traceback.format_exc())

		if self.frame.winfo_exists():
			self.frame.after_idle(self._plot_props)

	def _plot_props(self):
		"""绘制热力学性质图"""
		self.ax_gibbs.clear()
		self.ax_enthalpy.clear()
		self.ax_activity.clear()

		for model_key in self.model_vars:
			key = f"{model_key}_props"
			if key not in self.results_data:
				continue
			d = self.results_data[key]
			label = d['label']
			x_scan = d['x_scan']

			if d.get('gibbs') is not None:
				self.ax_gibbs.plot(x_scan, d['gibbs'] / 1000.0, label=label, linewidth=2)
			if d.get('enthalpy_mix') is not None:
				self.ax_enthalpy.plot(x_scan, d['enthalpy_mix'] / 1000.0, label=label, linewidth=2)
			if d.get('activity') is not None:
				for comp, act_values in d['activity'].items():
					self.ax_activity.plot(x_scan, act_values,
					                      label=f"{label} - a({comp})", linewidth=2)

		self.ax_gibbs.set_xlabel('Composition', fontsize=12)
		self.ax_gibbs.set_ylabel('Gibbs Energy (kJ/mol)', fontsize=12)
		self.ax_gibbs.set_title('Gibbs自由能对比', fontsize=14, fontweight='bold')
		self.ax_gibbs.legend()
		self.ax_gibbs.grid(True, alpha=0.3)

		self.ax_enthalpy.set_xlabel('Composition', fontsize=12)
		self.ax_enthalpy.set_ylabel('Enthalpy of Mixing (kJ/mol)', fontsize=12)
		self.ax_enthalpy.set_title('混合焓对比', fontsize=14, fontweight='bold')
		self.ax_enthalpy.legend()
		self.ax_enthalpy.grid(True, alpha=0.3)

		self.ax_activity.set_xlabel('Composition', fontsize=12)
		self.ax_activity.set_ylabel('Activity', fontsize=12)
		self.ax_activity.set_title('活度对比', fontsize=14, fontweight='bold')
		self.ax_activity.legend()
		self.ax_activity.grid(True, alpha=0.3)

		self.canvas_gibbs.draw()
		self.canvas_enthalpy.draw()
		self.canvas_activity.draw()

	def clear_results(self):
		"""清除结果"""
		self.ax_gibbs.clear()
		self.ax_enthalpy.clear()
		self.ax_activity.clear()
		self.canvas_gibbs.draw()
		self.canvas_enthalpy.draw()
		self.canvas_activity.draw()
		self.log_text.delete('1.0', tk.END)
		self.data_text.delete('1.0', tk.END)
		self.results_data.clear()
		self.log("结果已清除")

	def export_data(self):
		"""导出热力学性质数据"""
		if not self.results_data:
			messagebox.showinfo("提示", "暂无数据可导出")
			return

		file_path = filedialog.asksaveasfilename(
			title="保存数据", defaultextension=".csv",
			filetypes=[("CSV", "*.csv"), ("All", "*.*")])
		if not file_path:
			return

		try:
			with open(file_path, 'w', encoding='utf-8') as f:
				f.write(f"# 热力学性质计算结果\n")
				f.write(f"# {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

				for key, data in self.results_data.items():
					if not key.endswith('_props'):
						continue
					label = data.get('label', key)
					x_vals = data.get('x_scan', [])
					gibbs = data.get('gibbs')
					enthalpy = data.get('enthalpy_mix')
					activity = data.get('activity', {})

					headers = ["X"]
					if gibbs is not None:
						headers.append("Gibbs(J/mol)")
					if enthalpy is not None:
						headers.append("Enthalpy_Mix(J/mol)")
					for comp in activity.keys():
						headers.append(f"Activity_{comp}")

					f.write(f"\n# {label}\n")
					f.write(",".join(headers) + "\n")

					for i, x in enumerate(x_vals):
						row = [f"{x:.6f}"]
						if gibbs is not None:
							row.append(f"{gibbs[i]:.4f}" if i < len(gibbs) else "NaN")
						if enthalpy is not None:
							row.append(f"{enthalpy[i]:.4f}" if i < len(enthalpy) else "NaN")
						for comp_acts in activity.values():
							val = comp_acts[i] if i < len(comp_acts) else np.nan
							row.append(f"{val:.6f}")
						f.write(",".join(row) + "\n")

			self.log(f"数据已导出到: {file_path}")
		except Exception as e:
			messagebox.showerror("错误", f"导出失败:\n{e}")


# =============================================================================
# 伪二元相图标签页
# =============================================================================
class PseudoBinaryTab(CalculatorTab):
	"""伪二元相图标签页 - 单模型"""

	def __init__(self, parent_notebook, app):
		super().__init__(parent_notebook, app, "伪二元相图")

	def _create_params_bar(self, parent):
		frame = ttk.LabelFrame(parent, text="伪二元相图参数", padding=5)
		frame.pack(fill=tk.X)

		row0 = ttk.Frame(frame)
		row0.pack(fill=tk.X, pady=2)

		ttk.Label(row0, text="基础合金:").pack(side=tk.LEFT, padx=3)
		self.base_entry = ttk.Entry(row0, width=12)
		self.base_entry.pack(side=tk.LEFT, padx=3)
		self.base_entry.insert(0, "AL")

		ttk.Label(row0, text="变化组分:").pack(side=tk.LEFT, padx=3)
		self.vary_combobox = ttk.Combobox(row0, width=8)
		self.vary_combobox.pack(side=tk.LEFT, padx=3)

		ttk.Label(row0, text="范围:").pack(side=tk.LEFT, padx=3)
		self.scan_start = ttk.Entry(row0, width=5)
		self.scan_start.pack(side=tk.LEFT, padx=1)
		self.scan_start.insert(0, "0.0")
		ttk.Label(row0, text="~").pack(side=tk.LEFT)
		self.scan_end = ttk.Entry(row0, width=5)
		self.scan_end.pack(side=tk.LEFT, padx=1)
		self.scan_end.insert(0, "1.0")
		ttk.Label(row0, text="共").pack(side=tk.LEFT, padx=2)
		self.scan_points = ttk.Entry(row0, width=4)
		self.scan_points.pack(side=tk.LEFT, padx=1)
		self.scan_points.insert(0, "51")
		ttk.Label(row0, text="点").pack(side=tk.LEFT)

		row1 = ttk.Frame(frame)
		row1.pack(fill=tk.X, pady=2)

		ttk.Label(row1, text="温度(K):").pack(side=tk.LEFT, padx=3)
		self.temp_min = ttk.Entry(row1, width=6)
		self.temp_min.pack(side=tk.LEFT, padx=1)
		self.temp_min.insert(0, "873")
		ttk.Label(row1, text="~").pack(side=tk.LEFT)
		self.temp_max = ttk.Entry(row1, width=6)
		self.temp_max.pack(side=tk.LEFT, padx=1)
		self.temp_max.insert(0, "2073")
		ttk.Label(row1, text="步长").pack(side=tk.LEFT, padx=3)
		self.temp_step = ttk.Entry(row1, width=4)
		self.temp_step.pack(side=tk.LEFT, padx=1)
		self.temp_step.insert(0, "10")
		ttk.Label(row1, text="K").pack(side=tk.LEFT)

		ttk.Separator(row1, orient='vertical').pack(side=tk.LEFT, fill=tk.Y, padx=8)

		ttk.Label(row1, text="模型:").pack(side=tk.LEFT, padx=3)
		self.model_var = tk.StringVar(value='RKM')
		for key in ('RKM', 'Muggianu', 'Toop', 'UEM1'):
			ttk.Radiobutton(row1, text=key, variable=self.model_var,
			                value=key).pack(side=tk.LEFT, padx=3)

		self.uem1_liquid_only = tk.BooleanVar(value=False)
		self.cb_uem_liq = ttk.Checkbutton(row1, text="UEM Liq Only",
		                                   variable=self.uem1_liquid_only,
		                                   state='disabled')
		self.cb_uem_liq.pack(side=tk.LEFT, padx=3)
		self.model_var.trace_add('write', self._on_model_change)

		ttk.Button(row1, text="开始计算",
		           command=self.start_calculation).pack(side=tk.RIGHT, padx=5)

	def _on_model_change(self, *args):
		if self.model_var.get() == 'UEM1':
			self.cb_uem_liq.config(state='normal')
		else:
			self.cb_uem_liq.config(state='disabled')
			self.uem1_liquid_only.set(False)

	def update_comps(self):
		self.vary_combobox['values'] = self.app.available_comps
		if len(self.app.available_comps) > 1:
			self.vary_combobox.current(1)

	def start_calculation(self):
		try:
			if not self.app.dbe:
				messagebox.showerror("错误", "请先加载数据库！")
				return

			base_alloy_str = self.base_entry.get().strip()
			varying_comp = self.vary_combobox.get().strip().upper()
			scan_start = float(self.scan_start.get())
			scan_end = float(self.scan_end.get())
			scan_points = int(self.scan_points.get())
			temp_min = float(self.temp_min.get())
			temp_max = float(self.temp_max.get())
			temp_step = float(self.temp_step.get())

			base_comps_dict = self._parse_base_alloy(base_alloy_str)

			if not varying_comp:
				messagebox.showerror("错误", "请选择变化组分！")
				return
			if varying_comp not in self.app.available_comps:
				messagebox.showerror("错误", f"变化组分 {varying_comp} 不在数据库中！")
				return

			all_comps = list(set(list(base_comps_dict.keys()) + [varying_comp]))
			phase_groups = classify_phases(self.app.dbe, all_comps)
			candidate_phases = sum(phase_groups.values(), [])

			if not candidate_phases:
				messagebox.showerror("错误", "没有找到适用于当前组分的相！")
				return

			self.log(f"智能筛选: 基于组分 {all_comps} 找到 {len(candidate_phases)} 个候选相")

			# AI 推荐相
			recommended = None
			if hasattr(self.app, 'llm_service') and self.app.llm_service.is_available():
				recommended = self.app.llm_service.predict_phases(
					all_comps, candidate_phases, 'pseudo_binary')

			dialog = PhaseSelectorDialog(self.app.root, candidate_phases,
			                             title="伪二元相图 - 相选择",
			                             recommended_phases=recommended,
			                             phase_groups=phase_groups)
			if dialog.result is None:
				self.log("用户取消了计算")
				return

			final_phases = dialog.result
			model_key = self.model_var.get()
			model_spec = self._get_model_spec(model_key, final_phases)
			T_points = int((temp_max - temp_min) / temp_step) + 1

			threading.Thread(
				target=self._run_thread,
				args=(base_comps_dict, varying_comp, (scan_start, scan_end),
				      scan_points, (temp_min, temp_max), T_points,
				      model_spec, model_key, final_phases),
				daemon=True
			).start()

		except ValueError as e:
			messagebox.showerror("参数错误", str(e))
		except Exception as e:
			messagebox.showerror("错误", f"启动计算失败:\n{e}")

	def _run_thread(self, base_comps, varying_comp, x_range, x_points,
	                T_range, T_points, model_spec, model_key, phases):
		try:
			self.log(f"开始计算伪二元相图: {base_comps} + {varying_comp}")
			self.log(f"模型: {model_key}, 成分范围: {x_range}, 点数: {x_points}")

			calculator = PseudoBinaryCalculator(
				database=self.app.dbe, logger=self.log)

			result = calculator.calculate(
				base_composition=base_comps,
				varying_comp=varying_comp,
				x_range=x_range,
				x_points=x_points,
				T_range=T_range,
				components=self.app.available_comps,
				phases=phases,
				model_spec=model_spec,
				pdens=self.app.pdens,
				T_points=T_points,
				verbose=True
			)

			if not isinstance(result, CalculationResult):
				self.log("✗ 返回类型不是 CalculationResult")
				return

			if not result.success:
				self.log(f"✗ 计算失败: {result.message}")
				return

			self.calc_result = result
			self.log(f"✓ 伪二元相图计算完成 (mode={result.mode})")

			if self.frame.winfo_exists():
				self.frame.after(0, lambda: self._plot_result(result, varying_comp, model_key))

		except Exception as e:
			self.log(f"✗ 计算线程异常: {e}")
			self.log(traceback.format_exc())

	def _plot_result(self, result, varying_comp, model_key):
		self._rebuild_chart()

		mode = result.mode
		if mode == CalculationMode.BINPLOT.value:
			self._copy_binplot_to_ax(result, self.ax, model_key, varying_comp)
		elif mode == CalculationMode.PHASE_MAP.value:
			self._plot_phase_map_to_ax(result, self.ax, model_key, varying_comp)
		else:
			self._plot_phase_map_to_ax(result, self.ax, model_key, varying_comp)

		self.ax.grid(True, alpha=0.3)
		self.canvas.draw()

	def export_data(self):
		if self.calc_result is None:
			messagebox.showwarning("警告", "没有可导出的数据！")
			return

		file_path = filedialog.asksaveasfilename(
			title="保存数据", defaultextension=".csv",
			filetypes=[("CSV", "*.csv"), ("All", "*.*")])
		if not file_path:
			return

		try:
			result = self.calc_result
			with open(file_path, 'w', encoding='utf-8') as f:
				f.write(f"# 伪二元相图数据\n")
				f.write(f"# {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

				if result.mode == CalculationMode.BINPLOT.value:
					f.write("Line_ID,Label,X,Temperature_K\n")
					if result.figure and result.figure.axes:
						ax = result.figure.axes[0]
						for i, line in enumerate(ax.get_lines()):
							ll = line.get_label()
							if not ll or ll.startswith('_'):
								ll = f"Boundary_{i}"
							for x, y in zip(line.get_xdata(), line.get_ydata()):
								f.write(f"{i},{ll},{x:.6f},{y:.6f}\n")

				elif result.mode == CalculationMode.PHASE_MAP.value:
					f.write("X,Temperature_K,Phase\n")
					if result.phase_map is not None:
						x_grid = result.x_axis
						t_grid = result.t_grid
						phase_legend = getattr(result, 'phase_legend', None)
						for r in range(len(t_grid)):
							for c in range(len(x_grid)):
								val = result.phase_map[r, c]
								if phase_legend is not None:
									phase_str = phase_legend.get(int(val), "Error").replace(',', ';')
								else:
									phase_str = str(val).replace(',', ';')
								f.write(f"{x_grid[c]:.6f},{t_grid[r]:.6f},{phase_str}\n")

			self.log(f"数据已导出到: {file_path}")
		except Exception as e:
			messagebox.showerror("错误", f"导出失败:\n{e}")


# =============================================================================
# 三元等温相图标签页
# =============================================================================
class TernaryTab(CalculatorTab):
	"""三元等温相图标签页"""

	def __init__(self, parent_notebook, app):
		super().__init__(parent_notebook, app, "三元相图")

	def _create_params_bar(self, parent):
		frame = ttk.LabelFrame(parent, text="三元等温相图参数", padding=5)
		frame.pack(fill=tk.X)

		row0 = ttk.Frame(frame)
		row0.pack(fill=tk.X, pady=2)

		ttk.Label(row0, text="组分(3个,逗号分隔):").pack(side=tk.LEFT, padx=3)
		self.comps_entry = ttk.Entry(row0, width=15)
		self.comps_entry.pack(side=tk.LEFT, padx=3)
		self.comps_entry.insert(0, "AL,CR,NI")

		ttk.Label(row0, text="温度(K):").pack(side=tk.LEFT, padx=3)
		self.temp_entry = ttk.Entry(row0, width=8)
		self.temp_entry.pack(side=tk.LEFT, padx=3)
		self.temp_entry.insert(0, "1273")

		ttk.Label(row0, text="网格点数:").pack(side=tk.LEFT, padx=3)
		self.points_entry = ttk.Entry(row0, width=5)
		self.points_entry.pack(side=tk.LEFT, padx=3)
		self.points_entry.insert(0, "51")

		# 第二行：模型 + 开始计算
		row1 = ttk.Frame(frame)
		row1.pack(fill=tk.X, pady=2)

		ttk.Label(row1, text="模型:").pack(side=tk.LEFT, padx=3)
		self.model_var = tk.StringVar(value='RKM')
		for key in ('RKM', 'Muggianu', 'Toop', 'UEM1'):
			ttk.Radiobutton(row1, text=key, variable=self.model_var,
			                value=key).pack(side=tk.LEFT, padx=3)

		self.uem1_liquid_only = tk.BooleanVar(value=False)
		self.cb_uem_liq = ttk.Checkbutton(row1, text="UEM Liq Only",
		                                   variable=self.uem1_liquid_only,
		                                   state='disabled')
		self.cb_uem_liq.pack(side=tk.LEFT, padx=3)
		self.model_var.trace_add('write', self._on_model_change)

		ttk.Button(row1, text="开始计算",
		           command=self.start_calculation).pack(side=tk.RIGHT, padx=5)

	def _on_model_change(self, *args):
		if self.model_var.get() == 'UEM1':
			self.cb_uem_liq.config(state='normal')
		else:
			self.cb_uem_liq.config(state='disabled')
			self.uem1_liquid_only.set(False)

	def start_calculation(self):
		try:
			if not self.app.dbe:
				messagebox.showerror("错误", "请先加载数据库！")
				return

			comps_str = self.comps_entry.get().strip()
			temperature = float(self.temp_entry.get())
			grid_points = int(self.points_entry.get())

			comps = [c.strip().upper() for c in comps_str.split(',')]
			if len(comps) != 3:
				messagebox.showerror("错误", "必须指定3个组分！")
				return

			for comp in comps:
				if comp not in self.app.available_comps:
					messagebox.showerror("错误", f"组分 {comp} 不在数据库中！")
					return

			phase_groups = classify_phases(self.app.dbe, comps)
			candidate_phases = sum(phase_groups.values(), [])
			if not candidate_phases:
				messagebox.showerror("错误", "没有找到适用于当前组分的相！")
				return

			# AI 推荐相
			recommended = None
			if hasattr(self.app, 'llm_service') and self.app.llm_service.is_available():
				recommended = self.app.llm_service.predict_phases(
					comps, candidate_phases, 'ternary')

			dialog = PhaseSelectorDialog(self.app.root, candidate_phases,
			                             title="三元等温相图 - 相选择",
			                             recommended_phases=recommended,
			                             phase_groups=phase_groups)
			if dialog.result is None:
				self.log("用户取消了计算")
				return

			final_phases = dialog.result
			model_key = self.model_var.get()
			model_spec = self._get_model_spec(model_key, final_phases)

			threading.Thread(
				target=self._run_thread,
				args=(comps, temperature, grid_points, final_phases, model_spec, model_key),
				daemon=True
			).start()

		except ValueError as e:
			messagebox.showerror("参数错误", str(e))
		except Exception as e:
			messagebox.showerror("错误", f"启动计算失败:\n{e}")

	def _run_thread(self, comps, temperature, grid_points, phases, model_spec, model_key):
		try:
			self.log(f"开始计算三元等温相图: {'-'.join(comps)}")
			self.log(f"模型: {model_key}, 温度: {temperature} K, 网格: {grid_points}")

			calculator = TernaryCalculator(
				database=self.app.dbe, logger=self.log)

			result = calculator.calculate(
				components=comps + ['VA'],
				temperature=temperature,
				phases=phases,
				model_spec=model_spec,
				pdens=self.app.pdens,
				grid_points=grid_points,
				verbose=True
			)

			if not isinstance(result, CalculationResult) or not result.success:
				self.log(f"✗ 计算失败: {getattr(result, 'message', '未知错误')}")
				return

			self.calc_result = result
			self.log(f"✓ 三元相图计算完成 (mode={result.mode})")

			if self.frame.winfo_exists():
				self.frame.after(0, lambda: self._plot_result(result, comps, temperature, model_key))

		except Exception as e:
			self.log(f"✗ 计算线程异常: {e}")
			self.log(traceback.format_exc())

	def _plot_result(self, result, comps, temperature, model_key):
		source_fig = result.figure
		if source_fig is None or len(source_fig.axes) == 0:
			self.log("警告: ternary figure 为空")
			return

		try:
			# 三元相图使用三角投影，不能简单复制到普通 axes
			# 直接将 ternplot 生成的 figure 嵌入 canvas
			parent = self.canvas.get_tk_widget().master
			self.toolbar.destroy()
			self.canvas.get_tk_widget().destroy()

			# 设置标题
			source_ax = source_fig.axes[0]
			source_ax.set_title(
				f'三元等温相图: {"-".join(comps)} @ {temperature}K ({model_key})',
				fontsize=12, fontweight='bold')

			self.fig = source_fig
			self.ax = source_ax

			self.canvas = FigureCanvasTkAgg(self.fig, parent)
			self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
			self.toolbar = NavigationToolbar2Tk(self.canvas, parent)
			self.toolbar.update()
			self.canvas.draw()
			self.log("✓ ternary 图形绘制完成")

		except Exception as e:
			self.log(f"✗ ternary 绘制失败: {e}")
			self.log(traceback.format_exc())

	def export_data(self):
		if self.calc_result is None:
			messagebox.showwarning("警告", "没有可导出的数据！")
			return

		file_path = filedialog.asksaveasfilename(
			title="保存数据", defaultextension=".csv",
			filetypes=[("CSV", "*.csv"), ("All", "*.*")])
		if not file_path:
			return

		try:
			result = self.calc_result
			with open(file_path, 'w', encoding='utf-8') as f:
				f.write(f"# 三元相图数据\n")
				f.write(f"# {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
				f.write("Type,ID,X,Y\n")
				if result.figure and result.figure.axes:
					ax = result.figure.axes[0]
					for i, line in enumerate(ax.get_lines()):
						for x, y in zip(line.get_xdata(), line.get_ydata()):
							f.write(f"Line,{i+1},{x:.6f},{y:.6f}\n")
					for i, coll in enumerate(ax.collections):
						offsets = coll.get_offsets()
						for off in offsets:
							f.write(f"Scatter,{i+1},{off[0]:.6f},{off[1]:.6f}\n")
			self.log(f"数据已导出到: {file_path}")
		except Exception as e:
			messagebox.showerror("错误", f"导出失败:\n{e}")


# =============================================================================
# 三元液相面标签页
# =============================================================================
class LiquidusSurfaceTab(CalculatorTab):
	"""三元液相面标签页"""

	def __init__(self, parent_notebook, app):
		super().__init__(parent_notebook, app, "液相面")

	def _create_params_bar(self, parent):
		frame = ttk.LabelFrame(parent, text="三元液相面参数", padding=5)
		frame.pack(fill=tk.X)

		row0 = ttk.Frame(frame)
		row0.pack(fill=tk.X, pady=2)

		ttk.Label(row0, text="组分(3个):").pack(side=tk.LEFT, padx=3)
		self.comps_entry = ttk.Entry(row0, width=15)
		self.comps_entry.pack(side=tk.LEFT, padx=3)
		self.comps_entry.insert(0, "AL,CR,NI")

		ttk.Label(row0, text="网格点数:").pack(side=tk.LEFT, padx=3)
		self.points_entry = ttk.Entry(row0, width=5)
		self.points_entry.pack(side=tk.LEFT, padx=3)
		self.points_entry.insert(0, "31")

		ttk.Label(row0, text="温度(K):").pack(side=tk.LEFT, padx=3)
		self.temp_min = ttk.Entry(row0, width=6)
		self.temp_min.pack(side=tk.LEFT, padx=1)
		self.temp_min.insert(0, "873")
		ttk.Label(row0, text="~").pack(side=tk.LEFT)
		self.temp_max = ttk.Entry(row0, width=6)
		self.temp_max.pack(side=tk.LEFT, padx=1)
		self.temp_max.insert(0, "2073")

		ttk.Label(row0, text="绘图:").pack(side=tk.LEFT, padx=3)
		self.plot_type = tk.StringVar(value='contour')
		ttk.Radiobutton(row0, text="3D", variable=self.plot_type, value='3d').pack(side=tk.LEFT, padx=2)
		ttk.Radiobutton(row0, text="等高线", variable=self.plot_type, value='contour').pack(side=tk.LEFT, padx=2)

		# 第二行：模型 + 开始计算
		row1 = ttk.Frame(frame)
		row1.pack(fill=tk.X, pady=2)

		ttk.Label(row1, text="模型:").pack(side=tk.LEFT, padx=3)
		self.model_var = tk.StringVar(value='RKM')
		for key in ('RKM', 'Muggianu', 'Toop', 'UEM1'):
			ttk.Radiobutton(row1, text=key, variable=self.model_var,
			                value=key).pack(side=tk.LEFT, padx=3)

		self.uem1_liquid_only = tk.BooleanVar(value=False)
		self.cb_uem_liq = ttk.Checkbutton(row1, text="UEM Liq Only",
		                                   variable=self.uem1_liquid_only,
		                                   state='disabled')
		self.cb_uem_liq.pack(side=tk.LEFT, padx=3)
		self.model_var.trace_add('write', self._on_model_change)

		# AI 代理加速选项
		self.ai_accel_var = tk.BooleanVar(value=True)
		self.cb_ai_accel = ttk.Checkbutton(row1, text="AI 加速",
		                                    variable=self.ai_accel_var)
		if AI_SURROGATE_AVAILABLE:
			ttk.Separator(row1, orient='vertical').pack(side=tk.LEFT, fill='y', padx=3, pady=2)
			self.cb_ai_accel.pack(side=tk.LEFT, padx=3)

		ttk.Button(row1, text="开始计算",
		           command=self.start_calculation).pack(side=tk.RIGHT, padx=5)

	def _on_model_change(self, *args):
		if self.model_var.get() == 'UEM1':
			self.cb_uem_liq.config(state='normal')
		else:
			self.cb_uem_liq.config(state='disabled')
			self.uem1_liquid_only.set(False)

	def start_calculation(self):
		try:
			if not self.app.dbe:
				messagebox.showerror("错误", "请先加载数据库！")
				return

			comps_str = self.comps_entry.get().strip()
			grid_points = int(self.points_entry.get())
			temp_min = float(self.temp_min.get())
			temp_max = float(self.temp_max.get())

			comps = [c.strip().upper() for c in comps_str.split(',')]
			if len(comps) != 3:
				messagebox.showerror("错误", "必须指定3个组分！")
				return

			for comp in comps:
				if comp not in self.app.available_comps:
					messagebox.showerror("错误", f"组分 {comp} 不在数据库中！")
					return

			phase_groups = classify_phases(self.app.dbe, comps)
			candidate_phases = sum(phase_groups.values(), [])
			if not candidate_phases:
				messagebox.showerror("错误", "没有找到适用于当前组分的相！")
				return

			# AI 推荐相
			recommended = None
			if hasattr(self.app, 'llm_service') and self.app.llm_service.is_available():
				recommended = self.app.llm_service.predict_phases(
					comps, candidate_phases, 'liquidus_surface')

			dialog = PhaseSelectorDialog(self.app.root, candidate_phases,
			                             title="三元液相面 - 相选择",
			                             recommended_phases=recommended,
			                             phase_groups=phase_groups)
			if dialog.result is None:
				self.log("用户取消了计算")
				return

			final_phases = dialog.result
			model_key = self.model_var.get()
			model_spec = self._get_model_spec(model_key, final_phases)
			plot_type = self.plot_type.get()
			use_surrogate = (AI_SURROGATE_AVAILABLE and
			                 self.ai_accel_var.get())

			threading.Thread(
				target=self._run_thread,
				args=(comps, grid_points, (temp_min, temp_max), final_phases,
				      model_spec, model_key, plot_type, use_surrogate),
				daemon=True
			).start()

		except ValueError as e:
			messagebox.showerror("参数错误", str(e))
		except Exception as e:
			messagebox.showerror("错误", f"启动计算失败:\n{e}")

	def _run_thread(self, comps, grid_points, T_range, phases,
	                model_spec, model_key, plot_type, use_surrogate=False):
		try:
			self.log(f"开始计算三元液相面: {'-'.join(comps)}")
			self.log(f"模型: {model_key}, 网格: {grid_points}")

			if use_surrogate:
				# ---- AI 代理加速模式 ----
				self.log("🚀 AI 代理加速模式已启用")

				# 获取或创建共享缓存（挂在 app 上跨标签页复用）
				if not hasattr(self.app, '_surrogate_cache'):
					self.app._surrogate_cache = SurrogateCache()

				calculator = AdaptiveSurfaceCalculator(
					database=self.app.dbe,
					logger=self.log,
					cache=self.app._surrogate_cache,
				)

				# 进度回调：在主线程更新图表
				def on_progress(stage, result, info):
					anchor_n = info.get('anchor_count', 0)
					total_n = info.get('total_target', 0)
					elapsed = info.get('elapsed_seconds', 0)

					if stage == 'preview' and result is not None:
						self.calc_result = result
						if self.frame.winfo_exists():
							self.frame.after(0, lambda r=result: self._plot_result(
								r, comps, model_key + " [预览]", plot_type))
					elif stage == 'refined' and result is not None:
						speedup = info.get('estimated_speedup', 1.0)
						self.calc_result = result
						if self.frame.winfo_exists():
							self.frame.after(0, lambda r=result: self._plot_result(
								r, comps, model_key, plot_type))

				result = calculator.calculate(
					comp_b_range=(0, 1),
					comp_c_range=(0, 1),
					grid_points=grid_points,
					T_range=T_range,
					components=comps + ['VA'],
					phases=phases,
					model_spec=model_spec,
					pdens=self.app.pdens,
					verbose=True,
					progress_callback=on_progress,
				)
			else:
				# ---- 原始全量计算模式 ----
				calculator = LiquidusSurfaceCalculator(
					database=self.app.dbe, logger=self.log)

				result = calculator.calculate(
					comp_b_range=(0, 1),
					comp_c_range=(0, 1),
					grid_points=grid_points,
					T_range=T_range,
					components=comps + ['VA'],
					phases=phases,
					model_spec=model_spec,
					pdens=self.app.pdens,
					verbose=True
				)

			if not isinstance(result, CalculationResult) or not result.success:
				self.log(f"✗ 计算失败: {getattr(result, 'message', '未知错误')}")
				return

			self.calc_result = result
			self.log(f"✓ 三元液相面计算完成")

			# 非加速模式需要手动绘图；加速模式已在回调中绘制
			if not use_surrogate:
				if self.frame.winfo_exists():
					self.frame.after(0, lambda: self._plot_result(result, comps, model_key, plot_type))

		except Exception as e:
			self.log(f"✗ 计算线程异常: {e}")
			self.log(traceback.format_exc())

	# ---- 三角形坐标变换 ----
	@staticmethod
	def _ternary_xy(xb, xc):
		"""将 (xB, xC) 变换为等边三角形笛卡尔坐标"""
		x = xb + xc * 0.5
		y = xc * (np.sqrt(3) / 2)
		return x, y

	def _draw_ternary_frame(self, ax, comp_a, comp_b, comp_c):
		"""绘制等边三角形边框、刻度、组分标签"""
		H = np.sqrt(3) / 2
		# 边框
		ax.plot([0, 1], [0, 0], 'k-', linewidth=2)
		ax.plot([1, 0.5], [0, H], 'k-', linewidth=2)
		ax.plot([0.5, 0], [H, 0], 'k-', linewidth=2)

		# 顶点标签
		ax.text(-0.03, -0.04, comp_b, ha='center', va='top', fontsize=12, fontweight='bold')
		ax.text(1.03, -0.04, comp_c, ha='center', va='top', fontsize=12, fontweight='bold')
		ax.text(0.5, H + 0.04, comp_a, ha='center', va='bottom', fontsize=12, fontweight='bold')

		# 网格线与刻度
		for i in range(1, 10):
			f = i / 10
			# 平行于 B-C 底边（等 xA 线）
			x0, y0 = self._ternary_xy(0, f)
			x1, y1 = self._ternary_xy(1 - f, f)
			ax.plot([x0, x1], [y0, y1], color='gray', linewidth=0.4, alpha=0.4)
			# 平行于 A-C 右边（等 xB 线）
			x0, y0 = self._ternary_xy(f, 0)
			x1, y1 = self._ternary_xy(f, 1 - f)
			ax.plot([x0, x1], [y0, y1], color='gray', linewidth=0.4, alpha=0.4)
			# 平行于 A-B 左边（等 xC 线）
			x0, y0 = self._ternary_xy(0, f)
			x1, y1 = self._ternary_xy(1 - f, 0)
			ax.plot([x0, x1], [y0, y1], color='gray', linewidth=0.4, alpha=0.4)

		# 底边刻度标签 (xB)
		for i in range(0, 11, 2):
			f = i / 10
			ax.text(f, -0.03, f'{f:.1f}', ha='center', va='top', fontsize=7, color='gray')
		# 左边刻度标签 (xA)
		for i in range(0, 11, 2):
			f = i / 10
			x, y = self._ternary_xy(0, f)
			ax.text(x - 0.03, y, f'{f:.1f}', ha='right', va='center', fontsize=7,
			        color='gray', rotation=60)
		# 右边刻度标签 (xC)
		for i in range(0, 11, 2):
			f = i / 10
			x, y = self._ternary_xy(1 - f, f)
			ax.text(x + 0.03, y, f'{f:.1f}', ha='left', va='center', fontsize=7,
			        color='gray', rotation=-60)

		ax.set_aspect('equal')
		ax.set_xlim(-0.12, 1.12)
		ax.set_ylim(-0.10, H + 0.10)
		ax.axis('off')

	def _plot_result(self, result, comps, model_key, plot_type):
		surface_data = result.surface_data
		if surface_data is not None:
			xb = surface_data.get('xb_data', result.x_axis)
			xc = surface_data.get('xc_data', result.y_axis)
			t_data = surface_data.get('t_data', result.z_axis)
			comp_a = surface_data.get('comp_a', comps[0])
			comp_b = surface_data.get('comp_b', comps[1])
			comp_c = surface_data.get('comp_c', comps[2])
		else:
			xb = result.x_axis
			xc = result.y_axis
			t_data = result.z_axis
			comp_a = comps[0]
			comp_b = comps[1]
			comp_c = comps[2]

		if xb is None or xc is None or t_data is None:
			self.log("警告: surface 数据不完整")
			return

		valid_mask = ~np.isnan(t_data)
		xb_valid = xb[valid_mask]
		xc_valid = xc[valid_mask]
		t_valid = t_data[valid_mask]

		if len(t_valid) == 0:
			self.log("警告: 没有有效的液相面数据")
			return

		try:
			if plot_type == '3d':
				self._rebuild_chart(projection='3d')
				# 3D 也用三角形坐标 (x, y) + 温度 z
				x3d, y3d = self._ternary_xy(xb_valid, xc_valid)
				scatter = self.ax.scatter(x3d, y3d, t_valid,
				                          c=t_valid, cmap='coolwarm', s=15, alpha=0.8)
				self.ax.set_xlabel(f'X({comp_b})', fontsize=10)
				self.ax.set_ylabel(f'X({comp_c})', fontsize=10)
				self.ax.set_zlabel('液相线温度 (K)', fontsize=10)
				self.ax.set_title(f'三元液相面: {"-".join(comps)} ({model_key})',
				                  fontsize=12, fontweight='bold')
				self.fig.colorbar(scatter, ax=self.ax, label='温度 (K)', shrink=0.6)
			else:
				self._rebuild_chart()
				from scipy.interpolate import griddata

				# 在原始 (xb, xc) 空间插值，然后变换到三角坐标绘图
				grid_n = 80
				grid_b = np.linspace(0, 1, grid_n)
				grid_c = np.linspace(0, 1, grid_n)
				Xb, Xc = np.meshgrid(grid_b, grid_c)
				ternary_mask = (Xb + Xc <= 1.0 + 1e-9)

				T_grid = griddata((xb_valid, xc_valid), t_valid,
				                  (Xb, Xc), method='linear')
				T_grid[~ternary_mask] = np.nan

				# 变换到等边三角形坐标
				X_cart, Y_cart = self._ternary_xy(Xb, Xc)

				levels = 20
				contour = self.ax.contourf(X_cart, Y_cart, T_grid, levels=levels,
				                           cmap='coolwarm', alpha=0.85)
				self.ax.contour(X_cart, Y_cart, T_grid, levels=levels, colors='k',
				                linewidths=0.4, alpha=0.4)
				self.fig.colorbar(contour, ax=self.ax, label='液相线温度 (K)')

				# 绘制等边三角形边框 + 刻度
				self._draw_ternary_frame(self.ax, comp_a, comp_b, comp_c)
				self.ax.set_title(f'三元液相面: {"-".join(comps)} ({model_key})',
				                  fontsize=12, fontweight='bold')

			self.canvas.draw()
			self.log(f"✓ 液相面绘制完成 ({len(t_valid)} 个有效点)")

		except Exception as e:
			self.log(f"✗ 液相面绘制失败: {e}")
			self.log(traceback.format_exc())

	def export_data(self):
		if self.calc_result is None:
			messagebox.showwarning("警告", "没有可导出的数据！")
			return

		file_path = filedialog.asksaveasfilename(
			title="保存数据", defaultextension=".csv",
			filetypes=[("CSV", "*.csv"), ("All", "*.*")])
		if not file_path:
			return

		try:
			result = self.calc_result
			with open(file_path, 'w', encoding='utf-8') as f:
				f.write(f"# 三元液相面数据\n")
				f.write(f"# {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
				f.write("X_B,X_C,Temperature_K\n")
				sd = result.surface_data
				if sd:
					xb = sd.get('xb_data')
					xc = sd.get('xc_data')
					t_data = sd.get('t_data')
					if xb is not None:
						valid_mask = ~np.isnan(t_data)
						for i in np.where(valid_mask)[0]:
							f.write(f"{xb[i]:.6f},{xc[i]:.6f},{t_data[i]:.6f}\n")
			self.log(f"数据已导出到: {file_path}")
		except Exception as e:
			messagebox.showerror("错误", f"导出失败:\n{e}")


# =============================================================================
# 溶解度曲线标签页
# =============================================================================
class SolubilityTab(CalculatorTab):
	"""溶解度曲线标签页"""

	def __init__(self, parent_notebook, app):
		super().__init__(parent_notebook, app, "溶解度")

	def _create_params_bar(self, parent):
		frame = ttk.LabelFrame(parent, text="溶解度参数", padding=5)
		frame.pack(fill=tk.X)

		row0 = ttk.Frame(frame)
		row0.pack(fill=tk.X, pady=2)

		ttk.Label(row0, text="基础合金:").pack(side=tk.LEFT, padx=3)
		self.base_entry = ttk.Entry(row0, width=12)
		self.base_entry.pack(side=tk.LEFT, padx=3)
		self.base_entry.insert(0, "AL")

		ttk.Label(row0, text="溶解元素:").pack(side=tk.LEFT, padx=3)
		self.solute_combobox = ttk.Combobox(row0, width=8)
		self.solute_combobox.pack(side=tk.LEFT, padx=3)

		ttk.Label(row0, text="温度(K):").pack(side=tk.LEFT, padx=3)
		self.temp_start = ttk.Entry(row0, width=6)
		self.temp_start.pack(side=tk.LEFT, padx=1)
		self.temp_start.insert(0, "673")
		ttk.Label(row0, text="~").pack(side=tk.LEFT)
		self.temp_end = ttk.Entry(row0, width=6)
		self.temp_end.pack(side=tk.LEFT, padx=1)
		self.temp_end.insert(0, "2073")
		ttk.Label(row0, text="共").pack(side=tk.LEFT, padx=2)
		self.temp_points = ttk.Entry(row0, width=4)
		self.temp_points.pack(side=tk.LEFT, padx=1)
		self.temp_points.insert(0, "20")
		ttk.Label(row0, text="点").pack(side=tk.LEFT)

		# 第二行：模型 + 开始计算
		row1 = ttk.Frame(frame)
		row1.pack(fill=tk.X, pady=2)

		ttk.Label(row1, text="模型:").pack(side=tk.LEFT, padx=3)
		self.model_var = tk.StringVar(value='RKM')
		for key in ('RKM', 'Muggianu', 'Toop', 'UEM1'):
			ttk.Radiobutton(row1, text=key, variable=self.model_var,
			                value=key).pack(side=tk.LEFT, padx=3)

		self.uem1_liquid_only = tk.BooleanVar(value=False)
		self.cb_uem_liq = ttk.Checkbutton(row1, text="UEM Liq Only",
		                                   variable=self.uem1_liquid_only,
		                                   state='disabled')
		self.cb_uem_liq.pack(side=tk.LEFT, padx=3)
		self.model_var.trace_add('write', self._on_model_change)

		ttk.Button(row1, text="开始计算",
		           command=self.start_calculation).pack(side=tk.RIGHT, padx=5)

	def _on_model_change(self, *args):
		if self.model_var.get() == 'UEM1':
			self.cb_uem_liq.config(state='normal')
		else:
			self.cb_uem_liq.config(state='disabled')
			self.uem1_liquid_only.set(False)

	def update_comps(self):
		self.solute_combobox['values'] = self.app.available_comps
		if len(self.app.available_comps) > 2:
			self.solute_combobox.current(2)

	def start_calculation(self):
		try:
			if not self.app.dbe:
				messagebox.showerror("错误", "请先加载数据库！")
				return

			base_alloy_str = self.base_entry.get().strip()
			solute_comp = self.solute_combobox.get().strip().upper()
			temp_start = float(self.temp_start.get())
			temp_end = float(self.temp_end.get())
			temp_points = int(self.temp_points.get())

			base_comps_dict = self._parse_base_alloy(base_alloy_str)

			if not solute_comp:
				messagebox.showerror("错误", "请选择溶解元素！")
				return
			if solute_comp not in self.app.available_comps:
				messagebox.showerror("错误", f"溶解元素 {solute_comp} 不在数据库中！")
				return

			all_comps = sorted(list(set(list(base_comps_dict.keys()) + [solute_comp])))
			phase_groups = classify_phases(self.app.dbe, all_comps)
			candidate_phases = sum(phase_groups.values(), [])

			if not candidate_phases:
				messagebox.showerror("错误", "没有找到适用于当前组分的相！")
				return

			# AI 推荐相
			recommended = None
			if hasattr(self.app, 'llm_service') and self.app.llm_service.is_available():
				recommended = self.app.llm_service.predict_phases(
					all_comps, candidate_phases, 'solubility')

			dialog = PhaseSelectorDialog(self.app.root, candidate_phases,
			                             title="溶解度计算 - 相选择",
			                             recommended_phases=recommended,
			                             phase_groups=phase_groups)
			if dialog.result is None:
				self.log("用户取消了计算")
				return

			model_key = self.model_var.get()
			final_phases = dialog.result
			model_spec = self._get_model_spec(model_key, final_phases)

			threading.Thread(
				target=self._run_thread,
				args=(base_comps_dict, solute_comp, (temp_start, temp_end),
				      temp_points, model_spec, model_key, final_phases),
				daemon=True
			).start()

		except ValueError as e:
			messagebox.showerror("参数错误", str(e))
		except Exception as e:
			messagebox.showerror("错误", f"启动计算失败:\n{e}")

	def _run_thread(self, base_comps, solute, T_range, T_points,
	                model_spec, model_key, phases):
		try:
			self.log(f"开始计算溶解度曲线: 基础={base_comps}, 溶质={solute}")

			calculator = SolubilityCalculator(
				database=self.app.dbe,
				model_spec=model_spec,
				pdens=self.app.pdens,
				logger=self.log
			)

			result = calculator.calculate_solvus(
				base_composition=base_comps,
				solute=solute,
				T_range=T_range,
				T_points=T_points,
				verbose=True,
				phases=phases
			)

			if not isinstance(result, CalculationResult) or not result.success:
				self.log(f"✗ 计算失败: {getattr(result, 'message', '未知错误')}")
				return

			self.calc_result = result
			self.log(f"✓ 溶解度曲线计算完成")

			if self.frame.winfo_exists():
				self.frame.after(0, lambda: self._plot_result(result, base_comps, solute, model_key))

		except Exception as e:
			self.log(f"✗ 计算线程异常: {e}")
			self.log(traceback.format_exc())

	def _plot_result(self, result, base_comps, solute, model_key):
		self._rebuild_chart()

		temperatures = result.x_axis
		solubilities = result.y_axis

		if temperatures is None or solubilities is None:
			sol_data = result.solubility_data
			if sol_data is not None:
				temperatures = sol_data.get('temperatures')
				solubilities = sol_data.get('solubilities')

		if temperatures is None or solubilities is None:
			self.log("警告: 溶解度数据不完整")
			return

		try:
			temperatures = np.asarray(temperatures)
			solubilities = np.asarray(solubilities)

			if np.nanmax(solubilities) <= 1.0:
				sol_percent = solubilities * 100
				ylabel = f'{solute} 溶解度 (at%)'
			else:
				sol_percent = solubilities
				ylabel = f'{solute} 溶解度'

			valid_mask = ~np.isnan(sol_percent)
			if not np.any(valid_mask):
				self.log("警告: 没有有效的溶解度数据")
				return

			self.ax.plot(temperatures[valid_mask], sol_percent[valid_mask],
			             'b-o', linewidth=2, markersize=6, label=f'{solute} 溶解度')

			base_str = '+'.join([f'{k}{v:.1f}' if v != 1.0 else k
			                     for k, v in base_comps.items()])
			self.ax.set_xlabel('温度 (K)', fontsize=12)
			self.ax.set_ylabel(ylabel, fontsize=12)
			self.ax.set_title(f'溶解度曲线: {solute} in {base_str} ({model_key})',
			                  fontsize=12, fontweight='bold')
			self.ax.legend(loc='best')
			self.ax.grid(True, alpha=0.3)

			valid_sol = sol_percent[valid_mask]
			valid_temp = temperatures[valid_mask]
			info_text = (f'有效点数: {len(valid_sol)}\n'
			             f'温度: {np.min(valid_temp):.0f}-{np.max(valid_temp):.0f} K\n'
			             f'溶解度: {np.min(valid_sol):.4f}-{np.max(valid_sol):.4f}')
			self.ax.text(0.02, 0.98, info_text,
			             transform=self.ax.transAxes, fontsize=9,
			             verticalalignment='top',
			             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

			self.canvas.draw()
			self.log(f"✓ 溶解度曲线绘制完成 ({np.sum(valid_mask)} 个有效点)")

		except Exception as e:
			self.log(f"✗ 溶解度曲线绘制失败: {e}")
			self.log(traceback.format_exc())

	def export_data(self):
		if self.calc_result is None:
			messagebox.showwarning("警告", "没有可导出的数据！")
			return

		file_path = filedialog.asksaveasfilename(
			title="保存数据", defaultextension=".csv",
			filetypes=[("CSV", "*.csv"), ("All", "*.*")])
		if not file_path:
			return

		try:
			result = self.calc_result
			with open(file_path, 'w', encoding='utf-8') as f:
				f.write(f"# 溶解度曲线数据\n")
				f.write(f"# {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
				f.write("Temperature_K,Solubility\n")

				if result.solubility_data and 'temperatures' in result.solubility_data:
					temps = result.solubility_data['temperatures']
					sols = result.solubility_data['solubilities']
				else:
					temps = result.x_axis
					sols = result.y_axis

				if temps is not None and sols is not None:
					for t, s in zip(temps, sols):
						f.write(f"{t:.6f},{s:.6f}\n")

			self.log(f"数据已导出到: {file_path}")
		except Exception as e:
			messagebox.showerror("错误", f"导出失败:\n{e}")


# =============================================================================
# LLM 助手面板（嵌入侧边栏）
# =============================================================================
class LLMAssistantPanel:
	"""LLM 智能助手面板 - 嵌入侧边栏底部"""

	def __init__(self, parent, app):
		self.app = app
		self.frame = ttk.LabelFrame(parent, text="AI 助手", padding=5)
		self.frame.pack(fill=tk.BOTH, expand=True, pady=5, padx=5)

		self._create_widgets()

	def _create_widgets(self):
		# 状态指示
		status_frame = ttk.Frame(self.frame)
		status_frame.pack(fill=tk.X, pady=(0, 5))

		self.status_label = ttk.Label(status_frame, text="● 未配置",
		                              foreground="gray", font=('', 9))
		self.status_label.pack(side=tk.LEFT)

		ttk.Button(status_frame, text="设置", width=6,
		           command=self._open_settings).pack(side=tk.RIGHT)
		ttk.Button(status_frame, text="测试", width=6,
		           command=self._test_connection).pack(side=tk.RIGHT, padx=2)

		# 对话显示
		self.chat_text = scrolledtext.ScrolledText(self.frame, wrap=tk.WORD,
		                                           font=('', 9), height=10,
		                                           state=tk.DISABLED)
		self.chat_text.pack(fill=tk.BOTH, expand=True, pady=2)

		# 输入区
		input_frame = ttk.Frame(self.frame)
		input_frame.pack(fill=tk.X, pady=(2, 0))

		self.input_entry = ttk.Entry(input_frame, font=('', 10))
		self.input_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 3))
		self.input_entry.bind('<Return>', lambda e: self._send_message())

		ttk.Button(input_frame, text="发送", width=6,
		           command=self._send_message).pack(side=tk.RIGHT)

		# 快捷按钮
		btn_frame = ttk.Frame(self.frame)
		btn_frame.pack(fill=tk.X, pady=(3, 0))

		buttons = [
			("解析合金", self._action_parse_comp),
			("推荐相", self._action_predict_phases),
			("推荐参数", self._action_suggest_params),
			("推荐模型", self._action_recommend_model),
			("解读结果", self._action_interpret_result),
		]

		for text, cmd in buttons:
			ttk.Button(btn_frame, text=text, width=8,
			           command=cmd).pack(side=tk.LEFT, padx=1, expand=True, fill=tk.X)

	def update_status(self):
		"""更新连接状态"""
		if not LLM_AVAILABLE:
			self.status_label.config(text="● LLM模块未安装", foreground="gray")
			return

		if not hasattr(self.app, 'llm_service'):
			self.status_label.config(text="● 未初始化", foreground="gray")
			return

		if self.app.llm_service.is_available():
			cfg = self.app.llm_service.config
			self.status_label.config(
				text=f"● {cfg.provider}/{cfg.model}",
				foreground="green")
		elif self.app.llm_service.config.enabled:
			self.status_label.config(text="● 已启用(未连接)", foreground="orange")
		else:
			self.status_label.config(text="● 未启用", foreground="gray")

	def show_message(self, message, sender="AI"):
		"""显示消息"""
		timestamp = datetime.now().strftime("%H:%M:%S")

		def _update():
			try:
				self.chat_text.config(state=tk.NORMAL)
				self.chat_text.insert(tk.END, f"[{timestamp}] {sender}:\n{message}\n\n")
				self.chat_text.see(tk.END)
				self.chat_text.config(state=tk.DISABLED)
			except Exception:
				pass

		try:
			self.frame.after(0, _update)
		except Exception:
			pass

	def _send_message(self):
		"""发送消息"""
		msg = self.input_entry.get().strip()
		if not msg:
			return

		self.input_entry.delete(0, tk.END)
		self.show_message(msg, sender="你")

		if not hasattr(self.app, 'llm_service') or not self.app.llm_service.is_available():
			self.show_message("LLM 未配置或不可用，请先在设置中配置。")
			return

		def _do_chat():
			try:
				context = {}
				if self.app.dbe:
					context['components'] = self.app.available_comps
					context['database_loaded'] = True

				reply = self.app.llm_service.chat(msg, context)
				if reply:
					self.show_message(reply)
				else:
					self.show_message("未收到回复，请检查 LLM 配置。")
			except Exception as e:
				self.show_message(f"调用失败: {e}")

		threading.Thread(target=_do_chat, daemon=True).start()

	def _open_settings(self):
		"""打开 LLM 设置"""
		LLMSettingsDialog(self.app.root, self.app)
		self.update_status()

	def _test_connection(self):
		"""测试连接"""
		if not hasattr(self.app, 'llm_service'):
			self.show_message("LLM 未初始化")
			return

		def _do_test():
			success, msg = self.app.llm_service.test_connection()
			self.show_message(f"连接测试: {msg}")
			self.frame.after(0, self.update_status)

		self.show_message("正在测试连接...", sender="系统")
		threading.Thread(target=_do_test, daemon=True).start()

	def _get_llm(self):
		if not hasattr(self.app, 'llm_service') or not self.app.llm_service.is_available():
			self.show_message("LLM 未配置或不可用")
			return None
		return self.app.llm_service

	def _action_parse_comp(self):
		"""解析合金成分"""
		llm = self._get_llm()
		if not llm:
			return

		# 获取当前活动标签页的合金输入
		current_tab = self.app.get_current_tab()
		if current_tab and hasattr(current_tab, 'comps_entry'):
			text = current_tab.comps_entry.get().strip()
		elif current_tab and hasattr(current_tab, 'base_entry'):
			text = current_tab.base_entry.get().strip()
		else:
			self.show_message("请先在当前标签页输入合金成分")
			return

		if not text:
			self.show_message("合金输入框为空")
			return

		self.show_message(f"正在解析: {text}", sender="系统")

		def _do():
			result = llm.parse_composition(text, self.app.available_comps)
			if result:
				comp_str = ', '.join(f"{k}:{v}" for k, v in result.items())
				self.show_message(f"解析结果: {comp_str}")
			else:
				self.show_message("无法解析该成分表达式")

		threading.Thread(target=_do, daemon=True).start()

	def _action_predict_phases(self):
		"""推荐相"""
		llm = self._get_llm()
		if not llm:
			return

		selected_comps = self.app.get_selected_components()
		if not selected_comps:
			self.show_message("请先选择组分")
			return

		all_phases = [str(p) for p in self.app.dbe.phases.keys()] if self.app.dbe else []

		self.show_message(f"正在为 {selected_comps} 推荐相...", sender="系统")

		def _do():
			result = llm.predict_phases(selected_comps, all_phases, 'general')
			if result:
				self.show_message(f"推荐 {len(result)} 个相:\n{', '.join(result)}")
			else:
				self.show_message("无法获取推荐")

		threading.Thread(target=_do, daemon=True).start()

	def _action_suggest_params(self):
		"""推荐参数"""
		llm = self._get_llm()
		if not llm:
			return

		selected_comps = self.app.get_selected_components()
		current_tab = self.app.get_current_tab()
		calc_type = type(current_tab).__name__ if current_tab else 'general'

		self.show_message(f"正在推荐参数...", sender="系统")

		def _do():
			result = llm.suggest_parameters(selected_comps, calc_type)
			if result:
				parts = []
				if 'temp_min' in result:
					parts.append(f"温度: {result['temp_min']}-{result['temp_max']}K")
				if 'temp_step' in result:
					parts.append(f"步长: {result['temp_step']}K")
				if 'reasoning' in result:
					parts.append(f"理由: {result['reasoning']}")
				self.show_message("参数建议:\n" + "\n".join(parts))
			else:
				self.show_message("无法获取参数建议")

		threading.Thread(target=_do, daemon=True).start()

	def _action_recommend_model(self):
		"""推荐模型"""
		llm = self._get_llm()
		if not llm:
			return

		selected_comps = self.app.get_selected_components()
		current_tab = self.app.get_current_tab()
		calc_type = type(current_tab).__name__ if current_tab else 'general'

		self.show_message("正在推荐模型...", sender="系统")

		def _do():
			result = llm.recommend_model(selected_comps, calc_type)
			if result:
				self.show_message(
					f"推荐模型: {result.get('model', '未知')}\n"
					f"理由: {result.get('reasoning', '无')}")
			else:
				self.show_message("无法获取模型推荐")

		threading.Thread(target=_do, daemon=True).start()

	def _action_interpret_result(self):
		"""解读当前结果"""
		llm = self._get_llm()
		if not llm:
			return

		current_tab = self.app.get_current_tab()
		if not current_tab:
			self.show_message("没有活动标签页")
			return

		result = current_tab.calc_result
		if result is None:
			self.show_message("当前标签页没有计算结果")
			return

		self.show_message("正在分析结果...", sender="系统")

		def _do():
			context = {'calc_type': type(current_tab).__name__}
			result_dict = result.to_dict() if hasattr(result, 'to_dict') else {}
			interpretation = llm.interpret_result(result_dict, context)
			if interpretation:
				self.show_message(interpretation)
			else:
				self.show_message("无法获取结果解读")

		threading.Thread(target=_do, daemon=True).start()


# =============================================================================
# LLM 设置对话框
# =============================================================================
class LLMSettingsDialog(tk.Toplevel):
	"""LLM 设置对话框"""

	def __init__(self, parent, app):
		super().__init__(parent)
		self.app = app
		self.title("LLM 设置")
		self.geometry("500x480")
		self.resizable(False, False)
		self.transient(parent)
		self.grab_set()

		self._create_widgets()
		self._load_current()

		self.update_idletasks()
		x = parent.winfo_x() + (parent.winfo_width() - 500) // 2
		y = parent.winfo_y() + (parent.winfo_height() - 480) // 2
		self.geometry(f"+{x}+{y}")

	def _create_widgets(self):
		main = ttk.Frame(self, padding=15)
		main.pack(fill=tk.BOTH, expand=True)

		ttk.Label(main, text="LLM 配置", font=('', 14, 'bold')).pack(pady=(0, 10))

		# 启用开关
		self.enabled_var = tk.BooleanVar(value=False)
		ttk.Checkbutton(main, text="启用 LLM 辅助", variable=self.enabled_var,
		                command=self._on_toggle).pack(anchor=tk.W)

		ttk.Separator(main, orient='horizontal').pack(fill=tk.X, pady=8)

		# 提供商选择
		self.prov_frame = ttk.Frame(main)
		prov_frame = self.prov_frame
		prov_frame.pack(fill=tk.X, pady=3)
		ttk.Label(prov_frame, text="提供商:", width=10).pack(side=tk.LEFT)
		self.provider_var = tk.StringVar(value='ollama')
		self.provider_combo = ttk.Combobox(prov_frame, textvariable=self.provider_var,
		                                   values=['ollama', 'openai', 'anthropic', 'grok',
		                                           'deepseek', 'siliconflow', 'custom'],
		                                   width=15, state='readonly')
		self.provider_combo.pack(side=tk.LEFT, padx=5)
		self.provider_combo.bind('<<ComboboxSelected>>', self._on_provider_change)

		# API Key（本地模型如 Ollama 不需要，按提供商自动显隐）
		self.key_frame = ttk.Frame(main)
		self.key_frame.pack(fill=tk.X, pady=3)
		ttk.Label(self.key_frame, text="API Key:", width=10).pack(side=tk.LEFT)
		self.api_key_var = tk.StringVar()
		self.api_key_entry = ttk.Entry(self.key_frame, textvariable=self.api_key_var, show='*', width=35)
		self.api_key_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

		# Base URL
		url_frame = ttk.Frame(main)
		url_frame.pack(fill=tk.X, pady=3)
		ttk.Label(url_frame, text="Base URL:", width=10).pack(side=tk.LEFT)
		self.base_url_var = tk.StringVar()
		ttk.Entry(url_frame, textvariable=self.base_url_var, width=35).pack(
			side=tk.LEFT, padx=5, fill=tk.X, expand=True)

		# 模型（Combobox，可选择预设或自定义输入）
		model_frame = ttk.Frame(main)
		model_frame.pack(fill=tk.X, pady=3)
		ttk.Label(model_frame, text="模型:", width=10).pack(side=tk.LEFT)
		self.model_var = tk.StringVar()
		self.model_combo = ttk.Combobox(model_frame, textvariable=self.model_var, width=28)
		self.model_combo.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
		ttk.Button(model_frame, text="拉取模型",
		           command=self._fetch_models).pack(side=tk.RIGHT, padx=3)

		# 高级参数
		adv_frame = ttk.LabelFrame(main, text="高级参数", padding=5)
		adv_frame.pack(fill=tk.X, pady=5)

		adv_grid = ttk.Frame(adv_frame)
		adv_grid.pack(fill=tk.X)

		ttk.Label(adv_grid, text="Temperature:").grid(row=0, column=0, sticky=tk.W, padx=3)
		self.temp_var = tk.StringVar(value="0.3")
		ttk.Entry(adv_grid, textvariable=self.temp_var, width=8).grid(row=0, column=1, padx=3)

		ttk.Label(adv_grid, text="Max Tokens:").grid(row=0, column=2, sticky=tk.W, padx=3)
		self.max_tokens_var = tk.StringVar(value="2048")
		ttk.Entry(adv_grid, textvariable=self.max_tokens_var, width=8).grid(row=0, column=3, padx=3)

		ttk.Label(adv_grid, text="Timeout(s):").grid(row=1, column=0, sticky=tk.W, padx=3, pady=3)
		self.timeout_var = tk.StringVar(value="30")
		ttk.Entry(adv_grid, textvariable=self.timeout_var, width=8).grid(row=1, column=1, padx=3)

		# 按钮
		btn_frame = ttk.Frame(main)
		btn_frame.pack(fill=tk.X, pady=(10, 0))

		ttk.Button(btn_frame, text="测试连接", command=self._test).pack(side=tk.LEFT, padx=3)
		ttk.Button(btn_frame, text="保存", command=self._save).pack(side=tk.RIGHT, padx=3)
		ttk.Button(btn_frame, text="取消", command=self.destroy).pack(side=tk.RIGHT, padx=3)

		self.status_label = ttk.Label(main, text="", font=('', 9))
		self.status_label.pack(pady=5)

	def _load_current(self):
		"""加载当前配置"""
		if not LLM_AVAILABLE:
			self.status_label.config(text="LLM 模块未安装 (pip install httpx)", foreground="red")
			return

		if hasattr(self.app, 'llm_service'):
			cfg = self.app.llm_service.config
			self.enabled_var.set(cfg.enabled)
			self.provider_var.set(cfg.provider)
			self.api_key_var.set(cfg.api_key)
			self.base_url_var.set(cfg.base_url)
			self.model_var.set(cfg.model)
			self.temp_var.set(str(cfg.temperature))
			self.max_tokens_var.set(str(cfg.max_tokens))
			self.timeout_var.set(str(cfg.timeout))
			# 设置模型下拉列表
			if cfg.provider in PROVIDER_PRESETS:
				preset = PROVIDER_PRESETS[cfg.provider]
				self.model_combo['values'] = preset.get('models', [])
				# 本地模型隐藏 API Key
				if not preset.get('requires_key', True):
					self.key_frame.pack_forget()
			# Ollama: 自动拉取本地已安装模型
			if cfg.provider == 'ollama':
				self._fetch_ollama_models(cfg.base_url)

	def _on_provider_change(self, event=None):
		"""提供商切换时自动填充预设"""
		if not LLM_AVAILABLE:
			return
		provider = self.provider_var.get()
		if provider in PROVIDER_PRESETS:
			preset = PROVIDER_PRESETS[provider]
			self.base_url_var.set(preset['base_url'])
			self.model_var.set(preset['model'])
			self.model_combo['values'] = preset.get('models', [])
			# 需要 API Key 的厂商显示输入框，不需要的隐藏
			if preset.get('requires_key', True):
				self.key_frame.pack(fill=tk.X, pady=3, after=self.prov_frame)
			else:
				self.key_frame.pack_forget()
				self.api_key_var.set('')
			# Ollama: 自动拉取本地已安装模型
			if provider == 'ollama':
				self._fetch_ollama_models(preset['base_url'])

	def _fetch_models(self):
		"""拉取模型列表（按钮回调，适配所有厂商）"""
		provider = self.provider_var.get()
		base_url = self.base_url_var.get().strip()
		if not base_url:
			self.status_label.config(text="请先填写 Base URL", foreground="orange")
			return

		if provider == 'ollama':
			self._fetch_ollama_models(base_url)
		else:
			# 通用 OpenAI 兼容端点拉取
			api_key = self.api_key_var.get().strip()
			self._fetch_openai_compatible_models(base_url, api_key)

	def _fetch_openai_compatible_models(self, base_url, api_key):
		"""后台线程拉取 OpenAI 兼容端点的模型列表"""
		saved_model = self.model_var.get()
		self.status_label.config(text="正在获取模型列表...", foreground="gray")

		def _do_fetch():
			models, msg = fetch_openai_models(base_url, api_key)
			self.after(0, lambda: self._on_models_fetched(models, msg, saved_model))

		threading.Thread(target=_do_fetch, daemon=True).start()

	def _on_models_fetched(self, models, msg, saved_model):
		"""通用模型列表获取完毕，更新 Combobox"""
		if models:
			self.model_combo['values'] = models
			if saved_model in models:
				self.model_var.set(saved_model)
			else:
				self.model_var.set(models[0])
			self.status_label.config(text=msg, foreground="green")
		else:
			self.status_label.config(text=msg, foreground="orange")

	def _fetch_ollama_models(self, base_url):
		"""后台线程拉取 Ollama 本地已安装模型列表"""
		saved_model = self.model_var.get()
		self.status_label.config(text="正在获取本地模型列表...", foreground="gray")

		def _do_fetch():
			models, msg = fetch_ollama_models(base_url)
			# 回到主线程更新 UI
			self.after(0, lambda: self._on_ollama_models_fetched(
				models, msg, saved_model))

		threading.Thread(target=_do_fetch, daemon=True).start()

	def _on_ollama_models_fetched(self, models, msg, saved_model):
		"""Ollama 模型列表获取完毕，更新 Combobox"""
		if models:
			self.model_combo['values'] = models
			# 若之前选中的模型在列表中，保持选中；否则选第一个
			if saved_model in models:
				self.model_var.set(saved_model)
			else:
				self.model_var.set(models[0])
			self.status_label.config(text=msg, foreground="green")
		else:
			# 拉取失败，回退到预设列表
			preset_models = PROVIDER_PRESETS['ollama'].get('models', [])
			self.model_combo['values'] = preset_models
			self.status_label.config(text=msg, foreground="orange")

	def _on_toggle(self):
		pass

	def _test(self):
		"""测试连接"""
		if not LLM_AVAILABLE:
			self.status_label.config(text="LLM 模块未安装", foreground="red")
			return

		try:
			config = LLMConfig(
				provider=self.provider_var.get(),
				api_key=self.api_key_var.get(),
				base_url=self.base_url_var.get(),
				model=self.model_var.get(),
				temperature=float(self.temp_var.get()),
				max_tokens=int(self.max_tokens_var.get()),
				timeout=int(self.timeout_var.get()),
				enabled=True
			)
			service = LLMService(config)
			success, msg = service.test_connection()
			color = "green" if success else "red"
			self.status_label.config(text=msg, foreground=color)
		except Exception as e:
			self.status_label.config(text=f"测试失败: {e}", foreground="red")

	def _save(self):
		"""保存配置"""
		if not LLM_AVAILABLE:
			messagebox.showwarning("警告", "LLM 模块未安装")
			return

		try:
			config = LLMConfig(
				provider=self.provider_var.get(),
				api_key=self.api_key_var.get(),
				base_url=self.base_url_var.get(),
				model=self.model_var.get(),
				temperature=float(self.temp_var.get()),
				max_tokens=int(self.max_tokens_var.get()),
				timeout=int(self.timeout_var.get()),
				enabled=self.enabled_var.get()
			)

			if hasattr(self.app, 'llm_service'):
				self.app.llm_service.update_config(config)
			else:
				self.app.llm_service = LLMService(config)

			# 持久化
			self.app._save_config(llm=config.to_dict())

			if hasattr(self.app, 'llm_panel'):
				self.app.llm_panel.update_status()

			self.status_label.config(text="配置已保存", foreground="green")
			self.after(1000, self.destroy)

		except Exception as e:
			messagebox.showerror("错误", f"保存失败:\n{e}")


# =============================================================================
# 主应用类
# =============================================================================
class AlloyCalculatorApp:
	"""合金热力学计算工具 v3.0 - 单窗口多标签架构 + LLM 辅助"""

	SIDEBAR_WIDTH = 380

	def __init__(self, root):
		self.root = root
		self.root.title("合金热力学计算工具（UEM-Pycalphad）v3.0")
		self.root.geometry("1500x950")

		# 后端
		self.core = DatabaseManager()

		# 状态
		self.available_comps = []
		self.available_phases = []
		self.dbe = None
		self.pdens_window = None

		# 配置
		self.config_file = Path.home() / '.pycalphad_gui_config.json'
		config = self._load_config()
		self.last_system_directory = config.get('last_system_directory', os.getcwd())
		self.pdens = config.get('pdens', 2000)

		# 模型
		self.model_labels = {
			'RKM': 'R-K-M (Default)',
			'Muggianu': 'Muggianu',
			'Toop': 'Toop',
			'UEM1': 'UEM'
		}
		self.available_models = self.core.available_models_cls
		self.model_vars = {}
		self.uem1_liquid_only = tk.BooleanVar(value=False)
		self.comp_vars = {}

		# LLM
		self.llm_service = None
		if LLM_AVAILABLE:
			llm_config_data = config.get('llm', {})
			if llm_config_data:
				llm_config = LLMConfig.from_dict(llm_config_data)
				self.llm_service = LLMService(llm_config, logger=lambda m: None)
			else:
				self.llm_service = LLMService(LLMConfig(), logger=lambda m: None)

		# UI
		self._create_menu()
		self._create_main_layout()

		self.core.set_pdens(self.pdens)

	# =====================================================================
	# 配置
	# =====================================================================
	def _load_config(self):
		if self.config_file.exists():
			try:
				with open(self.config_file, 'r', encoding='utf-8') as f:
					config = json.load(f)
					for key in ['last_system_directory']:
						if key in config and not os.path.isdir(config[key]):
							config[key] = os.getcwd()
					return config
			except (json.JSONDecodeError, IOError) as e:
				print(f"加载配置文件失败: {e}")
		return {'last_system_directory': os.getcwd(), 'pdens': 2000}

	def _save_config(self, **kwargs):
		try:
			config = self._load_config()
			config.update(kwargs)
			with open(self.config_file, 'w', encoding='utf-8') as f:
				json.dump(config, f, indent=2, ensure_ascii=False)
		except IOError as e:
			print(f"保存配置文件失败: {e}")

	# =====================================================================
	# 菜单
	# =====================================================================
	def _create_menu(self):
		menubar = tk.Menu(self.root)
		self.root.config(menu=menubar)

		# 文件
		file_menu = tk.Menu(menubar, tearoff=0)
		menubar.add_cascade(label="文件(F)", menu=file_menu)
		file_menu.add_command(label="加载 TDB 数据库", command=self.load_system_database)
		file_menu.add_separator()
		file_menu.add_command(label="退出", command=self.root.quit)

		# 设置
		settings_menu = tk.Menu(menubar, tearoff=0)
		menubar.add_cascade(label="设置(S)", menu=settings_menu)
		settings_menu.add_command(label="点密度设置", command=self.open_pdens_settings)
		settings_menu.add_command(label="LLM 设置", command=self._open_llm_settings)

		# 帮助
		help_menu = tk.Menu(menubar, tearoff=0)
		menubar.add_cascade(label="帮助(H)", menu=help_menu)
		help_menu.add_command(label="关于", command=self._show_about)

	def _open_llm_settings(self):
		LLMSettingsDialog(self.root, self)
		if hasattr(self, 'llm_panel'):
			self.llm_panel.update_status()

	def _show_about(self):
		messagebox.showinfo("关于",
			"合金热力学计算工具 v3.0\n"
			"UEM-Pycalphad\n\n"
			"功能:\n"
			"- 液相线/固相线计算\n"
			"- 热力学性质计算\n"
			"- 伪二元相图\n"
			"- 三元等温相图\n"
			"- 三元液相面\n"
			"- 溶解度曲线\n"
			"- LLM 智能辅助")

	# =====================================================================
	# 主布局
	# =====================================================================
	def _create_main_layout(self):
		style = ttk.Style()
		style.configure('TLabelframe', font=('', 11, 'bold'))
		style.configure('TLabelframe.Label', font=('', 11, 'bold'))

		self.paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
		self.paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

		# 左侧边栏
		sidebar_container = ttk.Frame(self.paned, width=self.SIDEBAR_WIDTH)
		self.paned.add(sidebar_container, weight=0)
		self._create_sidebar(sidebar_container)

		# 右侧标签页
		notebook_frame = ttk.Frame(self.paned)
		self.paned.add(notebook_frame, weight=1)
		self._create_notebook(notebook_frame)

	def _create_sidebar(self, parent):
		canvas = tk.Canvas(parent, width=self.SIDEBAR_WIDTH - 20, highlightthickness=0)
		scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
		scrollable_frame = ttk.Frame(canvas)

		scrollable_frame.bind("<Configure>",
		                      lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

		window_id = canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")

		def _configure_canvas(event):
			canvas.itemconfig(window_id, width=event.width)

		canvas.bind("<Configure>", _configure_canvas)
		canvas.configure(yscrollcommand=scrollbar.set)

		canvas.pack(side="left", fill="both", expand=True)
		scrollbar.pack(side="right", fill="y")

		def _on_mousewheel(event):
			canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

		def _bind_mousewheel(event):
			canvas.bind_all("<MouseWheel>", _on_mousewheel)

		def _unbind_mousewheel(event):
			canvas.unbind_all("<MouseWheel>")

		canvas.bind("<Enter>", _bind_mousewheel)
		canvas.bind("<Leave>", _unbind_mousewheel)

		self._create_database_section(scrollable_frame)
		self._create_system_section(scrollable_frame)

		# LLM 助手面板
		self.llm_panel = LLMAssistantPanel(scrollable_frame, self)
		self.llm_panel.update_status()

	def _create_database_section(self, parent):
		frame = ttk.LabelFrame(parent, text="数据源", padding=8)
		frame.pack(fill=tk.X, pady=5, padx=5)

		ttk.Button(frame, text="加载 TDB 数据库",
		           command=self.load_system_database).pack(fill=tk.X, pady=3)

		status_frame = ttk.Frame(frame)
		status_frame.pack(fill=tk.X, pady=3)
		ttk.Label(status_frame, text="状态:", font=('', 10, 'bold')).pack(anchor=tk.W)
		self.system_db_label = ttk.Label(status_frame, text="未加载数据库",
		                                 foreground="red", wraplength=320, font=('', 10))
		self.system_db_label.pack(anchor=tk.W, pady=2)

	def _create_system_section(self, parent):
		frame = ttk.LabelFrame(parent, text="体系定义", padding=8)
		frame.pack(fill=tk.X, pady=5, padx=5)

		# 摘要行：显示已选组分数量 + 折叠/展开按钮
		self._comp_summary_frame = ttk.Frame(frame)
		self._comp_summary_frame.pack(fill=tk.X, pady=(0, 3))

		self._comp_expanded = False
		self._comp_toggle_btn = ttk.Button(
			self._comp_summary_frame, text="▶ 组分选择",
			command=self._toggle_comp_panel, width=14)
		self._comp_toggle_btn.pack(side=tk.LEFT)

		self._comp_summary_label = ttk.Label(
			self._comp_summary_frame, text="请先加载数据库",
			foreground="gray", font=('', 9))
		self._comp_summary_label.pack(side=tk.LEFT, padx=8)

		ttk.Button(self._comp_summary_frame, text="全选", width=5,
		           command=self._select_all_comps).pack(side=tk.RIGHT, padx=1)
		ttk.Button(self._comp_summary_frame, text="全不选", width=5,
		           command=self._deselect_all_comps).pack(side=tk.RIGHT, padx=1)

		# 可折叠的组分复选框面板（默认收起）
		self.comp_checkboxes_frame = ttk.Frame(frame)
		# 初始不 pack，保持折叠状态

	def _toggle_comp_panel(self):
		"""切换组分面板展开/收起"""
		self._comp_expanded = not self._comp_expanded
		if self._comp_expanded:
			self._comp_toggle_btn.config(text="▼ 组分选择")
			self.comp_checkboxes_frame.pack(fill=tk.X, pady=3,
			                               after=self._comp_summary_frame)
		else:
			self._comp_toggle_btn.config(text="▶ 组分选择")
			self.comp_checkboxes_frame.pack_forget()

	def _update_comp_summary(self):
		"""更新组分摘要标签"""
		if not self.comp_vars:
			self._comp_summary_label.config(text="请先加载数据库", foreground="gray")
			return
		total = len(self.comp_vars)
		selected = sum(1 for v in self.comp_vars.values() if v.get())
		self._comp_summary_label.config(
			text=f"已选 {selected}/{total}",
			foreground='#1565C0' if selected == total else '#E65100')

	def _create_models_section(self, parent):
		frame = ttk.LabelFrame(parent, text="模型选择 (液相线/热力学)", padding=8)
		frame.pack(fill=tk.X, pady=5, padx=5)

		checkbox_frame = ttk.Frame(frame)
		checkbox_frame.pack(fill=tk.X, pady=3)

		model_positions = {
			'RKM': (0, 0), 'Muggianu': (0, 1),
			'Toop': (1, 0), 'UEM1': (1, 1)
		}

		for model_key in self.available_models.keys():
			label = self.model_labels[model_key]
			var = tk.BooleanVar(value=(model_key == 'RKM'))
			self.model_vars[model_key] = var

			row, col = model_positions.get(model_key, (0, 0))
			cb = ttk.Checkbutton(checkbox_frame, text=label, variable=var)
			cb.grid(row=row, column=col, sticky=tk.W, padx=10, pady=2)

		self.cb_uem1_liq = ttk.Checkbutton(checkbox_frame, text="└─ Only Liquid",
		                                   variable=self.uem1_liquid_only)
		self.cb_uem1_liq.grid(row=2, column=1, sticky=tk.W, padx=15, pady=2)
		self.cb_uem1_liq.config(state='disabled')
		self.model_vars['UEM1'].trace_add('write', self._on_uem1_toggle)

		ttk.Label(frame, text="(伪二元/三元/液相面/溶解度标签页有独立模型选择)",
		          font=('', 8), foreground='gray').pack(anchor=tk.W)

	def _create_notebook(self, parent):
		self.notebook = ttk.Notebook(parent)
		self.notebook.pack(fill=tk.BOTH, expand=True)

		self.tab_liquidus = LiquidusSolidusTab(self.notebook, self)
		self.tab_props = ThermodynamicPropsTab(self.notebook, self)
		self.tab_pseudo_binary = PseudoBinaryTab(self.notebook, self)
		self.tab_ternary = TernaryTab(self.notebook, self)
		self.tab_surface = LiquidusSurfaceTab(self.notebook, self)
		self.tab_solubility = SolubilityTab(self.notebook, self)

		self.all_tabs = [
			self.tab_liquidus, self.tab_props,
			self.tab_pseudo_binary, self.tab_ternary,
			self.tab_surface, self.tab_solubility
		]

	# =====================================================================
	# 工具方法
	# =====================================================================
	def _on_uem1_toggle(self, *args):
		if self.model_vars['UEM1'].get():
			self.cb_uem1_liq.config(state='normal')
		else:
			self.cb_uem1_liq.config(state='disabled')
			self.uem1_liquid_only.set(False)

	def get_selected_models(self):
		selected_keys = []
		for model_key, var in self.model_vars.items():
			if var.get():
				if self.available_models[model_key] is None and model_key != 'RKM':
					pass
				else:
					selected_keys.append(model_key)

		if not selected_keys:
			messagebox.showwarning("警告", "请至少选择一个计算模型！")
			return []

		return selected_keys

	def get_selected_components(self):
		selected = []
		for comp, var in self.comp_vars.items():
			if var.get():
				selected.append(comp)
		return selected

	def get_current_tab(self):
		"""获取当前活动标签页"""
		try:
			idx = self.notebook.index(self.notebook.select())
			return self.all_tabs[idx]
		except Exception:
			return None

	# =====================================================================
	# 数据库加载
	# =====================================================================
	def load_system_database(self):
		file_paths = filedialog.askopenfilenames(
			title="选择系统TDB数据库文件（可多选）",
			initialdir=self.last_system_directory,
			filetypes=[("TDB文件", "*.tdb"), ("所有文件", "*.*")]
		)

		if not file_paths:
			return

		try:
			dir_path = os.path.dirname(file_paths[0])
			self._save_config(last_system_directory=dir_path)
			self.last_system_directory = dir_path

			result = self.core.load_database(file_paths)

			if isinstance(result, tuple) and len(result) == 2:
				loaded_files, db_obj = result
				self.dbe = db_obj
			else:
				loaded_files = result if isinstance(result, list) else file_paths
				self.dbe = self.core.dbe

			file_names = [os.path.basename(f) for f in loaded_files]
			self.system_db_label.config(
				text=f"已加载: {', '.join(file_names)}", foreground="green")

			self._update_database_ui()

		except Exception as e:
			messagebox.showerror("错误", f"数据库加载失败:\n{e}")

	def _update_database_ui(self):
		if self.dbe is not None:
			all_elements = set()
			for element in self.dbe.elements:
				if element.upper() != 'VA':
					all_elements.add(element.upper())
			self.available_comps = sorted(list(all_elements))
			self.available_phases = [str(phase) for phase in self.dbe.phases.keys()]
		else:
			self.available_comps = []
			self.available_phases = []

		self._update_component_checkboxes()

		# 更新所有标签页
		for tab in self.all_tabs:
			if hasattr(tab, 'update_comps'):
				tab.update_comps()

	def _update_component_checkboxes(self):
		for widget in self.comp_checkboxes_frame.winfo_children():
			widget.destroy()

		self.comp_vars.clear()

		if not self.available_comps:
			self._update_comp_summary()
			return

		num_cols = 4
		for i, comp in enumerate(self.available_comps):
			var = tk.BooleanVar(value=True)
			var.trace_add('write', lambda *_: self._update_comp_summary())
			self.comp_vars[comp] = var

			row = i // num_cols
			col = i % num_cols

			cb = ttk.Checkbutton(self.comp_checkboxes_frame, text=comp, variable=var, width=6)
			cb.grid(row=row, column=col, sticky=tk.W, padx=3, pady=2)

		# 更新摘要 + 自动展开（组分少于 8 个时自动展开）
		self._update_comp_summary()
		if len(self.available_comps) <= 8 and not self._comp_expanded:
			self._toggle_comp_panel()

	def _select_all_comps(self):
		for var in self.comp_vars.values():
			var.set(True)

	def _deselect_all_comps(self):
		for var in self.comp_vars.values():
			var.set(False)

	# =====================================================================
	# 点密度设置
	# =====================================================================
	def open_pdens_settings(self):
		if self.pdens_window is not None and self.pdens_window.winfo_exists():
			self.pdens_window.lift()
			return

		self.pdens_window = tk.Toplevel(self.root)
		self.pdens_window.title("点密度设置 (pdens)")
		self.pdens_window.geometry("400x350")
		self.pdens_window.resizable(False, False)

		main_frame = ttk.Frame(self.pdens_window, padding="15")
		main_frame.pack(fill=tk.BOTH, expand=True)

		ttk.Label(main_frame, text="平衡计算点密度设置",
		          font=('', 14, 'bold')).pack(pady=(0, 10))

		info_frame = ttk.LabelFrame(main_frame, text="说明", padding="8")
		info_frame.pack(fill=tk.X, pady=(0, 10))
		info_text = ("点密度 (pdens) 影响平衡计算的精度和速度\n"
		             "数值越大越精确但速度越慢")
		ttk.Label(info_frame, text=info_text, justify=tk.LEFT,
		          font=('', 9), foreground="#555").pack(anchor=tk.W)

		input_frame = ttk.Frame(main_frame)
		input_frame.pack(fill=tk.X, pady=(0, 10))
		ttk.Label(input_frame, text=f"当前值: {self.pdens}",
		          font=('', 10, 'bold')).pack(side=tk.LEFT, padx=(0, 10))
		entry_var = tk.StringVar(value=str(self.pdens))
		entry = ttk.Entry(input_frame, textvariable=entry_var, font=('', 11), width=15)
		entry.pack(side=tk.LEFT)

		preset_frame = ttk.LabelFrame(main_frame, text="预设值", padding="8")
		preset_frame.pack(fill=tk.X, pady=(0, 10))
		presets = [(500, "快速"), (1000, "一般"), (2000, "标准"), (3000, "高精度"), (5000, "超高")]
		for val, desc in presets:
			ttk.Button(preset_frame, text=f"{val} ({desc})",
			           command=lambda v=val: entry_var.set(str(v))).pack(side=tk.LEFT, padx=3)

		def save():
			try:
				v = int(entry_var.get())
				if v < 100 or v > 20000:
					messagebox.showwarning("警告", "建议值范围: 100-20000",
					                       parent=self.pdens_window)
				self.pdens = v
				self.core.set_pdens(v)
				self._save_config(pdens=v)
				self.pdens_window.destroy()
				self.pdens_window = None
			except ValueError:
				messagebox.showerror("错误", "请输入有效的整数", parent=self.pdens_window)

		def cancel():
			self.pdens_window.destroy()
			self.pdens_window = None

		btn_frame = ttk.Frame(main_frame)
		btn_frame.pack(fill=tk.X, pady=(10, 0))
		ttk.Button(btn_frame, text="保存", command=save, width=15).pack(side=tk.LEFT, padx=(0, 5))
		ttk.Button(btn_frame, text="取消", command=cancel, width=15).pack(side=tk.LEFT)

		entry.bind('<Return>', lambda e: save())
		self.pdens_window.bind('<Escape>', lambda e: cancel())
		entry.focus()
		entry.select_range(0, tk.END)


# =============================================================================
# 主程序入口
# =============================================================================
def main():
	"""主函数"""
	try:
		plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
		plt.rcParams['axes.unicode_minus'] = False
	except Exception as e:
		print(f"警告: 设置中文字体失败: {e}")

	root = tk.Tk()
	app = AlloyCalculatorApp(root)
	root.mainloop()


if __name__ == "__main__":
	import multiprocessing

	multiprocessing.freeze_support()
	main()
