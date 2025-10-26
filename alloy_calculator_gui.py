#!/usr/bin/env python
# -*- coding: utf-8 -*-


import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import numpy as np
import matplotlib

matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt
from pycalphad import Database, equilibrium, variables as v
from pycalphad.core.utils import get_pure_elements
from pycalphad import Model, calculate
import threading
from datetime import datetime
import os
import re
from pycalphad.advanced_uem_model import ModelMuggianu, ModelToop, ModelUEM1


# =============================================================================
# 主GUI类
# =============================================================================
class AlloyCalculatorGUI:
	"""
	合金热力学计算GUI主类

	主要功能：
	1. 数据库管理
	2. 液相线计算（多模型对比）
	3. 热力学性质计算（自由能、活度）
	4. 相图生成（伪二元相图）
	"""
	
	def __init__ (self, root):
		"""初始化GUI"""
		self.root = root
		self.root.title("合金热力学计算工具（UEM-Pycalphad）")
		self.root.geometry("1400x950")  # 增加高度以适应所有内容
		
		# 数据库相关
		self.dbe = None
		self.available_phases = []
		self.available_comps = []
		self.results_data = {}
		
		# 可用模型配置
		self._initialize_models()
		
		# 创建界面
		self.create_widgets()
	
	def _initialize_models (self):
		"""初始化可用的热力学模型"""
		self.available_models = {
			'RKM': Model,  # 默认Redlich-Kister-Muggianu模型
			'Muggianu': ModelMuggianu ,
			'Toop': ModelToop ,
			'UEM1': ModelUEM1
		}
		
		self.model_labels = {
			'RKM': 'R-K-M (Default)',
			'Muggianu': 'Muggianu',
			'Toop': 'Toop',
			'UEM1': 'UEM'
		}
		
		# UEM1特殊选项：是否只应用于液相
		self.uem1_liquid_only = tk.BooleanVar(value=False)
	
	# =========================================================================
	# GUI界面构建
	# =========================================================================
	def create_widgets (self):
		"""创建GUI界面组件"""
		# 配置全局样式
		style = ttk.Style()
		style.configure('TLabelframe', font=('', 11, 'bold'))
		style.configure('TLabelframe.Label', font=('', 11, 'bold'))

		# 左侧控制面板
		left_frame = ttk.Frame(self.root, width=500)
		left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=False, padx=5, pady=5)
		left_frame.pack_propagate(False)

		# 右侧结果显示区
		right_frame = ttk.Frame(self.root)
		right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)

		# 创建左侧控制组件
		self._create_control_panel(left_frame)

		# 创建右侧结果显示组件
		self._create_result_panel(right_frame)
	
	def _create_control_panel (self, parent):
		"""创建左侧控制面板"""
		
		# 1. 数据库加载区
		self._create_database_section(parent)
		
		# 2. 组分选择区
		self._create_component_section(parent)
		
		# 3. 成分扫描设置区
		self._create_composition_section(parent)
		
		# 4. 温度设置区
		self._create_temperature_section(parent)
		
		# 5. 模型选择区
		self._create_model_section(parent)
		
		# 6. 计算控制区
		self._create_control_section(parent)
	
	def _create_database_section (self, parent):
		"""数据库加载区域"""
		db_frame = ttk.LabelFrame(parent, text="1. 数据库", padding=10)
		db_frame.pack(fill=tk.X, pady=5)

		ttk.Button(db_frame, text="加载TDB数据库",
		           command=self.load_database).pack(fill=tk.X, pady=5)

		self.db_label = ttk.Label(db_frame, text="未加载",
		                          foreground="red", wraplength=400,
		                          font=('', 10))
		self.db_label.pack(pady=5)

		ttk.Label(db_frame, text="可用组分:",
		          font=('', 10, 'bold')).pack(anchor=tk.W, pady=(5, 2))
		self.comps_label = ttk.Label(db_frame, text="-",
		                             foreground="blue", wraplength=400,
		                             font=('', 10))
		self.comps_label.pack(anchor=tk.W, pady=(0, 5))

		ttk.Label(db_frame, text="可用相:",
		          font=('', 10, 'bold')).pack(anchor=tk.W, pady=(5, 2))
		self.phases_label = ttk.Label(db_frame, text="-",
		                              foreground="blue", wraplength=400,
		                              font=('', 10))
		self.phases_label.pack(anchor=tk.W)
	
	def _create_component_section (self, parent):
		"""组分选择区域"""
		comp_frame = ttk.LabelFrame(parent, text="2. 组分选择", padding=10)
		comp_frame.pack(fill=tk.X, pady=5)
		comp_frame.columnconfigure(1, weight=1)

		ttk.Label(comp_frame, text="研究组分 (逗号分隔):",
		          font=('', 10)).grid(row=0, column=0, sticky=tk.W, pady=3)

		self.comps_entry = ttk.Entry(comp_frame, width=30, font=('', 10))
		self.comps_entry.grid(row=0, column=1, sticky=tk.W + tk.E, pady=3, padx=(5, 0))
		self.comps_entry.insert(0, "AL,CR,NI")

		ttk.Label(comp_frame,
		          text="例: AL,CR (二元) 或 AL,CR,NI (三元)",
		          font=('', 9), foreground="gray").grid(
				row=1, column=0, columnspan=2, sticky=tk.W)
	
	def _create_composition_section (self, parent):
		"""成分扫描设置区域"""
		comp_frame = ttk.LabelFrame(parent, text="3. 成分扫描", padding=10)
		comp_frame.pack(fill=tk.X, pady=5)
		comp_frame.columnconfigure(1, weight=1)

		# 扫描组分（下拉框，从可用组分中选择）
		ttk.Label(comp_frame, text="扫描组分:",
		          font=('', 10)).grid(row=0, column=0, sticky=tk.W, pady=3)
		self.scan_comp_combobox = ttk.Combobox(comp_frame, width=18, font=('', 10), state='readonly')
		self.scan_comp_combobox.grid(row=0, column=1, sticky=tk.W + tk.E, pady=3, padx=(5, 0))

		# 扫描范围 - 分成三个独立输入框
		ttk.Label(comp_frame, text="扫描范围:",
		          font=('', 10)).grid(row=1, column=0, sticky=tk.W, pady=3)

		range_frame = ttk.Frame(comp_frame)
		range_frame.grid(row=1, column=1, sticky=tk.W + tk.E, pady=3, padx=(5, 0))

		ttk.Label(range_frame, text="从", font=('', 10)).pack(side=tk.LEFT, padx=2)
		self.scan_start_entry = ttk.Entry(range_frame, width=7, font=('', 10))
		self.scan_start_entry.pack(side=tk.LEFT, padx=2)
		self.scan_start_entry.insert(0, "0.1")

		ttk.Label(range_frame, text="到", font=('', 10)).pack(side=tk.LEFT, padx=2)
		self.scan_end_entry = ttk.Entry(range_frame, width=7, font=('', 10))
		self.scan_end_entry.pack(side=tk.LEFT, padx=2)
		self.scan_end_entry.insert(0, "0.9")

		ttk.Label(range_frame, text="共", font=('', 10)).pack(side=tk.LEFT, padx=2)
		self.scan_points_entry = ttk.Entry(range_frame, width=5, font=('', 10))
		self.scan_points_entry.pack(side=tk.LEFT, padx=2)
		self.scan_points_entry.insert(0, "10")

		ttk.Label(range_frame, text="点", font=('', 10)).pack(side=tk.LEFT, padx=2)

		# 保留旧的entry用于兼容性（从三个字段读取）
		self.scan_range_entry = None  # 标记为已废弃

		# 其他组分比例（不预加载，根据数据库动态调整）
		ttk.Label(comp_frame, text="其他组分比例 (冒号分隔):",
		          font=('', 10)).grid(row=2, column=0, sticky=tk.W, pady=3)
		self.other_ratio_entry = ttk.Entry(comp_frame, width=20, font=('', 10))
		self.other_ratio_entry.grid(row=2, column=1, sticky=tk.W + tk.E, pady=3, padx=(5, 0))

	def _create_temperature_section (self, parent):
		"""温度设置区域"""
		temp_frame = ttk.LabelFrame(parent, text="4. 温度范围 (K)", padding=10)
		temp_frame.pack(fill=tk.X, pady=5)
		temp_frame.columnconfigure(1, weight=1)

		ttk.Label(temp_frame, text="最低温度:",
		          font=('', 10)).grid(row=0, column=0, sticky=tk.W, pady=3)
		self.temp_min_entry = ttk.Entry(temp_frame, width=10, font=('', 10))
		self.temp_min_entry.grid(row=0, column=1, sticky=tk.W + tk.E, pady=3, padx=(5, 0))
		self.temp_min_entry.insert(0, "1400")

		ttk.Label(temp_frame, text="最高温度:",
		          font=('', 10)).grid(row=1, column=0, sticky=tk.W, pady=3)
		self.temp_max_entry = ttk.Entry(temp_frame, width=10, font=('', 10))
		self.temp_max_entry.grid(row=1, column=1, sticky=tk.W + tk.E, pady=3, padx=(5, 0))
		self.temp_max_entry.insert(0, "2200")

		ttk.Label(temp_frame, text="温度步长:",
		          font=('', 10)).grid(row=2, column=0, sticky=tk.W, pady=3)
		self.temp_step_entry = ttk.Entry(temp_frame, width=10, font=('', 10))
		self.temp_step_entry.grid(row=2, column=1, sticky=tk.W + tk.E, pady=3, padx=(5, 0))
		self.temp_step_entry.insert(0, "10")
	
	def _create_model_section (self, parent):
		"""模型选择区域"""
		model_frame = ttk.LabelFrame(parent, text="5. 模型选择（可多选对比）", padding=10)
		model_frame.pack(fill=tk.X, pady=5)

		self.model_vars = {}
		for model_key in self.available_models.keys():
			label = self.model_labels[model_key]
			var = tk.BooleanVar(value=(model_key == 'RKM'))  # 默认选RKM
			self.model_vars[model_key] = var

			# 创建带字体的Checkbutton
			cb_frame = ttk.Frame(model_frame)
			cb_frame.pack(anchor=tk.W, pady=3)

			cb = ttk.Checkbutton(cb_frame, text=label, variable=var)
			cb.pack(side=tk.LEFT)

			# 为checkbox添加字体（通过Style）
			# Note: ttk不直接支持font参数，但可以通过configure实现
			style = ttk.Style()
			style.configure('TCheckbutton', font=('', 10))

			# UEM1的特殊选项：仅用于液相（依赖于UEM1是否选中）
			if model_key == 'UEM1':
				self.cb_uem1_liq = ttk.Checkbutton(
						model_frame,
						text="└─ Only Liquid",
						variable=self.uem1_liquid_only)
				self.cb_uem1_liq.pack(anchor=tk.W, padx=20, pady=2)
				# 初始状态：禁用（因为UEM1默认未选中）
				self.cb_uem1_liq.config(state='disabled')

				# 添加回调：当UEM1状态改变时，更新"仅用于液相"的启用状态
				var.trace_add('write', self._on_uem1_toggle)
	
	def _on_uem1_toggle (self, *args):
		"""当UEM1选中状态改变时，更新"仅用于液相"选项的启用状态"""
		if self.model_vars['UEM1'].get():
			# UEM1被选中，启用"仅用于液相"选项
			self.cb_uem1_liq.config(state='normal')
		else:
			# UEM1未选中，禁用"仅用于液相"选项并清除勾选
			self.cb_uem1_liq.config(state='disabled')
			self.uem1_liquid_only.set(False)

	def _create_control_section (self, parent):
		"""计算控制区域"""
		control_frame = ttk.LabelFrame(parent, text="6. 计算控制", padding=10)
		control_frame.pack(fill=tk.X, pady=5)

		# 配置按钮样式
		style = ttk.Style()
		style.configure('TButton', font=('', 10), padding=5)

		# 配置列权重，让按钮均匀分布
		control_frame.columnconfigure(0, weight=1)
		control_frame.columnconfigure(1, weight=1)

		# 使用grid布局，2x2排列按钮
		ttk.Button(control_frame, text="液相线/固相线",
		           command=self.calculate_liquidus).grid(
				row=0, column=0, sticky=tk.W+tk.E, padx=3, pady=3)
		ttk.Button(control_frame, text="热力学性质",
		           command=self.calculate_properties).grid(
				row=0, column=1, sticky=tk.W+tk.E, padx=3, pady=3)
		ttk.Button(control_frame, text="伪二元相图",
		           command=self.calculate_phase_diagram).grid(
				row=1, column=0, sticky=tk.W+tk.E, padx=3, pady=3)
		ttk.Button(control_frame, text="清除结果",
		           command=self.clear_results).grid(
				row=1, column=1, sticky=tk.W+tk.E, padx=3, pady=3)

		# 进度提示（跨两列）
		self.progress_var = tk.StringVar(value="就绪")
		ttk.Label(control_frame, textvariable=self.progress_var,
		          foreground="green", font=('', 10, 'bold')).grid(
				row=2, column=0, columnspan=2, pady=5)
	
	def _create_result_panel (self, parent):
		"""创建右侧结果显示面板"""
		self.notebook = ttk.Notebook(parent)
		self.notebook.pack(fill=tk.BOTH, expand=True)
		
		# 标签页1: 液相线对比图
		self._create_liquidus_tab()
		
		# 标签页2: 热力学性质
		self._create_properties_tab()
		
		# 标签页3: 相图
		self._create_phase_diagram_tab()
		
		# 标签页4: 数据表格
		self._create_data_tab()
		
		# 标签页5: 日志
		self._create_log_tab()
	
	def _create_liquidus_tab (self):
		"""创建液相线对比标签页"""
		liquidus_tab = ttk.Frame(self.notebook)
		self.notebook.add(liquidus_tab, text="液相线/固相线")
		
		self.fig_liquidus = Figure(figsize=(8, 6))
		self.ax_liquidus = self.fig_liquidus.add_subplot(111)
		self.canvas_liquidus = FigureCanvasTkAgg(self.fig_liquidus, liquidus_tab)
		self.canvas_liquidus.get_tk_widget().pack(fill=tk.BOTH, expand=True)
		
		toolbar = NavigationToolbar2Tk(self.canvas_liquidus, liquidus_tab)
		toolbar.update()
	
	def _create_properties_tab (self):
		"""创建热力学性质标签页"""
		props_tab = ttk.Frame(self.notebook)
		self.notebook.add(props_tab, text="热力学性质")
		
		props_notebook = ttk.Notebook(props_tab)
		props_notebook.pack(fill=tk.BOTH, expand=True)
		
		# 子标签1: Gibbs自由能
		gibbs_frame = ttk.Frame(props_notebook)
		props_notebook.add(gibbs_frame, text="Gibbs自由能 (GM)")
		
		self.fig_gibbs = Figure(figsize=(8, 5))
		self.ax_gibbs = self.fig_gibbs.add_subplot(111)
		self.canvas_gibbs = FigureCanvasTkAgg(self.fig_gibbs, gibbs_frame)
		self.canvas_gibbs.get_tk_widget().pack(fill=tk.BOTH, expand=True)
		
		toolbar_gibbs = NavigationToolbar2Tk(self.canvas_gibbs, gibbs_frame)
		toolbar_gibbs.update()
		
		# 子标签2: 活度
		activity_frame = ttk.Frame(props_notebook)
		props_notebook.add(activity_frame, text="活度 (Activity)")
		
		self.fig_activity = Figure(figsize=(8, 5))
		self.ax_activity = self.fig_activity.add_subplot(111)
		self.canvas_activity = FigureCanvasTkAgg(self.fig_activity, activity_frame)
		self.canvas_activity.get_tk_widget().pack(fill=tk.BOTH, expand=True)
		
		toolbar_activity = NavigationToolbar2Tk(self.canvas_activity, activity_frame)
		toolbar_activity.update()
	
	def _create_phase_diagram_tab (self):
		"""创建相图标签页"""
		phase_tab = ttk.Frame(self.notebook)
		self.notebook.add(phase_tab, text="伪二元相图")
		
		self.fig_phase = Figure(figsize=(8, 6))
		self.canvas_phase = FigureCanvasTkAgg(self.fig_phase, phase_tab)
		self.canvas_phase.get_tk_widget().pack(fill=tk.BOTH, expand=True)
		
		toolbar = NavigationToolbar2Tk(self.canvas_phase, phase_tab)
		toolbar.update()
	
	def _create_data_tab (self):
		"""创建数据表格标签页"""
		data_tab = ttk.Frame(self.notebook)
		self.notebook.add(data_tab, text="数据")
		
		self.data_text = scrolledtext.ScrolledText(
				data_tab, wrap=tk.NONE, font=('Courier', 10))
		self.data_text.pack(fill=tk.BOTH, expand=True)
		
		ttk.Button(data_tab, text="导出数据",
		           command=self.export_data).pack(pady=5)
	
	def _create_log_tab (self):
		"""创建日志标签页"""
		log_tab = ttk.Frame(self.notebook)
		self.notebook.add(log_tab, text="日志")
		
		self.log_text = scrolledtext.ScrolledText(
				log_tab, wrap=tk.WORD, font=('Courier', 9))
		self.log_text.pack(fill=tk.BOTH, expand=True)
	
	# =========================================================================
	# 工具函数
	# =========================================================================
	def log (self, message):
		"""
		添加日志信息

		参数:
			message: 日志消息
		"""
		timestamp = datetime.now().strftime("%H:%M:%S")
		self.log_text.insert(tk.END, f"[{timestamp}] {message}\n")
		self.log_text.see(tk.END)
		self.root.update_idletasks()
	
	def get_selected_models (self):
		"""
		获取用户选中的模型

		返回:
			list: 选中的模型键名列表，如 ['RKM', 'Muggianu']
		"""
		selected_keys = []
		for model_key, var in self.model_vars.items():
			if var.get():
				if self.available_models[model_key] is None and model_key != 'RKM':
					self.log(f"警告: {self.model_labels[model_key]} 模型不可用，跳过")
				else:
					selected_keys.append(model_key)
		
		if not selected_keys:
			messagebox.showwarning("警告", "请至少选择一个可用的模型！")
			return None
		
		return selected_keys
	
	def get_model_spec (self, model_key):
		"""
		根据模型键获取传递给equilibrium的模型参数

		参数:
			model_key: 模型键名，如 'RKM', 'Muggianu'等

		返回:
			模型类或字典或None
		"""
		model_class = self.available_models.get(model_key)
		
		if model_key == 'RKM':
			return None  # 使用默认RKM
		elif model_key == 'UEM1' and model_class:
			if self.uem1_liquid_only.get():
				# 只替换液相
				return {ph: model_class if 'LIQUID' in ph.upper() else Model
				        for ph in self.available_phases}
			else:
				# 应用于所有相（利用自动回退）
				return model_class
		elif model_class:  # Muggianu, Toop
			return model_class
		else:
			return None
	
	# =========================================================================
	# 数据库管理
	# =========================================================================
	def load_database (self):
		"""
		加载TDB数据库文件 (已修正为支持多文件 - v3)
		"""
		file_paths = filedialog.askopenfilename(
				title="选择TDB数据库文件 (可多选, e.g., 'alcrni.tdb' + 'pure.tdb')",
				filetypes=[("TDB文件", "*.tdb"), ("所有文件", "*.*")],
				multiple=True
		)
		
		if not file_paths:  # file_paths 是一个元组
			return
		
		try:
			# (!!) 核心修正:
			# pycalphad的Database()构造函数不接受多参数。
			# 正确的方法是：读取所有TDB文件的内容，将它们合并(concatenate)成一个
			# 包含换行符的单一字符串，然后将该字符串传递给Database()。
			
			all_tdb_content = ""
			loaded_files_list = []
			
			for path in file_paths:
				try:
					# 确保使用 utf-8 或 'latin-1'，TDB文件编码可能不同
					with open(path, 'r', encoding='latin-1') as f:
						all_tdb_content += f.read() + "\n"  # 添加换行符确保文件间分离
					loaded_files_list.append(os.path.basename(path))
				except Exception as file_read_e:
					# 如果 latin-1 失败，尝试 utf-8
					try:
						with open(path, 'r', encoding='utf-8') as f:
							all_tdb_content += f.read() + "\n"
						loaded_files_list.append(os.path.basename(path))
					except Exception as e:
						self.log(f"读取文件失败: {path} - {e}")
						raise e  # 重新抛出异常以触发外层try-except
			
			# (!!) 修正: 将合并后的单一字符串传递给 Database
			self.dbe = Database(all_tdb_content)
			
			loaded_files_str = ", ".join(loaded_files_list)
			self.db_label.config(
					text=f"已加载: {loaded_files_str}",
					foreground="green")
			
			# 提取可用组分和相
			# (!! 以下代码与之前相同，无需修改 !!)
			all_elements = []
			for element in self.dbe.elements:
				if element != 'VA':
					elem_str = str(element).strip().upper()  # 转大写并去空格
					if elem_str:  # 确保非空
						all_elements.append(elem_str)
			
			self.log(f"从数据库提取的元素: {all_elements}")
			
			# 过滤出有效的元素符号
			self.available_comps = []
			for elem in all_elements:
				if len(elem) <= 2 and elem[0].isalpha():
					if len(elem) == 1:
						standardized = elem.upper()
					else:
						standardized = elem[0].upper() + elem[1].lower()
					self.available_comps.append(standardized)
			
			self.available_comps = sorted(list(set(self.available_comps)))
			self.available_phases = sorted(list(self.dbe.phases.keys()))
			
			# 更新界面显示
			if self.available_comps:
				self.comps_label.config(
						text=f"{', '.join(self.available_comps)}",
						foreground="blue")
			else:
				self.comps_label.config(
						text="警告: 未检测到有效组分！",
						foreground="red")
			
			self.phases_label.config(text=f"{', '.join(self.available_phases)}")
			
			self.log(f"成功加载数据库: {loaded_files_str}")
			self.log(f"可用组分 ({len(self.available_comps)}): {self.available_comps}")
			self.log(f"可用相 ({len(self.available_phases)}): {self.available_phases}")

			# 自动填充研究组分（从可用组分中提取）
			if self.available_comps:
				# 将可用组分填充到研究组分输入框
				self.comps_entry.delete(0, tk.END)
				self.comps_entry.insert(0, ",".join(self.available_comps))

				# 更新扫描组分下拉框的选项
				self.scan_comp_combobox['values'] = self.available_comps
				# 默认选择第一个组分
				if len(self.available_comps) > 0:
					self.scan_comp_combobox.current(0)

				self.log(f"已自动填充研究组分: {','.join(self.available_comps)}")

		except Exception as e:
			messagebox.showerror("错误", f"加载数据库失败:\n{e}")
			self.log(f"数据库加载失败: {e}")
			import traceback
			self.log(traceback.format_exc())
			self.dbe = None
			self.db_label.config(text="加载失败", foreground="red")
	# =========================================================================
	# 输入解析
	# =========================================================================
	def _parse_inputs (self):
		"""
		解析和验证用户输入 (已修正大小写问题)

		返回:
			dict: 包含所有解析后参数的字典，失败返回None
		"""
		if self.dbe is None:
			messagebox.showwarning("警告", "请先加载数据库！")
			return None
		
		try:
			# 解析研究组分
			comps_str = self.comps_entry.get().strip().upper()  # -> "AL,CR,NI"
			if not comps_str:
				raise ValueError("研究组分不能为空！")
			
			study_comps = [c.strip() for c in comps_str.split(',') if c.strip()]  # -> ['AL', 'CR', 'NI']
			
			# (!!) 修正: 我们的内部标准是 ALL-UPPERCASE
			available_comps_upper = [c.upper() for c in self.available_comps]
			
			invalid_comps = [c for c in study_comps
			                 if c not in available_comps_upper]
			if invalid_comps:
				raise ValueError(
						f"组分无效: {invalid_comps}\n"
						f"可用组分: {self.available_comps}")
			
			# (!!) 移除错误的混合大小写标准化
			# study_comps = [comp_map[c] for c in study_comps] # <--- 这行是错误的根源
			
			# study_comps 已经是 ['AL', 'CR', 'NI'] (全大写)
			
			if len(study_comps) < 2:
				raise ValueError("至少需要2个研究组分！")
			
			study_comps.append('VA')
			study_comps = sorted(list(set(study_comps)))  # -> ['AL', 'CR', 'NI', 'VA']
			
			# 解析扫描组分（从下拉框获取）
			scan_comp_input = self.scan_comp_combobox.get().strip().upper()
			if not scan_comp_input:
				raise ValueError("扫描组分不能为空！请从下拉框中选择")
			
			# (!!) 修正: 使用全大写进行验证
			if scan_comp_input not in available_comps_upper:
				raise ValueError(
						f"扫描组分 '{scan_comp_input}' 无效\n"
						f"可用组分: {self.available_comps}")
			scan_comp = scan_comp_input  # -> 'NI' (全大写)
			
			if scan_comp not in study_comps or scan_comp == 'VA':
				raise ValueError(
						f"扫描组分 '{scan_comp}' 必须是研究组分之一且不能是VA")
			
			# 解析扫描范围（从三个独立字段读取）
			try:
				start = float(self.scan_start_entry.get().strip())
				stop = float(self.scan_end_entry.get().strip())
				num = int(self.scan_points_entry.get().strip())
			except ValueError as e:
				raise ValueError(f"扫描范围输入无效: {e}")

			if not (0 < start < 1 and 0 < stop < 1):
				raise ValueError("扫描起点和终点必须在 (0, 1) 之间")
			if start >= stop:
				raise ValueError("扫描起点必须小于终点")
			if num < 2:
				raise ValueError("扫描点数必须至少为2")

			x_scan_range = np.linspace(start, stop, num)
			
			# 解析温度范围
			temp_min = float(self.temp_min_entry.get())
			temp_max = float(self.temp_max_entry.get())
			temp_step = float(self.temp_step_entry.get())
			if temp_min >= temp_max or temp_step <= 0:
				raise ValueError("温度范围或步长无效")
			
			# 解析其他组分比例
			# (!!) 修正: 这里的 study_comps 已经是全大写
			other_comps = [c for c in study_comps
			               if c != scan_comp and c != 'VA']  # -> ['AL', 'CR'] (全大写)
			
			if not other_comps:  # 二元系
				if len(study_comps) - 1 != 2:
					raise ValueError("二元系配置错误")
				ratios = []
			else:  # 多元系
				other_ratio_str = self.other_ratio_entry.get().strip()
				if not other_ratio_str:
					raise ValueError("其他组分比例不能为空（多元系）")
				
				ratio_parts = [p.strip() for p in other_ratio_str.split(':')]
				if len(ratio_parts) != len(other_comps):
					raise ValueError(
							f"比例数量 ({len(ratio_parts)}) 与"
							f"其他组分数量 ({len(other_comps)}) 不匹配")
				
				ratios = [float(r) for r in ratio_parts]
				if any(r < 0 for r in ratios) or sum(ratios) <= 0:
					raise ValueError("比例值必须为正数且总和大于0")
			
			# 返回解析结果
			return {
				"study_comps": study_comps,  # ['AL', 'CR', 'NI', 'VA']
				"scan_comp": scan_comp,  # 'NI'
				"x_scan_range": x_scan_range,
				"temp_min": temp_min,
				"temp_max": temp_max,
				"temp_step": temp_step,
				"other_comps": other_comps,  # ['AL', 'CR']
				"ratios": ratios,
				"db_phases": list(self.dbe.phases.keys())
			}
		
		except ValueError as ve:
			messagebox.showerror("输入错误", str(ve))
			self.log(f"输入错误: {ve}")
			return None
		except Exception as e:
			messagebox.showerror("错误", f"解析输入错误: {e}")
			self.log(f"解析输入错误: {e}")
			return None
	# =========================================================================
	# 模块1: 液相线温度计算
	# =========================================================================
	def calculate_liquidus (self):
		"""
		计算液相线/固相线温度（主入口）

		功能：
		- 液相线温度（Liquidus）：完全为液相的温度
		- 固相线温度（Solidus）：开始出现液相的温度
		支持多模型对比
		"""
		inputs = self._parse_inputs()
		if not inputs:
			return
		
		selected_model_keys = self.get_selected_models()
		if not selected_model_keys:
			return
		
		thread = threading.Thread(
				target=self._calculate_liquidus_thread,
				args=(selected_model_keys, inputs))
		thread.daemon = True
		thread.start()
	
	def _calculate_liquidus_thread (self, selected_model_keys, inputs):
		"""
		液相线/固相线计算线程

		计算逻辑：
		1. 遍历不同成分点
		2. 对每个成分点，扫描温度范围
		3. 找到液相线温度（全液相，液相分数≥99.5%）
		4. 找到固相线温度（开始出现液相，液相分数≥0.5%）
		"""
		try:
			self.results_data = {
				'x_scan': inputs['x_scan_range'],
				'scan_comp': inputs['scan_comp']
			}
			ratio_sum = sum(inputs['ratios']) if inputs['ratios'] else 1.0
			
			for model_key in selected_model_keys:
				model_label = self.model_labels[model_key]
				if model_key == 'UEM1' and self.uem1_liquid_only.get():
					model_label += " (Liq Only)"

				self.log(f"\n{'=' * 60}")
				self.log(f"[液相线/固相线计算] 模型: {model_label}")
				self.log(f"{'=' * 60}")
				self.progress_var.set(f"计算液相线/固相线: {model_label}")

				model_spec = self.get_model_spec(model_key)
				if model_spec is None and model_key != 'RKM':
					continue

				liquidus_temps = []
				solidus_temps = []

				for idx, x_scan in enumerate(inputs['x_scan_range']):
					# 计算当前成分
					composition = self._calculate_composition(
							x_scan, inputs, ratio_sum)

					# 查找液相线温度
					T_liq = self._find_liquidus_temperature(
							inputs, composition, model_spec)

					# 查找固相线温度
					T_sol = self._find_solidus_temperature(
							inputs, composition, model_spec)

					liquidus_temps.append(
							float(T_liq) if T_liq is not None
							                and not np.isnan(T_liq) else np.nan)

					solidus_temps.append(
							float(T_sol) if T_sol is not None
							                and not np.isnan(T_sol) else np.nan)

					# 记录日志
					comp_str = ", ".join([f"X({k})={v:.3f}"
					                      for k, v in composition.items()])
					liq_str = f"{float(T_liq):.1f} K" if T_liq is not None and not np.isnan(T_liq) else "未找到"
					sol_str = f"{float(T_sol):.1f} K" if T_sol is not None and not np.isnan(T_sol) else "未找到"
					self.log(f"  {comp_str} -> T_liq = {liq_str}, T_sol = {sol_str}")

				# 存储结果
				self.results_data[model_key] = {
					'liquidus': np.array(liquidus_temps),
					'solidus': np.array(solidus_temps),
					'type': 'liquidus_solidus',
					'label': model_label
				}
				self.log(f"{model_label} 液相线/固相线计算完成！\n")
			
			# 绘图和数据显示
			self._plot_liquidus_comparison()
			self._display_liquidus_data()
			self.progress_var.set("液相线计算完成！")
			self.log("\n所有液相线计算完成！")
		
		except Exception as e:
			self.log(f"液相线计算错误: {e}")
			import traceback
			self.log(traceback.format_exc())
			self.progress_var.set("计算失败")
	
	def _calculate_composition (self, x_scan, inputs, ratio_sum):
		"""
		计算给定扫描组分值时的合金成分

		参数:
			x_scan: 扫描组分的摩尔分数
			inputs: 输入参数字典
			ratio_sum: 其他组分比例之和

		返回:
			dict: 各组分摩尔分数字典
		"""
		composition = {inputs['scan_comp']: x_scan}
		remaining = 1.0 - x_scan
		
		if inputs['other_comps']:  # 多元系
			for i, comp in enumerate(inputs['other_comps']):
				composition[comp] = remaining * inputs['ratios'][i] / ratio_sum
		else:  # 二元系
			other_comp = [c for c in inputs['study_comps']
			              if c != inputs['scan_comp'] and c != 'VA'][0]
			composition[other_comp] = remaining
		
		return composition
	
	def _find_liquidus_temperature (self, inputs, composition, model_spec):
		"""
		查找给定成分的液相线温度

		方法：从高温向低温扫描，找到第一个全液相点

		参数:
			inputs: 输入参数
			composition: 成分字典
			model_spec: 模型规格

		返回:
			float: 液相线温度(K)，失败返回NaN
		"""
		T_range = (inputs['temp_min'], inputs['temp_max'], inputs['temp_step'])
		conditions = {v.T: T_range, v.P: 101325}
		
		# 设置成分条件
		scan_comp = inputs.get('scan_comp')
		comps_to_set = [
			c for c in inputs['study_comps']
			if c != 'VA' and c != scan_comp
		]
		
		
			
		comp_lookup = {k.upper(): value for k, value in composition.items()}
		for comp in comps_to_set:
			if comp in comp_lookup:
				conditions[v.X(comp)] = comp_lookup[comp]
			else:
				self.log(f"警告: 组分 {comp} 未在 'composition' 字典中找到。")
		
		try:
			active_elements = [c for c in inputs['study_comps'] if c != 'VA']
			eq = equilibrium(
					self.dbe, active_elements, inputs['db_phases'],
					model=model_spec, conditions=conditions,
					calc_opts={'pdens': 50})
			
			T_vals = eq.T.values.squeeze()
			phase_array = eq.Phase.values.squeeze()
			np_array = eq.NP.values.squeeze()
			
			# 处理单点情况
			if T_vals.ndim == 0:
				T_vals = np.array([T_vals.item()])
				phase_array = np.array([[phase_array.item()]])
				np_array = np.array([[np_array.item()]])
			
			liquidus_T = np.nan
			all_liquid_found = False
			
			# 从高温向低温扫描
			for t_idx in range(len(T_vals) - 1, -1, -1):
				is_liquid_present = False
				liquid_fraction = 0.0
				
				phases_at_T = phase_array[t_idx]
				fracs_at_T = np_array[t_idx]
				
				# 处理单个值情况
				if isinstance(phases_at_T, (str, bytes)):
					phases_at_T = [phases_at_T]
				if isinstance(fracs_at_T, (float, np.float64)):
					fracs_at_T = [fracs_at_T]
				
				# 检查液相
				for phase_name, frac in zip(phases_at_T, fracs_at_T):
					if phase_name == '':
						continue
					if isinstance(phase_name, bytes):
						phase_name = phase_name.decode('utf-8')
					
					if 'LIQUID' in phase_name.upper():
						is_liquid_present = True
						liquid_fraction = float(frac)
						break
				
				# 判断是否为全液相
				if is_liquid_present and liquid_fraction >= 0.995:
					liquidus_T = float(T_vals[t_idx])
					all_liquid_found = True
				elif all_liquid_found:
					break
			
			return liquidus_T
		
		except Exception as e:
			comp_str = ", ".join([f"X({k})={v:.3f}"
			                      for k, v in composition.items()])
			self.log(f"  查找液相线失败 (Comp={comp_str}): {e}")
			return np.nan

	def _find_solidus_temperature (self, inputs, composition, model_spec):
		"""
		查找给定成分的固相线温度

		方法：从低温向高温扫描，找到第一个开始出现液相的点

		参数:
			inputs: 输入参数
			composition: 成分字典
			model_spec: 模型规格

		返回:
			float: 固相线温度(K)，失败返回NaN
		"""
		T_range = (inputs['temp_min'], inputs['temp_max'], inputs['temp_step'])
		conditions = {v.T: T_range, v.P: 101325}

		# 设置成分条件
		scan_comp = inputs.get('scan_comp')
		comps_to_set = [
			c for c in inputs['study_comps']
			if c != 'VA' and c != scan_comp
		]

		comp_lookup = {k.upper(): value for k, value in composition.items()}
		for comp in comps_to_set:
			if comp in comp_lookup:
				conditions[v.X(comp)] = comp_lookup[comp]
			else:
				self.log(f"警告: 组分 {comp} 未在 'composition' 字典中找到。")

		try:
			active_elements = [c for c in inputs['study_comps'] if c != 'VA']
			eq = equilibrium(
					self.dbe, active_elements, inputs['db_phases'],
					model=model_spec, conditions=conditions,
					calc_opts={'pdens': 50})

			T_vals = eq.T.values.squeeze()
			phase_array = eq.Phase.values.squeeze()
			np_array = eq.NP.values.squeeze()

			# 处理单点情况
			if T_vals.ndim == 0:
				T_vals = np.array([T_vals.item()])
				phase_array = np.array([[phase_array.item()]])
				np_array = np.array([[np_array.item()]])

			solidus_T = np.nan

			# 从低温向高温扫描，找到第一个出现液相的点
			for t_idx in range(len(T_vals)):
				is_liquid_present = False
				liquid_fraction = 0.0

				phases_at_T = phase_array[t_idx]
				fracs_at_T = np_array[t_idx]

				# 处理单个值情况
				if isinstance(phases_at_T, (str, bytes)):
					phases_at_T = [phases_at_T]
				if isinstance(fracs_at_T, (float, np.float64)):
					fracs_at_T = [fracs_at_T]

				# 检查液相
				for phase_name, frac in zip(phases_at_T, fracs_at_T):
					if phase_name == '':
						continue
					if isinstance(phase_name, bytes):
						phase_name = phase_name.decode('utf-8')

					if 'LIQUID' in phase_name.upper():
						is_liquid_present = True
						liquid_fraction = float(frac)
						break

				# 判断是否开始出现液相（固相线定义：液相分数 > 0.5%）
				if is_liquid_present and liquid_fraction >= 0.005:
					solidus_T = float(T_vals[t_idx])
					break

			return solidus_T

		except Exception as e:
			comp_str = ", ".join([f"X({k})={v:.3f}"
			                      for k, v in composition.items()])
			self.log(f"  查找固相线失败 (Comp={comp_str}): {e}")
			return np.nan

	def _plot_liquidus_comparison (self):
		"""绘制液相线/固相线对比图"""
		self.ax_liquidus.clear()

		if not self.results_data or 'x_scan' not in self.results_data:
			return

		x_scan = self.results_data['x_scan']
		scan_comp = self.results_data['scan_comp']

		colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
		markers = ['o', 's', '^', 'D']
		idx = 0
		plotted_something = False

		# 按固定顺序绘制
		for model_key in ['RKM', 'UEM1', 'Muggianu', 'Toop']:
			if (model_key in self.results_data and
					self.results_data[model_key]['type'] == 'liquidus_solidus'):

				liquidus = self.results_data[model_key]['liquidus']
				solidus = self.results_data[model_key]['solidus']
				label = self.results_data[model_key]['label']

				# 绘制液相线（实线）
				valid_liq = ~np.isnan(liquidus)
				if np.any(valid_liq):
					self.ax_liquidus.plot(
							x_scan[valid_liq], liquidus[valid_liq],
							marker=markers[idx], linestyle='-', linewidth=2,
							markersize=5, label=f'{label} (Liquidus)',
							color=colors[idx], alpha=0.9)
					plotted_something = True

				# 绘制固相线（虚线）
				valid_sol = ~np.isnan(solidus)
				if np.any(valid_sol):
					self.ax_liquidus.plot(
							x_scan[valid_sol], solidus[valid_sol],
							marker=markers[idx], linestyle='--', linewidth=2,
							markersize=5, label=f'{label} (Solidus)',
							color=colors[idx], alpha=0.7)
					plotted_something = True

				idx += 1

		if plotted_something:
			self.ax_liquidus.set_xlabel(f'X({scan_comp})', fontsize=12)
			self.ax_liquidus.set_ylabel('Temperature (K)', fontsize=12)
			self.ax_liquidus.set_title(
					'Liquidus & Solidus Temperature Comparison',
					fontsize=14, fontweight='bold')
			self.ax_liquidus.legend(fontsize=9, loc='best', ncol=1)
			self.ax_liquidus.grid(True, alpha=0.3, linestyle='--')
		else:
			self.ax_liquidus.text(
					0.5, 0.5, 'No valid data to plot',
					ha='center', va='center')

		self.fig_liquidus.tight_layout()
		self.canvas_liquidus.draw()
	
	def _display_liquidus_data (self):
		"""显示液相线/固相线数据表格"""
		self.data_text.delete('1.0', tk.END)

		if not self.results_data or 'x_scan' not in self.results_data:
			self.data_text.insert(tk.END, "没有液相线/固相线数据\n")
			return

		x_scan = self.results_data['x_scan']
		scan_comp = self.results_data['scan_comp']
		header_parts = [f"X({scan_comp})"]
		model_keys_with_data = []

		# 确定有数据的模型
		for model_key in ['RKM', 'UEM1', 'Muggianu', 'Toop']:
			if (model_key in self.results_data and
					self.results_data[model_key]['type'] == 'liquidus_solidus'):
				model_keys_with_data.append(model_key)
				label = self.results_data[model_key]['label']
				header_parts.append(f"T_liq_{label}(K)")
				header_parts.append(f"T_sol_{label}(K)")

		self.data_text.insert(tk.END, "\t".join(header_parts) + "\n")
		self.data_text.insert(tk.END, "=" * 100 + "\n")

		# 数据行
		for i, x_val in enumerate(x_scan):
			line_parts = [f"{x_val:.4f}"]
			for model_key in model_keys_with_data:
				T_liq = self.results_data[model_key]['liquidus'][i]
				T_sol = self.results_data[model_key]['solidus'][i]
				line_parts.append(
						f"{T_liq:.2f}" if not np.isnan(T_liq) else "N/A")
				line_parts.append(
						f"{T_sol:.2f}" if not np.isnan(T_sol) else "N/A")
			self.data_text.insert(tk.END, "\t".join(line_parts) + "\n")
	
	# =========================================================================
	# 模块2: 热力学性质计算
	# =========================================================================
	def calculate_properties (self):
		"""
		计算热力学性质（主入口）

		功能：计算以下热力学性质
		1. Gibbs自由能 (GM)
		2. 化学势 (MU)
		3. 活度 (Activity = exp((μ-μ_ref)/RT))
		"""
		inputs = self._parse_inputs()
		if not inputs:
			return
		
		selected_model_keys = self.get_selected_models()
		if not selected_model_keys:
			return
		
		thread = threading.Thread(
				target=self._calculate_properties_thread,
				args=(selected_model_keys, inputs))
		thread.daemon = True
		thread.start()
	
	def _calculate_properties_thread (self, selected_model_keys, inputs):
		"""
		热力学性质计算线程

		计算步骤：
		1. 计算参考态化学势（纯组元在稳定相）
		2. 对每个成分点计算平衡
		3. 提取Gibbs自由能
		4. 提取化学势并计算活度
		"""
		try:
			temp_calc = inputs['temp_max']  # 使用最高温度
			self.log(f"\n{'=' * 60}")
			self.log("[热力学性质计算]")
			self.log(f"{'=' * 60}")
			
			self.log(f"计算温度: {temp_calc} K")
			liquid_phase_name = None
			for phase in inputs['db_phases']:
				if 'LIQUID' in phase.upper():
					liquid_phase_name = phase
					break
			
			if liquid_phase_name is None:
				self.log("错误: 数据库中未找到 LIQUID 相")
				raise ValueError("数据库中未找到 LIQUID 相")
			
			self.log(f"将计算 {liquid_phase_name} 相的性质")
			
			
			# --- 步骤1: 计算参考态化学势 ---
			ref_mus = self._calculate_reference_potentials(
					inputs['study_comps'], temp_calc, liquid_phase_name)
			
			RT = float(v.R) * temp_calc
			
			# 清除旧的性质数据
			keys_to_remove = [k for k in self.results_data if '_props' in k]
			for k in keys_to_remove:
				del self.results_data[k]
			
			# --- 步骤2: 对每个模型计算性质 ---
			ratio_sum = sum(inputs['ratios']) if inputs['ratios'] else 1.0
			
			for model_key in selected_model_keys:
				model_label = self.model_labels[model_key]
				if model_key == 'UEM1' and self.uem1_liquid_only.get():
					model_label += " (Liq Only)"
				
				self.progress_var.set(f"计算性质: {model_label}")
				self.log(f"\n计算 {model_label} 热力学性质 @ {temp_calc}K")
				
				model_spec = self.get_model_spec(model_key)
				if model_spec is None and model_key != 'RKM':
					continue
				
				# 准备数据存储
				active_comps_no_va = [c for c in inputs['study_comps']
				                      if c != 'VA']
				gibbs_list = []
				activity_data = {comp: [] for comp in active_comps_no_va}
				
				# --- 步骤3: 遍历成分点 ---
				for x_scan in inputs['x_scan_range']:
					composition = self._calculate_composition(
							x_scan, inputs, ratio_sum)
					
					# 设置平衡条件
					conditions = {v.T: temp_calc, v.P: 101325}
					comps_to_set = [c for c in inputs['study_comps']
					                if c != 'VA'][:-1]
					for comp in comps_to_set:
						if comp in composition:
							conditions[v.X(comp)] = composition[comp]
					
					try:
						# 计算平衡
						eq = equilibrium(
								self.dbe, inputs['study_comps'],
								[liquid_phase_name],  # <--- 修正: 强制只计算液相
								model=model_spec, conditions=conditions,
								calc_opts={'pdens': 50})
						
						# 提取Gibbs自由能
						G = float(eq.GM.squeeze().item())
						gibbs_list.append(G)
						
						# 提取化学势并计算活度
						for comp in activity_data.keys():
							if comp in ref_mus and not np.isnan(ref_mus[comp]):
								try:
									# 安全提取化学势
									mu_data = eq.MU.sel(component=comp)
									if mu_data.size == 0:
										raise ValueError(f"{comp} 无 MU 数据")

									mu_mix = float(mu_data.values.flatten()[0])

									# 验证有效性
									if np.isnan(mu_mix) or np.isinf(mu_mix):
										raise ValueError(f"{comp} MU 值无效")

									# 计算活度
									activity = np.exp((mu_mix - ref_mus[comp]) / RT)

									# 验证活度值
									if np.isnan(activity) or np.isinf(activity) or activity < 0:
										raise ValueError(f"{comp} 活度值无效: {activity}")

									activity_data[comp].append(activity)

								except (KeyError, Exception) as e:
									activity_data[comp].append(np.nan)
									# 详细错误日志（但不要太频繁）
									if x_scan < 0.15 or x_scan > 0.85:  # 只在范围边界记录
										self.log(f"    {comp} 活度计算失败 (X={x_scan:.3f}): {e}")
							else:
								activity_data[comp].append(np.nan)
					
					except Exception as e:
						gibbs_list.append(np.nan)
						for comp in activity_data.keys():
							activity_data[comp].append(np.nan)
						self.log(
								f"  性质计算失败 "
								f"(X({inputs['scan_comp']})={x_scan:.3f}): {e}")
				
				# 存储结果
				self.results_data[f'{model_key}_props'] = {
					'x_scan': inputs['x_scan_range'],
					'scan_comp': inputs['scan_comp'],
					'gibbs': np.array(gibbs_list),
					'activity': {k: np.array(v) for k, v in activity_data.items()},
					'temperature': temp_calc,
					'type': 'properties',
					'label': model_label
				}
				self.log(f"{model_label} 性质计算完成！")
			
			# 绘图和显示
			self._plot_gibbs()
			self._plot_activity()
			self._display_properties_data()
			
			self.progress_var.set("性质计算完成！")
			self.log("\n所有热力学性质计算完成！")
		
		except Exception as e:
			self.log(f"性质计算错误: {e}")
			import traceback
			self.log(traceback.format_exc())
			self.progress_var.set("计算失败")
	
	def _calculate_reference_potentials (self, study_comps, temperature, liquid_phase_name):
		"""
		计算参考态化学势（纯组元在稳定相）

		参数:
			study_comps: 研究组分列表
			temperature: 温度(K)

		返回:
			dict: {组分: 参考态化学势(J/mol)}
		"""
		self.log("计算参考态化学势...")
		ref_mus = {}

		active_comps = [c for c in study_comps if c != 'VA']

		for comp in active_comps:
			mu_ref = None
			method_used = None

			# 方法1: 尝试使用液相（必须用 equilibrium 才能得到 MU）
			try:
				ref_eq = equilibrium(
						self.dbe, [comp, 'VA'],
						liquid_phase_name,
						conditions={v.T: temperature, v.P: 101325, v.N: 1},
						model=Model)

				# 安全提取化学势
				mu_data = ref_eq.MU.sel(component=comp)
				if mu_data.size == 0:
					raise ValueError(f"{comp} 在 {liquid_phase_name} 中无 MU 数据")

				mu_ref = float(mu_data.values.flatten()[0])

				# 验证有效性
				if np.isnan(mu_ref) or np.isinf(mu_ref):
					raise ValueError(f"{comp} MU 值无效 (NaN/Inf)")

				method_used = f"liquid({liquid_phase_name})"

			except Exception as e:
				self.log(f"  {comp} 液相参考态失败: {e}")

				# 方法2: 使用稳定相作为备用方案
				try:
					self.log(f"  尝试使用 {comp} 的稳定相...")
					stable_eq = equilibrium(
							self.dbe, [comp, 'VA'],
							list(self.dbe.phases.keys()),
							conditions={v.T: temperature, v.P: 101325, v.N: 1})

					mu_data = stable_eq.MU.sel(component=comp)
					if mu_data.size == 0:
						raise ValueError(f"{comp} 稳定相中无 MU 数据")

					mu_ref = float(mu_data.values.flatten()[0])

					if np.isnan(mu_ref) or np.isinf(mu_ref):
						raise ValueError(f"{comp} MU 值无效 (NaN/Inf)")

					# 获取稳定相名称
					stable_phases = stable_eq.Phase.values.flatten()
					stable_phase = stable_phases[0] if len(stable_phases) > 0 else 'unknown'
					method_used = f"stable({stable_phase})"

				except Exception as e2:
					self.log(f"  {comp} 稳定相计算也失败: {e2}")
					mu_ref = None

			# 存储结果
			if mu_ref is not None:
				ref_mus[comp] = mu_ref
				self.log(f"  ✓ MU_ref({comp}) = {mu_ref:.2f} J/mol  [方法: {method_used}]")
			else:
				ref_mus[comp] = np.nan
				self.log(f"  ✗ {comp} 所有方法均失败，设置为 NaN")

		return ref_mus
	
	def _plot_gibbs (self):
		"""绘制Gibbs自由能图"""
		self.ax_gibbs.clear()
		
		colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
		markers = ['o', 's', '^', 'D']
		idx = 0
		plot_title = "Gibbs Free Energy"
		temp = None
		scan_comp = None
		plotted_something = False
		
		# 固定顺序绘制
		for model_key in ['RKM', 'UEM1', 'Muggianu', 'Toop']:
			prop_key = f'{model_key}_props'
			if (prop_key in self.results_data and
					self.results_data[prop_key]['type'] == 'properties'):
				
				data = self.results_data[prop_key]
				x_scan = data['x_scan']
				gibbs = data['gibbs']
				scan_comp = data['scan_comp']
				temp = data['temperature']
				label = data['label']
				valid = ~np.isnan(gibbs)
				
				if np.any(valid):
					self.ax_gibbs.plot(
							x_scan[valid], gibbs[valid] / 1000.0,  # 转kJ/mol
							marker=markers[idx], linestyle='-',
							linewidth=2, markersize=5,
							label=label, color=colors[idx], alpha=0.9)
					plotted_something = True
				idx += 1
		
		if plotted_something and temp is not None and scan_comp is not None:
			self.ax_gibbs.set_xlabel(f'X({scan_comp})', fontsize=12)
			self.ax_gibbs.set_ylabel('Gibbs Energy (kJ/mol)', fontsize=12)
			plot_title = f'Gibbs Free Energy at {temp:.0f} K'
			self.ax_gibbs.legend(fontsize=10, loc='best')
			self.ax_gibbs.grid(True, alpha=0.3, linestyle='--')
		else:
			self.ax_gibbs.text(
					0.5, 0.5, 'No valid data to plot',
					ha='center', va='center')
		
		self.ax_gibbs.set_title(plot_title, fontsize=14, fontweight='bold')
		self.fig_gibbs.tight_layout()
		self.canvas_gibbs.draw()
	
	def _plot_activity (self):
		"""绘制活度图"""
		self.ax_activity.clear()
		
		colors = plt.cm.tab10.colors
		linestyles = ['-', '--', '-.', ':']
		plot_title = "Activity"
		temp = None
		scan_comp = None
		plotted_something = False
		
		model_idx = 0
		# 固定顺序绘制
		for model_key in ['RKM', 'UEM1', 'Muggianu', 'Toop']:
			prop_key = f'{model_key}_props'
			if (prop_key in self.results_data and
					self.results_data[prop_key]['type'] == 'properties'):
				
				data = self.results_data[prop_key]
				x_scan = data['x_scan']
				activity_data = data['activity']
				scan_comp = data['scan_comp']
				temp = data['temperature']
				model_label_short = model_key
				
				comp_idx = 0
				for comp, acts_array in activity_data.items():
					valid = ~np.isnan(acts_array)
					if np.any(valid):
						self.ax_activity.plot(
								x_scan[valid], acts_array[valid],
								linestyle=linestyles[model_idx % len(linestyles)],
								linewidth=1.5, marker=None,
								label=f'$a_{{{comp}}}$ ({model_label_short})',
								color=colors[comp_idx % len(colors)])
						plotted_something = True
					comp_idx += 1
				model_idx += 1
		
		if plotted_something and temp is not None and scan_comp is not None:
			self.ax_activity.set_xlabel(f'X({scan_comp})', fontsize=12)
			self.ax_activity.set_ylabel('Activity', fontsize=12)
			plot_title = f'Activity at {temp:.0f} K'
			self.ax_activity.legend(fontsize=9, loc='best', ncol=2)
			self.ax_activity.grid(True, alpha=0.3, linestyle='--')
			self.ax_activity.set_ylim(bottom=0)
		else:
			# 显示详细的错误信息
			error_msg = "❌ 无有效活度数据\n\n"
			error_msg += "可能原因:\n"
			error_msg += "• 参考态化学势计算失败\n"
			error_msg += "• 混合相化学势提取失败\n"
			error_msg += "• 所有数据点计算失败\n\n"
			error_msg += "请查看 [日志] 标签页获取详细信息"

			self.ax_activity.text(
					0.5, 0.5, error_msg,
					ha='center', va='center', fontsize=10, color='#d62728',
					bbox=dict(boxstyle='round', facecolor='#ffe6e6', alpha=0.8, edgecolor='#d62728'))
		
		self.ax_activity.set_title(plot_title, fontsize=14, fontweight='bold')
		self.fig_activity.tight_layout()
		self.canvas_activity.draw()
	
	def _display_properties_data (self):
		"""显示热力学性质数据表格"""
		self.data_text.delete('1.0', tk.END)
		
		data_available = False
		scan_comp = None
		x_scan = None
		temp = None
		model_keys_with_data = []
		all_comps = set()
		
		# 收集数据
		for model_key in ['RKM', 'UEM1', 'Muggianu', 'Toop']:
			key = f'{model_key}_props'
			if (key in self.results_data and
					self.results_data[key]['type'] == 'properties'):
				data_available = True
				model_keys_with_data.append(model_key)
				if scan_comp is None:
					scan_comp = self.results_data[key]['scan_comp']
					x_scan = self.results_data[key]['x_scan']
					temp = self.results_data[key]['temperature']
				all_comps.update(self.results_data[key]['activity'].keys())
		
		if not data_available:
			self.data_text.insert(tk.END, "没有热力学性质数据\n")
			return
		
		sorted_comps = sorted(list(all_comps))
		
		# 表头
		header_parts = [f"T={temp:.0f}K", f"X({scan_comp})"]
		for model_key in model_keys_with_data:
			label = self.results_data[f'{model_key}_props']['label']
			header_parts.append(f"GM_{label}(kJ/mol)")
		for comp in sorted_comps:
			for model_key in model_keys_with_data:
				label = self.results_data[f'{model_key}_props']['label']
				header_parts.append(f"Act_{comp}_{label}")
		
		header = "\t".join(header_parts) + "\n"
		self.data_text.insert(tk.END, header)
		self.data_text.insert(tk.END, "=" * 120 + "\n")
		
		# 数据行
		for i, x_val in enumerate(x_scan):
			line_parts = ["", f"{x_val:.4f}"]
			
			# Gibbs自由能
			for model_key in model_keys_with_data:
				key = f'{model_key}_props'
				G = self.results_data[key]['gibbs'][i]
				line_parts.append(
						f"{G / 1000.0:.3f}" if not np.isnan(G) else "N/A")
			
			# 活度
			for comp in sorted_comps:
				for model_key in model_keys_with_data:
					key = f'{model_key}_props'
					if comp in self.results_data[key]['activity']:
						Act = self.results_data[key]['activity'][comp][i]
						line_parts.append(
								f"{Act:.4f}" if not np.isnan(Act) else "N/A")
					else:
						line_parts.append("-")
			
			self.data_text.insert(tk.END, "\t".join(line_parts) + "\n")
	
	# =========================================================================
	# 模块3: 伪二元相图计算
	# =========================================================================
	def calculate_phase_diagram (self):
		"""
		计算伪二元相图（主入口）

		说明：
		伪二元相图 - 固定其他组分比例，只改变一个扫描组分的摩尔分数
		例如：三元系Al-Cr-Ni，固定Al:Cr=1:1，扫描Ni的摩尔分数
		这实际上是在三元成分空间中的一条直线上的相图投影

		计算内容：
		- 温度范围：temp_min ~ temp_max
		- 成分范围：扫描组分的摩尔分数范围
		- 每个点计算平衡相
		"""
		inputs = self._parse_inputs()
		if not inputs:
			return
		
		selected_model_keys = self.get_selected_models()
		if not selected_model_keys:
			return
		
		# 只使用第一个选中的模型
		model_key = selected_model_keys[0]
		
		thread = threading.Thread(
				target=self._calculate_phase_diagram_thread,
				args=(model_key, inputs))
		thread.daemon = True
		thread.start()
	
	def _calculate_phase_diagram_thread (self, model_key, inputs):
		"""
		伪二元相图计算线程

		计算逻辑：
		1. 在温度-成分网格上遍历所有点
		2. 对每个点计算平衡
		3. 识别存在的相（分数>阈值）
		4. 绘制相区图
		"""
		try:
			model_label = self.model_labels[model_key]
			if model_key == 'UEM1' and self.uem1_liquid_only.get():
				model_label += " (Liq Only)"
			
			self.log(f"\n{'=' * 60}")
			self.log(f"[伪二元相图计算] 模型: {model_label}")
			self.log(f"{'=' * 60}")
			self.log(f"说明: 固定其他组分比例，扫描{inputs['scan_comp']}组分")
			self.progress_var.set(f"生成相图: {model_label}")
			
			# 温度和成分网格
			temp_num = 50  # 温度点数
			T_range = np.linspace(inputs['temp_min'], inputs['temp_max'], temp_num)
			
			model_spec = self.get_model_spec(model_key)
			if model_spec is None and model_key != 'RKM':
				raise ValueError(f"模型 {model_label} 不可用")
			
			# 相图数据存储：{(T, X): 相名称字符串}
			phase_data_dict = {}
			
			total_calcs = len(inputs['x_scan_range']) * len(T_range)
			calc_count = 0
			ratio_sum = sum(inputs['ratios']) if inputs['ratios'] else 1.0
			
			# 遍历成分
			for i, x_scan in enumerate(inputs['x_scan_range']):
				# 计算当前成分
				composition = self._calculate_composition(
						x_scan, inputs, ratio_sum)
				
				# 遍历温度
				for j, T in enumerate(T_range):
					conditions = {v.T: T, v.P: 101325}
					comps_to_set = [c for c in inputs['study_comps']
					                if c != 'VA'][:-1]
					for comp in comps_to_set:
						if comp in composition:
							conditions[v.X(comp)] = composition[comp]
					
					phase_string = 'ERROR'
					try:
						eq = equilibrium(
								self.dbe, inputs['study_comps'], inputs['db_phases'],
								model=model_spec, conditions=conditions,
								calc_opts={'pdens': 50})
						
						# 提取存在的相
						present_phases = eq.Phase.values.squeeze()
						phase_fracs = eq.NP.values.squeeze()
						
						# 处理单点结果
						if present_phases.ndim == 0:
							present_phases = np.array([present_phases.item()])
							phase_fracs = np.array([phase_fracs.item()])
						
						# 筛选显著相（分数>1e-5）
						valid_indices = (present_phases != '') & (phase_fracs > 1e-5)
						phases_found = present_phases[valid_indices]
						
						# 解码bytes
						phases_found_str = []
						for ph in phases_found:
							if isinstance(ph, bytes):
								phases_found_str.append(ph.decode('utf-8'))
							elif isinstance(ph, str):
								phases_found_str.append(ph)
						
						# 去重、排序、用'+'连接
						phase_string = '+'.join(sorted(list(set(phases_found_str))))
						if not phase_string:
							phase_string = 'NoPhase?'
					
					except Exception as e:
						phase_string = 'ERROR'
					
					phase_data_dict[(T, x_scan)] = phase_string
					
					calc_count += 1
					if calc_count % 50 == 0:
						progress = calc_count / total_calcs * 100
						self.progress_var.set(f"相图计算: {progress:.1f}%")
			
			# 绘制相图
			self._plot_phase_diagram(
					inputs['x_scan_range'], T_range, phase_data_dict,
					inputs['scan_comp'], model_label)
			
			# 保存相图
			timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
			safe_model_label = re.sub(r'[^\w\-]+', '_', model_label)
			safe_comps = "-".join([c for c in inputs['study_comps'] if c != 'VA'])
			filename = (f"phase_diagram_{safe_model_label}_{safe_comps}_"
			            f"{inputs['scan_comp']}_{timestamp}.png")
			
			self.fig_phase.savefig(filename, dpi=300, bbox_inches='tight')
			self.log(f"\n相图已保存: {filename}")
			self.progress_var.set("相图生成完成！")
		
		except Exception as e:
			self.log(f"相图计算错误: {e}")
			import traceback
			self.log(traceback.format_exc())
			self.progress_var.set("计算失败")
	
	def _plot_phase_diagram (self, x_range, T_range, phase_data_dict,
	                         scan_comp, model_label):
		"""
		绘制伪二元相图

		使用imshow和colorbar显示不同相区
		"""
		# 清除旧图
		self.fig_phase.clear()
		self.ax_phase = self.fig_phase.add_subplot(111)
		
		# 构建数值网格
		phase_numeric_grid = np.zeros((len(T_range), len(x_range)))
		unique_phase_strings = sorted(
				list(set(phase_data_dict.values()) - {'ERROR'}))
		phase_map = {name: i for i, name in enumerate(unique_phase_strings)}
		phase_map['ERROR'] = -1
		
		for j, T in enumerate(T_range):
			for i, x_val in enumerate(x_range):
				phase_string = phase_data_dict.get((T, x_val), 'ERROR')
				phase_numeric_grid[j, i] = phase_map.get(phase_string, -1)
		
		# 绘图
		cmap = plt.cm.get_cmap('tab20', len(unique_phase_strings))
		bounds = np.arange(len(unique_phase_strings) + 1) - 0.5
		
		im = self.ax_phase.imshow(
				phase_numeric_grid,
				extent=[x_range.min(), x_range.max(),
				        T_range.min(), T_range.max()],
				origin='lower', aspect='auto', cmap=cmap,
				interpolation='nearest',
				vmin=-0.5, vmax=len(unique_phase_strings) - 0.5)
		
		self.ax_phase.set_xlabel(f'X({scan_comp})', fontsize=12)
		self.ax_phase.set_ylabel('Temperature (K)', fontsize=12)
		self.ax_phase.set_title(
				f'Pseudo-Binary Phase Diagram ({model_label})',
				fontsize=14, fontweight='bold')
		
		# 添加颜色条
		cbar = self.fig_phase.colorbar(
				im, ax=self.ax_phase,
				boundaries=bounds,
				ticks=np.arange(len(unique_phase_strings)))
		cbar.ax.set_yticklabels(unique_phase_strings)
		cbar.set_label('Phase(s)', rotation=270, labelpad=15)
		
		self.fig_phase.tight_layout()
		self.canvas_phase.draw()
	
	# =========================================================================
	# 辅助功能
	# =========================================================================
	def export_data (self):
		"""导出数据到文件"""
		if not self.data_text.get('1.0', tk.END).strip():
			messagebox.showwarning("警告", "没有可导出的数据！")
			return
		
		file_path = filedialog.asksaveasfilename(
				title="保存数据",
				defaultextension=".txt",
				filetypes=[("制表符分隔文本", "*.txt"),
				           ("CSV文件", "*.csv"),
				           ("所有文件", "*.*")]
		)
		
		if file_path:
			try:
				content = self.data_text.get('1.0', tk.END)
				if file_path.lower().endswith(".csv"):
					content = content.replace('\t', ',')
				
				with open(file_path, 'w', encoding='utf-8') as f:
					f.write(content)
				
				messagebox.showinfo("成功", f"数据已导出到:\n{file_path}")
				self.log(f"数据已导出: {file_path}")
			
			except Exception as e:
				messagebox.showerror("错误", f"导出失败:\n{e}")
				self.log(f"导出数据失败: {e}")
	
	def clear_results (self):
		"""清除所有结果"""
		self.ax_liquidus.clear()
		self.ax_gibbs.clear()
		self.ax_activity.clear()
		
		self.fig_phase.clear()
		self.ax_phase = self.fig_phase.add_subplot(111)
		
		self.canvas_liquidus.draw()
		self.canvas_gibbs.draw()
		self.canvas_activity.draw()
		self.canvas_phase.draw()
		
		self.data_text.delete('1.0', tk.END)
		self.results_data = {}
		
		self.progress_var.set("就绪")
		self.log("\n结果已清除")


# =============================================================================
# 主程序入口
# =============================================================================
def main ():
	"""主函数"""
	root = tk.Tk()
	app = AlloyCalculatorGUI(root)
	root.mainloop()


if __name__ == "__main__":
	main()