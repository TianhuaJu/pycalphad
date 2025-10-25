#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
合金热力学计算与对比工具 - 简化版GUI
支持多模型对比、自动相识别、热力学性质计算、相图生成
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt
from pycalphad import Database, equilibrium, variables as v
from pycalphad.core.utils import instantiate_models
from pycalphad import Model
import threading
from datetime import datetime
import os


class AlloyCalculatorGUI:
    """合金热力学计算GUI - 简化版"""

    def __init__(self, root):
        self.root = root
        self.root.title("合金热力学计算与模型对比工具")
        self.root.geometry("1400x900")

        self.dbe = None
        self.available_phases = []
        self.available_comps = []
        self.results_data = {}  # 存储多个模型的结果

        # 可用模型
        self.available_models = {
            'RKM': Model,
            'Muggianu': None,
            'Toop': None,
            'UEM1': None
        }

        self.create_widgets()

    def create_widgets(self):
        """创建界面组件"""

        # 左侧控制面板
        left_frame = ttk.Frame(self.root, width=400)
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=False, padx=5, pady=5)

        # 右侧结果显示区
        right_frame = ttk.Frame(self.root)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)

        # === 左侧面板内容 ===

        # 1. 数据库加载
        db_frame = ttk.LabelFrame(left_frame, text="1. 数据库", padding=10)
        db_frame.pack(fill=tk.X, pady=5)

        ttk.Button(db_frame, text="加载TDB数据库", command=self.load_database).pack(fill=tk.X)
        self.db_label = ttk.Label(db_frame, text="未加载", foreground="red")
        self.db_label.pack(pady=5)

        # 显示可用组分和相
        ttk.Label(db_frame, text="数据库可用组分:", font=('', 8, 'bold')).pack(anchor=tk.W)
        self.comps_label = ttk.Label(db_frame, text="", foreground="blue", wraplength=350, font=('', 8))
        self.comps_label.pack(anchor=tk.W)

        ttk.Label(db_frame, text="数据库可用相:", font=('', 8, 'bold')).pack(anchor=tk.W)
        self.phases_label = ttk.Label(db_frame, text="", foreground="blue", wraplength=350, font=('', 8))
        self.phases_label.pack(anchor=tk.W)

        # 2. 组分选择
        comp_select_frame = ttk.LabelFrame(left_frame, text="2. 组分选择", padding=10)
        comp_select_frame.pack(fill=tk.X, pady=5)

        ttk.Label(comp_select_frame, text="研究组分列表 (逗号分隔):").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.comps_entry = ttk.Entry(comp_select_frame, width=20)
        self.comps_entry.grid(row=0, column=1, sticky=tk.W+tk.E, pady=2)
        self.comps_entry.insert(0, "AL,CU,Y")

        ttk.Label(comp_select_frame, text="说明: 可以指定子系统", font=('', 8), foreground="gray").grid(
            row=1, column=0, columnspan=2, sticky=tk.W, pady=2)
        ttk.Label(comp_select_frame, text="例如: AL,CU (二元) 或 AL,CU,Y (三元)", font=('', 8), foreground="gray").grid(
            row=2, column=0, columnspan=2, sticky=tk.W)

        # 3. 成分设置
        comp_frame = ttk.LabelFrame(left_frame, text="3. 成分扫描", padding=10)
        comp_frame.pack(fill=tk.X, pady=5)

        ttk.Label(comp_frame, text="扫描组分:").grid(row=0, column=0, sticky=tk.W)
        self.scan_comp_entry = ttk.Entry(comp_frame, width=15)
        self.scan_comp_entry.grid(row=0, column=1, sticky=tk.W)
        self.scan_comp_entry.insert(0, "Y")

        ttk.Label(comp_frame, text="扫描范围 (起,止,点数):").grid(row=1, column=0, sticky=tk.W)
        self.scan_range_entry = ttk.Entry(comp_frame, width=15)
        self.scan_range_entry.grid(row=1, column=1, sticky=tk.W)
        self.scan_range_entry.insert(0, "0.1,0.9,10")

        ttk.Label(comp_frame, text="其他组分比例:").grid(row=2, column=0, sticky=tk.W)
        self.other_ratio_entry = ttk.Entry(comp_frame, width=15)
        self.other_ratio_entry.grid(row=2, column=1, sticky=tk.W)
        self.other_ratio_entry.insert(0, "1:1")

        # 4. 温度设置
        temp_frame = ttk.LabelFrame(left_frame, text="4. 温度范围 (K)", padding=10)
        temp_frame.pack(fill=tk.X, pady=5)

        ttk.Label(temp_frame, text="最低温度:").grid(row=0, column=0, sticky=tk.W)
        self.temp_min_entry = ttk.Entry(temp_frame, width=10)
        self.temp_min_entry.grid(row=0, column=1, sticky=tk.W)
        self.temp_min_entry.insert(0, "300")

        ttk.Label(temp_frame, text="最高温度:").grid(row=1, column=0, sticky=tk.W)
        self.temp_max_entry = ttk.Entry(temp_frame, width=10)
        self.temp_max_entry.grid(row=1, column=1, sticky=tk.W)
        self.temp_max_entry.insert(0, "2000")

        ttk.Label(temp_frame, text="温度步长:").grid(row=2, column=0, sticky=tk.W)
        self.temp_step_entry = ttk.Entry(temp_frame, width=10)
        self.temp_step_entry.grid(row=2, column=1, sticky=tk.W)
        self.temp_step_entry.insert(0, "10")

        # 5. 模型选择（多选）
        model_frame = ttk.LabelFrame(left_frame, text="5. 模型选择（可多选对比）", padding=10)
        model_frame.pack(fill=tk.X, pady=5)

        self.model_vars = {}
        for model_name in ['RKM', 'Muggianu', 'Toop', 'UEM1']:
            var = tk.BooleanVar(value=(model_name == 'RKM'))
            self.model_vars[model_name] = var
            ttk.Checkbutton(model_frame, text=model_name, variable=var).pack(anchor=tk.W)

        # 6. 计算控制
        control_frame = ttk.LabelFrame(left_frame, text="6. 计算控制", padding=10)
        control_frame.pack(fill=tk.X, pady=5)

        ttk.Button(control_frame, text="计算液相线对比",
                  command=self.calculate_liquidus).pack(fill=tk.X, pady=2)

        ttk.Button(control_frame, text="计算热力学性质",
                  command=self.calculate_properties).pack(fill=tk.X, pady=2)

        ttk.Button(control_frame, text="生成相图",
                  command=self.calculate_phase_diagram).pack(fill=tk.X, pady=2)

        ttk.Button(control_frame, text="清除结果",
                  command=self.clear_results).pack(fill=tk.X, pady=2)

        # 进度显示
        self.progress_var = tk.StringVar(value="就绪")
        ttk.Label(control_frame, textvariable=self.progress_var,
                 foreground="green").pack(pady=5)

        # === 右侧结果显示 ===

        # 创建Notebook（多标签页）
        self.notebook = ttk.Notebook(right_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # 标签页1: 液相线对比图
        liquidus_tab = ttk.Frame(self.notebook)
        self.notebook.add(liquidus_tab, text="液相线对比")

        self.fig_liquidus = Figure(figsize=(8, 6))
        self.ax_liquidus = self.fig_liquidus.add_subplot(111)
        self.canvas_liquidus = FigureCanvasTkAgg(self.fig_liquidus, liquidus_tab)
        self.canvas_liquidus.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar_liquidus = NavigationToolbar2Tk(self.canvas_liquidus, liquidus_tab)
        toolbar_liquidus.update()

        # 标签页2: 热力学性质
        props_tab = ttk.Frame(self.notebook)
        self.notebook.add(props_tab, text="热力学性质")

        # 创建子标签页
        props_notebook = ttk.Notebook(props_tab)
        props_notebook.pack(fill=tk.BOTH, expand=True)

        # 子标签: 自由能
        gibbs_frame = ttk.Frame(props_notebook)
        props_notebook.add(gibbs_frame, text="自由能")
        self.fig_gibbs = Figure(figsize=(8, 5))
        self.ax_gibbs = self.fig_gibbs.add_subplot(111)
        self.canvas_gibbs = FigureCanvasTkAgg(self.fig_gibbs, gibbs_frame)
        self.canvas_gibbs.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # 子标签: 活度
        activity_frame = ttk.Frame(props_notebook)
        props_notebook.add(activity_frame, text="活度")
        self.fig_activity = Figure(figsize=(8, 5))
        self.ax_activity = self.fig_activity.add_subplot(111)
        self.canvas_activity = FigureCanvasTkAgg(self.fig_activity, activity_frame)
        self.canvas_activity.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # 标签页3: 相图
        phase_diagram_tab = ttk.Frame(self.notebook)
        self.notebook.add(phase_diagram_tab, text="相图")

        self.fig_phase = Figure(figsize=(8, 6))
        self.ax_phase = self.fig_phase.add_subplot(111)
        self.canvas_phase = FigureCanvasTkAgg(self.fig_phase, phase_diagram_tab)
        self.canvas_phase.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # 标签页4: 数据表格
        data_tab = ttk.Frame(self.notebook)
        self.notebook.add(data_tab, text="数据")

        self.data_text = scrolledtext.ScrolledText(data_tab, wrap=tk.WORD,
                                                   font=('Courier', 10))
        self.data_text.pack(fill=tk.BOTH, expand=True)

        ttk.Button(data_tab, text="导出数据", command=self.export_data).pack(pady=5)

        # 标签页5: 日志
        log_tab = ttk.Frame(self.notebook)
        self.notebook.add(log_tab, text="日志")

        self.log_text = scrolledtext.ScrolledText(log_tab, wrap=tk.WORD,
                                                  font=('Courier', 9))
        self.log_text.pack(fill=tk.BOTH, expand=True)

    def log(self, message):
        """添加日志"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_text.insert(tk.END, f"[{timestamp}] {message}\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()

    def load_database(self):
        """加载数据库"""
        file_path = filedialog.askopenfilename(
            title="选择TDB数据库文件",
            filetypes=[("TDB文件", "*.tdb"), ("所有文件", "*.*")]
        )

        if file_path:
            try:
                self.dbe = Database(file_path)
                self.db_label.config(text=f"已加载: {os.path.basename(file_path)}",
                                    foreground="green")

                # 提取可用组分和相
                # 过滤掉非正常元素符号（只保留1-2个大写字母的元素）
                import re
                all_elements = [str(c) for c in self.dbe.elements if c != 'VA']
                self.available_comps = sorted([
                    c for c in all_elements
                    if re.match(r'^[A-Z][A-Z]?$', c)  # 只保留正常元素符号
                ])
                self.available_phases = sorted(list(self.dbe.phases.keys()))

                self.comps_label.config(text=f"组分: {', '.join(self.available_comps)}")
                self.phases_label.config(text=f"相: {', '.join(self.available_phases)}")

                self.log(f"成功加载数据库: {file_path}")
                self.log(f"组分: {self.available_comps}")
                self.log(f"相: {self.available_phases}")

                # 动态导入高级模型
                try:
                    from pycalphad.advanced_uem_model import ModelMuggianu, ModelToop
                    from pycalphad.uem1_Model import uem1_model

                    self.available_models['Muggianu'] = ModelMuggianu
                    self.available_models['Toop'] = ModelToop
                    self.available_models['UEM1'] = uem1_model

                    self.log("成功加载所有模型: RKM, Muggianu, Toop, UEM1")

                except ImportError as e:
                    self.log(f"警告: 部分高级模型不可用 - {e}")
                    self.log("仅RKM模型可用")

            except Exception as e:
                messagebox.showerror("错误", f"加载数据库失败:\n{e}")
                self.log(f"数据库加载失败: {e}")

    def get_selected_models(self):
        """获取选中的模型"""
        selected = []
        for model_name, var in self.model_vars.items():
            if var.get():
                if self.available_models[model_name] is None:
                    self.log(f"警告: {model_name} 模型未加载，跳过")
                else:
                    selected.append(model_name)
        return selected

    def calculate_liquidus(self):
        """计算液相线（多模型对比）"""
        if self.dbe is None:
            messagebox.showwarning("警告", "请先加载数据库！")
            return

        selected_models = self.get_selected_models()
        if not selected_models:
            messagebox.showwarning("警告", "请至少选择一个模型！")
            return

        # 在后台线程中运行
        thread = threading.Thread(target=self._calculate_liquidus_thread,
                                 args=(selected_models,))
        thread.daemon = True
        thread.start()

    def _calculate_liquidus_thread(self, selected_models):
        """液相线计算线程"""
        try:
            # 解析用户指定的研究组分
            comps_str = self.comps_entry.get().strip().upper()
            study_comps = [c.strip() for c in comps_str.split(',') if c.strip()]

            # 验证组分是否都在数据库中
            invalid_comps = [c for c in study_comps if c not in self.available_comps]
            if invalid_comps:
                raise ValueError(
                    f"以下组分不在数据库中: {invalid_comps}\n"
                    f"数据库可用组分: {self.available_comps}"
                )

            if len(study_comps) < 2:
                raise ValueError("至少需要2个组分！")

            self.log(f"研究组分: {study_comps}")

            # 解析其他输入
            scan_comp = self.scan_comp_entry.get().strip().upper()
            scan_range_str = self.scan_range_entry.get().strip()
            start, stop, num = [float(x) for x in scan_range_str.split(',')]

            temp_min = float(self.temp_min_entry.get())
            temp_max = float(self.temp_max_entry.get())
            temp_step = float(self.temp_step_entry.get())

            other_ratio_str = self.other_ratio_entry.get().strip()

            x_scan_range = np.linspace(start, stop, int(num))

            # 检查扫描组分是否有效
            if scan_comp not in study_comps:
                raise ValueError(f"扫描组分 '{scan_comp}' 不在研究组分列表中: {study_comps}")

            # 检查比例数量
            other_comps = [c for c in study_comps if c != scan_comp]
            ratios = [float(r) for r in other_ratio_str.split(':')]

            if len(ratios) != len(other_comps):
                raise ValueError(
                    f"比例数量不匹配！\n"
                    f"扫描组分: {scan_comp}\n"
                    f"其他组分: {other_comps} ({len(other_comps)}个)\n"
                    f"输入比例: {other_ratio_str} ({len(ratios)}个值)\n\n"
                    f"请输入 {len(other_comps)} 个比例值，用冒号分隔。\n"
                    f"例如: {':'.join(['1'] * len(other_comps))}"
                )

            # 对每个模型计算
            self.results_data = {'x_scan': x_scan_range, 'scan_comp': scan_comp}

            for model_name in selected_models:
                self.log(f"\n{'='*60}")
                self.log(f"开始计算: {model_name}")
                self.log(f"{'='*60}")

                self.progress_var.set(f"计算中: {model_name}")

                liquidus_temps = []
                ratio_sum = sum(ratios)

                for idx, x_scan in enumerate(x_scan_range):
                    # 计算成分
                    remaining = 1.0 - x_scan

                    composition = {scan_comp: x_scan}
                    for i, comp in enumerate(other_comps):
                        composition[comp] = remaining * ratios[i] / ratio_sum

                    # 构建模型字典（液相用选定模型，其他相用RKM）
                    model_class = self.available_models[model_name]
                    models = {}
                    for phase in self.available_phases:
                        if 'LIQUID' in phase.upper():
                            models[phase] = model_class
                        else:
                            models[phase] = Model

                    # 计算液相线
                    T_liq = self.find_liquidus(
                        self.dbe, study_comps, self.available_phases,
                        composition, temp_min, temp_max, temp_step, models
                    )

                    liquidus_temps.append(T_liq)

                    comp_str = ", ".join([f"X({k})={v:.3f}" for k, v in composition.items()])
                    if T_liq is not None:
                        self.log(f"  {comp_str} -> T_liq = {T_liq:.1f} K")
                    else:
                        self.log(f"  {comp_str} -> 未找到")

                self.results_data[model_name] = {
                    'liquidus': np.array(liquidus_temps),
                    'type': 'liquidus'
                }

                self.log(f"{model_name} 计算完成！\n")

            # 绘图
            self.plot_liquidus_comparison()
            self.display_liquidus_data()

            self.progress_var.set("计算完成！")
            self.log("\n所有模型计算完成！")

        except Exception as e:
            self.log(f"错误: {e}")
            import traceback
            self.log(traceback.format_exc())
            self.progress_var.set("计算失败")

    def find_liquidus(self, dbe, comps, phases, composition,
                     temp_min, temp_max, temp_step, models):
        """查找液相线温度（修复版）"""

        T_range = (temp_min, temp_max, temp_step)

        conditions = {
            v.T: T_range,
            v.P: 101325,
        }

        # 只设置N-1个组分
        comps_to_set = comps[:-1]
        for comp in comps_to_set:
            if comp in composition:
                conditions[v.X(comp)] = composition[comp]

        try:
            eq = equilibrium(dbe, comps, phases, model=models, conditions=conditions)

            T_vals = eq.T.values
            phase_array = eq.Phase.values
            np_array = eq.NP.values

            liquidus_T = None

            # 修复后的逻辑：从低温向高温扫描，找到最后一个100%液相的点
            for t_idx in range(len(T_vals)):
                liquid_np = 0.0

                # 查找LIQUID相的分数
                for v_idx in range(phase_array.shape[-1]):
                    phase_name = phase_array[0, 0, t_idx, 0, 0, v_idx]
                    if isinstance(phase_name, bytes):
                        phase_name = phase_name.decode('utf-8')

                    if 'LIQUID' in phase_name.upper():
                        liquid_np = np_array[0, 0, t_idx, 0, 0, v_idx]
                        break

                # 如果液相分数>=99.5%，更新液相线温度
                if liquid_np >= 0.995:
                    liquidus_T = T_vals[t_idx]
                # 如果液相分数<99.5%且之前找到过100%液相，说明已经过了液相线
                elif liquidus_T is not None:
                    break

            return liquidus_T

        except Exception as e:
            self.log(f"    计算错误: {e}")
            return None

    def plot_liquidus_comparison(self):
        """绘制液相线对比图"""
        self.ax_liquidus.clear()

        x_scan = self.results_data['x_scan']
        scan_comp = self.results_data['scan_comp']

        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
        markers = ['o', 's', '^', 'v']

        idx = 0
        for model_name in ['RKM', 'Muggianu', 'Toop', 'UEM1']:
            if model_name in self.results_data and self.results_data[model_name]['type'] == 'liquidus':
                liquidus = self.results_data[model_name]['liquidus']

                # 过滤None和NaN
                liquidus_float = np.array([x if x is not None else np.nan for x in liquidus], dtype=float)
                valid = ~np.isnan(liquidus_float)

                if np.any(valid):
                    self.ax_liquidus.plot(x_scan[valid], liquidus_float[valid],
                                         marker=markers[idx], linestyle='-', linewidth=2,
                                         markersize=6, label=model_name,
                                         color=colors[idx])
                idx += 1

        self.ax_liquidus.set_xlabel(f'X({scan_comp})', fontsize=12, fontweight='bold')
        self.ax_liquidus.set_ylabel('Liquidus Temperature (K)', fontsize=12, fontweight='bold')
        self.ax_liquidus.set_title('Liquidus Temperature Comparison', fontsize=14, fontweight='bold')
        self.ax_liquidus.legend(fontsize=11, loc='best')
        self.ax_liquidus.grid(True, alpha=0.3, linestyle='--')

        self.fig_liquidus.tight_layout()
        self.canvas_liquidus.draw()

    def display_liquidus_data(self):
        """显示液相线数据"""
        self.data_text.delete('1.0', tk.END)

        x_scan = self.results_data['x_scan']
        scan_comp = self.results_data['scan_comp']

        # 表头
        header = f"X({scan_comp})\t"
        for model_name in ['RKM', 'Muggianu', 'Toop', 'UEM1']:
            if model_name in self.results_data and self.results_data[model_name]['type'] == 'liquidus':
                header += f"T_liq_{model_name} (K)\t"
        header += "\n"
        self.data_text.insert(tk.END, header)
        self.data_text.insert(tk.END, "=" * 100 + "\n")

        # 数据行
        for i, x_val in enumerate(x_scan):
            line = f"{x_val:.4f}\t"

            for model_name in ['RKM', 'Muggianu', 'Toop', 'UEM1']:
                if model_name in self.results_data and self.results_data[model_name]['type'] == 'liquidus':
                    T_liq = self.results_data[model_name]['liquidus'][i]
                    if T_liq is not None and not np.isnan(T_liq):
                        line += f"{T_liq:.2f}\t"
                    else:
                        line += "N/A\t"

            line += "\n"
            self.data_text.insert(tk.END, line)

    def calculate_properties(self):
        """计算热力学性质"""
        if self.dbe is None:
            messagebox.showwarning("警告", "请先加载数据库！")
            return

        selected_models = self.get_selected_models()
        if not selected_models:
            messagebox.showwarning("警告", "请至少选择一个模型！")
            return

        thread = threading.Thread(target=self._calculate_properties_thread,
                                 args=(selected_models,))
        thread.daemon = True
        thread.start()

    def _calculate_properties_thread(self, selected_models):
        """热力学性质计算线程"""
        try:
            self.log(f"\n{'='*60}")
            self.log("开始计算热力学性质")
            self.log(f"{'='*60}")

            # 解析用户指定的研究组分
            comps_str = self.comps_entry.get().strip().upper()
            study_comps = [c.strip() for c in comps_str.split(',') if c.strip()]

            # 验证组分是否都在数据库中
            invalid_comps = [c for c in study_comps if c not in self.available_comps]
            if invalid_comps:
                raise ValueError(
                    f"以下组分不在数据库中: {invalid_comps}\n"
                    f"数据库可用组分: {self.available_comps}"
                )

            if len(study_comps) < 2:
                raise ValueError("至少需要2个组分！")

            self.log(f"研究组分: {study_comps}")

            # 解析其他输入
            scan_comp = self.scan_comp_entry.get().strip().upper()
            scan_range_str = self.scan_range_entry.get().strip()
            start, stop, num = [float(x) for x in scan_range_str.split(',')]

            temp_calc = float(self.temp_max_entry.get())  # 使用最高温度作为计算温度
            other_ratio_str = self.other_ratio_entry.get().strip()

            x_scan_range = np.linspace(start, stop, int(num))

            # 检查扫描组分是否有效
            if scan_comp not in study_comps:
                raise ValueError(f"扫描组分 '{scan_comp}' 不在研究组分列表中: {study_comps}")

            # 检查比例数量
            other_comps = [c for c in study_comps if c != scan_comp]
            ratios = [float(r) for r in other_ratio_str.split(':')]

            if len(ratios) != len(other_comps):
                raise ValueError(
                    f"比例数量不匹配！\n"
                    f"扫描组分: {scan_comp}\n"
                    f"其他组分: {other_comps} ({len(other_comps)}个)\n"
                    f"输入比例: {other_ratio_str} ({len(ratios)}个值)\n\n"
                    f"请输入 {len(other_comps)} 个比例值，用冒号分隔。\n"
                    f"例如: {':'.join(['1'] * len(other_comps))}"
                )

            # 存储性质数据
            for model_name in selected_models:
                self.progress_var.set(f"计算性质: {model_name}")
                self.log(f"\n计算 {model_name} 模型的热力学性质...")

                gibbs_list = []
                activity_data = {comp: [] for comp in study_comps}
                ratio_sum = sum(ratios)

                for x_scan in x_scan_range:
                    # 计算成分
                    remaining = 1.0 - x_scan

                    composition = {scan_comp: x_scan}
                    for i, comp in enumerate(other_comps):
                        composition[comp] = remaining * ratios[i] / ratio_sum

                    # 构建模型
                    model_class = self.available_models[model_name]
                    models = {}
                    for phase in self.available_phases:
                        if 'LIQUID' in phase.upper():
                            models[phase] = model_class
                        else:
                            models[phase] = Model

                    # 平衡计算
                    conditions = {v.T: temp_calc, v.P: 101325}
                    comps_to_set = study_comps[:-1]
                    for comp in comps_to_set:
                        if comp in composition:
                            conditions[v.X(comp)] = composition[comp]

                    try:
                        eq = equilibrium(self.dbe, study_comps,
                                       self.available_phases, model=models,
                                       conditions=conditions)

                        # 提取Gibbs自由能
                        G = eq.GM.values[0, 0, 0, 0, 0]
                        gibbs_list.append(G)

                        # 提取活度
                        for comp in activity_data.keys():
                            try:
                                act = eq.ACR(comp).values[0, 0, 0, 0, 0]
                                activity_data[comp].append(act)
                            except:
                                activity_data[comp].append(np.nan)

                    except Exception as e:
                        gibbs_list.append(np.nan)
                        for comp in activity_data.keys():
                            activity_data[comp].append(np.nan)
                        self.log(f"  计算失败 (X({scan_comp})={x_scan:.3f}): {e}")

                self.results_data[f'{model_name}_gibbs'] = {
                    'x_scan': x_scan_range,
                    'scan_comp': scan_comp,
                    'gibbs': np.array(gibbs_list),
                    'activity': activity_data,
                    'temperature': temp_calc,
                    'type': 'properties'
                }

                self.log(f"{model_name} 性质计算完成！")

            # 绘图
            self.plot_gibbs()
            self.plot_activity()

            self.progress_var.set("性质计算完成！")
            self.log("\n热力学性质计算完成！")

        except Exception as e:
            self.log(f"错误: {e}")
            import traceback
            self.log(traceback.format_exc())
            self.progress_var.set("计算失败")

    def plot_gibbs(self):
        """绘制Gibbs自由能图"""
        self.ax_gibbs.clear()

        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
        markers = ['o', 's', '^', 'v']

        idx = 0
        for model_name in ['RKM', 'Muggianu', 'Toop', 'UEM1']:
            key = f'{model_name}_gibbs'
            if key in self.results_data:
                data = self.results_data[key]
                x_scan = data['x_scan']
                gibbs = data['gibbs']
                scan_comp = data['scan_comp']
                temp = data['temperature']

                valid = ~np.isnan(gibbs)

                if np.any(valid):
                    self.ax_gibbs.plot(x_scan[valid], gibbs[valid],
                                      marker=markers[idx], linestyle='-',
                                      linewidth=2, markersize=6,
                                      label=model_name, color=colors[idx])
                idx += 1

        if 'scan_comp' in self.results_data:
            scan_comp = self.results_data['scan_comp']
        else:
            # 从gibbs数据中获取
            for key in self.results_data:
                if '_gibbs' in key:
                    scan_comp = self.results_data[key]['scan_comp']
                    temp = self.results_data[key]['temperature']
                    break

        self.ax_gibbs.set_xlabel(f'X({scan_comp})', fontsize=12, fontweight='bold')
        self.ax_gibbs.set_ylabel('Gibbs Energy (J/mol)', fontsize=12, fontweight='bold')
        self.ax_gibbs.set_title(f'Gibbs Free Energy at {temp:.0f} K',
                               fontsize=14, fontweight='bold')
        self.ax_gibbs.legend(fontsize=10, loc='best')
        self.ax_gibbs.grid(True, alpha=0.3, linestyle='--')

        self.fig_gibbs.tight_layout()
        self.canvas_gibbs.draw()

    def plot_activity(self):
        """绘制活度图"""
        self.ax_activity.clear()

        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
        linestyles = ['-', '--', '-.', ':']

        # 只绘制第一个模型的活度（多个组分）
        for model_name in ['RKM', 'Muggianu', 'Toop', 'UEM1']:
            key = f'{model_name}_gibbs'
            if key in self.results_data:
                data = self.results_data[key]
                x_scan = data['x_scan']
                activity_data = data['activity']
                scan_comp = data['scan_comp']
                temp = data['temperature']

                for idx, (comp, acts) in enumerate(activity_data.items()):
                    acts_array = np.array(acts)
                    valid = ~np.isnan(acts_array)

                    if np.any(valid):
                        self.ax_activity.plot(x_scan[valid], acts_array[valid],
                                            linestyle=linestyles[idx % 4],
                                            linewidth=2, marker='o', markersize=4,
                                            label=f'{comp} ({model_name})',
                                            color=colors[idx % 4])

                # 只绘制第一个找到的模型
                break

        self.ax_activity.set_xlabel(f'X({scan_comp})', fontsize=12, fontweight='bold')
        self.ax_activity.set_ylabel('Activity', fontsize=12, fontweight='bold')
        self.ax_activity.set_title(f'Activity at {temp:.0f} K',
                                  fontsize=14, fontweight='bold')
        self.ax_activity.legend(fontsize=10, loc='best')
        self.ax_activity.grid(True, alpha=0.3, linestyle='--')

        self.fig_activity.tight_layout()
        self.canvas_activity.draw()

    def calculate_phase_diagram(self):
        """生成相图并保存"""
        if self.dbe is None:
            messagebox.showwarning("警告", "请先加载数据库！")
            return

        selected_models = self.get_selected_models()
        if not selected_models:
            messagebox.showwarning("警告", "请至少选择一个模型！")
            return

        # 只使用第一个选中的模型
        model_name = selected_models[0]

        thread = threading.Thread(target=self._calculate_phase_diagram_thread,
                                 args=(model_name,))
        thread.daemon = True
        thread.start()

    def _calculate_phase_diagram_thread(self, model_name):
        """相图计算线程"""
        try:
            self.log(f"\n{'='*60}")
            self.log(f"开始计算相图 (模型: {model_name})")
            self.log(f"{'='*60}")

            # 解析用户指定的研究组分
            comps_str = self.comps_entry.get().strip().upper()
            study_comps = [c.strip() for c in comps_str.split(',') if c.strip()]

            # 验证组分是否都在数据库中
            invalid_comps = [c for c in study_comps if c not in self.available_comps]
            if invalid_comps:
                raise ValueError(
                    f"以下组分不在数据库中: {invalid_comps}\n"
                    f"数据库可用组分: {self.available_comps}"
                )

            if len(study_comps) < 2:
                raise ValueError("至少需要2个组分！")

            self.log(f"研究组分: {study_comps}")

            # 解析其他输入
            scan_comp = self.scan_comp_entry.get().strip().upper()
            scan_range_str = self.scan_range_entry.get().strip()
            start, stop, num = [float(x) for x in scan_range_str.split(',')]

            temp_min = float(self.temp_min_entry.get())
            temp_max = float(self.temp_max_entry.get())
            temp_num = 50  # 温度点数

            other_ratio_str = self.other_ratio_entry.get().strip()

            x_scan_range = np.linspace(start, stop, int(num))
            T_range = np.linspace(temp_min, temp_max, temp_num)

            # 检查扫描组分是否有效
            if scan_comp not in study_comps:
                raise ValueError(f"扫描组分 '{scan_comp}' 不在研究组分列表中: {study_comps}")

            # 检查比例数量
            other_comps = [c for c in study_comps if c != scan_comp]
            ratios = [float(r) for r in other_ratio_str.split(':')]

            if len(ratios) != len(other_comps):
                raise ValueError(
                    f"比例数量不匹配！\n"
                    f"扫描组分: {scan_comp}\n"
                    f"其他组分: {other_comps} ({len(other_comps)}个)\n"
                    f"输入比例: {other_ratio_str} ({len(ratios)}个值)\n\n"
                    f"请输入 {len(other_comps)} 个比例值，用冒号分隔。\n"
                    f"例如: {':'.join(['1'] * len(other_comps))}"
                )

            self.progress_var.set(f"生成相图: {model_name}")

            # 构建模型
            model_class = self.available_models[model_name]
            models = {}
            for phase in self.available_phases:
                if 'LIQUID' in phase.upper():
                    models[phase] = model_class
                else:
                    models[phase] = Model

            # 创建网格数据
            phase_data = np.zeros((len(T_range), len(x_scan_range)), dtype=object)

            total_calcs = len(x_scan_range) * len(T_range)
            calc_count = 0
            ratio_sum = sum(ratios)

            for i, x_scan in enumerate(x_scan_range):
                # 计算成分
                remaining = 1.0 - x_scan

                composition = {scan_comp: x_scan}
                for j, comp in enumerate(other_comps):
                    composition[comp] = remaining * ratios[j] / ratio_sum

                # 对每个温度计算
                for j, T in enumerate(T_range):
                    conditions = {v.T: T, v.P: 101325}
                    comps_to_set = study_comps[:-1]
                    for comp in comps_to_set:
                        if comp in composition:
                            conditions[v.X(comp)] = composition[comp]

                    try:
                        eq = equilibrium(self.dbe, study_comps,
                                       self.available_phases, model=models,
                                       conditions=conditions)

                        # 找到主相（相分数最大的相）
                        phase_array = eq.Phase.values
                        np_array = eq.NP.values

                        max_np = 0.0
                        dominant_phase = ''

                        for v_idx in range(phase_array.shape[-1]):
                            phase_name = phase_array[0, 0, 0, 0, 0, v_idx]
                            if isinstance(phase_name, bytes):
                                phase_name = phase_name.decode('utf-8')
                            phase_frac = np_array[0, 0, 0, 0, 0, v_idx]

                            if phase_frac > max_np:
                                max_np = phase_frac
                                dominant_phase = phase_name

                        phase_data[j, i] = dominant_phase

                    except Exception as e:
                        phase_data[j, i] = 'ERROR'

                    calc_count += 1
                    if calc_count % 10 == 0:
                        progress = calc_count / total_calcs * 100
                        self.progress_var.set(f"相图计算: {progress:.1f}%")

            # 绘制相图
            self.plot_phase_diagram(x_scan_range, T_range, phase_data, scan_comp, model_name)

            # 保存相图
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"phase_diagram_{model_name}_{timestamp}.png"
            self.fig_phase.savefig(filename, dpi=300, bbox_inches='tight')

            self.log(f"\n相图已保存: {filename}")
            self.progress_var.set("相图生成完成！")

        except Exception as e:
            self.log(f"错误: {e}")
            import traceback
            self.log(traceback.format_exc())
            self.progress_var.set("计算失败")

    def plot_phase_diagram(self, x_range, T_range, phase_data, scan_comp, model_name):
        """绘制相图"""
        self.ax_phase.clear()

        # 为每个相分配颜色
        unique_phases = list(set(phase_data.flatten()))
        color_map = {}
        cmap = plt.cm.get_cmap('tab10')
        for i, phase in enumerate(unique_phases):
            color_map[phase] = cmap(i / len(unique_phases))

        # 转换为数值数据
        phase_numeric = np.zeros_like(phase_data, dtype=float)
        for i, phase in enumerate(unique_phases):
            phase_numeric[phase_data == phase] = i

        # 绘制
        X, Y = np.meshgrid(x_range, T_range)
        im = self.ax_phase.pcolormesh(X, Y, phase_numeric, shading='auto',
                                     cmap='tab10', vmin=0, vmax=len(unique_phases))

        self.ax_phase.set_xlabel(f'X({scan_comp})', fontsize=12, fontweight='bold')
        self.ax_phase.set_ylabel('Temperature (K)', fontsize=12, fontweight='bold')
        self.ax_phase.set_title(f'Phase Diagram ({model_name})',
                               fontsize=14, fontweight='bold')

        # 添加图例
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=color_map[phase], label=phase)
                          for phase in unique_phases if phase != 'ERROR']
        self.ax_phase.legend(handles=legend_elements, loc='upper left',
                            fontsize=9, framealpha=0.9)

        self.ax_phase.grid(True, alpha=0.2, linestyle='--')

        self.fig_phase.tight_layout()
        self.canvas_phase.draw()

    def export_data(self):
        """导出数据"""
        if not self.results_data:
            messagebox.showwarning("警告", "没有可导出的数据！")
            return

        file_path = filedialog.asksaveasfilename(
            title="保存数据",
            defaultextension=".txt",
            filetypes=[("文本文件", "*.txt"), ("CSV文件", "*.csv")]
        )

        if file_path:
            try:
                with open(file_path, 'w', encoding='utf-8') as f:
                    f.write(self.data_text.get('1.0', tk.END))

                messagebox.showinfo("成功", f"数据已导出到:\n{file_path}")
                self.log(f"数据已导出: {file_path}")

            except Exception as e:
                messagebox.showerror("错误", f"导出失败:\n{e}")

    def clear_results(self):
        """清除结果"""
        self.ax_liquidus.clear()
        self.ax_gibbs.clear()
        self.ax_activity.clear()
        self.ax_phase.clear()

        self.canvas_liquidus.draw()
        self.canvas_gibbs.draw()
        self.canvas_activity.draw()
        self.canvas_phase.draw()

        self.data_text.delete('1.0', tk.END)
        self.results_data = {}

        self.progress_var.set("就绪")
        self.log("\n结果已清除")


def main():
    """主函数"""
    root = tk.Tk()
    app = AlloyCalculatorGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
