#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
合金液相线和固相线计算图形界面

功能：
1. 加载热力学数据库
2. 选择计算模型（RKM, Muggianu, Toop, UEM1）
3. 设定合金成分
4. 设定温度范围
5. 计算液相线/固相线
6. 图形化显示结果
7. 导出数据
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np
from pycalphad import Database, equilibrium, variables as v
from pycalphad.model import Model
import os
import threading


class PhaseDiagramGUI:
    """合金相图计算GUI"""

    def __init__(self, root):
        self.root = root
        self.root.title("合金液相线/固相线计算工具")
        self.root.geometry("1400x900")

        # 数据存储
        self.dbe = None
        self.dbe_path = None
        self.available_models = {
            'RKM (Redlich-Kister-Muggianu)': Model,
            'Muggianu': None,  # 将在加载时动态导入
            'Toop': None,
            'UEM1': None
        }
        self.calculation_running = False
        self.results_data = None

        # 创建界面
        self.create_widgets()

    def create_widgets(self):
        """创建所有界面组件"""

        # ============ 左侧控制面板 ============
        left_frame = ttk.Frame(self.root, padding="10")
        left_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)

        # 1. 数据库加载区域
        db_frame = ttk.LabelFrame(left_frame, text="数据库", padding="10")
        db_frame.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=5)

        ttk.Button(db_frame, text="加载数据库", command=self.load_database).grid(
            row=0, column=0, padx=5, pady=5, sticky=tk.W)

        self.db_label = ttk.Label(db_frame, text="未加载", foreground="red")
        self.db_label.grid(row=0, column=1, padx=5, pady=5, sticky=tk.W)

        # 2. 模型选择区域
        model_frame = ttk.LabelFrame(left_frame, text="模型选择", padding="10")
        model_frame.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=5)

        ttk.Label(model_frame, text="液相模型:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.model_var = tk.StringVar(value='RKM (Redlich-Kister-Muggianu)')
        self.model_combo = ttk.Combobox(model_frame, textvariable=self.model_var,
                                        values=list(self.available_models.keys()),
                                        state='readonly', width=30)
        self.model_combo.grid(row=0, column=1, pady=2, padx=5)

        # 3. 成分设定区域
        comp_frame = ttk.LabelFrame(left_frame, text="成分设定", padding="10")
        comp_frame.grid(row=2, column=0, sticky=(tk.W, tk.E), pady=5)

        ttk.Label(comp_frame, text="组分列表 (逗号分隔):").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.comps_entry = ttk.Entry(comp_frame, width=30)
        self.comps_entry.insert(0, "AL,CR,NI")
        self.comps_entry.grid(row=0, column=1, pady=2, padx=5)

        ttk.Label(comp_frame, text="扫描组分:").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.scan_comp_entry = ttk.Entry(comp_frame, width=30)
        self.scan_comp_entry.insert(0, "NI")
        self.scan_comp_entry.grid(row=1, column=1, pady=2, padx=5)

        ttk.Label(comp_frame, text="扫描范围 (起,止,点数):").grid(row=2, column=0, sticky=tk.W, pady=2)
        self.scan_range_entry = ttk.Entry(comp_frame, width=30)
        self.scan_range_entry.insert(0, "0.1,0.9,10")
        self.scan_range_entry.grid(row=2, column=1, pady=2, padx=5)

        ttk.Label(comp_frame, text="其他组分比例:").grid(row=3, column=0, sticky=tk.W, pady=2)
        self.other_ratio_entry = ttk.Entry(comp_frame, width=30)
        self.other_ratio_entry.insert(0, "1:1")
        self.other_ratio_entry.grid(row=3, column=1, pady=2, padx=5)

        # 4. 温度范围设定
        temp_frame = ttk.LabelFrame(left_frame, text="温度范围 (K)", padding="10")
        temp_frame.grid(row=3, column=0, sticky=(tk.W, tk.E), pady=5)

        ttk.Label(temp_frame, text="最低温度:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.temp_min_entry = ttk.Entry(temp_frame, width=15)
        self.temp_min_entry.insert(0, "1400")
        self.temp_min_entry.grid(row=0, column=1, pady=2, padx=5)

        ttk.Label(temp_frame, text="最高温度:").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.temp_max_entry = ttk.Entry(temp_frame, width=15)
        self.temp_max_entry.insert(0, "2400")
        self.temp_max_entry.grid(row=1, column=1, pady=2, padx=5)

        ttk.Label(temp_frame, text="温度步长:").grid(row=2, column=0, sticky=tk.W, pady=2)
        self.temp_step_entry = ttk.Entry(temp_frame, width=15)
        self.temp_step_entry.insert(0, "10")
        self.temp_step_entry.grid(row=2, column=1, pady=2, padx=5)

        # 5. 相选择区域
        phase_frame = ttk.LabelFrame(left_frame, text="相选择", padding="10")
        phase_frame.grid(row=4, column=0, sticky=(tk.W, tk.E), pady=5)

        ttk.Label(phase_frame, text="相列表 (逗号分隔):").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.phases_entry = ttk.Entry(phase_frame, width=30)
        self.phases_entry.insert(0, "LIQUID,FCC_A1,BCC_A2,B2,L12_FCC")
        self.phases_entry.grid(row=0, column=1, pady=2, padx=5)

        # 6. 计算按钮区域
        button_frame = ttk.Frame(left_frame, padding="10")
        button_frame.grid(row=5, column=0, sticky=(tk.W, tk.E), pady=10)

        self.calc_button = ttk.Button(button_frame, text="开始计算",
                                      command=self.start_calculation)
        self.calc_button.grid(row=0, column=0, padx=5)

        ttk.Button(button_frame, text="导出数据",
                  command=self.export_data).grid(row=0, column=1, padx=5)

        ttk.Button(button_frame, text="清除结果",
                  command=self.clear_results).grid(row=0, column=2, padx=5)

        # 进度显示
        self.progress_var = tk.StringVar(value="就绪")
        self.progress_label = ttk.Label(button_frame, textvariable=self.progress_var)
        self.progress_label.grid(row=1, column=0, columnspan=3, pady=5)

        # ============ 右侧结果显示区域 ============
        right_frame = ttk.Frame(self.root, padding="10")
        right_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)

        # 创建Notebook用于多标签页显示
        self.notebook = ttk.Notebook(right_frame)
        self.notebook.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # 图形标签页
        self.plot_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.plot_frame, text="图形结果")

        # 创建matplotlib图形
        self.fig = Figure(figsize=(10, 6))
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # 添加matplotlib工具栏
        toolbar = NavigationToolbar2Tk(self.canvas, self.plot_frame)
        toolbar.update()

        # 数据标签页
        self.data_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.data_frame, text="数据表格")

        self.data_text = scrolledtext.ScrolledText(self.data_frame, width=80, height=30)
        self.data_text.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # 日志标签页
        self.log_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.log_frame, text="计算日志")

        self.log_text = scrolledtext.ScrolledText(self.log_frame, width=80, height=30)
        self.log_text.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # 配置网格权重
        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)
        right_frame.columnconfigure(0, weight=1)
        right_frame.rowconfigure(0, weight=1)

    def log(self, message):
        """添加日志消息"""
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()

    def load_database(self):
        """加载热力学数据库"""
        file_path = filedialog.askopenfilename(
            title="选择数据库文件",
            filetypes=[("TDB文件", "*.tdb"), ("所有文件", "*.*")],
            initialdir="examples"
        )

        if file_path:
            try:
                self.dbe = Database(file_path)
                self.dbe_path = file_path
                filename = os.path.basename(file_path)
                self.db_label.config(text=f"已加载: {filename}", foreground="green")
                self.log(f"成功加载数据库: {file_path}")

                # 尝试导入高级模型
                try:
                    from pycalphad.advanced_uem_model import ModelMuggianu, ModelToop
                    from pycalphad.uem1_Model import uem1_model
                    self.available_models['Muggianu'] = ModelMuggianu
                    self.available_models['Toop'] = ModelToop
                    self.available_models['UEM1'] = uem1_model
                    self.log("成功加载所有模型: RKM, Muggianu, Toop, UEM1")
                except ImportError as e:
                    self.log(f"警告: 部分高级模型加载失败: {e}")
                    self.log("仅支持 RKM 模型")

            except Exception as e:
                messagebox.showerror("错误", f"加载数据库失败:\n{e}")
                self.log(f"错误: {e}")

    def start_calculation(self):
        """开始计算（在新线程中运行）"""
        if self.dbe is None:
            messagebox.showwarning("警告", "请先加载数据库！")
            return

        if self.calculation_running:
            messagebox.showinfo("提示", "计算正在进行中，请稍候...")
            return

        # 在新线程中运行计算
        calc_thread = threading.Thread(target=self.calculate)
        calc_thread.daemon = True
        calc_thread.start()

    def calculate(self):
        """执行液相线/固相线计算"""
        try:
            self.calculation_running = True
            self.calc_button.config(state='disabled')
            self.progress_var.set("计算中...")

            # 获取参数
            comps_str = self.comps_entry.get().strip().upper()
            comps = [c.strip() for c in comps_str.split(',')]

            phases_str = self.phases_entry.get().strip().upper()
            phases = [p.strip() for p in phases_str.split(',')]

            scan_comp = self.scan_comp_entry.get().strip().upper()
            scan_range_str = self.scan_range_entry.get().strip()
            start, stop, num = [float(x) for x in scan_range_str.split(',')]

            temp_min = float(self.temp_min_entry.get())
            temp_max = float(self.temp_max_entry.get())
            temp_step = float(self.temp_step_entry.get())

            other_ratio_str = self.other_ratio_entry.get().strip()

            self.log(f"\n{'='*60}")
            self.log("开始计算液相线")
            self.log(f"{'='*60}")
            self.log(f"组分: {comps}")
            self.log(f"相: {phases}")
            self.log(f"扫描组分: {scan_comp}")
            self.log(f"扫描范围: {start} ~ {stop}, {int(num)}个点")
            self.log(f"温度范围: {temp_min} ~ {temp_max} K, 步长 {temp_step} K")

            # 生成扫描成分数组
            x_scan_range = np.linspace(start, stop, int(num))

            # 选择模型
            model_name = self.model_var.get()
            model_class = self.available_models[model_name]

            if model_class is None:
                raise ValueError(f"模型 {model_name} 未加载")

            # 构建模型字典
            if model_name == 'RKM (Redlich-Kister-Muggianu)':
                models = {ph: Model for ph in phases}
            else:
                models = {ph: model_class if ph == 'LIQUID' else Model for ph in phases}

            self.log(f"使用模型: {model_name}")

            # 计算液相线
            liquidus_temps = []

            for idx, x_scan in enumerate(x_scan_range):
                self.progress_var.set(f"计算中... {idx+1}/{len(x_scan_range)}")

                # 根据比例计算其他组分
                other_comps = [c for c in comps if c != scan_comp]
                remaining = 1.0 - x_scan

                # 解析比例
                ratios = [float(r) for r in other_ratio_str.split(':')]
                ratio_sum = sum(ratios)

                composition = {scan_comp: x_scan}
                for i, comp in enumerate(other_comps):
                    composition[comp] = remaining * ratios[i] / ratio_sum

                # 计算液相线
                T_liquidus = self.find_liquidus(
                    self.dbe, comps, phases, composition,
                    temp_min, temp_max, temp_step, models
                )

                liquidus_temps.append(T_liquidus)

                comp_str = ", ".join([f"X({k})={v:.3f}" for k, v in composition.items()])
                if T_liquidus is not None:
                    self.log(f"  {comp_str} -> T_liquidus = {T_liquidus:.1f} K")
                else:
                    self.log(f"  {comp_str} -> 未找到液相线")

            # 保存结果
            self.results_data = {
                'x_scan': x_scan_range,
                'scan_comp': scan_comp,
                'liquidus': np.array(liquidus_temps),
                'comps': comps,
                'model': model_name
            }

            # 绘图
            self.plot_results()

            # 显示数据
            self.display_data()

            self.log(f"\n{'='*60}")
            self.log("计算完成！")
            self.log(f"{'='*60}")
            self.progress_var.set("计算完成")

        except Exception as e:
            self.log(f"\n错误: {e}")
            messagebox.showerror("计算错误", str(e))
            self.progress_var.set("计算失败")

        finally:
            self.calculation_running = False
            self.calc_button.config(state='normal')

    def find_liquidus(self, dbe, comps, phases, composition,
                     temp_min, temp_max, temp_step, models):
        """查找液相线温度"""

        T_range = (temp_min, temp_max, temp_step)

        conditions = {
            v.T: T_range,
            v.P: 101325,
        }

        # 对于N组分系统，只能指定N-1个独立的摩尔分数
        # 最后一个组分的摩尔分数由总和=1自动确定
        comps_to_set = comps[:-1]  # 取前N-1个组分
        for comp in comps_to_set:
            if comp in composition:
                conditions[v.X(comp)] = composition[comp]

        try:
            eq = equilibrium(dbe, comps, phases, model=models, conditions=conditions)

            T_vals = eq.T.values
            phase_array = eq.Phase.values
            np_array = eq.NP.values

            liquidus_T = None

            # 从高温到低温扫描，找到100%液相的最低温度
            for t_idx in range(len(T_vals) - 1, -1, -1):
                liquid_np = None
                for v_idx in range(phase_array.shape[-1]):
                    phase_name = phase_array[0, 0, t_idx, 0, 0, v_idx]
                    if phase_name == 'LIQUID':
                        liquid_np = np_array[0, 0, t_idx, 0, 0, v_idx]
                        break

                if liquid_np is not None and liquid_np >= 0.995:
                    liquidus_T = T_vals[t_idx]
                elif liquid_np is not None and liquid_np < 0.995:
                    break

            return liquidus_T

        except Exception as e:
            self.log(f"    计算错误: {e}")
            return None

    def plot_results(self):
        """绘制结果图"""
        if self.results_data is None:
            return

        x_scan = self.results_data['x_scan']
        scan_comp = self.results_data['scan_comp']
        liquidus = self.results_data['liquidus']
        model = self.results_data['model']

        self.ax.clear()

        # 过滤None和NaN值
        # 首先将None转换为NaN，然后过滤
        liquidus_float = np.array([x if x is not None else np.nan for x in liquidus], dtype=float)
        valid = ~np.isnan(liquidus_float)

        self.ax.plot(x_scan[valid], liquidus_float[valid], 'o-',
                    linewidth=2, markersize=6, label=model)

        self.ax.set_xlabel(f'X({scan_comp})', fontsize=12)
        self.ax.set_ylabel('Liquidus Temperature (K)', fontsize=12)
        self.ax.set_title(f'Liquidus Temperature vs {scan_comp} Composition',
                         fontsize=14, fontweight='bold')
        self.ax.legend(fontsize=11)
        self.ax.grid(True, alpha=0.3)

        self.fig.tight_layout()
        self.canvas.draw()

    def display_data(self):
        """显示数据表格"""
        if self.results_data is None:
            return

        x_scan = self.results_data['x_scan']
        liquidus = self.results_data['liquidus']
        scan_comp = self.results_data['scan_comp']
        comps = self.results_data['comps']

        self.data_text.delete('1.0', tk.END)

        # 表头
        header = f"X({scan_comp})\t"
        for comp in comps:
            if comp != scan_comp:
                header += f"X({comp})\t"
        header += "T_liquidus (K)\n"
        self.data_text.insert(tk.END, header)
        self.data_text.insert(tk.END, "=" * 80 + "\n")

        # 数据行
        other_comps = [c for c in comps if c != scan_comp]
        other_ratio_str = self.other_ratio_entry.get().strip()
        ratios = [float(r) for r in other_ratio_str.split(':')]
        ratio_sum = sum(ratios)

        for x_scan_val, T_liq in zip(x_scan, liquidus):
            line = f"{x_scan_val:.4f}\t"

            remaining = 1.0 - x_scan_val
            for i, comp in enumerate(other_comps):
                x_comp = remaining * ratios[i] / ratio_sum
                line += f"{x_comp:.4f}\t"

            if T_liq is not None and not np.isnan(T_liq):
                line += f"{T_liq:.2f}\n"
            else:
                line += "N/A\n"

            self.data_text.insert(tk.END, line)

    def export_data(self):
        """导出数据到文件"""
        if self.results_data is None:
            messagebox.showwarning("警告", "没有可导出的数据！")
            return

        file_path = filedialog.asksaveasfilename(
            title="保存数据",
            defaultextension=".txt",
            filetypes=[("文本文件", "*.txt"), ("CSV文件", "*.csv"), ("所有文件", "*.*")]
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
        self.ax.clear()
        self.canvas.draw()
        self.data_text.delete('1.0', tk.END)
        self.results_data = None
        self.progress_var.set("就绪")
        self.log("\n结果已清除")


def main():
    """主函数"""
    root = tk.Tk()
    app = PhaseDiagramGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
