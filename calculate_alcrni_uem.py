#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
对比测试脚本：
对比使用 uem1_model（所有相）和传统 Model 计算 Al-Cr-Ni 系统液相线温度的差异。
条件：xAl/xCr = 1/1
"""

# --- 导入 ---
import numpy as np
import matplotlib.pyplot as plt
from pycalphad import Database, equilibrium, variables as v, Model
from pycalphad.core.errors import DofError
import logging
from pycalphad.uem1_Model import uem1_model
from pycalphad.advanced_uem_model import ModelUEM1
# 禁用 pycalphad 的冗余日志
logging.getLogger('pycalphad').setLevel(logging.WARNING)


# ---------------------------------------------------------------------------
# --- 液相线计算辅助函数 ---
# ---------------------------------------------------------------------------

def calculate_liquidus_line(dbe, comps, phases, x_ni_range, model_dict, label):
	"""
	计算液相线温度（xAl/xCr = 1/1）

	参数：
	- dbe: Database 对象
	- comps: 组元列表
	- phases: 相列表
	- x_ni_range: Ni 含量数组
	- model_dict: 模型字典（可以是 uem1_model 或 None）
	- label: 标签（用于打印）

	返回：
	- liquidus_temps: 液相线温度列表
	"""
	liquidus_temps = []

	print(f"\n开始计算 {label} 模型的液相线...")
	print("=" * 70)

	for x_ni in x_ni_range:
		x_al = (1.0 - x_ni) / 2.0
		x_cr = (1.0 - x_ni) / 2.0

		# 跳过极端成分
		if x_al < 0.01 or x_cr < 0.01:
			liquidus_temps.append(np.nan)
			continue

		# 设置温度范围，寻找液相线
		T_range = (1400, 2400, 10)  # 10K 步长，覆盖可能的液相线范围

		conditions = {
			v.T: T_range,
			v.P: 101325,
			v.X('AL'): x_al,
			v.X('CR'): x_cr,
		}

		try:
			# 计算平衡
			eq = equilibrium(dbe, comps, phases, model=model_dict, conditions=conditions)

			# 寻找液相线温度（LIQUID 相分数从 1 降到 <1 的温度）
			# 即：从高温到低温，找到100%液相的最低温度
			T_vals = eq.T.values
			phase_array = eq.Phase.values
			np_array = eq.NP.values

			liquidus_T = None

			# 从高温到低温扫描，找到液相线（100%液相的最低温度）
			for t_idx in range(len(T_vals) - 1, -1, -1):
				# 检查 LIQUID 相分数
				liquid_np = None
				for v_idx in range(phase_array.shape[-1]):
					phase_name = phase_array[0, 0, t_idx, 0, 0, v_idx]
					if phase_name == 'LIQUID':
						liquid_np = np_array[0, 0, t_idx, 0, 0, v_idx]
						break

				# 如果 LIQUID 相分数 >= 0.995（基本是纯液相）
				if liquid_np is not None and liquid_np >= 0.995:
					liquidus_T = T_vals[t_idx]
					# 继续向下找，直到找到更低的纯液相温度
				elif liquid_np is not None and liquid_np < 0.995:
					# 已经进入两相区，前面找到的温度就是液相线
					break

			if liquidus_T is not None:
				liquidus_temps.append(liquidus_T)
				print(f"X(NI)={x_ni:.3f}, X(AL)={x_al:.3f}, X(CR)={x_cr:.3f} -> T_liquidus={liquidus_T:.1f} K ({label})")
			else:
				liquidus_temps.append(np.nan)
				print(f"X(NI)={x_ni:.3f}, X(AL)={x_al:.3f}, X(CR)={x_cr:.3f} -> 未找到液相线 ({label})")

		except Exception as e:
			liquidus_temps.append(np.nan)
			print(f"X(NI)={x_ni:.3f} -> 计算错误: {e}")

	return liquidus_temps


# ---------------------------------------------------------------------------
# --- 主测试函数 ---
# ---------------------------------------------------------------------------

def main():
	"""
	对比使用 uem1_model（所有相）和传统 Model 计算液相线的差异
	"""

	print("=" * 70)
	print("Al-Cr-Ni 体系液相线温度计算（uem1_model vs Model 对比）")
	print("条件: xAl/xCr = 1/1")
	print("说明: uem1_model 应用于所有相，Model 使用默认的修正 Muggianu 方法")
	print("=" * 70)

	# 1. 加载数据库
	try:
		dbe = Database('examples/alcrni.tdb')
		print("\n数据库加载成功: examples/alcrni.tdb")
	except Exception as e:
		print(f"错误: 无法加载数据库: {e}")
		return

	comps = ['AL', 'CR', 'NI']
	phases = ['LIQUID', 'FCC_A1', 'BCC_A2', 'B2', 'L12_FCC']

	# 2. 设置成分范围
	num_points = 10
	x_ni_range = np.linspace(0.1, 0.9, num_points)

	# 3. 构建 uem1_model 字典（所有相使用 UEM1）
	print("\n正在为所有相构建 uem1_model 实例...")
	
			
	uem_models = {ph: uem1_model if ph == 'LIQUID' else Model for ph in phases}
	

	print(f"已成功为 {len(uem_models)} 个相构建 uem1_model: {list(uem_models.keys())}")

	# 4. 计算液相线（传统 Model）
	liquidus_model = calculate_liquidus_line(
		dbe, comps, phases, x_ni_range,
		model_dict=None,  # None 表示使用默认 Model
		label="传统 Model"
	)

	# 5. 计算液相线（uem1_model - 所有相）
	liquidus_uem = calculate_liquidus_line(
		dbe, comps, phases, x_ni_range,
		model_dict=ModelUEM1,  # 所有相使用 uem1_model
		label="uem1_model（所有相）"
	)

	# 6. 绘制对比图
	print("\n" + "=" * 70)
	print("绘制液相线对比图")
	print("=" * 70)

	# 过滤掉 NaN 值
	valid_indices = [i for i in range(len(liquidus_model))
	                 if not np.isnan(liquidus_model[i]) and not np.isnan(liquidus_uem[i])]

	if len(valid_indices) > 0:
		x_ni_valid = x_ni_range[valid_indices]
		liquidus_model_valid = np.array(liquidus_model)[valid_indices]
		liquidus_uem_valid = np.array(liquidus_uem)[valid_indices]

		# --- 解决中文乱码问题 ---
		try:
			plt.rcParams['font.sans-serif'] = ['SimHei']  # 指定默认字体为黑体
			plt.rcParams['axes.unicode_minus'] = False  # 解决保存图像时负号'-'显示为方块的问题
			print("已尝试设置中文字体 'SimHei'。")
		except Exception as e:
			print(f"警告: 设置 'SimHei' 字体失败，中文可能仍会显示为乱码: {e}")
		# --- 结束 ---

		plt.figure(figsize=(10, 6))
		plt.plot(x_ni_valid, liquidus_model_valid, 'b-o', label='传统 Model（修正 Muggianu）', linewidth=2)
		plt.plot(x_ni_valid, liquidus_uem_valid, 'r-s', label='uem1_model（所有相）', linewidth=2)
		plt.xlabel('Mole Fraction of Ni ($x_{Ni}$)', fontsize=12)
		plt.ylabel('Liquidus Temperature (K)', fontsize=12)
		plt.title('Al-Cr-Ni Liquidus Comparison (xAl/xCr = 1/1)\nuem1_model vs Traditional Model', fontsize=14)
		plt.legend(fontsize=11)
		plt.grid(True, alpha=0.3)
		plt.xlim(0, 1)

		plot_filename = 'alcrni_liquidus_UEM_vs_Model.png'
		plt.savefig(plot_filename, dpi=150, bbox_inches='tight')
		print(f"图片已保存到: {plot_filename}")

		# 7. 保存数据文件
		data_filename = 'liquidus_uem_vs_model.txt'
		with open(data_filename, 'w', encoding='utf-8') as f:
			f.write("# Al-Cr-Ni Liquidus Temperature Data (xAl/xCr = 1/1)\n")
			f.write("# Comparison: uem1_model (all phases) vs Traditional Model\n")
			f.write("# X(Ni)\tX(Al)\tX(Cr)\tT_liquidus_Model(K)\tT_liquidus_UEM(K)\tDifference(K)\n")
			for i in valid_indices:
				x_ni = x_ni_range[i]
				x_al = (1.0 - x_ni) / 2.0
				x_cr = (1.0 - x_ni) / 2.0
				t_model = liquidus_model[i]
				t_uem = liquidus_uem[i]
				diff = t_uem - t_model
				f.write(f"{x_ni:.4f}\t{x_al:.4f}\t{x_cr:.4f}\t{t_model:.2f}\t{t_uem:.2f}\t{diff:.2f}\n")

		print(f"数据已保存到: {data_filename}")

		# 8. 统计信息
		print("\n" + "=" * 70)
		print("统计信息")
		print("=" * 70)

		differences = liquidus_uem_valid - liquidus_model_valid
		print(f"液相线温度差异 (uem1_model - Model):")
		print(f"  平均值: {np.mean(differences):.2f} K")
		print(f"  最大值: {np.max(differences):.2f} K")
		print(f"  最小值: {np.min(differences):.2f} K")
		print(f"  标准差: {np.std(differences):.2f} K")
		print(f"\n注意: uem1_model 应用于所有相（包括固相），")
		print(f"      这会改变固相稳定性，导致液相线位置明显不同。")
		print(f"      如果只想对比液相模型的影响，请使用 calculate_liquidus_alcrni.py")

	else:
		print("没有有效的数据点可供对比。")

	print("\n" + "=" * 70)
	print("测试完成！")
	print("=" * 70)


if __name__ == "__main__":
	main()