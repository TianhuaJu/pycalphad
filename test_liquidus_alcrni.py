"""
液相线温度计算测试 - Al-Cr-Ni 体系
计算条件: xAl/xCr = 1/1, 液相线温度随 Ni 含量变化
"""
import numpy as np
import matplotlib.pyplot as plt
from pycalphad import Database, equilibrium, variables as v
from pycalphad.model import Model
from pycalphad.model_uem_integrated import ModelUEM1


def calculate_liquidus_line(dbe, comps, phases, x_ni_range, model_dict, label):
	"""
	计算液相线温度

	Parameters
	----------
	dbe : Database
		热力学数据库
	comps : list
		组分列表
	phases : list
		相列表
	x_ni_range : array
		Ni 含量范围
	model_dict : dict
		模型字典
	label : str
		标签（用于绘图）

	Returns
	-------
	liquidus_temperatures : array
		液相线温度数组
	"""
	liquidus_temps = []

	for x_ni in x_ni_range:
		# xAl/xCr = 1/1, 即 xAl = xCr = (1 - x_ni) / 2
		x_al = (1.0 - x_ni) / 2.0
		x_cr = (1.0 - x_ni) / 2.0

		# 如果 Ni 含量过高，跳过
		if x_al < 0.01 or x_cr < 0.01:
			liquidus_temps.append(np.nan)
			continue

		# 设置温度范围，寻找液相线
		T_range = (1000, 2000, 20)  # 20K 步长，加快计算

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
			T_vals = eq.T.values
			phase_array = eq.Phase.values
			np_array = eq.NP.values

			liquidus_T = None

			# 从高温到低温扫描
			for t_idx in range(len(T_vals) - 1, -1, -1):
				# 检查所有 vertex
				for v_idx in range(phase_array.shape[-1]):
					phase_name = phase_array[0, 0, t_idx, 0, 0, v_idx]
					if phase_name == 'LIQUID':
						np_val = np_array[0, 0, t_idx, 0, 0, v_idx]
						if np_val > 0.995:  # 基本全是 LIQUID
							liquidus_T = T_vals[t_idx]
							break

				# 如果找到了液相线，停止搜索
				if liquidus_T is not None:
					break

			if liquidus_T is not None:
				liquidus_temps.append(liquidus_T)
				print(f"X(NI)={x_ni:.3f}, X(AL)={x_al:.3f}, X(CR)={x_cr:.3f} -> T_liquidus={liquidus_T:.1f} K ({label})")
			else:
				liquidus_temps.append(np.nan)
				print(f"X(NI)={x_ni:.3f}, X(AL)={x_al:.3f}, X(CR)={x_cr:.3f} -> 未找到液相线 ({label})")

		except Exception as e:
			print(f"X(NI)={x_ni:.3f} 计算出错: {e}")
			liquidus_temps.append(np.nan)

	return np.array(liquidus_temps)


def main():
	"""主测试函数"""
	print("=" * 70)
	print("Al-Cr-Ni 体系液相线温度计算")
	print("条件: xAl/xCr = 1/1")
	print("=" * 70)
	print()

	# 加载数据库
	dbe = Database('examples/alcrni.tdb')
	comps = ['AL', 'CR', 'NI']
	# 只使用主要的相以加快计算
	phases = ['LIQUID', 'FCC_A1', 'BCC_A2']

	# 定义 Ni 含量范围 - 减少点数以加快计算
	x_ni_range = np.linspace(0.1, 0.70, 7)  # 只计算 7 个点

	# 准备两种模型
	models_muggianu = {ph: Model for ph in phases}
	models_uem1 = {ph: ModelUEM1 for ph in phases}

	print("=" * 70)
	print("开始计算 Muggianu 模型的液相线")
	print("=" * 70)
	liquidus_muggianu = calculate_liquidus_line(
		dbe, comps, phases, x_ni_range, models_muggianu, "Muggianu"
	)

	print()
	print("=" * 70)
	print("开始计算 UEM1 模型的液相线")
	print("=" * 70)
	liquidus_uem1 = calculate_liquidus_line(
		dbe, comps, phases, x_ni_range, models_uem1, "UEM1"
	)

	# 绘制结果
	print()
	print("=" * 70)
	print("绘制液相线对比图")
	print("=" * 70)

	plt.figure(figsize=(10, 6))

	# 过滤掉 NaN 值
	valid_muggianu = ~np.isnan(liquidus_muggianu)
	valid_uem1 = ~np.isnan(liquidus_uem1)

	plt.plot(x_ni_range[valid_muggianu], liquidus_muggianu[valid_muggianu],
	         'o-', label='Muggianu', linewidth=2, markersize=6)
	plt.plot(x_ni_range[valid_uem1], liquidus_uem1[valid_uem1],
	         's-', label='UEM1', linewidth=2, markersize=6)

	plt.xlabel('X(Ni)', fontsize=12)
	plt.ylabel('Liquidus Temperature (K)', fontsize=12)
	plt.title('Al-Cr-Ni System Liquidus (xAl/xCr = 1/1)', fontsize=14, fontweight='bold')
	plt.legend(fontsize=11)
	plt.grid(True, alpha=0.3)
	plt.tight_layout()

	# 保存图片
	output_file = 'liquidus_alcrni_comparison.png'
	plt.savefig(output_file, dpi=150)
	print(f"图片已保存到: {output_file}")

	# 显示图片（如果在交互环境中）
	try:
		plt.show()
	except:
		pass

	# 打印统计信息
	print()
	print("=" * 70)
	print("统计信息")
	print("=" * 70)

	# 计算差异
	diff = liquidus_muggianu - liquidus_uem1
	valid_diff = ~np.isnan(diff)

	if np.any(valid_diff):
		print(f"液相线温度差异 (Muggianu - UEM1):")
		print(f"  平均值: {np.nanmean(diff):.2f} K")
		print(f"  最大值: {np.nanmax(diff):.2f} K")
		print(f"  最小值: {np.nanmin(diff):.2f} K")
		print(f"  标准差: {np.nanstd(diff):.2f} K")

	# 保存数据到文件
	output_data = 'liquidus_data.txt'
	with open(output_data, 'w') as f:
		f.write("# Al-Cr-Ni Liquidus Temperature Data (xAl/xCr = 1/1)\n")
		f.write("# X(Ni)\tX(Al)\tX(Cr)\tT_liquidus_Muggianu(K)\tT_liquidus_UEM1(K)\tDifference(K)\n")
		for i, x_ni in enumerate(x_ni_range):
			x_al = (1.0 - x_ni) / 2.0
			x_cr = (1.0 - x_ni) / 2.0
			t_mug = liquidus_muggianu[i]
			t_uem = liquidus_uem1[i]
			diff_val = t_mug - t_uem if not np.isnan(t_mug) and not np.isnan(t_uem) else np.nan
			f.write(f"{x_ni:.4f}\t{x_al:.4f}\t{x_cr:.4f}\t{t_mug:.2f}\t{t_uem:.2f}\t{diff_val:.2f}\n")

	print(f"\n数据已保存到: {output_data}")

	print()
	print("=" * 70)
	print("测试完成！")
	print("=" * 70)


if __name__ == '__main__':
	main()
