#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
测试脚本：
使用自定义的 UEMModel (Unified Extrapolation Model)
计算 Al-Cr-Ni 系统在 xAl/xCr=1 条件下的液相线温度。
"""

# --- 导入 ---
import numpy as np
import matplotlib.pyplot as plt
from pycalphad import Database, equilibrium, variables as v
from pycalphad.core.errors import DofError
import logging
from pycalphad.UEMModel import UEMModel

# 导入 pycalphad 的基类 Model 和 UEM 所需的库
try:
	from pycalphad.model import Model
	from itertools import combinations
	from tinydb import where
	from symengine import exp, Add, Piecewise, S, Mul, Pow
except ImportError:
	print("错误: 无法导入 'pycalphad.model'。")
	print("请确保已安装 pycalphad (pip install pycalphad)")
	exit()

# 禁用 pycalphad 的冗余日志
logging.getLogger('pycalphad').setLevel(logging.WARNING)

# ---------------------------------------------------------------------------
# --- 主测试函数 (与您的文件 相同) ---
# ---------------------------------------------------------------------------

def calculate_liquidus_with_uem ():
	"""
	使用 UEMModel 计算 Al-Cr-Ni 液相线的主函数
	"""
	
	# 1. 加载数据库
	print("正在加载 'examples/alcrni.tdb' 数据库...")
	try:
		dbf = Database('examples/alcrni.tdb')
		print("数据库加载成功。")
	except Exception as e:
		print(f"错误: 无法加载数据库: {e}")
		print("请确保已安装 'pycalphad[examples]' (例如: pip install pycalphad[examples])")
		return
	
	comps = ['AL', 'CR', 'NI', 'VA']
	all_phase_names = list(dbf.phases.keys())
	
	# 2. **关键步骤: 为所有相实例化 UEMModel**
	print("正在为所有相关相构建 UEMModel 实例...")
	uem_models = {}
	for phase_name in all_phase_names:
		try:
			# 尝试为该相构建 UEMModel
			# (现在 UEMModel 已在上方定义)
			uem_models[phase_name] = UEMModel(dbf, comps, phase_name)
		except DofError:
			# 如果该相不能由 'AL', 'CR', 'NI' 构成，
			# DofError 会被触发。我们忽略这个相。
			pass
		except Exception as e:
			print(f"  警告: 构建相 {phase_name} 时出错: {e}")
	
	print(f"已成功为 {len(uem_models)} 个相构建 UEMModel。")
	print(f"将使用 UEMModel 计算的相: {list(uem_models.keys())}")
	
	# 3. 设置成分条件
	# 减少点数以加快计算速度（从40减少到10）
	num_points = 10
	xNi_values = np.linspace(0.1, 0.9, num_points)
	xAl_values = (1.0 - xNi_values) / 2.0
	
	liquidus_temps = []
	calculated_xNi = []
	
	print(f"\n开始使用 UEMModel 计算 {num_points} 个成分点的液相线温度...")
	print("这可能需要几分钟时间...")
	
	# 4. 循环计算每个成分点
	for xni, xal in zip(xNi_values, xAl_values):
		
		conditions = {
			v.P: 101325,
			v.T: (1400, 2400, 10),  # 温度范围 (K)，扩展范围以覆盖高熔点成分
			v.X('NI'): xni,
			v.X('AL'): xal
		}
		
		try:
			# **运行平衡计算，并传入 uem_models 字典**
			eq_result = equilibrium(dbf, comps, list(uem_models.keys()),
			                        conditions,
			                        model=uem_models,  # <-- 关键！
			                        output='NP')

			# 5. 提取液相线温度
			# 从高温到低温扫描，找到100%液相的最低温度
			T_vals = eq_result.T.values
			phase_array = eq_result.Phase.values
			np_array = eq_result.NP.values

			liquidus_T = None

			# 从高温到低温扫描
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
				calculated_xNi.append(xni)
				print(f"  [点 {len(calculated_xNi)}/{num_points}] 成功: xNi={xni:.3f}, T_liq={liquidus_T:.2f} K")
			else:
				print(f"  [点 {len(calculated_xNi) + 1}/{num_points}] 警告: xNi={xni:.3f} - 在指定温度范围内未找到100%液相区。")
		
		except Exception as e:
			print(f"  [点 {len(calculated_xNi) + 1}/{num_points}] 错误: xNi={xni:.3f} 计算失败: {e}")
	
	print("\n计算完成。")
	
	# 6. 绘图
	if liquidus_temps:
		print("正在绘制结果图...")
		plt.figure(figsize=(10, 6))
		plt.plot(calculated_xNi, liquidus_temps, 'r.-', label='Liquidus (UEMModel)')
		plt.title(f'Al-Cr-Ni Liquidus (Section xAl/xCr = 1)\nCALCULATED WITH UEMMODEL')
		plt.xlabel('Mole Fraction of Ni ($x_{Ni}$)')
		plt.ylabel('Liquidus Temperature (K)')
		plt.xlim(0, 1)
		plt.grid(True)
		plt.legend()
		
		plot_filename = 'alcrni_liquidus_UEM_plot.png'
		plt.savefig(plot_filename)
		print(f"图像已保存到: {plot_filename}")
	
	else:
		print("没有可供绘图的数据。")


if __name__ == "__main__":
	calculate_liquidus_with_uem()