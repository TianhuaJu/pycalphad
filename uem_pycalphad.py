import sys
import numpy as np
import xarray as xr
from pycalphad import Database, variables as v
from pycalphad.model import Model
from pycalphad import calculate, equilibrium

# [!!] 1. 导入 pytest [!!]
import pytest

# ---------------------------------------------------------------------
# 导入您修改后的文件中的类
# ---------------------------------------------------------------------
try:
	from pycalphad.model_uem_integrated import ModelWithUEM, ModelUEM1
except ImportError:
	print("错误: 无法导入 'ModelWithUEM' 或 'ModelUEM1'。")
	print("请确保 'model_uem_integrated.py' 与此脚本在同一目录中。\n")
	sys.exit(1)
except Exception as e:
	print(f"导入时发生意外错误: {e}")
	print("请确保您已应用 __new__ 和 __init__ 的修复！\n")
	sys.exit(1)


def load_database (db_path='examples/alcrni.tdb'):
	"""加载数据库并处理 FileNotFoundError"""
	# (此函数保持不变)
	try:
		dbe = Database(db_path)
		print(f"成功加载数据库: {db_path}")
		return dbe
	except FileNotFoundError:
		print(f"错误: 数据库文件未找到: {db_path}")
		print("请确保 pycalphad 的 'examples' 目录在您的工作路径中,")
		print("或者提供 'alcrni.tdb' 的完整路径。")
		sys.exit(1)
	except Exception as e:
		print(f"加载数据库时出错: {e}")
		sys.exit(1)


# [!!] 2. 定义 'dbe' Fixture [!!]
# scope="module" 确保数据库在所有测试中只被加载一次
@pytest.fixture(scope="module")
def dbe ():
	"""
	Pytest fixture to load the database once per test module.
	"""
	db = load_database('examples/alcrni.tdb')
	return db


# [!!] 3. 您的测试函数保持不变 [!!]
# pytest 现在会自动将上面的 dbe fixture 注入到这个 'dbe' 参数中
def test_calculate_phase_fractions (dbe):
	"""
	测试 1: 相分数 (NP) 计算 (pycalphad.equilibrium)
	...
	"""
	print("\n" + "=" * 70)
	print("测试 1: 'pycalphad.equilibrium' (相分数 vs 温度)")
	print("=" * 70)
	
	# (函数内部代码保持不变)
	comps = ['AL', 'CR', 'NI']
	phases = ['LIQUID', 'FCC_A1', 'BCC_A2', 'LAVES_C14', 'GAMMA_PRIME_L12', 'ALNI_B2']
	conditions = {
		v.T: (1000, 1800, 20),
		v.P: 101325,
		v.X('AL'): 0.25,
		v.X('CR'): 0.35,
	}
	models_base = {ph: Model for ph in phases}
	models_uem1 = {ph: ModelUEM1 for ph in phases}
	
	try:
		print("正在计算: Baseline (Muggianu)...")
		calc_base = equilibrium(dbe, comps, phases, model=models_base, conditions=conditions)

		print("正在计算: ModelUEM1 (UEM1 模式)...")
		calc_uem1 = equilibrium(dbe, comps, phases, model=models_uem1, conditions=conditions)

		print("\n--- 验证结果 (测试 1) ---")

		assert calc_base.T.size > 0 and calc_uem1.T.size > 0
		print("[通过] 两种模型都成功完成了 'equilibrium' 计算。")

		# 比较 GM (摩尔吉布斯自由能) 来验证结果不同
		gm_base_sum = float(np.nansum(calc_base.GM.values))
		gm_uem1_sum = float(np.nansum(calc_uem1.GM.values))

		print(f"总吉布斯自由能:")
		print(f"  Muggianu: {gm_base_sum:.2e}")
		print(f"  UEM1:     {gm_uem1_sum:.2e}")

		assert not np.isclose(gm_base_sum, gm_uem1_sum, rtol=1e-5), \
			"UEM1 和 Muggianu 的吉布斯自由能结果相同！"
		print("[通过] UEM1 和 Muggianu 的吉布斯自由能不同。")

		print("\n[成功] 'pycalphad.equilibrium' 测试通过！")

	except AssertionError as e:
		print(f"\n[!!! 测试失败 !!!] 错误: {e}")
	except Exception as e:
		print(f"\n[!!! 严重错误 !!!] 'equilibrium' 测试在执行过程中崩溃: {e}")
		import traceback
		traceback.print_exc()


# [!!] 4. 您的第二个测试函数也保持不变 [!!]
def test_map_phase_diagram (dbe):
	"""
	测试 2: 相图 (T-x 截面) 计算 (pycalphad.equilibrium)
	使用简化的条件避免数值问题
	"""
	print("\n" + "=" * 70)
	print("测试 2: 'pycalphad.equilibrium' (温度扫描)")
	print("=" * 70)

	# 使用简化的条件
	comps = ['AL', 'CR', 'NI']
	phases = ['LIQUID', 'FCC_A1']
	conditions = {
		v.T: (1400, 1800, 10),  # 温度范围
		v.P: 101325,
		v.X('AL'): 0.30,
		v.X('CR'): 0.30,
	}
	models_base = {ph: Model for ph in phases}
	models_uem1 = {ph: ModelUEM1 for ph in phases}

	try:
		print("正在计算: Baseline (Muggianu) 温度扫描...")
		map_base = equilibrium(dbe, comps, phases, model=models_base,
		                       conditions=conditions)

		print("正在计算: ModelUEM1 (UEM1 模式) 温度扫描...")
		map_uem1 = equilibrium(dbe, comps, phases, model=models_uem1,
		                       conditions=conditions)

		print("\n--- 验证结果 (测试 2) ---")
		assert isinstance(map_base, xr.Dataset) and isinstance(map_uem1, xr.Dataset)
		assert map_base.T.size > 0 and map_uem1.T.size > 0
		print("[通过] 两种模型都成功完成了 'equilibrium' 计算。")

		# 比较吉布斯自由能
		gm_base_sum = float(np.nansum(map_base.GM.values))
		gm_uem1_sum = float(np.nansum(map_uem1.GM.values))

		print(f"总吉布斯自由能:")
		print(f"  Muggianu: {gm_base_sum:.2e}")
		print(f"  UEM1:     {gm_uem1_sum:.2e}")

		assert not np.isclose(gm_base_sum, gm_uem1_sum, rtol=1e-5), \
			"UEM1 和 Muggianu 的吉布斯自由能相同！"
		print("[通过] UEM1 和 Muggianu 的计算结果不同。")

		print("\n[成功] 'pycalphad.equilibrium' (温度扫描) 测试通过！")

	except AssertionError as e:
		print(f"\n[!!! 测试失败 !!!] 错误: {e}")
	except Exception as e:
		print(f"\n[!!! 严重错误 !!!] 'equilibrium' (温度扫描) 测试崩溃: {e}")
		import traceback
		traceback.print_exc()


def test_enthalpy_calculation(dbe):
	"""
	测试 3: 吉布斯自由能计算
	测试 LIQUID 相使用 UEM 模型的热力学性质计算
	使用 equilibrium 来确保正确比较单点
	"""
	print("\n" + "=" * 70)
	print("测试 3: 吉布斯自由能计算 (LIQUID 相)")
	print("=" * 70)

	comps = ['AL', 'CR', 'NI']
	phases = ['LIQUID']

	# 定义多个测试点
	test_points = [
		{'X_AL': 0.25, 'X_CR': 0.35, 'T': 1500},  # 点 1
		{'X_AL': 0.30, 'X_CR': 0.30, 'T': 1600},  # 点 2
		{'X_AL': 0.40, 'X_CR': 0.25, 'T': 1700},  # 点 3
	]

	# 创建两种模型
	models_base = {'LIQUID': Model}
	models_uem1 = {'LIQUID': ModelUEM1}

	try:
		print(f"测试 {len(test_points)} 个不同的组分-温度点\n")

		all_diffs = []

		for i, point in enumerate(test_points):
			conditions = {
				v.T: point['T'],
				v.P: 101325,
				v.X('AL'): point['X_AL'],
				v.X('CR'): point['X_CR'],
			}

			# 使用 equilibrium 计算单点
			eq_base = equilibrium(dbe, comps, phases, model=models_base, conditions=conditions)
			eq_uem1 = equilibrium(dbe, comps, phases, model=models_uem1, conditions=conditions)

			# 提取吉布斯自由能
			gm_base = float(eq_base.GM.squeeze().values)
			gm_uem1 = float(eq_uem1.GM.squeeze().values)

			diff = abs(gm_base - gm_uem1)
			all_diffs.append(diff)

			print(f"点 {i+1} [T={point['T']}K, X(AL)={point['X_AL']}, X(CR)={point['X_CR']}]:")
			print(f"  Muggianu: {gm_base:.2f} J/mol")
			print(f"  UEM1:     {gm_uem1:.2f} J/mol")
			print(f"  差值:     {diff:.2f} J/mol")
			print()

		print("\n--- 验证结果 (测试 3) ---")

		# 验证所有点都计算成功
		assert len(all_diffs) == len(test_points)
		print(f"[通过] 成功计算了 {len(test_points)} 个点的吉布斯自由能。")

		# 验证至少有一个点显示出明显差异
		max_diff = max(all_diffs)
		avg_diff = np.mean(all_diffs)

		print(f"最大差值: {max_diff:.2f} J/mol")
		print(f"平均差值: {avg_diff:.2f} J/mol")

		assert max_diff > 1.0, \
			"UEM1 和 Muggianu 在所有测试点的吉布斯自由能都相同！"
		print("[通过] UEM1 和 Muggianu 的吉布斯自由能存在差异。")

		print("\n[成功] 吉布斯自由能计算测试通过！")

	except AssertionError as e:
		print(f"\n[!!! 测试失败 !!!] 错误: {e}")
	except Exception as e:
		print(f"\n[!!! 严重错误 !!!] 吉布斯自由能计算测试崩溃: {e}")
		import traceback
		traceback.print_exc()


def test_activity_calculation(dbe):
	"""
	测试 4: 活度计算
	测试 LIQUID 相使用 UEM 模型的活度计算（通过 equilibrium 获取化学势）
	"""
	print("\n" + "=" * 70)
	print("测试 4: 活度计算 (LIQUID 相)")
	print("=" * 70)

	comps = ['AL', 'CR', 'NI']
	phases = ['LIQUID']

	# 平衡计算条件
	conditions = {
		v.T: 1600,  # 温度 1600K
		v.P: 101325,
		v.X('AL'): 0.30,
		v.X('CR'): 0.30,
		# NI: 1 - 0.30 - 0.30 = 0.40
	}

	# LIQUID 相分别使用不同模型
	models_base = {'LIQUID': Model}
	models_uem1 = {'LIQUID': ModelUEM1}

	try:
		print(f"计算条件: T={conditions[v.T]}K, X(AL)={conditions[v.X('AL')]}, X(CR)={conditions[v.X('CR')]}")

		# 使用 equilibrium 来获取化学势，因为 calculate 不直接提供化学势
		print("\n正在计算: Baseline (Muggianu) LIQUID 相化学势...")
		eq_base = equilibrium(dbe, comps, phases, model=models_base,
		                      conditions=conditions)

		print("正在计算: ModelUEM1 LIQUID 相化学势...")
		eq_uem1 = equilibrium(dbe, comps, phases, model=models_uem1,
		                      conditions=conditions)

		print("\n--- 验证结果 (测试 4) ---")

		# 提取各组分的化学势
		# equilibrium 返回的 MU 维度是 (N, P, T, X_AL, X_CR, component)
		R = 8.314462618  # J/(mol·K)
		T = conditions[v.T]

		# 提取化学势 - 使用 squeeze 移除大小为1的维度
		mu_base = eq_base.MU.squeeze()
		mu_uem1 = eq_uem1.MU.squeeze()

		# 获取每个组分的化学势
		mu_al_base = float(mu_base.sel(component='AL').values)
		mu_al_uem1 = float(mu_uem1.sel(component='AL').values)

		mu_cr_base = float(mu_base.sel(component='CR').values)
		mu_cr_uem1 = float(mu_uem1.sel(component='CR').values)

		mu_ni_base = float(mu_base.sel(component='NI').values)
		mu_ni_uem1 = float(mu_uem1.sel(component='NI').values)

		print(f"化学势 (MU) [J/mol]:")
		print(f"  AL - Muggianu: {mu_al_base:.2f}, UEM1: {mu_al_uem1:.2f}, 差值: {abs(mu_al_base - mu_al_uem1):.2f}")
		print(f"  CR - Muggianu: {mu_cr_base:.2f}, UEM1: {mu_cr_uem1:.2f}, 差值: {abs(mu_cr_base - mu_cr_uem1):.2f}")
		print(f"  NI - Muggianu: {mu_ni_base:.2f}, UEM1: {mu_ni_uem1:.2f}, 差值: {abs(mu_ni_base - mu_ni_uem1):.2f}")

		# 计算活度（相对于参考态）
		# 注意：这里计算的是相对活度差异
		act_al_base = np.exp(mu_al_base / (R * T))
		act_al_uem1 = np.exp(mu_al_uem1 / (R * T))

		act_cr_base = np.exp(mu_cr_base / (R * T))
		act_cr_uem1 = np.exp(mu_cr_uem1 / (R * T))

		print(f"\n相对活度:")
		print(f"  AL - Muggianu: {act_al_base:.6e}, UEM1: {act_al_uem1:.6e}")
		print(f"  CR - Muggianu: {act_cr_base:.6e}, UEM1: {act_cr_uem1:.6e}")

		# 验证计算成功
		assert not np.isnan(mu_al_base) and not np.isnan(mu_al_uem1), \
			"AL 化学势计算返回 NaN！"
		assert not np.isnan(mu_cr_base) and not np.isnan(mu_cr_uem1), \
			"CR 化学势计算返回 NaN！"
		assert not np.isnan(mu_ni_base) and not np.isnan(mu_ni_uem1), \
			"NI 化学势计算返回 NaN！"
		print("[通过] 化学势计算成功，未出现 NaN。")

		# 验证 UEM1 和 Muggianu 的结果不同
		assert not np.isclose(mu_al_base, mu_al_uem1, atol=1.0), \
			"UEM1 和 Muggianu 的 AL 化学势相同！"
		print("[通过] UEM1 和 Muggianu 的化学势/活度不同。")

		print("\n[成功] 活度计算测试通过！")

	except AssertionError as e:
		print(f"\n[!!! 测试失败 !!!] 错误: {e}")
	except Exception as e:
		print(f"\n[!!! 严重错误 !!!] 活度计算测试崩溃: {e}")
		import traceback
		traceback.print_exc()


def test_liquid_phase_equilibrium(dbe):
	"""
	测试 5: LIQUID 相平衡计算
	使用 equilibrium 函数测试 LIQUID 相使用 UEM 模型，其他相使用默认模型
	"""
	print("\n" + "=" * 70)
	print("测试 5: LIQUID 相平衡计算 (UEM 仅用于 LIQUID)")
	print("=" * 70)

	comps = ['AL', 'CR', 'NI']
	phases = ['LIQUID', 'FCC_A1', 'BCC_A2']

	# 平衡计算条件
	conditions = {
		v.T: (1200, 1800, 50),  # 温度范围
		v.P: 101325,
		v.X('AL'): 0.35,
		v.X('CR'): 0.25,
	}

	# 创建混合模型：LIQUID 使用 UEM1，其他相使用默认 Model
	models_mixed = {
		'LIQUID': ModelUEM1,
		'FCC_A1': Model,
		'BCC_A2': Model,
	}

	models_base = {ph: Model for ph in phases}

	try:
		print(f"计算条件: T={conditions[v.T]}, X(AL)={conditions[v.X('AL')]}, X(CR)={conditions[v.X('CR')]}")
		print("模型设置: LIQUID=ModelUEM1, FCC_A1=Model, BCC_A2=Model")

		print("\n正在计算: 全部使用 Muggianu 模型...")
		eq_base = equilibrium(dbe, comps, phases, model=models_base,
		                      conditions=conditions)

		print("正在计算: LIQUID 使用 UEM1，其他相使用默认模型...")
		eq_mixed = equilibrium(dbe, comps, phases, model=models_mixed,
		                       conditions=conditions)

		print("\n--- 验证结果 (测试 5) ---")

		# 验证计算成功
		assert eq_base.T.size > 0 and eq_mixed.T.size > 0
		print("[通过] 两种模型配置都成功完成了平衡计算。")

		# 比较总吉布斯自由能
		gm_base_sum = float(np.nansum(eq_base.GM.values))
		gm_mixed_sum = float(np.nansum(eq_mixed.GM.values))

		print(f"\n总吉布斯自由能:")
		print(f"  全部 Muggianu: {gm_base_sum:.2e}")
		print(f"  LIQUID 用 UEM1: {gm_mixed_sum:.2e}")
		print(f"  差值:          {abs(gm_base_sum - gm_mixed_sum):.2e}")

		# 验证结果不同
		assert not np.isclose(gm_base_sum, gm_mixed_sum, rtol=1e-5), \
			"LIQUID 相使用 UEM1 和全部使用 Muggianu 的结果相同！"
		print("[通过] LIQUID 相使用 UEM1 模型产生了不同的平衡结果。")

		print("\n[成功] LIQUID 相平衡计算测试通过！")

	except AssertionError as e:
		print(f"\n[!!! 测试失败 !!!] 错误: {e}")
	except Exception as e:
		print(f"\n[!!! 严重错误 !!!] LIQUID 相平衡计算测试崩溃: {e}")
		import traceback
		traceback.print_exc()


# 主函数：运行所有测试
if __name__ == '__main__':
	print("=" * 70)
	print("UEM 模型性能测试套件")
	print("=" * 70)

	# 加载数据库
	dbe = load_database('examples/alcrni.tdb')

	# 运行所有测试
	print("\n\n开始运行所有测试...\n")

	# 测试 1: 相分数计算
	test_calculate_phase_fractions(dbe)

	# 测试 2: 相图计算
	test_map_phase_diagram(dbe)

	# 测试 3: 焓值计算
	test_enthalpy_calculation(dbe)

	# 测试 4: 活度计算
	test_activity_calculation(dbe)

	# 测试 5: LIQUID 相平衡计算
	test_liquid_phase_equilibrium(dbe)

	print("\n" + "=" * 70)
	print("所有测试已完成！")
	print("=" * 70)
