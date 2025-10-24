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

# 导入 pycalphad 的基类 Model 和 UEM 所需的库
try:
	from pycalphad.models.model import Model
	from itertools import combinations
	from tinydb import where
	from symengine import exp, Add, Piecewise, S, Mul, Pow
except ImportError:
	print("错误: 无法导入 'pycalphad.models.model'。")
	print("请确保已安装 pycalphad (pip install pycalphad)")
	exit()

# 禁用 pycalphad 的冗余日志
logging.getLogger('pycalphad').setLevel(logging.WARNING)


# ---------------------------------------------------------------------------
# --- UEMModel 类的完整定义 (从上一问中复制) ---
# ---------------------------------------------------------------------------
class UEMModel(Model):
	"""
	继承自 pycalphad.Model，
	使用 UEM (Unified Extrapolation Model) 逻辑
	重写 excess_mixing_energy。

	该模型假设多元系的过剩能完全由二元系外推得到，
	并将忽略数据库中所有的三元及更高阶交互参数。
	"""
	
	def _build_binary_L_cache (self, dbe):
		"""
		辅助函数：搜索并缓存所有单亚晶格二元交互参数。
		"""
		self._binary_L_param_cache = {}
		param_search = dbe.search
		phase = dbe.phases[self.phase_name]
		
		#
		param_query = (
				(where('phase_name') == self.phase_name) &
				((where('parameter_type') == 'G') | (where('parameter_type') == 'L')) &
				(where('constituent_array').test(self._interaction_test))  #
		)
		params = param_search(param_query)
		
		for param in params:
			mixing_subl_indices = []
			is_binary = True
			target_subl = -1
			comps = []
			
			for subl_idx, const_array in enumerate(param['constituent_array']):
				if len(const_array) > 2:
					is_binary = False
					break
				if len(const_array) == 2:
					mixing_subl_indices.append(subl_idx)
					target_subl = subl_idx
					comps = tuple(sorted([v.Species(c) for c in const_array]))
				elif len(const_array) == 1 and const_array[0] == v.Species('*'):
					pass
				elif len(const_array) == 1:
					pass
				else:
					is_binary = False
					break
			
			#
			if not is_binary or len(mixing_subl_indices) != 1:
				continue
			
			key = (target_subl, comps)
			if key not in self._binary_L_param_cache:
				self._binary_L_param_cache[key] = {'all': [], 'odd': []}
			
			L_expr = param['parameter']
			order = param.get('parameter_order', 0)
			
			self._binary_L_param_cache[key]['all'].append((L_expr, order))
			
			if order % 2 != 0:
				self._binary_L_param_cache[key]['odd'].append(L_expr)
	
	def _get_uem_d_term (self, dbe, comp_k, comp_i, subl_index):
		"""
		根据UEM公式计算 d_ki 项。
		d_ki = (2/RT) * sum(L_ki^{(v)}) (v 为奇数)
		"""
		
		if not hasattr(self, '_binary_L_param_cache'):
			self._build_binary_L_cache(dbe)
		
		key = (subl_index, tuple(sorted((comp_k, comp_i))))
		odd_L_terms = self._binary_L_param_cache.get(key, {}).get('odd', [])
		
		if not odd_L_terms:
			return S.Zero
		
		d_ki = (S(2) / (v.R * v.T)) * Add(*odd_L_terms)
		return d_ki
	
	def excess_mixing_energy (self, dbe):
		"""
		重写基类的 excess_mixing_energy。

		该实现完全基于 UEM 公式：
		G^E = sum [ (x_i * x_j) * rk_sum(X_i, X_j) ]

		注意：这会忽略所有三元及更高阶的参数！
		"""
		
		if not hasattr(self, '_binary_L_param_cache'):
			self._build_binary_L_cache(dbe)
		
		total_excess_energy = S.Zero
		phase = dbe.phases[self.phase_name]
		
		#
		for subl_index, sublattice_comps in enumerate(self.constituents):
			#
			active_comps = sorted(list(sublattice_comps.intersection(self.components)))
			
			if len(active_comps) < 2:
				continue
			
			#
			site_fracs = {comp: v.SiteFraction(self.phase_name, subl_index, comp) for comp in active_comps}
			
			sublattice_total_G_E = S.Zero
			
			for comp_i, comp_j in combinations(active_comps, 2):
				x_i = site_fracs[comp_i]
				x_j = site_fracs[comp_j]
				
				G_ij_binary_expr_rk_sum = S.Zero
				key = (subl_index, tuple(sorted((comp_i, comp_j))))
				binary_params = self._binary_L_param_cache.get(key, {}).get('all', [])
				
				if not binary_params:
					continue
				
				p_i_species, p_j_species = key[1]
				p_i = site_fracs[p_i_species]
				p_j = site_fracs[p_j_species]
				
				for L_expr, order in binary_params:
					#
					G_ij_binary_expr_rk_sum += L_expr * (p_i - p_j) ** order
				
				X_ij_i_num = x_i
				X_ij_j_num = x_j
				X_denom = x_i + x_j
				
				other_comps = [c for c in active_comps if c != comp_i and c != comp_j]
				
				if not other_comps:
					sublattice_total_G_E += p_i * p_j * G_ij_binary_expr_rk_sum
					continue
				
				for comp_k in other_comps:
					x_k = site_fracs[comp_k]
					
					d_ki = self._get_uem_d_term(dbe, comp_k, comp_i, subl_index)
					d_kj = self._get_uem_d_term(dbe, comp_k, comp_j, subl_index)
					
					d_sum = d_ki + d_kj
					d_sum_safe = Piecewise((d_sum, d_sum != 0), (1, True))
					
					alpha_i_k = (d_kj / d_sum_safe) * exp(-d_ki)
					alpha_j_k = (d_ki / d_sum_safe) * exp(-d_kj)
					
					X_ij_i_num += alpha_i_k * x_k
					X_ij_j_num += alpha_j_k * x_k
					X_denom += (alpha_i_k + alpha_j_k) * x_k
				
				X_denom_safe = Piecewise((X_denom, X_denom != 0), (1, True))
				X_ij_i = X_ij_i_num / X_denom_safe
				X_ij_j = X_ij_j_num / X_denom_safe
				
				substitution_dict = {}
				if p_i_species == comp_i:
					substitution_dict[p_i] = X_ij_i
					substitution_dict[p_j] = X_ij_j
				else:
					substitution_dict[p_i] = X_ij_j
					substitution_dict[p_j] = X_ij_i
				
				rk_sum_modified = G_ij_binary_expr_rk_sum.xreplace(substitution_dict)
				
				final_term = (x_i * x_j) * rk_sum_modified
				sublattice_total_G_E += final_term
			
			total_excess_energy += sublattice_total_G_E
		
		#
		return total_excess_energy / self._site_ratio_normalization


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
	num_points = 40
	xNi_values = np.linspace(0.01, 0.99, num_points)
	xAl_values = (1.0 - xNi_values) / 2.0
	
	liquidus_temps = []
	calculated_xNi = []
	
	print(f"\n开始使用 UEMModel 计算 {num_points} 个成分点的液相线温度...")
	print("这可能需要几分钟时间...")
	
	# 4. 循环计算每个成分点
	for xni, xal in zip(xNi_values, xAl_values):
		
		conditions = {
			v.P: 101325,
			v.T: (1300, 3000, 10),  # 温度范围 (K)
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
			liquid_phase_fraction = eq_result.NP.sel(phase='LIQUID')
			all_liquid_temps = eq_result.T.where(liquid_phase_fraction > 0.999, drop=True)
			
			if all_liquid_temps.size > 0:
				liq_t = float(all_liquid_temps.min())
				liquidus_temps.append(liq_t)
				calculated_xNi.append(xni)
				print(f"  [点 {len(calculated_xNi)}/{num_points}] 成功: xNi={xni:.3f}, T_liq={liq_t:.2f} K")
			else:
				print(
						f"  [点 {len(calculated_xNi) + 1}/{num_points}] 警告: xNi={xni:.3f} - 在指定温度范围内未找到100%液相区。")
		
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