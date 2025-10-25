# ---------------------------------------------------------------------------
# --- UEM Model Implementation (Start) ---
# ---------------------------------------------------------------------------
#
# 您需要将此代码粘贴到 model.py 文件的末尾，
# 确保在文件顶部导入了 necessary的模块:
#
# from symengine import exp, log, Abs, Add, And, Float, Mul, Piecewise, Pow, S, sin, StrictGreaterThan, Symbol, zoo, oo
# import pycalphad.variables as v
# from tinydb import where
# from itertools import combinations # <-- 确保导入
#
# ---------------------------------------------------------------------------

from itertools import combinations
from tinydb import where
from symengine import exp, Add, Piecewise, S, Mul, Pow, Not, StrictGreaterThan
import pycalphad.variables as v
from pycalphad.model import  Model


class uem1_model(Model):
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

		缓存格式:
		self._binary_L_param_cache[ (subl_index, (comp_i, comp_j_sorted)) ] = {
			'all': [(L0_expr, 0), (L1_expr, 1), ...],
			'odd': [L1_expr, L3_expr, ...],
			'first_comp_in_db': comp_i  # 记录数据库中的第一个组分
		}
		"""
		self._binary_L_param_cache = {}
		param_search = dbe.search
		phase = dbe.phases[self.phase_name]

		# 查询所有 'G' 或 'L' 交互参数
		#
		param_query = (
				(where('phase_name') == self.phase_name) &
				((where('parameter_type') == 'G') | (where('parameter_type') == 'L')) &
				(where('constituent_array').test(self._interaction_test))
		)
		params = param_search(param_query)

		for param in params:
			# 检查这是否是一个单亚晶格二元交互
			mixing_subl_indices = []
			is_binary = True
			target_subl = -1
			comps_original = None

			for subl_idx, const_array in enumerate(param['constituent_array']):
				if len(const_array) > 2:  # 超过2元，不是二元交互
					is_binary = False
					break
				if len(const_array) == 2:
					mixing_subl_indices.append(subl_idx)
					target_subl = subl_idx
					# 保留原始顺序（不使用 sorted）
					comps_original = tuple([v.Species(c) for c in const_array])
				elif len(const_array) == 1 and const_array[0] == v.Species('*'):
					pass  # 允许通配符
				elif len(const_array) == 1:
					pass  # 允许单个固定组元
				else:
					is_binary = False
					break

			# 我们只关心在单个亚晶格上混合的二元交互
			if not is_binary or len(mixing_subl_indices) != 1 or comps_original is None:
				continue

			# 使用规范化的 key（sorted）用于查找，但记录原始顺序
			key_canonical = (target_subl, tuple(sorted(comps_original)))

			if key_canonical not in self._binary_L_param_cache:
				self._binary_L_param_cache[key_canonical] = {
					'all': [],
					'odd': [],
					'first_comp_in_db': comps_original[0]  # 记录数据库中的第一个组分
				}

			L_expr = param['parameter']
			order = param.get('parameter_order', 0)

			self._binary_L_param_cache[key_canonical]['all'].append((L_expr, order))

			if order % 2 != 0:  # 奇数阶
				self._binary_L_param_cache[key_canonical]['odd'].append(L_expr)
	
	def _get_uem_d_term (self, dbe, comp_k, comp_i, subl_index):
		"""
		根据UEM公式计算 d_ki 项。
		d_ki = (1/RT) * |g_i^∞(in k) - g_k^∞(in i)|

		对于Redlich-Kister模型，这简化为:
		d_ki = (2/RT) * |sum(L_ki^{(v)})| (v 为奇数)

		关键：
		1. L1(A,B) = -L1(B,A)，需要根据数据库顺序调整符号
		2. 必须取绝对值，确保 d_ki >= 0
		"""

		# 确保缓存已建立
		if not hasattr(self, '_binary_L_param_cache'):
			self._build_binary_L_cache(dbe)

		# 使用规范化的 key 查找
		key_canonical = (subl_index, tuple(sorted((comp_k, comp_i))))

		cache_entry = self._binary_L_param_cache.get(key_canonical, None)
		if cache_entry is None:
			return S.Zero

		odd_L_terms = cache_entry.get('odd', [])
		if not odd_L_terms:
			return S.Zero

		# 根据数据库中的顺序判断是否需要符号修正
		first_comp_in_db = cache_entry['first_comp_in_db']

		# 判断逻辑：
		# - 如果数据库中是 (k, i) 顺序，我们需要 (k, i)，sign = +1
		# - 如果数据库中是 (i, k) 顺序，我们需要 (k, i)，sign = -1
		#   因为 L1(k, i) = -L1(i, k)
		if first_comp_in_db == comp_k:
			sign_correction = 1
		elif first_comp_in_db == comp_i:
			sign_correction = -1
		else:
			# 不应该发生这种情况
			raise ValueError(f"Unexpected first component: {first_comp_in_db}")
		
		term_sum = sign_correction * Add(*odd_L_terms)
		d_ki_signed = Piecewise(
				(Mul(S(2) / (v.R * v.T), term_sum), StrictGreaterThan(v.T, S.Zero)),  # T > 0
				(S.Zero, True)  # T <= 0 (或 T 不是符号)
		)
		d_ki = Piecewise(
				(d_ki_signed, Not(StrictGreaterThan(S.Zero, d_ki_signed))),  # 条件: d_ki_signed >= 0
				(-d_ki_signed, True)  # 其他情况: d_ki_signed < 0
		)
		
		return d_ki
	
	def excess_mixing_energy (self, dbe):
		"""
		重写基类的 excess_mixing_energy。

		该实现完全基于 UEM 公式：
		G^E = sum_{i<j} [ (x_i x_j) / (X_ij^i X_ij^j) ] * G_ij^E(X_ij^i, X_ij^j)

		注意：这里完全忽略所有三元及更高阶的参数！
		"""
		
		# 1. 建立二元 L 参数的缓存
		if not hasattr(self, '_binary_L_param_cache'):
			self._build_binary_L_cache(dbe)
		
		total_excess_energy = S.Zero
		phase = dbe.phases[self.phase_name]
		
		# 2. 遍历所有亚晶格
		for subl_index, sublattice_comps in enumerate(self.constituents):
			active_comps = sorted(list(sublattice_comps.intersection(self.components)))
			
			if len(active_comps) < 2:
				continue  # 该亚晶格没有混合
			
			# 获取该亚晶格的位组分数 (Site Fractions)
			site_fracs = {comp: v.SiteFraction(self.phase_name, subl_index, comp) for comp in active_comps}
			
			sublattice_total_G_E = S.Zero
			
			# 3. 遍历所有二元组 (i, j)
			for comp_i, comp_j in combinations(active_comps, 2):
				x_i = site_fracs[comp_i]
				x_j = site_fracs[comp_j]
				
				# 4. 构建 G_ij^E(x_i, x_j) 二元 R-K 表达式
				G_ij_binary_expr_rk_sum = S.Zero
				key = (subl_index, tuple(sorted((comp_i, comp_j))))
				binary_params = self._binary_L_param_cache.get(key, {}).get('all', [])
				
				if not binary_params:
					continue  # 没有 (i, j) 交互
				
				# 获取参数定义时的组元顺序
				p_i_species, p_j_species = key[1]
				p_i = site_fracs[p_i_species]
				p_j = site_fracs[p_j_species]
				
				for L_expr, order in binary_params:
					G_ij_binary_expr_rk_sum += L_expr * (p_i - p_j) ** order
				
				# 完整的 G_ij^E = p_i * p_j * G_ij_binary_expr_rk_sum
				
				# 5. 计算所有其他组元 k 的 d_ki 和 d_kj
				X_ij_i_num = x_i
				X_ij_j_num = x_j
				X_denom = x_i + x_j
				
				other_comps = [c for c in active_comps if c != comp_i and c != comp_j]
				
				if not other_comps:
					# 纯二元系，无需外推
					sublattice_total_G_E += p_i * p_j * G_ij_binary_expr_rk_sum
					continue
				
				for comp_k in other_comps:
					x_k = site_fracs[comp_k]
					
					# 6. 计算 d_ki, d_kj (注意组元顺序)
					d_ki = self._get_uem_d_term(dbe, comp_k, comp_i, subl_index)
					d_kj = self._get_uem_d_term(dbe, comp_k, comp_j, subl_index)
					
					# 7. 计算 alpha 系数
					d_sum = d_ki + d_kj
					# 使用 Piecewise 避免 0/0 (如果 d_sum=0, alpha=0)
					d_sum_safe = Piecewise((d_sum, d_sum != 0), (1, True))
					
					alpha_i_k = (d_kj / d_sum_safe) * exp(-d_ki)
					alpha_j_k = (d_ki / d_sum_safe) * exp(-d_kj)
					
					# 8. 累加等效组分
					X_ij_i_num += alpha_i_k * x_k
					X_ij_j_num += alpha_j_k * x_k
					X_denom += (alpha_i_k + alpha_j_k) * x_k
				
				# 9. 计算 X_ij^i 和 X_ij^j
				X_denom_safe = Piecewise((X_denom, X_denom != 0), (1, True))
				X_ij_i = X_ij_i_num / X_denom_safe
				X_ij_j = X_ij_j_num / X_denom_safe
				
				# 10. 代入 G_ij^E 表达式
				# 我们需要替换 G_ij_binary_expr_rk_sum 中的 p_i 和 p_j
				substitution_dict = {}
				# 检查 p_i 是 comp_i 还是 comp_j
				if p_i_species == comp_i:
					substitution_dict[p_i] = X_ij_i
					substitution_dict[p_j] = X_ij_j
				else:
					substitution_dict[p_i] = X_ij_j
					substitution_dict[p_j] = X_ij_i
				
				rk_sum_modified = G_ij_binary_expr_rk_sum.xreplace(substitution_dict)
				
				# 11. 根据UEM公式计算最终贡献
				# G^E = sum [ (x_i * x_j) / (X_i * X_j) ] * [ X_i * X_j * rk_sum(X_i, X_j) ]
				# 这可以简化为:
				# G^E = sum [ (x_i * x_j) * rk_sum(X_i, X_j) ]
				# 这种形式也可以避免 (X_i * X_j) = 0 时的除零问题
				
				final_term = (x_i * x_j) * rk_sum_modified
				sublattice_total_G_E += final_term
			
			total_excess_energy += sublattice_total_G_E
		
		# 12. 返回总能量，并按照基类的要求进行归一化
		#
		return total_excess_energy / self._site_ratio_normalization

# ---------------------------------------------------------------------------
# --- UEM Model Implementation (End) ---
# ---------------------------------------------------------------------------