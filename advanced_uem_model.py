"""
UEM Integration for pycalphad Model Class (Optimized Version)
==============================================================

性能优化版本，主要改进：
1. 添加缓存机制，避免重复计算
2. 批量预处理二元参数
3. 优化符号表达式构建
4. 减少Piecewise嵌套

[优化版本 2025-12-13]
基于原始版本进行性能优化
"""

from symengine import (Add, Mul, Pow, S, exp, log, Symbol,
                       StrictGreaterThan, Piecewise, Abs, Float, Or, And)
from tinydb import where
import itertools
from functools import lru_cache
import pycalphad.variables as v
from pycalphad.model import Model


class ModelWithUEM(Model):
    """
    增强的Model类, 支持UEM外推方法（性能优化版）
    
    优化特性：
    - 二元参数缓存
    - 贡献系数缓存
    - 批量预处理
    - 表达式优化
    
    支持的外推方法 (extrapolation_method):
    - 'muggianu': Muggianu对称模型 (α=0.5, 0.5)
    - 'kohler':   Kohler对称模型 (α=0, 0)
    - 'toop':     Toop-Kohler非对称模型 (需指定asymmetric_component)
    - 'uem1':     UEM1 — 基于无限稀释性质差 D=(2/RT)|Σ L^{2m+1}|
    - 'uem2_n':   UEM2 — 基于混合焓积分 [TODO: 需实现独立D项公式]
    - 'uem_adv':  UEM-Adv高级模型 [TODO: 需实现独立D项公式]
    """
    
    # ========================================================================
    # 类级别缓存（跨实例共享）
    # ========================================================================
    _binary_params_cache = {}
    _contribution_coeff_cache = {}
    
    # 自定义参数名列表，__new__中需要从kwargs过滤掉
    _custom_kwargs = {'extrapolation_method', 'uem_variant', 'delta_function',
                      'use_cache', 'debug', 'asymmetric_component'}
    
    def __new__(cls, *args, **kwargs):
        """
        重写__new__以兼容pycalphad>=0.11的_dispatch_on机制。
        _dispatch_on只接受(dbe, comps, phase_name, parameters)，
        自定义关键字参数必须在传递前过滤掉。
        """
        filtered_kwargs = {k: v for k, v in kwargs.items()
                           if k not in cls._custom_kwargs}
        return super().__new__(cls, *args, **filtered_kwargs)
    
    def __init__(self, dbe, comps, phase_name, parameters=None,
                 extrapolation_method=None, uem_variant='uem1',
                 delta_function=None, use_cache=True, debug=False,
                 asymmetric_component=None):
        
        # 保存 dbe 供后续使用
        self.dbe = dbe
        self.use_cache = use_cache
        self.debug = debug
        
        # 确定使用的外推方法
        phase = dbe.phases[phase_name.upper()]
        if extrapolation_method is None:
            extrapolation_method = phase.model_hints.get('extrapolation_method', 'muggianu')
        
        self.extrapolation_method = extrapolation_method.lower()
        self.uem_variant = uem_variant
        self.delta_function = delta_function
        
        # Toop不对称组元: 默认取字母序第一个非VA组元
        self.asymmetric_component = asymmetric_component
        
        # 实例级别缓存（用于当前计算）
        self._instance_cache = {
            'binary_params': {},
            'inf_dilution': {},
            'alpha': {},
            'd_terms': {}
        }
        
        # 调用父类初始化
        super().__init__(dbe, comps, phase_name, parameters)
    
    def _get_cache_key(self, *args):
        """生成缓存键"""
        return (self.phase_name, self.extrapolation_method) + tuple(str(a) for a in args)
    
    # ========================================================================
    # 优化1: 批量预加载二元参数
    # ========================================================================
    def _preload_binary_params(self, dbe, subl_idx, active_components):
        """
        批量预加载所有二元参数，避免重复查询数据库
        """
        cache_key = self._get_cache_key('preload', subl_idx, tuple(sorted(str(c) for c in active_components)))
        
        if self.use_cache and cache_key in self._instance_cache['binary_params']:
            return self._instance_cache['binary_params'][cache_key]
        
        phase_name_raw = self.phase_name
        phase_name_clean = phase_name_raw.split(':')[0]
        
        # 一次性查询所有二元参数
        def _match_phase_name(p_name):
            return p_name == phase_name_raw or p_name == phase_name_clean
        
        def _is_any_binary_param(const_array):
            if not self._array_validity(const_array):
                return False
            if len(const_array) != len(self.constituents):
                return False
            if len(const_array[subl_idx]) != 2:
                return False
            # 检查是否是活性组元的二元组合
            binary_comps = set(const_array[subl_idx])
            if not binary_comps.issubset(set(active_components)):
                return False
            for k, subl in enumerate(const_array):
                if k != subl_idx and len(subl) != 1:
                    return False
            return True
        
        all_binary_query = (
            (where('phase_name').test(_match_phase_name)) &
            (where('parameter_type').one_of(['G', 'L'])) &
            (where('constituent_array').test(_is_any_binary_param))
        )
        
        all_binary_params = dbe.search(all_binary_query)
        
        # 按组元对组织参数
        params_by_pair = {}
        for param in all_binary_params:
            pair = frozenset(param['constituent_array'][subl_idx])
            if pair not in params_by_pair:
                params_by_pair[pair] = []
            params_by_pair[pair].append(param)
        
        # 处理每个组元对
        result = {}
        for pair, params in params_by_pair.items():
            comp_list = list(pair)
            if len(comp_list) != 2:
                continue
            
            comp1, comp2 = comp_list[0], comp_list[1]
            params.sort(key=lambda p: p['parameter_order'])
            
            if params:
                max_order = params[-1]['parameter_order']
                param_dict = {p['parameter_order']: p['parameter'] for p in params}
                L_params = [param_dict.get(order, S.Zero) for order in range(max_order + 1)]
                first_comp = params[0]['constituent_array'][subl_idx][0]
            else:
                L_params = []
                first_comp = None
            
            result[pair] = (L_params, first_comp)
        
        if self.use_cache:
            self._instance_cache['binary_params'][cache_key] = result
        
        return result
    
    def _get_binary_L_params_cached(self, comp1, comp2, dbe, subl_idx, preloaded_params):
        """
        从预加载的参数中获取二元L参数（优化版）
        """
        pair = frozenset([comp1, comp2])
        if pair in preloaded_params:
            return preloaded_params[pair]
        return ([], None)
    
    # ========================================================================
    # 优化2: 无限稀释能量计算（带缓存）
    # ========================================================================
    def _get_infinite_dilution_energies_cached(self, comp_i, comp_k, L_params_ik, first_comp_in_db):
        """
        计算无限稀释能量（优化版，使用缓存）
        """
        if not L_params_ik:
            return (S.Zero, S.Zero)
        
        # 使用元组作为缓存键（L_params可能包含符号，需要特殊处理）
        cache_key = (str(comp_i), str(comp_k), str(first_comp_in_db), len(L_params_ik))
        
        if self.use_cache and cache_key in self._instance_cache['inf_dilution']:
            return self._instance_cache['inf_dilution'][cache_key]
        
        # 优化：使用列表推导式替代多次Add调用
        sum_L_m = sum(L_params_ik, S.Zero)
        sum_L_m_alternating = sum(
            (param * ((-1) ** order) for order, param in enumerate(L_params_ik)),
            S.Zero
        )
        
        if first_comp_in_db == comp_i:
            g_k_inf_in_i = sum_L_m
            g_i_inf_in_k = sum_L_m_alternating
        elif first_comp_in_db == comp_k:
            g_i_inf_in_k = sum_L_m
            g_k_inf_in_i = sum_L_m_alternating
        else:
            # 默认处理
            g_k_inf_in_i = sum_L_m
            g_i_inf_in_k = sum_L_m_alternating
        
        result = (g_i_inf_in_k, g_k_inf_in_i)
        
        if self.use_cache:
            self._instance_cache['inf_dilution'][cache_key] = result
        
        return result
    
    # ========================================================================
    # 优化3: D项计算（减少重复）
    # ========================================================================
    def _compute_d_terms_batch(self, active_components, dbe, subl_idx, preloaded_params):
        """
        批量计算所有组元对的D项
        """
        cache_key = self._get_cache_key('d_terms', subl_idx, tuple(sorted(str(c) for c in active_components)))
        
        if self.use_cache and cache_key in self._instance_cache['d_terms']:
            return self._instance_cache['d_terms'][cache_key]
        
        R = v.R
        T = v.T
        
        d_terms = {}
        
        for comp_a, comp_b in itertools.permutations(active_components, 2):
            pair_key = (str(comp_a), str(comp_b))
            
            L_params, first_comp = self._get_binary_L_params_cached(
                comp_a, comp_b, dbe, subl_idx, preloaded_params
            )
            
            if not L_params:
                d_terms[pair_key] = S.Zero
                continue
            
            g_a_in_b, g_b_in_a = self._get_infinite_dilution_energies_cached(
                comp_a, comp_b, L_params, first_comp
            )
            
            d_num = g_a_in_b - g_b_in_a
            
            # 简化的安全除法（减少Piecewise复杂度）
            term = d_num / (R * T)
            d_terms[pair_key] = term
        
        if self.use_cache:
            self._instance_cache['d_terms'][cache_key] = d_terms
        
        return d_terms
    
    # ========================================================================
    # 优化4: 贡献系数计算（批量+缓存）
    # ========================================================================
    def _calculate_alpha_batch(self, active_components, dbe, subl_idx, preloaded_params):
        """
        批量计算所有需要的贡献系数
        """
        if self.extrapolation_method == 'muggianu':
            # Muggianu: 所有alpha都是0.5
            alpha_dict = {}
            for k in active_components:
                for i, j in itertools.combinations(active_components, 2):
                    if k == i or k == j:
                        continue
                    alpha_dict[(str(k), str(i), str(j))] = (S.Half, S.Half)
            return alpha_dict
        
        elif self.extrapolation_method == 'kohler':
            # Kohler: 所有alpha都是0
            alpha_dict = {}
            for k in active_components:
                for i, j in itertools.combinations(active_components, 2):
                    if k == i or k == j:
                        continue
                    alpha_dict[(str(k), str(i), str(j))] = (S.Zero, S.Zero)
            return alpha_dict
        
        elif self.extrapolation_method in ('toop', 'toop-kohler'):
            # Toop-Kohler非对称模型
            # 确定不对称组元: 用户指定 > 默认字母序第一个非VA组元
            asym = self.asymmetric_component
            if asym is None:
                # 默认: 字母序第一个
                asym_str = sorted(str(c) for c in active_components)[0]
            else:
                asym_str = str(asym)
            
            alpha_dict = {}
            for k in active_components:
                for i, j in itertools.combinations(active_components, 2):
                    if k == i or k == j:
                        continue
                    
                    key = (str(k), str(i), str(j))
                    i_str, j_str = str(i), str(j)
                    
                    if i_str == asym_str:
                        # 不对称组元是i: α_ki=0, α_kj=1
                        # 效果: X_ij_i = x_i, X_ij_j = 1-x_i (Muggianu对i的投影)
                        alpha_dict[key] = (S.Zero, S.One)
                    elif j_str == asym_str:
                        # 不对称组元是j: α_ki=1, α_kj=0
                        # 效果: X_ij_j = x_j, X_ij_i = 1-x_j
                        alpha_dict[key] = (S.One, S.Zero)
                    else:
                        # 不对称组元不在binary (i,j)中: Kohler投影
                        alpha_dict[key] = (S.Zero, S.Zero)
            
            return alpha_dict
        
        # UEM方法：需要计算
        cache_key = self._get_cache_key('alpha_batch', subl_idx,
                                        tuple(sorted(str(c) for c in active_components)))
        
        if self.use_cache and cache_key in self._instance_cache['alpha']:
            return self._instance_cache['alpha'][cache_key]
        
        # 预计算所有D项
        d_terms = self._compute_d_terms_batch(active_components, dbe, subl_idx, preloaded_params)
        
        alpha_dict = {}
        T = v.T
        
        # 预定义安全常量
        EPSILON = Float(1e-9)
        MAX_EXP = Float(50.0)
        
        for k in active_components:
            for i, j in itertools.combinations(active_components, 2):
                if k == i or k == j:
                    continue
                
                key = (str(k), str(i), str(j))
                
                # 获取D项
                d_ki_signed = d_terms.get((str(i), str(k)), S.Zero)
                d_kj_signed = d_terms.get((str(j), str(k)), S.Zero)

                # 使用Piecewise替代Abs()以避免求导问题
                d_ki = Piecewise(
                    (-d_ki_signed, d_ki_signed < 0),
                    (d_ki_signed, True)
                )
                d_kj = Piecewise(
                    (-d_kj_signed, d_kj_signed < 0),
                    (d_kj_signed, True)
                )

                denominator = d_ki + d_kj

                # 优化：使用简化的Piecewise结构
                # 当denominator接近0时，返回0.5
                is_valid = StrictGreaterThan(denominator, EPSILON)

                # 安全的指数计算
                exp_ki = Piecewise(
                    (S.Zero, StrictGreaterThan(d_ki, MAX_EXP)),
                    (exp(-d_ki), True)
                )
                exp_kj = Piecewise(
                    (S.Zero, StrictGreaterThan(d_kj, MAX_EXP)),
                    (exp(-d_kj), True)
                )
                
                # 计算alpha
                alpha_ki = Piecewise(
                    ((d_kj / denominator) * exp_ki, is_valid),
                    (S.Half, True)
                )
                alpha_kj = Piecewise(
                    ((d_ki / denominator) * exp_kj, is_valid),
                    (S.Half, True)
                )
                
                alpha_dict[key] = (alpha_ki, alpha_kj)
        
        if self.use_cache:
            self._instance_cache['alpha'][cache_key] = alpha_dict
        
        return alpha_dict
    
    # ========================================================================
    # 优化5: 二元过剩能计算
    # ========================================================================
    def _get_binary_GE_rk_optimized(self, i, j, X_ij_i, X_ij_j, L_params_list, first_comp_in_db):
        """
        计算R-K形式的二元过剩Gibbs能（优化版）
        """
        if not L_params_list:
            return S.Zero
        
        # 确定差值项的符号
        if first_comp_in_db == i:
            diff_term = X_ij_i - X_ij_j
        elif first_comp_in_db == j:
            diff_term = X_ij_j - X_ij_i
        else:
            diff_term = X_ij_i - X_ij_j
        
        # 优化：使用Horner法则计算多项式（减少幂运算）
        # L0 + L1*d + L2*d^2 + ... = L0 + d*(L1 + d*(L2 + ...))
        if len(L_params_list) == 1:
            rk_sum = L_params_list[0]
        else:
            # Horner法则
            rk_sum = L_params_list[-1]
            for k in range(len(L_params_list) - 2, -1, -1):
                rk_sum = L_params_list[k] + diff_term * rk_sum
        
        return X_ij_i * X_ij_j * rk_sum
    
    # ========================================================================
    # 优化6: 主计算方法
    # ========================================================================
    def excess_mixing_energy(self, dbe=None):
        """
        计算过剩混合能（性能优化版）
        """
        if dbe is None:
            dbe = self.dbe
        
        total_GE = S.Zero
        
        # 预定义常量
        EPSILON = Float(1e-9)
        
        for subl_idx, subl_constituents in enumerate(self.constituents):
            active_components = sorted(
                list(subl_constituents.intersection(self.components)),
                key=str
            )
            
            if len(active_components) < 2:
                continue
            
            # 优化：预加载所有二元参数
            preloaded_params = self._preload_binary_params(dbe, subl_idx, active_components)
            
            # 优化：批量计算所有贡献系数
            alpha_dict = self._calculate_alpha_batch(
                active_components, dbe, subl_idx, preloaded_params
            )
            
            # 构建位置分数字典
            x = {
                comp: v.SiteFraction(self.phase_name, subl_idx, comp)
                for comp in active_components
            }
            
            sublattice_GE = S.Zero
            
            for i, j in itertools.combinations(active_components, 2):
                # 计算有效摩尔分数
                numerator_i = x[i]
                denominator = x[i] + x[j]
                
                # 累加第三组元的贡献
                for k in active_components:
                    if k == i or k == j:
                        continue
                    
                    key = (str(k), str(i), str(j))
                    alpha_ki, alpha_kj = alpha_dict.get(key, (S.Half, S.Half))
                    
                    numerator_i = numerator_i + alpha_ki * x[k]
                    denominator = denominator + (alpha_ki + alpha_kj) * x[k]
                
                # 计算有效摩尔分数（带安全保护）
                # 使用Or条件替代Abs()以避免求导问题: |x| > eps <=> x > eps OR x < -eps
                is_denom_valid = Or(denominator > EPSILON, denominator < -EPSILON)
                X_ij_i = Piecewise(
                    (numerator_i / denominator, is_denom_valid),
                    (S.Half, True)
                )
                X_ij_j = S.One - X_ij_i

                # 获取二元参数
                pair = frozenset([i, j])
                L_params, first_comp = preloaded_params.get(pair, ([], None))

                # 计算二元过剩能
                G_ij_E = self._get_binary_GE_rk_optimized(
                    i, j, X_ij_i, X_ij_j, L_params, first_comp
                )

                if G_ij_E == S.Zero:
                    continue

                # 缩放因子
                scaling_num = x[i] * x[j]
                scaling_den = X_ij_i * X_ij_j

                is_scaling_valid = Or(scaling_den > EPSILON, scaling_den < -EPSILON)
                G_E_term = Piecewise(
                    ((scaling_num / scaling_den) * G_ij_E, is_scaling_valid),
                    (S.Zero, True)
                )
                
                sublattice_GE = sublattice_GE + G_E_term
            
            total_GE = total_GE + sublattice_GE
        
        return total_GE / self._site_ratio_normalization
    
    # ========================================================================
    # 缓存管理方法
    # ========================================================================
    def clear_instance_cache(self):
        """清除实例级别缓存"""
        self._instance_cache = {
            'binary_params': {},
            'inf_dilution': {},
            'alpha': {},
            'd_terms': {}
        }
    
    @classmethod
    def clear_class_cache(cls):
        """清除类级别缓存"""
        cls._binary_params_cache.clear()
        cls._contribution_coeff_cache.clear()
    
    def get_cache_stats(self):
        """获取缓存统计信息"""
        return {
            'binary_params': len(self._instance_cache['binary_params']),
            'inf_dilution': len(self._instance_cache['inf_dilution']),
            'alpha': len(self._instance_cache['alpha']),
            'd_terms': len(self._instance_cache['d_terms'])
        }

    def get_contribution_coefficients(self, active_components, subl_idx=0, temperature=1000.0, logger=None):
        """
        获取贡献系数的数值（用于调试和输出）

        Parameters
        ----------
        active_components : list
            活性组分列表
        subl_idx : int
            亚晶格索引（默认0）
        temperature : float
            温度（K），用于计算数值
        logger : callable
            日志输出函数

        Returns
        -------
        dict : {(k, i, j): (alpha_ki_value, alpha_kj_value)}
        """
        import sympy as sp
        from pycalphad import variables as v

        # 使用logger或print
        log = logger if logger else print

        # 预加载参数
        preloaded_params = self._preload_binary_params(self.dbe, subl_idx, active_components)

        # 计算alpha_dict（符号表达式）
        alpha_dict = self._calculate_alpha_batch(active_components, self.dbe, subl_idx, preloaded_params)

        # 转换为数值
        result = {}
        for key, (alpha_ki, alpha_kj) in alpha_dict.items():
            alpha_ki_expr = None
            alpha_kj_expr = None
            try:
                # 将符号表达式代入温度值
                # 确保temperature是Python float而不是numpy类型
                T_val = float(temperature)

                # 使用字典方式进行替换，可能更稳定
                from pycalphad import variables as v
                alpha_ki_expr = alpha_ki.subs({v.T: T_val})
                alpha_kj_expr = alpha_kj.subs({v.T: T_val})

                # 转换为数值
                import numpy as np

                if hasattr(alpha_ki_expr, 'evalf'):
                    # 符号表达式，使用evalf
                    alpha_ki_val = float(alpha_ki_expr.evalf())
                    alpha_kj_val = float(alpha_kj_expr.evalf())
                elif isinstance(alpha_ki_expr, (int, float, np.number)):
                    # 已经是数值
                    alpha_ki_val = float(alpha_ki_expr)
                    alpha_kj_val = float(alpha_kj_expr)
                else:
                    # 尝试直接转换
                    alpha_ki_val = float(alpha_ki_expr)
                    alpha_kj_val = float(alpha_kj_expr)

                result[key] = (alpha_ki_val, alpha_kj_val)

            except Exception as e:
                # 如果转换失败，记录详细错误并使用默认值
                import traceback
                log(f"  [警告] 计算失败 {key}:")
                log(f"    错误类型: {type(e).__name__}")
                log(f"    错误信息: {e}")
                if alpha_ki_expr is not None:
                    log(f"    alpha_ki_expr类型: {type(alpha_ki_expr)}")
                else:
                    log(f"    subs()调用失败，alpha_ki_expr未定义")
                # 打印堆栈跟踪的最后几行
                tb_lines = traceback.format_exc().split('\n')
                for line in tb_lines[-4:-1]:
                    log(f"      {line}")
                result[key] = (0.5, 0.5)

        return result


# ============================================================================
# 便捷类（不同UEM变体）
# ============================================================================

class ModelUEM1(ModelWithUEM):
    """UEM1变体 - 基于无限稀释性质（推荐）"""
    def __init__(self, dbe, comps, phase_name, parameters=None):
        super().__init__(dbe, comps, phase_name, parameters,
                         extrapolation_method='uem1', uem_variant='uem1')


class ModelUEMAdv(ModelWithUEM):
    """UEM-Adv变体 - 高级几何方法
    
    WARNING: 当前实现与UEM1相同，尚未实现独立的D项公式。
    TODO: 需要实现基于几何平均的D项计算方法。
    """
    def __init__(self, dbe, comps, phase_name, parameters=None):
        import warnings
        warnings.warn(
            "ModelUEMAdv当前使用UEM1的D项公式，尚未实现独立算法。"
            "计算结果与ModelUEM1完全相同。",
            UserWarning, stacklevel=2
        )
        super().__init__(dbe, comps, phase_name, parameters,
                         extrapolation_method='uem_adv', uem_variant='uem_adv')


class ModelUEM2N(ModelWithUEM):
    """UEM2-N变体 - 基于混合焓积分
    
    WARNING: 当前实现与UEM1相同，尚未实现独立的D项公式。
    TODO: 需要实现基于积分型性质差的D项计算:
          D_{k-i}^{UEM2} 应使用二元过剩性质的积分平均，
          而非UEM1的无限稀释性质差。
    """
    def __init__(self, dbe, comps, phase_name, parameters=None):
        import warnings
        warnings.warn(
            "ModelUEM2N当前使用UEM1的D项公式，尚未实现独立算法。"
            "计算结果与ModelUEM1完全相同。",
            UserWarning, stacklevel=2
        )
        super().__init__(dbe, comps, phase_name, parameters,
                         extrapolation_method='uem2_n', uem_variant='uem2_n')


class ModelMuggianu(ModelWithUEM):
    """Muggianu模型（传统对称模型）"""
    def __init__(self, dbe, comps, phase_name, parameters=None):
        super().__init__(dbe, comps, phase_name, parameters,
                         extrapolation_method='muggianu', uem_variant='muggianu')


class ModelKohler(ModelWithUEM):
    """Kohler模型（传统对称模型）"""
    def __init__(self, dbe, comps, phase_name, parameters=None):
        super().__init__(dbe, comps, phase_name, parameters,
                         extrapolation_method='kohler', uem_variant='kohler')


class ModelToop(ModelWithUEM):
    """Toop-Kohler模型（传统非对称模型）
    
    Parameters
    ----------
    asymmetric_component : str, optional
        不对称组元名称 (如 'AL')。未指定时默认取字母序第一个非VA组元。
        Toop模型中，包含不对称组元的二元子系统用Muggianu投影，
        不包含不对称组元的二元子系统用Kohler投影。
    """
    def __init__(self, dbe, comps, phase_name, parameters=None,
                 asymmetric_component=None):
        super().__init__(dbe, comps, phase_name, parameters,
                         extrapolation_method='toop', uem_variant='toop-kohler',
                         asymmetric_component=asymmetric_component)


# 向后兼容别名
ModelUEM = ModelUEM1
ModelUEMRedlichKister = ModelWithUEM


# ============================================================================
# 性能测试工具
# ============================================================================

def benchmark_model(dbe, comps, phase_name, n_iterations=10):
    """
    性能基准测试
    
    使用示例:
        from pycalphad import Database
        dbe = Database('your_database.tdb')
        results = benchmark_model(dbe, ['AL', 'FE', 'NI', 'VA'], 'LIQUID')
    """
    import time
    
    results = {}
    
    # 测试原始Model
    try:
        start = time.perf_counter()
        for _ in range(n_iterations):
            model_orig = Model(dbe, comps, phase_name)
        end = time.perf_counter()
        results['Original Model'] = (end - start) / n_iterations
    except Exception as e:
        results['Original Model'] = f"Error: {e}"
    
    # 测试优化后的UEM Model
    try:
        start = time.perf_counter()
        for _ in range(n_iterations):
            model_uem = ModelWithUEM(dbe, comps, phase_name, use_cache=True)
        end = time.perf_counter()
        results['UEM Model (cached)'] = (end - start) / n_iterations
    except Exception as e:
        results['UEM Model (cached)'] = f"Error: {e}"
    
    # 测试无缓存版本
    try:
        start = time.perf_counter()
        for _ in range(n_iterations):
            model_uem_nc = ModelWithUEM(dbe, comps, phase_name, use_cache=False)
        end = time.perf_counter()
        results['UEM Model (no cache)'] = (end - start) / n_iterations
    except Exception as e:
        results['UEM Model (no cache)'] = f"Error: {e}"
    
    return results


# ============================================================================
# 全局启用UEM
# ============================================================================

def enable_uem_globally():
    """
    全局启用UEM支持（替换pycalphad.model.Model）
    """
    import pycalphad.model
    pycalphad.model.Model = ModelWithUEM
    print("UEM support enabled globally (optimized version)")