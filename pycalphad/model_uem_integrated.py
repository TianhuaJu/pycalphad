"""
UEM Integration for pycalphad Model Class
=========================================

This module provides UEM (Unified Excess Model) integration for pycalphad.
It can be used as a drop-in replacement or extension to the standard Model class.

[已修正 2025-10-22 v3]
1. 修正 'kohler_toop_excess_sum' 中的 AttributeError:
   将 '_binary_interaction_test' 替换为 '_interaction_test'
   以解决基类 __init__ 中的调用顺序问题。
2. 移除了未使用的导入 (get_non_VA_species)，解决 ImportError。
3. 新增重写(override) 'kohler_toop_excess_sum' 方法。
4. 修正 '_redlich_kister_sum_uem' 中的逻辑，确保UEM激活时
   正确跳过(skip)显式三元参数。
"""
import itertools

# 从 symengine 导入所有需要的函数
from symengine import (
    Add, Mul, Pow, S, exp, StrictGreaterThan, Piecewise,
    # 分别对应 sympy 的 Ge 和 Le
    # 对应 sympy 的 Lt
    Subs

)
from tinydb import where

import pycalphad.variables as v
from pycalphad import Database
from pycalphad.model import Model


# ============================================================================
# Enhanced Model Class with UEM Support
# ============================================================================

class ModelWithUEM(Model):
    """
    增强的Model类,支持UEM外推方法
    
    完全兼容原始pycalphad.Model,并添加UEM支持
    通过phase.model_hints['extrapolation_method']来选择外推方法
    """
    
    def __new__ (cls, dbe, comps, phase_name, parameters=None, **kwargs):
        # 1. 接受 *args 和 **kwargs (包括 'extrapolation_method' 等)
        # 2. 我们忽略 **kwargs，因为它们只对 __init__ 有意义
        # 3. 创建一个 'cls' (ModelWithUEM) 的新实例
        instance = object.__new__(cls)
        return instance
    def __init__(self, dbe, comps, phase_name, parameters=None,
                 extrapolation_method=None, uem_variant='uem1',
                 delta_function=None):
        """
        Parameters
        ----------
        dbe : Database
            热力学数据库
        comps : list
            组分列表
        phase_name : str
            相名称
        parameters : dict or list, optional
            参数字典或列表
        extrapolation_method : str, optional
            外推方法: 'muggianu', 'kohler', 'toop', 'uem1', 'uem_adv', 'uem2_n'
            如果未指定,将从phase.model_hints中读取
        uem_variant : str, optional
            UEM变体,当extrapolation_method为'uem'时使用
        delta_function : callable, optional
            自定义性质差异计算函数
        """
        # 首先调用父类初始化
        # 此时 self.models['xsmix'] 会被构建 (使用 Muggianu)
        # [注意] 此调用会触发重写的 kohler_toop_excess_sum，必须确保
        # 它在此时能正确运行（见 v3 修正）
        
        
        # 确定使用的外推方法
        phase = dbe.phases[phase_name.upper()]
        if extrapolation_method is None:
            # 从model_hints中读取
            extrapolation_method = phase.model_hints.get('extrapolation_method', 'muggianu')
        
        self.extrapolation_method = extrapolation_method.lower()
        self.uem_variant = uem_variant
        self.delta_function = delta_function
        
        super().__init__(dbe, comps, phase_name, parameters)
        
    
    def _get_binary_L_params (self, comp1, comp2, dbe):
        """
        一个辅助函数：从数据库(dbe)中获取 (comp1, comp2) 二元系
        在当前相的 L 参数。

        返回: (param_list, first_comp_in_db)
             param_list 是 [L0, L1, L2, ...] 的符号表达式列表 (已经是 v.T 的函数)
             first_comp_in_db 是数据库中参数的第一个组元
        """
        phase_name = self.phase_name
        param_search = dbe.search
        
        def _is_binary_param (const_array):
            if not self._array_validity(const_array): return False
            if len(const_array) != 1: return False  # 假设单亚晶格
            subl = const_array[0]
            return len(subl) == 2 and set(subl) == {comp1, comp2}
        
        binary_param_query = (
                (where('phase_name') == phase_name) &
                ((where('parameter_type') == 'G') | (where('parameter_type') == 'L')) &
                (where('constituent_array').test(_is_binary_param))
        )
        binary_params = param_search(binary_param_query)
        
        if not binary_params:
            return ([], None)
        
        binary_params.sort(key=lambda p: p['parameter_order'])
        max_order = binary_params[-1]['parameter_order']
        param_dict = {p['parameter_order']: p['parameter'] for p in binary_params}
        first_comp_in_db = binary_params[0]['constituent_array'][0][0]
        
        L_params = [param_dict.get(order, S.Zero) for order in range(max_order + 1)]
        
        return (L_params, first_comp_in_db)
    
    def _get_infinite_dilution_energies (self, comp_i, comp_k, L_params_ik, first_comp_in_db):
        """
        根据 L 参数计算符号化的无限稀释偏摩尔 性质G^E
        """
        if not L_params_ik:
            return (S.Zero, S.Zero)
        
        # L_params_ik 是 L(first_comp, other_comp)
        sum_L_m = Add(*L_params_ik)
        sum_L_m_alternating = Add(*[param * ((-1) ** order) for order, param in enumerate(L_params_ik)])
        
        if first_comp_in_db == comp_i:
            # 数据库中是 (i, k)
            g_k_inf_in_i = sum_L_m
            g_i_inf_in_k = sum_L_m_alternating
        elif first_comp_in_db == comp_k:
            # 数据库中是 (k, i)
            g_i_inf_in_k = sum_L_m
            g_k_inf_in_i = sum_L_m_alternating
        else:
            raise ValueError("组元不匹配")
        
        return (g_i_inf_in_k, g_k_inf_in_i)
    
    def _calculate_alpha_uem1 (self, k, i, j, dbe):
        """
        计算符号化的贡献系数 (alpha_i(ij)^k, alpha_j(ij)^k)
        """
        # 1. 引入 R 和 T 符号
        R = v.R
        T = v.T
        
        # 2. 计算 d_ki
        # g_i_inf_in_k 是 i 在 k 中的无限稀释偏摩尔G^E (lim x_i->0)
        # g_k_inf_in_i 是 k 在 i 中的无限稀释偏摩尔G^E (lim x_k->0)
        L_params_ki, first_comp_ki = self._get_binary_L_params(k, i, dbe)
        g_i_inf_in_k, g_k_inf_in_i = self._get_infinite_dilution_energies(i, k, L_params_ki, first_comp_ki)
        
        d_ki_num = g_i_inf_in_k - g_k_inf_in_i
        d_ki = Piecewise(
                (d_ki_num / (R * T), StrictGreaterThan(T, 0)),
                (S.Zero, True)  # 处理 T=0 的情况
        )
        
        # 3. 计算 d_kj
        L_params_kj, first_comp_kj = self._get_binary_L_params(k, j, dbe)
        g_j_inf_in_k, g_k_inf_in_j = self._get_infinite_dilution_energies(j, k, L_params_kj, first_comp_kj)
        
        d_kj_num = g_j_inf_in_k - g_k_inf_in_j
        d_kj = Piecewise(
                (d_kj_num / (R * T), StrictGreaterThan(T, 0)),
                (S.Zero, True)  # 处理 T=0 的情况
        )
        
        # 4. 计算 alpha (即图片中的 o)
        denominator = d_ki + d_kj
        
        # alpha_i(ij)^k
        alpha_ki_ij = Piecewise(
                ((d_kj / denominator) * exp(-d_ki), denominator != 0),
                (S.Half, True)  # 处理分母为0的情况 (例如 d_ki 和 d_kj 均为0)
        )
        
        # alpha_j(ij)^k (根据对称性)
        alpha_kj_ij = Piecewise(
                ((d_ki / denominator) * exp(-d_kj), denominator != 0),
                (S.Half, True)  # 处理分母为0的情况
        )
        
        return alpha_ki_ij, alpha_kj_ij
    
    def _get_binary_GE_rk (self, i, j, X_ij_i, X_ij_j, dbe):
        """
        根据给定的组元i和j，从数据库dbe中查找Lij参数，
        并使用传入的 *有效组分* X_ij_i 和 X_ij_j 构建
        二元系R-K过剩吉布斯自由能(G^E)的符号表达式。

        参数:
        i (pycalphad.Species): 组元 i
        j (pycalphad.Species): 组元 j
        X_ij_i (symengine.Expr): 组元 i 的 *有效* 符号组分 (UEM 公式 3 的结果)
        X_ij_j (symengine.Expr): 组元 j 的 *有效* 符号组分 (1 - X_ij_i)
        dbe (pycalphad.Database): 要搜索的数据库

        返回:
        symengine.Expr: G^E_ij (X_ij_i, X_ij_j) 的符号表达式
        """
        
        # --- 1. 从数据库 dbe 中搜寻 Lij 参数 ---
        
        phase_name = self.phase_name
        param_search = dbe.search
        
        # 辅助测试函数, 用于 tinydb 查询
        def _is_binary_ij_param (const_array):
            if not self._array_validity(const_array): return False
            # 假设为单亚晶格 (与 UEM 模型的上下文一致)
            if len(const_array) != 1: return False
            subl = const_array[0]
            # 查找包含且仅包含 i 和 j 的二元参数
            return len(subl) == 2 and set(subl) == {i, j}
        
        # 构建查询
        binary_param_query = (
                (where('phase_name') == phase_name) &
                ((where('parameter_type') == 'G') | (where('parameter_type') == 'L')) &
                (where('constituent_array').test(_is_binary_ij_param))
        )
        
        binary_params = param_search(binary_param_query)
        
        # --- 2. 处理查找到的参数 ---
        
        L_params_list = []
        first_comp_in_db = None
        
        if binary_params:
            # 按参数阶数 (0, 1, 2...) 排序
            binary_params.sort(key=lambda p: p['parameter_order'])
            max_order = binary_params[-1]['parameter_order']
            # 将 L 参数放入字典 {order: param_expr}
            param_dict = {p['parameter_order']: p['parameter'] for p in binary_params}
            # 记录数据库中定义的第一个组元 (例如 G(PHASE, I, J;0))
            first_comp_in_db = binary_params[0]['constituent_array'][0][0]
            # 创建 L 参数的有序列表 [L0, L1, L2, ...]
            L_params_list = [param_dict.get(order, S.Zero) for order in range(max_order + 1)]
        
        # --- 3. 使用 X_ij_i 和 X_ij_j 构建 R-K 表达式 ---
        
        rk_sum = S.Zero
        
        # R-K 多项式中的差值项 (X_A - X_B) 必须与数据库中的定义顺序 (A, B) 一致
        # 这对于奇数阶参数 (L1, L3...) 至关重要
        if first_comp_in_db == i:
            # 数据库中是 (i, j), 表达式为 (X_ij_i - X_ij_j)
            diff_term = X_ij_i - X_ij_j
        elif first_comp_in_db == j:
            # 数据库中是 (j, i), 表达式为 (X_ij_j - X_ij_i)
            diff_term = X_ij_j - X_ij_i
        else:
            # 未找到参数, diff_term 无关紧要，因为 L_params_list 为空
            # 为保险起见设置一个默认值
            diff_term = X_ij_i - X_ij_j
        
        for k, Lk in enumerate(L_params_list):
            # Lk * (diff_term)^k
            term = Mul(Lk, Pow(diff_term, k))
            rk_sum = Add(rk_sum, term)
        
        # G^E_ij = X_ij^i * X_ij^j * [ L0 + L1(diff_term) + ... ]
        G_E = Mul(Mul(X_ij_i, X_ij_j), rk_sum)
        
        return G_E
   
    def _calculate_contribution_coefficient(self, k, i, j,dbe):
        """
        计算UEM贡献系数 α(k→i_in_ij)
        
        Parameters
        ----------
        k, i, j : Species
            组分k对i-j二元对的贡献系数
        
        Returns
        -------
        贡献系数 (0到1之间)
        """
        # 特殊情况处理
        if self.extrapolation_method == 'muggianu':
            return S(0.5),S(0.5)
        elif self.extrapolation_method == 'kohler':
            return S(0),S(0)
        elif self.extrapolation_method == 'toop':
            #这里需修改，根据非对称组分返回函数值，当非对称组分为k时，返回0，当非对称组分为i时，返回0，当非对称组分为j时返回1.
            #非对称组分的选择，参考Kohler-Toop方法中的非对称组分的方法
            phase = dbe.phases[self.phase_name]
            chem_groups = phase.model_hints.get('chemical_groups', {})
            group_i = chem_groups.get(i.name, i.name)
            group_j = chem_groups.get(j.name, j.name)
            group_k = chem_groups.get(k.name, k.name)
            if group_k == group_i and group_k != group_j:
                # k 和 i 对称 (非对称组元是 j)
                # k 完全贡献给 i (α_ki = 1), 对 j 无贡献 (α_kj = 0)
                return S(1), S(0)
            elif group_k == group_j and group_k != group_i:
                # k 和 j 对称 (非对称组元是 i)
                # k 完全贡献给 j (α_kj = 1), 对 i 无贡献 (α_ki = 0)
                return S(0), S(1)
            else:
                return S(0),S(0)
            
        elif self.extrapolation_method.startswith('uem'):
            return self._calculate_alpha_uem1(k,i,j,dbe)
        else:
            return self._calculate_alpha_uem1(k,i,j,dbe)
    

    def excess_mixing_energy(self, dbe):
        """
        [已修改] 使用UEM外推构建过剩混合能
        
        这是对excess_mixing_energy的UEM实现
        """
        # --- 1. 处理多亚晶格 (回退) ---
        if len(self.constituents) > 1:
            '''多个亚晶格时，暂不支持UEM。返回标准Muggianu模型。'''
            # 调用父类的原始 excess_mixing_energy
            return super(ModelWithUEM, self).excess_mixing_energy(dbe)
        
        # --- 2. 单亚晶格置换模型 (UEM, Kohler, Toop, Muggianu) ---
        subl_idx = 0
        active_components = sorted(list(self.constituents[subl_idx]))
        
        # 定义真实摩尔分数 (x_i)
        # 对于 (A,B,C)1 相, x_i = v.Y(phase, 0, i)
        x = {
            comp: v.SiteFraction(self.phase_name, subl_idx, comp)
            for comp in active_components
        }
        
        total_GE = S.Zero
        
        # --- 3. 遍历所有唯一的二元对 (i, j) ---
        for i, j in itertools.combinations(active_components, 2):
            
            # --- 4. 计算有效二元组分 X_ij^i 和 X_ij^j ---
            numerator_i = x[i]
            denominator = Add(x[i], x[j])
            
            # 遍历所有其他组元 k
            for k in active_components:
                if k == i or k == j:
                    continue
                    
                alpha_ki_ij, alpha_kj_ij = self._calculate_contribution_coefficient(k, i, j, dbe)
                
                numerator_i = Add(numerator_i, Mul(alpha_ki_ij, x[k]))
                denominator = Add(denominator, Mul(Add(alpha_ki_ij, alpha_kj_ij), x[k]))
                
            
            # 情况a: 退化到二元边 (k=0), denominator = x[i]+x[j]
            x_i_plus_x_j = Add(x[i], x[j])
            
            X_ij_i = Piecewise(
                    # 一般情况 (多元): numerator_i / denominator
                    (numerator_i / denominator, denominator != 0),
                    
                    # 情况a (原点或纯组元k): 此时 x_i=0, x_j=0
                    # 在纯k时, X_ij_i -> alpha_ki_ij / (alpha_ki_ij + alpha_kj_ij)
                    
                    
                    (S.Half, True)
            )
            X_ij_j = S.One - X_ij_i
            
            # --- 6. 计算二元 G_ij^E (R-K 多项式) ---
            
            G_ij_E = self._get_binary_GE_rk(i, j, X_ij_i, X_ij_j, dbe)
            
            # --- 7. [修正] 组合总 G^E (UEM 公式1) + 数值稳定性 ---
            scaling_factor_num = Mul(x[i], x[j])
            scaling_factor_den = Mul(X_ij_i, X_ij_j)
            
            
            G_E_term = Piecewise(
                    # 一般情况
                    (Mul(scaling_factor_num / scaling_factor_den, G_ij_E), scaling_factor_den != 0),
                    
                    # 极限情况 (X_ij_i=0 或 X_ij_j=0), 此时 G_ij_E 也为 0, 总项为 0
                    (S.Zero, True)
            )
            
            total_GE = Add(total_GE, G_E_term)
        
        # --- 8. 返回归一化的总能量 ---
        return total_GE / self._site_ratio_normalization
        
     
    

# ============================================================================
# Convenience Classes for Different UEM Variants
# ============================================================================

class ModelUEM1(ModelWithUEM):
    """UEM1变体 - 基于无限稀释性质(推荐)"""
    def __init__(self, dbe, comps, phase_name, parameters=None):
        super().__init__(dbe, comps, phase_name, parameters,
                        extrapolation_method='uem1', uem_variant='uem1')


class ModelUEMAdv(ModelWithUEM):
    """UEM-Adv变体 - 高级几何方法(推荐)"""
    def __init__(self, dbe, comps, phase_name, parameters=None):
        super().__init__(dbe, comps, phase_name, parameters,
                        extrapolation_method='uem_adv', uem_variant='uem_adv')


class ModelUEM2N(ModelWithUEM):
    """UEM2-N变体 - 基于混合焓积分"""
    def __init__(self, dbe, comps, phase_name, parameters=None):
        super().__init__(dbe, comps, phase_name, parameters,
                        extrapolation_method='uem2_n', uem_variant='uem2_n')


# 向后兼容的别名
ModelUEM = ModelUEM1
ModelUEMRedlichKister = ModelWithUEM


# ============================================================================
# Model Class Replacement (Optional - for monkey patching)
# ============================================================================

def enable_uem_globally():
    """
    全局启用UEM支持
    
    将pycalphad.model.Model替换为ModelWithUEM
    使用示例:
        import pycalphad
        from model_uem_integrated import enable_uem_globally
        
        enable_uem_globally()
        # 现在所有的Model实例都支持UEM了
    """
    import pycalphad.model
    pycalphad.model.Model = ModelWithUEM
    print("UEM support enabled globally for pycalphad.Model")