"""
多元系热力学性质计算模块

该模块提供基于R-K多项式和UEM外推的多元系统热力学性质计算。

主要功能:
---------
- 混合焓计算
- 过剩Gibbs能计算
- 活度系数计算
- 活度计算

计算方法:
---------
1. 对每个二元对(i,j)，计算其在多元系中的有效摩尔分数
2. 使用外推模型(UEM等)计算第三组分的贡献系数
3. 基于有效摩尔分数和贡献系数，外推得到多元系性质
"""

import math
from typing import Callable, Dict, Optional
from itertools import combinations
from sympy import symbols, diff, N, Symbol

from pycalphad.models.rk_uem.rk_binary import RKBinaryPolynomial


# 气体常数
R_GAS_CONSTANT = 8.314  # J/(mol·K)


class ThermodynamicCalculator:
    """
    多元系热力学性质计算器

    基于R-K多项式二元数据和外推模型计算多元系统的热力学性质。

    Parameters
    ----------
    database_path : str
        R-K参数数据库路径

    Attributes
    ----------
    database_path : str
        数据库路径
    warnings : list
        计算过程中的警告信息列表
    """

    def __init__(self, database_path: str):
        """
        初始化热力学计算器

        Parameters
        ----------
        database_path : str
            R-K参数数据库路径
        """
        self.database_path = database_path
        self.warnings = []

    @staticmethod
    def _safe_float(expr) -> float:
        """
        安全地将表达式转换为浮点数

        Parameters
        ----------
        expr : various
            需要转换的表达式

        Returns
        -------
        float
            转换后的浮点值
        """
        try:
            return float(expr)
        except (TypeError, ValueError):
            try:
                return float(N(expr))
            except (TypeError, ValueError, ZeroDivisionError):
                return 0.0

    def _validate_and_normalize_composition(self, composition: Dict[str, float]) -> Dict[str, float]:
        """
        验证并归一化组成

        Parameters
        ----------
        composition : dict
            组成字典 {组分: 摩尔分数}

        Returns
        -------
        dict
            归一化的组成字典
        """
        if not composition:
            raise ValueError("组成字典不能为空")

        # 检查负值
        for comp, frac in composition.items():
            if frac < 0:
                raise ValueError(f"组分 {comp} 的摩尔分数不能为负值: {frac}")

        # 归一化
        total = sum(composition.values())
        if total <= 1e-10:
            raise ValueError("总摩尔分数必须大于0")

        if abs(total - 1.0) > 1e-6:
            return {comp: frac / total for comp, frac in composition.items()}
        else:
            return composition.copy()

    def get_specified_binary_composition(
        self,
        composition: Dict[str, float],
        binary_pair: tuple,
        temperature: float,
        extrapolation_model: Callable,
    ) -> Optional[tuple]:
        """
        获取指定二元对在多元系中的有效摩尔分数

        通过外推模型计算第三组分对二元对的贡献，得到有效摩尔分数。

        Parameters
        ----------
        composition : dict
            多元系组成
        binary_pair : tuple
            二元对 (comp_i, comp_j)
        temperature : float
            温度 (K)
        extrapolation_model : callable
            外推模型函数

        Returns
        -------
        tuple or None
            (RKBinaryPolynomial对象, X_i_effective, X_j_effective)
        """
        comp_i, comp_j = binary_pair
        composition = self._validate_and_normalize_composition(composition)

        # 检查组分是否存在
        if comp_i not in composition or comp_j not in composition:
            raise ValueError(f"组分 {comp_i} 或 {comp_j} 不在组成中")

        if comp_i == comp_j:
            raise ValueError(f"不能指定相同的组分: {comp_i}")

        x_i = composition[comp_i]
        x_j = composition[comp_j]

        # 创建二元系对象
        binary_sys = RKBinaryPolynomial((comp_i, comp_j), self.database_path)

        # 如果就是二元系，直接返回
        if len(composition) == 2:
            return binary_sys, x_i, x_j

        # 计算第三组分的贡献
        sum_contrib_to_i = 0.0
        sum_contrib_to_j = 0.0

        for comp_k, x_k in composition.items():
            if comp_k not in (comp_i, comp_j):
                try:
                    # 计算贡献系数
                    r_ki = extrapolation_model(comp_k, comp_i, comp_j, temperature, self.database_path)
                    r_kj = extrapolation_model(comp_k, comp_j, comp_i, temperature, self.database_path)

                    sum_contrib_to_i += r_ki * x_k
                    sum_contrib_to_j += r_kj * x_k

                except Exception as e:
                    self.warnings.append(f"计算组分 {comp_k} 的贡献系数时出错: {e}")
                    continue

        # 计算有效摩尔分数
        x_eff_i = x_i + sum_contrib_to_i
        x_eff_j = x_j + sum_contrib_to_j
        total_eff = x_eff_i + x_eff_j

        if total_eff == 0:
            # 符号计算情况
            is_symbolic = any(isinstance(v, Symbol) for v in composition.values())
            if is_symbolic:
                return binary_sys, 0, 0
            raise ValueError(f"二元对 {comp_i}-{comp_j} 的有效摩尔分数为零")

        # 归一化有效摩尔分数
        X_i_eff = x_eff_i / total_eff
        X_j_eff = x_eff_j / total_eff

        return binary_sys, X_i_eff, X_j_eff

    def get_mixing_enthalpy(
        self,
        composition: Dict[str, float],
        temperature: float,
        extrapolation_model: Callable,
    ) -> float:
        """
        计算多元系混合焓

        Parameters
        ----------
        composition : dict
            组成字典 {组分: 摩尔分数}
        temperature : float
            温度 (K)
        extrapolation_model : callable
            外推模型函数

        Returns
        -------
        float
            混合焓 (J/mol)
        """
        return self._calculate_excess_property(
            'mixing_enthalpy', composition, temperature, extrapolation_model
        )

    def get_excess_gibbs(
        self,
        composition: Dict[str, float],
        temperature: float,
        extrapolation_model: Callable,
    ) -> float:
        """
        计算多元系过剩Gibbs能

        Parameters
        ----------
        composition : dict
            组成字典 {组分: 摩尔分数}
        temperature : float
            温度 (K)
        extrapolation_model : callable
            外推模型函数

        Returns
        -------
        float
            过剩Gibbs能 (J/mol)
        """
        return self._calculate_excess_property(
            'excess_gibbs_energy', composition, temperature, extrapolation_model
        )

    def _calculate_excess_property(
        self,
        property_name: str,
        composition: Dict[str, float],
        temperature: float,
        extrapolation_model: Callable,
    ) -> float:
        """
        计算过剩性质的通用函数

        Parameters
        ----------
        property_name : str
            性质名称 ('mixing_enthalpy' 或 'excess_gibbs_energy')
        composition : dict
            组成字典
        temperature : float
            温度 (K)
        extrapolation_model : callable
            外推模型函数

        Returns
        -------
        float
            过剩性质值
        """
        self.warnings.clear()
        composition = self._validate_and_normalize_composition(composition)

        total_property = 0.0
        components = list(composition.keys())
        n = len(components)

        # 检查是否为符号输入
        is_symbolic = any(isinstance(v, Symbol) for v in composition.values())

        # 遍历所有二元对
        for i in range(n):
            for j in range(i + 1, n):
                comp_i, comp_j = components[i], components[j]
                x_i, x_j = composition[comp_i], composition[comp_j]

                # 获取有效摩尔分数
                result = self.get_specified_binary_composition(
                    composition, (comp_i, comp_j), temperature, extrapolation_model
                )

                if result is None:
                    continue

                binary_sys, x_i_eff, x_j_eff = result
                binary_property = 0.0

                if is_symbolic:
                    # 符号计算路径
                    if property_name == 'excess_gibbs_energy':
                        params = binary_sys.aij_values if binary_sys.aij_values else binary_sys.hij_values
                    else:
                        params = binary_sys.hij_values

                    if params:
                        binary_sys._update_parameters_for_temperature(temperature)
                        std_c1, std_c2 = binary_sys._standardized_components
                        x_std1, x_std2 = (x_i_eff, x_j_eff) if comp_i == std_c1 else (x_j_eff, x_i_eff)
                        x_diff = x_std1 - x_std2
                        rk_sum = sum(p * (x_diff ** k) for k, p in enumerate(params))
                        binary_property = x_i_eff * x_j_eff * rk_sum

                else:
                    # 数值计算路径
                    if x_i_eff * x_j_eff != 0:
                        if property_name == 'excess_gibbs_energy':
                            val = binary_sys.excess_gibbs_energy(comp_i, x_i_eff, temperature)
                            if val is None:
                                # 回退到混合焓
                                val = binary_sys.mixing_enthalpy(comp_i, x_i_eff, temperature)
                                if val is not None:
                                    binary_property = val
                                    self.warnings.append(
                                        f"系统 {binary_sys.name}: gE数据缺失, 使用H_mix近似"
                                    )
                            else:
                                binary_property = val
                        else:
                            val = binary_sys.mixing_enthalpy(comp_i, x_i_eff, temperature)
                            if val is not None:
                                binary_property = val

                # 计算权重并累加
                if binary_property != 0:
                    denominator = x_i_eff * x_j_eff
                    if denominator != 0:
                        weight = (x_i * x_j) / denominator
                        total_property += weight * binary_property

        return total_property

    def calculate_activity_coefficient(
        self,
        composition: Dict[str, float],
        solute: str,
        solvent: str,
        temperature: float,
        extrapolation_model: Callable,
    ) -> float:
        """
        计算活度系数 ln(γi)

        基于热力学关系式: ln(γi) = [GE + ∂GE/∂xi - Σj(xj·∂GE/∂xj)] / RT

        Parameters
        ----------
        composition : dict
            组成字典
        solute : str
            溶质组分
        solvent : str
            溶剂组分（参考组分）
        temperature : float
            温度 (K)
        extrapolation_model : callable
            外推模型函数

        Returns
        -------
        float
            活度系数的自然对数 ln(γi)
        """
        # 输入验证
        if solute not in composition:
            raise ValueError(f"溶质 {solute} 不在组成字典中")
        if solvent not in composition:
            raise ValueError(f"溶剂 {solvent} 不在组成字典中")
        if temperature <= 0:
            raise ValueError("温度必须大于0")

        # 归一化组成
        normalized_comp = self._validate_and_normalize_composition(composition)

        # 单组分情况
        if len(normalized_comp) == 1:
            return 0.0

        # 创建符号变量
        independent_components = [c for c in normalized_comp if c != solvent]
        if not independent_components:
            return 0.0

        symbols_map = {name: symbols(name, real=True, positive=True) for name in independent_components}

        # 溶剂的符号表达式
        solvent_expression = 1 - sum(symbols_map.values())
        symbolic_comp = symbols_map.copy()
        symbolic_comp[solvent] = solvent_expression

        try:
            # 构造过剩Gibbs函数
            gmEx = self.get_excess_gibbs(symbolic_comp, temperature, extrapolation_model)

            # 检查是否为理想溶液
            if isinstance(gmEx, (int, float)) and gmEx == 0:
                return 0.0

            # 计算偏导数
            derivatives = {}
            for name, sym in symbols_map.items():
                try:
                    derivatives[name] = diff(gmEx, sym)
                except Exception as e:
                    self.warnings.append(f"计算组分 {name} 的偏导数失败: {e}")
                    return 0.0

            # 在实际组成点评估
            substitution = {symbols_map[name]: normalized_comp[name] for name in independent_components}

            # 评估过剩Gibbs能
            try:
                gmE_value = gmEx.subs(substitution)
                gmE_numerical = self._safe_float(gmE_value)
            except Exception as e:
                self.warnings.append(f"评估过剩Gibbs能失败: {e}")
                return 0.0

            # 计算加权偏导数之和
            sum_weighted_derivatives = 0.0
            for name in independent_components:
                try:
                    derivative_value = derivatives[name].subs(substitution)
                    derivative_numerical = self._safe_float(derivative_value)
                    sum_weighted_derivatives += normalized_comp[name] * derivative_numerical
                except Exception as e:
                    self.warnings.append(f"处理组分 {name} 的加权偏导数失败: {e}")
                    continue

            # 根据溶质类型计算活度系数
            if solute == solvent:
                ln_gamma_RT = gmE_numerical - sum_weighted_derivatives
            else:
                if solute in derivatives:
                    try:
                        dGE_dsolute = derivatives[solute].subs(substitution)
                        dGE_dsolute_numerical = self._safe_float(dGE_dsolute)
                        ln_gamma_RT = gmE_numerical + dGE_dsolute_numerical - sum_weighted_derivatives
                    except Exception as e:
                        self.warnings.append(f"计算溶质 {solute} 的偏导数失败: {e}")
                        return 0.0
                else:
                    ln_gamma_RT = gmE_numerical - sum_weighted_derivatives

            # 转换为ln(γ)
            ln_gamma = ln_gamma_RT / (R_GAS_CONSTANT * temperature)

            return self._safe_float(ln_gamma)

        except Exception as e:
            self.warnings.append(f"活度系数计算失败: {e}")
            return 0.0

    def calculate_all_activity_coefficients(
        self,
        composition: Dict[str, float],
        solvent: str,
        temperature: float,
        extrapolation_model: Callable,
    ) -> Dict[str, float]:
        """
        计算所有组分的活度系数

        Parameters
        ----------
        composition : dict
            组成字典
        solvent : str
            溶剂组分（参考组分）
        temperature : float
            温度 (K)
        extrapolation_model : callable
            外推模型函数

        Returns
        -------
        dict
            {组分: ln(γi)}
        """
        result = {}
        for component in composition.keys():
            try:
                ln_gamma = self.calculate_activity_coefficient(
                    composition, component, solvent, temperature, extrapolation_model
                )
                result[component] = ln_gamma
            except Exception as e:
                self.warnings.append(f"计算组分 {component} 的活度系数失败: {e}")
                result[component] = 0.0
        return result

    def calculate_activity(
        self,
        composition: Dict[str, float],
        solute: str,
        solvent: str,
        temperature: float,
        extrapolation_model: Callable,
    ) -> float:
        """
        计算组分的活度

        活度 = 摩尔分数 × 活度系数

        Parameters
        ----------
        composition : dict
            组成字典
        solute : str
            溶质组分
        solvent : str
            溶剂组分（参考组分）
        temperature : float
            温度 (K)
        extrapolation_model : callable
            外推模型函数

        Returns
        -------
        float
            活度
        """
        try:
            normalized_comp = self._validate_and_normalize_composition(composition)

            ln_gamma = self.calculate_activity_coefficient(
                composition, solute, solvent, temperature, extrapolation_model
            )
            x_i = normalized_comp[solute]
            gamma = math.exp(ln_gamma)
            return x_i * gamma

        except Exception as e:
            self.warnings.append(f"活度计算失败: {e}")
            # 返回理想溶液的活度
            return composition.get(solute, 0.0) / sum(composition.values())
