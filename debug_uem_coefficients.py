#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
UEM 模型调试工具 - 打印贡献系数和中间计算值

用于调试 UEM 模型中化学势（活度）异常的问题。
打印所有关键的中间变量：d_ki, d_kj, alpha 系数, X_ij 等。
"""

import numpy as np
from pycalphad import Database, equilibrium, variables as v
from pycalphad.UEMModel import UEMModel
from symengine import sympify
import pycalphad.variables as v


class UEMModelDebug(UEMModel):
    """
    继承 UEMModel，添加调试输出功能
    """

    def __init__(self, dbe, comps, phase_name, verbose=True, debug_composition=None):
        """
        参数：
        - verbose: 是否打印调试信息
        - debug_composition: 需要调试的特定成分点（字典格式），例如 {'AL': 0.3, 'CR': 0.3, 'NI': 0.4}
        """
        super().__init__(dbe, comps, phase_name)
        self.verbose = verbose
        self.debug_composition = debug_composition
        self.debug_data = []  # 存储调试数据

    def excess_mixing_energy(self, dbe):
        """
        重写 excess_mixing_energy，添加调试输出
        """
        from itertools import combinations
        from symengine import exp, Add, Piecewise, S

        # 1. 建立二元 L 参数的缓存
        if not hasattr(self, '_binary_L_param_cache'):
            self._build_binary_L_cache(dbe)

        total_excess_energy = S.Zero
        phase = dbe.phases[self.phase_name]

        if self.verbose:
            print("\n" + "=" * 80)
            print(f"UEM 调试输出 - 相: {self.phase_name}")
            print("=" * 80)

        # 2. 遍历所有亚晶格
        for subl_index, sublattice_comps in enumerate(self.constituents):
            active_comps = sorted(list(sublattice_comps.intersection(self.components)))

            if len(active_comps) < 2:
                continue  # 该亚晶格没有混合

            if self.verbose:
                print(f"\n亚晶格 {subl_index}: 活性组分 = {active_comps}")

            # 获取该亚晶格的位组分数 (Site Fractions)
            site_fracs = {comp: v.SiteFraction(self.phase_name, subl_index, comp) for comp in active_comps}

            sublattice_total_G_E = S.Zero

            # 3. 遍历所有二元组 (i, j)
            for comp_i, comp_j in combinations(active_comps, 2):
                if self.verbose:
                    print(f"\n  {'=' * 76}")
                    print(f"  二元对: {comp_i} - {comp_j}")
                    print(f"  {'=' * 76}")

                x_i = site_fracs[comp_i]
                x_j = site_fracs[comp_j]

                # 4. 构建 G_ij^E(x_i, x_j) 二元 R-K 表达式
                G_ij_binary_expr_rk_sum = S.Zero
                key = (subl_index, tuple(sorted((comp_i, comp_j))))
                binary_params = self._binary_L_param_cache.get(key, {}).get('all', [])

                if not binary_params:
                    if self.verbose:
                        print(f"    ⚠️  没有 {comp_i}-{comp_j} 交互参数")
                    continue  # 没有 (i, j) 交互

                # 获取参数定义时的组元顺序
                p_i_species, p_j_species = key[1]
                p_i = site_fracs[p_i_species]
                p_j = site_fracs[p_j_species]

                for L_expr, order in binary_params:
                    G_ij_binary_expr_rk_sum += L_expr * (p_i - p_j) ** order

                # 5. 计算所有其他组元 k 的 d_ki 和 d_kj
                X_ij_i_num = x_i
                X_ij_j_num = x_j
                X_denom = x_i + x_j

                other_comps = [c for c in active_comps if c != comp_i and c != comp_j]

                if not other_comps:
                    # 纯二元系，无需外推
                    sublattice_total_G_E += p_i * p_j * G_ij_binary_expr_rk_sum
                    if self.verbose:
                        print(f"    ℹ️  纯二元系，无需 UEM 外推")
                    continue

                if self.verbose:
                    print(f"    第三组分: {other_comps}")
                    print(f"    {'─' * 72}")

                for comp_k in other_comps:
                    x_k = site_fracs[comp_k]

                    # 6. 计算 d_ki, d_kj (注意组元顺序)
                    d_ki = self._get_uem_d_term(dbe, comp_k, comp_i, subl_index)
                    d_kj = self._get_uem_d_term(dbe, comp_k, comp_j, subl_index)

                    if self.verbose:
                        print(f"\n    组分 k = {comp_k}:")
                        print(f"      d({comp_k},{comp_i}) = {d_ki}")
                        print(f"      d({comp_k},{comp_j}) = {d_kj}")

                    # 7. 计算 alpha 系数
                    d_sum = d_ki + d_kj
                    d_sum_safe = Piecewise((d_sum, d_sum != 0), (1, True))

                    alpha_i_k = (d_kj / d_sum_safe) * exp(-d_ki)
                    alpha_j_k = (d_ki / d_sum_safe) * exp(-d_kj)

                    if self.verbose:
                        print(f"      d_sum = d({comp_k},{comp_i}) + d({comp_k},{comp_j}) = {d_sum}")
                        print(f"      α({comp_i}|{comp_k}) = (d_kj / d_sum) * exp(-d_ki) = {alpha_i_k}")
                        print(f"      α({comp_j}|{comp_k}) = (d_ki / d_sum) * exp(-d_kj) = {alpha_j_k}")

                    # 存储调试数据
                    self.debug_data.append({
                        'binary_pair': (str(comp_i), str(comp_j)),
                        'third_comp': str(comp_k),
                        'd_ki': str(d_ki),
                        'd_kj': str(d_kj),
                        'd_sum': str(d_sum),
                        'alpha_i_k': str(alpha_i_k),
                        'alpha_j_k': str(alpha_j_k)
                    })

                    # 8. 累加等效组分
                    X_ij_i_num += alpha_i_k * x_k
                    X_ij_j_num += alpha_j_k * x_k
                    X_denom += (alpha_i_k + alpha_j_k) * x_k

                # 9. 计算 X_ij^i 和 X_ij^j
                X_denom_safe = Piecewise((X_denom, X_denom != 0), (1, True))
                X_ij_i = X_ij_i_num / X_denom_safe
                X_ij_j = X_ij_j_num / X_denom_safe

                if self.verbose:
                    print(f"\n    等效组分:")
                    print(f"      X_ij^{comp_i} = {X_ij_i}")
                    print(f"      X_ij^{comp_j} = {X_ij_j}")

                # 10. 代入 G_ij^E 表达式
                substitution_dict = {}
                if p_i_species == comp_i:
                    substitution_dict[p_i] = X_ij_i
                    substitution_dict[p_j] = X_ij_j
                else:
                    substitution_dict[p_i] = X_ij_j
                    substitution_dict[p_j] = X_ij_i

                rk_sum_modified = G_ij_binary_expr_rk_sum.xreplace(substitution_dict)

                # 11. 根据UEM公式计算最终贡献
                final_term = (x_i * x_j) * rk_sum_modified
                sublattice_total_G_E += final_term

                if self.verbose:
                    print(f"    最终贡献: (x_{comp_i} * x_{comp_j}) * G_ij^E(X_ij) = {final_term}")

            total_excess_energy += sublattice_total_G_E

        if self.verbose:
            print("\n" + "=" * 80)
            print(f"总过剩能: {total_excess_energy}")
            print("=" * 80 + "\n")

        # 12. 返回总能量，并按照基类的要求进行归一化
        return total_excess_energy / self._site_ratio_normalization


def evaluate_uem_at_point(dbe, phase_name, comps, composition, temperature, pressure=101325):
    """
    在特定成分点计算 UEM 模型，并打印详细的调试信息

    参数：
    - dbe: Database 对象
    - phase_name: 相名称，例如 'LIQUID'
    - comps: 组分列表，例如 ['AL', 'CR', 'NI']
    - composition: 成分字典，例如 {'AL': 0.35, 'CR': 0.35, 'NI': 0.30}
    - temperature: 温度 (K)
    - pressure: 压力 (Pa)，默认 101325

    返回：
    - excess_energy_value: 过剩能数值 (J/mol)
    - debug_data: 调试数据列表
    """

    print("\n" + "=" * 90)
    print("UEM 模型调试 - 在特定成分点评估")
    print("=" * 90)
    print(f"相: {phase_name}")
    print(f"温度: {temperature} K")
    print(f"压力: {pressure} Pa")
    print(f"成分: {composition}")
    print("=" * 90)

    # 创建调试模型
    model_debug = UEMModelDebug(dbe, comps + ['VA'], phase_name, verbose=True)

    # 获取过剩能表达式
    excess_expr = model_debug.excess_mixing_energy(dbe)

    # 准备替换字典（将符号替换为数值）
    subs_dict = {
        v.T: temperature,
        v.P: pressure,
        v.R: 8.31446
    }

    # 添加成分替换
    for comp, x_val in composition.items():
        subs_dict[v.SiteFraction(phase_name, 0, v.Species(comp))] = x_val

    # 计算数值
    try:
        excess_energy_value = float(excess_expr.xreplace(subs_dict))

        print("\n" + "=" * 90)
        print("数值结果")
        print("=" * 90)
        print(f"过剩能 G^E = {excess_energy_value:.2f} J/mol")
        print("=" * 90 + "\n")

        return excess_energy_value, model_debug.debug_data

    except Exception as e:
        print(f"\n⚠️  计算数值时出错: {e}")
        print(f"表达式: {excess_expr}")
        return None, model_debug.debug_data


def scan_composition_range(dbe, phase_name, comps, base_comp1, base_comp2,
                           scan_comp, num_points=10, temperature=1700):
    """
    扫描成分范围，识别异常的化学势/活度

    参数：
    - base_comp1, base_comp2: 保持比例的两个组分（例如 'AL', 'CR'）
    - scan_comp: 扫描的组分（例如 'NI'）
    - num_points: 扫描点数
    """

    print("\n" + "=" * 90)
    print("UEM 模型成分扫描")
    print("=" * 90)
    print(f"固定比例: {base_comp1}/{base_comp2} = 1/1")
    print(f"扫描组分: {scan_comp}")
    print(f"温度: {temperature} K")
    print("=" * 90 + "\n")

    scan_range = np.linspace(0.1, 0.8, num_points)

    results = []

    for x_scan in scan_range:
        x_base1 = (1.0 - x_scan) / 2.0
        x_base2 = (1.0 - x_scan) / 2.0

        composition = {
            base_comp1: x_base1,
            base_comp2: x_base2,
            scan_comp: x_scan
        }

        print(f"\n{'━' * 90}")
        print(f"测试点: X({scan_comp})={x_scan:.3f}, X({base_comp1})={x_base1:.3f}, X({base_comp2})={x_base2:.3f}")
        print(f"{'━' * 90}")

        g_excess, debug_data = evaluate_uem_at_point(
            dbe, phase_name, comps, composition, temperature
        )

        results.append({
            'composition': composition.copy(),
            'g_excess': g_excess,
            'debug_data': debug_data
        })

        # 检查异常值
        if g_excess is not None and abs(g_excess) > 100000:
            print(f"\n⚠️⚠️⚠️  警告：过剩能异常大！G^E = {g_excess:.2f} J/mol")
            print("详细调试数据:")
            for data in debug_data:
                print(f"  二元对 {data['binary_pair']}, 第三组分 {data['third_comp']}:")
                print(f"    d_ki = {data['d_ki']}")
                print(f"    d_kj = {data['d_kj']}")
                print(f"    α_i_k = {data['alpha_i_k']}")
                print(f"    α_j_k = {data['alpha_j_k']}")

    print("\n" + "=" * 90)
    print("成分扫描完成")
    print("=" * 90)

    return results


# ============================================================================
# 主测试函数
# ============================================================================

def main():
    """主测试函数"""

    print("\n" + "=" * 90)
    print("UEM 贡献系数调试工具")
    print("=" * 90)
    print()

    # 加载数据库
    dbe = Database('examples/alcrni.tdb')
    comps = ['AL', 'CR', 'NI']
    phase_name = 'LIQUID'

    # 测试1: 单个成分点的详细调试
    print("\n【测试1】单个成分点的详细调试")
    print("-" * 90)

    composition_test = {'AL': 0.35, 'CR': 0.35, 'NI': 0.30}
    temperature_test = 1700

    g_excess, debug_data = evaluate_uem_at_point(
        dbe, phase_name, comps, composition_test, temperature_test
    )

    # 测试2: 成分扫描（可选）
    print("\n【测试2】成分扫描（识别异常点）")
    print("-" * 90)

    user_input = input("是否进行成分扫描？(y/n): ").strip().lower()

    if user_input == 'y':
        results = scan_composition_range(
            dbe, phase_name, comps,
            base_comp1='AL',
            base_comp2='CR',
            scan_comp='NI',
            num_points=5,
            temperature=1700
        )

        # 汇总结果
        print("\n" + "=" * 90)
        print("汇总结果")
        print("=" * 90)
        print(f"{'X(NI)':<10} {'X(AL)':<10} {'X(CR)':<10} {'G_excess (J/mol)':<20}")
        print("-" * 90)

        for res in results:
            comp = res['composition']
            g = res['g_excess']
            if g is not None:
                flag = "  ⚠️  异常" if abs(g) > 100000 else ""
                print(f"{comp['NI']:<10.3f} {comp['AL']:<10.3f} {comp['CR']:<10.3f} {g:<20.2f}{flag}")

        print("=" * 90)

    print("\n测试完成！")


if __name__ == "__main__":
    main()
