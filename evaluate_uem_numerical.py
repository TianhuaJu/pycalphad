#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
UEM 数值评估工具 - 直接计算并打印贡献系数的数值

更简洁的版本，直接计算数值而不是符号表达式
"""

import numpy as np
from pycalphad import Database
import pycalphad.variables as v
from itertools import combinations


def calculate_d_term_numerical(dbe, phase_name, comp_k, comp_i, subl_index, temperature):
    """
    计算 d_ki 的数值

    d_ki = (2/RT) * sum(L_ki^{(odd)})
    """
    from tinydb import where

    R = 8.31446  # J/(mol·K)

    # 确保组分是 Species 对象
    if not isinstance(comp_k, v.Species):
        comp_k = v.Species(comp_k)
    if not isinstance(comp_i, v.Species):
        comp_i = v.Species(comp_i)

    # 查找二元参数
    param_query = (
        (where('phase_name') == phase_name) &
        ((where('parameter_type') == 'G') | (where('parameter_type') == 'L'))
    )
    params = dbe.search(param_query)

    # 筛选出 (comp_k, comp_i) 或 (comp_i, comp_k) 的参数
    odd_L_values = []
    first_comp_in_db = None

    for param in params:
        const_array = param.get('constituent_array', [])

        # 检查是否是单亚晶格二元交互
        if len(const_array) > subl_index:
            subl_comps = const_array[subl_index]

            if len(subl_comps) == 2:
                # 将元组转换为 Species 对象
                comp_a = subl_comps[0] if isinstance(subl_comps[0], v.Species) else v.Species(subl_comps[0])
                comp_b = subl_comps[1] if isinstance(subl_comps[1], v.Species) else v.Species(subl_comps[1])

                # 检查是否匹配 (comp_k, comp_i)
                if {comp_a, comp_b} == {comp_k, comp_i}:
                    order = param.get('parameter_order', 0)

                    # 记录第一个组分（用于符号修正）
                    if first_comp_in_db is None:
                        first_comp_in_db = comp_a

                    if order % 2 != 0:  # 奇数阶
                        L_expr = param['parameter']
                        # 计算数值
                        try:
                            L_value = float(L_expr.xreplace({v.T: temperature}))
                            odd_L_values.append(L_value)
                        except:
                            # 如果无法转换为float，尝试 evalf
                            try:
                                L_value = float(L_expr.subs({v.T: temperature}))
                                odd_L_values.append(L_value)
                            except:
                                pass

    if not odd_L_values:
        return 0.0, first_comp_in_db

    # 符号修正
    sign_correction = 1
    if first_comp_in_db is not None:
        if first_comp_in_db == comp_k:
            sign_correction = 1
        elif first_comp_in_db == comp_i:
            sign_correction = -1

    # ✅ 关键修正：取绝对值！
    # d_ki = (1/RT) * |g_i^∞(in k) - g_k^∞(in i)|
    # 对于 R-K 模型: d_ki = (2/RT) * |sum(L_odd)|
    d_ki = abs((2.0 / (R * temperature)) * sign_correction * sum(odd_L_values))

    return d_ki, first_comp_in_db


def evaluate_uem_coefficients(dbe, phase_name, composition, temperature, verbose=True):
    """
    在指定成分点数值计算 UEM 贡献系数

    参数:
    - dbe: Database 对象
    - phase_name: 相名称
    - composition: 成分字典，例如 {'AL': 0.35, 'CR': 0.35, 'NI': 0.30}
    - temperature: 温度 (K)
    - verbose: 是否打印详细信息

    返回:
    - coefficients: 贡献系数数据列表
    """

    if verbose:
        print("\n" + "=" * 90)
        print("UEM 贡献系数数值计算")
        print("=" * 90)
        print(f"相: {phase_name}")
        print(f"温度: {temperature} K")
        print(f"成分: {composition}")
        print("=" * 90)

    comps = list(composition.keys())
    subl_index = 0  # 假设单亚晶格

    coefficients = []

    # 遍历所有二元对
    for comp_i, comp_j in combinations(comps, 2):
        x_i = composition[comp_i]
        x_j = composition[comp_j]

        if verbose:
            print(f"\n{'─' * 90}")
            print(f"二元对: {comp_i} - {comp_j}")
            print(f"  x_{comp_i} = {x_i:.4f}, x_{comp_j} = {x_j:.4f}")
            print(f"{'─' * 90}")

        # 找到第三组分
        other_comps = [c for c in comps if c != comp_i and c != comp_j]

        if not other_comps:
            if verbose:
                print("  ℹ️  纯二元系，无第三组分")
            continue

        for comp_k in other_comps:
            x_k = composition[comp_k]

            if verbose:
                print(f"\n  第三组分: {comp_k}, x_{comp_k} = {x_k:.4f}")

            # 计算 d_ki 和 d_kj
            d_ki, first_comp_ki = calculate_d_term_numerical(
                dbe, phase_name, comp_k, comp_i, subl_index, temperature
            )

            d_kj, first_comp_kj = calculate_d_term_numerical(
                dbe, phase_name, comp_k, comp_j, subl_index, temperature
            )

            if verbose:
                print(f"    d_{comp_k}{comp_i} = {d_ki:.6f}")
                print(f"    d_{comp_k}{comp_j} = {d_kj:.6f}")

            # 计算 d_sum
            d_sum = d_ki + d_kj

            if verbose:
                print(f"    d_sum = {d_sum:.6f}")

            # 计算 alpha 系数
            if abs(d_sum) < 1e-10:
                alpha_i_k = 0.0
                alpha_j_k = 0.0
                if verbose:
                    print(f"    ⚠️  d_sum ≈ 0, 设置 α = 0")
            else:
                alpha_i_k = (d_kj / d_sum) * np.exp(-d_ki)
                alpha_j_k = (d_ki / d_sum) * np.exp(-d_kj)

            if verbose:
                print(f"    α_{comp_i}^{comp_k} = (d_kj / d_sum) * exp(-d_ki) = {alpha_i_k:.6f}")
                print(f"    α_{comp_j}^{comp_k} = (d_ki / d_sum) * exp(-d_kj) = {alpha_j_k:.6f}")

            # 检查异常值
            if abs(alpha_i_k) > 100 or abs(alpha_j_k) > 100:
                if verbose:
                    print(f"    ⚠️⚠️⚠️  警告：贡献系数异常大！")

            if np.isnan(alpha_i_k) or np.isnan(alpha_j_k):
                if verbose:
                    print(f"    ⚠️⚠️⚠️  错误：贡献系数为 NaN！")

            # 存储结果
            coefficients.append({
                'binary_pair': (comp_i, comp_j),
                'third_comp': comp_k,
                'x_i': x_i,
                'x_j': x_j,
                'x_k': x_k,
                'd_ki': d_ki,
                'd_kj': d_kj,
                'd_sum': d_sum,
                'alpha_i_k': alpha_i_k,
                'alpha_j_k': alpha_j_k,
            })

            # 计算等效组分
            X_ij_i_num = x_i + alpha_i_k * x_k
            X_ij_j_num = x_j + alpha_j_k * x_k
            X_denom = x_i + x_j + (alpha_i_k + alpha_j_k) * x_k

            if abs(X_denom) < 1e-10:
                X_ij_i = 0.0
                X_ij_j = 0.0
            else:
                X_ij_i = X_ij_i_num / X_denom
                X_ij_j = X_ij_j_num / X_denom

            if verbose:
                print(f"\n  等效组分:")
                print(f"    X_{comp_i}{comp_j}^{comp_i} = {X_ij_i:.6f}")
                print(f"    X_{comp_i}{comp_j}^{comp_j} = {X_ij_j:.6f}")

            # 检查等效组分是否合理（应该在 0-1 之间）
            if X_ij_i < -0.01 or X_ij_i > 1.01 or X_ij_j < -0.01 or X_ij_j > 1.01:
                if verbose:
                    print(f"    ⚠️⚠️⚠️  警告：等效组分超出 [0, 1] 范围！")

            coefficients[-1]['X_ij_i'] = X_ij_i
            coefficients[-1]['X_ij_j'] = X_ij_j

    if verbose:
        print("\n" + "=" * 90)
        print("计算完成")
        print("=" * 90 + "\n")

    return coefficients


def analyze_problematic_points(dbe, phase_name, temperature_range, composition_grid):
    """
    分析多个成分点，识别异常的贡献系数

    参数:
    - dbe: Database 对象
    - phase_name: 相名称
    - temperature_range: 温度列表
    - composition_grid: 成分点列表
    """

    print("\n" + "=" * 90)
    print("UEM 异常点分析")
    print("=" * 90)

    problematic_points = []

    for temp in temperature_range:
        for comp in composition_grid:
            print(f"\n检查点: T={temp}K, 成分={comp}")

            coeffs = evaluate_uem_coefficients(dbe, phase_name, comp, temp, verbose=False)

            # 检查是否有异常值
            has_problem = False
            for c in coeffs:
                if abs(c['alpha_i_k']) > 100 or abs(c['alpha_j_k']) > 100:
                    has_problem = True
                if np.isnan(c['alpha_i_k']) or np.isnan(c['alpha_j_k']):
                    has_problem = True
                if c['X_ij_i'] < -0.01 or c['X_ij_i'] > 1.01:
                    has_problem = True
                if c['X_ij_j'] < -0.01 or c['X_ij_j'] > 1.01:
                    has_problem = True

            if has_problem:
                print(f"  ⚠️  发现异常！")
                problematic_points.append({
                    'temperature': temp,
                    'composition': comp,
                    'coefficients': coeffs
                })

                # 打印详细信息
                evaluate_uem_coefficients(dbe, phase_name, comp, temp, verbose=True)
            else:
                print(f"  ✓ 正常")

    print("\n" + "=" * 90)
    print(f"分析完成：共发现 {len(problematic_points)} 个异常点")
    print("=" * 90 + "\n")

    return problematic_points


# ============================================================================
# 主程序
# ============================================================================

def main():
    """主程序"""

    # 加载数据库
    dbe = Database('examples/alcrni.tdb')
    phase_name = 'LIQUID'

    # 测试1: 单点详细分析
    print("\n【测试1】单点详细分析")
    print("=" * 90)

    composition = {'AL': 0.35, 'CR': 0.35, 'NI': 0.30}
    temperature = 1700

    coeffs = evaluate_uem_coefficients(dbe, phase_name, composition, temperature, verbose=True)

    # 测试2: 用户自定义成分点
    print("\n【测试2】自定义成分点分析")
    print("=" * 90)

    user_input = input("\n是否测试自定义成分点？(y/n): ").strip().lower()

    if user_input == 'y':
        try:
            x_al = float(input("请输入 X(AL) (0-1): "))
            x_cr = float(input("请输入 X(CR) (0-1): "))
            x_ni = float(input("请输入 X(NI) (0-1): "))
            temp = float(input("请输入温度 (K): "))

            if abs(x_al + x_cr + x_ni - 1.0) > 0.001:
                print(f"⚠️  警告：成分总和 = {x_al + x_cr + x_ni:.4f} ≠ 1.0")

            custom_comp = {'AL': x_al, 'CR': x_cr, 'NI': x_ni}
            coeffs = evaluate_uem_coefficients(dbe, phase_name, custom_comp, temp, verbose=True)

        except ValueError as e:
            print(f"输入错误: {e}")

    # 测试3: 批量扫描
    print("\n【测试3】批量扫描识别异常点")
    print("=" * 90)

    user_input = input("\n是否进行批量扫描？(y/n): ").strip().lower()

    if user_input == 'y':
        # 生成成分网格
        x_ni_vals = np.linspace(0.2, 0.8, 5)
        comp_grid = []

        for x_ni in x_ni_vals:
            x_al = (1.0 - x_ni) / 2.0
            x_cr = (1.0 - x_ni) / 2.0
            comp_grid.append({'AL': x_al, 'CR': x_cr, 'NI': x_ni})

        temp_range = [1600, 1700, 1800]

        problems = analyze_problematic_points(dbe, phase_name, temp_range, comp_grid)

        if problems:
            print("\n异常点汇总:")
            print("=" * 90)
            for p in problems:
                print(f"\nT={p['temperature']}K, 成分={p['composition']}")
                for c in p['coefficients']:
                    if abs(c['alpha_i_k']) > 100 or abs(c['alpha_j_k']) > 100:
                        print(f"  ⚠️  {c['binary_pair']} - {c['third_comp']}: "
                              f"α_i={c['alpha_i_k']:.2f}, α_j={c['alpha_j_k']:.2f}")

    print("\n测试完成！")


if __name__ == "__main__":
    main()
