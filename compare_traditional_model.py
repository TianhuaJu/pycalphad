#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
对比测试：验证两个脚本中传统模型计算的差异

关键差异：
1. calculate_alcrni_uem.py: model_dict=None (默认模型)
2. calculate_liquidus_alcrni.py: models_muggianu = {ph: Model for ph in phases}
3. calculate_alcrni_uem.py: 温度范围 (1400, 2400, 10)
4. calculate_liquidus_alcrni.py: 温度范围 (1400, 2200, 10)
"""

import numpy as np
from pycalphad import Database, equilibrium, variables as v
from pycalphad.model import Model


def test_single_point(dbe, comps, phases, x_ni, T_range, model_dict, label):
    """
    测试单个成分点的液相线温度
    """
    x_al = (1.0 - x_ni) / 2.0
    x_cr = (1.0 - x_ni) / 2.0

    conditions = {
        v.T: T_range,
        v.P: 101325,
        v.X('AL'): x_al,
        v.X('CR'): x_cr,
    }

    try:
        eq = equilibrium(dbe, comps, phases, model=model_dict, conditions=conditions)

        T_vals = eq.T.values
        phase_array = eq.Phase.values
        np_array = eq.NP.values

        liquidus_T = None

        # 从高温到低温扫描
        for t_idx in range(len(T_vals) - 1, -1, -1):
            liquid_np = None
            for v_idx in range(phase_array.shape[-1]):
                phase_name = phase_array[0, 0, t_idx, 0, 0, v_idx]
                if phase_name == 'LIQUID':
                    liquid_np = np_array[0, 0, t_idx, 0, 0, v_idx]
                    break

            if liquid_np is not None and liquid_np >= 0.995:
                liquidus_T = T_vals[t_idx]
            elif liquid_np is not None and liquid_np < 0.995:
                break

        if liquidus_T is not None:
            print(f"{label}: X(NI)={x_ni:.3f} -> T_liquidus={liquidus_T:.1f} K")
            return liquidus_T
        else:
            print(f"{label}: X(NI)={x_ni:.3f} -> 未找到液相线")
            return None

    except Exception as e:
        print(f"{label}: X(NI)={x_ni:.3f} -> 计算错误: {e}")
        return None


def main():
    """
    对比两种传统模型的使用方式
    """
    print("=" * 80)
    print("传统 Muggianu 模型对比测试")
    print("=" * 80)
    print()

    # 加载数据库
    dbe = Database('examples/alcrni.tdb')
    comps = ['AL', 'CR', 'NI']
    phases = ['LIQUID', 'FCC_A1', 'BCC_A2', 'B2', 'L12_FCC']

    # 测试几个成分点
    test_points = [0.1, 0.3, 0.5, 0.7, 0.9]

    print("\n差异1: 模型字典的使用方式")
    print("-" * 80)
    print("方式A (calculate_alcrni_uem.py):     model_dict=None")
    print("方式B (calculate_liquidus_alcrni.py): models = {ph: Model for ph in phases}")
    print("-" * 80)
    print()

    # 使用相同的温度范围测试
    T_range_common = (1400, 2200, 10)

    results_none = []
    results_dict = []
    models_explicit = {ph: Model for ph in phases}

    for x_ni in test_points:
        print(f"\n{'=' * 80}")
        print(f"测试点: X(NI) = {x_ni:.3f}")
        print(f"{'=' * 80}")

        # 方式A: model_dict=None
        t_none = test_single_point(
            dbe, comps, phases, x_ni, T_range_common,
            model_dict=None,
            label="方式A (model=None)"
        )
        results_none.append(t_none)

        # 方式B: 明确指定 Model 字典
        t_dict = test_single_point(
            dbe, comps, phases, x_ni, T_range_common,
            model_dict=models_explicit,
            label="方式B (model={ph: Model})"
        )
        results_dict.append(t_dict)

        # 显示差异
        if t_none is not None and t_dict is not None:
            diff = t_dict - t_none
            print(f"  → 差异 (B - A): {diff:.2f} K")
            if abs(diff) > 1.0:
                print(f"  ⚠️  警告：差异超过 1K！")
        print()

    # 统计结果
    print("\n" + "=" * 80)
    print("统计结果")
    print("=" * 80)

    differences = []
    for i, x_ni in enumerate(test_points):
        if results_none[i] is not None and results_dict[i] is not None:
            diff = results_dict[i] - results_none[i]
            differences.append(diff)
            print(f"X(NI)={x_ni:.3f}: 方式A={results_none[i]:.1f}K, "
                  f"方式B={results_dict[i]:.1f}K, 差异={diff:.2f}K")

    if len(differences) > 0:
        differences = np.array(differences)
        print(f"\n差异统计:")
        print(f"  平均值: {np.mean(differences):.4f} K")
        print(f"  最大值: {np.max(differences):.4f} K")
        print(f"  最小值: {np.min(differences):.4f} K")
        print(f"  标准差: {np.std(differences):.4f} K")

        if np.max(np.abs(differences)) < 0.01:
            print("\n✅ 结论: model=None 和 model={ph: Model} 给出完全相同的结果")
        else:
            print(f"\n⚠️  结论: 存在差异！最大差异 {np.max(np.abs(differences)):.2f} K")

    # 测试温度范围的影响
    print("\n\n" + "=" * 80)
    print("差异2: 温度范围的影响")
    print("-" * 80)
    print("范围A (calculate_alcrni_uem.py):     (1400, 2400, 10)")
    print("范围B (calculate_liquidus_alcrni.py): (1400, 2200, 10)")
    print("-" * 80)
    print()

    # 选择一个可能受影响的成分点（低 Ni，高 Cr-Al，液相线可能较高）
    x_ni_test = 0.1

    print(f"测试成分: X(NI)={x_ni_test:.3f}, X(AL)={0.45:.3f}, X(CR)={0.45:.3f}")
    print()

    T_range_wide = (1400, 2400, 10)
    T_range_narrow = (1400, 2200, 10)

    t_wide = test_single_point(
        dbe, comps, phases, x_ni_test, T_range_wide,
        model_dict=None,
        label="范围A (1400-2400K)"
    )

    t_narrow = test_single_point(
        dbe, comps, phases, x_ni_test, T_range_narrow,
        model_dict=None,
        label="范围B (1400-2200K)"
    )

    if t_wide is not None and t_narrow is not None:
        diff = t_wide - t_narrow
        print(f"\n差异: {diff:.2f} K")
        if abs(diff) < 0.01:
            print("✅ 温度范围不影响此成分点的液相线计算")
        else:
            print(f"⚠️  温度范围影响液相线计算！差异 {diff:.2f} K")
    elif t_wide is not None and t_narrow is None:
        print(f"\n⚠️  范围B 无法找到液相线，但范围A 找到了 {t_wide:.1f} K")
        print(f"    这意味着液相线高于 2200K，需要更宽的温度范围！")
    elif t_wide == t_narrow:
        print("\n✅ 两种温度范围给出相同结果")

    print("\n" + "=" * 80)
    print("测试完成！")
    print("=" * 80)


if __name__ == "__main__":
    main()
