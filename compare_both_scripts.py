#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
精确对比两个测试脚本的结果

在完全相同的成分点上运行两个脚本的计算逻辑，
对比传统 Model 和 UEM 模型的结果差异。
"""

import numpy as np
from pycalphad import Database, equilibrium, variables as v
from pycalphad.model import Model
from pycalphad.uem1_Model import uem1_model
from pycalphad.advanced_uem_model import ModelUEM1


def calculate_liquidus_single(dbe, comps, phases, x_ni, T_range, model_dict, label):
    """计算单个成分点的液相线温度"""
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

        return liquidus_T

    except Exception as e:
        print(f"{label} - X(NI)={x_ni:.3f} 计算错误: {e}")
        return None


def main():
    """主对比函数"""
    print("=" * 90)
    print("两个测试脚本的精确对比")
    print("=" * 90)
    print()

    # 加载数据库
    dbe = Database('examples/alcrni.tdb')
    comps = ['AL', 'CR', 'NI']
    phases = ['LIQUID', 'FCC_A1', 'BCC_A2', 'B2', 'L12_FCC']

    # 使用相同的成分点
    x_ni_range = np.linspace(0.1, 0.9, 10)

    # 相同的温度范围
    T_range = (1400, 2400, 10)

    print("测试条件:")
    print(f"  成分范围: X(NI) = 0.1 to 0.9 ({len(x_ni_range)} 个点)")
    print(f"  温度范围: {T_range[0]}-{T_range[1]} K，步长 {T_range[2]} K")
    print(f"  相: {', '.join(phases)}")
    print()

    # 准备模型
    print("准备模型...")

    # 1. 传统 Model (None)
    model_traditional_none = None

    # 2. 传统 Model ({ph: Model})
    model_traditional_dict = {ph: Model for ph in phases}

    # 3. uem1_model - 所有相 (calculate_alcrni_uem.py 的方式)
    uem_models_all = {}
    for phase_name in phases:
        try:
            uem_models_all[phase_name] = uem1_model(dbe, comps + ['VA'], phase_name)
        except:
            pass

    # 4. ModelUEM1 - 只有液相 (calculate_liquidus_alcrni.py 的方式)
    models_uem1_liquid_only = {ph: ModelUEM1 if ph == 'LIQUID' else Model for ph in phases}

    print(f"  ✓ 传统 Model (None)")
    print(f"  ✓ 传统 Model ({{ph: Model}})")
    print(f"  ✓ uem1_model - 所有相 ({len(uem_models_all)} 个相)")
    print(f"  ✓ ModelUEM1 - 只有液相")
    print()

    # 存储结果
    results = {
        '传统Model(None)': [],
        '传统Model(Dict)': [],
        'uem1_model(所有相)': [],
        'ModelUEM1(液相)': []
    }

    print("=" * 90)
    print("开始计算...")
    print("=" * 90)
    print()

    for x_ni in x_ni_range:
        x_al = (1.0 - x_ni) / 2.0
        x_cr = (1.0 - x_ni) / 2.0

        print(f"X(NI)={x_ni:.3f}, X(AL)={x_al:.3f}, X(CR)={x_cr:.3f}:")

        # 计算各个模型
        t1 = calculate_liquidus_single(dbe, comps, phases, x_ni, T_range,
                                      model_traditional_none, "传统(None)")
        results['传统Model(None)'].append(t1)

        t2 = calculate_liquidus_single(dbe, comps, phases, x_ni, T_range,
                                      model_traditional_dict, "传统(Dict)")
        results['传统Model(Dict)'].append(t2)

        t3 = calculate_liquidus_single(dbe, comps, phases, x_ni, T_range,
                                      uem_models_all, "UEM(所有相)")
        results['uem1_model(所有相)'].append(t3)

        t4 = calculate_liquidus_single(dbe, comps, phases, x_ni, T_range,
                                      models_uem1_liquid_only, "UEM1(液相)")
        results['ModelUEM1(液相)'].append(t4)

        # 打印结果
        print(f"  传统Model(None):     {t1:.1f} K" if t1 else "  传统Model(None):     未找到")
        print(f"  传统Model(Dict):     {t2:.1f} K" if t2 else "  传统Model(Dict):     未找到")
        print(f"  uem1_model(所有相):    {t3:.1f} K" if t3 else "  uem1_model(所有相):    未找到")
        print(f"  ModelUEM1(液相):     {t4:.1f} K" if t4 else "  ModelUEM1(液相):     未找到")

        # 计算差异
        if t1 and t2:
            diff_trad = t2 - t1
            if abs(diff_trad) > 0.1:
                print(f"  ⚠️  传统模型差异: {diff_trad:.2f} K")

        if t3 and t4:
            diff_uem = t3 - t4
            print(f"  UEM差异(所有相-液相): {diff_uem:.2f} K")

        print()

    # 统计分析
    print("=" * 90)
    print("统计分析")
    print("=" * 90)
    print()

    # 对比1: 传统 Model (None) vs (Dict)
    print("【对比1】传统 Model: None vs {ph: Model}")
    print("-" * 90)
    diff_traditional = []
    for i in range(len(x_ni_range)):
        t1 = results['传统Model(None)'][i]
        t2 = results['传统Model(Dict)'][i]
        if t1 and t2:
            diff = t2 - t1
            diff_traditional.append(diff)
            if abs(diff) > 0.1:
                print(f"  X(NI)={x_ni_range[i]:.3f}: {t1:.1f} K vs {t2:.1f} K, 差异 {diff:.2f} K")

    if len(diff_traditional) > 0:
        diff_traditional = np.array(diff_traditional)
        print(f"\n统计:")
        print(f"  平均差异: {np.mean(np.abs(diff_traditional)):.4f} K")
        print(f"  最大差异: {np.max(np.abs(diff_traditional)):.4f} K")
        if np.max(np.abs(diff_traditional)) < 0.01:
            print("  ✅ 结论: 两种方式完全一致")
        else:
            print("  ⚠️  结论: 存在差异")
    print()

    # 对比2: uem1_model (所有相) vs ModelUEM1 (液相)
    print("【对比2】UEM模型: uem1_model(所有相) vs ModelUEM1(只液相)")
    print("-" * 90)
    diff_uem = []
    for i in range(len(x_ni_range)):
        t3 = results['uem1_model(所有相)'][i]
        t4 = results['ModelUEM1(液相)'][i]
        if t3 and t4:
            diff = t3 - t4
            diff_uem.append(diff)
            print(f"  X(NI)={x_ni_range[i]:.3f}: {t3:.1f} K vs {t4:.1f} K, 差异 {diff:.2f} K")

    if len(diff_uem) > 0:
        diff_uem = np.array(diff_uem)
        print(f"\n统计:")
        print(f"  平均差异: {np.mean(diff_uem):.2f} K")
        print(f"  最大差异: {np.max(diff_uem):.2f} K")
        print(f"  最小差异: {np.min(diff_uem):.2f} K")
        print(f"  标准差:   {np.std(diff_uem):.2f} K")
        print("\n说明:")
        print("  uem1_model(所有相) 将 UEM 应用于所有相（包括固相），")
        print("  这会改变固相的稳定性，从而显著影响液相线位置。")
        print("  ModelUEM1(只液相) 仅将 UEM 应用于液相，固相使用传统 Muggianu，")
        print("  这样可以隔离液相模型的影响。")
    print()

    # 对比3: 传统 vs UEM (液相)
    print("【对比3】传统 Muggianu vs UEM (只液相)")
    print("-" * 90)
    diff_mug_uem1 = []
    for i in range(len(x_ni_range)):
        t1 = results['传统Model(None)'][i]
        t4 = results['ModelUEM1(液相)'][i]
        if t1 and t4:
            diff = t1 - t4
            diff_mug_uem1.append(diff)
            print(f"  X(NI)={x_ni_range[i]:.3f}: Muggianu={t1:.1f} K, UEM1={t4:.1f} K, 差异={diff:.2f} K")

    if len(diff_mug_uem1) > 0:
        diff_mug_uem1 = np.array(diff_mug_uem1)
        print(f"\n统计:")
        print(f"  平均差异: {np.mean(diff_mug_uem1):.2f} K")
        print(f"  最大差异: {np.max(diff_mug_uem1):.2f} K")
        print(f"  最小差异: {np.min(diff_mug_uem1):.2f} K")
    print()

    print("=" * 90)
    print("测试完成！")
    print("=" * 90)


if __name__ == "__main__":
    main()
