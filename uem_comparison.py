#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
对比测试：比较两个不同的 UEM 实现
========================================

检查 pycalphad.uem1_model.uem1_model 和 pycalphad.model_uem_integrated.ModelUEM1
是否产生相同的结果。
"""

import numpy as np
from pycalphad import Database, equilibrium, variables as v
from pycalphad.model import Model

# 导入两个不同的 UEM 实现
from pycalphad.uem1_Model import uem1_model
from pycalphad.advanced_uem_model import ModelUEM1

def calcu_single_point():
    """测试单个成分点的液相线计算"""

    print("=" * 70)
    print("UEM 实现对比测试")
    print("=" * 70)
    print()

    # 加载数据库
    dbe = Database('examples/alcrni.tdb')
    comps = ['AL', 'CR', 'NI']
    phases = ['LIQUID', 'FCC_A1', 'BCC_A2']

    # 测试成分
    x_ni = 0.3
    x_al = 0.35
    x_cr = 0.35

    print(f"测试成分: X(NI)={x_ni}, X(AL)={x_al}, X(CR)={x_cr}")
    print()

    conditions = {
        v.T: (1400, 2400, 10),  # 扩展温度范围到 2400K
        v.P: 101325,
        v.X('AL'): x_al,
        v.X('CR'): x_cr,
    }

    # 方案 1: 使用 uem1_model（所有相）
    print("方案 1: 使用 pycalphad.uem1_model.uem1_model（所有相）")
    print("-" * 70)
    models_uem = {}
    for phase_name in phases:
        try:
            models_uem[phase_name] = uem1_model(dbe, comps + ['VA'], phase_name)
        except Exception as e:
            print(f"  警告: 构建 {phase_name} 的 uem1_model 失败: {e}")

    print(f"  成功构建: {list(models_uem.keys())}")

    eq1 = equilibrium(dbe, comps, phases, model=models_uem, conditions=conditions)

    # 提取液相线
    liquidus_1 = None
    for t_idx in range(len(eq1.T.values) - 1, -1, -1):
        liquid_np = None
        for v_idx in range(eq1.Phase.shape[-1]):
            phase_name = eq1.Phase.values[0, 0, t_idx, 0, 0, v_idx]
            if phase_name == 'LIQUID':
                liquid_np = eq1.NP.values[0, 0, t_idx, 0, 0, v_idx]
                break
        if liquid_np is not None and liquid_np >= 0.995:
            liquidus_1 = eq1.T.values[t_idx]
        elif liquid_np is not None and liquid_np < 0.995:
            break

    print(f"  液相线温度: {liquidus_1:.1f} K" if liquidus_1 else "  液相线温度: 未找到")
    print()

    # 方案 2: 使用 ModelUEM1（所有相）
    print("方案 2: 使用 pycalphad.model_uem_integrated.ModelUEM1（所有相）")
    print("-" * 70)
    models_uem1 = {}
    for phase_name in phases:
        try:
            models_uem1[phase_name] = ModelUEM1(dbe, comps + ['VA'], phase_name)
        except Exception as e:
            print(f"  警告: 构建 {phase_name} 的 ModelUEM1 失败: {e}")

    print(f"  成功构建: {list(models_uem1.keys())}")

    eq2 = equilibrium(dbe, comps, phases, model=models_uem1, conditions=conditions)

    # 提取液相线
    liquidus_2 = None
    for t_idx in range(len(eq2.T.values) - 1, -1, -1):
        liquid_np = None
        for v_idx in range(eq2.Phase.shape[-1]):
            phase_name = eq2.Phase.values[0, 0, t_idx, 0, 0, v_idx]
            if phase_name == 'LIQUID':
                liquid_np = eq2.NP.values[0, 0, t_idx, 0, 0, v_idx]
                break
        if liquid_np is not None and liquid_np >= 0.995:
            liquidus_2 = eq2.T.values[t_idx]
        elif liquid_np is not None and liquid_np < 0.995:
            break

    print(f"  液相线温度: {liquidus_2:.1f} K" if liquidus_2 else "  液相线温度: 未找到")
    print()

    # 方案 3: 使用传统 Model（作为基准）
    print("方案 3: 使用传统 Model（基准）")
    print("-" * 70)

    eq3 = equilibrium(dbe, comps, phases, model=None, conditions=conditions)

    # 提取液相线
    liquidus_3 = None
    for t_idx in range(len(eq3.T.values) - 1, -1, -1):
        liquid_np = None
        for v_idx in range(eq3.Phase.shape[-1]):
            phase_name = eq3.Phase.values[0, 0, t_idx, 0, 0, v_idx]
            if phase_name == 'LIQUID':
                liquid_np = eq3.NP.values[0, 0, t_idx, 0, 0, v_idx]
                break
        if liquid_np is not None and liquid_np >= 0.995:
            liquidus_3 = eq3.T.values[t_idx]
        elif liquid_np is not None and liquid_np < 0.995:
            break

    print(f"  液相线温度: {liquidus_3:.1f} K" if liquidus_3 else "  液相线温度: 未找到")
    print()

    # 对比结果
    print("=" * 70)
    print("结果对比")
    print("=" * 70)
    if liquidus_1 and liquidus_2 and liquidus_3:
        print(f"uem1_model 液相线:        {liquidus_1:.1f} K")
        print(f"ModelUEM1 液相线:       {liquidus_2:.1f} K")
        print(f"传统 Model 液相线:      {liquidus_3:.1f} K")
        print()
        print(f"uem1_model - ModelUEM1:   {liquidus_1 - liquidus_2:+.1f} K")
        print(f"uem1_model - Model:       {liquidus_1 - liquidus_3:+.1f} K")
        print(f"ModelUEM1 - Model:      {liquidus_2 - liquidus_3:+.1f} K")
        print()

        if abs(liquidus_1 - liquidus_2) < 1.0:
            print("✅ 两个 UEM 实现产生相同结果（差异 < 1K）")
        else:
            print("❌ 两个 UEM 实现产生不同结果！")
            print(f"   差异: {abs(liquidus_1 - liquidus_2):.1f} K")
    else:
        print("❌ 至少有一个方案未找到液相线")

    print("=" * 70)


if __name__ == "__main__":
    calcu_single_point()
