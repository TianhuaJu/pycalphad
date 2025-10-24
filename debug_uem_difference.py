#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
调试 UEMModel vs ModelUEM1 的差异
打印出计算的每个步骤，找出差异点
"""

import numpy as np
from pycalphad import Database, variables as v
from pycalphad.model import Model
from symengine import symbols

# 导入两个实现
from pycalphad.UEMModel import UEMModel
from pycalphad.model_uem_integrated import ModelUEM1

def test_binary_limit():
    """
    测试二元边界（应该回归到原始二元参数）
    """
    print("=" * 80)
    print("测试1：二元边界（Al-Cr 二元，xNi=0）")
    print("=" * 80)
    print()

    dbe = Database('examples/alcrni.tdb')
    comps = ['AL', 'CR', 'NI']
    phase_name = 'LIQUID'

    # 创建两个模型实例
    print("创建 UEMModel 实例...")
    model_uem = UEMModel(dbe, comps + ['VA'], phase_name)

    print("创建 ModelUEM1 实例...")
    model_uem1 = ModelUEM1(dbe, comps + ['VA'], phase_name)

    print()
    print("计算二元边界（xAl=0.5, xCr=0.5, xNi=0）的过剩能...")
    print()

    # 获取过剩能表达式（这是符号表达式）
    GE_uem = model_uem.excess_mixing_energy(dbe)
    GE_uem1 = model_uem1.excess_mixing_energy(dbe)

    print("UEMModel 过剩能表达式（前100个字符）:")
    print(str(GE_uem)[:100] + "...")
    print()

    print("ModelUEM1 过剩能表达式（前100个字符）:")
    print(str(GE_uem1)[:100] + "...")
    print()

    # 代入数值
    subs_dict = {
        v.SiteFraction('LIQUID', 0, v.Species('AL')): 0.5,
        v.SiteFraction('LIQUID', 0, v.Species('CR')): 0.5,
        v.SiteFraction('LIQUID', 0, v.Species('NI')): 0.0,
        v.T: 1800,
        v.R: 8.314
    }

    GE_uem_val = float(GE_uem.subs(subs_dict))
    GE_uem1_val = float(GE_uem1.subs(subs_dict))

    print(f"UEMModel 过剩能（数值）:  {GE_uem_val:.2f} J/mol")
    print(f"ModelUEM1 过剩能（数值）: {GE_uem1_val:.2f} J/mol")
    print(f"差异:                     {abs(GE_uem_val - GE_uem1_val):.2f} J/mol")
    print()

    # 也用传统 Model 测试
    model_std = Model(dbe, comps + ['VA'], phase_name)
    GE_std = model_std.excess_mixing_energy(dbe)
    GE_std_val = float(GE_std.subs(subs_dict))

    print(f"传统 Model 过剩能（数值）:{GE_std_val:.2f} J/mol")
    print()

    if abs(GE_uem_val - GE_std_val) < 1.0 and abs(GE_uem1_val - GE_std_val) < 1.0:
        print("✅ 二元边界测试通过：所有模型在二元边界回归到相同值")
    else:
        print("❌ 二元边界测试失败：UEM 在二元边界应该等于原始参数")

    print("=" * 80)
    print()


def test_ternary_point():
    """
    测试三元点
    """
    print("=" * 80)
    print("测试2：三元点（Al=0.35, Cr=0.35, Ni=0.3）")
    print("=" * 80)
    print()

    dbe = Database('examples/alcrni.tdb')
    comps = ['AL', 'CR', 'NI']
    phase_name = 'LIQUID'

    model_uem = UEMModel(dbe, comps + ['VA'], phase_name)
    model_uem1 = ModelUEM1(dbe, comps + ['VA'], phase_name)

    GE_uem = model_uem.excess_mixing_energy(dbe)
    GE_uem1 = model_uem1.excess_mixing_energy(dbe)

    # 代入数值
    subs_dict = {
        v.SiteFraction('LIQUID', 0, v.Species('AL')): 0.35,
        v.SiteFraction('LIQUID', 0, v.Species('CR')): 0.35,
        v.SiteFraction('LIQUID', 0, v.Species('NI')): 0.30,
        v.T: 1800,
        v.R: 8.314
    }

    GE_uem_val = float(GE_uem.subs(subs_dict))
    GE_uem1_val = float(GE_uem1.subs(subs_dict))

    print(f"UEMModel 过剩能:   {GE_uem_val:.2f} J/mol")
    print(f"ModelUEM1 过剩能:  {GE_uem1_val:.2f} J/mol")
    print(f"差异:              {abs(GE_uem_val - GE_uem1_val):.2f} J/mol")
    print()

    # 传统 Model
    model_std = Model(dbe, comps + ['VA'], phase_name)
    GE_std = model_std.excess_mixing_energy(dbe)
    GE_std_val = float(GE_std.subs(subs_dict))

    print(f"传统 Model 过剩能: {GE_std_val:.2f} J/mol")
    print()

    print(f"UEMModel vs 传统:  {GE_uem_val - GE_std_val:+.2f} J/mol")
    print(f"ModelUEM1 vs 传统: {GE_uem1_val - GE_std_val:+.2f} J/mol")

    print()
    print("=" * 80)
    print()


if __name__ == "__main__":
    test_binary_limit()
    test_ternary_point()

    print()
    print("=" * 80)
    print("总结")
    print("=" * 80)
    print()
    print("如果二元边界测试通过，说明两个实现在二元情况下是正确的。")
    print("如果三元点差异很大，说明三元外推的算法实现有差异。")
    print()
    print("关键检查点：")
    print("1. d_ki 的计算是否相同")
    print("2. alpha 的计算是否相同")
    print("3. X_ij 的计算是否相同")
    print("4. R-K 多项式的构建是否相同")
    print()
