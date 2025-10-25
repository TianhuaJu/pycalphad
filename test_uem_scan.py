#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
快速扫描 UEM 贡献系数，查找异常点
"""

import numpy as np
from pycalphad import Database
from evaluate_uem_numerical import evaluate_uem_coefficients

# 加载数据库
dbe = Database('examples/alcrni.tdb')
phase_name = 'LIQUID'

print("=" * 90)
print("UEM 贡献系数扫描 - 查找异常点")
print("=" * 90)
print()

# 扫描成分范围
x_ni_vals = np.linspace(0.1, 0.9, 9)
temperature = 1700

print(f"温度: {temperature} K")
print(f"扫描范围: X(NI) = 0.1 ~ 0.9，固定 X(AL) = X(CR)")
print("=" * 90)
print()

results = []

for x_ni in x_ni_vals:
    x_al = (1.0 - x_ni) / 2.0
    x_cr = (1.0 - x_ni) / 2.0

    composition = {'AL': x_al, 'CR': x_cr, 'NI': x_ni}

    print(f"{'─' * 90}")
    print(f"成分: X(NI)={x_ni:.2f}, X(AL)={x_al:.2f}, X(CR)={x_cr:.2f}")
    print(f"{'─' * 90}")

    coeffs = evaluate_uem_coefficients(dbe, phase_name, composition, temperature, verbose=False)

    # 检查异常
    has_issue = False
    issues = []

    for c in coeffs:
        # 检查负的 alpha
        if c['alpha_i_k'] < 0:
            issues.append(f"  ⚠️  α_{c['binary_pair'][0]}^{c['third_comp']} = {c['alpha_i_k']:.6f} < 0")
            has_issue = True
        if c['alpha_j_k'] < 0:
            issues.append(f"  ⚠️  α_{c['binary_pair'][1]}^{c['third_comp']} = {c['alpha_j_k']:.6f} < 0")
            has_issue = True

        # 检查极大的 alpha
        if abs(c['alpha_i_k']) > 10:
            issues.append(f"  ⚠️  α_{c['binary_pair'][0]}^{c['third_comp']} = {c['alpha_i_k']:.6f} 极大！")
            has_issue = True
        if abs(c['alpha_j_k']) > 10:
            issues.append(f"  ⚠️  α_{c['binary_pair'][1]}^{c['third_comp']} = {c['alpha_j_k']:.6f} 极大！")
            has_issue = True

        # 检查 NaN
        if np.isnan(c['alpha_i_k']) or np.isnan(c['alpha_j_k']):
            issues.append(f"  ⚠️⚠️⚠️  {c['binary_pair']} - {c['third_comp']}: NaN!")
            has_issue = True

        # 检查等效组分超出范围
        if c['X_ij_i'] < -0.01 or c['X_ij_i'] > 1.01:
            issues.append(f"  ⚠️  X_{c['binary_pair'][0]}{c['binary_pair'][1]}^{c['binary_pair'][0]} = {c['X_ij_i']:.4f} 超出 [0,1]")
            has_issue = True
        if c['X_ij_j'] < -0.01 or c['X_ij_j'] > 1.01:
            issues.append(f"  ⚠️  X_{c['binary_pair'][0]}{c['binary_pair'][1]}^{c['binary_pair'][1]} = {c['X_ij_j']:.4f} 超出 [0,1]")
            has_issue = True

    if has_issue:
        print("  发现异常：")
        for issue in issues:
            print(issue)

        # 打印详细信息
        print("\n  详细数据：")
        for c in coeffs:
            print(f"    {c['binary_pair']} - {c['third_comp']}:")
            print(f"      d_ki = {c['d_ki']:.6f}, d_kj = {c['d_kj']:.6f}, d_sum = {c['d_sum']:.6f}")
            print(f"      α_i = {c['alpha_i_k']:.6f}, α_j = {c['alpha_j_k']:.6f}")
            print(f"      X_i = {c['X_ij_i']:.6f}, X_j = {c['X_ij_j']:.6f}")

        results.append({
            'composition': composition,
            'coeffs': coeffs,
            'issues': issues
        })
    else:
        print("  ✓ 正常")

    print()

# 汇总
print("\n" + "=" * 90)
print("汇总结果")
print("=" * 90)

if results:
    print(f"\n共发现 {len(results)} 个异常点：\n")
    for r in results:
        comp = r['composition']
        print(f"X(NI)={comp['NI']:.2f}, X(AL)={comp['AL']:.2f}, X(CR)={comp['CR']:.2f}")
        for issue in r['issues']:
            print(issue)
        print()
else:
    print("\n✓ 未发现异常点！所有成分点的贡献系数都在合理范围内。")

print("=" * 90)
print("扫描完成！")
print("=" * 90)
