"""
正确的UEM测试脚本示例 - 强制使用本地pycalphad

这个脚本展示如何正确使用calculate()函数测试UEM模型。
"""
import sys
import os

# 强制使用本地pycalphad（在最前面！）
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

print(f"Python路径优先级:")
print(f"  1. {sys.path[0]}")
print(f"  2. {sys.path[1] if len(sys.path) > 1 else 'N/A'}")

import numpy as np
from pycalphad import Database, calculate, variables as v
from pycalphad.models.model_uem import ModelUEM

# 确认使用的是本地版本
import pycalphad
print(f"\n正在使用的pycalphad位置:")
print(f"  {pycalphad.__file__}")

if 'site-packages' in pycalphad.__file__:
    print("\n⚠️ 警告: 正在使用site-packages中安装的版本!")
    print("请运行: pip install -e . 以使用本地开发版本")
    sys.exit(1)
else:
    print("✓ 正在使用本地开发版本")

# 加载数据库
dbf = Database('pycalphad/tests/databases/alcrni.tdb')

print("\n" + "=" * 60)
print("UEM模型测试脚本")
print("=" * 60)

# ============================================================================
# 测试1: 二元系统（UEM应等于标准模型）
# ============================================================================
print("\n【测试1】二元系统 Al-Ni (UEM应与标准模型相同)")
print("-" * 60)

comps = ['AL', 'NI', 'VA']
phases = ['LIQUID']

# 正确用法：使用条件字典，键为状态变量对象
conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.5, v.N: 1}

# 计算标准模型（Muggianu）
res_std = calculate(dbf, comps, phases, conds)

# 计算UEM模型
res_uem = calculate(dbf, comps, phases, conds, model=ModelUEM)

gm_std = float(res_std.GM.values.flatten()[0])
gm_uem = float(res_uem.GM.values.flatten()[0])
diff = abs(gm_uem - gm_std)

print(f"标准模型 GM: {gm_std:.2f} J/mol")
print(f"UEM模型  GM: {gm_uem:.2f} J/mol")
print(f"差异:        {diff:.6f} J/mol")
print(f"✓ 二元系统UEM = 标准模型（差异 < 1e-6）" if diff < 1e-6 else "✗ 失败")

# ============================================================================
# 测试2: 三元系统（UEM应不同于Muggianu）
# ============================================================================
print("\n【测试2】三元系统 Al-Cr-Ni (UEM应与Muggianu不同)")
print("-" * 60)

comps = ['AL', 'CR', 'NI', 'VA']
phases = ['LIQUID']

conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.33, v.X('CR'): 0.33, v.N: 1}

res_std = calculate(dbf, comps, phases, conds)
res_uem = calculate(dbf, comps, phases, conds, model=ModelUEM)

gm_std = float(res_std.GM.values.flatten()[0])
gm_uem = float(res_uem.GM.values.flatten()[0])
abs_diff = abs(gm_uem - gm_std)
rel_diff = abs_diff / abs(gm_std) * 100

print(f"Muggianu GM: {gm_std:.2f} J/mol")
print(f"UEM GM:      {gm_uem:.2f} J/mol")
print(f"绝对差异:    {abs_diff:.2f} J/mol")
print(f"相对差异:    {rel_diff:.2f}%")
print(f"✓ 三元系统UEM ≠ Muggianu（差异1-15%）" if abs_diff > 100 else "✗ 差异过小")

# ============================================================================
# 测试3: 多个组分测试
# ============================================================================
print("\n【测试3】不同Al组分的二元系统测试")
print("-" * 60)

for x_al in [0.2, 0.5, 0.8]:
    conds = {v.T: 1800, v.P: 101325, v.X('AL'): x_al, v.N: 1}
    res_std = calculate(dbf, ['AL', 'NI', 'VA'], ['LIQUID'], conds)
    res_uem = calculate(dbf, ['AL', 'NI', 'VA'], ['LIQUID'], conds, model=ModelUEM)

    gm_std = float(res_std.GM.values.flatten()[0])
    gm_uem = float(res_uem.GM.values.flatten()[0])
    diff = abs(gm_uem - gm_std)
    status = "✓" if diff < 1e-6 else "✗"
    print(f"X(AL)={x_al:.1f}: 差异 = {diff:.8f} J/mol {status}")

# ============================================================================
# 测试4: 纯组分极限
# ============================================================================
print("\n【测试4】纯组分极限稳定性")
print("-" * 60)

comps = ['AL', 'CR', 'NI', 'VA']
phases = ['LIQUID']

# 纯AL
conds = {v.T: 1800, v.P: 101325, v.X('AL'): 1.0, v.X('CR'): 0.0, v.N: 1}
res = calculate(dbf, comps, phases, conds, model=ModelUEM)
gm_al = float(res.GM.values.flatten()[0])
print(f"纯AL: GM = {gm_al:.2f} J/mol {'✓ 有限' if np.isfinite(gm_al) else '✗ 无穷'}")

# 纯CR
conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.0, v.X('CR'): 1.0, v.N: 1}
res = calculate(dbf, comps, phases, conds, model=ModelUEM)
gm_cr = float(res.GM.values.flatten()[0])
print(f"纯CR: GM = {gm_cr:.2f} J/mol {'✓ 有限' if np.isfinite(gm_cr) else '✗ 无穷'}")

# 纯NI
conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.0, v.X('CR'): 0.0, v.N: 1}
res = calculate(dbf, comps, phases, conds, model=ModelUEM)
gm_ni = float(res.GM.values.flatten()[0])
print(f"纯NI: GM = {gm_ni:.2f} J/mol {'✓ 有限' if np.isfinite(gm_ni) else '✗ 无穷'}")

print("\n" + "=" * 60)
print("测试完成！")
print("=" * 60)

# ============================================================================
# 关键点说明
# ============================================================================
print("\n【关键语法说明】")
print("-" * 60)
print("✗ 错误用法：calculate(dbf, comps, phases, T=1800, X_AL=0.5)")
print("✓ 正确用法：calculate(dbf, comps, phases, {v.T: 1800, v.X('AL'): 0.5})")
print("\n组分摩尔分数必须使用 v.X('元素符号') 而不是 X_元素=值")
print("=" * 60)
