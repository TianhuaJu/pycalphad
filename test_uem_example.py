"""
正确的UEM测试脚本示例

这个脚本展示如何正确使用calculate()函数测试UEM模型。
"""

from pycalphad import Database, calculate, variables as v
from pycalphad.models.model_uem import ModelUEM

# 加载数据库
dbf = Database('pycalphad/tests/databases/alcrni.tdb')

print("=" * 60)
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

print(f"标准模型 GM: {res_std.GM.values[0]:.2f} J/mol")
print(f"UEM模型  GM: {res_uem.GM.values[0]:.2f} J/mol")
print(f"差异:        {abs(res_uem.GM.values[0] - res_std.GM.values[0]):.6f} J/mol")
print(f"✓ 二元系统UEM = 标准模型（差异 < 1e-6）" if abs(res_uem.GM.values[0] - res_std.GM.values[0]) < 1e-6 else "✗ 失败")

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

print(f"Muggianu GM: {res_std.GM.values[0]:.2f} J/mol")
print(f"UEM GM:      {res_uem.GM.values[0]:.2f} J/mol")
print(f"绝对差异:    {abs(res_uem.GM.values[0] - res_std.GM.values[0]):.2f} J/mol")
print(f"相对差异:    {abs(res_uem.GM.values[0] - res_std.GM.values[0]) / abs(res_std.GM.values[0]) * 100:.2f}%")
print(f"✓ 三元系统UEM ≠ Muggianu（差异1-15%）" if abs(res_uem.GM.values[0] - res_std.GM.values[0]) > 100 else "✗ 差异过小")

# ============================================================================
# 测试3: 多个组分测试
# ============================================================================
print("\n【测试3】不同Al组分的二元系统测试")
print("-" * 60)

for x_al in [0.2, 0.5, 0.8]:
    conds = {v.T: 1800, v.P: 101325, v.X('AL'): x_al, v.N: 1}
    res_std = calculate(dbf, ['AL', 'NI', 'VA'], ['LIQUID'], conds)
    res_uem = calculate(dbf, ['AL', 'NI', 'VA'], ['LIQUID'], conds, model=ModelUEM)

    diff = abs(res_uem.GM.values[0] - res_std.GM.values[0])
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
print(f"纯AL: GM = {res.GM.values[0]:.2f} J/mol {'✓ 有限' if np.isfinite(res.GM.values[0]) else '✗ 无穷'}")

# 纯CR
conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.0, v.X('CR'): 1.0, v.N: 1}
res = calculate(dbf, comps, phases, conds, model=ModelUEM)
print(f"纯CR: GM = {res.GM.values[0]:.2f} J/mol {'✓ 有限' if np.isfinite(res.GM.values[0]) else '✗ 无穷'}")

# 纯NI
conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.0, v.X('CR'): 0.0, v.N: 1}
res = calculate(dbf, comps, phases, conds, model=ModelUEM)
print(f"纯NI: GM = {res.GM.values[0]:.2f} J/mol {'✓ 有限' if np.isfinite(res.GM.values[0]) else '✗ 无穷'}")

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
