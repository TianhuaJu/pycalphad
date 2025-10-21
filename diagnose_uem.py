"""诊断UEM模块加载问题"""
import sys
import os

print("=" * 70)
print("UEM模块诊断")
print("=" * 70)

# 1. 检查Python路径
print("\n【1】Python搜索路径:")
print("-" * 70)
for i, path in enumerate(sys.path, 1):
    print(f"{i}. {path}")

# 2. 检查模块是否已加载
print("\n【2】已加载的pycalphad模块:")
print("-" * 70)
pycalphad_modules = [key for key in sys.modules.keys() if 'pycalphad' in key]
for mod in pycalphad_modules:
    print(f"  - {mod}")

if pycalphad_modules:
    print(f"\n⚠️ 发现 {len(pycalphad_modules)} 个已加载的pycalphad模块!")
    print("这些模块可能使用了旧的缓存代码。")

# 3. 强制重新导入
print("\n【3】强制清除并重新导入:")
print("-" * 70)

# 清除所有pycalphad模块
modules_to_remove = list(pycalphad_modules)
for mod in modules_to_remove:
    del sys.modules[mod]
    print(f"  清除: {mod}")

# 重新导入
print("\n重新导入pycalphad...")
from pycalphad import Database, Model, variables as v
from pycalphad.models.model_uem import ModelUEM
import pycalphad.models.uem_symbolic as uem

print("✓ 导入成功")

# 4. 检查源代码
print("\n【4】检查UEM源代码:")
print("-" * 70)
uem_file = uem.__file__
print(f"uem_symbolic.py 路径: {uem_file}")

# 检查关键行
with open(uem_file, 'r', encoding='utf-8') as f:
    lines = f.readlines()

# 查找第118行（应该是使用SiteFraction的那一行）
for i in range(115, 120):
    if i < len(lines):
        print(f"{i+1:3d}: {lines[i].rstrip()}")

# 检查是否使用了SiteFraction
line_118 = lines[117] if len(lines) > 117 else ""
if 'SiteFraction' in line_118:
    print("\n✓ 第118行正确使用了SiteFraction")
elif 'v.X(' in line_118:
    print("\n✗ 第118行仍在使用v.X() - 代码未更新!")
else:
    print(f"\n⚠️ 第118行内容异常: {line_118.strip()}")

# 5. 测试UEM模型创建
print("\n【5】测试UEM模型:")
print("-" * 70)

dbf = Database('pycalphad/tests/databases/alcrni.tdb')
comps = ['AL', 'NI', 'VA']
phase_name = 'LIQUID'

try:
    mod_uem = ModelUEM(dbf, comps, phase_name)
    print("✓ UEM模型创建成功")

    # 获取表达式
    xsmix = mod_uem.excess_mixing_energy(dbf)
    xsmix_str = str(xsmix)

    # 检查变量类型
    if 'SiteFraction' in str(type(xsmix)) or 'Y(' in xsmix_str:
        print("✓ excess_mixing_energy使用了SiteFraction")
    elif 'X(' in xsmix_str and 'SiteFraction' not in xsmix_str:
        print("✗ excess_mixing_energy仍在使用X() - 缓存问题!")
        print(f"表达式: {xsmix_str[:200]}...")
    else:
        print(f"⚠️ 表达式内容: {xsmix_str[:200]}...")

    # 检查是否为零
    if xsmix == 0 or xsmix_str == '0':
        print("\n⚠️ 警告: excess_mixing_energy返回零!")
        print("可能原因:")
        print("  1. 二元参数未找到")
        print("  2. 表达式列表为空")
        print("  3. 数据库加载问题")

    # 尝试数值代入
    print("\n【6】数值代入测试:")
    print("-" * 70)
    subs_dict = {
        v.T: 1800,
        v.P: 101325,
        v.SiteFraction(phase_name, 0, 'AL'): 0.5,
        v.SiteFraction(phase_name, 0, 'NI'): 0.5,
        v.N: 1
    }

    gm_expr = mod_uem.GM
    gm_val = gm_expr.subs(subs_dict)

    print(f"GM (代入后): {gm_val}")

    try:
        gm_float = float(gm_val)
        print(f"GM (float): {gm_float:.2f} J/mol")

        import math
        if math.isnan(gm_float):
            print("✗ 结果是NaN!")
        else:
            print("✓ 结果是有效数值")
    except Exception as e:
        print(f"✗ 无法转换为float: {e}")

except Exception as e:
    print(f"✗ 错误: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 70)
print("诊断完成")
print("=" * 70)
