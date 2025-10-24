"""诊断三元系统UEM"""
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from pycalphad import Database, Model, variables as v
from pycalphad.models.model_uem import ModelUEM
import pycalphad.models.uem_symbolic as uem

print("=" * 70)
print("三元系统UEM诊断")
print("=" * 70)

dbf = Database('pycalphad/tests/databases/alcrni.tdb')

print("\n【测试1】三元系统 Al-Cr-Ni")
print("-" * 70)

comps = ['AL', 'CR', 'NI', 'VA']
phase_name = 'LIQUID'

try:
    mod_uem = ModelUEM(dbf, comps, phase_name)
    print("✓ UEM模型创建成功")

    # 检查组分数量
    comps_no_va = [str(c) for c in mod_uem.components if str(c) != 'VA']
    print(f"不含VA的组分: {comps_no_va}")
    print(f"组分数量: {len(comps_no_va)}")

    if len(comps_no_va) == 2:
        print("⚠️ 会触发二元特殊处理")
    else:
        print("✓ 会使用UEM三元算法")

    # 获取excess表达式
    print("\n调用get_uem1_excess_gibbs_expr...")
    xsmix = mod_uem.excess_mixing_energy(dbf)
    print(f"excess_mixing_energy类型: {type(xsmix)}")
    print(f"表达式长度: {len(str(xsmix))} 字符")

    if 'sympy' in str(type(xsmix)):
        print("✗ 警告：仍是sympy对象!")
    else:
        print("✓ 是symengine对象")

    # 获取GM
    gm_expr = mod_uem.GM
    print(f"\nGM表达式类型: {type(gm_expr)}")

    # 数值代入
    print("\n数值代入测试...")
    subs_dict = {
        v.T: 1800,
        v.P: 101325,
        v.SiteFraction(phase_name, 0, 'AL'): 0.33,
        v.SiteFraction(phase_name, 0, 'CR'): 0.33,
        v.SiteFraction(phase_name, 0, 'NI'): 0.34,
        v.N: 1
    }

    gm_val = gm_expr.subs(subs_dict)
    print(f"GM (代入后) 类型: {type(gm_val)}")

    # 检查是否还有符号
    gm_str = str(gm_val)
    if 'Y(LIQUID' in gm_str:
        print("✗ 代入后仍包含Y符号!")
        print(f"前200字符: {gm_str[:200]}")
    else:
        print("✓ 所有符号已替换")

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

print("\n【测试2】纯组分AL（三元系统中）")
print("-" * 70)

try:
    mod_uem = ModelUEM(dbf, comps, phase_name)

    subs_dict = {
        v.T: 1800,
        v.P: 101325,
        v.SiteFraction(phase_name, 0, 'AL'): 1.0,
        v.SiteFraction(phase_name, 0, 'CR'): 0.0,
        v.SiteFraction(phase_name, 0, 'NI'): 0.0,
        v.N: 1
    }

    gm_val = mod_uem.GM.subs(subs_dict)
    print(f"GM (代入后) 类型: {type(gm_val)}")

    gm_str = str(gm_val)
    if 'Y(LIQUID' in gm_str:
        print(f"✗ 仍包含Y符号: {gm_str[:200]}")

    try:
        gm_float = float(gm_val)
        print(f"GM (float): {gm_float:.2f} J/mol")
    except Exception as e:
        print(f"✗ 无法转换为float: {e}")

except Exception as e:
    print(f"✗ 错误: {e}")
    import traceback
    traceback.print_exc()

print("\n【测试3】直接测试_binary_excess")
print("-" * 70)

try:
    # 测试新的_binary_excess函数
    x_al = v.SiteFraction(phase_name, 0, 'AL')
    x_ni = v.SiteFraction(phase_name, 0, 'NI')

    print("调用_binary_excess('AL', 'NI')...")
    binary_expr = uem._binary_excess(dbf, 'AL', 'NI', phase_name, x_al, x_ni)

    print(f"binary_excess类型: {type(binary_expr)}")
    print(f"表达式长度: {len(str(binary_expr))} 字符")

    if 'sympy' in str(type(binary_expr)):
        print("✗ 返回sympy对象!")
    else:
        print("✓ 返回symengine对象")

    # 数值代入
    subs = {x_al: 0.5, x_ni: 0.5, v.T: 1800}
    result = binary_expr.subs(subs)

    try:
        result_float = float(result)
        print(f"数值结果: {result_float:.2f} J/mol")
    except Exception as e:
        print(f"✗ 无法转换: {e}")

except Exception as e:
    print(f"✗ 错误: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 70)
print("诊断完成")
print("=" * 70)
