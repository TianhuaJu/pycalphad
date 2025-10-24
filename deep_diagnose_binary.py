"""深度诊断：检查UEM二元特殊处理"""
import sys
import os

# 强制使用本地pycalphad
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from pycalphad import Database, Model, variables as v
from pycalphad.models.model_uem import ModelUEM
import pycalphad.models.uem_symbolic as uem

print("=" * 70)
print("UEM二元系统深度诊断")
print("=" * 70)

# 加载数据库
dbf = Database('pycalphad/tests/databases/alcrni.tdb')

print("\n【步骤1】检查组分列表")
print("-" * 70)
comps_with_va = ['AL', 'NI', 'VA']
comps_no_va = [c for c in comps_with_va if c != 'VA']
print(f"包含VA的组分: {comps_with_va}")
print(f"不含VA的组分: {comps_no_va}")
print(f"不含VA的组分数量: {len(comps_no_va)}")

if len(comps_no_va) == 2:
    print("✓ 应该触发二元特殊处理")
else:
    print("✗ 不会触发二元特殊处理")

print("\n【步骤2】直接调用get_uem1_excess_gibbs_expr")
print("-" * 70)
phase_name = 'LIQUID'

try:
    expr = uem.get_uem1_excess_gibbs_expr(dbf, comps_no_va, phase_name, v.T)
    print(f"返回表达式类型: {type(expr)}")
    print(f"表达式长度: {len(str(expr))} 字符")
    print(f"表达式前500字符:\n{str(expr)[:500]}")

    # 检查是否包含复杂的权重
    expr_str = str(expr)
    if '(Y(LIQUID,0,AL) + Y(LIQUID,0,NI))**2' in expr_str:
        print("\n✗ 错误：表达式仍包含 (Y(AL) + Y(NI))**2 权重!")
        print("这意味着二元特殊处理没有生效!")
    else:
        print("\n✓ 表达式没有复杂的权重项")

    # 检查符号
    from sympy import preorder_traversal
    symbols_found = set()
    for term in preorder_traversal(expr):
        if hasattr(term, 'name'):
            symbols_found.add(str(term))

    print(f"\n表达式中的符号:")
    for sym in sorted(symbols_found):
        print(f"  - {sym}")

except Exception as e:
    print(f"✗ 错误: {e}")
    import traceback
    traceback.print_exc()

print("\n【步骤3】创建ModelUEM并获取GM表达式")
print("-" * 70)

try:
    mod_uem = ModelUEM(dbf, comps_with_va, phase_name)
    print("✓ UEM模型创建成功")

    gm_expr = mod_uem.GM
    print(f"GM表达式类型: {type(gm_expr)}")
    print(f"GM表达式长度: {len(str(gm_expr))} 字符")

except Exception as e:
    print(f"✗ 错误: {e}")
    import traceback
    traceback.print_exc()

print("\n【步骤4】数值代入测试")
print("-" * 70)

try:
    # 使用SiteFraction代入
    subs_dict = {
        v.T: 1800,
        v.P: 101325,
        v.SiteFraction(phase_name, 0, 'AL'): 0.5,
        v.SiteFraction(phase_name, 0, 'NI'): 0.5,
        v.N: 1
    }

    print(f"代入字典:")
    for key, val in subs_dict.items():
        print(f"  {key} = {val}")

    gm_val = gm_expr.subs(subs_dict)
    print(f"\nGM (代入后) 类型: {type(gm_val)}")
    print(f"GM (代入后) 值: {gm_val}")

    # 尝试转换为浮点数
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
        print(f"gm_val内容: {gm_val}")

        # 检查是否还有未求值的符号
        if 'Y(LIQUID' in str(gm_val):
            print("\n⚠️ 代入后仍包含 Y(LIQUID,...) 符号!")
            print("这意味着符号没有被正确替换")

except Exception as e:
    print(f"✗ 错误: {e}")
    import traceback
    traceback.print_exc()

print("\n【步骤5】与标准模型对比")
print("-" * 70)

try:
    mod_std = Model(dbf, comps_with_va, phase_name)
    gm_std_expr = mod_std.GM

    subs_dict = {
        v.T: 1800,
        v.P: 101325,
        v.SiteFraction(phase_name, 0, 'AL'): 0.5,
        v.SiteFraction(phase_name, 0, 'NI'): 0.5,
        v.N: 1
    }

    gm_std_val = gm_std_expr.subs(subs_dict)
    gm_std_float = float(gm_std_val)

    print(f"标准模型 GM: {gm_std_float:.2f} J/mol")
    print("✓ 标准模型工作正常")

except Exception as e:
    print(f"✗ 标准模型错误: {e}")

print("\n" + "=" * 70)
print("诊断完成")
print("=" * 70)
