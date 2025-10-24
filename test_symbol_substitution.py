"""测试符号替换问题"""
import sys
import os

# 强制使用本地pycalphad
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from pycalphad import Database, Model, variables as v
from pycalphad.models.model_uem import ModelUEM

print("=" * 70)
print("符号替换测试")
print("=" * 70)

dbf = Database('pycalphad/tests/databases/alcrni.tdb')
comps = ['AL', 'NI', 'VA']
phase_name = 'LIQUID'

print("\n【测试1】检查_site_ratio_normalization")
print("-" * 70)

mod_std = Model(dbf, comps, phase_name)
mod_uem = ModelUEM(dbf, comps, phase_name)

print(f"标准Model._site_ratio_normalization: {mod_std._site_ratio_normalization}")
print(f"UEM Model._site_ratio_normalization: {mod_uem._site_ratio_normalization}")

print(f"\n标准Model._site_ratio_normalization类型: {type(mod_std._site_ratio_normalization)}")
print(f"UEM Model._site_ratio_normalization类型: {type(mod_uem._site_ratio_normalization)}")

print("\n【测试2】检查excess_mixing_energy（除以归一化之前）")
print("-" * 70)

xsmix_std = mod_std.redlich_kister_sum(
    dbf.phases[phase_name],
    dbf.search,
    dbf._parameters.search(
        (v.where('phase_name') == phase_name) &
        ((v.where('parameter_type') == 'G') | (v.where('parameter_type') == 'L')) &
        (v.where('constituent_array').test(mod_std._interaction_test))
    )
)

import pycalphad.models.uem_symbolic as uem
comps_no_va = [str(c) for c in mod_uem.components if str(c) != 'VA']
xsmix_uem_raw = uem.get_uem1_excess_gibbs_expr(dbf, comps_no_va, phase_name, v.T)

print(f"标准 excess (除以归一化之前) 类型: {type(xsmix_std)}")
print(f"UEM excess (除以归一化之前) 类型: {type(xsmix_uem_raw)}")

print(f"\n标准 excess 前200字符: {str(xsmix_std)[:200]}")
print(f"\nUEM excess 前200字符: {str(xsmix_uem_raw)[:200]}")

print("\n【测试3】检查符号类型")
print("-" * 70)

# 提取表达式中的符号
from sympy import preorder_traversal

def get_symbols(expr):
    symbols = set()
    for term in preorder_traversal(expr):
        if hasattr(term, 'is_Symbol') and term.is_Symbol:
            symbols.add(term)
    return symbols

std_symbols = get_symbols(xsmix_std)
uem_symbols = get_symbols(xsmix_uem_raw)

print(f"标准模型符号数量: {len(std_symbols)}")
print(f"UEM模型符号数量: {len(uem_symbols)}")

print(f"\n标准模型符号:")
for sym in sorted(std_symbols, key=str):
    print(f"  {sym} - 类型: {type(sym).__module__}.{type(sym).__name__}")

print(f"\nUEM模型符号:")
for sym in sorted(uem_symbols, key=str):
    print(f"  {sym} - 类型: {type(sym).__module__}.{type(sym).__name__}")

print("\n【测试4】直接测试符号替换")
print("-" * 70)

# 创建测试符号
test_symbol1 = v.SiteFraction(phase_name, 0, 'AL')
test_symbol2 = v.SiteFraction(phase_name, 0, 'AL')

print(f"test_symbol1: {test_symbol1}")
print(f"test_symbol2: {test_symbol2}")
print(f"是同一对象: {test_symbol1 is test_symbol2}")
print(f"是否相等: {test_symbol1 == test_symbol2}")
print(f"hash相等: {hash(test_symbol1) == hash(test_symbol2)}")

# 测试替换
from sympy import symbols
test_expr = test_symbol1 * test_symbol1
print(f"\ntest_expr: {test_expr}")
print(f"test_expr类型: {type(test_expr)}")

subs_result = test_expr.subs({test_symbol2: 0.5})
print(f"用test_symbol2替换后: {subs_result}")
print(f"结果类型: {type(subs_result)}")

try:
    result_float = float(subs_result)
    print(f"转换为float: {result_float}")
except Exception as e:
    print(f"无法转换为float: {e}")

print("\n" + "=" * 70)
