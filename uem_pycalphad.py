# UEM 集成测试脚本
import time
from pycalphad import Database, calculate, variables as v
from pycalphad.models.model_uem import ModelUEM, ModelMuggianu, ModelKohler

dbf = Database('examples/alcrni.tdb')
comps = ['AL', 'NI', 'CR', 'VA']
phases = ['LIQUID']

# UEM1 计算
t0 = time.perf_counter()
result_uem = calculate(dbf, comps, phases, model=ModelUEM, T=1800, P=101325)
t_uem = time.perf_counter() - t0
print(f"UEM1:     min={result_uem.GM.values.min():.1f}  max={result_uem.GM.values.max():.1f}  ({t_uem:.2f}s)")

# Muggianu (与原始 pycalphad 等价)
t0 = time.perf_counter()
result_mug = calculate(dbf, comps, phases, model=ModelMuggianu, T=1800, P=101325)
t_mug = time.perf_counter() - t0
print(f"Muggianu: min={result_mug.GM.values.min():.1f}  max={result_mug.GM.values.max():.1f}  ({t_mug:.2f}s)")

# 原始 R-K-M
t0 = time.perf_counter()
result_rkm = calculate(dbf, comps, phases, T=1800, P=101325)
t_rkm = time.perf_counter() - t0
print(f"R-K-M:    min={result_rkm.GM.values.min():.1f}  max={result_rkm.GM.values.max():.1f}  ({t_rkm:.2f}s)")
