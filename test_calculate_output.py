#!/usr/bin/env python3
"""
Test what calculate() returns
"""
import sys
sys.path.insert(0, '.')

from pycalphad import Database, calculate, variables as v
from pycalphad.models.model_uem import ModelUEM

dbf = Database('pycalphad/tests/databases/alcrni.tdb')
comps = ['AL', 'CR', 'NI', 'VA']
phases = ['LIQUID']
conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.33, v.X('CR'): 0.33, v.N: 1}

print("=" * 70)
print("Calculate() Output Test")
print("=" * 70)

# Calculate with standard model
res_std = calculate(dbf, comps, phases, conds)
print(f"\\nStandard result dataset:")
print(f"  Variables: {list(res_std.data_vars)}")
print(f"  Coords: {list(res_std.coords)}")
print(f"  GM shape: {res_std.GM.shape}")
print(f"  GM values: {res_std.GM.values}")
print(f"  GM value [0]: {float(res_std.GM.values.flatten()[0]):.6f} J/mol")

# Calculate with UEM model
res_uem = calculate(dbf, comps, phases, conds, model=ModelUEM)
print(f"\\nUEM result dataset:")
print(f"  Variables: {list(res_uem.data_vars)}")
print(f"  Coords: {list(res_uem.coords)}")
print(f"  GM shape: {res_uem.GM.shape}")
print(f"  GM values: {res_uem.GM.values}")
print(f"  GM value [0]: {float(res_uem.GM.values.flatten()[0]):.6f} J/mol")

# Check if there are other properties
print(f"\\n{'-' * 70}")
print("Checking all data variables")
print("-" * 70)
for var in res_std.data_vars:
    std_val = float(res_std[var].values.flatten()[0])
    uem_val = float(res_uem[var].values.flatten()[0])
    diff = abs(uem_val - std_val)
    if diff > 1e-6:
        print(f"{var:20s}: std={std_val:12.2f}, uem={uem_val:12.2f}, diff={diff:10.2f}")
    else:
        print(f"{var:20s}: {std_val:12.2f} (identical)")

print("\\n" + "=" * 70)
