#!/usr/bin/env python3
"""
Search ALL points for our composition
"""
import sys
sys.path.insert(0, '.')
import numpy as np

from pycalphad import Database, calculate, variables as v
from pycalphad.models.model_uem import ModelUEM

dbf = Database('pycalphad/tests/databases/alcrni.tdb')
comps = ['AL', 'CR', 'NI', 'VA']
phases = ['LIQUID']
conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.33, v.X('CR'): 0.33, v.N: 1}

print("Calculating...")
res_std = calculate(dbf, comps, phases, conds)
res_uem = calculate(dbf, comps, phases, conds, model=ModelUEM)

X_data = res_std.X.values
GM_std = res_std.GM.values.flatten()
GM_uem = res_uem.GM.values.flatten()

AL_idx = [i for i, c in enumerate(res_std.coords['component'].values) if str(c) == 'AL'][0]
CR_idx = [i for i, c in enumerate(res_std.coords['component'].values) if str(c) == 'CR'][0]
NI_idx = [i for i, c in enumerate(res_std.coords['component'].values) if str(c) == 'NI'][0]

print(f"Searching {len(GM_std)} points...")

# Find ALL points where models differ significantly
differences = []
for i in range(len(GM_std)):
    if np.isnan(GM_std[i]) or np.isnan(GM_uem[i]):
        continue

    x_al = X_data[0, i, AL_idx]
    x_cr = X_data[0, i, CR_idx]
    x_ni = X_data[0, i, NI_idx]

    diff = GM_uem[i] - GM_std[i]

    if abs(diff) > 100:  # More than 100 J/mol difference
        differences.append((i, x_al, x_cr, x_ni, GM_std[i], GM_uem[i], diff))

print(f"\\nFound {len(differences)} points with |diff| > 100 J/mol")

if differences:
    # Sort by closeness to target composition
    differences.sort(key=lambda x: abs(x[1] - 0.33) + abs(x[2] - 0.33))

    print("\\nTop 10 closest to X(AL)=0.33, X(CR)=0.33:")
    print("-" * 90)
    print(f"{'Idx':>6s} {'X(AL)':>8s} {'X(CR)':>8s} {'X(NI)':>8s} {'GM_std':>12s} {'GM_uem':>12s} {'Diff':>10s}")
    print("-" * 90)
    for i, x_al, x_cr, x_ni, gm_s, gm_u, diff in differences[:10]:
        print(f"{i:6d} {x_al:8.4f} {x_cr:8.4f} {x_ni:8.4f} {gm_s:12.2f} {gm_u:12.2f} {diff:10.2f}")
else:
    print("\\nNo points found with significant differences!")
    print("Checking maximum difference:")
    diffs = np.abs(GM_uem - GM_std)
    valid_diffs = diffs[~np.isnan(diffs)]
    if len(valid_diffs) > 0:
        max_diff = np.max(valid_diffs)
        max_idx = np.argmax(diffs)
        print(f"Maximum difference: {max_diff:.6f} J/mol at index {max_idx}")
