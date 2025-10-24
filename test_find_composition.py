#!/usr/bin/env python3
"""
Find the point that matches our composition
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

print("=" * 70)
print("Finding Composition Point")
print("=" * 70)

# Calculate with both models
res_std = calculate(dbf, comps, phases, conds)
res_uem = calculate(dbf, comps, phases, conds, model=ModelUEM)

# Get composition data
X_data_std = res_std.X.values
Y_data_std = res_std.Y.values

print(f"\\nDataset shape: {X_data_std.shape}")
print(f"Number of points: {X_data_std.shape[1]}")

# Find points close to X(AL)=0.33, X(CR)=0.33
# X is indexed by [N_point, component]
# Components order: AL, CR, NI (excluding VA)
AL_idx = [i for i, c in enumerate(res_std.coords['component'].values) if str(c) == 'AL'][0]
CR_idx = [i for i, c in enumerate(res_std.coords['component'].values) if str(c) == 'CR'][0]
NI_idx = [i for i, c in enumerate(res_std.coords['component'].values) if str(c) == 'NI'][0]

print(f"\\nComponent indices: AL={AL_idx}, CR={CR_idx}, NI={NI_idx}")

# Find points where X(AL) ≈ 0.33 and X(CR) ≈ 0.33
tolerance = 0.01
for i in range(min(100, X_data_std.shape[1])):  # Check first 100 points
    x_al = X_data_std[0, i, AL_idx] if X_data_std.ndim == 3 else X_data_std[i, AL_idx]
    x_cr = X_data_std[0, i, CR_idx] if X_data_std.ndim == 3 else X_data_std[i, CR_idx]
    x_ni = X_data_std[0, i, NI_idx] if X_data_std.ndim == 3 else X_data_std[i, NI_idx]

    if abs(x_al - 0.33) < tolerance and abs(x_cr - 0.33) < tolerance:
        gm_std = res_std.GM.values.flatten()[i]
        gm_uem = res_uem.GM.values.flatten()[i]

        print(f"\\nFound matching point at index {i}:")
        print(f"  X(AL) = {x_al:.6f}")
        print(f"  X(CR) = {x_cr:.6f}")
        print(f"  X(NI) = {x_ni:.6f}")
        print(f"  GM (std) = {gm_std:.6f} J/mol")
        print(f"  GM (uem) = {gm_uem:.6f} J/mol")
        print(f"  Difference = {gm_uem - gm_std:.6f} J/mol")
        break
else:
    print("\\nNo exact match found in first 100 points")
    print("Showing first few points:")
    for i in range(min(5, X_data_std.shape[1])):
        x_al = X_data_std[0, i, AL_idx] if X_data_std.ndim == 3 else X_data_std[i, AL_idx]
        x_cr = X_data_std[0, i, CR_idx] if X_data_std.ndim == 3 else X_data_std[i, CR_idx]
        gm_std = res_std.GM.values.flatten()[i]
        gm_uem = res_uem.GM.values.flatten()[i]
        print(f"  Point {i}: X(AL)={x_al:.3f}, X(CR)={x_cr:.3f}, GM_std={gm_std:.2f}, GM_uem={gm_uem:.2f}, diff={gm_uem-gm_std:.2f}")

print("\\n" + "=" * 70)
