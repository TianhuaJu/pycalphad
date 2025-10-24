#!/usr/bin/env python3
"""
Simple diagnostic test for ternary UEM
"""
import sys
import os

# Force local pycalphad
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from pycalphad import Database, calculate, variables as v
from pycalphad.models.model_uem import ModelUEM

# Load database
dbf = Database('pycalphad/tests/databases/alcrni.tdb')

print("=" * 70)
print("Ternary UEM Diagnostic Test")
print("=" * 70)

# Ternary test
comps = ['AL', 'CR', 'NI', 'VA']
phases = ['LIQUID']
conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.33, v.X('CR'): 0.33, v.N: 1}

print(f"\nComponents: {comps[:3]}")
print(f"Phase: {phases[0]}")
print(f"T = {conds[v.T]} K")
print(f"X(AL) = {conds[v.X('AL')]}")
print(f"X(CR) = {conds[v.X('CR')]}")
print(f"X(NI) = {1 - conds[v.X('AL')] - conds[v.X('CR')]:.2f} (implicit)")

print("\n" + "-" * 70)
print("Calculating with Muggianu (default)")
print("-" * 70)
res_std = calculate(dbf, comps, phases, conds)
gm_std = float(res_std.GM.values.flatten()[0])
print(f"GM (Muggianu): {gm_std:.6f} J/mol")

print("\n" + "-" * 70)
print("Calculating with UEM")
print("-" * 70)
res_uem = calculate(dbf, comps, phases, conds, model=ModelUEM)
gm_uem = float(res_uem.GM.values.flatten()[0])
print(f"GM (UEM):      {gm_uem:.6f} J/mol")

print("\n" + "-" * 70)
print("Comparison")
print("-" * 70)
abs_diff = abs(gm_uem - gm_std)
rel_diff = abs_diff / abs(gm_std) * 100 if gm_std != 0 else 0
print(f"Absolute difference: {abs_diff:.2f} J/mol")
print(f"Relative difference: {rel_diff:.4f}%")

if abs_diff > 100:
    print("✓ UEM differs from Muggianu (expected)")
elif abs_diff < 1e-6:
    print("⚠️ UEM identical to Muggianu (unexpected for ternary)")
    print("\nDEBUG: Checking model construction...")

    # Create model to check internal state
    from pycalphad import Model
    std_model = Model(dbf, comps, 'LIQUID')
    uem_model = ModelUEM(dbf, comps, 'LIQUID')

    # Get component list
    std_comps = [str(c) for c in std_model.components if str(c) != 'VA']
    uem_comps = [str(c) for c in uem_model.components if str(c) != 'VA']

    print(f"Standard model components: {std_comps}")
    print(f"UEM model components: {uem_comps}")
    print(f"Number of components: {len(uem_comps)}")

    if len(uem_comps) == 2:
        print("⚠️ WARNING: Model sees only 2 components (binary case)")
    elif len(uem_comps) == 3:
        print("✓ Model correctly sees 3 components (ternary case)")
else:
    print(f"? Small difference detected ({abs_diff:.2f} J/mol)")

print("\n" + "=" * 70)
