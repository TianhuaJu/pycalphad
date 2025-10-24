#!/usr/bin/env python3
"""
Debug UEM expression building
"""
import sys
import os

# Force local pycalphad
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from pycalphad import Database, Model, variables as v
from pycalphad.models.model_uem import ModelUEM

# Load database
dbf = Database('pycalphad/tests/databases/alcrni.tdb')

print("=" * 70)
print("UEM Expression Building Debug")
print("=" * 70)

comps = ['AL', 'CR', 'NI', 'VA']
phase = 'LIQUID'

print(f"\nComponents: {comps[:3]}")
print(f"Phase: {phase}")

# Create models
print("\n" + "-" * 70)
print("Building Models")
print("-" * 70)

std_model = Model(dbf, comps, phase)
uem_model = ModelUEM(dbf, comps, phase)

print(f"Standard model type: {type(std_model).__name__}")
print(f"UEM model type: {type(uem_model).__name__}")

# Get excess mixing energy expressions
print("\n" + "-" * 70)
print("Excess Mixing Energy Expressions")
print("-" * 70)

std_excess = std_model.excess_mixing_energy(dbf)
uem_excess = uem_model.excess_mixing_energy(dbf)

print(f"\nStandard excess type: {type(std_excess)}")
print(f"UEM excess type: {type(uem_excess)}")

print(f"\nStandard excess (first 200 chars):\n{str(std_excess)[:200]}...")
print(f"\nUEM excess (first 200 chars):\n{str(uem_excess)[:200]}...")

# Check if they're identical
if str(std_excess) == str(uem_excess):
    print("\n⚠️ WARNING: Expressions are IDENTICAL!")
    print("This explains why results are the same.")
else:
    print("\n✓ Expressions are DIFFERENT (as expected)")
    print(f"Standard expression length: {len(str(std_excess))} chars")
    print(f"UEM expression length: {len(str(uem_excess))} chars")

# Substitute values to check numerical evaluation
print("\n" + "-" * 70)
print("Numerical Evaluation at T=1800, X(AL)=0.33, X(CR)=0.33")
print("-" * 70)

# Create substitution dictionary
subs_dict = {
    v.T: 1800,
    v.SiteFraction('LIQUID', 0, 'AL'): 0.33,
    v.SiteFraction('LIQUID', 0, 'CR'): 0.33,
    v.SiteFraction('LIQUID', 0, 'NI'): 0.34
}

try:
    std_val = float(std_excess.subs(subs_dict))
    uem_val = float(uem_excess.subs(subs_dict))

    print(f"Standard excess: {std_val:.2f} J/mol")
    print(f"UEM excess: {uem_val:.2f} J/mol")
    print(f"Difference: {abs(uem_val - std_val):.2f} J/mol")
except Exception as e:
    print(f"Error evaluating: {e}")

# Check contribution ratios
print("\n" + "-" * 70)
print("UEM Contribution Ratios")
print("-" * 70)

import pycalphad.models.uem_symbolic as uem

T = v.T
for k in ['AL', 'CR', 'NI']:
    for i in ['AL', 'CR', 'NI']:
        for j in ['AL', 'CR', 'NI']:
            if k != i and k != j and i != j:
                r_ki = uem.uem1_contribution_ratio(dbf, k, i, j, phase, T)
                r_ki_val = float(r_ki.subs({T: 1800})) if hasattr(r_ki, 'subs') else float(r_ki)
                print(f"r_{k}{i} (in {i}-{j}): {r_ki_val:.6f}")

print("\n" + "=" * 70)
