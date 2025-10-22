#!/usr/bin/env python3
"""
Test with logging enabled
"""
import sys
import os
import logging

# Force local pycalphad
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

# Enable logging
logging.basicConfig(level=logging.INFO, format='%(name)s - %(levelname)s - %(message)s')

from pycalphad import Database, calculate, variables as v
from pycalphad.models.model_uem import ModelUEM

print("=" * 70)
print("UEM Test with Logging")
print("=" * 70)

# Load database
dbf = Database('pycalphad/tests/databases/alcrni.tdb')

# Ternary test
comps = ['AL', 'CR', 'NI', 'VA']
phases = ['LIQUID']
conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.33, v.X('CR'): 0.33, v.N: 1}

print("\n" + "-" * 70)
print("Calculating with Standard Model")
print("-" * 70)
res_std = calculate(dbf, comps, phases, conds)
gm_std = float(res_std.GM.values.flatten()[0])
print(f"GM (Standard): {gm_std:.2f} J/mol")

print("\n" + "-" * 70)
print("Calculating with UEM Model")
print("-" * 70)
res_uem = calculate(dbf, comps, phases, conds, model=ModelUEM)
gm_uem = float(res_uem.GM.values.flatten()[0])
print(f"GM (UEM):      {gm_uem:.2f} J/mol")

print("\n" + "-" * 70)
print("Result")
print("-" * 70)
diff = abs(gm_uem - gm_std)
print(f"Difference: {diff:.2f} J/mol")

if diff > 100:
    print("✓ Models differ as expected")
elif diff < 1e-6:
    print("✗ Models identical (unexpected)")
else:
    print(f"? Small difference ({diff:.2f} J/mol)")

print("\n" + "=" * 70)
