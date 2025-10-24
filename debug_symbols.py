#!/usr/bin/env python3
"""
Debug what symbols are in the expressions
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
print("Symbol Debugging")
print("=" * 70)

comps = ['AL', 'CR', 'NI', 'VA']
phase = 'LIQUID'

# Create models
std_model = Model(dbf, comps, phase)
uem_model = ModelUEM(dbf, comps, phase)

print("\n" + "-" * 70)
print("Standard Model GM Symbols")
print("-" * 70)

std_gm = std_model.GM
print(f"GM type: {type(std_gm)}")
print(f"GM expression (first 500 chars):\n{str(std_gm)[:500]}...")

# Get free symbols
std_symbols = std_gm.free_symbols
print(f"\nFree symbols in standard GM:")
for sym in sorted(std_symbols, key=str):
    print(f"  {sym}")

print("\n" + "-" * 70)
print("UEM Model GM Symbols")
print("-" * 70)

uem_gm = uem_model.GM
print(f"GM type: {type(uem_gm)}")
print(f"GM expression (first 500 chars):\n{str(uem_gm)[:500]}...")

# Get free symbols
uem_symbols = uem_gm.free_symbols
print(f"\nFree symbols in UEM GM:")
for sym in sorted(uem_symbols, key=str):
    print(f"  {sym}")

print("\n" + "=" * 70)
