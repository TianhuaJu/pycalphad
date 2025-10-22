#!/usr/bin/env python3
"""
Test GM evaluation
"""
import sys
sys.path.insert(0, '.')

from pycalphad import Database, Model, variables as v
from pycalphad.models.model_uem import ModelUEM

dbf = Database('pycalphad/tests/databases/alcrni.tdb')
comps = ['AL', 'CR', 'NI', 'VA']

# Create models
std_model = Model(dbf, comps, 'LIQUID')
uem_model = ModelUEM(dbf, comps, 'LIQUID')

# Substitution dictionary
subs_dict = {
    v.T: 1800,
    v.P: 101325,
    v.SiteFraction('LIQUID', 0, 'AL'): 0.33,
    v.SiteFraction('LIQUID', 0, 'CR'): 0.33,
    v.SiteFraction('LIQUID', 0, 'NI'): 0.34
}

# Evaluate GM
std_gm_val = std_model.GM.subs(subs_dict)
uem_gm_val = uem_model.GM.subs(subs_dict)

print("=" * 70)
print("GM Evaluation Test")
print("=" * 70)
print(f"\\nConditions:")
print(f"  T = {subs_dict[v.T]} K")
print(f"  P = {subs_dict[v.P]} Pa")
print(f"  X(AL) = {subs_dict[v.SiteFraction('LIQUID', 0, 'AL')]}")
print(f"  X(CR) = {subs_dict[v.SiteFraction('LIQUID', 0, 'CR')]}")
print(f"  X(NI) = {subs_dict[v.SiteFraction('LIQUID', 0, 'NI')]}")

print(f"\\n{'-' * 70}")
print("Results")
print("-" * 70)
print(f"Standard GM: {float(std_gm_val):.6f} J/mol")
print(f"UEM GM:      {float(uem_gm_val):.6f} J/mol")
print(f"Difference:  {float(uem_gm_val - std_gm_val):.6f} J/mol")

diff = abs(float(uem_gm_val - std_gm_val))
if diff > 100:
    print(f"\\n✓ Models differ by {diff:.2f} J/mol (expected)")
elif diff < 1e-6:
    print(f"\\n✗ Models identical (unexpected)")
else:
    print(f"\\n? Small difference: {diff:.2f} J/mol")

print("\\n" + "=" * 70)
