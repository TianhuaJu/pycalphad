#!/usr/bin/env python3
"""
Test site ratio normalization
"""
import sys
import os

# Force local pycalphad
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from pycalphad import Database, Model
from pycalphad.models.model_uem import ModelUEM

# Load database
dbf = Database('pycalphad/tests/databases/alcrni.tdb')

comps = ['AL', 'CR', 'NI', 'VA']
phase = 'LIQUID'

print("=" * 70)
print("Site Ratio Normalization Test")
print("=" * 70)

# Create models
std_model = Model(dbf, comps, phase)
uem_model = ModelUEM(dbf, comps, phase)

print(f"\nStandard model site ratio normalization:")
print(f"  {std_model._site_ratio_normalization}")
print(f"  Type: {type(std_model._site_ratio_normalization)}")

print(f"\nUEM model site ratio normalization:")
print(f"  {uem_model._site_ratio_normalization}")
print(f"  Type: {type(uem_model._site_ratio_normalization)}")

print(f"\nAre they equal? {std_model._site_ratio_normalization == uem_model._site_ratio_normalization}")

# Try getting excess before and after normalization
import pycalphad.models.uem_symbolic as uem
from pycalphad import variables as v

comps_no_va = ['AL', 'CR', 'NI']

print("\n" + "-" * 70)
print("Getting UEM expression")
print("-" * 70)

try:
    raw_expr = uem.get_uem1_excess_gibbs_expr(dbf, comps_no_va, phase, v.T)
    print(f"Raw UEM expression length: {len(str(raw_expr))} chars")
    print(f"Raw UEM expression type: {type(raw_expr)}")

    norm_expr = raw_expr / uem_model._site_ratio_normalization
    print(f"\nNormalized expression length: {len(str(norm_expr))} chars")
    print(f"Normalized expression type: {type(norm_expr)}")

    # Try substituting values
    from pycalphad import variables as v
    subs_dict = {
        v.T: 1800,
        v.SiteFraction('LIQUID', 0, 'AL'): 0.33,
        v.SiteFraction('LIQUID', 0, 'CR'): 0.33,
        v.SiteFraction('LIQUID', 0, 'NI'): 0.34
    }

    raw_val = raw_expr.subs(subs_dict)
    print(f"\nRaw UEM value: {float(raw_val):.2f} J/mol")

    norm_val = norm_expr.subs(subs_dict)
    print(f"Normalized UEM value: {float(norm_val):.2f} J/mol")

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 70)
