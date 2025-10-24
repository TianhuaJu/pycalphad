#!/usr/bin/env python3
"""
Debug total energy calculation
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
print("Total Energy Calculation Debug")
print("=" * 70)

comps = ['AL', 'CR', 'NI', 'VA']
phase = 'LIQUID'

# Create models
std_model = Model(dbf, comps, phase)
uem_model = ModelUEM(dbf, comps, phase)

# Substitution dictionary
subs_dict = {
    v.T: 1800,
    v.P: 101325,
    v.SiteFraction('LIQUID', 0, 'AL'): 0.33,
    v.SiteFraction('LIQUID', 0, 'CR'): 0.33,
    v.SiteFraction('LIQUID', 0, 'NI'): 0.34
}

print("\n" + "-" * 70)
print("Standard Model Energy Components")
print("-" * 70)

# Standard model components
ref_std = std_model.reference_energy(dbf).subs(subs_dict)
idmix_std = std_model.ideal_mixing_energy(dbf).subs(subs_dict)
xsmix_std = std_model.excess_mixing_energy(dbf).subs(subs_dict)
total_std = std_model.GM.subs(subs_dict)

print(f"Reference energy:      {float(ref_std):12.2f} J/mol")
print(f"Ideal mixing energy:   {float(idmix_std):12.2f} J/mol")
print(f"Excess mixing energy:  {float(xsmix_std):12.2f} J/mol")
print(f"Total GM:              {float(total_std):12.2f} J/mol")

print("\n" + "-" * 70)
print("UEM Model Energy Components")
print("-" * 70)

# UEM model components
ref_uem = uem_model.reference_energy(dbf).subs(subs_dict)
idmix_uem = uem_model.ideal_mixing_energy(dbf).subs(subs_dict)
xsmix_uem = uem_model.excess_mixing_energy(dbf).subs(subs_dict)
total_uem = uem_model.GM.subs(subs_dict)

print(f"Reference energy:      {float(ref_uem):12.2f} J/mol")
print(f"Ideal mixing energy:   {float(idmix_uem):12.2f} J/mol")
print(f"Excess mixing energy:  {float(xsmix_uem):12.2f} J/mol")
print(f"Total GM:              {float(total_uem):12.2f} J/mol")

print("\n" + "-" * 70)
print("Differences")
print("-" * 70)

print(f"Δ Reference:           {float(ref_uem - ref_std):12.2f} J/mol")
print(f"Δ Ideal mixing:        {float(idmix_uem - idmix_std):12.2f} J/mol")
print(f"Δ Excess mixing:       {float(xsmix_uem - xsmix_std):12.2f} J/mol")
print(f"Δ Total GM:            {float(total_uem - total_std):12.2f} J/mol")

print("\n" + "-" * 70)
print("Site Ratio Normalization")
print("-" * 70)

print(f"Standard model: {std_model._site_ratio_normalization}")
print(f"UEM model: {uem_model._site_ratio_normalization}")

# Check raw excess before normalization
raw_xsmix_std = std_model.excess_mixing_energy(dbf)
raw_xsmix_uem = uem_model.excess_mixing_energy(dbf)

print(f"\nRaw excess (std) length: {len(str(raw_xsmix_std))} chars")
print(f"Raw excess (UEM) length: {len(str(raw_xsmix_uem))} chars")

print("\n" + "=" * 70)
