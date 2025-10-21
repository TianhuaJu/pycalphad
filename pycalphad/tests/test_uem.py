"""
Tests for Unified Extrapolation Model (UEM).
"""
from pycalphad import Database, Model, calculate, variables as v
from pycalphad.models.model_uem import ModelUEM
from pycalphad.tests.fixtures import select_database, load_database
import numpy as np
import pytest


@select_database("alcrni.tdb")
def test_uem_binary_equivalence(load_database):
    """UEM should give identical results to standard Model for binary systems."""
    dbf = load_database()
    comps = ['AL', 'NI', 'VA']
    phases = ['LIQUID']

    # Test at several compositions
    for x_al in [0.2, 0.5, 0.8]:
        conds = {v.T: 1800, v.P: 101325, v.X('AL'): x_al, v.N: 1}
        res_std = calculate(dbf, comps, phases, conds)
        res_uem = calculate(dbf, comps, phases, conds, model=ModelUEM)

        # Should be identical for binary (within numerical tolerance)
        assert np.allclose(res_std.GM.values, res_uem.GM.values, atol=1e-6)


@select_database("alcrni.tdb")
def test_uem_ternary_finite(load_database):
    """UEM should give finite results for ternary systems."""
    dbf = load_database()
    comps = ['AL', 'CR', 'NI', 'VA']
    phases = ['LIQUID']

    conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.33, v.X('CR'): 0.33, v.N: 1}
    res = calculate(dbf, comps, phases, conds, model=ModelUEM)

    # Should be finite
    assert np.all(np.isfinite(res.GM.values))
    assert np.all(np.isfinite(res.HM.values))
    assert np.all(np.isfinite(res.SM.values))


@select_database("alcrni.tdb")
def test_uem_vs_muggianu_differ(load_database):
    """UEM and Muggianu should differ for ternary systems."""
    dbf = load_database()
    comps = ['AL', 'CR', 'NI', 'VA']
    phases = ['LIQUID']

    conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.33, v.X('CR'): 0.33, v.N: 1}
    res_std = calculate(dbf, comps, phases, conds)
    res_uem = calculate(dbf, comps, phases, conds, model=ModelUEM)

    # Should be different (this is expected and correct)
    # We just check they're both valid and not identical
    assert np.all(np.isfinite(res_std.GM.values))
    assert np.all(np.isfinite(res_uem.GM.values))
    # Typically differ by 1-15%
    diff = np.abs(res_uem.GM.values - res_std.GM.values)
    assert np.any(diff > 100)  # At least 100 J/mol difference expected


@select_database("alcrni.tdb")
def test_uem_pure_components(load_database):
    """UEM should be stable at pure component limits."""
    dbf = load_database()
    comps = ['AL', 'CR', 'NI', 'VA']
    phases = ['LIQUID']

    # Pure AL
    conds = {v.T: 1800, v.P: 101325, v.X('AL'): 1.0, v.X('CR'): 0.0, v.N: 1}
    res = calculate(dbf, comps, phases, conds, model=ModelUEM)
    assert np.all(np.isfinite(res.GM.values))

    # Pure CR
    conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.0, v.X('CR'): 1.0, v.N: 1}
    res = calculate(dbf, comps, phases, conds, model=ModelUEM)
    assert np.all(np.isfinite(res.GM.values))


def test_uem_symmetric_binary():
    """UEM property difference should be zero for symmetric binary (L0 only)."""
    from pycalphad.models.uem_symbolic import uem1_delta_expr

    dbf = Database()
    dbf.add_structure_entry('A', 'LIQUID', 0)
    dbf.add_structure_entry('B', 'LIQUID', 0)
    dbf.add_parameter(phase_name='LIQUID', parameter_type='G',
                     constituents=[['A', 'B']], parameter_order=0,
                     parameter=-10000.0, diffusing_species=None)

    # For symmetric interaction (L0 only), delta should be zero
    delta = uem1_delta_expr(dbf, 'A', 'B', 'LIQUID', v.T)

    # Evaluate at T=1000
    from sympy import lambdify
    delta_func = lambdify(v.T, delta, 'numpy')
    delta_val = float(delta_func(1000))

    assert abs(delta_val) < 1e-10  # Should be essentially zero


def test_uem_asymmetric_binary():
    """UEM property difference should be non-zero for asymmetric binary."""
    from pycalphad.models.uem_symbolic import uem1_delta_expr

    dbf = Database()
    dbf.add_structure_entry('A', 'LIQUID', 0)
    dbf.add_structure_entry('B', 'LIQUID', 0)
    dbf.add_parameter(phase_name='LIQUID', parameter_type='G',
                     constituents=[['A', 'B']], parameter_order=0,
                     parameter=-10000.0, diffusing_species=None)
    dbf.add_parameter(phase_name='LIQUID', parameter_type='G',
                     constituents=[['A', 'B']], parameter_order=1,
                     parameter=-3000.0, diffusing_species=None)

    # For asymmetric interaction (L0 + L1), delta should be non-zero
    delta = uem1_delta_expr(dbf, 'A', 'B', 'LIQUID', v.T)

    from sympy import lambdify
    delta_func = lambdify(v.T, delta, 'numpy')
    delta_val = float(delta_func(1000))

    assert delta_val > 1e-6  # Should be clearly non-zero
