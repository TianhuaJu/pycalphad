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
        # Use equal_nan=True to handle NaN values in sampled points
        assert np.allclose(res_std.GM.values, res_uem.GM.values, atol=1e-6, equal_nan=True)


@select_database("alcrni.tdb")
def test_uem_ternary_finite(load_database):
    """UEM should give finite results for ternary systems."""
    dbf = load_database()
    comps = ['AL', 'CR', 'NI', 'VA']
    phases = ['LIQUID']

    conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.33, v.X('CR'): 0.33, v.N: 1}
    res = calculate(dbf, comps, phases, conds, model=ModelUEM)

    # Check that there are SOME finite values (not all NaN)
    # Some points may be NaN due to invalid compositions in the sampling
    assert np.any(np.isfinite(res.GM.values))
    # Check that finite values exist for all properties
    gm_finite = res.GM.values[np.isfinite(res.GM.values)]
    assert len(gm_finite) > 100  # Should have many valid points


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
    # Filter to only compare finite values (some points may be NaN)
    valid_mask = np.isfinite(res_std.GM.values) & np.isfinite(res_uem.GM.values)
    assert np.sum(valid_mask) > 100  # Should have many valid points

    # Typically differ by 1-15%
    diff = np.abs(res_uem.GM.values[valid_mask] - res_std.GM.values[valid_mask])
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
    # Check that we have valid finite values (some sampled points may be NaN)
    assert np.any(np.isfinite(res.GM.values))

    # Pure CR
    conds = {v.T: 1800, v.P: 101325, v.X('AL'): 0.0, v.X('CR'): 1.0, v.N: 1}
    res = calculate(dbf, comps, phases, conds, model=ModelUEM)
    assert np.any(np.isfinite(res.GM.values))


def test_uem_symmetric_binary():
    """UEM property difference should be zero for symmetric binary (L0 only)."""
    from pycalphad.models.uem_symbolic import uem1_delta_expr

    # Create a simple binary database with symmetric interaction (L0 only)
    tdb_str = """
    ELEMENT A LIQUID 0 0 0 !
    ELEMENT B LIQUID 0 0 0 !

    PHASE LIQUID % 1 1 !
    CONSTITUENT LIQUID : A,B : !

    PARAMETER G(LIQUID,A;0) 298.15 0; 6000 N !
    PARAMETER G(LIQUID,B;0) 298.15 0; 6000 N !
    PARAMETER G(LIQUID,A,B;0) 298.15 -10000; 6000 N !
    """

    dbf = Database(tdb_str)

    # For symmetric interaction (L0 only), delta should be zero
    delta = uem1_delta_expr(dbf, 'A', 'B', 'LIQUID', v.T)

    # Evaluate at T=1000
    from symengine import lambdify
    delta_func = lambdify([v.T], delta, 'llvm')
    delta_val = float(delta_func(1000))

    assert abs(delta_val) < 1e-10  # Should be essentially zero


def test_uem_asymmetric_binary():
    """UEM property difference should be non-zero for asymmetric binary."""
    from pycalphad.models.uem_symbolic import uem1_delta_expr

    # Create a binary database with asymmetric interaction (L0 + L1)
    tdb_str = """
    ELEMENT A LIQUID 0 0 0 !
    ELEMENT B LIQUID 0 0 0 !

    PHASE LIQUID % 1 1 !
    CONSTITUENT LIQUID : A,B : !

    PARAMETER G(LIQUID,A;0) 298.15 0; 6000 N !
    PARAMETER G(LIQUID,B;0) 298.15 0; 6000 N !
    PARAMETER G(LIQUID,A,B;0) 298.15 -10000; 6000 N !
    PARAMETER G(LIQUID,A,B;1) 298.15 -3000; 6000 N !
    """

    dbf = Database(tdb_str)

    # For asymmetric interaction (L0 + L1), delta should be non-zero
    delta = uem1_delta_expr(dbf, 'A', 'B', 'LIQUID', v.T)

    from symengine import lambdify
    delta_func = lambdify([v.T], delta, 'llvm')
    delta_val = float(delta_func(1000))

    assert delta_val > 1e-6  # Should be clearly non-zero
