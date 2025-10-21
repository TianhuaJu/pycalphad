"""
Unified Extrapolation Model (UEM) for multicomponent thermodynamics.

Implements property-difference-based extrapolation from binary to multicomponent systems.
Based on Ju et al. (2024) Thermochimica Acta 740, 179824.
"""
from sympy import Symbol, Add, Mul, Pow, Abs, exp, simplify, S, nan
from pycalphad import variables as v
from pycalphad.variables import R
from tinydb import where
import logging

logger = logging.getLogger(__name__)


def _is_binary(arr, comp1, comp2):
    """Check if constituent array contains only comp1 and comp2."""
    comps = set(str(s) for subl in arr for s in subl)
    return comp1 in comps and comp2 in comps and len(comps) == 2


def uem1_delta_expr(dbe, comp1, comp2, phase_name, T):
    """
    Calculate property difference delta_ij from binary parameters.

    delta_ij = |dG_ex/dx|_{x=0} - |dG_ex/dx|_{x=1}| / (R*T)

    Based on Ju et al. (2024) non-interactive property difference method.
    """
    x = Symbol('x')
    G_ex = S.Zero

    # Get binary interaction parameters
    query = ((where('phase_name') == phase_name) &
             (where('parameter_type') == 'G') &
             (where('constituent_array').test(lambda arr: _is_binary(arr, comp1, comp2))))
    params = dbe.search(query)

    if not params:
        return S.Zero

    # Build Redlich-Kister expression: G_ex = x*(1-x)*sum(L_n*(2x-1)^n)
    for p in params:
        order = p['parameter_order']
        param = p['parameter']
        G_ex += x * (1 - x) * param * (2*x - 1)**order

    # Calculate boundary derivatives
    dGdx_0 = G_ex.diff(x).subs(x, 0)
    dGdx_1 = G_ex.subs(x, 1-x).diff(x).subs(x, 0)

    # Property difference
    delta = Abs(dGdx_0 - dGdx_1) / (R * T)
    return simplify(delta)


def uem1_contribution_ratio(dbe, k, i, j, phase_name, T):
    """
    Calculate contribution coefficient r_ki.

    r_ki = (delta_kj / (delta_ki + delta_kj)) * exp(-delta_ki)

    Determines how much component k contributes to i in the i-j pair.
    """
    delta_ki = uem1_delta_expr(dbe, k, i, phase_name, T)
    delta_kj = uem1_delta_expr(dbe, k, j, phase_name, T)

    # Handle edge cases
    if delta_ki == S.Zero and delta_kj == S.Zero:
        return S.Half

    delta_sum = delta_ki + delta_kj
    if delta_sum == S.Zero:
        # If they sum to zero, return 0.5 (neutral contribution)
        return S.Half

    return simplify((delta_kj / delta_sum) * exp(-delta_ki))


def _binary_excess(dbe, comp_i, comp_j, phase_name, x_eff_i, x_eff_j):
    """Build binary excess energy with effective mole fractions."""
    G_ex = S.Zero

    query = ((where('phase_name') == phase_name) &
             (where('parameter_type') == 'G') &
             (where('constituent_array').test(lambda arr: _is_binary(arr, comp_i, comp_j))))
    params = dbe.search(query)

    for p in params:
        order = p['parameter_order']
        param = p['parameter']
        G_ex += x_eff_i * x_eff_j * param * (x_eff_i - x_eff_j)**order

    return G_ex


def get_uem1_excess_gibbs_expr(dbe, comps, phase_name, T, subl_index=0):
    """
    Calculate UEM excess Gibbs energy for multicomponent system.

    Algorithm:
    1. For each binary pair (i,j):
       - Calculate effective mole fractions: x_eff_i = x_i + sum(r_ki * x_k)
       - Normalize: X_ij = x_eff_i / (x_eff_i + x_eff_j)
       - Build binary excess: G_ex_ij(X_ij, X_ji)
       - Weight and add to total

    Parameters
    ----------
    dbe : Database
    comps : list of str
        Component names (excluding VA)
    phase_name : str
    T : StateVariable
    subl_index : int
        Sublattice index (default 0 for single sublattice phases like LIQUID)

    Returns
    -------
    SymPy expression
        Total excess Gibbs energy (J/mol)
    """
    # Use site fractions instead of mole fractions
    x = {comp: v.SiteFraction(phase_name, subl_index, comp) for comp in comps}
    expr_list = []

    # For binary systems, UEM reduces to standard Redlich-Kister
    if len(comps) == 2:
        comp_i, comp_j = comps[0], comps[1]
        G_ex = _binary_excess(dbe, comp_i, comp_j, phase_name, x[comp_i], x[comp_j])
        return G_ex

    # Iterate binary pairs for ternary and higher systems
    for i in range(len(comps)):
        for j in range(i + 1, len(comps)):
            comp_i, comp_j = comps[i], comps[j]

            # Calculate effective mole fractions
            x_eff_i = x[comp_i]
            x_eff_j = x[comp_j]

            for k in comps:
                if k not in [comp_i, comp_j]:
                    r_ki = uem1_contribution_ratio(dbe, k, comp_i, comp_j, phase_name, T)
                    r_kj = uem1_contribution_ratio(dbe, k, comp_j, comp_i, phase_name, T)
                    x_eff_i += r_ki * x[k]
                    x_eff_j += r_kj * x[k]

            # Normalize to i-j subsystem
            denom = x_eff_i + x_eff_j

            # Skip if denominator would be zero
            if denom == S.Zero:
                continue

            Xi_ij = x_eff_i / denom
            Xj_ij = x_eff_j / denom

            # Binary excess with effective mole fractions
            G_ex_ij = _binary_excess(dbe, comp_i, comp_j, phase_name, Xi_ij, Xj_ij)

            # Weight by (x_i*x_j) / (Xi_ij*Xj_ij)
            # This simplifies to: (x_i * x_j) * denom^2 / (x_eff_i * x_eff_j)
            if x_eff_i != S.Zero and x_eff_j != S.Zero:
                weight = (x[comp_i] * x[comp_j] * denom * denom) / (x_eff_i * x_eff_j)
                expr_list.append(weight * G_ex_ij)

    if not expr_list:
        return S.Zero

    return Add(*expr_list)
