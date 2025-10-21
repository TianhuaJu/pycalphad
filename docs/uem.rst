Unified Extrapolation Model (UEM)
===================================

The UEM provides an alternative extrapolation method for predicting multicomponent
solution properties from binary data, based on property differences rather than
geometric averaging.

Overview
--------

**Binary description**: Same Redlich-Kister polynomials as standard Model

**Multicomponent extrapolation**:
  - Standard: Muggianu (geometric symmetric averaging)
  - UEM: Property-difference-based effective mole fractions

**Key principle**: Components with similar properties contribute more to each other's
effective mole fractions in binary subsystems.

Usage
-----

.. code-block:: python

    from pycalphad import Database, calculate
    from pycalphad.models.model_uem import ModelUEM

    dbf = Database('mydb.tdb')

    # Standard Muggianu
    result_std = calculate(dbf, ['AL', 'CR', 'NI', 'VA'], ['LIQUID'],
                          T=1800, P=101325, X_AL=0.33, X_CR=0.33)

    # UEM extrapolation
    result_uem = calculate(dbf, ['AL', 'CR', 'NI', 'VA'], ['LIQUID'],
                          model=ModelUEM,
                          T=1800, P=101325, X_AL=0.33, X_CR=0.33)

Expected Behavior
-----------------

- **Binary systems**: UEM = Standard (identical results)
- **Ternary+ systems**: UEM ≠ Muggianu (differences expected and correct)

The UEM typically predicts different values (1-15% difference) for multicomponent
systems, especially those with highly asymmetric binary subsystems.

Algorithm
---------

For each binary pair (i,j):

1. Calculate property difference: ``δ_ij = |dG_ex/dx|_{x=0} - dG_ex/dx|_{x=1}| / (RT)``
2. Calculate contribution coefficients: ``r_ki = (δ_kj / (δ_ki + δ_kj)) * exp(-δ_ki)``
3. Build effective mole fractions: ``x_eff_i = x_i + Σ_k(r_ki * x_k)``
4. Normalize and construct binary excess with effective fractions
5. Sum weighted contributions from all pairs

References
----------

Primary reference:

Ju, T., Huang, Z., Ding, X., Yan, X. & Liao, C. (2024).
"A Unified Extrapolation thermodynamic model for multicomponent solutions
based on binary data." *Thermochimica Acta*, 740, 179824.
https://doi.org/10.1016/j.tca.2024.179824

Earlier work:

- Ju, T., et al. (2020). *Fluid Phase Equilibria*, 507, 112416.
- Ju, T., et al. (2020). *J. Molecular Liquids*, 298, 111951.

See Also
--------

- :ref:`models`
- :class:`pycalphad.Model`
- :class:`pycalphad.models.model_uem.ModelUEM`
