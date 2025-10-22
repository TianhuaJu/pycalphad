"""
Redlich-Kister-UEM: Property-difference-based multicomponent extrapolation.

This implements the UEM extrapolation method as an alternative to Muggianu/Kohler/Toop
for predicting multicomponent solution properties from binary data.

References
----------
Ju, T., Huang, Z., Ding, X., Yan, X. & Liao, C. (2024).
Thermochimica Acta, 740, 179824. DOI: 10.1016/j.tca.2024.179824
"""
from pycalphad import Model, variables as v
from sympy import S
import pycalphad.models.uem_symbolic as uem
import logging

logger = logging.getLogger(__name__)


class ModelUEM(Model):
    """
    UEM extrapolation model for solution phases.

    Uses property-difference-based effective mole fractions instead of
    geometric averaging (Muggianu) for multicomponent extrapolation.

    For binary systems, gives identical results to standard Redlich-Kister.
    For ternary+, provides alternative predictions based on component similarity.

    Examples
    --------
    >>> from pycalphad import Database, calculate
    >>> from pycalphad.models.model_uem import ModelUEM
    >>> dbf = Database('mydb.tdb')
    >>> result = calculate(dbf, ['AL', 'CR', 'NI', 'VA'], ['LIQUID'],
    ...                    model=ModelUEM, T=1800, P=101325, X_AL=0.33, X_CR=0.33)
    """

    contributions = [
        ('ref', 'reference_energy'),
        ('idmix', 'ideal_mixing_energy'),
        ('xsmix', 'excess_mixing_energy')
    ]

    def excess_mixing_energy(self, dbe):
        """
        UEM-based excess mixing energy.

        Uses effective mole fractions calculated from property differences
        instead of geometric averaging.

        For binary systems, UEM is mathematically equivalent to standard
        Redlich-Kister, so we use the parent class implementation for efficiency.
        """
        comps = [str(c) for c in self.components if str(c) != 'VA']

        if len(comps) < 2:
            return S.Zero

        # For binary systems, UEM = standard RK (use parent class for efficiency)
        if len(comps) == 2:
            return super(ModelUEM, self).excess_mixing_energy(dbe)

        # For ternary+, use UEM extrapolation
        expr = uem.get_uem1_excess_gibbs_expr(dbe, comps, self.phase_name, v.T)
        return expr / self._site_ratio_normalization
