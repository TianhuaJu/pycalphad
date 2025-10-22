"""
R-K-UEM (Redlich-Kister Unified Excess Model) 模块

该模块实现了基于R-K多项式的二元系统描述和UEM外推方法的多元系统热力学计算。

主要组件:
-----------
- rk_binary: R-K二元系多项式模块
- extrapolation: 外推模型模块（UEM, Kohler, Muggianu, Toop等）
- thermodynamics: 多元系热力学性质计算模块

UEM (Unified Excess Model) 是一种通过贡献系数统一传统外推模型的方法，
可以作为Kohler、Muggianu、Toop等模型的通用替代。
"""

from pycalphad.models.rk_uem.rk_binary import RKBinaryPolynomial
from pycalphad.models.rk_uem.extrapolation import (
    UEM1, UEM2_N, UEM2_O, UEM2_Adv,
    Kohler, Muggianu, Toop_Kohler, Toop_Muggianu,
    GSM
)
from pycalphad.models.rk_uem.thermodynamics import ThermodynamicCalculator

__all__ = [
    'RKBinaryPolynomial',
    'ThermodynamicCalculator',
    'UEM1', 'UEM2_N', 'UEM2_O', 'UEM2_Adv',
    'Kohler', 'Muggianu', 'Toop_Kohler', 'Toop_Muggianu',
    'GSM'
]

__version__ = '0.1.0'
