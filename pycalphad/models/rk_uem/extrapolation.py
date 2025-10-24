"""
外推模型模块

该模块实现了多种外推模型，用于从二元系性质推断多元系性质。

包含的模型:
-----------
- UEM1: 非交互作用条件下的统一外推模型
- UEM2_O, UEM2_N, UEM2_Adv: 考虑交互作用的UEM变体
- Kohler: Kohler模型
- Muggianu: Muggianu模型
- Toop_Kohler: Toop-Kohler模型
- Toop_Muggianu: Toop-Muggianu模型
- GSM: 几何相似模型

UEM的核心思想:
通过计算贡献系数将第三组分的影响分配到i-j二元对上，
统一了传统的Kohler、Muggianu、Toop等外推方法。
"""

import math
from typing import Callable
from scipy import integrate

from pycalphad.models.rk_uem.rk_binary import RKBinaryPolynomial


def exp(x: float) -> float:
    """指数函数，带有数值保护"""
    try:
        return math.exp(min(x, 700))  # 防止溢出
    except (OverflowError, ValueError):
        return math.exp(700) if x > 0 else 0.0


def safe_float(expr) -> float:
    """
    安全地将表达式转换为浮点数

    Parameters
    ----------
    expr : various
        需要转换的表达式

    Returns
    -------
    float
        转换后的浮点值
    """
    try:
        return float(expr)
    except (TypeError, ValueError):
        try:
            from sympy import N
            return float(N(expr))
        except (TypeError, ValueError, ZeroDivisionError):
            return 0.0


def _asymmetric_component_choose(i: str, j: str, k: str, T: float, database_path: str) -> str:
    """
    选择非对称组分（用于Toop模型）

    Parameters
    ----------
    i, j, k : str
        三个组分
    T : float
        温度 (K)
    database_path : str
        数据库路径

    Returns
    -------
    str
        非对称组分
    """
    bij = RKBinaryPolynomial((i, j), database_path)
    bik = RKBinaryPolynomial((i, k), database_path)
    bjk = RKBinaryPolynomial((j, k), database_path)

    a = bij.mixing_enthalpy(component=i, x_component=0.5, temperature=T) or 0.0
    b = bik.mixing_enthalpy(component=i, x_component=0.5, temperature=T) or 0.0
    c = bjk.mixing_enthalpy(component=j, x_component=0.5, temperature=T) or 0.0

    if (a > 0 and b > 0 and c > 0) or (a < 0 and b < 0 and c < 0):
        if a > 0:
            return k if a == min(a, b, c) else (j if b == min(a, b, c) else i)
        else:
            return k if a == max(a, b, c) else (j if b == max(a, b, c) else i)
    else:
        if a * b > 0:
            return i
        if a * c > 0:
            return j
        if b * c > 0:
            return k

    return i  # 默认返回


def _deviation_function(k: str, i: str, j: str, T: float, database_path: str) -> float:
    """
    计算最小二元偏差项（用于GSM模型）

    Parameters
    ----------
    k, i, j : str
        组分
    T : float
        温度 (K)
    database_path : str
        数据库路径

    Returns
    -------
    float
        偏差函数值
    """
    bij = RKBinaryPolynomial((i, j), database_path)
    bik = RKBinaryPolynomial((i, k), database_path)

    func_ij = lambda x: bij.mixing_enthalpy(component=i, x_component=x, temperature=T) or 0.0
    func_ik = lambda x: bik.mixing_enthalpy(component=i, x_component=x, temperature=T) or 0.0

    integrand = lambda x: pow(func_ij(x) - func_ik(x), 2)

    try:
        yeta = integrate.quad(integrand, 0, 1)[0]
        return yeta
    except Exception:
        return 0.0


def _get_graph_center(k: str, i: str, T: float, database_path: str) -> tuple:
    """
    计算函数与坐标轴围成图形的中心坐标（用于UEM2_Adv）

    Parameters
    ----------
    k, i : str
        组分
    T : float
        温度 (K)
    database_path : str
        数据库路径

    Returns
    -------
    tuple
        (x_center, y_center)
    """
    bki = RKBinaryPolynomial((k, i), database_path)
    f = lambda x: bki.mixing_enthalpy(component=k, x_component=x, temperature=T) or 0.0

    try:
        xf = lambda x: x * f(x)
        ff = lambda x: f(x) * f(x)

        x_ = integrate.quad(xf, 0, 1, limit=200, epsabs=1e-6, epsrel=1e-6)[0]
        A = integrate.quad(f, 0, 1, limit=200, epsabs=1e-6, epsrel=1e-6)[0]
        y_ = integrate.quad(ff, 0, 1, limit=200, epsabs=1e-6, epsrel=1e-6)[0]

        if abs(A) < 1e-10:
            return 0.0, 0.0

        x = x_ / A
        y = y_ / (2 * A)

        return x - 0.5, y
    except Exception:
        return 0.0, 0.0


def _get_dki_by_uem2adv(k: str, i: str, j: str, T: float, database_path: str) -> float:
    """
    计算UEM2-Adv中的组分间性质差

    Parameters
    ----------
    k, i, j : str
        组分
    T : float
        温度 (K)
    database_path : str
        数据库路径

    Returns
    -------
    float
        性质差
    """
    xkj, ykj = _get_graph_center(k, j, T, database_path)
    xij, yij = _get_graph_center(i, j, T, database_path)

    theta10 = math.atan2(yij, xij)
    theta20 = math.atan2(ykj, xkj)

    a = math.sqrt(xij * xij + yij * yij)
    b = math.sqrt(xkj * xkj + ykj * ykj)

    # 角度象限判断
    if abs(abs(theta10 - theta20) - math.pi/2) < 0.1:
        delta_theta_corr = math.pi / 2
    else:
        delta_theta_corr = 0

    theta_sum_sq = theta10 * theta10 + theta20 * theta20
    if abs(theta_sum_sq) < 1e-10:
        delta_theta = 0.0
    else:
        delta_theta = abs((math.pi / 2.0 * (theta10 * theta10 - theta20 * theta20) + delta_theta_corr) / theta_sum_sq)

    d_ki = delta_theta * delta_theta

    return d_ki


# ============================================================================
# 外推模型函数
# ============================================================================

def Kohler(k: str, i: str, j: str, T: float = 298.0, database_path: str = None) -> float:
    """
    Kohler模型贡献系数

    在Kohler模型中，第三组分的贡献在两个组分间平均分配。

    Parameters
    ----------
    k : str
        第三组分
    i, j : str
        二元对组分
    T : float
        温度 (K)
    database_path : str
        数据库路径（此模型不需要）

    Returns
    -------
    float
        贡献系数 (总是0.0，表示完全对称)
    """
    return 0.0


def Muggianu(k: str, i: str, j: str, T: float = 298.0, database_path: str = None) -> float:
    """
    Muggianu模型贡献系数

    在Muggianu模型中，第三组分的贡献均等分配。

    Parameters
    ----------
    k : str
        第三组分
    i, j : str
        二元对组分
    T : float
        温度 (K)
    database_path : str
        数据库路径（此模型不需要）

    Returns
    -------
    float
        贡献系数 (总是0.5)
    """
    return 0.5


def Toop_Kohler(k: str, i: str, j: str, T: float = 298.0, database_path: str = None) -> float:
    """
    Toop-Kohler模型贡献系数

    基于非对称组分选择的Toop模型变体。

    Parameters
    ----------
    k : str
        第三组分
    i, j : str
        二元对组分
    T : float
        温度 (K)
    database_path : str
        数据库路径

    Returns
    -------
    float
        贡献系数
    """
    if not database_path:
        return 0.5

    asmym = _asymmetric_component_choose(k, i, j, T, database_path)

    if asmym == k or asmym == i:
        return 0.0
    else:
        return 1.0


def Toop_Muggianu(k: str, i: str, j: str, T: float = 298.0, database_path: str = None) -> float:
    """
    Toop-Muggianu模型贡献系数

    Parameters
    ----------
    k : str
        第三组分
    i, j : str
        二元对组分
    T : float
        温度 (K)
    database_path : str
        数据库路径

    Returns
    -------
    float
        贡献系数
    """
    if not database_path:
        return 0.5

    asmym = _asymmetric_component_choose(k, i, j, T, database_path)

    if asmym == k:
        return 0.5
    elif asmym == i:
        return 0.0
    else:
        return 1.0


def UEM1(k: str, i: str, j: str, T: float, database_path: str) -> float:
    """
    UEM1 - 非交互作用条件下的统一外推模型

    基于无限稀释性质差计算贡献系数。

    Parameters
    ----------
    k : str
        第三组分
    i, j : str
        二元对组分
    T : float
        温度 (K)
    database_path : str
        数据库路径

    Returns
    -------
    float
        贡献系数
    """
    if not database_path:
        return 0.5

    try:
        bki = RKBinaryPolynomial((k, i), database_path)
        bik = RKBinaryPolynomial((i, k), database_path)
        bjk = RKBinaryPolynomial((j, k), database_path)
        bkj = RKBinaryPolynomial((k, j), database_path)

        RT = 8.314 * T

        wik = (bik.infinite_dilution_property(component=i, temperature=T) or 0.0) / RT
        wki = (bki.infinite_dilution_property(component=k, temperature=T) or 0.0) / RT
        wjk = (bjk.infinite_dilution_property(component=j, temperature=T) or 0.0) / RT
        wkj = (bkj.infinite_dilution_property(component=k, temperature=T) or 0.0) / RT

        df_ki = abs(wik - wki)
        df_kj = abs(wjk - wkj)

        if df_ki == 0 and df_kj == 0:
            return 0.5
        else:
            return (df_kj / (df_ki + df_kj)) * exp(-df_ki)

    except Exception:
        return 0.5


def UEM2_O(k: str, i: str, j: str, T: float, database_path: str) -> float:
    """
    UEM2-O - 考虑交互作用的UEM（原始版本）

    基于偏摩尔性质积分计算贡献系数。

    Parameters
    ----------
    k : str
        第三组分
    i, j : str
        二元对组分
    T : float
        温度 (K)
    database_path : str
        数据库路径

    Returns
    -------
    float
        贡献系数
    """
    if not database_path:
        return 0.5

    try:
        bkj = RKBinaryPolynomial((k, j), database_path)
        bji = RKBinaryPolynomial((j, i), database_path)
        bki = RKBinaryPolynomial((k, i), database_path)
        bij = RKBinaryPolynomial((i, j), database_path)

        RT = 8.314 * T

        def safe_integrand(binary_sys, comp):
            def integrand(x):
                try:
                    result = binary_sys.partial_molar_property(component=comp, x_component=x, temperature=T)
                    if math.isinf(result) or math.isnan(result):
                        return 0.0
                    return result / RT
                except Exception:
                    return 0.0
            return integrand

        Wki = integrate.quad(safe_integrand(bki, k), 0, 1, limit=200, epsabs=1e-6, epsrel=1e-6)[0]
        Wkj = integrate.quad(safe_integrand(bkj, k), 0, 1, limit=200, epsabs=1e-6, epsrel=1e-6)[0]
        Wij = integrate.quad(safe_integrand(bij, i), 0, 1, limit=200, epsabs=1e-6, epsrel=1e-6)[0]
        Wji = integrate.quad(safe_integrand(bji, j), 0, 1, limit=200, epsabs=1e-6, epsrel=1e-6)[0]

        dki = abs(Wkj - Wij)
        dkj = abs(Wki - Wji)

        if dki == 0 and dkj == 0:
            return 0.5

        return (dkj / (dki + dkj)) * exp(-dki)

    except Exception:
        return 0.5


def UEM2_N(k: str, i: str, j: str, T: float, database_path: str) -> float:
    """
    UEM2-N - 基于混合焓面积的UEM

    Parameters
    ----------
    k : str
        第三组分
    i, j : str
        二元对组分
    T : float
        温度 (K)
    database_path : str
        数据库路径

    Returns
    -------
    float
        贡献系数
    """
    if not database_path:
        return 0.5

    try:
        bij = RKBinaryPolynomial((i, j), database_path)
        bkj = RKBinaryPolynomial((k, j), database_path)
        bki = RKBinaryPolynomial((k, i), database_path)
        bji = RKBinaryPolynomial((j, i), database_path)

        RT = 8.314 * T

        def safe_mixing_enthalpy(binary_sys, comp):
            def integrand(x):
                try:
                    result = binary_sys.mixing_enthalpy(component=comp, x_component=x, temperature=T)
                    if result is None or math.isinf(result) or math.isnan(result):
                        return 0.0
                    return result * 2 / RT
                except Exception:
                    return 0.0
            return integrand

        Aij = integrate.quad(safe_mixing_enthalpy(bij, i), 0, 1, limit=200, epsabs=1e-6, epsrel=1e-6)[0]
        Akj = integrate.quad(safe_mixing_enthalpy(bkj, k), 0, 1, limit=200, epsabs=1e-6, epsrel=1e-6)[0]
        Aki = integrate.quad(safe_mixing_enthalpy(bki, k), 0, 1, limit=200, epsabs=1e-6, epsrel=1e-6)[0]
        Aji = integrate.quad(safe_mixing_enthalpy(bji, j), 0, 1, limit=200, epsabs=1e-6, epsrel=1e-6)[0]

        D_ki = abs((Aij - Akj) / (Aij + Akj)) if (Aij + Akj) != 0 else 0.0
        D_kj = abs((Aji - Aki) / (Aji + Aki)) if (Aji + Aki) != 0 else 0.0

        if D_ki == 0 and D_kj == 0:
            return 0.5

        return (D_kj / (D_ki + D_kj)) * exp(-D_ki)

    except Exception:
        return 0.5


def UEM2_Adv(k: str, i: str, j: str, T: float, database_path: str) -> float:
    """
    UEM2-Adv - 基于几何性质差的高级UEM

    考虑混合焓曲线的几何形状差异。

    Parameters
    ----------
    k : str
        第三组分
    i, j : str
        二元对组分
    T : float
        温度 (K)
    database_path : str
        数据库路径

    Returns
    -------
    float
        贡献系数
    """
    if not database_path:
        return 0.5

    try:
        dki = _get_dki_by_uem2adv(k, i, j, T, database_path)
        dkj = _get_dki_by_uem2adv(k, j, i, T, database_path)

        if dki == 0 and dkj == 0:
            return 0.5

        return (dkj / (dki + dkj)) * exp(-dki)

    except Exception:
        return 0.5


def GSM(k: str, i: str, j: str, T: float, database_path: str) -> float:
    """
    GSM - 几何相似模型

    基于二元系混合焓曲线的偏差函数计算贡献系数。

    Parameters
    ----------
    k : str
        第三组分
    i, j : str
        二元对组分
    T : float
        温度 (K)
    database_path : str
        数据库路径

    Returns
    -------
    float
        贡献系数
    """
    if not database_path:
        return 0.5

    try:
        n1 = _deviation_function(k=k, i=i, j=j, T=T, database_path=database_path)
        n2 = _deviation_function(k=k, i=j, j=i, T=T, database_path=database_path)

        if n1 == 0 and n2 == 0:
            return 0.5
        else:
            return n1 / (n1 + n2)

    except Exception:
        return 0.5
