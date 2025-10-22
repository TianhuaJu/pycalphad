"""
R-K (Redlich-Kister) 二元系多项式模块

该模块提供基于Redlich-Kister多项式的二元系统热力学性质计算。

R-K多项式形式:
F = x_i * x_j * Σ(A_ij * (x_i - x_j)^k)

特性:
- 支持任意组分标识
- 自动识别二元系的等价性 (i-j 与 j-i 为同一系统)
- 动态温度依赖的参数计算
- 严格处理R-K多项式的数学对称性

注意: 由于 (x_j - x_i) = -(x_i - x_j)，为保持函数一致性:
- 偶数项参数: A0, A2, A4... 符号相同
- 奇数项参数: A1, A3, A5... 符号相反
"""

import math
import re
from typing import Tuple, Optional, Dict
from contextlib import contextmanager
import sqlite3
import os


def ln(x: float) -> float:
    """自然对数函数"""
    return math.log(x)


def sqrt(x: float) -> float:
    """平方根函数"""
    return math.sqrt(x)


class RKBinaryPolynomial:
    """
    R-K二元系统热力学计算类

    支持计算:
    - 过剩Gibbs自由能
    - 混合焓
    - 偏摩尔性质
    - 无限稀释性质

    Parameters
    ----------
    composition : tuple of (str, str)
        二元系统组成，如 ("Fe", "Cr"), ("A", "B") 等
    database_path : str, optional
        R-K参数数据库路径

    Attributes
    ----------
    component1 : str
        第一个组分
    component2 : str
        第二个组分
    system_id : str
        标准化的系统标识符
    """

    # 常数
    R = 8.314  # 气体常数 J/(mol·K)
    DEFAULT_T = 298.15  # 默认温度 K

    # 预编译正则表达式用于参数解析
    _PAT_A = re.compile(r'^(\-?\d+\.?\d*)$')
    _PAT_B = re.compile(r'^(\-?\d+\.?\d*)(\*)(T)$')
    _PAT_C = re.compile(r'^(\-?\d+\.?\d*)(\*)(T)(\*{2})(\-?\d+\.?\d*)$')
    _PAT_CL = re.compile(r'^(\-?\d+\.?\d*)(\*)(T)(\*lnT)$')

    def __init__(self, composition: Tuple[str, str], database_path: str = None):
        """
        初始化R-K二元系统

        Parameters
        ----------
        composition : tuple
            二元系统组成 (组分1, 组分2)
        database_path : str, optional
            R-K参数数据库路径
        """
        if len(composition) != 2:
            raise ValueError("组成必须包含恰好两个元素")

        # 保持用户输入的顺序
        self.component1, self.component2 = composition

        # 标准化组分顺序（用于数据库查询）
        self._standardized_components = tuple(sorted(composition))

        # 生成系统标识
        self.system_id = f"{self._standardized_components[0]}-{self._standardized_components[1]}"
        self.name = f"{self.component1}-{self.component2}"

        # 判断是否需要反转参数符号
        self._is_reversed = (self.component1, self.component2) != self._standardized_components

        # 数据库路径
        self.database_path = database_path

        # 存储原始字符串参数
        self._aij_strings: Tuple[str, ...] = ()
        self._hij_strings: Tuple[str, ...] = ()

        # 当前温度下的计算参数
        self.aij_values: Tuple[float, ...] = ()
        self.hij_values: Tuple[float, ...] = ()

        # 从数据库加载参数
        if database_path:
            self._load_parameters()

    @contextmanager
    def _get_db_connection(self):
        """数据库连接上下文管理器"""
        conn = None
        try:
            if not os.path.exists(self.database_path):
                raise RuntimeError(f"数据库文件不存在: {self.database_path}")

            conn = sqlite3.connect(self.database_path)
            yield conn
        except sqlite3.Error as e:
            raise RuntimeError(f"数据库操作失败: {e}")
        finally:
            if conn:
                conn.close()

    def _load_parameters(self) -> None:
        """从数据库加载R-K参数"""
        self._aij_strings = self._query_parameters('R_K_gE')
        self._hij_strings = self._query_parameters('R_K_Hmix')

    def _query_parameters(self, table: str) -> Tuple[str, ...]:
        """
        从数据库查询参数

        Parameters
        ----------
        table : str
            表名 ('R_K_gE' 或 'R_K_Hmix')

        Returns
        -------
        tuple
            参数字符串元组
        """
        std_symbol = f"{self._standardized_components[0]}-{self._standardized_components[1]}"
        reverse_symbol = f"{self._standardized_components[1]}-{self._standardized_components[0]}"

        with self._get_db_connection() as conn:
            cursor = conn.cursor()
            query = f"SELECT A0,A1,A2,A3,A4,A5 FROM {table} WHERE Symbol = ?"

            # 尝试标准顺序
            cursor.execute(query, (std_symbol,))
            data = cursor.fetchall()

            if data:
                return self._process_query_result(data, reverse_odd=False)

            # 尝试反向顺序
            cursor.execute(query, (reverse_symbol,))
            data = cursor.fetchall()

            if data:
                return self._process_query_result(data, reverse_odd=True)

            return ()

    def _process_query_result(self, data: list, reverse_odd: bool) -> Tuple[str, ...]:
        """
        处理数据库查询结果，考虑R-K多项式的对称性

        由于 (x_j - x_i) = -(x_i - x_j)，为保持函数值相同:
        - 偶数项系数相同
        - 奇数项系数相反

        Parameters
        ----------
        data : list
            查询结果
        reverse_odd : bool
            是否对奇数项取负号

        Returns
        -------
        tuple
            参数字符串元组
        """
        parameters = ()
        for row in data:
            for i, param_str in enumerate(row):
                if param_str is not None:
                    if reverse_odd and i % 2 == 1:
                        parameters += (f"-({param_str})",)
                    else:
                        parameters += (param_str,)
                else:
                    parameters += ("0",)
        return parameters

    def _calculate_parameter_at_temperature(self, param_str: str, temperature: float) -> float:
        """
        计算指定温度下的参数值

        支持的表达式格式:
        - 常数: a
        - 线性: b*T
        - 幂次: c*T**n
        - 对数: c*T*lnT

        Parameters
        ----------
        param_str : str
            参数字符串表达式
        temperature : float
            温度 (K)

        Returns
        -------
        float
            计算出的参数值
        """
        if not param_str or param_str == "0":
            return 0.0

        # 处理负号包围的表达式
        if param_str.startswith("-(") and param_str.endswith(")"):
            return -self._calculate_parameter_at_temperature(param_str[2:-1], temperature)

        str_parts = param_str.split('+')
        total = 0.0

        for part in str_parts:
            part = part.strip()
            total += self._evaluate_parameter_part(part, temperature)

        return total

    def _evaluate_parameter_part(self, part: str, temperature: float) -> float:
        """评估参数表达式的单个部分"""
        # 常数项
        match = self._PAT_A.match(part)
        if match:
            return float(match.group(1))

        # 线性项
        match = self._PAT_B.match(part)
        if match:
            return float(match.group(1)) * temperature

        # 幂次项
        match = self._PAT_C.match(part)
        if match:
            coeff = float(match.group(1))
            exponent = float(match.group(5))
            return coeff * (temperature ** exponent)

        # 对数项
        match = self._PAT_CL.match(part)
        if match:
            coeff = float(match.group(1))
            return coeff * temperature * ln(temperature)

        return 0.0

    def _update_parameters_for_temperature(self, temperature: float) -> None:
        """更新指定温度下的参数值"""
        self.aij_values = tuple(
            self._calculate_parameter_at_temperature(param_str, temperature)
            for param_str in self._aij_strings
        )

        self.hij_values = tuple(
            self._calculate_parameter_at_temperature(param_str, temperature)
            for param_str in self._hij_strings
        )

    def excess_gibbs_energy(self, component: str, x_component: float,
                           temperature: float = DEFAULT_T) -> Optional[float]:
        """
        计算过剩Gibbs自由能

        Parameters
        ----------
        component : str
            组分标识
        x_component : float
            该组分的摩尔分数
        temperature : float
            温度 (K)

        Returns
        -------
        float or None
            过剩Gibbs自由能 (J/mol)，如果无数据返回None
        """
        if component not in (self.component1, self.component2):
            raise ValueError(f"组分 '{component}' 不在二元系 {self.name} 中")

        if not 0 <= x_component <= 1:
            raise ValueError("摩尔分数必须在0到1之间")

        self._update_parameters_for_temperature(temperature)

        if not self.aij_values:
            return None

        # 设置摩尔分数
        if component == self.component1:
            x1, x2 = x_component, 1 - x_component
        else:
            x1, x2 = 1 - x_component, x_component

        # 计算标准化顺序下的摩尔分数差
        if self._is_reversed:
            x_std1, x_std2 = x2, x1
        else:
            x_std1, x_std2 = x1, x2

        x_diff = x_std1 - x_std2
        rk_sum = sum(aij * (x_diff ** i) for i, aij in enumerate(self.aij_values))

        return x1 * x2 * rk_sum

    def mixing_enthalpy(self, component: str, x_component: float,
                       temperature: float = DEFAULT_T) -> Optional[float]:
        """
        计算混合焓

        Parameters
        ----------
        component : str
            组分标识
        x_component : float
            该组分的摩尔分数
        temperature : float
            温度 (K)

        Returns
        -------
        float or None
            混合焓 (J/mol)，如果无数据返回None
        """
        if component not in (self.component1, self.component2):
            raise ValueError(f"组分 '{component}' 不在二元系 {self.name} 中")

        if not 0 <= x_component <= 1:
            raise ValueError("摩尔分数必须在0到1之间")

        self._update_parameters_for_temperature(temperature)

        if not self.hij_values:
            return None

        # 设置摩尔分数
        if component == self.component1:
            x1, x2 = x_component, 1 - x_component
        else:
            x1, x2 = 1 - x_component, x_component

        # 计算标准化顺序下的摩尔分数差
        if self._is_reversed:
            x_std1, x_std2 = x2, x1
        else:
            x_std1, x_std2 = x1, x2

        x_diff = x_std1 - x_std2
        rk_sum = sum(hij * (x_diff ** i) for i, hij in enumerate(self.hij_values))

        return x1 * x2 * rk_sum

    def infinite_dilution_property(self, component: str, temperature: float = DEFAULT_T) -> float:
        """
        计算指定组分在无限稀释时的偏摩尔性质

        Parameters
        ----------
        component : str
            组分标识
        temperature : float
            温度 (K)

        Returns
        -------
        float
            无限稀释偏摩尔性质
        """
        if component not in (self.component1, self.component2):
            raise ValueError(f"组分 '{component}' 不在二元系 {self.name} 中")

        self._update_parameters_for_temperature(temperature)

        if not self.hij_values:
            return 0.0

        # 根据R-K多项式特性计算无限稀释性质
        if component == self.component1:
            if self._is_reversed:
                return sum(self.hij_values)
            else:
                return sum((-1) ** i * hij for i, hij in enumerate(self.hij_values))
        else:
            if self._is_reversed:
                return sum((-1) ** i * hij for i, hij in enumerate(self.hij_values))
            else:
                return sum(self.hij_values)

    def partial_molar_property(self, component: str, x_component: float,
                              temperature: float = DEFAULT_T) -> float:
        """
        计算指定组分的偏摩尔性质

        Parameters
        ----------
        component : str
            组分标识
        x_component : float
            该组分的摩尔分数
        temperature : float
            温度 (K)

        Returns
        -------
        float
            偏摩尔性质
        """
        if component not in (self.component1, self.component2):
            raise ValueError(f"组分 '{component}' 不在二元系 {self.name} 中")

        if not 0 <= x_component <= 1:
            raise ValueError("摩尔分数必须在0到1之间")

        self._update_parameters_for_temperature(temperature)

        if not self.hij_values:
            return 0.0

        # 设置摩尔分数
        if component == self.component1:
            x1, x2 = x_component, 1 - x_component
        else:
            x1, x2 = 1 - x_component, x_component

        # 计算标准化顺序下的摩尔分数和差值
        if self._is_reversed:
            x_std1, x_std2 = x2, x1
        else:
            x_std1, x_std2 = x1, x2

        x_diff = x_std1 - x_std2

        # 计算R-K多项式的各项和
        sum0 = sum(hij * (x_diff ** k) for k, hij in enumerate(self.hij_values))
        sum1 = sum(hij * (x_diff ** (k + 1)) for k, hij in enumerate(self.hij_values))
        sum2 = sum(hij * (x_diff ** (k - 1)) for k, hij in enumerate(self.hij_values) if k >= 1)

        # 根据组分计算偏摩尔性质
        if component == self.component1:
            if self._is_reversed:
                result = x1 * x2 * sum0 + (1 - x2) * (-2 * x1 * x2 * sum2 + sum1)
            else:
                result = x1 * x2 * sum0 + (1 - x1) * (2 * x1 * x2 * sum2 - sum1)
        else:
            if self._is_reversed:
                result = x1 * x2 * sum0 + (1 - x1) * (2 * x1 * x2 * sum2 - sum1)
            else:
                result = x1 * x2 * sum0 + (1 - x2) * (-2 * x1 * x2 * sum2 + sum1)

        return result

    def get_component_info(self) -> dict:
        """获取组分信息"""
        return {
            'original_composition': (self.component1, self.component2),
            'standardized_composition': self._standardized_components,
            'system_id': self.system_id,
            'name': self.name,
            'is_reversed': self._is_reversed
        }

    def __str__(self) -> str:
        return f"RKBinaryPolynomial({self.name})"

    def __repr__(self) -> str:
        return f"RKBinaryPolynomial(composition=({self.component1!r}, {self.component2!r}))"
