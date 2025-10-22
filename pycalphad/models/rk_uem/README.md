# R-K-UEM 模型模块

## 概述

R-K-UEM (Redlich-Kister Unified Excess Model) 模块为pycalphad提供了基于Redlich-Kister多项式和统一外推模型(UEM)的多元合金系统热力学计算功能。

### 主要特性

- **R-K多项式二元系统描述**: 支持基于R-K多项式的二元系统热力学性质计算
- **UEM外推方法**: 通过贡献系数统一传统外推模型（Kohler、Muggianu、Toop等）
- **多元系统计算**: 支持任意多元系统的混合焓、过剩Gibbs能、活度系数等性质计算
- **灵活的外推模型**: 提供多种外推模型选择，包括UEM1、UEM2系列、传统模型等

### UEM的优势

UEM（Unified Excess Model）是一种创新的外推方法，它通过计算**贡献系数**来分配第三组分对二元对的影响，从而避免了传统Toop模型中寻找非对称组分的复杂过程。

**核心思想**:
- 对于三元系统i-j-k，计算第三组分k对i-j二元对的贡献系数 α(k→i) 和 α(k→j)
- 贡献系数基于二元系统的性质差异计算
- 统一了Kohler (α=0), Muggianu (α=0.5)等传统模型

## 模块结构

```
pycalphad/models/rk_uem/
├── __init__.py           # 模块初始化，导出主要类和函数
├── rk_binary.py          # R-K二元系多项式类
├── extrapolation.py      # 外推模型函数（UEM, Kohler, Muggianu等）
├── thermodynamics.py     # 多元系热力学计算类
└── README.md            # 本文档
```

## 快速开始

### 1. 二元系统计算

```python
from pycalphad.models.rk_uem import RKBinaryPolynomial

# 创建二元系统对象
binary = RKBinaryPolynomial(("Fe", "Cr"), database_path="path/to/database.db")

# 计算混合焓
H_mix = binary.mixing_enthalpy("Fe", x_component=0.5, temperature=1800.0)
print(f"混合焓: {H_mix} J/mol")

# 计算过剩Gibbs能
G_ex = binary.excess_gibbs_energy("Fe", x_component=0.5, temperature=1800.0)
print(f"过剩Gibbs能: {G_ex} J/mol")

# 计算无限稀释性质
W_inf = binary.infinite_dilution_property("Fe", temperature=1800.0)
print(f"无限稀释性质: {W_inf} J/mol")
```

### 2. 多元系统计算

```python
from pycalphad.models.rk_uem import ThermodynamicCalculator, UEM1

# 创建热力学计算器
calc = ThermodynamicCalculator(database_path="path/to/database.db")

# 定义组成
composition = {
    "Fe": 0.5,
    "Cr": 0.3,
    "Ni": 0.2
}

temperature = 1800.0  # K

# 计算混合焓
H_mix = calc.get_mixing_enthalpy(composition, temperature, UEM1)
print(f"混合焓: {H_mix} J/mol")

# 计算过剩Gibbs能
G_ex = calc.get_excess_gibbs(composition, temperature, UEM1)
print(f"过剩Gibbs能: {G_ex} J/mol")
```

### 3. 活度系数计算

```python
from pycalphad.models.rk_uem import ThermodynamicCalculator, UEM1
import math

calc = ThermodynamicCalculator(database_path="path/to/database.db")

composition = {"Fe": 0.6, "Cr": 0.25, "Ni": 0.15}
temperature = 1800.0
solvent = "Fe"  # 参考组分

# 计算单个组分的活度系数
ln_gamma_Cr = calc.calculate_activity_coefficient(
    composition, solute="Cr", solvent=solvent,
    temperature=temperature, extrapolation_model=UEM1
)
gamma_Cr = math.exp(ln_gamma_Cr)
print(f"Cr的活度系数: γ = {gamma_Cr}")

# 计算所有组分的活度系数
ln_gammas = calc.calculate_all_activity_coefficients(
    composition, solvent, temperature, UEM1
)
for comp, ln_gamma in ln_gammas.items():
    print(f"{comp}: ln(γ) = {ln_gamma:.4f}")

# 计算活度
activity_Cr = calc.calculate_activity(
    composition, "Cr", solvent, temperature, UEM1
)
print(f"Cr的活度: a = {activity_Cr}")
```

### 4. 比较不同外推模型

```python
from pycalphad.models.rk_uem import (
    ThermodynamicCalculator,
    UEM1, UEM2_N, Kohler, Muggianu, Toop_Kohler
)

calc = ThermodynamicCalculator(database_path="path/to/database.db")
composition = {"Fe": 0.5, "Cr": 0.3, "Ni": 0.2}
T = 1800.0

models = [
    ("UEM1", UEM1),
    ("UEM2-N", UEM2_N),
    ("Kohler", Kohler),
    ("Muggianu", Muggianu),
    ("Toop-Kohler", Toop_Kohler),
]

print("不同外推模型的混合焓比较:")
for name, model in models:
    H_mix = calc.get_mixing_enthalpy(composition, T, model)
    print(f"{name:15s}: {H_mix:>10.2f} J/mol")
```

## 外推模型说明

### UEM系列

1. **UEM1** - 非交互作用条件下的UEM
   - 基于无限稀释性质差
   - 计算速度快，适用于弱交互作用系统

2. **UEM2-O** - 考虑交互作用的UEM（原始版本）
   - 基于偏摩尔性质积分
   - 更精确但计算成本较高

3. **UEM2-N** - 基于混合焓面积的UEM
   - 基于归一化的混合焓积分
   - 平衡了精度和效率

4. **UEM2-Adv** - 高级UEM
   - 考虑混合焓曲线的几何形状
   - 最高精度但计算最复杂

### 传统模型

1. **Kohler** - Kohler模型
   - 完全对称外推
   - 贡献系数 α = 0

2. **Muggianu** - Muggianu模型
   - 均等分配外推
   - 贡献系数 α = 0.5

3. **Toop-Kohler** - Toop-Kohler模型
   - 基于非对称组分选择
   - 需要识别非对称组分

4. **Toop-Muggianu** - Toop-Muggianu模型
   - Toop模型的Muggianu变体

5. **GSM** - 几何相似模型
   - 基于二元系曲线的偏差函数

## 数据库要求

### 数据库格式

R-K参数数据库应为SQLite格式，包含以下表:

#### R_K_gE表 (过剩Gibbs能参数)
```sql
CREATE TABLE R_K_gE (
    Symbol TEXT,      -- 系统标识，如 "Fe-Cr"
    A0 TEXT,          -- R-K参数 A0
    A1 TEXT,          -- R-K参数 A1
    A2 TEXT,          -- R-K参数 A2
    A3 TEXT,          -- R-K参数 A3
    A4 TEXT,          -- R-K参数 A4
    A5 TEXT           -- R-K参数 A5
);
```

#### R_K_Hmix表 (混合焓参数)
```sql
CREATE TABLE R_K_Hmix (
    Symbol TEXT,      -- 系统标识，如 "Fe-Cr"
    A0 TEXT,          -- R-K参数 A0
    A1 TEXT,          -- R-K参数 A1
    A2 TEXT,          -- R-K参数 A2
    A3 TEXT,          -- R-K参数 A3
    A4 TEXT,          -- R-K参数 A4
    A5 TEXT           -- R-K参数 A5
);
```

### 参数表达式格式

参数支持以下格式的温度依赖表达式:
- 常数: `"12345.6"`
- 线性: `"12.5*T"`
- 幂次: `"1.5*T**2"`
- 对数: `"10.0*T*lnT"`
- 组合: `"10000+5.2*T"`

示例:
```sql
INSERT INTO R_K_Hmix VALUES (
    'Fe-Cr',
    '-17737+7.996546*T',  -- A0
    '1331',               -- A1
    NULL,                 -- A2
    NULL,                 -- A3
    NULL,                 -- A4
    NULL                  -- A5
);
```

## API 参考

### RKBinaryPolynomial

二元系统R-K多项式类。

**构造函数**:
```python
RKBinaryPolynomial(composition: Tuple[str, str], database_path: str = None)
```

**主要方法**:
- `mixing_enthalpy(component, x_component, temperature)` - 计算混合焓
- `excess_gibbs_energy(component, x_component, temperature)` - 计算过剩Gibbs能
- `partial_molar_property(component, x_component, temperature)` - 计算偏摩尔性质
- `infinite_dilution_property(component, temperature)` - 计算无限稀释性质

### ThermodynamicCalculator

多元系统热力学计算器。

**构造函数**:
```python
ThermodynamicCalculator(database_path: str)
```

**主要方法**:
- `get_mixing_enthalpy(composition, temperature, extrapolation_model)` - 计算混合焓
- `get_excess_gibbs(composition, temperature, extrapolation_model)` - 计算过剩Gibbs能
- `calculate_activity_coefficient(composition, solute, solvent, temperature, extrapolation_model)` - 计算活度系数
- `calculate_activity(composition, solute, solvent, temperature, extrapolation_model)` - 计算活度
- `calculate_all_activity_coefficients(composition, solvent, temperature, extrapolation_model)` - 计算所有组分的活度系数

### 外推模型函数

所有外推模型函数具有统一的接口:
```python
def extrapolation_model(k: str, i: str, j: str, T: float, database_path: str) -> float:
    """
    计算组分k对i-j二元对的贡献系数

    Parameters
    ----------
    k : str
        第三组分
    i, j : str
        二元对的两个组分
    T : float
        温度 (K)
    database_path : str
        数据库路径

    Returns
    -------
    float
        贡献系数 (0到1之间)
    """
```

**可用的外推模型**:
- `UEM1`, `UEM2_O`, `UEM2_N`, `UEM2_Adv`
- `Kohler`, `Muggianu`
- `Toop_Kohler`, `Toop_Muggianu`
- `GSM`

## 理论背景

### R-K多项式

对于二元系统i-j，性质F表示为:

$$F = x_i x_j \sum_{k=0}^{n} A_k (x_i - x_j)^k$$

其中:
- $x_i, x_j$ 是摩尔分数
- $A_k$ 是R-K参数
- $n$ 是多项式阶数（通常为5）

### UEM外推

对于三元系统i-j-k，i-j二元对的有效摩尔分数为:

$$X_i^{eff} = x_i + \alpha_{ki}^{ij} x_k$$
$$X_j^{eff} = x_j + \alpha_{kj}^{ij} x_k$$

其中 $\alpha_{ki}^{ij}$ 是贡献系数，计算公式为:

$$\alpha_{ki}^{ij} = \frac{\Delta_{kj}}{\Delta_{ki} + \Delta_{kj}} \exp(-\Delta_{ki})$$

$\Delta$ 表示组分间的性质差异。

## 运行示例

```bash
# 查看API用法演示
python examples/rk_uem_example.py

# 使用实际数据库运行示例
python examples/rk_uem_example.py /path/to/database.db
```

## 注意事项

1. **数据库路径**: 所有计算都需要提供R-K参数数据库
2. **组成归一化**: 模块会自动归一化输入的组成
3. **温度单位**: 所有温度单位为开尔文(K)
4. **能量单位**: 所有能量单位为焦耳每摩尔(J/mol)
5. **数值稳定性**: 模块内置了数值保护机制，避免除零和溢出

## 参考文献

1. Zhang et al., "General formalism for new generation geometrical model: Application to the thermodynamics of liquid mixtures", Calphad, 2010
2. UEM模型相关文献（请补充）

## 许可证

遵循pycalphad项目的许可证。

## 贡献

欢迎贡献代码和报告问题！

## 联系方式

如有问题或建议，请通过pycalphad项目的GitHub仓库联系。
