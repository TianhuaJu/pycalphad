# Pycalphad数据库机制与多数据库使用指南

## 📚 核心概念

### 1. 参考态与模型的关系

> **重要结论：参考态是模型无关的**

**为什么？**
- 参考态是纯组元的化学势（例如纯Fe、纯Ni）
- 纯组元性质由数据库中的纯组元参数决定（如GHSERFE、GHSERNI等）
- **模型（RKM、UEM等）只影响多元系的交互项**
- 模型不改变纯组元的本征性质

**代码体现：**
```python
# ❌ 错误理解（已修正）
ref_eq = equilibrium(dbe, [comp, 'VA'], phase,
                    conditions={...},
                    model=ModelUEM1)  # 不应该传递模型！

# ✓ 正确理解
ref_eq = equilibrium(dbe, [comp, 'VA'], phase,
                    conditions={...})
                    # 参考态不传model参数
```

### 2. 活度计算原理

活度公式：
```
a_i = exp((μ_i - μ_i^ref) / RT)
```

其中：
- `μ_i`: 混合态中组分i的化学势（受模型影响）
- `μ_i^ref`: 纯组元i的参考态化学势（与模型无关）
- `R`: 气体常数 8.315 J/(mol·K)
- `T`: 温度 (K)

**关键点：**
- **混合态**计算需要使用自定义模型（UEM等）
- **参考态**计算不使用模型（纯组元性质）
- 两者的差值反映了组元在混合物中的偏离程度

## 🗄️ Database类结构

### 核心数据存储

```python
from pycalphad import Database

dbe = Database('system.tdb')

# 数据存储位置：
dbe.symbols      # dict: FUNCTION定义 {'GHSERFE': <SymEngine表达式>}
dbe._parameters  # TinyDB: 所有参数 (G, L, TC, BMAGN等)
dbe.elements     # set: 元素名称
dbe.species      # set: Species对象
dbe.phases       # dict: Phase对象
dbe.refstates    # dict: 参考态定义
```

### 参数查询机制

```python
from tinydb import where

# 查询纯组元G参数
pure_params = dbe.search(
    (where('parameter_type') == 'G') &
    (where('parameter_order') == 0) &
    (where('phase_name') == 'LIQUID')
)

# 查询交互参数
interaction_params = dbe.search(
    (where('parameter_type') == 'L') |
    (where('parameter_order') > 0)
)
```

## 🔀 多数据库使用策略

### 当前GUI实现（字符串连接法）

**位置：** `alloy_calculator_gui.py` 第538-558行

```python
def load_database(self):
    # 用户可以选择多个TDB文件
    file_paths = filedialog.askopenfilename(multiple=True)

    # 连接所有文件内容
    all_tdb_content = ""
    for path in file_paths:
        with open(path, 'r', encoding='latin-1') as f:
            all_tdb_content += f.read() + "\n"

    # 创建单一数据库对象
    self.dbe = Database(all_tdb_content)
```

**优点：**
- 简单直接
- pycalphad原生支持
- 自动处理参数合并

**缺点：**
- 如果有重复FUNCTION定义会报错
- 无法精确控制参数优先级
- 后加载的参数会覆盖先加载的参数

### 高级用法：选择性数据库合并

#### 应用场景

> **用户需求：**
> "能否自定义纯组元性质的读取数据库，比如unary50.tdb读取纯组元的数据，交互数据或其它化合相从其他数据库读取。这样做的目的是将来想扩展pycalphad，计算更多组元（比如四元、五元等）热力学性质时方便。"

#### 解决方案：手动数据库合并

```python
from pycalphad import Database
from copy import deepcopy
from tinydb import where

def create_merged_database(unary_tdb_path, system_tdb_path):
    """
    创建合并数据库：纯组元数据来自unary50.tdb，交互参数来自系统数据库

    参数：
        unary_tdb_path: 纯组元数据库路径（如unary50.tdb）
        system_tdb_path: 系统数据库路径（如Al-Cr-Ni-Fe.tdb）

    返回：
        merged_db: 合并后的数据库对象
    """
    # 1. 加载两个数据库
    unary_db = Database(unary_tdb_path)
    system_db = Database(system_tdb_path)

    # 2. 以系统数据库为基础
    merged_db = deepcopy(system_db)

    # 3. 覆盖纯组元FUNCTION（GHSER等）
    for symbol_name, symbol_expr in unary_db.symbols.items():
        if symbol_name.startswith('GHSER'):
            merged_db.symbols[symbol_name] = symbol_expr
            print(f"✓ 覆盖纯组元函数: {symbol_name}")

    # 4. 更新参考态定义
    for elem, refstate in unary_db.refstates.items():
        merged_db.refstates[elem] = refstate

    # 5. 选择性合并纯组元G参数
    def is_pure_endmember(param):
        """判断是否为纯组元端元参数"""
        return (param['parameter_type'] == 'G' and
                param['parameter_order'] == 0 and
                all(len(subl) == 1 for subl in param['constituent_array']))

    # 移除系统数据库中的纯组元参数
    merged_db._parameters.remove(where('parameter_order') == 0)

    # 插入unary数据库中的纯组元参数
    pure_params = [p for p in unary_db._parameters.all()
                   if is_pure_endmember(p)]
    merged_db._parameters.insert_multiple(pure_params)

    print(f"✓ 添加了 {len(pure_params)} 个纯组元参数")

    # 6. 合并元素和物种
    merged_db.elements.update(unary_db.elements)
    merged_db.species.update(unary_db.species)

    return merged_db

# 使用示例
merged_db = create_merged_database(
    'databases/unary50.tdb',
    'databases/Al-Cr-Ni-Fe-Co.tdb'
)

# 用于计算
from pycalphad import equilibrium
import pycalphad.variables as v

result = equilibrium(
    merged_db,
    ['AL', 'CR', 'NI', 'FE', 'CO'],
    ['FCC_A1', 'BCC_A2', 'LIQUID'],
    {v.T: 1200, v.P: 101325,
     v.X('AL'): 0.2, v.X('CR'): 0.2,
     v.X('NI'): 0.2, v.X('FE'): 0.2}
)
```

### 参数优先级控制

```python
def merge_with_priority(primary_db, secondary_db, override_symbols=True):
    """
    合并数据库并控制优先级

    参数：
        primary_db: 主数据库（高优先级）
        secondary_db: 次数据库（低优先级）
        override_symbols: 是否用primary覆盖secondary的符号

    返回：
        merged_db: 合并后的数据库
    """
    from copy import deepcopy

    merged_db = deepcopy(secondary_db)

    # 符号覆盖（如果启用）
    if override_symbols:
        for name, expr in primary_db.symbols.items():
            merged_db.symbols[name] = expr
    else:
        # 只添加不存在的符号
        for name, expr in primary_db.symbols.items():
            if name not in merged_db.symbols:
                merged_db.symbols[name] = expr

    # 参数合并（检查是否存在相同参数）
    from tinydb import where

    for param in primary_db._parameters.all():
        # 检查是否已存在相同参数
        existing = merged_db._parameters.search(
            (where('phase_name') == param['phase_name']) &
            (where('constituent_array') == param['constituent_array']) &
            (where('parameter_type') == param['parameter_type']) &
            (where('parameter_order') == param['parameter_order'])
        )

        if existing:
            # 移除旧参数，插入新参数
            merged_db._parameters.remove(
                (where('phase_name') == param['phase_name']) &
                (where('constituent_array') == param['constituent_array']) &
                (where('parameter_type') == param['parameter_type']) &
                (where('parameter_order') == param['parameter_order'])
            )

        merged_db._parameters.insert(param)

    # 合并其他属性
    merged_db.elements.update(primary_db.elements)
    merged_db.species.update(primary_db.species)
    merged_db.refstates.update(primary_db.refstates)

    for phase_name, phase in primary_db.phases.items():
        if phase_name not in merged_db.phases:
            merged_db.phases[phase_name] = phase

    return merged_db
```

## 🔍 纯组元性质读取机制

### Model.reference_energy() 方法

**位置：** `pycalphad/model.py` 第866-881行

```python
def reference_energy(self, dbe):
    """
    返回端元能量的加权平均（符号形式）
    """
    # 查询纯组元G参数
    pure_param_query = (
        (where('phase_name') == self.phase_name) &
        (where('parameter_order') == 0) &
        (where('parameter_type') == "G") &
        (where('constituent_array').test(self._purity_test))
    )

    phase = dbe.phases[self.phase_name]
    param_search = dbe.search

    # 使用Redlich-Kister求和
    pure_energy_term = self.redlich_kister_sum(
        phase, param_search, pure_param_query
    )

    return pure_energy_term / self._site_ratio_normalization
```

**关键点：**
1. 只查询 `parameter_order == 0` 的参数（端元参数）
2. 使用 `_purity_test` 过滤纯组元（每个亚晶格只有一个组分）
3. **不涉及任何模型特定的计算**
4. 结果仅依赖于数据库中的纯组元参数

### 符号替换机制

```python
# TDB文件中的定义
FUNCTION GHSERFE  298.15
  +1225.7+124.134*T-23.5143*T*LN(T)-.00439752*T**2-5.8927E-08*T**3
  +77359*T**(-1); 1811.00 Y
  -25383.581+299.31255*T-46*T*LN(T)+2.29603E+31*T**(-9); 6000.00 N !

# 数据库加载后
dbe.symbols['GHSERFE'] = <SymEngine表达式>

# 在Model中使用
# PARAMETER G(BCC_A2,FE:VA;0) 298.15 GHSERFE; 6000 N !
# 计算时会自动替换：
G = GHSERFE  # 被替换为实际的数学表达式
```

## 💡 最佳实践

### 1. 单一系统计算

如果只计算固定的系统（如Al-Cr-Ni），**使用字符串连接法**：

```python
# 简单且有效
dbe = Database('alcrni.tdb')
```

### 2. 多系统复用纯组元数据

如果要计算多个系统（四元、五元等），**使用数据库合并法**：

```python
# 加载一次unary数据库
unary_db = Database('unary50.tdb')

# 为不同系统创建合并数据库
alcrni_db = create_merged_database(unary_db, 'AlCrNi.tdb')
alcrnife_db = create_merged_database(unary_db, 'AlCrNiFe.tdb')
alcrnifeco_db = create_merged_database(unary_db, 'AlCrNiFeCo.tdb')
```

### 3. GUI集成建议

```python
class AlloyCalculatorGUI:
    def __init__(self):
        self.unary_db = None  # 可选的纯组元数据库
        self.system_db = None  # 系统数据库
        self.dbe = None  # 最终使用的数据库

    def load_unary_database(self):
        """加载纯组元数据库（可选）"""
        path = filedialog.askopenfilename(
            title="选择纯组元数据库 (可选，如unary50.tdb)"
        )
        if path:
            self.unary_db = Database(path)
            self.log("✓ 已加载纯组元数据库")

    def load_system_database(self):
        """加载系统数据库"""
        paths = filedialog.askopenfilename(
            title="选择系统数据库文件",
            multiple=True
        )

        if paths:
            # 合并系统数据库文件
            all_content = ""
            for path in paths:
                with open(path, 'r', encoding='latin-1') as f:
                    all_content += f.read() + "\n"

            self.system_db = Database(all_content)

            # 如果有纯组元数据库，合并它们
            if self.unary_db:
                self.dbe = merge_databases(self.unary_db, self.system_db)
                self.log("✓ 已合并纯组元数据库和系统数据库")
            else:
                self.dbe = self.system_db
                self.log("✓ 已加载系统数据库")
```

## 📊 验证数据库合并

```python
def verify_database_merge(merged_db, unary_db, system_db):
    """验证数据库合并是否成功"""

    print("=== 数据库合并验证 ===\n")

    # 1. 检查纯组元函数
    print("1. 纯组元FUNCTION检查:")
    ghser_funcs = [s for s in merged_db.symbols.keys()
                   if s.startswith('GHSER')]
    print(f"   找到 {len(ghser_funcs)} 个GHSER函数")
    for func in ghser_funcs[:5]:  # 显示前5个
        print(f"   - {func}")

    # 2. 检查元素
    print(f"\n2. 元素列表:")
    print(f"   unary数据库: {len(unary_db.elements)} 个元素")
    print(f"   系统数据库: {len(system_db.elements)} 个元素")
    print(f"   合并数据库: {len(merged_db.elements)} 个元素")

    # 3. 检查相
    print(f"\n3. 相列表:")
    print(f"   unary数据库: {len(unary_db.phases)} 个相")
    print(f"   系统数据库: {len(system_db.phases)} 个相")
    print(f"   合并数据库: {len(merged_db.phases)} 个相")

    # 4. 检查参数数量
    from tinydb import where

    pure_params = merged_db.search(
        (where('parameter_order') == 0) &
        (where('parameter_type') == 'G')
    )
    interaction_params = merged_db.search(
        (where('parameter_order') > 0) |
        (where('parameter_type') == 'L')
    )

    print(f"\n4. 参数统计:")
    print(f"   纯组元参数: {len(pure_params)}")
    print(f"   交互参数: {len(interaction_params)}")

    # 5. 测试计算
    print(f"\n5. 测试计算:")
    try:
        from pycalphad import equilibrium
        import pycalphad.variables as v

        # 选择系统中的前两个元素进行测试
        test_comps = list(merged_db.elements)[:2] + ['VA']
        test_phases = list(merged_db.phases.keys())[:3]

        result = equilibrium(
            merged_db, test_comps, test_phases,
            {v.T: 1000, v.P: 101325, v.X(test_comps[0]): 0.5}
        )
        print(f"   ✓ 成功计算 {test_comps[0]}-{test_comps[1]} 系统")
    except Exception as e:
        print(f"   ✗ 计算失败: {e}")

# 使用示例
verify_database_merge(merged_db, unary_db, system_db)
```

## 🚀 扩展到多元系统

### 四元系统示例

```python
# 加载数据库
merged_db = create_merged_database(
    'unary50.tdb',
    'Al-Cr-Ni-Fe.tdb'
)

# 计算四元相图截面
from pycalphad import equilibrium
import pycalphad.variables as v
import numpy as np

# 固定Fe和Ni的比例，扫描Al和Cr
results = []

for x_al in np.linspace(0, 1, 50):
    for x_cr in np.linspace(0, 1-x_al, 50):
        x_fe = (1 - x_al - x_cr) * 0.5
        x_ni = (1 - x_al - x_cr) * 0.5

        eq = equilibrium(
            merged_db,
            ['AL', 'CR', 'FE', 'NI', 'VA'],
            ['FCC_A1', 'BCC_A2', 'LIQUID'],
            {v.T: 1200, v.P: 101325,
             v.X('AL'): x_al,
             v.X('CR'): x_cr,
             v.X('FE'): x_fe}
        )
        results.append(eq)
```

### 五元系统示例

```python
# 高熵合金计算
merged_db = create_merged_database(
    'unary50.tdb',
    'AlCrFeCoNi.tdb'
)

# 等摩尔比高熵合金
eq = equilibrium(
    merged_db,
    ['AL', 'CR', 'FE', 'CO', 'NI', 'VA'],
    ['FCC_A1', 'BCC_A2'],
    {v.T: (300, 2000, 10),  # 温度扫描
     v.P: 101325,
     v.X('AL'): 0.2,
     v.X('CR'): 0.2,
     v.X('FE'): 0.2,
     v.X('CO'): 0.2}
    # NI自动为1 - 0.2*4 = 0.2
)
```

## 📚 参考资料

### 核心文件位置

- **Database类**: `pycalphad/io/database.py`
- **TDB解析**: `pycalphad/io/tdb.py`
- **Model类**: `pycalphad/model.py`
- **参考态处理**: `pycalphad/model.py` 第866-881行

### 相关文档

- pycalphad官方文档: https://pycalphad.org/docs/latest/
- Database API: https://pycalphad.org/docs/latest/api/pycalphad.io.database.html
- Model API: https://pycalphad.org/docs/latest/api/pycalphad.model.html

## 📝 总结

### 关键认识

1. **参考态是模型无关的** - 纯组元性质不受模型影响
2. **模型只影响交互项** - UEM、RKM等模型只改变多元系混合能
3. **数据库可以合并** - 通过字符串连接或手动合并
4. **优先级可控** - 可以选择性地覆盖参数

### 实现建议

1. **简单系统** - 使用字符串连接法加载多个TDB文件
2. **复杂系统** - 使用数据库合并法，分离纯组元数据和交互数据
3. **多元扩展** - 使用unary50.tdb作为纯组元数据的统一来源
4. **验证重要** - 始终验证合并后的数据库是否正确

通过正确理解参考态和模型的关系，以及灵活运用数据库合并技术，可以高效地进行多元系统的热力学计算！
