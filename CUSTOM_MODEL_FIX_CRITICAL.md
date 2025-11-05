# 自定义模型传递机制 - 关键修正文档

## ⚠️ 核心问题

### 用户反馈
> "本项目开发的目的是采用自定义的UEM模型替代pycalphad中默认的模型计算液相（或固溶相）的热力学性质，因此，在计算相图时，**不能完全退回到默认模型**，而一定要保证，计算时指定相是采用了用户指定的模型。"

### 问题发现
在之前的实现中，`_calculate_using_binplot()` 方法**没有传递自定义模型参数**，导致：
- 官方binplot使用默认RKM模型
- **用户选择的UEM等自定义模型被忽略**
- **违背了项目的核心目的**

## ✅ 解决方案

### 关键修改

#### 1. 修改 `_calculate_using_binplot()` 签名
```python
# 旧版（错误）
def _calculate_using_binplot(self, model_key, inputs, model_label):
    # 没有model_spec参数！

# 新版（正确）
def _calculate_using_binplot(self, model_key, inputs, model_label, model_spec):
    # ✓ 接收model_spec参数
```

#### 2. 通过eq_kwargs传递模型
```python
# 构建eq_kwargs以传递自定义模型
eq_kwargs = {}
if model_spec is not None:
    eq_kwargs['model'] = model_spec
    self.log(f"✓ 使用自定义模型: {model_label}")
else:
    self.log("使用默认RKM模型")

# 使用binplot绘制，关键：通过eq_kwargs传递模型
ax = binplot(
    self.dbe,
    inputs['study_comps'],
    inputs['db_phases'],
    conditions,
    eq_kwargs=eq_kwargs,  # ⭐关键：传递模型参数到equilibrium
    plot_kwargs={'ax': self.fig_phase.gca()}
)
```

#### 3. 更新调用逻辑
```python
# 旧版（错误）- 只有RKM才用binplot
if is_true_binary and model_key == 'RKM':
    self._calculate_using_binplot(model_key, inputs, model_label)  # 缺少model_spec!

# 新版（正确）- 所有真二元都可用binplot，并传递模型
if is_true_binary:
    self._calculate_using_binplot(model_key, inputs, model_label, model_spec)
```

## 🔍 模型传递路径

### 完整调用链

```
1. 用户界面选择
   └─> 模型: UEM

2. get_model_spec('UEM1')
   └─> 返回: ModelUEM1 类

3. _calculate_phase_diagram_thread(..., model_spec=ModelUEM1)
   └─> 检测系统类型

4a. 真二元路径 (binplot)
    └─> _calculate_using_binplot(..., model_spec=ModelUEM1)
        └─> eq_kwargs = {'model': ModelUEM1}
            └─> binplot(..., eq_kwargs=eq_kwargs)
                └─> map_binary(..., eq_kwargs)
                    └─> equilibrium(..., model=ModelUEM1)
                        └─> ✓ 使用UEM模型计算！

4b. 伪二元路径 (equilibrium)
    └─> _calculate_using_equilibrium(..., model_spec=ModelUEM1)
        └─> equilibrium(..., model=ModelUEM1)
            └─> ✓ 使用UEM模型计算！
```

## 📚 pycalphad模型系统

### 模型参数传递方式

#### 方式1: 直接传递 (equilibrium)
```python
eq = equilibrium(
    dbe, components, phases, conditions,
    model=ModelUEM1  # 直接传递模型类
)
```

#### 方式2: 通过eq_kwargs (binplot/ternplot)
```python
binplot(
    dbe, components, phases, conditions,
    eq_kwargs={'model': ModelUEM1}  # 间接传递
)
# binplot内部调用equilibrium时会使用该模型
```

#### 方式3: 字典方式（按相指定）
```python
# 只对液相使用UEM，其他用默认
model_dict = {
    'LIQUID': ModelUEM1,
    'FCC_A1': Model,  # 默认
    'BCC_A2': Model
}

eq = equilibrium(
    dbe, components, phases, conditions,
    model=model_dict
)
```

### get_model_spec() 返回值

```python
def get_model_spec(self, model_key):
    if model_key == 'RKM':
        return None  # 使用默认

    elif model_key == 'UEM1':
        if self.uem1_liquid_only.get():
            # 只替换液相
            return {
                'LIQUID': ModelUEM1,
                'FCC_A1': Model,
                ...
            }
        else:
            # 所有相都用UEM
            return ModelUEM1

    elif model_key in ['Muggianu', 'Toop']:
        return ModelMuggianu / ModelToop

    else:
        return None
```

## ✓ 验证方法

### 1. 日志验证
查看GUI日志输出：
```
[相图计算] 模型: UEM
使用官方binplot方法计算真二元相图
✓ 使用自定义模型: UEM
✓ 二元相图绘制完成（官方binplot+自定义模型）
```

### 2. 代码追踪
在关键位置添加断点：
```python
# pycalphad/plot/binary/map.py:68
models = eq_kwargs.get('model')
print(f"模型: {models}")  # 应该显示 ModelUEM1

# pycalphad/core/equilibrium.py
# 查看model参数是否被使用
```

### 3. 结果验证
- UEM模型的热力学性质应与RKM有区别
- 相图边界位置应反映UEM的计算结果
- 对比RKM和UEM的相图应有明显差异

## 🎯 所有计算路径的模型传递检查表

### ✅ 已验证的计算路径

| 功能 | 方法 | 模型传递 | 状态 |
|------|------|----------|------|
| 液相线计算 | `calculate_liquidus_thread()` | `model=model_spec` | ✅ 正确 |
| Gibbs自由能 | `calculate_gibbs_thread()` | `model=model_spec` | ✅ 正确 |
| 活度计算 | `calculate_activity_thread()` | `model=model_spec` | ✅ 正确 |
| **二元相图(binplot)** | `_calculate_using_binplot()` | `eq_kwargs={'model': model_spec}` | ✅ **已修正** |
| 伪二元相图(equilibrium) | `_calculate_using_equilibrium()` | `model=model_spec` | ✅ 正确 |

### 未来需要的计算路径

| 功能 | 方法 | 模型传递方式 |
|------|------|--------------|
| 三元相图 | `ternplot()` | `eq_kwargs={'model': model_spec}` |
| 等温截面 | `equilibrium()` | `model=model_spec` |
| 步冷曲线 | `equilibrium()` | `model=model_spec` |

## 🔬 测试用例

### 测试用例1: 真二元 + UEM
```python
系统: Al-Ni
模型: UEM1
温度: 300-2000K
成分: X(NI) = 0-1

预期行为:
1. 检测为真二元 ✓
2. 使用 _calculate_using_binplot() ✓
3. 日志显示"使用自定义模型: UEM" ✓
4. binplot调用时传递 eq_kwargs={'model': ModelUEM1} ✓
5. 相图反映UEM热力学性质 ✓
```

### 测试用例2: 伪二元 + UEM
```python
系统: Al-Cr-Ni (固定Al:Cr=1:1)
模型: UEM1
温度: 800-1800K
成分: X(NI) = 0-1

预期行为:
1. 检测为伪二元 ✓
2. 使用 _calculate_using_equilibrium() ✓
3. 每次equilibrium调用传递 model=ModelUEM1 ✓
4. 相图反映UEM热力学性质 ✓
```

### 测试用例3: 真二元 + RKM
```python
系统: Al-Ni
模型: RKM (默认)
温度: 300-2000K
成分: X(NI) = 0-1

预期行为:
1. 检测为真二元 ✓
2. 使用 _calculate_using_binplot() ✓
3. model_spec = None ✓
4. eq_kwargs = {} (空字典) ✓
5. binplot使用默认模型 ✓
```

## 📖 关键代码位置

### 主要文件
- `alloy_calculator_gui.py:1697-1755` - _calculate_using_binplot()
- `alloy_calculator_gui.py:1757-1847` - _calculate_using_equilibrium()
- `alloy_calculator_gui.py:479-504` - get_model_spec()

### pycalphad源码
- `pycalphad/plot/binary/plot.py:75` - binplot()
- `pycalphad/plot/binary/map.py:18` - map_binary()
- `pycalphad/plot/binary/map.py:68` - 读取eq_kwargs['model']
- `pycalphad/core/equilibrium.py` - equilibrium()

## 💡 重要教训

### 1. 理解项目核心目的至关重要
> 如果不理解"用UEM模型替代默认模型"是核心目的，
> 就可能在优化时破坏这个功能。

### 2. 模型参数必须贯穿整个调用链
```
GUI选择 → get_model_spec → 计算方法 → pycalphad API → equilibrium
          ↑___________必须在每一步传递___________↓
```

### 3. 官方wrapper函数的模型传递
- binplot, ternplot等高级函数
- 需要通过 `eq_kwargs` 传递模型
- 不能直接传递 `model` 参数

### 4. 验证日志的重要性
```python
# 添加清晰的日志
self.log(f"✓ 使用自定义模型: {model_label}")
```
让用户能够确认模型确实被使用。

## 🚀 后续建议

### 1. 添加模型使用统计
```python
# 在日志中显示
self.log(f"模型使用统计:")
self.log(f"  - 调用equilibrium次数: {n}")
self.log(f"  - 使用自定义模型: {model_label}")
self.log(f"  - 受影响的相: {phases}")
```

### 2. 添加模型对比功能
```python
# 同时计算RKM和UEM，显示差异
result_rkm = equilibrium(..., model=None)
result_uem = equilibrium(..., model=ModelUEM1)
delta = result_uem.GM - result_rkm.GM
self.log(f"模型差异: ΔG = {delta}")
```

### 3. 支持更多相图类型
- 三元相图 (ternplot)
- 等温截面
- 垂直截面

### 4. 模型验证工具
```python
def verify_model_usage(self, result, expected_model):
    """验证计算结果确实使用了指定模型"""
    # 通过某些特征值判断
    pass
```

## 📝 总结

### 修复前（错误）
```python
# binplot调用没有传递模型
ax = binplot(dbe, comps, phases, conds, plot_kwargs={...})
# ❌ 使用默认RKM模型，UEM被忽略！
```

### 修复后（正确）
```python
# binplot通过eq_kwargs传递模型
eq_kwargs = {'model': ModelUEM1} if model_spec else {}
ax = binplot(dbe, comps, phases, conds,
            eq_kwargs=eq_kwargs,  # ✓ 传递模型
            plot_kwargs={...})
# ✓ UEM模型被正确使用！
```

### 影响范围
- **所有真二元相图计算** - 现在正确使用自定义模型
- **binplot路径** - 从忽略模型到正确传递
- **项目核心功能** - 恢复正常

### 重要性评级
⭐⭐⭐⭐⭐ **最高优先级修复**

这是项目的**核心功能修复**，没有这个修改，项目的主要目的
（使用自定义UEM模型）在真二元相图计算中完全失效。

感谢用户的及时发现和反馈！
