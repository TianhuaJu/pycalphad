# UEMModel.py 修复方案详解

## 问题根源

```python
# ❌ 错误代码（第 70 行）
comps = tuple(sorted([v.Species(c) for c in const_array]))

# ❌ 错误代码（第 108 行）
key = (subl_index, tuple(sorted((comp_k, comp_i))))
```

**为什么这样是错的？**

Redlich-Kister 参数有符号约定：
- **偶数阶**（L0, L2, L4...）：**对称** → L0(A,B) = L0(B,A)
- **奇数阶**（L1, L3, L5...）：**反对称** → L1(A,B) = -L1(B,A)

使用 `sorted()` 会丢失原始顺序，导致奇数阶参数符号错误！

---

## 正确的修复方案（推荐 v2 版本）

### 1. 修复 `_build_binary_L_cache` 函数

```python
def _build_binary_L_cache(self, dbe):
    """
    构建二元 Redlich-Kister 参数的缓存。

    关键修改：记录数据库中的原始顺序（first_comp_in_db）
    """
    self._binary_L_param_cache = {}
    param_search = dbe.search
    phase = dbe.phases[self.phase_name]

    param_query = (
        (where('phase_name') == self.phase_name) &
        ((where('parameter_type') == 'G') | (where('parameter_type') == 'L')) &
        (where('constituent_array').test(self._interaction_test))
    )
    params = param_search(param_query)

    for param in params:
        mixing_subl_indices = []
        is_binary = True
        target_subl = -1
        comps_original = None

        for subl_idx, const_array in enumerate(param['constituent_array']):
            if len(const_array) > 2:
                is_binary = False
                break
            if len(const_array) == 2:
                mixing_subl_indices.append(subl_idx)
                target_subl = subl_idx
                # ✅ 关键修改1：保留原始顺序
                comps_original = tuple([v.Species(c) for c in const_array])
            elif len(const_array) == 1 and const_array[0] == v.Species('*'):
                pass
            elif len(const_array) == 1:
                pass
            else:
                is_binary = False
                break

        if not is_binary or len(mixing_subl_indices) != 1 or comps_original is None:
            continue

        # 使用规范化的 key（sorted）用于查找，但记录原始顺序
        key_canonical = (target_subl, tuple(sorted(comps_original)))

        if key_canonical not in self._binary_L_param_cache:
            self._binary_L_param_cache[key_canonical] = {
                'all': [],
                'odd': [],
                'first_comp_in_db': comps_original[0]  # ✅ 关键修改2：记录第一个组分
            }

        L_expr = param['parameter']
        order = param.get('parameter_order', 0)

        self._binary_L_param_cache[key_canonical]['all'].append((L_expr, order))

        if order % 2 != 0:  # 奇数阶
            self._binary_L_param_cache[key_canonical]['odd'].append(L_expr)
```

**关键修改：**
1. **第 70 行附近**：不再对 `comps` 排序，保留 `comps_original`
2. **新增字段**：在缓存中添加 `'first_comp_in_db': comps_original[0]`

---

### 2. 修复 `_get_uem_d_term` 函数

```python
def _get_uem_d_term(self, dbe, comp_k, comp_i, subl_index):
    """
    根据UEM公式计算 d_ki 项。
    d_ki = (1/RT) * (lim(dG_ki/dx_i) - lim(dG_ki/dx_k))

    对于Redlich-Kister模型，这简化为:
    d_ki = (2/RT) * sum(L_ki^{(v)}) (v 为奇数)

    关键：根据数据库中的顺序判断是否需要符号修正
    """

    # 确保缓存已建立
    if not hasattr(self, '_binary_L_param_cache'):
        self._build_binary_L_cache(dbe)

    # 使用规范化的 key 查找（与缓存构建时一致）
    key_canonical = (subl_index, tuple(sorted((comp_k, comp_i))))

    cache_entry = self._binary_L_param_cache.get(key_canonical, None)
    if cache_entry is None:
        return S.Zero

    odd_L_terms = cache_entry.get('odd', [])
    if not odd_L_terms:
        return S.Zero

    # ✅ 关键修改：根据 first_comp_in_db 判断符号修正
    first_comp_in_db = cache_entry['first_comp_in_db']

    # 判断逻辑：
    # - 如果数据库中是 (k, i) 顺序，我们需要 (k, i)，sign = +1
    # - 如果数据库中是 (i, k) 顺序，我们需要 (k, i)，sign = -1
    if first_comp_in_db == comp_k:
        sign_correction = 1
    elif first_comp_in_db == comp_i:
        sign_correction = -1
    else:
        # 不应该发生这种情况
        raise ValueError(f"Unexpected first component: {first_comp_in_db}")

    # 应用符号修正
    d_ki = (S(2) / (v.R * v.T)) * sign_correction * Add(*odd_L_terms)
    return d_ki
```

**关键修改：**
1. 获取 `first_comp_in_db`
2. 比较 `first_comp_in_db` 与 `comp_k`
3. 如果不同，则需要符号修正（`sign_correction = -1`）

---

## 修复原理图解

### 示例：Al-Cr 二元系统

数据库中定义：
```
L0(Al,Cr) = -30000 J/mol
L1(Al,Cr) = +10000 J/mol
```

### 情况1：调用 `_get_uem_d_term(dbe, Al, Cr, 0)`

```
请求顺序：(Al, Cr)
数据库顺序：(Al, Cr)  ← first_comp_in_db = Al
first_comp_in_db == comp_k  ✅
sign_correction = +1

d_Al_Cr = (2/RT) * (+1) * L1(Al,Cr)
        = (2/RT) * (+10000)
        = +20000/RT
```

### 情况2：调用 `_get_uem_d_term(dbe, Cr, Al, 0)`

```
请求顺序：(Cr, Al)
数据库顺序：(Al, Cr)  ← first_comp_in_db = Al
first_comp_in_db == comp_i  ❌
sign_correction = -1

d_Cr_Al = (2/RT) * (-1) * L1(Al,Cr)
        = (2/RT) * (-10000)
        = -20000/RT
```

**验证正确性：**
```
理论值：L1(Cr,Al) = -L1(Al,Cr) = -10000
d_Cr_Al = (2/RT) * L1(Cr,Al) = (2/RT) * (-10000) ✅
```

---

## 对比：旧代码 vs 新代码

| 方面 | ❌ 旧代码 | ✅ 新代码 |
|------|----------|----------|
| **第70行** | `tuple(sorted(...))` | `tuple([...])` 保留顺序 |
| **第108行** | `tuple(sorted(...))` | `tuple(sorted(...))` + `first_comp_in_db` |
| **缓存内容** | `{'all': [], 'odd': []}` | `{'all': [], 'odd': [], 'first_comp_in_db': ...}` |
| **符号处理** | ❌ 忽略顺序，总是正号 | ✅ 根据顺序动态调整 |
| **二元边界** | ✅ 正确（偶数阶为主） | ✅ 正确 |
| **三元外推** | ❌ 错误（+32351 J/mol） | ✅ 正确（-19395 J/mol） |

---

## 完整的修复代码

完整的修复代码已保存在 `UEMModel_FIXED.py` 文件中。

**应用修复的步骤：**

1. 打开 `pycalphad/UEMModel.py`

2. 替换 `_build_binary_L_cache` 函数（约第 29-93 行）

3. 替换 `_get_uem_d_term` 函数（约第 95-116 行）

4. 运行测试验证：
   ```bash
   python debug_uem_difference.py
   ```

5. 预期结果：
   ```
   UEMModel 过剩能:   -19394.61 J/mol  ✅
   ModelUEM1 过剩能:  -19394.61 J/mol  ✅
   差异:              0.00 J/mol       ✅
   ```

---

## 为什么推荐 v2 版本？

v2 版本相比 v1 版本的优势：

1. **更简洁**：不需要维护双向引用（原始key + 反向key）
2. **更清晰**：逻辑直观，通过 `first_comp_in_db` 一目了然
3. **更高效**：只存储一份数据，通过符号修正实现双向查找
4. **更易维护**：代码量更少，bug 更少

参考 `ModelUEM1` 的实现，它使用的也是类似的思路。

---

## 总结

**核心原则：保留参数的原始定义顺序，在使用时动态调整符号。**

这就像一个坐标变换：
- 数据库给你一个向量 **v** = (L0, L1, L2, L3)
- 如果你需要反方向的向量，不能简单地使用 **v**
- 而应该计算 **v'** = (L0, -L1, L2, -L3)
- 偶数分量不变，奇数分量取反

修复后，UEMModel 将给出与 ModelUEM1 相同的正确结果！
