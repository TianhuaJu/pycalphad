# UEMModel.py 的修复方案
#
# 问题：使用 sorted() 导致丢失参数顺序，造成奇数阶 R-K 参数符号错误
# 解决方案：保留原始顺序，并根据调用顺序动态调整符号

def _build_binary_L_cache(self, dbe):
    """
    构建二元 Redlich-Kister 参数的缓存。

    修复：不再使用 sorted()，而是保留数据库中的原始顺序。
    对于每个二元交互，在缓存中同时存储两个顺序的 key。
    """
    self._binary_L_param_cache = {}
    param_search = dbe.search
    phase = dbe.phases[self.phase_name]

    # 查询所有 'G' 或 'L' 交互参数
    param_query = (
        (where('phase_name') == self.phase_name) &
        ((where('parameter_type') == 'G') | (where('parameter_type') == 'L')) &
        (where('constituent_array').test(self._interaction_test))
    )
    params = param_search(param_query)

    for param in params:
        # 检查这是否是一个单亚晶格二元交互
        mixing_subl_indices = []
        is_binary = True
        target_subl = -1
        comps = []

        for subl_idx, const_array in enumerate(param['constituent_array']):
            if len(const_array) > 2:  # 超过2元，不是二元交互
                is_binary = False
                break
            if len(const_array) == 2:
                mixing_subl_indices.append(subl_idx)
                target_subl = subl_idx
                # ❌ 旧代码：comps = tuple(sorted([v.Species(c) for c in const_array]))
                # ✅ 新代码：保留原始顺序
                comps = tuple([v.Species(c) for c in const_array])
            elif len(const_array) == 1 and const_array[0] == v.Species('*'):
                pass  # 允许通配符
            elif len(const_array) == 1:
                pass  # 允许单个固定组元
            else:
                is_binary = False
                break

        # 我们只关心在单个亚晶格上混合的二元交互
        if not is_binary or len(mixing_subl_indices) != 1:
            continue

        # 保存原始顺序的 key
        key_original = (target_subl, comps)
        # 同时保存反向顺序的 key
        key_reversed = (target_subl, (comps[1], comps[0]))

        L_expr = param['parameter']
        order = param.get('parameter_order', 0)

        # 为原始顺序存储参数
        if key_original not in self._binary_L_param_cache:
            self._binary_L_param_cache[key_original] = {
                'all': [],
                'odd': [],
                'first_comp': comps[0]  # 记录第一个组分
            }

        self._binary_L_param_cache[key_original]['all'].append((L_expr, order))

        if order % 2 != 0:  # 奇数阶
            self._binary_L_param_cache[key_original]['odd'].append(L_expr)

        # 为反向顺序创建引用（指向原始顺序）
        if key_reversed not in self._binary_L_param_cache:
            self._binary_L_param_cache[key_reversed] = {
                'all': [],
                'odd': [],
                'first_comp': comps[0],  # 仍然记录原始第一个组分
                'is_reversed': True,  # 标记这是反向引用
                'original_key': key_original  # 指向原始 key
            }


def _get_uem_d_term(self, dbe, comp_k, comp_i, subl_index):
    """
    根据UEM公式计算 d_ki 项。
    d_ki = (1/RT) * (lim(dG_ki/dx_i) - lim(dG_ki/dx_k))

    对于Redlich-Kister模型，这简化为:
    d_ki = (2/RT) * sum(L_ki^{(v)}) (v 为奇数)

    修复：考虑参数定义顺序，对奇数阶参数进行符号修正。

    关键：L1(A,B) = -L1(B,A)，所以需要根据调用顺序判断符号
    """

    # 确保缓存已建立
    if not hasattr(self, '_binary_L_param_cache'):
        self._build_binary_L_cache(dbe)

    # ❌ 旧代码：key = (subl_index, tuple(sorted((comp_k, comp_i))))
    # ✅ 新代码：尝试查找两个可能的顺序

    key_requested = (subl_index, (comp_k, comp_i))
    key_reversed = (subl_index, (comp_i, comp_k))

    cache_entry = None
    sign_correction = 1  # 默认不需要符号修正

    if key_requested in self._binary_L_param_cache:
        cache_entry = self._binary_L_param_cache[key_requested]
        # 检查是否是反向引用
        if cache_entry.get('is_reversed', False):
            # 我们请求的是 (k, i)，但数据库中是 (i, k)
            sign_correction = -1
            # 从原始 key 获取实际的奇数阶参数
            original_key = cache_entry['original_key']
            odd_L_terms = self._binary_L_param_cache[original_key].get('odd', [])
        else:
            # 顺序一致
            odd_L_terms = cache_entry.get('odd', [])
    elif key_reversed in self._binary_L_param_cache:
        cache_entry = self._binary_L_param_cache[key_reversed]
        # 我们请求的是 (k, i)，但找到的是 (i, k)
        if cache_entry.get('is_reversed', False):
            # 缓存中的 (i, k) 本身是反向引用，实际原始是 (k, i)
            sign_correction = 1
            original_key = cache_entry['original_key']
            odd_L_terms = self._binary_L_param_cache[original_key].get('odd', [])
        else:
            # 缓存中的 (i, k) 是原始顺序，我们需要 (k, i)，所以要取反
            sign_correction = -1
            odd_L_terms = cache_entry.get('odd', [])
    else:
        # 找不到任何参数
        return S.Zero

    if not odd_L_terms:
        return S.Zero

    # 应用符号修正
    # d_ki = (2/RT) * sign * sum(L_ki^{(odd)})
    d_ki = (S(2) / (v.R * v.T)) * sign_correction * Add(*odd_L_terms)
    return d_ki


# ============================================================================
# 更简洁的替代方案（推荐）
# ============================================================================
# 基于 ModelUEM1 的实现思路，更简单直接

def _build_binary_L_cache_v2(self, dbe):
    """
    构建二元 Redlich-Kister 参数的缓存（简化版本）。

    策略：为每个二元对，无论顺序如何，都统一存储到规范化的 key 下，
    但记录原始的第一个组分。
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
                # 保留原始顺序
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

        # 使用规范化的 key（按字母序排序），但记录原始顺序
        key_canonical = (target_subl, tuple(sorted(comps_original)))

        if key_canonical not in self._binary_L_param_cache:
            self._binary_L_param_cache[key_canonical] = {
                'all': [],
                'odd': [],
                'first_comp_in_db': comps_original[0]  # 记录数据库中的第一个组分
            }

        L_expr = param['parameter']
        order = param.get('parameter_order', 0)

        self._binary_L_param_cache[key_canonical]['all'].append((L_expr, order))

        if order % 2 != 0:  # 奇数阶
            self._binary_L_param_cache[key_canonical]['odd'].append(L_expr)


def _get_uem_d_term_v2(self, dbe, comp_k, comp_i, subl_index):
    """
    根据UEM公式计算 d_ki 项（简化版本）。

    关键：通过 first_comp_in_db 判断是否需要符号修正
    """

    if not hasattr(self, '_binary_L_param_cache'):
        self._build_binary_L_cache(dbe)

    # 使用规范化的 key 查找
    key_canonical = (subl_index, tuple(sorted((comp_k, comp_i))))

    cache_entry = self._binary_L_param_cache.get(key_canonical, None)
    if cache_entry is None:
        return S.Zero

    odd_L_terms = cache_entry.get('odd', [])
    if not odd_L_terms:
        return S.Zero

    # 判断是否需要符号修正
    first_comp_in_db = cache_entry['first_comp_in_db']

    # 如果数据库中第一个组分是 comp_k，则顺序一致，sign = 1
    # 如果数据库中第一个组分是 comp_i，则顺序相反，sign = -1
    if first_comp_in_db == comp_k:
        sign_correction = 1
    elif first_comp_in_db == comp_i:
        sign_correction = -1
    else:
        # 不应该发生
        raise ValueError(f"Unexpected first component in database: {first_comp_in_db}")

    # 应用符号修正
    d_ki = (S(2) / (v.R * v.T)) * sign_correction * Add(*odd_L_terms)
    return d_ki


# ============================================================================
# 总结
# ============================================================================
"""
推荐使用 v2 版本（简化版本），因为：

1. 逻辑更清晰：统一使用规范化的 key（sorted），但记录原始顺序信息
2. 代码更简洁：不需要维护双向引用
3. 易于理解：直接通过 first_comp_in_db 判断符号

修复步骤：
1. 将 _build_binary_L_cache 替换为 _build_binary_L_cache_v2
2. 将 _get_uem_d_term 替换为 _get_uem_d_term_v2
3. 在缓存中添加 'first_comp_in_db' 字段记录原始顺序
4. 在获取 d_ki 时根据 first_comp_in_db 判断符号

关键原理：
- L0(A,B) = L0(B,A)  （偶数阶对称）
- L1(A,B) = -L1(B,A) （奇数阶反对称）
- d_ki = (2/RT) * (L1 + L3 + L5 + ...) [以 (k,i) 顺序]
- 如果数据库是 (i,k) 顺序，则需要取反：d_ki = -(2/RT) * (L1_db + L3_db + ...)
"""
