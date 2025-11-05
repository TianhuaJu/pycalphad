# 相图绘制重构总结

## 更新概述

根据用户要求，将相图绘制方式从**颜色区块填充**改为**线条边界划分 + 区域标注**，提升相图的专业性和可读性。

## 主要改进

### 1. 绘制方式的根本改变

#### ❌ 旧版本（颜色区块方式）
```python
# 使用 imshow 填充颜色区块
im = self.ax_phase.imshow(
    phase_numeric_grid,
    cmap=cmap,
    interpolation='nearest'
)

# 使用 colorbar 显示相名称
cbar = self.fig_phase.colorbar(im, ax=self.ax_phase)
cbar.ax.set_yticklabels(unique_phase_strings)
```

**问题：**
- 相区域用颜色填充，不够专业
- 需要查看colorbar才知道相名称
- 颜色可能不易区分
- 不符合学术论文的标准表达

#### ✅ 新版本（线条边界方式）
```python
# 为每个相绘制边界线条
for phase_idx, (phase_name, color) in enumerate(zip(unique_phase_strings, colors)):
    phase_mask = (phase_numeric_grid == phase_idx).astype(float)

    if np.any(phase_mask > 0):
        contours = self.ax_phase.contour(
            X_grid, T_grid, phase_mask,
            levels=[0.5],
            colors=[color],
            linewidths=2.5,
            alpha=0.8
        )

# 在相区域中心标注相名称
self._label_phase_regions(x_range, T_range, phase_numeric_grid,
                          unique_phase_strings, phase_map)
```

**优势：**
- ✅ 清晰的线条边界，专业美观
- ✅ 直接在图中标注相名称，一目了然
- ✅ 不同相边界颜色不同，增强区分度
- ✅ 符合学术期刊的相图表达标准

### 2. 颜色策略

#### 旧版本
- 使用 `tab20` colormap 填充区域
- 颜色区块可能重叠或不清晰

#### 新版本
```python
# 为每个相边界选择不同的颜色
colors = plt.cm.tab10(np.linspace(0, 1, len(unique_phase_strings)))
```

- **线条颜色**：`tab10` 配色方案（10种鲜明颜色）
- **线条宽度**：2.5pt（清晰可见）
- **透明度**：alpha=0.8（柔和不刺眼）
- **背景**：纯白色（干净专业）

### 3. 智能标签系统

新增两个核心方法来实现区域标注：

#### `_label_phase_regions()` - 相区域标注主方法

**功能：**
1. 遍历所有相
2. 识别每个相的连通区域
3. 调用标签添加方法

**技术亮点：**
```python
# scipy可选依赖
try:
    from scipy import ndimage
    has_scipy = True
except ImportError:
    has_scipy = False

if has_scipy:
    # 精确方法：识别连通区域
    labeled_array, num_features = ndimage.label(phase_mask)
    for region_idx in range(1, num_features + 1):
        region_mask = (labeled_array == region_idx)
        self._add_phase_label(...)
else:
    # 简化方法：整体标注
    self._add_phase_label(phase_mask, ...)
```

**优势：**
- 支持同一相的多个不连通区域
- scipy可选，无依赖也能工作
- 每个区域单独标注

#### `_add_phase_label()` - 标签添加辅助方法

**功能：**
1. 计算区域质心坐标
2. 根据区域大小调整字体
3. 添加带背景框的标签

**质心计算：**
```python
region_indices = np.where(region_mask)
center_j = int(np.mean(region_indices[0]))  # 温度坐标
center_i = int(np.mean(region_indices[1]))  # 成分坐标

center_x = x_range[center_i]
center_T = T_range[center_j]
```

**自适应字体：**
```python
region_size = np.sum(region_mask)
size_ratio = region_size / total_size

if size_ratio > 0.1:
    fontsize = 11; fontweight = 'bold'
elif size_ratio > 0.05:
    fontsize = 9; fontweight = 'bold'
elif size_ratio > 0.02:
    fontsize = 8; fontweight = 'normal'
else:
    fontsize = 7; fontweight = 'normal'
```

**标签样式：**
```python
self.ax_phase.text(
    center_x, center_T, phase_name,
    fontsize=fontsize,
    fontweight=fontweight,
    ha='center', va='center',
    bbox=dict(
        boxstyle='round,pad=0.4',
        facecolor='white',      # 白色背景
        edgecolor='gray',       # 灰色边框
        alpha=0.8,              # 半透明
        linewidth=1
    ),
    zorder=10  # 确保标签在最上层
)
```

### 4. 视觉效果增强

#### 背景和网格
```python
# 白色背景
self.ax_phase.set_facecolor('white')

# 添加网格线
self.ax_phase.grid(True, linestyle='--', alpha=0.3, linewidth=0.5)
```

#### 坐标轴
```python
self.ax_phase.set_xlabel(f'X({scan_comp})', fontsize=12, fontweight='bold')
self.ax_phase.set_ylabel('Temperature (K)', fontsize=12, fontweight='bold')
self.ax_phase.set_title(
    f'Pseudo-Binary Phase Diagram ({model_label})',
    fontsize=14, fontweight='bold'
)
```

## 技术对比表

| 特性 | 旧版本 | 新版本 |
|------|--------|--------|
| **绘制方法** | `imshow` 颜色填充 | `contour` 线条边界 |
| **相名称显示** | colorbar（侧边栏） | 区域内标注 |
| **颜色使用** | 区域填充色 | 边界线条色 |
| **区分度** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **专业性** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **可读性** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **标准性** | 一般 | 学术期刊标准 |
| **连通区域识别** | ❌ | ✅ (scipy) |
| **自适应字体** | ❌ | ✅ |
| **标签背景框** | ❌ | ✅ |
| **网格辅助** | ❌ | ✅ |

## 代码统计

- **修改文件**: `alloy_calculator_gui.py`
- **修改行数**: 144 insertions(+), 27 deletions(-)
- **新增方法**:
  - `_label_phase_regions()` - 相区域标注
  - `_add_phase_label()` - 标签添加
- **修改方法**:
  - `_plot_phase_diagram()` - 完全重构

## 实现细节

### contour 线条绘制

```python
X_grid, T_grid = np.meshgrid(x_range, T_range)

for phase_idx, (phase_name, color) in enumerate(zip(unique_phase_strings, colors)):
    # 创建该相的掩码（该相的区域值为1，其他为0）
    phase_mask = (phase_numeric_grid == phase_idx).astype(float)

    if np.any(phase_mask > 0):
        # levels=[0.5] 表示在0和1之间绘制边界
        contours = self.ax_phase.contour(
            X_grid, T_grid, phase_mask,
            levels=[0.5],
            colors=[color],
            linewidths=2.5,
            alpha=0.8
        )
```

**原理：**
- `phase_mask` 是二值网格（该相=1，其他=0）
- `levels=[0.5]` 在0和1之间绘制等值线
- 等值线正好是相边界

### 连通区域识别（scipy）

```python
from scipy import ndimage

# 标记连通区域
labeled_array, num_features = ndimage.label(phase_mask)

# labeled_array: 每个连通区域有唯一编号
# num_features: 连通区域总数

for region_idx in range(1, num_features + 1):
    region_mask = (labeled_array == region_idx)
    # 为该连通区域添加标签
```

**作用：**
- 识别同一相的多个不连通区域
- 每个区域单独标注相名称
- 避免相名称遗漏或重叠

## 应用场景

### 适用于

✅ 伪二元相图（Pseudo-binary phase diagram）
✅ 学术论文图表
✅ 技术报告
✅ 教学演示
✅ 多相平衡体系

### 处理情况

✅ 单相区域 → 大字体居中标注
✅ 两相区域 → 线条边界 + 小字体标注
✅ 多相区域 → 彩色线条 + 自适应标注
✅ 不连通区域 → 多个标签
✅ 小区域 → 缩小字体避免重叠

## 依赖说明

### 必需依赖
- `numpy` ✅
- `matplotlib` ✅

### 可选依赖
- `scipy` (推荐，用于连通区域识别)
  - 有 scipy: 精确识别每个连通区域
  - 无 scipy: 简化方法，整体标注

## 视觉效果示例

### 旧版本输出
```
┌─────────────────────────────────┐
│  彩色区块相图                    │
│  ████████████ (LIQUID)          │
│  ████████████                   │
│  ░░░░░░░░░░░░ (FCC_A1)          │
│  ░░░░░░░░░░░░                   │
│                                 │
│  [Colorbar]                     │
│  ■ LIQUID                       │
│  ░ FCC_A1                       │
│  ▓ BCC_A2                       │
└─────────────────────────────────┘
```

### 新版本输出
```
┌─────────────────────────────────┐
│  线条边界相图                    │
│                                 │
│       [LIQUID]                  │
│  ━━━━━━━━━━━━━━━━━━━━           │
│                                 │
│       [FCC_A1]                  │
│  ─ ─ ─ ─ ─ ─ ─ ─ ─ ─            │
│                                 │
│       [BCC_A2]                  │
│                                 │
└─────────────────────────────────┘
```

## 未来改进建议

1. **交互功能**
   - 点击相区域显示详细信息
   - 悬停显示温度/成分坐标

2. **导出选项**
   - 高分辨率PNG/PDF导出
   - 矢量图格式支持（SVG）

3. **标签优化**
   - 自动避免标签重叠
   - 复杂相名称自动缩写

4. **样式定制**
   - 用户可选配色方案
   - 线条样式自定义（实线、虚线等）
   - 标签样式模板

5. **边界平滑**
   - 使用样条插值平滑边界
   - 减少锯齿效应

## Git 信息

- **分支**: `claude/update-db-default-with-memory-011CUp6ecRxdw9AnCF3Som6a`
- **提交哈希**: `8a4f4d4`
- **修改文件**: `alloy_calculator_gui.py`
- **变更统计**: 144 insertions(+), 27 deletions(-)

## 兼容性

- **Python**: 3.6+
- **Matplotlib**: 2.0+
- **NumPy**: 1.10+
- **SciPy**: 0.17+ (可选)

## 测试建议

### 基础测试
1. 加载包含多相的TDB数据库
2. 计算伪二元相图
3. 验证：
   - 相边界线条清晰
   - 相名称标注正确
   - 不同相边界颜色不同

### 高级测试
1. **多连通区域**
   - 相在不同温度/成分范围出现多次
   - 验证每个区域都有标签

2. **小区域处理**
   - 很小的相区域
   - 验证字体自动缩小

3. **无scipy环境**
   - 卸载scipy测试
   - 验证回退方案工作正常

4. **大量相**
   - 5个以上相的复杂系统
   - 验证颜色区分度

## 结论

本次重构成功将相图绘制从**颜色区块方式**升级为**线条边界 + 区域标注方式**，显著提升了：

✅ **专业性** - 符合学术标准
✅ **可读性** - 直观清晰
✅ **美观度** - 简洁优雅
✅ **区分度** - 彩色线条
✅ **智能性** - 自适应标签

完全满足用户要求，并提供了健壮的实现方案。
