# 数据库加载器更新总结

## 更新概述

基于 `claude/add-uem-integrated-model-011CUQHL2w2AVzHnBV4zr7dm` 分支，对 `alloy_calculator_gui.py` 进行了改进，实现了以下功能：

1. **默认目录改为工作目录** - 首次使用时默认打开当前工作目录
2. **目录记忆功能** - 自动记住上次打开的目录，下次启动时直接使用

## 主要修改

### 1. 添加必要的导入

```python
import json
from pathlib import Path
```

### 2. 初始化配置文件支持

在 `__init__` 方法中添加：
```python
# 配置文件路径 - 用于记住上次打开的目录
self.config_file = Path.home() / '.pycalphad_gui_config.json'
self.last_directory = self._load_last_directory()
```

配置文件位置：`~/.pycalphad_gui_config.json`

### 3. 实现目录加载功能

新增 `_load_last_directory()` 方法：
- **首次使用**：返回当前工作目录 `os.getcwd()`
- **后续使用**：从配置文件读取上次使用的目录
- **验证机制**：如果保存的目录不再存在，自动回退到工作目录
- **错误处理**：JSON解析失败时回退到工作目录

```python
def _load_last_directory(self):
    """加载上次打开的目录，默认为当前工作目录"""
    if self.config_file.exists():
        try:
            with open(self.config_file, 'r', encoding='utf-8') as f:
                config = json.load(f)
                last_dir = config.get('last_directory', os.getcwd())
                # 验证目录是否仍然存在
                if os.path.isdir(last_dir):
                    return last_dir
        except (json.JSONDecodeError, IOError) as e:
            print(f"加载配置文件失败: {e}")

    # 默认返回当前工作目录
    return os.getcwd()
```

### 4. 实现目录保存功能

新增 `_save_last_directory()` 方法：
- 将目录路径保存到JSON配置文件
- 支持中文路径（`ensure_ascii=False`）
- 错误处理：保存失败时打印警告

```python
def _save_last_directory(self, directory):
    """保存当前目录到配置文件"""
    try:
        config = {'last_directory': directory}
        with open(self.config_file, 'w', encoding='utf-8') as f:
            json.dump(config, f, indent=2, ensure_ascii=False)
    except IOError as e:
        print(f"保存配置文件失败: {e}")
```

### 5. 优化 `load_database()` 方法

**移除了旧的逻辑：**
```python
# 旧代码（已删除）
examples_dir = os.path.join(os.getcwd(), 'examples')

if os.path.isdir(examples_dir):
    default_dir = examples_dir
else:
    default_dir = os.getcwd()
```

**新的实现：**
```python
# 使用记忆的目录（首次为工作目录，后续为上次打开的目录）
default_dir = self.last_directory

file_paths = filedialog.askopenfilename(
    title="选择TDB数据库文件 (可多选, e.g., 'alcrni.tdb' + 'pure.tdb')",
    initialdir=default_dir,  # 使用记忆的目录
    filetypes=[("TDB文件", "*.tdb"), ("所有文件", "*.*")],
    multiple=True
)

if not file_paths:
    return

# 保存新的目录到配置文件（使用第一个文件的目录）
if file_paths:
    new_directory = os.path.dirname(file_paths[0])
    self.last_directory = new_directory
    self._save_last_directory(new_directory)
```

## 工作流程

```
1. 启动GUI
   ↓
2. 读取配置文件 ~/.pycalphad_gui_config.json
   ↓
3. 如果配置存在且目录有效 → 使用保存的目录
   否则 → 使用当前工作目录 os.getcwd()
   ↓
4. 用户点击"加载TDB数据库"
   ↓
5. 打开文件对话框（initialdir = 记忆的目录）
   ↓
6. 用户选择文件
   ↓
7. 提取文件所在目录
   ↓
8. 更新 self.last_directory
   ↓
9. 保存到配置文件
   ↓
10. 下次启动时使用该目录
```

## 配置文件格式

```json
{
  "last_directory": "/path/to/last/opened/directory"
}
```

示例：
```json
{
  "last_directory": "/home/user/pycalphad/examples"
}
```

## 功能对比

| 特性 | 旧版本 | 新版本 |
|------|--------|--------|
| 默认目录 | examples 子目录优先 | 工作目录或上次目录 |
| 目录记忆 | ❌ 无 | ✅ 持久化保存 |
| 配置文件 | ❌ 无 | ✅ ~/.pycalphad_gui_config.json |
| 目录验证 | ❌ 无 | ✅ 自动验证并回退 |
| 中文路径 | ❌ 可能有问题 | ✅ 完全支持 |
| 用户友好性 | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |

## 优势

1. **更符合用户习惯**
   - 首次使用时从工作目录开始，更直观
   - 记住常用目录，减少导航时间

2. **灵活性**
   - 不再强制使用 examples 目录
   - 支持多个项目，每次自动回到上次位置

3. **持久化**
   - 配置跨会话保存
   - 重启应用后保持用户偏好

4. **健壮性**
   - 目录验证机制
   - 完善的错误处理
   - 配置文件损坏时自动回退

5. **国际化支持**
   - 支持中文路径
   - UTF-8 编码

## Git 信息

- **基础分支**: `claude/add-uem-integrated-model-011CUQHL2w2AVzHnBV4zr7dm`
- **新分支**: `claude/update-db-default-with-memory-011CUp6ecRxdw9AnCF3Som6a`
- **提交哈希**: `0a1e96b`
- **修改文件**: `alloy_calculator_gui.py`
- **修改统计**: 1 file changed, 50 insertions(+), 13 deletions(-)

## 创建 Pull Request

```
https://github.com/TianhuaJu/pycalphad/pull/new/claude/update-db-default-with-memory-011CUp6ecRxdw9AnCF3Som6a
```

## 测试建议

1. **首次启动测试**
   - 删除配置文件（如果存在）
   - 启动GUI
   - 验证文件对话框默认目录是工作目录

2. **记忆功能测试**
   - 从某个目录加载数据库文件
   - 关闭并重启GUI
   - 再次点击加载，验证对话框打开在上次的目录

3. **目录验证测试**
   - 手动编辑配置文件，设置一个不存在的目录
   - 启动GUI
   - 验证回退到工作目录

4. **中文路径测试**
   - 在包含中文的路径下加载文件
   - 验证保存和读取都正常工作

## 兼容性说明

- **Python 版本**: 兼容 Python 3.6+
- **依赖**: 仅使用标准库（json, pathlib）
- **向后兼容**: 完全兼容现有代码，无破坏性更改
- **配置迁移**: 旧版本用户首次启动时会自动创建配置文件

## 后续改进建议

1. 添加"最近打开"列表，快速访问多个常用目录
2. 添加收藏夹功能，保存多个项目目录
3. 在GUI上显示当前默认目录
4. 添加"重置为工作目录"按钮
5. 支持更多配置选项（如窗口大小、位置等）
