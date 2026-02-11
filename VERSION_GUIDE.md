# skyborn-cdo 版本管理指南

## 版本号格式

```
CDO_MAJOR.CDO_MINOR.CDO_PATCH.WRAPPER_VERSION
例如: 2.5.3.1
```

- **前三位** (2.5.3): 跟随上游 CDO 版本
- **最后一位** (.1, .2, .3...): skyborn-cdo 包装器的版本号

## 何时更新版本号

### 更新最后一位（包装器版本）- 仅 skyborn-cdo 代码变更，CDO 版本不变

**场景：**
- ✅ 修复 Bug（如 Windows 进程挂起问题）
- ✅ 添加新功能（如通配符支持、help 系统）
- ✅ 性能优化
- ✅ 文档改进
- ✅ **添加新 Python 版本支持**（如 Python 3.14）
- ✅ 改进错误处理
- ✅ API 增强

**操作：**
```
2.5.3.1 → 2.5.3.2 → 2.5.3.3 ...
```

**需要修改的文件：**
1. `src/skyborn_cdo/__init__.py` 中的 `__version__`
2. `pyproject.toml` 中的 `version`

### 更新前三位（CDO 版本升级）

**场景：**
- 🔄 升级到 CDO 2.5.4
- 🔄 升级到 CDO 2.6.0
- 🔄 升级到 CDO 3.0.0

**操作：**
```
2.5.3.x → 2.5.4.1  (CDO 2.5.3 → 2.5.4, 重置包装版本为 .1)
2.5.4.x → 2.6.0.1  (CDO 2.5.4 → 2.6.0)
```

**需要修改的文件：**
1. `src/skyborn_cdo/__init__.py` 中的 `__version__` 和 `__cdo_version__`
2. `pyproject.toml` 中的 `version`
3. 重新编译 CDO 二进制文件

## 版本演进示例

### 场景 1: 首次发布
```
2.5.3.1  ← 首次发布到 PyPI（2026-02-12）
```

### 场景 2: 添加 Python 3.14 支持
```
2.5.3.1  ← 当前版本
↓
2.5.3.2  ← 添加 Python 3.14 wheel 构建（只更新包装器版本）
```

**修改：**
- `pyproject.toml`: 添加 `"Programming Language :: Python :: 3.14"` 分类器
- `.github/workflows/build_wheels.yml`: 添加 Python 3.14 构建
- 更新版本号到 `2.5.3.2`

### 场景 3: 修复 Bug
```
2.5.3.2  ← 当前版本
↓
2.5.3.3  ← 修复错误处理 Bug
```

### 场景 4: 添加新功能
```
2.5.3.3  ← 当前版本
↓
2.5.3.4  ← 添加进度回调功能
```

### 场景 5: 升级 CDO 版本
```
2.5.3.4  ← 基于 CDO 2.5.3 的最后版本
↓
2.5.4.1  ← 升级到 CDO 2.5.4（重置包装版本号）
```

**修改：**
- 重新编译 CDO 2.5.4 二进制
- 更新 `__version__ = "2.5.4.1"`
- 更新 `__cdo_version__ = "2.5.4"`
- 测试所有功能

## 发布清单

每次发布到 PyPI 前检查：

- [ ] 更新 `src/skyborn_cdo/__init__.py` 中的 `__version__`
- [ ] 更新 `pyproject.toml` 中的 `version`
- [ ] 更新 `README.md`（如果有功能变更）
- [ ] 运行完整测试套件
- [ ] 检查 CHANGELOG（如果维护的话）
- [ ] Git 提交并打标签：`git tag v2.5.3.1`
- [ ] 构建 wheels：`python -m build`
- [ ] 上传到 PyPI：`twine upload dist/*`

## 版本号比较（PEP 440）

Python 会按以下规则比较版本号：

```
2.5.3.1 < 2.5.3.2 < 2.5.3.10 < 2.5.4.1 < 2.6.0.1
```

✅ 这确保了包管理器能正确识别新版本。

## 常见问题

**Q: 为什么不从 0.1.0 开始？**  
A: 因为 skyborn-cdo 是 CDO 的包装器，版本号跟随 CDO 能让用户立即知道绑定的 CDO 版本。

**Q: 如果只更新文档需要发布新版本吗？**  
A: 看情况。如果是重大文档更新（如添加大量示例），可以发布补丁版本（如 2.5.3.2）。小改动可以在下次版本一起发布。

**Q: 如何处理 breaking changes？**  
A: 如果是 API 不兼容变更，考虑：
- 在 CHANGELOG 中明确标注
- 在 README 中添加迁移指南
- 如果变更很大，可以考虑主版本号跳跃（但这应该很少见，因为主要跟随 CDO）

**Q: 如果 CDO 有安全补丁怎么办？**  
A: 立即升级 CDO 并发布新版本（如 2.5.3.5 包含 CDO 2.5.3 的安全补丁）。

## 自动化建议

可以考虑使用工具自动同步版本号：

```bash
# 使用 bump2version 或 bumpver
pip install bumpver

# 配置 .bumpversion.cfg 或 pyproject.toml
# 然后可以：
bumpver update --patch  # 2.5.3.1 → 2.5.3.2
```

或者创建一个简单的 Python 脚本：

```python
# scripts/bump_version.py
import re

def bump_wrapper_version():
    # 读取当前版本
    with open("src/skyborn_cdo/__init__.py") as f:
        content = f.read()
    
    # 提取版本
    match = re.search(r'__version__ = "(\d+)\.(\d+)\.(\d+)\.(\d+)"', content)
    major, minor, patch, wrapper = match.groups()
    
    # 递增包装版本
    new_wrapper = int(wrapper) + 1
    new_version = f"{major}.{minor}.{patch}.{new_wrapper}"
    
    # 更新所有文件
    # ... 实现自动更新逻辑
    
    print(f"Version bumped: {match.group(1)} → {new_version}")
```

## 总结

**当前版本：** `2.5.3.1` （首次 PyPI 发布）

**未来更新：**
- 小更新/Bug 修复/Python 版本支持 → 增加最后一位：`.2`, `.3`, `.4` ...
- CDO 版本升级 → 更新前三位并重置最后一位：`2.5.4.1`, `2.6.0.1` ...

这个方案简单、清晰，符合 Python 包管理规范（PEP 440）。
