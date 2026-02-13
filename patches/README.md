# Windows 兼容性补丁系统

## 📖 概述

本项目使用**智能代码匹配**系统（而非传统的行号式 patch 文件）对 CDO 源码进行 Windows 编译适配。这种方式**显著降低维护成本**，CDO 版本更新后大概率无需修改补丁逻辑。

---

## 🔧 补丁脚本：`scripts/patch_cdo_windows.py`

### 工作原理

- **基于代码特征匹配**：使用正则表达式和字符串模式查找需要修改的代码
- **自动容错**：如果某个模式未找到会给出警告，而非直接失败
- **版本适应性强**：CDO 代码结构小幅变化时仍能正常工作

### 使用方法

```bash
# 1. 验证补丁（不修改文件，检查所有补丁点是否存在）
python scripts/patch_cdo_windows.py verify

# 2. 应用补丁（修改源码）
python scripts/patch_cdo_windows.py apply

# 3. 恢复原始代码（撤销所有修改）
python scripts/patch_cdo_windows.py restore

# 可选：指定 CDO 源码目录
python scripts/patch_cdo_windows.py apply --cdo-src /path/to/cdo
```

---

## 📝 补丁内容说明

| 文件 | 修改内容 | 原因 |
|------|---------|------|
| **src/cdo.cc** | 添加 `io.h` / `windows.h` 头文件 | Windows 需要 `_isatty` / `_fileno` |
| | 修改 `cdo_init_is_tty()` 为 Windows 实现 | fstat/S_ISCHR 在 MinGW 中不稳定 |
| | 在 `clear_processes()` 前添加 `fflush(stdout)` | 避免 Windows CI 退出时输出丢失 |
| **src/cdo_getopt.cc** | 屏蔽 `sys/ioctl.h` | Windows 不支持 |
| | 屏蔽 `unistd.h` | 使用 Windows API 替代 |
| **其他 6 个文件** | 统一屏蔽 `unistd.h` | MinGW 需条件编译 |

---

## 🚀 构建流程中的集成

### CI 构建 (``.github/workflows/build.yml`)

```bash
# scripts/build_cdo_windows.sh 中自动调用
python scripts/patch_cdo_windows.py apply --cdo-src vendor/cdo
# ... configure && make ...
# 构建完成后 vendor/ 目录自动被 Git 清除，无需手动恢复
```

### 本地构建 (`scripts/build_local_windows.sh`)

```bash
# 应用补丁 → 编译 → 自动恢复
python scripts/patch_cdo_windows.py apply
# ... compile ...
python scripts/patch_cdo_windows.py restore  # 保持 Git 仓库干净
```

---

## ✨ 相比传统 patch 文件的优势

| 传统 .patch 文件 | 智能 Python 脚本 |
|----------------|----------------|
| ❌ 依赖固定行号 | ✅ 基于代码特征匹配 |
| ❌ CDO 更新后必定失败 | ✅ 小幅变化仍可工作 |
| ❌ 错误信息模糊 | ✅ 明确显示每个修改点状态 |
| ❌ 修改逻辑隐藏在 diff 格式中 | ✅ 代码清晰可读 |
| ❌ 难以调试 | ✅ 易于修改和扩展 |

---

## 🔄 CDO 版本升级时的处理流程

### 典型场景：CDO 2.5.4 → 2.5.5

1. **更新 vendor/cdo 源码**
   ```bash
   cd vendor/cdo
   # ... 同步上游代码或解压新版本 ...
   ```

2. **验证补丁是否仍然有效**
   ```bash
   python scripts/patch_cdo_windows.py verify
   ```

3. **根据验证结果采取行动**

   - **情况 A：所有补丁点均找到** → ✅ 无需修改
   - **情况 B：部分补丁点未找到** → ⚠️ 需调整正则表达式
   - **情况 C：CDO 已自行修复 Windows 兼容性** → 🎉 删除对应补丁逻辑

4. **如需调整，编辑 `patch_cdo_windows.py`**
   - 打开脚本，找到 `patches = [...]` 部分
   - 修改对应的正则表达式或替换文本
   - 重新验证

5. **提交更新**
   ```bash
   git add scripts/patch_cdo_windows.py
   git commit -m "chore: update Windows patches for CDO 2.5.5"
   ```

---

## 🛠️ 添加新的补丁逻辑

如果发现新的 Windows 兼容性问题，在 `patch_cdo_windows.py` 中添加：

```python
patches = [
    # 现有补丁...
    
    ("src/new_file.cc", [
        ("描述这个补丁的作用",
         re.compile(r'原始代码的正则表达式'),
         r'替换后的代码'),
    ]),
]
```

---

## 📚 参考信息

- **补丁脚本**：`scripts/patch_cdo_windows.py`
- **CI 构建脚本**：`scripts/build_cdo_windows.sh`
- **本地构建脚本**：`scripts/build_local_windows.sh`
- **历史 patch 文件（已弃用）**：`patches/windows-compat.patch`

---

**最后更新**：2026-02-13  
**适用 CDO 版本**：2.5.4
