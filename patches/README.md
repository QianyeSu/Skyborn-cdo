# Windows Compatibility Patch System

## 📖 Overview

This project uses a **smart code pattern matching system** (rather than traditional line-number-based patch files) to adapt CDO source code for Windows compilation. This approach **significantly reduces maintenance costs** and is likely to work without modification even after CDO version updates.

---

## 🔧 Patch Script: `scripts/patch_cdo_windows.py`

### How it Works

- **Code pattern matching**: Uses regular expressions and string patterns to find code that needs modification
- **Automatic error handling**: Issues warnings if a pattern is not found instead of failing outright
- **Version resilience**: Works reliably even with minor CDO code structure changes

### Usage

```bash
# 1. Verify patches (don't modify files, check if all patch points exist)
python scripts/patch_cdo_windows.py verify

# 2. Apply patches (modify source code)
python scripts/patch_cdo_windows.py apply

# 3. Restore original code (undo all modifications)
python scripts/patch_cdo_windows.py restore

# Optional: specify CDO source directory
python scripts/patch_cdo_windows.py apply --cdo-src /path/to/cdo
```

---

## 📝 Patch Details

| File | Modification | Reason |
|------|---------|------|
| **src/cdo.cc** | Add `io.h` / `windows.h` headers | Windows requires `_isatty` / `_fileno` |
| | Modify `cdo_init_is_tty()` to Windows implementation | fstat/S_ISCHR unreliable in MinGW |
| | Add `fflush(stdout)` before `clear_processes()` | Prevent output loss on Windows CI exit |
| **src/cdo_getopt.cc** | Guard `sys/ioctl.h` | Not available on Windows |
| | Guard `unistd.h` | Replaced with Windows API |
| **Other 6 files** | Guard `unistd.h` consistently | Needs conditional compilation for MinGW |

---

## 🚀 Integration with Build Process

### CI Build (``.github/workflows/build.yml`)

```bash
# Automatically called in scripts/build_cdo_windows.sh
python scripts/patch_cdo_windows.py apply --cdo-src vendor/cdo
# ... configure && make ...
# vendor/ directory is automatically cleaned by Git after build, no manual restore needed
```

### Local Build (`scripts/build_local_windows.sh`)

```bash
# Apply patches → Compile → Auto-restore
python scripts/patch_cdo_windows.py apply
# ... compile ...
python scripts/patch_cdo_windows.py restore  # Keep Git repository clean
```

---

## ✨ Advantages Over Traditional Patch Files

| Traditional .patch File | Smart Python Script |
|----------------|----------------|
| ❌ Depends on fixed line numbers | ✅ Based on code pattern matching |
| ❌ Breaks after any CDO update | ✅ Still works with minor code changes |
| ❌ Cryptic error messages | ✅ Clear status for each patch point |
| ❌ Logic hidden in diff format | ✅ Code is clear and readable |
| ❌ Difficult to debug | ✅ Easy to modify and extend |

---

## 🔄 CDO Version Upgrade Workflow

### Typical scenario: CDO 2.5.4 → 2.5.5

1. **Update vendor/cdo source**
   ```bash
   cd vendor/cdo
   # ... synchronize upstream code or extract new version ...
   ```

2. **Verify if patches still work**
   ```bash
   python scripts/patch_cdo_windows.py verify
   ```

3. **Act based on verification results**

   - **Case A: All patch points found** → ✅ No modification needed
   - **Case B: Some patch points not found** → ⚠️ Adjust regex patterns
   - **Case C: CDO already fixed Windows compatibility** → 🎉 Remove corresponding patch logic

4. **If adjustment needed, edit `patch_cdo_windows.py`**
   - Open the script and find the `patches = [...]` section
   - Modify the corresponding regex patterns or replacement text
   - Re-verify

5. **Commit the update**
   ```bash
   git add scripts/patch_cdo_windows.py
   git commit -m "chore: update Windows patches for CDO 2.5.5"
   ```

---

## 🛠️ Adding New Patch Logic

If you discover a new Windows compatibility issue, add it to `patch_cdo_windows.py`:

```python
patches = [
    # Existing patches...
    
    ("src/new_file.cc", [
        ("Description of what this patch does",
         re.compile(r'regex pattern for original code'),
         r'replacement code'),
    ]),
]
```

---

## 📚 Reference Information

- **Patch script**: `scripts/patch_cdo_windows.py`
- **CI build script**: `scripts/build_cdo_windows.sh`
- **Local build script**: `scripts/build_local_windows.sh`
- **Legacy patch file (deprecated)**: `patches/windows-compat.patch`

---

**Last Updated**: 2026-02-13  
**Compatible with CDO**: 2.5.4
