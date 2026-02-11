# skyborn-cdo Version Management Guide

## Version Number Format

```
CDO_MAJOR.CDO_MINOR.CDO_PATCH.WRAPPER_VERSION
Example: 2.5.3.1
```

- **First three digits** (2.5.3): Follow upstream CDO version
- **Last digit** (.1, .2, .3...): skyborn-cdo wrapper version number

## When to Update Version Number

### Update the last digit (wrapper version) - Only skyborn-cdo code changes, CDO version unchanged

**Scenarios:**
- ✅ Bug fixes (e.g., Windows process hanging issue)
- ✅ New features (e.g., wildcard support, help system)
- ✅ Performance optimizations
- ✅ Documentation improvements
- ✅ **Add new Python version support** (e.g., Python 3.14)
- ✅ Error handling improvements
- ✅ API enhancements

**Operation:**
```
2.5.3.1 → 2.5.3.2 → 2.5.3.3 ...
```

**Files to modify:**
1. `__version__` in `src/skyborn_cdo/__init__.py`
2. `version` in `pyproject.toml`

### Update first three digits (CDO version upgrade)

**Scenarios:**
- 🔄 Upgrade to CDO 2.5.4
- 🔄 Upgrade to CDO 2.6.0
- 🔄 Upgrade to CDO 3.0.0

**Operation:**
```
2.5.3.x → 2.5.4.1  (CDO 2.5.3 → 2.5.4, reset wrapper version to .1)
2.5.4.x → 2.6.0.1  (CDO 2.5.4 → 2.6.0)
```

**Files to modify:**
1. `__version__` and `__cdo_version__` in `src/skyborn_cdo/__init__.py`
2. `version` in `pyproject.toml`
3. Recompile CDO binary files

## Version Evolution Examples

### Scenario 1: First Release
```
2.5.3.1  ← First PyPI release (2026-02-12)
```

### Scenario 2: Add Python 3.14 Support
```
2.5.3.1  ← Current version
↓
2.5.3.2  ← Add Python 3.14 wheel builds (only update wrapper version)
```

**Modifications:**
- `pyproject.toml`: Add `"Programming Language :: Python :: 3.14"` classifier
- `.github/workflows/build_wheels.yml`: Add Python 3.14 builds
- Update version to `2.5.3.2`

### Scenario 3: Bug Fix
```
2.5.3.2  ← Current version
↓
2.5.3.3  ← Fix error handling bug
```

### Scenario 4: Add New Feature
```
2.5.3.3  ← Current version
↓
2.5.3.4  ← Add progress callback feature
```

### Scenario 5: Upgrade CDO Version
```
2.5.3.4  ← Last version based on CDO 2.5.3
↓
2.5.4.1  ← Upgrade to CDO 2.5.4 (reset wrapper version)
```

**Modifications:**
- Recompile CDO 2.5.4 binaries
- Update `__version__ = "2.5.4.1"`
- Update `__cdo_version__ = "2.5.4"`
- Test all functionality

## Release Checklist

Before each PyPI release, verify:

- [ ] Update `__version__` in `src/skyborn_cdo/__init__.py`
- [ ] Update `version` in `pyproject.toml`
- [ ] Update `README.md` (if feature changes)
- [ ] Run complete test suite
- [ ] Check CHANGELOG (if maintaining one)
- [ ] Git commit and tag: `git tag v2.5.3.1`
- [ ] Build wheels: `python -m build`
- [ ] Upload to PyPI: `twine upload dist/*`

## Version Number Comparison (PEP 440)

Python compares version numbers using these rules:

```
2.5.3.1 < 2.5.3.2 < 2.5.3.10 < 2.5.4.1 < 2.6.0.1
```

✅ This ensures package managers correctly identify newer versions.

## FAQ

**Q: Why not start from 0.1.0?**  
A: Because skyborn-cdo is a CDO wrapper, version numbers that follow CDO let users immediately know which CDO version is bundled.

**Q: Do we need to publish a new version for documentation-only updates?**  
A: It depends. For major documentation updates (e.g., adding many examples), consider releasing a patch version (e.g., 2.5.3.2). Minor changes can be batched in the next release.

**Q: How to handle breaking changes?**  
A: If API-incompatible changes occur, consider:
- Clearly mark in CHANGELOG
- Add migration guide in README
- For major changes, consider a major version bump (though this should be rare since we follow CDO versioning)

**Q: What if CDO issues a security patch?**  
A: Immediately upgrade CDO and release a new version (e.g., 2.5.3.5 containing CDO 2.5.3 security patches).

## Automation Suggestions

Consider using tools to automatically sync version numbers:

```bash
# Using bump2version or bumpver
pip install bumpver

# Configure .bumpversion.cfg or pyproject.toml
# Then:
bumpver update --patch  # 2.5.3.1 → 2.5.3.2
```

Or create a simple Python script:

```python
# scripts/bump_version.py
import re

def bump_wrapper_version():
    # Read current version
    with open("src/skyborn_cdo/__init__.py") as f:
        content = f.read()
    
    # Extract version
    match = re.search(r'__version__ = "(\d+)\.(\d+)\.(\d+)\.(\d+)"', content)
    major, minor, patch, wrapper = match.groups()
    
    # Increment wrapper version
    new_wrapper = int(wrapper) + 1
    new_version = f"{major}.{minor}.{patch}.{new_wrapper}"
    
    # Update all files
    # ... implement automatic update logic
    
    print(f"Version bumped: {match.group(1)} → {new_version}")
```

## Summary

**Current version:** `2.5.3.1` (First PyPI release)

**Future updates:**
- Minor fixes/Bug fixes/Python version support → Increment last digit: `.2`, `.3`, `.4` ...
- CDO version upgrade → Update first three digits and reset last digit: `2.5.4.1`, `2.6.0.1` ...

This approach is simple, clear, and complies with Python packaging standards (PEP 440).
