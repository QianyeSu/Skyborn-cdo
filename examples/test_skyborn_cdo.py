"""
skyborn-cdo Local Test Script
=============================
Run in the skyborn_dev environment:
    conda activate skyborn_dev
    set PATH=C:\\msys64\\mingw64\\bin;%PATH%
    python examples/test_skyborn_cdo.py
"""

import os
import sys
import glob
import tempfile

# ─── Configuration ───────────────────────────────────────
# Ensure MinGW DLL is in PATH (required for local development)
MINGW_BIN = r"C:\msys64\mingw64\bin"
if os.path.isdir(MINGW_BIN) and MINGW_BIN not in os.environ.get("PATH", ""):
    os.environ["PATH"] = MINGW_BIN + os.pathsep + os.environ["PATH"]

GPM_DATA_DIR = r"I:\GPM_Month\Data"
OUTPUT_DIR = r"I:\GPM_Month"


def separator(title: str):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")


def test_import():
    """Test 1: Package Import"""
    separator("Test 1: Import skyborn_cdo")
    from skyborn_cdo import Cdo
    cdo = Cdo()
    print(f"  ✓ Import successful: {repr(cdo)}")
    return cdo


def test_version(cdo):
    """Test 2: CDO Version"""
    separator("Test 2: CDO Version Information")
    version = cdo.version()
    print(f"  Version output:\n{version}")
    assert "2.5.3" in version, f"Version mismatch: {version}"
    print(f"  ✓ CDO 2.5.3 confirmed")


def test_operators(cdo):
    """Test 3: Operator List"""
    separator("Test 3: Operator List")
    ops = cdo.operators()
    print(f"  Available operators: {len(ops)}")
    key_ops = ["sinfo", "mergetime", "remapbil", "sellonlatbox",
               "timavg", "yearmonmean", "fldmean", "griddes"]
    for op in key_ops:
        status = "✓" if op in ops else "✗"
        print(f"  {status} {op}")
    print(f"  ✓ Operator list loaded successfully")


def test_sinfon(cdo):
    """Test 4: sinfon - Read NetCDF File Information"""
    separator("Test 4: sinfon — Read File Information")
    infiles = sorted(glob.glob(os.path.join(GPM_DATA_DIR, "GPM_Precip_*.nc")))
    if not infiles:
        print(f"  ⚠ No GPM files found: {GPM_DATA_DIR}")
        return None

    print(f"  Found {len(infiles)} GPM files")
    test_file = infiles[0]
    print(f"  Test file: {os.path.basename(test_file)}")

    info = cdo.sinfon(input=test_file)
    if isinstance(info, list):
        info = "\n".join(info)
    print(f"  sinfon output:\n{info[:500]}")
    print(f"  ✓ sinfon successful")
    return infiles


def test_griddes(cdo, infiles):
    """Test 5: griddes - Get Grid Description"""
    separator("Test 5: griddes — Grid Description")
    grid = cdo.griddes(input=infiles[0])
    if isinstance(grid, list):
        grid = "\n".join(grid)
    print(f"  Grid information:\n{grid[:400]}")
    print(f"  ✓ griddes successful")


def test_remapbil(cdo, infiles):
    """Test 6: remapbil - Remap to 1° Grid"""
    separator("Test 6: remapbil — Remap to 1° Grid")
    outfile = os.path.join(OUTPUT_DIR, "test_remap_1deg.nc")
    if os.path.exists(outfile):
        os.remove(outfile)

    cdo.remapbil("r360x180", input=infiles[0], output=outfile)
    size_mb = os.path.getsize(outfile) / (1024 * 1024)
    print(f"  Input: {os.path.basename(infiles[0])}")
    print(f"  Output: {outfile}")
    print(f"  Size: {size_mb:.2f} MB")

    # Verify output grid
    grid = cdo.griddes(input=outfile)
    if isinstance(grid, list):
        grid = "\n".join(grid)
    assert "xsize     = 360" in grid, "Remapped grid X is incorrect"
    assert "ysize     = 180" in grid, "Remapped grid Y is incorrect"
    print(f"  ✓ Remapping to r360x180 successful (360x180)")

    os.remove(outfile)
    return True


def test_mergetime(cdo, infiles):
    """Test 7: mergetime - Merge Multiple Time Steps"""
    separator("Test 7: mergetime — Merge 3 Months")
    files_3months = infiles[:3]
    outfile = os.path.join(OUTPUT_DIR, "test_merge_3months.nc")
    if os.path.exists(outfile):
        os.remove(outfile)

    print(f"  Input files:")
    for f in files_3months:
        print(f"    {os.path.basename(f)}")

    cdo.mergetime(input=files_3months, output=outfile)
    size_mb = os.path.getsize(outfile) / (1024 * 1024)
    print(f"  Output: {outfile}")
    print(f"  Size: {size_mb:.2f} MB")
    print(f"  ✓ mergetime successful")

    os.remove(outfile)
    return True


def test_sellonlatbox(cdo, infiles):
    """Test 8: sellonlatbox - Crop China Region"""
    separator("Test 8: sellonlatbox — Crop China Region")
    outfile = os.path.join(OUTPUT_DIR, "test_china_region.nc")
    if os.path.exists(outfile):
        os.remove(outfile)

    # China's approximate boundaries: 73-135E, 18-54N
    cdo.sellonlatbox("73,135,18,54", input=infiles[0], output=outfile)
    size_mb = os.path.getsize(outfile) / (1024 * 1024)
    print(f"  Region: 73-135°E, 18-54°N (China)")
    print(f"  Output: {outfile}")
    print(f"  Size: {size_mb:.2f} MB")
    print(f"  ✓ sellonlatbox crop successful")

    os.remove(outfile)
    return True


def test_pipeline(cdo, infiles):
    """Test 9: Pipeline Combination - mergetime + remapbil"""
    separator("Test 9: Pipeline — mergetime + remapbil")
    files_3months = infiles[:3]
    outfile = os.path.join(OUTPUT_DIR, "test_pipeline_merge_remap.nc")
    if os.path.exists(outfile):
        os.remove(outfile)

    # CDO supports -op1 -op2 input pipeline chaining
    cdo.remapbil(
        "r360x180",
        input="-mergetime " + " ".join(files_3months),
        output=outfile
    )
    size_mb = os.path.getsize(outfile) / (1024 * 1024)
    print(f"  Pipeline: remapbil(r360x180, mergetime(3 files))")
    print(f"  Output: {outfile}")
    print(f"  Size: {size_mb:.2f} MB")
    print(f"  ✓ Pipeline operation successful")

    os.remove(outfile)
    return True


def test_tempfile_output(cdo, infiles):
    """Test 10: Using Temporary File as Output"""
    separator("Test 10: Temporary File Output")
    with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as tmp:
        tmpfile = tmp.name
    try:
        cdo.timmean(input=infiles[0], output=tmpfile)
        size_kb = os.path.getsize(tmpfile) / 1024
        print(f"  timmean output: {tmpfile}")
        print(f"  Temporary file size: {size_kb:.1f} KB")
        print(f"  ✓ Temporary file output successful")
    finally:
        if os.path.exists(tmpfile):
            os.remove(tmpfile)


def main():
    print("=" * 60)
    print("  skyborn-cdo Comprehensive Test Suite")
    print(f"  Python {sys.version}")
    print("=" * 60)

    passed = 0
    failed = 0
    total = 0

    try:
        # Basic tests (no data files required)
        cdo = test_import()
        total += 1
        passed += 1

        test_version(cdo)
        total += 1
        passed += 1

        test_operators(cdo)
        total += 1
        passed += 1

        # Data tests (require GPM files)
        infiles = test_sinfon(cdo)
        total += 1
        passed += 1

        if infiles:
            test_griddes(cdo, infiles)
            total += 1
            passed += 1

            test_remapbil(cdo, infiles)
            total += 1
            passed += 1

            test_mergetime(cdo, infiles)
            total += 1
            passed += 1

            test_sellonlatbox(cdo, infiles)
            total += 1
            passed += 1

            test_pipeline(cdo, infiles)
            total += 1
            passed += 1

            test_tempfile_output(cdo, infiles)
            total += 1
            passed += 1

    except Exception as e:
        failed += 1
        total += 1
        print(f"\n  ✗ Test failed: {e}")
        import traceback
        traceback.print_exc()

    separator("Test Results")
    print(f"  Passed: {passed}/{total}")
    if failed:
        print(f"  Failed: {failed}/{total}")
    else:
        print(f"  All tests passed! ✓")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
