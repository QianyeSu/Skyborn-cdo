#!/usr/bin/env python
"""
CI wheel smoke test — verifies CDO binary + ecCodes GRIB functionality.

This script is executed by cibuildwheel after each wheel is installed.
It validates:
  1. CDO binary loads and reports its version
  2. CDO can create data (NetCDF output, no ecCodes needed)
  3. CDO can convert NetCDF -> GRIB (ecCodes WRITE / encode)
  4. CDO can read GRIB back to NetCDF (ecCodes READ / decode)

If step 3/4 fails, ecCodes definitions are missing or broken.
"""

import os
import sys
import tempfile


def main():
    from skyborn_cdo import Cdo, get_cdo_path

    cdo = Cdo()
    tmpdir = tempfile.mkdtemp()
    passed = 0
    failed = 0

    # ------------------------------------------------------------------
    # Diagnostics: show package layout for debugging
    # ------------------------------------------------------------------
    import skyborn_cdo as _pkg

    pkg_dir = os.path.dirname(_pkg.__file__)
    print(f"Package dir: {pkg_dir}")
    try:
        cdo_bin = get_cdo_path()
        print(f"CDO binary:  {cdo_bin}")
    except FileNotFoundError as e:
        print(f"CDO binary:  NOT FOUND - {e}")

    for subdir in ["bin", "lib", "share"]:
        full = os.path.join(pkg_dir, subdir)
        if os.path.isdir(full):
            files = os.listdir(full)
            shown = files[:15]
            suffix = f"... +{len(files) - 15} more" if len(files) > 15 else ""
            print(f"  {subdir}/ ({len(files)} items): {shown}{suffix}")

    # Check auditwheel .libs/ directory
    libs_dir = os.path.join(os.path.dirname(pkg_dir), "skyborn_cdo.libs")
    if os.path.isdir(libs_dir):
        files = os.listdir(libs_dir)
        print(f"  skyborn_cdo.libs/ ({len(files)} items): {files[:10]}")
    else:
        print("  skyborn_cdo.libs/ does not exist")

    # On Linux, show RPATH and ldd for CDO binary
    if sys.platform == "linux":
        import subprocess as _sp

        try:
            cdo_bin = get_cdo_path()
            r = _sp.run(
                ["readelf", "-d", cdo_bin],
                capture_output=True, text=True, timeout=5,
            )
            for line in r.stdout.splitlines():
                if "RPATH" in line or "RUNPATH" in line:
                    print(f"  RPATH: {line.strip()}")
            r = _sp.run(
                ["ldd", cdo_bin],
                capture_output=True, text=True, timeout=5,
            )
            print(f"  ldd output ({len(r.stdout.splitlines())} libs):")
            for line in r.stdout.splitlines()[:30]:
                print(f"    {line.strip()}")
            if "not found" in r.stdout:
                print("  WARNING: some libraries not found!")
        except Exception as e:
            print(f"  Diagnostics error: {e}")

    print()  # blank line before tests

    # ------------------------------------------------------------------
    # Test 1: CDO version
    # ------------------------------------------------------------------
    version = cdo.version()
    print(f"CDO version: {version}")
    passed += 1

    # ------------------------------------------------------------------
    # Test 2: Create data in NetCDF format (doesn't need ecCodes)
    # ------------------------------------------------------------------
    nc_file = os.path.join(tmpdir, "topo.nc")
    try:
        cdo.topo(output=nc_file)
        size = os.path.getsize(nc_file)
        assert size > 0, "NetCDF file is empty"
        print(f"  PASS: NetCDF write ({size} bytes)")
        passed += 1
    except Exception as e:
        print(f"  FAIL: NetCDF write - {e}", file=sys.stderr)
        failed += 1

    # ------------------------------------------------------------------
    # Test 3: Convert NetCDF -> GRIB (tests ecCodes GRIB encoding)
    # ------------------------------------------------------------------
    grb_file = os.path.join(tmpdir, "topo.grb")
    try:
        cdo.copy(input=nc_file, output=grb_file, options="-f grb")
        size = os.path.getsize(grb_file)
        assert size > 0, "GRIB file is empty"
        print(f"  PASS: GRIB write / ecCodes encode ({size} bytes)")
        passed += 1
    except Exception as e:
        print(f"  FAIL: GRIB write / ecCodes encode - {e}", file=sys.stderr)
        failed += 1

    # ------------------------------------------------------------------
    # Test 4: Read GRIB -> NetCDF (tests ecCodes GRIB decoding)
    # ------------------------------------------------------------------
    nc2_file = os.path.join(tmpdir, "roundtrip.nc")
    try:
        if os.path.exists(grb_file) and os.path.getsize(grb_file) > 0:
            cdo.copy(input=grb_file, output=nc2_file)
            size = os.path.getsize(nc2_file)
            assert size > 0, "Roundtrip NetCDF file is empty"
            print(f"  PASS: GRIB read / ecCodes decode ({size} bytes)")
            passed += 1
        else:
            print("  SKIP: GRIB read (no GRIB file from previous step)")
    except Exception as e:
        print(f"  FAIL: GRIB read / ecCodes decode - {e}", file=sys.stderr)
        failed += 1

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print(f"\nResults: {passed} passed, {failed} failed")
    if failed > 0:
        sys.exit(1)
    print("All wheel smoke tests passed!")


if __name__ == "__main__":
    main()
