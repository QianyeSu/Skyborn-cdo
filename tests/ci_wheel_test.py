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
    from skyborn_cdo import Cdo

    cdo = Cdo()
    tmpdir = tempfile.mkdtemp()
    passed = 0
    failed = 0

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
