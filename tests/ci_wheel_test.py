#!/usr/bin/env python
"""
CI wheel smoke test — verifies CDO binary and common operators.

Executed by cibuildwheel after each wheel is installed.
Validates:
  1. CDO binary loads and reports its version
  2. CDO can create data (NetCDF)
  3. CDO can convert NetCDF -> GRIB (ecCodes encode)
  4. CDO can read GRIB -> NetCDF (ecCodes decode)
  5. CDO mergetime operator
  6. CDO remapbil (bilinear interpolation)
  7. CDO fldmean (field statistics)
"""

import os
import sys
import tempfile


# ======================================================================
# Diagnostics — print package layout and dynamic linking info
# ======================================================================

def diagnostics():
    import skyborn_cdo as _pkg
    from skyborn_cdo import get_cdo_path

    pkg_dir = os.path.dirname(_pkg.__file__)
    print(f"Package dir: {pkg_dir}")
    print(f"Python:      {sys.version}")
    print(f"Platform:    {sys.platform}")

    try:
        cdo_bin = get_cdo_path()
        print(f"CDO binary:  {cdo_bin}")
    except FileNotFoundError as e:
        print(f"CDO binary:  NOT FOUND — {e}")
        return

    for subdir in ["bin", "lib", "share"]:
        full = os.path.join(pkg_dir, subdir)
        if os.path.isdir(full):
            items = os.listdir(full)
            shown = items[:15]
            extra = f" ... +{len(items) - 15} more" if len(items) > 15 else ""
            print(f"  {subdir}/ ({len(items)} items): {shown}{extra}")

    # auditwheel .libs/ directory
    libs_dir = os.path.join(os.path.dirname(pkg_dir), "skyborn_cdo.libs")
    if os.path.isdir(libs_dir):
        items = os.listdir(libs_dir)
        print(f"  skyborn_cdo.libs/ ({len(items)} items): {items[:10]}")
    else:
        print("  skyborn_cdo.libs/ does not exist (good — all libs in lib/)")

    if sys.platform != "linux":
        print()
        return

    # ---- Linux-specific diagnostics ----
    import subprocess as _sp
    cdo_bin = get_cdo_path()

    # file type
    try:
        r = _sp.run(["file", cdo_bin], capture_output=True, text=True, timeout=5)
        print(f"  file: {r.stdout.strip()}")
    except Exception:
        pass

    # ELF interpreter
    try:
        r = _sp.run(["readelf", "-l", cdo_bin], capture_output=True, text=True, timeout=5)
        for line in r.stdout.splitlines():
            if "interpreter" in line.lower():
                print(f"  {line.strip()}")
    except Exception:
        pass

    # RPATH / RUNPATH
    try:
        r = _sp.run(["readelf", "-d", cdo_bin], capture_output=True, text=True, timeout=5)
        for line in r.stdout.splitlines():
            if "RPATH" in line or "RUNPATH" in line:
                print(f"  {line.strip()}")
    except Exception:
        pass

    # ldd — capture BOTH stdout and stderr
    try:
        r = _sp.run(["ldd", cdo_bin], capture_output=True, text=True, timeout=10)
        stdout_lines = r.stdout.strip().splitlines() if r.stdout.strip() else []
        print(f"  ldd ({len(stdout_lines)} libs):")
        for line in stdout_lines[:30]:
            print(f"    {line.strip()}")
        not_found = [l for l in stdout_lines if "not found" in l]
        if not_found:
            print(f"  *** WARNING: {len(not_found)} libraries NOT FOUND ***")
        if r.stderr.strip():
            print(f"  ldd stderr: {r.stderr.strip()[:300]}")
        if not stdout_lines and not r.stderr.strip():
            print("  *** ldd produced NO output — binary may be corrupted ***")
    except Exception as e:
        print(f"  ldd error: {e}")

    # LD_LIBRARY_PATH that will be used at runtime
    from skyborn_cdo._cdo_binary import get_bundled_env
    env = get_bundled_env()
    print(f"  LD_LIBRARY_PATH: {env.get('LD_LIBRARY_PATH', '(not set)')[:200]}")

    # LD_DEBUG=libs — trace library loading for CDO --version
    try:
        env_dbg = env.copy()
        env_dbg["LD_DEBUG"] = "libs"
        r = _sp.run(
            [cdo_bin, "--version"],
            capture_output=True, text=True, timeout=10,
            env=env_dbg, stdin=_sp.DEVNULL,
        )
        dbg = r.stderr.splitlines()
        print(f"  LD_DEBUG ({len(dbg)} lines, rc={r.returncode}):")
        for line in dbg[:15]:
            print(f"    {line}")
        # Highlight errors
        for i, line in enumerate(dbg):
            if i >= 15 and ("error" in line.lower() or "not found" in line.lower()):
                print(f"    ...[{i}] {line}")
    except Exception as e:
        print(f"  LD_DEBUG error: {e}")

    print()  # blank line before tests


# ======================================================================
# Test runner
# ======================================================================

def main():
    from skyborn_cdo import Cdo

    diagnostics()

    cdo = Cdo()
    tmpdir = tempfile.mkdtemp()
    passed = 0
    failed = 0

    def run_test(name, fn):
        nonlocal passed, failed
        try:
            fn()
            passed += 1
        except Exception as e:
            print(f"  FAIL: {name} — {e}", file=sys.stderr)
            failed += 1

    # ---- Test 1: CDO version ----
    def test_version():
        v = cdo.version()
        print(f"CDO version: {v}")
        if not v or "Error" in str(v):
            raise RuntimeError(f"CDO version returned empty or error: '{v}'")
    run_test("CDO version", test_version)

    # ---- Test 2: NetCDF write (topo) ----
    nc_file = os.path.join(tmpdir, "topo.nc")
    def test_netcdf_write():
        cdo.topo(output=nc_file)
        sz = os.path.getsize(nc_file)
        assert sz > 0, "NetCDF file is empty"
        print(f"  PASS: NetCDF write ({sz:,} bytes)")
    run_test("NetCDF write", test_netcdf_write)

    # ---- Test 3: GRIB encode (NetCDF -> GRIB) ----
    grb_file = os.path.join(tmpdir, "topo.grb")
    def test_grib_encode():
        if not os.path.exists(nc_file) or os.path.getsize(nc_file) == 0:
            raise RuntimeError("SKIP — no NetCDF source")
        cdo.copy(input=nc_file, output=grb_file, options="-f grb")
        sz = os.path.getsize(grb_file)
        assert sz > 0, "GRIB file is empty"
        print(f"  PASS: GRIB encode ({sz:,} bytes)")
    run_test("GRIB encode (ecCodes write)", test_grib_encode)

    # ---- Test 4: GRIB decode (GRIB -> NetCDF) ----
    nc2_file = os.path.join(tmpdir, "roundtrip.nc")
    def test_grib_decode():
        if not os.path.exists(grb_file) or os.path.getsize(grb_file) == 0:
            raise RuntimeError("SKIP — no GRIB source")
        cdo.copy(input=grb_file, output=nc2_file)
        sz = os.path.getsize(nc2_file)
        assert sz > 0, "Roundtrip NetCDF is empty"
        print(f"  PASS: GRIB decode ({sz:,} bytes)")
    run_test("GRIB decode (ecCodes read)", test_grib_decode)

    # ---- Test 5: mergetime ----
    def test_mergetime():
        if not os.path.exists(nc_file) or os.path.getsize(nc_file) == 0:
            raise RuntimeError("SKIP — no source file")
        t1 = os.path.join(tmpdir, "t1.nc")
        t2 = os.path.join(tmpdir, "t2.nc")
        merged = os.path.join(tmpdir, "merged.nc")
        cdo.settaxis("2020-01-01,12:00:00,1day", input=nc_file, output=t1)
        cdo.settaxis("2020-01-02,12:00:00,1day", input=nc_file, output=t2)
        cdo.mergetime(input=f"{t1} {t2}", output=merged)
        sz = os.path.getsize(merged)
        assert sz > 0, "Merged file is empty"
        nsteps = cdo.ntime(input=merged)
        n = int(str(nsteps).strip()) if nsteps else 0
        assert n == 2, f"Expected 2 timesteps, got {n}"
        print(f"  PASS: mergetime ({sz:,} bytes, {n} timesteps)")
    run_test("mergetime", test_mergetime)

    # ---- Test 6: remapbil (bilinear interpolation) ----
    def test_remapbil():
        if not os.path.exists(nc_file) or os.path.getsize(nc_file) == 0:
            raise RuntimeError("SKIP — no source file")
        remap_out = os.path.join(tmpdir, "remapped.nc")
        cdo.remapbil("r10x10", input=nc_file, output=remap_out)
        sz = os.path.getsize(remap_out)
        assert sz > 0, "Remapped file is empty"
        print(f"  PASS: remapbil / bilinear interpolation ({sz:,} bytes)")
    run_test("remapbil (bilinear interpolation)", test_remapbil)

    # ---- Test 7: fldmean (field mean) ----
    def test_fldmean():
        if not os.path.exists(nc_file) or os.path.getsize(nc_file) == 0:
            raise RuntimeError("SKIP — no source file")
        fm_out = os.path.join(tmpdir, "fldmean.nc")
        cdo.fldmean(input=nc_file, output=fm_out)
        sz = os.path.getsize(fm_out)
        assert sz > 0, "Fldmean file is empty"
        print(f"  PASS: fldmean / field statistics ({sz:,} bytes)")
    run_test("fldmean (field statistics)", test_fldmean)

    # ---- Summary ----
    total = passed + failed
    print(f"\nResults: {passed}/{total} passed, {failed} failed")
    if failed > 0:
        sys.exit(1)
    print("All wheel smoke tests passed!")


if __name__ == "__main__":
    main()
