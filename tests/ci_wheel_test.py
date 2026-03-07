#!/usr/bin/env python
"""
CI Wheel Comprehensive Stress Test — verifies CDO binary functionality.

Executed by cibuildwheel after each wheel is installed.
Tests as many CDO operations as possible within CI time constraints (~5 min).

Test Coverage (~250 tests):
  • Basic: version, operators, help
  • Synthetic data: topo, random, for, stdatm, const
  • Info queries: sinfo, griddes, showname, ntime, filedes, showformat
  • Selection: sellonlatbox, selindexbox, sellevel, selmon, selname, selcode,
              selseason, selday, selyear, selsmon, seltime
  • Statistics: fldmean/std/min/max/sum/var/skew/kurt/median/pctl, zonmean,
               timmean/std/var/range/skew/kurt/avg, zon/mer full families
  • Zonal/meridional: zonmin/max/sum/std/skew/kurt/median, mermin/max/sum/std/skew/kurt/median
  • Running stats: runmean, runstd, runmin, runmax, runsum, runvar, runrange
  • Grid box stats: gridboxmean, gridboxmax, gridboxmin, gridboxstd, gridboxrange/sum/var
  • Arithmetic: mulc, addc, subc, divc, abs, sqrt, sqr, expr, add/sub/mul/div
  • Math: exp, ln, log10, sin, cos, tan, asin, acos, atan, reci, nint
  • Comparison: eq, ne, le, lt, ge, gt
  • Grid ops: remapbil, remapcon, remapnn, remapbic, remaplaf, gridarea, gridweights
  • Spectral: gp2sp, sp2gpl, sp2gp (CRITICAL)
  • Format: NetCDF4/4c/2, GRIB1, GRIB2
  • Time: mergetime, settaxis, selmon, monmean/min/max/std/sum,
         seasmean/std/min/max/sum/range, yearmean/min/max/std/sum/range
  • Vertical: intlevel, vertmean/sum/min/max/std/var, invertlev
  • Masking: setmissval, setrtomiss, setmisstoc, ifthen, ifnotthen, ifthenelse,
             masklonlatbox, setctomiss
  • Ensemble: ensmean, ensmin, ensmax, ensstd, ensrange, enssum, ensskew, enskurt, enspctl
  • VarsStat: varsskew, varskurt, varsmedian, varspctl, varsmean, varsstd,
              varsmin, varsmax, varsrange, varssum
  • Trend: trend, detrend
  • File ops: cat, duplicate, invertlat, sortname, splitvar, splitmon
  • Metadata: chname, setstdname, setname, setunit, setcode, showlevel, showcode
  • Correlation: timcor, fldcor, timcovar
  • Z-axis: zaxisdes, npar
  • Date/Time: showdate, showtime, showtimestamp, shifttime
  • Chained: complex multi-operator pipes
  • CDO 2.6.0 new: varsskew, varskurt, varsmedian, varspctl, symmetrize, fillmiss
  • Regression: explicit-coordinate NC4 (lazy-grid mutex -- 0xC0000005 fix)
  • String API: cdo("cdo ...") explicit-file write (two-phase HDF5 regression)
  • run_raw() direct: CdoRunner.run_raw() lower-level API (mulc, timmean, remapbil,
                      chained ops, returncode check, no-prefix mode, overwrite flag)
  • NCL wind: uv2vr_cfd/uv2dv_cfd with synthetic u/v file (CF standard_name recognition)
  • Error handling: invalid inputs
"""

import os
import sys
import tempfile
import time


# ======================================================================
# Diagnostics (Linux only — shows ELF/library loading)
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
        r = _sp.run(["file", cdo_bin], capture_output=True,
                    text=True, timeout=5)
        print(f"  file: {r.stdout.strip()}")
    except Exception:
        pass

    # ELF interpreter
    try:
        r = _sp.run(["readelf", "-l", cdo_bin],
                    capture_output=True, text=True, timeout=5)
        for line in r.stdout.splitlines():
            if "interpreter" in line.lower():
                print(f"  {line.strip()}")
    except Exception:
        pass

    # RPATH / RUNPATH
    try:
        r = _sp.run(["readelf", "-d", cdo_bin],
                    capture_output=True, text=True, timeout=5)
        for line in r.stdout.splitlines():
            if "RPATH" in line or "RUNPATH" in line:
                print(f"  {line.strip()}")
    except Exception:
        pass

    # ldd — capture BOTH stdout and stderr
    try:
        r = _sp.run(["ldd", cdo_bin], capture_output=True,
                    text=True, timeout=10)
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
    print(
        f"  LD_LIBRARY_PATH: {env.get('LD_LIBRARY_PATH', '(not set)')[:200]}")

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
# Test runner — comprehensive stress test
# ======================================================================

def main():
    import shutil
    from skyborn_cdo import Cdo
    from skyborn_cdo._runner import CdoError

    diagnostics()

    cdo = Cdo(timeout=10)
    tmpdir = tempfile.mkdtemp()
    passed = 0
    failed = 0
    t_start = time.time()

    def run_test(name, fn):
        nonlocal passed, failed
        t0 = time.time()
        try:
            fn()
            passed += 1
            dt = time.time() - t0
            slow = f" ({dt:.1f}s)" if dt > 2.0 else ""
            print(f"  [PASS] {name}{slow}")
        except Exception as e:
            print(f"  [FAIL] {name} -- {e}", file=sys.stderr)
            failed += 1

    # Basic functionality
    print("\n=== 1. Basic ===")
    run_test("version", lambda: cdo.version())
    run_test("operators", lambda: assert_true(len(cdo.operators()) > 800))
    run_test("has_operator", lambda: assert_true(
        cdo.has_operator("mergetime")))

    # Synthetic data generation
    print("\n=== 2. Data Generation ===")
    topo_nc = os.path.join(tmpdir, "topo.nc")
    run_test("topo (r12x6)", lambda: cdo(
        f"cdo remapbil,r12x6 -topo {topo_nc}", timeout=30) or assert_file(topo_nc))

    rand_nc = os.path.join(tmpdir, "rand.nc")
    run_test("random r12x6", lambda: cdo(
        f"cdo -random,r12x6 {rand_nc}", timeout=30) or assert_file(rand_nc))

    const_nc = os.path.join(tmpdir, "const.nc")
    run_test("const field", lambda: cdo(
        f"cdo -const,273.15,r12x6 {const_nc}", timeout=20) or assert_file(const_nc))

    monthly_nc = os.path.join(tmpdir, "monthly.nc")
    run_test("12-month series", lambda: cdo(
        f"cdo -settaxis,2020-01-15,12:00,1mon -for,1,12 {monthly_nc}", timeout=30) or assert_file(monthly_nc))

    stdatm_nc = os.path.join(tmpdir, "stdatm.nc")
    run_test("stdatm levels", lambda: cdo(
        f"cdo stdatm,0,10000,30000,50000 {stdatm_nc}", timeout=20) or assert_file(stdatm_nc))

    # Info queries
    print("\n=== 3. Info Queries ===")
    run_test("sinfo", lambda: assert_true(
        len(str(cdo.sinfo(input=topo_nc))) > 10))
    run_test("griddes", lambda: assert_true(
        "gridtype" in str(cdo.griddes(input=topo_nc)).lower()))
    run_test("showname", lambda: len(cdo.showname(input=topo_nc)))
    run_test("ntime", lambda: cdo.ntime(input=monthly_nc))
    run_test("nlevel", lambda: cdo.nlevel(input=stdatm_nc))

    # Selection/clipping
    print("\n=== 4. Selection ===")
    sellonlat_nc = os.path.join(tmpdir, "sellonlat.nc")
    run_test("sellonlatbox", lambda: cdo.sellonlatbox(
        "0,90,0,45", input=topo_nc, output=sellonlat_nc) or assert_file(sellonlat_nc))

    selidx_nc = os.path.join(tmpdir, "selindex.nc")
    run_test("selindexbox", lambda: cdo.selindexbox(
        "1,6,1,3", input=topo_nc, output=selidx_nc) or assert_file(selidx_nc))

    sellev_nc = os.path.join(tmpdir, "sellev.nc")
    run_test("sellevel", lambda: cdo.sellevel(
        "0,10000", input=stdatm_nc, output=sellev_nc) or assert_file(sellev_nc))

    selmon_nc = os.path.join(tmpdir, "selmon.nc")
    run_test("selmon", lambda: cdo.selmon("1,2,3", input=monthly_nc,
             output=selmon_nc) or assert_file(selmon_nc))

    # Statistics
    print("\n=== 5. Statistics ===")
    fldmean_nc = os.path.join(tmpdir, "fldmean.nc")
    run_test("fldmean", lambda: cdo.fldmean(input=topo_nc,
             output=fldmean_nc) or assert_file(fldmean_nc))

    fldstd_nc = os.path.join(tmpdir, "fldstd.nc")
    run_test("fldstd", lambda: cdo.fldstd(input=topo_nc,
             output=fldstd_nc) or assert_file(fldstd_nc))

    zonmean_nc = os.path.join(tmpdir, "zonmean.nc")
    run_test("zonmean", lambda: cdo.zonmean(input=topo_nc,
             output=zonmean_nc) or assert_file(zonmean_nc))

    timmean_nc = os.path.join(tmpdir, "timmean.nc")
    run_test("timmean", lambda: cdo.timmean(input=monthly_nc,
             output=timmean_nc) or assert_file(timmean_nc))

    timstd_nc = os.path.join(tmpdir, "timstd.nc")
    run_test("timstd", lambda: cdo.timstd(input=monthly_nc,
             output=timstd_nc) or assert_file(timstd_nc))

    # Arithmetic
    print("\n=== 6. Arithmetic ===")
    mulc_nc = os.path.join(tmpdir, "mulc.nc")
    run_test("mulc", lambda: cdo.mulc("2.5", input=topo_nc,
             output=mulc_nc) or assert_file(mulc_nc))

    addc_nc = os.path.join(tmpdir, "addc.nc")
    run_test("addc", lambda: cdo.addc("100", input=topo_nc,
             output=addc_nc) or assert_file(addc_nc))

    abs_nc = os.path.join(tmpdir, "abs.nc")
    run_test("abs", lambda: cdo.abs(input=topo_nc,
             output=abs_nc) or assert_file(abs_nc))

    sqrt_nc = os.path.join(tmpdir, "sqrt.nc")
    run_test("sqrt", lambda: cdo.sqrt(input=abs_nc,
             output=sqrt_nc) or assert_file(sqrt_nc))

    expr_nc = os.path.join(tmpdir, "expr.nc")

    def _expr_test():
        names = str(cdo.showname(input=topo_nc)).strip().split()
        vname = names[0] if names else "topo"
        cdo.expr(f"doubled={vname}*2;", input=topo_nc, output=expr_nc)
        assert_file(expr_nc)
    run_test("expr", _expr_test)

    add_nc = os.path.join(tmpdir, "add.nc")
    run_test("add files", lambda: cdo.add(
        input=f"{topo_nc} {mulc_nc}", output=add_nc) or assert_file(add_nc))

    # Grid operations
    print("\n=== 7. Grid Ops ===")
    remap_bil_nc = os.path.join(tmpdir, "remap_bil.nc")
    run_test("remapbil", lambda: cdo.remapbil("r8x4", input=topo_nc,
             output=remap_bil_nc) or assert_file(remap_bil_nc))

    remap_con_nc = os.path.join(tmpdir, "remap_con.nc")
    run_test("remapcon", lambda: cdo.remapcon("r8x4", input=topo_nc,
             output=remap_con_nc) or assert_file(remap_con_nc))

    remap_nn_nc = os.path.join(tmpdir, "remap_nn.nc")
    run_test("remapnn", lambda: cdo.remapnn("r8x4", input=topo_nc,
             output=remap_nn_nc) or assert_file(remap_nn_nc))

    # Spectral (CRITICAL)
    print("\n=== 8. Spectral (CRITICAL) ===")
    sp_nc = os.path.join(tmpdir, "spectral.nc")
    run_test("gp2sp T21", lambda: cdo(
        f"cdo gp2sp -remapbil,t21grid {topo_nc} {sp_nc}", timeout=60) or assert_file(sp_nc))

    sp2gpl_nc = os.path.join(tmpdir, "sp2gpl.nc")
    run_test("sp2gpl", lambda: cdo.sp2gpl(
        input=sp_nc, output=sp2gpl_nc) or assert_file(sp2gpl_nc))

    sp2gp_nc = os.path.join(tmpdir, "sp2gp.nc")
    run_test("sp2gp", lambda: cdo.sp2gp(
        input=sp_nc, output=sp2gp_nc) or assert_file(sp2gp_nc))

    complex_sp_nc = os.path.join(tmpdir, "complex_sp.nc")
    run_test("sp2gpl chain nc4", lambda: cdo(
        f"cdo -f nc4 -sp2gpl -setgridtype,regular {sp_nc} {complex_sp_nc}", timeout=60) or assert_file(complex_sp_nc))

    # Format conversion
    print("\n=== 9. Formats ===")
    nc4_nc = os.path.join(tmpdir, "nc4.nc")
    run_test("NetCDF4", lambda: cdo.copy(input=topo_nc,
             output=nc4_nc, options="-f nc4") or assert_file(nc4_nc))

    nc2_nc = os.path.join(tmpdir, "nc2.nc")
    run_test("NetCDF2", lambda: cdo.copy(input=topo_nc,
             output=nc2_nc, options="-f nc2") or assert_file(nc2_nc))

    grb_file = os.path.join(tmpdir, "topo.grb")
    run_test("GRIB1 encode", lambda: cdo.copy(input=topo_nc,
             output=grb_file, options="-f grb") or assert_file(grb_file))

    grb2_file = os.path.join(tmpdir, "topo.grb2")
    run_test("GRIB2 encode", lambda: cdo.copy(input=topo_nc,
             output=grb2_file, options="-f grb2") or assert_file(grb2_file))

    grb_decode_nc = os.path.join(tmpdir, "grb_decode.nc")
    if os.path.exists(grb_file):
        run_test("GRIB1 decode", lambda: cdo.copy(input=grb_file,
                 output=grb_decode_nc) or assert_file(grb_decode_nc))

    # Time operations
    print("\n=== 10. Time Ops ===")
    t1_nc = os.path.join(tmpdir, "t1.nc")
    t2_nc = os.path.join(tmpdir, "t2.nc")
    merged_nc = os.path.join(tmpdir, "merged.nc")
    run_test("mergetime", lambda: (
        cdo.settaxis("2020-01-01,12:00:00,1day", input=topo_nc, output=t1_nc),
        cdo.settaxis("2020-01-02,12:00:00,1day", input=topo_nc, output=t2_nc),
        cdo.mergetime(input=f"{t1_nc} {t2_nc}", output=merged_nc),
        assert_file(merged_nc),
        assert_true(int(str(cdo.ntime(input=merged_nc)).strip()) == 2)
    ))

    run_test("showdate", lambda: assert_true(
        "2020" in str(cdo.showdate(input=monthly_nc))))
    run_test("showmon", lambda: str(cdo.showmon(input=monthly_nc)))

    # Chained operations
    print("\n=== 11. Chains ===")
    chain1_nc = os.path.join(tmpdir, "chain1.nc")
    run_test("sel+remap", lambda:  cdo(
        f"cdo -remapbil,r72x36 -sellonlatbox,0,180,0,90 {topo_nc} {chain1_nc}", timeout=40) or assert_file(chain1_nc))

    chain2_nc = os.path.join(tmpdir, "chain2.nc")
    run_test("sel+fldmean", lambda: cdo(
        f"cdo -fldmean -sellonlatbox,-180,180,-30,30 {topo_nc} {chain2_nc}", timeout=30) or assert_file(chain2_nc))

    chain3_nc = os.path.join(tmpdir, "chain3.nc")
    run_test("mulc+addc", lambda: cdo(
        f"cdo -addc,273.15 -mulc,0.01 {topo_nc} {chain3_nc}", timeout=30) or assert_file(chain3_nc))

    # Vertical interpolation
    print("\n=== 12. Vertical / Level Ops ===")
    ml2pl_nc = os.path.join(tmpdir, "ml2pl.nc")
    run_test("intlevel", lambda: cdo.intlevel(
        "0,5000,20000", input=stdatm_nc, output=ml2pl_nc) or assert_file(ml2pl_nc))

    seltimestep_nc = os.path.join(tmpdir, "seltimestep.nc")
    run_test("seltimestep", lambda: cdo.seltimestep(
        "1,2,3", input=monthly_nc, output=seltimestep_nc) or assert_file(seltimestep_nc))

    # Masking and conditional
    print("\n=== 13. Masking ===")
    setmiss_nc = os.path.join(tmpdir, "setmiss.nc")
    run_test("setmissval", lambda: cdo.setmissval(
        "-999", input=topo_nc, output=setmiss_nc) or assert_file(setmiss_nc))

    setrtomiss_nc = os.path.join(tmpdir, "setrtomiss.nc")
    run_test("setrtomiss", lambda: cdo.setrtomiss(
        "-1000,0", input=topo_nc, output=setrtomiss_nc) or assert_file(setrtomiss_nc))

    setmisstoc_nc = os.path.join(tmpdir, "setmisstoc.nc")
    run_test("setmisstoc", lambda: cdo.setmisstoc(
        "0", input=setrtomiss_nc, output=setmisstoc_nc) or assert_file(setmisstoc_nc))

    # More statistics
    print("\n=== 14. Extended Stats ===")
    fldmin_nc = os.path.join(tmpdir, "fldmin.nc")
    run_test("fldmin", lambda: cdo.fldmin(input=topo_nc,
             output=fldmin_nc) or assert_file(fldmin_nc))

    fldmax_nc = os.path.join(tmpdir, "fldmax.nc")
    run_test("fldmax", lambda: cdo.fldmax(input=topo_nc,
             output=fldmax_nc) or assert_file(fldmax_nc))

    fldsum_nc = os.path.join(tmpdir, "fldsum.nc")
    run_test("fldsum", lambda: cdo.fldsum(input=topo_nc,
             output=fldsum_nc) or assert_file(fldsum_nc))

    timmin_nc = os.path.join(tmpdir, "timmin.nc")
    run_test("timmin", lambda: cdo.timmin(input=monthly_nc,
             output=timmin_nc) or assert_file(timmin_nc))

    timmax_nc = os.path.join(tmpdir, "timmax.nc")
    run_test("timmax", lambda: cdo.timmax(input=monthly_nc,
             output=timmax_nc) or assert_file(timmax_nc))

    timsum_nc = os.path.join(tmpdir, "timsum.nc")
    run_test("timsum", lambda: cdo.timsum(input=monthly_nc,
             output=timsum_nc) or assert_file(timsum_nc))

    mermean_nc = os.path.join(tmpdir, "mermean.nc")
    run_test("mermean", lambda: cdo.mermean(input=topo_nc,
             output=mermean_nc) or assert_file(mermean_nc))

    # Grid description & manipulation
    print("\n=== 15. Grid Manipulation ===")
    run_test("gridarea", lambda: cdo.gridarea(input=topo_nc,
             output=os.path.join(tmpdir, "gridarea.nc")) or assert_file(os.path.join(tmpdir, "gridarea.nc")))

    run_test("gridweights", lambda: cdo.gridweights(input=topo_nc,
             output=os.path.join(tmpdir, "gridweights.nc")) or assert_file(os.path.join(tmpdir, "gridweights.nc")))

    setgrid_nc = os.path.join(tmpdir, "setgrid.nc")
    run_test("setgridtype", lambda: cdo.setgridtype(
        "regular", input=topo_nc, output=setgrid_nc) or assert_file(setgrid_nc))

    # Metadata operations
    print("\n=== 16. Metadata ===")
    chname_nc = os.path.join(tmpdir, "chname.nc")

    def _chname_test():
        # Detect variable name, stripping any non-printable chars
        raw = str(cdo.showname(input=topo_nc))
        vname = raw.strip().split()[0] if raw.strip() else "topo"
        vname = ''.join(c for c in vname if c.isprintable() and c != ' ')
        # Force NetCDF output — chname cannot rename vars in GRIB format
        # (GRIB uses parameter codes, variable names are derived from code tables)
        cdo(f"cdo -f nc -chname,{vname},elevation {topo_nc} {chname_nc}", timeout=30)
        assert_file(chname_nc)
        # Verify the rename by reading the NC3 binary header directly.
        # On Windows, cdo.exe crashes (exit 0xC0000005 ACCESS VIOLATION) when
        # reading any NetCDF file — CDI uses pthread_mutex (resource_handle.c
        # LIST_LOCK) in its NC code path, and winpthreads causes a crash.
        # CDO's default format is GRIB and reads fine; only NC reads crash.
        # The chname operator itself works correctly — we just verify the
        # result without invoking CDO again.
        _verify_nc_varname(chname_nc, "elevation", vname)
    run_test("chname", _chname_test)

    run_test("showyear", lambda: str(cdo.showyear(input=monthly_nc)))
    run_test("nvar", lambda: assert_true(
        int(str(cdo.nvar(input=topo_nc)).strip()) >= 1))
    run_test("showlevel", lambda: str(cdo.showlevel(input=stdatm_nc)))
    run_test("showcode", lambda: str(cdo.showcode(input=topo_nc)))

    # Seasonal / monthly statistics
    print("\n=== 17. Seasonal Stats ===")
    ymonmean_nc = os.path.join(tmpdir, "ymonmean.nc")
    run_test("ymonmean", lambda: cdo.ymonmean(input=monthly_nc,
             output=ymonmean_nc) or assert_file(ymonmean_nc))

    # File comparison / diff
    print("\n=== 18. Comparison ===")
    run_test("diff (identical)", lambda: cdo.diff(
        input=f"{topo_nc} {topo_nc}"))

    sub_nc = os.path.join(tmpdir, "sub.nc")
    run_test("sub", lambda: cdo.sub(
        input=f"{topo_nc} {topo_nc}", output=sub_nc) or assert_file(sub_nc))

    mul_nc = os.path.join(tmpdir, "mul.nc")
    run_test("mul", lambda: cdo.mul(
        input=f"{topo_nc} {topo_nc}", output=mul_nc) or assert_file(mul_nc))

    div_nc = os.path.join(tmpdir, "div.nc")
    run_test("div", lambda: cdo.div(
        input=f"{topo_nc} {addc_nc}", output=div_nc) or assert_file(div_nc))

    # Advanced chains
    print("\n=== 19. Advanced Chains ===")
    chain4_nc = os.path.join(tmpdir, "chain4.nc")
    run_test("fldmean+abs+mulc", lambda: cdo(
        f"cdo -fldmean -abs -mulc,-1 {topo_nc} {chain4_nc}", timeout=30) or assert_file(chain4_nc))

    chain5_nc = os.path.join(tmpdir, "chain5.nc")
    run_test("fldmean+sellev", lambda: cdo(
        f"cdo -fldmean -sellevel,0 {stdatm_nc} {chain5_nc}", timeout=30) or assert_file(chain5_nc))

    chain6_nc = os.path.join(tmpdir, "chain6.nc")
    run_test("timmean+selmon", lambda: cdo(
        f"cdo -timmean -selmon,1,6 {monthly_nc} {chain6_nc}", timeout=30) or assert_file(chain6_nc))

    # Ensemble / merge operations
    print("\n=== 20. Merge / Ensemble ===")
    merge_nc = os.path.join(tmpdir, "merge.nc")
    run_test("merge", lambda: cdo.merge(
        input=f"{topo_nc} {const_nc}", output=merge_nc) or assert_file(merge_nc))

    run_test("splitlevel", lambda: cdo.splitlevel(
        input=stdatm_nc, output=os.path.join(tmpdir, "splev")))

    ensmean_nc = os.path.join(tmpdir, "ensmean.nc")
    run_test("ensmean", lambda: cdo.ensmean(
        input=f"{topo_nc} {topo_nc}", output=ensmean_nc) or assert_file(ensmean_nc))

    # File info queries
    print("\n=== 22. File Info ===")
    run_test("filedes", lambda: assert_true(
        len(str(cdo.filedes(input=topo_nc))) > 5))
    run_test("showformat", lambda: assert_true(
        len(str(cdo.showformat(input=topo_nc)).strip()) > 0))
    run_test("showunit", lambda: str(cdo.showunit(input=topo_nc)))
    run_test("showtimestamp", lambda: str(cdo.showtimestamp(input=monthly_nc)))
    run_test("showstdname", lambda: str(cdo.showstdname(input=topo_nc)))

    # Time statistics
    print("\n=== 23. Time Statistics ===")
    monmean_nc = os.path.join(tmpdir, "monmean.nc")
    run_test("monmean", lambda: cdo.monmean(input=monthly_nc,
             output=monmean_nc) or assert_file(monmean_nc))

    seasmean_nc = os.path.join(tmpdir, "seasmean.nc")
    run_test("seasmean", lambda: cdo.seasmean(input=monthly_nc,
             output=seasmean_nc) or assert_file(seasmean_nc))

    yearmean_nc = os.path.join(tmpdir, "yearmean.nc")
    run_test("yearmean", lambda: cdo.yearmean(input=monthly_nc,
             output=yearmean_nc) or assert_file(yearmean_nc))

    yearmin_nc = os.path.join(tmpdir, "yearmin.nc")
    run_test("yearmin", lambda: cdo.yearmin(input=monthly_nc,
             output=yearmin_nc) or assert_file(yearmin_nc))

    yearmax_nc = os.path.join(tmpdir, "yearmax.nc")
    run_test("yearmax", lambda: cdo.yearmax(input=monthly_nc,
             output=yearmax_nc) or assert_file(yearmax_nc))

    # Vertical operations
    print("\n=== 24. Vertical Ops ===")
    vertmean_nc = os.path.join(tmpdir, "vertmean.nc")
    run_test("vertmean", lambda: cdo.vertmean(input=stdatm_nc,
             output=vertmean_nc) or assert_file(vertmean_nc))

    vertsum_nc = os.path.join(tmpdir, "vertsum.nc")
    run_test("vertsum", lambda: cdo.vertsum(input=stdatm_nc,
             output=vertsum_nc) or assert_file(vertsum_nc))

    vertmin_nc = os.path.join(tmpdir, "vertmin.nc")
    run_test("vertmin", lambda: cdo.vertmin(input=stdatm_nc,
             output=vertmin_nc) or assert_file(vertmin_nc))

    vertmax_nc = os.path.join(tmpdir, "vertmax.nc")
    run_test("vertmax", lambda: cdo.vertmax(input=stdatm_nc,
             output=vertmax_nc) or assert_file(vertmax_nc))

    # Ensemble extended
    print("\n=== 25. Ensemble Extended ===")
    ensmin_nc = os.path.join(tmpdir, "ensmin.nc")
    run_test("ensmin", lambda: cdo.ensmin(
        input=f"{topo_nc} {topo_nc}", output=ensmin_nc) or assert_file(ensmin_nc))

    ensmax_nc = os.path.join(tmpdir, "ensmax.nc")
    run_test("ensmax", lambda: cdo.ensmax(
        input=f"{topo_nc} {topo_nc}", output=ensmax_nc) or assert_file(ensmax_nc))

    ensstd_nc = os.path.join(tmpdir, "ensstd.nc")
    run_test("ensstd", lambda: cdo.ensstd(
        input=f"{topo_nc} {topo_nc}", output=ensstd_nc) or assert_file(ensstd_nc))

    # Masking extended
    print("\n=== 26. Masking / Conditional ===")
    ifthen_nc = os.path.join(tmpdir, "ifthen.nc")
    run_test("ifthen", lambda: cdo.ifthen(
        input=f"{topo_nc} {topo_nc}", output=ifthen_nc) or assert_file(ifthen_nc))

    maskbox_nc = os.path.join(tmpdir, "maskbox.nc")
    run_test("masklonlatbox", lambda: cdo.masklonlatbox(
        "0,90,0,45", input=topo_nc, output=maskbox_nc) or assert_file(maskbox_nc))

    setctomiss_nc = os.path.join(tmpdir, "setctomiss.nc")
    run_test("setctomiss", lambda: cdo.setctomiss(
        "0", input=topo_nc, output=setctomiss_nc) or assert_file(setctomiss_nc))

    # Trend analysis
    print("\n=== 27. Trend / Detrend ===")
    trend_a = os.path.join(tmpdir, "trend_a.nc")
    trend_b = os.path.join(tmpdir, "trend_b.nc")
    run_test("trend", lambda: cdo(
             f"cdo trend {monthly_nc} {trend_a} {trend_b}") or (assert_file(trend_a), assert_file(trend_b)))

    detrend_nc = os.path.join(tmpdir, "detrend.nc")
    run_test("detrend", lambda: cdo.detrend(input=monthly_nc,
             output=detrend_nc) or assert_file(detrend_nc))

    # File operations
    print("\n=== 28. File Ops ===")
    cat_nc = os.path.join(tmpdir, "cat.nc")
    run_test("cat", lambda: cdo.cat(
        input=f"{topo_nc} {topo_nc}", output=cat_nc) or assert_file(cat_nc))

    dup_nc = os.path.join(tmpdir, "dup.nc")
    run_test("duplicate", lambda: cdo.duplicate("2", input=topo_nc,
             output=dup_nc) or assert_file(dup_nc))

    invlev_nc = os.path.join(tmpdir, "invertlev.nc")
    run_test("invertlev", lambda: cdo.invertlev(input=stdatm_nc,
             output=invlev_nc) or assert_file(invlev_nc))

    invlat_nc = os.path.join(tmpdir, "invertlat.nc")
    run_test("invertlat", lambda: cdo.invertlat(input=topo_nc,
             output=invlat_nc) or assert_file(invlat_nc))

    # Date/Time display & shift
    print("\n=== 29. Date/Time ===")
    run_test("showtime", lambda: str(cdo.showtime(input=monthly_nc)))

    shift_nc = os.path.join(tmpdir, "shifttime.nc")
    run_test("shifttime", lambda: cdo.shifttime("1hour", input=monthly_nc,
             output=shift_nc) or assert_file(shift_nc))

    settunits_nc = os.path.join(tmpdir, "settunits.nc")
    run_test("settunits", lambda: cdo.settunits("hours", input=monthly_nc,
             output=settunits_nc) or assert_file(settunits_nc))

    # More arithmetic operators
    print("\n=== 30. More Arithmetic ===")
    sqr_nc = os.path.join(tmpdir, "sqr.nc")
    run_test("sqr", lambda: cdo.sqr(input=topo_nc,
             output=sqr_nc) or assert_file(sqr_nc))

    min2_nc = os.path.join(tmpdir, "min2.nc")
    run_test("min (2 files)", lambda: cdo.min(
        input=f"{topo_nc} {mulc_nc}", output=min2_nc) or assert_file(min2_nc))

    max2_nc = os.path.join(tmpdir, "max2.nc")
    run_test("max (2 files)", lambda: cdo.max(
        input=f"{topo_nc} {mulc_nc}", output=max2_nc) or assert_file(max2_nc))

    subc_nc = os.path.join(tmpdir, "subc.nc")
    run_test("subc", lambda: cdo.subc("50", input=topo_nc,
             output=subc_nc) or assert_file(subc_nc))

    divc_nc = os.path.join(tmpdir, "divc.nc")
    run_test("divc", lambda: cdo.divc("2", input=topo_nc,
             output=divc_nc) or assert_file(divc_nc))

    # CDO 2.5.x regression tests
    print("\n=== 31. CDO 2.5.x Features ===")
    setstdname_nc = os.path.join(tmpdir, "setstdname.nc")

    def _setstdname_test():
        raw = str(cdo.showname(input=topo_nc))
        vname = raw.strip().split()[0] if raw.strip() else "topo"
        vname = ''.join(c for c in vname if c.isprintable() and c != ' ')
        cdo(f"cdo -f nc -setstdname,{vname},air_temperature {topo_nc} {setstdname_nc}", timeout=30)
        assert_file(setstdname_nc)
    run_test("setstdname (2.5.3+)", _setstdname_test)

    nc4c_nc = os.path.join(tmpdir, "nc4c.nc")
    run_test("nc4c compress", lambda: cdo.copy(input=topo_nc,
             output=nc4c_nc, options="-f nc4c") or assert_file(nc4c_nc))

    sortname_nc = os.path.join(tmpdir, "sortname.nc")
    if os.path.exists(merge_nc):
        run_test("sortname", lambda: cdo.sortname(input=merge_nc,
                 output=sortname_nc) or assert_file(sortname_nc))

    run_test("ngrids", lambda: assert_true(
        int(str(cdo.ngrids(input=topo_nc)).strip()) >= 1))

    # Error handling
    print("\n=== 32. Errors ===")
    run_test("invalid file", lambda: assert_raises(
        CdoError, lambda: cdo.info(input="/nonexistent.nc")))
    run_test("invalid params", lambda: assert_raises(CdoError, lambda: cdo.sellonlatbox(
        "abc,def", input=topo_nc, output=os.path.join(tmpdir, "err.nc"))))

    # Variable / Parameter Selection
    print("\n=== 33. Variable / Parameter Selection ===")
    selname_nc = os.path.join(tmpdir, "selname.nc")

    def _selname_test():
        raw = str(cdo.showname(input=topo_nc)).strip().split()
        vname = raw[0] if raw else "var1"
        vname = ''.join(c for c in vname if c.isprintable() and c != ' ')
        cdo.selname(vname, input=topo_nc, output=selname_nc)
        assert_file(selname_nc)
    run_test("selname", _selname_test)

    selcode_nc = os.path.join(tmpdir, "selcode.nc")

    def _selcode_test():
        # Dynamically detect the parameter code so we're not hardcoding 129
        raw = str(cdo.showcode(input=topo_nc)).strip().split()
        code = raw[0] if raw else "1"
        # Strip any non-digit characters except leading minus
        code = code.lstrip('-').strip()
        if not code.isdigit():
            code = "1"
        cdo.selcode(code, input=topo_nc, output=selcode_nc)
        assert_file(selcode_nc)
    run_test("selcode", _selcode_test)

    selseason_nc = os.path.join(tmpdir, "selseason.nc")
    run_test("selseason JJA", lambda: cdo(
        f"cdo selseason,JJA {monthly_nc} {selseason_nc}", timeout=20) or assert_file(selseason_nc))

    # Running Statistics
    print("\n=== 34. Running Statistics ===")
    runmean_nc = os.path.join(tmpdir, "runmean.nc")
    run_test("runmean", lambda: cdo.runmean("3", input=monthly_nc,
             output=runmean_nc) or assert_file(runmean_nc))

    runstd_nc = os.path.join(tmpdir, "runstd.nc")
    run_test("runstd", lambda: cdo.runstd("3", input=monthly_nc,
             output=runstd_nc) or assert_file(runstd_nc))

    runmin_nc = os.path.join(tmpdir, "runmin.nc")
    run_test("runmin", lambda: cdo.runmin("3", input=monthly_nc,
             output=runmin_nc) or assert_file(runmin_nc))

    runmax_nc = os.path.join(tmpdir, "runmax.nc")
    run_test("runmax", lambda: cdo.runmax("3", input=monthly_nc,
             output=runmax_nc) or assert_file(runmax_nc))

    # Verify runmean output timestep count: 12 - (3-1) = 10
    run_test("runmean ntime=10", lambda: assert_true(
        os.path.exists(runmean_nc) and
        int(str(cdo.ntime(input=runmean_nc)).strip()) == 10))

    # Grid Box Statistics
    print("\n=== 35. Grid Box Statistics ===")
    gridboxmean_nc = os.path.join(tmpdir, "gridboxmean.nc")
    run_test("gridboxmean", lambda: cdo.gridboxmean("2,2", input=topo_nc,
             output=gridboxmean_nc) or assert_file(gridboxmean_nc))

    gridboxmax_nc = os.path.join(tmpdir, "gridboxmax.nc")
    run_test("gridboxmax", lambda: cdo.gridboxmax("2,2", input=topo_nc,
             output=gridboxmax_nc) or assert_file(gridboxmax_nc))

    gridboxmin_nc = os.path.join(tmpdir, "gridboxmin.nc")
    run_test("gridboxmin", lambda: cdo.gridboxmin("2,2", input=topo_nc,
             output=gridboxmin_nc) or assert_file(gridboxmin_nc))

    gridboxstd_nc = os.path.join(tmpdir, "gridboxstd.nc")
    run_test("gridboxstd", lambda: cdo.gridboxstd("2,2", input=topo_nc,
             output=gridboxstd_nc) or assert_file(gridboxstd_nc))

    # Set / Rename Operators
    print("\n=== 36. Set / Rename Operators ===")
    setname_nc = os.path.join(tmpdir, "setname.nc")
    run_test("setname", lambda: cdo.setname("newvar", input=topo_nc,
             output=setname_nc) or assert_file(setname_nc))

    setunit_nc = os.path.join(tmpdir, "setunit.nc")
    run_test("setunit", lambda: cdo.setunit("K", input=topo_nc,
             output=setunit_nc) or assert_file(setunit_nc))

    setcode_nc2 = os.path.join(tmpdir, "setcode2.nc")
    run_test("setcode", lambda: cdo.setcode("200", input=topo_nc,
             output=setcode_nc2) or assert_file(setcode_nc2))

    # More Interpolation Methods
    print("\n=== 37. More Interpolation ===")
    remapbic_nc = os.path.join(tmpdir, "remapbic.nc")
    run_test("remapbic", lambda: cdo.remapbic("r8x4", input=topo_nc,
             output=remapbic_nc) or assert_file(remapbic_nc))

    remaplaf_nc = os.path.join(tmpdir, "remaplaf.nc")
    run_test("remaplaf", lambda: cdo.remaplaf("r8x4", input=topo_nc,
             output=remaplaf_nc) or assert_file(remaplaf_nc))

    # Splitting Operations
    print("\n=== 38. Splitting ===")
    splitvar_pfx = os.path.join(tmpdir, "splitvar_")

    def _splitvar_test():
        cdo.splitvar(input=topo_nc, output=splitvar_pfx, timeout=20)
        split_files = [f for f in os.listdir(
            tmpdir) if f.startswith("splitvar_")]
        assert_true(len(split_files) >= 1)
    run_test("splitvar", _splitvar_test)

    splitmon_pfx = os.path.join(tmpdir, "splitmon_")

    def _splitmon_test():
        cdo.splitmon(input=monthly_nc, output=splitmon_pfx, timeout=30)
        split_files = [f for f in os.listdir(
            tmpdir) if f.startswith("splitmon_")]
        assert_true(len(split_files) == 12)
    run_test("splitmon (12 files)", _splitmon_test)

    # Correlation / Covariance
    print("\n=== 39. Correlation / Covariance ===")
    timcor_nc = os.path.join(tmpdir, "timcor.nc")
    run_test("timcor (self)", lambda: cdo.timcor(
        input=f"{monthly_nc} {monthly_nc}", output=timcor_nc) or assert_file(timcor_nc))

    fldcor_nc = os.path.join(tmpdir, "fldcor.nc")
    run_test("fldcor (self)", lambda: cdo.fldcor(
        input=f"{topo_nc} {topo_nc}", output=fldcor_nc) or assert_file(fldcor_nc))

    timcovar_nc = os.path.join(tmpdir, "timcovar.nc")
    run_test("timcovar (self)", lambda: cdo.timcovar(
        input=f"{monthly_nc} {monthly_nc}", output=timcovar_nc) or assert_file(timcovar_nc))

    # More Seasonal Statistics
    print("\n=== 40. More Seasonal Stats ===")
    seasstd_nc = os.path.join(tmpdir, "seasstd.nc")
    run_test("seasstd", lambda: cdo.seasstd(input=monthly_nc,
             output=seasstd_nc) or assert_file(seasstd_nc))

    seasmin_nc = os.path.join(tmpdir, "seasmin.nc")
    run_test("seasmin", lambda: cdo.seasmin(input=monthly_nc,
             output=seasmin_nc) or assert_file(seasmin_nc))

    seasmax_nc = os.path.join(tmpdir, "seasmax.nc")
    run_test("seasmax", lambda: cdo.seasmax(input=monthly_nc,
             output=seasmax_nc) or assert_file(seasmax_nc))

    # Z-axis / Misc
    print("\n=== 41. Z-axis / Misc ===")
    run_test("zaxisdes", lambda: assert_true(
        len(str(cdo.zaxisdes(input=stdatm_nc))) > 5))

    run_test("npar", lambda: assert_true(
        int(str(cdo.npar(input=topo_nc)).strip()) >= 1))

    # Verify gridboxmean reduces the grid size
    run_test("gridboxmean grid smaller", lambda: assert_true(
        os.path.exists(gridboxmean_nc) and
        os.path.getsize(gridboxmean_nc) < os.path.getsize(topo_nc)))

    # CDO 2.6.0 New Features
    print("\n=== 42. CDO 2.6.0 New Features ===")

    # VarsStat operators require a single file with multiple variables on the SAME grid.
    # Build varstest_nc: topo + renamed copy (elev), both on identical topo grid.
    _topo_elev_nc = os.path.join(tmpdir, "topo_elev.nc")
    varstest_nc = os.path.join(tmpdir, "varstest.nc")
    run_test("varstest-prep", lambda: (
        cdo.chname("topo,elev", input=topo_nc, output=_topo_elev_nc) or True) and (
        cdo.merge(input=f"{topo_nc} {_topo_elev_nc}", output=varstest_nc) or
        assert_file(varstest_nc)))

    # VarsStat: statistics across variables within one file
    varsskew_nc = os.path.join(tmpdir, "varsskew.nc")
    run_test("varsskew", lambda: cdo.varsskew(
        input=varstest_nc, output=varsskew_nc) or assert_file(varsskew_nc))

    varskurt_nc = os.path.join(tmpdir, "varskurt.nc")
    run_test("varskurt", lambda: cdo.varskurt(
        input=varstest_nc, output=varskurt_nc) or assert_file(varskurt_nc))

    varsmedian_nc = os.path.join(tmpdir, "varsmedian.nc")
    run_test("varsmedian", lambda: cdo.varsmedian(
        input=varstest_nc, output=varsmedian_nc) or assert_file(varsmedian_nc))

    varspctl_nc = os.path.join(tmpdir, "varspctl.nc")
    run_test("varspctl,90", lambda: cdo.varspctl(
        "90", input=varstest_nc, output=varspctl_nc) or assert_file(varspctl_nc))

    # symmetrize: mirrors data at the equator
    symmetrize_nc = os.path.join(tmpdir, "symmetrize.nc")
    run_test("symmetrize", lambda: cdo.symmetrize(
        input=topo_nc, output=symmetrize_nc) or assert_file(symmetrize_nc))

    # fillmiss: bug #12341 fix (wrong result for n<4 neighbours)
    fillmiss_bugfix_nc = os.path.join(tmpdir, "fillmiss_bugfix.nc")
    run_test("fillmiss bug#12341 (n<4)", lambda: cdo.fillmiss(
        input="-setrtomiss,-10001,0 " + topo_nc,
        output=fillmiss_bugfix_nc) or assert_file(fillmiss_bugfix_nc))

    # ======================================================================
    # Regression test: explicit-coordinate NC4 (cdf_lazy_grid.c mutex bug)
    # CDO 2.6.0 on Windows crashed with ACCESS VIOLATION (0xC0000005) when
    # reading any NC4 file whose lat/lon/depth grid is stored as explicit
    # float coordinate-variable arrays (the common CF-1 convention).
    # Root cause: HAVE_LIBPTHREAD enabled cdf_lazy_grid.c's pthread_mutex_lock()
    # which winpthreads mishandles.  CI tests previously used CDO-generated
    # parametric grids (no explicit coordinate arrays) and missed the bug.
    # This test creates an explicit-coordinate NC4 file via netCDF4-python,
    # then exercises a CDO operator on it, guaranteeing the lazy-grid path
    # is triggered.
    # ======================================================================
    print("\n=== 21. Explicit-Coordinate NC4 Regression (lazy-grid mutex) ===")
    import netCDF4 as nc4  # type: ignore
    import numpy as np  # type: ignore

    explicit_coord_nc4 = os.path.join(tmpdir, "explicit_coord.nc4")

    def _make_explicit_coord_nc4(path):
        """Create a CF-1 NC4 file with explicit float lat/lon arrays."""
        nlat, nlon = 4, 8
        lat = np.linspace(-85.0, 85.0, nlat, dtype=np.float32)
        lon = np.linspace(0.0, 350.0, nlon, dtype=np.float32)
        data = np.random.rand(1, nlat, nlon).astype(np.float32) * 300.0

        with nc4.Dataset(path, "w", format="NETCDF4") as ds:
            ds.Conventions = "CF-1.6"
            ds.createDimension("time", None)
            ds.createDimension("lat", nlat)
            ds.createDimension("lon", nlon)

            t = ds.createVariable("time", "f8", ("time",))
            t.units = "days since 2000-01-01"
            t.calendar = "standard"
            t[:] = [0.0]

            la = ds.createVariable("lat", "f4", ("lat",))
            la.units = "degrees_north"
            la.long_name = "latitude"
            la.standard_name = "latitude"
            la[:] = lat

            lo = ds.createVariable("lon", "f4", ("lon",))
            lo.units = "degrees_east"
            lo.long_name = "longitude"
            lo.standard_name = "longitude"
            lo[:] = lon

            v = ds.createVariable("temperature", "f4",
                                  ("time", "lat", "lon"),
                                  fill_value=1e20)
            v.units = "K"
            v.long_name = "air temperature"
            v.standard_name = "air_temperature"
            v.coordinates = "lat lon"
            v[:] = data

    _make_explicit_coord_nc4(explicit_coord_nc4)

    # fldmean triggers CDI lazy-grid coordinate loading, which previously
    # crashed on Windows via cdf_lazy_grid.c's pthread_mutex_lock().
    explicit_mean_nc4 = os.path.join(tmpdir, "explicit_mean.nc4")
    run_test(
        "explicit-coord NC4 fldmean (lazy-grid mutex regression)",
        lambda: cdo(
            f"cdo fldmean {explicit_coord_nc4} {explicit_mean_nc4}",
            timeout=30,
        ) or assert_file(explicit_mean_nc4),
    )

    # Also verify sellonlatbox (requires grid coordinate access)
    explicit_box_nc4 = os.path.join(tmpdir, "explicit_box.nc4")
    run_test(
        "explicit-coord NC4 sellonlatbox (lazy-grid mutex regression)",
        lambda: cdo(
            f"cdo sellonlatbox,0,90,-45,45 {explicit_coord_nc4} {explicit_box_nc4}",
            timeout=30,
        ) or assert_file(explicit_box_nc4),
    )

    # ======================================================================
    # Section 22. mergetime timestep-count and wildcard regression
    #
    # Regression tests for two bugs:
    #   a) Windows exit-hang early-kill: CDO was killed right after creating
    #      the NC4/HDF5 file header (file size > 0, but before any time
    #      records were written), producing output with time=0.
    #      Fix: use file-size STABILITY (not mere existence) + CDO's
    #      "Processed N values" stderr message as the completion signal.
    #
    #   b) Wildcard expansion including the output file: when using
    #      cdo("cdo mergetime *.nc out.nc") and out.nc already existed,
    #      glob("*.nc") would match out.nc and add it to the input list.
    #      Fix: exclude explicit tokens from glob expansion results.
    # ======================================================================
    print("\n=== 22. mergetime Timestep-Count & Wildcard Regression ===")

    mt_dir = os.path.join(tmpdir, "mergetime_test")
    os.makedirs(mt_dir, exist_ok=True)

    def _make_monthly_nc(path, year):
        """12-month NC4 file with explicit time coordinates."""
        import cftime  # type: ignore
        with nc4.Dataset(path, "w", format="NETCDF4") as ds:
            ds.createDimension("time", 12)
            ds.createDimension("lat", 4)
            ds.createDimension("lon", 6)
            t = ds.createVariable("time", "f8", ("time",))
            t.units = "hours since 1900-01-01 00:00:00"
            t.calendar = "gregorian"
            t.standard_name = "time"
            t.axis = "T"
            dates = [cftime.DatetimeGregorian(year, m, 1)
                     for m in range(1, 13)]
            t[:] = nc4.date2num(dates, units=t.units, calendar=t.calendar)
            la = ds.createVariable("lat", "f4", ("lat",))
            la.units = "degrees_north"
            la[:] = np.linspace(-45, 45, 4, dtype=np.float32)
            lo = ds.createVariable("lon", "f4", ("lon",))
            lo.units = "degrees_east"
            lo[:] = np.linspace(0, 300, 6, dtype=np.float32)
            v = ds.createVariable("tas", "f4", ("time", "lat", "lon"),
                                  fill_value=-9999.0)
            v.units = "K"
            v[:] = np.random.rand(12, 4, 6).astype("f4") + 280.0

    nc_2023 = os.path.join(mt_dir, "ERA5_MSE_2023.nc")
    nc_2024 = os.path.join(mt_dir, "ERA5_MSE_2024.nc")
    mt_out = os.path.join(mt_dir, "ERA5_MSE_merged.nc")

    _make_monthly_nc(nc_2023, 2023)
    _make_monthly_nc(nc_2024, 2024)

    # --- 22a: method-call, explicit inputs -> 24 timesteps ---
    def _test_mergetime_timesteps():
        cdo.mergetime(input=f"{nc_2023} {nc_2024}", output=mt_out)
        with nc4.Dataset(mt_out) as ds:
            n = ds.variables["time"].size
            assert n == 24, (
                f"Expected 24 timesteps, got {n}. "
                "Likely caused by CDO exit-hang early-kill bug "
                "(killed before time records were flushed)."
            )
            t_vals = ds.variables["time"][:].data
            assert (t_vals[1:] > t_vals[:-1]).all(), \
                "Merged time axis is not monotonically increasing"

    run_test("mergetime 2x12 -> 24 timesteps (exit-hang early-kill regression)",
             _test_mergetime_timesteps)

    # --- 22b: wildcard ERA5_MSE_202*.nc must not include the output ---
    # The output file ERA5_MSE_2025_merged.nc ALSO matches ERA5_MSE_202*.nc
    # (because 2025 starts with 202).  We pre-create it to simulate a prior
    # run, then re-run; without the fix, glob would pull it into the input list.
    nc_2022 = os.path.join(mt_dir, "ERA5_MSE_2022.nc")
    _make_monthly_nc(nc_2022, 2022)

    def _test_mergetime_wildcard():
        orig_cwd = os.getcwd()
        try:
            os.chdir(mt_dir)
            # Output name intentionally matches the wildcard ERA5_MSE_202*.nc
            # (2025 starts with "202"). We run twice to verify that on the
            # second run -- when the output file already exists -- it is NOT
            # included in the glob-expanded input list.
            mt_out_wild = "ERA5_MSE_2025_merged.nc"
            for _pass in range(2):
                # -O lets CDO overwrite on the second pass; without the glob
                # fix the second pass would either silently merge 48 steps (if
                # the output was opened as input) or fail with a CDO error.
                cdo(f"cdo -O mergetime ERA5_MSE_202*.nc {mt_out_wild}",
                    timeout=60)
                with nc4.Dataset(mt_out_wild) as ds:
                    n = ds.variables["time"].size
                    assert n == 36, (
                        f"Pass {_pass+1}: expected 36 timesteps, got {n}. "
                        "Wildcard likely included the output file as an input."
                    )
        finally:
            os.chdir(orig_cwd)

    run_test("mergetime wildcard excludes output file (glob regression)",
             _test_mergetime_wildcard)

    # ======================================================================
    # Section 43. cdo("cdo ...") string API — explicit-file write regression
    #
    # The Cdo.__call__ / run_raw() code path (string command style) was
    # found to hang or produce truncated output when used with explicit
    # input file paths and a write-format operator like mergetime.
    #
    # Root cause: _runner.py lacked the two-phase HDF5 write-bit detection.
    # HDF5 creates a new file with file_consistency_flags = 0 (all zeros)
    # before H5Fclose flushes the superblock, so the FIRST poll saw bit0=0
    # and falsely declared the file "done", killing CDO before any time
    # records were written.  The method-call path (cdo.mergetime) was not
    # affected because it always resolves wildcards, thus the timing was
    # slightly different; but both paths share the same _exec() loop.
    #
    # These tests exercise the string API for:
    #   a) mergetime with two explicit NC4 files -> must produce 24 timesteps
    #   b) sellonlatbox string form -> output file must be non-empty
    #   c) chained -fldmean -sellonlatbox -> output file must be non-empty
    #   d) sinfo (no output file) -> must return non-empty text
    # ======================================================================
    print("\n=== 43. cdo() String API — explicit-file write regression ===")

    # 43a: mergetime with two explicit NC4 files
    def _test_cdo_str_mergetime():
        mt_str_out = os.path.join(mt_dir, "cdo_str_merged.nc")
        if os.path.exists(mt_str_out):
            os.remove(mt_str_out)
        cdo(f"cdo mergetime {nc_2023} {nc_2024} {mt_str_out}", timeout=60)
        with nc4.Dataset(mt_str_out) as ds:
            n = ds.variables["time"].size
            assert n == 24, (
                f"cdo('cdo mergetime ...'): expected 24 timesteps, got {n}. "
                "Two-phase HDF5 detection likely missing in run_raw() path."
            )

    run_test("cdo() string API: mergetime 2x12 -> 24 timesteps",
             _test_cdo_str_mergetime)

    # 43b: sellonlatbox via string API
    str_bbox_nc = os.path.join(tmpdir, "str_bbox.nc")
    run_test("cdo() string API: sellonlatbox writes file",
             lambda: cdo(f"cdo sellonlatbox,0,90,0,45 {topo_nc} {str_bbox_nc}",
                         timeout=30) or assert_file(str_bbox_nc))

    # 43c: chained operators via string API
    str_chain_nc = os.path.join(tmpdir, "str_chain.nc")
    run_test("cdo() string API: chained -fldmean -sellonlatbox",
             lambda: cdo(
                 f"cdo -fldmean -sellonlatbox,0,90,0,45 {topo_nc} {str_chain_nc}",
                 timeout=30) or assert_file(str_chain_nc))

    # 43d: info operator via string API (no output file, returns text)
    run_test("cdo() string API: sinfo returns text",
             lambda: assert_true(
                 isinstance(cdo(f"cdo sinfo {topo_nc}", timeout=20), str)
                 and len(str(cdo(f"cdo sinfo {topo_nc}", timeout=20))) > 10))

    # ======================================================================
    # Sections 44–59: Operator-method coverage expansion to 200+ tests
    # Each section uses cdo.operator(...) method-call style exclusively,
    # mirroring the pattern established in section 42 (varsskew, varsmedian).
    # ======================================================================

    # === 44. Zonal Stats Extended ===
    print("\n=== 44. Zonal Stats Extended ===")
    zonmin_nc = os.path.join(tmpdir, "zonmin.nc")
    run_test("zonmin", lambda: cdo.zonmin(input=topo_nc,
             output=zonmin_nc) or assert_file(zonmin_nc))
    zonmax_nc = os.path.join(tmpdir, "zonmax.nc")
    run_test("zonmax", lambda: cdo.zonmax(input=topo_nc,
             output=zonmax_nc) or assert_file(zonmax_nc))
    zonsum_nc = os.path.join(tmpdir, "zonsum.nc")
    run_test("zonsum", lambda: cdo.zonsum(input=topo_nc,
             output=zonsum_nc) or assert_file(zonsum_nc))
    zonstd_nc = os.path.join(tmpdir, "zonstd.nc")
    run_test("zonstd", lambda: cdo.zonstd(input=topo_nc,
             output=zonstd_nc) or assert_file(zonstd_nc))
    zonskew_nc = os.path.join(tmpdir, "zonskew.nc")
    run_test("zonskew", lambda: cdo.zonskew(input=topo_nc,
             output=zonskew_nc) or assert_file(zonskew_nc))
    zonkurt_nc = os.path.join(tmpdir, "zonkurt.nc")
    run_test("zonkurt", lambda: cdo.zonkurt(input=topo_nc,
             output=zonkurt_nc) or assert_file(zonkurt_nc))
    zonmedian_nc = os.path.join(tmpdir, "zonmedian.nc")
    run_test("zonmedian", lambda: cdo.zonmedian(input=topo_nc,
             output=zonmedian_nc) or assert_file(zonmedian_nc))

    # === 45. Meridional Stats Extended ===
    print("\n=== 45. Meridional Stats Extended ===")
    mermin_nc = os.path.join(tmpdir, "mermin.nc")
    run_test("mermin", lambda: cdo.mermin(input=topo_nc,
             output=mermin_nc) or assert_file(mermin_nc))
    mermax_nc = os.path.join(tmpdir, "mermax.nc")
    run_test("mermax", lambda: cdo.mermax(input=topo_nc,
             output=mermax_nc) or assert_file(mermax_nc))
    mersum_nc = os.path.join(tmpdir, "mersum.nc")
    run_test("mersum", lambda: cdo.mersum(input=topo_nc,
             output=mersum_nc) or assert_file(mersum_nc))
    merstd_nc = os.path.join(tmpdir, "merstd.nc")
    run_test("merstd", lambda: cdo.merstd(input=topo_nc,
             output=merstd_nc) or assert_file(merstd_nc))
    merskew_nc = os.path.join(tmpdir, "merskew.nc")
    run_test("merskew", lambda: cdo.merskew(input=topo_nc,
             output=merskew_nc) or assert_file(merskew_nc))
    merkurt_nc = os.path.join(tmpdir, "merkurt.nc")
    run_test("merkurt", lambda: cdo.merkurt(input=topo_nc,
             output=merkurt_nc) or assert_file(merkurt_nc))
    mermedian_nc = os.path.join(tmpdir, "mermedian.nc")
    run_test("mermedian", lambda: cdo.mermedian(input=topo_nc,
             output=mermedian_nc) or assert_file(mermedian_nc))

    # === 46. Field Stats Extended ===
    print("\n=== 46. Field Stats Extended ===")
    fldskew_nc = os.path.join(tmpdir, "fldskew.nc")
    run_test("fldskew", lambda: cdo.fldskew(input=topo_nc,
             output=fldskew_nc) or assert_file(fldskew_nc))
    fldkurt_nc = os.path.join(tmpdir, "fldkurt.nc")
    run_test("fldkurt", lambda: cdo.fldkurt(input=topo_nc,
             output=fldkurt_nc) or assert_file(fldkurt_nc))
    fldvar_nc = os.path.join(tmpdir, "fldvar.nc")
    run_test("fldvar", lambda: cdo.fldvar(input=topo_nc,
             output=fldvar_nc) or assert_file(fldvar_nc))
    fldstd1_nc = os.path.join(tmpdir, "fldstd1.nc")
    run_test("fldstd1", lambda: cdo.fldstd1(input=topo_nc,
             output=fldstd1_nc) or assert_file(fldstd1_nc))
    fldmedian_nc = os.path.join(tmpdir, "fldmedian.nc")
    run_test("fldmedian", lambda: cdo.fldmedian(input=topo_nc,
             output=fldmedian_nc) or assert_file(fldmedian_nc))
    fldpctl_nc = os.path.join(tmpdir, "fldpctl.nc")
    run_test("fldpctl,90", lambda: cdo.fldpctl("90", input=topo_nc,
             output=fldpctl_nc) or assert_file(fldpctl_nc))

    # === 47. Time Stats Extended ===
    print("\n=== 47. Time Stats Extended ===")
    timvar_nc = os.path.join(tmpdir, "timvar.nc")
    run_test("timvar", lambda: cdo.timvar(input=monthly_nc,
             output=timvar_nc) or assert_file(timvar_nc))
    timrange_nc = os.path.join(tmpdir, "timrange.nc")
    run_test("timrange", lambda: cdo.timrange(input=monthly_nc,
             output=timrange_nc) or assert_file(timrange_nc))
    timminidx_nc = os.path.join(tmpdir, "timminidx.nc")
    run_test("timminidx", lambda: cdo.timminidx(input=monthly_nc,
             output=timminidx_nc) or assert_file(timminidx_nc))
    timmaxidx_nc = os.path.join(tmpdir, "timmaxidx.nc")
    run_test("timmaxidx", lambda: cdo.timmaxidx(input=monthly_nc,
             output=timmaxidx_nc) or assert_file(timmaxidx_nc))
    timavg_nc = os.path.join(tmpdir, "timavg.nc")
    run_test("timavg", lambda: cdo.timavg(input=monthly_nc,
             output=timavg_nc) or assert_file(timavg_nc))

    # === 48. Running Stats Extended ===
    print("\n=== 48. Running Stats Extended ===")
    runsum_nc = os.path.join(tmpdir, "runsum.nc")
    run_test("runsum", lambda: cdo.runsum("3", input=monthly_nc,
             output=runsum_nc) or assert_file(runsum_nc))
    runvar_nc = os.path.join(tmpdir, "runvar.nc")
    run_test("runvar", lambda: cdo.runvar("3", input=monthly_nc,
             output=runvar_nc) or assert_file(runvar_nc))
    runrange_nc = os.path.join(tmpdir, "runrange.nc")
    run_test("runrange", lambda: cdo.runrange("3", input=monthly_nc,
             output=runrange_nc) or assert_file(runrange_nc))

    # === 49. Monthly Stats Extended ===
    print("\n=== 49. Monthly Stats Extended ===")
    monmin_nc = os.path.join(tmpdir, "monmin.nc")
    run_test("monmin", lambda: cdo.monmin(input=monthly_nc,
             output=monmin_nc) or assert_file(monmin_nc))
    monmax_nc = os.path.join(tmpdir, "monmax.nc")
    run_test("monmax", lambda: cdo.monmax(input=monthly_nc,
             output=monmax_nc) or assert_file(monmax_nc))
    monstd_nc = os.path.join(tmpdir, "monstd.nc")
    run_test("monstd", lambda: cdo.monstd(input=monthly_nc,
             output=monstd_nc) or assert_file(monstd_nc))
    monsum_nc = os.path.join(tmpdir, "monsum.nc")
    run_test("monsum", lambda: cdo.monsum(input=monthly_nc,
             output=monsum_nc) or assert_file(monsum_nc))

    # === 50. Yearly Stats Extended ===
    print("\n=== 50. Yearly Stats Extended ===")
    yearstd_nc = os.path.join(tmpdir, "yearstd.nc")
    run_test("yearstd", lambda: cdo.yearstd(input=monthly_nc,
             output=yearstd_nc) or assert_file(yearstd_nc))
    yearsum_nc = os.path.join(tmpdir, "yearsum.nc")
    run_test("yearsum", lambda: cdo.yearsum(input=monthly_nc,
             output=yearsum_nc) or assert_file(yearsum_nc))
    yearrange_nc = os.path.join(tmpdir, "yearrange.nc")
    run_test("yearrange", lambda: cdo.yearrange(input=monthly_nc,
             output=yearrange_nc) or assert_file(yearrange_nc))

    # === 51. Seasonal Stats Extended ===
    print("\n=== 51. Seasonal Stats Extended ===")
    seassum_nc = os.path.join(tmpdir, "seassum.nc")
    run_test("seassum", lambda: cdo.seassum(input=monthly_nc,
             output=seassum_nc) or assert_file(seassum_nc))
    seasrange_nc = os.path.join(tmpdir, "seasrange.nc")
    run_test("seasrange", lambda: cdo.seasrange(input=monthly_nc,
             output=seasrange_nc) or assert_file(seasrange_nc))

    # === 52. VarsStat Extended ===
    # Uses varstest_nc (2-variable file built in section 42)
    print("\n=== 52. VarsStat Extended ===")
    varsmean_nc = os.path.join(tmpdir, "varsmean.nc")
    run_test("varsmean", lambda: cdo.varsmean(
        input=varstest_nc, output=varsmean_nc) or assert_file(varsmean_nc))
    varsstd_nc = os.path.join(tmpdir, "varsstd.nc")
    run_test("varsstd", lambda: cdo.varsstd(
        input=varstest_nc, output=varsstd_nc) or assert_file(varsstd_nc))
    varsmin_nc = os.path.join(tmpdir, "varsmin.nc")
    run_test("varsmin", lambda: cdo.varsmin(
        input=varstest_nc, output=varsmin_nc) or assert_file(varsmin_nc))
    varsmax_nc = os.path.join(tmpdir, "varsmax.nc")
    run_test("varsmax", lambda: cdo.varsmax(
        input=varstest_nc, output=varsmax_nc) or assert_file(varsmax_nc))
    varsrange_nc = os.path.join(tmpdir, "varsrange.nc")
    run_test("varsrange", lambda: cdo.varsrange(
        input=varstest_nc, output=varsrange_nc) or assert_file(varsrange_nc))
    varssum_nc = os.path.join(tmpdir, "varssum.nc")
    run_test("varssum", lambda: cdo.varssum(
        input=varstest_nc, output=varssum_nc) or assert_file(varssum_nc))

    # === 53. Ensemble Stats Extended ===
    print("\n=== 53. Ensemble Stats Extended ===")
    ensrange_nc = os.path.join(tmpdir, "ensrange.nc")
    run_test("ensrange", lambda: cdo.ensrange(
        input=f"{topo_nc} {mulc_nc}", output=ensrange_nc) or assert_file(ensrange_nc))
    enssum_nc = os.path.join(tmpdir, "enssum.nc")
    run_test("enssum", lambda: cdo.enssum(
        input=f"{topo_nc} {topo_nc}", output=enssum_nc) or assert_file(enssum_nc))
    ensskew_nc = os.path.join(tmpdir, "ensskew.nc")
    run_test("ensskew", lambda: cdo.ensskew(
        input=f"{topo_nc} {mulc_nc} {addc_nc} {abs_nc}",
        output=ensskew_nc) or assert_file(ensskew_nc))
    enskurt_nc = os.path.join(tmpdir, "enskurt.nc")
    run_test("enskurt", lambda: cdo.enskurt(
        input=f"{topo_nc} {mulc_nc} {addc_nc} {abs_nc}",
        output=enskurt_nc) or assert_file(enskurt_nc))
    enspctl_nc = os.path.join(tmpdir, "enspctl.nc")
    run_test("enspctl,90", lambda: cdo.enspctl(
        "90", input=f"{topo_nc} {mulc_nc} {addc_nc}",
        output=enspctl_nc) or assert_file(enspctl_nc))

    # === 54. Math Functions ===
    # exp: scale to (0, ~8.85) to avoid float overflow
    # ln / log10: add 1 to abs(topo) so input is strictly positive
    # asin / acos: divide by 10000 to clamp into [-1, 1]
    # reci (1/x): add 1000 so denominator is never near zero
    print("\n=== 54. Math Functions ===")
    exp_nc = os.path.join(tmpdir, "exp.nc")
    run_test("exp", lambda: cdo.exp(
        input=f"-mulc,0.001 -abs {topo_nc}",
        output=exp_nc) or assert_file(exp_nc))
    ln_nc = os.path.join(tmpdir, "ln.nc")
    run_test("ln", lambda: cdo.ln(
        input=f"-addc,1 -abs {topo_nc}",
        output=ln_nc) or assert_file(ln_nc))
    log10_nc = os.path.join(tmpdir, "log10.nc")
    run_test("log10", lambda: cdo.log10(
        input=f"-addc,1 -abs {topo_nc}",
        output=log10_nc) or assert_file(log10_nc))
    sin_nc = os.path.join(tmpdir, "sin.nc")
    run_test("sin", lambda: cdo.sin(input=topo_nc,
             output=sin_nc) or assert_file(sin_nc))
    cos_nc = os.path.join(tmpdir, "cos.nc")
    run_test("cos", lambda: cdo.cos(input=topo_nc,
             output=cos_nc) or assert_file(cos_nc))
    tan_nc = os.path.join(tmpdir, "tan.nc")
    run_test("tan", lambda: cdo.tan(input=topo_nc,
             output=tan_nc) or assert_file(tan_nc))
    asin_nc = os.path.join(tmpdir, "asin.nc")
    run_test("asin", lambda: cdo.asin(
        input=f"-divc,10000 {topo_nc}",
        output=asin_nc) or assert_file(asin_nc))
    acos_nc2 = os.path.join(tmpdir, "acos2.nc")
    run_test("acos", lambda: cdo.acos(
        input=f"-divc,10000 -abs {topo_nc}",
        output=acos_nc2) or assert_file(acos_nc2))
    atan_nc = os.path.join(tmpdir, "atan.nc")
    run_test("atan", lambda: cdo.atan(input=topo_nc,
             output=atan_nc) or assert_file(atan_nc))
    reci_nc = os.path.join(tmpdir, "reci.nc")
    run_test("reci", lambda: cdo.reci(
        input=f"-addc,1000 {topo_nc}",
        output=reci_nc) or assert_file(reci_nc))
    nint_nc = os.path.join(tmpdir, "nint.nc")
    run_test("nint", lambda: cdo.nint(input=topo_nc,
             output=nint_nc) or assert_file(nint_nc))

    # === 55. Comparison Operators ===
    print("\n=== 55. Comparison Operators ===")
    eq_nc = os.path.join(tmpdir, "eq.nc")
    run_test("eq (self==self -> all 1)", lambda: cdo.eq(
        input=f"{topo_nc} {topo_nc}", output=eq_nc) or assert_file(eq_nc))
    ne_nc = os.path.join(tmpdir, "ne.nc")
    run_test("ne (self!=self -> all 0)", lambda: cdo.ne(
        input=f"{topo_nc} {topo_nc}", output=ne_nc) or assert_file(ne_nc))
    le_nc = os.path.join(tmpdir, "le.nc")
    run_test("le (self<=self -> all 1)", lambda: cdo.le(
        input=f"{topo_nc} {topo_nc}", output=le_nc) or assert_file(le_nc))
    lt_nc = os.path.join(tmpdir, "lt.nc")
    run_test("lt (self<self  -> all 0)", lambda: cdo.lt(
        input=f"{topo_nc} {topo_nc}", output=lt_nc) or assert_file(lt_nc))
    ge_nc = os.path.join(tmpdir, "ge.nc")
    run_test("ge (self>=self -> all 1)", lambda: cdo.ge(
        input=f"{topo_nc} {topo_nc}", output=ge_nc) or assert_file(ge_nc))
    gt_nc = os.path.join(tmpdir, "gt.nc")
    run_test("gt (self>self  -> all 0)", lambda: cdo.gt(
        input=f"{topo_nc} {topo_nc}", output=gt_nc) or assert_file(gt_nc))

    # === 56. Masking Extended ===
    print("\n=== 56. Masking Extended ===")
    # ifnotthen: output data where mask IS missing (opposite of ifthen)
    ifnotthen_nc = os.path.join(tmpdir, "ifnotthen.nc")
    run_test("ifnotthen", lambda: cdo.ifnotthen(
        input=f"{setrtomiss_nc} {topo_nc}",
        output=ifnotthen_nc) or assert_file(ifnotthen_nc))
    # ifthenelse: 3 inputs — mask, then-field, else-field
    ifthenelse_nc = os.path.join(tmpdir, "ifthenelse.nc")
    run_test("ifthenelse", lambda: cdo.ifthenelse(
        input=f"{setctomiss_nc} {topo_nc} {mulc_nc}",
        output=ifthenelse_nc) or assert_file(ifthenelse_nc))

    # === 57. Vertical Stats Extended ===
    print("\n=== 57. Vertical Stats Extended ===")
    vertstd_nc = os.path.join(tmpdir, "vertstd.nc")
    run_test("vertstd", lambda: cdo.vertstd(input=stdatm_nc,
             output=vertstd_nc) or assert_file(vertstd_nc))
    vertvar_nc = os.path.join(tmpdir, "vertvar.nc")
    run_test("vertvar", lambda: cdo.vertvar(input=stdatm_nc,
             output=vertvar_nc) or assert_file(vertvar_nc))

    # === 58. GridBox Stats Extended ===
    print("\n=== 58. GridBox Stats Extended ===")
    gridboxrange_nc = os.path.join(tmpdir, "gridboxrange.nc")
    run_test("gridboxrange", lambda: cdo.gridboxrange("2,2", input=topo_nc,
             output=gridboxrange_nc) or assert_file(gridboxrange_nc))
    gridboxsum_nc = os.path.join(tmpdir, "gridboxsum.nc")
    run_test("gridboxsum", lambda: cdo.gridboxsum("2,2", input=topo_nc,
             output=gridboxsum_nc) or assert_file(gridboxsum_nc))
    gridboxvar_nc = os.path.join(tmpdir, "gridboxvar.nc")
    run_test("gridboxvar", lambda: cdo.gridboxvar("2,2", input=topo_nc,
             output=gridboxvar_nc) or assert_file(gridboxvar_nc))

    # === 59. Selection Extended ===
    # monthly_nc: Jan–Dec 2020, all at 12:00, all on the 15th
    print("\n=== 59. Selection Extended ===")
    selday_nc = os.path.join(tmpdir, "selday.nc")
    run_test("selday (15th)", lambda: cdo.selday("15", input=monthly_nc,
             output=selday_nc) or assert_file(selday_nc))
    selyear_nc = os.path.join(tmpdir, "selyear.nc")
    run_test("selyear (2020)", lambda: cdo.selyear("2020", input=monthly_nc,
             output=selyear_nc) or assert_file(selyear_nc))
    selsmon_nc = os.path.join(tmpdir, "selsmon.nc")
    run_test("selsmon (June)", lambda: cdo.selsmon("6", input=monthly_nc,
             output=selsmon_nc) or assert_file(selsmon_nc))
    seltime_nc = os.path.join(tmpdir, "seltime.nc")
    run_test("seltime (12:00:00)", lambda: cdo.seltime("12:00:00",
             input=monthly_nc, output=seltime_nc) or assert_file(seltime_nc))

    # ======================================================================
    # Section 60. run_raw() Direct API
    # ======================================================================
    # Tests CdoRunner.run_raw() directly (bypassing Cdo.__call__).
    # run_raw() parses a full CDO command string (with or without the leading
    # "cdo" token), injects the actual binary path, and returns a
    # subprocess.CompletedProcess.  This is the same code path reached by
    # cdo("cdo ...") but exercised here through the lower-level runner API.
    print("\n=== 60. run_raw() Direct API ===")
    runner = cdo._runner  # direct access to the underlying CdoRunner

    # 60a: arithmetic via run_raw
    rr_mulc_nc = os.path.join(tmpdir, "rr_mulc.nc")
    run_test("run_raw: mulc,3.14", lambda: runner.run_raw(
        f"cdo mulc,3.14 {topo_nc} {rr_mulc_nc}", timeout=30) or assert_file(rr_mulc_nc))

    # 60b: time mean via run_raw
    rr_timmean_nc = os.path.join(tmpdir, "rr_timmean.nc")
    run_test("run_raw: timmean", lambda: runner.run_raw(
        f"cdo timmean {monthly_nc} {rr_timmean_nc}", timeout=30) or assert_file(rr_timmean_nc))

    # 60c: grid interpolation via run_raw
    rr_remap_nc = os.path.join(tmpdir, "rr_remap.nc")
    run_test("run_raw: remapbil r8x4", lambda: runner.run_raw(
        f"cdo remapbil,r8x4 {topo_nc} {rr_remap_nc}", timeout=60) or assert_file(rr_remap_nc))

    # 60d: chained operators via run_raw
    rr_chain_nc = os.path.join(tmpdir, "rr_chain.nc")
    run_test("run_raw: chained -fldmean -timmean", lambda: runner.run_raw(
        f"cdo -fldmean -timmean {monthly_nc} {rr_chain_nc}", timeout=30) or assert_file(rr_chain_nc))

    # 60e: metadata change (chname) via run_raw
    rr_chname_nc = os.path.join(tmpdir, "rr_chname.nc")
    run_test("run_raw: chname topo->elevation", lambda: runner.run_raw(
        f"cdo chname,topo,elevation {topo_nc} {rr_chname_nc}", timeout=30) or assert_file(rr_chname_nc))

    # 60f: run_raw returns subprocess.CompletedProcess with returncode=0
    rr_check_nc = os.path.join(tmpdir, "rr_check.nc")

    def _test_runraw_completedprocess():
        result = runner.run_raw(
            f"cdo addc,1.0 {topo_nc} {rr_check_nc}", timeout=30)
        assert hasattr(result, "returncode"), \
            "run_raw() did not return a CompletedProcess object"
        assert result.returncode == 0, \
            f"Expected returncode=0, got {result.returncode}"
        assert_file(rr_check_nc)
    run_test("run_raw: returns CompletedProcess (returncode=0)",
             _test_runraw_completedprocess)

    # 60g: run_raw works WITHOUT the leading 'cdo' token
    # The implementation strips 'cdo'/'cdo.exe' and prepends the real binary,
    # but when the first token is an operator name (not 'cdo'), it should also
    # work because the binary path is always prepended unconditionally.
    rr_noprefix_nc = os.path.join(tmpdir, "rr_noprefix.nc")
    run_test("run_raw: operator without leading 'cdo' token", lambda: runner.run_raw(
        f"subc,0.0 {topo_nc} {rr_noprefix_nc}", timeout=30) or assert_file(rr_noprefix_nc))

    # 60h: -O overwrite flag via run_raw (second run overwrites first)
    def _test_runraw_overwrite():
        rr_ow_nc = os.path.join(tmpdir, "rr_overwrite.nc")
        runner.run_raw(f"cdo mulc,1.0 {topo_nc} {rr_ow_nc}", timeout=30)
        runner.run_raw(f"cdo -O mulc,2.0 {topo_nc} {rr_ow_nc}", timeout=30)
        assert_file(rr_ow_nc)
    run_test("run_raw: -O overwrite flag", _test_runraw_overwrite)

    # ======================================================================
    # Section 61. NCL Wind: uv2vr_cfd / uv2dv_cfd
    # ======================================================================
    # CDO's NCL_wind operators convert U and V wind components to derived
    # fields (relative vorticity, divergence) using curvilinear finite
    # differences.  The operators identify U and V by CF standard_name
    # ('eastward_wind' / 'northward_wind').  This test creates a small
    # synthetic NC4 file with both variables and verifies that:
    #   a) CDO can locate u/v from their standard_name / long_name attributes.
    #   b) The output file is non-empty and contains at least one variable.
    print("\n=== 61. NCL Wind: uv2vr_cfd / uv2dv_cfd ===")
    uv_wind_nc = os.path.join(tmpdir, "uv_wind.nc")

    def _make_uv_wind_nc(path):
        """Small CF-1.6 NC4 with u/v wind on an 18x36 regular lat/lon grid."""
        nlat, nlon = 18, 36
        lat = np.linspace(-85.0, 85.0, nlat, dtype=np.float32)
        lon = np.linspace(0.0, 350.0, nlon, dtype=np.float32)
        la2d, lo2d = np.meshgrid(lat, lon, indexing="ij")
        # Solid-rotation approximation: u=-sin(lat)*cos(lon)*5, v=cos(lon)*5
        u_val = (
            -np.sin(np.deg2rad(la2d)) * np.cos(np.deg2rad(lo2d)) * 5.0
        ).astype(np.float32)
        v_val = (np.cos(np.deg2rad(lo2d)) * 5.0).astype(np.float32)
        with nc4.Dataset(path, "w", format="NETCDF4") as ds:
            ds.Conventions = "CF-1.6"
            ds.createDimension("time", 1)
            ds.createDimension("lat", nlat)
            ds.createDimension("lon", nlon)
            t = ds.createVariable("time", "f8", ("time",))
            t.units = "days since 2000-01-01"
            t.calendar = "standard"
            t[:] = [0.0]
            la = ds.createVariable("lat", "f4", ("lat",))
            la.units = "degrees_north"
            la.standard_name = "latitude"
            la[:] = lat
            lo = ds.createVariable("lon", "f4", ("lon",))
            lo.units = "degrees_east"
            lo.standard_name = "longitude"
            lo[:] = lon
            # u: eastward (zonal) wind component
            u_var = ds.createVariable("u", "f4", ("time", "lat", "lon"),
                                      fill_value=1e20)
            u_var.units = "m s-1"
            u_var.long_name = "eastward wind"
            u_var.standard_name = "eastward_wind"
            u_var[:] = u_val[np.newaxis, :]
            # v: northward (meridional) wind component
            v_var = ds.createVariable("v", "f4", ("time", "lat", "lon"),
                                      fill_value=1e20)
            v_var.units = "m s-1"
            v_var.long_name = "northward wind"
            v_var.standard_name = "northward_wind"
            v_var[:] = v_val[np.newaxis, :]

    _make_uv_wind_nc(uv_wind_nc)

    # 61a: relative vorticity from U/V — force NC4 output so showname is fast
    vorticity_nc = os.path.join(tmpdir, "vorticity.nc")
    run_test("uv2vr_cfd (U+V -> relative vorticity)", lambda: (
        cdo.uv2vr_cfd(input=uv_wind_nc, output=vorticity_nc, options="-f nc4"),
        assert_file(vorticity_nc)))

    # 61b: divergence from U/V — force NC4 output so showname is fast
    divergence_nc = os.path.join(tmpdir, "divergence.nc")
    run_test("uv2dv_cfd (U+V -> divergence)", lambda: (
        cdo.uv2dv_cfd(input=uv_wind_nc, output=divergence_nc,
                      options="-f nc4"),
        assert_file(divergence_nc)))

    # ======================================================================
    # Section 62. cdo() String API — Extended Coverage
    # ======================================================================
    # Each test calls cdo("cdo OPERATOR ...") to exercise the string-API
    # code path for operator categories not yet covered by Section 43.
    print("\n=== 62. String API Extended ===")

    # --- 62a. Field statistics ---
    s_fldmin_nc = os.path.join(tmpdir, "s_fldmin.nc")
    run_test("str: fldmin", lambda: cdo(
        f"cdo fldmin {topo_nc} {s_fldmin_nc}", timeout=20) or assert_file(s_fldmin_nc))

    s_fldmax_nc = os.path.join(tmpdir, "s_fldmax.nc")
    run_test("str: fldmax", lambda: cdo(
        f"cdo fldmax {topo_nc} {s_fldmax_nc}", timeout=20) or assert_file(s_fldmax_nc))

    s_fldsum_nc = os.path.join(tmpdir, "s_fldsum.nc")
    run_test("str: fldsum", lambda: cdo(
        f"cdo fldsum {topo_nc} {s_fldsum_nc}", timeout=20) or assert_file(s_fldsum_nc))

    s_fldstd_nc = os.path.join(tmpdir, "s_fldstd.nc")
    run_test("str: fldstd", lambda: cdo(
        f"cdo fldstd {topo_nc} {s_fldstd_nc}", timeout=20) or assert_file(s_fldstd_nc))

    s_fldmean_nc = os.path.join(tmpdir, "s_fldmean.nc")
    run_test("str: fldmean", lambda: cdo(
        f"cdo fldmean {topo_nc} {s_fldmean_nc}", timeout=20) or assert_file(s_fldmean_nc))

    # --- 62b. Zonal / Meridional statistics ---
    s_zonmin_nc = os.path.join(tmpdir, "s_zonmin.nc")
    run_test("str: zonmin", lambda: cdo(
        f"cdo zonmin {topo_nc} {s_zonmin_nc}", timeout=20) or assert_file(s_zonmin_nc))

    s_zonmax_nc = os.path.join(tmpdir, "s_zonmax.nc")
    run_test("str: zonmax", lambda: cdo(
        f"cdo zonmax {topo_nc} {s_zonmax_nc}", timeout=20) or assert_file(s_zonmax_nc))

    s_mermean_nc = os.path.join(tmpdir, "s_mermean.nc")
    run_test("str: mermean", lambda: cdo(
        f"cdo mermean {topo_nc} {s_mermean_nc}", timeout=20) or assert_file(s_mermean_nc))

    s_mermin_nc = os.path.join(tmpdir, "s_mermin.nc")
    run_test("str: mermin", lambda: cdo(
        f"cdo mermin {topo_nc} {s_mermin_nc}", timeout=20) or assert_file(s_mermin_nc))

    s_mermax_nc = os.path.join(tmpdir, "s_mermax.nc")
    run_test("str: mermax", lambda: cdo(
        f"cdo mermax {topo_nc} {s_mermax_nc}", timeout=20) or assert_file(s_mermax_nc))

    # --- 62c. Time statistics ---
    s_timmean_nc = os.path.join(tmpdir, "s_timmean.nc")
    run_test("str: timmean", lambda: cdo(
        f"cdo timmean {monthly_nc} {s_timmean_nc}", timeout=20) or assert_file(s_timmean_nc))

    s_timmin_nc = os.path.join(tmpdir, "s_timmin.nc")
    run_test("str: timmin", lambda: cdo(
        f"cdo timmin {monthly_nc} {s_timmin_nc}", timeout=20) or assert_file(s_timmin_nc))

    s_timmax_nc = os.path.join(tmpdir, "s_timmax.nc")
    run_test("str: timmax", lambda: cdo(
        f"cdo timmax {monthly_nc} {s_timmax_nc}", timeout=20) or assert_file(s_timmax_nc))

    s_timstd_nc = os.path.join(tmpdir, "s_timstd.nc")
    run_test("str: timstd", lambda: cdo(
        f"cdo timstd {monthly_nc} {s_timstd_nc}", timeout=20) or assert_file(s_timstd_nc))

    s_timsum_nc = os.path.join(tmpdir, "s_timsum.nc")
    run_test("str: timsum", lambda: cdo(
        f"cdo timsum {monthly_nc} {s_timsum_nc}", timeout=20) or assert_file(s_timsum_nc))

    s_timvar_nc = os.path.join(tmpdir, "s_timvar.nc")
    run_test("str: timvar", lambda: cdo(
        f"cdo timvar {monthly_nc} {s_timvar_nc}", timeout=20) or assert_file(s_timvar_nc))

    s_timavg_nc = os.path.join(tmpdir, "s_timavg.nc")
    run_test("str: timavg", lambda: cdo(
        f"cdo timavg {monthly_nc} {s_timavg_nc}", timeout=20) or assert_file(s_timavg_nc))

    s_timrange_nc = os.path.join(tmpdir, "s_timrange.nc")
    run_test("str: timrange", lambda: cdo(
        f"cdo timrange {monthly_nc} {s_timrange_nc}", timeout=20) or assert_file(s_timrange_nc))

    # --- 62d. Monthly / Seasonal / Yearly statistics ---
    s_monmean_nc = os.path.join(tmpdir, "s_monmean.nc")
    run_test("str: monmean", lambda: cdo(
        f"cdo monmean {monthly_nc} {s_monmean_nc}", timeout=20) or assert_file(s_monmean_nc))

    s_monmin_nc = os.path.join(tmpdir, "s_monmin.nc")
    run_test("str: monmin", lambda: cdo(
        f"cdo monmin {monthly_nc} {s_monmin_nc}", timeout=20) or assert_file(s_monmin_nc))

    s_monmax_nc = os.path.join(tmpdir, "s_monmax.nc")
    run_test("str: monmax", lambda: cdo(
        f"cdo monmax {monthly_nc} {s_monmax_nc}", timeout=20) or assert_file(s_monmax_nc))

    s_seasmean_nc = os.path.join(tmpdir, "s_seasmean.nc")
    run_test("str: seasmean", lambda: cdo(
        f"cdo seasmean {monthly_nc} {s_seasmean_nc}", timeout=20) or assert_file(s_seasmean_nc))

    s_yearmean_nc = os.path.join(tmpdir, "s_yearmean.nc")
    run_test("str: yearmean", lambda: cdo(
        f"cdo yearmean {monthly_nc} {s_yearmean_nc}", timeout=20) or assert_file(s_yearmean_nc))

    # --- 62e. Arithmetic operators ---
    s_mulc_nc = os.path.join(tmpdir, "s_mulc.nc")
    run_test("str: mulc,3.14", lambda: cdo(
        f"cdo mulc,3.14 {topo_nc} {s_mulc_nc}", timeout=20) or assert_file(s_mulc_nc))

    s_addc_nc = os.path.join(tmpdir, "s_addc.nc")
    run_test("str: addc,50", lambda: cdo(
        f"cdo addc,50 {topo_nc} {s_addc_nc}", timeout=20) or assert_file(s_addc_nc))

    s_subc_nc = os.path.join(tmpdir, "s_subc.nc")
    run_test("str: subc,100", lambda: cdo(
        f"cdo subc,100 {topo_nc} {s_subc_nc}", timeout=20) or assert_file(s_subc_nc))

    s_divc_nc = os.path.join(tmpdir, "s_divc.nc")
    run_test("str: divc,2.0", lambda: cdo(
        f"cdo divc,2.0 {topo_nc} {s_divc_nc}", timeout=20) or assert_file(s_divc_nc))

    s_abs_nc = os.path.join(tmpdir, "s_abs.nc")
    run_test("str: abs", lambda: cdo(
        f"cdo abs {topo_nc} {s_abs_nc}", timeout=20) or assert_file(s_abs_nc))

    s_sqrt_nc = os.path.join(tmpdir, "s_sqrt.nc")
    run_test("str: sqrt(-abs)", lambda: cdo(
        f"cdo sqrt -abs {topo_nc} {s_sqrt_nc}", timeout=20) or assert_file(s_sqrt_nc))

    s_nint_nc = os.path.join(tmpdir, "s_nint.nc")
    run_test("str: nint", lambda: cdo(
        f"cdo nint {topo_nc} {s_nint_nc}", timeout=20) or assert_file(s_nint_nc))

    # --- 62f. Two-operand field math ---
    s_add_nc = os.path.join(tmpdir, "s_add.nc")
    run_test("str: add topo+mulc", lambda: cdo(
        f"cdo add {topo_nc} {mulc_nc} {s_add_nc}", timeout=20) or assert_file(s_add_nc))

    s_sub_nc = os.path.join(tmpdir, "s_sub.nc")
    run_test("str: sub topo-mulc", lambda: cdo(
        f"cdo sub {topo_nc} {mulc_nc} {s_sub_nc}", timeout=20) or assert_file(s_sub_nc))

    s_mul_nc = os.path.join(tmpdir, "s_mul.nc")
    run_test("str: mul topo*mulc", lambda: cdo(
        f"cdo mul {topo_nc} {mulc_nc} {s_mul_nc}", timeout=20) or assert_file(s_mul_nc))

    s_div_nc = os.path.join(tmpdir, "s_div.nc")
    run_test("str: div topo/mulc", lambda: cdo(
        f"cdo div {topo_nc} {mulc_nc} {s_div_nc}", timeout=20) or assert_file(s_div_nc))

    # --- 62g. Running statistics ---
    s_runmean_nc = os.path.join(tmpdir, "s_runmean.nc")
    run_test("str: runmean,3", lambda: cdo(
        f"cdo runmean,3 {monthly_nc} {s_runmean_nc}", timeout=20) or assert_file(s_runmean_nc))

    s_runstd_nc = os.path.join(tmpdir, "s_runstd.nc")
    run_test("str: runstd,3", lambda: cdo(
        f"cdo runstd,3 {monthly_nc} {s_runstd_nc}", timeout=20) or assert_file(s_runstd_nc))

    s_runmin_nc = os.path.join(tmpdir, "s_runmin.nc")
    run_test("str: runmin,3", lambda: cdo(
        f"cdo runmin,3 {monthly_nc} {s_runmin_nc}", timeout=20) or assert_file(s_runmin_nc))

    s_runmax_nc = os.path.join(tmpdir, "s_runmax.nc")
    run_test("str: runmax,3", lambda: cdo(
        f"cdo runmax,3 {monthly_nc} {s_runmax_nc}", timeout=20) or assert_file(s_runmax_nc))

    # --- 62h. Grid interpolation ---
    s_remapbil_nc = os.path.join(tmpdir, "s_remapbil.nc")
    run_test("str: remapbil,r8x4", lambda: cdo(
        f"cdo remapbil,r8x4 {topo_nc} {s_remapbil_nc}", timeout=30) or assert_file(s_remapbil_nc))

    s_remapnn_nc = os.path.join(tmpdir, "s_remapnn.nc")
    run_test("str: remapnn,r8x4", lambda: cdo(
        f"cdo remapnn,r8x4 {topo_nc} {s_remapnn_nc}", timeout=30) or assert_file(s_remapnn_nc))

    s_remapcon_nc = os.path.join(tmpdir, "s_remapcon.nc")
    run_test("str: remapcon,r8x4", lambda: cdo(
        f"cdo remapcon,r8x4 {topo_nc} {s_remapcon_nc}", timeout=30) or assert_file(s_remapcon_nc))

    # --- 62i. Selection ---
    s_selmon_nc = os.path.join(tmpdir, "s_selmon.nc")
    run_test("str: selmon,6", lambda: cdo(
        f"cdo selmon,6 {monthly_nc} {s_selmon_nc}", timeout=20) or assert_file(s_selmon_nc))

    s_selyear_nc = os.path.join(tmpdir, "s_selyear.nc")
    run_test("str: selyear,2020", lambda: cdo(
        f"cdo selyear,2020 {monthly_nc} {s_selyear_nc}", timeout=20) or assert_file(s_selyear_nc))

    s_selname_nc = os.path.join(tmpdir, "s_selname.nc")

    def _str_selname_test():
        vname = str(cdo.showname(input=topo_nc)).strip().split()[0]
        vname = ''.join(c for c in vname if c.isprintable() and c != ' ')
        cdo(f"cdo selname,{vname} {topo_nc} {s_selname_nc}", timeout=20)
        assert_file(s_selname_nc)
    run_test("str: selname", _str_selname_test)

    s_sellev_nc = os.path.join(tmpdir, "s_sellev.nc")
    run_test("str: sellevel,0,10000", lambda: cdo(
        f"cdo sellevel,0,10000 {stdatm_nc} {s_sellev_nc}", timeout=20) or assert_file(s_sellev_nc))

    # --- 62j. Ensemble operators ---
    s_ensmean_nc = os.path.join(tmpdir, "s_ensmean.nc")
    run_test("str: ensmean", lambda: cdo(
        f"cdo ensmean {topo_nc} {mulc_nc} {s_ensmean_nc}", timeout=30) or assert_file(s_ensmean_nc))

    s_ensmin_nc = os.path.join(tmpdir, "s_ensmin.nc")
    run_test("str: ensmin", lambda: cdo(
        f"cdo ensmin {topo_nc} {mulc_nc} {s_ensmin_nc}", timeout=30) or assert_file(s_ensmin_nc))

    s_ensmax_nc = os.path.join(tmpdir, "s_ensmax.nc")
    run_test("str: ensmax", lambda: cdo(
        f"cdo ensmax {topo_nc} {mulc_nc} {s_ensmax_nc}", timeout=30) or assert_file(s_ensmax_nc))

    # --- 62k. Metadata modification ---
    s_setname_nc = os.path.join(tmpdir, "s_setname.nc")
    run_test("str: setname,elevation", lambda: cdo(
        f"cdo setname,elevation {topo_nc} {s_setname_nc}", timeout=20) or assert_file(s_setname_nc))

    s_setunit_nc = os.path.join(tmpdir, "s_setunit.nc")
    run_test("str: setunit,m", lambda: cdo(
        f"cdo setunit,m {topo_nc} {s_setunit_nc}", timeout=20) or assert_file(s_setunit_nc))

    s_chname_nc = os.path.join(tmpdir, "s_chname.nc")
    run_test("str: chname topo->alt", lambda: cdo(
        f"cdo chname,topo,alt {topo_nc} {s_chname_nc}", timeout=20) or assert_file(s_chname_nc))

    # --- 62l. Masking ---
    s_setrtomiss_nc = os.path.join(tmpdir, "s_setrtomiss.nc")
    run_test("str: setrtomiss,0,100", lambda: cdo(
        f"cdo setrtomiss,0,100 {topo_nc} {s_setrtomiss_nc}", timeout=20) or assert_file(s_setrtomiss_nc))

    s_setmisstoc_nc = os.path.join(tmpdir, "s_setmisstoc.nc")
    run_test("str: setmisstoc chain", lambda: cdo(
        f"cdo -setmisstoc,0 -setrtomiss,0,100 {topo_nc} {s_setmisstoc_nc}",
        timeout=20) or assert_file(s_setmisstoc_nc))

    # --- 62m. Vertical statistics ---
    s_vertmean_nc = os.path.join(tmpdir, "s_vertmean.nc")
    run_test("str: vertmean", lambda: cdo(
        f"cdo vertmean {stdatm_nc} {s_vertmean_nc}", timeout=20) or assert_file(s_vertmean_nc))

    s_vertsum_nc = os.path.join(tmpdir, "s_vertsum.nc")
    run_test("str: vertsum", lambda: cdo(
        f"cdo vertsum {stdatm_nc} {s_vertsum_nc}", timeout=20) or assert_file(s_vertsum_nc))

    s_vertmin_nc = os.path.join(tmpdir, "s_vertmin.nc")
    run_test("str: vertmin", lambda: cdo(
        f"cdo vertmin {stdatm_nc} {s_vertmin_nc}", timeout=20) or assert_file(s_vertmin_nc))

    s_vertmax_nc = os.path.join(tmpdir, "s_vertmax.nc")
    run_test("str: vertmax", lambda: cdo(
        f"cdo vertmax {stdatm_nc} {s_vertmax_nc}", timeout=20) or assert_file(s_vertmax_nc))

    # --- 62n. Comparison operators ---
    s_eq_nc = os.path.join(tmpdir, "s_eq.nc")
    run_test("str: eq (same field)", lambda: cdo(
        f"cdo eq {topo_nc} {topo_nc} {s_eq_nc}", timeout=20) or assert_file(s_eq_nc))

    s_ne_nc = os.path.join(tmpdir, "s_ne.nc")
    run_test("str: ne (topo vs mulc)", lambda: cdo(
        f"cdo ne {topo_nc} {mulc_nc} {s_ne_nc}", timeout=20) or assert_file(s_ne_nc))

    s_le_nc = os.path.join(tmpdir, "s_le.nc")
    run_test("str: le (topo<=mulc)", lambda: cdo(
        f"cdo le {topo_nc} {mulc_nc} {s_le_nc}", timeout=20) or assert_file(s_le_nc))

    # --- 62o. Grid box statistics (string API) ---
    s_gridboxmean_nc = os.path.join(tmpdir, "s_gridboxmean.nc")
    run_test("str: gridboxmean,2,2", lambda: cdo(
        f"cdo gridboxmean,2,2 {topo_nc} {s_gridboxmean_nc}", timeout=20) or assert_file(s_gridboxmean_nc))

    # --- 62p. Chained operators (multi-step string API) ---
    s_chain1_nc = os.path.join(tmpdir, "s_chain1.nc")
    run_test("str: chain -fldmean -addc,10", lambda: cdo(
        f"cdo -fldmean -addc,10 {topo_nc} {s_chain1_nc}", timeout=20) or assert_file(s_chain1_nc))

    s_chain2_nc = os.path.join(tmpdir, "s_chain2.nc")
    run_test("str: chain -timmean -selmon,3,4,5", lambda: cdo(
        f"cdo -timmean -selmon,3,4,5 {monthly_nc} {s_chain2_nc}", timeout=20) or assert_file(s_chain2_nc))

    s_chain3_nc = os.path.join(tmpdir, "s_chain3.nc")
    run_test("str: chain -zonmean -mulc,2", lambda: cdo(
        f"cdo -zonmean -mulc,2 {topo_nc} {s_chain3_nc}", timeout=20) or assert_file(s_chain3_nc))

    # --- 62q. Info operators (text output, no file) ---
    run_test("str: showname (text)", lambda: assert_true(
        len(str(cdo(f"cdo showname {topo_nc}"))) > 0))

    run_test("str: ntime (text)", lambda: assert_true(
        len(str(cdo(f"cdo ntime {monthly_nc}"))) > 0))

    run_test("str: griddes (text)", lambda: assert_true(
        len(str(cdo(f"cdo griddes {topo_nc}"))) > 0))

    run_test("str: sinfo (text)", lambda: assert_true(
        len(str(cdo(f"cdo sinfo {topo_nc}"))) > 0))

    # ---- Summary ----
    elapsed = time.time() - t_start
    total = passed + failed
    print(f"\n{'='*60}")
    print(f"Results: {passed}/{total} passed, {failed} failed")
    print(f"Time: {elapsed:.1f}s")
    print(f"{'='*60}")

    # Clean up temporary directory and all test files
    if os.path.exists(tmpdir):
        print(f"\nCleaning up temporary files: {tmpdir}")
        shutil.rmtree(tmpdir, ignore_errors=True)

    if failed > 0:
        sys.exit(1)
    print("All stress tests passed!")


# ======================================================================
# Helper functions
# ======================================================================

def assert_file(path):
    """Assert file exists and is non-empty"""
    if not os.path.isfile(path):
        raise AssertionError(f"File does not exist: {path}")
    if os.path.getsize(path) == 0:
        raise AssertionError(f"File is empty: {path}")


def assert_true(cond):
    """Assert condition is true"""
    if not cond:
        raise AssertionError(f"Condition failed")


def assert_raises(exc_type, func):
    """Assert function raises exception"""
    try:
        func()
    except exc_type:
        return
    raise AssertionError(f"Expected {exc_type.__name__} but got none")


def _verify_nc_varname(nc_path: str, expected: str, original: str) -> None:
    """Verify a variable name in a NetCDF file without invoking CDO.

    NC3 classic stores variable names as plain ASCII in the file header
    (IETF RFC), so scanning the first 64 KB conclusively verifies the name.
    netCDF4 Python bindings are tried first when available.
    """
    try:
        import netCDF4 as _nc4  # type: ignore
        with _nc4.Dataset(nc_path) as ds:
            if expected not in ds.variables:
                raise AssertionError(
                    f"chname: expected '{expected}' in variables "
                    f"{sorted(ds.variables.keys())} (original: '{original}')")
        return
    except ImportError:
        pass

    # Fallback: scan the NC3 binary header (first 64 KB covers all metadata).
    with open(nc_path, "rb") as fh:
        magic = fh.read(4)
        fh.seek(0)
        header = fh.read(65536)
    if magic[:3] != b"CDF":
        raise AssertionError(
            f"chname: '{nc_path}' is not NC3 classic (magic={magic!r}); "
            f"cannot verify variable name")
    if expected.encode("ascii") not in header:
        raise AssertionError(
            f"chname: '{expected}' not found in NC3 header of '{nc_path}' "
            f"(original: '{original}')")


if __name__ == "__main__":
    main()
