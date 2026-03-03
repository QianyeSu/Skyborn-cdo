#!/usr/bin/env bash
# =============================================================================
# build_cdo_windows.sh - Compile CDO from source on Windows (MSYS2/MinGW64)
# =============================================================================
# Runs after build_deps_windows.sh has installed all dependencies.
# Applies Windows compatibility patches to CDO source before building.
#
# Environment variables:
#   CDO_SOURCE_DIR     - Path to CDO source (default: vendor/cdo)
#   CDO_DEPS_PREFIX    - Where dependencies are (default: /mingw64)
#   CDO_INSTALL_PREFIX - Where CDO will be installed (default: /opt/cdo-install)
#   PARALLEL_JOBS      - Number of parallel make jobs (default: nproc)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"

CDO_SOURCE="${CDO_SOURCE_DIR:-${PROJECT_DIR}/vendor/cdo}"
DEPS_PREFIX="${CDO_DEPS_PREFIX:-/mingw64}"
INSTALL_PREFIX="${CDO_INSTALL_PREFIX:-/opt/cdo-install}"
JOBS="${PARALLEL_JOBS:-$(nproc)}"

echo "============================================"
echo "Building CDO for Windows (MSYS2/MinGW64)"
echo "  Source:     ${CDO_SOURCE}"
echo "  Deps:       ${DEPS_PREFIX}"
echo "  Install to: ${INSTALL_PREFIX}"
echo "============================================"

cd "${CDO_SOURCE}"

# Apply Windows compatibility patches using intelligent Python script
echo "[skyborn-cdo] Applying Windows compatibility patches..."
python "${PROJECT_DIR}/scripts/patch_cdo_windows.py" apply --cdo-src "${CDO_SOURCE}"
if [ $? -ne 0 ]; then
    echo "[skyborn-cdo] ERROR: Patches failed to apply!"
    exit 1
fi

# Verify critical patches were applied
echo "[skyborn-cdo] Verifying patches..."
if ! grep -q '#ifdef __cplusplus' "${CDO_SOURCE}/src/process.h"; then
    echo "[skyborn-cdo] ERROR: __cplusplus guard not applied to process.h!"
    head -25 "${CDO_SOURCE}/src/process.h"
    exit 1
fi
if ! grep -q 'std::to_string' "${CDO_SOURCE}/src/mpmo_color.h"; then
    echo "[skyborn-cdo] ERROR: mpmo_color.h rewrite not applied!"
    head -90 "${CDO_SOURCE}/src/mpmo_color.h"
    exit 1
fi
if grep -q '#ifdef HAVE_LIBPTHREAD' "${CDO_SOURCE}/src/processManager.h"; then
    echo "[skyborn-cdo] processManager.h pthread guard: applied"
fi
echo "[skyborn-cdo] Patches verified OK"

# Defensive fix: ensure fftw3.h / <mutex> / fftwMutex are included unconditionally
# in the affected CDO source files.  The patch_cdo_windows.py script handles this
# via regex, but this inline block acts as a belt-and-suspenders fallback in case
# the regex doesn't match (e.g. different CDO snapshot with slightly different text).
for _fftw_file in \
    "src/operators/Filter.cc" \
    "src/operators/Fourier.cc" \
    "src/cdo_fctrans.cc"; do
  if [[ -f "${CDO_SOURCE}/${_fftw_file}" ]]; then
    python - "${CDO_SOURCE}/${_fftw_file}" <<'PY'
import sys, re
fpath = sys.argv[1]
text = open(fpath, encoding="utf-8", errors="ignore").read()
original = text
# Remove the first #ifdef HAVE_LIBFFTW3 block that contains fftw3.h/mutex/fftwMutex.
# Replace with unconditional includes so fftw_complex is always visible.
text = re.sub(
    r'#ifdef HAVE_LIBFFTW3\s*\n#include <fftw3\.h>\s*\n(?:#include <mutex>\s*\n)?(?:static std::mutex fftwMutex;\s*\n)?#endif',
    '#include <fftw3.h>\n#include <mutex>\nstatic std::mutex fftwMutex;',
    text, flags=re.MULTILINE)
# Also handle the split-block form used in cdo_fctrans.cc:
# Block 1: just #include <fftw3.h>
text = re.sub(
    r'#ifdef HAVE_LIBFFTW3\s*\n(#include <fftw3\.h>)\s*\n#endif',
    r'\1', text, flags=re.MULTILINE)
# Block 2: mutex + fftwMutex
text = re.sub(
    r'#ifdef HAVE_LIBFFTW3\s*\n(#include <mutex>\s*\nstatic std::mutex fftwMutex;)\s*\n#endif',
    r'\1', text, flags=re.MULTILINE)
if text != original:
    open(fpath, 'w', encoding="utf-8", newline="\n").write(text)
    print(f"[skyborn-cdo] Defensive fftw3 fix applied: {fpath}")
PY
  fi
done

# Defensive fix: rename include guards in both table.h files to prevent
# the TABLE_H guard clash that hides cdo::define_table from Setpartab.cc.
# patch_cdo_windows.py handles this, but the inline fallback ensures it
# works even if the patcher pattern doesn't match a slightly different snapshot.
for _tbl_file in \
    "src/table.h" \
    "libcdi/src/table.h"; do
  if [[ -f "${CDO_SOURCE}/${_tbl_file}" ]]; then
    python - "${CDO_SOURCE}/${_tbl_file}" "${_tbl_file}" <<'PY'
import sys
fpath, rel = sys.argv[1], sys.argv[2]
text = open(fpath, encoding="utf-8", errors="ignore").read()
original = text
if rel == "src/table.h":
    # CDO's own table.h: use a unique guard so it is never skipped
    text = text.replace("#ifndef TABLE_H\n#define TABLE_H",
                        "#ifndef CDO_SRC_TABLE_H\n#define CDO_SRC_TABLE_H")
else:
    # libcdi's auto-generated table.h: rename guard to avoid polluting TABLE_H
    text = text.replace("#ifndef TABLE_H\n#define TABLE_H",
                        "#ifndef CDI_TABLE_H\n#define CDI_TABLE_H")
if text != original:
    open(fpath, 'w', encoding="utf-8", newline="\n").write(text)
    print(f"[skyborn-cdo] Defensive table.h guard fix applied: {fpath}")
PY
  fi
done

# Prevent make from trying to regenerate autotools files.
# The vendored source includes pre-generated configure/Makefile.in/aclocal.m4,
# but git checkout sets all timestamps to the same time, which can cause make
# to think the generated files are stale and try to re-run aclocal/autoconf.
# Touch source files first, then generated files 1s later to ensure correct ordering.
echo "[skyborn-cdo] Fixing autotools timestamps..."
find . -name 'configure.ac' -exec touch {} +
find . -name 'Makefile.am' -exec touch {} +
sleep 1
find . -name aclocal.m4 -exec touch {} +
find . -name configure -exec touch {} +
find . -name Makefile.in -exec touch {} +
find . -name config.h.in -exec touch {} +

# Ensure configure and autotools helper scripts are executable (they may be stored as 644 in git)
find . -name configure -exec chmod +x {} +
find . \( -name 'config.sub' -o -name 'config.guess' -o -name 'install-sh' \
    -o -name 'missing' -o -name 'compile' -o -name 'depcomp' \
    -o -name 'ltmain.sh' -o -name 'test-driver' \) -exec chmod +x {} +

# If configure doesn't exist, run autoreconf
if [[ ! -f configure ]]; then
    echo "[skyborn-cdo] Running autoreconf..."
    autoreconf -fvi
fi

export PKG_CONFIG_PATH="${DEPS_PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export PATH="${DEPS_PREFIX}/bin:${PATH}"

echo "[skyborn-cdo] Configuring CDO for Windows..."

# NOTE: Do NOT use --host=x86_64-w64-mingw32 here!
# In MSYS2/MinGW64, the compiler is already the MinGW cross-compiler and
# produces native Windows executables.  Adding --host makes autotools treat
# this as a cross-compilation, which disables AC_RUN_IFELSE runtime tests.
# Those skipped tests cause CDO's NetCDF format-selection code path (-f nc*)
# to be mis-configured, leading to hangs on Windows.
./configure \
    --prefix="${INSTALL_PREFIX}" \
    --with-netcdf="${DEPS_PREFIX}" \
    --with-hdf5="${DEPS_PREFIX}" \
    --with-eccodes="${DEPS_PREFIX}" \
    --with-fftw3 \
    --with-proj="${DEPS_PREFIX}" \
    --with-udunits2="${DEPS_PREFIX}" \
    --with-szlib="${DEPS_PREFIX}" \
    --with-threads=no \
    --disable-fortran \
    --disable-across \
    --disable-custom-modules \
    --enable-cgribex \
    CFLAGS="-O2 -I${DEPS_PREFIX}/include" \
    CXXFLAGS="-D_USE_MATH_DEFINES -O2 -std=c++20 -Wno-template-body -I${DEPS_PREFIX}/include" \
    CPPFLAGS="-I${DEPS_PREFIX}/include" \
    LDFLAGS="-L${DEPS_PREFIX}/lib" \
    LIBS="-lz -lm -lws2_32 -lrpcrt4"

# GCC 15 + MinGW: avoid header shadowing of system <process.h> by
# local src/process.h (same basename). In the generated src/Makefile,
# DEFAULT_INCLUDES typically contains -I. and -I$(srcdir), which lets
# <process.h> resolve to CDO's C++ header during libstdc++ gthread include.
# Force quote-only include search for src headers to prevent this collision.
if [[ -f src/Makefile ]]; then
    echo "[skyborn-cdo] Patching src/Makefile include search to avoid process.h collision..."
    sed -i 's|^DEFAULT_INCLUDES = .*|DEFAULT_INCLUDES = -iquote $(srcdir) -idirafter $(srcdir)|' src/Makefile
fi

echo "[skyborn-cdo] Building CDO..."
# Disable libcdi tests that use POSIX-only functions (srand48, lrand48, etc.)
# not available on MinGW/Windows
if [[ -f libcdi/Makefile ]]; then
    sed -i 's/ tests$//; s/ tests / /g' libcdi/Makefile
fi
make -j"${JOBS}"

echo "[skyborn-cdo] Installing CDO..."
make install

echo "============================================"
echo "CDO Windows build complete!"
echo "============================================"

if [[ -f "${INSTALL_PREFIX}/bin/cdo.exe" ]]; then
    echo "Binary: ${INSTALL_PREFIX}/bin/cdo.exe"
    echo "Size:   $(du -h "${INSTALL_PREFIX}/bin/cdo.exe" | cut -f1)"
    echo ""
    echo "--- DLL dependencies ---"
    ldd "${INSTALL_PREFIX}/bin/cdo.exe" || true
fi
