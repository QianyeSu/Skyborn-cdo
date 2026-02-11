#!/bin/bash
# ============================================================================
# fix_linux_deps.sh — Bundle missing runtime libraries for manylinux wheel
# ============================================================================
# Copies non-baseline system libraries that CDO and its deps need into the
# package lib/ directory.  This prevents auditwheel from vendoring them into
# .libs/ with hashed filenames, which can corrupt DT_NEEDED resolution for
# the CDO executable (auditwheel patches .so RPATH/DT_NEEDED but may corrupt
# ELF executables when combined with prior patchelf calls).
#
# Usage: bash fix_linux_deps.sh <package_dir> [<cdo_bin>]
#   <package_dir>: path to skyborn_cdo/ package (contains bin/ and lib/)
#   <cdo_bin>:     optional, path to original CDO binary for ldd analysis
#                  (defaults to <package_dir>/bin/cdo)
# ============================================================================

set -euo pipefail

PKG_DIR="${1:?Usage: $0 <package_dir> [<cdo_bin>]}"
CDO_BIN="${2:-$PKG_DIR/bin/cdo}"
LIB_DIR="$PKG_DIR/lib"

if [ ! -f "$CDO_BIN" ]; then
    echo "ERROR: CDO binary not found at $CDO_BIN"
    exit 1
fi

mkdir -p "$LIB_DIR"

# ---- manylinux_2_28 baseline libraries (PEP 600) — MUST NOT be bundled ----
# These are guaranteed to exist on any manylinux_2_28-compatible system.
BASELINE='linux-vdso|ld-linux|/ld[.-]|/libc\.|/libm\.|/libpthread|/libdl\.|/librt\.|/libstdc\+\+|/libgcc_s|/libresolv|/libz\.|/libutil|/libnsl|/libcrypt|/libBrokenLocale|/libSegFault|/libthread_db|/libanl|/libmvec|/libnss|/libmcheck'

copy_missing_libs() {
    local binary="$1"
    local label="$2"
    ldd "$binary" 2>/dev/null | \
        grep '=>' | \
        grep -vE "$BASELINE" | \
        awk '$3 ~ /^\// {print $3}' | sort -u | \
        while read -r lib; do
            if [ -f "$lib" ]; then
                local base
                base=$(basename "$lib")
                if [ ! -e "$LIB_DIR/$base" ]; then
                    echo "  + $base  <-  $lib  ($label)"
                    cp -L "$lib" "$LIB_DIR/" 2>/dev/null || true
                fi
            fi
        done
}

echo "=== Collecting CDO runtime dependencies ==="
copy_missing_libs "$CDO_BIN" "cdo"

echo "=== Collecting .so transitive dependencies ==="
find "$LIB_DIR" -name '*.so*' -type f 2>/dev/null | while read -r so; do
    copy_missing_libs "$so" "$(basename "$so")"
done

# Second pass: some deps may have pulled in new .so files that have
# their own additional dependencies.
echo "=== Second pass (transitive) ==="
find "$LIB_DIR" -name '*.so*' -type f 2>/dev/null | while read -r so; do
    copy_missing_libs "$so" "$(basename "$so") [2nd]"
done

NFILES=$(find "$LIB_DIR" -name '*.so*' | wc -l)
echo "=== Library directory: $NFILES shared objects ==="
