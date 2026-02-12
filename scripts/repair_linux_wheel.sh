#!/usr/bin/env bash
# ===========================================================================
# repair_linux_wheel.sh — auditwheel repair with CDO binary protection
# ===========================================================================
# auditwheel + patchelf can corrupt ELF executables (as opposed to shared
# libraries).  The CDO binary often has many PT_NOTE segments and debug
# info that confuse patchelf, leading to a broken binary (ldd shows nothing,
# running it yields SIGSEGV / exit -11).
#
# Strategy:
#   1. Strip the CDO binary *before* auditwheel runs (removes excess notes
#      and debug info that trigger patchelf bugs).
#   2. Run auditwheel repair normally.
#   3. Verify the repaired binary is functional (ldd, --version).
#   4. If verification fails, replace the corrupted binary with the
#      pre-repair copy and manually set RPATH via patchelf.
# ===========================================================================
set -euo pipefail

DEST_DIR="$1"   # {dest_dir}
WHEEL="$2"      # {wheel}

echo "=== repair_linux_wheel.sh ==="
echo "  wheel:    $WHEEL"
echo "  dest_dir: $DEST_DIR"

# -------------------------------------------------------------------
# 1.  Strip the CDO binary inside the wheel before auditwheel touches it
# -------------------------------------------------------------------
# Unpack wheel to a temp directory, strip the binary, repack.
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

unzip -q "$WHEEL" -d "$WORK/unpacked"

CDO_BIN=$(find "$WORK/unpacked" -type f -name "cdo" -path "*/bin/cdo" | head -1)
if [ -n "$CDO_BIN" ]; then
    echo "  Found CDO binary: $CDO_BIN"

    # Save a pristine copy before any tool touches it
    cp "$CDO_BIN" "$WORK/cdo_pristine"

    # Full strip — removes debug info, notes, symbol tables
    strip -s "$CDO_BIN" 2>/dev/null || strip --strip-all "$CDO_BIN" 2>/dev/null || true

    # Also strip excess ELF notes that confuse patchelf
    # (objcopy --remove-section works even when strip alone isn't enough)
    objcopy --remove-section=.note.gnu.property \
            --remove-section=.note.gnu.build-id \
            --remove-section=.note.ABI-tag \
            "$CDO_BIN" 2>/dev/null || true

    echo "  Stripped CDO binary ($(stat -c%s "$CDO_BIN") bytes)"

    # Repack into a new wheel
    STRIPPED_WHEEL="$WORK/$(basename "$WHEEL")"
    (cd "$WORK/unpacked" && zip -q -r "$STRIPPED_WHEEL" .)
    WHEEL="$STRIPPED_WHEEL"
    echo "  Repacked stripped wheel"
else
    echo "  WARNING: CDO binary not found in wheel, proceeding without strip"
fi

# -------------------------------------------------------------------
# 2.  Run auditwheel repair
# -------------------------------------------------------------------
echo "  Running auditwheel repair..."
LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}" \
    auditwheel repair -w "$DEST_DIR" "$WHEEL"

# -------------------------------------------------------------------
# 3.  Verify the repaired wheel
# -------------------------------------------------------------------
echo "  Verifying repaired wheel..."
REPAIRED_WHEEL=$(ls -1 "$DEST_DIR"/*.whl | head -1)
VERIFY=$(mktemp -d)

unzip -q "$REPAIRED_WHEEL" -d "$VERIFY"

REPAIRED_CDO=$(find "$VERIFY" -type f -name "cdo" -path "*/bin/cdo" | head -1)
if [ -z "$REPAIRED_CDO" ]; then
    echo "  ERROR: CDO binary not found in repaired wheel!"
    rm -rf "$VERIFY"
    exit 1
fi

# Determine library paths inside the wheel for verification
PKG_LIB=$(dirname "$(dirname "$REPAIRED_CDO")")/lib
PKG_PARENT=$(dirname "$(dirname "$REPAIRED_CDO")")
AUDITWHEEL_LIBS=$(dirname "$PKG_PARENT")/skyborn_cdo.libs

VERIFY_LD_PATH="$PKG_LIB"
[ -d "$AUDITWHEEL_LIBS" ] && VERIFY_LD_PATH="$VERIFY_LD_PATH:$AUDITWHEEL_LIBS"

# Test 1: ldd should produce output
LDD_OUT=$(LD_LIBRARY_PATH="$VERIFY_LD_PATH" ldd "$REPAIRED_CDO" 2>&1 || true)
LDD_LINES=$(echo "$LDD_OUT" | grep -c "=>" || true)

# Test 2: --version should work (exit 0)
VERSION_RC=0
LD_LIBRARY_PATH="$VERIFY_LD_PATH" "$REPAIRED_CDO" --version >/dev/null 2>&1 || VERSION_RC=$?

echo "  ldd lines: $LDD_LINES, --version rc: $VERSION_RC"

if [ "$LDD_LINES" -gt 0 ] && [ "$VERSION_RC" -eq 0 ]; then
    echo "  Verification PASSED — repaired binary is functional"
    rm -rf "$VERIFY"
    exit 0
fi

# -------------------------------------------------------------------
# 4.  Repair failed — restore pristine binary + manual RPATH
# -------------------------------------------------------------------
echo "  WARNING: Repaired binary appears broken, restoring pristine copy..."

if [ ! -f "$WORK/cdo_pristine" ]; then
    echo "  ERROR: No pristine binary saved, cannot recover"
    rm -rf "$VERIFY"
    exit 1
fi

# Calculate the RPATH the binary needs inside the installed wheel
# bin/cdo -> needs $ORIGIN/../lib and $ORIGIN/../../skyborn_cdo.libs
RPATH='$ORIGIN/../lib:$ORIGIN/../../skyborn_cdo.libs'

# Restore pristine binary
cp "$WORK/cdo_pristine" "$REPAIRED_CDO"
chmod +x "$REPAIRED_CDO"

# Set RPATH with patchelf (on the pristine copy, which is safe)
if command -v patchelf >/dev/null 2>&1; then
    patchelf --set-rpath "$RPATH" "$REPAIRED_CDO" 2>/dev/null || true
    echo "  Set RPATH to: $RPATH"
fi

# Verify again
VERSION_RC2=0
LD_LIBRARY_PATH="$VERIFY_LD_PATH" "$REPAIRED_CDO" --version >/dev/null 2>&1 || VERSION_RC2=$?
echo "  Post-restore --version rc: $VERSION_RC2"

# Repack the fixed wheel
rm "$REPAIRED_WHEEL"
(cd "$VERIFY" && zip -q -r "$REPAIRED_WHEEL" .)
echo "  Repacked wheel with restored binary"

rm -rf "$VERIFY"

if [ "$VERSION_RC2" -ne 0 ]; then
    echo "  WARNING: Binary still fails after restore (rc=$VERSION_RC2)"
    echo "  The wheel may not work correctly at runtime."
fi

exit 0
