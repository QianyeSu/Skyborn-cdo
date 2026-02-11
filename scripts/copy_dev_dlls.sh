#!/bin/bash
# Copy CDO DLL dependencies from MinGW64 to src/skyborn_cdo/bin for development
set -e

CDO_BIN="/f/skyborn-cdo-package/src/skyborn_cdo/bin"
MINGW_BIN="/mingw64/bin"

echo "=== Iteratively resolving ALL DLL dependencies ==="

for pass in 1 2 3 4 5; do
    # Collect all missing DLLs from all exe/dll files in CDO_BIN
    MISSING=""
    for f in "$CDO_BIN"/*.exe "$CDO_BIN"/*.dll; do
        [ -f "$f" ] || continue
        FOUND=$(ldd "$f" 2>/dev/null | grep "not found" | awk '{print $1}')
        if [ -n "$FOUND" ]; then
            MISSING="$MISSING $FOUND"
        fi
    done
    # Deduplicate
    MISSING=$(echo "$MISSING" | tr ' ' '\n' | sort -u | grep -v '^$')

    if [ -z "$MISSING" ]; then
        echo "Pass $pass: All dependencies resolved!"
        break
    fi

    echo "Pass $pass: Resolving missing DLLs..."
    COPIED=0
    for dll in $MISSING; do
        src="$MINGW_BIN/$dll"
        if [ -f "$src" ] && [ ! -f "$CDO_BIN/$dll" ]; then
            echo "  Copying $dll"
            cp "$src" "$CDO_BIN/"
            COPIED=$((COPIED + 1))
        fi
    done
    if [ "$COPIED" -eq 0 ]; then
        echo "  No new DLLs to copy. Remaining missing:"
        echo "$MISSING"
        break
    fi
done

echo ""
echo "=== Testing cdo.exe ==="
"$CDO_BIN/cdo.exe" --version 2>&1 | head -5 || echo "cdo.exe failed with exit code $?"
