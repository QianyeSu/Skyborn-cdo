#!/bin/bash
# Local Windows CDO build using MSYS2 MinGW-w64
# Run this from MSYS2 MinGW64 shell:
#   /c/msys64/mingw64.exe  (or from PowerShell: C:\msys64\usr\bin\bash.exe -lc "...")
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CDO_SRC="$PROJECT_DIR/vendor/cdo"
BUILD_DIR="$PROJECT_DIR/build_win64"
PREFIX="$PROJECT_DIR/install_win64"

echo "=== CDO Windows (MinGW-w64) Build ==="
echo "Source:  $CDO_SRC"
echo "Build:   $BUILD_DIR"
echo "Install: $PREFIX"

# Clean previous build
rm -rf "$BUILD_DIR" "$PREFIX"
mkdir -p "$BUILD_DIR" "$PREFIX"

cd "$CDO_SRC"

# Apply Windows compatibility patches using intelligent Python script
echo ""
echo "=== Applying Windows compatibility patches ==="
python "${PROJECT_DIR}/scripts/patch_cdo_windows.py" apply --cdo-src "$CDO_SRC"
if [ $? -ne 0 ]; then
    echo "Warning: Some patches failed to apply"
fi

# Defensive fix for snapshots where Filter.cc uses std::scoped_lock without
# stable mutex include/fftwMutex declaration.
if [ -f "$CDO_SRC/src/operators/Filter.cc" ]; then
    python - <<'PY'
from pathlib import Path
path = Path('vendor/cdo/src/operators/Filter.cc')
text = path.read_text(encoding='utf-8', errors='ignore')
changed = False
if 'std::scoped_lock' in text and '#include <mutex>' not in text:
    if '#include "field_functions.h"\n' in text:
        text = text.replace('#include "field_functions.h"\n', '#include "field_functions.h"\n#include <mutex>\n', 1)
        changed = True
    elif '#include <fftw3.h>\n' in text:
        text = text.replace('#include <fftw3.h>\n', '#include <fftw3.h>\n#include <mutex>\n', 1)
        changed = True

if 'std::scoped_lock' in text and 'fftwMutex' in text and 'static std::mutex fftwMutex;' not in text:
    anchor = '#include "field_functions.h"\n'
    if anchor in text:
        text = text.replace(anchor, anchor + '\nstatic std::mutex fftwMutex;\n', 1)
        changed = True

if changed:
    path.write_text(text, encoding='utf-8', newline='\n')
PY
fi

if [ -f "$CDO_SRC/src/operators/Fourier.cc" ]; then
    python - <<'PY'
from pathlib import Path
path = Path('vendor/cdo/src/operators/Fourier.cc')
text = path.read_text(encoding='utf-8', errors='ignore')
changed = False
if 'std::scoped_lock' in text and '#include <mutex>' not in text:
    if '#include "field_functions.h"\n' in text:
        text = text.replace('#include "field_functions.h"\n', '#include "field_functions.h"\n#include <mutex>\n', 1)
        changed = True
    elif '#include <fftw3.h>\n' in text:
        text = text.replace('#include <fftw3.h>\n', '#include <fftw3.h>\n#include <mutex>\n', 1)
        changed = True

if 'std::scoped_lock' in text and 'fftwMutex' in text and 'static std::mutex fftwMutex;' not in text:
    anchor = '#include "field_functions.h"\n'
    if anchor in text:
        text = text.replace(anchor, anchor + '\nstatic std::mutex fftwMutex;\n', 1)
        changed = True

if changed:
    path.write_text(text, encoding='utf-8', newline='\n')
PY
fi

# Check if configure exists (pre-generated)
if [ ! -f configure ]; then
    echo "Running autoreconf..."
    autoreconf -fiv
fi

# Configure CDO
cd "$BUILD_DIR"
echo ""
echo "=== Configuring CDO ==="
"$CDO_SRC/configure" \
    --prefix="$PREFIX" \
    --with-netcdf=/mingw64 \
    --with-hdf5=/mingw64 \
    --with-eccodes=/mingw64 \
    --with-fftw3 \
    --with-proj=/mingw64 \
    --with-udunits2=/mingw64 \
    --disable-custom-modules \
    --with-threads=yes \
    CFLAGS="-O2 -I/mingw64/include" \
    CXXFLAGS="-O2 -std=c++20 -I/mingw64/include" \
    CPPFLAGS="-I/mingw64/include" \
    LDFLAGS="-L/mingw64/lib" \
    LIBS="-lws2_32"

# GCC 15 + MinGW: prevent local src/process.h from shadowing system <process.h>
# during libstdc++ gthread includes.
if [ -f "$BUILD_DIR/src/Makefile" ]; then
    sed -i 's|^DEFAULT_INCLUDES = .*|DEFAULT_INCLUDES = -iquote $(srcdir) -idirafter $(srcdir)|' "$BUILD_DIR/src/Makefile"
fi

echo ""
echo "=== Building CDO ==="
make -j$(nproc)

echo ""
echo "=== Installing CDO ==="
make install

echo ""
echo "=== Build Complete ==="
echo "CDO binary: $PREFIX/bin/cdo.exe"
ls -la "$PREFIX/bin/"

echo ""
echo "=== Testing CDO ==="
"$PREFIX/bin/cdo.exe" --version || true

# Restore vendor directory to clean state (undo patch modifications)
echo ""
echo "=== Restoring vendor directory ==="
python "${PROJECT_DIR}/scripts/patch_cdo_windows.py" restore --cdo-src "$CDO_SRC"
