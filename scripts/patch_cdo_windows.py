#!/usr/bin/env python3
"""
Smart Windows Compatibility Patcher for CDO Source Code

Uses code pattern matching (instead of fixed line numbers) to modify CDO source
for Windows compilation. More resilient to CDO version updates.

Usage:
    python patch_cdo_windows.py apply [--cdo-src PATH]    # Apply patches
    python patch_cdo_windows.py restore [--cdo-src PATH]  # Restore original
    python patch_cdo_windows.py verify [--cdo-src PATH]   # Verify (dry-run)
"""

import re
import sys
import argparse
from pathlib import Path
from typing import Tuple


class WindowsPatcher:
    """CDO Windows compatibility patch manager"""

    def __init__(self, cdo_src_dir: Path):
        self.cdo_src = Path(cdo_src_dir).resolve()
        self.backup_dir = self.cdo_src / ".patch_backup"

    def _backup_map_path(self) -> Path:
        return self.backup_dir / "_path_map.txt"

    def _save_backup_mapping(self, rel_path: str, backup_name: str):
        """Append a rel_path -> backup_name mapping to the map file"""
        with open(self._backup_map_path(), 'a', encoding='utf-8') as f:
            f.write(f"{backup_name}\t{rel_path}\n")

    def _load_backup_mappings(self) -> dict:
        """Load backup_name -> rel_path mappings from the map file"""
        mappings = {}
        map_file = self._backup_map_path()
        if map_file.exists():
            for line in map_file.read_text(encoding='utf-8').splitlines():
                if '\t' in line:
                    backup_name, rel_path = line.split('\t', 1)
                    mappings[backup_name] = rel_path
        return mappings

    def patch_file(self, rel_path: str, patches: list, dry_run: bool = False) -> Tuple[bool, int]:
        """Apply patches to a single file. Returns (success, count)"""
        file_path = self.cdo_src / rel_path

        if not file_path.exists():
            print(f"[X] {rel_path}: File not found")
            return False, 0

        print(f"[*] {rel_path}")

        try:
            content = file_path.read_text(encoding='utf-8', errors='ignore')
        except Exception as e:
            print(f"   [X] Read failed: {e}")
            return False, 0

        original = content
        applied_count = 0

        for desc, pattern, replacement in patches:
            if isinstance(pattern, str):
                # Simple string replacement
                if pattern in content:
                    content = content.replace(pattern, replacement, 1)
                    print(f"   [+] {desc}")
                    applied_count += 1
                else:
                    print(f"   [ ] Not found: {desc}")
            else:
                # Regex replacement
                new_content, count = pattern.subn(
                    replacement, content, count=1)
                if count > 0:
                    content = new_content
                    print(f"   [+] {desc}")
                    applied_count += 1
                else:
                    print(f"   [ ] Not found: {desc}")

        # Write modifications
        if content != original and not dry_run:
            # Backup original file
            backup_name = rel_path.replace('/', '__').replace('\\', '__')
            backup_path = self.backup_dir / backup_name
            backup_path.parent.mkdir(parents=True, exist_ok=True)
            backup_path.write_text(original, encoding='utf-8', newline='\n')
            self._save_backup_mapping(rel_path, backup_name)

            # Write modified content
            file_path.write_text(content, encoding='utf-8', newline='\n')

        return applied_count > 0, applied_count

    def apply_all(self, dry_run: bool = False) -> bool:
        """Apply all patches"""
        print(f"CDO source directory: {self.cdo_src}\n")

        if not dry_run:
            self.backup_dir.mkdir(exist_ok=True)

        total_files = 0
        total_patches = 0

        # =================================================================
        # Patch definitions: based on code pattern matching
        # =================================================================

        patches = [
            # --- src/cdo.cc ---
            ("src/cdo.cc", [
                ("Add Windows headers (io.h, windows.h)",
                 re.compile(r'^(#include\s+<unistd\.h>.*?)$', re.MULTILINE),
                 r'#ifdef _WIN32\n#include <io.h>\n#include <windows.h>\n#else\n\1\n#endif'),

                ("cdo_init_is_tty: Windows implementation",
                 re.compile(
                     r'(static\s+void\s+cdo_init_is_tty\s*\(\s*\)\s*\{\s*\n)'
                     r'(\s+struct stat statbuf;\s*\n'
                     r'\s+fstat\(0, &statbuf\);\s*\n'
                     r'\s+if \(S_ISCHR\(statbuf\.st_mode\)\)[^}]+\}\s*\n'
                     r'\s+fstat\(1, &statbuf\);\s*\n'
                     r'[^}]+stdoutIsTerminal[^}]+\}\s*\n'
                     r'\s+fstat\(2, &statbuf\);\s*\n'
                     r'[^}]+stderrIsTerminal[^}]+\}\s*\n)'
                     r'(\})',
                     re.MULTILINE
                 ),
                 r'\1#ifdef _WIN32\n'
                 r'  cdo::stdinIsTerminal = _isatty(_fileno(stdin));\n'
                 r'  cdo::stdoutIsTerminal = _isatty(_fileno(stdout));\n'
                 r'  cdo::stderrIsTerminal = _isatty(_fileno(stderr));\n'
                 r'#else\n'
                 r'\2#endif\n'
                 r'\3'),

                ("Add fflush before clear_processes",
                 re.compile(
                     r'(\s+)(g_processManager\.clear_processes\s*\(\s*\)\s*;)'),
                 r'\1fflush(stdout);\n\1\2'),
            ]),

            # --- src/cdo_getopt.cc ---
            ("src/cdo_getopt.cc", [
                ("Add windows.h for console API on Windows",
                 re.compile(r'^(#include\s+<map>)$', re.MULTILINE),
                 r'\1\n\n#ifdef _WIN32\n#include <windows.h>\n#endif'),

                ("Guard sys/ioctl.h (not available on Windows)",
                 re.compile(r'^(#include\s+<sys/ioctl\.h>)$', re.MULTILINE),
                 r'#ifndef _WIN32\n\1\n#endif'),

                ("Guard unistd.h",
                 re.compile(r'^(#include\s+<unistd\.h>)$', re.MULTILINE),
                 r'#ifndef _WIN32\n\1\n#endif'),
            ]),

            # --- src/dcw_reader.cc: Windows access()/R_OK compatibility ---
            ("src/dcw_reader.cc", [
                ("Add io.h/access fallback on Windows",
                 re.compile(r'^(#include\s+<cstring>)$', re.MULTILINE),
                 r'\1\n#ifdef _WIN32\n#include <io.h>\n#ifndef R_OK\n#define R_OK 4\n#endif\n#ifndef access\n#define access _access\n#endif\n#endif'),
            ]),

            # --- src/process.h: guard C++ content and pthread usage ---
            # 1. Guard C++ content from C compiler: MinGW's unistd.h includes
            #    "process.h" (Windows process API), which resolves to CDO's
            #    src/process.h via -I src. CDO's process.h is C++ only.
            # 2. Guard pthread_t start_thread() with HAVE_LIBPTHREAD.
            #    In MSYS2 MINGW64, pthreads is provided by winpthreads (POSIX
            #    thread model). This guard ensures graceful degradation if
            #    the build is done on a system without pthreads.
            ("src/process.h", [
                ("Guard C++ content from C compiler",
                 re.compile(
                     r'(#ifndef PROCESS_H\s*\n'
                     r'#define PROCESS_H\s*\n)',
                     re.MULTILINE
                 ),
                 r'\1\n#ifdef __cplusplus\n'),

                ("Guard pthread_t start_thread() with HAVE_LIBPTHREAD",
                 re.compile(
                     r'^(\s+)(pthread_t start_thread\(\);)$',
                     re.MULTILINE
                 ),
                 r'#ifdef HAVE_LIBPTHREAD\n\1\2\n#endif'),

                ("Close C++ guard at end",
                 "#endif /* PROCESS_H */",
                 "#endif /* __cplusplus */\n#endif /* PROCESS_H */"),
            ]),

            # --- src/processManager.h: guard pthread includes ---
            # processManager.h unconditionally includes <pthread.h>.
            # In MSYS2 MINGW64, POSIX threads (winpthreads) are used -- NOT MCF.
            # MCF is only the default in UCRT64/CLANG64 environments.
            # We move <pthread.h> to BEFORE the C++ STL includes as a defensive
            # measure (belt-and-suspenders), in case the include order matters
            # for edge cases in future GCC versions.
            ("src/processManager.h", [
                ("Add pthread.h before C++ STL headers (safer include order)",
                 "// Stdlib includes\n#include <map>",
                 "#ifdef HAVE_LIBPTHREAD\n"
                 "// Must be included before C++ STL to ensure correct threading setup\n"
                 "#include <pthread.h>\n"
                 "#endif\n\n"
                 "// Stdlib includes\n#include <map>"),

                ("Guard existing pthread.h include with HAVE_LIBPTHREAD",
                 re.compile(r'^(#include <pthread\.h>)$', re.MULTILINE),
                 r'#ifdef HAVE_LIBPTHREAD\n\1\n#endif'),

                ("Guard m_threadIDs vector",
                 re.compile(
                     r'^(\s+)(std::vector<pthread_t> m_threadIDs;)$',
                     re.MULTILINE
                 ),
                 r'#ifdef HAVE_LIBPTHREAD\n\1\2\n#endif'),
            ]),

            # --- src/processManager.cc: guard pthread function calls ---
            ("src/processManager.cc", [
                ("Guard start_thread() call in run_processes",
                 '      m_threadIDs.push_back(idProcessPair.second->start_thread());',
                 '#ifdef HAVE_LIBPTHREAD\n      m_threadIDs.push_back(idProcessPair.second->start_thread());\n#endif'),

                ("Guard pthread_self() push in run_processes",
                 '  m_threadIDs.push_back(pthread_self());',
                 '#ifdef HAVE_LIBPTHREAD\n  m_threadIDs.push_back(pthread_self());\n#endif'),

                ("Guard kill_processes pthread calls",
                 re.compile(
                     r'(ProcessManager::kill_processes\(\)\s*\n\{)\s*\n'
                     r'(  for \(auto threadID.*?)\n(\})',
                     re.DOTALL
                 ),
                 r'\1\n#ifdef HAVE_LIBPTHREAD\n\2\n#endif\n\3'),
            ]),

            # --- src/process.cc: guard start_thread() implementation ---
            ("src/process.cc", [
                ("Guard start_thread() implementation with HAVE_LIBPTHREAD",
                 re.compile(
                     r'^(pthread_t\nProcess::start_thread\(\))',
                     re.MULTILINE
                 ),
                 r'#ifdef HAVE_LIBPTHREAD\n\1'),

                ("Close HAVE_LIBPTHREAD guard after start_thread()",
                 re.compile(
                     r'(  return thrID;\n\})',
                     re.MULTILINE
                 ),
                 r'\1\n#endif /* HAVE_LIBPTHREAD */'),
            ]),

            # --- src/varray.h: add missing <iostream> for std::cerr ---
            # GCC 15 no longer transitively includes <iostream> headers; any
            # file that writes to std::cerr must include it explicitly.
            ("src/varray.h", [
                ("Add iostream include for std::cerr",
                 re.compile(
                     r'^(#include "compare\.h")$',
                     re.MULTILINE
                 ),
                 r'\1\n#include <iostream>'),
            ]),

            # --- src/field.h: add missing <stdexcept> for std::runtime_error ---
            ("src/field.h", [
                ("Add stdexcept header for std::runtime_error",
                 re.compile(
                     r'(#include\s+"cdo_vlist\.h"\s*\n)',
                     re.MULTILINE
                 ),
                 r'\1#include <stdexcept>\n'),
            ]),

            # --- src/mpmo_color.h: rewrite text_color() to avoid sstream ---
            # GCC 15.2.0 has a bug where basic_stringstream's constructor
            # triggers a "has no member named init" error. The root cause is
            # std::locale being incomplete when sstream templates are defined.
            # Fix: replace std::stringstream with std::string + std::to_string().
            ("src/mpmo_color.h", [
                ("Remove sstream include (replaced by string-based impl)",
                 '#include <sstream>\n',
                 '// #include <sstream>  // removed: GCC 15 sstream bug workaround\n'),

                ("Rewrite text_color() without stringstream",
                 re.compile(
                     r'static inline std::string\s*\n'
                     r'text_color\(TextColor foreground = NO_COLOR, TextMode mode = MODELESS, TextColor background = NO_COLOR\)\s*\n'
                     r'\{\s*\n'
                     r'.*?return s\.str\(\);\s*\n'
                     r'\}',
                     re.DOTALL
                 ),
                 # Use a lambda to avoid regex backreference interpretation of \0 in \033
                 lambda m: (
                     'static inline std::string\n'
                     'text_color(TextColor foreground = NO_COLOR, TextMode mode = MODELESS, TextColor background = NO_COLOR)\n'
                     '{\n'
                     '  if (!color_enabled()) return "";\n'
                     '\n'
                     '  std::string s = "\\033[";\n'
                     '  bool tty = true;\n'
                     '\n'
                     '  if (!tty)\n'
                     '    return s + "m";\n'
                     '\n'
                     '  if (!foreground && !background)\n'
                     '    s += "0";\n'
                     '\n'
                     '  if (foreground)\n'
                     '  {\n'
                     '    s += std::to_string(static_cast<int>(foreground));\n'
                     '    if (background) s += ";";\n'
                     '  }\n'
                     '  if (background)\n'
                     '  {\n'
                     '    s += std::to_string(10 + static_cast<int>(background));\n'
                     '    if (mode) s += ";";\n'
                     '  }\n'
                     '  else if (mode) { s += ";"; }\n'
                     '\n'
                     '  if (mode)\n'
                     '    s += std::to_string(static_cast<int>(mode));\n'
                     '\n'
                     '  s += "m";\n'
                     '  return s;\n'
                     '}'
                 )),
            ]),

            # --- GCC 15 sstream workaround: add <locale> before <sstream> ---
            # GCC 15.2.0 bug: std::stringstream needs std::locale complete, not just
            # forward-declared. Add #include <locale> before #include <sstream> in
            # files that use stringstream.

            ("src/util_string.h", [
                ("Add locale before sstream for GCC 15",
                 re.compile(
                     r'(#include <sstream>)',
                     re.MULTILINE
                 ),
                 r'#include <locale>\n\1'),
            ]),

            ("src/util_string.cc", [
                ("Add locale before sstream for GCC 15",
                 re.compile(
                     r'(#include <sstream>)',
                     re.MULTILINE
                 ),
                 r'#include <locale>\n\1'),
            ]),

            ("src/process.cc", [
                ("Add locale before sstream for GCC 15",
                 re.compile(
                     r'(#include <sstream>)',
                     re.MULTILINE
                 ),
                 r'#include <locale>\n\1'),
            ]),

            ("src/cdo_getopt.cc", [
                ("Add locale before sstream for GCC 15",
                 re.compile(
                     r'(#include <sstream>)',
                     re.MULTILINE
                 ),
                 r'#include <locale>\n\1'),
            ]),

            ("src/operators/Smooth.cc", [
                ("Add locale before sstream for GCC 15",
                 re.compile(
                     r'(#include <sstream>)',
                     re.MULTILINE
                 ),
                 r'#include <locale>\n\1'),
            ]),

            # --- Missing <sstream> includes for files using std::stringstream ---
            ("src/process_int.cc", [
                ("Add sstream include for std::stringstream",
                 re.compile(r'^(#include <cassert>)$', re.MULTILINE),
                 r'\1\n#include <sstream>'),
            ]),

            ("src/processManager.cc", [
                ("Add sstream include for std::stringstream",
                 re.compile(r'^(#include <mutex>)$', re.MULTILINE),
                 r'\1\n#include <sstream>'),
            ]),

            ("src/operators/Bitrounding.cc", [
                ("Add sstream include for std::stringstream",
                 re.compile(r'^(#include <cdi\.h>)$', re.MULTILINE),
                 r'\1\n#include <sstream>'),
            ]),

            ("src/operators/CMOR.cc", [
                ("Add sstream include for std::stringstream",
                 re.compile(r'^(#include <iomanip>)$', re.MULTILINE),
                 r'\1\n#include <sstream>'),
            ]),

            ("src/operators/Command.cc", [
                ("Add sstream include for std::stringstream",
                 re.compile(r'^(#include <iterator>)$', re.MULTILINE),
                 r'\1\n#include <sstream>'),
            ]),

            ("src/operators/Expr.cc", [
                ("Add sstream include for std::stringstream",
                 re.compile(r'^(#include <cassert>)$', re.MULTILINE),
                 r'\1\n#include <sstream>'),
            ]),

            ("src/operators/Getgridcell.cc", [
                ("Add sstream include for std::stringstream",
                 re.compile(r'^(#include <cdi\.h>)$', re.MULTILINE),
                 r'\1\n#include <sstream>'),
            ]),

            ("src/operators/Info.cc", [
                ("Add sstream include for std::stringstream",
                 re.compile(r'^(#include <numbers>)$', re.MULTILINE),
                 r'\1\n#include <sstream>'),
            ]),

            ("src/operators/Pack.cc", [
                ("Add sstream include for std::stringstream",
                 re.compile(r'^(#include <fstream>)$', re.MULTILINE),
                 r'\1\n#include <sstream>'),
            ]),

            ("src/operators/Remapgrid.cc", [
                ("Add sstream include for std::stringstream",
                 re.compile(r'^(#include <algorithm>)$', re.MULTILINE),
                 r'\1\n#include <sstream>'),
            ]),

            ("src/operators/Remapweights.cc", [
                ("Add sstream include for std::stringstream",
                 re.compile(r'^(#include <thread>)$', re.MULTILINE),
                 r'\1\n#include <sstream>'),
            ]),

            ("src/operators/Filter.cc", [
                # Move fftw3.h + mutex + fftwMutex OUTSIDE the HAVE_LIBFFTW3 guard.
                # This ensures fftw_complex/fftw_plan are always defined regardless
                # of whether configure detects HAVE_LIBFFTW3.  Our build always has
                # fftw3 available (--with-fftw3 is passed to configure), so the
                # unconditional include is safe.
                ("Unconditionally include fftw3.h, mutex and fftwMutex",
                 re.compile(
                     r'#ifdef HAVE_LIBFFTW3\s*\n#include <fftw3\.h>\s*\n(?:#include <mutex>\s*\n)?(?:static std::mutex fftwMutex;\s*\n)?#endif',
                     re.MULTILINE
                 ),
                 '#include <fftw3.h>\n#include <mutex>\nstatic std::mutex fftwMutex;'),
            ]),

            ("src/operators/Fourier.cc", [
                # Same treatment as Filter.cc: unconditionally include fftw3.h,
                # <mutex> and declare fftwMutex so their types/symbols are always
                # visible regardless of HAVE_LIBFFTW3.  The struct members and
                # function bodies that actually use FFTW are still guarded by
                # their own #ifdef HAVE_LIBFFTW3 blocks and will be compiled only
                # when the library was detected -- no runtime impact.
                ("Unconditionally include fftw3.h, mutex and fftwMutex",
                 re.compile(
                     r'#ifdef HAVE_LIBFFTW3\s*\n#include <fftw3\.h>\s*\n(?:#include <mutex>\s*\n)?(?:static std::mutex fftwMutex;\s*\n)?#endif',
                     re.MULTILINE
                 ),
                 '#include <fftw3.h>\n#include <mutex>\nstatic std::mutex fftwMutex;'),
            ]),

            # --- src/cdo_fctrans.cc: fftw3.h and mutex blocks unconditional ---
            # cdo_fctrans.cc has TWO separate HAVE_LIBFFTW3 guards:
            #   1) #ifdef HAVE_LIBFFTW3 / #include <fftw3.h> / #endif
            #   2) #ifdef HAVE_LIBFFTW3 / #include <mutex> / static...fftwMutex / #endif
            # Struct members (fftw_complex etc.) are inside a THIRD guard.
            # If configure fails to link fftw3 and HAVE_LIBFFTW3 stays undefined
            # via config.h we'd be safe (all code skipped), but if some external
            # HAVE_LIBFFTW3 define fires without the types being available we'd
            # get compile errors.  Unconditionally include to be bulletproof.
            ("src/cdo_fctrans.cc", [
                ("Unconditionally include fftw3.h",
                 "#ifdef HAVE_LIBFFTW3\n#include <fftw3.h>\n#endif",
                 "#include <fftw3.h>"),
                ("Unconditionally include mutex and fftwMutex",
                 "#ifdef HAVE_LIBFFTW3\n#include <mutex>\nstatic std::mutex fftwMutex;\n#endif",
                 "#include <mutex>\nstatic std::mutex fftwMutex;"),
            ]),

            # --- HDF5: unconditionally include hdf5.h ---
            # griddes_h5.cc and operators/Importcmsaf.cc wrap #include "hdf5.h"
            # in #ifdef HAVE_LIBHDF5, but their function signatures use hid_t /
            # herr_t types directly outside the guard.  Our build always links
            # HDF5 (--with-hdf5 is passed to configure), so the unconditional
            # include is safe and ensures the types are visible everywhere.
            ("src/griddes_h5.cc", [
                ("Unconditionally include hdf5.h",
                 "#ifdef HAVE_LIBHDF5\n#include \"hdf5.h\"\n#endif",
                 "#include \"hdf5.h\""),
            ]),

            ("src/operators/Importcmsaf.cc", [
                ("Unconditionally include hdf5.h",
                 "#ifdef HAVE_LIBHDF5\n#include \"hdf5.h\"\n#endif",
                 "#include \"hdf5.h\""),
            ]),

            # --- libcdi/src/table.h: fix guard clash and add tablepar.h ---
            # libcdi/src/table.h uses param_type without including tablepar.h.
            # It also uses the SAME include guard "TABLE_H" as CDO's src/table.h
            # which declares cdo::define_table.
            # When the libcdi version is processed first it sets TABLE_H, making
            # CDO's src/table.h a no-op -> cdo::define_table unavailable.
            # Fix: (1) rename guard to CDI_TABLE_H so both files are processed,
            #      (2) add #include "tablepar.h" so param_type/tableLink are visible.
            ("libcdi/src/table.h", [
                ("Rename include guard to avoid clash with src/table.h",
                 "/* Automatically generated, do not edit! */\n#ifndef TABLE_H\n#define TABLE_H",
                 "/* Automatically generated, do not edit! */\n#ifndef CDI_TABLE_H\n#define CDI_TABLE_H"),
                ("Add tablepar.h include to define param_type",
                 "#ifndef CDI_TABLE_H\n#define CDI_TABLE_H\n\n// clang-format off",
                 "#ifndef CDI_TABLE_H\n#define CDI_TABLE_H\n\n#include \"tablepar.h\"\n\n// clang-format off"),
            ]),

            # --- src/table.h: rename include guard to CDO_SRC_TABLE_H ---
            # Belt-and-suspenders: even if something outside our control sets
            # TABLE_H (e.g. a system-installed CDI header), CDO's own table.h
            # must still be processed so that cdo::define_table is declared.
            # Using a unique guard name ensures it is always included.
            ("src/table.h", [
                ("Rename include guard to CDO_SRC_TABLE_H",
                 "#ifndef TABLE_H\n#define TABLE_H\n\n#include <string>",
                 "#ifndef CDO_SRC_TABLE_H\n#define CDO_SRC_TABLE_H\n\n#include <string>"),
            ]),

            # --- src/operators/Setpartab.cc: fix table.h include resolution ---
            # The operator is compiled from src/operators/, and the include
            # search path has -I../libcdi/src listed before -idirafter .
            # So `#include "table.h"` resolves to libcdi/src/table.h (which
            # only contains data tables, not cdo::define_table) instead of
            # src/table.h.  Using the relative path "../table.h" unambiguously
            # resolves to src/table.h from the operators subdirectory.
            ("src/operators/Setpartab.cc", [
                ("Fix table.h include to use relative path to src/table.h",
                 '#include "table.h"',
                 '#include "../table.h"'),
            ]),

            # --- src/mpim_grid/grid_proj.cc: unconditionally include proj.h ---
            # PJ, proj_create, proj_destroy, proj_errno, proj_errno_string and
            # PJ_DEFAULT_CTX are used throughout the file but their declarations
            # are wrapped in #ifdef HAVE_LIBPROJ.  On Windows, PROJ is always
            # available as a build dependency, so include the header
            # unconditionally so the types are visible in every code path.
            ("src/mpim_grid/grid_proj.cc", [
                ("Unconditionally include proj.h",
                 "#ifdef HAVE_LIBPROJ\n#include \"proj.h\"\n#endif",
                 "#include \"proj.h\""),
            ]),

            # --- libcdi/src/util.c: provide cdiCreateUUID for Windows ---
            # The entire cdiCreateUUID implementation is wrapped in
            # `#ifndef _WIN32`, so on Windows the symbol is never defined,
            # causing an undefined-reference linker error.
            # The build already links -lrpcrt4, so we can use UuidCreate()
            # from <rpc.h> to generate a proper UUID v4.
            ("libcdi/src/util.c", [
                ("Add Windows cdiCreateUUID implementation using UuidCreate",
                 "#endif\n}\n#endif\n#endif\n\n/*\n * Local Variables:",
                 "#endif\n}\n#endif\n#endif\n\n"
                 "#ifdef _WIN32\n"
                 "#include <rpc.h>\n"
                 "void\n"
                 "cdiCreateUUID(unsigned char *uuid)\n"
                 "{\n"
                 "  UUID winUUID;\n"
                 "  UuidCreate(&winUUID);\n"
                 "  memcpy(uuid, &winUUID, CDI_UUID_SIZE);\n"
                 "}\n"
                 "#endif\n"
                 "\n/*\n * Local Variables:"),
            ]),

            # --- libcdi/configure: bypass POSIX.1-2001 check ---
            # MinGW does not define _POSIX_VERSION in <unistd.h>, but libcdi
            # is still buildable.  Force the check result to "yes".
            ("libcdi/configure", [
                ("Bypass POSIX.1-2001 conformance check",
                 "e) acx_cv_cc_posix_support2001=no ;;",
                 "e) acx_cv_cc_posix_support2001=yes ;;"),
            ]),

            # --- libcdi/src/input_file.c: implement pread for Windows ---
            # MinGW doesn't have pread(). The existing "#define pread read" is
            # incorrect (wrong argument count). Replace with proper implementation.
            ("libcdi/src/input_file.c", [
                ("Implement pread for Windows",
                 re.compile(
                     r'// On Windows, define ssize_t and pread manually\s*\n'
                     r'#ifdef _WIN32\s*\n'
                     r'#define ssize_t __int64\s*\n'
                     r'#define pread read\s*\n'
                     r'#include <io\.h>\s*\n'
                     r'#else\s*\n'
                     r'#include <unistd\.h>\s*\n'
                     r'#endif',
                     re.MULTILINE
                 ),
                 '''// On Windows, implement pread using _lseeki64 + _read
#ifdef _WIN32
#include <io.h>
typedef __int64 ssize_t;
typedef __int64 off_t;

static ssize_t pread_windows(int fd, void *buf, size_t count, off_t offset) {
    off_t current = _lseeki64(fd, 0, SEEK_CUR);
    if (current == -1) return -1;
    if (_lseeki64(fd, offset, SEEK_SET) == -1) return -1;
    int result = _read(fd, buf, (unsigned int)count);
    _lseeki64(fd, current, SEEK_SET);
    return result;
}
#define pread pread_windows
#else
#include <unistd.h>
#endif'''),
            ]),

            # --- libcdi/src/gribapi_utilities.c: implement setenv/unsetenv for Windows ---
            # MinGW doesn't have POSIX setenv()/unsetenv(). Implement using _putenv_s().
            ("libcdi/src/gribapi_utilities.c", [
                ("Implement setenv/unsetenv for Windows",
                 re.compile(
                     r'(#include <time\.h>)\s*\n'
                     r'\s*\n'
                     r'(#define FAIL_ON_GRIB_ERROR)',
                     re.MULTILINE
                 ),
                 r'''\1

// Windows compatibility: implement setenv/unsetenv
#ifdef _WIN32
#include <stdlib.h>
static int setenv(const char *name, const char *value, int overwrite) {
    if (!overwrite && getenv(name)) return 0;
    return _putenv_s(name, value);
}
static int unsetenv(const char *name) {
    return _putenv_s(name, "");
}
#endif

\2'''),
            ]),

            # --- src/cdo_settings.cc: disable CDI I/O threading on Windows ---
            # set_cdi_options() unconditionally calls cdiDefGlobal("THREADSAFE", 1)
            # whenever cdoLockIO == false (the default).  This enables CDI's
            # internal pthread_mutex protection for file I/O (CDI_IO_Mutex and
            # lazy-grid loadSerialize mutexes in cdf_lazy_grid.c).
            #
            # On Windows (MSYS2/MinGW64 with winpthreads), this causes an
            # ACCESS VIOLATION (0xC0000005) when reading NetCDF files created
            # by piped operators (e.g. chname, setname).  The winpthreads
            # pthread_mutex implementation has compatibility issues with CDI's
            # lazy-grid initialisation on Windows.
            #
            # CDO's *process* threading (pipe chains, parallel operators) is
            # a separate mechanism controlled by --with-threads=yes and does
            # NOT depend on CDI_Threadsafe.  Disabling CDI_Threadsafe on
            # Windows does not affect multi-operator pipelines.
            ("src/cdo_settings.cc", [
                ("Guard CDI THREADSAFE enable with #ifndef _WIN32",
                 "  if (Threading::cdoLockIO == false) cdiDefGlobal(\"THREADSAFE\", 1);",
                 "#ifndef _WIN32\n"
                 "  // CDI I/O threading is disabled on Windows: winpthreads\n"
                 "  // mutexes in CDI's lazy-grid code cause ACCESS VIOLATION\n"
                 "  // (0xC0000005) when reading files created by pipe operators.\n"
                 "  // CDO's own process threading (--with-threads=yes) is\n"
                 "  // unaffected by this setting.\n"
                 "  if (Threading::cdoLockIO == false) cdiDefGlobal(\"THREADSAFE\", 1);\n"
                 "#endif"),
            ]),

            # --- libcdi/src/cdf_lazy_grid.c: disable winpthreads mutex -------
            # CDI's lazy-grid code uses pthread_mutex_lock() / pthread_once()
            # to serialise deferred loading of coordinate variable arrays from
            # NC4 files (lat/lon/depth values stored as actual float variables).
            #
            # On Windows (MSYS2/MinGW64 + winpthreads), these mutex calls cause
            # an ACCESS VIOLATION (0xC0000005) when CDO opens any NetCDF file
            # that has explicit coordinate variables — e.g. real-world climate
            # datasets with lat/lon/depth arrays (the crash path NOT fixed by
            # the cdo_settings.cc THREADSAFE guard above).
            #
            # With --host cross-compile (CDO 2.5.3) HAVE_LIBPTHREAD was never
            # defined, so the branch was dead.  With native MSYS2 compilation
            # (CDO 2.6.0+) HAVE_LIBPTHREAD is defined and the mutex activates.
            #
            # Fix: undefine HAVE_LIBPTHREAD at file scope in cdf_lazy_grid.c
            # only.  CDO single-operator processes are inherently single-threaded
            # per file; the grid-load serialisation is not needed.  CDO's own
            # process-level parallelism (--with-threads=yes pipe chains) runs
            # at the operator level and is entirely unaffected.
            ("libcdi/src/cdf_lazy_grid.c", [
                ("Disable winpthreads mutex in lazy-grid (ACCESS VIOLATION fix)",
                 "#ifdef HAVE_CONFIG_H\n"
                 "#include \"config.h\"\n"
                 "#endif\n"
                 "\n"
                 "#ifdef HAVE_LIBNETCDF",
                 "#ifdef HAVE_CONFIG_H\n"
                 "#include \"config.h\"\n"
                 "#endif\n"
                 "\n"
                 "// Windows: winpthreads pthread_mutex_lock inside CDI's lazy-grid\n"
                 "// coordinate-loading code causes ACCESS VIOLATION (0xC0000005)\n"
                 "// when reading any NC4 file with explicit lat/lon/depth variables.\n"
                 "// Single-operator CDO processes don't require this thread-safety;\n"
                 "// disabling HAVE_LIBPTHREAD here makes this translation unit use\n"
                 "// the no-op mutex path without affecting CDO's process threading.\n"
                 "#ifdef _WIN32\n"
                 "#undef HAVE_LIBPTHREAD\n"
                 "#endif\n"
                 "\n"
                 "#ifdef HAVE_LIBNETCDF"),
            ]),

            # --- libcdi/src/resource_handle.c: disable winpthreads mutex ----
            # CDI's resource handle (reshLock/reshUnlock) calls
            # pthread_mutex_lock() and pthread_once() to protect the internal
            # resource list.  On Windows (MSYS2/MinGW64 + winpthreads) these
            # calls cause ACCESS VIOLATION (0xC0000005) in ptaxisCopy(), which
            # is called from cdf_read_timesteps() -> cdfInqContents() when CDO
            # opens any NC4 file that has a time dimension.
            #
            # Fix: same #undef HAVE_LIBPTHREAD pattern used for cdf_lazy_grid.c.
            # Single-operator CDO processes are effectively single-threaded per
            # open file; the reshLock serialisation is not needed on Windows.
            ("libcdi/src/resource_handle.c", [
                ("Disable winpthreads mutex in resource handle (ACCESS VIOLATION fix)",
                 "#ifdef HAVE_CONFIG_H\n"
                 "#include \"config.h\"\n"
                 "#endif\n"
                 "\n"
                 "#ifndef _XOPEN_SOURCE",
                 "#ifdef HAVE_CONFIG_H\n"
                 "#include \"config.h\"\n"
                 "#endif\n"
                 "\n"
                 "// Windows: winpthreads pthread_mutex_lock inside CDI's resource handle\n"
                 "// (reshLock/reshUnlock called from ptaxisCopy) causes ACCESS VIOLATION\n"
                 "// (0xC0000005) when reading NC4 files with a time dimension.\n"
                 "// Single-operator CDO processes don't require this thread-safety;\n"
                 "// disabling HAVE_LIBPTHREAD makes reshLock/reshUnlock become no-ops\n"
                 "// without affecting CDO's process-level parallelism.\n"
                 "#ifdef _WIN32\n"
                 "#undef HAVE_LIBPTHREAD\n"
                 "#endif\n"
                 "\n"
                 "#ifndef _XOPEN_SOURCE"),
            ]),

            # --- libcdi/src/stream_cdf_i.c: config.h must precede cdi.h -----
            # On Windows (LLP64 ABI), off_t is 'long' = 4 bytes UNLESS
            # _FILE_OFFSET_BITS=64 is defined before <sys/types.h> is first
            # included.  config.h (generated by autoconf) defines
            # _FILE_OFFSET_BITS=64, but stream_cdf_i.c was including cdi.h
            # (which transitively pulls in <sys/types.h>) BEFORE config.h.
            #
            # tsteps.c correctly includes config.h first (it even has a
            # comment warning about this exact pitfall).  Because of the
            # mismatch, tsteps_t::taxis sat at different offsets in the two
            # translation units (off_t=8 in tsteps.c vs off_t=4 in
            # stream_cdf_i.c), corrupting taxis->name with garbage and
            # causing ACCESS VIOLATION in ptaxisCopy().
            #
            # Fix: wrap config.h in the standard HAVE_CONFIG_H guard and move
            # it to the very top of stream_cdf_i.c, before any other includes.
            ("libcdi/src/stream_cdf_i.c", [
                ("Move config.h before cdi.h to fix off_t size (ACCESS VIOLATION fix)",
                 "#include \"cdi.h\"\n"
                 "#include \"cdi_limits.h\"\n"
                 "#ifdef HAVE_CONFIG_H\n"
                 "#include \"config.h\"\n"
                 "#endif",
                 "/* config.h MUST be first: on Windows (LLP64), off_t = long = 4 bytes unless\n"
                 " * _FILE_OFFSET_BITS=64 is defined (from config.h) before <sys/types.h>.\n"
                 " * cdi.h pulls in <sys/types.h> which fixes the size of off_t. Including\n"
                 " * config.h after cdi.h results in a struct layout mismatch for tsteps_t\n"
                 " * vs. tsteps.c (config.h first), causing ACCESS_VIOLATION in ptaxisCopy. */\n"
                 "#ifdef HAVE_CONFIG_H\n"
                 "#include \"config.h\"\n"
                 "#endif\n"
                 "#include \"cdi.h\"\n"
                 "#include \"cdi_limits.h\""),
            ]),

            # Keep unistd.h includes on MinGW: it is available and required by
            # multiple files for POSIX-like APIs (access, R_OK, etc.).
        ]

        # Apply all patches
        for rel_path, file_patches in patches:
            success, count = self.patch_file(rel_path, file_patches, dry_run)
            if success:
                total_files += 1
                total_patches += count
            print()

        # Summary
        mode = "Verify mode - no files modified" if dry_run else f"{total_files} files modified"
        print(f"{'='*70}")
        print(f"Completed: {mode}, {total_patches} patches applied")
        print(f"{'='*70}")

        return total_files > 0

    def restore_all(self) -> bool:
        """Restore all files from backups"""
        if not self.backup_dir.exists():
            print("No backups found - nothing to restore")
            return True

        print(f"Restoring from backup: {self.backup_dir}\n")

        # Load the mapping file to get correct original paths
        mappings = self._load_backup_mappings()

        restored = 0
        for backup_file in self.backup_dir.iterdir():
            if backup_file.is_file() and backup_file.name != "_path_map.txt":
                # Use mapping to recover original path; fall back to double-underscore split
                rel_path = mappings.get(backup_file.name)
                if not rel_path:
                    rel_path = backup_file.name.replace('__', '/')
                original_file = self.cdo_src / rel_path

                try:
                    content = backup_file.read_text(encoding='utf-8')
                    original_file.write_text(
                        content, encoding='utf-8', newline='\n')
                    print(f"[+] Restored: {rel_path}")
                    restored += 1
                except Exception as e:
                    print(f"[X] Restore failed {rel_path}: {e}")

        # Clean up backup directory
        try:
            for f in self.backup_dir.iterdir():
                f.unlink()
            self.backup_dir.rmdir()
            print(
                f"\n[OK] {restored} files restored, backup directory removed")
        except Exception as e:
            print(
                f"\n[!] {restored} files restored, but failed to remove backup dir: {e}")

        return True


def main():
    parser = argparse.ArgumentParser(description="CDO Windows Smart Patcher")
    parser.add_argument("action", choices=["apply", "restore", "verify"],
                        help="Action: apply/restore/verify")
    parser.add_argument("--cdo-src", type=Path,
                        help="CDO source directory (default: ../vendor/cdo)")

    args = parser.parse_args()

    # Determine CDO source directory
    if args.cdo_src:
        cdo_src = args.cdo_src
    else:
        script_dir = Path(__file__).parent
        cdo_src = script_dir.parent / "vendor" / "cdo"

    if not cdo_src.exists():
        print(f"[X] CDO source directory not found: {cdo_src}")
        print("   Please specify --cdo-src or ensure vendor/cdo exists")
        return 1

    patcher = WindowsPatcher(cdo_src)

    if args.action == "apply":
        success = patcher.apply_all(dry_run=False)
    elif args.action == "restore":
        success = patcher.restore_all()
    elif args.action == "verify":
        success = patcher.apply_all(dry_run=True)

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
