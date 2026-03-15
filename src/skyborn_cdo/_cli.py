"""
CLI entry point for skyborn-cdo.

Allows running CDO commands directly via: skyborn-cdo <cdo args>
or checking the installation via: skyborn-cdo --info
"""

import glob
import os
import sys

from skyborn_cdo._cdo_binary import get_bundled_env, get_cdo_path, get_cdo_version


def main():
    """Entry point for the skyborn-cdo console script."""
    args = sys.argv[1:]

    if not args or args[0] in ("--info", "-i"):
        _print_info()
        return

    # --help / -h with no further args → Python help; with args → passthrough to CDO
    if args[0] in ("--help", "-h"):
        if len(args) == 1:
            _print_help()
            return
        # e.g. "skyborn-cdo -h sellonlatbox" → forward to CDO

    # Convenience operator help style:
    #   skyborn-cdo mergetime --help
    #   skyborn-cdo mergetime --h
    # Rewrite to canonical CDO style: -h <operator>
    if len(args) >= 2 and not args[0].startswith("-") and args[1] in ("--help", "--h"):
        args = ["-h", args[0]]

    # PowerShell does not expand *.nc for native commands.  Mirror the Python
    # API behavior so `skyborn-cdo mergetime data_202*.nc out.nc` works on
    # Windows exactly like the shell-expanded form does on POSIX shells.
    if os.name == "nt" and any(c in " ".join(args) for c in ("*", "?", "[")):
        args = _expand_cli_wildcards(args)

    import subprocess

    cdo_path = get_cdo_path()
    env = get_bundled_env()

    # For help commands (-h operator), capture and print output
    # because CDO writes help to stderr and may return non-zero.
    if args[0] in ("-h",) and len(args) >= 2:
        _creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
        result = subprocess.run(
            [cdo_path] + args,
            capture_output=True,
            text=True,
            env=env,
            creationflags=_creationflags,
        )
        if result.stdout.strip():
            print(result.stdout.strip())
        if result.stderr.strip():
            print(result.stderr.strip())
        sys.exit(0)

    if os.name == "nt":
        # On Windows, CDO (MinGW build) hangs at process exit when any
        # NC4/HDF5 file has been processed (HDF5 DllMain cleanup deadlock).
        # subprocess.run() waits for the process to exit, so it hangs too.
        #
        # _run_windows() streams CDO's stdout/stderr to the terminal in
        # real-time (same as Linux) and kills CDO as soon as the output
        # file is fully written, exactly like using CDO on Linux.
        output_file = _extract_output_file(args)
        rc = _run_windows([cdo_path] + args, env, output_file)
        sys.exit(rc)
    else:
        result = subprocess.run([cdo_path] + args, env=env)
        sys.exit(result.returncode)


# ---------------------------------------------------------------------------
# Windows helpers
# ---------------------------------------------------------------------------

# File extensions that CDO uses for output files.
_CDO_FILE_EXTS = frozenset([
    ".nc", ".nc4", ".nc3", ".nc2", ".ncf",
    ".grib", ".grb", ".grb2",
    ".bin", ".srv", ".ext", ".ctl",
])


def _extract_output_file(args):
    """Return the likely CDO output file from *args*, or ``None``.

    CDO convention: the *last* positional argument is always the output
    file.  We only consider arguments that look like file paths (known
    climate-data extension, or contains a path separator).  If only one
    such argument is found it is an input-only call (``cdo info in.nc``)
    and we return ``None`` so the caller uses activity-based detection.
    """
    file_args = []
    for a in args:
        if a.startswith("-"):
            continue
        ext = os.path.splitext(a)[1].lower()
        if ext in _CDO_FILE_EXTS or os.sep in a or "/" in a or "\\" in a:
            file_args.append(a)
    return file_args[-1] if len(file_args) >= 2 else None


def _expand_cli_wildcards(args):
    """Expand wildcard file arguments for Windows shells.

    PowerShell does not expand globs for external commands, so CDO receives the
    literal pattern and fails with "Open failed on >*.nc<".  This mirrors the
    Python API's wildcard handling: expand only arguments that contain glob
    characters and leave unmatched patterns untouched so CDO can report them.
    """
    explicit_files = {
        a for a in args
        if not a.startswith("-") and not any(c in a for c in ("*", "?", "[", "]"))
    }
    expanded_args = []
    for arg in args:
        if any(c in arg for c in ("*", "?", "[", "]")):
            matches = sorted(m for m in glob.glob(
                arg) if m not in explicit_files)
            if matches:
                expanded_args.extend(matches)
            else:
                expanded_args.append(arg)
        else:
            expanded_args.append(arg)
    return expanded_args


def _run_windows(cmd, env, output_file=None):
    """Run CDO on Windows with real-time streaming and exit-hang workaround.

    CDO (MinGW build) hangs at process exit when NC4/HDF5 is involved.
    This function:

    * Pipes CDO stdout → sys.stdout and stderr → sys.stderr **in real
      time**, so the user sees CDO's progress messages as they are printed,
      exactly as on Linux.
    * Detects that CDO has finished its work by reading the HDF5 superblock
      (instant) or file-size stability (fallback for non-HDF5).  For
      info/show operators that write to stdout instead of a file, completion
      is detected by a short quiet period after output stops.
    * Kills the hung CDO process and returns exit code 0.

    For CDO *failures* (bad input, unknown operator, …) CDO exits cleanly
    with a non-zero code before any hang; those cases fall through to a
    normal ``proc.wait()`` and the real return code is preserved.
    """
    import subprocess
    import threading
    import time

    from skyborn_cdo._runner import _hdf5_file_is_closed

    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdin=subprocess.DEVNULL,
            env=env,
            creationflags=subprocess.CREATE_NO_WINDOW,
        )
    except FileNotFoundError:
        sys.stderr.write(f"skyborn-cdo: CDO binary not found: {cmd[0]}\n")
        return 1

    # Stream stdout/stderr to the terminal in real time.
    # We also track total bytes + last-activity timestamp for detecting
    # when info operators (no output file) have finished.
    _bytes_written = [0]
    _last_activity = [time.monotonic()]
    _lock = threading.Lock()

    def _stream(pipe, dest):
        try:
            fd = pipe.fileno()
            while True:
                chunk = os.read(fd, 4096)
                if not chunk:
                    break
                dest.buffer.write(chunk)
                dest.buffer.flush()
                with _lock:
                    _bytes_written[0] += len(chunk)
                    _last_activity[0] = time.monotonic()
        except (OSError, ValueError, AttributeError):
            pass
        finally:
            try:
                pipe.close()
            except OSError:
                pass

    t_out = threading.Thread(target=_stream, args=(
        proc.stdout, sys.stdout), daemon=True)
    t_err = threading.Thread(target=_stream, args=(
        proc.stderr, sys.stderr), daemon=True)
    t_out.start()
    t_err.start()

    # ------------------------------------------------------------------
    # Completion-detection state
    # ------------------------------------------------------------------
    # HDF5 two-phase detection:
    #   Phase A – wait for write-bit to be SET   (CDO opened the file)
    #   Phase B – wait for write-bit to be CLEARED (H5Fclose() returned)
    #
    # For a NEW output file:      initial_hdf5 = None  → write_started = True
    #   None → False (phase A done) → True (done!)
    #
    # For a PRE-EXISTING output file (overwrite with -O):
    #                             initial_hdf5 = True → write_started = False
    #   True (wait) → False (phase A done, write_started → True) → True (done!)
    #
    initial_hdf5 = _hdf5_file_is_closed(output_file) if output_file else None
    # True for new files; False for pre-existing
    write_started = initial_hdf5 is not True

    # Size-stability fallback for non-HDF5 output (GRIB, NC3 classic, …)
    _last_fsize = -1
    _last_fsize_t = 0.0
    _SIZE_STABLE = 1.0  # seconds

    # Quiet-period for info / show operators (output goes to stdout, no file)
    _QUIET_SECS = 0.5

    _POLL = 0.3       # proc.wait() timeout per iteration
    _GRACE = 0.1      # extra grace after detection before kill
    _elapsed = 0.0
    _detected_at = None

    while True:
        try:
            proc.wait(timeout=_POLL)
            # CDO exited cleanly (success or error) – drain and return.
            t_out.join(timeout=1)
            t_err.join(timeout=1)
            return proc.returncode
        except subprocess.TimeoutExpired:
            _elapsed += _POLL

        if _detected_at is None:
            done = False

            if output_file:
                closed = _hdf5_file_is_closed(output_file)
                if closed is not None:
                    # HDF5 file: two-phase write-bit detection.
                    if not write_started:
                        if closed is False:        # write bit set → phase A done
                            write_started = True
                    else:
                        if closed is True:         # write bit cleared → done!
                            done = True
                else:
                    # Not HDF5 (GRIB / NC3 classic): size-stability fallback.
                    try:
                        cur = (os.path.getsize(output_file)
                               if os.path.isfile(output_file) else -1)
                        if cur != _last_fsize:
                            _last_fsize = cur
                            _last_fsize_t = _elapsed
                        elif cur > 0 and (_elapsed - _last_fsize_t) >= _SIZE_STABLE:
                            done = True
                    except OSError:
                        pass
            else:
                # No output file (info / show operators).
                # CDO writes output, then hangs.  Detect via quiet period.
                with _lock:
                    total = _bytes_written[0]
                    since = time.monotonic() - _last_activity[0]
                if total > 0 and since >= _QUIET_SECS:
                    done = True

            if done:
                _detected_at = _elapsed

        if _detected_at is not None and (_elapsed - _detected_at) >= _GRACE:
            # CDO finished its work but is stuck in exit cleanup.  Kill it.
            try:
                proc.kill()
            except OSError:
                pass
            try:
                subprocess.Popen(
                    ["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
            except OSError:
                pass
            try:
                proc.wait(timeout=2)
            except subprocess.TimeoutExpired:
                pass
            t_out.join(timeout=1)
            t_err.join(timeout=1)
            return 0

    return 0  # unreachable, satisfies linters


def _print_info():
    """Print skyborn-cdo installation info."""
    import skyborn_cdo

    print(f"skyborn-cdo version: {skyborn_cdo.__version__}")
    print(f"CDO version target:  {skyborn_cdo.__cdo_version__}")
    try:
        cdo_path = get_cdo_path()
        print(f"CDO binary path:     {cdo_path}")
        print()
        version = get_cdo_version(cdo_path)
        print(version)
    except FileNotFoundError as e:
        print(f"\nCDO binary NOT FOUND: {e}")


def _print_help():
    print("skyborn-cdo: Pre-compiled CDO (Climate Data Operators) for Python")
    print()
    print("Usage:")
    print("  skyborn-cdo --info              Show CDO binary info and version")
    print("  skyborn-cdo --help              Show this help message")
    print("  skyborn-cdo -h <operator>       Show help for a specific CDO operator")
    print("  skyborn-cdo --operators         List all available CDO operators")
    print("  skyborn-cdo <cdo-args>          Pass arguments directly to CDO")
    print()
    print("Operator Help Examples:")
    print("  skyborn-cdo -h sellonlatbox     Show help for sellonlatbox")
    print("  skyborn-cdo -h mergetime        Show help for mergetime")
    print("  skyborn-cdo -h remapbil         Show help for remapbil")
    print()
    print("Command Examples:")
    print("  skyborn-cdo -O mergetime in1.nc in2.nc out.nc")
    print("  skyborn-cdo -O -f nc4 sellonlatbox,0,30,0,30 input.nc output.nc")
    print("  skyborn-cdo -O -f nc4 -fldmean -sellonlatbox,70,140,10,55 input.nc output.nc")
    print()
    print("Python API:")
    print("  from skyborn_cdo import Cdo")
    print('  cdo = Cdo(options="-O")')
    print('  cdo("mergetime in1.nc in2.nc out.nc")')
    print('  cdo.mergetime(input="in1.nc in2.nc", output="out.nc")')
    print('  print(cdo.help("sellonlatbox"))  # operator help in Python')


if __name__ == "__main__":
    main()
