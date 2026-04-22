"""
CLI entry point for skyborn-cdo.

Allows running CDO commands directly via: skyborn-cdo <cdo args>
or checking the installation via: skyborn-cdo --info
"""

import glob
import locale
import os
import sys

from skyborn_cdo._cdo_binary import get_bundled_env, get_cdo_path, get_cdo_version


def main():
    """Entry point for the skyborn-cdo console script."""
    args = _normalize_cli_args(sys.argv[1:])

    if not args or args[0] in ("--info", "-i"):
        _print_info()
        return

    # --help / -h with no further args: Python help; with args: passthrough to CDO
    if args[0] in ("--help", "-h"):
        if len(args) == 1:
            _print_help()
            return
        # e.g. "skyborn-cdo -h sellonlatbox": forward to CDO

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
            env=env,
            creationflags=_creationflags,
        )
        stdout = _decode_cli_output(result.stdout)
        stderr = _decode_cli_output(result.stderr)
        if stdout.strip():
            print(stdout.strip())
        if stderr.strip():
            print(stderr.strip())
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


def _normalize_cli_args(args):
    """Normalize convenience help forms to the canonical ``-h <operator>``."""
    if len(args) >= 2 and not args[0].startswith("-") and args[1] in ("--help", "--h", "-h"):
        return ["-h", args[0]]
    return args


def _decode_cli_output(data):
    """Decode subprocess output without assuming the console code page."""
    if data is None:
        return ""
    if isinstance(data, str):
        return data

    # MinGW-built CDO on Windows may emit UTF-8 help text even when the
    # active console encoding is GBK. Decode bytes ourselves so help output
    # never crashes inside subprocess' text-mode reader thread.
    encodings = []
    if os.name == "nt":
        encodings.append("utf-8")
    encodings.append(locale.getpreferredencoding(False) or "utf-8")
    encodings.extend(["utf-8", "latin-1"])

    seen = set()
    for enc in encodings:
        if not enc or enc in seen:
            continue
        seen.add(enc)
        try:
            return data.decode(enc)
        except UnicodeDecodeError:
            continue

    return data.decode("utf-8", errors="replace")


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

    * Pipes CDO stdout to sys.stdout and stderr to sys.stderr **in real
      time**, so the user sees CDO's progress messages as they are printed,
      exactly as on Linux.
    * Detects that CDO has finished its work by reading the HDF5 superblock
      (instant) or checking whether the output file can be opened
      exclusively after CDO releases its handle.  For
      info/show operators that write to stdout instead of a file, completion
      is detected by a short quiet period after output stops.
    * Kills the hung CDO process and returns exit code 0.

    For CDO *failures* (bad input, unknown operator, etc.) CDO exits cleanly
    with a non-zero code before any hang; those cases fall through to a
    normal ``proc.wait()`` and the real return code is preserved.
    """
    import subprocess
    import threading
    import time

    from skyborn_cdo._runner import (
        CdoRunner,
        _win_default_deadline,
        _win_init_completion_state,
        _win_timeout_completion_succeeded,
        _win_update_completion_state,
    )

    state = _win_init_completion_state(output_file, now=time.monotonic())

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
    _bytes_written = [0]
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

    # Quiet-period for info / show operators (output goes to stdout, no file)
    _QUIET_SECS = 0.5

    _POLL = 0.3       # proc.wait() timeout per iteration
    _GRACE = 0.1      # extra grace after detection before kill
    _elapsed = 0.0
    _detected_at = None
    _deadline = _win_default_deadline(output_file)
    env_timeout = os.environ.get("SKYBORN_CDO_TIMEOUT")
    if env_timeout:
        try:
            parsed_timeout = int(env_timeout)
            if parsed_timeout > 0:
                _deadline = parsed_timeout
            else:
                sys.stderr.write(
                    "skyborn-cdo: ignoring invalid SKYBORN_CDO_TIMEOUT value; "
                    "expected a positive integer number of seconds.\n"
                )
        except ValueError:
            sys.stderr.write(
                "skyborn-cdo: ignoring invalid SKYBORN_CDO_TIMEOUT value; "
                "expected a positive integer number of seconds.\n"
            )
    timed_out = False

    while True:
        try:
            proc.wait(timeout=_POLL)
            # CDO exited cleanly (success or error); drain and return.
            t_out.join(timeout=1)
            t_err.join(timeout=1)
            return proc.returncode
        except subprocess.TimeoutExpired:
            _elapsed += _POLL

        if _detected_at is None:
            with _lock:
                total = _bytes_written[0]
            done = _win_update_completion_state(
                state,
                stdout_size=total,
                now=time.monotonic(),
                quiet_secs=_QUIET_SECS,
            )
            if done:
                _detected_at = _elapsed

        if _detected_at is not None and (_elapsed - _detected_at) >= _GRACE:
            # CDO finished its work but is stuck in exit cleanup.  Kill it.
            CdoRunner._kill_proc_tree(proc)
            try:
                proc.wait(timeout=2)
            except subprocess.TimeoutExpired:
                pass
            t_out.join(timeout=1)
            t_err.join(timeout=1)
            return 0

        if _deadline and _elapsed >= _deadline:
            timed_out = True
            CdoRunner._kill_proc_tree(proc)
            try:
                proc.wait(timeout=2)
            except subprocess.TimeoutExpired:
                pass
            break

    if timed_out:
        t_out.join(timeout=1)
        t_err.join(timeout=1)
        completed = _win_timeout_completion_succeeded(state)
        if completed:
            return 0

        sys.stderr.write(
            f"skyborn-cdo: command timed out after {_elapsed:.0f}s before "
            f"the Windows completion safeguard could confirm the result.\n"
        )
        sys.stderr.write(
            "Set SKYBORN_CDO_TIMEOUT=<seconds> or use the Python API "
            "timeout=... override for longer jobs.\n"
        )
        return 1

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
    print("  cdo --info                      Windows-only alias for skyborn-cdo")
    print("  cdo-win --info                  Windows-only alias for skyborn-cdo")
    print("  skyborn-cdo --help              Show this help message")
    print("  skyborn-cdo -h <operator>       Show help for a specific CDO operator")
    print("  skyborn-cdo <operator> --help   Show help for a specific CDO operator")
    print("  skyborn-cdo --operators         List all available CDO operators")
    print("  skyborn-cdo <cdo-args>          Pass arguments directly to CDO")
    print()
    print("Operator Help Examples:")
    print("  skyborn-cdo -h sellonlatbox     Show help for sellonlatbox")
    print("  skyborn-cdo -h mergetime        Show help for mergetime")
    print("  skyborn-cdo mergetime --help    Alternate operator help syntax")
    print("  skyborn-cdo -h remapbil         Show help for remapbil")
    print()
    print("Command Examples:")
    print("  skyborn-cdo -O mergetime in1.nc in2.nc out.nc")
    print("  skyborn-cdo -O -f nc4 sellonlatbox,0,30,0,30 input.nc output.nc")
    print("  skyborn-cdo -O -f nc4 -fldmean -sellonlatbox,70,140,10,55 input.nc output.nc")
    print("  cdo -O sellonlatbox,70,140,15,55 input.nc output.nc")
    print("  cdo-win -O sellonlatbox,70,140,15,55 input.nc output.nc")
    print()
    print("Python API:")
    print("  from skyborn_cdo import Cdo")
    print('  cdo = Cdo(options="-O")')
    print('  cdo("mergetime in1.nc in2.nc out.nc")')
    print('  cdo.mergetime(input="in1.nc in2.nc", output="out.nc")')
    print('  print(cdo.help("sellonlatbox"))  # operator help in Python')
    print()
    print("Windows:")
    print("  Set SKYBORN_CDO_TIMEOUT=<seconds> to raise the default safeguard limit")


if __name__ == "__main__":
    main()
