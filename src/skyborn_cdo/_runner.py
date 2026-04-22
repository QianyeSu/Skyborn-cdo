"""
Low-level subprocess runner for CDO commands.
"""

import os
import shlex
import subprocess
import tempfile
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union

if os.name == "nt":
    import ctypes
    from ctypes import wintypes

# HDF5 file signature (first 8 bytes of every .nc4 / HDF5 file).
_HDF5_SIG = b"\x89HDF\r\n\x1a\n"
_WIN_OUTPUTFILE_DEADLINE = 3600
_WIN_NO_OUTFILE_DEADLINE = 30

if os.name == "nt":
    _GENERIC_READ = 0x80000000
    _OPEN_EXISTING = 3
    _ERROR_SHARING_VIOLATION = 32
    _ERROR_LOCK_VIOLATION = 33
    _INVALID_HANDLE_VALUE = wintypes.HANDLE(-1).value
    _kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    _CreateFileW = _kernel32.CreateFileW
    _CreateFileW.argtypes = [
        wintypes.LPCWSTR,
        wintypes.DWORD,
        wintypes.DWORD,
        wintypes.LPVOID,
        wintypes.DWORD,
        wintypes.DWORD,
        wintypes.HANDLE,
    ]
    _CreateFileW.restype = wintypes.HANDLE
    _CloseHandle = _kernel32.CloseHandle
    _CloseHandle.argtypes = [wintypes.HANDLE]
    _CloseHandle.restype = wintypes.BOOL


def _hdf5_file_is_closed(path: str):
    """Detect whether an HDF5/NC4 file has been properly closed by CDO.

    HDF5 stores a ``file_consistency_flags`` byte at offset 11 of the
    superblock (version 2+).  Bit 0 is set when the file is opened for
    writing and *cleared* by ``H5Fclose()``.  Reading that single byte
    gives an instant, reliable completion signal without any arbitrary
    stability wait.

    Returns
    -------
    True   HDF5 file exists and has been properly closed (bit 0 == 0).
    False  HDF5 superblock was read successfully and bit 0 == 1,
             meaning CDO still has the file open for writing.  The caller
             should trust this signal and NOT fall back to weaker completion
             probes while the file is still open for writing.
    None   Cannot determine: file does not exist, is not an HDF5 file,
             uses a superblock version < 2, or an OSError prevented reading
             the header (e.g. transient Windows sharing violation).  In all
             these cases the caller should use another completion probe.
    """
    try:
        if not os.path.isfile(path):
            return None
        with open(path, "rb") as fh:
            hdr = fh.read(12)
        if len(hdr) < 12 or hdr[:8] != _HDF5_SIG:
            return None   # not HDF5, caller uses fallback
        if hdr[8] < 2:
            return None   # superblock v0/v1: no flags at byte 11
        # True = write-access bit cleared = closed
        return (hdr[11] & 0x01) == 0
    except OSError:
        # Cannot read the superblock (e.g. transient Windows sharing
        # violation).  Return None so the caller can fall back to another
        # completion probe instead of spinning forever on a False sentinel.
        return None


def _win_file_is_exclusive_ready(path: str):
    """Return whether *path* can be opened exclusively on Windows.

    When CDO still holds the output file open, ``CreateFileW(..., share=0)``
    fails with a sharing violation.  Once CDO has closed the handle, the same
    probe succeeds even if the process itself is still hung in exit cleanup.
    """
    if os.name != "nt":
        return None

    try:
        if not os.path.isfile(path):
            return None
        handle = _CreateFileW(
            path,
            _GENERIC_READ,
            0,
            None,
            _OPEN_EXISTING,
            0,
            None,
        )
        if handle == _INVALID_HANDLE_VALUE:
            err = ctypes.get_last_error()
            if err in (_ERROR_SHARING_VIOLATION, _ERROR_LOCK_VIOLATION):
                return False
            return None
        _CloseHandle(handle)
        return True
    except OSError:
        return None


def _win_default_deadline(output_file: Optional[str]) -> int:
    """Return the default Windows safeguard timeout in seconds."""
    return _WIN_OUTPUTFILE_DEADLINE if output_file else _WIN_NO_OUTFILE_DEADLINE


@dataclass(frozen=True)
class _WinOutputFingerprint:
    """Minimal output-file identity used to detect real work on Windows."""

    existed: bool
    size: int = 0
    mtime_ns: int = 0


@dataclass
class _WinCompletionState:
    """Shared Windows completion state used by both runner and CLI."""

    output_file: Optional[str]
    initial_output: Optional[_WinOutputFingerprint]
    write_started: bool
    output_verified: bool
    quiet_period_observed: bool
    last_stdout_size: int
    last_stdout_t: float


def _win_capture_output_fingerprint(output_file: Optional[str]) -> Optional[_WinOutputFingerprint]:
    """Return a stable fingerprint for *output_file*, or ``None`` for stdout-only calls."""
    if not output_file:
        return None

    try:
        stat = os.stat(output_file)
    except OSError:
        return _WinOutputFingerprint(existed=False)

    mtime_ns = getattr(stat, "st_mtime_ns", int(stat.st_mtime * 1_000_000_000))
    return _WinOutputFingerprint(
        existed=True,
        size=stat.st_size,
        mtime_ns=mtime_ns,
    )


def _win_init_completion_state(
    output_file: Optional[str],
    *,
    now: Optional[float] = None,
) -> _WinCompletionState:
    """Capture the pre-launch Windows completion baseline."""
    if now is None:
        now = time.monotonic()

    return _WinCompletionState(
        output_file=output_file,
        initial_output=_win_capture_output_fingerprint(output_file),
        write_started=False,
        output_verified=False,
        quiet_period_observed=False,
        last_stdout_size=0,
        last_stdout_t=now,
    )


def _win_output_fingerprint_changed(
    initial: Optional[_WinOutputFingerprint],
    current: Optional[_WinOutputFingerprint],
) -> bool:
    """Return whether a pre-existing output file changed identity during this run."""
    return bool(
        initial
        and current
        and initial.existed
        and current.existed
        and (current.size != initial.size or current.mtime_ns != initial.mtime_ns)
    )


def _win_output_created_nonempty(
    initial: Optional[_WinOutputFingerprint],
    current: Optional[_WinOutputFingerprint],
) -> bool:
    """Return whether this run produced a new non-empty output file."""
    return bool(
        initial
        and current
        and not initial.existed
        and current.existed
        and current.size > 0
    )


def _win_update_output_completion_state(
    state: _WinCompletionState,
    *,
    process_alive: bool = True,
) -> bool:
    """Update output-file completion state and return whether completion is verified."""
    if not state.output_file:
        return False

    current = _win_capture_output_fingerprint(state.output_file)
    hdf5_closed = _hdf5_file_is_closed(state.output_file)
    exclusive_ready = _win_file_is_exclusive_ready(state.output_file)
    fingerprint_changed = _win_output_fingerprint_changed(state.initial_output, current)
    created_nonempty = _win_output_created_nonempty(state.initial_output, current)

    if not state.write_started and (
        hdf5_closed is False
        or (process_alive and exclusive_ready is False)
        or fingerprint_changed
        or created_nonempty
    ):
        state.write_started = True

    # A readable HDF5 "still open for write" bit is a stronger signal than
    # the Windows handle probe. Do not accept exclusive-open success while
    # the superblock still says the file was not cleanly closed.
    if hdf5_closed is False:
        file_ready = False
    elif hdf5_closed is True:
        file_ready = True
    elif process_alive and exclusive_ready is True:
        file_ready = True
    else:
        file_ready = False

    if state.write_started and file_ready and (fingerprint_changed or created_nonempty):
        state.output_verified = True

    return state.output_verified


def _win_update_stdout_completion_state(
    state: _WinCompletionState,
    *,
    stdout_size: int,
    now: float,
    quiet_secs: float,
) -> bool:
    """Update stdout-only completion state and return whether quiet completion was observed."""
    if stdout_size != state.last_stdout_size:
        state.last_stdout_size = stdout_size
        state.last_stdout_t = now
        return False

    if stdout_size > 0 and (now - state.last_stdout_t) >= quiet_secs:
        state.quiet_period_observed = True

    return state.quiet_period_observed


def _win_update_completion_state(
    state: _WinCompletionState,
    *,
    stdout_size: int,
    now: float,
    quiet_secs: float,
) -> bool:
    """Shared Windows completion update for output-file and stdout-only commands."""
    if state.output_file:
        return _win_update_output_completion_state(state, process_alive=True)

    return _win_update_stdout_completion_state(
        state,
        stdout_size=stdout_size,
        now=now,
        quiet_secs=quiet_secs,
    )


def _win_timeout_completion_succeeded(state: _WinCompletionState) -> bool:
    """Return whether a hard-timeout can still be treated as success."""
    if state.output_file:
        if state.output_verified:
            return True
        return _win_update_output_completion_state(state, process_alive=False)

    return state.quiet_period_observed


class CdoError(Exception):
    """Exception raised when a CDO command fails."""

    def __init__(self, message: str, returncode: int = -1, stderr: str = "", cmd: str = ""):
        self.returncode = returncode
        self.stderr = stderr
        self.cmd = cmd
        super().__init__(message)


class CdoRunner:
    """
    Low-level CDO command executor.

    Wraps ``subprocess.run`` to invoke the CDO binary with proper
    environment setup and error handling.
    """

    def __init__(self, cdo_path: str, env: Optional[dict] = None, debug: bool = False):
        """
        Parameters
        ----------
        cdo_path : str
            Absolute path to the CDO executable.
        env : dict, optional
            Environment variables for the CDO process.
        debug : bool
            If True, print commands and stderr to stdout.
        """
        self.cdo_path = cdo_path
        self.env = env or os.environ.copy()
        self.debug = debug

    # -----------------------------------------------------------------
    # Process management helpers
    # -----------------------------------------------------------------

    @staticmethod
    def _kill_proc_tree(proc: subprocess.Popen) -> None:
        """Kill a process and all its children.

        On Windows ``proc.kill()`` only terminates the main process;
        child processes (e.g. HDF5/NetCDF helpers) can survive and hold
        file locks or block pipes.  ``taskkill /F /T`` kills the whole
        process tree but can be slow (~5 s).  We do proc.kill() first
        for speed, then background taskkill as cleanup.
        """
        if os.name == "nt":
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
        else:
            try:
                proc.kill()
            except OSError:
                pass

    # -----------------------------------------------------------------
    # Core execution with Windows exit-hang workaround
    # -----------------------------------------------------------------

    def _exec(self, cmd: List[str], timeout: Optional[int],
              cmd_label: Optional[str] = None,
              output_file: Optional[str] = None) -> subprocess.CompletedProcess:
        """Run *cmd* and return a CompletedProcess.

        On all platforms the CDO binary is executed with PIPE on
        stdout / stderr.  On Windows certain builds of CDO hang during
        process exit (after data processing has already completed) when
        the output format is set via ``-f nc*``.  Because the process
        never exits, the pipe EOF is never reached and
        ``communicate()`` blocks forever.

        To work around this, on Windows we read stdout/stderr in daemon
        threads (so reads are non-blocking relative to ``wait()``) and
        use ``proc.wait()`` instead of ``proc.communicate()``.  Completion
        is detected by inspecting the HDF5/NC4 output file's superblock:
        ``H5Fclose()`` clears bit 0 of ``file_consistency_flags`` (offset 11
        of a version-2+ superblock), giving an instant signal with no
        arbitrary wait.  When the work is done we kill the hung process and
        return a successful result.

        Parameters
        ----------
        output_file : str, optional
            Path to the expected output file.  Used on Windows to detect
            that CDO completed its work even when the process hangs at
            exit (stderr may be empty due to C-runtime buffering).
        """
        label = cmd_label or " ".join(cmd)
        _creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0

        if os.name != "nt":
            # --- POSIX: straightforward communicate() ----------------
            try:
                proc = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    stdin=subprocess.DEVNULL,
                    text=True,
                    env=self.env,
                )
                stdout, stderr = proc.communicate(timeout=timeout)
                return subprocess.CompletedProcess(
                    cmd, proc.returncode, stdout, stderr)
            except subprocess.TimeoutExpired:
                self._kill_proc_tree(proc)
                try:
                    proc.communicate(timeout=5)
                except (subprocess.TimeoutExpired, OSError):
                    pass
                raise CdoError(
                    f"CDO command timed out after {timeout}s: {label}",
                    returncode=-1, stderr="", cmd=label,
                )
            except FileNotFoundError:
                raise CdoError(
                    f"CDO binary not found at: {self.cdo_path}",
                    returncode=-1, cmd=label,
                )

        # --- Windows: threaded-pipe strategy --------------------------
        #
        # We still use PIPE (not file redirect) for stdout/stderr and
        # drain them via daemon threads so that ``proc.wait()`` can time
        # out independently of the pipe EOF.
        _STDOUT_STABLE_SECS = 0.8  # stdout-only operators: growth window
        _state = _win_init_completion_state(output_file, now=time.monotonic())

        try:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                stdin=subprocess.DEVNULL,
                env=self.env,
                creationflags=_creationflags,
            )
        except FileNotFoundError:
            raise CdoError(
                f"CDO binary not found at: {self.cdo_path}",
                returncode=-1, cmd=label,
            )

        # Read pipes in daemon threads so proc.wait() can detect the
        # timeout independently of whether the pipes have reached EOF.
        #
        # IMPORTANT: use os.read(fd, n) rather than pipe.read(n).
        # BufferedReader.read(4096) on a pipe can wait for the full 4096-byte
        # request or EOF, which means small outputs such as ``sinfo`` may never
        # appear in stdout_chunks if CDO hangs during exit before EOF arrives.
        # os.read() returns as soon as any bytes are available.
        stdout_chunks: list = []
        stderr_chunks: list = []

        def _reader(pipe, buf):
            try:
                fd = pipe.fileno()
                while True:
                    chunk = os.read(fd, 4096)
                    if not chunk:
                        break
                    buf.append(chunk)
            except (OSError, ValueError):
                pass
            finally:
                try:
                    pipe.close()
                except OSError:
                    pass

        tout = threading.Thread(target=_reader,
                                args=(proc.stdout, stdout_chunks),
                                daemon=True)
        terr = threading.Thread(target=_reader,
                                args=(proc.stderr, stderr_chunks),
                                daemon=True)
        tout.start()
        terr.start()

        # --- Polling loop -------------------------------------------
        #
        # Instead of blocking on proc.wait(timeout) for the full
        # duration, we poll every _POLL seconds.  Most CDO operations
        # finish in < 0.1 s; the exit-hang means the *process* stays
        # alive but the work is done.
        #
        # Completion detection strategy (in priority order):
        #   1. HDF5/NC4 output: read the superblock file_consistency_flags
        #      byte (offset 11, version 2+).  HDF5 sets bit 0 on open and
        #      clears it on H5Fclose().  When the bit is 0 the file is
        #      completely written: instant detection, zero wait.
        #   2. Windows file-handle probe: once CDO closes the output file,
        #      it becomes exclusively openable even if the process is still
        #      stuck in exit cleanup. This avoids treating a stable 30 KB
        #      NetCDF4 header as a finished multi-GB output file.
        #   3. stdout STABILITY fallback for stdout-only operators like
        #      info/sinfo/showname.  Once stdout stops growing for a short
        #      window, the useful payload has arrived and the hung process can
        #      be killed.
        #
        # Why NOT rely on "Processed N values" stderr message:
        #   CDO 2.x prints this to stdout (not stderr).  On Windows the
        #   MinGW C-runtime stdio buffer is never flushed to the pipe
        #   because the exit-hang (inside HDF5 DllMain cleanup) occurs
        #   before stdio streams are flushed by exit().  The message
        #   therefore never arrives on the pipe, regardless of buffer size.
        #
        _POLL = 0.3            # seconds between polls
        _GRACE_AFTER = 0.1     # short grace after close detected, then kill
        _elapsed = 0.0
        _deadline = timeout if timeout else 0

        # On Windows, the wrapper needs a safeguard deadline because the
        # completion probe may never resolve if the process hangs in exit
        # cleanup or another process interferes with the output file handle.
        # (e.g. inside HDF5 H5Fclose, before atexit flushes stdio), the
        # buffer is never written to the pipe and stdout_chunks stays empty
        # forever; the loop would spin indefinitely without a deadline.
        # Apply a conservative fallback so the caller gets a CdoError
        # instead of an infinite hang.  Users who need a longer window can
        # pass timeout= explicitly (which sets _deadline above).
        if os.name == "nt" and _deadline == 0:
            _deadline = _win_default_deadline(output_file)
        _detected_at = None    # wall-clock time when completion first seen
        hung = False

        while True:
            try:
                proc.wait(timeout=_POLL)
                break  # process exited normally
            except subprocess.TimeoutExpired:
                _elapsed += _POLL

            _now = time.monotonic()

            # -- Quick-check: did CDO already finish its work? --------
            if _detected_at is None:
                _done = _win_update_completion_state(
                    _state,
                    stdout_size=sum(len(chunk) for chunk in stdout_chunks),
                    now=_now,
                    quiet_secs=_STDOUT_STABLE_SECS,
                )
                if _done:
                    _detected_at = _elapsed

            # If completion was detected, allow a short grace period
            # for the process to exit cleanly, then kill it.
            if _detected_at is not None:
                if (_elapsed - _detected_at) >= _GRACE_AFTER:
                    hung = True
                    self._kill_proc_tree(proc)
                    try:
                        proc.wait(timeout=2)
                    except subprocess.TimeoutExpired:
                        pass
                    break

            # Hard timeout: no completion detected
            if _deadline and _elapsed >= _deadline:
                hung = True
                self._kill_proc_tree(proc)
                try:
                    proc.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    pass
                break

        # After the process is dead the pipes will reach EOF and the
        # reader threads will finish.
        tout.join(timeout=2)
        terr.join(timeout=2)

        stdout = b"".join(stdout_chunks).decode("utf-8", errors="replace")
        stderr = b"".join(stderr_chunks).decode("utf-8", errors="replace")

        if hung:
            # CDO may have finished processing but hangs during exit
            # cleanup (known issue with certain MinGW/HDF5 builds).
            #
            # Two scenarios reach here:
            #   a) The polling loop detected completion (output file or
            #      stdout data) and killed the zombie; _detected_at
            #      is set, so completion was definitely observed.
            #   b) Hard timeout with no early detection; check
            #      stdout/stderr/output-file one final time.
            if _detected_at is not None:
                completed = True
            else:
                completed = _win_timeout_completion_succeeded(_state)

            if completed:
                if self.debug:
                    print("[skyborn-cdo] Process hung at exit after successful "
                          f"completion, killed ({_elapsed:.1f}s).")
                return subprocess.CompletedProcess(cmd, 0, stdout, stderr)

            raise CdoError(
                f"CDO command timed out after {_elapsed:.0f}s: {label}\n"
                f"  Hint: CDO may have hung inside HDF5/NetCDF4 file operations "
                f"before the Windows completion safeguard could confirm the "
                f"result. If the operation legitimately takes longer, pass "
                f"timeout=<seconds> to override the default limit.",
                returncode=-1, stderr=stderr, cmd=label,
            )

        return subprocess.CompletedProcess(
            cmd, proc.returncode, stdout, stderr)

    # -----------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------

    def run(
        self,
        args: List[str],
        input_files: Optional[List[str]] = None,
        output_file: Optional[str] = None,
        options: Optional[List[str]] = None,
        return_output: bool = False,
        timeout: Optional[int] = None,
    ) -> Union[str, int]:
        """
        Execute a CDO command.

        Parameters
        ----------
        args : list of str
            CDO operator and its arguments, e.g. ["-mergetime"].
        input_files : list of str, optional
            Input file paths.
        output_file : str, optional
            Output file path.
        options : list of str, optional
            Global CDO options like ["-O", "-s"].
        return_output : bool
            If True, return stdout content instead of returncode.
        timeout : int, optional
            Timeout in seconds.

        Returns
        -------
        str or int
            stdout content if ``return_output`` is True, else returncode.

        Raises
        ------
        CdoError
            If CDO exits with a non-zero return code.
        """
        cmd = [self.cdo_path]

        if options:
            cmd.extend(options)

        cmd.extend(args)

        if input_files:
            cmd.extend(input_files)

        if output_file:
            cmd.append(output_file)

        if self.debug:
            print(f"[skyborn-cdo] Running: {' '.join(cmd)}")

        result = self._exec(cmd, timeout, output_file=output_file)

        if self.debug and result.stderr:
            print(f"[skyborn-cdo] stderr: {result.stderr}")

        if result.returncode != 0:
            raise CdoError(
                f"CDO command failed (exit code {result.returncode}):\n"
                f"  Command: {' '.join(cmd)}\n"
                f"  Error: {result.stderr.strip()}",
                returncode=result.returncode,
                stderr=result.stderr,
                cmd=" ".join(cmd),
            )

        if return_output:
            return result.stdout

        return result.returncode

    def run_raw(self, cmd_string: str, timeout: Optional[int] = None) -> subprocess.CompletedProcess:
        """
        Execute a raw CDO command string.

        The command string is parsed with ``shlex.split`` and the first
        token (``cdo``) is replaced with the actual CDO binary path.

        Parameters
        ----------
        cmd_string : str
            Full CDO command string, e.g. "cdo mergetime in1.nc in2.nc out.nc"
            or "cdo -O -s mergetime in*.nc out.nc". The leading ``cdo`` is optional.

        Returns
        -------
        subprocess.CompletedProcess
            The completed process object.

        Raises
        ------
        CdoError
            If CDO exits with non-zero return code.
        """
        if os.name == 'nt':
            parts = shlex.split(cmd_string, posix=False)
            parts = [p.strip('"').strip("'") for p in parts]
        else:
            parts = shlex.split(cmd_string)

        # Strip leading 'cdo' if present
        if parts and parts[0].lower() in ("cdo", "cdo.exe"):
            parts = parts[1:]

        cmd = [self.cdo_path] + parts

        if self.debug:
            print(f"[skyborn-cdo] Running: {' '.join(cmd)}")

        # Guess the output file: last argument that looks like a file path
        # (contains a dot or a path separator) among ALL such arguments.
        # We require at least 2 file-like args so that a command whose only
        # file argument is an INPUT file (e.g. "sinfo data.nc" or
        # "showname data.nc") never has that input file mistaken for an
        # output file.  Misidentifying an existing input HDF5 file as the
        # output file would cause _hdf5_file_is_closed() to return True
        # immediately, killing CDO before it does any work.
        _file_like_args = [
            p for p in cmd[1:]
            if not p.startswith("-") and ("." in p or os.sep in p or "/" in p)
        ]
        _outf = _file_like_args[-1] if len(_file_like_args) >= 2 else None

        result = self._exec(cmd, timeout, cmd_label=cmd_string,
                            output_file=_outf)

        if result.returncode != 0:
            raise CdoError(
                f"CDO command failed (exit code {result.returncode}):\n"
                f"  Command: {cmd_string}\n"
                f"  Error: {result.stderr.strip()}",
                returncode=result.returncode,
                stderr=result.stderr,
                cmd=cmd_string,
            )

        return result
