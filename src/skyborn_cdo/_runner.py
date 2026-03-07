"""
Low-level subprocess runner for CDO commands.
"""

import os
import re
import shlex
import subprocess
import tempfile
import threading
import time
from pathlib import Path
from typing import List, Optional, Union

# Pattern CDO prints to stdout when processing completes successfully.
# NOTE: CDO 2.x writes this to stdout (not stderr), and on Windows the
# C-runtime stdio buffer is never flushed to the pipe before the exit-hang,
# so this regex is used only in the post-kill final check, not in the
# live polling loop where the HDF5 superblock method is used instead.
_CDO_DONE_RE = re.compile(r"Processed \d+ values? from \d+ variable")

# HDF5 file signature (first 8 bytes of every .nc4 / HDF5 file).
_HDF5_SIG = b"\x89HDF\r\n\x1a\n"


def _hdf5_file_is_closed(path: str):
    """Detect whether an HDF5/NC4 file has been properly closed by CDO.

    HDF5 stores a ``file_consistency_flags`` byte at offset 11 of the
    superblock (version 2+).  Bit 0 is set when the file is opened for
    writing and *cleared* by ``H5Fclose()``.  Reading that single byte
    gives an instant, reliable completion signal without any arbitrary
    stability wait.

    Returns
    -------
    True   – HDF5 file exists *and* has been properly closed (bit 0 == 0).
    False  – HDF5 file exists but is still open / too small to read yet.
    None   – file does not exist yet, or is not an HDF5 file (caller
             should fall back to size-stability detection).
    """
    try:
        if not os.path.isfile(path):
            return None
        with open(path, "rb") as fh:
            hdr = fh.read(12)
        if len(hdr) < 12 or hdr[:8] != _HDF5_SIG:
            return None   # not HDF5 → caller uses fallback
        if hdr[8] < 2:
            return None   # superblock v0/v1: no flags at byte 11
        return (hdr[11] & 0x01) == 0  # True = write-access bit cleared = closed
    except OSError:
        return False


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
    # Core execution – with Windows exit-hang workaround
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
        stdout_chunks: list = []
        stderr_chunks: list = []

        def _reader(pipe, buf):
            try:
                while True:
                    chunk = pipe.read(4096)
                    if not chunk:
                        break
                    buf.append(chunk)
            except (OSError, ValueError):
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
        #      completely written – instant detection, zero wait.
        #   2. stdout_chunks non-empty (info/read-only operators like
        #      showname) → output is the result, no file needed.
        #   3. Non-HDF5 output: file-size STABILITY fallback (rarely
        #      needed – CDO exits cleanly for GRIB/NC2 output).
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
        _SIZE_STABLE_SECS = 1.0  # fallback (non-HDF5 only): stable-size window
        _elapsed = 0.0
        _deadline = timeout if timeout else 0
        _detected_at = None    # timestamp when completion first seen
        _last_fsize = -1       # last observed output-file size
        _last_fsize_changed = 0.0  # time when file size last changed
        hung = False

        while True:
            try:
                proc.wait(timeout=_POLL)
                break  # process exited normally
            except subprocess.TimeoutExpired:
                _elapsed += _POLL

            # -- Quick-check: did CDO already finish its work? --------
            if _detected_at is None:
                _done = False
                if output_file:
                    # 1) NC4/HDF5: check superblock file_consistency_flags.
                    #    H5Fclose() clears bit 0 → instant "done" signal.
                    _closed = _hdf5_file_is_closed(output_file)
                    if _closed is True:
                        _done = True
                    elif _closed is None:
                        # Not HDF5 (GRIB/NC2): fall back to size-stability.
                        # CDO exits cleanly for these formats so this path
                        # is almost never reached on a healthy system.
                        try:
                            _cur_size = (
                                os.path.getsize(output_file)
                                if os.path.isfile(output_file) else -1
                            )
                            if _cur_size != _last_fsize:
                                _last_fsize = _cur_size
                                _last_fsize_changed = _elapsed
                            elif (_cur_size > 0
                                  and (_elapsed - _last_fsize_changed)
                                  >= _SIZE_STABLE_SECS):
                                _done = True
                        except OSError:
                            pass
                # 2) stdout data available (info / print operators)
                if not _done and stdout_chunks:
                    _done = True
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

            # Hard timeout — no completion detected
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
            #      stdout data) and killed the zombie → _detected_at
            #      is set → definitely completed.
            #   b) Hard timeout with no early detection → check
            #      stdout/stderr/output-file one final time.
            if _detected_at is not None:
                completed = True
            else:
                completed = (_CDO_DONE_RE.search(stdout)
                             or _CDO_DONE_RE.search(stderr))
                if not completed and output_file:
                    try:
                        completed = os.path.isfile(output_file) and \
                            os.path.getsize(output_file) > 0
                    except OSError:
                        completed = False
                if not completed and stdout.strip():
                    completed = True

            if completed:
                if self.debug:
                    print("[skyborn-cdo] Process hung at exit after successful "
                          f"completion – killed ({_elapsed:.1f}s).")
                return subprocess.CompletedProcess(cmd, 0, stdout, stderr)

            raise CdoError(
                f"CDO command timed out after {timeout}s: {label}",
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

        # Guess the output file: last non-option argument.
        _outf = None
        for _p in reversed(cmd[1:]):
            if not _p.startswith("-"):
                _outf = _p
                break

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
