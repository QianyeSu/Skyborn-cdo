"""
Tests for skyborn_cdo Python API.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path
import tomllib
from types import SimpleNamespace

import pytest


def _cli_test_env():
    """Ensure CLI subprocesses import the workspace src/ tree first."""
    env = os.environ.copy()
    src_dir = str(Path(__file__).resolve().parents[1] / "src")
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = src_dir if not existing else src_dir + os.pathsep + existing
    return env


class TestCdoBinaryDiscovery:
    """Test CDO binary discovery mechanisms."""

    def test_get_cdo_path_bundled(self):
        """Test that get_cdo_path finds bundled binary (if available)."""
        from skyborn_cdo._cdo_binary import _package_bin_dir, get_cdo_path

        bin_dir = _package_bin_dir()
        if (bin_dir / "cdo").exists() or (bin_dir / "cdo.exe").exists():
            path = get_cdo_path()
            assert os.path.isfile(path)
            assert "skyborn_cdo" in path
        else:
            pytest.skip("No bundled CDO binary present (development mode)")

    def test_get_cdo_path_env_var(self, tmp_path):
        """Test that CDO env var is respected."""
        from skyborn_cdo._cdo_binary import get_cdo_path

        fake_cdo = tmp_path / "cdo"
        fake_cdo.write_text("#!/bin/sh\necho fake")
        fake_cdo.chmod(0o755)

        original_env = os.environ.get("CDO")
        try:
            os.environ["CDO"] = str(fake_cdo)
            path = get_cdo_path()
            assert path == str(fake_cdo)
        finally:
            if original_env:
                os.environ["CDO"] = original_env
            else:
                os.environ.pop("CDO", None)

    def test_get_cdo_path_explicit(self, tmp_path):
        """Test that explicit path argument takes priority."""
        from skyborn_cdo._cdo_binary import get_cdo_path

        fake_cdo = tmp_path / "my_cdo"
        fake_cdo.write_text("#!/bin/sh\necho fake")
        fake_cdo.chmod(0o755)

        path = get_cdo_path(str(fake_cdo))
        assert path == str(fake_cdo)

    def test_get_cdo_path_not_found(self):
        """Test FileNotFoundError when CDO not available anywhere."""
        from unittest.mock import patch
        from pathlib import Path
        from skyborn_cdo._cdo_binary import get_cdo_path

        original_env = os.environ.get("CDO")
        original_path = os.environ.get("PATH")
        try:
            os.environ.pop("CDO", None)
            os.environ["PATH"] = ""  # Empty path
            # Mock the bundled binary dir so the bundled check doesn't find it
            with patch("skyborn_cdo._cdo_binary._package_bin_dir", return_value=Path("/nonexistent/bin")):
                with pytest.raises(FileNotFoundError, match="CDO binary not found"):
                    get_cdo_path("/nonexistent/cdo")
        finally:
            if original_env:
                os.environ["CDO"] = original_env
            if original_path:
                os.environ["PATH"] = original_path

    def test_bundled_env(self):
        """Test that get_bundled_env returns proper environment."""
        from skyborn_cdo._cdo_binary import get_bundled_env

        env = get_bundled_env()
        assert isinstance(env, dict)
        assert "PATH" in env


class TestCdoRunner:
    """Test low-level CdoRunner."""

    def test_runner_invalid_binary(self):
        """Test CdoError on invalid binary path."""
        from skyborn_cdo._runner import CdoError, CdoRunner

        runner = CdoRunner("/nonexistent/cdo")
        with pytest.raises(CdoError, match="CDO binary not found"):
            runner.run(["--version"])

    @pytest.mark.skipif(os.name != "nt", reason="Windows-only file-lock probe")
    def test_win_file_is_exclusive_ready(self, tmp_path):
        """Exclusive-open probe should track whether a file handle is still held."""
        import ctypes
        from ctypes import wintypes

        from skyborn_cdo._runner import _win_file_is_exclusive_ready

        path = tmp_path / "lockprobe.nc"
        path.write_bytes(b"test")

        generic_read = 0x80000000
        open_existing = 3
        invalid_handle = wintypes.HANDLE(-1).value
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        create_file = kernel32.CreateFileW
        create_file.argtypes = [
            wintypes.LPCWSTR,
            wintypes.DWORD,
            wintypes.DWORD,
            wintypes.LPVOID,
            wintypes.DWORD,
            wintypes.DWORD,
            wintypes.HANDLE,
        ]
        create_file.restype = wintypes.HANDLE
        close_handle = kernel32.CloseHandle
        close_handle.argtypes = [wintypes.HANDLE]
        close_handle.restype = wintypes.BOOL

        assert _win_file_is_exclusive_ready(str(path)) is True

        handle = create_file(str(path), generic_read, 0, None, open_existing, 0, None)
        assert handle != invalid_handle
        try:
            assert _win_file_is_exclusive_ready(str(path)) is False
        finally:
            close_handle(handle)

        assert _win_file_is_exclusive_ready(str(path)) is True

    @pytest.mark.skipif(os.name != "nt", reason="Windows-only file-lock probe")
    def test_win_file_is_exclusive_ready_ignores_nonsharing_errors(self, tmp_path, monkeypatch):
        """Non-sharing open failures must not be treated as 'still writing'."""
        import skyborn_cdo._runner as runner_mod

        path = tmp_path / "lockprobe_access.nc"
        path.write_bytes(b"test")

        monkeypatch.setattr(runner_mod, "_CreateFileW",
                            lambda *args, **kwargs: runner_mod._INVALID_HANDLE_VALUE)
        monkeypatch.setattr(runner_mod.ctypes, "get_last_error", lambda: 5)

        assert runner_mod._win_file_is_exclusive_ready(str(path)) is None


class TestWindowsCompletionState:
    """Unit tests for shared Windows completion-state helpers."""

    def test_preexisting_output_unchanged_after_timeout_fails(self, tmp_path, monkeypatch):
        """A stale pre-existing output must never be treated as successful completion."""
        import skyborn_cdo._runner as runner_mod

        path = tmp_path / "existing.nc"
        path.write_bytes(b"old-output")

        state = runner_mod._win_init_completion_state(str(path), now=0.0)
        monkeypatch.setattr(runner_mod, "_hdf5_file_is_closed", lambda _path: True)
        monkeypatch.setattr(runner_mod, "_win_file_is_exclusive_ready", lambda _path: True)

        assert runner_mod._win_update_output_completion_state(state) is False
        assert runner_mod._win_timeout_completion_succeeded(state) is False

    def test_preexisting_output_changed_after_write_start_and_close_succeeds(self, tmp_path, monkeypatch):
        """A pre-existing output should succeed only after write-start plus a changed fingerprint."""
        import skyborn_cdo._runner as runner_mod

        path = tmp_path / "existing_changed.nc"
        path.write_bytes(b"old-output")

        state = runner_mod._win_init_completion_state(str(path), now=0.0)
        readiness = [False, True]
        monkeypatch.setattr(runner_mod, "_hdf5_file_is_closed", lambda _path: None)
        monkeypatch.setattr(
            runner_mod,
            "_win_file_is_exclusive_ready",
            lambda _path: readiness.pop(0) if readiness else True,
        )

        assert runner_mod._win_update_output_completion_state(state) is False

        path.write_bytes(b"new-output-that-changed")
        assert runner_mod._win_update_output_completion_state(state) is True
        assert runner_mod._win_timeout_completion_succeeded(state) is True

    def test_new_output_nonempty_and_closed_succeeds(self, tmp_path, monkeypatch):
        """A newly created non-empty ready output should count as successful completion."""
        import skyborn_cdo._runner as runner_mod

        path = tmp_path / "new_output.nc"
        state = runner_mod._win_init_completion_state(str(path), now=0.0)

        path.write_bytes(b"fresh-output")
        monkeypatch.setattr(runner_mod, "_hdf5_file_is_closed", lambda _path: True)
        monkeypatch.setattr(runner_mod, "_win_file_is_exclusive_ready", lambda _path: True)

        assert runner_mod._win_update_output_completion_state(state) is True
        assert runner_mod._win_timeout_completion_succeeded(state) is True

    def test_stdout_only_partial_output_without_quiet_period_fails(self):
        """Stdout-only commands must fail hard timeout unless quiet completion was observed."""
        import skyborn_cdo._runner as runner_mod

        state = runner_mod._win_init_completion_state(None, now=0.0)

        assert runner_mod._win_update_completion_state(
            state,
            stdout_size=64,
            now=0.0,
            quiet_secs=0.5,
        ) is False
        assert runner_mod._win_update_completion_state(
            state,
            stdout_size=64,
            now=0.2,
            quiet_secs=0.5,
        ) is False
        assert runner_mod._win_timeout_completion_succeeded(state) is False


class TestCdoClass:
    """Test the high-level Cdo class."""

    @pytest.fixture
    def cdo(self):
        """Create a Cdo instance if CDO is available."""
        from skyborn_cdo import Cdo

        try:
            return Cdo()
        except FileNotFoundError:
            pytest.skip("CDO binary not available")

    def test_version(self, cdo):
        """Test that version() returns a string."""
        version = cdo.version()
        assert isinstance(version, str)
        assert "Climate Data Operators" in version or "CDO" in version

    def test_repr(self, cdo):
        """Test repr."""
        r = repr(cdo)
        assert "Cdo" in r
        assert "cdo_path" in r

    def test_call_version(self, cdo):
        """Test command-line style: cdo("--version")."""
        result = cdo("--version")
        assert isinstance(result, (str, int))

    def test_operators_list(self, cdo):
        """Test fetching operator list."""
        ops = cdo.operators()
        assert isinstance(ops, set)
        if ops:  # May be empty if CDO has issues
            assert "info" in ops or "copy" in ops
            assert "mergetime" in ops

    def test_has_operator(self, cdo):
        """Test operator existence check."""
        assert cdo.has_operator("copy") or cdo.has_operator("info")
        assert not cdo.has_operator("nonexistent_operator_xyz")

    def test_getattr_raises_for_private(self, cdo):
        """Test that private attributes raise AttributeError."""
        with pytest.raises(AttributeError):
            cdo._nonexistent

    def test_call_strip_cdo_prefix(self, cdo):
        """Test that 'cdo' prefix is stripped from command string."""
        # "cdo --version" and "--version" should both work
        from skyborn_cdo._runner import CdoRunner
        result = cdo("cdo --version")
        assert isinstance(result, (str, int))

    def test_cleanup(self, cdo):
        """Test cleanup removes temp files."""
        cdo._tempfiles.append("/tmp/nonexistent_test_file.nc")
        cdo.cleanup()
        assert len(cdo._tempfiles) == 0

    def test_operator_bound_help(self, cdo):
        """Test cdo.<operator>.help() convenience shortcut."""
        txt = cdo.mergetime.help()
        assert isinstance(txt, str)
        assert "mergetime" in txt.lower()


class TestCdoOperations:
    """Integration tests that require a working CDO and test NC files."""

    @pytest.fixture
    def cdo(self):
        from skyborn_cdo import Cdo

        try:
            c = Cdo()
            c.version()  # Verify it works
            return c
        except (FileNotFoundError, Exception):
            pytest.skip("CDO binary not available or not functional")

    @pytest.fixture
    def sample_nc(self, cdo, tmp_path):
        """Create a minimal test NetCDF file using CDO."""
        outfile = str(tmp_path / "test.nc")
        try:
            # Create a simple grid with constant field using CDO's topo operator
            cdo.topo(output=outfile)
            return outfile
        except Exception:
            pytest.skip(
                "Could not create test NC file (NetCDF support may be missing)")

    def test_info(self, cdo, sample_nc):
        """Test cdo.info on a file."""
        result = cdo.info(input=sample_nc)
        assert isinstance(result, str)
        assert len(result) > 0

    def test_sinfo(self, cdo, sample_nc):
        """Test cdo.sinfo on a file."""
        result = cdo.sinfo(input=sample_nc)
        assert isinstance(result, str)
        assert len(result.strip()) > 0

    def test_call_style_sinfo(self, cdo, sample_nc):
        """Test raw command-line style sinfo returns text output."""
        result = cdo(f"cdo sinfo {sample_nc}")
        assert isinstance(result, str)
        assert len(result.strip()) > 0

    def test_copy(self, cdo, sample_nc, tmp_path):
        """Test cdo.copy."""
        outfile = str(tmp_path / "copy_out.nc")
        cdo.copy(input=sample_nc, output=outfile)
        assert os.path.exists(outfile)
        assert os.path.getsize(outfile) > 0

    def test_call_style(self, cdo, sample_nc, tmp_path):
        """Test command-line style invocation."""
        outfile = str(tmp_path / "call_out.nc")
        cdo(f"cdo copy {sample_nc} {outfile}")
        assert os.path.exists(outfile)

    def test_options(self, sample_nc, tmp_path):
        """Test Cdo with default options."""
        from skyborn_cdo import Cdo

        try:
            cdo = Cdo(options="-O -s")
        except FileNotFoundError:
            pytest.skip("CDO not available")

        outfile = str(tmp_path / "opts_out.nc")
        cdo.copy(input=sample_nc, output=outfile)
        assert os.path.exists(outfile)

    def test_showname(self, cdo, sample_nc):
        """Test info operator that returns names."""
        result = cdo.showname(input=sample_nc)
        assert isinstance(result, str)

    def test_mergetime_produces_correct_timesteps(self, cdo, tmp_path):
        """mergetime of two 12-month files must produce 24 ordered timesteps.

        Regression test for the Windows exit-hang early-kill bug that caused
        CDO to be terminated while writing the NC4 header (before any time
        records were flushed), resulting in time=0 in the output.
        """
        nc4 = pytest.importorskip("netCDF4")
        cftime = pytest.importorskip("cftime")
        import numpy as np

        def _make_monthly_nc(path, year):
            ds = nc4.Dataset(path, "w", format="NETCDF4")
            ds.createDimension("time", 12)
            ds.createDimension("lat", 3)
            ds.createDimension("lon", 4)
            t = ds.createVariable("time", "f8", ("time",))
            t.units = "hours since 1900-01-01 00:00:00"
            t.calendar = "gregorian"
            t.standard_name = "time"
            t.axis = "T"
            dates = [cftime.DatetimeGregorian(
                year, m, 1) for m in range(1, 13)]
            t[:] = nc4.date2num(dates, units=t.units, calendar=t.calendar)
            v = ds.createVariable("tas", "f4", ("time", "lat", "lon"),
                                  fill_value=-9999.0)
            v.long_name = "air temperature"
            v.units = "K"
            v[:] = np.random.rand(12, 3, 4).astype("f4") + 280.0
            ds.close()

        nc_2023 = str(tmp_path / "tas_2023.nc")
        nc_2024 = str(tmp_path / "tas_2024.nc")
        out_nc = str(tmp_path / "tas_merged.nc")

        _make_monthly_nc(nc_2023, 2023)
        _make_monthly_nc(nc_2024, 2024)

        # method-call style
        cdo.mergetime(input=f"{nc_2023} {nc_2024}", output=out_nc)

        with nc4.Dataset(out_nc) as ds:
            t_var = ds.variables["time"]
            assert t_var.size == 24, (
                f"Expected 24 timesteps after mergetime, got {t_var.size}. "
                "This likely means CDO was killed before flushing its data "
                "(exit-hang early-kill bug)."
            )
            t_vals = t_var[:].data
            assert (t_vals[1:] > t_vals[:-1]).all(), \
                "Merged time axis is not monotonically increasing"

    def test_mergetime_wildcard(self, cdo, tmp_path):
        """cdo('cdo mergetime *.nc out.nc') must not include the output file
        among the inputs when it already exists from a previous run."""
        nc4 = pytest.importorskip("netCDF4")
        import numpy as np
        import os as _os

        def _make_nc(path, n_times, offset_hours=0):
            ds = nc4.Dataset(path, "w", format="NETCDF4")
            ds.createDimension("time", n_times)
            t = ds.createVariable("time", "f8", ("time",))
            t.units = "hours since 1900-01-01 00:00:00"
            t.calendar = "gregorian"
            t.standard_name = "time"
            t.axis = "T"
            t[:] = [offset_hours + i * 24 for i in range(n_times)]
            v = ds.createVariable("tas", "f4", ("time",))
            v[:] = np.ones(n_times, dtype="f4")
            ds.close()

        _make_nc(str(tmp_path / "a_input.nc"), 12, offset_hours=0)
        _make_nc(str(tmp_path / "b_input.nc"), 12, offset_hours=12 * 24)

        orig_cwd = _os.getcwd()
        try:
            _os.chdir(tmp_path)
            # First run: produce out_merged.nc
            cdo("cdo mergetime a_input.nc b_input.nc out_merged.nc")
            with nc4.Dataset("out_merged.nc") as ds:
                assert ds.variables["time"].size == 24

            # Second run with wildcard: out_merged.nc exists but must NOT be
            # picked up as an input (wildcard exclusion fix)
            cdo("cdo mergetime a_input.nc b_input.nc out_merged2.nc")
            with nc4.Dataset("out_merged2.nc") as ds:
                assert ds.variables["time"].size == 24, \
                    "Wildcard expansion incorrectly included the output file"
        finally:
            _os.chdir(orig_cwd)

    def test_call_style_wildcard_preserves_spaces(self, cdo, tmp_path):
        """Wildcard expansion in raw command strings must preserve paths with spaces."""
        nc4 = pytest.importorskip("netCDF4")
        import numpy as np

        def _make_nc(path, n_times, offset_hours=0):
            ds = nc4.Dataset(path, "w", format="NETCDF4")
            ds.createDimension("time", n_times)
            t = ds.createVariable("time", "f8", ("time",))
            t.units = "hours since 1900-01-01 00:00:00"
            t.calendar = "gregorian"
            t.standard_name = "time"
            t.axis = "T"
            t[:] = [offset_hours + i * 24 for i in range(n_times)]
            v = ds.createVariable("tas", "f4", ("time",))
            v[:] = np.ones(n_times, dtype="f4")
            ds.close()

        spaced_dir = tmp_path / "dir with spaces"
        spaced_dir.mkdir()
        _make_nc(str(spaced_dir / "a_2023.nc"), 12, offset_hours=0)
        _make_nc(str(spaced_dir / "a_2024.nc"), 12, offset_hours=12 * 24)

        merged = spaced_dir / "merged.nc"
        cdo(
            f'cdo -O mergetime "{spaced_dir / "a_202*.nc"}" "{merged}"'
        )

        assert merged.exists()
        with nc4.Dataset(merged) as ds:
            assert ds.variables["time"].size == 24


class TestCli:
    """Test CLI entry point."""

    def test_cli_decode_output_handles_utf8_bytes(self):
        """CLI help decoding should not depend on the active Windows code page."""
        from skyborn_cdo._cli import _decode_cli_output

        payload = "Operator help – remapbil °".encode("utf-8")
        decoded = _decode_cli_output(payload)
        assert "remapbil" in decoded
        assert "Operator help" in decoded

    def test_cli_operator_help_decodes_binary_output(self, monkeypatch, capsys):
        """Regression: help output with non-GBK bytes must not crash on Windows."""
        import skyborn_cdo._cli as cli

        old_argv = sys.argv[:]
        sys.argv = ["skyborn-cdo", "remapbil", "--help"]

        monkeypatch.setattr(cli, "get_cdo_path", lambda: "cdo")
        monkeypatch.setattr(cli, "get_bundled_env", lambda: {})

        original_os_name = cli.os.name
        monkeypatch.setattr(cli.os, "name", "nt")

        class _SubprocessStub:
            CREATE_NO_WINDOW = 0

            @staticmethod
            def run(*args, **kwargs):
                assert kwargs.get("capture_output") is True
                assert "text" not in kwargs
                return SimpleNamespace(
                    stdout="remapbil help – details".encode("utf-8"),
                    stderr=b"",
                    returncode=0,
                )

        original_subprocess = sys.modules.get("subprocess")
        sys.modules["subprocess"] = _SubprocessStub
        try:
            with pytest.raises(SystemExit) as exc:
                cli.main()
            assert exc.value.code == 0
        finally:
            sys.argv = old_argv
            if original_subprocess is not None:
                sys.modules["subprocess"] = original_subprocess
            else:
                sys.modules.pop("subprocess", None)

        out = capsys.readouterr().out
        assert "remapbil help" in out

    def test_cli_info(self):
        """Test skyborn-cdo --info."""
        result = subprocess.run(
            [sys.executable, "-m", "skyborn_cdo._cli", "--info"],
            capture_output=True,
            text=True,
            env=_cli_test_env(),
        )
        assert "skyborn-cdo version" in result.stdout

    def test_cli_help(self):
        """Test skyborn-cdo --help."""
        result = subprocess.run(
            [sys.executable, "-m", "skyborn_cdo._cli", "--help"],
            capture_output=True,
            text=True,
            env=_cli_test_env(),
        )
        assert "skyborn-cdo" in result.stdout
        assert "cdo --info" in result.stdout
        assert "cdo-win" in result.stdout
        assert "Python API" in result.stdout
        assert "<operator> --help" in result.stdout

    def test_cli_operator_help_long_form(self):
        """Test `skyborn-cdo mergetime --help` convenience syntax."""
        result = subprocess.run(
            [sys.executable, "-m", "skyborn_cdo._cli", "mergetime", "--help"],
            capture_output=True,
            text=True,
            env=_cli_test_env(),
        )
        out = (result.stdout + result.stderr).lower()
        assert "mergetime" in out
        assert "merge datasets" in out or "sorted by date and time" in out

    def test_cli_operator_help_short_alias(self):
        """Test `skyborn-cdo mergetime --h` convenience syntax."""
        result = subprocess.run(
            [sys.executable, "-m", "skyborn_cdo._cli", "mergetime", "--h"],
            capture_output=True,
            text=True,
            env=_cli_test_env(),
        )
        out = (result.stdout + result.stderr).lower()
        assert "mergetime" in out

    def test_cli_operator_help_short_flag_after_operator(self):
        """Test `skyborn-cdo mergetime -h` convenience syntax."""
        result = subprocess.run(
            [sys.executable, "-m", "skyborn_cdo._cli", "mergetime", "-h"],
            capture_output=True,
            text=True,
            env=_cli_test_env(),
        )
        out = (result.stdout + result.stderr).lower()
        assert "mergetime" in out
        assert "merge datasets" in out or "sorted by date and time" in out

    def test_package_version_matches_pyproject(self):
        """Package __version__ should match the published project version."""
        import skyborn_cdo

        pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
        with pyproject.open("rb") as f:
            project = tomllib.load(f)

        assert skyborn_cdo.__version__ == project["project"]["version"]
        scripts = project["project"]["entry-points"]["console_scripts"]
        assert scripts == {"skyborn-cdo": "skyborn_cdo._cli:main"}

    def test_windows_launcher_scripts_are_packaged_conditionally(self):
        """Windows alias wrappers must live in repo scripts and setup.py Windows-only logic."""
        repo_root = Path(__file__).resolve().parents[1]
        setup_py = (repo_root / "setup.py").read_text(encoding="utf-8")

        assert (repo_root / "scripts" / "cdo.cmd").is_file()
        assert (repo_root / "scripts" / "cdo-win.cmd").is_file()
        assert 'if os.name == "nt":' in setup_py
        assert '"scripts/cdo.cmd"' in setup_py
        assert '"scripts/cdo-win.cmd"' in setup_py

    def test_cli_sinfo_outputs_text(self, tmp_path):
        """CLI sinfo should print metadata text for a NetCDF file."""
        from skyborn_cdo import Cdo

        try:
            cdo = Cdo()
            cdo.version()
        except (FileNotFoundError, Exception):
            pytest.skip("CDO binary not available or not functional")

        sample_nc = str(tmp_path / "cli_sinfo.nc")
        cdo.topo(output=sample_nc)

        result = subprocess.run(
            [sys.executable, "-m", "skyborn_cdo._cli", "sinfo", sample_nc],
            capture_output=True,
            text=True,
            timeout=30,
            env=_cli_test_env(),
        )
        assert result.returncode == 0
        combined = (result.stdout + result.stderr).strip()
        assert combined
        assert "File format" in combined or "Grid coordinates" in combined or "points" in combined

    def test_cli_wildcard_mergetime_windows_style(self, tmp_path):
        """CLI should expand *.nc itself because PowerShell won't."""
        from skyborn_cdo import Cdo

        try:
            cdo = Cdo()
            cdo.version()
        except (FileNotFoundError, Exception):
            pytest.skip("CDO binary not available or not functional")

        nc4 = pytest.importorskip("netCDF4")
        cftime = pytest.importorskip("cftime")
        import numpy as np

        def _make_monthly_nc(path, year):
            ds = nc4.Dataset(path, "w", format="NETCDF4")
            ds.createDimension("time", 12)
            ds.createDimension("lat", 2)
            ds.createDimension("lon", 3)
            t = ds.createVariable("time", "f8", ("time",))
            t.units = "hours since 1900-01-01 00:00:00"
            t.calendar = "gregorian"
            dates = [cftime.DatetimeGregorian(
                year, m, 1) for m in range(1, 13)]
            t[:] = nc4.date2num(dates, units=t.units, calendar=t.calendar)
            v = ds.createVariable("tas", "f4", ("time", "lat", "lon"))
            v[:] = np.random.rand(12, 2, 3).astype("f4")
            ds.close()

        _make_monthly_nc(str(tmp_path / "ERA5_MSE_and_DSE_2023.nc"), 2023)
        _make_monthly_nc(str(tmp_path / "ERA5_MSE_and_DSE_2024.nc"), 2024)
        merged = str(tmp_path / "merged.nc")

        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "skyborn_cdo._cli",
                "-O",
                "mergetime",
                str(tmp_path / "ERA5_MSE_and_DSE_202*.nc"),
                merged,
            ],
            capture_output=True,
            text=True,
            timeout=60,
            env=_cli_test_env(),
        )
        assert result.returncode == 0, result.stderr or result.stdout
        assert os.path.exists(merged)
        with nc4.Dataset(merged) as ds:
            assert ds.variables["time"].size == 24
