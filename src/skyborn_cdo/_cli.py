"""
CLI entry point for skyborn-cdo.

Allows running CDO commands directly via: skyborn-cdo <cdo args>
or checking the installation via: skyborn-cdo --info
"""

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

    # Pass-through to CDO
    import os
    import subprocess

    cdo_path = get_cdo_path()
    env = get_bundled_env()

    _creationflags = subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0

    # For help commands (-h operator), capture and print output
    # because CDO writes help to stderr and may return non-zero.
    if args[0] in ("-h",) and len(args) >= 2:
        result = subprocess.run(
            [cdo_path] + args,
            capture_output=True,
            text=True,
            env=env,
            creationflags=_creationflags,
        )
        output = result.stdout.strip()
        err = result.stderr.strip()
        if output:
            print(output)
        if err:
            print(err)
        sys.exit(0)

    result = subprocess.run(
        [cdo_path] + args,
        env=env,
        creationflags=_creationflags,
    )
    sys.exit(result.returncode)


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
