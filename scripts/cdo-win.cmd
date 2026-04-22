@echo off
setlocal
if exist "%~dp0skyborn-cdo.exe" (
    "%~dp0skyborn-cdo.exe" %*
    set "RC=%ERRORLEVEL%"
    exit /b %RC%
)
"%~dp0python.exe" -m skyborn_cdo %*
set "RC=%ERRORLEVEL%"
exit /b %RC%
