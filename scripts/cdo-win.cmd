@echo off
if exist "%~dp0skyborn-cdo.exe" goto run_launcher
"%~dp0python.exe" -m skyborn_cdo %*
exit /b %ERRORLEVEL%

:run_launcher
"%~dp0skyborn-cdo.exe" %*
exit /b %ERRORLEVEL%
