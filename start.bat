@echo off
rem ===========================================================================
rem PRAMANA - one command to start everything:
rem   finds Python, builds the solver (incremental), prepares the benchmark data
rem   if missing, and opens the PRAMANA terminal console.
rem Usage: start.bat            (any extra arguments are passed to the console)
rem ===========================================================================
setlocal EnableDelayedExpansion
cd /d "%~dp0"
chcp 65001 >nul
title PRAMANA

for /f %%a in ('echo prompt $E^| cmd') do set "ESC=%%a"
set "O=%ESC%[38;2;255;140;26m"
set "Y=%ESC%[38;2;255;214;10m"
set "R=%ESC%[38;2;239;68;68m"
set "D=%ESC%[38;2;135;135;135m"
set "N=%ESC%[0m"

echo.
echo  %O%PRAMANA%N% %D%- certified LP / MILP / QP solver%N%
echo.

rem ---- 1. Python ---------------------------------------------------------------
set PYEXE=
set PYARGS=
if defined PYTHON (
  "%PYTHON%" -c "import sys" >nul 2>nul && set "PYEXE=%PYTHON%"
)
if not defined PYEXE (
  python -c "import sys; assert sys.version_info >= (3, 8)" >nul 2>nul && set PYEXE=python
)
if not defined PYEXE (
  for %%v in (3.13 3.12 3.11 3.10 3.9 3) do (
    if not defined PYEXE (
      py -%%v -c "import sys" >nul 2>nul && (set PYEXE=py& set PYARGS=-%%v)
    )
  )
)
if not defined PYEXE (
  echo  %R%[xx] Python 3.8+ not found.%N% Install it from python.org, or set PYTHON=C:\path\to\python.exe
  exit /b 1
)
echo  %Y%[ok]%N% python      %D%!PYEXE! !PYARGS!%N%

rem ---- 2. Build (incremental; a no-op when nothing changed) --------------------
if not exist build mkdir build
call "%~dp0build.bat" > build\start_build.log 2>&1
if errorlevel 1 (
  if exist build\pramana.exe (
    echo  %O%[--]%N% build       %D%could not rebuild ^(no Visual Studio?^) - using the existing build\pramana.exe%N%
  ) else (
    echo  %R%[xx] build failed%N% - see build\start_build.log. Needs Visual Studio 2022+ with "Desktop development with C++".
    exit /b 1
  )
) else (
  echo  %Y%[ok]%N% build       %D%build\pramana.exe%N%
)

rem ---- 3. Benchmark data ---------------------------------------------------------
if not exist data\netlib\afiro.mps.gz (
  echo  %O%[..]%N% data        %D%downloading Netlib / MIPLIB / Maros-Meszaros ^(first run only^)%N%
  "!PYEXE!" !PYARGS! bench\fetch_data.py netlib netlib-infeas miplib3 maros || echo  %O%[--]%N% data download failed - public sets unavailable, generated models still work
)
if not exist data\gen\plan_10x12.mps "!PYEXE!" !PYARGS! gen\refinery.py suite --out data\gen >nul
if not exist data\adversarial\scaled_afiro.mps "!PYEXE!" !PYARGS! gen\adversarial.py --out data\adversarial >nul
echo  %Y%[ok]%N% data        %D%benchmark and industrial models ready%N%

rem ---- 4. Console ---------------------------------------------------------------
set "PYTHONPATH=%~dp0python;%PYTHONPATH%"
"!PYEXE!" !PYARGS! -m pramana.tui %*
exit /b %ERRORLEVEL%
