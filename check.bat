@echo off
rem ===========================================================================
rem PRAMANA self-check: builds (if needed), runs all tests, the from-scratch
rem audit and a set of certified solves (LP / infeasible / unbounded / badly
rem scaled / MILP / QP / refinery parametric), each re-verified independently.
rem Usage: check.bat            (about 1-3 minutes after the first build)
rem Exit code 0 = every check passed.
rem ===========================================================================
setlocal EnableDelayedExpansion
cd /d "%~dp0"
set EXE=build\pramana.exe
set OUT=build\check
set FAILS=0
set PASSES=0

rem ---- Python (only needed for tests, verifier, audit and data generation)
if "%PYTHON%"=="" (
  where python >nul 2>nul && set PYTHON=python
  if "!PYTHON!"=="" where py >nul 2>nul && set PYTHON=py
)
if "%PYTHON%"=="" (
  echo WARNING: no Python found - Python-based checks are skipped. Set PYTHON=path\to\python.exe
)

echo.
echo === [1/6] Build
call "%~dp0build.bat"
if errorlevel 1 (
  if not exist %EXE% (echo FAIL: build & exit /b 1)
  echo WARNING: could not rebuild - using the existing %EXE%
)
if not exist %OUT% mkdir %OUT%

echo.
echo === [2/6] Benchmark data
if not exist data\netlib\afiro.mps.gz (
  if "%PYTHON%"=="" (echo FAIL: data missing and no Python to fetch it & exit /b 1)
  "%PYTHON%" bench\fetch_data.py netlib netlib-infeas miplib3 maros || (echo FAIL: fetch_data - internet needed & exit /b 1)
)
if not exist data\gen\plan_10x12.mps "%PYTHON%" gen\refinery.py suite --out data\gen
if not exist data\adversarial\scaled_afiro.mps "%PYTHON%" gen\adversarial.py --out data\adversarial
echo data OK

echo.
echo === [3/6] Unit tests (C++)
build\pramana_tests.exe
call :result "C++ unit tests" %ERRORLEVEL% 0

if not "%PYTHON%"=="" (
  echo.
  echo === [3b] End-to-end tests: Python bindings, exact verifier, generators
  "%PYTHON%" -m pytest tests\python -q
  call :result "Python end-to-end tests" !ERRORLEVEL! 0

  echo.
  echo === [4/6] From-scratch audit: no solver library linked or included
  "%PYTHON%" bench\audit_deps.py > %OUT%\audit.txt
  call :result "from-scratch dependency audit, docs\FROM_SCRATCH_AUDIT.md" !ERRORLEVEL! 0
)

echo.
echo === [5/6] Certified solves (exit code 0 = proven AND certificate accepted)
call :solve "LP  Netlib AFIRO"                          data\netlib\afiro.mps.gz
call :solve "LP  Netlib PILOT87, ill-conditioned"      data\netlib\pilot87.mps.gz
call :solve "LP  Netlib DEGEN3, highly degenerate"     data\netlib\degen3.mps.gz
call :solve "LP  infeasible KLEIN3, Farkas proof"      data\netlib_infeas\klein3.mps
call :solve "LP  unbounded, primal ray"                tests\data\unbounded.mps
call :solve "LP  badly scaled AFIRO, 1e+-6 scaling"           data\adversarial\scaled_afiro.mps
call :solve "LP  Klee-Minty n=20"                       data\adversarial\kleeminty_20.mps
call :solve "LP  refinery planning 10x12"               data\gen\plan_10x12.mps
call :solve "MILP MIPLIB P0033"                         data\miplib3\p0033.mps.gz
call :solve "MILP MIPLIB P0282"                         data\miplib3\p0282.mps.gz
call :solve "MILP unit commitment 10x24"                data\gen\uc_10x24.mps
call :solve "MILP crude unloading, tight formulation"  data\gen\unload_tight.mps
call :solve "QP  Maros-Meszaros HS21"                   data\maros\HS21.qps
call :solve "QP  economic dispatch 50x24"               data\gen\dispatch_50x24.qps

echo.
echo === [6/6] Independent re-verification of saved results
%EXE% data\netlib\afiro.mps.gz --log 0 --json %OUT%\afiro.json --vectors >nul
%EXE% verify data\netlib\afiro.mps.gz %OUT%\afiro.json
call :result "C++ certifier re-check, afiro" %ERRORLEVEL% 0
if not "%PYTHON%"=="" (
  pushd python
  "%PYTHON%" -m pramana.verify ..\data\netlib\afiro.mps.gz ..\%OUT%\afiro.json --exact-basis
  call :result "exact-rational verifier + exact basis, afiro" !ERRORLEVEL! 0
  popd
)
%EXE% data\netlib_infeas\klein3.mps --log 0 --json %OUT%\klein3.json --vectors >nul
if not "%PYTHON%"=="" (
  pushd python
  "%PYTHON%" -m pramana.verify ..\data\netlib_infeas\klein3.mps ..\%OUT%\klein3.json
  call :result "exact-rational Farkas check, klein3" !ERRORLEVEL! 0
  popd
)
%EXE% parametric data\gen\plan_10x12.mps --col BUY_CR03_0 --kind upper --from 0 --to 120 --log 0 --json %OUT%\param.json > %OUT%\param.txt
call :result "certified parametric crude valuation" %ERRORLEVEL% 0

echo.
echo === GPU (informational, not a pass/fail check)
%EXE% info | findstr /C:"gpu" /C:"compute" /C:"kernels"

echo.
echo ===========================================================================
echo  PASSED %PASSES%   FAILED %FAILS%
echo ===========================================================================
if %FAILS% gtr 0 exit /b 1
exit /b 0

:solve
%EXE% %2 --log 0 --time 120 > %OUT%\last.txt
set RC=%ERRORLEVEL%
for /f "tokens=2" %%s in ('findstr /B /C:"Status" %OUT%\last.txt') do set ST=%%s
for /f "tokens=2" %%o in ('findstr /B /C:"Objective" %OUT%\last.txt') do set OBJ=%%o
findstr /B /C:"Certificate      ACCEPTED" %OUT%\last.txt >nul
if errorlevel 1 set RC=9
call :result "%~1 : !ST! !OBJ!" %RC% 0
set OBJ=
exit /b 0

:result
if "%~2"=="%~3" (
  echo   PASS  %~1
  set /a PASSES+=1
) else (
  echo   FAIL  %~1   [exit code %~2]
  set /a FAILS+=1
)
exit /b 0
