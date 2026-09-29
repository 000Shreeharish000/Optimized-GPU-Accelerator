@echo off
rem PRAMANA end-to-end demo (see docs/DEMO.md). Run after build.bat and bench\fetch_data.py.
setlocal
cd /d "%~dp0"
set EXE=build\pramana.exe
if "%PYTHON%"=="" set PYTHON=python

echo.
echo === 1. From-scratch audit (imported DLLs, source scan, lines of code)
%PYTHON% bench\audit_deps.py | findstr /C:"pramana.exe" /C:"Result:" /C:"total"

echo.
echo === 2. GPU and build report
%EXE% info | findstr /C:"gpu" /C:"kernels" /C:"compute"

echo.
echo === 3. Ill-conditioned Netlib LP (PILOT87) solved and certified
%EXE% data\netlib\pilot87.mps.gz --log 0 --json build\demo_pilot87.json --vectors

echo.
echo === 4. Independent exact-rational re-verification of that result (own parser, Fractions)
pushd python
%PYTHON% -m pramana.verify ..\data\netlib\pilot87.mps.gz ..\build\demo_pilot87.json | findstr /C:"verdict" /C:"verified"
popd

echo.
echo === 5. Badly scaled model (rows/cols scaled by 10^+-6): true optimum -464.7531 with certificate
%EXE% data\adversarial\scaled_afiro.mps --log 0 | findstr /C:"Status" /C:"Objective" /C:"safe_dual"

echo.
echo === 6. Infeasible LP: verified Farkas certificate
%EXE% data\netlib_infeas\klein3.mps --log 0 | findstr /C:"Status" /C:"farkas"

echo.
echo === 7. Race: CPU dual simplex vs GPU PDHG on a refinery planning LP
%EXE% data\gen\plan_20x12.mps --algo race --log 1 | findstr /C:"engine" /C:"Engine" /C:"Status" /C:"Time"

echo.
echo === 8. MILP: unit commitment, cuts on (then off, 60 s limit)
%EXE% data\gen\uc_10x24.mps --log 0 | findstr /C:"Status" /C:"Objective" /C:"Time"
%EXE% data\gen\uc_10x24.mps --log 0 --no-cuts --time 60 | findstr /C:"Status" /C:"Time"

echo.
echo === 9. MILP racing on 4 cores (bell5)
%EXE% data\miplib3\bell5.mps.gz --log 0 --mip-threads 4 | findstr /C:"Status" /C:"Objective" /C:"Time"

echo.
echo === 10. Crude valuation: exact certified value curve + strategy comparison
%EXE% parametric data\gen\plan_10x12.mps --col BUY_CR03_0 --kind upper --from 0 --to 120 --json build\demo_param.json --log 0
%EXE% family data\gen\plan_10x12.mps --col BUY_CR03_0 --kind upper --from 0 --to 120 --cases 32 --log 0

echo.
echo === Done. Open web\index.html and drop build\demo_pilot87.json or build\demo_param.json onto it.
