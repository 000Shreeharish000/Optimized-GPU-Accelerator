@echo off
rem Build PRAMANA with MSVC + CMake + Ninja (all shipped with Visual Studio).
rem Usage: build.bat [Release|Debug]
setlocal
set CONFIG=%1
if "%CONFIG%"=="" set CONFIG=Release
if "%PRAMANA_VCVARS%"=="" set PRAMANA_VCVARS=C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat
where cl >nul 2>nul
if errorlevel 1 call "%PRAMANA_VCVARS%" >nul
cd /d "%~dp0"
if not exist build\build.ninja cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=%CONFIG% || exit /b 1
cmake --build build || exit /b 1
echo Build OK: build\pramana.exe  build\libpramana.dll  build\pramana_tests.exe
