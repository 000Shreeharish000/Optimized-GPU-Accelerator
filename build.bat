@echo off
rem Build PRAMANA with MSVC + CMake + Ninja (CMake and Ninja ship with Visual Studio / Build Tools).
rem Usage: build.bat [Release|Debug]
rem Finds Visual Studio automatically (vswhere); override with PRAMANA_VCVARS=<path to vcvars64.bat>.
setlocal
set CONFIG=%1
if "%CONFIG%"=="" set CONFIG=Release
cd /d "%~dp0"

where cl >nul 2>nul
if not errorlevel 1 goto :build
if not "%PRAMANA_VCVARS%"=="" goto :vcvars
set VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe
if not exist "%VSWHERE%" (
  echo ERROR: Visual Studio / Build Tools with "Desktop development with C++" not found.
  echo        Install it, or set PRAMANA_VCVARS to your vcvars64.bat.
  exit /b 1
)
for /f "usebackq delims=" %%i in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set VSDIR=%%i
set PRAMANA_VCVARS=%VSDIR%\VC\Auxiliary\Build\vcvars64.bat
:vcvars
if not exist "%PRAMANA_VCVARS%" (
  echo ERROR: vcvars64.bat not found at "%PRAMANA_VCVARS%"
  exit /b 1
)
call "%PRAMANA_VCVARS%" >nul

:build
if not exist build\build.ninja cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=%CONFIG% || exit /b 1
cmake --build build || exit /b 1
echo Build OK: build\pramana.exe  build\libpramana.dll  build\pramana_tests.exe
