@echo off
call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>nul
cd /d "%~dp0"
cl /nologo /O2 /w /D_CRT_SECURE_NO_WARNINGS emps.c /Fe:emps.exe
