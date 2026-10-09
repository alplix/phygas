@echo off
rem Phygas CUDA GPU app for Windows x86_64 (MSVC + nvcc, the same toolchain and BOINC MSVC libraries as
rem bitboinc's other CUDA builds). Output: build\gpuwin\phygas_cuda_v%V%_windows_x86_64.exe
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul 2>nul
set NVCC="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.3\bin\nvcc.exe"
set V=2
if not "%PHYGAS_VERSION%"=="" set V=%PHYGAS_VERSION%
set SRC=%~dp0..\src
set BOINC=C:\Users\Alp\Documents\bitboinc-app\gpu_test_local\boinc
set GENCODES=-gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_89,code=sm_89 -gencode arch=compute_90,code=sm_90 -gencode arch=compute_100,code=sm_100 -gencode arch=compute_120,code=sm_120 -gencode arch=compute_120,code=compute_120
set OUT=%~dp0..\build\gpuwin
if not exist %OUT% mkdir %OUT%
cd /d %OUT%
for %%f in (coords potential montecarlo phygas) do cl /nologo /O2 /std:c11 /MT /I %SRC% /c %SRC%\%%f.c /Fo%%f.obj || exit /b 1
%NVCC% -O3 -std=c++17 --fmad=false %GENCODES% -I %SRC% -c %SRC%\phygas_cuda.cu -o phygas_cuda.obj -Xcompiler /MT || exit /b 1
%NVCC% -O2 -std=c++17 -DPG_GPU -DPG_CUDA -I %SRC% -I %BOINC%\api -I %BOINC%\lib -c %SRC%\main_phygas_boinc.cpp -o main.obj -Xcompiler /MT || exit /b 1
%NVCC% -o phygas_cuda_v%V%_windows_x86_64.exe main.obj phygas_cuda.obj coords.obj potential.obj montecarlo.obj phygas.obj %BOINC%\win_build\Build\x64\Release\libboincapi.lib %BOINC%\win_build\Build\x64\Release\libboinc.lib -Xcompiler /MT -Xlinker advapi32.lib -Xlinker shell32.lib -Xlinker user32.lib -Xlinker ole32.lib -Xlinker psapi.lib -Xlinker /SUBSYSTEM:CONSOLE || exit /b 1
echo === DONE ===
