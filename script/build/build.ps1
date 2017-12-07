#!/usr/bin/env powershell

#$env:CC = "clang-cl.exe"
#$env:CXX = "clang-cl.exe"

Remove-Item build -Force -Recurse -ErrorAction SilentlyContinue
New-Item -Path .\build -ItemType directory -Force
Set-Location build

cmake -G "Ninja" "-DCMAKE_TOOLCHAIN_FILE=$env:WORKSPACE\vcpkg\scripts\buildsystems\vcpkg.cmake" "-DVCPKG_CHAINLOAD_TOOLCHAIN_FILE=$env:WORKSPACE\sysconfig\cmake\physycom_toolchain.cmake" "-DCMAKE_BUILD_TYPE=Release" ..
cmake --build .
Set-Location ..
