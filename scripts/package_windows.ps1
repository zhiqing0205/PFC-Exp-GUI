param(
  [string]$BuildDir = "build",
  [string]$Config = "Release",
  [string]$OutDir = "dist"
)

$ErrorActionPreference = "Stop"

cmake -S . -B $BuildDir -DBUILD_SIM=OFF
cmake --build $BuildDir --config $Config

$exe = Get-ChildItem -Path $BuildDir -Recurse -Filter "qt_gui_client.exe" | Select-Object -First 1
if (-not $exe) {
  throw "qt_gui_client.exe not found under $BuildDir. Check your generator output."
}

$stage = Join-Path $BuildDir "package\qt-experiment-runner"
New-Item -ItemType Directory -Force -Path $stage | Out-Null
Copy-Item $exe.FullName (Join-Path $stage "qt_gui_client.exe") -Force

$windeployqt = Get-Command windeployqt.exe -ErrorAction SilentlyContinue
if (-not $windeployqt) {
  throw "windeployqt.exe not found in PATH. Install Qt and add <Qt>/bin to PATH."
}

& $windeployqt.Path (Join-Path $stage "qt_gui_client.exe") --qmldir (Join-Path $PSScriptRoot "..\src\gui\qml") --compiler-runtime

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$zip = Join-Path $OutDir "qt-experiment-runner_windows.zip"
if (Test-Path $zip) { Remove-Item $zip -Force }
Compress-Archive -Path (Join-Path $stage "*") -DestinationPath $zip

Write-Host "Done: $zip"

