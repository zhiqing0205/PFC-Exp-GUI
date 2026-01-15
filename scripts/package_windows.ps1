param(
  [string]$BuildDir = "build",
  [string]$Config = "Release",
  [string]$OutDir = "dist"
)

$ErrorActionPreference = "Stop"

$cmakeArgs = @("-DBUILD_SIM=OFF")
if ($env:QT_EXPERIMENT_RUNNER_VERSION) {
  $cmakeArgs += "-DQT_EXPERIMENT_RUNNER_VERSION=$($env:QT_EXPERIMENT_RUNNER_VERSION)"
}

cmake -S . -B $BuildDir @cmakeArgs
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
$arch = $env:PROCESSOR_ARCHITECTURE
if ($arch -eq "AMD64") { $arch = "x64" }
elseif ($arch -eq "ARM64") { $arch = "arm64" }
else { $arch = $arch.ToLowerInvariant() }

$zip = Join-Path $OutDir "qt-experiment-runner_windows_$arch.zip"
if (Test-Path $zip) { Remove-Item $zip -Force }
Compress-Archive -Path (Join-Path $stage "*") -DestinationPath $zip

Write-Host "Done: $zip"
