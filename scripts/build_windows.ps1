param(
  [string]$BuildDir = "build",
  [string]$Config = "Release"
)

$ErrorActionPreference = "Stop"

cmake -S . -B $BuildDir -DBUILD_SIM=OFF
cmake --build $BuildDir --config $Config

Write-Host "Done. Look for qt_gui_client(.exe) under $BuildDir\\bin or the generator output directory."

