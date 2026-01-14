#!/usr/bin/env bash
set -euo pipefail

gen_args=()
if command -v ninja >/dev/null 2>&1; then
  gen_args=(-G Ninja)
fi

cmake -S . -B build "${gen_args[@]}" -DCMAKE_BUILD_TYPE=Release -DBUILD_SIM=OFF
cmake --build build

APP="build/bin/qt_gui_client.app"
if [[ ! -d "$APP" ]]; then
  # fallback search
  APP="$(find build -maxdepth 4 -type d -name 'qt_gui_client.app' | head -n 1 || true)"
fi
if [[ -z "${APP:-}" || ! -d "$APP" ]]; then
  echo "qt_gui_client.app not found. Ensure you built on macOS and qt_gui_client is a bundle." >&2
  exit 1
fi

if ! command -v macdeployqt >/dev/null 2>&1; then
  echo "macdeployqt not found in PATH. Ensure Qt is installed and <Qt>/bin is in PATH." >&2
  exit 1
fi

DMG_NAME="$(basename "$APP" .app).dmg"
rm -f "./$DMG_NAME"

macdeployqt "$APP" -qmldir=src/gui/qml -dmg

DMG=""
if [[ -f "./$DMG_NAME" ]]; then
  DMG="./$DMG_NAME"
else
  APP_DIR="$(dirname "$APP")"
  DMG="$(ls -t "./"*.dmg "$APP_DIR"/*.dmg 2>/dev/null | head -n 1 || true)"
fi
if [[ -z "${DMG:-}" || ! -f "$DMG" ]]; then
  echo "macdeployqt succeeded but no .dmg was found (expected: ./$DMG_NAME)" >&2
  exit 1
fi

mkdir -p dist
rm -f dist/*.dmg
cp -f "$DMG" "dist/qt-experiment-runner_macos.dmg"

echo "Done: dist/qt-experiment-runner_macos.dmg"
