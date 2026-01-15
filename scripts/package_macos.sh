#!/usr/bin/env bash
set -euo pipefail

set -x

gen_args=()
if command -v ninja >/dev/null 2>&1; then
  gen_args=(-G Ninja)
fi

cmake_args=(-DCMAKE_BUILD_TYPE=Release -DBUILD_SIM=OFF)
if [[ -n "${QT_EXPERIMENT_RUNNER_VERSION:-}" ]]; then
  cmake_args+=(-DQT_EXPERIMENT_RUNNER_VERSION="${QT_EXPERIMENT_RUNNER_VERSION}")
fi

cmake -S . -B build "${gen_args[@]}" "${cmake_args[@]}"
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

macdeployqt_bin="$(command -v macdeployqt || true)"
if [[ -z "${macdeployqt_bin:-}" ]]; then
  for base in "${QT_ROOT_DIR:-}" "${QTDIR:-}" "${Qt6_Dir:-}" "${Qt6_DIR:-}" "${Qt_DIR:-}"; do
    if [[ -n "${base:-}" && -x "$base/bin/macdeployqt" ]]; then
      macdeployqt_bin="$base/bin/macdeployqt"
      break
    fi
  done
fi
if [[ -z "${macdeployqt_bin:-}" ]]; then
  echo "macdeployqt not found. PATH=$PATH" >&2
  echo "Hint: ensure Qt bin dir is in PATH, or set QT_ROOT_DIR to your Qt install." >&2
  exit 1
fi

mkdir -p dist
rm -f dist/*.dmg

"$macdeployqt_bin" "$APP" -qmldir=src/gui/qml

if ! command -v hdiutil >/dev/null 2>&1; then
  echo "hdiutil not found; cannot create .dmg on this environment." >&2
  exit 1
fi

STAGE_DIR="build/dmg_stage"
rm -rf "$STAGE_DIR"
mkdir -p "$STAGE_DIR"
cp -R "$APP" "$STAGE_DIR/"

arch="$(uname -m)"
out="dist/qt-experiment-runner_macos_${arch}.dmg"
hdiutil create -volname "Qt Experiment Runner" -srcfolder "$STAGE_DIR" -ov -format UDZO "$out"

echo "Done: $out"
