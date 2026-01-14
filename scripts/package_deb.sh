#!/usr/bin/env bash
set -euo pipefail

gen_args=()
if command -v ninja >/dev/null 2>&1; then
  gen_args=(-G Ninja)
fi

cmake -S . -B build "${gen_args[@]}" -DCMAKE_BUILD_TYPE=Release
cmake --build build

rm -rf build/stage
cmake --install build --prefix build/stage

cpack --config build/CPackConfig.cmake -G DEB

mkdir -p dist
rm -f dist/*.deb
mv -f qt-experiment-runner_*.deb dist/

echo "Done. Look for *.deb under ./dist/."
