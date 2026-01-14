#!/usr/bin/env bash
set -euo pipefail

sudo apt-get update

sudo apt-get install -y \
  build-essential \
  cmake \
  ninja-build \
  pkg-config \
  qt6-base-dev \
  qt6-declarative-dev \
  libqt6quickcontrols2-6 \
  qml6-module-qtquick \
  qml6-module-qtquick-controls \
  qml6-module-qtquick-layouts \
  qml6-module-qtquick-dialogs \
  openmpi-bin \
  libopenmpi-dev \
  libfftw3-dev \
  libfftw3-mpi-dev

echo "Done. You can now build with:"
echo "  cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release"
echo "  cmake --build build"
