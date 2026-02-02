#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/package-linux.sh --version <tag> --deb-arch <amd64|arm64> --appimage-arch <x86_64|aarch64> [--build-dir build] [--out-dir dist]

Inputs:
  --version         Git tag name (e.g. 0.3.10 or v0.3.10)
  --deb-arch        Debian arch for .deb (amd64 / arm64)
  --appimage-arch   AppImage arch (x86_64 / aarch64)
  --build-dir       CMake build dir (default: build)
  --out-dir         Output dir for final artifacts (default: dist)
EOF
}

VERSION=""
DEB_ARCH=""
APPIMAGE_ARCH=""
BUILD_DIR="build"
OUT_DIR="dist"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --version) VERSION="${2:-}"; shift 2 ;;
    --deb-arch) DEB_ARCH="${2:-}"; shift 2 ;;
    --appimage-arch) APPIMAGE_ARCH="${2:-}"; shift 2 ;;
    --build-dir) BUILD_DIR="${2:-}"; shift 2 ;;
    --out-dir) OUT_DIR="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "${VERSION}" || -z "${DEB_ARCH}" || -z "${APPIMAGE_ARCH}" ]]; then
  echo "Missing required args." >&2
  usage
  exit 2
fi

APP_NAME="MID Nano"
APP_ID="mid-nano"

mkdir -p "${OUT_DIR}"
OUT_DIR="$(cd "${OUT_DIR}" && pwd)"
BUILD_DIR="$(cd "${BUILD_DIR}" && pwd)"

BIN_GUI="${BUILD_DIR}/bin/${APP_NAME}"
BIN_CLI="${BUILD_DIR}/bin/pfc-exp-cli"

if [[ ! -f "${BIN_GUI}" ]]; then
  echo "GUI binary not found: ${BIN_GUI}" >&2
  exit 1
fi
if [[ ! -f "${BIN_CLI}" ]]; then
  echo "CLI binary not found: ${BIN_CLI}" >&2
  exit 1
fi

DEB_VERSION="${VERSION#v}"

tmp="$(mktemp -d)"
cleanup() { rm -rf "${tmp}"; }
trap cleanup EXIT

###############################################################################
# .deb (system-deps package)
###############################################################################

deb_root="${tmp}/debroot"
mkdir -p "${deb_root}/DEBIAN" \
         "${deb_root}/opt/${APP_ID}" \
         "${deb_root}/usr/bin" \
         "${deb_root}/usr/share/applications" \
         "${deb_root}/usr/share/icons/hicolor/scalable/apps"

cp "${BIN_GUI}" "${deb_root}/opt/${APP_ID}/${APP_NAME}"
cp "${BIN_CLI}" "${deb_root}/opt/${APP_ID}/pfc-exp-cli"
chmod +x "${deb_root}/opt/${APP_ID}/${APP_NAME}" "${deb_root}/opt/${APP_ID}/pfc-exp-cli"

cat > "${deb_root}/usr/bin/${APP_ID}" <<EOF
#!/usr/bin/env sh
set -eu
DIR="/opt/${APP_ID}"
exec "\${DIR}/${APP_NAME}" "\$@"
EOF
chmod +x "${deb_root}/usr/bin/${APP_ID}"

install -m 0644 "packaging/linux/mid-nano.desktop" "${deb_root}/usr/share/applications/${APP_ID}.desktop"
install -m 0644 "packaging/mid-nano.svg" "${deb_root}/usr/share/icons/hicolor/scalable/apps/${APP_ID}.svg"

cat > "${deb_root}/DEBIAN/control" <<EOF
Package: ${APP_ID}
Version: ${DEB_VERSION}
Section: science
Priority: optional
Architecture: ${DEB_ARCH}
Maintainer: MID Nano CI <noreply@example.com>
Homepage: https://github.com/zhiqing0205/PFC-Exp-GUI
Depends: libc6, libstdc++6, libgcc-s1, libqt6core6, libqt6gui6, libqt6widgets6, libfftw3-3
Description: MID Nano (Qt6) experiment GUI
 A small GUI to configure experiment parameters and launch pfc-exp-cli.
EOF

deb_out="${OUT_DIR}/${APP_NAME}-${VERSION}-linux-${DEB_ARCH}.deb"
dpkg-deb --build --root-owner-group "${deb_root}" "${deb_out}"

###############################################################################
# .AppImage (bundled Qt)
###############################################################################

if ! command -v linuxdeploy >/dev/null 2>&1; then
  echo "linuxdeploy is required in PATH for AppImage packaging." >&2
  exit 1
fi

appdir="${tmp}/AppDir"
mkdir -p "${appdir}/usr/bin" \
         "${appdir}/usr/share/applications" \
         "${appdir}/usr/share/icons/hicolor/scalable/apps"

cp "${BIN_GUI}" "${appdir}/usr/bin/${APP_NAME}"
cp "${BIN_CLI}" "${appdir}/usr/bin/pfc-exp-cli"
chmod +x "${appdir}/usr/bin/${APP_NAME}" "${appdir}/usr/bin/pfc-exp-cli"

cat > "${appdir}/usr/bin/${APP_ID}" <<EOF
#!/usr/bin/env sh
set -eu
DIR="\$(CDPATH= cd -- "\$(dirname -- "\$0")" && pwd)"
exec "\${DIR}/${APP_NAME}" "\$@"
EOF
chmod +x "${appdir}/usr/bin/${APP_ID}"

install -m 0644 "packaging/linux/mid-nano.desktop" "${appdir}/usr/share/applications/${APP_ID}.desktop"
install -m 0644 "packaging/mid-nano.svg" "${appdir}/usr/share/icons/hicolor/scalable/apps/${APP_ID}.svg"

appimage_work="${tmp}/appimage"
mkdir -p "${appimage_work}"

pushd "${appimage_work}" >/dev/null
export APPIMAGE_EXTRACT_AND_RUN=1
export VERSION="${DEB_VERSION}"

linuxdeploy \
  --appdir "${appdir}" \
  --executable "${appdir}/usr/bin/${APP_NAME}" \
  --executable "${appdir}/usr/bin/pfc-exp-cli" \
  --desktop-file "${appdir}/usr/share/applications/${APP_ID}.desktop" \
  --icon-file "${appdir}/usr/share/icons/hicolor/scalable/apps/${APP_ID}.svg" \
  --plugin qt \
  --output appimage

produced="$(ls -1 ./*.AppImage 2>/dev/null | head -n 1 || true)"
if [[ -z "${produced}" ]]; then
  echo "Failed to produce AppImage." >&2
  ls -la
  exit 1
fi

case "${APPIMAGE_ARCH}" in
  x86_64) appimage_suffix="linux-x86_64.AppImage" ;;
  aarch64) appimage_suffix="linux-arm64.AppImage" ;;
  *) echo "Unknown appimage arch: ${APPIMAGE_ARCH}" >&2; exit 2 ;;
esac

appimage_out="${OUT_DIR}/${APP_NAME}-${VERSION}-${appimage_suffix}"
mv "${produced}" "${appimage_out}"
chmod +x "${appimage_out}"
popd >/dev/null
