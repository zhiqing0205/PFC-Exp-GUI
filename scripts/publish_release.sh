#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${GH_TOKEN:-}" ]]; then
  echo "GH_TOKEN is required." >&2
  exit 1
fi

repo="${GITHUB_REPOSITORY:?}"
ref="${GITHUB_REF:?}"
ref_name="${GITHUB_REF_NAME:?}"
sha="${GITHUB_SHA:?}"

mode="${RELEASE_MODE:-}"
if [[ -z "${mode:-}" ]]; then
  if [[ "$ref" == refs/tags/* ]]; then
    mode="tag"
  else
    mode="nightly"
  fi
fi

assets_dir="release_assets"
rm -rf "$assets_dir"
mkdir -p "$assets_dir"

package_files=()
while IFS= read -r f; do
  package_files+=("$f")
done < <(find dist -type f \( -name '*.deb' -o -name '*.zip' -o -name '*.dmg' \) -print)

if [[ "${#package_files[@]}" -eq 0 ]]; then
  echo "No package files found under dist/." >&2
  exit 1
fi

base_prefix="qt-experiment-runner"
if [[ "$mode" == "nightly" ]]; then
  base_prefix="qt-experiment-runner-nightly"
fi

for f in "${package_files[@]}"; do
  base="$(basename "$f")"
  case "$f" in
    */linux-deb/*)
      if [[ "$base" == *.deb ]]; then
        stem="${base%.deb}"
        arch="${stem##*_}"
        cp -f "$f" "$assets_dir/${base_prefix}-linux-${arch}.deb"
      else
        cp -f "$f" "$assets_dir/$base"
      fi
      ;;
    */windows-zip/*)
      if [[ "$base" == *.zip ]]; then
        stem="${base%.zip}"
        arch="${stem##*_}"
        if [[ "$arch" == "windows" || -z "${arch:-}" ]]; then
          arch="x64"
        fi
        cp -f "$f" "$assets_dir/${base_prefix}-windows-${arch}.zip"
      else
        cp -f "$f" "$assets_dir/$base"
      fi
      ;;
    */macos-dmg/*)
      if [[ "$base" == *.dmg ]]; then
        stem="${base%.dmg}"
        arch="${stem##*_}"
        if [[ "$arch" == "macos" || -z "${arch:-}" ]]; then
          arch="arm64"
        fi
        cp -f "$f" "$assets_dir/${base_prefix}-macos-${arch}.dmg"
      else
        cp -f "$f" "$assets_dir/$base"
      fi
      ;;
    *)
      cp -f "$f" "$assets_dir/$base"
      ;;
  esac
done

(
  cd "$assets_dir"
  shopt -s nullglob
  checksum_files=( *.deb *.zip *.dmg )
  if [[ "${#checksum_files[@]}" -eq 0 ]]; then
    echo "No package files to checksum under ${assets_dir}/." >&2
    exit 1
  fi
  sha256sum "${checksum_files[@]}" > SHA256SUMS.txt
)

if [[ "$mode" == "tag" ]]; then
  tag="$ref_name"
  if gh release view "$tag" >/dev/null 2>&1; then
    gh release upload "$tag" "$assets_dir"/* --clobber
  else
    gh release create "$tag" "$assets_dir"/* --title "$tag" --generate-notes
  fi
else
  tag="nightly"
  gh release delete "$tag" --yes --cleanup-tag >/dev/null 2>&1 || true
  gh release create "$tag" "$assets_dir"/* --title "Nightly" --notes "Automated build for ${sha}" --prerelease --target "$sha"
fi

export RELEASE_TAG="$tag"

python3 - <<'PY'
import glob
import hashlib
import json
import os
import subprocess

repo = os.environ["GITHUB_REPOSITORY"]
tag = os.environ["RELEASE_TAG"]

release_raw = subprocess.check_output(["gh", "api", f"repos/{repo}/releases/tags/{tag}"])
release = json.loads(release_raw.decode("utf-8"))
assets = {asset.get("name"): asset.get("id") for asset in release.get("assets", [])}

for path in sorted(glob.glob("release_assets/*")):
  if not os.path.isfile(path):
    continue
  name = os.path.basename(path)
  if name == "SHA256SUMS.txt":
    continue

  hasher = hashlib.sha256()
  with open(path, "rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      hasher.update(chunk)
  digest = hasher.hexdigest()

  asset_id = assets.get(name)
  if not asset_id:
    raise SystemExit(f"Asset not found in release: {name}")

  subprocess.check_call([
    "gh",
    "api",
    "--method",
    "PATCH",
    f"repos/{repo}/releases/assets/{asset_id}",
    "-f",
    f"label=sha256:{digest}",
  ])
PY
