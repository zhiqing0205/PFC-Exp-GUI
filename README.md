# MID Nano

使用 Qt（Qt6 Widgets）制作的实验参数配置 GUI，用来生成/组合实验参数并调用 C++ 实验程序运行（`pfc-exp-cli`）。

## 组成

- `MID Nano`：图形界面（单次实验 + 批量扫参）
- `pfc-exp-cli`：实验程序命令行版本（参数通过命令行传入，输出写入指定目录）

## 功能

- 单次实验：配置 `u0`（密度/序参量均值）、`con0`（浓度）、`steps`（迭代步数）等参数并运行。
- 批量实验：对 `u0 / con0 / steps` 设置区间与步长，按笛卡尔积组合后顺序运行多组实验。
- 每次运行都会在输出目录写入 `params.json`（GUI 生成）与 `run_config.txt`（CLI 写入）。

## 构建（CMake）

依赖：

- CMake ≥ 3.20
- C++17 编译器
- Qt6（Widgets）
- FFTW3（`fftw3`）

### Linux（Ubuntu/Debian）

```bash
sudo apt-get update
sudo apt-get install -y qt6-base-dev qt6-base-dev-tools libfftw3-dev ninja-build cmake g++

cmake -S . -B build -G Ninja
cmake --build build
```

产物位于 `build/bin/`：

- `build/bin/MID Nano`
- `build/bin/pfc-exp-cli`

### Windows（建议使用 Release）

你可以直接使用 GitHub Actions 构建好的发布产物（见下文）。如果要本地构建：

1. 安装 Qt6（MSVC 64-bit 版本）。
2. 安装 FFTW3（推荐 vcpkg）。

示例（vcpkg）：

```powershell
vcpkg install fftw3:x64-windows
cmake -S . -B build -A x64 `
  -DCMAKE_TOOLCHAIN_FILE=C:\path\to\vcpkg\scripts\buildsystems\vcpkg.cmake `
  -DCMAKE_PREFIX_PATH="C:\Qt\6.x.x\msvc2019_64"
cmake --build build --config Release
```

## 使用

### GUI

运行 `MID Nano`（确保同目录下存在 `pfc-exp-cli`），在界面中填写参数与输出目录：

- Single Run：单次运行
- Batch Sweep：扫参批量运行（`u0 / con0 / steps`）

### CLI

```bash
./pfc-exp-cli --help
./pfc-exp-cli --steps 200 --mod 25 --seed 20200604 --u0 0.05 --con0 0.2 --outdir ./out/run1
```

## 网格尺寸（编译期）

当前实验代码使用编译期网格（默认 `256×256×1`）。可通过 CMake 缓存变量修改：

```bash
cmake -S . -B build -G Ninja -DPFC_GRID_L=512 -DPFC_GRID_M=512 -DPFC_GRID_N=1
cmake --build build
```

## GitHub Actions（自动构建与发布）

仓库内提供工作流：在 push tag（例如 `0.1.10`，也兼容 `v0.1.10`）时自动构建 Windows/macOS/Linux 并上传 release 产物（可适当简化为 zip/tar.gz）。

更多实验程序原理说明见 `docs/misfit_uniform_time_seed_fftw_early_overview.md`。
