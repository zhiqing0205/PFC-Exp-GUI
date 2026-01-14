# Qt GUI 实验参数配置与批量运行（misfit_uniform_time_seed_fftw_early）

本项目包含两部分：
- `misfit_sim`：MPI + FFTW 的 C++ 实验程序（由 `misfit_uniform_time_seed_fftw_early.cpp` 构建），已支持命令行参数传入与 `--outdir` 输出目录。
- `qt_gui_client`：Qt6(QML) 图形界面，用于配置参数、设置范围/步长做组合扫描，并按队列逐个启动实验。

另附：`misfit_uniform_time_seed_fftw_early_overview.md`（原程序概览）。

## 1. Ubuntu 24.04 安装依赖

在 Ubuntu 24.04 上可直接安装：

```bash
./scripts/install_ubuntu_24.04.sh
```

（脚本使用 `apt-get`，需要 sudo 权限与网络访问。）

## 2. 编译（推荐 CMake + Ninja）

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

编译产物在 `build/bin/`：
- `build/bin/misfit_sim`
- `build/bin/qt_gui_client`

如果只想编译 GUI：

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_SIM=OFF
cmake --build build
```

## 3. 运行 GUI

```bash
./build/bin/qt_gui_client
```

GUI 中：
- `实验程序 (misfit_sim)` 可以留空（默认使用与 GUI 同目录的 `misfit_sim`）。
- 输出会在 “父目录/会话前缀_时间戳/每个 run 子目录” 下生成。
- 扫描参数启用后会对多个参数做组合（笛卡尔积），按队列顺序依次运行。

每个 run 的目录下会生成：
- `run_config.txt`（实验程序写出）
- `gui_process.log`（GUI 捕获的标准输出/错误）

## 4. 运行命令行（单次实验）

无需 GUI 也可直接运行（示例：单进程）：

```bash
./build/bin/misfit_sim --outdir runs/manual_01 --u0 0.05 --con0 0.2 --steps 5002 --dt 0.05 --dx 0.125
```

使用 MPI（示例：4 进程）：

```bash
mpirun -np 4 ./build/bin/misfit_sim --outdir runs/manual_mpi_01 --u0 0.05 --con0 0.2 --steps 5002
```

查看支持的参数：

```bash
./build/bin/misfit_sim --help
```

## 5. Windows / macOS / Linux 构建说明

### Linux
- 直接参考 “Ubuntu 24.04 安装依赖” + CMake 编译即可。
- 生成 `.deb`（Ubuntu/Debian）：
  ```bash
  ./scripts/package_deb.sh
  ```

### macOS（Homebrew）
```bash
brew install cmake ninja qt open-mpi fftw pkg-config
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

生成 `.dmg`（需 Qt 自带 `macdeployqt` 在 PATH 中）：
```bash
./scripts/package_macos.sh
```

> 注：`misfit_sim` 依赖 `fftw3_mpi`；如果你本地 FFTW 未启用 MPI，需要用支持 MPI 的 FFTW 包或自行编译 FFTW。

### Windows
- GUI 推荐使用 Qt 官方安装器安装 Qt6 + CMake（或 Visual Studio + CMake）。
- `misfit_sim`（MPI+FFTW）在 Windows 原生环境下依赖链较重，建议：
  - 在 WSL2/Ubuntu 内编译运行实验程序；
  - Windows 上仅构建 GUI：`-DBUILD_SIM=OFF`，并在 GUI 中手动选择 WSL 侧或远端的实验程序（或改为生成命令行脚本使用）。

生成可分发的 Windows 包（需 Qt 自带 `windeployqt.exe` 在 PATH 中）：
```powershell
.\scripts\package_windows.ps1
```

## 6. 常见注意事项
- 该实验程序网格很大（默认 1920×1920×1）且数组很多，单次运行可能占用较高内存；建议在具备足够内存/多核的环境运行。
- 批量扫描组合数过大时，GUI 会限制（默认最多 5000 组）以避免误触产生极长队列，可通过缩小范围或增大步长控制规模。

## 7. 推一次自动产出安装包（GitHub Actions）

仓库内已提供工作流：`.github/workflows/build-and-package.yml`，触发条件：
- `push` 到 `main/master`：自动构建并生成 artifacts（在 Actions 页面下载）
- `push` tag `v*`（例如 `v0.1.0`）：自动创建 GitHub Release，并把安装包作为 Release assets 上传

产物（默认策略）：
- Linux：生成 `.deb`（包含 `qt_gui_client` + `misfit_sim`）
- Windows：生成可直接解压运行的 `.zip`（GUI 客户端；实验程序建议在 WSL/远端/同机 Linux 构建后在 GUI 里选择）
- macOS：生成 `.dmg`（GUI 客户端）

你需要做的步骤：
1) 把当前目录初始化为 git 仓库并推到 GitHub（如果还没做）
2) 确认仓库已开启 Actions
3) 正常 `git push`：到 GitHub → Actions → 选择本工作流 → 下载 artifacts
4) 想要自动生成 Release：打 tag 并 push
   ```bash
   git tag v0.1.0
   git push origin v0.1.0
   ```
