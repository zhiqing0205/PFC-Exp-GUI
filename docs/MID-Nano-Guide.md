# MID Nano 使用手册（草案）

（本手册面向当前仓库版本：以 GitHub Releases / Tags 为准）

2026 年 2 月

---

## 声明

MID Nano 是一个用于配置实验参数、运行 `pfc-exp-cli` 并查看输出结果的桌面工具（Qt6 Widgets）。本软件主要面向 Windows 使用场景，同时提供 macOS / Linux 构建产物，便于跨平台演示与小规模测试。

本手册为“操作说明 + 截图占位”版本：文中的截图位置暂留空，后续可直接插入真实截图，并按建议标注进行标注说明。

---

## 目录

### 1. 前言
#### 1.1 软件简介
#### 1.2 版本信息与更新方式

### 2. 系统要求
#### 2.1 操作系统兼容性
#### 2.2 硬件建议
#### 2.3 依赖与网络要求

### 3. 安装与卸载
#### 3.1 下载与校验（Release）
#### 3.2 Windows / macOS / Linux 运行说明
#### 3.3 卸载方法

### 4. 用户界面概览
#### 4.1 主窗口布局
#### 4.2 标签页说明（Experiment / Visualizer）
#### 4.3 日志（Log）与进度条（Progress）

### 5. 核心功能操作指南
#### 5.1 实验（Experiment：固定参数）
#### 5.2 实验（Experiment：扫参组合）
#### 5.3 输出查看（Visualizer）
#### 5.4 输出目录与文件说明

### 6. 故障排除（FAQ）
#### 6.1 CLI 找不到 / 无法启动
#### 6.2 macOS 被系统终止（SIGKILL 9）
#### 6.3 Windows 缺 DLL（0xC0000135）
#### 6.4 常见参数与输出问题

### 7. 附录
#### 7.1 参数速查（CLI）
#### 7.2 快捷操作与交互说明

---

## 1. 前言

### 1.1 软件简介

MID Nano 的目标是把“改参数 → 跑一次 → 看结果”的流程做得更直观：

- 在 GUI 中配置实验参数（如 `u0 / con0 / sig / dt / dx / steps / mod / seed`）。
- `Experiment` 标签页支持固定参数运行，也支持扫参组合运行（Fixed / Range / List 三种模式）。
- 自动把每次运行写入输出目录，并生成 `params.json / run_config.txt / checkpoint_timestamps.txt / Phimax_*.txt` 等文件。
- 在 `Visualizer` 标签页中直接查看输出（表格 / 热图），并支持导出 PNG。

### 1.2 版本信息与更新方式

- 版本号以 Git 标签为准，例如 `0.4.1`。
- 推荐从 GitHub Releases 下载对应版本的压缩包；更新时替换整个目录即可。

> （截图占位：Release 页面/Tag 列表）
>
> 建议标注：
> ① Releases 入口；② 版本号（Tag）；③ 各平台下载产物；④ 校验信息（如有）。

---

## 2. 系统要求

### 2.1 操作系统兼容性

- Windows：Windows 10/11 64 位（主要使用目标）。
- macOS：支持运行（用于演示/测试）；首次运行可能需要处理隔离属性（见 FAQ）。
- Linux：支持运行（用于演示/测试）；当前发布版仅提供 x86_64 产物，arm64 暂不提供。

### 2.2 硬件建议

本项目以“小规模实验/演示”为主，不要求高端硬件。建议：

| 项目 | 最低配置 | 推荐配置 |
|---|---:|---:|
| CPU | 双核 | 四核及以上 |
| 内存 | 4 GB | 8–16 GB |
| 磁盘 | 1 GB 可用空间 | 5 GB 以上（用于保存输出） |

### 2.3 依赖与网络要求

- 使用发布版：无需额外依赖，解压即可运行。
- 从源码构建：需要 Qt6、CMake、编译器、FFTW3（Windows 推荐 vcpkg）。
- 网络：仅在下载发布版、拉取依赖或使用 CI 构建时需要。

---

## 3. 安装与卸载

### 3.1 下载与校验（Release）

从 GitHub Releases 下载对应平台产物（文件名以版本号为准），例如：

```
MID Nano-0.4.1-win-x64.exe
MID Nano-0.4.1-win-x64.zip
MID Nano-0.4.1-mac-arm64.dmg
MID Nano-0.4.1-mac-arm64.zip
MID Nano-0.4.1-mac-x64.dmg
MID Nano-0.4.1-mac-x64.zip
MID Nano-0.4.1-linux-amd64.deb
MID Nano-0.4.1-linux-x86_64.AppImage
```

> （截图占位：下载产物列表）
>
> 建议标注：
> ① 各平台产物（zip/dmg/deb/AppImage）；② 同一版本号（Tag）对应多种产物；③ 校验信息（如有）。

### 3.2 Windows / macOS / Linux 运行说明

**Windows**

推荐两种方式（二选一）：

1) 安装版（`*-win-x64.exe`）

1. 双击安装包，按提示安装。
2. 从开始菜单启动 `MID Nano`。

2) 便携版（`*-win-x64.zip`）

1. 解压 zip 到任意目录（建议放在非系统目录，如 `D:\Apps\MID Nano\`）。
2. 双击运行 `MID Nano.exe`。

无论哪种方式，请确保 `pfc-exp-cli.exe` 与 GUI 在同一目录（发布版已包含）。

> （截图占位：Windows 解压后的目录结构）
>
> 建议标注：① `MID Nano.exe`；② `pfc-exp-cli.exe`；③ 依赖 DLL（如有）。

**macOS**

推荐两种方式（二选一）：

1) DMG（`*.dmg`）

1. 双击打开 DMG。
2. 将 `MID Nano.app` 拖拽到 `/Applications`。

2) ZIP（`*.zip`）

1. 解压 zip。
2. 将 `MID Nano.app` 放到 `/Applications`（或任意目录）运行。

若首次运行被阻止，可尝试执行：

- `xattr -dr com.apple.quarantine "/Applications/MID Nano.app"`

> （截图占位：macOS Finder 中的 app）
>
> 建议标注：① app 名称；② 右键打开提示（如出现）；③ 终端命令示例。

**Linux**

推荐两种方式（二选一）：

1) DEB（`*.deb`）

- Debian/Ubuntu：
  - `sudo apt-get install -y "./MID Nano-0.4.1-linux-amd64.deb"`

2) AppImage（`*.AppImage`）

- `chmod +x "MID Nano-*.AppImage"`
- `./"MID Nano-0.4.1-linux-x86_64.AppImage"`

> （截图占位：Linux 终端运行示例）
>
> 建议标注：① 安装/运行命令；② 权限处理；③ 启动后的主界面。

### 3.3 卸载方法

不同平台卸载方式如下：

- Windows 安装版：运行开始菜单中的卸载项，或执行安装目录内的 `Uninstall.exe`。
- Windows 便携版：直接删除解压目录即可。
- macOS：删除 `MID Nano.app`。
- Linux（DEB）：`sudo apt remove mid-nano`（或 `sudo dpkg -r mid-nano`）。
- Linux（AppImage）：删除 `.AppImage` 文件即可。

输出目录默认在“文档/MID Nano/outputs”，如需清理请手动删除输出目录。

---

## 4. 用户界面概览

### 4.1 主窗口布局

主窗口由三部分组成：顶部为快捷按钮区（Upload License / About），中间为功能标签页，下方为日志区域（可拖动分隔条调整高度）。

> （截图占位：主窗口全景）
>
> 建议标注：
> ① 顶部快捷按钮；② 标签页；③ 参数区；④ 输出区；⑤ 进度条；⑥ Log（可调高度）。

### 4.2 标签页说明

**Experiment（实验/扫参）**：每个参数可选 Fixed / Range / List；当全部为 Fixed 时等价于单次实验，当任意参数使用 Range/List 且产生多值时会自动组合生成多组任务依次运行。

**Visualizer（结果查看）**：查看 `*.txt` 输出；对 `Phimax_*.txt` 提供热图渲染与 PNG 导出。

其余标签页（Manufacturing / Transformation / Stress–Strain）为后续规划入口，目前为占位说明。

### 4.3 日志（Log）与进度条（Progress）

- `Log` 会打印每次运行的输出目录、执行命令、CLI 输出与错误提示。
- Experiment 会显示“当前 run 的 step 进度”；当总运行数 > 1 时，还会显示“整体 run 进度”（按队列完成比例更新）。

### 4.4 关于（About）与 License

- `About`：显示项目名称、版本号、作者、仓库链接（GitHub/Gitee），以及依赖信息（Qt/FFTW）与构建信息（Build type / Git commit / build time / 运行环境信息）。
- `Upload License…`：选择一个 license 文件并保存到本机应用数据目录（仅本地保存，不会上传网络）。该文件用于后续功能（如授权校验、功能解锁等）的对接。

> （截图占位：About 弹窗）
>
> 建议标注：① 项目名称；② 版本号；③ 仓库链接；④ 作者信息。

---

## 5. 核心功能操作指南

### 5.1 实验（Experiment：固定参数）

1. 在 `Parameters` 中设置参数（建议先从小步数开始验证流程）。
2. 在 `Output` 中选择输出根目录，并可选填写 `Experiment folder name`。
3. 点击 `Run` 开始运行；需要中止时点击 `Stop`。
4. 点击 `Open Output` 可直接打开当前输出目录。

> （截图占位：Experiment（固定参数）参数区）
>
> 建议标注：① 关键参数（u0/con0/steps/mod/seed）；② Run/Stop；③ Step 进度条；④ Preview（1 run）。

### 5.2 实验（Experiment：扫参组合）

Experiment 的“参数模式”支持三种：

- **Fixed**：固定单个值。
- **Range**：通过 Start/End/Step 生成序列（右侧会显示生成的列表预览）。
- **List**：直接输入逗号分隔的数值列表。

GUI 会对所有参数的序列做笛卡尔积组合，生成任务队列并依次运行。

> （截图占位：Experiment（扫参组合）参数表）
>
> 建议标注：① Param 列；② Mode 下拉框；③ Range/List 输入与预览；④ Preview（总运行数）；⑤ 整体 run 进度条。

### 5.3 输出查看（Visualizer）

1. 打开 `Visualizer` 标签页后，会默认切到“最近一次运行”的输出目录并加载 `*.txt`。
2. 左侧选择文件：
   - `run_config.txt`：以 `key/value` 表格展示参数；
   - `checkpoint_timestamps.txt`：展示时间戳、可读时间、累计耗时与间隔；
   - `Phimax_*.txt` 等：启用 `Plot` 热图，支持缩放/拖动与导出 PNG。

> （截图占位：Visualizer（Table））
>
> 建议标注：① 文件列表；② Table 列名；③ 特殊文件的展示效果（run_config、timestamps）。
>
> （截图占位：Visualizer（Plot 热图））
>
> 建议标注：① Value 选择；② 平滑参数；③ Fit；④ Export PNG；⑤ 颜色分布区域。

### 5.4 输出目录与文件说明

每次运行会生成一个 run 目录，通常包含：

- `params.json`：GUI 写入的参数快照（用于复现）。
- `run_config.txt`：CLI 写入的运行配置（用于核对）。
- `checkpoint_timestamps.txt`：检查点时间戳（用于耗时分析）。
- `Phimax_*.txt` / `Phimaxq4q6_*.txt`：快照数据（用于热图可视化）。
- 其他程序输出文件（如 `*.vti / *.pvti` 等，视 CLI 输出而定）。

---

## 6. 故障排除（FAQ）

### 6.1 CLI 找不到 / 无法启动

现象：GUI 弹窗提示找不到 `pfc-exp-cli`，或点击 Run 后立即失败。

排查：

- 确认 `pfc-exp-cli`（Windows 为 `pfc-exp-cli.exe`）与 GUI 在同一目录。
- 查看 `Log` 是否有 “CLI not found” 或缺 DLL 的提示。

### 6.2 macOS 被系统终止（SIGKILL 9）

现象：日志显示 `exitStatus=CrashExit exitCode=9`。

处理：

- 尝试执行：`xattr -dr com.apple.quarantine "/Applications/MID Nano.app"`
- 若仍异常，检查 `Contents/Frameworks` 依赖是否完整。

### 6.3 Windows 缺 DLL（0xC0000135）

现象：日志显示 `0xC0000135`。

处理：

- 该错误通常表示依赖 DLL 缺失（例如 FFTW）。
- 请确认使用的是发布版产物，或确保 `*.dll` 与 `pfc-exp-cli.exe` 位于同目录。

### 6.4 常见参数与输出问题

- 建议先用较小 `steps` 做流程验证，再逐步增加规模。
- 若输出目录为空或权限不足，请更换到有写入权限的路径（如用户文档目录）。

---

## 7. 附录

### 7.1 参数速查（CLI）

常用参数（以 `--help` 输出为准）：

- `--u0`：密度/序参量均值
- `--con0`：浓度均值
- `--sig`：核函数参数
- `--dt` / `--dx`：时间/空间步长
- `--steps`：迭代步数
- `--mod`：输出/分析间隔
- `--seed`：随机种子
- `--outdir`：输出目录

### 7.2 快捷操作与交互说明

- `Log` 区域高度：拖动中间分隔条调整。
- `Plot` 热图缩放：滚轮缩放。
- `Plot` 热图平移：鼠标左键按住拖动。
- 目录快速定位：`Open Output` 打开当前输出目录。
- 快捷入口：顶部 `About` 查看版本与链接；`Upload License…` 保存 license 文件。
