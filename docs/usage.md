# 使用说明

## 1. 实验（Experiment）

Experiment 采用“逐参数模式”，每个参数都可以选择一种生成方式，并按笛卡尔积组合顺序运行多组实验：

- `Fixed`：固定一个值
- `Range`：通过 `Start / End / Step` 生成序列（右侧会显示生成后的列表预览）
- `List`：直接输入逗号分隔的数值列表（例如 `0.05,0.10,0.15`）

运行步骤：

1. 选择 `Base output dir`（输出根目录）。
2. 选择 `Experiment folder name`（可留空，自动使用时间戳）。
3. 配置参数（可全固定，也可扫参组合）。
4. 点击 `Run`。

输出目录：

- 当总运行数 = 1：输出写入 `Base output dir/Experiment folder name/`。
- 当总运行数 > 1：输出写入 `Base output dir/Experiment folder name/run_0001/` 等子目录。

示例：

- `u0`: Range 0.00 → 0.20，步长 0.05（5 组）
- `con0`: List 0.10,0.20,0.30（3 组）
- `steps`: Fixed 200（1 组）

总运行数 = 5 × 3 × 1 = 15（界面会显示预览与总数）。

## 2. 输出内容

每个 run 目录中至少包含：

- `params.json`：GUI 写入的参数快照（便于复现实验）
- `run_config.txt`：CLI 写入的运行参数
- `phi_*.vti / con_*.vti / *.pvti` 等（实验程序输出）

## 3. 可视化（Visualizer）

GUI 提供 `Visualizer` 标签页，用于快速查看输出目录中的 `*.txt`。打开该标签页时，会默认切换到当前运行的输出目录（Current Output），左侧自动加载该目录下的 `*.txt` 文件列表。

- 普通 `*.txt`（如能量、时间戳等日志类输出）：默认使用 `Table` 查看原始数值列。
- `Phimax_*.txt` / `Phimaxq4q6_*.txt` 等二维快照：会启用 `Plot`（热图）渲染，可滚轮缩放、鼠标拖动，支持 `Fit` 与 `Export PNG…` 导出图片。

## 4. macOS 常见问题

- 如果运行时日志显示 `exitCode=9 (SIGKILL)`，通常是 Gatekeeper/隔离属性（quarantine）或依赖 dylib 缺失导致。可以尝试：
  - `xattr -dr com.apple.quarantine "/Applications/MID Nano.app"`
