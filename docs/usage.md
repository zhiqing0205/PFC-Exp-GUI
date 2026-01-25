# 使用说明

## 1. 单次实验（Single Run）

1. 选择 `Base output dir`（输出根目录）。
2. 选择 `Run folder name`（可留空，自动使用时间戳）。
3. 设置关键参数（示例）：
   - `u0`：密度/序参量均值
   - `con0`：浓度均值
   - `steps`：迭代步数
   - `seed`：随机数种子（用于初始晶粒/噪声等随机过程，便于复现）
4. 点击 `Run`。

输出将写入 `Base output dir/Run folder name/`。

## 2. 批量扫参（Batch Sweep）

Batch Sweep 支持对 `u0 / con0 / steps` 设置区间与步长，并按笛卡尔积组合顺序运行：

- 取消勾选 `Sweep`：使用 `Single` 的固定值。
- 勾选 `Sweep`：使用 `Start/End/Step` 生成序列。

例如：

- `u0`: 0.00 → 0.20，步长 0.05（5 组）
- `con0`: 0.10 → 0.30，步长 0.05（5 组）
- `steps`: 100 → 300，步长 100（3 组）

总运行数 = 5 × 5 × 3 = 75。

## 3. 输出内容

每个 run 目录中至少包含：

- `params.json`：GUI 写入的参数快照（便于复现实验）
- `run_config.txt`：CLI 写入的运行参数
- `phi_*.vti / con_*.vti / *.pvti` 等（实验程序输出）

## 4. TXT 可视化（TXT Visualizer）

GUI 提供 `TXT Visualizer` 标签页，用于快速查看输出目录中的 `*.txt`（例如 `energy.txt`、`phimax.txt`、`checkpoint_timestamps.txt`）：

1. 打开 `TXT Visualizer`。
2. 在 `Run directory` 选择某个 run 的输出目录。
3. 左侧选择 `*.txt` 文件；右侧 `Plot` 可选择要绘制的列（默认使用第 0 列为 x）。

## 5. macOS 常见问题

- 如果运行时日志显示 `exitCode=9 (SIGKILL)`，通常是 Gatekeeper/隔离属性（quarantine）或依赖 dylib 缺失导致。可以尝试：
  - `xattr -dr com.apple.quarantine /Applications/PFC-Exp-GUI.app`
