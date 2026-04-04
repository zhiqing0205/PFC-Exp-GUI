# 压缩上下文摘要（2026-04-04）

## 当前结论

- 项目是 Qt GUI（`MID Nano`）+ C++ CLI（`pfc-exp-cli`）
- GUI 通过命令行参数调用 CLI，不是把实验代码直接嵌在 GUI 里执行
- 统一入口：`src/experiment/main_cli.cpp`
- 统一模型：`misfit / cvd / elastic`

## 模型关系

- `misfit`：主仿真模型，老成员，早就在项目里
- `cvd`：新增并接入 GUI 的主仿真模型
- `elastic`：新增并接入 GUI 的后处理分析器

## 最关键判断

- GUI/CLI 层面三者是并列入口
- 算法/数据流层面：
  - `misfit` / `cvd` 负责产出场数据
  - `elastic` 负责读取已有 `phi_*.vti / con_*.vti` 做分析

## 实际依赖关系

- `misfit` 会持续写 `phi_*.vti / con_*.vti`
- `cvd` 会持续写 `phi_*.vti / con_*.vti`
- `elastic` 会读取这些 `.vti`
- 因此现在两条链都能走：
  - `cvd -> elastic`
  - `misfit -> elastic`

## `ref/` 目录状态

- `ref/` 是本地未跟踪目录，不是 Git 已跟踪源码
- 真正参与编译的是 `src/experiment/*.cpp`
- `ref/*.cpp` 更像参考原始实现

## 关键提交

- `1f61f46`：接入 `cvd` 和 `elastic`，并做 unified CLI + GUI 三模型页

## 关键文件

- 统一入口：`src/experiment/main_cli.cpp`
- 公共基础：`src/experiment/pfc_common.h`
- 老模型：`src/experiment/misfit_uniform_time_seed_fftw_early.cpp`
- 新模型：`src/experiment/cvd.cpp`
- 后处理：`src/experiment/elastic.cpp`
- GUI 模型页：`src/gui/MainWindow.cpp`

## 当前已落盘文档

- 完整分析：`docs/ref-models-analysis-2026-04-04.md`
- 本压缩摘要：`docs/context-compressed-2026-04-04.md`

## 后续推荐继续分析顺序

1. `misfit` 主循环
2. `cvd` 相对 `misfit` 的关键差异
3. `elastic` 的后处理链路
4. GUI 的传参和输出目录组织

---

## 2026-04-04 当日晚些时候的 GUI 调整结果

- `Experiment` 顶层语义改成 `Simulation`
- `Simulation` 页现在只保留：
  - `Misfit`
  - `CVD`
- 原先放在 `Experiment` 里的 `Elastic` 已移出，改为独立的 `Elastic Analysis` 页
- `Elastic Analysis` 现在显式区分：
  - `input dir`：必须包含 `phi_*.vti / con_*.vti`
  - `output dir`：默认与 input 相同，也可单独指定
- GUI 会扫描输入目录里的 VTI checkpoint，并尝试从 `run_config.txt` 自动回填参数
- Visualizer 增加了跳转到 `Elastic Analysis` 的入口
- CLI `elastic` 也已补强：
  - 真正使用 `--input-dir`
  - 会扫描实际存在的 checkpoint step，而不是盲目依赖手填 step 编号

### 现在的推荐逻辑流

- `Misfit / CVD`：前向仿真
- `Elastic Analysis`：读取已有 VTI 做应变/结构后处理
- `Visualizer`：查看 `Phimax_*.txt / strain_*.txt / q4q6_*.txt`

### 当前可用链路

- `Misfit -> Elastic Analysis -> Visualizer`
- `CVD -> Elastic Analysis -> Visualizer`
