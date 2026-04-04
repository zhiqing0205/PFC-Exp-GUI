# ref 模型关系分析（2026-04-04）

## 结论摘要

这三个文件不是完全独立、互不相关的三个程序；它们更像是同一套 PFC/MPI/FFTW 代码骨架的三个分支：

- `misfit`：主仿真模型，做失配应变相关的 PFC 演化
- `cvd`：另一个主仿真模型，做 CVD / 两相生长类演化
- `elastic`：后处理分析器，读取前面仿真产出的 `phi_*.vti / con_*.vti`，计算 `Phimax / strain / q4q6`

因此：

- 在 GUI / CLI 入口层面：它们被当成并列模型
- 在算法 / 数据流层面：`elastic` 更像下游分析程序，不完全和 `misfit/cvd` 同级

---

## 1. 项目整体结构

这个项目是一个 Qt GUI + C++ CLI 实验程序：

- Qt GUI：`MID Nano`
- CLI：`pfc-exp-cli`
- GUI 不直接在窗口内跑实验，而是拼接命令行参数并调用 CLI，再做结果可视化

关键位置：

- `README.md`
- `CMakeLists.txt`
- `src/experiment/main_cli.cpp`
- `src/gui/MainWindow.cpp`

---

## 2. 最近提交与实际接入情况

关键提交：

- `1f61f46`  
  `Add multi-model support: CVD and Elastic models with unified CLI`

这个提交新增/接入了：

- `src/experiment/cvd.cpp`
- `src/experiment/elastic.cpp`
- `src/experiment/main_cli.cpp`
- `src/experiment/pfc_common.h`

同时把 GUI 扩成了 3 个模型子页：

- Misfit
- CVD
- Elastic

### 关于 misfit

`misfit` 并不是这次新加的；它在多模型提交之前就已经在项目中。

### 关于 `ref/`

当前 `ref/` 目录是本地未跟踪目录，不是 Git 已跟踪源码目录：

- `git status --short` 显示 `?? ref/`
- `git ls-files` 中没有 `ref/*.cpp`

也就是说：

- 真正参与编译的是 `src/experiment/*.cpp`
- `ref/*.cpp` 更像参考原始版本 / 本地备份源码

---

## 3. GUI / CLI 里的并列关系

### 统一 CLI 入口

`src/experiment/main_cli.cpp` 使用：

- `--model misfit`
- `--model cvd`
- `--model elastic`

进行分发。

### GUI

`src/gui/MainWindow.cpp` 中也明确把 Experiment 页面拆成 3 个子标签页：

- Misfit
- CVD
- Elastic

所以从产品/UI 视角，这三者是并列入口。

---

## 4. 三个模型分别做什么

## 4.1 `ref/misfit.cpp`

这是主仿真程序。

特征：

- 使用 `phi` 与 `con` 两个场做 PFC 演化
- 有完整的 FFTW/MPI 时步推进
- 有噪声项 `Gauss(...)`
- 会写 `phi_*.vti / con_*.vti`
- 也带峰值/结构分析函数，比如 `phimax_atmposi`、`cryst_order`

可以理解为：

> 负责生成晶体场随时间演化的主求解器

---

## 4.2 `ref/cvd.cpp`

这也是主仿真程序，而且与 `misfit` 代码血缘很近。

共享的大骨架包括：

- MPI + FFTW
- `phi/con` 双场
- `C2AK / C2BK`
- `UPDATEINN3 / UPDATEINFMIX / UPDATEINFCMIX`
- VTK 输出

但有几个重要差异：

- 初始化函数换成 `INITIAL_C2(...)`
- 参数组不同
- 主方程主路径走 `MAINEQUATIONC2MIX(...)`
- 噪声路径大部分注释掉了
- 周期性输出 `phi_*.vti / con_*.vti`

因此它不是 `misfit` 的下游，而是：

> 同一套求解框架下的另一个主仿真模型

---

## 4.3 `ref/elastic.cpp`

这个不一样。

虽然代码表面上仍保留了不少与主求解器相似的函数声明和全局参数，但主流程不再做场演化，而是：

- 按步数循环
- 每隔 `mod`
  - 读取 `phi_<step>_<rank>.vti`
  - 读取 `con_<step>_<rank>.vti`
  - 提取峰值位置
  - 输出 `Phimax_<step>.txt`
  - 做应变分析，输出 `strain_<step>.txt`
  - 还定义了晶体序参数分析 `q4q6_*.txt`

所以它本质上是：

> 读取已有场数据并做后处理分析的程序

不是主仿真器。

---

## 5. 在 `src/experiment` 中的改造版情况

## 5.1 `src/experiment/cvd.cpp`

文件头已经明确标注：

- `CVD (Controlled Vapor Deposition) PFC simulation model`
- `Adapted from ref/cvd.cpp for unified CLI framework`

特点：

- 被改成 `run_cvd(...)`
- 接入统一参数解析
- 支持 `--outdir`
- 启动后 `chdir` 到输出目录
- 主循环中持续写 `phi_*.vti / con_*.vti`

因此它是可直接由 GUI 调用的主仿真模型。

---

## 5.2 `src/experiment/elastic.cpp`

文件头已经明确标注：

- `Elastic strain analysis PFC model`
- `Post-processing: reads pre-computed phi/con fields, computes strain tensors and crystallographic order`

主流程特征：

- `run_elastic(...)`
- 循环里按步数去读 `phi_<step>_<rank>.vti / con_<step>_<rank>.vti`
- 调 `phimax_atmposi(...)`
- 输出 `Phimax_*.txt`
- 调 `elastic3D(...)`
- 输出 `strain_*.txt`

这里的关键结论是：

> `elastic` 在当前代码中是后处理分析器

---

## 5.3 `src/experiment/misfit_uniform_time_seed_fftw_early.cpp`

这是项目原有主模型，当前统一入口中对应 `run_misfit(...)`。

特点：

- 仍然是主仿真程序
- 会做时步推进
- 当前 GUI 接入版里主要输出 `Phimax_*.txt`
- 初始化时会写一次 `phi/con` VTI
- 但循环里原本持续写 `phi/con` VTI 的代码现在被注释掉了

因此：

- 它仍然是主求解器
- 但它当前不像 `cvd` 那样持续产出供 `elastic` 消费的 VTI checkpoint

---

## 6. 三者的“关系”怎么理解最准确

## 6.1 从代码血缘看

三者明显属于同一家族：

- 都是 MPI + FFTW + PFC 双场
- 共享大量函数命名和全局参数风格
- `cvd`、`elastic` 看起来就是在老 PFC 代码上改出来的分支

所以不是三个毫无关系的程序。

---

## 6.2 从 GUI / CLI 入口看

它们被统一包装成并列模型：

- `misfit`
- `cvd`
- `elastic`

这只是产品层面的并列。

---

## 6.3 从数据流看

真正的依赖关系更像这样：

- `misfit` / `cvd`：生产场数据
- `elastic`：读取已有场数据并做分析

---

## 7. 基于上述关系，GUI 已做的调整（2026-04-04）

为了让产品结构和真实数据流更一致，已经按下面的思路调整 GUI / CLI：

### 7.1 GUI 层

- 顶层 `Experiment` 语义改成了 `Simulation`
- `Simulation` 里只保留两个主模型页：
  - `Misfit`
  - `CVD`
- `Elastic` 不再作为仿真模型页混在一起，而是迁移为独立的 `Elastic Analysis` 页

这样用户在界面层面也能直接看出：

- `Misfit/CVD` 是“生成数据”的
- `Elastic` 是“消费已有数据做分析”的

### 7.2 Elastic Analysis 页现在的逻辑

新的 `Elastic Analysis` 页显式拆开了输入输出：

- `input dir`
  - 要求包含 `phi_*.vti / con_*.vti`
- `output dir`
  - 默认直接写回 input 目录
  - 也支持单独指定

同时 GUI 会：

- 扫描 input 目录里真实存在的 VTI checkpoint
- 尝试从 `run_config.txt` 自动回填 `u0 / con0 / sig / dt / dx / steps / mod / seed`
- 如果目录里根本没有成对的 `phi/con .vti`，会直接给出提示

### 7.3 CLI elastic 的补强

为了让 GUI 上的输入输出逻辑真正成立，也同步修了 `src/experiment/elastic.cpp`：

- `--input-dir` 现在真的参与读文件，不再只是“解析了但没用”
- 读入前会把 input 目录解析成绝对路径，避免 `--outdir` 切换 cwd 后找错文件
- 会优先扫描 input 目录里实际存在的 checkpoint step，再逐个分析
  - 这样就不会再因为仿真侧 checkpoint 编号风格不同而读不到文件
- 读文件失败时会返回非零退出码，GUI 能正确判定失败

### 7.4 当前更合理的使用链路

现在更推荐把整个产品理解成：

1. `Simulation / Misfit`
2. `Simulation / CVD`
3. `Elastic Analysis`
4. `Visualizer`

其中当前最通顺的默认链路依然是：

- `CVD -> Elastic Analysis -> Visualizer`

因为当前 GUI 接入版 `misfit` 仍没有像 `cvd` 那样稳定地持续输出供 elastic 消费的 VTI checkpoint。

也就是说：

> `elastic` 是下游分析器

---

## 6.4 最准确的关系描述

最准确的描述应为：

> `misfit` 和 `cvd` 是同一套 PFC 求解框架下的两个并列主模型；`elastic` 不是同类主求解器，而是读取前者输出场文件做应变/结构分析的后处理器。  
> 在 GUI 中它被“产品化地并列”了，但在算法依赖关系上它其实是下游。

---

## 7. 当前实现里的一个重要实际问题

`elastic` 虽然支持 `--input-dir` 参数，但当前实现存在一个实际情况：

- GUI 构造命令行时没有传 `--input-dir`
- `PFC_SetupOutputDir(...)` 会直接切换当前工作目录到 `--outdir`
- `elastic` 读文件时直接 `fopen("phi_<step>_<rank>.vti")`

这意味着当前 `elastic` 的真实运行假设是：

> 它运行所在的输出目录里，已经存在现成的 `phi/con` VTI 文件

否则会直接打不开。

---

## 8. `misfit -> elastic` 与 `cvd -> elastic` 的现实可用性

### `cvd -> elastic`

链路基本是通的，因为 `cvd` 会持续写 `phi_*.vti / con_*.vti`。

### 当前 GUI 版 `misfit -> elastic`

链路不太顺，因为当前 `misfit` 在循环中把后续 checkpoint VTI 输出注释掉了，只保留了初始 VTI 与 `Phimax_*.txt` 相关输出。

所以现实上更像：

- CVD 产出 checkpoint
- Elastic 做下游分析

---

## 9. `elastic` 里还有哪些分析输出

`elastic.cpp` 里除了 `strain_*.txt` 外，还定义了晶体序参数分析：

- `cryst_order(...)`
- 输出 `q4q6_*.txt`

但当前主流程里实际调用的是：

- `elastic3D(...)`

没有看到主流程真正调用 `cryst_order(...)`。

因此当前更确定的主输出是：

- `Phimax_*.txt`
- `strain_*.txt`

GUI 可视化侧则已经支持：

- `Phimax_*`
- `strain_*`
- `q4q6_*`

---

## 10. 后续继续分析时建议的切入顺序

如果下一步要继续看代码，建议按这个顺序：

1. `misfit` 主循环  
   看它每一步的物理项与输出
2. `cvd` 相比 `misfit` 的差异  
   重点看初始化、主方程、参数默认值
3. `elastic` 的后处理流程  
   看 VTI -> 峰值提取 -> 应变计算
4. GUI 调度链  
   看每个 tab 如何传参、输出目录如何组织

---

## 11. 一句话版本

一句话总结：

> `misfit` 和 `cvd` 是并列的主仿真模型；`elastic` 代码血缘上属于同一套体系，但功能上更像读取它们输出做后处理分析的下游程序。
