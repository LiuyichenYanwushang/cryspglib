# R6：参数化 k 的分导（能力契约与推进计划）

R5 已把**普通离散标量**分导在固定语料上闭合（`366,260/366,260` 完整分解，见
[subduction-audit.md](subduction-audit.md)）。R6 处理的是 R5 分母之外的那一族：
pinned `isotropy_w_subduce_*` 的 **5,756 行**，它们的母群 irrep 活在**参数化 k 域**
（73 个源，全部是过 Γ 的直线 `k = Γ + t·v`）。

## 0. R5 基线冻结（R6 的起点）

**代码基线**（本文件与 R6.0 首个提交的父提交）：

```
c501ac0  Follow-up review: close the two remaining gate holes in the generation pipeline
```

R5 收口链：`dbf6d5c`（R5 报告）→ `ff2e4a4`（三门口禁判词）→ `1d08b3c`、`e9e8e8e`
（复核 A 修正）→ `6c3f335`（复核 A/B 处置）→ `c37d3b9`（第三方复核处置）→
`c501ac0`（跟进门禁）。

**R5 验收证据的适用范围（诚实口径，勿写成"c501ac0 通过全表审计"）**：

| 提交 | 已实跑 | 说明 |
|---|---|---|
| `c37d3b9` | **全表三门口禁**：517.1 s、exit 0、判词 `VERDICT complete scope=global full_decomposition=complete`；`embedded_records=15239/15239`、`positive_stored_compared=94271/94271`、`probe_full_success=366260/366260`、`identity_only=0`、`absent_zero=271989`、`frobenius_evaluated=1895/1895`、`w_computed=5756/5756`、`hard_failures=0`；同批 Rust/Python/example 测试与 clippy | 运行时 Rust 源码在这一版 |
| `c501ac0` | Python 门禁回归（`test_build_table` 7、`test_subduction_settings` 10、reviewer 的 `gate_regression_tests.py` 4/4） | **只改脚本与文档**；运行时部分沿用 `c37d3b9` 未变更基线的审计证据，未重跑全表审计 |

**R6 起点必须复现的运行时数据校验值**（SHA-256，`git show <commit>:<path>` 可核）：

```
30dc93c63988b8cf44521355cc60435ffc081f04d90c24acd66e31c26331a165  src/irrep/w_little_characters_data.rs
9f8635c4cf6dd769446f51d3bb7436b18d7e6064062da878568553983da011e5  src/irrep/subduction_settings_data.rs
1b7285390c986e8bb1c1ba4b660d5d479d0df015b5f156328faebcc57edbce63  src/irrep/generated_data.rs
568667bfc8027095537d642297b319c872d00016b868143c666f90d5931d9f7b  isotropy_subgroup/iso.zip
65512dabbbf0934690aaaa25448d59b4d8b4f5bdda5c5c5ec1834f50f6ed016a  isotropy_subgroup/data_isotropy.txt
```

**历史重生成（开放事项）**：`src/irrep/subduction_settings_data.rs` 的**逐字节**重生成
需要未跟踪的 `target/task9/` 证据链（最终 child shift 是多次增量修补的合并），本仓库
不声称它可从 tracked 输入重现；已核的是内容、覆盖、身份、幺模性与引擎嵌入
（见 [subduction-audit.md](subduction-audit.md) 的 Provenance status）。
**边界**：R6 **不修改**这张表，因此该事项不阻断 R6；一旦 R6 需要修改或重新生成正式
冻结表，数据来源与确定性重建就是那次数据变更的**前置条件**，不能一边改旧表一边继续
依赖未跟踪的增量文件。

## 1. 两个能力承诺：A 与 B（R6 最重要的验收边界）

**能力 A：给定明确参数值，完成该点的完整分导分解。**
输入 = `(isotropy 记录, 母群参数化源, 有理参数 t)`；输出 = 该 `t` 处的**完整分解**
（每个子群 star 的 q、源行 k、目标 irrep 身份、维数、重数），或一个**明确错误**。
A 只承诺"这一个参数点算对"，不承诺任何参数区间上的结论。

**能力 B：对整个参数族给出分区适用的结论。**
要求说清哪些参数区间共用同一套处理、哪些特殊取值必须转入别的分支、每条结论的适用
条件（例如"对 `t ∉ S`（有限例外集）分解形式为 …"）。

**推进顺序：先 A 后 B。** 尤其禁止把"若干参数采样点通过"写成"整个参数族覆盖完整"；
R6.2 的覆盖说明必须逐条区分**已计算 / 有独立对照 / 仅内部一致性 / 不支持**。

## 2. 范围（本轮不做的事）

* 只做 **ordinary scalar**：不加入双群/spinor、不加入磁共表示（那是 R9–R11）；
* 不做参数族的符号化处理（B 阶段只给分区与例外集，符号化推导不在本轮）；
* 不改 `subduction_settings_data.rs`、不改 pinned w 特征标表（`w_little_characters_data.rs`）。

## 3. 阶段与验收

### R6.0 能力契约与测试样例（**已完成**，2026-09-23）

交付：`docs/subduction-conventions.md` §16（契约）、`OFFICIAL_LINE_PARAMETER` /
`official_line_parameter()` / `subduce_line_at_parameter` / `LineSubduction`。
每个"支持"声明都有测试，不支持一律显式报错（见 §16 的失败语义两条负例）。

交付：`docs/subduction-conventions.md` 新增一节（输入坐标/参数表示/支持的表示类型/
输出含义/错误状态），以及代码侧的显式入口与错误变体。验收：

1. 契约里每一个"支持"声明都有对应测试；**不支持不能降级成零结果**（缺数据 →
   `MissingChildStarData`，源不匹配 → `LineSourceMismatch`，参数非法 → 明确错误）；
2. 明确写出参数帧：`t` 乘的是**母群 conventional 倒格基**下的方向（与 pinned
   `little_k` 的换算见 §4），波矢在引擎内部换算到母群 primitive 倒格基；
3. 每个错误变体至少一条负例，且断言错误类型（不只看 `is_err()`）。

### R6.1 具体参数点的端到端分解（**已完成**，2026-09-23）

交付：线源完整分解入口 + 轨道化折叠修复 + 5 条常驻单元测试；审计的 w 门禁由
Γ-only 路径升级为**完整分解**。实测：`t = 1/4` 下 5,756/5,756 条 pinned 行完整分解
且恒等重数等于 pinned 频率（604.8 s，三门口禁 exit 0）；一般参数、特殊值两侧、
独立于多重度求解器的臂数对照、外来源与越界点的负例全部通过（明细见 §16 的验收表）。
一般参数上**没有外部 oracle**：那里的证据级别是内部一致性加一条独立的几何计数，
唯一的外部对照仍是 `t = 1/4` 的 pinned 频率（与全表审计），R6.2 必须如实分列。

交付：`subduce_line_at_parameter(...)` 的完整分解（复用 `build_block`/`reconstruct`），
先在少量有理参数点上跑通，再扩大。验收（三组，见 §5 的测试重点）：

1. **已知点回归**：`t = 1/4`（官方参数约定，`OFFICIAL_LINE_PARAMETER`）下，完整分解
   算出的**恒等重数**必须等于 pinned 频率行（独立来源：官方程序 `SHOW FREQUENCY`），
   覆盖全部可得上下文；
2. **一般位置与特殊值两侧**：一般 `t`（折叠 q 落一般位置）、`t` 使折叠 q 命中
   stored 子群 k、`t` 使折叠 q 落特殊点（非平凡小余群），以及这些特殊值**两侧**的
   精确有理点，都要有端到端样例；
3. **代数一致性**：维数守恒（`Σ mult × dim × star = 母群 full-star 维数`）、重数为
   非负整数、逐子群代表元重构逐项通过（引擎门禁，出错即 Err）。

### R6.2 参数族覆盖与审计（**下一步**）

交付：可审计的覆盖说明（区间 + 例外集 + 每条结论的证据级别），以及对 5,756 行的
重述口径：**哪些行的恒等频率在任意 `t` 下成立、哪些只在特定 `t`（如 1/4）成立**。
验收：覆盖说明里的每个数字都能由仓库内脚本/测试复算；四种证据级别分开列。

## 4. 帧与参数约定（R6.0 必须写死的东西）

**整套量都在母群 conventional 倒格坐标下**（R6.1 实测更正：本节早先写的"pinned
`little_k` 在 primitive 倒格基、引擎要把 `k_conv` 换算到 primitive 再折叠"是**错的**，
代码从来没有做过那次换算）：

* 冻结 little 表的 `direction` 是**官方打印的 conventional 方向**（SG 196 `DT` 的
  `("0","2","0")`），冻结的 little 群操作/平移也自述为 conventional 帧；
* pinned `little_k` 同样是 conventional：SG 196 的 `X1 = (0,1,0)`、`L1 = (1,1,1)/2`
  （primitive 读法应为 `X = (1/2,0,1/2)`）。`exact_primitive_basis` 只提供**格子**，
  不改变坐标系；
* 引擎入口接受 `(表, 有理 t)`：`k = t · direction`，**同一向量**直接进
  `fold_wave_vector(T, k)`（`T` = 嵌入的精确仿射变换），臂集合、字符求值与折叠共用
  这一个向量，没有任何二次换算；
* **没有任何约化**（R6.2 最终读法；R6.1 曾在 `canonical_wave_vector` 里把 `k` 约化进
  母群倒格基本胞，**该修复已被撤销**，见 §16）：冻结的 `D` 是在 `k = Γ` 处解出的纯 Γ
  点字符，参数通过 Bloch 因子 `exp(2πi t (v·T))` 进入，所以 `k(t) = t·direction`
  **原样**参与字符求值与折叠。把 `k` 约化而保留 `D` 等于换一条带来算；一个参数步是
  标签的 monodromy 位移（`(k + K, M_K(α)) ~ (k, α)`，`cryspglib::irrep::line_monodromy`），
  不是 gauge。
* 交叉证据：`10038 DT1` 在 `t = 1/4` 折出 6 个臂（2 + 4），little 群冻结为
  `{E, C2y}`（正是 conventional 帧下 `(0,2,0)` 的稳定子，12/2 = 6），且
  `parent_dimension = little_dim × arms = 6 = pinned DT1 dim`，三者同时吻合。
* `t` 与 `t + Δ` 何时是同一个表示：当 `Δ · direction` 是母群倒格矢量（含 centering
  消光）时相同。R6.1 已用**可执行的等价性测试**钉住这一条
  （`a_reciprocal_vector_shift_of_the_parameter_changes_nothing`，覆盖 40 条曾受
  gauge 影响的 pinned 行在 `t = 5/4, 9/4, 13/4` 上与 `t = 1/4` 逐项相同），而不是只
  写口号；一般位移规则的完整刻画仍属 R6.2。

## 5. 测试重点（reviewer 建议，写入验收）

1. **参数与分支**：一般参数点、特殊取值、特殊值两侧的精确有理点；不能只测容易回落到
   旧固定点数据的 `t = 1/4`。
2. **数学一致性与独立对照**：维数守恒、重数非负整数、限制后特征标重构**之外**，
   保留少量可手算样例与独立 oracle 样例（pinned 频率行 = 官方程序输出）；
   不以内部自洽代替外部正确性证据。
3. **约定与失败语义**：等价坐标/原点变换、分数 setting（ordinal 26 一类）、缺数据或
   不支持的输入；守住"计算失败不是重数为零"这条已经在 R5 修过的边界。

## 6. Review 安排

按 reviewer 建议：不再为"证明没有 bug"无限期冻结开发；下一轮 adversarial review
放在 **R6.0 契约确定 + R6.1 第一条端到端路径完成之后**，围绕新增风险（参数帧、
分支、失败语义、独立对照）检查。
