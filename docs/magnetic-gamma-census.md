# R28 受限 Γ 磁分导覆盖诊断（census）

> **历史统计边界（R35）**：以下是 R28 运行记录，不是当前拒绝计数。生产输入契约已加入
> [精确原点对齐](magnetic-gamma-input.md#verified-origin-only-source-alignment-r35)，
> 980 条合格 Γ 源行全部通过源帧校验（152 条使用非零原点）。这不等于全分导成功率；
> 新版完整 census 的结果单独记入[当前账本](../AGENTS.md)，不能沿用下文 13,264 的旧拒绝数。
> **R36** 又增加了[经幺正子群的精确桥接](magnetic-unitary-bridge.md)：保留原先的固定见证，
> 只对缺口使用普通 H 嵌入与 M→H frame 构造额外候选并完整回验。因此以下 7,970/8,751
> 也是旧版记录数，不能作为新版覆盖率；R24 自身的候选族及其认证口径没有改动。
> **R37** 在保留这些选择之后，再尝试完整源格的 primitive 归一化与同格的其它普通 H
> 记录。它们只提供经验证的几何见证，不声称证明 OP stabilizer；新统计仍见当前账本。
> **R38** 增加独立的字面记录格障碍证书：56 条不可能在该字面格上 exact 嵌入，
> 当轮另 1,078 条缺口仍未解决。它不改变正见证或 TSV，也不否证其它 cell/axis convention。
> **R39** 再以成对辅助群与反幺正平移提议 map，最后用原始 M 精确验收；旧路线优先级不变。
> 候选识别未决独立保留，不作为不存在性证明。
> **R40** 在此后枚举 Type III 的有限分级点扩张，同样用原始平移完整验收；
> 不改变旧见证或表示域。最新覆盖、识别器限制与冻结数据见当前账本。
> **R41/R42** 修复共享定心/选轴，不扩表示域。R42 fixed 为 16,591，150,100 次探针成功，
> 仍缺 130 条几何记录。旧 149,956 个成功探针均保留，但 SG20 T1 /UNI23 的 8 行选中
> 另一合法嵌入（4 行目标标签改变），不是整个旧 TSV 不变。路线优先级不变，而同一路线
> 候选增加可改变全序最小者；[新旧明确嵌入均有逐操作回归](<../tests/magnetic_antitranslation_bridge.rs>)。
> 自动几何选择仍不证明 OP stabilizer，也不应将不同嵌入的目标标签直接当成同一结果比较。
> **R43** 修复显式操作查询的约化/坐标还原，候选未决 6→0；fixed 16,591、成功探针
> 150,100、剩余 130 均不变，151,068 行 TSV 与 R42 逐字节相同。首-H 单独报告新增
> 29 条见证，但原已由其它候选覆盖，不能计作普查新增。仍不承诺任意病态基的识别完备性。

本文件是**只读诊断**的记录，不含生产算法改动。程序
[examples/magnetic_gamma_subduction_census.rs](../examples/magnetic_gamma_subduction_census.rs)
把已接受的三个入口指向全语料：

* R25 输入契约 `irrep::magnetic_input::verify_scalar_gamma_input`；
* R26/R27 受限分导 `irrep::magnetic_subduction::subduce_scalar_gamma`；
* R24 嵌入搜索 `irrep::magnetic_embedding::parent_embedding_verdict`（只用来**固定**每个
  record 的一个 embedding：**复用 R24 现有候选族，本轮不新增、不扩大 setting 搜索**）。

基线：`HEAD = ec1a530`（未提交的 R27）＋本轮的**两个新文件**。R25/R26/R27 生产源码、
现有测试、现有文档与旧 WIP 均未改动。R28 未 commit/push。

## 1. 分母：谁是 population，谁被排除

| 层 | 定义 | 计数 |
|---|---|---|
| record | 230 个母群磁 isotropy 表的行（`irreps_of(sg)` × `magnetic_subgroups()`） | **16,721** |
| record（独立口径） | `query::magnetic_isotropy_subgroups_of(sg)` 的行数之和 | **16,721**（一致，gate 逐 SG 断言） |
| Γ source | k 矢量 numerators 全 0 的 irrep | **2,059** |
| ├ 非 spinor、dim = 1、字符实 ±1 | 进入 probe 域 | **980** |
| ├ 排除：spinor | 具名 `spinor_source` | 746 |
| ├ 排除：dim > 1 | 具名 `dimension_gt_one` | 333 |
| ├ 排除：字符行不可读 | 具名 `character_unavailable_*` | 0 |
| ├ 排除：复字符 / 网格外 | 具名 `non_real_character` / `character_off_grid` | 0 / 0 |
| 非 Γ 的 irrep（**含 spinor**） | 按定义出域（仅记账） | 6,329 |

`dim` 直方图（非 spinor Γ 源）：`{1: 980, 2: 246, 3: 87}`。注意 980 与 R25 契约文档提到的
“980 个 scalar Γ 行”数量吻合（本轮测到字符行不可读、复字符、网格外各为 0，所以按
“非 spinor + dim 1 + Γ + 字符行可读且实 ±1”筛出的数量就是全部 dim 1 Γ 源的数量）；两种
读法是否逐行同一集合，本轮没有独立逐行核对。`k_label == "GM"` 的标签读法也恰好得到
2,059，与数据读法（numerators 全 0）一致，但筛选本身只用数据读法（契约自己的判据）。

**pair** 是**每个母群内部**的乘积和 `Σ_sg records(sg) × sources(sg)`（不是两个总量相乘）：

* 全部 pair（record × 任何 Γ 源）：**187,905**；
* 其中源合格：**75,534**；
* 源域外 pair：**112,371**。

## 2. 固定 embedding：一个 record 一个，且与源/结果无关

对每个 record 调用 `parent_embedding_verdict`，只保留
`realises_record_lattice == true` 的条目，按**全序键**取第一个，策略固定为
`RequireExact`：

1. 首选 `map = I`（表格帧，与已接受的 R25/R26 见证同帧）；
2. 其余按 `(map, shift)` 的精确 `(分子, 分母)` 项排序。

键只读 record 本身，**不看 source、不看分导结果**，因此不可能“按成功挑候选”。本轮不枚举
多个 embedding，故这是**条件分母**：成功率是“在 R24 族的第一个 realising embedding 且
`RequireExact` 下”的，不是“任何 embedding 下”。

| record 有无固定 embedding | 计数 | 说明 |
|---|---|---|
| 有 | **7,970** | 其中 2,051 用恒等 map，1,647 恒等 map 且 shift = 0 |
| 无 | **8,751** | 不 probe；TSV 中 `map`/`shift` 为 `-`，stage = `no_fixed_embedding` |
| ├ `no_realising_candidate` | 635 | R24 返回了包含性条目，但没有一个精确实现记录格 |
| ├ `certified_absence` | 8,116 | 候选族 `{S·B}` 内每个 map 都被精确条件否证 |
| ├ `undecided_residual` | **0** | 算术上限/分量上限造成的**未决**，绝不并入“无解” |
| └ record error | 0 | geometry / operations / verdict 错误 |

`certified_absence` 的作用域只是 R24 的候选族（48 个有符号置换 × 5 个基），**不是**
“任何 setting 下都不存在嵌入”，也不是证明记录不属于母群。

## 3. probe：真实双 η 调用与唯一 stage

每个合格 pair 真实调用两次（`η = +1` 与 `η = -1`）：
`verify_scalar_gamma_input` → （成功时）`subduce_scalar_gamma`。每个 probe 只归**一个首次
失败 stage**，并保留完整具名 variant 与单行 `detail`：

* `verify_input`：R25 契约拒绝；
* `subduce`：契约接受但 R26 分导拒绝；
* `success_check`：两者成功，但本 census 的独立复核不成立（**必须为 0**）；
* `success`：`mult`、维数恒等式（`dim(source) = multiplicity·dim(target)`）、有限非零
  单位 `U`，以及逐类重建 `A(g)·(U z) = U·(D(g) z)`（anti 时左因子共轭）都通过。

复核的两个判据都是**正谓词**，但形状不同：`is_finite_unit` 拒绝 NaN、无穷**与零**（交错因子
`U` 必须是非零单位），`within_tolerance` 要求残差有限且 `≤ 1e-9`——**合法的零残差通过**，
NaN/无穷才失败。两者都不会反过来从否定式比较里漏过 NaN。

## 4. 全语料结果（230 SG，debug，8 线程）

两次全量：gate 重做前 **354.5 s**，重做后 **344.9 s**。两次的 TSV 第 1–16 列逐字节相同
（`md5 = 79d6ee48…`），唯一差异是 `certified_absence` 的 `detail` 措辞按导师要求补上“只对
R24 候选族成立”的作用域；统计路径未变。

| 量 | 计数 |
|---|---|
| probes（= 2 × pairs_probed） | **70,374** = 2 × 35,187 |
| `success` | **57,110**（即 28,555 个 pair 两 η 都成功） |
| `verify_input` 拒绝 | **13,264**（6,632 个 pair 两 η 都拒绝） |
| ├ `CharacterRowSettingMismatch` | 13,264（**唯一**出现的输入拒绝变体） |
| `subduce` 拒绝 | **0** |
| `success_check` 失败 | **0** |
| 无固定 embedding 的 pair（TSV 单列） | 40,347 |
| TSV 数据行 | **151,068** = 2 × (35,187 + 40,347) |
| η 状态/target 标签/operator 表不一致 | 0 |
| η 之间 intertwiner 不同 | 25,853 pair（反幺正号的自然结果，非不一致） |

**没有**出现 `LatticePolicyUnsatisfied`（因为候选已按 `realises_record_lattice` 过滤）、
`EmbeddingRejected`、`Arithmetic`、`UnresolvedTargetData`、`Target*`、`NoMatchingTarget`、
`AmbiguousTarget` 或 `DimensionNotConserved`。算术未决为 0，未被当作无解。

### 4.1 具体阻塞：`CharacterRowSettingMismatch`

失败集合非常集中（可逐行在 TSV 里追溯）：

* 失败出现在 **24 个空间群**：48, 50, 59, 68, 70, 85, 86, 88, 125, 126, 129, 130, 133, 134,
  137, 138, 141, 142, 201, 203, 222, 224, 227, 228；
* 至少失败一次的 **152** 个 `(sg, source)` 源行；其中 **24** 个源行在任何被 probe 的
  record 上都不成功，其余 128 个源行只对**部分** record 失败（同一源行在别的 record 的
  映射下落在该行操作集内）；
* 1,159 个 record 参与过失败；SG 48 / 50 / 59 的 probe 全失败（成功 0），SG 89 的 224 个
  probe 全成功；
* 与 R25 契约文档记载的“980 个 scalar Γ 行中 152 行、24 个空间群”[设置不一致] 数字**一致**。
  口径并不相同：文档按**逐操作比对**母群操作表与源字符行，本处按
  `record × 源` 的**可映射性**。这里是**不同测量单位的数值对照**，不是独立算法验证：本 census
  的每个 probe 都直接调用同一个 R25 输入契约（`verify_scalar_gamma_input`），与该历史测量
  不构成独立复核；152/24 只是两把不同的尺子量到同一个规模。

例（TSV 原样字段）：`48  GM1+  direction_index 0  UNI 355  GM1+  +1  verify_input
CharacterRowSettingMismatch  map=1/1,…  shift=0/1,…  target=-`，
`detail = SG 48 (UNI 355): the source character row and the parent operation table do not
describe the same setting (rotation [[-1,0,0],[0,-1,0],[0,0,1]] has no row entry modulo the
parent lattice)`。

**这不等于数据缺失**：`CharacterRowSettingMismatch` 只陈述“母群操作表与源字符行不在同一
setting”，本轮没有搜索对齐变换，也没有证明对齐不存在（与 R25 的措辞边界一致）。

## 5. 门禁（`--gate`）与它的独立性

门禁的核心不是自洽计数，而是**读实际发出的 TSV 行**：

1. **行身份 validator**（`validate_tsv` / `validate_global_tsv`）逐 SG 把 TSV 文本解析回
   `(sg, condensing irrep, direction_index, UNI, source, η)` 键，与 **inventory 推出的期望
   集合**（records × 合格源 × {+1,−1}）对照：重复键、未知键、缺失键、字段数、stage 名、
   `no_fixed_embedding` 的 reason 名、以及**每个 record 的 `(map, shift)` 跨 source 与两个 η
   必须是同一个**，任何一条不成立都让 gate 非零。文本与 counters 对账（成功/各阶段/行数），
   删行或同长度重复 η 会同时触发 missing 和 duplicated。
2. **遍历完整**：逐 SG 的 record 数等于独立 `magnetic_isotropy_subgroups_of` 计数
   （16,721 == 16,721）。
3. **partition**：records = 有固定 + 无固定；无固定再按四种 reason 分区；源排除之和 =
   Γ 源 − 合格源；pair 恒等式在**每个母群内部**成立（乘积和，不用总量相乘）。
4. **锚点**：R25/R26/R27 全部见证按其**自己的** pinned embedding 重跑（SG 2 UNI 4 的 Type I、
   SG 2 UNI 6 的 Type III、SG 1 X1/UNI 3 的 `map=2I` Type IV、SG 3 UNI 10 的独立目标身份、
   SG 10 UNI 53、SG 12 UNI 56、SG 5 UNI 22、SG 11 UNI 59 的 `GM1+`/`GM1-`），并断言
   SG 1 X1/UNI 3 在 `RequireExact` 下以 `LatticePolicyUnsatisfied` 具名拒绝；violations = 0。
5. **η 一致**：同一 pair 的两个 η 必须在**状态与具名 variant** 上一致（0 不符）。
6. **6 + 3 = 9 个 example 单测**：库存/独立查询一致（含 16,721）、源筛选网格、候选键全序与
   顺序无关、锚点、2D 边界具名，加上 3 个新测试 —— TSV validator 的 6 种变异（删行 /
   同长度重复 η / 未知 record / 未知 source / 单 η map 不一致 / counter 与文本不符）、
   validator 对无固定 embedding 行的覆盖（含未知 reason）、两个正谓词各自的边界
   （`is_finite_unit` 拒绝 NaN/无穷/零；`within_tolerance` 接受合法零残差、拒绝 NaN/无穷）。

复现（本机，显式 target，离线）：

```text
CARGO_TARGET_DIR=/home/liuyichen/TB_rs/cryspglib/target CARGO_NET_OFFLINE=true \
  cargo build --example magnetic_gamma_subduction_census
RAYON_NUM_THREADS=8 ./target/debug/examples/magnetic_gamma_subduction_census --gate \
  --tsv target/logs/r28-census-worker/census.tsv > target/logs/r28-census-worker/census.txt 2>&1
cargo test --example magnetic_gamma_subduction_census
```

线程无关性：同一份源码下 `--smoke` 的 8 线程与 4 线程 TSV **逐字节相同**
（`smoke2.tsv` == `smoke3-4threads.tsv`，SHA `945d934c…`）。

## 6. 二维 Γ 下卡的来源可用性证据（不实现）

三条**真实调用**的边界证据（不是从维度推断）：

| 源 | SG | 标签 | dim | compound | matrices_len | characters_len | 契约回答 |
|---|---|---|---|---|---|---|---|
| 扫描中最早的 dim ≥ 2 Γ 源 | 75 | `GM3GM4` | 2 | **true**（复合/实化行） | 16 | 4 | `NotScalarSource` |
| 最早的**普通**（非复合）dim ≥ 2 Γ 源 | 89 | `GM5` | 2 | false | 32 = 8 ops × 2² | 8 | `NotScalarSource` |
| 显式点名的第二个普通父群交叉核对 | 156 | `GM3` | 2 | false | 24 = 6 ops × 2² | 6 | `NotScalarSource` |

**不能**把第一条读成“普通复不可约 dim 2 已具备”：`GM3GM4` 是复合标签，其 2 维来自实化/
复合行；第二、三条才是普通单值维度 2 源，且它们确实带完整矩阵存储。三条都按
“(sg, dim) 全序”之外的独立路径取到：第 2 条由扫描得到（SG 89 早于 SG 156），第 3 条是
显式命名的父群（同一个表，独立读出），所以它们同时覆盖了“扫描结果”和“点名核对”。

两侧存量都真实存在：

* 源侧：非 spinor Γ 源 dim 直方图 `{1: 980, 2: 246, 3: 87}`；
* 目标侧（在 Γ 读磁 summary）：SG 75 的 24 行/10 个 UNI 中 9 个 UNI 的 Γ corep 含 dim ≥ 2，
  直方图 `{1: 23, 2: 18}`，其中 1 个是“单一非 spinor 源”的 corep；SG 89 的 52 行/25 个 UNI
  中 21 个 UNI 含 dim ≥ 2，直方图 `{1: 92, 2: 35}`，其中 6 个是单一非 spinor 源；SG 156 的
  52 行/19 个 UNI 中 14 个 UNI 含 dim ≥ 2，直方图 `{1: 52, 2: 26}`，其中 6 个是单一非
  spinor 源。

代码层面的卡点（可追溯）：

* [magnetic_input.rs:894–900](../src/irrep/magnetic_input.rs#L894-L900)：`dim != 1` ⇒
  `NotScalarSource`（spinor 在 887–893、非 Γ 在 901–909 各自具名）；
* [magnetic_subduction.rs:745–762](../src/irrep/magnetic_subduction.rs#L745-L762)：
  `resolve_source` 要求 `dim == 1` **且** `k_vector().denominator == 1`；后者是比输入契约
  更严的一道（本轮 980 个合格源上没有触发，因为它们的 kd 均为 1）；
* [magnetic_subduction.rs:765–771](../src/irrep/magnetic_subduction.rs#L765-L771)：
  `is_one_dimensional_single_source` 只把 `dim == 1`、`CorepType::A`、单一 1 维非 spinor
  源的 corep 当候选；
* [magnetic_subduction.rs:1258](../src/irrep/magnetic_subduction.rs#L1258)：
  `solve_intertwiner` 解的是**标量** `U`。

结论的边界：`NotScalarSource` 只证明当前契约在输入侧拒绝 dim ≥ 2 源；矩阵存储非空只证明
数据存在。本轮**没有**证明二维矩阵法（相位、重数、逐类重建、setting）已经具备，也没有
实现它。

## 7. 覆盖限制（不作夸大）

* **条件分母**：只在“R24 realising 候选族里排序第一的那个固定 embedding + `RequireExact`”
  下统计；不是任何 embedding / 任何 setting 的成功率。
* **8,751 / 16,721（52.3%）的 record 根本没有被 probe**（无 realising 候选），其中
  `certified_absence` 只是候选族内的结论。
* 域外**未覆盖**：746 个 spinor Γ 源、333 个 dim ≥ 2 Γ 源、6,329 个非 Γ irrep；以及 grey
  记录不被本输入来源接受这一既有边界。
* 失败只保留具名 variant：`CharacterRowSettingMismatch` 是拒绝名，不是根因证明；本轮无
  setting 搜索/逆变换，也不声称“不存在对齐”。
* 无 full-star、无 B/C 型、无高维分导、无新生产算法；旧极端 i128 格原语边界不受影响。

## 8. 交付物与版本

* 新程序：[examples/magnetic_gamma_subduction_census.rs](../examples/magnetic_gamma_subduction_census.rs)
  SHA-256 `9a9840fdf0a54506ea3b943106af7507e2b0417685dfbb3095033f788b6b6090`；
* 新报告：本文件；
* worker 日志：`target/logs/r28-census-worker/`
  * **最终全量**：`census.txt` SHA-256
    `ea63880193e826881382bd78701eac6025573a23eb87103dae3be6f9e074340c`、
    `census.tsv` SHA-256 `1e9e5194628ea93a06f4ca75a2f9e95bfa641b486b5af84477e6116897cb3f6c`
    （30,540,820 B，151,069 行 = 1 表头 + 151,068 数据行，gate exit 0，332.3 s）；
  * **执行版本边界**：上面这次全量是用“纯注释/措辞修正之前”的那一版源码构建的二进制跑的；
    之后只改了 example 的注释与本文措辞（没有统计、选择或算法路径改动），example 的 SHA
    由 `7b0866c1…` 变为 `9a9840fd…`。当前源码已重新 `cargo clippy -- -D warnings` 干净、
    `cargo test --example` 9/0，并用 `--smoke --gate` 复核（`smoke4-wording.tsv`）：TSV 与旧
    smoke 逐字节相同（SHA `945d934c…`）。**完整 230 全量的独立复跑由 Leader 用当前源码执行。**
  * 第二次全量（加“点名交叉核对父群”的边界行之前）：`census-v2.*`，TSV 与**最终全量逐字节
    相同**（同一 SHA）；
  * 第一次全量（gate 重做前，保留为历史）：`census-v1.txt` / `census-v1.tsv` SHA-256
    `cc81ad1c40d5edb69e7c88e0c9cc10fbfbe9d8c75873577c5fda3bd0bc49af7f`；与最终全量的
    **第 1–16 列逐字节相同**（`md5 = 79d6ee48…`），唯一差异是 `certified_absence` 的 detail
    措辞补上了候选族 scope；
  * `smoke2.tsv`（8 线程）与 `smoke3-4threads.tsv`（4 线程）SHA-256 都是
    `945d934ce32925ca2c6549b9ab0f000665b4183aaee97a973a8da9588512c483` ⇒ 逐字节相同；
  * `hashes-v1.txt` 记录历史全量的哈希。
* 版本：`git rev-parse HEAD` = `ec1a5306d5e314075482ea00fa0c47fefac868b8`；`git status --porcelain`
  中除本文件与 example 两个新文件外，其余（`AGENTS.md`、`docs/magnetic-gamma-subduction.md`、
  `src/irrep/magnetic_subduction.rs`、`tests/magnetic_subduction.rs` 的 M，以及
  `src/irrep/magnetic_embedding/joint_congruence.rs`、
  `tests/magnetic_subduction_provenance_contract.rs` 的 `??`）都是本轮开始前就存在的
  R25/R26/R27 状态，本轮未触碰。R28 未 commit/push。
