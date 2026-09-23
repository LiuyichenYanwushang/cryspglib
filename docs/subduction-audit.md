# 完整分导的全表审计（任务 9）

任务 9 的工具分别检查几何来源和生产 API 的计算覆盖。**收集到官方基矢、恒等项
吻合、完整分解成功是三个不同的结论**；任何一个都不能替代另外两个。

（R2 之前）14,713 个恒等-only probe 的逐星缺口清点见
[subduction-gap-census.md](subduction-gap-census.md)：21,136 个缺失星，
按子群/k-star/setting 归并为 989 组；R2 补齐子群 #1 后余 12,878 条，R4 批次 1
补齐平凡小余群的构造目标（批次 1）、一维投影 catalogue（批次 2a）与二维投影表
（批次 2b）后**余 0 条**：全表 366,260 个 probe 全部完整分解（详见下文两个门禁的
当前结果）。

## 生产 API 审计

在工作区根目录运行：

```bash
CARGO_TARGET_DIR=$PWD/cryspglib/target cargo run --release -p cryspglib \
  --example audit_irrep_subduction -- --output /tmp/subduction.tsv
```

`--parent 221` 或 `--ordinal 12400` 限定诊断范围。退出码 2 表示所选门禁覆盖不完整，
退出码 1 表示发现不一致（计算错误、频率冲突、计数守恒破坏优先于一切门禁），
退出码 0 只表示在**本次报告的范围内**未发现不一致。三个门禁互相独立：

- `--require-complete` 要求所选范围内**普通恒等分导表**闭合：15,239 条记录全部冻结
  embedding、94,271 条存储正项全部复现、366,260 个标量 probe 全部有精确结果
  （完整分解**或**恒等内容）、几何与 Frobenius 检查无未计算项。
- `--require-full-decomposition` 要求所选范围内每个普通 probe 都是**完整分解**；
  恒等-only 结果在这个门禁下**计为不完整**（它在 `--require-complete` 下算通过，
  因为回答恒等重数是另一个问题）。**当前状态**：R4 三个批次闭合后全局为
  `identity_only=0`、`incomplete=0`，该门禁 **exit 0**（判词见本文 R5 报告）。
  **历史**：R4 之前它曾 exit 2 并报告 14,713 个缺口，摘要行形如
  `full_decomposition: ... identity_only=14713 incomplete=14713 global=covered`；
  这段历史口径保留在 R4 各批次小节与第六轮表里，不再代表当前全局状态。
  这条门禁的存在就是为了防止再次把"恒等项通过"当成"完整分解完成"。
  范围语义：`--parent`/`--ordinal` 限定运行时只报告该范围的完整分解状态，
  判词带 `global_coverage=not_established`；局部全过不会被表述成全局覆盖完成。
  一个**选不中任何记录**的范围（例如 `--parent 2 --ordinal 0`，ordinal 0 属于
  SG 1）不是"空范围内通过"，而是 `empty scope` 错误、非零退出：零 probe 的报告
  证明不了任何事，CI 里按 (parent, ordinal) 循环时一个笔误不会冒充成功。
- `--require-w-complete` 另外要求 5,756 条 `other_wave_vector_subduction` 也被
  **引擎**计算。**现已闭合**：审计只走块路线
  （`line_trivial_content_via_blocks`：把直线源的臂按 `LINE_PARAMETER = 1/4` 折叠、
  归组后交给既有的 `build_block` 折叠/解块流程），逐行与 pinned 频率比较，
  实测 `w_scope: rows=5756 computed=5756 uncomputed=0 mismatched=0 engine_errors=0`，
  两个旧门禁同时开启仍退出 0、判词 `VERDICT complete scope=global`。
  w 门禁只覆盖**恒等重数**；这些源的完整分解是独立范围（R6）。
  注意审计**不按期望答案选择算法**：不一致一律计入 `w_frequency_mismatch`，
  块路线返回 `Err` 一律计入 `w_engine_error`；两个计数都已并入
  `hard_failures()`，因此在任何开关组合下都退出 1（永久负例见 example 的
  `a_frequency_mismatch_fails_under_every_flag_combination` 与
  `a_w_engine_error_fails_under_every_flag_combination`）。
  （历史：第十七轮已把 73 个源的 k 域解出，缺口曾是 little 群特征标；
  后续由 `w_little_characters_data` 的 73/73 冻结表补齐。）

### 其它波矢行的独立官方对照

这 5,756 行现在由引擎计算，且它们的**数值**另有可复算的官方对照：
`scripts/verify_w_subduction_oracle.py` 对每个 `(母群 SG, irrep)` 运行随包 `iso` 的
`DISPLAY ISOTROPY` + `SHOW FREQUENCY`，把频率表按 `(子群号, Dir 标签)` 与 pinned 的
`isotropy_subduce_*` ∪ `isotropy_w_subduce_*` 双向比较（语义与命令见
`docs/isotropy-data-semantics.md` §4）。实测：**28 个 (SG, irrep) 组、1,150 条记录、
oracle 5,756 条 w 行对 pinned 5,756 条 w 行、0 不匹配**，其中 144 条记录是
“oracle 侧没有 w 条目”的负向核对。离线回归
`scripts/test_verify_w_subduction_oracle.py` 14 项（含折行续行解析）。

因此仓库对 pinned 表的覆盖分两条轨道，发布时必须按轨道声明：

| 范围 | 证据 | 状态 |
|---|---|---|
| 普通离散标量分导（15,239 记录 / 94,271 正项 / 366,260 probe / Γ Frobenius 1,895） | **cryspglib 引擎计算** + 几何与零项检查 | **范围内完全闭合**：366,260/366,260 个 probe 都是**完整分解**（`identity_only = missing = error = uncomputed = 0`），恒等正项 94,271/0 不匹配；正式报告见本文「R5：普通离散标量覆盖闭合」一节 |
| 其它波矢 w 行（1,006 记录 / 5,756 行） | **cryspglib 引擎计算**（`line_trivial_content_via_blocks` + 冻结 little 特征标 73/73，单一算法、不按答案选择）+ 官方 `iso` live oracle 逐行复核 | **5,756/5,756 计算且与 pinned 相同**，`mismatched=0`、`engine_errors=0`；`--require-w-complete` 退出 0，判词 `VERDICT complete scope=global`。w API 目前只返回**恒等重数**（不含完整分解） |

这 5,756 行的源现在**全部有冻结特征标**：
`src/irrep/w_little_characters_data.rs` 覆盖 73/73（`W_LITTLE_CHARACTERS_UNRESOLVED`
为空），`scripts/check_other_wave_vector_rows.py` 打印
`rows_with_frozen_table=5756 rows_blocked=0`。其中 65 个源由 Γ 兼容行唯一确定，
SG 202/203/209/210 的 `DT3`/`DT4`（348 行）由小群配对路线闭合（Γ 方程定出配对和、
pinned 源次序给出拆分，并用归档 CIR 在 X 点的字符逐操作交叉核对；详见
`docs/isotropy-data-semantics.md` §4 与 `docs/task9-remaining-work.md`）。
锚点回归：`tests/w_little_characters.rs`（6 项）。此后的轮次已在引擎侧完成
Mackey/特征标求和（`line_trivial_content_via_blocks` 走既有 `build_block` 折叠/解块
流程），把 5,756 行的频率全部算出并与 pinned 表逐行比较；下文的 Python 原型记录
保留为历史推导。

Mackey/特征标求和的 Python 原型（`target/task9/explore/proto_freq.py`，未入库）已把公式
跑通一半：对 `(母群 SG, irrep, Dir)` 用程序打印的**子群操作**（`VALUE IRREP` +
`VALUE DIRECTION <lab>` + `SHOW EL`，在母群帧里）与 little 群 `{R : Rv = v}` 的陪集，
按 `(1/|H|) Σ_{h∈H} Σ_{s∈G_L\G, s⁻¹hs∈G_L} χ_{W'}(s⁻¹hs)` 求值。SG 196 的
`6D1`（P1 子群）三源全部命中（`DT1=6, DT2=6, SM1=12`，正是 `dim V'` 锚点），
但 `4D1`/`C5`/`C11` 偏高（如 `4D1` 得 12 而非 4）。

已定位到根因方向：原型的群元规范化把平移**按 Z³ 取模**，而母群格子是**带心的**
（SG 196 的 cF 在 conventional 坐标下不是 Z³）；子群超胞平移（如 `(2,0,0)`）因此
被错误折叠成单位元，陪集/成员判定全部失真。引擎侧有正确的 `Lattice` /
`parent_primitive_basis` 层，Rust 实现必须按**母群格子**而不是 Z³ 规范化。

第二轮原型把成员判据改成**只看旋转**（`G_L = {(R,t) : Rv = v}` 本来就只约束 R ✓，
这才是对的），并把平移留给相位，于是 `6D1` 仍然全中，但 `4D1` 变成 6/6/12
（pinned 4/4/4），而且 `SHOW EL` 对该子群只打印**一个**操作（`|H| = 1`）——
原因是它给的是子群操作**模母群格子**，而 P1 超胞子群的真实平移陪集（4 个）不在其中。
于是公式要按有限群投影写完整：

```text
mult(trivial_H, V'|_H)
  = 1/(|P_H| · n) · Σ_{R ∈ P_H} Σ_{T ∈ L_H/L_parent} χ_{V'}((R,t_R)·(E,T))
  = 1/(|P_H| · n) · Σ_R Σ_T χ_{V'}(R, t_R + T)
χ_{V'}(g) = Σ_{s ∈ G_L\G, R_g 固定该臂} χ_{W'}(s⁻¹ g s)
```

其中 `n = |L_H/L_parent|` 是子群胞相对母群格子的平移陪集数（P1 情形就是记录里的
`Size` ✓），`P_H` 是子群点群（`SHOW EL` 打印的操作的旋转部分 ✓），
`t_R` 取程序打印的代表平移 ✓。对每个 `T` 求和正是「折叠到子群 Γ」的判据来源：
`χ_{V'}(E,T) = Σ_arms exp(-2πi k·s⁻¹T)`，只有 `k·T ∈ Z` 的臂不互相抵消。
`6D1` 的锚点自洽：`n = 6`、`|P_H| = 1`、`Σ_T χ = 72` ⇒ `72/6 = 12 = dim V'` ✓。
`L_H/L_parent` 的陪集代表可从官方 `SHOW BASIS` / pinned isotropy basis 取。

**与引擎现有实现的等价关系（第二十三轮）**：对 `T` 的平均并不需要显式枚举子群胞平移
——它对每个臂给出 `Σ_{T ∈ L_H/L_parent} exp(-2πi k·(s⁻¹T))`，正好是「该臂折叠到子群
Γ」的判据（只有 `k·T ∈ Z ∀T` 时不为零）。所以

```text
frequency = Σ_{折叠到子群 Γ 的母群臂的 H-轨道} m,
m = mult(trivial_{H ∩ sG_Ls⁻¹}, W'^s)
```

与 `trivial_content_with_embedding`（第六轮）已经实现的那条精确路径**同一件事**：
引擎只需要把「臂的字符来源」从 CIR 矩阵换成冻结的 little 表（外加 Bloch 相位），
折叠加 Γ 判定、子群 k 数据缺失的处理都可以原样复用。这是 Rust 侧最小改动的入口。

**（第二十四轮）不要另写折叠判定**：原型用官方 `SHOW BASIS` 打印的子群基矢
（SG 196 ord 10039 的 `4D1`：`(0,0,2),(-2,0,0),(0,-1,0)`，|det| = 4 ✓）自算
「臂 k 是否落在子群倒格」时，DT/SM 两个方向都得到 **12**（全部母群点群操作），
而 pinned 是 **4**；差别来自臂集合（应为 little 群陪集）与子群倒格的帧/中心化约定。
结论：Rust 实现**必须复用引擎里已经验证过的折叠与 Γ 判定**
（`trivial_content_with_embedding` 用的那一套，含 primitive/conventional 帧与
centering 消光），只替换臂字符来源；另行手写折叠判定会重复已经踩过的帧错误。

`--require-complete` 只表述第一条轨道（`VERDICT complete scope=global`），不会把
第二条轨道算作已完成。

每条子群记录遍历该母群的全部标量源表示，调用实际的
`subduce_full_star_with_embedding`。恒等项为零的几何证明作为交叉检查保留，但
**不能跳过完整分解调用**。历史反例（现已闭合）：ordinal 13345 的 W1–W5 当初既没有
恒等项、又缺少完整分解所需的子群 k 数据，把它们仅记成零项就会错误宣告完整覆盖；
R4 批次 1 给它们补上了构造目标（实测 31/31 完整分解，见
[subduction-r4-batches.md](subduction-r4-batches.md)）。

### 恒等内容：第二条精确路径

完整分解缺子群 k 数据时，审计改用
`trivial_content_with_embedding` 只回答**恒等重数**。这不是放宽：折叠到非 Γ 子群
k 点的块**不可能**含子群恒等表示——子群格平移 `t` 在任何 `q` 表示上作用为
`exp(-2πi q·t)`，而恒等表示作用为 1，因此共有不可约成分要求 `q` 落在子群倒格
（含 centering 消光）里。这些块被精确跳过，Γ 块仍走完整分解同一套 `build_block`
与冻结 child origin，因此结果精确。该路径的行标为 `identity_only`，与完整结果一样
逐条与存储表比较；TSV 的 detail 保留完整分解的原始错误。

报告分别列出：

- embedding 成功、歧义、无有效 setting 及其它错误；
- 每个源表示的完整结果（`full_success`）、恒等-only 结果（`identity_only`）、
  数据缺失、计算错误和 embedding 不可用；
- 94,271 条存储正项的比较状态，以及未存储的零项（`absent_zero` 再分引擎零项与
  几何零项）；
- Γ Frobenius 检查实际完成的条数、未计算条数；
- 5,756 条 `other_wave_vector_subduction` 的解析状态（是否属于冻结源、父群是否
  一致、频率非零）与 k 参数缺失。

每个总数都必须由互斥类别之和复原。重复源表示不能直接累加频率；同一源的不同
频率必须报错。当前 pinned 原始表独立清点得到 94,271 个不同的
`(isotropy ordinal, source irrep)`，没有重复行。

成功的完整分解经过维数、整数重数和逐操作特征标重建检查。这些计算检查及
恒等频率比较仍不足以固定非恒等 irrep 标签；标签依赖下面的官方坐标约定证据。

## setting 的来源与冻结

```bash
python3 scripts/audit_subduction_settings.py \
  --output /tmp/subduction-settings.jsonl --workers 8
python3 scripts/generate_subduction_settings.py --check
```

采集器从校验 SHA-256 的 `iso.zip` 使用官方程序及数据，每个 irrep 查询放在独立
临时目录，固定 `SET I ALL OR 1`。输出保留实际 ordinal、源 irrep 身份、方向标签、
官方基矢及 origin、查询失败原因和统计分母。官方空表、格式错误、origin 不一致
与不能得到整数 unimodular 换基的情况分别记录。产物是候选证据，不直接启用引擎。

冻结项按 ordinal 选择。换基由官方 conventional 基矢 `B` 精确反算：

```text
U = W P_parent (P_sub B)^-1
T = B^T
```

不能用“从母群里找到了一个同构的子群”替代这条有向基矢关系；也不能从操作集
匹配中任意挑选 child origin。冻结的 child origin 约定须保留来源，并逐操作验证。
冻结项若验证失败必须报错，不能回退到另一个候选。

`tests/subduction_settings.rs` 固定了 #43 的坐标轴陷阱：交换两个子群轴保持完整
操作集，却交换 GM3/GM4。源特征标在 x 法向镜面上的迹分别为 +1/−1，因此这种
错误可以通过所有群闭包检查，却必须被非恒等标签回归拦住。另有真正的 shear
矩阵见证，保证冻结机制不会被无意限制成 signed permutations。

## 全表结果（2026-09-22 第六轮，普通恒等分导表闭合）

> **历史小节**：下表是第六轮（R4 之前）的结果，`351,547 完整 + 14,713 恒等-only`
> 的分法只描述当时的状态；R4 三个批次之后 `identity_only = 0`，当前口径见本文末节
> 「R5：普通离散标量覆盖闭合」。

`audit_irrep_subduction --require-complete` 全表运行 452 s，退出码 0，判词
`VERDICT complete scope=global`：

| 项目 | 实际结果 |
|---|---|
| 冻结 embedding | 15,239 / 15,239（0 拒绝、0 歧义） |
| 所有标量 probe 请求 | 366,260 |
| probe 结果 | 351,547 完整分解 + 14,713 恒等-only = **366,260 全部有精确结果**；0 缺数据、0 计算错误、0 embedding 不可用 |
| 恒等正项 | **94,271 / 94,271 通过**；0 不匹配、0 不可用 |
| 未存储项 | 271,989 全部为 0（引擎零项 114,770、几何零项 157,219）；0 假阳性 |
| Γ Frobenius | 1,895 / 1,895 通过（strict 1,673、DistinctComponentSum 222） |
| 生产自检 | 维数、整数重数、逐操作重建、CIR 来源不匹配全部为 0 |
| 其它波矢记录 | 1,006 条记录 / 5,756 行全部解析到冻结源、父群一致、频率非零；**引擎 5,756/5,756 全部计算**（`--require-w-complete` 退出 0，`mismatched=0`、`engine_errors=0`）；数值同时由官方 `iso` live oracle 逐行对照，5,756/5,756 一致 |
| 不一致及计数错误 | 0（`hard_failures=0`、`accounting_violations=0`） |

同一张表在三个门禁下的当前结果（R0 起必须分别报告）：

| 门禁 | 退出码 | 判词/摘要 |
|---|---:|---|
| `--require-complete --require-w-complete` | 0 | `VERDICT complete scope=global gates=--require-complete,--require-w-complete full_decomposition=not_gated` |
| `--require-full-decomposition` | **0** | `full_decomposition: scope=global probes=366260 full_success=366260 identity_only=0 missing=0 error=0 uncomputed=0 incomplete=0 global=covered`，`VERDICT complete ... full_decomposition=complete` |
| 三个门禁同时 | **0** | `VERDICT complete scope=global gates=--require-complete,--require-full-decomposition,--require-w-complete full_decomposition=complete` |

即：**三个门禁同时 exit 0**。普通离散标量的完整分解在 R4 三个批次后是
366,260/366,260 = 100%，恒等分导表与 w 行的恒等重数各自闭合。缺口清零前本表第三行
曾是"完整分解缺口优先、exit 2"，那段历史口径保留在 R4 各批次小节里。

### R1/R2：表外目标由「子群操作 + 精确 q」现场构造

完整分解的目标不再必须绑定一条静态 `IrrepRecord`：

- 目标来源分两类，`FullStarTarget` / `SubductionTarget` 的 `ml` / `bc` / `row_ml` /
  `irnumber` 都是 `Option`：**存储成分**带冻结 CIR 身份；**构造成分**
  （`SubductionComponent::Constructed { q, index }`）只带精确 q 与构造表内序号，
  `None` 就是"没有来源标签"——不借 Γ 标签，也不填假 CIR 号。CDML/BC 命名与表示
  计算分开，构造成分在拿到有来源的标签映射前始终显示为未命名。
- 子群 #1 的任意 q：小群就是平移群，一维表示是 Bloch 相位
  `D_q(T_L) = exp(+2 pi i q.L)`，符号与存储行的 `bloch_phase` 一致；点先按子群倒格
  约化，因此 q 与 q+G 是同一个身份。存储行优先：只有该点确实没有 pinned 行时才构造。
- 复用同一条解块与重建流程：单位性、正交性、整数重数、维数和与逐操作重建仍由既有
  `solve_prepared_character_block` / `reconstruct` 检查；构造行没有跳过任何门禁。
- 实测（R2 交付时，`--require-complete` 全表，`hard_failures=0`、exit 0）：
  `full_success=353,382 identity_only=12,878`，恒等正项 94,271/94,271、Γ Frobenius
  1,895/1,895、w 行 5,756/5,756 全部不变。child-#1 的 **1,125 条记录 / 20,099 个
  probe** 全部完整分解，构造目标恰为清点里的 **3,992** 个缺失星（永久回归
  `tests/subduction_star_decomposition.rs::child_p1_records_decompose_every_scalar_probe_without_pinned_data`）。
- 余下 probe 仍返回 `MissingChildStarData`：它们的来源分类是 R3，可复用生成是 R4。
  逐子群的 **probe / 缺失星 / 全部折叠星**三个计数见
  [subduction-gap-census.md](subduction-gap-census.md) 末节，排优先级按 probe。
- **R4 批次 1（2026-09-22）**：构造判据由「子群 #1」改为「该星的小余群平凡」，
  并补上非平凡子群点群的**子群自星诱导**（`ConstructedStar`）。全表
  `full_success=357,033 identity_only=9,227`（+3,651 / −3,651）、`hard_failures=0`。
- **R4 批次 2b（2026-09-22）**：最后 58 组（52 组非退化 C2×C2 唯一二维不可约表示、
  6 组 D3 的规范化普通表示）由 `catalogue::projective_targets` 回答，带结构门禁与
  正交性门禁。全表 `full_success=366,260 identity_only=0 error=0`、两个门禁 exit 0、
  `VERDICT complete`；缺口集合为空，清点工具接受空输入并报零计数。
- **R4 批次 2a（2026-09-22）**：小余群非平凡但**一维投影特征标 catalogue 完整**
  （精确 cocycle 求解，解数 = `|P_q|`）的星由 `ConstructedLittleRep::Projective`
  现场回答。全表 `full_success=366,039 identity_only=221`（+9,006 / −9,006）、
  `error=0`、`hard_failures=0`；剩余 221 个 probe / 58 组需要二维投影表示（批次 2b）。
  范围、文件所有权、负例与一个已修 bug（原始/约化 q 混用导致复数重数）见
  [subduction-r4-batches.md](subduction-r4-batches.md)。
- 复核修复（`d140c10` 复核）：Γ 便捷入口的普通目标曾填占位 CIR 号 `Some(0)`，与
  full-star 入口对同一目标不一致；现在读 `record.source_identity()` 的真实编号，永久
  测试 `gamma_and_full_star_entries_agree_on_every_target_identity` 逐项比较两个入口
  的身份。`FullStarBlock::constructed_multiplicity` 改为只接受完整目标身份
  （`Constructed { q, index }`，q 按子群倒格规范化），不再按原始折叠坐标查询而返回
  假零；ordinal 1045 的两个构造成分身份 `(3/4,3/4,0)` 与 `(1/4,1/4,0)`、各重数 2
  已钉住。

关键的恒等-only probe（R2 前 14,713 条、R2 后 12,878 条、R4 批 1 后 9,227 条、
批 2a 后 221 条）中，**160 条是存储正项**（例如 SG 196 W1→#24 的
`W1`，存储频率 1）：本轮之前它们因另一条折叠星缺子群 k 数据而无法计算，现在由
恒等-only 路径逐条复现，0 不匹配、0 假阳性。永久测试
`identity_only_content_answers_probes_without_full_child_data` 与
`identity_only_content_agrees_with_the_full_decomposition_and_covers_the_pinned_set`
固定了：缺数据 probe 全部被精确回答（pinned 15 个，R4 批 1 关掉 5 个、批 2a 关掉
其余 10 个），另 2,075 个 probe 上恒等-only 结果与完整验证过的完整星分解完全一致。

官方采集器完成 4,777 个查询，输出 **0-based ordinal 0..15238** 的全部记录。
13,861 条候选通过 origin 精确相等与整数 unimodular 换基初检（其中 95 条 U
不是 signed permutation）；238 条基矢关系不符、945 条其余 origin 不符、195 条
官方空表，共 15,239 条。基矢不符项中有 112 条同时 origin 不符，因此不区分优先级
的 origin 不符总数是 1,057。这些候选数不代表引擎覆盖；未验证完整操作及坐标
约定的候选不会自动冻结。

所以两条轨道现在都在引擎侧范围内闭合：普通恒等分导表（15,239 条记录 / 94,271 条
正项 / 366,260 个 probe）无未支持项；`other_wave_vector_subduction` 的 5,756 行也
全部由引擎算出——k 域经 `data_little.txt` 的 `little_k` 与官方 `DISPLAY KPOINT`
双向确认，73 个源 irrep 的 little 群特征标由 `w_little_characters_data` 的 73/73
冻结表提供（不再依赖未解码的 `little_irr_full_matrices`）。
`--require-w-complete` 与 `w_scope` 行显式报告这条轨道，`mismatched` 与
`engine_errors` 都计入 `hard_failures()`。

后续扩充必须重新运行审计并更新实际覆盖。磁群、spinor 和离散子群 irrep 数据
未提供的 k 不会因这些工具而自动获得支持。

## R5：普通离散标量覆盖闭合（正式报告，2026-09-22）

**结论**：在固定语料上，普通离散标量的完整分解**已闭合**。三个门禁
（`--require-complete`、`--require-w-complete`、`--require-full-decomposition`）
同时运行时 **exit 0**、判词 `VERDICT complete scope=global
full_decomposition=complete`。本节是 R5 卡要求的"覆盖报告"；全表复现约 520–540 s
（本仓库记录的两次：521.7 s 与独立复核 541.0 s，同一份代码、数字逐项一致）。

**注意作用域**：`--parent`/`--ordinal` 的局部运行即使带上全部门禁也仍 exit 0，
判词里是 `global=not_established`；只有不带 scope 的全表运行才证明全局覆盖，
CI 里按 ordinal 循环不能替代它。

### 语料与分母（先声明范围，再谈百分比）

| 项 | 数量 | 说明 |
|---|---:|---|
| 冻结 isotropy 记录 | 15,239 | 全部 embedding 可用（`embedding_ok=15239`） |
| probe（普通标量） | **366,260** | 每条记录在普通标量探针上的展开，**分母就是这个数** |
| 其中完整分解 | **366,260（100%）** | `identity_only = missing = error = uncomputed = 0` |
| 存储恒等正项 | 94,271 | 全部比对通过，`mismatch = 0` |
| Γ Frobenius | 1,895 | 全部通过，`failures = 0` |
| w 频率行（别的波矢） | 5,756 | 恒等重数全部算出并匹配，`mismatched = engine_errors = 0` |
| spinor 记录 | 3,611 | **不在分母内**（普通标量分导明确拒绝） |

`identity_only` 路径（第六轮的恒等-only 入口）仍然存在，但在这份语料上已经**没有
probe 需要它**（`probe_identity_only = 0`），所以审计本身不再复算任何恒等-only 行；
这条入口的精确性现在由 `tests/subduction_identity_regressions.rs` 在 2,075 个 probe
上逐条对照完整分解（458 个正项）来保证。

### 每条完整结果被验证什么（不是只比恒等列）

先分清"谁来强制"和"谁能独立反驳"（独立复核 A 指出这里原先的说法过强，本节已改）：

| 不变式 | 强制点（真实门禁） | 审计字段 | 当前值 | 审计层的独立性 |
|---|---|---|---:|---|
| `Σ multiplicity × child_little_dim × child_star_size = parent_full_dim` | 引擎 `TotalDimensionMismatch` / `StarDimensionMismatch`（`subduction_star_decompose.rs` 构造块时即 Err） | `production_dim_mismatch` | 0 | **同义反复**：审计比较的是引擎同一批对象算出的量 |
| 重数为非负整数 | 引擎 `NonIntegralMultiplicity`（`subduction.rs`） | `production_integrality_mismatch` | 0 | 该字段实际查的是**维数配平**，不是整数性（口径已修正） |
| 逐操作重建子导字符 | 引擎重建检查（同一容差 `SUBDUCTION_TOLERANCE`） | `production_recon_mismatch` | 0 | 审计拿到的就是引擎已比较过的同一对向量 |
| 目标对上冻结 CIR 来源身份 | 引擎 `TargetSourceMismatch` | `target_source_unmatched` | 0 | 现有语料不可达（4105 个 CIR 号无重复），只可能由代码 bug 触发 |
| 来源号与标签读数一致 | 引擎 `TargetSourceMismatch`（by_label vs total） | `label_source_disagreement` | 0 | 现有语料不可达（同 SG 内 0 个重复 ml、0 条 compound 含 Γ 一维成分） |

它们的**共同效果**是实的：任何一项非零都会进入 `hard_failures`，并在**所有**门禁组合下
让运行 exit 1（负例测试 `every_counted_production_violation_is_a_hard_failure` 钉住接线；
复核 B 指出该测试只设置计数器、不驱动生产自增点，名字与注释已按此改口）。
但要诚实说明：这三项"审计侧"检查并不是对引擎的独立复算（引擎在更早的层已经用同一批
对象、同一容差比较过并以 Err 收口），另两项在现有语料上不可达。因此**真正独立的证据是
引擎的 Err 路径本身**（`probes.error` → hard failure → exit 1）以及下面这些外部 fixture，
而不是这张表里的零值。

**对 R4 新增的「构造目标」这五项还要再退一步（复核 B）**：构造目标没有冻结 CIR 来源
（`irnumber = None`），所以"目标对上来源身份"与"来源号/标签一致"两项对它是**恒真**的；
又由于每个空间群都有平凡 Γ 行，子群 Γ 星永远命中 stored 行，恒等频率比较**从不**看到
构造星（复核 B 在 6 个见证上下文里数到 44 个构造块，`gamma_constructed = 0`）。构造路径
因此只剩引擎自身的维数/整数性/重建检查，加下一条的外部对照。

另有独立的来源矩阵 fixture 继续参与门禁（`tests/subduction_star_source.rs` 的
1,232 次字符比较、`tests/subduction_compound_stars.rs` 的复成分与 k/-k 对照），
它们不经过诱导算法，因此不是自证。**非平凡标签可以置换**：恒等频率吻合不能替代
规范标签检查，所以 `label_source_disagreement` 与 `target_source_unmatched` 单独计入
`hard_failures`，而不是只看恒等列。

### 构造目标的外部对照（唯一一条，计数已钉死）

`the_catalogue_reproduces_pinned_little_group_characters` 在**离散** pinned k 点上
（归档有字符表、引擎已自行验证过的小群行）逐操作比较 catalogue 与 pinned 字符：
**1,328 条记录 / 7,578 个操作 / 94 条走二维表 / 1,660 条延后**，四个数在测试里都是
`assert_eq!` 钉值（复核 B 之前只有下限断言，静默缩水不会被发现）。它不经过缺口数据，
所以不循环；范围只有"离散点上的非平凡小余群"，不覆盖参数化折叠点本身。

### 缺口清点（R5 卡要求的口径修正）

清点工具 `examples/census_subduction_gaps.rs` 仍然是"从**新审计**生成剩余缺口"，
replay 校验没有被关闭：它现在拒绝非 `identity_only` 行、仍校验首个缺失星与引擎错误
一致。覆盖闭合后审计里没有 `identity_only` 行，工具按**空 manifest + 零计数**返回
（回归测试 `a_closed_audit_replays_as_an_empty_manifest`），而不是静默失败或伪造缺口。
清点结果（批 2a 后）为 `records=83 probes=221 stars=331 missing_stars=326`，
批 2b 后清零。

**复核 A/B 之后加固（2026-09-23）**：空 manifest 的两种来源（审计真的闭合 vs 传错/
截断/格式漂移的文件）此前不可区分（A 提出，B 用 `identity-only;` 漂移实跑复现）。
现在 `read_requests` 用**精确 token**（而不是子串）识别 marker，摘要行给出
`audit_rows`、`probe_rows`、`identity_only_rows`、`unanswered_probe_rows` 与 `closed`
（不再打印恒为 0 的 `replay_errors`/`dimension_errors` 常量），没有识别到 marker 时
打警告，并新增 `--require-empty`：要求审计**至少有一条 probe 行、每条 probe 行都被
引擎回答、且没有 identity-only 行**，否则 exit 1。在 R5 审计上：

```text
$ census_subduction_gaps target/r5_audit.tsv --require-empty
records=0 probes=0 stars=0 missing_stars=0 constructed_stars=0 reachable_stars=0
audit_rows=389150 probe_rows=366260 identity_only_rows=0 unanswered_probe_rows=0 closed=true   （exit 0）
```

残留限制写在工具文档里：**被改名**的 marker 若其余字段仍合法，本工具无法与"已闭合"
区分，所以它打印两个计数而不是一个光秃秃的 0。replay 的 6 条拒绝路径现在各有一条
负例（`the_replay_rejections_are_all_reachable`），marker 漂移另有一条
（`a_renamed_identity_only_marker_is_not_silently_dropped`）。

### 复现

```bash
CARGO_TARGET_DIR=$PWD/target cargo run --release -p cryspglib \
  --example audit_irrep_subduction -- \
  --require-complete --require-w-complete --require-full-decomposition \
  --output target/r5_audit.tsv
```

实测（2026-09-22，约 520 s，exit 0）：

```text
coverage: embedded_records=15239/15239 positive_stored_compared=94271/94271
  probe_full_success=366260/366260 probe_identity_only=0/366260 probe_answered=366260/366260
  absent_zero=271989 frobenius_evaluated=1895/1895 w_computed=5756/5756
production_checks: dimension_mismatch=0 integrality_mismatch=0 reconstruction_mismatch=0
  target_source_unmatched=0 label_source_disagreement=0
w_scope: rows=5756 computed=5756 uncomputed=0 mismatched=0 engine_errors=0
VERDICT complete scope=global
  gates=--require-complete,--require-full-decomposition,--require-w-complete
  full_decomposition=complete
```

### 边界（这些不在 100% 里）

* spinor / 双群分导：明确拒绝，单独扩展；
* 参数化 k 源的**完整**分解（R6）：当前只算恒等重数（5,756 行），而且那 5,756 行是
  **冻结在官方程序的参数约定 `t = 1/4`**（`subduction_star_decompose.rs` 的
  `LINE_PARAMETER`）上的频率，不代表任意 t 都成立；
* 磁共表示分导（R9–R11）：磁表只提供候选记录查询，不是经过验证的磁嵌入；
  230 个 Type-II/grey UNI 在磁表中没有记录，是数据来源边界；
* 官方方向 descriptor（R7）与正式 API（R8）：`dim = 2`/`dim ≥ 4` 仍是内部记法，
  `irrep::subduction` 仍是 `#[doc(hidden)]`。

因此"100%"只等于**这份固定 ordinary scalar 语料的完整分解**，不是项目完成度。

### 对抗性复核（2026-09-23，两个独立 reviewer）

复核 A（覆盖声明与门禁真实性）与复核 B（实现与回归质量）都是只读审查，结论均为
**可接受、无阻断项**；B 用不依赖本 crate 的路径（自建 cocycle/扭群代数、直接读
`iso.zip` 重算分母 366,260、独立重跑三门口禁 535 s）复现了本节所有头条数字。已处理：

| 发现 | 处置 |
|---|---|
| A-P1：`13345 W1–W5`「永久反例」的说法过期 | 文档改写（`e9e8e8e`）；该上下文现已完全回答 |
| A-P2：五项生产检查里有三项同义反复、两项在语料上不可达 | 本节独立性表逐项标注（不再称"load-bearing"） |
| A-P2：`CLAUDE.md` 基线数字过期 | 基线块更新（18/4、161、三门禁、520–540 s） |
| B-P1：C2×C2 族"唯一二维"的**证明前提** `omega(g,g)=±1` 为假（实测含 1/6…5/6，child #43 有 `omega(g,g)=i`） | 结论不变、论证改写（`subduction_catalogue.rs` 与批次卡）；`1d1b95e` 提交信息以文档为准 |
| B-P2：构造目标在五项检查里结构性只剩一项有效；恒等比较永不看到构造星 | 本节新增「构造目标的外部对照」小节，并说明原因 |
| B-P2：构造路径唯一的 pinned 对照只有下限断言 | 四个计数改为 `assert_eq!` 钉值（1,328 / 7,578 / 94 / 1,660） |
| B-P2：端到端 fail-closed 负例被三个批次替换殆尽 | 补 `an_out_of_scope_co_group_still_reports_missing_child_star_data`（真实嵌入 + 手工构造越界星）与 solver 三条负方向（含 `MAX_GRID`） |
| B-P2：census 空输入放宽后无法区分"闭合"与"格式漂移" | 精确 token + 计数 + 警告 + `--require-empty` + 6 条拒绝路径负例（见上） |
| B-P2：`every_production_check_violation_is_a_hard_failure` 只测接线 | 改名 `every_counted_production_violation_is_a_hard_failure` 并在注释里写明它不驱动生产自增点 |
| B-P3：`ConstructedStar::dimension()` 硬编码 1（二维族会答错） | 维数改为显式携带并加回归 |
| B-P3：`subduction_catalogue.rs` 的重复 doc 行、`4 != order` 死代码、多处过期的 `#[allow(dead_code)]`（实际都在调用链上） | 全部删除，严格 clippy 仍零警告 |
| B-P3：census 错误路径 `{other:?}` 单行 >6 KB | 改为报告块数与错误本身 |
| B-P3：SG 38/40 的 pinned 行对小群外操作存字面 0 | 写入 `subduction-conventions.md` §14 的存储行契约（"row = 每个列出操作的特征标"不成立） |
| B-P2（潜在，未触发）：`stored_child_components_at` 对无法展开的记录直接 `?` | 保留大声报错并写明理由（伪装成"缺数据"更糟）；语料上 8,388 条记录 0 条不可展开 |

**复核 A 列出的"无法验证"项，本轮主线程已补两项**（记录在
`target/r5_review_a_unverified.md`，输入全部是仓库内冻结文件 + tracked 脚本）：

* R3/R4 中间离线数字：重跑 `scripts/classify_subduction_gap_sources.py` 得到
  `899 = 215 解析 + 684 参数化`、`684 组 / 123 子群 / 81,576 矩阵元`、批 1 的
  `215 组 / 41 子群 / 4,212 星 / 3,651 + 119 probe`、批 2a 的
  `58 = 52 + 6 / 221 / 83 / 29 / 331 = 326 + 2 + 3`，全部与文档一致；且重跑产物与冻结的
  `target/r4_groups.tsv` **SHA-256 完全相同**（`acec85f0…`）。
* 73 张 w 特征标表：重跑 `scripts/freeze_w_little_characters.py`（调用随包官方 `iso`）
  得到 `sources solved: 73 | failures: 0`，生成的 Rust 表与提交的
  `src/irrep/w_little_characters_data.rs` **逐字节相同**。
* 仍未验证（如实保留）：几何 oracle 只有 22 组 (SG, irrep) / 62 行抽样（占 15,239 行的
  0.31%），扩大它需要为 4,777 个普通 irrep 各起一次官方 `iso`，是独立工作量，与本节的
  100% 结论无关，**不因为上面两项已复现而改写**。

### 第三方复核（2026-09-23，公开仓库浏览 + 隔离复现）

第三个独立 reviewer 从公开仓库审阅了 `0386cfa`..`ff2e4a4` 区间（本节的 R5 收口提交
在 `ff2e4a4` 之后，部分发现针对的是生成管线而不是 R5 数字），给出 5 条发现 + 3 条
提醒，**没有**质疑 R5 的 366,260/366,260 与恒等正项 94,271/0 不匹配。逐条处置
（全部在主线程复现后修）：

| 发现 | 复现结果（主线程） | 处置 |
|---|---|---|
| P1：`scripts/generate_subduction_settings.py` 仍是 69+6 条、五字段的旧生成器，`--write` 会覆盖 15,239 条六字段的正式模块 | **属实，且已有后果**：`scripts/test_subduction_settings.py` 因此**当前就是红的**（`parse_committed` 解析不了六字段，2 个 setUpClass error） | 旧生成器的 `parse_committed` 改为六字段；`--write` 改为**拒绝执行**（指向 `scripts/task9/build_table.py`）；`--check` 改为校验它自己的 75 条遗留记录 + 模块覆盖全部 15,239 个 ordinal；该测试重写后 10 passed；数据模块头与生成器头同步改指 task9 管线 |
| P1：嵌入验证缺 `L_H ⊆ L_G`，`det U = ±1` 不够 | **属实**：`probe_embedding` 对 SG 1 自嵌入 + `U = diag(2,1/2,1)` 返回 `Ok`（已用公共 API 复现） | `validate_candidate` 增加包含性检查；全表 15,239 行回归 + 负例各一条；审计 embedding 仍 15,239/15,239（未误伤任何冻结行） |
| P2：`IrrepSubduction`/`FullStarSubduction` 只存 setting 分子，丢失分母（ordinal 26 分母为 2） | **属实**：两处只复制 `embedding.setting()` | 两个结果类型新增 `setting_denominator()`；ordinal 26 上的回归钉住与嵌入一致 |
| P2：`probe_subduction_settings.rs` 把计算失败写成 `trivial=0`，与"真的 0"混淆 | **属实**（`Err(_) => 0`） | 三态分开（`<n>` / `?` / `error` + stderr 原因），profile 失败输出 `profile=error`；两条单测 |
| P2：`scripts/task9/build_table.py` 先写文件后统计、不去重、空输入也能 exit 0 | **属实** | 重写为"先验证后写"：重复 ordinal 报错、必须覆盖既有模块的 ordinal 集合（或 `--expected`）、未覆盖记录按 status 报告、`os.replace` 原子替换、新增 `--check` 与显式 `--partial`；新增 `scripts/test_build_table.py`（5 条回归）；`--shifts/--derived` 接受管线实际写出的三种形状并可重复 |
| 提醒：`--parent 2 --ordinal 0` 选不中记录却报 `clean` | **属实**：`VERDICT clean ... probe_full_success=0/0`、exit 0 | 审计新增 `empty scope` 错误（非零退出）+ 回归；合法 scope 仍 exit 0 |
| 提醒：`scripts/task9/README.md` 路径多退一级、12 个脚本硬编码 `/home/liuyichen/...` | **属实** | README 路径改为从 `target/task9` 出发的 `../../…`；12 个脚本改为从 `__file__` 推导 `REPO`；补上只在 `target/task9/` 里存在、被三个 tracked 脚本 import 的 `derive_shift.py`（现在已入库），全部脚本可 import |
| 提醒：审计文档仍把"exit 2 / identity_only=14713"写成当前状态 | 属实 | 改为显式历史 + 当前状态（R4 后 `identity_only=0`、该门禁 exit 0），并写明 `empty scope` 的语义 |

**跟进复核（同一 reviewer，针对 `c37d3b9`）指出两处仍可绕过的门禁，均已补**：

1. `build_table.py` 在 `--out` **尚不存在**时没有预期全集，完整性检查被跳过，
   空输入会写出零条目表并 exit 0（`test_an_empty_input_cannot_write_a_table` 只覆盖了
   "已有非空基线"这一半）。现在预期全集只能来自 `--expected`、已存在的 `--out`、或
   **tracked 的 `--baseline` 模块**（默认 `src/irrep/subduction_settings_data.rs`）；
   三者不一致即报错，全部缺失时拒绝写入（除非显式 `--partial`）。reviewer 给出的最小
   复现命令现在 exit 1、不创建输出；`scripts/test_build_table.py` 加到 7 条
   （新增"新输出路径必须知道全集"与"out 与 baseline 必须一致"）。
2. `generate_subduction_settings.py --check` 的覆盖检查是**数量相等 + 逐条身份**，
   重复 ordinal（如 `[0,1,1]` 对 3 条记录）能顶替缺失项通过，而 `compare()` 以
   ordinal 为键会把重复合并。现在抽出共享实现 `check_ordinal_coverage`（要求
   ordinal 序列严格等于 `range(len(records))`，报出缺失/重复/越界），
   **在线 `--check` 与离线测试调用同一个函数**，避免"测试严格、生产宽松"的分叉。
   reviewer 的四个回归用例（`gate_regression_tests.py`）现在 4/4 通过。

**关于正式生成入口的统一**（发现 1 的建议）：现在只有一个 writer——
`scripts/task9/build_table.py`；`--check` 是它的门禁。诚实说明一条边界：
**历史那份表的逐字节重生成无法只靠 tracked 输入完成**，因为最后的 child-shift 集合是
多次增量修补文件（`recorded_shift.json`、`final_shifts.json`、`one27_out.json` …，
合计 2,703 条非零 shift）的合并，这些证据在未跟踪的 `target/task9/` 里；因此本仓库
**不声称**该表的字节级可复现，只声称：模块内容离线可核（每个 pinned ordinal 一条、
身份相符、`|det U| = 分母³`、shift 最简）、引擎可嵌入全部 15,239 行（审计
`embedded_records=15239/15239` + 全表 `L_H ⊆ L_G` 回归）、以及 75 条遗留记录**在线**
经官方 oracle 复推一致（`generate_subduction_settings.py --check` 刚刚实跑通过：
`checked 75 legacy entries ... the committed module covers all 15239 pinned records`）。
未来重建的口径见 `scripts/task9/README.md` 的「Provenance status」小节。
