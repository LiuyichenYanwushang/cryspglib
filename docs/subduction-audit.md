# 完整分导的全表审计（任务 9）

任务 9 的工具分别检查几何来源和生产 API 的计算覆盖。**收集到官方基矢、恒等项
吻合、完整分解成功是三个不同的结论**；任何一个都不能替代另外两个。

14,713 个恒等-only probe 的逐星缺口清点见
[subduction-gap-census.md](subduction-gap-census.md)：21,136 个缺失星，
按子群/k-star/setting 归并为 989 组。清点不改变下述生产覆盖数字。

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
  因为回答恒等重数是另一个问题）。当前全局基线在新门禁下 **exit 2**，摘要行
  `full_decomposition: ... identity_only=14713 incomplete=14713 global=covered`
  明确报告 14,713 个缺口；旧的 `--require-complete`（以及 w 门禁）仍 exit 0。
  这条门禁的存在就是为了防止再次把"恒等项通过"当成"完整分解完成"。
  范围语义：`--parent`/`--ordinal` 限定运行时只报告该范围的完整分解状态，
  判词带 `global_coverage=not_established`；局部全过不会被表述成全局覆盖完成。
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
| 普通恒等分导（15,239 记录 / 94,271 正项 / 366,260 probe / Γ Frobenius 1,895） | **cryspglib 引擎计算** + 几何与零项检查 | 范围内闭合，0 未支持。366,260 个 probe 结果分两类：**351,547 个完整分解 + 14,713 个仅恒等重数**（例：ordinal 13345 `W1` 的完整分解仍 `MissingChildStarData`，其恒等重数由 `trivial_content_with_embedding` 精确回答） |
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
**不能跳过完整分解调用**。ordinal 13345 的 W1–W5 仍是永久反例：它们没有恒等项，
却缺少完整分解所需的子群 k 数据；把它们仅记成零项会错误宣告完整覆盖。

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
| `--require-full-decomposition` | **2** | `full_decomposition: scope=global ... identity_only=14713 incomplete=14713 global=covered`，`VERDICT incomplete ... gates=--require-full-decomposition` |
| 三个门禁同时 | **2** | 同上；完整分解缺口优先于其它门禁的通过 |

即：恒等分导表与 w 行的**恒等重数**已经闭合（旧两个门禁 exit 0），但**普通离散
标量的完整分解**仍是 351,547/366,260 = 95.98%，缺口 14,713 条，R1 起的里程碑按
`docs/subduction-next-milestones.md` 逐批补齐；在缺口清零前，新门禁一直退出 2，
不得用恒等项通过代替完整分解验收。

关键的 14,713 条恒等-only probe 中，**160 条是存储正项**（例如 SG 196 W1→#24 的
`W1`，存储频率 1）：本轮之前它们因另一条折叠星缺子群 k 数据而无法计算，现在由
恒等-only 路径逐条复现，0 不匹配、0 假阳性。永久测试
`identity_only_content_answers_probes_without_full_child_data` 与
`identity_only_content_agrees_with_the_full_decomposition_and_covers_the_pinned_set`
固定了：15 个 pinned 缺数据 probe 全部被精确回答，另 2,060 个 probe 上恒等-only
结果与完整验证过的完整星分解完全一致。

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
