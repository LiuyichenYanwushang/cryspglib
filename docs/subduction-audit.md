# 完整分导的全表审计（任务 9）

任务 9 的工具分别检查几何来源和生产 API 的计算覆盖。**收集到官方基矢、恒等项
吻合、完整分解成功是三个不同的结论**；任何一个都不能替代另外两个。

## 生产 API 审计

在工作区根目录运行：

```bash
CARGO_TARGET_DIR=$PWD/cryspglib/target cargo run --release -p cryspglib \
  --example audit_irrep_subduction -- --output /tmp/subduction.tsv
```

`--parent 221` 或 `--ordinal 12400` 限定诊断范围。退出码 2 表示所选门禁覆盖不完整，
退出码 1 表示发现不一致，退出码 0 只表示未发现不一致。两个门禁分开：

- `--require-complete` 要求所选范围内**普通恒等分导表**闭合：15,239 条记录全部冻结
  embedding、94,271 条存储正项全部复现、366,260 个标量 probe 全部有精确结果
  （完整分解或恒等内容）、几何与 Frobenius 检查无未计算项。
- `--require-w-complete` 另外要求 5,756 条 `other_wave_vector_subduction` 也被
  **引擎**计算。目前 pinned 数据回答不了：这 73 个“别的波矢”irrep 在
  `data_irreps.txt` 里只有 `irrep_w_label/_space_group/_dimension/_type` 四张表
  （无 k、无特征标）。第十七轮已把这 73 个源的 k 域全部解出（参数化点/线/面/一般
  位置；方向记录在母群 primitive 倒格基里，官方 `DISPLAY KPOINT` 用 conventional
  帧打印），所以 **k 参数已经不是缺口**；缺的是 little 群特征标
  （`little_irr_full_matrices` 的 2.22M 整数编码尚未解出）。该开关因此仍是引擎侧
  门禁，全表运行退出 2，并打印
  `w_scope: rows=5756 computed=0 uncomputed=5756 reason=...`。

### 其它波矢行的独立官方对照

引擎侧尚未计算这 5,756 行，但它们的**数值**已有可复算的官方对照：
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
| 普通恒等分导（15,239 记录 / 94,271 正项 / 366,260 probe / Γ Frobenius 1,895） | **cryspglib 引擎计算** + 几何与零项检查 | 范围内闭合，0 未支持 |
| 其它波矢 w 行（1,006 记录 / 5,756 行） | 官方 `iso` live oracle 逐行（键 = 子群号 + Dir） | 逐行一致；**引擎未计算**，`--require-w-complete` 退出 2 |

这 5,756 行里 **348 行（194 条记录，6.0%）**引用 SG 202/203/209/210 的 `DT3`/`DT4`
—— 这 8 个源的特征标目前解不出来（Γ 兼容表只能定它们的和，线上其它特殊点的
`SHOW CHARACTER` 又被星污染，pinned 通道也区分不了它们：两个块逐行相同）。
其余 **5,408 行（93.99%）**的源已有冻结特征标
（`src/irrep/w_little_characters_data.rs`，65/73，锚点回归 `tests/w_little_characters.rs`），
下一步是引擎侧的 Mackey/特征标求和。

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
| 其它波矢记录 | 1,006 条记录 / 5,756 行全部解析到冻结源、父群一致、频率非零；k 域已解出但**引擎 0 行可计算**（`--require-w-complete` 退出 2）；数值由官方 `iso` live oracle 逐行对照，5,756/5,756 一致 |
| 不一致及计数错误 | 0（`hard_failures=0`、`accounting_violations=0`） |

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

所以普通恒等分导表（15,239 条记录 / 94,271 条正项 / 366,260 个 probe）**引擎侧
范围内闭合**：范围内未支持项为 0。范围之外的另一条轨道是
`other_wave_vector_subduction` 的 5,756 行：它们的 k 域已经解出（参数化 k 域，
`data_little.txt` 的 `little_k` + 官方 `DISPLAY KPOINT` 双向确认），数值也与官方
`iso` 的 live oracle 逐行一致，但**引擎还不能算**——73 个源 irrep 的 little 群
特征标不在 pinned irrep 表里，`little_irr_full_matrices` 的编码未解出。
`--require-w-complete` 与 `w_scope` 行持续显式报告这一条，不会被静默算作已完成。

后续扩充必须重新运行审计并更新实际覆盖。磁群、spinor 和离散子群 irrep 数据
未提供的 k 不会因这些工具而自动获得支持。
