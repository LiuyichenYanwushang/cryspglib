# R4 批次卡（按可验证的目标族扩展）

R4 的规则来自 [subduction-next-milestones.md](subduction-next-milestones.md) §R4：一次只做
一个明确目标族，每批先离线证明特征标、代表元与相位，再接生产求解，并给出该批
**原缺口 probe 的完整入口成功数 / 仍缺数据数 / 新错误数**；既有成功集合不得退化。
批次顺序以 R3 的可行性为准（[subduction-gap-sources.md](subduction-gap-sources.md)）。

全表基线（R2/R3 后）：`full_success = 353,382`、`identity_only = 12,878`
（`probe = 366,260`）、`hard_failures = 0`。

## 批次 1：平凡小余群的构造目标（215 组解析路线）

**状态：已交付（2026-09-22）。** 实测 `full_success = 357,033`
（+3,651）、`identity_only = 9,227`（−3,651）、`hard_failures = 0`；
存储恒等正项 94,271/0 不匹配、Γ Frobenius 1,895/1,895、w 行 5,756/0 错误不变。
重跑清点后剩余 9,227 个 probe、12,932 个缺失星（原 12,878 / 17,144），
R3 分类器在新缺口上只剩 **684 个 `parameterized_source`**（0 个解析、0 个无源），
即批次 1 恰好关掉了 R3 判定为解析的那一批。

### 范围（R3 输入，可复算）

R3 的 899 个缺口组里 `analytic_general_position = 215`（41 个子群），它们的**每个**
星臂小余群阶都是 1：小群只剩平移，因此 q 上的目标就是一维 Bloch 相位
`D(T_L) = exp(+2πi q·L)`，不需要任何归档字符数据。

| 项目 | 数量 |
|---|---:|
| 组 | 215（41 个子群） |
| manifest 行（缺失星） | 4,212 |
| 只含解析星、本批可完全闭合的 probe | 3,651 |
| 同时含参数化星的混合 probe（本批只闭合其解析星） | 119 |

命令：

```bash
python3 scripts/classify_subduction_gap_sources.py target/r12_gaps.tsv > target/r3_groups.tsv
```

### 机制

R2 只对子群 #1 构造目标（`child_sg != 1` 直接返回空，保持 `MissingChildStarData`）。
本批把判据从「子群号」换成**几何事实**：子群自身 data-Hall 操作里，固定在 q（模子群
原胞倒格）的旋转只有恒等 → 构造一维 Bloch 相位目标；否则仍返回空，保持
`MissingChildStarData`，由批次 2 处理。

两条**星级**规则是本批实现中实测踩出来的，写在这里以免重犯：

1. **构造是星级回退，不是点级回退**：一个物理子群 irrep 活在**整个星**上。星里只要
   任何一条臂命中了 pinned 行，就必须只用 pinned 行（否则同一个 irrep 会以两个身份
   出现，求解器的正交性检查直接报错：ordinal 1007 的 `V1`/`L1` 就是见证）；
2. **一个星只构造一次**，落在该星规范化排序的第一个点上：平凡小余群在整个星上只有
   一个小群 irrep，逐点各构造一个会得到两个相同字符的行（ordinal 1007 的
   `(0,1/2,1/2)` 与 `(0,3/2,1/2)` 就是见证）。

两处必须同时成立才能算「目标完整」：

1. **小群块（`H_q`）**：`transported_values` 用构造成的相位表示在子群操作上取值（已有）；
2. **子群全星**：重建阶段把目标按**子群自己的星**诱导后逐子群陪集代表元比较。子群点群
   非平凡时星大小 > 1，因此构造目标也必须有诱导星（本批新增 `ConstructedStar`），
   不能像 #1 那样把「小群特征标 = 全星特征标」。

### 文件所有权

- `src/irrep/subduction_star.rs`：`ConstructedLittleRep` 移入此处；新增
  `LittleCharacter`（存储行 / 构造表示的统一字符来源）与 `ConstructedStar`；
  `arm_character`、`induced_component_character`、`induced_character` 改为走该来源，
  存储行路径行为不变。
- `src/irrep/subduction_star_decompose.rs`：`constructed_child_components_at` 的判据、
  `ChildStarEvaluator::Constructed` 持有构造星、单元测试。
- `tests/subduction_constructed_stars.rs`（新）：端到端与负例。
- `docs/subduction-r4-batches.md`、`docs/subduction-gap-sources.md`、`CLAUDE.md`。

### 验收（实测）

- `cargo test --release -p cryspglib --tests --doc` 全绿（lib 392、integration 160、
  doctest 27；audit example 17、census example 2）；严格 all-target clippy 通过。
- 全表审计（`--require-complete --require-full-decomposition`，518 s）：
  `full_success=357033 identity_only=9227`、`hard_failures=0`、
  恒等正项 `94271/94271 passed, mismatch=0`、Γ Frobenius `1895/1895`、
  w 行 `5756 computed / 0 mismatched / 0 engine_errors`、
  `production_checks` 五项全 0、`accounting_violations=0`。
  `--require-complete` 仍 exit 0；`--require-full-decomposition` exit 2
  （`incomplete=9227`，R5 的目标）。
- 新增回归：
  - `tests/subduction_constructed_stars.rs`（4 项）：ordinal 13345（SG 225 → #8）
    31/31 完整且恒等重数与 pinned 频率一致、两臂星的诱导（`block_dimension =
    star_size × little_dimension`）、ordinal 13346 的 10 个参数化星仍是
    `MissingChildStarData`、以及「有 pinned 行的星绝不再构造」；
  - `src/irrep/subduction_star_decompose.rs`：
    `a_constructed_star_induces_over_the_child_star`（子群 #2，手算
    `χ(E)=2`、`χ(T_(1,0,0))=0`、`χ(T_(0,1,0))=−1`、`χ(−I)=0`）；
  - `examples/census_subduction_gaps.rs`：清点的可达性模型改为「无 pinned 行的星，
    平凡小余群 → `constructed_trivial_co_group`，否则 → `missing_discrete_scalar_data`」，
    两条测试分别钉住 13345（0/3/0）与 13346（0/1/4）；
  - `examples/audit_irrep_subduction.rs`：13345 由「五个 W 恒等-only」改为
    「31/31 完整、全门禁 clean」，缺口的微型基线移到 13346；
  - `tests/subduction_identity_regressions.rs`：冻结上下文里 15 个缺数据 probe → 10 个，
    比较数 2,060 → 2,065、非 Γ 1,522 → 1,527。
- 负例：非平凡小余群且无 pinned 行的折叠星**仍然**是 `MissingChildStarData`
  （ordinal 13346 的 `B` 线四星，R3 的 `parameterized_source`），批次边界未被越过。

### 尚未做（批次 2 入口）

- 684 个 `parameterized_source` 组（123 个子群、9,227 个 probe）：需要把归档参数域
  代入求解参数、取出该 q 上的小群字符/矩阵并验证适用域与目标完整性；
  R3 已确认 684/684 的 PIR 矩阵块可读（81,576 个矩阵元）。
- 参数化来源可复用则冻结最小数据，运行时保持纯 Rust（R4 规则）。
- `little_k`（w 源数组）与 PIR `k_arms` 是不同数组，各自的坐标约定要分别验证，
  不能互相套用。

## 批次 2 侦察（2026-09-22，只读，未改引擎）

目标是把 R3 剩下的 **684 组 / 123 个子群 / 9,227 个 probe** 从"来源候选"变成可求值的
目标。侦察结论（复算脚本 `target/r4_batch2_{recon,structure,classes}.py`，全部只读）：

1. **归档记录的小群覆盖目标小群**：684/684 组的 `matched_*` 记录的 `operations` 旋转集
   都包含 q 处**实际**小余群的全部旋转（即记录的通用域小群不比 q 处小），所以
   "记录小群 vs 特殊点小群"的失配不是普遍现象。
2. **目标小余群很小**：|P_q| ∈ {2,3,4,6}（2:452、4:222、3:4、6:6）；除 6 组为
   非交换 D3 外全部交换。
3. **catalogue 由精确 cocycle 决定**：因子系统 `omega_ij = exp(2πi q·L_ij)`（`L` 取
   未约化乘积）给出 ω-正则类数 k，配合 `Σ dim² = |P_q|` 唯一确定小群不可约表示的
   维数多重集：|P|=2→{1,1}、|P|=3→{1,1,1}、|P|=4 且 k=4→{1,1,1,1}、
   |P|=4 且 k=1→{2}、|P|=6→{2,1,1}。

| \|P_q\| | #ω-正则类 | 星大小 | 组数 | 需要的小群 irrep |
|---|---:|---:|---:|---|
| 2 | 2 | 1/2/3/4/6 | 138/252/4/54/4 = 452 | 两个一维（χ(g)=±1 或 ±i） |
| 3 | 3 | 1/2 | 2/2 = 4 | 三个一维 |
| 4 | 4 | 1/2/4 | 58/90/22 = 170 | 四个一维 |
| 4 | 1 | 1/2/4 | 18/32/2 = 52 | **一个二维** |
| 6 | 3 | 1/2 | 4/2 = 6 | 一个二维 + 两个一维 |

4. **为什么走"精确 cocycle 求小群 irrep"而不是冻结归档参数表**：归档记录按 **k 域**
   索引，而目标是域上的**特殊参数点**——通用小群可能比 q 处小，域的星在 q 处会合并
   （child #3 的 `V1VA1`、child #25 的 `LD1LE4`：记录 dim 2 = 通用星 2 × 小维 1，
   但在 q 处星塌成 1、小群升到 4）。要复用归档字符就必须逐域重建"域 irrep → 点 irrep"
   的对应，而 R3 已明确把这条对应留作未验证项。目标 catalogue 本身则**完全由 q 处的
   精确 cocycle 决定**，不需要任何归档字符数据。

### 批次 2a：一维投影特征标 catalogue（已交付，2026-09-22 第三轮）

**状态：已交付。** 全表审计（`--require-complete --require-full-decomposition`，560 s）：

| 项目 | 批次 1 后 | 批次 2a 后 |
|---|---:|---:|
| `full_success` | 357,033 | **366,039** |
| `identity_only` | 9,227 | **221** |
| `error` / `hard_failures` | 0 / 0 | **0 / 0** |
| 恒等正项 | 94,271/0 不匹配 | 94,271/0 不匹配 |
| Γ Frobenius | 1,895/1,895 | 1,895/1,895 |
| w 行 | 5,756 / 0 错误 | 5,756 / 0 错误 |
| `production_checks` | 全 0 | 全 0 |

覆盖率 366,039/366,260 = **99.94%**；`--require-complete` 仍 exit 0，
`--require-full-decomposition` exit 2（`incomplete=221`）。清点后剩余
**221 个 probe / 83 条记录 / 29 个子群 / 326 个缺失星**，重新分类后**全部**是
`parameterized_source`：**58 组**（星阶 4 的 52 组 + 星阶 6 的 6 组），与侦察的 2b
预测完全一致。

**实现**（`src/irrep/subduction_catalogue.rs` + `subduction_star.rs` 的
`ConstructedLittleRep::Projective` + `constructed_child_components_at` 的新分支）：

* catalogue 的因子系统、常量与取值共用**同一个约化后的点**；
* 只有解数恰为 `|P_q|` 时才使用；否则保持 `MissingChildStarData`（fail closed）；
* 一维 character 的身份是 `Constructed { q, index }`，`index` 按解序编号。

**本轮修掉的 bug（教训）**：第一版接线把常量建在**原始**折叠点上、却用**约化**点求值，
30 个 probe 因此得到复数重数（`1±i`，ordinal 14090，SG 226 W5 → #98）。原因：重建出的
小群操作**不是**其代表元的格平移（平移可以带 1/4），两半相位必须共用同一个 `q`；
混用会让相位差一个非整数。修复后 14090 变成 25/25 完整、`VERDICT clean`。
值得记下的是引擎当时是 fail-closed（报错而不是给错值），所以错误可见、没有污染结果。

**证据**：

* `the_catalogue_reproduces_pinned_little_group_characters`：在**离散** pinned k 点上
  对照引擎已验证的小群字符行，**1,176 条记录 / 6,318 个操作全部命中**
  （另 1,812 条属更高维批次，跳过），2.8 s；
* `a_two_fold_co_group_has_two_characters`：手算 C2（φ=0 → ψ∈{0,1/2}；φ=1/2 → ψ=±1/4；
  φ=1/3 → ψ=1/6，钉住"规范可以比 cocycle 更细"）；
* `tests/subduction_constructed_stars.rs`：13346 的五个 `W` probe 全部由 catalogue 回答
  且恒等频率与 pinned 表一致；3988（二维小群）仍 `MissingChildStarData`；
* 审计微型基线由 13345/13346（均已闭合）移到 3988；恒等回归的 2,075 个 probe 现在
  **全部**走完整入口（缺数据集合为空）；settings 一组的 120 个 probe 全部完整。

### 批次 2b：二维投影表示（已交付，2026-09-22 第四轮）

**状态：已交付；普通离散标量覆盖闭合。** 全表审计（`--require-complete
--require-full-decomposition`，约 9 分钟）：

| 项目 | 批次 2a 后 | 批次 2b 后 |
|---|---:|---:|
| `full_success` | 366,039 | **366,260 / 366,260（100%）** |
| `identity_only` | 221 | **0** |
| `error` / `hard_failures` | 0 / 0 | **0 / 0** |

判词 `VERDICT complete scope=global gates=--require-complete,--require-full-decomposition
full_decomposition=complete`、**exit 0**；恒等正项 94,271/0 不匹配、Γ Frobenius
1,895/1,895、w 行 5,756/0 错误、`production_checks` 全 0 —— 即
[subduction-next-milestones.md](subduction-next-milestones.md) §R5 的验收清单全部满足
（由 R4 的三个批次达成）。

**结构（离线先钉死，脚本可复算）**：58 组的 326 个缺失星只有两个族：

| 族 | 组数 | 结构证据 | 投影不可约表示 |
|---|---:|---|---|
| 非退化 C2×C2 | 52 | 元素阶 (1,2,2,2)、交换、交换子配对 `beta_ij=phi_ij-phi_ji` 的根只有单位元 | **一个二维**，字符 `(2,0,0,0)` |
| D3 且 cocycle 是上边界 | 6 | 元素阶 (1,3,3,2,2,2)、非交换、配对恒为 0、ψ 解恰 2 个（`|Hom(D3,U(1))|`） | 规范化的普通表示 `{1,1,2}` |

第一条不是经验规则而是可证，但**证明不能**用「`g^2=e` ⇒ `omega(g,g)=±1`」——
`g^2 = e` 只对旋转部分成立，cocycle 记的是代表元乘积带出的格平移相位，所以
`omega(g,g)` 是任意非零单位根（归档的 58 组里实测出现 1/6、1/4、1/3、2/3、3/4、5/6；
child #43 的 `q=(0,1,1/2)` 上 `omega(g,g)=i`）。正确的论证：`u_g^2 = omega(g,g)·1`
且 `omega(g,g) ≠ 0`，故 `u_g` 的特征值是相异的 `±sqrt(omega(g,g))`；`u_g` 又不可能
是标量（标量会迫使 `omega(g,h)=omega(h,g)` 对所有 `h` 成立，即 `g` 与全群正交，与
配对非退化矛盾），两个特征空间都是一维，迹为 0。第二条用 `M(D3)=0`（每个 cocycle
都是上边界）：解出规范 ψ 后，三个目标 = ψ 规范 × 普通不可约表示；取哪个 ψ 解不影响
**集合**（两个解相差 sign 特征标，而 sign ⊗ 普通不可约只是置换它们）。
（复核 B 指出原论证的前提为假；`1d1b95e` 的提交信息沿用旧说法，以本段为准。）

**失败关闭**：`projective_targets` 只覆盖上面两族，且每族都有自己的结构门禁与
正交性门禁（每个字符 `(1/|P_q|)Σ|χ|² = 1`、`Σ dim² = |P_q|`）；其它结构（例如立方
Γ 点、阶 48）返回空表，生产入口继续报 `MissingChildStarData`（单测钉住）。

**证据**：

* `the_projective_tables_cover_only_the_two_gated_families`：child #43 `(0,1,1/2)` 得
  唯一二维目标且**整行**字符钉为 `(2,0,0,0)`、child #160 `(0,0,3/4)` 得 `{1,1,2}`、
  三行整表钉为 `(1,1,1,1,1,1)` / `(1,1,1,-1,-1,-1)` / `(2,-1,-1,0,0,0)` 且 `Σdim²=6`、
  child #221 Γ 为空；
* 扩展后的 `the_catalogue_reproduces_pinned_little_group_characters`：**1,328 条 pinned
  记录 / 7,578 个操作全部命中**，其中 **94 条走二维投影表**（另 1,660 条属范围外）；
  这四个数现在是 `assert_eq!` 钉值（复核 B：原先只有 `>= 100` / `>= 300` 下限），
  因为这条对照是构造目标**唯一**的外部证据；
* 端到端：ordinal 3988（C2×C2 族）12/12、12041（D3 族）19/19，均 `VERDICT clean`；
* 全表审计 366,260/366,260、`identity_only=0`、`error=0`、两个门禁 exit 0。

**边界说明**：缺口清零后 `examples/census_subduction_gaps.rs` 的输入为空；该工具现在
对空审计返回空 manifest 与零计数（不再报错），R5 的收口由完整门禁证明。也因此
"仍然缺数据"的端到端负例在真实数据里已不存在：fail-closed 由
`the_one_dimensional_solver_returns_nothing_instead_of_a_subset`（非上边界 / `MAX_ORDER`
/ `MAX_GRID` 三条都返回空表）与
`an_out_of_scope_co_group_still_reports_missing_child_star_data`（真实嵌入 + 手工构造
的越界折叠星，`build_block` 必须报 `MissingChildStarData`）两条测试保证 —— 这是覆盖
闭合后的正常状态，不应误读为门禁被删除。复核 B 之前，后一条的端到端断言在三个批次
里被逐步换成了正例，已按本节补回。

### 批次 2a 的证据计划（实现前先离线）

1. **精确性**：字符指数为有理数；逐对验证 `χ(g)χ(h) = ω(g,h)χ(gh)`、正交性
   `(1/|P|)Σ|χ|² = 1`、解集大小 = `|Hom(P_q,U(1))|`。
2. **独立对照**：在**离散** pinned k 点上（引擎已有 4,777 条普通标量行、含大量非平凡
   小余群），用同一 cocycle 求解器算出 catalogue，与 pinned 行的字符逐项比较——
   这是算法自己的证据集，不依赖待验证的缺口行；目标规模上千个 k 点。
3. **归档候选对照**：684 组的候选记录/成分数与算出的 catalogue 维数多重集比较（批次 2a
   的 626 组应一致）；不一致的组保留为缺口并单独记录，不猜。
4. **引擎门禁**：每个 probe 的维数守恒、正交性、逐操作重建；全表审计的恒等项交叉检查
   （`isotropy_subduce_*` 频率）与 `exact_target` 计数。
