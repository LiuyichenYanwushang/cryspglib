# R31 有限 Type-C Γ 可行性门禁（固定三源）

本文件记录第 31 轮**有限可行性门禁**的实际结果与边界。它不是 Type-C 生产实现，
没有修改任何生产接口、算法、依赖或旧未跟踪 WIP。固定的父群 source 是
`SG 89 GM5`、`SG 156 GM3`、`SG 191 GM5+`（均为 CDML、非 spinor、二维 Γ）。

- 基线：`cf0a4baab95ce999cb96084c9571b604d0596807`（R30 已发布版本）
- 新增文件：`examples/magnetic_type_c_gamma_probe.rs`、本文件
- 运行产物：`target/logs/r31-type-c-worker/`（不在提交范围）
- 未 commit/push；旧未跟踪 `src/irrep/magnetic_embedding/joint_congruence.rs`
  （SHA `d61d4eb7…`）原样保留、未接入、未编辑。

## 1. 结论

在**固定有限域**内找到了 **2 个实际见证**：父群二维 Γ source 的 unitary
restriction 恰为两个实 `±1`、一维、非 spinor、互不等价的 H irrep，且反幺正
coset 按 a-扭曲配对交换这两个 H 分量。两个见证都通过：

- R24 既有嵌入族的显式 `(map, shift)`：均为单位映射、零 shift；
- `verify_matrix_gamma_input` 在 `RequireExact` 下对 `η=+1` 与 `η=-1` 都成功；
- 逐类 unitary restriction `trace A(h) = χ1(h)+χ2(h)`；
- 逐类 a-扭曲配对 `χ1(a0⁻¹ h a0) = χ2(h)`、`χ2(a0⁻¹ h a0) = χ1(h)`；
- 矩阵级反幺正交换：在 H 的公共实本征基下，每个反幺正类的两个对角元为零、
  两个非对角块非零（两个 η 都检查）。

| 父群 source | record（condensing，row） | UNI | H（unitary subgroup） | Type-C corep | H 源 | map / shift | 结论 |
|---|---|---|---|---|---|---|---|
| SG 89 `GM5` | `GM2`，direction 0 | 749 | SG 16（P222） | `GM3 + GM4`，dim 2 | `GM3`（BC Γ4）、`GM4`（BC Γ2），均 dim 1 实 ±1 | `I` / `0` | witness；η=+1、−1 均通过；restriction classes=4；pairing classes=4（discriminating=2）；exchange anti=4，unitary=4 |
| SG 89 `GM5` | `GM4`，direction 0 | 751 | SG 21（C222） | `GM3 + GM4`，dim 2 | `GM3`（BC Γ4）、`GM4`（BC Γ2），均 dim 1 实 ±1 | `I` / `0` | witness；η=+1、−1 均通过；restriction classes=4；pairing classes=4（discriminating=2）；exchange anti=4，unitary=4 |

两个见证都是 SG 89 的 **Type-C 记录**，不是旧 R29 Type-A 见证
（`UNI 747`、`101`、`10`、`753`）的重复。`GM3 + GM4` 是 corep 的显示标签，
只用于报告；身份由 typed H source（`sg:ml`）和 records 决定。

## 2. 分层计数（实际运行）

固定域：3 个父群；遍历每个父群的全部标量 irrep 磁 isotropy rows；source 与
condensing 身份独立（同一 row 的 condensing 可以等于或不等于固定 source）。

| 层 | 计数 |
|---|---|
| 父群 rows | SG 89: 52；SG 156: 52；SG 191: 216；合计 **320**；与 `query::magnetic_isotropy_subgroups_of` 的 320 行一致 |
| spinor 父群记录 | **48** 条被具名跳过（生成表访问器本身拒绝 spinor isotropy 数据；它们没有本域的磁 rows） |
| distinct UNI | **110**（SG 89: 25；SG 156: 19；SG 191: 66） |
| UNI 级 Γ 读取 | summary failed **0**；Γ missing **0**；Γ unresolved **0**；**非 Γ** unresolved **0**；Γ character columns 合计 **1666** |
| Γ coreps | 合计 **777**：Type A **739**、Type B **2**、Type C **36** |
| Type-C dim 2 | 唯一 corep **36**；按 row 展开 **119** 次（同一 UNI 被多行引用时不合并记录） |
| 实 H 对筛选 | 唯一元数据候选 **2**；row 展开候选 **2**；excluded **34**；unresolved **0** |
| excluded 明细 | spinor H source **30**；复合（CIR）分量非实 **4**；没有 dimension 或 standalone 非实排除 |
| 第 4 层 | 候选 **2**，witness **2**，incomplete **0** |
| 状态 / exit | `complete_all_candidates_are_witnesses` / **0** |

SG 156 与 SG 191 在本固定域内 **0 个实 H 对候选**：SG 156 的 28 次 row 展开
Type-C dim 2 全部被具名排除；SG 191 的 86 次同样。这不是“这些群没有 Type-C”，
只是本固定筛选层的结果；没有任何后续层被声称对这些行执行过。

## 3. 方法与来源

### 3.1 rows、UNI 与 Γ

- rows 由 `query::irreps_of(sg)` 的每个非 spinor 记录的 `magnetic_subgroups()`
  枚举；另与公共行查询逐父群对账。
- 每个 distinct UNI 只读一次
  `magnetic_summary::magnetic_irrep_summary_by_uni_partial(uni)`。使用 partial
  形式是为了把 **Γ 未决**与**无关 k 点未决**分开计数，而不是因为非 Γ 失败
  就把 Γ 报成零。实际运行中两者均为 0。
- Γ 由 `kpoint.coords == (0,0,0,1)` 识别，不用显示标签代替坐标。
- 没有找到现成的、同时保留 compound 分支与 Type-C dedup 的公共“仅 Γ”summary
  入口；因此使用 partial summary 后只读取 Γ k-point。这是有意的取舍：不复制
  summary 的 dedup/compound 逻辑，也不因为全 k 计算而把非 Γ 结果混入 Γ。
- corep 的 `source_irreps` 是 H（unitary subgroup）的 typed source 身份；
  corep 显示 label 不作为身份。

### 3.2 两个实一维 H 源的核验

对每个 Type-C dim 2 corep 的每个 H source：

1. `source_irreps` 必须恰好 2 个；各 source 非 spinor、dimension 1；
2. 在 H 的 `query::find_irreps(sg, ml, CDML)` 中解析；若该 label 是个别命名的
   compound（CIR）分量而没有独立记录，则从所属 compound 的
   `compound_selected_arm_view()` 取该分量的 typed row（不是回退、不是改名）；
3. `ordinary_scalar_selected_arm_block_trace()` 的每个 character 走**有限浮点
   容差**判定（`GRID_TOLERANCE = 1e-9`）：值必须有限；imag 在容差内为 0，
   且 re 在容差内为 `+1` 或 `-1`，才进入实 `±1` 域。NaN、`±Inf`、imag NaN
   记为 `NonFinite`（**未决**，不是域外排除）；有限但非实记为 `NonReal`
   排除；有限实数但不是 `±1` 记为 `OffGrid` 排除；identity 必须为 `+1`。
   这是对生成 `f64` 值的容差分类/重构，不是对其二进制表示的精确证明；
4. 两个 source 的 typed operation/character rows 必须逐个 operation 相等或
   不等价；相等则记为 `EquivalentPair`（不可决，不当成“无候选”），不等价才
   进入候选。

实际结果：30 个 Type-C dim 2 corep 因第一个 H source 为 spinor 被
`Spinor(0)` 排除；4 个因 compound 分量 character 非实被 `NonReal(0)` 排除；
当前真实数据没有触发 `NonFinite`，也没有 `Dimension`、`MissingHSource` 或
`EquivalentPair` 未决。

### 3.3 固定嵌入

对每个过滤后的候选，使用与 R28 相同的 R24 口径：`parent_embedding_verdict`
的既有候选族 `{S·B}`，在 `realises_record_lattice` 为真的条目中按
source-independent 全序键取第一个，`RequireExact`，不搜索新的 setting。
两个见证都取到 `kind=Plain`、`map=I`、`shift=0`。没有候选需要
`CertifiedAbsence`/`NoRealising` 分支；报告仍保留这些分支。

### 3.4 source 模型与双 η

`verify_matrix_gamma_input` 用上述显式 `(map, shift)` 对固定父群 source 与
record 做输入认证，并构造 parent grey action：

- source 的 operation-aware 实矩阵表必须先通过形状、正交性、trace/typed row、
  完整有序对群律和实不可约性认证；
- 对 M 的每个类，用 R25/R27 的 class binding 绑定到 source 的存储操作；
- 反幺正类乘以 `η=±1`；
- 全 ordered-pair corepresentation law 被重新检查。

两个见证在 `η=+1` 与 `η=-1` 下都成功。没有构造 Type-C 目标，也没有调用
`subduce_matrix_gamma`（其 target 仍限 Type A）。

### 3.5 unitary restriction

对每个见证，用 summary 的 Γ character column（R27 raw first-Hall index）建立
`raw index → M class` 对照，再对每个 **unitary** class 比较：

```text
trace(A(h)) = χ1(h) + χ2(h)
```

其中 `A(h)` 是 `MatrixGammaModel` 在该 class 的二维 real action。检查同时使用
两路来源：summary 的 Γ character column 和两个 H typed row 直接求出的
`χ1(h)+χ2(h)`；两路必须一致，且都必须等于 `trace(A(h))`。source 是二维、两个
H source 都是一维且互不等价，因此 character 等式与维数一起推出
`A|_H ≅ χ1 ⊕ χ2`；不是只比较两个标签。见证的 unitary classes 各为 4。

### 3.6 a-扭曲配对与矩阵级交换

- 取一个反幺正列作为 `a0`，在 raw M 操作中求逆并计算 `σ(h)=a0⁻¹ h a0`；
  用 canonical translation 归到 M class，再映射回 summary 的列。
- 对每个 unitary column，从两个 H typed row 取 `χ1,χ2`，检查
  `χ1(σ(h)) = χ2(h)` 且 `χ2(σ(h)) = χ1(h)`。因为本域字符为实 `±1`，
  `conj` 是恒等，所以这就是 `Δ^a(h)=conjΔ(a0⁻¹ h a0)` 在该域的具体形式。
  两个见证各检查 4 个 unitary columns，其中 2 个是区分两个 H 分量的列；
  其余列上两侧字符相同，是恒等验证而非额外证据。
- 矩阵级：取一个 unitary class，其二维 action 的特征值为 `+1/-1`（trace 0、
  det −1），由它给出 H 的公共实本征基；要求每个 unitary class 在该基下对角，
  每个反幺正类在该基下两个对角元为零、两个非对角块非零。两个 η 都检查，
  见证各为 anti=4、unitary=4。

## 4. 错误、未执行与有限负结果

- **exit 0**：本次运行。所有适用层完成，2/2 候选是 witness；没有把候选数当作
  可达证明，也没有把零当作全域否证。
- **exit 1**：自检或守恒恒等式失败（程序/读数问题，不是物理负结果）。
- **exit 2**：固定域不完整：UNI summary 失败、Γ 缺失、Γ source 未决、
  real-H-pair 筛选因 H 数据不可用或 character 非有限而无法决、候选没有
  可用/可决的 R24 嵌入、或 source model 无法认证。此类运行绝不打印“0 候选”。
- 明确的排除（spinor H、有限非实 character、有限非 `±1`、dimension ≠ 1）
  与未决（缺 H row、character 不可用、`NonFinite`、equivalent pair、无
  separating class 等）在计数和输出中分开；`PairOutcome` 不把两者混成
  同一个零。

## 5. 自检

`--self-check` 使用与实际筛选相同的 `screen_pair` 函数，至少覆盖：

- 正控制：Type-C dim 2 + 两个非 spinor dim 1 实 `±1` 不等价 H 行 → candidate；
- Type-A dim 1 → `NotTypeCDim2`，不进入候选；
- 有限复 character（imag ≠ 0）→ `NonReal` 排除；
- 有限实数 0.5 → `OffGrid` 排除；
- NaN、`+Inf`、`-Inf`、imag NaN 经真实
  `real_one_dim_row_from_pairs → grid_sign → screen_pair` 数字路径 → 
  `Character{problem: NonFinite}` **未决**，不是 exclusion；
- spinor → `Spinor` 排除；
- dimension 2 → `Dimension` 排除；
- 两个相等实行 → `EquivalentPair` 未决；
- 缺 character row → `Character` 未决（不是零候选）；
- 只有一个 source 的 Type C → `SourceCount` 未决。

实际运行 `self_check violations=0`，exit 0。

## 6. 复现命令与退出

```text
CARGO_TARGET_DIR=/home/liuyichen/TB_rs/cryspglib/target \
CARGO_BUILD_JOBS=8 cargo build --offline --example magnetic_type_c_gamma_probe

RAYON_NUM_THREADS=8 ./target/debug/examples/magnetic_type_c_gamma_probe --self-check
RAYON_NUM_THREADS=8 ./target/debug/examples/magnetic_type_c_gamma_probe \
  --tsv target/logs/r31-type-c-worker/candidates.tsv
```

最终实际命令与计数：

- build：exit 0；
- self-check：exit 0，`self_check violations=0`；
- 全门禁：exit 0，`status complete_all_candidates_are_witnesses`；
- `cargo clippy --offline --example magnetic_type_c_gamma_probe -- -D warnings`：
  exit 0（只有 workspace 级 Cargo profile 前导 warning，与 example 无关）；
- `rustfmt --edition 2024 --check examples/magnetic_type_c_gamma_probe.rs`：
  exit 0；
- `RAYON_NUM_THREADS=1` 重跑：exit 0；`candidates-1.tsv` 与 8 线程的
  `candidates.tsv` 逐字节相同，去掉 `tsv` 路径行后的 stdout 也相同；
- 调试构建单次扫描约 1 秒量级（8 线程实测 WALL=1.12 s）；TSV 行数 2
  （不含表头）。

`--tsv` 的稳定字段（表头）：

```text
parent_sg source_sg source condensing direction_index uni unitary_sg corep_label
h1_sg h1 h1_bc h2_sg h2 h2_bc h1_signs h2_signs corep_type corep_dim
stage reason embedding_kind map shift eta_plus eta_minus restriction pairing
exchange_plus exchange_minus
```

`map` 为 9 个 `num/den` 行主序值（`;`、`,` 分隔），`shift` 为 3 个
`num/den` 值；`stage` 为 `witness` / `not_established` / `unexecuted`。

## 7. 未验证 / 明确不在范围

- 没有实现生产 Type-C 目标、完整 target corep 矩阵、intertwiner 或全分导；
  `subduction`/`magnetic_subduction` 的 Type-A-only 状态未改。
- 只检查 R24 既有嵌入族 `{S·B}` 中 **第一个** exact-realising 条目；不是所有
  setting、不是所有 map，也没有 setting 搜索。R24 的 `certified_absence`
  只表示该候选族内无成员，不是“记录不属于任何子群”。
- SG 156、SG 191 的零是**本固定 real-H-pair 筛选层**的有限结果，不是全域
  Type-C 不可能性；没有对它们继续构造或求解未出现的候选。
- 30 个 spinor-H Type-C dim 2 和 4 个复合非实分量按声明的实 `±1` 域排除；
  这不是说 spinor/复 Type C 不存在。
- 48 个 spinor 父群记录因当前生成表访问器拒绝 spinor isotropy 数据而跳过；
  没有声称它们没有磁子群，只报告本 API 边界。
- Type-C 分类本身取自现有 magnetic summary 的 Wigner 结果；本轮没有独立重算
  Wigner `W` 判据，也没有外部 corep oracle。门禁验证的是**源 restriction 与
  a-扭曲配对**，不是重新证明 Type-C 分类。
- H 的实 `±1` 字符来自 generated typed row，并用 `1e-9` 容差做**有限浮点**
  分类/重构；这不是对 f64 字符值的精确证明，也没有外部字符表对照。
  source 矩阵认证沿用 R29 的 Γ-equivalent 认证口径，不是对原始 PIR
  provenance 的重新证明。
- 没有 non-Γ source、full-star、general B/C、高维、complex scalar、setting
  solver、band assignment 或外部 oracle 对照。
- 旧 `canonical_translation`/`in_lattice` 的极端 i128 输入边界未在本 probe
  中系统测试；本运行只走到实际语料。
- 两个 witness 的 H source 字符在 typed row 上为实 `±1`，但未对每个 operation
  打印完整字符表；打印的是维数与 `+`/`-` 计数，完整逐操作配对由
  `pairing=ok classes=4 discriminating=2` 和矩阵交换计数支持。

## 8. 本次改动边界

- 只新增 example 与本 doc；未编辑现有生产源、测试、依赖、Cargo 配置或旧 WIP。
- 未 commit/push。
- 运行日志：`target/logs/r31-type-c-worker/`。其中早期 `run*.out` 与
  `metadata*` 是开发迭代记录；最终证据是 `final.out`、`candidates.tsv`、
  `self-check.out`、`clippy.out`。
- 最终文件哈希（SHA-256）：
  - `examples/magnetic_type_c_gamma_probe.rs`
    `132b7e24f8fa98e8d5d6103f0b1342e8cdbaa734b31a0e5c47ee483ab2763b1a`
  - `final.out`
    `4cbbb33a1e6c23dca0153a730a73b1adac52a30e9e8771f602adc61a7e913fd1`
  - `candidates.tsv`
    `44d838e2745efedc38ec5ce27d669b21272be436614a3e91b8701e6b41b4af89`
  - `self-check.out`
    `bc45bf9327dfa8c5a82560cb504ccb15841be9ab54d2a5f75f4a31ff31426a30`
  - `clippy.out`
    `5761ca6a98e389bab2581ac3181a6eda6b0e2060c02c110067ae8fb17552abb4`
- 工作树：只有两个新增未跟踪文件与旧未跟踪 WIP
  `src/irrep/magnetic_embedding/joint_congruence.rs`；没有 tracked 文件改动。
- 开发过程中一次 `cargo fmt -- <file>` 误把整个 crate 的 tracked 源文件格式化了；
  已用 `git restore --worktree -- .` 还原，并在还原后的基线上重新构建、运行、
  自检和 clippy。以上最终数字均来自还原后的基线。
