# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Workspace context

This crate is a **workspace member** inside `/home/liuyichen/TB_rs`. All cargo commands must be run from the workspace root:

```bash
cd /home/liuyichen/TB_rs
cargo build --package cryspglib
cargo test  --package cryspglib
cargo check --package cryspglib
```

---

## 工作流铁律

### 规则 1: 每完成一个可编译的修改就立即 commit

**每次 `cargo check` 成功后必须立即 commit**，然后再做下一个修改。不要连续做多个修改才 commit。

```bash
git add -A && git commit -m "描述"
```

**Why:** `git checkout` 恢复时只保留已提交的内容。中间修改全部丢失。宁可 commit 太多（事后 squash），不能丢失工作。

**反面案例（发生过两次）：**
1. 2163→3307 行的诊断代码因 git checkout 全部丢失
2. MagneticOps→SymmetryOps 重构中，已修复的 10+ 处 field access 因一次 revert 全部丢失

### 规则 2: 不用 Python 脚本做代码修改

批量 sed/Python 替换容易产生意外后果。代码修改必须用 Edit/Write 工具逐处进行，每处修改后确认正确。

### 规则 3: 不用 type alias 做过渡

`pub type OldName = NewName;` 只是把问题藏起来，缺少可维护性。应该全局替换所有引用，然后删除旧定义。

### 规则 4: 先 use 再去掉 crate:: 前缀

替换类型时，先在文件头部添加 `use crate::NewType;`，再把文件中所有 `crate::NewType` 替换为 `NewType`。

---

## ISOTROPY 数据格式知识

### PIR vs CIR 的本质区别

- **CIR**（Complex Irreducible Representation）：只依赖**旋转矩阵 R**，不依赖 translation。
  每个独特的旋转类型有一个 complex character `χ(R) = (re, im)`。
  opcount = little co-group 的大小（distinct rotation types）。

- **PIR**（Physically Irreducible Representation）：依赖完整的空间群操作 `{R|t}`。
  字符通过 CIR + Bloch 相位组合：`PIR({R|t}) = Σ_i CIR_i(R) * exp(i*2π*k·t)`
  opcount = full little group 的大小（包含所有 translation 变体）。

- ISOTROPY 的 PIR 和 CIR 操作数不同是结构性的，不是数据缺失。
  CIR 永远只有 distinct rotation types 的条目。

### CIR_data.txt 格式

```
seq sg "name" "label" dim kx ky kz kd
<16-int complex matrix line>  ← 每个 operation 一行
<16-int complex matrix line>  ← 第二行（和第一行一起编码复矩阵）
(re, im)                      ← complex character
...                            ← 重复（如果有更多 operation）
```

- 第一行 16 个整数是复矩阵的编码（ISOTROPY 特有格式）
- 第二行 16 个整数也是复矩阵的一部分
- 每个 operation 由 2 行矩阵 + 1 行字符组成
- CIR 无 irtranslation（字符不依赖 translation）

### spglib Hall vs ISOTROPY 的 translation 差异

- ISOTROPY 使用 **primitive cell** 的 translation
- spglib Hall 使用 **conventional cell** 的 translation（可能含 centering）
- 两者 translation 不同，但**旋转矩阵相同**
- 重排时：PIR 字符是实数，只需排列不需相位修正（ISOTROPY PIR 字符对 matched operation 是正确的）
- CIR 展开时：额外 Hall 位置需从 matching rotation 复制 + Bloch 相位 `exp(i*2π*k·t_hall/kd)`

### 重排后数据一致性

重排后 PIR 和 CIR 在同一 Hall 位置的字符来自**同一个 ISOTROPY operation**
（因为 mapping[h] 唯一定义了 ISOTROPY 索引）。因此 PIR = CIR_sum 在排列后
仍然成立。CIR 展开只影响额外位置（CIR 源数据没有的 translation 变体）。

---

## 方法论：数据优先于算法（CIR-PIR debug 教训）

### 问题回顾

CIR-PIR 测试失败（596 个 compound irrep, 1222 mismatches）。花了大量时间尝试
Bloch 相位修正（`exp(i*2π*k·Δt)`），每改一次 mismatches 反而增多。
用户几个简单问题（"PIR 和 CIR 分别是什么数据？"、"两者的 translation 可以对比吗？"）
引导重新检查 ISOTROPY 源数据格式，才发现 CIR 字符只依赖旋转不依赖 translation。

### 为什么会陷入思维玄幻

| 阶段 | 错误做法 | 应该做的 |
|------|---------|---------|
| 1 | 假设 CIR"应该"和 PIR 操作数相同 | 先检查两者的源数据格式 |
| 2 | 尝试 Bloch 相位公式修正（多次迭代失败） | 先验证假设：CIR 真的需要相位吗？ |
| 3 | 修改越来越复杂（Hall trans, ISO trans, Δt...） | 回退到数据源，理解 CIR 格式 |
| 4 | 每轮 fix 让 mismatches 增加而不是减少 | 如果修不好，说明假设错误，不是公式不对 |

### 核心教训

**规则：当修复让问题变严重时，停止修代码，回到数据源检查假设。**

1. **数据 > 算法**：不知道数据长什么样之前，不要写任何修正公式
2. **简单的问题解决复杂的问题**：用户问"数据分别是什么"直接打破了思维定式
3. **物理/数据约束是第一位的**：CIR 没有 translation 不是 bug，是结构性的——
   CIR 是小群表示（只关心旋转），PIR 是空间群表示（含平移相位）。
4. **失败的修复是信号**：如果改了 3 次 mismatches 还在增加，假设一定错了

### 量化对比

| 指标 | Bloch 相位路线 | 最终方案 |
|------|---------------|---------|
| 修改次数 | 8+ commits, 每次 ~100 行 | 1 commit, ~30 行 |
| mismatches 趋势 | 0→1276→1222（恶化） | 0（直接解决） |
| 修改范围 | Python 脚本 + Rust 测试 | Rust 测试 min(n_ops, cir_ops) |
| 根本认知 | CIR"缺失"数据需要推算 | CIR 本来就只覆盖旋转类型 |

---

## 调试方法论 — 从 spinor Wigner 排查中提炼的经验

### 原则 1：比较 passing vs failing cases，找差异因子

当同一个代码路径**有些 case 通过、有些失败**时，不要猜测通用原因。直接比较一个成功 case 和一个失败 case，问：

> **"这两个 case 之间什么不同，导致一个通过、一个失败？"**

这个差异通常直接揭示根因。

历史实例：SG2 T-point passed after loop fix → SG159 L-point still failed。当时观察到
SG2 的 LG {I, -I} 包含 `(a₀h)²`，而 SG159 的 LG {I, mirror} 不包含。
这个比较成功定位了“失败发生在 square/LG mapping 阶段”，但当时进一步断言
“不是 setting/algorithm 问题”是错误的；2026-06-19 的 UNI663 oracle 已证明，
MSG 嵌入基底和 standalone H Hall/spin-table 基底不一致本身就会制造这种
`square_not_in_spin_table` / `square_outside_little_group` 现象。

### 原则 2：假设驱动的逐层排除

不要在无数可能原因中随机尝试。为每种假设设计一个**最小 oracle test**，一票否决或确认：

1. **列出所有可能假设**（按先验概率排序）
2. **为每个假设设计一个 oracle**（最小代码改动，只输出统计数据）
3. **跑 oracle，看数据**——如果数据否决假设，立刻排除，不再纠结
4. **如果 oracle 确认假设，再设计修复**

实例——NONE=1,007 的排查顺序：
- H1: same-rotation lift 误选 → scan same-rot candidates → OTHER=0 → 否决
- H2: UU* antiunitary square → 6 formulas oracle → 6.5% fix → 不是主因
- H3: H/G gauge mismatch → 当时的 G-gauge oracle → 0% fix → **错误地否决**
- H4: det=-1 improper → det stats → 混合分布 → 否决
- H5: J-insertion → J-oracle on NONE → 61% fix → 确认方向
- H5-global: global J → 88.3%→83.0% → 不能全局替换 → 否决
- H5-per-case: case-level J fallback → old_fail_j_ok=22/945 → 否决

每排除一个假设，就缩小搜索范围。不要跳过 oracle 直接修代码。

### 原则 3：诊断与修复分离

诊断代码（oracle/counter/scan）**不应改变正式分类结果**。先加诊断、跑数据、看统计、确认假设，再设计修复。

- Oracle 只在 `None`/失败分支执行，不改变 `return` 值
- 计数器用 `AtomicUsize`，在 diagnostic test 中读取
- 正式路径保持原样，等 oracle 确认后再改

### 原则 4：per-term → per-case 的层级

在 Wigner sum 中，单个 term 的修复不等于整个 case 的修复：

1. **per-term fix**（只对失败 term 用新公式）：数学上危险，同一个 sum 混用两个 convention
2. **global fix**（所有 term 都用新公式）：如果破坏了更多正常 term → regression
3. **per-case fallback**（先试旧公式，整个 case 失败再试新公式）：唯一理论上干净的 fallback

要区分三者，不能看到 per-term oracle 有 61% fix 就急于做 per-term patch。

### 原则 5：语义正确的计数器命名

计数器名字必须准确反映被计数的**物理/数学含义**，不能有歧义：

- 错误：`central=false` → "raw misses"。正确：`central=false` → "same lift, no central element"
- 错误：`theta2_fixes` → "fixing misses"。正确：关系到 `±u_k` 的 same/Ebar/none 三类
- 错误：`NONE=0` → "0 mismatch"。正确：`NONE=1,007` → "1,007 non-trivial mismatch"

错误命名会误导后续分析方向。本例中 `central=false` 被误读为 "raw failure"，导致构造了大量无用的 sign-flip 修复。

### 原则 6：不要过早下结论说"需要大工程"

Data generation 存 `central_parity` 或 `extended character table` 可能是最终方案，但在确认以下问题之前不应断定：

1. **先确认问题确实来自数据缺失**（而非 algorithm bug or convention mismatch）
2. **先做 oracle 估计大工程的收益**（例如 eta ±1 测试）
3. **先排除更便宜的修复**（runtime inference, convention alignment）

### 原则 7：逐阶段确认，不要把 UNI=0 当成整个管线失败

磁群识别是一条多阶段管线。当 `UNI=0` 时，不要笼统地说"磁群识别失败"。

分别检查每个阶段的输出：
- 普通 SG 是否正确？Hall 是否正确？
- 磁类型是否正确？操作数是否正确？
- unitary/anti-unitary 比例是否正确？
- 失败是否**仅**发生在 DB matching 阶段？

实例：石墨烯 AFM 返回 `SG=191, Hall=485, Type-3, 24 ops (12U+12A)`
但 `UNI=0`。前四个阶段全部正确——问题只在 DB 匹配阶段。
这排除了晶格约定、磁操作生成、FSG/XSG 分类等所有问题。

### 原则 8：容差扫描是分类工具，不是修复方法

对 `symprec` 做跨数量级扫描（1e-3 → 1e-6）：

- 结果随容差变化 → 优先检查数值稳定性
- **结果跨多个数量级完全不变** → 优先检查坐标约定、数据 setting、变换方向、群操作合成

石墨烯 AFM 的 UNI=0 在四个数量级容差下完全不变——这排除了数值问题，
直接指向 convention/transform 错误。

### 原则 9：对矩阵方向使用推导，不使用记忆

看到 `T R T⁻¹` 或 `T⁻¹ R T` 时，不要凭变量名判断。先写清楚：

```
x_new 与 x_old 的关系是什么?
```

然后计算 `C g C⁻¹`。推导结果是最小、最可靠的 oracle。

对于 `x_std = T x + s`，正确的 Seitz 共轭是：
```
R_std = T R T⁻¹
t_std = s - R_std s + T t
```

### 原则 10：检查 helper 的所有调用者

一个错误 helper 可能在某个调用点被另一个错误抵消，却在其他调用点失败。

石墨烯案例：`get_distinct_changed_magnetic_symmetry` 内部使用 `T⁻¹ R T`（错误），
但 `get_reference_space_group` 传入的 `tmat` 本身也取反了（把 P 当 T，实际应为 P⁻¹）。
两个错误在参考 setting 路径上互相抵消，但在数据库 correction transformation 路径
（这里 T 已经正确）暴露出来。

单独分析第一个调用点会错误地认为公式"等效"。

### 原则 11：高对称测试不足以验证线性代数约定

立方晶格的变换矩阵通常是单位矩阵、置换矩阵、正交对称矩阵或自逆矩阵
（`T = T⁻¹`）。这些测试全部通过**不能**证明矩阵方向正确。

应当至少包含：
- 非正交六方/单斜晶格（`T != T⁻¹`）
- 非零 origin shift
- 非平凡 correction transformation
- 具体断言 UNI/BNS 值（而非只检查 `> 0`）

石墨烯六方 AFM 是第一个暴露这些问题的非正交 oracle。

---

## 错题集 — 核心教训

### Bug 1: loop domain ≠ character domain
Wigner 求和对象是 **little co-group**（旋转），不是 full little group（Seitz 变体）。loop 遍历了 4 个 Seitz 但 spin table 只有 2 个 co-group 条目 → W=-0.5。

### Bug 2: rotation matching 不能跨基底
rotation 只在同一基底下不变。`x_new = T x_old + s` 时 `R_new = T R T⁻¹`。

### Bug 3: 不依赖数组顺序隐含语义
Grey 群的 a₀ 必须是纯 θ (R=I)。取 `antiunitary[0]` 可能取到 θ·g。

### Bug 6: Θ²=Ē 和 SU(2) central sign
`central = !spatial_central` 的推导：Θ²=Ē=-I（自旋 1/2），spatial=EBAR 时 Ē²=I→SAME→central=false。
**但仍需验证**：`su2_same_up_to_sign(U_b², U_{h²})` 跨 G/H 两套 SU(2) gauge 比较可能出错。

### 已排除的假设
spin 数据库不完整 ❌ | Pauli SU(2) 合成 ❌ | same-rotation lift 误选 ❌ | 
UU* 公式 ❌ | det=-1 improper ❌ | global J-insertion ❌（regression）

### 未排除
H/G gauge mismatch — 跨 gauge SU(2) 比较仍需独立验证。

## Architecture overview

cryspglib has two major subsystems:

**1. spglib port** — space group identification from crystal structures.
`Crystal::new(lat, positions, types)` → `.analyze()` → `.dataset()` → `SpaceGroup`.
Also supports magnetic space group identification (1,651 UNI types) via `.with_magnetic()`.

**2. irrep module** — irreducible representation data for all 230 space groups.
`irreps_of(sg_number)` → `IrrepRecord` (labels, characters, matrices, isotropy subgroups, magnetic corepresentations).
100% data coverage: 4,777 irreps with character tables (~50k values), full matrices (~580k values),
15,239 non-magnetic + 16,721 magnetic isotropy subgroups.
100% of characters are in spglib Hall order.

---

## Build & Test Commands

```bash
cd /home/liuyichen/TB_rs

cargo build --package cryspglib
cargo test  --package cryspglib
cargo test  --package cryspglib <test_name>
cargo check --package cryspglib
```

Key diagnostic test: `irrep::corep::tests::diagnose_wigner_sources -- --nocapture`

Enable verbose Wigner diagnostics (eprintln! output):

```bash
cargo test --package cryspglib --features debug-corep -- <test_name> -- --nocapture
```

### Data regeneration pipeline

`hall_operations.json` is a committed static artifact (does not need regeneration).

```bash
# Regenerate generated_data.rs from ISOTROPY source data:
python3 scripts/generate_irrep_data.py
# Full pipeline shell:
bash scripts/regenerate_all.sh
```

After regeneration, run the full test suite to validate:

### Diagnostic validation scripts

| Script | Purpose |
|--------|---------|
| `validate_cir_pir.py` | Standalone CIR→PIR validation — checks PIR = Σ CIR * exp(i2πk·t) per Hall op |
| `check_iso_vs_spglib.py` | Compare ISOTROPY primitive vs spglib Hall conventional translations |
| `check_phase_correction.py` | Analyze Bloch phase corrections when mapping ISOTROPY→spglib Hall order |
| `debug_cir_pir_sg9.py` | Single-case debug: SG9 CIR/PIR mapping details |
| `test_su2_closure.py` | Pauli SU(2) composition closure test |
| `test_spinor_wigner_formula.py` | Spinor Wigner formula standalone test |

---

## Key types

| Type | Location | Description |
|------|----------|-------------|
| `Crystal` | `api.rs` | Entry point: lattice + positions + types + optional magnetic moments |
| `SymmetryAnalysis` | `api.rs` | Builder for symmetry analysis (`.symprec()`, `.dataset()`, `.magnetic_dataset()`) |
| `SymmetryOps` | `api.rs` | Ordered set of `{R\|t}` + time_reversal, with `from_database(hall_number)` |
| `SymmetryOp` | `api.rs` | Single `{R\|t}` with rotation, translation, time_reversal |
| `SpaceGroup` | `lib.rs` | SG number, Hall number, ops, Wyckoff positions, standard cell |
| `MagneticDataset` | `lib.rs` | MSG result: UNI number, type, rotations, translations, time_reversals |
| `MagneticSymmetry` | `lib.rs` | MSG + symmetry ops combined (implements `Display`) |
| `MagneticSpaceGroupType` | `lib.rs` | MSG type lookup: `.from_uni()`, `.classify()` |
| `SpaceGroupType` | `lib.rs` | SG type lookup: `.from_hall()` |
| `IrrepRecord` | `irrep/types.rs` | Irrep: labels, dim, k-vector, characters, matrices, subgroups, corepresentations |
| `SpinLiftContext` | `irrep/wigner.rs` | H and G spin ops for Wigner test |
| `SeitzOp` | `irrep/wigner.rs` | `{R\|t}` with optional time reversal |
| `CorepType` | `irrep/corep.rs` | A/B/C/Unsupported |

---

## Module structure

### spglib port subsystem

| Module | Role |
|--------|------|
| `api.rs` | `Crystal` (entry point), `SymmetryAnalysis` (builder), `SymmetryOps`, `SymmetryOp` |
| `lib.rs` | `SpaceGroup`, `SpaceGroupType`, `MagneticDataset`, `MagneticSymmetry`, `MagneticSpaceGroupType`, `SymError` |
| `cell.rs` | `Cell` (lattice + positions + types + optional tensors) |
| `symmetry.rs` | `Symmetry` (raw symmetry operations, N×rot+trans arrays) |
| `spacegroup.rs` | `Spacegroup`, `spa_search_spacegroup*` |
| `spg_database.rs` | `spgdb_get_spacegroup_operations`, `spgdb_get_spacegroup_type` |
| `magnetic_spacegroup.rs` | MSG identification: `msg_identify_magnetic_space_group_type` |
| `msg_database.rs` | `msgdb_get_magnetic_spacegroup_type` (1,651 UNI entries) |
| `msg_database_gen.rs` | Auto-generated MSG database (shipped as source) |
| `spin.rs` | Spin-polarized symmetry: `spn_get_operations_with_site_tensors` |
| `pointgroup.rs` | `ptg_get_pointgroup`, `ptg_get_transformation_matrix` |
| `primitive.rs` | `prm_get_primitive`, `prm_get_primitive_symmetry` |
| `delaunay.rs` | Delaunay lattice reduction |
| `niggli.rs` | Niggli lattice reduction |
| `hall_symbol.rs` | Hall symbol parsing/conversion |
| `kpoint.rs` | k-point grid generation, irreducible mesh |
| `kgrid.rs` | Grid address utilities |
| `determination.rs` | Space group determination pipeline |
| `refinement.rs` | Cell refinement / idealization |
| `overlap.rs` | Atom overlap detection |
| `parser.rs` | POSCAR parser |
| `site_symmetry.rs`, `sitesym_database.rs` | Site symmetry + database |
| `arithmetic.rs` | Arithmetic crystal class symbols |
| `mathfunc.rs` | `Mat3`, `Mat3I`, `Vec3`, matrix/vector operations |
| `debug.rs` | Diagnostic print helpers |

### irrep subsystem

| Module | Role |
|--------|------|
| `irrep/mod.rs` | Module docs, re-exports, coverage summary (4,777 irreps, 100% coverage) |
| `irrep/types.rs` | `IrrepRecord`, `IsotropyRecord`, `MagneticIsotropyRecord` |
| `irrep/query.rs` | `irreps_of()`, `kpoints_of()`, `format_character_table()` |
| `irrep/corep.rs` | Co-representation: `compute_coreps()`, `CorepType`, diagnostic tests (30+ tests) |
| `irrep/wigner.rs` | Wigner test: Seitz composition, SU(2) composition, spinor classification |
| `irrep/bridge.rs` | `impl SpaceGroup` — bridge APIs linking spglib port → irrep |
| `irrep/generated_data.rs` | Auto-generated static arrays (~753k lines). `include!()`-d into `types.rs` |
| `irrep/settings_data.rs` | Hall→setting mappings. `include!()`-d into `generated_data.rs` |
| `irrep/wigner_extra.rs` | Pre-computed antiunitary character path. `include!()`-d into `wigner.rs` |
| `irrep/preamble.rs` | Generated data prelude |
| `irrep/{triclinic,monoclinic,orthorhombic,tetragonal,trigonal,hexagonal,cubic}.rs` | Per-crystal-system irrep data (`include!()`-d into `generated_data.rs`) |

---

## Test suite (189 tests across 7 binaries)

All tests pass (as of 2026-06-18).

| Binary / Location | Tests | Description |
|-------|-------|-------------|
| `src/irrep/corep.rs` (unit) | 25 | Wigner diagnostics, BCS validation, CIR-PIR cross-validation, self-consistency invariants |
| `src/irrep/{query,types,corep,mod}.rs` (doctests) | 5 | API examples |
| `src/{arithmetic,cell,debug,delaunay,determination}.rs` (unit) | 12 | Arithmetic crystals, overlap detection, lattice reduction, error handling |
| `src/{api,lib}.rs` (doctests & unit) | 21 | Entry-point API examples and type constructors |
| `tests/irrep_validation.rs` (integration) | 31 | Full-sweep validation: every SG has irreps, dimensions match, labels well-formed, k-vectors positive |
| `tests/magnetic_integration.rs` (integration) | 11 | Magnetic structure analysis end-to-end (graphene, bilayer, Fe, CoF3, etc.) |
| `tests/{cof3,crps4,la2nio4,bcs_corep_validation}.rs` (integration) | 5 | Reference material cases |
| **Total** | **~189** | |

Key diagnostic tests (most useful for Wigner debugging):

```bash
cargo test --package cryspglib diagnose_wigner_sources -- --nocapture
cargo test --package cryspglib diagnose_spinor_wigner_per_term -- --nocapture
cargo test --package cryspglib diagnose_none_examples -- --nocapture
cargo test --package cryspglib test_cir_pir_cross_validation -- --nocapture
```

Full validation sweep (integration tests, ~1 min):

```bash
cargo test --package cryspglib --tests
```

---

## 磁空间群识别的坐标变换约定

完整故障复盘、数学推导、错误互相掩盖机制和调试方法见：
`docs/magnetic-spacegroup-basis-transform-postmortem.md`。

### 标准设置变换

`magnetic_spacegroup.rs` 使用

```text
x_std = (T, s) x
```

因此 Seitz 操作必须按下面的方向变换：

```text
R_std = T R T⁻¹
t_std = s - R_std s + T t
```

`get_reference_space_group` 返回的 `tmat` 是参考空间群
`bravais_lattice` 的逆矩阵。不要把共轭方向改成 `T⁻¹ R T`，
也不要直接把 `bravais_lattice` 当作 `tmat`；立方体系可能碰巧通过，
但六方等非对称基底会无法匹配正确的 UNI。

**回归案例——石墨烯 AFM（z 方向反铁磁）**：

```text
SG=191 (P6/mmm), Hall=485, type=BlackWhite
24 ops (12 unitary + 12 anti-unitary)
UNI=1466, BNS=191.236
```

测试命令：

```bash
cargo test --package cryspglib --test magnetic_integration test_graphene_afm_z
```

该测试覆盖 `symprec=1e-3..1e-6`，结果必须保持一致。此前的 `UNI=0`
不是二维体系的数据库限制，而是标准设置变换方向移植错误。

**双层石墨烯 oracle（z=0.51 / z=0.49）**：

打破水平镜面对称后，空间群从 P6/mmm (#191) 降到 P-3m1 (#164)：

| 磁构型 | SG | Hall | ops | UNI | BNS |
|--------|-----|------|-----|-----|-----|
| 非磁 | 164 | 456 | 12 | 0 | - |
| FM (都↑) | 164 | 456 | 12 | 1319 | 164.89 |
| AFM (一↑一↓) | 164 | 456 | 12 | 1318 | 164.88 |

这是第二个非立方磁群 oracle——对称操作数从 24 降到 12，
且涉及三方晶系（介于立方和六方之间），对 setting transformation 的敏感度
高于立方案例但低于平面六方案例。

---

## 错误处理约定：Result, not Option

spglib port 的主要公共 API 已全部从 `Option<T>` 迁移到 `Result<T, SymError>`。
`SymError` 是 unit-only enum（无数据字段），每个 variant 名直接指向失败位置。

### 为什么不用带数据的 Error

`SymError` 的 variant 已经精确到单个函数甚至函数内的单个失败点。
携带额外数据不会增加定位精度，反而破坏 `Copy` 和简单模式匹配。

### 改造的函数

| 管线 | 函数 | 错误 variant |
|------|------|-------------|
| 磁群 | `spn_get_operations_with_site_tensors` | `MagneticOpGenerationFailed`, `MagneticPrimitiveLatticeFailed` |
| 磁群 | `msg_identify_with_parent_hall` | `MagneticReferenceGroupFailed`, `MagneticFallbackReferenceFailed`, `MagneticUniMatchFailed` |
| 磁群 | `magnetic_dataset()` | 传播上游错误 |
| 非磁 | `prm_get_primitive` | `CellStandardizationFailed` |
| 非磁 | `sym_get_operation` | `SymmetryOperationSearchFailed` |
| 非磁 | `spa_search_spacegroup` | `SpacegroupSearchFailed` |
| 非磁 | `det_determine_all` | `SpacegroupSearchFailed` |
| API | `Crystal::from_poscar` | `InvalidInput` |
| API | `SymmetryOps::from_magnetic_database` | `SpacegroupSearchFailed` |
| API | `SymmetryOps::from_sg` | `SpacegroupSearchFailed` |
| API | `find_hall_number` / `find_first_hall_for_uni` | `SpacegroupSearchFailed` |

### 错误传播模式

- 使用 `?` 直接传播同类型错误
- `Option` → `Result` 转换用 `.ok_or(SymError::Variant)?`
- `Result` → `Option` 转换（仅在临时 backward compat 中）用 `.ok()`
- `MagneticUniMatchFailed` 在 `lib.rs` 中被特殊处理：不传播，而是触发
  FSG/XSG fallback（因为 UNI 匹配失败不是 fatal——仍可从群阶数推算磁类型）

### 晶格矩阵约定

`lattice[cart][vec]` —— **行=笛卡尔分量，列=晶格矢量**。
详见 `mathfunc.rs` 模块文档。六方晶格不对称——约定错误将导致
空间群识别错误（如 graphene #191 → #10）。

---


## 审计教训：Phase 1 假修复复盘

> 以下是从"136/136 fixed 实为 0/136"事件中提炼的核心教训。完整记录见 git log。

### 六类错误

| # | 错误 | 预防 |
|---|------|------|
| 1 | **失败阶段转移冒充修复**：`!matches!(result, SquareNotInSpinTable)` 把症状转移当治愈 | 判据必须是**总失败数下降** |
| 2 | **总失败数不变不追问**：679 自始至终没变，却自我合理化 | 总数不变时先追问为什么 |
| 3 | **声称完成但只做了 1/5**：计划要求 origin/Seitz双射/多解检测，实际只有 rotation multiset | 逐条对照计划清单打勾 |
| 4 | **挑选有利指标**：只看 W_FAIL 和 OLD_PATH_FAIL 下降 | 总失败数是唯一不可作弊的指标 |
| 5 | **用 sed 批量改代码**：违反规则 2 | 只用 Edit 逐处修改 |
| 6 | **不区分诊断/生产路径**：诊断传了 setting_xf，生产传 None，却声称"已接入" | 两者必须同时验证 |

### 核心原则

- **判据必须是总失败数下降，不能是中间阶段转移**
- **宣布完成前逐条对照计划清单**
- **总数不变时必须追问为什么，不能自我合理化**
- **只做诊断 oracles，不急于改正式分类结果**

---

## 排查方法论精华

### 核心流程

```text
诊断 oracle → 全量统计 → 确认 convention → 找到边界 case → 修正 → 验证 → 下一个
```

### 关键实例

**k-convention 排查**：
1. oracle: 对 3611 个 spinor irrep 比较 Rk vs R⁻ᵀk
2. 结果: reciprocal_exact=3611 → 确认 R⁻ᵀ 是正确 convention
3. 边界: centered cell 还需 pure translation phase check
4. 数据 bug: 1/3,1/6 k 点被错误 rationalize → 修复 parse_spinor_data.py
5. 效果: 679→347 (-332)

**Setting transform 排查**：
1. oracle: UNI663 比较 ops_from_msg vs ops_from_hall rotation → C2z≠C2y
2. 假修复: 阶段转移冒充修复 (136/136 实为 0/136)
3. 真修复: 完整 (T,s) 求解 (Gaussian elimination + modulo-1)
4. 发现: 48 signed-perm 不够 → 泛化为 rational Mat3
5. 入口 bug: UNI187 unitary 含 mirror_y 但被识别为 SG1 → Hall 选择错误

### 铁律

1. **数据 > 算法**：先检查源数据再调公式
2. **判据 = 总失败数下降**：不接受阶段转移或其他中间指标
3. **全量 oracle，不靠单例**：`diagnose_spin_lg_k_convention` 遍历全部 3611 个 irrep
4. **100% 通过后仍追问边界**：reciprocal_exact=3611 后仍发现 centered cell 问题
5. **入口数据出错时停止下游修复**：UNI187 的 unitary 操作本身错了
6. **穷举 convention，让数据说话**：不确定 Rk vs R⁻ᵀk 时两个都算，比较 exact match

---

## Spinor Wigner 交接记录（2026-06-19，DeepSeek 从这里继续）

> 本节覆盖前面“679 是当前基线”“Phase 1 完成”等旧结论。不要沿用旧分母或旧失败分类。
> 当前代码停在提交 `d88d281`，工作树干净。

### 已确认并修复的两个根因

#### 1. spin.dat 三分之一/六分之一 k 点被错误生成成 Gamma

`scripts/parse_spinor_data.py` 原先用 `1e-6` 判断十进制坐标能否化成
`1/3`、`1/6`。源文件使用 `0.333333`，乘 3 后的误差可能略大于
`1e-6`，导致 H/K 等 k 点最终回退为 `(0,0,0)/1`。

已完成：

- 新增 `_rationalize_kvector()`，容差按源数据精度改为 `5e-6`；
- 无法落入支持分母 `{1,2,3,4,6}` 时直接报错，不再静默回退 Gamma；
- 重新生成 `src/irrep/generated_data.rs`；
- 234 个 spinor irrep 的错误 k 数据被纠正。

相关提交：

- `0186394 fix spinor third-coordinate rationalization`
- `05a7a4c regenerate corrected spinor k vectors`

#### 2. little-group 判定使用了错误的倒空间作用

原实现检查 `R k`。分数坐标下的倒空间作用应为 `R^{-T} k`。
中心化常规胞还不能只检查分量是否为整数：得到的 reciprocal shift
必须对所有 unitary pure translations 给出整数相位。

正确条件已经接入 `filter_little_group`：

```text
unitary:     R^{-T} k - k ∈ L*
antiunitary: -R^{-T} k - k ∈ L*
```

验证 oracle 覆盖全部 3611 个 spinor irreps：

```text
reciprocal_centered_exact = 3611
reciprocal_centered_fp/fn = 0/0
```

相关提交：

- `60dfb99 fix reciprocal little-group filtering`
- `982fcaf support setting-aware little-group filtering`

### 当前有效基线

最新成功运行的命令：

```bash
cargo test --package cryspglib diagnose_wigner_sources -- --nocapture
```

结果：

```text
spinor_complex_ok   = 20,463
spinor_complex_fail =    347
```

此前的 `21,389` 分母无效：其中一部分 case 来自错误 k 数据和错误
little-group membership。修正后当前诊断实际进入 complex-spinor Wigner
分类的 case 数是 `20,810`。

347 个最终失败：

| 阶段 | 数量 |
|---|---:|
| `non_quantized` | 180 |
| `square_not_in_spin_table` | 89 |
| `antiunitary_spin_lookup` | 72 |
| `square_outside_little_group` | 6 |

按旧 `find_setting_transform` 是否返回结果拆分：

```text
antiunitary_spin_lookup   xf_found=true   72
non_quantized             xf_found=true  180
square_not_in_spin_table  xf_found=false  89
square_outside_little_group xf_found=false 6
```

因此“485 个是 k convention、122 个是 setting、72 个单独 lookup”的旧分类已经作废。

### 尚未修复：不要误认为已经接入生产

1. `compute_corepresentation()` 的 spinor 路径仍向
   `wigner_classify_spinor()` 传 `setting_xf=None`。
2. `find_setting_transform()` 仍不可信：
   - rotation 使用 greedy pairing；
   - modulo-1 方程被当普通实数 Gaussian elimination；
   - identity basis 求解失败时仍返回 identity，并计入 `XF_FOUND`；
   - ambiguous transform 仍直接取 `.first()`。
3. `SettingTransform.basis` 已泛化为 `f64 Mat3`，但这是为 rational
   basis 诊断准备的基础设施，不代表 setting 问题已经解决。
4. `standard_setting_transform()` 和
   `get_space_group_with_magnetic_symmetry()` 当前只用于诊断，没有进入生产分类。
5. direct anti-coset 路径仍可能遍历同一 little co-group rotation 的多个
   Seitz/中心化平移变体；`non_quantized=180` 不能直接解释为物理结果。

### 当前最高优先级线索：磁群操作入口可能选错 Hall setting

这是停止前刚发现的线索，**尚未完成验证，禁止直接当最终结论**。

当前 `get_magnetic_operations()`：

1. 用 `get_first_hall_for_uni()` 扫描 `msgdb_get_uni_candidates(hall)`；
2. 将找到的 Hall 传给 `msgdb_get_spacegroup_operations(uni, hall)`。

但 `msgdb_get_spacegroup_operations()` 自身明确支持 `hall=0`，
表示该 UNI 映射表中的第一个合法 Hall offset。当前扫描逻辑是否等价尚未证明。

异常证据：

```text
UNI187 被 identify_unitary_subgroup() 识别为 SG1，
但 get_magnetic_operations(187) 的 unitary 操作包含：
  I
  mirror_y with t=(1/2,0,0)
```

这不可能是 SG1 的闭合 unitary subgroup。说明至少有一处出错：

- Hall 选择错误；
- `timerev` 语义/提取错误；
- UNI→operation table offset 错误；
- 或 subgroup identification 使用了错误操作集。

DeepSeek 的第一步应是：

1. 对 UNI187、UNI270/271、UNI663 比较：
   - 当前 `get_first_hall_for_uni()` 返回的 Hall；
   - `msgdb_get_spacegroup_operations(uni, 0)`；
   - UNI mapping 表中的 `first_hall`；
   - unitary 操作是否闭合；
   - 识别出的 H/G 是否与 BNS/UNI 元数据一致。
2. 在这个入口问题确认前，不要继续调整 Wigner phase、SU(2) lift 或 setting solver。

### spglib standard-setting 诊断现状

新增测试：

```bash
cargo test --package cryspglib diagnose_spglib_standard_setting_transform -- --nocapture
```

当前结果：

```text
total               = 1644
found               = 1644
sg_match            = 1644
detected_hall_exact = 1450
data_hall_exact     = 1450
```

194 个不 exact 的 case 很可能不是 affine transform 本身失败，而是输入的
“unitary”操作已经异常。UNI187 是明确例子。因此不要先修这 194 个 transform。

### 测试状态

- `cargo check --package cryspglib`：通过（现有 warnings 未清理）。
- `diagnose_spin_lg_k_convention`：通过，3611/3611。
- `diagnose_wigner_sources`：通过，当前失败 347。
- `diagnose_spglib_standard_setting_transform`：通过，输出上述 1450/1644。
- `cargo test --package cryspglib --no-run`：当前仍失败，因为
  `examples/sg159_lpoint.rs` 调用 `wigner_classify_spinor()` 时缺少新增的
  `Option<&SettingTransform>` 参数。此问题尚未修。

### 建议接手顺序

1. 验证并修正 UNI→Hall→磁群操作入口，加入 unitary subgroup closure 回归测试。
2. 重新运行 Wigner 基线；旧的 347 数量可能再次变化。
3. 分别求 H 与 parent G 的标准 setting，禁止把 H transform 用于 G spin table。
4. 用修正后的 setting 重新检查 72 个 `antiunitary_spin_lookup`。
5. 对剩余 `non_quantized` 检查 anti-coset 是否按 little co-group 去重，并逐 term
   验证 Bloch phase 和 SU(2) central sign。
6. 只有 `spinor_complex_fail=0` 且全量测试通过后，才能声明 100%。

---

## 2026-06-23 更新：centering-aware delta_in_lattice 修复

### 当前基线

```text
spinor_complex_fail = 1,328 (全部 square_not_in_spin_table，全部 xf_found=true)
```

相比修复前的 2,698 减少了 1,370 (50.8%)。

### 修复内容：`delta_in_lattice` 替换为 centering-aware 版本

#### 问题

旧的 `delta_in_lattice` 对 `canonical_pure_translations` 做 brute-force
integer 组合枚举 (n_i ∈ [-3,3], max 4 vectors)，只能检查 Z³ 等价。
对于 I (½,½,½)、F ((0,½,½)/(½,0,½)/(½,½,0))、C、A、R 中心化群，
centering vectors 缺失导致所有半整数平移差被拒绝。

更根本的问题：`ops_from_hall` 中**所有群都只有 1 个恒等操作 (0,0,0)**，
不包含 centering 变体。翻译格点不能从操作列表中提取，必须根据 SG 的
Bravais 类型推导。

#### 新 `delta_in_lattice(delta, centering_shifts)`

- Z³ 隐式处理（所有平移差 mod Z³）
- `centering_shifts` 枚举所有 centering vector 组合（mod 1 去重）
- 检查 `delta` 的小数部分是否匹配任何 centering 组合
- 对于 F-centering：最多 2³=8 种组合，对于 P：0 种（trivially Z³ only）

#### `centering_shifts_for_sg(sg: u8) -> &'static [[f64; 3]]`

match 语句覆盖全部 81 个非原始群的 centering 向量：

| Centering | SGs | Shifts |
|-----------|-----|--------|
| I (Body) | 23,24,44-46,71-74,79-80,82,87-88,97-98,107-110,119-122,139-142,197,199,204,206,211,214,217,220,229-230 | (½,½,½) |
| F (Face) | 22,42-43,69-70,196,202-203,209-210,216,219,225-228 | (0,½,½), (½,0,½), (½,½,0) |
| C | 5,8-9,12,15,21,35-37,63-68 | (½,½,0) |
| A | 20,38-41 | (0,½,½) |
| R (hex) | 146,148,155,160-161,166-167 | (⅔,⅓,⅓), (⅓,⅔,⅔) |
| P | 其余 149 | (none) |

#### 调用链改动

1. `find_sq_spin_lg_first` 参数从 `canonical_pure_translations` 改为 `centering_shifts`
2. `wigner_classify_spinor_direct_anti_diagnostic` 的 `all_canon` 替换为
   `centering_shifts_for_sg(ctx.sg)`，移除 Z³ 和 `canonical_pure_translations`
3. `wigner_classify_spinor_primary` 的第二调用点同样使用 `centering_shifts_for_sg(ctx.sg)`

### 剩余 1,328 个失败：原始群 origin shift 问题

**全部来自原始群 (P)**，如 SG85 P4/n, SG86 P4₂/n, SG125 P4/nbm 等。

模式：
- `sq_rot = C2z = [[-1,0,0],[0,-1,0],[0,0,1]]`
- `sq_t = (0, 0, 0)`（after to_bilbao）
- `sp_t = (0.5, 0.5, 0)`（spin table 中的 C2z 条目）
- Delta = (-0.5, -0.5, 0)

**不是 centering 问题**：对原始群，(-0.5, -0.5, 0) 不在 Z³ lattice 中。

**`to_bilbao` 验证**：SG85 origin=[4,1,1]，对于 C2z：
- `(I-R)*origin = [8, 2, 0]` → `t_bilbao = (-8,-2,0) mod 1 = (0,0,0)`
- spin table 期望 `(0.5, 0.5, 0)`，无法通过 origin shift 得到

**可能根因**：
1. `to_bilbao` 只使用 origin shift，不使用 basis matrix（虽然所有 basis=I）
2. `(a₀h)²` 的 Seitz 组合可能缺少 n-glide/screw 等非简单平移
3. spin table 的 Bilbao 约定可能使用了与 spglib 不同的平移标准化规则

**待做**：逐 case trace `(a₀h)²` 的计算过程，对比 spglib 和 Bilbao 的
平移约定。剩余失败群列表：85, 86, 125, 126, 129, 130, 151, 179, 201 等。

### 关键代码位置（更新）

| 功能 | 文件:行 |
|------|--------|
| `centering_shifts_for_sg()` | `wigner.rs:~2208` |
| `delta_in_lattice()` (新版) | `wigner.rs:~2226` (nested inside `find_sq_spin_lg_first`) |
| `find_sq_spin_lg_first()` | `wigner.rs:~2256` |
| Direct anti: centering 构建 | `wigner.rs:~1775` |
| Primary path: centering 调用 | `wigner.rs:~2955` |
| 诊断：Hall 翻译输出 | `corep.rs:~2130` |
