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

### 规则 1: cargo check 成功后立即 commit

每次 `cargo check` 成功后必须立即 `git add -A && git commit`。宁可多几个小 commit，不能累积未提交改动。

**Why:** 未提交代码跨 session 容易因 git checkout 等操作永久丢失（发生过：2163→3307 行的诊断代码全部丢失）。

### 规则 2: 不用 Python 脚本做代码修改

批量 sed/Python 替换容易产生意外后果。代码修改必须用 Edit/Write 工具逐处进行，每处修改后确认正确。

### 规则 3: 不用 type alias 做过渡

`pub type OldName = NewName;` 只是把问题藏起来，缺少可维护性。应该全局替换所有引用，然后删除旧定义。

### 规则 4: 先 use 再去掉 crate:: 前缀

替换类型时，先在文件头部添加 `use crate::NewType;`，再把文件中所有 `crate::NewType` 替换为 `NewType`。

---

## 调试方法论 — 从 spinor Wigner 排查中提炼的经验

### 原则 1：比较 passing vs failing cases，找差异因子

当同一个代码路径**有些 case 通过、有些失败**时，不要猜测通用原因。直接比较一个成功 case 和一个失败 case，问：

> **"这两个 case 之间什么不同，导致一个通过、一个失败？"**

这个差异通常直接揭示根因。

实例：SG2 T-point passed after loop fix → SG159 L-point still failed. Difference: SG2's LG {I, -I} always contains (a₀h)²; SG159's LG {I, mirror} doesn't. Root cause: spinor little group coverage, not algorithm correctness.

### 原则 2：假设驱动的逐层排除

不要在无数可能原因中随机尝试。为每种假设设计一个**最小 oracle test**，一票否决或确认：

1. **列出所有可能假设**（按先验概率排序）
2. **为每个假设设计一个 oracle**（最小代码改动，只输出统计数据）
3. **跑 oracle，看数据**——如果数据否决假设，立刻排除，不再纠结
4. **如果 oracle 确认假设，再设计修复**

实例——NONE=1,007 的排查顺序：
- H1: same-rotation lift 误选 → scan same-rot candidates → OTHER=0 → 否决
- H2: UU* antiunitary square → 6 formulas oracle → 6.5% fix → 不是主因
- H3: H/G gauge mismatch → G-gauge oracle → 0% fix → 否决
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

---

## 错题集 — Spinor Wigner SU(2) 调试记录

### Bug 1: 归一化分母 mismatch（SG2 T-point W=-0.5）

**现象**：SG2 T2/T3 spinor irreps 返回 `None`。per-term 全部通过但 W=-0.5 不量子化。

**根因**：Wigner 公式 `W = Σ χ̃(a₀h) / |H₀|` 的求和对象是 **little co-group**（不同旋转），不是 full little group（所有 Seitz 翻译变体）。旧 loop 遍历了 4 个 Seitz 变体但 spinor 表只有 2 个条目。

**修复**：Loop 改为遍历 `spin_lg_op_indices[0..n_lg_ops]`（co-group 规范代表），归一化分母 = `n_lg_ops`。

**教训**：loop domain 和 character domain 必须一致。

### Bug 2: sq 匹配目标错误（高分群 0% pass）

**现象**：Loop fix 后 SG180/SG148/SG179 等高分群全部 0% pass。

**根因**：`(a₀h)²` 的翻译来自磁群，`h_spin_seitz` 只有规范翻译。Full Seitz matching 因翻译不匹配失败。

**修复**：sq 用 **rotation-only matching**（rotation 在不同 setting 下不变，translation 受 origin 影响）。

**教训**：匹配时区分 invariant（rotation）和 variant（translation）。

### Bug 3: a₀ 选择错误（grey 群用了 θ·g 而非 θ）

**现象**：SG159 L-point 产生 C3 旋转，不在 LG 中。

**根因**：代码取 `antiunitary[0]` 碰巧取到 θ·g（mirror）。对 grey 群必须取纯 θ (R=I)。

**修复**：`select_spinor_a0` helper 显式找 R=I 的反酉操作。

**教训**：不依赖数组顺序隐含的语义。Grey 群的 a₀ 必须是 θ。

### Bug 4: 示例 ctx.g 设置错误（假阳性）

**现象**：SG159.63 L2 Wigner=None

**真因**：示例代码 `ctx.g = h_spin`（SG143, 3 ops）而非 `spin_ops_for_sg(159)`（6 ops）。

**教训**：BlackWhite 群 G≠H，必须用 G 的 spin ops。排查时先确认数据存在再怀疑算法。

### Bug 5: LG-first sq matching 修复（v5, +11 cases）

**现象**：`h_spin_seitz.iter().position(|s| s.rot == sq.rot)` 可能先匹配到不在 `spin_lg_op_indices` 中的 candidate。

**修复**：`find_sq_spin_lg_first()` — LG 内 full Seitz → LG 内 unique rotation → 全局 rotation。

**效果**：88.1% → 88.3%（+11 ok, -11 fail）

### Bug 6: Θ²=Ē 中心元和 antiunitary square convention

**现象**：grey group h=I 时 SU(2) 无法检测 Θ²=Ē。`(JU)(JU)*` 修复 61% NONE 但 global 替换产生 regression。

**排查过程**（完整假设排除链）：
- Paoli SU(2) closure: 47,486/0 matched → SU(2) 合成本身正确
- Same-rotation scan: OTHER=0 → lift 选择正确
- 6 antiunitary square formulas: UU* only 6.5% → 不是简单 square formula
- H/G gauge mismatch: G-gauge oracle 0% → 不是 gauge mixing
- det distribution: 混合 60/40 → 不是 improper rotation
- J-insertion on NONE: 61% fix → J 是关键
- Global J: 88.3%→83.0% → 不能全局替换
- Case-level J fallback: only 22/945 fix → 不能作为 fallback

**当前结论**：J-insertion `(JU)(JU)*` 确认了 direction（antiunitary square 需要显式 Θ=JK），但不能全局应用。923 个 both_fail + 1,547 个 both_ok_diff_type 需要更深层诊断。

### ✅ 已排除的问题

1. **spin 数据库不完整** ❌
2. **a₀ improper rotation 缺 SU(2) lift** ❌
3. **`(a₀h)²` 超出 little group（grey 群）** ❌
4. **Seitz 翻译变体 double counting** ❌
5. **sq 匹配因翻译不匹配失败** ❌
6. **Pauli SU(2) 合成约定错误** ❌（closure test 验证）
7. **same-rotation lift 误选** ❌（OTHER=0）
8. **UU* antiunitary square formula** ❌（6.5% fix）
9. **H/G gauge mismatch** ❌（G-gauge 0%）
10. **det=-1 improper 独占** ❌（混合分布）
11. **global J-insertion** ❌（regression）

### SU(2) 覆盖率演进

| 版本 | ok | fail | rate | 关键变化 |
|------|-----|------|------|----------|
| v3 (θ fix) | 7,108 | 956 | 88.1% | Bug 1+2+3 全部修复 |
| v5 (LG-first sq match) | 7,119 | 945 | 88.3% | `find_sq_spin_lg_first` |
| v6 (eta_ebar) | 7,119 | 945 | 88.3% | always missing, no gain |
| J-left (global) | 6,690 | 1,374 | 83.0% | reverted |
| J-left (per-case) | — | 22/945 | — | too few, reverted |

### 当前失败分类 (945 total, 2026-05-12)

| 类别 | 数量 | 说明 |
|------|------|------|
| `sq_in_lg_but_su2_fail` | 883 (93%) | sq 在 LG 内，SU(2) 或 W 失败 |
| `sq_outside_lg` | 54 (6%) | sq 在数据库但不在 LG |
| `sq_not_in_spin` | 8 (1%) | sq rotation 不在 H spin ops |

### SU(2) central-detection relation (per-term)

| 关系 | 数量 | 占比 | 含义 |
|------|------|------|------|
| SAME | 31,357 | 24% | u_sq ≈ u_k (same lift) |
| EBAR | 98,356 | 75% | u_sq ≈ -u_k (Ebar lift) |
| NONE | 1,007 | 1% | unrelated → function returns None |

NONE 1,007 的 sub-category: 全部 no other same-rotation candidate (not a lift-selection issue).

### 当前 wigner.rs 工具函数

| 函数 | 用途 |
|------|------|
| `find_sq_spin_lg_first()` | LG-first sq matching |
| `infer_eta_ebar()` | Central parity inference (always missing) |
| `conj_pauli()`, `neg_pauli()` | Pauli operations |
| `antiunitary_square_pauli()` | J-left antiunitary square |
| `SquareKernel` enum | Pluggable square kernel (OldU2 / JLeft) |
| `find_spin_in_db()` | Spin DB lookup with -R fallback |
| `su2_compose()`, `su2_same_up_to_sign()` | Core SU(2) operations |

---

## Architecture overview

cryspglib has two major subsystems:

**1. spglib port** — space group identification from crystal structures.
`Crystal::new(lat, positions, types)` → `.analyze()` → `.dataset()` → `SpaceGroup`.

**2. irrep module** — irreducible representation data for all 230 space groups.
`irreps_of(sg_number)` → `IrrepRecord` (labels, characters, matrices, isotropy subgroups, magnetic corepresentations).

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

### Data regeneration pipeline

```bash
cargo test --package cryspglib --lib \
  irrep::corep::tests::export_hall_operations -- --nocapture
python3 scripts/generate_irrep_data.py
# or: bash scripts/regenerate_all.sh
```

---

## Key types

| Type | Location | Description |
|------|----------|-------------|
| `Crystal` | `api.rs` | Entry point: lattice + positions + types |
| `SymmetryOps` | `api.rs` | Ordered set of `{R|t}` + time_reversal |
| `SpaceGroup` | `lib.rs` | SG number, Hall number, ops |
| `IrrepRecord` | `irrep/types.rs` | irrep: labels, dim, k-vector, characters, matrices, subgroups, corepresentations |
| `SpinLiftContext` | `irrep/wigner.rs` | H and G spin ops for Wigner test |
| `SeitzOp` | `irrep/wigner.rs` | `{R|t}` with optional time reversal |
| `CorepType` | `irrep/corep.rs` | A/B/C/Unsupported |

---

## Module structure

| Module | Role |
|--------|------|
| `api.rs` | `Crystal`, `SymmetryOps`, `find_hall_number` |
| `irrep/types.rs` | `IrrepRecord`, `IsotropyRecord` |
| `irrep/query.rs` | `irreps_of()`, `kpoints_of()`, `format_character_table()` |
| `irrep/corep.rs` | Co-representation: `compute_coreps()`, `CorepType`, diagnostic tests |
| `irrep/wigner.rs` | Wigner test: Seitz composition, SU(2) composition, spinor classification |
| `irrep/wigner_extra.rs` | Pre-computed antiunitary character path (included into wigner.rs) |
| `irrep/bridge.rs` | `impl SpaceGroup` — bridge APIs |
| `irrep/generated_data.rs` | Auto-generated static arrays (~20 MB) |

---

## Test suite (177 tests)

| Layer | Count | Description |
|-------|-------|-------------|
| Wigner algorithm | 11 | Seitz composition, k filtering, Type A/B/C |
| Diagnostic | 15+ | `diagnose_wigner_sources`, `diagnose_per_term_su2_trace`, etc. |
| Full self-consistency | 8 | χ(E)=dim, k-point grouping, isotropy subgroup validity |
| BCS validation | 5 | Known reference cases |
| Integration | 6 | bcs_corep_validation, irrep_validation, etc. |
