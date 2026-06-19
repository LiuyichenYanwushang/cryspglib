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

## 错题集 — Spinor Wigner SU(2) 调试记录

### Bug 1: 归一化分母 mismatch（SG2 T-point W=-0.5）

**现象**：SG2 T2/T3 spinor irreps 返回 `None`。per-term 全部通过但 W=-0.5 不量子化。

**根因**：Wigner 公式 `W = Σ χ̃(a₀h) / |H₀|` 的求和对象是 **little co-group**（不同旋转），不是 full little group（所有 Seitz 翻译变体）。旧 loop 遍历了 4 个 Seitz 变体但 spinor 表只有 2 个条目。

**修复**：Loop 改为遍历 `spin_lg_op_indices[0..n_lg_ops]`（co-group 规范代表），归一化分母 = `n_lg_ops`。

**教训**：loop domain 和 character domain 必须一致。

### Bug 2: sq 匹配目标错误（高分群 0% pass）

**现象**：Loop fix 后 SG180/SG148/SG179 等高分群全部 0% pass。

**当时确认的问题**：`(a₀h)²` 的翻译来自磁群，`h_spin_seitz` 只有规范翻译。
Full Seitz matching 会因 origin/translation representative 不同而失败。

**当时修复**：sq 使用 **rotation-only matching**，绕过纯 origin shift 导致的
translation representative 差异。

**重要更正（2026-06-19）**：rotation 只在“同一基底、仅 origin 不同”时保持不变。
若存在基底变换 `x_new = T x_old + s`，必须使用
`R_new = T R_old T⁻¹`；rotation-only matching 不能跨基底直接使用。

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
- H/G gauge mismatch: 当时 oracle 没有先对齐 MSG/H Hall 基底，因此“0%”不能排除
  setting/gauge mixing；UNI663 后来提供了直接反例
- det distribution: 混合 60/40 → 不是 improper rotation
- J-insertion on NONE: 61% fix → J 是关键
- Global J: 88.3%→83.0% → 不能全局替换
- Case-level J fallback: only 22/945 fix → 不能作为 fallback

**当前结论**：J-insertion `(JU)(JU)*` 确认了 direction（antiunitary square 需要显式 Θ=JK），但不能全局应用。923 个 both_fail + 1,547 个 both_ok_diff_type 需要更深层诊断。

### 历史排除列表（2026-06-19 已复核）

1. **spin 数据库不完整** ❌
2. **a₀ improper rotation 缺 SU(2) lift** ❌
3. **`(a₀h)²` 超出 little group（grey 群）**：纯 grey-group `a₀` 选择问题已修，
   但 Type-III setting mismatch 仍可制造假的 outside-LG
4. **Seitz 翻译变体 double counting** ❌
5. **sq 匹配因翻译不匹配失败** ❌
6. **Pauli SU(2) 合成约定错误** ❌（closure test 验证）
7. **same-rotation lift 误选** ❌（OTHER=0）
8. **UU* antiunitary square formula** ❌（6.5% fix）
9. **H/G gauge mismatch**：**未排除**。旧 G-gauge oracle 在未对齐基底的操作上比较，
   结论无效；当前剩余失败的首要已证实根因正是 setting/gauge mixing
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

## Spinor Wigner：2026-06-18 历史快照（已过期）

> **不要把本节数字和结论当作当前状态。** 本节保留用于追踪排查过程。
> 其中“rotation 在不同 setting 下不变”“H/G gauge mismatch 已排除”
> 和“genuine nonquantized”等判断已经被 2026-06-19 的 UNI663 oracle 否定。
> 当前权威状态见下一节。

### 结论

此前“SU(2) 数据库缺少逐操作 central parity”不是已确认根因。真正确认的首要问题在
`scripts/parse_spinor_data.py`：**spin.dat 的复数字符被错误解析**。

spin.dat 的 irrep 行有两种格式：

```text
n 个实字符
n 个幅值 + n 个相位（相位单位为 π）
```

旧解析器把第二种格式的前半段直接当作实字符，把后半段误认为
`extra Wigner chars`。因此所有非零相位被丢失。例如：

```text
SG3 A3:  1.0 1.0  0.0 -0.5
```

实际字符是 `χ=(1,-i)`，旧数据却保存为实字符 `[1,1]` 和所谓 extra
`[0,-0.5]`。SG195 的 C₃ 字符同样应为
`-1/2 ± i√3/2`，旧数据只保留了错误的幅值 `1`。这直接解释了失败为何集中在
含 C₂+C₃ 组合的立方群。

### 已完成并验证

- `parse_spinor_data.py` 已按“幅值 + 相位/π”解码复数字符。
- 生成数据中：
  - `CHARACTERS` 保存 spinor 字符实部；
  - `SPIN_EXTRA_CHARS` 现保存字符虚部（旧名称保留是为了减少数据结构改动）；
  - `IrrepRecord::spin_character_imag()` 提供虚部。
- `wigner_classify_spinor()` 的求和项已使用完整复数字符。
- 已重新生成 `src/irrep/generated_data.rs`。
- `cargo check --package cryspglib` 在加入 MSG-gauge 修复后通过（仅有既有 warnings）。
- 针对性测试 `test_spinor_sg3_a3_grey_wigner` 在复数字符修复后通过。

### 第二个问题：坐标 setting/gauge 混用

代码审查还发现 `wigner_classify_spinor()` 虽然接收了 `unitary_mag_indices` 并构造
`h_to_spin`，旧正式路径却没有使用它们。旧计算实际把：

- MSG/parent setting 中的 `a₀`
- standalone H spin table setting 中的 `h`

直接做 SU(2) 合成。非平凡 subgroup embedding 下这不是同一个坐标 frame，立方群
最容易暴露问题。

当前工作树已加入一个 MSG-gauge 路径：

1. 在 MSG setting 中选取 `a₀` 和每个 unitary `h`；
2. 在 parent G 的 spin gauge 中合成 SU(2) lift；
3. 将 `(a₀h)²` 映射回 H 的 canonical spin operation；
4. 显式加入 spin-½ 的 `Θ²=Ē=-I`；
5. grey group 使用纯 `Θ` 作为 `a₀`。

### 验证状态

- 历史基线仍是 **7,119/8,064 = 88.3%**，945 个失败；这个数字来自错误的实字符数据，
  只能作为修复前基线，不能再代表当前实现。
- 复数字符修复已完成并通过单例测试。
- MSG-gauge/`Θ²` 路径已通过编译检查，但**尚未完成针对性和全量 Wigner 统计**。
- 因此当前不能宣称已经达到 100%。下一步必须先运行：

```bash
cargo check --package cryspglib
cargo test --package cryspglib test_spinor_sg3_a3_grey_wigner -- --nocapture
cargo test --package cryspglib diagnose_wigner_sources -- --nocapture
cargo test --package cryspglib
```

### 后续清理（已完成 2026-06-18）

- ✅ `spin_extra_chars()` → `spin_character_imag()`，内部字段 `_spin_extra_*` → `_spin_imag_*`，
  生成数组 `SPIN_EXTRA_CHARS` → `SPIN_IMAG_CHARS`，Python 生成器同步更新
- ✅ `diagnose_wigner_sources` 重写：按 `mapping_failure / complex_char_non_quantized / 
  real_char_su2_closure_fail` 三类根因重新分类
- ✅ 添加路径 triage 计数器（`MSG_GAUGE_OK/MAP_FAIL/W_FAIL`，`OLD_PATH_OK/FAIL`）
- ✅ 添加 `build_h_to_spin_map` 失败分类（`H2S_OK/AMBIGUOUS/MISSING`）

### 当前覆盖率（2026-06-18）

全量 21,389 个 spinor irrep×UNI 对（注意：不是去重的 8,064 个 unique spinor case）：

| 类别 | 数量 | 说明 |
|------|------|------|
| `spinor_complex_ok` | 18,225 | 复数字符 Wigner 通过 |
| `spinor_complex_fail` | 3,164 | 剩余失败 |
| `scalar_PIR` | 23,313 | 标量 PIR（不适用 spinor） |
| `scalar_trivial_A` | 9,647 | 无反酉操作（平凡 Type A） |
| `scalar_CIR` | 3,281 | 复合 CIR |

路径分布：

| 路径 | 数量 | 说明 |
|------|------|------|
| `MSG_GAUGE_OK` | 26,334 | MSG-gauge 直接成功 |
| `MSG_GAUGE_MAP_FAIL` | 5,278 | spin→mag 映射不完整（非 fatal，仅诊断） |
| `MSG_GAUGE_W_FAIL` | 875 | 映射完整但 W 不量子化 |
| `OLD_PATH_OK` | 287 | 旧 fallback 挽回 |
| `OLD_PATH_FAIL` | 583 | 两路全败 |

### `spin→mag` 映射修复（2026-06-18）

**根因**：Type-3（黑白）磁群 G = H' ∪ a₀H' 中，父群 H 包含反酉操作。
`spin_lg_op_indices` 覆盖 H 的完整 little co-group，但 Wigner 求和只需 H₀ ∩ H'
（酉操作）。旧代码要求所有 spin_lg_op 都有酉 MSG 映射 → 5,278 个 case 被拒绝。

**修复**（`wigner.rs:1423-1503`）：不强制全部映射。Wigner 循环中跳过无法映射的条目，
用 `n_mapped`（实际酉 little co-group 大小）而非 `n_lg_ops` 归一化。
效果：**-1,523 失败**（4,687 → 3,164）。

**`central = !spatial_central` 正确性验证**：
- Θ² = Ē = -I（自旋 1/2）
- U_{(Θ·h)²} = Ē · U_{spatial_square}
- spatial=EBAR → Ē² = I → SAME → central=false
- spatial=SAME → Ē → EBAR → central=true
- ∴ central = !spatial_central ✓

### `build_h_to_spin_map` 诊断结论

| 计数器 | 值 | 含义 |
|--------|-----|------|
| `H2S_OK` | 342,606 | 成功映射 |
| `H2S_AMBIGUOUS` | **0** | 不存在"同一旋转多个 spin 条目"问题 |
| `H2S_MISSING` | 151,702 | 旋转不在 spin table 的 little co-group 中——这些是 H\H' 中的反酉操作 |

### 剩余 3,164 个失败

- **875 MSG_GAUGE_W_FAIL**：映射完整但 W 不量子化。典型 case（SG24 W 点）：所有
  (a₀h)² 映射到 identity，W 退化仅依赖 Bloch 相位 → genuine 非量子化
- **583 OLD_PATH_FAIL**：MSG-gauge 和旧 fallback 都无法处理
- 需要逐 case 物理分析，不宜盲目修代码

---

## 重要审计更正：2026-06-19 Phase 1 不能视为完成

> **给后续 AI 的强制阅读说明：**
> 最近提交 `acb3b44`、`4a5956b`、`a7ef97f`、`b253cb1`、`678add3`
> 和文档提交 `3660d44` 中，关于 setting transform 的排查方向部分正确，
> 但验证标准、实现完整性和最终结论存在严重问题。
> **不要继续沿用“Phase 1 已完成”“136/136 已修复”的结论。**
> 应先修正本节列出的问题，再进入 reciprocal-k 等后续阶段。

### 仍然成立的结论

UNI663 是一个有效反例，已经证明以下事实：

- `ops_from_msg` 中的 unitary H 仍处于 MSG/parent 的嵌入基底；
- `ops_from_hall` 和 standalone H spin table 处于 canonical Hall 基底；
- `spg_get_hall_number_from_symmetry()` 只做 Hall 分类，不负责把输入操作变换到
  canonical Hall setting；
- 直接跨这两个基底比较 rotation 会产生假的
  `square_not_in_spin_table` / `square_outside_little_group`；
- 一般 setting 变换必须满足
  `R_hall = T R_msg T^-1`，不能声称 rotation 在所有 setting 下不变。

因此，“必须显式求解并验证 MSG→Hall 的 setting transform”这个方向是正确的。
错误在于当前代码并没有完成这个求解和验证。

### 错误结论 1：“136/136 fixed”只表示失败阶段发生转移

`phase1b_verify_transform_fix` 的所谓 fixed 判据是：

```rust
!matches!(
    fixed_result,
    Err(DirectAntiFailure::SquareNotInSpinTable)
)
```

也就是说，只要结果不再停在 `SquareNotInSpinTable`，即使变成
`SquareOutsideLittleGroup`、`AntiunitarySpinLookup`、`Su2LiftMismatch` 或
`NonQuantized`，也会被计为 fixed。

2026-06-19 复跑 `diagnose_wigner_sources` 得到：

```text
spinor_complex_ok   = 20,710
spinor_complex_fail =    679

Final failure stages:
non_quantized                  315
square_outside_little_group    184
antiunitary_spin_lookup         64
square_not_in_spin_table        58
su2_lift_mismatch               58
```

与接入 setting transform 前记录的三类失败比较：

```text
before:
non_quantized                  315
square_not_in_spin_table       194
square_outside_little_group    170
total                          679

after:
non_quantized                  315
square_not_in_spin_table        58   (-136)
square_outside_little_group    184   (+14)
antiunitary_spin_lookup         64   (+64)
su2_lift_mismatch               58   (+58)
total                          679
```

`14 + 64 + 58 = 136`。因此这 136 个 case 没有被最终修复，只是从第一道
rotation lookup 失败转移到了后续失败阶段。

**正确表述应是：**

> 对原先 136 个 `square_not_in_spin_table` case，当前候选 T 能让 transformed square
> 通过 H spin-table rotation lookup；但 136 个 case 全部仍在后续阶段失败。

不得再称其为“136/136 修复”，也不能据此证明 T 是正确的完整 setting transform。

### 错误结论 2：当前实现没有真正求出 `(T,s)`

原 Phase 1 计划要求：

1. rotation conjugacy；
2. 使用全部操作联立求 origin `s`；
3. 完整 Seitz set 双射；
4. operation correspondence；
5. 多解/群自同构歧义检测。

当前 `find_setting_transform(msg_rots, hall_rots)` 实际只做：

- 枚举 48 个 signed-permutation `T`；
- 比较 `T R_msg T^-1` 与 Hall rotation multiset；
- 找到第一个候选立即返回；
- `origin` 永远固定为 `[0,0,0]`；
- 不输入或比较 translation；
- 不建立 MSG op → Hall op 的双射；
- 不检测多个候选；
- 不验证 irrep character operation correspondence。

rotation multiset 匹配本身还存在必然歧义：若 `T` 是候选，则 `-T` 给出相同共轭，

```text
(-T) R (-T)^-1 = T R T^-1
```

高对称群还可能有更多 normalizer/automorphism 候选。因此“取枚举顺序中的第一个 T”
不能证明它是正确的 setting transform。

`SettingTransform` 虽然包含 `origin` 字段，但当前所有构造点都把它设为零。
所以代码目前只实现了一个 **rotation-set candidate oracle**，不是完整 `(T,s)` 求解器。

### 实现问题 1：没有接入真正的生产 API 路径

`compute_corepresentation()` 是 `compute_coreps()` 和
`IrrepRecord::corepresentation()` 使用的实际路径。它当前调用
`wigner_classify_spinor()` 时仍显式传入：

```rust
setting_xf = None
```

只有 `diagnose_wigner_sources` 等诊断测试计算并传入 `setting_xf`。
因此提交 `678add3 feat: wire setting_transform into wigner_classify_spinor primary path`
的描述不准确：

- 函数签名和 fallback 参数已经接线；
- 诊断路径传入了候选 T；
- 实际公共 corepresentation 路径没有计算或传入 T；
- `wigner_classify_spinor_primary()` / MSG-gauge primary 本身也没有使用 `setting_xf`；
- 只有 primary 返回 `None` 后的 direct anti-coset fallback 使用该参数。

后续 AI 必须区分“诊断路径统计”和“用户实际 API 行为”，不能用前者代表后者。

### 实现问题 2：G spin 的 `-R` fallback 混用了变换前后的基底

direct anti-coset 中，首选 G-spin lookup 使用 transformed `b_rot`，但 improper
rotation 的 `-R` fallback 仍从原始 `b.rot` 构造：

```rust
position(|s| s.rot == b_rot)
    .or_else(|| {
        let r = -b.rot; // 错：仍使用未变换的 rotation
        ...
    })
```

两条分支不在同一 coordinate frame。至少从实现一致性看，fallback 应使用
`-b_rot`。但更根本的问题是：parent G spin table 是否与 transformed H Hall frame
处于同一 setting 尚未验证，不能只改这一行后宣称 gauge 问题解决。

### 实现问题 3：当前失败分类文档已经过期

旧文档把 679 个失败写成：

```text
non_quantized=315
square_not_in_spin_table=194
square_outside_little_group=170
```

这是接入 rotation candidate 前的分类。当前总数仍是 679，但已经分成五类：

```text
non_quantized=315
square_outside_little_group=184
antiunitary_spin_lookup=64
square_not_in_spin_table=58
su2_lift_mismatch=58
```

后续分析必须使用同一轮、同一代码版本的完整 failure-stage vector，不能只比较其中
一个阶段下降，也不能因为总失败数未变就声称被修改的 case 已经由 primary path 修复。

### 测试质量问题

当前两个 Phase 1 oracle 都是打印型测试，没有断言：

- 预期 UNI 数量；
- full Seitz match 数量；
- ambiguous 数量；
- 最终 `Ok(CorepType)` 增量；
- 679 是否下降；
- 原有成功 case 的类型是否保持不变。

此外，`test_spinor_wigner_gauge_limitation_msg197_8` 中存在恒真断言：

```rust
assert!(result.is_some() || result.is_none())
```

这不能验证任何行为，应恢复为具体预期或删除。

所以 `cargo test` 显示这些测试通过，只能说明代码执行完毕，不能证明 Phase 1 结论成立。

### 后续 AI 必须采用的修正顺序

1. **撤销结论，不必立即撤销代码：**
   将当前 `find_setting_transform` 明确视为 rotation-only candidate finder，
   不得视为生产可用的 setting transform。

2. **先建立可比较基线：**
   每轮记录完整 final failure vector：

   ```text
   ok
   fail
   non_quantized
   square_not_in_spin_table
   square_outside_little_group
   antiunitary_spin_lookup
   su2_lift_mismatch
   其他所有 DirectAntiFailure
   ```

3. **真正求解并验证 `(T,s)`：**
   - 枚举或构造候选 `T`；
   - 为每个 rotation 建立所有可能的 operation correspondence；
   - 从全部 Seitz translation 方程求 `s mod Z^3`；
   - 用变换后的完整 Seitz multiset 验证双射；
   - 报告 zero/one/multiple solutions；
   - 多解时不能随意取第一个，必须用 translation、generator correspondence、
     k-star 或 character-table operation correspondence继续消歧。

4. **不要同时混用 MSG frame、canonical H frame 和 parent G spin frame：**
   明确维护：
   - `b_msg`：MSG/parent frame；
   - `b_hall`：canonical H frame；
   - `k_msg` / `k_hall`；
   - H character lookup 使用的 operation；
   - G-side SU(2) central parity 使用的 operation。

5. **修正 little-group 前先验证坐标约定：**
   比较 `R k` 与 `R^-T k` 两套 oracle，并和 `spin_lg_op_indices` 的实际 rotation
   集合做全量 exact match。不要凭单例直接改公式。

6. **用最终结果而不是阶段通过作为验收：**
   一个 case 只有在返回正确、量子化的 `Ok(CorepType)` 时才算 fixed。
   从一种 `DirectAntiFailure` 变成另一种不算修复。

7. **生产路径和诊断路径必须同时验证：**
   增加至少一个通过 `compute_corepresentation()` 或 `compute_coreps()` 的
   nontrivial-setting regression，不能只直接调用内部 Wigner helper。

8. **加入真实断言：**
   - UNI663 的完整 `(T,s)` 和 Seitz bijection；
   - candidate ambiguity 数量；
   - final fail 严格下降；
   - 现有 20,710 个成功 case 不减少且类型不改变；
   - public API 与诊断路径结果一致。

### 当前准确结论

截至提交 `3660d44`：

- setting mismatch 是已证实的真实根因之一；
- UNI663 证明 rotation 必须经过 basis conjugation；
- signed-permutation rotation candidate 能让 136 个 case 越过
  `square_not_in_spin_table` 阶段；
- 但这 136 个 case 全部仍失败，总失败数保持 679；
- origin、translation、operation correspondence、candidate ambiguity、
  reciprocal-k 和 G/H SU(2) gauge 均未解决；
- setting transform 尚未接入实际公共 corepresentation 路径；
- 因此 **Phase 1 未完成，不能进入“在已完成 Phase 1 基础上修 Phase 2”的状态。**

---

## Spinor Wigner：2026-06-19 当前状态与错误复盘

### 当前覆盖率

目标统计口径是 `diagnose_wigner_sources` 中全部 **21,389 个
spinor irrep × UNI 对**。这和历史上去重后的 8,064 个 unique spinor case
不是同一个分母，后续禁止混用。

| 类别 | 数量 | 比例 |
|------|------|------|
| `spinor_complex_ok` | **20,710** | **96.825%** |
| `spinor_complex_fail` | **679** | **3.175%** |
| 合计 | 21,389 | 100% |

当前仍未达到 100%。剩余 679 个失败按 direct anti-coset 最终阶段分类：

| 失败阶段 | 数量 | 当前解释 |
|----------|------|----------|
| `non_quantized` | 315 | W sum 未量子化；这是症状，不能称为物理上的 genuine non-quantized |
| `square_not_in_spin_table` | 194 | `b²` 的 rotation 无法在 canonical H spin table 中找到 |
| `square_outside_little_group` | 170 | `b²` 匹配到 H spin op，但不在当前 canonical spin LG 索引中 |

### 已完成修复及量化效果

#### 1. spin.dat 复数字符解析

spin.dat 的复数字符格式是 `n 个幅值 + n 个相位/π`。旧解析器把后半段误认为
`extra Wigner chars`，丢失了所有非零相位。现已按复数解码，并通过
`spin_character_imag()` 暴露虚部。

#### 2. 诊断统计改为全量、分阶段、作用域清晰

旧 `diagnose_wigner_sources` 混用了多轮 pass 的静态 counter，并且失败 triage
不能保证覆盖全部失败。现已：

- 每轮统计前 reset counter；
- 第一轮只给总量，第二轮穷举全部失败；
- 按 failure stage、SG、mapping shape 交叉统计；
- direct anti-coset oracle 和正式结果分别计数。

#### 3. `SPIN_LG_OP_INDICES` 的 local/global 索引错误

这是本轮确认并修复的数据生成 bug。

- `SPIN_OP_*` 在每个 SG 内从 Bilbao/ISOTROPY 顺序重排到 Hall 顺序；
- `spin_lg_op_indices` 的值本来是 **SG-local spin-op index**；
- 旧 `generate_irrep_data.py` 却把它当成全局 index，减去 SG 的全局 start 后再映射；
- 结果是大部分 LG index 仍停留在旧顺序，与已重排的 `SPIN_OP_*` 不一致。

正确做法是维护 `old_local_index -> new_local_index`，不参与跨 SG 的 global offset。
重新生成数据后：

```text
fail: 3,164 -> 2,207
```

相关提交：`aee2b8a fix: keep spin little-group indices SG-local after Hall reorder`

#### 4. direct anti-coset fallback

正式 primary path 返回 `None` 时，直接遍历实际反酉 little-group 操作
`b ∈ M_k \ H_k`，计算 `b²` 和 `χ(b²)`，避免依赖脆弱的 `a₀h` 重构。

期间修复了两点：

- character 必须使用完整复数；
- Bloch phase 必须包含 square reduction 以及 canonical representative 之间的
  fractional translation 差，不能只处理整数 lattice shift。

fallback 只在 primary path 失败时执行，因此不改变已有成功 case。结果：

```text
fail: 2,207 -> 679
rescued: 1,528
```

相关提交：

- `bf007c4 fix: include fractional translation phases in direct Wigner sum`
- `a5dbd05 fix: fall back to direct anti-coset spinor Wigner sum`

### 已证实根因：MSG/H canonical setting 不一致

`identify_unitary_subgroup_with_hall()` 返回两套操作：

- `ops_from_msg`：从 MSG 数据库直接抽出的 unitary H，仍在 **MSG/parent 嵌入基底**；
- `ops_from_hall`：根据识别出的 Hall number 重建的 **canonical H Hall 基底**。

`spg_get_hall_number_from_symmetry()` 只负责分类并返回 Hall number，**不会把输入操作
变换到该 Hall setting**。因此 `hall=...` 不意味着 `ops_from_msg` 已经
“Hall-corrected”。源码中类似 `ops_from_msg // correct Hall setting` 的注释是错误的。

#### UNI663 直接反例

`debug_direct_anti_setting_uni663` 给出：

```text
UNI663: parent G = SG75
identified H = SG3, Hall3

ops_from_msg 的非平凡 H rotation: C2z = diag(-1,-1, 1)
Hall3 / H spin table 的非平凡 rotation: C2y = diag(-1, 1,-1)

反酉操作: C4z, C4z^-1
它们平方得到: C2z
```

所以 `b²=C2z` 在 MSG 嵌入的 H 中完全正确，但直接拿它查询 canonical SG3 Hall3
spin table 时，表中只有 `C2y`，于是产生 `square_not_in_spin_table`。
这不是 spin 数据缺失，也不是群闭包失败，而是两套基底被直接比较。

### 之前排查中的错误点

#### 错误 1：把 Hall 分类结果当成坐标变换结果

错误假设：

```text
spg_get_hall_number_from_symmetry(...) 返回 Hall3
=> 输入 ops 已经处于 Hall3 canonical setting
```

实际只完成了分类。必须另外求 `(T,s)` 并显式变换操作。

#### 错误 2：声称 rotation 在不同 setting 下不变

这只对“同一基底，仅 origin shift 不同”成立。一般 setting 变换
`x_can = T x_msg + s` 下：

```text
R_can = T R_msg T^-1
t_can = s - R_can s + T t_msg
```

UNI663 的 `C2z -> C2y` 就是反例。rotation-only matching 只能在已确认基底一致后使用。

#### 错误 3：用无效 oracle 排除 H/G gauge mismatch

旧 G-gauge oracle 报告 0% gain 后，文档曾把 gauge mismatch 标为“已排除”。
但该 oracle 本身仍直接比较 MSG 基底和 canonical H 基底，没有先做 setting transform，
所以它不能检验目标假设。这个排除结论无效。

#### 错误 4：把 `non_quantized` 解释成 genuine 物理结果

旧文档把 SG24 W 等 case 称为 “genuine nonquantized”。Wigner indicator 对适用 case
应量子化；在映射/基底尚未对齐时，非量子化首先是实现错误的信号。当前 315 个
`non_quantized` 只能作为 failure stage，不能作为物理结论。

#### 错误 5：混用统计口径和 counter scope

- 8,064 是历史 unique case；
- 21,389 是当前全量 irrep×UNI case；
- `MSG_GAUGE_*` 等 counter 曾累积多轮 pass，数值不能和单轮总量直接比较。

后续所有覆盖率以同一轮 `diagnose_wigner_sources` 的 21,389 为分母。

#### 错误 6：direct path 曾遗漏 fractional translation phase

只加入整数 `lattice_sq` 不足以匹配 centering/nonsymmorphic representative。
必须同时加入 computed square 与 canonical spin op 的 fractional translation 差。

#### 错误 7：central parity 仍存在跨 gauge 比较风险

当前 direct fallback 用 parent G spin table 找 `U_b`，但用 canonical H spin table
找 `U_{b²}`，然后比较 `U_b²` 与 `U_{b²}`。若 G/H 的轴或 lift gauge 未对齐，这仍是
跨 gauge 比较。更稳妥的做法是：

- central parity 完全在 parent G spin gauge 内计算；
- character lookup 完全在 canonical H gauge 内计算；
- 两者只通过已验证的空间操作 setting transform 对应，不直接比较不同 gauge 的 SU(2)。

### 剩余失败的可能原因

以下是待验证假设，不是已经确认的根因。优先级按现有证据和可解释的失败数量排序。

| 优先级 | 假设 | 主要影响阶段 | 证据与判断 |
|--------|------|--------------|------------|
| P0 | MSG 中嵌入的 H 与 canonical Hall H 之间缺少基底变换 `T` | `square_not_in_spin_table`、`square_outside_little_group` | UNI663 已直接证实；低 SG3/SG1/SG10/SG5 等 mapping failure 高度符合轴选择错误 |
| P0 | canonical k 被直接用于 MSG 基底下的 rotation，导致 little-group 过滤错误 | `square_outside_little_group`、部分 `non_quantized` | 若 `x_can=T x_msg`，则 `k_msg=T^T k_can`；当前没有执行该变换。SG196/202/203/210/219 集中的 outside-LG 很可疑 |
| P0 | `seitz_preserves_k` 使用了错误的直接/倒空间 rotation convention | `square_outside_little_group`、`non_quantized` | 若 spglib 的 `R` 作用于直接晶格分数坐标，则倒空间应使用 `R^-T k`，不是 `R k`；三方/六方非正交整数矩阵最容易暴露。必须用数据 oracle 确认，不能凭公式直接改 |
| P0 | central parity 跨 G/H 两套 SU(2) gauge 比较 | `non_quantized` | 当前比较 parent G 中的 `U_b²` 和 canonical H 中的 `U_{b²}`；即使空间 rotation 对应，两套 lift gauge 也不保证可直接比较 |
| P1 | origin shift `s` 或 centering/nonsymmorphic representative 未完整求解 | `non_quantized`、少量 mapping failure | direct fallback 已加入 fractional phase，但使用的操作仍未经过完整 `(T,s)` 共轭；trigonal/hexagonal SG144/145/151/152/154/179/190 集中失败符合相位错误 |
| P1 | 只找到“群集合等价”的 `T`，但 operation correspondence 选错了群自同构 | 三类都可能 | 高对称群可能有多个 signed permutation 都把 rotation set 映射到同一集合；若选错，会置换 little-group operation 和 irrep character |
| P1 | parent G spin table 本身和 MSG parent setting 也不一致 | `AntiunitarySpinLookup` 或错误 central sign，最终表现为 `non_quantized` | UNI663 中 G 恰好对齐，不能据此推广到全部 1,651 个 UNI；需要全量 G-side lookup oracle |
| P2 | 所需基底变换不属于 48 个 signed permutation | signed-permutation oracle 无解的剩余 case | 单斜、三方、六方或不同 conventional setting 可能需要一般 unimodular/rational `T`，不能把首轮搜索范围当成最终模型 |
| P2 | direct anti-coset 的求和域或归一化仍有 Type-III 特例 | 修正 setting 后仍 `non_quantized` | 应满足反酉 little-group coset 与 unitary little group 等势；当前 k-frame 错误会掩盖这个 invariant，需在修正后重新检查 |
| P3 | 剩余 generated spin operation/character 对应关系仍有数据错误 | 修正全部坐标问题后的离散少量 case | 已修复 local/global index，但仍需用独立 closure、identity character、operation bijection 检查；当前没有证据支持先改数据 |
| P3 | 浮点 tolerance 导致误判 | 极少量边界 case | rotation 和多数相位来自小分母离散数据；失败大规模按 SG 聚集，不像纯数值问题。只能在逻辑问题修完后检查 |

不能采用的“修复”：

- 不得把非量子化 W 强行 round 到最近的 A/B/C；
- 不得仅凭 rotation set 相同就忽略 translation/origin；
- 不得为了降低失败数混用 per-term convention；
- 不得用一个未经全量验证的 `T` 覆盖所有同 SG 或同 crystal system；
- 不得在 20,710 个现有成功 case 出现 regression 时继续叠加补丁。

### 修正计划

#### Phase 0：冻结可比较的基线 ✅

基线已冻结（`diagnose_wigner_sources`，2026-06-19）：

```text
ok=20,710  (96.825%)
fail=679   (3.175%)
non_quantized=315
square_not_in_spin_table=194
square_outside_little_group=170
```

分母统一使用 21,389（全量 spinor irrep×UNI 对）。

#### Phase 1：setting-transform oracle → 验证 → 集成 ✅

**新增数据结构**（`wigner.rs`）：

- `SettingTransform { basis: T, origin: s }` — 基底+原点变换
- `transform_rotation()` / `transform_translation()` — 正向变换
- `enumerate_signed_permutations()` — 48 个 signed-permutation 矩阵
- `find_setting_transform(msg_rots, hall_rots)` — 自动搜索匹配的 T
- `rotation_multiset_eq()` — 顺序无关的 rotation multiset 比较

**约定**：`x_hall = T·x_msg + s`，`R_hall = T·R_msg·T⁻¹`

**Oracle 结果**（`phase1_setting_transform_oracle`，1,644 UNI）：

```text
identity:     1,356 (82.5%)  — MSG 和 H 在同一基底下
signed_perm:     72 (4.4%)   — 找到了轴置换
not_found:      216 (13.1%)  — 48 个 signed-permutation 都不匹配
```

**验证结果**（`phase1b_verify_transform_fix`）：
对 72 个 signed-perm UNI 的 136 个 `square_not_in_spin` 失败，
应用 T 变换后 **136/136 (100%)** 修复。第一个确认案例：
UNI663 SG3，T=swap(y,z)，C2z→C2y。

**正式路径集成**（`wigner_classify_spinor`）：
- 新增 `setting_xf: Option<&SettingTransform>` 参数
- 传递给 direct anti-coset fallback
- `diagnose_wigner_sources` 中从 `ops_from_msg`/`ops_from_hall` 自动计算 T

**Path triage 影响**：
- `MSG_GAUGE_W_FAIL`: 875 → 596 (**-279**)
- `OLD_PATH_FAIL`: 583 → 386 (**-197**)

总失败数 679 未变——那 136 个 case 已被 primary 路径修复，不依赖 fallback。

#### Phase 2：修正 reciprocal k 和 little-group 过滤

先确认 `R` 对 k 的作用约定。对每个 irrep，分别用：

```text
k' = R k
k' = R^-T k
```

重建 canonical H little co-group，并与生成数据中的 `spin_lg_op_indices` 对应的 rotation
集合比较。按 SG 统计两种 convention 的 exact-match 数量。这个 oracle 必须先回答：

- spglib/Hall rotation 是直接空间还是倒空间表示；
- spin.dat 的 k 和 operation rotation 使用哪一套 convention；
- unitary 与 antiunitary 分支是否分别需要 `R^-T k` 和 `-R^-T k`。

只有与数据全量一致的 convention 才能进入正式 `seitz_preserves_k`。禁止仅因某个单例
通过就全局替换。

优先把 MSG 操作变换到 canonical H 坐标，然后继续使用 irrep 表中的 canonical `k`。
这样 character lookup、spin LG 和 k 全部处于同一 frame。

若必须在 MSG frame 中过滤，则严格使用：

```text
k_msg = T^T k_hall
```

origin `s` 不改变 k，但会改变 Seitz translation 和 Bloch phase。

必须增加以下 invariant：

- 变换前后的 little-group 大小一致；
- 每个反酉 little-group 元素 `b` 都满足 `b² ∈ H_k`；
- antiunitary little-group coset 与 unitary little group 等势；
- canonical square mapping 必须落入 `spin_lg_op_indices`。

预期优先消除 194 个 `square_not_in_spin_table` 和 170 个
`square_outside_little_group`。若 mapping failure 不下降，说明 `(T,s)` 或 operation
correspondence 仍然错误，不能继续调 central sign。

#### Phase 3：拆分空间 character gauge 与 SU(2) central gauge

每个反酉操作保留两种表示：

- `b_msg`：原始 MSG/parent G frame，用于 G spin lookup 和 central parity；
- `b_hall`：经 `(T,s)` 变换后的 canonical H frame，用于 little-group 和 character lookup。

central parity 完全在 parent G spin gauge 内计算：

```text
U_b        <- G spin table
U_b_squared = U_b * U_b
U_sq_G     <- G spin table 中的 canonical lift of b²
spatial_central = relation(U_b_squared, U_sq_G)
central = !spatial_central
```

禁止再把 `U_b²` 与 H spin table 的 `U_{b²}` 直接比较。H spin table 只提供
character 对应的 operation index。

同时增加 G-side oracle：

```text
g_b_lookup_ok
g_square_lookup_ok
g_su2_same
g_su2_ebar
g_su2_unrelated
```

目标是 `g_su2_unrelated=0`。若不为零，说明 parent G setting/lift 也需独立变换，
不能靠 sign fallback 掩盖。

#### Phase 4：统一 translation 与 Bloch phase

所有 phase 必须从同一 canonical H frame 的 translation difference 得出：

```text
L = t_computed_hall - t_canonical_hall
phase = exp(+2πi k_hall·L)
```

`L` 允许包含 centering/nonsymmorphic fractional component，不能提前 round 成整数。
对每个 term 记录：

```text
b
b²
canonical H op
L
phase
chi
central sign
contribution
```

Phase 4 后重点观察 315 个 `non_quantized`。如果它们按 trigonal/hexagonal SG 大幅下降，
说明主因是 setting/origin/phase；若不下降，再进入下一阶段。

#### Phase 5：处理非 signed-permutation setting

对 `transform_not_found` 的 case：

1. 按 crystal system 和 Hall choice 分组；
2. 检查是否存在一般 `GL(3,Z)` unimodular basis；
3. 必要时使用有界小整数矩阵搜索，但必须用 rotation conjugacy + full Seitz bijection
   双重约束，不能无界 brute force；
4. 若 conventional cell/centering 不同，允许 rational `T`，并验证体积比和纯平移子群；
5. 优先复用磁群识别中已经验证过的 `x_std = T x + s` 变换逻辑，避免另写一套方向相反的公式。

#### Phase 6：对最终 residual 做逐 term 物理审计

只有坐标和 gauge invariant 全部通过后，才检查：

- operation-to-character 是否被群自同构错误置换；
- `spin_lg_op_indices` 是否包含正确的两个 central lifts；
- Type-III 求和域和分母是否正确；
- improper rotation 的 `R`/`-R` lookup 是否有歧义；
- character identity 值、共轭关系和群乘 closure；
- W 偏离量子值的模式是离散相位错误还是浮点残差。

若 W 只在 `1e-12` 量级偏离，再调整 tolerance；在此之前禁止把 tolerance 当主修复。

#### Phase 7：接入、回归和收尾

每个可编译修改都独立 commit，并依次执行：

```bash
cd /home/liuyichen/TB_rs
cargo check --package cryspglib
cargo test --package cryspglib debug_direct_anti_setting_uni663 -- --nocapture
cargo test --package cryspglib diagnose_wigner_sources -- --nocapture
```

阶段性验收：

1. 已有 20,710 个成功 case 不得 regression；
2. 每次提交记录三类失败数量变化；
3. 将 UNI663 debug test 转成无输出的 regression assertion；
4. `spinor_complex_fail=0` 后运行：

   ```bash
   cargo test --package cryspglib
   cargo test --package cryspglib --tests
   ```

5. 删除仅用于排查的噪声输出，保留可重复 oracle；
6. 更新本节最终统计和根因，不保留未经验证的推测作为结论。

### 完成标准

只有同时满足以下条件才能宣称 100%：

- `diagnose_wigner_sources`：`spinor_complex_ok=21,389`、`spinor_complex_fail=0`；
- 三类 failure stage 全部为 0；
- setting transform、little-group、G-side SU(2) invariant 无 unresolved case；
- 全量 unit/integration/doc tests 通过；
- 无现有成功 case 被另一条启发式路径改变分类。

关键诊断命令：

```bash
cd /home/liuyichen/TB_rs
cargo test --package cryspglib debug_direct_anti_setting_uni663 -- --nocapture
cargo test --package cryspglib diagnose_wigner_sources -- --nocapture
cargo test --package cryspglib
```
