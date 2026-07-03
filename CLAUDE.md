# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Irrep 终极目标与详细实施计划

最终目标：给定一个磁空间群（结构识别得到的 UNI、直接输入 UNI、或 BNS/OG 标签），程序应能一次性回答：

1. 这个磁群有哪些高对称 k 点/线/面，以及每个标签的标准名称和坐标。
2. 这个磁群在这些 k 点上的磁共表示（corepresentation）：来源 H-irrep、corep type A/B/C、维数、字符、antiunitary 完整性、Hall/setting 约定。
3. 每个磁共表示对应的可能 isotropy subgroup：普通 isotropy subgroup 和 magnetic isotropy subgroup 都要能追溯到 source irrep/corep、k 点、方向、domain/arm 信息。

### 当前完成度（2026-07-03）

已完成或基本可用：

- `src/irrep/query.rs`: `irreps_of(sg)`, `kpoints_of(sg)`, `IrrepRecord::k_label()` 已能列出 230 个 SG 的 ISOTROPY 高对称 k 标签、坐标和 irreps。
- `src/irrep/types.rs` + `generated_data.rs`: `IrrepRecord::subgroups()` 和 `IrrepRecord::magnetic_subgroups()` 已保存普通/磁 isotropy subgroup 数据。
- `src/irrep/corep.rs`: `compute_corepresentation()` / `compute_coreps(bns, k_label)` 已实现 scalar PIR、scalar CIR、spinor SU(2) Wigner 分类。
- 最新 spinor Wigner 全扫诊断已清零失败：`cargo test --package cryspglib diagnose_wigner_sources -- --nocapture` 中 `spinor_complex_ok = 21216`，无 `spinor_complex_fail`。
- `magnetic_isotropy_coreps_of_irrep()` / `magnetic_isotropy_coreps_of_sg_k()` 已有 corep 与 magnetic isotropy 的早期桥接，但还不是最终用户 API。

未完成的关键缺口：

- 缺少一个面向用户的统一入口：当前能力分散在 `query.rs`, `corep.rs`, `api.rs`, `magnetic_spacegroup.rs`。
- 磁群 k 点语义需要固定：第一版采用 unitary subgroup H 的 ISOTROPY k 点作为标准列表，同时显式返回 magnetic little group 的 unitary/antiunitary 阶数；后续再处理 antiunitary 合并 k-star 的展示策略。
- Type-C corep 需要 pairing/dedup：不能把一对互为 antiunitary 共轭的 H irreps 重复展示成两个磁 coreps。
- corep 到 isotropy subgroup 的映射要明确为“候选”：第一版把 source H-irrep 的 ordinary/magnetic isotropy subgroup 合并挂到 corep；后续再验证 Type-C/compound/spinor 的物理筛选规则。

### 目标 API（新增）

新增文件：`src/irrep/magnetic_summary.rs`。

对外入口：

```rust
pub fn magnetic_irrep_summary(input: MagneticIrrepInput)
    -> Result<MagneticIrrepSummary, MagneticIrrepError>;

pub fn magnetic_irrep_summary_by_uni(uni: usize)
    -> Result<MagneticIrrepSummary, MagneticIrrepError>;

pub fn magnetic_irrep_summary_by_bns(bns: &str)
    -> Result<MagneticIrrepSummary, MagneticIrrepError>;

pub fn magnetic_irrep_summary_from_ops(
    uni: usize,
    mag_ops: &crate::SymmetryOps,
) -> Result<MagneticIrrepSummary, MagneticIrrepError>;
```

输入类型：

```rust
pub enum MagneticIrrepInput<'a> {
    Uni(usize),
    Bns(&'a str),
    Operations { uni: usize, ops: &'a crate::SymmetryOps },
}
```

核心返回结构：

```rust
pub struct MagneticIrrepSummary {
    pub uni: usize,
    pub bns_label: String,
    pub magnetic_type: crate::MagneticType,
    pub parent_sg: u8,
    pub unitary_sg: u8,
    pub unitary_hall: usize,
    pub kpoints: Vec<MagneticKPointSummary>,
}

pub struct MagneticKPointSummary {
    pub label: String,
    pub coords: (i8, i8, i8, i8),
    pub little_group_order: usize,
    pub unitary_order: usize,
    pub antiunitary_order: usize,
    pub coreps: Vec<MagneticCorepSummary>,
}

pub struct MagneticCorepSummary {
    pub label: String,
    pub source_irreps: Vec<SourceIrrepSummary>,
    pub corep_type: crate::irrep::corep::CorepType,
    pub source: crate::irrep::corep::WignerSource,
    pub dim: usize,
    pub characters: Vec<f64>,
    pub timerev: Vec<bool>,
    pub completeness: crate::irrep::corep::CharacterCompleteness,
    pub isotropy_candidates: Vec<CorepIsotropyCandidate>,
}

pub struct SourceIrrepSummary {
    pub sg: u8,
    pub ml: &'static str,
    pub bc: &'static str,
    pub dim: u8,
    pub spinor: bool,
}

pub struct CorepIsotropyCandidate {
    pub source_ml: &'static str,
    pub ordinary: Vec<crate::irrep::types::IsotropyRecord>,
    pub magnetic: Vec<crate::irrep::types::MagneticIsotropyRecord>,
    pub relation: IsotropyCandidateRelation,
}

pub enum IsotropyCandidateRelation {
    DirectSourceIrrep,
    TypeCPairedSource,
    CompoundSource,
    SpinorNoIsotropyData,
}
```

错误类型：

```rust
pub enum MagneticIrrepError {
    InvalidUni(usize),
    UnknownBns(String),
    MissingMagneticOperations(usize),
    MissingUnitarySubgroup(usize),
    MissingIrrepData { sg: u8 },
    CorepComputationFailed { uni: usize, sg: u8, k_label: String },
}
```

Re-export 路径：

- `src/irrep/mod.rs`: `pub mod magnetic_summary;`
- `src/irrep/mod.rs`: re-export 常用类型，或要求用户显式 `use cryspglib::irrep::magnetic_summary::*;`
- `src/lib.rs`: 暂不顶层 re-export，等 API 稳定后再决定是否暴露到 crate root。

### 实现顺序与路径

#### Phase 1: 只读 summary API 骨架

路径：`src/irrep/magnetic_summary.rs`, `src/irrep/mod.rs`。

实现方法：

1. 新增上面的 public structs/enums。
2. `magnetic_irrep_summary_by_uni(uni)`:
   - 校验 `1 <= uni <= 1651`。
   - `let mag_ops = SymmetryOps::from_magnetic_database(uni)?`。
   - 调用 `magnetic_irrep_summary_from_ops(uni, &mag_ops)`。
3. `magnetic_irrep_summary_by_bns(bns)`:
   - 复用 `corep.rs` 现有 BNS→UNI helper；若 helper 当前私有，移动/改为 `pub(crate)`。
   - 转入 `magnetic_irrep_summary_by_uni(uni)`。
4. `magnetic_irrep_summary_from_ops(uni, mag_ops)`:
   - 调用 `identify_unitary_subgroup_with_hall(uni)` 取得 H 信息。
   - 记录 `unitary_sg`, `unitary_hall`, `msg_to_data`。
   - 遍历 `query::kpoints_of(unitary_sg)` 生成 `MagneticKPointSummary`。

验收测试：

- `magnetic_summary_by_uni_smoke`: UNI 1599 或 BNS `221.97` 返回非空 kpoints。
- `magnetic_summary_by_bns_matches_uni`: BNS 和 UNI 两条路径返回相同 `uni/unitary_sg/kpoints.len()`。

#### Phase 2: magnetic little group 元数据

路径：`src/irrep/magnetic_summary.rs`，复用 `src/irrep/wigner.rs`。

实现方法：

1. 提取 helper:

```rust
fn canonical_pure_translations(h_ops: &crate::SymmetryOps) -> Vec<[f64; 3]>;

fn magnetic_little_group_indices(
    k: (i8, i8, i8, i8),
    mag_ops: &crate::SymmetryOps,
    setting_xf: Option<&crate::irrep::wigner::SettingTransform>,
    canonical_translations: &[[f64; 3]],
) -> Vec<usize>;
```

2. 对每个 k 点调用 `filter_little_group_with_transform`。
3. 填充 `little_group_order`, `unitary_order`, `antiunitary_order`。
4. 第一版 k 点列表保持 H 的 ISOTROPY k 点，不试图减少/合并 magnetic star；API 文档明确这一点。

验收测试：

- grey group: 每个有 antiunitary 的 k 点 `antiunitary_order > 0`。
- nonmag/type-I equivalent: `antiunitary_order == 0` 时 corep type 应走 trivial A。

#### Phase 3: corep 计算接入 summary

路径：`src/irrep/magnetic_summary.rs`, 必要时调整 `src/irrep/corep.rs` helper 可见性。

实现方法：

1. 对每个 `KPointSummary.irreps` 找到 H 的 `IrrepRecord`。
2. 调用：

```rust
corep::compute_corepresentation(ir, uni, mag_ops)
```

3. 把 `Corepresentation` 转成 `MagneticCorepSummary`。
4. `label` 第一版使用 source ML label；Type-C dedup 完成后改为组合 label。
5. 保留 `CharacterCompleteness`，不隐藏 pending/unsupported。

验收测试：

- `128.406` at `Z`: 至少返回 Type-C 和 Type-A coreps，维数符合现有 corep tests。
- `221.97` at `GM`: 返回非空 coreps，identity character 等于 dim。

#### Phase 4: Type-C pairing/dedup

路径：`src/irrep/magnetic_summary.rs`，必要时在 `src/irrep/corep.rs` 增加 `pub(crate)` helper。

新增内部类型：

```rust
struct CorepDedupKey {
    corep_type: CorepType,
    dim: usize,
    rounded_characters: Vec<i64>,
    timerev: Vec<bool>,
}
```

实现方法：

1. 第一版 dedup 使用 `corep_type + dim + rounded characters + timerev`，字符按 `1e-8` 量化。
2. Type-C 合并时 `source_irreps` 追加两个 H-irrep。
3. `label` 规则：
   - 单 source: `"GM4-"`。
   - Type-C pair: `"Z1Z4 + Z2Z3"` 或按 ML 排序 join。
   - compound source 保留原 compound ML label。
4. 后续增强：用 antiunitary conjugation 显式找 partner，而不是只靠 character key。

验收测试：

- `test_type_c_coreps_are_deduplicated` 的语义迁移到 summary API：Type-C pair 不重复。
- 每个 `MagneticCorepSummary.source_irreps` 非空；Type-C 至少能出现两个 source 的 case。

#### Phase 5: isotropy candidates 挂接

路径：`src/irrep/magnetic_summary.rs`。

实现方法：

1. 对每个 `source_irrep` 收集：
   - `ir.subgroups()` → ordinary candidates。
   - `ir.magnetic_subgroups()` → magnetic candidates。
2. 普通 subgroup 去重 key:

```rust
(sg, symbol, direction, domains, arms)
```

3. 磁 subgroup 去重 key:

```rust
(mag_sg, bns_label, direction)
```

4. `relation` 规则：
   - 单 source scalar: `DirectSourceIrrep`
   - Type-C 合并 source: `TypeCPairedSource`
   - `cir_component_count() > 0`: `CompoundSource`
   - spinor 且没有 isotropy 数据: `SpinorNoIsotropyData`
5. 第一版只声明 candidates，不声明这些 subgroup 已按 corep order parameter 方向完成物理筛选。

验收测试：

- SG 221 GM4- 路径能返回包含 ordinary/magnetic subgroup 的 candidates。
- Type-C 合并后 candidates 去重稳定，不因 source 顺序变化而重复。

#### Phase 6: 格式化与示例

路径：`src/irrep/magnetic_summary.rs`, `README.md` 或 `examples/`。

新增格式化 API：

```rust
pub fn format_magnetic_irrep_summary(summary: &MagneticIrrepSummary) -> String;
pub fn format_magnetic_kpoint_summary(kpoint: &MagneticKPointSummary) -> String;
```

示例文件：

- `examples/magnetic_irrep_summary.rs`

示例目标：

```rust
let summary = magnetic_irrep_summary_by_bns("221.97")?;
for kp in &summary.kpoints {
    println!("{} {:?}", kp.label, kp.coords);
    for c in &kp.coreps {
        println!("  {} {:?} dim={}", c.label, c.corep_type, c.dim);
    }
}
```

验收测试：

- example 能编译运行。
- 格式化输出不依赖 HashMap 随机顺序。

### 推荐提交顺序

1. `feat: add magnetic irrep summary types`
2. `feat: summarize magnetic k-points by UNI`
3. `feat: attach corepresentations to magnetic summary`
4. `feat: deduplicate type-c magnetic coreps`
5. `feat: attach isotropy candidates to magnetic coreps`
6. `docs: add magnetic irrep summary example`

每个提交前运行：

```bash
cd /home/liuyichen/TB_rs
cargo check --package cryspglib
```

关键节点额外运行：

```bash
cargo test --package cryspglib diagnose_wigner_sources -- --nocapture
cargo test --package cryspglib test_type_c_coreps_are_deduplicated -- --nocapture
cargo test --package cryspglib --tests
```

### 非目标（先不要做）

- 不在第一版里重新生成 ISOTROPY 数据。
- 不在第一版里发明新的磁 k 点命名系统；先复用 H 的 k labels，并暴露 magnetic little group 元数据。
- 不在第一版里声称 isotropy candidates 已完成 order-parameter 方向的唯一筛选；先返回可追溯候选。
- 不把 summary API 暴露到 crate root；先稳定在 `cryspglib::irrep::magnetic_summary`。

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

### 当前主要问题
Spinor Wigner 的 `square_not_in_spin_table` 主线已经在 2026-07-03 清零。
当前 irrep 方向的主要问题是产品化入口：把 H-irrep/k-point/corep/isotropy
这些已经分散可用的数据组织成稳定的 `magnetic_summary` API。

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

## Test suite (~205 tests across 8 binaries)

Core irrep diagnostics pass as of 2026-07-03.  The full spinor Wigner sweep reports
`spinor_complex_ok = 21216` with no `spinor_complex_fail`.

| Binary / Location | Tests | Description |
|-------|-------|-------------|
| `src/irrep/corep.rs` (unit) | 132 | Wigner diagnostics, BCS validation, CIR-PIR cross-validation, setting transform oracles, H2S triage |
| `src/irrep/query.rs` (unit) | 5 | API examples, character table formatting |
| `src/{api,lib,arithmetic,cell,debug,delaunay,determination}.rs` (unit) | ~26 | Entry-point API, arithmetic crystals, overlap detection, lattice reduction, error handling |
| `tests/irrep_validation.rs` (integration) | 31 | Full-sweep validation: every SG has irreps, dimensions match, labels well-formed, k-vectors positive |
| `tests/magnetic_integration.rs` (integration) | 11 | Magnetic structure analysis end-to-end (graphene, bilayer, Fe, CoF3, etc.) |
| `tests/{cof3,crps4,la2nio4,bcs_corep_validation}.rs` (integration) | 7 | Reference material cases |
| **Total** | **~205** | |

Key diagnostic tests (most useful for Wigner debugging):

```bash
# Primary: full-sweep Wigner failure diagnosis (~217s)
cargo test --package cryspglib diagnose_wigner_sources -- --nocapture

# Setting transform oracle (identity/signed_perm/unimodular/none x ok/other/square_not_in_spin)
cargo test --package cryspglib phase1_setting_transform_oracle -- --nocapture

# Per-term and per-case failure analysis
cargo test --package cryspglib diagnose_spinor_wigner_per_term -- --nocapture
cargo test --package cryspglib diagnose_nonquantized_per_term -- --nocapture
cargo test --package cryspglib diagnose_none_examples -- --nocapture

# CIR-PIR consistency and data integrity
cargo test --package cryspglib test_cir_pir_cross_validation -- --nocapture

# Magnetic entry point diagnostics
cargo test --package cryspglib diagnose_magnetic_entry_hall_anomalies -- --nocapture
cargo test --package cryspglib diagnose_spglib_standard_setting_transform -- --nocapture

# k-convention oracle (3611 irreps)
cargo test --package cryspglib diagnose_spin_lg_k_convention -- --nocapture
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
