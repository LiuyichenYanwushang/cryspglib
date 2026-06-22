# 工作总结与后续计划（供 Codex Review）

## 一、Baseline 演变

| 阶段 | ok | fail | 改动 |
|------|-----|------|------|
| new.md 初始基线 | 20,921 | 301 | 诚实基线 (strict dispatcher) |
| `enumerate_unimodular_bases()` | 21,015 | 201 | SG13/14 Hall73→72/82→81 (-100) |
| 容差 1e-6→1e-4 | 21,079 | 137 | W≈0 浮点假阳性 (-64) |
| **当前 (step 1-7 后)** | **21,127** | **89** | G→H spin parity oracle 修复 SG92/95 (-48) |

## 二、所有修改的提交

### 扩增 candidate pool (`src/irrep/wigner.rs`)
- `35be90f`：`enumerate_unimodular_bases()` — 生成所有 entries∈{-1,0,1}, det=±1 的 3×3 矩阵，替换 `enumerate_signed_permutations()` 用于 `find_setting_transform`
- `d67f0fe`：`OnceLock<Vec<Mat3I>>` 缓存 unimodular bases，避免每次调用重新生成 6,960 个矩阵

### LiftRelation enum + 显式 -U² (`src/irrep/wigner.rs`)
- `52c9800`：新增 `LiftRelation` enum (`Same`/`EBar`) 替代 `su2_same_up_to_sign` 的模糊 bool 语义。新增 `su2_lift_relation()` 为主 API，`su2_same_up_to_sign()` 保留为 legacy wrapper
- `bde55f8`：反酉平方显式写为 `neg_pauli(&su2_compose(&u_b, &u_b))` = -U²，而非依赖 spin table gauge 隐式检测 Θ² = -1

### G→H spin frame parity (`src/irrep/wigner.rs`)
- `bde55f8`：`compute_signed_perm_spin_parity()` — 对 signed-permutation setting transform，计算轴向向量变换 Q = det(P)·P 的 SU(2) 表示 U_Q，将 G spin table lift 变换到 H frame 后比较，输出 parity ε(R) = ±1
- `bfac77d`：扩展 parity oracle 支持 det(P)=1（原本仅 det=-1）
- parity 修正逻辑位置：`wigner_classify_spinor_direct_anti_diagnostic` 中 `central` 计算之后 (`:1737-1809`)

### 容差与分类规则 (`src/irrep/wigner.rs`)
- `0065378`：容差从 1e-4 收紧至 1e-5
- `0065378`：删除 W=±dim 接受分支，仅保留 W=0, ±1（arXiv:2211.10740）

### Per-term trace 诊断 (`src/irrep/corep.rs`)
- `47ddfbc`：`PerTermTrace` struct — 记录每个 b 的完整数据流（旋转、SU(2) lift、parity、Bloch phase、contribution）
- `bbdbc09`：`diagnose_nonquantized_per_term` test — 按 (SG, k-point) 聚类展示 per-term Wigner sum，含 SU(2) 数据

## 三、codex review 后的修复路线与执行状态

codex review 指出根因不是"G spin table gauge 不一致"（(-U)²=U²，符号不影响平方），而是**缺少 G→H spin frame parity**。

### 建议执行顺序与完成状态

| # | 任务 | 状态 | 说明 |
|---|------|------|------|
| 1 | 撤回 4460be0，改用 LiftRelation enum | ✅ | `52c9800` |
| 2 | 反酉平方显式写成 -U² | ✅ | `bde55f8` |
| 3 | signed-permutation 轴向向量 quaternion oracle，验证 SG92/95 | ✅ | `bde55f8` — **SG92(24)+SG95(24)=48 eliminated** |
| 4 | 通用 G→H Z₂ cocycle/1-cochain parity（非 signed-permutation） | ❌ | 剩余 89 个失败需要此步 |
| 5 | 容差收紧为 1e-5 | ✅ | `0065378` |
| 6 | 缓存 unimodular 候选和 Hall-pair 结果 | ✅ | `d67f0fe` — unimodular 已缓存 |
| 7 | 删除 ±dim 接受规则 | ✅ | `0065378` |
| 8 | 重新统计剩余失败 | ✅ | 89 个 |

## 四、SG92 具体修复验证

UNI808: parent G = SG96, H = SG92, MSG→H basis P = diag(1,-1,1), det(P) = -1

轴向向量变换：Q = det(P)·P = diag(-1,1,-1)（180° 绕 y 轴）

对 C2z (b²=180°绕z)：U_G(C2z) = -k → Q 变换后 → +k。H canonical lift = -k。ε(C2z) = -1。

应用 ε = -1 后，8 个 term 全部贡献 -1 → W = -8/8 = -1 → Type B ✓

SG95 (P4₃2₂) 同理：P 也是 diag(1,-1,1) 或等价 signed-permutation。

## 五、当前剩余 89 个 failure 的分布

全部在体心/面心/高对称 cubic 群，具有非 signed-permutation 的 setting transform（含 1/2, 1/4 等剪切分量）：

| SG | 数量 | 说明 |
|----|------|------|
| 88 | 8 | I4₁/a，体心四方 |
| 122 | 8 | I-42d，体心四方 |
| 141 | 10 | I4₁/amd，体心四方 |
| 142 | 10 | I4₁/acd，体心四方 |
| 199 | 4 | I2₁3，体心立方 |
| 201 | 4 | Pn-3，原始立方 |
| 203 | 8 | Fd3，面心立方 |
| 206 | 4 | Ia3，体心立方 |
| 220 | 5 | I-43d，体心立方 |
| 222 | 2 | Pn3n，原始立方 |
| 224 | 6 | Pn3m，原始立方 |
| 227 | 10 | Fd3m，面心立方 |
| 228 | 10 | Fd3c，面心立方 |

## 六、关键代码位置

### 生产路径

| 功能 | 文件:行 |
|------|--------|
| `LiftRelation` enum | `wigner.rs:2076` |
| `su2_lift_relation()` | `wigner.rs:2086` |
| `su2_same_up_to_sign()` (legacy) | `wigner.rs:2106` |
| `compute_signed_perm_spin_parity()` | `wigner.rs:226` |
| G→H parity 应用 | `wigner.rs:1805-1840` |
| 显式 -U² | `wigner.rs:1774` |
| 容差 1e-5 | `wigner.rs:1864` |
| `enumerate_unimodular_bases()` | `wigner.rs:1197` |
| `UNIMODULAR_BASES` static | `wigner.rs:1194` |
| `PerTermTrace` struct | `wigner.rs:1517` |
| `wigner_classify_spinor_direct_anti_diagnostic` | `wigner.rs:1634` |

### 诊断测试

| 测试 | 文件:行 |
|------|--------|
| `diagnose_wigner_sources` | `corep.rs:2078` |
| `diagnose_nonquantized_per_term` | `corep.rs:3274` |
| 运行命令 | `cargo test --package cryspglib --release diagnose_wigner_sources -- --nocapture` |

## 七、尚未完成的工作

1. **通用 G→H Z₂ cocycle/1-cochain parity**（codex step 4）：
   - 对非 signed-permutation 的 rational basis，无法直接用 U_Q 变换 quaternion
   - codex 建议：分别构造 G restricted-to-H 和 H spin table 的 Z₂ multiplication cocycle，求解连接两套 section 的 1-cochain ε(h)，用 ε(b²) 修正 central parity
   - 位置：`wigner.rs:1805-1840`（当前仅处理 signed-permutation）

2. **antiunitary_square_pauli() 顺序错误**（codex 指出）：
   - 当前使用 `J U`（位置 `wigner.rs:234`），应为 `J U*`
   - 影响：legacy path 可能受影响（但当前仅 direct path 用于生产）

3. **Hall-pair transform 结果缓存**（codex 建议）：
   - `find_setting_transform` 的 Hall-pair 结果可缓存以减少重复搜索

4. **其他同级 `let tol = 1e-6` 实例**（`wigner.rs:808,884,2358,2780`）：
   - 这些在 legacy path 中，暂不影响生产路径但需统一

5. **`neg_rot` 未使用变量**（`wigner.rs:261`）：
   - 在 `compute_signed_perm_spin_parity` 中声明但未使用

## 八、后续计划

### 第一步：通用 G→H parity（预计影响 89 个剩余失败中的大部分）

对非 signed-permutation rational basis P：
1. 构造 G restricted-to-H spin table 的 Z₂ multiplication cocycle μ_G(h₁,h₂)
2. 构造 H spin table 的 Z₂ multiplication cocycle μ_H(h₁,h₂)
3. 求解 1-cochain ε: G→{±1} 满足 μ_H = μ_G · δε
4. 对每个 b² term：ε(b²) 修正 central parity

### 第二步：清理遗留问题
- 修复 `antiunitary_square_pauli` 顺序
- 统一所有 `tol = 1e-6` 实例
- 移除 `neg_rot` 未使用变量

### 验收标准
- 总样本数不变
- mapping failure = 0（已达成）
- non_quantized → 0 或确认为物理上真正非量子化
- 反酉 LG 非空时必须量子化为 0, ±1
