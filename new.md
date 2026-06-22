# 工作总结与后续计划（供 Codex Review）

## 一、Baseline 演变

| 阶段 | ok | fail | 改动 |
|------|-----|------|------|
| 初始 (strict dispatcher) | 20,921 | 301 | 诚实基线 |
| `enumerate_unimodular_bases()` | 21,015 | 201 | SG13/14 Hall73→72/82→81 (-100) |
| 容差 1e-6→1e-4 | 21,079 | 137 | W≈0 浮点假阳性 (-64) |
| **当前 (codex 审核修正后)** | **21,127** | **89** | SG92/95 仍在 ok 中，但 parity 实现有 bug 待验证 |

### codex 审查结论

SG92/95 的 48 个"修复"目前判定为**多个错误互相抵消**，不足以证明物理正确性。在对 parity oracle 的错误进行修正后，baseline 未变（89），但 SG92/95 还在 ok 中的原因待独立验证。

## 二、codex 指出的 4 个关键 bug 及修复

### Bug 1: central 判断反了 (`wigner.rs:1874`)
- **问题**: 显式 `-U²` 后 `Same` 表示无 Ē → `central=false`，但代码写成了 `spatial_central == Same`
- **修复**: `let mut central = spatial_central == LiftRelation::EBar` (34a57b2)

### Bug 2: frame conjugation 用错函数 (`wigner.rs:333`)
- **问题**: `U_Q · U_G · conj_pauli(U_Q)` 中 `conj_pauli` 给出 `[u0,-u1,u2,-u3]`（Pauli 复共轭），而非四元数逆 `[u0,-u1,-u2,-u3]`
- **修复**: 新增 `quat_conj()` 替代 `conj_pauli()` (34a57b2)

### Bug 3: parity oracle 收到错误坐标帧的 rotation (`wigner.rs:1917`)
- **问题**: `compute_signed_perm_spin_parity` 声明需要 G frame 的 `sq_rot: b² rotation in G/MSG frame`，但传入的是已变换到 H frame 的 `sq.rot`
- **修复**: 传入 G frame 的 `b_sq_rot_ms_g` (34a57b2)

### Bug 4: `is_signed_perm` 判断不充分 (`wigner.rs:1892`)
- **问题**: 仅检查 entry ∈ {-1,0,1}，剪切矩阵也会通过
- **修复**: 增加检验每行/列恰好一个非零、P·Pᵀ=I (34a57b2)

## 三、codex 指出的其他问题

| # | 问题 | 位置 | 状态 |
|---|------|------|------|
| 5 | 180° quaternion 提取只找对角线+1 | `wigner.rs:307` | 未修 — 无法处理绕 (1,1,0) 等轴的 signed perm |
| 6 | quaternion 向量符号需与 spin table Pauli convention 回归验证 | — | 未做 |
| 7 | `antiunitary_square_pauli()` 顺序错误 (J U → J U*) | `wigner.rs:234` | 未修 — 仅 legacy path 受影响 |
| 8 | `wigner.rs:1651` 注释写 χ((a₀b)²) 实际循环计算 χ(b²) | — | 未修 |
| 9 | new.md 中很多行号已过期 | — | 本次更新 |
| 10 | `neg_rot` 未使用变量 | `wigner.rs:261` | 未修 |

### 正当的修改（codex 确认）

- `OnceLock` 缓存 unimodular bases
- `LiftRelation` enum
- `1e-5` 容差
- 删除 `±dim` 接受规则

## 四、当前剩余 89 个 failure 的分布

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

**重要**: 之前 `new.md:66` 声称"剩余全部是非 signed-permutation"，但由于 bug 4（`is_signed_perm` 判断不充分），这个结论尚未得到可靠验证。修正 `is_signed_perm` 后需要重新确认哪些 group 的 transform 真正是 signed-permutation。

## 五、关键代码位置（修正后）

### 生产路径

| 功能 | 文件:行 |
|------|--------|
| `LiftRelation` enum | `wigner.rs:2082` |
| `su2_lift_relation()` | `wigner.rs:2092` |
| `su2_same_up_to_sign()` (legacy) | `wigner.rs:2111` |
| `quat_conj()` | `wigner.rs:230` |
| `compute_signed_perm_spin_parity()` | `wigner.rs:240` |
| G→H parity 应用 | `wigner.rs:1886-1940` |
| 显式 -U²: `neg_pauli(&su2_compose(&u_b, &u_b))` | `wigner.rs:1851` |
| central 判断: `spatial_central == EBar` | `wigner.rs:1881` |
| 容差 1e-5 | `wigner.rs:1879` |
| `enumerate_unimodular_bases()` | `wigner.rs:1203` |
| `UNIMODULAR_BASES` cache | `wigner.rs:1200` |
| `is_signed_perm` 验证 (修正后) | `wigner.rs:1892-1913` |
| `PerTermTrace` struct | `wigner.rs:1527` |
| `wigner_classify_spinor_direct_anti_diagnostic` | `wigner.rs:1640` |

### 诊断测试

| 测试 | 文件:行 |
|------|--------|
| `diagnose_wigner_sources` | `corep.rs:2078` |
| `diagnose_nonquantized_per_term` | `corep.rs:3274` |
| 运行命令 | `cargo test --package cryspglib --release diagnose_wigner_sources -- --nocapture` |

## 六、codex 建议的后续计划

### 第一步：添加硬性测试验证 signed-permutation parity

- ε(I) = +1
- U_Q · U_I · U_Q⁻¹ = U_I (identity invariant)
- G rotation 变换后必须等于 H rotation
- 所有真正 signed-permutation 的 multiplication parity 满足群乘法

### 第二步：修正 signed-permutation parity 实现

- 修复 180° quaternion 提取（支持绕非坐标轴的 signed perm）
- 回归验证 quaternion 向量符号与 spin table Pauli convention
- 统一处理 antiunitary_square_pauli 的 J·U → J·U* 或统一用 -U²

### 第三步：重新统计 SG92/95

预计 baseline 可能回退（这些 case 目前靠错误抵消通过），这是必要的诚实结果。

### 第四步：通用 Z₂ cocycle/1-cochain parity

- 对非 signed-permutation rational basis，构造 G restricted-to-H 和 H spin table 的 Z₂ multiplication cocycle
- 求解连接两套 section 的 1-cochain ε(h)
- **重要**: ε 不唯一，不同解相差 H→Z₂ 群同态。需要 generator anchor、固定 lift convention 或 Cartesian spin-frame 关系保持确定性

### 验收标准
- 总样本数不变
- mapping failure = 0（已达成）
- non_quantized → 0 或确认为物理上真正非量子化
- 反酉 LG 非空时必须量子化为 0, ±1
