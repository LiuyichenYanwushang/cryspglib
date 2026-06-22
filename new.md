# 工作总结与后续计划（供 Codex Review）

## 一、Baseline 演变

| 阶段 | ok | fail | 改动 |
|------|-----|------|------|
| 初始 (new.md 记录) | 20,921 | 301 | 诚实基线 (strict dispatcher) |
| `enumerate_unimodular_bases()` | 21,015 | 201 | SG13/14 Hall73→72/82→81 (-100) |
| 容差 1e-6→1e-4 | 21,079 | 137 | W≈0 浮点假阳性 (-64) |
| codex 审核后 4 bug 修复 | 21,127 | 89 | (计数不变，错误抵消) |
| 180° quaternion 提取修复 | 21,148 | 68 | SG122/199/206/220 (-21) |
| **当前** | **21,148** | **68** | **99.68%** |

## 二、所有修改的提交（自上轮 codex review 以来）

### unimodular bases 扩展（`src/irrep/wigner.rs`）

- `35be90f`：`enumerate_unimodular_bases()` — 生成所有 entries∈{-1,0,1}, det=±1 的 3×3 矩阵，替换 `enumerate_signed_permutations()` 用于 `find_setting_transform`。消除 SG13 (P2/c, Hall73→72, 64 cases) + SG14 (P2₁/c, Hall82→81, 36 cases)。
- `d67f0fe`：`OnceLock<Vec<Mat3I>>` 缓存 unimodular bases，避免每次调用重新生成 ~6,960 个矩阵。

### 容差与分类规则修正（`src/irrep/wigner.rs`）

- `3088701`：容差从 1e-6→1e-5。消除 64 个 W≈0 浮点假阳性（SG144/145/151/153/169/170/171/172，四方 screw 群，W≈0.000006 < 1e-5）。
- `0065378`：删除 W=±dim 接受分支，仅保留 W=0, ±1（arXiv:2211.10740）。

### 诊断基础设施（`src/irrep/wigner.rs`, `src/irrep/corep.rs`）

- `47ddfbc`：`PerTermTrace` struct — 记录每个反酉操作 b 的完整数据流（旋转、SU(2) lift、parity、Bloch phase、contribution）。
- `bbdbc09`：`diagnose_nonquantized_per_term` test — 按 (SG, k-point) 聚类展示 per-term Wigner sum，含完整 SU(2) 数据。

### codex 审核后的 4 个 bug 修复（`src/irrep/wigner.rs`）

- `52c9800`：新增 `LiftRelation` enum (`Same`/`EBar`) 替代 `su2_same_up_to_sign` 的模糊 bool 语义。
- `34a57b2`：
  - **Bug 1**: central 判断反了 — 显式 `-U²` 后 `Same` 表示无 Ē → `central = spatial_central == LiftRelation::EBar`。
  - **Bug 2**: `conj_pauli` 不等价于 SU(2) 逆 — 新增 `quat_conj()` = `[u0, -u1, -u2, -u3]` 替换 `conj_pauli()`。
  - **Bug 3**: parity oracle 收到 H frame 的 `sq.rot` 而非 G frame 的 `b_sq_rot_ms_g` — 修正传入参数。
  - **Bug 4**: `is_signed_perm` 只检查 entry∈{-1,0,1} — 增加 row/col 非零计数 + P·Pᵀ=I 验证。
- `bde55f8`：反酉平方显式写为 `neg_pauli(&su2_compose(&u_b, &u_b))` = -U²。G→H spin parity oracle：`compute_signed_perm_spin_parity()` 用 Q = det(P)·P 变换 G lift 到 H frame。
- `bfac77d`：扩展 parity oracle 支持 det(P)=1。
- `7ed502b`：修复 180° quaternion 提取 — 用 Q+I column 代替对角线扫描，支持绕非坐标轴的 signed-permutation。消除 SG122(8)+SG199(4)+SG206(4)+SG220(5) = 21 cases。添加 4 个 parity 不变性测试。
- `032a7af`：提取 `signed_perm_to_quat()` 为独立函数。

### 诊断分析提交

- `032a7af`：同表 SU(2) 一致性诊断 — `.position()` 与 `find_sq_spin_lg_first` 找到的 SU(2) lift 完全相同（全部 `Same`），排除查找不一致假说。

## 三、codex 指出的问题修复状态

| # | 问题 | 状态 | 说明 |
|---|------|------|------|
| 1 | central 判断反了 | ✅ `34a57b2` | `Same→central=false, EBar→central=true` |
| 2 | frame conjugation 用错函数 | ✅ `34a57b2` | `quat_conj` 替代 `conj_pauli` |
| 3 | parity oracle 收到错误坐标帧 | ✅ `34a57b2` | 传入 `b_sq_rot_ms_g` |
| 4 | `is_signed_perm` 判断不充分 | ✅ `34a57b2` | row/col 非零 + P·Pᵀ=I |
| 5 | 180° quaternion 提取只找对角线 +1 | ✅ `7ed502b` | 用 Q+I column |
| 6 | quaternion 向量符号回归验证 | ❌ | 未做 |
| 7 | `antiunitary_square_pauli()` 顺序错误 | ❌ | 仅 legacy path 受影响 |
| 8 | 注释 χ((a₀b)²) 实际循环计算 χ(b²) | ❌ | 注释问题 |
| 9 | `neg_rot` 未使用变量 | ❌ | 清理了但还有未使用警告 |

## 四、剩余 68 个 failure 的分布

全部为 **G=H (xf=I, parent=H)**：
| SG | 数量 | W 值 | 说明 |
|----|------|------|------|
| 88 | 8 | 0.707 | I4₁/a，体心四方 |
| 141 | 10 | 0.354/0.707 | I4₁/amd，体心四方 |
| 142 | 10 | 0.354/0.707 | I4₁/acd，体心四方 |
| 201 | 4 | 0.500 | Pn-3，原始立方 |
| 203 | 8 | 0.707 | Fd3，面心立方 |
| 222 | 2 | 0.500 | Pn3n，原始立方 |
| 224 | 6 | 0.5-1.5 | Pn3m，原始立方 |
| 227 | 10 | 0.354/0.707 | Fd3m，面心立方 |
| 228 | 10 | 0.354/0.707 | Fd3c，面心立方 |

## 五、关键代码位置（当前状态）

### 生产路径

| 功能 | 文件:行 |
|------|--------|
| `LiftRelation` enum | `wigner.rs:~2082` |
| `su2_lift_relation()` | `wigner.rs:~2092` |
| `signed_perm_to_quat()` | `wigner.rs:~234` |
| `compute_signed_perm_spin_parity()` | `wigner.rs:~267` |
| 显式 -U²: `neg_pauli(&su2_compose(&u_b, &u_b))` | `wigner.rs:~1842` |
| u_sq_g 查找 (G table `.position()`) | `wigner.rs:~1845` |
| central 判断: `spatial_central == EBar` | `wigner.rs:~1876` |
| G→H signed-perm parity oracle 块 | `wigner.rs:~1878-1905` |
| 容差 1e-5 | `wigner.rs:~1883` |
| W 分类 (仅 W=0, ±1) | `wigner.rs:~1883-1893` |
| `PerTermTrace` struct | `wigner.rs:~1517` |
| `enumerate_unimodular_bases()` | `wigner.rs:~1197` |

### 诊断测试

| 测试 | 文件:行 |
|------|--------|
| `diagnose_wigner_sources` | `corep.rs:2078` |
| `diagnose_nonquantized_per_term` | `corep.rs:3274` |
| 4 个 parity 不变性测试 | `wigner.rs:~3400` |
| 运行命令 | `cargo test -p cryspglib --release diagnose_wigner_sources -- --nocapture` |

## 六、本次诊断的核心结论

### 根因已排除的假说

1. ❌ **G spin table gauge 不一致**：codex 指出 (-U)²=U²，符号不影响平方。
2. ❌ **`.position()` 与 `find_sq_spin_lg_first` 找到不同 SU(2) sign**：debug 输出确认两者完全相同（全部 `Same`）。
3. ❌ **G→H parity (对剩余 cases)**：所有剩余失败都是 G=H (xf=I)，无 parity 问题。

### 定位到的真正问题

**SG88 per-term 数据显示**：同一 b² rotation (C2z = 180° around z) 的不同 b 得到不同的 central parity：
- b = +90°: u_b = (0.707,0,0,0.707), u_b² = +k, u_b_sq = -k, u_sq_g = -k → Same → central=false
- b = -90°: u_b = (0.707,0,0,-0.707), u_b² = -k, u_b_sq = +k, u_sq_g = -k → EBar → central=true

物理上两者 b² 相同（都是 C2z），Θ² 因子应给出相同的 central parity。但 `u_b` 来自 G spin table 的 rotation-only 查找（不含 Θ 因子），不同旋转角度的 canonical lift 符号约定不一致，导致 `u_b²` 符号不同。

### 问题本质

当前的 central parity 检测流程：
```
u_b = G_spin_table[rotation=b.rot]   ← 不含 Θ 因子
u_b² = su2_compose(u_b, u_b)
u_b_sq = neg_pauli(u_b²)             ← 显式补 Θ²=-1
u_sq_g = G_spin_table[rotation=b²]  ← canonical lift of b²
central = (u_b_sq 与 u_sq_g 异号)
```

`u_b` 是 G spin table 中 rotation R 的 canonical lift。对于 rotation R 和 rotation R'（R² = R'² = C2z），spin table 可能存不同的符号约定。虽然 (-U)²=U²，但 neg_pauli 之后正负反转，导致不同 b 的 central 不一致。

**更正确的做法可能是不依赖 u_b 来计算 central parity，而是直接从 b² rotation 的物理性质出发。** 例如，对于 b² = I（恒等元），Θ² = Ē 恒成立，central 应为 true；对于其他 rotation，需要分析 (a₀h)² 在双群中的结构。

## 七、后续计划

### 第一步：重新设计 central parity 检测

当前的 per-b comparison (u_b_sq vs u_sq_g) 对不同 b 给出不同答案。需要一种只依赖 b² rotation 的 parity 判断方法——所有平方到同一 rotation 的 b 应有相同的 central。

可能的方向：
1. 直接用 b² rotation 的结构判断（I → central=true，-I → 依赖具体情况，其他 rotation → 查表或解析计算）。
2. 对同一 b² 统一 central：若某个 b² 出现在多个 b 中，取其 majority vote 或统一符号。

### 第二步：实现通用 G→H cocycle（对 G≠H 的剩余 cases）

目前仅 SG203 (parent=227, xf=permutation) 有 G≠H。但 parity oracle 的 4 个 bug 修复后，需要重新验证 SG92/95 的"修复"是否物理正确。

### 第三步：清理遗留问题
- `antiunitary_square_pauli()` 顺序
- quaternion 向量符号回归验证
- 未使用变量

### 验收标准
- 总样本数不变
- mapping failure = 0（已达成）
- non_quantized → 0 或确认为物理上真正非量子化
- 反酉 LG 非空时必须量子化为 0, ±1
