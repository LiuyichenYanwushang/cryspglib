# 工作总结与后续计划（供 Codex Review）

## 一、Baseline 变化

| 指标 | 本次开始 | 当前 | 变化 |
|------|---------|------|------|
| ok | 20,611 | **20,921** | +310 |
| fail | 147 | **301** | +154 (诚实基线) |
| `square_not_in_spin_table` | 26 | **0** | 消除 |
| `square_outside_little_group` | 9 | **0** | 消除 |
| `non_quantized` | 112 | **301** | +189 (诚实基线暴露) |
| 全部测试 | 199 | **200** | +1 |

fail 增加是因为 strict dispatcher 暴露了 198 个之前被 legacy MSG-gauge path 掩盖的 case。

## 二、所有修改的提交（15 个）

### spglib Hall 搜索修复（`src/spacegroup.rs`）
- `028fd97`：`spa_search_spacegroup_with_symmetry` 搜索全部 530 Hall（之前只搜 230 个首选 Hall）
- `a3589e0`：拒绝 tolerance retry 产生的 Hall1 假阳性

### SettingTransform 防御（`src/irrep/wigner.rs`）
- `ee403dd`：`transform_rotation()`/`transform_translation()` → `Option` 返回；新增 `transform_seitz()`
- `1cf4748`：`filter_little_group_with_transform` 使用 transform 前原子验证全部 ops
- `3d24083`：`SettingTransform::then()` 组合方法 (P_total = P_2·P_1, s_total = P_2·s_1 + s_2)
- `1cb3706`：Hall→data-Hall 对齐
- `4fec731`：零原点 Seitz 验证 fallback（修复 Hall 210→209, 213→212 的 0 candidates）
- `43fc047`：`standard_setting_transform` fallback

### 帧统一（`src/irrep/corep.rs`）
- `fef7573`+`4aafedf`：canonical H 平移子群（先错用 `ops_from_msg` → 修正为 `ops_from_hall`）
- `8424b57`：`msg_to_data` 在 `identify_unitary_subgroup_with_hall` 中一次计算，`compute_corepresentation` 直接使用
- `9632324`：**Codex frame pipeline 统一**（核心提交，详见第三节）

### 路径修复
- `224cdc0`：b²∉H₀ skip → **被 `9632324` 撤回**
- `6c84bc5`：移除 G→H SU(2) cross-gauge fallback
- `175d775`：direct anti-coset 路径优先于 msg-gauge primary
- `23aa76a`：**strict dispatcher** — 仅 `MissingSpinData` 允许 legacy fallback，建立诚实基线

## 三、Codex Frame Pipeline 统一（`9632324`）引入的关键模式

### 1. `transform_embeds_ops()` — 允许 subset 嵌入
```rust
fn transform_embeds_ops(xf, source, target) -> bool
```
验证 source ops 变换后是 target ops 的**子集**（不要求逐项等长）。
中心化群 MSG 只有 8 个代表元，canonical Hall 有 16 个。
位置：`corep.rs:532`

### 2. `transform_applies_to_all_ops()` — 原子验证
所有 op 的 rotation 变换必须产生整数矩阵。
位置：`corep.rs:554`

### 3. 保留 `msg_to_detected`，不再丢弃 `_xf`
```rust
if let Some((std_sg, std_hall, xf)) = standard_setting_transform(...) {
    (std_sg, std_hall, oh, Some(xf))  // xf saved, not discarded
}
```
位置：`corep.rs:613`

### 4. 正确组合顺序
```rust
let msg_to_data = msg_to_detected.then(detected_to_data);
```
位置：`corep.rs:668`

### 5. 候选筛选用 `find()` 而非 `first()`
```rust
for detected_to_data in &xfs {
    if !transform_embeds_ops(detected_to_data, &ops_from_hall, &data_ops) { continue; }
    // ...
}
```
位置：`corep.rs:662`

### 6. 诊断与生产路径统一
诊断不再独立计算 `setting_xf`，直接使用 `h_info.msg_to_data` 和 canonical translations。
位置：`corep.rs:2118`

### 7. 撤回 b²∉H₀ skip
恢复 `SquareOutsideLittleGroup` 硬错误。frame 已统一，invariant 不再被违反。
位置：`wigner.rs:1580`

### 8. 新增回归测试
`test_centered_type3_msg_to_data_transform_is_validated` — 覆盖 SG22/42/43/69/70 的中心化 Type III 群。
位置：`corep.rs:1110`

## 四、当前状态（诚实基线）

严格 dispatcher 禁止 NonQuantized/mapping/SU2 错误落入 legacy fallback 后：
- **20,921 ok / 301 fail**
- 全部 301 个为 direct-path `non_quantized`
- Mapping failure = 0
- 200 测试通过

### 301 个 failure 的 SG 分布

| SG | 数量 | 备注 |
|----|------|------|
| 13 | 64 | P2/c, Hall 73→72: 0 transform candidates |
| 14 | 36 | P2₁/c, Hall 82→81: 0 transform candidates |
| 92 | 24 | P4₁2₁2 |
| 95 | 24 | P4₃2₂ |
| 144 | 12 | |
| 145 | 12 | |
| 169 | 12 | |
| 170 | 12 | |
| 141 | 10 | |
| 142 | 10 | |
| 227 | 10 | Fd3m |
| 228 | 10 | Fd3c |
| 88 | 8 | |
| 122 | 8 | |
| 203 | 8 | |
| 其他 | <8 each | |

SG13+14 占 100/301 (33%)。

## 五、已尝试但撤回的修改

### G→H gauge parity（已撤回）

尝试添加 Z₂ gauge mapping: ε(h) = U_h^(G) / U_h^(H) ∈ {+1,-1}

失败原因：比较了**不同帧中的不同旋转**的 SU(2) lifts：
- `u_sq_g` 来自 `b.rot²`（G 帧，未变换）
- `h_sq_lift` 来自 `sq.rot`（data-Hall 帧，变换后）

当 `setting_xf` 改变旋转轴（signed-permutation）时，这两个旋转是**不同的**，比较它们的 SU(2) lifts 没有意义。

修复方向：要么全在 H 帧内计算 `central`，要么全在 G 帧且显式变换到 H。

## 六、尚未完成的工作

1. **h_seitz 仍来自 MSG frame**（`corep.rs:210`）：未切换到 canonical `ops_from_hall`
2. **SG13/14 Hall 对**：`find_setting_transform` 返回 0 candidates。需要通用仿射求解器（非 signed-permutation）
3. **LG filter 混帧 fallback**（`wigner.rs:497`）：transform 失败时回退 raw MSG 帧 + canonical translations
4. **旧 `standard_setting_transform` fallback 不安全**（`corep.rs:711`）：未验证 target Hall，使用 `.next()`
5. **G→H gauge parity**：需要在**同一帧内**计算 `central`，而非跨帧比较
6. **301 个 non_quantized per-term trace**：按 (SG, k, b² rotation, W value) 聚类

## 七、后续计划

### 第一步：per-term Wigner trace for 301 non_quantized

在 direct anti-coset 路径添加诊断，对每个 failing case 记录：
- SG, k-point, irrep label, dim
- Per-term: b index, b.rot, b².rot, sq_spin_idx, local_idx
- spin character χ₀, Bloch phase, term contribution
- Cluster by W value pattern (e.g., W=0.5 vs W=±0.25 vs W=irrational)

按 pattern 分组后，对每个 pattern 的 representative case 深入分析。

### 第二步：修复 identified pattern

根据第一步的模式针对性修复，而非盲目全局修改。

### 第三步：修复 SG13/14 Hall 对（100 cases）

需要通用仿射求解器（rational basis，非仅 signed-permutation）。
优先级高因为占 1/3 的失败。

### 第四步：h_seitz canonical frame + G→H gauge parity

在帧完全统一后，重新实现 gauge parity（此时所有旋转在同一帧，比较有意义）。

### 第五步：清理不安全 fallback

- 移除 LG filter 混帧 fallback → Result 返回
- 移除旧 `standard_setting_transform` 不安全 fallback

### 验收标准
- 总样本数不变
- mapping failure = 0（已达成）
- non_quantized → 0 或确认为物理上真正的非量子化
- 反酉 LG 非空时必须量子化为 0, ±1
