# 工作总结与后续计划（供 Codex Review）

## 一、Baseline 演变

| 阶段 | fail | 改动 |
|------|------|------|
| codex 审核后基线 | 68 | 详见上轮 new.md |
| unimodular bases 扩展 | ~200 | SG13/14 Hall73→72/82→81 |
| 容差修正 + 分类规则 | ~140 | W≈0 假阳性、仅 W=0/±1 |
| codex 4 bug 修复 | 89 | LiftRelation, quat_conj, parity frame, is_signed_perm |
| 180° quaternion 修复 | 68 | SG122/199/206/220 |
| strict Seitz matching (codex Phase 1) | **2,698** | 删除 rotation-only fallback，暴露翻译格点问题 |
| **centering-aware delta_in_lattice** | **1,328** | F/I/C/A/R 中心化群全部修复 (-1,370, 51%) |
| **当前（诊断中）** | **1,328** | 全部剩余 = 原始群 origin shift 问题 |

## 二、codex 审核后的所有修改

### Phase 0-3（上轮 new.md 已覆盖）
- `52c9800`: `LiftRelation` enum (`Same`/`EBar`)
- `34a57b2`: 4 bug 修复 (central 判断, conj_pauli→quat_conj, parity oracle frame, is_signed_perm)
- `bde55f8`: 显式 -U²: `neg_pauli(&su2_compose(&u_b, &u_b))`
- `7ed502b`: 180° quaternion 修复 (Q+I column)
- `032a7af`: `signed_perm_to_quat()` 独立函数

### Phase 1: Strict Seitz Matching（本轮核心工作）

#### 移除 rotation-only fallback
- `find_sq_spin_lg_first` 改为严格 Seitz 匹配（旋转 + 平移 mod lattice）
- 暴露 **2,698** 翻译匹配失败
- 根因：`ops_from_hall` 只有 1 个恒等操作 (0,0,0)，不包含 centering 变体

#### 诊断：Hall 翻译格点数据
- 添加 `DIAG HallTrans` 输出，确认**所有 Hall 群都只有 1 个 I=(0,0,0)**
- 翻译格点**不能**从 `ops_from_hall` 中提取，必须从 Bravais 格点类型推导

#### `delta_in_lattice` 重写（`81a4d84`）

旧版：brute-force 枚举 n_i ∈ [-3,3] × 4 vectors → 只能检查 Z³

新版：
```rust
fn delta_in_lattice(delta: &[f64; 3], centering_shifts: &[[f64; 3]]) -> bool {
    // 1. 枚举 centering_shifts 所有组合 (mod 1 去重)
    // 2. 检查 delta 小数部分是否匹配任何组合
    // Z³ 隐式处理
}
```

- F-centering: 最多 2³=8 种组合
- P: 0 种 → trivially Z³ only

#### `centering_shifts_for_sg(sg: u8)` 新增

81 个非原始群的 centering vectors，按 SG 号 match：

| Centering | SGs | Shifts |
|-----------|-----|--------|
| I (Body) | 23,24,44-46,71-74,79-80,82,87-88,97-98,107-110,119-122,139-142,197,199,204,206,211,214,217,220,229-230 | (½,½,½) |
| F (Face) | 22,42-43,69-70,196,202-203,209-210,216,219,225-228 | (0,½,½), (½,0,½), (½,½,0) |
| C | 5,8-9,12,15,21,35-37,63-68 | (½,½,0) |
| A | 20,38-41 | (0,½,½) |
| R (hex) | 146,148,155,160-161,166-167 | (⅔,⅓,⅓), (⅓,⅔,⅔) |
| P | 其余 149 | (none) |

#### 调用链改动
1. `find_sq_spin_lg_first`: 参数 `canonical_pure_translations` → `centering_shifts`
2. `wigner_classify_spinor_direct_anti_diagnostic`:
   - `all_canon` (Z³ + centering) → `centering_shifts_for_sg(ctx.sg)` (仅 centering)
   - `canonical_pure_translations` 参数标记 `_unused`
3. `wigner_classify_spinor_primary`: 第二调用点同样使用 `centering_shifts_for_sg(ctx.sg)`

#### 效果
**2,698 → 1,328** (-1,370, -50.8%)。全部 F/I/C/A/R 中心化群通过。

---

## 三、剩余 1,328 个失败：原始群 origin shift 问题

### 失败分布
全部为 **primitive (P)** 群：

| SG | 数量 | 群符号 |
|----|------|--------|
| 85 | 80 | P4/n |
| 86 | 80 | P4₂/n |
| 201 | 56 | Pn-3 |
| 151 | 48 | P3₁12 |
| 179 | 45 | P6₅22 |
| 125 | 40 | P4/nbm |
| 126 | 40 | P4/nnc |
| 129 | 40 | P4/nmm |
| 130 | 40 | P4/ncc |
| 其他 | ~ | P6₁22, P6₄22, 等 |

### 失败模式（以 SG85 为例）

**SQ_FAIL_DETAIL 诊断输出：**
```
SG85:
  b_bilbao: rot=C4z+  trans=(0.5, 0.5, 0)
  sq:       rot=C2z   trans=(0, 0, 0)
  centering_shifts=[]
  sg_setting_origin=[4.0, 1.0, 1.0]
```

**SG85 spin table (Bilbao 约定)：**
```
[1] C4z+:  (0.5, 0,   0)
[2] C2z:   (0.5, 0.5, 0)
[3] C4z-:  (0,   0.5, 0)
```

### 根因分析

1. **b_bilbao 有 C4z+ at (0.5, 0.5, 0)**，但 spin table 有 C4z+ at (0.5, 0, 0)
2. **平移差 Δb = (0, 0.5, 0)**——不是 Z³ lattice vector
3. `to_bilbao(C4z+, (0.5, 0.5, 0))` 返回 `(0.5, 0.5, 0)`（未改变）
4. `origin=[4,1,1]` mod 1 = `(0,0,0)` → `(I-C4z+)*origin` 全是整数 → **to_bilbao 对小数平移是 no-op**

5. **验证**：spin table 自己的 C4z+ at (0.5, 0, 0) 平方：
   ```
   {C4z+|(0.5,0,0)}² = {C2z | C4z+*(0.5,0,0) + (0.5,0,0)}
                      = {C2z | (0,0.5,0) + (0.5,0,0)}
                      = {C2z | (0.5,0.5,0)}  ← 匹配 spin table!
   ```
   内部一致。问题出在 **b_bilbao 的平移本身就不对**。

6. **SG_SETTING_ORIGIN 数据问题**：
   - `extract_sg_settings.py` 从 ISOTROPY `iso.zip` 提取 origin 值
   - 注释说 "Convert origin to rational form: fractions with denominator 1,2,4" ——但**代码中实际未做转换**
   - 值如 `[4, 1, 1]` 写为 `{v:1.1f}`，直接存储为 f64
   - mod 1 后全是 0，无法产生非零平移修正

### 可能的修复方向

**方向 A：修复 origin 数据生成**
- 理解 ISOTROPY `isotropy_origin` 的真实含义（可能是 Cartesian 坐标或缩放后的值）
- 正确转换为 fractional 坐标
- 需要深入理解 ISOTROPY 数据格式
- 优点：根本性修复
- 缺点：需要理解外部数据格式，可能涉及 cell parameters

**方向 B：Rotation-only match + Bloch 相位修正**
- 对原始群，当无法完全 Seitz 匹配时，接受 rotation match
- 通过 Bloch phase 修正补偿平移差：
  ```
  χ_corrected = χ_spin × exp(i·2π·k·(t_sq − t_spin))
  ```
- 优点：物理正确，phase correction 已存在于旧代码
- 缺点：恢复 rotation-only fallback（codex 曾建议删除）

**方向 C：接受 delta 在合理容差内的任何匹配**
- 将 delta 比较容差从 Z³ 放宽到 1/2 和 1/3 分母
- 但这是 hack，不解决根本问题

---

## 四、关键代码位置（当前状态）

### 新增/修改

| 功能 | 文件:行 | 说明 |
|------|--------|------|
| `centering_shifts_for_sg()` | `wigner.rs:~2208` | SG→centering vectors 查找 (match 81 非 P 群) |
| `delta_in_lattice()` (新版) | `wigner.rs:~2266` | centering-aware 格点检查 (nested in `find_sq_spin_lg_first`) |
| `find_sq_spin_lg_first()` | `wigner.rs:~2256` | 参数改为 `centering_shifts` |
| Direct anti: centering 构建 | `wigner.rs:~1775` | `centering_shifts_for_sg(ctx.sg)` 替代 `all_canon` |
| Primary path: centering 调用 | `wigner.rs:~2955` | `centering_shifts_for_sg(ctx.sg)` |
| SQ_FAIL_DETAIL 诊断 | `wigner.rs:~1810` | 打印 b→b_bilbao→sq→spin table 全链路 |
| HallTrans 诊断 | `corep.rs:~2130` | 打印 Hall 恒等操作翻译数 |

### 未改动（上次 codex 遗留）

| # | 问题 | 状态 |
|---|------|------|
| 6 | quaternion 向量符号回归验证 | ❌ 未做 |
| 7 | `antiunitary_square_pauli()` 顺序错误 | ❌ 仅 legacy path |
| 8 | 注释 χ((a₀b)²) 实际循环计算 χ(b²) | ❌ 注释问题 |
| 9 | `neg_rot` 等未使用变量 | ❌ 未清理 |

---

## 五、测试命令

```bash
# 当前基线诊断
cargo test --package cryspglib --release diagnose_wigner_sources -- --nocapture

# 查看 SQ_FAIL 详细诊断（前 5 个 primitive 失败）
cargo test --package cryspglib --release diagnose_wigner_sources -- --nocapture 2>&1 | grep -A10 "SQ_FAIL_DETAIL"

# 编译检查
cargo check --package cryspglib
```

---

## 六、后续计划

### 第一步：修复剩余 1,328 个原始群失败

优先尝试**方向 B**（rotation-only match + Bloch phase correction）：
1. 在 `find_sq_spin_lg_first` 中，当 `centering_shifts` 为空且 rotation match 成功时，不返回 None
2. caller 中计算 Bloch 相位修正 `exp(i·2π·k·delta)`
3. 应用到 Wigner sum 的字符贡献

如果方向 B 效果不佳，深入**方向 A**（修复 origin 数据）。

### 第二步：接入生产路径
- `compute_corepresentation()` 的 spinor 路径仍传 `setting_xf=None`
- 需将修正后的 setting 接入生产管线

### 第三步：清理遗留问题
- `antiunitary_square_pauli()` 顺序
- quaternion 向量符号回归验证
- 未使用变量清理

### 验收标准
- `square_not_in_spin_table` → 0 或确认为数据约定差异
- 全量测试通过
- 总样本数不变
