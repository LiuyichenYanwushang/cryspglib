# Wigner 共表示算法 — 错题本

> 当前 `wigner.rs` 的框架（M = H ∪ a₀H → Wigner test → type A/B/C）是正确的，但实现细节对 nonsymmorphic + BZ 边界 k 点 + double group 的组合不够。以下是所有已知问题，按严重程度排序。

---

## P0 — 致命：只用 rotation 匹配 b²，没有 Seitz 平移

**问题**：`wigner_classify` 只计算 `(R_{g₀}R_h)²` 的旋转部分，然后按 rotation 查找 H 中的操作。忽略了平移。

**正确做法**：Seitz 乘法是

```
{R₁|t₁}{R₂|t₂} = {R₁R₂ | t₁ + R₁t₂}
```

对于 `b = a₀h`（空间部分 g₀h），

```
b² = (g₀h)² = { (R₀R_h)² | t_c + (R₀R_h)t_c }
```

其中 `t_c = t₀ + R₀t_h`。

**影响**：对 128.406 的 Z=(0,0,1/2) 点：
- 找错 b² 对应的 unitary operation
- 漏掉 `e^{iπ·Δt_z}` 的 Bloch 相位因子
- Wigner sum 符号可能被反转，Type A/B/C 全判错

**需要的数据结构**：
```rust
struct SeitzOp {
    rot: Mat3I,        // 3×3 integer rotation
    trans_num: [i32; 3],  // translation numerator (fractional)
    trans_den: i32,        // translation denominator
    timerev: bool,
}
```

---

## P0 — 致命：`filter_little_group` 对 antiunitary 操作写错了 k 作用

**问题**：对所有操作都用 `Rk ≡ k` 判断。

**正确做法**：
- Unitary: `Rk ≡ k (mod reciprocal lattice)`
- Antiunitary: `-Rk ≡ k (mod reciprocal lattice)`

**代码**：
```rust
let sign: i32 = if ops.timerev[i] { -1 } else { 1 };
(sign * rx - kx_i) % kd_i == 0
    && (sign * ry - ky_i) % kd_i == 0
    && (sign * rz - kz_i) % kd_i == 0
```

**注意**：对 Z=(0,0,1/2)，因为 -Z ≡ Z，这个 bug 被偶然掩盖。

---

## P0 — 致命：`find_rotation_in_ops` 只比较 rotation，不比较 translation

**问题**：nonsymmorphic 群中存在多个 `{R|t}` 共享同一个 R 但 t 不同。在非 Γ 点它们的 character 不同（差 Bloch 相位）。

**正确做法**：匹配完整 Seitz operation，并返回 lattice shift：
```rust
struct MatchResult {
    op_index: usize,           // which H operation
    lattice_shift: [i32; 3],   // L = t_computed - t_stored
}
```

然后 character 乘以 Bloch 相位：
```
χ_effective = χ_stored × exp(i·k·L)
```

其中 `k = 2π·(kx,ky,kz)/kd`。

**Bilbao convention**：`{1|t₁,t₂,t₃} ↦ exp(iπ·t₃)` 对 Z=(0,0,1/2)。

---

## P1 — Type B 维数和 character 全错

**问题**：代码把 Type A 和 Type B 同样处理（dim=d, χ̃=χ）。

**正确做法**：

| Type | W | dim | χ̃(h) unitary | χ̃(a₀h) antiunitary |
|------|---|-----|-------------|-------------------|
| A | +1 | d | χ_i(h) | 需要 intertwiner U |
| B | -1 | **2d** | χ_i(h) + χ_i^A(h) ≈ 2χ_i(h) | **0** |
| C | 0 | 2d | χ_i(h) + χ_partner(h) | **0** |

Type B 是 Kramers 加倍，维数必须是 2d。

---

## P1 — Type C 的 character 公式不对

**问题**：代码写 `2.0 * h_chars[h_idx]`，文档写 `2·Re[χ(h)]`。

**正确做法**：Type C 的 corep 是 D̃ = Δ_i ⊕ Δ_j，其中 Δ_j ∼ Δ_i^{a₀}。

```
χ̃(h) = χ_i(h) + χ_j(h)
```

其中 `χ_j(h) = χ_i(a₀⁻¹ha₀)*`。

只有当 `Δ_i^{a₀} = Δ_i*`（即 antiunitary 把 irrep 映到普通复共轭）时才有 `2·Re[χ(h)]`。Bilbao 128.406 的 Z₁Z₂ 说明 antiunitary 操作把 irrep **配对**，不是简单取复共轭。

**需要的数据结构**：
```rust
fn build_corep_characters(
    irrep_i_chars: &[Complex64],
    partner_chars: &[Complex64],    // for type C
    a0_conjugation_map: &[usize],   // h ↦ a₀⁻¹ha₀ in H
    ...
)
```

---

## P1 — `f64` character 不够，必须用 `Complex64`

**问题**：Bilbao 矩阵中大量出现 `i, -i, e^{iπ/4}, e^{-i3π/4}` 等复数。

**影响**：double group / nonsymmorphic 小群的自然表示是复数的。用 `f64` 存 character 会丢失虚部信息，Wigner sum 会出错。

**建议**：
```rust
use num_complex::Complex64;
// h_chars: &[Complex64]
// w_sum: Complex64
// 最后检查 w.im.abs() < tol
```

---

## P1 — 找不到 b² 时不能静默加 0

**问题**：
```rust
if let Some(sq_h_idx) = find_rotation_in_ops(...) { ... }
// else: silently contributes 0
```

**正确做法**：`b² ∈ H_k` 必须成立。找不到说明 Seitz 乘法、setting、或 operation list 有问题。应该返回 `Result<CorepType, Error>` 而不是静默吞掉。

---

## P2 — antiunitary 操作的 character 不能用 H 的 character 近似

**问题**：文档写 "antiunitary ops without H mapping: character left as 0.0"，但在 Type A 中 antiunitary 操作可以有非零 trace。

**正确做法**：
- Type A: 需要 intertwiner U 满足 `U·Δ(a₀⁻¹ha₀)*·U⁻¹ = Δ(h)`，然后 χ̃(a₀h) = Tr[U·Δ(h)*]
- Type B/C: antiunitary sector 是 off-diagonal block → trace = 0 ✓

**分两层目标**：
1. 只做 Wigner 类型判断 → 只需要正确计算 W（不需要 antiunitary character）
2. 生成完整 character table → 需要构造 intertwiner U

---

## P3 — 缺少 Bloch 相位约定

Bilbao 对 Z=(0,0,1/2) 用的约定类似 `{1|t₁,t₂,t₃} ↦ exp(iπ·t₃)`。

需要在全局声明 Bloch 相位 convention，并一致应用于：
- 小群过滤（操作是否真的在 k 点小群中）
- 特征标查表（lattice shift 带来的附加相位）
- Wigner sum 中的 χ(b²)

---

## 修复优先级

1. **Seitz 乘法** — 定义 `SeitzOp` 结构体，实现完整 `{R|t}` 复合
2. **Seitz 匹配** — 替换 `find_rotation_in_ops`，返回 `(op_index, lattice_shift)`
3. **antiunitary 小群过滤** — `k ↦ -Rk` 
4. **Complex64 character** — `f64` → `Complex64`
5. **Type B 维数加倍** — dim = 2d, χ̃(h) = 2χ(h), χ̃(a₀h) = 0
6. **Type C partner irrep** — 计算 Δ_i^{a₀}，找到配对 irrep
7. **错误处理** — b² 找不到时报错，不静默
8. **Bloch 相位** — 统一的 convention + e^{ik·L} 因子
9. **Intertwiner U** — 完整 matrix table 需要

---

## 一句话总结

> Wigner test 框架对，但实现只适合 symmorphic + Γ 点 + spinless。对 128.406 Z 点必须加 Seitz 平移、Bloch 相位、antiunitary 的 -k 作用、Complex64 character、Type B 加倍。
