# 无 SOC 共线自旋群：有限幺正商群

本接口独立于 MSG，不改动 `SymmetryOp.time_reversal`、`OperationKind`、旧磁操作搜索或磁群数据库。
`SpinChannelAction::Swap` 表示独立的自旋 π 旋转代表元，**不是时间反演**。所有新操作均幺正：
交换 up/dn 不共轭轨道矩阵，动量为 `W^{-T} k`，没有额外负号。

## 范围与非目标

- 三维周期输入胞、无 SOC、非零共线磁序、已经自旋块对角的两个通道。
- 空间操作必须正规化**输入胞的整数平移格**：`W ∈ GL(3,Z)`，商掉的是输入胞的 `Z³`。
  包括胞内分数平移；不按旋转矩阵去掉这些操作。任意超胞中不一定能恢复原胞的全部点操作。
- Preserve 选自旋恒等代表，Swap 选绕垂直磁轴的固定轴旋转 π。绕磁轴的连续旋转无需对两个
  对角自旋块数值采样。此商群的 Swap²=Preserve **不是**自旋子 π 旋转平方为 +I 的断言。
- 不提供自旋群编号、数据库、不可约表示、非共线扩展、SU(2) 提升或中心相位。
- 输入只有几何与磁矩，不能据此判断实际 Hamiltonian 是否有 SOC、外磁场、轨道磁通，或满足
  任何反幺正约束。这些前提由调用者排除或另行验证；不能把几何证书当作 Hamiltonian 证书。
- 层群（`Crystal::with_layer`）、全零/不可辨识磁序明确拒绝，不暗中选磁轴或复制两个通道。

## 调用

```rust
use cryspglib::{
    Crystal, CollinearMoments, CollinearSpinTolerances, CollinearSpinOperation,
    SpinChannelAction, detect_collinear_spin_group, validate_collinear_spin_group,
};

let crystal = Crystal::new(
    [[2., 0., 0.], [0., 3., 0.], [0., 0., 4.]],
    vec![[0., 0., 0.], [0.5, 0., 0.]], vec![26, 26],
)?;
let moments = [2., -2.];
let tolerance = CollinearSpinTolerances::default();
let discovered = detect_collinear_spin_group(
    &crystal, CollinearMoments::Signed(&moments), tolerance,
)?;

// 显式指定完整的目标子群，而不是只给生成元。
let target = [
    CollinearSpinOperation::identity(),
    CollinearSpinOperation {
        rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        translation: [0.5, 0., 0.],
        channel_action: SpinChannelAction::Swap,
    },
];
let group = validate_collinear_spin_group(
    &crystal, CollinearMoments::Signed(&moments), &target, tolerance,
)?;
let square = group.multiplication_table()[1][1];
assert_eq!(square.operation_index, group.identity_index());
assert_eq!(square.cell_shift, [1, 0, 0]);
# Ok::<(), Box<dyn std::error::Error>>(())
```

磁矩也可传 `CollinearMoments::Vectors(&[[mx,my,mz], ...])`。取最大范数磁矩的方向作为共同轴，
将轴的首个非零分量规范为正；所有原子先检查垂直残差，磁矩大小保留实际向量范数，符号由轴向投影确定。
因此独立的共线容差不能把大小不匹配藏进横向分量。结果保留 `spin_axis()`、`signed_moments()` 和
`collinearity_residuals()`。超过共线容差报 `NonCollinear`，不会默默丢弃横向分量。允许容差内的近共线输入，
但这些检查不等于完整向量在某个指定自旋旋转下的 covariance 证明。

磁矩必须显式传入；不隐式使用或修改 `Crystal` 上附带的磁矩。要用已有向量时，可传
`CollinearMoments::Vectors(crystal.moments().ok_or(...)? )`。

## 容差与搜索

`CollinearSpinTolerances` 三个独立绝对容差（默认均 `1e-5`）：

| 字段 | 单位与用途 |
|---|---|
| `position` | Cartesian 长度；原子几何、Seitz 模格比较与晶格度量验证，有限且 >0 |
| `moment` | 输入磁矩单位；比较带符号大小，有限且 >=0 |
| `collinearity` | Cartesian 垂直磁矩范数，输入磁矩单位，有限且 >=0 |

晶格使用库的一贯列矢量约定 `lattice[cart][vec]`。晶格度量 `G=AᵀA` 的检查为逐项
`|(WᵀGW-G)ij| <= position*(sqrt(Gii)+sqrt(Gjj)) + position²`。
空间 W 不作用于自旋轴，不使用 `det(Q) Q m`。

发现入口复用非磁晶格旋转搜索（自动角容差），以最稀有 species 为锚逐个试探同 species 平移；
全部原子都必须有唯一的几何像并构成双射，再分别检查 `m_target = +m_source` 与 `-m_source`。
磁矩大小不均匀不一定非法：发现会过滤不相容操作，显式验证则报告 `MomentMismatch`。
零矩原子仍参与结构及磁矩映射。若所有 `|m| <= moment`，报 `AllZeroMoments`；向量轴也须可辨识。

每次周期像搜索用 `|d_j-ell_j| <= position * ||row_j(A^{-1})||` 枚举有限整数区间，
最终用 Cartesian 距离判定。**不是**分数坐标分量 round，也不是固定 27 个邻居。
有多个周期像或多个同 species 原子时明确拒绝，不按磁矩选一个几何映射。
第一版每次最多枚举 4096 个整数像，超出报 `ImageSearchLimit`；应先约化过度偏斜基底或减小容差，
而非偷偷裁掉区间。逆矩阵另验证 `B A = I-E` 的行和误差界（含浮点 dot-product 误差）小于 1/4，
用 Neumann 界向外包围真实逆矩阵的行范数；不把有限但错误的 cofactor inverse 当作奇异胞通过证明。

坐标/平移与待浮点解析整数超过 `2^52` 的表示边界明确拒绝。原始位置和平移先拆成整数与小数部分：
整数输运用 `i128`，几何残差只对居中的小数计算，避免大整数代表吞掉小残差；最终 shift 检查 `i64` 范围。
旋转乘法用 `i128` 累加并检查最终 `i32` 表示，不允许 wrap/saturating cast。结构候选搜索中的
Delaunay 基变换和旋转变换也用宽整数 determinant 并拒绝超范围转换；不能等旧内核溢出以后再验证。
无法用 `i32` 完整表达候选旋转时，整个搜索失败，不返回饱和转换后的子集。

小数 affine、Cartesian 距离、晶格度量和 Seitz 组合用 `mul_add` / TwoSum 跟踪保守舍入界，
周期像区间按误差向外扩张。向量的磁轴归一化、投影/垂直残差、范数和跨原子大小比较同样保留
舍入不确定性；两个 `hypot` 数值恰相等不算大小相等证明。无法解析的向量边界（包括零容差）明确拒绝，
发现入口也传播此错误而不是删掉相应操作。只有距离的上界不超过 `position` 才接受；无法分辨的临界带或
严重消去返回 `NumericalRange`，不以 float 恰好为零充当几何证明。输出 δ/ε 是受控舍入下的
f64 数值残差，等式按所述重建容差使用，不称精确实数等式。仅包含 subnormal 向量的
不可可靠归一化磁轴（即使用户设 `moment=0`）也报 `NumericalRange`，不发布非单位轴。

发现后的**整个**候选集合复用显式验证器。容差可能选出不闭合集合，此时返回失败操作，
不删操作、不补生成元、不自动降低容差，也不取所谓最大闭合子群。搜索是数值搜索，
不是所有真实参数/所有超胞的解析完备性证明。

## 输出与群律证书

`ValidatedCollinearSpinGroup` 字段私有，仅有只读访问器。显式验证保留输入的操作顺序、
原始平移整数部分及原始原子坐标代表，不标准化/理想化结构，也不把磁矩在结构原胞中平均。
每个 `MappedCollinearSpinOperation.atom_mapping[i]` 提供：

```text
source i -> target
cell_shift: [i64; 3]
position_residual: fractional [f64; 3]
moment_residual: signed f64
```

满足

```text
W_g r_i + t_g = r_{pi_g(i)} + ell_g(i) + delta_g(i)
mu_g(i) = m_{pi_g(i)} - eta_g m_i
||A delta_g(i)|| <= position; |mu_g(i)| <= moment
```

乘法表的行 g、列 h 表示先 h 后 g，记录所选代表 p、整数 q 与实际 Seitz 残差 ε：

```text
W_p = W_g W_h; eta_p = eta_g eta_h
t_g + W_g t_h = t_p + q(g,h) + epsilon(g,h)
pi_p(i) = pi_g(pi_h(i))
ell_g(pi_h(i)) + W_g ell_h(i) = ell_p(i) + q(g,h)
delta_g(pi_h(i)) + W_g delta_h(i) = delta_p(i) + epsilon(g,h)
```

最后一式按 `position` 检查浮点重建误差，ε 的 Cartesian 范数也必须在 `position` 内；
不将 δ 或 ε 吸收到整数 shift。特别是半平移 Swap 的平方虽然回到单位商类，但 q 为一整格。

验证覆盖单位元及其原子作用、重复模整数、每个有序乘积唯一命中、两侧单位/逆元、
精确乘法表结合律、置换与整数 shift 复合，以及三元组整数 cocycle：

```text
q(g,h) + q(gh,k) = W_g q(h,k) + q(g,hk)
```

所有失败均具型返回，包括具体操作/原子索引或缺失乘积。`NotClosed.product` 的平移是模格的
居中诊断代表，不是可能不可用 f64 还原的 raw 大整数平移；成功产品的完整整数回执始终保留。
`NotClosed` 不隐藏成一个成功子集。
认证成本为 `O(|G|^3 + |G|^2 N + |G| N^2)`（每个几何匹配另受像枚举预算约束）；
不为一般超大超胞引入额外索引/并行框架。非有限值、奇异胞、非法长度/容差在搜索之前拒绝。

## Rustb 消费边界

cryspglib 只提供群与几何。Rustb 必须另外给出、认证完整的轨道 sewing action `P_g`，特别是：

- up/dn 的共同 Wannier 规范；相同轨道标签不保证相同规范；
- 跨通道 action 的复合（含 q 的 Bloch 相位），以及中心映射和实际轨道 action 的一致性；
- 全实空间支撑上的 covariance、Hermiticity、投影幂等性。

有了该轨道表示后，第一版**幺正**投影才是

```text
H'_sigma = sum_g P_g[ H_{eta_g sigma} ] / |G|
H0 = (H_up + H_down)/2; Hz = (H_up - H_down)/2
H0' = sum_g P_g[H0]/|G|
Hz' = sum_g eta_g P_g[Hz]/|G|
```

`Swap` 只选另一通道，不添加复共轭。几何相容性不证明实际复杂 hopping 具备此对称性。
谱关系是 `E_up(k)=E_down(W^{-T}k)`，不是一般 k 点的同点简并。
本次没有修改 Rustb，也没有实现其局域规范、轨道表示或 Hamiltonian 联合投影。

## 可重复检查

```bash
cargo test -p cryspglib --test collinear_spin
cargo test -p cryspglib --lib collinear_spin::tests
cargo test -p cryspglib --doc
cargo clippy -p cryspglib --all-targets -- -D warnings
```

公开回归含 FM 的独立轴、半平移 AFM、带非磁 marker 的 C4 altermagnet、复数 hopping 的
**下游 toy oracle**（幺正与反幺正结果不同）、螺旋/换胞代表、偏斜胞唯一周期像、真实几何和磁矩
残差，以及非法输入/不闭合拒绝。toy H 不冒充 cryspglib 或 Rustb 的生产投影实现。
