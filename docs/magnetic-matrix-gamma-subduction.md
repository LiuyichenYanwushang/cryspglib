# R29 二维实 Γ 磁矩阵分导（`irrep::magnetic_matrix_subduction`）

本文件记录第 29 轮新增入口的行为、认证链、真实见证与边界。它**只**描述新模块；
R25（`magnetic_input`）、R26/R27（`magnetic_subduction`）的公开行为未改动（见 §6）。

## 1. 这个入口回答什么

[`verify_matrix_gamma_input`](../src/irrep/magnetic_matrix_subduction.rs) 接受
**非 spinor、非复合、`k = 0`、真正二维的实 Γ 母群 irrep**（源）与一个磁 isotropy
record（凝聚身份）＋调用者显式 `(map, shift)` ＋可选 `η`，验证嵌入、认证源矩阵并构造
父群 grey 作用；[`subduce_matrix_gamma`](../src/irrep/magnetic_matrix_subduction.rs)
把该限制分解到磁群自身的 Γ corepresentation，给出**矩阵 intertwiner** 与全空间重建。

源与凝聚是**两个独立身份**（同 R25）：`MatrixGammaRequest { source, record, embedding,
time_reversal }`。输出 `MatrixSubduction` 借用已认证的 `MatrixGammaModel`，字段私有，
唯一构造路径就是上述两个函数。

## 2. 为什么不直接把 `matrices()` 当 operation-aware 行

[`IrrepRecord::matrices`](../src/irrep/types.rs) 自己的文档（1178–1186 行）声明它是
legacy raw/source 存储，与 `pir_rotations`/`pir_translations` 的对齐**未经证明**，
长度相同不是证明。模块因此对每个被接受的 irrep 表**运行期**执行：

1. 形状/有限性：`m` 个操作行对 `m × d²` 个矩阵条目与 `m` 个字符；
2. 逐操作 `Mᵀ M = I`——Euclidean metric 与后面的 `U†U = I` 归一化**因此可用**，是被
   检查出来的，不是默认的；
3. `trace(Mᵢ)` 同时等于 legacy 字符行 `characters()` **和**生产 corep 路径读的 typed
   selected-arm Γ 行（`ordinary_scalar_selected_arm_block_trace`，底层是
   `SCALAR_LITTLE_CHARS`）：两套存储互为对照；
4. 在该操作表自身的格上做**完整有序对矩阵群律**——这是判别性检查：一个保 χ 的两行交换
   在这里失败而在第 3 步通过（见 §5 的 fixture）；
5. 实交换子维数 `= 1`（实不可约）。注意实交换子用实数未知量求解，不用
   `hom_rows`（后者把复未知量的实部/虚部当独立量，属另一件事）。

得到的语义只能是**经 χ（两套存储面）＋完整群律＋实不可约性认证的 Γ 等价 representation
basis**，它**不是**“原始 PIR provenance 已被证明”：`matrices()` 的告警涵盖 source 代表元
**以及 basis/gauge**，而 `k = 0` 让“换格代表元的 Bloch 相位为 1”只对代表元那一半成立，
并不能证明 basis/gauge 的对齐。模块不依赖该相位论证：绑定与等价的成立靠上面运行期的
χ 双面一致、群律与 class 绑定**独立建立**，Γ 是这个主张的**适用条件**（非 Γ 源被具名拒绝），
不是存储语义的证明。报告与文档都不得把结论升级为 PIR provenance 保证。

## 3. operation 绑定与 grey 作用

每个磁类用与 R25 **同一** parent frame 与同一 class mapping 映回父群
（`magnetic_input` 的 `ParentFrame`/`mapped_classes`，本轮仅加 `pub(crate)` 可见性），再按
`rotation` 相同、平移相差父群格矢唯一绑定到存储操作行；没有匹配或匹配到不同矩阵的多行
都是具名拒绝（`SourceOperationUnmatched` / `SourceOperationAmbiguous`）。中心化群同一
陪集的多余行在 Γ 必须携带同一矩阵，否则 `LawProductAmbiguous`。

父群 grey 作用为

```text
A(p)    = D(p)                  （幺正类）
A(p, Θ) = η · D(p) ∘ conj       （反幺正类，η = ±1）
```

构造后**逐有序对**检查 `A(g)·[conj(A(h)) if t(g)] = A(gh)`（`SourceLawViolated` 携带
失败对与两侧矩阵），所以“实表示在限制群上成立”与 η 都是被检查的，不是假设的。

## 4. 目标与分解

目标身份**独立**取自 `magnetic_summary` 的 Γ corep 元数据（label、维数、Wigner 型、
来源 H irrep），列操作经 R27 的 raw first-Hall provenance（
`magnetic_subduction::proven_gamma_class_columns`，本轮只加薄 `pub(crate)` 包装）解析。

* 幺正部分：由**该 H irrep 自己的**已认证矩阵，经 summary 用的 H data-Hall 变换、
  `find_seitz`、`build_h_to_irrep_op_map` 绑定；逐幺正类与 summary 自己的字符值交叉核对
  （不一致即 `TargetCorepInconsistent`）。
* 反幺正部分：解 `U_h X = X U_{σ(h)}` 的核（`σ` 是幺正类被规范反幺正代表共轭，kernel
  base 记为 `B`），由 `D(a0²) = λ·B²` 的**正有限**比例 `λ = D(a0²)/B²` 取定标
  `X = √λ·B`（不是 `B/√λ`；比例与定标后整矩阵的平方逐项验证，`λ≤0`、非有限或零槽
  不匹配都拒绝），
  定义 `D(a) = X · U_{a0⁻¹ a} · K`（`a` 是反幺正类，`K` 是复共轭因子）；**不读** summary 的
  pending 反幺正列。随后验平方、共轭相容与全 pair law。
* intertwiner：解实线性（反幺正分支带 conj）方程组的核，每个解按 `U†U = I` 归一化，只保留
  与已选 copy **交叉块为零**（`Uᵢ†Uⱼ = 0`）的 copy，直到覆盖 2 维。
* 成功条件：`Σ mᵢ dᵢ = 2` **且** `V†V = I` **且** `Σ U U† = I₂`（`SourceSpaceNotIsometric`
  与 `IncompleteDecomposition` 是两个不同的具名失败）；再逐类、逐 copy 做复数向量重建
  `S(g)∘φ = φ∘T(g)`，`φ(z) = U z`，反幺正分支共轭的是**像向量 `U z`**。

目标家族：Wigner **type A**、单一非 spinor 源、维数 1 或 2、且其 H irrep 通过同一认证。
type B/C、spinor 源、维数 > 2 的候选不进家族，**在加入 `examined` 之前就被跳过**（因此
`NoMatchingTarget.examined` 只列家族内候选）；已在家族内但认证/绑定/law 失败的候选是具名
硬拒绝，不会静默跳过。

## 5. 真实见证（冻结 record，全部 `RequireExact`）

源认证本身实测：SG 89 `GM5`（8 ops / 32 entries）与 SG 156 `GM3`（6 ops / 24 entries）
的 `trace==characters`、typed 行一致、逐 op 正交、PIR 表全 pair 群律（64/64、36/36）、实
交换子维数 1 全部通过；`canonical_hall_ops` 与 `symmetry_operations_of` 都**严格恒等**映到
PIR 操作表。这是样本级认证，不是全域覆盖率。

| # | 请求（源 / 凝聚 / row / map） | UNI | 结果 |
|---|---|---|---|
| 1 | SG89 `GM5` / `GM1` row0 / `I` | 747 | 目标 `GM5`（H = SG89）dim 2 mult 1；M = 母群本身，目标幺正部分必然等于源，两个 η 结果相同 |
| 2 | SG89 `GM5` / `M1` row0 / `[[-1,-1,0],[-1,1,0],[0,0,-1]]` | 753 | 目标 `GM5` dim 2 mult 1，**8 个反幺正算子**由 M 乘法构造；两个 η 都成功且 intertwiner 不同（η=+ ≈ `(1/√2)[[-i,-i],[-i,i]]`，η=− ≈ `(1/√2)[[-1,-1],[-1,1]]`） |
| 3 | SG89 `GM5` / `GM5` row0 / `I` | 101 | **2 → 1 + 1**：`GM1`(H=SG3) 与 `GM2`(H=SG3) 各 dim 1 mult 1；两个 η 的两个 intertwiner 都改变（η=+ 为 `i·e₂`/`i·e₁`，η=− 为实向量） |
| 4 | SG89 `GM5` / `GM5` row2 / `I` | 10 | **2 × 1**：`GM1`(H=P1) dim 1 **mult 2**，两份 copy 正交归一且张满 ℂ² |
| 5 | SG156 `GM3` / `GM3` row2 / `I` | 1 | 交叉成功：`GM1`(H=P1) dim 1 mult 2，copy = `e₁, e₂` |
| 6 | SG156 `GM3` / `GM1` row0 / `I` | 1279 | **完整三角**（P3m1 本身，6 类，无 anti）目标 `GM3`(H=SG156) dim 2 mult 1；源律只有数值相等（±√3/2 乘积），逐位比较曾误拒 1.11e-16 残差 |
| 7 | SG191 `GM5+` / `GM1-` row0 / `I` | 1471 | **非 echo 的 dim-2 anti 目标**：目标 `GM5` 由 **H = SG177（P622）**诱导（与源母群不同），24 类 12 反；coset 平方比例实测 **1/4**（`X = B/2`），两个 η 的 intertwiner 不同 |

见证 1 是“M = 母群”的 Type I 普通二维目标：目标与源在该 record 必然相同（数学后果），
因此它**只**覆盖“Type I 普通 dim-2 目标可被独立认证并重建”，不作为“目标独立于源”的证据。
见证 2 是非恒等嵌入（map 行列式 2）的**反平移**磁群；见证 6 覆盖完整三角（数值相等而非
逐位）；见证 **7 是真正的非 echo dim-2 目标**：目标诱导自另一个空间群 H=SG177，12 个反幺正
算子与 **非单位比例 `λ = D(a0²)/B² = 1/4`**（即 `B² = 4·D(a0²)`、`X = B/2`）的 coset
定标都被实际执行（旧的 `B/√λ = 2B` 约定在这里会把平方差 16 倍，因此它是可鉴别的反例，
不只是潜在 bug）；见证 3/4/5 覆盖 2→1+1、2×1 与跨空间群成功。

针对 `√λ` 定标的**非单位比例小 fixture**（模块单测）：
`B = [[1,-√3],[√3,1]]`，`B² = 4·R120`，因此 `λ = D(a0²)/B² = R120/(4·R120) = 1/4`，
正确 `X = √λ·B = B/2`（`X² = R120` 且正交），而错误的 `B/√λ = 2B` 的平方为 `16·R120`
被判失败。另有 `λ≤0`、非有限、零槽不匹配与 NaN/Inf/zero-U 的正有限守卫 fixture（合法零
残差必须通过）。为找定向反例还实际跑过一次额外探索：SG168–194、恒等 map、保持
`RequireExact`，**未扩大 setting 候选族**（40 个 dim-2 anti 成功，其中 uni 1469/1471 实测
`λ = 1/4`）；它不是全域搜索或覆盖率声明。

## 6. 边界与未验证项（如实）

* **没有 setting 搜索**。只使用调用者给出的嵌入；源/凝聚身份不硬编码任何 SG/标签。
* 记录 **UNI 139**（SG89 `GM5` row1）与 **UNI 39/41**（SG156 `GM3` row0/1）在 R24 候选族
  内 `certified_absence`，我另试 `I/Wᵀ/W⁻¹/(Wᵀ)⁻¹` 均 `SourceLatticeNotInParent`；
  **UNI 816**（`I_c422`，目标 H = SG97 的 2 维 `GM5`）同样在文档化 map 下
  `SourceLatticeNotInParent`（`Wᵀ·(1/2,1/2,1/2) = (-1/2,0,1/2)`）。这些**特定 SG89/SG156
  记录**因此没有见证；但“H ≠ 源母群”的 2 维 anti 目标并非不可达——见证 7（SG191/H=SG177）
  就是这样的记录，本节上一条只适用于列出的那几个 record，不是对该域的否证。
* `SourceMatrixUncertified` 有专门的模块单测（保 χ 破群律的置换 fixture：SG89 rows 1↔3、
  SG156 rows 1↔2，两个矩阵的迹都不变而群律失败），但**没有**通过公开 API 的动态见证
  （未找到自然失败的非复合二维 Γ 行）。
* 未覆盖：spinor、compound、非 Γ、复数/高维源、type B/C、full-star、setting 自动求解；
  旧极端 `i128` 原语 panic 边界未修。
* 数值语义：矩阵**系数**比较一律用正有限近似 predicate（`|a-b| ≤ 1e-9` 且两侧有限；
  NaN/∞ 永不通过），维度、anti flag、index/column 关联仍是严格结构相等；合法零残差通过。
  该 predicate 在源/target 全 pair law、正交性、trace/字符对照、intertwiner 归一化与
  Gram、copy 交叉块、`V†V`/`ΣUU†`、复数向量重建上统一使用。
* 无公开动态见证的具名分支（不得当全覆盖）：`SourceMatrixUncertified`、
  `SourceSpaceNotIsometric`、`IncompleteDecomposition`、`TargetMatrixUncertified`、
  `TargetSourceUnresolved`、`TargetOperationUnmatched`、`LawProductAmbiguous`、
  `IdentityMissing/IdentityNotIdentity`、`TypedRow*`、`ClassMapping`/`ParentFrameUnavailable`、
  `Arithmetic`、`AmbiguousLabel`；部分由模块单测或间接路径覆盖。
* 本文件中的数字都是这些冻结 record 上的执行结果，**不是**全域成功率或覆盖证明。
* 源码/数据/上游语义变化都必须重新认证。

R25/R26/R27 的既有公开入口、错误顺序与行为未改；本轮只新增模块，并在
`magnetic_input.rs` 加少量 `pub(crate)` 可见性、在 `magnetic_subduction.rs` 加一个薄
`pub(crate)` 包装复用 R27 的 raw provenance。未 commit / 未 push。
