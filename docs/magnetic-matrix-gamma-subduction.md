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

## 7. R32 增量：受限双 H 的 Type-C 目标（`subduce_matrix_gamma_type_c`）

第 32 轮在同一模块新增**受限** pair 入口，**不改** R29 的
`MatrixGammaModel`、`MatrixSubductionPiece`、`MatrixSubductionError` 或
`subduce_matrix_gamma`：

```text
pub fn subduce_matrix_gamma_type_c(
    model: &MatrixGammaModel,
) -> Result<MatrixTypeCSubduction<'_>, MatrixTypeCError>
```

输入模型仍只由 `verify_matrix_gamma_input` 构造；新入口不重复认证 source。目标域为：
磁群 Γ corep **恰为 Wigner type C、维数 2**，其 H  unitary restriction 恰为**两个互不等价、
非 spinor、非复合、维数 1、实 `±1`** 的 H irrep，反幺正 coset 交换两者。新结果
`MatrixTypeCSubduction` + `MatrixTypeCPiece` 用 `h_sources() -> [&'static str; 2]` 保留两个
真实 H 身份；`MatrixTypeCError` 是独立错误 enum，旧公开 enum 没有新增变体。

### 7.1 目标构造与重数（与 Type A 的区别）

- 每个 unitary class 的 H 字符由**两个独立 H record** 的已认证一维矩阵表给出
  （`certify_irrep_matrices(record, 1)`，其 typed selected-arm 行与 legacy 行仍按 R29 口径
  交叉核对），并用 summary 的 Γ 字符列独立核对 `trace(A(h)) = χ1(h)+χ2(h)`。
- `a0` 是首个反幺正 class；实际计算 `σ(h) = a0⁻¹ h a0`（model class 合成 + 格归约），并检查
  `χ1(σ(h)) = χ2(h)`、`χ2(σ(h)) = χ1(h)` 以及 `χ1(a0²) = χ2(a0²) = s ∈ {±1}`。
- 目标表是 `D(h) = diag(χ1(h), χ2(h))` 与规范 coset 矩阵
  `X = [[0, 1], [s, 0]]`（`X² = s·I = D(a0²)`，因此 **`s = -1` 的 `-I` 平方也是合法分支**），
  反幺正类为 `D(a0 h) = X·D(h)·K`。随后对全部有序对重新检查
  `D(g)·conj(D(h)) = D(gh)`。
- 该 type-C corep 的复线性实 Hom 空间是 **2 维**：`I` 与 `diag(i, -i)` 都满足幺正与
  反幺正半线性自 intertwiner 方程。**这不是两份 copies**：`matching_copies` 只保留复列空间
  两两正交（`Uᵢ†Uⱼ = 0`）且受源维数上限约束的 copy，因此这里 `multiplicity = 1`。
  生产代码不把 `kernel.len()` 当重数，也不使用实核维数。
- 成功条件：`multiplicity · dimension = 2`，唯一 copy 的 `U†U = I₂`、`|det U| = 1`
  （满复秩），并用非实探针逐类验证 `S(g)∘φ = φ∘T(g)`，`φ(z) = U z`，反幺正分支共轭像向量
  `U z`。`X` 与 `-X` 是同一目标的合法 coset gauge（整体相位可吸收），不作为歧义。

### 7.2 真实见证（全部 `RequireExact`、`map = I`、`shift = 0`、双 η）

| 源 / 凝聚 / row | UNI | H | H 对（BC） | s = χ(a0²) | 结果 |
|---|---|---|---|---|---|
| SG89 `GM5` / `GM2` row0 | 749 | SG16 P222 | `GM3`(Γ4) + `GM4`(Γ2) | **+1** | 目标 label `GM3 + GM4`，dim 2 mult 1，4 个反幺正类；双 η 成功 |
| SG89 `GM5` / `GM4` row0 | 751 | SG21 C222 | `GM3`(Γ4) + `GM4`(Γ2) | **-1** | 同上；`X² = D(a0²) = -I` 的 `-I` 分支被真实数据动态见证 |

两个见证的 unitary restriction、`a0` 扭曲配对、完整有序对 law、非零复 `U`、`U†U = I₂`、
`|det U| = 1`、covered 2 与复数向量重建都由实际运行执行。UNI749 的 `s=+1` 与 UNI751 的
`s=-1` 由模块单测直接读取 `χ(a0²)` 并核对 `X² = sI`。旧入口对 UNI749 仍返回
`NoMatchingTarget { examined = ["GM1", "GM2"] }`（既有 R29 测试保持），不会静默丢掉第二个 H。

### 7.3 错误与边界（如实）

- 域外排除以 `rejected` 具名列入 `NoMatchingTarget`：spinor H、维数不等于 1、
  复合 H、非 Γ H、有限非实字符、有限但非 `±1`、source 数不是 2。真实数据中的两个见证
  不受这些排除影响。
- 未决/不可用数据是独立错误，**不是**“无候选”：缺 H row / 多个 H row、typed 行不可用、
  `NonFinite`（NaN/∞/imag NaN）、两个 H 字符逐操作等价（`EquivalentHSources`）、
  unresolved summary corep、无 Γ、provenance 失败、H 映射失败、以及
  `source.sg` 与 `context.h_sg` 不一致或 summary/record 一维声明矛盾。
- 模块单测动态见证：同一 H pair 的未决拒绝、spinor/维数排除、跨 H 空间组身份拒绝、
  单个 anti 系数篡改（law 与 reconstruction 都拒绝）、twist 字符篡改（twist 守卫拒绝）、
  换另一反幺正代表元仍构造等价目标、错误 raw class 触发 `TargetProvenance`、
  NaN/±Inf/imag NaN → `NonFinite`、有限非实 → `NonReal`、0.5 → `OffGrid`。
- **未动态见证**：复合 H 与非 Γ H 在固定见证域内没有真实 pair 候选，因此只由代码路径和
  排除名单覆盖；`record.dim != 1` 的 summary/record 声明矛盾也没有真实数据负例。
  `s=-1` 有真实见证（UNI751）；`s=+1` 由 UNI749 见证。
- 不扩 type A/B、complex H、non-Γ 源、spinor、full-star、setting solver；无外部 corep
  oracle；Type C 分类仍取自 summary 的 Wigner 结果，本入口没有重算 Wigner `W`。

### 7.4 回归

新增 `tests/magnetic_type_c_subduction.rs` 6 项：两个见证双 η 的低层 law/重建、
高层双 CDML/BC 身份、旧 `subduce` 对 Type-C record 的具名拒绝、旧 Type A 见证不变、
scalar 源拒绝、Type A record 的 `NoMatchingTarget`。模块内新增 10 项单测（见 §7.3；含实际见证 a-扭曲配对逐类检查与 `s=+1/-1` 两个真实分支）。
R32 当时未改变旧 R29 公开行为；R47 的源原点增量见下节。

## 8. R47：二维源的共同原点对齐

在原矩阵认证之后、绑定磁类之前，[矩阵输入](<../src/irrep/magnetic_matrix_subduction.rs>)
复用[标量输入](<../src/irrep/magnetic_input.rs>)的完整 data-Hall 原点证明：源 Seitz 行必须
属于声明的数据 Hall，完整心格双向一致、母群旋转覆盖齐全，再联立全部母群操作求一个
共同 δ，逐操作验证 `τ_parent = τ_source + (I−R)δ mod L_parent`。
不是只在 M 的已用操作上配旋转，也不是每个操作单独择原点；没有新的 setting solver。

`MatrixGammaModel::source_origin_shift()` 保存源到母群帧的见证，与磁嵌入的 `shift()`
分开。Γ 下只变换操作坐标，不改矩阵、行顺序、η、源/凝聚身份或 intertwiner 的定义。
原认证链与 basis/gauge 限制完全保留，不能据此声称原始 PIR 存储已获 provenance 保证。
标量路径共用同一薄封装，计算与原点选择未变。Type-A 和既有受限 Type-C 共用这个模型；
没有扩大目标家族、复合/三维源、spinor、非 Γ 或 full-star 的适用域。

新增错误 `MatrixSubductionError::SourceOriginAlignment { sg, label, reason }`，
`reason` 保留 boxed `MagneticInputError` 并作为 error source；调用者的穷尽 match 需增加
此分支。数据/共同原点不匹配与 Arithmetic 未决不能混为同一种无解证明。

[公开几何/矩阵回归](<../tests/magnetic_matrix_subduction.rs>)固定 14 个母群的 28 个源/记录
组合，以独立 1/48-cell 整数坐标与完整心矢对照共同原点、母群双向覆盖及每个矩阵绑定，
再检验复 intertwiner 重建；[高层回归](<../tests/magnetic_subduction_api.rs>)覆盖双 η 和
CDML/BC 标签。模块测试遍历现有域内二维源，另拒绝一个矩阵群律和字符均不变、但
Seitz 行整体偏离声明 data-Hall 的合成输入。既有“逐操作可放置但无共同原点”、漏旋转、
错误心格、篡改来源的负例现在直接经过共用封装。
这些是固定源/嵌入检查，不是全部磁 record×matrix-source 的 census 或项目整体完成。

## 9. R48：先排除幺正限制正交的目标

原循环会构造所有元数据合格的目标，即使目标在源限制中根本不存在。真实例子是
SG207 GM3 沿 GM4 的 UNI750：源的 H=SG75 限制为 GM1⊕GM2，但扫描随后试图
构造不相关的复数 GM3/GM4 H 分量，因没有独立数据库行而中止整个分导。

现在在原 R27 操作/类来源证明之后、构造目标矩阵之前，仅对完整幺正子群计算
`⟨χ_source|H, χ_target|H⟩`；为零才跳过。Type-A 和既有受限 Type-C 共用此门槛，
每个 M 的幺正平移类只计一次，不用反幺正 trace，也不把 Hom 实核维数当重数。
所有重复 raw 列必须给出同一有限字符；空列、缺值、非有限值、H 身份不符、单位元
维数不符、异常字符范数（A=1、C=2）、非整数或负重叠都保留为具名数据错误。

正重叠不产生目标或重数，也不意味着数据已齐备：原 H 矩阵认证、coset 构造、
复 Hom、满秩等距与逐操作复重建仍全部执行。确实需要的复数 H 分量仍具名拒绝，
没有把“没有独立记录”笼统改成跳过；例如同一 UNI750 的 SG89 GM5 仍在此域外。
输入域、目标家族、嵌入选择与原有公开返回类型均不变。Type-C 对零重叠的拒绝理由
记为 `orthogonal-unitary-restriction`；这不是一般 B/C、复源或高维算法。

## 10. R49：实/复一维 Type-A 目标的显式伴随入口

[目标实现](<../src/irrep/magnetic_matrix_subduction/scalar_targets.rs>)提供
`subduce_matrix_gamma_scalar_targets(&MatrixGammaModel)`；[高层入口](<../src/irrep/magnetic_subduction_api.rs>)
`subduce_scalar_targets(&MagneticSubductionRequest)`复用原有请求解析、四类身份与输入认证。
源仍是现有**非 spinor、非 compound、实二维 Γ**模型；只接受能够完全分解到
**一维非 spinor Type-A Γ**目标的情况，不承诺所有已接受的二维源都能在这个入口分解。
旧 `subduce`、实矩阵目标算子和受限 Type-C API 不变；不在错误后静默换路由。
新报告的 `scope()` / `SCALAR_TARGET_SCOPE`描述本入口，不能用模型的旧 scope 代替。

### 10.1 字符来源不是伪造的独立 irrep

H 字符必须来自唯一的独立 scalar record，或 compound metadata 具名、由
`compound_selected_arm_view`实际公开的**权威 CIR 分量**。不把缺失行变成“无目标”，
不选第一个歧义来源，也不把未存储的共轭行冒充另一条权威 CIR 行。
`HCharacterStorage`保留 owner 的标签、component 位置及 CIR irnumber；
`piece.source().labels()`仅返回该分量自身的真实标签映射，不能借用 owner 的 compound BC 名。

一维情况下字符本身就是矩阵系数，因此允许直接使用复数 χ；这**不是**从高维字符恢复
矩阵的算法。完整 typed H 行逐项验证有限性与单位模，使用 χ 的忠实实二维化
`[[Re χ, −Im χ], [Im χ, Re χ]]`仅复用已有的完整有序对群律检查器。
检查包括心格冗余行的一致性，不从中宣称原 PIR basis/gauge 获证。
随后要求整个 H data-Hall Seitz 表与 typed 行严格双射，再经已有 H setting 映射绑定
各 M 幺正类，并逐一对照全部对应 summary 列。

R27 raw first-Hall 索引证明和 R48 零重叠排除门槛仍在前面执行；summary 的显示坐标
属于 H data-Hall，**不是** M 的原始坐标。独立回归也通过 `magnetic_operation_index`
关联原始带色 Seitz 操作，不把显示平移直接拿去配 M 类。

### 10.2 复数陪集、Hom 与完整重建

取第一个反幺正类 a₀，规范选择 `T(a₀)=K`，则

```text
T(h)     = χ(h)
T(a₀ h) = conj(χ(h)) K
```

右侧的共轭不可省略。复用原 scalar coset 构造器，检验全部半线性乘法；这也验证
`χ(a₀²)=+1`及扭曲共轭相容性。负平方不能用任意相位修成一维 Type A。
不读取反幺正 trace，不从字符重叠推断重数；算术错误仍保留为 Arithmetic。

输出复用已有 `TargetOperator`（复系数和反幺正标志），而非向实矩阵接口塞入丢失虚部
的数据。复 Hom 在原有实线性系统中增加目标虚部的交叉项；源仍是原实模型。
核、归一化、独立 copy 选择、全源 Gram/投影恒等与逐类复向量重建都复用原组件。
必须实证完整覆盖二维源，缺覆盖、非等距、重建失败仍是错误；实核维数不是重数。

[公开回归](<../tests/magnetic_scalar_targets.rs>)固定 SG89/91/97/141 的见证，含 screw、
心格、非零源原点、双 η、CDML/BC 与独立操作重建；同时保留实标量拆分/重数控制、
旧入口对确需复目标的拒绝、二维目标及非二维源的范围拒绝。模块负例覆盖缺失/错误 H
身份声明、缺值/非有限 summary、H frame 不匹配、所有反幺正代表及不允许的负平方。

这不扩展一般 B/C、二维复目标、复合/三维/复源、spinor、non-Γ/full-star、一般
setting/normalizer 或 OP-stabilizer 契约；全项目仍未完成。

## 11. R50：权威复共轭 H 配对的实二维载体

[独立低层入口](<../src/irrep/magnetic_matrix_subduction/conjugate_pairs.rs>)
subduce_matrix_gamma_conjugate_pair 与[高层入口](<../src/irrep/magnetic_subduction_api.rs>)
subduce_conjugate_pair 仍使用同一个已认证实二维、非 compound、非 spinor Γ 源。
只接纳两条已存储的一维非实 H 字符 χ、χ*；每条均通过 §10 的完整 typed 行认证与
H data-Hall 绑定，不从 seed 捏造另一条数据。报告保留两条 H 身份与各自的 HCharacterStorage，
标签必须用各分量的 labels()，不能借用 owner 的 compound BC 名。

ConjugatePairKind 严格分开两种情况：

- MagneticTypeC：存在反幺正陪集，原 summary 是二维 Type C。
- UnitaryPair：没有反幺正操作，原 summary 将两条 H 源并列成 Type-A 标签。
  此时输出的是**两个不等价复普通不可约表示的配对载体**，不是二维复不可约表示，更不是 C。

两者的 H 字符范数都是2。R48 门槛仅为后一种显式双源、纯幺正 packet 增加此规则；
旧单 H 源 A 仍为1，旧 C 仍为2，所有旧入口的候选过滤与结果类型不变。
没有把零重叠变成重数证明，也没有新增自动 fallback。

在规范实二维化 R(z)=[[Re z,−Im z],[Im z,Re z]] 基中取 T(h)=R(χ(h))。
若有反幺正代表 a₀，取 λ=sqrt(χ(a₀²))，并定义 T(a₀h)=R(λχ(h))K。
**此处不能共轭 χ(h)**：被 K 共轭的是实矩阵 R(χ(h))，不是其复参数。
全部有序对半线性群律验证 χ 对 a₀ 共轭的不变性、平方与陪集相容性，不只核验平方。
两个平方根相差的符号可由整体 i 相位的 intertwiner 吸收；不任意修改 H 字符。
无反幺正陪集时仅使用 R(χ)。因此可直接复用 TargetMatrixOperator 与原实目标 Hom，
不需要通用复矩阵框架或新求解器。

每个成功必须有完整的2×2复 intertwiner、全源 Gram/投影恒等与逐操作复向量重建。
单个磁 C 见证的 Hom 实核维数为2，纯幺正配对见证为4，但两者的完整载体重数都为1；
从不把核维数当重数。全目标群律检查器也保留 compose 的 Arithmetic 错误，
不再把溢出说成“乘积不属于群”；此修复有真实前失败/后通过的合成边界回归。

[公开回归](<../tests/magnetic_conjugate_pairs.rs>)包括 SG100/149/141 的固定非平凡嵌入、
原始反幺正平移、心格、非零源原点、双η和双标签。独立关联 raw M 索引与 H 显示列，
核验每个 CIR 字符、全部乘法、复基向量/混合向量重建及 Gram/投影；旧入口仍拒绝这些见证。
模块控制遍历三重/四重的所有反幺正代表，包括非实和负平方；错误 A 扭曲、缺长度、
非有限值、错误代表、篡改单个系数与算术溢出均拒绝。

该入口不接纳两个不同的实 H 字符、复一维 A 目标或单个实二维 H 目标；它们仍由原入口处理。
不扩大源域，也不宣称完成一般 B/C、高维/复合源、spinor、非Γ/full-star、一般 setting、
normalizer 或 OP-stabilizer。完整固定语料覆盖与项目整体完成是不同的口径。

## 12. 存储行认证的算术边界（R51）

源矩阵、实 H 目标以及 scalar/CIR H 行共用的存储行群律 checker，现在将 exact compose
失败保留为 MatrixCertificationProblem::Arithmetic { left, right, reason }。
这是公开错误枚举的新增变体，下游穷尽 match 需要处理它；SourceMatrixUncertified /
TargetMatrixUncertified 外层来源语义不变，内层保留原始算术原因，不能据此否证表示群律。
真正无乘积行、等价重复行矩阵冲突和矩阵乘法违例仍分别报告原有的错误。
合成极值回归不表示静态官方数据曾触发溢出；不扩大矩阵源/目标或数值可表示范围。
纯幺正配对另有 norm=2 的合法控制和错误norm/缺列/非有限字符负例；单源 A 的 norm=1
控制保留。公开配对测试从权威 H 行和原始 M 乘积独立核对反幺正系数 λχ(h)，并明确
区分错误的 λconj(χ(h)) 公式，而不只依赖完整群律/重建最终断言。

## 13. 显式二维 compound Γ 源（P2a / R53）

[新源认证](<../src/irrep/magnetic_matrix_subduction/compound_sources.rs>)
verify_compound_matrix_gamma_input 只接纳非spinor、物理维数2、Γ、DistinctComponentSum，
并要求两条不同的权威非实共轭一维CIR。旧 verify_matrix_gamma_input 和四个旧高层入口
仍拒绝compound；不自动换源或换目标。新高层 subduce_compound_scalar_targets 与
subduce_compound_conjugate_pair 分别只复用原scalar-A与非实共轭pair目标引擎。
低层显式认证的模型也可交给原实A/受限实H-C引擎；每条路线必须给出完整分解，不能拼接
不同路线的重叠成功。没有新增B目标、一般维数框架、setting搜索或依赖。

物理矩阵仍独立检查有限性、形状、正交、单位、raw/typed迹及全部乘法。仅此窄源分支
要求实commutant维数2；普通源与实H目标仍要求1，绝不全局把门槛改成2。每条具名CIR
另以忠实实二维化复用原scalar-H的完整群律认证；各自有限、单位模、维数1，完整Seitz
旋转**和平移**顺序与物理表逐项相同。两行须点点共轭且确实非实，之和给出物理矩阵迹。
只验证总迹不足：常驻反例将两行的同一非实位置同时共轭，保留总复迹/单位模/共轭关系，
但两条CIR各自违反群律。全局共轭反而仍是表示，因此具名身份必须来自不可变权威构造器，
不能从这些代数检查推断。算术错误仍保留Arithmetic，不改说没有乘积或没有表示。

之后复用完整data-Hall、双向心格、旋转覆盖与**一个共同源原点**证明，接原有色群律、Hom、
copy归一化、全源Gram/投影和逐操作复重建。不是只对所选磁群的部分操作找逐行原点。
源物理载体的反幺正作用仍是实矩阵 ηD(p)K；不能让两条非实CIR各自冒充可独立继承ηK。

MatrixGammaModel::source_components() 在旧模型返回None，新模型保留两条完整CIR行及
component/label/irnumber；这些操作坐标是**共同原点平移前的存储源帧**，不是目标H帧。
source_origin_shift 与嵌入shift仍分开。physical owner保留自己的真实CDML/BC，组件不借用
owner的合并BC：例如SG75 GM3GM4的owner BC为Γ2Γ4，不能给GM3/GM4伪造独立BC。
报告scope按显式源家族选择；旧模型/旧报告的原scope字符串不变。新增公开错误
SourceNotCompound、CompoundSourceInvalid，下游穷尽match需处理。这里认证的是发出的
Γ操作表示，不是原始档案矩阵字节、basis或gauge来源的更强证明。

[公开回归](<../tests/magnetic_compound_sources.rs>)对全部10个非零源原点compound，双η、
CDML/规范化BC、scalar与unitary-pair共80次高层调用，独立逐复基/混合向量重建、Gram、
全投影及所有有色乘法；另有SG75磁C、UNI663的B目标缺口、旧四入口拒绝和域外控制。
SG85/86/88/201/203覆盖四分之一/八分之一源原点及C/F等心格。旧高低层BC raw-LaTeX
别名解析差异是另一个既有接口边界；此处的双标签证据使用label(convention)返回的规范拼写。

固定R47几何下的新域单列为69源、3,778记录/源对、双η **7,556次候选**。探索完整测量：
输入全部认证；scalar-A **5,050**（也可由原实A目标引擎完成）、pair **2,466**（磁C1,896、
纯幺正packet570），联合 **7,516**，不相加重复A路线。另 **40** 次/20对仍无支持的目标，
涉及14个已清点非spinor-B UNI；这不是分导成功，也不是数据不存在。UNI663见证现在可
认证源，但本轮不实现T(a₀)=J K的B目标；窄B留给P3a。最终冻结数字与旧域逐字对照见账本。
