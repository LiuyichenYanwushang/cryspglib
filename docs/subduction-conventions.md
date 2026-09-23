# 完整分导的坐标、表示与来源约定（任务 1）

本文件是 `src/irrep/subduction.rs`（任务 3–12）的契约。它固定坐标帧、表示语义、
来源身份、精度与首版不支持边界；**不含数值算法**。配套任务卡见
`docs/full-irrep-subduction-plan.md`。下面所有数字都取自 pinned
`isotropy_subgroup/iso.zip`（`data_isotropy.txt`）与随包 `iso` 9.6.1 在
`SET I ALL OR 1` 下的输出，逐条可复现。

## 1. 两个独立输入

| 角色 | 输入 | 用途 | 不能混用 |
|---|---|---|---|
| 凝聚表示 | `(parent_sg, condensing_ml, IsotropyDirection)` | 决定子群 `H` 及其**具体嵌入** | 它的方向决定嵌入，不只是决定"哪一族子群" |
| 被分导表示 | `(parent_sg, probe_ml)` | 要分解的对象 | 它不改变子群、不影响嵌入 |

固定子群的键是**等距记录序号** `IsotropySubgroup::ordinal`（0-based，`data_isotropy.txt`
记录序），不是 `(sg_H, 方向标签)`：同一 `sg_H` 可以由不同方向、不同嵌入得到
（见 §2 例 B 与例 C）。任何以 `(sg_H, label)` 为键的缓存都必须在嵌入层再次校验。

## 2. 坐标约定（列向量）

以**已验证**的仿射变换为准（子群 conventional 帧 → 母群 conventional 帧）：

```text
x_G = T x_H + o
R_G = T R_H T^-1
t_G = T t_H + o - R_G o
k_H = T^T k_G
```

- 子群 conventional 基矢按**行**存为 `B`，则 `T = B^T`。推导：`x_G = Σ_i x_H,i·B_i = B^T x_H`。
- `k_H = T^T k_G` 由 Bloch 相位配对不变性 `k_G·x_G = k_H·x_H` 得到；origin 不改变 k，
  只改变操作的平移代表元。
- 操作一律用 Seitz 记号 `{R|t}`，`R` 行主序整数，`t` 为**该帧**的分数平移。

`SHOW XYZ` 打印的正是这个 `x_H = T^-1 (x_G − o)`（**不是** `B^-1`，两者只在 `B`
对称时相同），可作为 `T`/`o` 的函数级 oracle。

存储数据到 `B` 的链路（`W` 是 `data_isotropy.txt` 的 `isotropy_basis`，行向量，
母群 **primitive** 帧；`P_*` 来自 `irrep::isotropy::parent_primitive_basis(*)`）：

```text
W · P_parent      = 子群 primitive 基（母群 conventional 坐标）
P_sub · B         = 同一个格的另一组基
⇒ W · P_parent    = U · (P_sub · B)          U ∈ GL(3, Z), det U = ±1
⇒ B               = P_sub^-1 · U^-1 · W · P_parent
```

**格相等或体积相等不足以证明 `U = I`**：`U` 是"有向对应"，必须显式求出并验证
（§4）。`U` 的搜索范围只取小的 signed permutation/shear 候选，并且每个候选都必须
通过 §4 的完整验证；若有多个候选通过（只差子群自同构、会置换 irrep 标签），
按任务 4 走冻结元数据路径，不取第一项。

### 手算例 A：221 `GM4+` / `P1` → #83 `P4/m`（首个黄金用例）

存储记录（ordinal 12400）：`W = [[0,0,1],[0,-1,0],[1,0,0]]`，`origin = [0,0,0,1]`。
官方 `SHOW BASIS`/`SHOW ORIGIN`：`(0,0,1),(0,-1,0),(1,0,0)` / `(0,0,0)`；
`SHOW ELEMENTS`：8 个陪集代表
`(E|0,0,0), (C2x|0,0,0), (C4x+|0,0,0), (C4x-|0,0,0), (I|0,0,0), (SGx|0,0,0), (S4x-|0,0,0), (S4x+|0,0,0)`。

母群 221 与子群 #83 都是 P 点阵，`P_parent = P_sub = I`，官方基与 `W` 逐项相同
⇒ `U = I`，`B = W`，`T = B^T = W = [[0,0,1],[0,-1,0],[1,0,0]]`（对称且 `T^-1 = T`），
`o = (0,0,0)`。`k_G = (0,0,0)` ⇒ `k_H = T^T k_G = (0,0,0)`。

操作映射核对（真算）：#83 的规范 setting（P4/m，唯一轴 c）里四重轴是 `C4z+`，
即 `R_H = [[0,-1,0],[1,0,0],[0,0,1]]`；代入 `R_G = T R_H T^-1`（`T^-1 = T`）：

```text
T R_H   = [[0,0,1],[-1,0,0],[0,-1,0]]
R_G     = (T R_H) T = [[1,0,0],[0,0,-1],[0,1,0]] = C4x+   ← 官方打印的第 3 个代表
```

即"子群自己 c 轴上的四重轴"在母群帧中是 x 上的四重轴（`B` 的第 3 行 `(1,0,0)`
正是母群 x 轴）——`T` 必须显式携带方向信息，不能只看格。

**陷阱一（已核对）**：`SHOW ELEMENTS` 打印的 `C4x+` 等符号是**母群帧**里的名字，
已经含有 `T` 的效果；把它再代进 `R_G = T R_H T^-1` 是错的。验证时必须用子群
**自身 setting** 的操作做映射（判据见 §4）。判别依据：例 B 打印 `C2a`/`SGda`
（立方母群的 ⟨110⟩ 对角轴名），而 #12 自身 setting 的同一操作叫 `C2y`/`SGy`。

**陷阱二（已核对，2026-09-21）**：官方打印（以及 `data_space.txt` 的 `ipoint_op`）
的旋转矩阵是**行作用** `x' = x M`，而本引擎（和 `SymmetryOps`）是列作用
`x' = R x`，解码时**必须转置** `R = M^T`。判别证据：`ipoint_op` 里
`C4z+ = [[0,1,0],[-1,0,0],[0,0,1]]` 与它旁边标注的 `4[001]`（右手 +90°）只有按行
作用才一致；SG 167 `GM3+` P1 → #15 `C2/c` 是 fixture 里唯一"转置不闭合"的记录，
不转置时打印旋转不保持子群格子、转置后保持。`verify_isotropy_operations.py` 用
`LEGEND_ANCHORS`、母群帧格子保持、以及"逆映射回子群基必须得到整数旋转 + 保持子群胞
+ 平移回到 1/12 源网格"三条断言钉死这件事，fixture 存的是列作用矩阵。

### 手算例 B：221 `GM4+` / `P2` → #12 `C2/m`（非对称换基、`U ≠ I`）

存储（ordinal 12401）：`W = [[0,0,1],[1,0,0],[0,1,0]]`，origin `(0,0,0)`。
官方：`B = (1,-1,0),(1,1,0),(0,0,1)`（`det B = 2 = Z_sub·Size/Z_parent = 2·1/1`），
origin `(0,0,0)`。#12 是 C 心，`P_sub = [[1/2,1/2,0],[-1/2,1/2,0],[0,0,1]]`，于是

```text
P_sub · B = [[1,0,0],[0,1,0],[0,0,1]]        (逐行验算；等价地 B^-1 = P_sub)
U         = W · (P_sub · B)^-1 = W           (det U = 1)
T         = B^T = [[1,1,0],[-1,1,0],[0,0,1]]
```

所以对中心化子群，`B` 既不是 `W` 也不是 `W·P_parent`；`T` 由 §2 的链路给出，
`U` 必须显式记录。官方 `SHOW ELEMENTS` 该行是
`(E|0,0,0), (C2a|0,0,0), (I|0,0,0), (SGda|0,0,0)`（4 个陪集代表 = 2/m 的点群商）。
用子群自身 setting 核对同一映射：`R_H = diag(-1,1,-1)`（`C2y`）与 `diag(1,-1,1)`
（`SGy`），`T^-1 = [[1/2,-1/2,0],[1/2,1/2,0],[0,0,1]]`，

```text
T R_H T^-1 = [[0,1,0],[1,0,0],[0,0,-1]]   = C2a    (绕母群 [110] 的二重轴)
T R_H T^-1 = [[0,-1,0],[-1,0,0],[0,0,1]]  = SGda   (⊥ [110] 的镜面)
```
（前者取 `R_H = C2y`，后者取 `R_H = SGy`。）

与官方打印逐项一致。注意 `T` 的行列式为 2（`det B = 2`），**不是** unimodular；
`o = (0,0,0)` 且 C 心的平移 `(1/2,1/2,0)_H` 经 `T` 映到 `(1,0,0)`（母群格矢量，
与 `Size = 1` 一致）。

k 换算示意（同一 `B`）：若 `k_G = (1/2,1/2,0)`（母群 M 点），则
`k_H = B k_G = (0,1,0)`。别急着约化——C 心的子群倒格在 conventional 坐标里要求
`h1 + h2` 为偶数，所以 `(0,1,0)` 是真正的非 Γ k（单斜 Y 点），不是 (`1/2,1/2,0)` 的
零化。等价判定必须带子群 centring 消光，这正是任务 7 的范围。

### 手算例 C：221 `GM3+` / `P1` → #123、`C1` → #47（`W`/`o` 相同、嵌入不同）

两条记录的 `W = I`、`origin = (0,0,0)`，但官方 `SHOW ELEMENTS` 给出 **16** 与 **8**
个陪集代表（#123 `P4/mmm` 的点群商 16，#47 `Pmmm` 的 8）。**`(W, origin)` 不含点群
信息**，嵌入必须由子群自身的规范操作集经 `T/o` 搬入母群帧后验证得到。

## 3. L_G 与 L_H

- `L_G` = 母群平移格；`L_H` = 子群平移格（`W` 张成，母群 primitive 坐标）。二者满足
  `L_H ⊆ L_G`，指数 `= |det W| = Size`。
- 判断"某个操作属于母群"用 `L_G`；判断"两条操作在 H 内是同一元素"、做闭包/逆元去重
  必须用 `L_H`。**模 `L_G` 相同不蕴含模 `L_H` 相同**（超胞/klassengleiche 情形），
  任务 3 的测试必须钉死这一点。
- 归约 `v → v + shift·L_H` 必须把 `shift`（整数组合）返回给调用者：非 Γ 的 Bloch 相位
  与"操作代表元变更"都需要它；只返回代表元、丢弃 `shift` 的设计不成立。
  `wigner.rs::ExactSeitzReduction.lattice_shift` 是同一形态的既有先例，但它只针对
  canonical Hall 表、且分母固定 12（§7），不能直接复用。
- **行基格子的坐标映射是 `(L^T)^-1`，不是 `L^-1`**（`L` 的行是基矢，点写成
  `L^T · 坐标`）。两者只在 `L` 对称时相同，所以只用对角/对称格子的测试查不出这个
  转置错误；`subduction.rs` 的 `Lattice` 存的就是 `(L^T)^-1`，并有非对称格子的回归。

## 4. 嵌入的验证判据（任务 4 的前置契约）

候选 `(T, o)` 只有同时满足下列全部条件才能产生 `SubgroupEmbedding`：

1. 子群规范操作集（`SG_DATA_HALL[sg_H]` 的严格读取，§7）经 `(T,o)` 搬入母群帧后，
   每个完整 Seitz 操作都落在母群操作集内（模 `L_G`，逐操作配对，不是只比旋转）；
2. 该集合在模 `L_H` 下闭包、含逆元；映射后的**去重代表元数**必须等于由记录算出的
   期望数——中心化胞的 conventional 操作数 ≠ 点群商阶数，所以期望数要按"子群自身
   规范操作集取模 `T_H` 后的大小"来算，并由任务 2/4 的 fixture 逐个钉住。已核对的
   见证：例 A 8（#83，4/m 商 8）、例 B 4（#12，2/m 商 4）、例 C 16（#123）与 8（#47）、
   #139 `M1-` P1 16（#126）。`Size > 1`（超胞嵌入）的期望数由任务 4 用 fixture 与
   子群号识别共同确定，不得在实现里用点群阶直接顶替；
   **代表元必须由映射后的原始操作直接对 `L_H` 去重得到，不能先把它们对 `L_G` 约化**：
   `L_H ⊆ L_G` 只保证模 `L_H` 相同 ⇒ 模 `L_G` 相同，反向不成立。SG 139 `M1-`
   P1 → #126 是钉死的见证：反演 `(-I | 5/2,5/2,5/2)` 对 `L_G`（I 心）约化为零平移，
   但对 `L_H` 的代表元是 `(1/2,1/2,1/2)`。旧实现仍返回 16 个代表元，
   其中 8 个平移类错误，随后分导报 `OperationNotInCharacterRow`。
   `SubgroupEmbedding::operations()` 保留对 `L_G` 约化的母群成员视图；
   `representatives()` 保留对 `L_H` 的代表元，母群特征标与子群拉回均使用它；
3. 子群号识别结果等于 `record.sg`；
4. 若该 ordinal 在任务 2 的 fixture 集内，则与官方 `SHOW ELEMENTS` 逐操作相等；
5. `IsotropySubgroup` 的 `parent_sg/ordinal/record` 上下文自洽（字段是 public，
   伪造记录必须被拒绝，不能成为绕过验证的输入通道）。

`SHOW ELEMENTS` 打印的是**母群帧**中的陪集代表，且**打印什么符号取决于子群是否
在母群帧中命名得出来**：旋转用母群轴命名（例 A 的 `C4x±`、`SGx`、`S4x±`，
例 B 的 `C2a`、`SGda` 这类立方 ⟨110⟩ 对角名），平移在母群 conventional 坐标，
而且**不约化**（例 A/B 全为 `(0,0,0)`；#139 `M1-` 打印的 Origin 是 `(1,1,1)`、
代表元含 `(C2x|0,2,2)`、`(I|5/2,5/2,5/2)`。在该例的 `L_H = Z³` 下，
这两条平移分别归约为零与 `(1/2,1/2,1/2)`；在 I 心母群格 `L_G` 下两者都归约为零）。
因此第 4 条的比较对象是"子群自身 setting 的操作经 `(T,o)` 搬入母群
帧后的操作多重集"，**按 `L_H` 约化后比较**（子群元素身份由 `L_H` 决定；SG 139 的
反演因此必须显示 `(1/2,1/2,1/2)` 而不是 `0`），**不是**把打印符号再套一次
`R_G = T R_H T^-1`（§2 陷阱）。任务 2 的 fixture 必须同时钉住打印串与它对应的
`(R_G, t_G)`（符号→矩阵表随 fixture 冻结），并把这条写成断言而不是注释。
fixture 比较器 `tests/irrep_subduction.rs::keys` 即按 `subgroup_lattice()` 约化。

## 5. 表示语义：Γ / selected-arm / full-star，以及维数来源

| 名称 | 含义 | 本引擎中的维数来源 |
|---|---|---|
| Γ | `k = 0`；小群 = 全点群 | `IrrepRecord::dim`（物理维数）、`CharacterRow::dimension()`（该行空间的维数） |
| selected-arm | 只取星的一条臂，在臂的小群上解释 | `RepresentationSpaceKind::SelectedArmBlockTrace`；行仍按**完整 PIR 操作宇宙**索引，调用者必须自己限制到臂小群集合 |
| full-star | 整个星（诱导表示） | 诱导表示维数 = `little_dim × star_size`；普通标量的求值、折叠和独立分解入口见 §10–11 |

**首版（任务 5）只做 Γ、标量、普通（非 compound、非 spinor）请求**；其余输入返回
`Unsupported`/`MissingIrrepData`，不返回部分项。`SelectedArmBlockTrace` 的行不是
full-star 行：任务 7 的 selected-arm 结果必须用**独立名称**，任务 8 不得用
"full-star trace ÷ 臂数"伪造它。

## 6. 复/实语义与 compound 计数

- 首版分解的是**复不可约成分**（复化后的限制表示）；`IrrepRecord::dim` 仍是物理
  （可能为实）表示的维数，两者不得互相冒充。物理实表示的重新分组另行展示。
- `CompoundMetadata::semantics` 决定行的组装方式，**必须**分别处理：
  - `ConjugateRealification`：`χ = 2·Re(χ_CIR)`，该行范数为 2，**不能**喂给复不可约
    特征标正交性内积；
  - `DistinctComponentSum`：`χ = Σ χ_CIR`，按 constituent 分别分解。
- **禁止**用 ML 标签的拼接长度推断 compound 的重复计数；目标项必须能追溯到
  `CompoundMetadata::cir_irnumbers`/`cir_labels` 这样的**冻结 CIR 来源身份**。
- 目标候选行必须是完整、互异的复不可约集合；先用 Gram 矩阵验证正交性（任务 6），
  再套内积。
- **不得用 `block_trace` 反过来"验证"它自己的成分之和**：`compound_selected_arm_view()`
  的 `block_trace` 就是由同一对 constituent 组装的，比较它是循环论证。
  已删除该 gate 与只为它存在的 `InconsistentCompoundRow`。独立证据在
  `tests/compound_subduction_regressions.rs`：SG 83 的轴向/极向 2D 迹与 SG 23
  `W1W1` 的平移本征值用物理模型（`Rxx+Ryy`、`det(R)·(Rxx+Ryy)`、Bloch 相位）钉住，
  不经过 CIR 行求和。
- **拉回到子群帧必须保留平移**：`representatives()` 的母群帧代表元经 `unmap_operation`
  映回子群 setting 时不先对子群格约化，随后**确定性地**撤销冻结的 `child_shift`
  （`t' = t + δ - Rδ` 的逆），不再"在两个 origin 里取第一个通过检查的"。
  `character_of` 已经按"旋转＋平移模子群格＋Bloch 相位修正"配对，因此保留的
  平移既决定相位、也不引入歧义。

## 7. 来源身份、严格 Hall 来源与精度

- **来源身份 ≠ 显示标签**。`ml`/`bc` 是显示标签（112 个母群 irrep 的旧拼写官方二进制
  不接受）；稳定身份来自 CIR/PIR 来源编号与冻结 provenance 串。输出条目必须带来源
  身份，不能合成标签。
- **`bridge::canonical_hall_ops()` 不是严格 API**：它有 `hall == 0 → from_sg` 与
  `from_database` 失败 → `from_sg` **两处** first-Hall 回退。新引擎用私有
  `strict_sg_hall_ops(sg)`（只读 `SG_DATA_HALL[sg]`，为 0 或加载失败即报
  `StrictHallUnavailable`），公共 bridge 行为不变。
- **精度**：`wigner.rs::ExactSeitzOp` 把平移固定为分母 12 的网格
  （`ExactSeitzError::TranslationOffTwelfthGrid`），而 pinned origin 的分母包含
  **16**、变换过程还会产生新分母（如 1/8、1/16 的组合）。因此引擎使用**本模块局部**
  的 checked 有理类型，不改全库代数：

```text
Rat { num: i128, den: i128 }   // den > 0、gcd(num,den)=1，构造时归一
checked_add / checked_sub / checked_mul / checked_div / is_integer / to_i32_checked
Mat3R / Vec3R                   // 3×3 / 3 向量；inverse3 -> Result<_, SingularTransform>
```

  所有 `f64`/`i32` → 有理的转换走 checked 构造；**不得**用 `as` 强转、不得把分母
  截到 12、不得把"取模后丢弃平移"当作归约。整数溢出必须报错。
- `KVector { numerators: [i8;3], denominator: i8 }` 是窄整数类型；任务 7 的折叠 k
  若超出其范围必须报错或另用宽表示，**不得截断**。

## 8. 实际 API 清单（签名 / 帧 / 维数 / 精度）

| API（已核对签名） | 帧 | 维数 / 精度 | 限制 |
|---|---|---|---|
| `query::irreps_of(sg: u8) -> &'static [IrrepRecord]` | ISO/data-Hall setting | `dim: u8`；`kx,ky,kz: i8`、`kd: i8` | 含 spinor 记录；compound 用拼接 ML 标签 |
| `IrrepRecord::ordinary_scalar_selected_arm_block_trace() -> Result<CharacterRow, CharacterViewError>` | 行内 `operations` 与字符同帧（data-Hall/PIR 宇宙） | `CharacterRow::dimension()`；`Complex64` | `spinor` 或 compound 返回 `NotApplicable` |
| `IrrepRecord::compound_selected_arm_view() -> Result<CompoundSelectedArmCharacter, CharacterViewError>` | 同上 | `ConjugateRealification{seed, block_trace}` / `DistinctComponentSum{first, second, block_trace}`；constituent 各自 `dimension` | 必须按 `semantics` 分支，不得把 `block_trace` 当单个复 irrep |
| `IrrepRecord::spinor_selected_arm_view() -> Result<SpinCharacterRow, CharacterViewError>` | 同上 + SU(2) lift | `SpinSeitzOperation{seitz, pauli}` | 任务 11 之后单独处理，首版不用 |
| `CharacterRow::{representation_space, dimension, len, values, get, entry, operations, operation}` | 行自带 | `Complex64` + `SeitzOperation{rotation:[i32;9], translation:[f64;3]}` | 构造函数已校验 identity 唯一/维数/有限性；平移是 f64，精确化需走 §7 的 checked 转换 |
| `irrep::isotropy::isotropy_subgroup_for_direction(sg, ml, IsotropyDirection) -> Result<IsotropySubgroup, IsotropyError>` | 记录为 pinned 数据 | `record.basis: [[i32;3];3]`、`record.origin: [i32;4]` | 字段 public；嵌入层必须重新校验（§4.5） |
| `irrep::isotropy::{parent_primitive_basis, basis_in_parent_conventional, origin_shift_in_parent_conventional, subgroup_size}` | primitive → conventional 换算 | f64 / `u32` | 只承诺帧换算与 Size；不代表官方 Basis/Origin 列（见 `docs/isotropy-data-semantics.md` §3） |
| `SymmetryOps::from_database(hall: usize) -> Result<Self, SymError>` | 指定 Hall setting | 精确整数 R + f64 t | 新引擎的严格来源，经 `SG_DATA_HALL[sg]` |
| `bridge::canonical_hall_ops(sg) -> Result<SymmetryOps, SymError>` | data-Hall，**有回退** | — | 不得当作严格 API（§7） |
| `wigner::{ExactSeitzOp, exact_seitz_table, ExactSeitzReduction}` | canonical Hall 表 | 分母 12 网格、`lattice_shift: [i32;3]` | 不能承载 1/16 origin；按 L_H 泛化的部分在 `subduction.rs` 内实现 |
| `query::k_vectors_agree(a: KVector, b: KVector) -> bool` | — | i8 窄整数 | `d = 0` 返回 false；折叠 k 的等价判定要按子群倒格做（任务 7） |

## 9. 当前不支持边界（必须报错，不得部分返回）

- 现有 `subduce_irrep_with_embedding` 的多臂分解仍返回 `UnsupportedMultiArmStar`；
  标量完整星使用 §11–12 的独立入口，包含 compound 和 realification 的 k/-k。
  尚不支持 spinor/双群（任务 11）；磁子群与磁共表示（任务 10）；
  参数化 k（不在 `query::irreps_of(sg_H)` 离散表中的折叠 k）；无法唯一确定嵌入的
  ordinal（任务 4 冻结路径之外的）。
- 非 Γ（`k ≠ 0`）与 compound 的成分展开（任务 7）已落地；已实测覆盖：四个采样
  (母群, 子群) 对（16 `R1` P1、221 `GM4+` P1、221 `GM4+` P2、139 `M1-` P1）的
  全部标量探针（含 Γ）共
  folded 92 / 多臂 57 / 缺数据 0，其中 SG 139 该对贡献 folded 20（全部 Γ/M 单臂）
  与 17 个真多臂；扩充 fixture 后，Γ 侧清点为 1895 条记录，243 条命中已冻结的
  子群号，其中 14 条钉住、199 条多候选歧义、30 条 setting 在搜索空间外。
  历史 CLAUDE 里的 160/164 不出自任何当前被 pinned
  的检查，不要引用；独立的 stored-frequency gate（59 条记录 / 538 次比较）也全绿。
- 任务 6 的 SG 23 `W1W1` realification 测试使用**显式的四元平移商** 0, L, 2L, 3L
  （`L = (1/2,1/2,1/2)`，`k = (1/2,1/2,1/2)`）上验证 seed 与共轭成分的展开
  （seed 特征标 `[1,-i,-1,+i]`、共轭 `[1,+i,-1,-i]`，内积 (1,0)/(0,1)）；
  该测试本身不验证 k/-k 的 full-star 分组；任务 8c 的独立验收见 §12。
- 对应错误：`UnsupportedMultiArmStar`、`MissingIrrepData`、`AmbiguousEmbedding`、
  `StrictHallUnavailable`、`SingularTransform`、`NonIntegralMultiplicity`、
  `CharacterMismatch`、`DimensionSumMismatch`（任务 3/5 落地时确定到具体变体）。
- 任何"返回空分解冒充成功"、"取候选第一项"、"放宽误差/加豁免"、"用 legacy
  `characters()`/`matrices()` 与 Hall 操作配对"的做法都在本契约下不合格。

## 10. 任务 8a：普通标量完整星求值与折叠几何

暂存 API 位于 `irrep::subduction::star`。`OrdinaryStar::new(probe)` 只接受生成表
内的普通标量记录；compound 与 spinor 显式拒绝。它以严格 data-Hall 操作的完整
Seitz 元素 `g_i` 为 transporter，按母群倒格枚举 `k_i = R_i^-T k`。对任意母群
操作 `h`，完整特征标为

```text
χ_full(h) = Σ_{i: h k_i ≡ k_i (mod L_G*)} χ_seed(g_i^-1 h g_i).
```

共轭操作的平移保留到 selected-arm 行求值结束；不固定的臂贡献零。任意输入操作
须通过母群成员检查。显式 transporter 构造也必须覆盖完整星，不能靠缺失臂得到
较小但看似成功的表示。恒等特征标同时核对 `selected_dim × arm_count` 与源记录
维数。全 4105 条普通标量记录通过的是**构造、星完整性和恒等维数**检查，不能写成
全部操作的 full-star 分解已经验证。

`folded_stars(embedding)` 将每个臂按 `T^T` 折叠，按子群倒格合并同点并保留所有
母群臂索引，再按子群旋转分成 stars。每个 star 的不同 q 点须携带相同臂数；
`FoldedStar::block_dimension()` 是**整个该子群 star 所承载的母群子空间维数**，
不是子群某个小表示的维数。此阶段没有算子群 irrep 重数。

独立字符证据来自 `scripts/generate_subduction_star_fixtures.py`：校验 pinned
`CIR_data.zip` SHA256 后，直接对原始完整矩阵及各臂对角块取迹，生成
`tests/data/subduction_star_cir.rs`。10 个固定 CIR 记录（SG 92/139/198/221/225）
覆盖 308 个原始操作，含螺旋和中心化胞；Rust 测试先核对实际操作与星臂所在帧，
再核对原始及格平移后的字符，共 1232 次比较。平移后的独立期望为各个原始对角块
分别乘 `exp(+2πi k_i·L)` 后相加，绝不为整个 full-star trace 乘一个统一相位。

## 11. 任务 8b：普通标量完整星分解

暂存入口为
`irrep::subduction::star::decompose::subduce_full_star_with_embedding(&subgroup, &embedding, probe)`。
它重新验证缓存嵌入与输入记录的上下文，返回 `FullStarSubduction`。任务 8b 接受
普通标量母群；任务 8c 扩展到 compound（§12）。原单臂入口的支持范围不变。

每个输出 `FullStarBlock` 对应一个子群 k-star，保留精确的折叠代表点 `q()`、
所匹配数据行的 `stored_k()`、`star_size()` 和非零目标项。每项目标的 `dimension`
是**复小表示维数**；`irnumber` 是实际 CIR 来源号，compound 子群行
`DistinctComponentSum` 的两个成分分别报告。`little_dimension()` 是该 q 的总表示
空间维数，`block_dimension()` 是整个子群星承载的维数，两者相差 `star_size()`。

算法在整个折叠子群星里寻找有数据的代表臂，再在 `H_q` 上仅对折叠到该 q 的母群
臂计算字符。`H_q` 可以置换这些母群臂；被置换的臂在迹中贡献零，不能把完整迹
除以臂数。子群操作取模 `L_H` 的陪集代表，拉回 Hall 帧时完整保留平移及冻结
原点修正。共享字符求解器检查 Gram、整数重数、q 块维数和逐操作重建；随后检查

```text
Σ_targets multiplicity × child_little_dim × child_star_size = parent_full_dimension.
```

最后从**子群自身**的 Hall 操作和源行重新诱导各目标的完整星字符，在全部子群
陪集代表上与母群完整星限制比较。`reconstruction()` 同时返回这两组字符；这一步
不复用母群 q 块作为所谓的目标重建。

独立回归 `tests/subduction_star_decomposition.rs` 固定了 13 个完整分解，期望来自
归档 CIR 各臂对角块的单独内积计算；又直接用原始矩阵迹核对 92 个嵌入操作上的
完整重建。它覆盖 221 → #83/#12、139 → #126、225 → #8 的不同方向，包含重数 2、
二维小表示、Size=2、非零子群原点修正及 compound 子群成分。例如：

- 221 `GM4+` P1，探针 `X1+`：#83 的 `X1+`（star size 2）与 `Z1+`（size 1），各一次。
- 221 `GM4+` P2，探针 `X5+`：#12 的 `A1+`、`A2+` 各一次，`V1+` 两次；`2×1×2 + 1 + 1 = 6`。
- 139 `M1-` P1，探针 `N1+`：#126 的二维 `R1` 一次，star size 2，完整维数 4。

`tests/subduction_identity_regressions.rs` 还独立核对存储的恒等项频率：59 个冻结
嵌入下，2060 次普通探针对照（458 正项、1602 零项），其中 1522 次为非 Γ 探针。
缺数据的 15 个组合也作为完整清单钉住：ordinal 13345/13346/13351 的 `W1`–`W5`。
它们未计入成功分解；新增缺失会使回归失败。

缺少任一子群星的数据返回 `MissingChildStarData`，不返回部分分解。spinor、磁共表示
仍不支持。逐记录嵌入覆盖由任务 9 处理，不能把本阶段样例通过写成全表验收。

## 12. 任务 8c：compound 与 realification 完整星

`irrep::subduction::star::scalar_star::ScalarStar` 接受生成表内的标量记录；原
`OrdinaryStar` 保持普通记录限定。`DistinctComponentSum` 保留两个 CIR 来源；
`ConjugateRealification` 保留 seed 在 k 的完整星及其共轭在 -k 的完整星。
共轭作用于**整个字符求值，包括 Bloch 相位**。两个成分即使等价，母群也必须保留
两份；重合的星臂在折叠时携带成分身份，不因 q 相同而丢失维数。
`ScalarStar`、`FoldedStar`、`FullStarBlock` 的 `arm_count()` 都按成分计数；
子群不同 q 点的数量是 `star_size()`，两者不能混用。

子群目标目录枚举复不可约成分，允许 -k 没有独立数据行的情况。若两成分在同一
子群星上，先用完整 Seitz 操作把行搬到共同 q，再比较 `H_q` 上的字符。只对同一
realification 的 seed/共轭配对作等价检查；单位范数且逐操作相同才合并为一个目标，
正交则分别保留，其余情形报错。不同 CIR 来源不因标签或 k-star 相同而合并。
永久反例还固定了数值边界：逐项约 `1e-4` 的相位误差可能只造成约 `2.5e-9` 的
内积误差，因此“内积距 1 小于 `1e-7`”不能代替逐操作比较；Gram 对角元也使用
共享求解器相同的阈值，不对它开平方后再比较。

返回值的 `parent_source_identity()` 保存完整来源身份；`parent_irnumber()` 对普通
记录返回 `Some(id)`，对 compound 返回 `None`。目标同时保留真实 CIR 号和
`SubductionComponent`，区分 seed 与共轭。`FullStarBlock::stored_k()` 现在是所选
复成分的**有效波矢**，可为源行 k 的负值；它与精确折叠代表点 `q()` 同子群倒格类。

源数据审计 `scripts/generate_subduction_compound_fixtures.py` 校验 CIR archive 的
SHA256，并按冻结元数据中的来源号核对全部 672 条 compound 的臂表达式：519 条
distinct、112 条 k/-k 异星 realification、41 条同星 realification。两个 distinct
来源的完整臂表达式相同，所有 compound 成分的小表示维数相等。这是来源及几何审计，
不是全表非恒等字符验收。41 条同星来源的原始完整矩阵字符另在每个操作及有限平移
相位类上与其共轭比较，共 1484 次；SG 23 的异星负例保证此检查包含中心化平移。

`tests/subduction_compound_stars.rs` 从原始 CIR 矩阵独立核对 17 个来源、12 个物理
compound，含每个来源臂的格平移相位。分别检查各复成分很关键：SG 23 在 I 心平移
处 seed 为 -i、共轭为 +i，二者总和为零，单看总和无法发现相位符号交换。
完整分解钉值包括：

- #19 `R1R1` 自限制：CIR 559，复维数 2，重数 2。
- #23 `W1W1` 自限制：CIR 716 的 seed 与共轭分属两个星，各重数 1。
- #45 `W1W1` / `W2W2` 自限制：分别为 CIR 1805 / 1806，小维数 1、星大小 2、重数 2。
- #83 `GM3+GM4+` 自限制：CIR 4077、4078 各重数 1。
- 167 `GM3+` P1 → #15，探针 `T1T2`：二维 `M1`（CIR 363）重数 2。

上述源矩阵验收共有 928 次复成分字符、464 次物理和、28 次分解重建比较。
四个自身嵌入还遍历全部标量探针，验证自限制保留每个复来源：55 个普通、23 个
compound 成功，79 个 spinor 不支持且未计入成功。加上既有五个上下文的
161 个普通和 1 个 compound，共 240 次完整分解；这些上下文均无缺失子群星。
另在 69 条冻结嵌入记录中找到 78 次 compound 恒等项对照，全部是零项，其中 64 次
为非 Γ 探针，无缺失；它不提供正恒等项覆盖，也不代替全表 94271 条验收。

新增四个自身嵌入来自实际 `SHOW ELEMENTS`，全部 basis=I、origin=0、Size=1；操作
oracle 扩为 14 用例 / 90 代表。没有用“母群等于子群”绕过 setting 校验。全表逐记录
setting、离散表缺失的子群 k、磁群与 spinor 仍按各自后续任务处理。

## 13. 任务 9：恒等内容（identity-only）精确路径

`subduce_full_star_with_embedding` 的契约不变：任何一个折叠子群星缺少随包离散
k/irrep 数据就报 `MissingChildStarData`，**不部分返回**，并且仍然要求全部块的
维数覆盖整个母群星（`TotalDimensionMismatch`）。恒等分导表只需要一个更弱但同样
精确的量，因此任务 9 增加独立入口：

```rust
pub struct TrivialContent {
    pub total: u32,        // 按子群恒等行的冻结 CIR 来源号累计的重数（存储表口径）
    pub by_label: u32,     // 按子群 ML 标签累计的同一和；两者不等即 TargetSourceMismatch
    pub gamma_stars: usize,
    pub skipped_stars: usize,
}
pub fn trivial_content_with_embedding(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    probe: &'static IrrepRecord,
) -> Result<TrivialContent, FullStarError>;
```

**为什么跳过非 Γ 块是精确的**：子群格平移 `t` 在折叠波矢 `q` 的任何表示上作用为
标量 `exp(-2πi q·t)`，母群星在对应陪集上的特征标因此是该相位乘以 `q = 0` 处的值；
子群恒等表示在每个平移上作用为 1。所以两者共有不可约成分要求
`exp(-2πi q·t) = 1` 对所有子群格平移 `t` 成立，即 `q = 0` 模**子群**倒格（含
centering 消光，判定用 `Lattice::contains`）。落在其它 `q` 的折叠星对恒等重数的
贡献恒为 0，不需要任何子群数据，跳过它们不是近似。Γ 块本身仍由完整分解同一个
`build_block` 阶段构造：同样的字符配对、同样的冻结 child origin shift。

契约细节：

- 上下文（subgroup/embedding/probe）先用与完整入口相同的
  `validate_subduction_context` 重新校验；spinor probe 同样以
  `UnsupportedCharacterSpace` 拒绝。
- 子群必须存在唯一的“恒等 Γ 行”（一维、k = 0、在其自身存储 setting 下全为 +1）；
  该行由 `trivial_child_record` 在随包表里查找，找不到或有多个时报
  `MissingChildTrivialIrrep`，绝不返回 0。`every_space_group_has_one_trivial_gamma_row`
  对 230 个空间群逐一钉住该行存在。
- `by_label` 与 `total` 必须相等，否则报 `TargetSourceMismatch`（不是静默取其一）。
- 返回值只对**恒等重数**成立；它不给出完整分导分解，调用方不得把它当作
  `FullStarSubduction` 使用。

永久回归（`tests/subduction_identity_regressions.rs`）：

- `identity_only_content_answers_probes_without_full_child_data`：SG 196 `W1` P2 → #24，
  探针 `W1`（k = (1/2,1,0)）的完整分解由 `MissingChildStarData` 拒绝，恒等-only
  给出 `(total, by_label, gamma_stars, skipped_stars) = (1, 1, 1, 2)`，等于存储表
  频率 1；同几何的 `W2` 未被存储表列出，恒等-only 给出 0；`L1` 折叠出的星没有 Γ
  点，恒等-only 与完整分解都为 0，且 `gamma_stars + skipped_stars` 等于折叠星总数。
- `identity_only_content_agrees_with_the_full_decomposition_and_covers_the_pinned_set`：
  在任务 8 的 59 个冻结上下文里，2,060 个 probe 上恒等-only 与完整验证过的完整星
  分解**逐条相等**，另外 15 个（ordinal 13345/13346/13351 的 `W1`–`W5`）完整分解
  无数据、由恒等-only 精确回答并与存储表一致。
- `src/irrep/subduction_star_decompose.rs` 单元测试另钉住 SG 221 `GM4+` P1 → #83 的
  全部标量 probe 上两条路径一致，以及 SG 196 的上述四元组。

全表审计（`examples/audit_irrep_subduction.rs`）据此把探针结果分为
`full_success` 与 `identity_only` 两类，二者都与存储表逐条比较（历史：第六轮时
14,713 条恒等-only 中 160 条是存储正项，全部复现，0 不匹配、0 假阳性；R4 三个批次
把 `identity_only` 清零后，这条恒等-only 路径仍保留在审计里，只是不再有 probe
需要它，`probe_identity_only = 0`）。范围闭合的判词由
`--require-complete` 给出；`other_wave_vector_subduction` 的 5,756 行因 pinned
irrep 表没有其 73 个源 irrep 的 k 矢量与特征标行，单独由 `--require-w-complete`
与 `w_scope` 行报告，未计入该范围。它们的**小群表**在归档 `data_little.txt` 中完整
存在（73/73 源、维数相符），而且**解码后可见它们全部是参数化直线**：
`little_k` 每个 (Bravais 格, k 槽) 存 4 组 `(x,y,z,d)` = 基点 + 至多三个自由方向，
73/73 个源的基点都是 Γ、自由参数恰为一个（cF 的 `(1,0,1)`/`(1,1,2)`，cI 的
`(1,-1,1)`/`(0,0,1)` 等）。因此不存在单一数值 k 可以去折叠，剩下的工作是解码这些
直线上的小群特征标，再按线（而不是点）对照 5,756 个存储频率。
