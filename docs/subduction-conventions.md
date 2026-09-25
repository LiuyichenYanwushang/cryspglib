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

## 14. 存储行与构造目标的两条契约（对抗性复核 B，2026-09-23）

**存储行只在记录自己的小群上有定义。** 归档字符行是「该 irrep 在其 k 点的小群
（及其代表元）上的取值」，表里列出的**其它**操作可以是占位值：SG 38/40 有 4 条
pinned 行对小群之外的 `-I`、`m_y` 存的是字面 0。任何代码都不得把这种 0 当特征标
读；生产路径只在 `H_q` 内（或经完整星诱导出的臂上）取值，所以当前无影响，但
「row = 每个列出操作的特征标」**不是**契约。要用行外的值，必须先证明该行确实
给该操作定义了值。

**构造目标没有归档来源，因此它的独立证据只有一条。** 构造目标（`irnumber = None`）
不参与「来源身份」与「来源号/标签一致」两项审计检查（对它是恒真的），而子群 Γ 星
永远命中 stored 行（每个空间群都有平凡 Γ 行），于是恒等频率比较**从不**看到构造星。
所以构造路径的证据是：catalogue 与 pinned 离散字符的逐操作对照
（`the_catalogue_reproduces_pinned_little_group_characters`，计数钉死为
2,661 条 pinned 行 / 22,302 个操作 / 719 条走高维目录 / 293 条源行维数为 2
（其中 121 条 D4）/ 1,282 条延后）、各族结构门禁与正交门禁、以及引擎自身的
Gram/维数/重建检查。这些事实写在 R5 报告的独立性表里。

**fail-closed 的分层与它们的测试。** 判据链是
`has_trivial_little_co_group`（纯几何）→ 一维 solver 只在解数 `== |P_q|` 时返回
→ C2×C2、D3、D4 三族结构门禁 + 正交门禁 → 其余一律空表 → `select_representative` 报
`MissingChildStarData`。边界由以下正反回归钉住：

- `the_one_dimensional_solver_returns_nothing_instead_of_a_subset`：非上边界 cocycle
  返回空（不是子集）、`MAX_ORDER`（48 阶 Γ 点）返回空；
- `generator_order_search_handles_large_cocycle_denominators`：分母为 512 的合成 coboundary
  返回完整四个解，证明搜索成本由小余群生成元阶决定，不随相位分母增长；
- `an_out_of_scope_co_group_still_reports_missing_child_star_data`：在真实嵌入
  （225 `X1+` P3 → #221）上把折叠星放到越界点 (0,0,1/2)，`build_block` 必须报
  `MissingChildStarData`——覆盖闭合后 pinned 语料里已没有能触发该分支的 probe，
  所以这条负例是手工构造的，且**必须**保留；
- `constructed_targets_have_their_own_identity_and_no_borrowed_labels`：构造目标的
  身份只由精确点与构造小群给出，不借用任何 stored 标签。

`ConstructedStar::dimension()` 现在携带真实的小群维数（一维目标为 1，二维目标为 2），
不再硬编码 1；`stored_child_components_at` 遇到无法展开的子群记录**大声报错**而不
`continue`——把生成 bug 伪装成「缺数据」比报错更糟。这两条都由复核 B 指出。

## 15. 嵌入验证与结果元信息的两条契约（第三方复核，2026-09-23）

**嵌入的必要条件是 `L_H ⊆ L_G`。** `validate_candidate` 现在先检查候选映射推出的
子群格 `U^-1 · W · P_parent` 的每一行都落在母群格 `Lattice::new(P_parent)` 里，
再检查有限操作映射、代表元数目与模 `L_H` 的闭合性。原因：`det U = ±1` **不**蕴含
包含关系——SG 1 自嵌入取 `U = diag(2, 1/2, 1)`（即 `diag(4,1,2)/2`）时，
`det U = 1`、唯一有限操作 `(I|0)` 映射到自身、模 `L_H` 只有一个陪集，所有有限检查
都通过，但它把子群平移 `(1,0,0)` 送到 `(1/2,0,0) ∉ Z³`，根本不是空间群嵌入；
修复前 `probe_embedding` 会把它当 `Ok` 返回（已复现）。搜索路径只枚举幺模候选，
天然满足该条件；这条检查真正守护的是**分数 `U`** 的候选与冻结行。
回归：`a_setting_whose_lattice_is_not_a_parent_sublattice_is_rejected`（负例）
与 `every_frozen_setting_keeps_the_subgroup_lattice_inside_the_parent`（全表
15,239 行逐一重算包含关系，并断言分数行仍 ≥ 30 条）。

**结果对象必须带 setting 的分母。** `U = setting / setting_denominator` 是有理矩阵，
冻结表里有 30 条分母为 2、3 条分母为 -2（最简写法下）的记录（ordinal 26 是第一个
分母 2 的行：`U = [[1,2,1],[-1,2,-1],[-1,0,1]] / 2`）。`SubgroupEmbedding` 一直同时
暴露分子与分母，但 `IrrepSubduction` 与 `FullStarSubduction` 之前只复制分子，调用方
拿到 `setting()` 会按整数矩阵做坐标变换。两个结果类型现在都有
`setting_denominator()`；回归 `a_fractional_setting_reaches_the_results_with_its_denominator`
在 ordinal 26 上钉住 `setting()/setting_denominator()` 与嵌入一致（并说明当前语料里
没有 Γ 归属的分数行，所以 `IrrepSubduction` 一侧只钉接线 + 幺模 Γ 上下文）。

**空范围不是通过。** 审计的 `--parent`/`--ordinal` 组合若选不中任何非 spinor
isotropy 记录（如 `--parent 2 --ordinal 0`），现在直接报 `empty scope` 错误并非零
退出，不再打印 `VERDICT clean` 与 `probe_full_success=0/0`（修复前即如此，已复现）。
回归：`an_empty_scope_is_rejected_instead_of_reporting_clean`，同时确认
`--parent 1 --ordinal 0` 与 `--parent 1` 仍然正常。

**生成管线的门禁与运行时无关，但同样按"失败关闭"要求：** 全表只有一个 writer
（`scripts/task9/build_table.py`），它必须先确定"预期 ordinal 全集"才允许写：
`--expected` 显式列表优先，否则已存在的 `--out` 与 tracked 的 `--baseline` 模块
（默认 `src/irrep/subduction_settings_data.rs`）必须一致，两者都不可用时拒绝写入
（除非显式 `--partial`）——否则"输出路径打错/在新目录重建"会把"分不清完整与截断"
变成一次成功写入。旧生成器 `scripts/generate_subduction_settings.py` 的 `--check`
与离线测试共用 `check_ordinal_coverage`：ordinal 序列必须严格等于
`range(len(records))`，只比条数会被重复项顶替缺失项。

**工具不再吞错误。** `examples/probe_subduction_settings.rs` 之前把
`subduce_full_star_with_embedding` 的任何错误写成 `trivial=0`，而 0 在这条流水线里
是**合法结果**（表示约定落在共轭分支上），两者混淆会把引擎缺数据误判成物理结论。
现在三态分开：`trivial=<n>` 是算出的重数、`trivial=?` 是结构上不可用、`trivial=error`
是计算失败（原因打到 stderr），`--profile` 同理输出 `profile=error` 而不是一个
"少了若干项、可能恰好匹配目标"的更短 profile。


## 16. R6.0/R6.1 契约：参数化 k 的完整分解（显式参数值）

**能力 A（本轮交付）**：给定 `(isotropy 记录, 母群参数化源, 有理参数 t)`，算**这一个
参数点**的完整分解；不承诺任何参数区间上的结论。能力 B（整个参数族的覆盖说明）留到
R6.2，见 [subduction-r6-plan.md](subduction-r6-plan.md) §1。

```rust
pub const OFFICIAL_LINE_PARAMETER: (i128, i128) = (1, 4);
pub fn official_line_parameter() -> Result<Rat, SubductionError>;
pub fn subduce_line_at_parameter(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    table: &'static LittleCharacterTable,
    parameter: Rat,
) -> Result<LineSubduction, FullStarError>;
```

输入/输出契约：

* `table` 必须属于 `subgroup.parent_sg`（否则 `LineSourceMismatch`）；记录与嵌入必须
  互相一致（与离散入口同一套 `StaleIsotropyRecord`/`EmbeddingContextMismatch` 校验）。
* `parameter = t`：波矢是 `t · table.direction`，方向用的是**冻结表自己存的**那份
  （官方打印的 conventional 方向，例如 SG 196 `DT` 的 `("0","2","0")`）；臂、字符与
  折叠共用这一个向量，帧不再做二次换算。`t = 1/4` 是官方参数约定
  （`OFFICIAL_LINE_PARAMETER`），不是引擎选择。
* **帧：整套量都在母群 conventional 倒格坐标下**（实测钉死，勿按"primitive 折算"实现）：
  `direction`、pinned `little_k`（SG 196 `X1 = (0,1,0)`、`L1 = (1,1,1)/2`）与冻结
  little 群操作都是 conventional 帧；`exact_primitive_basis` 只提供**格子**，不改变
  坐标系；`fold_wave_vector` 直接作用在 `t·direction` 上，没有任何二次换算。
* **波矢不约化 + monodromy 契约（R6.2 最终读法；R6.1 的 `canonical_wave_vector`
  已撤销）**：冻结的 `D` 是 **`k = Γ` 处解出**的纯 Γ 点字符，参数通过 Bloch 因子
  进入：`χ_α^t(R, T) = D_α(R) · exp(2πi t (v·T))`。所以引擎必须用**原始**
  `k(t) = t·v`，把 `k` 约化进母群基本胞而保留 `D` 等于去算**另一条带**。
  一个参数步 `t → t + n` 乘上 reciprocal-shift twist
  `Φ_{nv}(R) = exp(2πi n (v·T_R))`，它在线小群上是**真正的一维特征标**，把标签送到

  ```text
  (k + K, α) ~ (k, M_K(α))，  decompose(α, t + n) == decompose(M_{n v}(α), t)（逐块逐目标）
  ```

  `M_K` 由**字符指纹**算出（绝不按标签名硬编码），实现在
  [`crate::irrep::line_monodromy`]（`monodromy` / `complex_conjugation` / `orbit`）；相位
  以有理数模 1 计算，再与冻结 Gaussian 整数特征标精确比较，不走浮点容差。当前
  冻结值是 `(实部, 虚部)` 整数对，因此非零字符项只直接支持 Gaussian 单位相位（四分之一圈）；
  其他有效根相位若不能匹配冻结表，会返回 `Missing`，表示当前冻结格式里没有可用像。
  沿源自身方向的非恒等映射为：SG 203/210 交换 `DT1↔DT2`、`DT3↔DT4`；SG 227/228
  交换 `DT1↔DT3`、`DT2↔DT4`；SG 209 的沿线映射为恒等，`SM` 源沿自身方向也不移动。
  SG 209/210 的 `DT3`/`DT4` 仍互为复共轭，但复共轭是另一种映射，不能与沿线位移混为一谈。
  非恒等沿线映射在 pinned 表上产生 **816** 条可比较的频率对；参数输运中同标签读法
  失效的 **40 行**由 `a_parameter_step_is_not_the_identity_on_a_nontrivial_monodromy`
  逐条见证。R6.1 把 SG 210 的 `Z1@1/4 → Z2@5/4` 读成 gauge 滑移，实际是正确位移：
  该行的 `M_v(DT3) = DT4`，而 pinned `DT4@1/4` 的分解正是 `Z2`。
  `monodromy(parent, K)` 是通用的倒格字符扭曲映射；只有当 `K = delta * v` 与某个源的
  冻结方向 `v` 平行时，它才表示该源参数 `t -> t + delta` 的沿线位移。SG 203 的
  `K=(2,0,0)` 见证是跨线扭曲示例，不是参数步进。函数先精确验证 `K ∈ G*_parent`；
  对每个线源再逐对检查小群
  操作相位是否乘法：对 `g=(R_g,t_g)`、`h=(R_h,t_h)`，必须有
  `((I−R_g^T)K)·t_h ∈ Z`。这是 `exp(2πi K·t_g)` 成为该线小群一维特征标的精确条件；
  `R_g^{-T}K=K` 只是更强的充分条件。不能通过精确检查时返回
  `LabelImage::UnsupportedShift`，不猜相位公式的像。对每个冻结源，令
  `d = gcd(coords_{G*}(v))`，则 `v/d` 是沿该线的
  最小倒格平移；全表 73/73 源的像唯一，且完整分解满足
  `decompose(α, 1/4 + 1/d) = decompose(M_{v/d}(α), 1/4)`；实测这 73 个源都 `d=1`，
  所以没有扩展参数集合，确认的仍是 `t = 1/4 + n`。
* **已有官方锚点且经输运验证的域**：逐源为 `t = 1/4 + n/d`，其中 `v/d` 是该母群
  倒格上的沿线最小正步长。它满足 `(t − 1/4)·v ∈ G*_parent`；完整分解等于相应
  monodromy 像在 anchor 的分解。pinned 频率在 monodromy 轨道上不变（816 条非平凡像
  全表钉住）。除此之外的退化或非平移等价参数仍没有同等级 oracle。`t = 3/4` 一类
  （共轭 coset 再走半步）若不满足该源的倒格等价条件，就既不等于 pinned，也不由
  该表的标签描述，引擎仍会给出一致的字符（hard failure 0）但**没有 oracle**，
  报告只给实测值，不作结论。
* 输出 `LineSubduction`：`parent_dimension = little_dim × arms`、每个折叠子群星一个
  `FullStarBlock`（`q`、`stored_k`、`star_size`、`arm_count`、`little_dimension`、
  `block_dimension`、目标表）、`reconstruction()`、`setting()/setting_denominator()`、
  `parameter()`、`wave_vector()`，以及 `trivial_content()`——恒等重数的两个读数
  （冻结 CIR 来源号与行标签）必须一致，否则 `TargetSourceMismatch`；子群没有唯一
  恒等 Γ 行时 `MissingChildTrivialIrrep`，**绝不以 0 代替**。这两个变体是**表损坏
  防御分支**：两个读数同源（同一 `ChildComponent` 的 `ml` 与 `irnumber`），唯一恒等
  Γ 行又对 230 个 SG 钉死过，所以公网 API 在 pinned 数据上**不可达**，也没有负例；
  它们与下面两条可达的失败语义不是同一等级，不要并列成"都已验证"。
* 不变式（引擎内强制，失败即 `Err`）：`Σ mult × dim × star = parent_dimension`
  （`TotalDimensionMismatch`）、逐子群代表元的完整星重构（`ReconstructionMismatch`）、
  每个 q 块的 `χ(E) = block_dimension`（`QBlockIdentityMismatch`）。
* 失败语义与离散路径相同：折叠点子群表里没有、且小余群不在已构造的三族内 →
  `MissingChildStarData`。ordinal 13543 的 D4 缺口已由正向端到端测试覆盖；真正的越界
  负例仍由 order-16 的手工构造星 `an_out_of_scope_co_group_still_reports_missing_child_star_data`
  钉住，不能返回部分分解或零。

**R6.1 修掉的 R5 遗留结构错误（Γ-only 路径看不见）**：旧的 `line_folded_stars` 把每个
约化 `q` 各当作一个子群星，而 `FoldedStar` 的语义是**子群点群下的轨道**。Γ-only 路径
只建 Γ 块、从不做重构，所以这个错误一直被掩盖；第一次在一般参数上做完整分解时，重构
在恒等元上给出 12（或 8）而不是 6，直接失败。现在改为复用离散路径的共享
`fold_arms`（轨道划分 + 每个轨道内臂数一致性检查），线源与离散源的折叠几何彻底统一；
R5 的 Γ-only 入口保留为 [`line_trivial_content_via_blocks`]，但审计的 w 门禁已经改用
**完整分解**（见下）。

**R6.1 验收（本轮实测）**：

| 组 | 内容 | 结果 |
|---|---|---|
| 已知点 | `t = 1/4` 下**全部 5,756 条 pinned w 行**的完整分解，恒等重数 == pinned 频率 | 审计 `w_computed=5756/5756`、`mismatched=0`、`engine_errors=0`、`hard_failures=0`（R6.1 升级当轮 604.8 s、规范波矢修复后复跑 557.0 s，三门口禁 exit 0）；另有 SG 196 的 106 行作为单元测试常驻 |
| 一般位置 | child #1（`10038` `W1` `4D1`）`t = 1/7`：6 个构造块、每块 1 维、恒等重数 0 | 与**独立于多重度求解器的几何计数**一致（0 条臂折到子群 Γ；臂集合/帧/Γ 判定与引擎共用，所以这是第二读数而非外部 oracle） |
| 特殊值两侧 | 同一记录 `t = 1/4`（2 块：`Z1`×2 + `GM1`×4、恒等重数 4 == pinned）对 `t = 1/6`、`t = 1/3`（各 6 块、恒等重数 0） | 两侧都完整、都等于几何计数（这一对的 `0 == 0` 只是弱断言，非零锚点是 `t = 1/4` 的 pinned 值）；证明"参数变了结论就变" |
| 约定（全表） | `t = 5/4` 的完整分解 == `M_v(α)` 在 `t = 1/4` 的完整分解 | **5,756/5,756**，`comparison_skipped=0`（非零时硬失败），同标签读法在 40 行上不同（审计 `w_parameter_shift` 硬门禁 + `tests/line_monodromy.rs`） |
| 约定（字符层） | `χ(5/4, α) == χ(1/4, M_v(α))` 逐代表元 | **5,756/5,756**（`line_family_coverage --gate`，容差 1e-9） |
| monodromy 代数 | `M_0 = id`、`M_-K = M_K^-1`、`M_{2K} = M_K²`、`C·M_K = M_-K·C` | 11 个母群、73 个源；M₀ 覆盖断言 73/73，逆/复合/共轭关系各 136 项；这些方向映射也包含合法跨线扭曲（`the_monodromy_contract_identities_hold`） |
| 沿线标签 | 每个源沿自身冻结方向的标签像 | 73/73 唯一且逐标签钉值；非恒等映射限于 SG 203/210/227/228 的 DT 源（`along_line_images_match_the_frozen_parent_labels`） |
| pinned 不变性 | pinned 频率在沿线 monodromy 轨道上不变 | **816** 条非平凡像逐条相等，测试精确断言 816（`the_pinned_frequencies_are_constant_on_monodromy_orbits`） |
| 失败语义 | 别家的 `table` → `LineSourceMismatch`；越界一般点 → `MissingChildStarData` | 两条负例都断言具体变体 |

**证据级别（防止把内部自洽写成外部正确）**：

| 声明 | 证据级别 |
|---|---|
| `t = 1/4` 的**恒等重数** == pinned 频率（全 5,756 行） | **外部**：pinned 行来自官方程序 `SHOW FREQUENCY`，另经 `verify_w_subduction_oracle.py` 的 live oracle 双向比较 |
| `t = 1/4` 的**完整分解**（维数守恒、逐代表元重构） | **内部一致性**：引擎自身的两个不变式；没有独立的完整分解 oracle（官方不打印载荷） |
| 一般 `t` 的分解 | **内部一致性 + 一条独立几何计数**（只替换多重度求解器）；无外部 oracle |
| 一般参数的非平凡小余群 | D4 有字符目录回归；ordinal 13543/14106/13691 三个见证在 42 个有理样本点均到达构造二维目标；另对 5,756 行全扫 42 点。其它非平凡小余群仍按已证明的 family gate 支持或 fail-closed；一般参数没有外部完整分解 oracle |
| `MissingChildStarData` / `LineSourceMismatch` | 有常驻负例（断言具体变体） |
| `MissingChildTrivialIrrep` / `TargetSourceMismatch` | **表损坏防御分支，公网 API 在 pinned 数据上不可达**（见上方失败语义条目），没有也无法写负例 |

**能力 A 的参数定义域（诚实边界，全部由 5,756 行全表清点或指定记录实测）**：

* **42 个小分母有理采样点**：旧的 54/5,756（`t = 1/7,1/6,1/3,2/7`）与 30/5,756
  （`t = 3/8`）个 `MissingChildStarData` 已通过 D4 构造族闭合。`--projective-sample-sweep`
  现扫描 `0<t<1` 中所有约分后分母 `3..12` 且不在 `(1/4)Z` 的 42 个参数；对全部 5,756
  行逐点重算，每点均为 `full=5756/5756, missing=0, other_errors=0, content_mismatches=0`。
  常驻 `sampled_d4_gap_witnesses_decompose_with_the_constructed_two_dimensional_target` 测试
  对三条见证在全部 42 点逐一要求构造二维目标、零恒等重数与逐操作重建。旧失败星是阶 8
  的 D4/C4v little co-group，缺少标准二维普通 irrep。该有限扫描不代表所有 rational `t`。
* **有限小余群的一维求解范围**：`subduction_catalogue::one_dimensional_characters`
  现在按每个生成元的有限阶枚举相位根，再逐式验证完整 projective character 方程；搜索
  不再随 cocycle 分母增长。合成分母 512 的 coboundary 回归返回完整四个解；真实见证
  ordinal 10030 `DT1`（child #18，`|P_q| = 2`）在 `t = 1/1000000` 现在完整分解成功，
  维数、重建与恒等内容均有断言。剩余边界来自表示族而非网格大小：一维 solver 只搜索
  `|P_q| ≤ 8`；阶 8 的阿贝尔 coboundary 小余群也由通用一维路径处理。高维射影目标覆盖
  经结构门禁证明的非退化 C2×C2、coboundary D3 与 coboundary D4 家族；其它小群、非
  coboundary D4 仍显式返回 `MissingChildStarData`。
  有理中间量超过 `i128` 时也会显式报算术错误。
* **`t = 0`（以及任何使 `t·v` 的稳定子严格大于冻结小群的 `t`）**：此时臂集合仍由
  **direction** 生成，引擎回答的是"由冻结 little 群表示诱导出的形式表示"，**不是**
  DT 线 irrep 在 Γ 的分导（后者不存在）。实测：`10030 DT1 t=0 → Ok(3) blocks=1`、
  `13543 DT5 t=0 → Ok(0) blocks=1`、`10038 DT1 t=0 → Ok(6)`（对照 `t = 1/4` 分别是
  `1`、`1`、`4`）。这是一个**形式值**，不是错误也不是物理结论；R6.2 必须显式分区，
  并决定是否改成显式错误。
* **计时**（同机同二进制，5,756 行顺序计时）：完整分解 `10.97 s` 对 R5 的 Γ-only
  `6.13 s` + `trivial_content()` `0.01 s` ⇒ R6.1 升级的真实代价约 **+4.8 s**；
  审计总时长的差异（517.1 / 557.0 / 604.8 s，reviewer 在并发负载下测得 697.6 s）是
  机器负载，不是本次升级的成本，不要把两组墙钟并列当成 A/B。
* **一般参数的性能长尾**：`t = 1/4` 走存储路径，全表约 11–14 s；一般参数（大分母 +
  非平凡小余群）单次调用可到 p99 ≈ 40 s、最大 53.6 s（reviewer C 实测 446 行 ×
  `t = 1/97` 合计约 724 s）。因此全表逐点扫描一般参数不可行，R6.2 的覆盖说明用
  "精确支持集 + 抽样佐证"，不假装逐点算过。

**能力 B（R6.2 参数族覆盖说明）**：全部数字由
`examples/line_family_coverage.rs --gate` 重算（判词与退出码可进 CI），正式报告见
[subduction-r6-coverage.md](subduction-r6-coverage.md)。结论分五档：

* 支持集恰好是**四分之一网格**：每个臂折到子群 Γ 的 `t` 集合是 `(1/4)Z`
  （5,756/5,756 行，精确有理枚举；103,608 个网格外 (行, 参数) 组合里折到 Γ 的臂数
  全为 0）。网格外（非退化点）`content(t) = 0` 是**证明级**结论。
* `t = 1/4` == pinned（**外部**官方输出）；`t = 5/4`（以及任何
  `(t - 1/4)·v ∈ G*_parent`）的**完整分解逐块逐目标**等于 `M_v^{n}(α)` 在 `t = 1/4`
  的分解——审计的 `w_parameter_shift` 硬门禁（5,756 行；`comparison_skipped=0`，
  同标签读法在 40 行上不同，作为"门禁非空转"的见证计数）。字符层同表核对
  （`line_family_coverage --gate`：`χ(5/4) == χ(M_v(α), 1/4)`，5,756/5,756）。
* **共轭/位移参数（`k(t) ≡ -k(1/4)`）的记账：三阶段，已闭合**（勿把中间阶段的数字
  当独立错误计数）：
  * **阶段 1（R6.1）**：假设"`t` 与 `t + 1` 是同一个 parent irrep"，据此把 `k` 约化
    （`canonical_wave_vector`）并让审计比较同标签的两侧。该假设对含分数平移的线小群
    **过强**：`v` 是母群倒格矢，而 `Φ_v = exp(2πi v·T)` 在这些群上非平凡。
  * **阶段 2（`efe8abb`/`315e49`）**：在**仍约化**的读法下记账：`t = 3/4` 上 58 行
    fail-closed（非整数重数）、50 行违反空间群共轭律，字符层 223 行不共轭（实表
    5,510 可比）。当时把 58 定位成 "transport/assembly 缺陷"。
  * **阶段 3（本轮，判决）**：撤销约化后，`t = 1/4, 3/4, 5/4` 三个参数的
    **hard failure 全为 0**（各 5,756 行），ledger 的独立复算里共轭律违反为 0
    （`line_transport_ledger --batch`），字符层位移比较 5,756/5,756 相等。即阶段 2 的
    58/50/223 **是规范波矢的产物，不是 frozen 表或求解器的缺陷**：把冻结表配到它
    不描述的波矢上，得到的自然不是表示。
  * **操作层面**：`line_family_coverage --gate` 的 `CONJUGATE_GAP_CONTENT/ERRORS` 由
    `187/58` 改为 `0/0`，字符层由"`χ(3/4)` 对 `conj(χ(1/4))`"（缺 twist，223 差异）
    改为"`χ(5/4)` 对 `χ(M_v(α), 1/4)`"（0 差异）。复表（`DT3`/`DT4`，188 行）继续用
    伙伴源 pinned 作 oracle，0 不一致、0 错误。
* **仍然成立的边界**：`t = 3/4` 一类**不在已验证域**。那里的字符仍是一致的
  （共轭律 0 违反、范数整数），但它既不是 pinned 的位移像、也不由任何冻结标签描述，
  所以**没有 oracle**：报告只给实测值。`t = 0`、`t = 1/2` 同理（稳定子严格大于冻结
  小群时引擎算的是形式表示）。
* `t = 0`、`t = 1/2` 一类（`t·v` 的稳定子严格大于冻结小群）若能完整计算，会在
  `LineSubduction::parameter_kind()` 标为 `ParameterKind::Formal`：结果是从冻结线小群
  得到的形式诱导，不是该点增强小群的线 irrep 分导，也不作为物理结论。判定条件是存在
  不同方向臂 `a` 使 `t·(a-v) ∈ L*_parent`；缺少子群目标字符时仍显式报错。永久类型测试
  钉住 SG 196 `DT1` 的 `t=0,1/2` 与 `1/4,1/6,1/3`；覆盖门禁尚未统计全表的类型标签。
* 一般参数下**非 Γ 目标**没有外部对照（官方 pinned 只有恒等频率），只有引擎自检。

**已验证域之外仍未承诺**：`t ∉ 1/4 + Z` 的参数（`1/2`、`3/4`、`1/6` …）没有 oracle，
引擎的答案是一致但**未经验证**的形式值；`t ∉ (1/4)Z` 的非退化点上 `content = 0` 有
证明级的几何论证（见 coverage 报告命题 1/2），但"某个具体非零/非 pinned 值"不作
结论。参数化源的"整族结论"、例外集与覆盖说明见
[subduction-r6-coverage.md](subduction-r6-coverage.md)。
