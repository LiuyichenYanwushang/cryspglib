# 磁空间群 UNI=0 故障复盘：标准设置变换方向错误

## 1. 问题摘要

目标测试：

```bash
cargo test --package cryspglib --test magnetic_integration test_graphene_afm_z
```

测试体系是石墨烯蜂窝晶格，两个子晶格具有相反的 z 方向磁矩。修复前，
程序能够正确识别：

```text
普通空间群: SG 191, P6/mmm
Hall 编号: 485
磁群类型: Type-3 BlackWhite
磁操作数: 24 = 12 unitary + 12 anti-unitary
```

但是磁空间群数据库匹配失败：

```text
UNI=0
BNS=""
```

而且从 `symprec=1e-3` 到 `1e-6`，结果始终不变。

修复后，所有这些容差都稳定返回：

```text
UNI=1466
BNS=191.236
```

根因不是二维材料的物理特殊性，也不是磁操作生成不完整，更不是容差不够。
根因是磁空间群标准设置变换的坐标约定在 Rust 移植中被写反了。

涉及两个紧密关联的错误：

1. `get_reference_space_group` 把 `ref_sg.bravais_lattice` 本身当成了
   坐标变换矩阵 `tmat`，实际上应该使用它的逆矩阵。
2. `get_distinct_changed_magnetic_symmetry` 使用了逆向共轭
   `T⁻¹ R T`，但调用方约定是 `x_std = T x + s`，因此正确公式应为
   `T R T⁻¹`。

这两个错误在部分路径上会互相抵消，导致立方体系的既有测试继续通过；
但该 helper 还用于数据库候选的 correction transformation，在那里没有
第二个错误可以抵消，所以非正交六方基底暴露了问题。

相关修复提交：

```text
8806eb9 fix: align magnetic symmetry basis transforms
be8992f test: lock graphene AFM magnetic group
```

---

## 2. 为什么首先判断这不是磁操作生成错误

故障输出已经提供了几个非常强的约束：

### 2.1 普通空间群正确

程序识别出：

```text
SG=191
Hall=485
```

这说明：

- 晶格矩阵布局 `lattice[cart][vec]` 正确；
- 原子位置可以被普通空间群搜索正确处理；
- 六方晶格没有因为矩阵转置而退化到低对称空间群。

如果晶格矩阵本身写错，石墨烯会先在普通空间群阶段失败，而不是只在 UNI
匹配阶段失败。

### 2.2 磁群类型正确

程序返回：

```text
MagneticType::BlackWhite
```

Type-3 分类来自 FSG/XSG 的阶数关系。它正确意味着：

- 忽略 time reversal 得到的 family space group 数量合理；
- 只保留 unitary 操作得到的 maximal subgroup 数量合理；
- unitary 和 anti-unitary 的陪集结构基本正确。

### 2.3 磁操作数和时间反演分布正确

程序生成：

```text
24 operations
12 unitary
12 anti-unitary
```

这说明 `operations_with_site_tensors` 已经正确处理了：

- 两个碳原子的子晶格交换；
- z 方向磁矩作为 axial vector 的变换；
- 需要附加 time reversal 的空间操作。

因此没有理由先去修改磁矩比较、原子映射或 axial-vector 公式。

### 2.4 改变 symprec 完全无效

以下四个容差给出完全相同的结果：

```text
1e-3, 1e-4, 1e-5, 1e-6
```

这基本排除了：

- 浮点误差刚好越过阈值；
- 原子重合判断不稳定；
- translation 只差一个很小数值；
- 数据库匹配需要更宽容差。

如果输出在多个数量级的容差下完全稳定，继续调容差通常是在掩盖结构性错误。

---

## 3. 失败链路的定位

磁空间群识别的主要链路是：

```text
Crystal::magnetic_dataset
  -> operations_with_site_tensors
  -> identify_with_parent_hall
  -> reduce_to_primitive_magsym
  -> get_reference_space_group
  -> get_changed_magnetic_symmetry
  -> get_uni_candidates
  -> get_distinct_changed_magnetic_symmetry
  -> is_subset
```

Graphene AFM 已经顺利通过了前面的空间群识别、磁操作生成和磁类型分类。
失败发生在 Hall 485 的 Type-3 UNI 候选匹配阶段。

Hall 485 的相关候选是 UNI 1465–1471。数据库记录能够正常读取，但经过当前
setting/origin correction 后，没有任何候选通过 `is_subset`。

这将搜索范围缩小为：

1. 数据库操作是否错误；
2. 输入磁操作是否没有变换到数据库使用的 setting；
3. correction transformation 的共轭方向是否错误；
4. origin shift 的符号或作用顺序是否错误。

数据库在大量其他测试中工作正常，而输入磁操作的数量和群结构也正确，所以
最值得优先检查的是“输入操作到标准设置”的坐标变换。

---

## 4. 关键诊断方法：比较 passing case 与 failing case

既有磁群测试主要是：

- Pm-3m；
- Im-3m；
- Fm-3m；
- 沿高对称方向降到四方或三方的立方晶格案例。

这些案例大量使用：

- 单位矩阵；
- 轴置换矩阵；
- 正交矩阵；
- 对称矩阵；
- 自逆矩阵。

而 graphene 使用非正交六方基底：

```text
a = (1, 0, 0)
b = (1/2, sqrt(3)/2, 0)
c = (0, 0, 2)
```

其基变换矩阵通常不满足：

```text
T = T⁻¹
```

也不与所有旋转矩阵对易。

因此出现了一个非常明确的差异因子：

> 立方测试通过、六方测试失败，优先检查 basis/setting transformation，
> 不要优先怀疑与晶系无关的磁矩分类算法。

这个判断把问题从整个磁性识别管线缩小到了
`magnetic_spacegroup.rs` 中的坐标变换代码。

---

## 5. 数学 oracle：从坐标变换重新推导 Seitz 共轭

不能只凭记忆判断 `T R T⁻¹` 还是 `T⁻¹ R T`。最可靠的方法是从调用方
声明的坐标关系直接推导。

代码和上游 spglib 的约定是：

```text
x_std = T x + s
```

原 setting 中的 Seitz 操作是：

```text
g(x) = R x + t
```

定义坐标变换：

```text
C(x) = T x + s
C⁻¹(x_std) = T⁻¹ (x_std - s)
```

标准 setting 中的操作必须是：

```text
g_std = C g C⁻¹
```

逐步展开：

```text
g_std(x_std)
  = C(g(T⁻¹(x_std - s)))
  = C(R T⁻¹(x_std - s) + t)
  = T R T⁻¹ x_std - T R T⁻¹ s + T t + s
```

因此：

```text
R_std = T R T⁻¹
t_std = s - R_std s + T t
```

这就是实现必须满足的数学 oracle。任何与此不一致的代码都在实现反方向变换。

随后将 Rust 实现逐函数对照上游 spglib：

```text
https://github.com/spglib/spglib/blob/develop/src/magnetic_spacegroup.c
```

重点比较了：

```text
get_reference_space_group
get_changed_magnetic_symmetry
get_distinct_changed_magnetic_symmetry
```

上游源码的注释、矩阵乘法顺序和上述独立推导一致，因此可以把它作为第二个
oracle。这里不是“照抄上游就算正确”，而是先由坐标定义推导公式，再确认上游
实现与推导相符；两条独立证据同时指向 Rust 移植的矩阵方向错误。

---

## 6. 发现的第一个错误：`tmat` 取反了

文件：

```text
src/magnetic_spacegroup.rs
```

函数：

```text
get_reference_space_group
```

### 6.1 旧代码

```rust
let tmat = ref_sg.bravais_lattice;
let shift = ref_sg.origin_shift;
```

但 `ref_sg.bravais_lattice` 在这里表示的是晶格关系中的矩阵 `P`：

```text
(a_std, b_std, c_std) = (a, b, c) P
```

分数坐标与晶格基矢做逆变换，所以对应的坐标关系是：

```text
x_std = P⁻¹ x + p
```

因此：

```text
T = P⁻¹
```

而不是 `P`。

### 6.2 缺失的真实晶格度量规整

旧代码还直接使用了在 unit lattice 搜索得到的 `bravais_lattice`。
上游 spglib 会先将该矩阵映射到输入晶格的真实笛卡尔度量中，调用
`ref_find_similar_bravais_lattice`，然后再映射回分数基底。

这一步对于非正交晶格尤其重要。修复后流程为：

```rust
let lattice_inv = mat_inverse_matrix_d3(lattice, 0.0).ok()?;

ref_sg.bravais_lattice =
    mat_multiply_matrix_d3(lattice, &ref_sg.bravais_lattice);

ref_find_similar_bravais_lattice(ref_sg, symprec);

ref_sg.bravais_lattice =
    mat_multiply_matrix_d3(&lattice_inv, &ref_sg.bravais_lattice);

let tmat =
    mat_inverse_matrix_d3(&ref_sg.bravais_lattice, 0.0).ok()?;
```

同时 `get_reference_space_group` 增加 `lattice` 参数，因为没有真实输入晶格，
无法完成这一步度量规整：

```rust
fn get_reference_space_group(
    lattice: &Mat3,
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
)
```

---

## 7. 发现的第二个错误：Seitz 共轭方向写反

文件：

```text
src/magnetic_spacegroup.rs
```

函数：

```text
get_distinct_changed_magnetic_symmetry
```

### 7.1 旧代码的旋转变换

```rust
let tmp = mat_multiply_matrix_d3(&inv_tmat, &rot_f64);
let r_new = mat_multiply_matrix_d3(&tmp, tmat);
```

即：

```text
R_new = T⁻¹ R T
```

这对应的是反方向坐标变换，不符合函数调用方的：

```text
x_std = T x + s
```

### 7.2 旧代码的平移变换

```rust
t_new = T⁻¹ [t + (R - I)s]
```

它同样属于逆向变换，而且 shift 是与旧旋转 `R` 组合，而正确公式需要使用
变换后的 `R_std`：

```text
t_std = s - R_std s + T t
```

### 7.3 修复后的实现

```rust
let tmp = mat_multiply_matrix_d3(tmat, &rot_f64);
let r_new = mat_multiply_matrix_d3(&tmp, &inv_tmat);

let rotated_shift =
    mat_multiply_matrix_vector_id3(&rot_i, shift);
let transformed_trans =
    mat_multiply_matrix_vector_d3(tmat, &sym_msg.trans[i]);

for j in 0..3 {
    t_new[j] = mat_dmod1(
        shift[j] - rotated_shift[j] + transformed_trans[j],
    );
}
```

对应：

```text
R_std = T R T⁻¹
t_std = s - R_std s + T t
```

最后用 `mat_dmod1` 将平移归一化到晶胞周期等价类。

---

## 8. 为什么两个错误没有让所有测试都失败

这是本次故障最值得保留的调试经验。

旧代码同时存在：

```text
错误 A: 把 T 写成 P，而真实 T 应是 P⁻¹
错误 B: helper 计算 T⁻¹ R T，而不是 T R T⁻¹
```

把错误 A 代入错误 B：

```text
T_old⁻¹ R T_old
  = P⁻¹ R P
```

而真实变换为：

```text
T_real R T_real⁻¹
  = P⁻¹ R P
```

所以在旋转部分，这两个错误可能刚好抵消。

这解释了为什么：

- 很多已有测试继续通过；
- 仅凭“立方测试全部通过”不能证明变换约定正确；
- 单独观察初始 reference-setting 变换可能看不出问题。

但是 `get_distinct_changed_magnetic_symmetry` 不只用于初始 reference setting。
它还用于：

```text
get_std_transformations
  -> correction transformation
  -> get_distinct_changed_magnetic_symmetry
```

数据库返回的 correction transformation 已经是调用方所需的 `T`。这里不存在
“把 `T` 写成 `T⁻¹`”的配对错误，因此 helper 内部的逆向共轭不会被抵消。

对于 identity correction：

```text
T = I
```

正向和逆向公式相同，所以仍可能通过。

对于非平凡 correction，尤其是六方非正交 setting：

```text
T R T⁻¹ != T⁻¹ R T
```

于是 transformed magnetic operations 无法与数据库中的 UNI 操作对齐，
`is_subset` 对所有候选返回 false。

这是一类典型的“两个 bug 在部分路径上互相掩盖，换一条调用路径后暴露”的问题。

---

## 9. 为什么不应该把问题归因于二维材料

修复前曾有一个错误结论：

> 二维/准二维结构与 3D MSG 数据库 setting 不兼容，所以 UNI=0 是已知限制。

这个结论的问题是，它没有解释以下事实：

1. 普通 SG、Hall、磁类型、操作数全部正确；
2. 失败对 `symprec` 不敏感；
3. 数据库中明确存在同 Hall、同 Type-3 的候选；
4. 失败点恰好位于 setting correction 后的集合比较；
5. 当前实现与上游的矩阵方向不一致。

一旦修正坐标变换，同一个二维结构立即稳定识别为 UNI 1466。

因此更可靠的原则是：

> 在宣布“数据库限制”或“物理特殊性”之前，先证明坐标变换、数据 setting
> 和操作共轭与权威实现一致。

“复杂体系失败”不自动意味着“复杂物理不受支持”；它也可能只是最先暴露了
一个被高对称测试掩盖的基础线性代数错误。

---

## 10. 修复后的回归测试

`tests/magnetic_integration.rs::test_graphene_afm_z` 不再只打印结果。
它对四个容差分别断言：

```rust
assert_eq!(r.spacegroup_number, 191);
assert_eq!(r.hall_number, 485);
assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
assert_eq!(r.uni_number, 1466);
assert_eq!(r.bns_number.trim(), "191.236");
assert_eq!(r.num_operations, 24);
```

并检查：

```text
12 unitary
12 anti-unitary
```

容差覆盖：

```text
1e-3
1e-4
1e-5
1e-6
```

这种断言比只检查 `uni_number > 0` 更强，因为它能捕获：

- 匹配到了错误候选；
- BNS/UNI 映射回归；
- time reversal 分布变化；
- 操作数异常；
- 对容差敏感的不稳定识别。

---

## 11. 回归验证结果

目标测试：

```bash
cargo test --package cryspglib \
  --test magnetic_integration test_graphene_afm_z
```

结果：

```text
1 passed
```

磁性集成测试：

```bash
cargo test --package cryspglib --test magnetic_integration
```

结果：

```text
10 passed
```

CoF3 磁群回归：

```bash
cargo test --package cryspglib --test cof3
```

结果：

```text
2 passed
```

完整 crate 测试：

```bash
cargo test --package cryspglib
```

结果：

```text
unit tests: 116 passed, 1 ignored
integration tests: all passed
doctests: 26 passed
```

重要的既有磁群结果保持不变：

```text
Fe SC [001]        -> UNI 1005, BNS 123.345
Fe SC [100]        -> UNI 1005, BNS 123.345
Fe BCC AFM [111]   -> UNI 1331, BNS 166.101
FCC FM [001]       -> UNI 1005, BNS 123.345
FCC FM [111]       -> UNI 1331, BNS 166.101
CoF3 magnetic      -> UNI 1333, BNS 167.103
Graphene AFM-z     -> UNI 1466, BNS 191.236
```

---

## 12. 可复用的调试方法

### 12.1 先按阶段确认哪些信息已经正确

不要把 `UNI=0` 当成“整个磁性识别都失败”。

分别检查：

```text
普通 SG 是否正确？
Hall 是否正确？
磁类型是否正确？
操作数是否正确？
unitary/anti-unitary 比例是否正确？
失败是否仅发生在 DB matching？
```

这些信息可以把搜索范围从整个系统缩小到单个阶段。

### 12.2 容差扫描是分类工具，不是修复方法

如果结果随容差变化，优先检查数值稳定性。

如果结果跨多个数量级完全不变，优先检查：

- 坐标约定；
- 数据 setting；
- 变换方向；
- 群操作合成；
- 数据库映射。

### 12.3 比较 passing 与 failing 的差异因子

本例：

```text
passing: 大量立方、正交、自逆变换
failing: 六方、非正交、T != T⁻¹
```

差异直接指向 basis transformation。

### 12.4 对矩阵方向使用推导，不使用记忆

看到以下任一形式时：

```text
T R T⁻¹
T⁻¹ R T
s - R' s + T t
T⁻¹(t + (R-I)s)
```

不要凭变量名判断。先写清：

```text
x_new 与 x_old 的关系是什么？
```

然后计算：

```text
C g C⁻¹
```

推导结果是最小、最可靠的 oracle。

### 12.5 检查 helper 的所有调用者

一个错误 helper 可能在某个调用点被另一个错误抵消，却在其他调用点失败。

本例中同一个 helper 同时服务于：

1. reference-setting transformation；
2. database correction transformation。

只分析第一个调用点会错误地认为公式“等效”。

### 12.6 高对称测试不足以验证线性代数约定

应当至少包含：

- 正交但非单位的变换；
- 非正交六方/单斜晶格；
- 非零 origin shift；
- 非平凡 alternative/correction transformation；
- `T != T⁻¹` 的案例。

只有立方测试很容易让转置、逆矩阵和共轭方向错误长期潜伏。

---

## 13. 后续修改时的检查清单

修改磁空间群 setting 相关代码前，逐项确认：

1. 当前坐标约定是否明确写成 `x_new = T x_old + s`？
2. 旋转是否使用 `R_new = T R_old T⁻¹`？
3. 平移是否使用 `t_new = s - R_new s + T t_old`？
4. 晶格基矢变换与分数坐标变换是否正确取逆？
5. `ref_sg.bravais_lattice` 表示的是晶格变换还是坐标变换？
6. 是否在真实输入晶格度量下调用了 lattice refinement？
7. helper 是否还被 correction/origin-alternative 路径调用？
8. 测试是否包含非正交且 `T != T⁻¹` 的晶格？
9. 测试是否断言具体 UNI/BNS，而不只是 `UNI > 0`？
10. 是否运行了已有立方、三方和真实材料磁群回归？

只要这十项中有一项无法回答，就不应断言 setting transformation 已经正确。
