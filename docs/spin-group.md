# 无 SOC 自旋群：共线有限群、SSG 数据库与标准胞识别

本模块独立于磁空间群（MSG）：不复用、也不写回 `SymmetryOp.time_reversal`、旧磁操作搜索或
MSG 数据库。自旋通道是独立的 `CollinearAction`（`Preserve`/`Swap`），反幺正性是独立的
`OperationKind`；两者不可互换。共线情形只需要一个有符号标量序参量，不需要有向磁轴。

当前版本实现**共线**部分；共面、非共面接口保留但未实现。旧 `collinear_spin` 模块已从公共
API 删除，作为新模块的内部引擎保留；公开入口统一在 `cryspglib::spin_group` 下。

## 范围与非目标

- 三维周期输入胞、无 SOC、非零共线磁序、已经自旋块对角的两个通道。
- 空间操作正规化输入胞的整数平移格：`W ∈ GL(3,Z)`，商掉输入胞的 `Z³`；保留胞内分数平移。
  任意超胞中不一定能恢复原胞的全部点操作。
- 共线序参量写成 `m_i = σ_i n`。`Preserve`/`Swap` 是作用在 σ 上的净符号；连续绕轴旋转、
  含轴镜面与反幺正自旋伙伴构成规范核，不被采样为有限操作行。
- 完整群包含幺正半群 `G` 和反幺正陪集 `G·A0`，`A0=[C2_perp T||E]`。`Swap` 的幺正 π 旋转
  代表不共轭自旋块、不产生动量负号；`A0` 反幺正、`k → -k`，但不与 `Swap` 混同。
- 不实现：共面/非共面、不可约表示、SU(2) 提升/中心相位、一般 setting 搜索、自旋子 sewing、
  Hamiltonian 投影。几何证书不是 Hamiltonian 证书；调用者必须自行排除 SOC、外磁场、轨道
  时间反演破缺等前提。
- 层群、全零磁矩、不可辨识磁轴、非共线输入明确拒绝，不静默选轴或补通道。

## 公共 API

```rust
use cryspglib::spin_group::{
    CollinearMoments, SpinTolerances, detect_collinear_spin_group,
    validate_collinear_spin_group,
};
use cryspglib::Crystal;

let crystal = Crystal::new(
    [[2., 0., 0.], [0., 3., 0.], [0., 0., 4.]],
    vec![[0., 0., 0.], [0.5, 0., 0.]],
    vec![26, 26],
)?.with_magnetic(vec![[0., 0., 2.], [0., 0., -2.]])?;

let moments = [2., -2.];
let tolerances = SpinTolerances::default();
let group = detect_collinear_spin_group(
    &crystal,
    CollinearMoments::Signed(&moments),
    tolerances,
)?;

assert_eq!(group.kind(), cryspglib::spin_group::SpinGroupKind::Collinear);
for operation in group.operations() {
    // operation.rotation / operation.translation / operation.action / operation.kind
    let _ = operation;
}
# Ok::<(), Box<dyn std::error::Error>>(())
```

主要类型：

- `SpinTolerances { position, moment, direction }`：三个独立绝对容差，默认均 `1e-5`。
- `SpinOperation<A> { rotation, translation, action, kind }`：`rotation`/`translation` 是输入胞
  分数坐标下的空间 Seitz 部分；`action` 是序参量符号；`kind` 是幺正/反幺正。
- `SpinGroup<A,O,M>`：只读的完整有限群。`operations()`、`atom_mappings(g)`、`product(g,h)`、
  `inverse_index(g)`、`identity_index()`、`tolerances()`、`order()`。
- `CollinearSpinGroup = SpinGroup<CollinearAction, CollinearOrder, f64>`；
  `CollinearOrder` 提供 `signed_moments()`、`axis()`、`collinearity_residuals()`。
- `SpinAtomMapping`：source atom → target、整数 cell shift、分数位置残差、磁矩残差。
- `SpinGroupProduct`：群乘积的 operation index、整数 cell shift、平移/作用残差。

`validate_collinear_spin_group(crystal, moments, target_operations, tolerances)` 接受显式的
**完整目标操作集**，逐项验证单位、重复、全部有序乘积、两侧逆、单位/逆、结合律、原子置换、
整数平移运输与 cocycle；不增删/不闭合/不降低容差。检测后的完整候选集也走同一验证器。

`CollinearMoments` 支持 `Signed(&[f64])` 和 `Vectors(&[[f64;3]])`。向量输入先检查共线残差，
保留实际范数与符号；不能把大小不匹配藏在横向分量里。零矩原子仍参与结构映射。

## 有限群证书

对每个操作 `g` 和原子 `i`：

```text
W_g r_i + t_g = r_{π_g(i)} + ℓ_g(i) + δ_g(i)
m_{π_g(i)} = η_g · m_i + μ_g(i)
||A δ_g(i)|| <= position,  |μ_g(i)| <= moment
```

乘积表记录先作用 h 再作用 g 的代表 p、整数 q 与 Seitz 残差 ε：

```text
W_p = W_g W_h,  η_p = η_g η_h,  kind_p = kind_g XOR kind_h
t_g + W_g t_h = t_p + q(g,h) + ε(g,h)
π_p(i) = π_g(π_h(i))
ℓ_g(π_h(i)) + W_g ℓ_h(i) = ℓ_p(i) + q(g,h)
```

验证器另外检查 `q(g,h)+q(gh,k) = W_g q(h,k)+q(g,hk)`。`Preserve`/`Swap` 的 `Swap²=Preserve`
只说明有限商类回到单位，**不**断言自旋子 π 旋转平方为 `+I`；SU(2) 中心相位不在本证书内。

## 共线 SSG 数据库

`spin_group::database` 带一个钉死的 FindSpinGroup 共线子集：

- 记录数 **1,421**；release `id=2`，`source_sha256=fcb74e80…`，
  manifest SHA-256 `4216cdee…`，identifier-list SHA-256 `29fa9858…`。
- 每条记录存**商签名**：空间操作 `{W|t}` 加通道符号。连续自旋旋转、含轴镜面与反幺正自旋伙伴
  是规范核，不重复成行。
- 公开编号为 `G0.L0.k_index.map_index.L`（四索引加 family 字母）。`k_index` 与 `map_index`
  对应导出器的 `k_index` 与 `equivalent map index`；原始 `t_index` 不在公开编号中。
- 查询：`collinear_spin_space_groups()`、`collinear_spin_space_group(number)`、
  `collinear_spin_space_group_by_identifier("221.221.1.1.L")`；记录暴露
  `identifier`、`international_symbol`、`international_notation_latex`、`group_type`、
  `basis`、`spin_point_group`、`tspin`、`gs`、`gso`、`num_all_elements` 与签名操作。
- 生成器 `scripts/generate_ssg_collinear_data.py` 只读 `target/` 下的钉死导出；生成文件
  `src/spin_group/database/ssg_collinear_data.rs` 提交入库，数据来源与 SHA 写在其头部。

## 共线 SSG 识别（标准胞精确）

```rust
use cryspglib::spin_group::{
    CollinearMoments, SpinTolerances, detect_collinear_spin_group,
    identify_collinear_spin_space_group,
};

let group = detect_collinear_spin_group(&crystal, CollinearMoments::Signed(&[2., -2.]), SpinTolerances::default())?;
let identification = identify_collinear_spin_space_group(&group)?;
assert_eq!(identification.identifier(), "47.47.2.1.L");
assert_eq!(identification.international_symbol(), "P_a1m1m1m∞m1");
# Ok::<(), Box<dyn std::error::Error>>(())
```

- **标准胞精确**：输入必须已经用数据库记录的基（`G0` 或 `L0`）和常规胞描述；不做 basis/origin/
  supercell 搜索。改变操作数据的非标准描述返回 `NoExactMatch`；只改变规范原点、不改变操作数据
  的变换仍识别同一编号。
- **无有向轴**：共线自旋空间在非相对论极限与实空间解耦；全局自旋旋转不改变编号。有限签名只用
  `(空间操作, 通道)`，`kind` 不进入签名（同一空间操作的幺正/反幺正伴带同一通道）。
- 比较容差：分数平移默认 `DEFAULT_COLLINEAR_IDENTIFICATION_TOLERANCE = 1e-8`，可用
  `identify_collinear_spin_space_group_with_tolerance` /
  `identify_collinear_signature_with_tolerance` 显式覆盖；旋转与通道精确比较。
- 低层入口 `identify_collinear_signature(&[SpinSpaceGroupOperation])` 可直接匹配商签名。
- 错误：`EmptyGroup`、`InconsistentChannel`、`NonFiniteTranslation`、`InvalidTolerance`、
  `NoExactMatch`（非标准胞/不在钉死 release）、`Ambiguous`（防御性；当前 release 签名唯一）。
- 边界：没有 `NeedsStandardization`/一般 setting 归一化；需要非标准胞识别时先做 basis/origin
  变换，或等待后续 setting 搜索。共面/非共面编号未实现。

## 数据来源与许可

数据库内容来自 FindSpinGroup（<https://findspingroup.com>、
<https://database.findspingroup.com>，Apache-2.0 代码与公开数据库），分类计数对应
Xiao et al., PRX 14, 031037 (2024) 的 1,421 共线 SSG。生成的 Rust 表只包含识别所需的商签名与
元数据；原始 JSON 导出不提交入库。引用时请同时引用数据库与分类论文。

## 测试与验收

Rust 验收覆盖：数据表计数/唯一性/来源 SHA、伪签名负控制、全部 1,421 条签名自识别、
共线 FM/AFM/反铁磁/三斜等检测与验证、全局自旋旋转与原点规范、容差边界、非标准胞拒绝、
完整群律与映射残差。每阶段另有独立子代理编写的 Rust 测试入库。

```bash
cargo test -p cryspglib --test spin_group_database
cargo test -p cryspglib --test spin_group_database_independent
cargo test -p cryspglib --test spin_group_identification
cargo test -p cryspglib --test spin_group_identification_independent
cargo test -p cryspglib --test spin_group_collinear
cargo clippy -p cryspglib --all-targets -- -D warnings
```
