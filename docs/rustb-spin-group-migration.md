# Rustb 迁移指南：旧 `collinear_spin` → `spin_group`

面向 Rustb 的 AI/维护者。cryspglib 的旧公共 API `cryspglib::collinear_spin` 与其根导出已删除；
新的稳定入口是 `cryspglib::spin_group`。Rustb 只需要**共线**部分；共面/非共面未实现。

## 0. 一句话变化

旧 API 的 `ValidatedCollinearSpinGroup` 是**幺正子群**，需要显式 `.with_spin_only_time_reversal()`
才得到含 `A0=[C2_perp T||E]` 的完整群；新 API 的 `detect_collinear_spin_group` 直接返回
`CollinearSpinGroup`（完整有限群，幺正半群 + 反幺正陪集），通道是 `CollinearAction`，反幺正性是
独立的 `OperationKind`，两者不再混在一个 `SpinChannelAction`/`time_reversal` 字段里。

## 1. Cargo.toml

开发期继续用路径依赖即可：

```toml
cryspglib = { path = "../cryspglib", version = "0.2.1", optional = true, default-features = false }
```

cryspglib 分支验收并发布/推送后，改为固定的 git revision 或 crates.io 版本；不要使用未验收的
本地路径作为发布依赖。

## 2. 类型与入口映射

| 旧 | 新 |
|---|---|
| `cryspglib::CollinearMoments` | `cryspglib::spin_group::CollinearMoments` |
| `cryspglib::CollinearSpinTolerances` | `cryspglib::spin_group::SpinTolerances` |
| `cryspglib::SpinChannelAction` | `cryspglib::spin_group::CollinearAction` |
| `cryspglib::ValidatedCollinearSpinGroup` | `cryspglib::spin_group::CollinearSpinGroup` |
| `cryspglib::CollinearSpinTimeReversalGroup` | `cryspglib::spin_group::CollinearSpinGroup`（同一个完整群类型） |
| `cryspglib::detect_collinear_spin_group` | `cryspglib::spin_group::detect_collinear_spin_group` |
| `cryspglib::validate_collinear_spin_group` | `cryspglib::spin_group::validate_collinear_spin_group` |
| `operation.operation.rotation` | `operation.rotation` |
| `operation.operation.translation` | `operation.translation` |
| `operation.operation.channel_action` | `operation.action` |
| `operation.kind.is_antiunitary()` | 不变（类型仍是 `cryspglib::OperationKind`，由 crate 根导出） |
| `group.with_spin_only_time_reversal()` | 不需要：新检测已含反幺正半群 |
| 旧完整群 `CollinearSpinTimeReversalGroup::operations()` 返回 `impl Iterator<...Ref>`；旧幺正子群 `ValidatedCollinearSpinGroup::operations()` 已是切片 | 新 `CollinearSpinGroup::operations()` 统一返回 `&[SpinOperation<CollinearAction>]`，直接 `.iter()` |
| `group.spin_only_time_reversal_index()` | 不再需要；用 `operations().iter().position(|op| op.kind.is_antiunitary())` 找代表 |
| `CollinearSpinTolerances { position, moment, collinearity }` | `SpinTolerances { position, moment, direction }`；当前共线引擎把 `direction` 按 Cartesian 垂直磁矩残差使用（与旧 `collinearity` 同量纲，不是角度），所以数值可 1:1 对应；仍建议先核对语义再用 `SpinTolerances::default()` |

`SpinTolerances` 默认三个字段均 `1e-5`。当前共线检测中 `direction` 是垂直磁矩残差
（Cartesian 范数，与 `moment` 同量纲），不是弧度；不要按角度换算。磁矩输入
`CollinearMoments::Signed`/`Vectors` 语义不变；
`CollinearOrder` 可通过 `group.order()` 读取 `signed_moments()`、`axis()`、`collinearity_residuals()`。

## 3. centre.rs 的具体改法

`src/hamiltonian_symmetry/centre.rs` 中：

1. `PairOperation.channel` 的类型改为 `cryspglib::spin_group::CollinearAction`（若 Rustb 重新导出，
   就改为 Rustb 的新 re-export）。所有 `match`/赋值里的 `SpinChannelAction::Preserve`、`::Swap`
   替换为 `CollinearAction::Preserve`、`::Swap`。
2. 反幺正 pair 操作生成函数改成只吃完整群：

```rust
fn pair_operations_from_spin_group(
    group: &cryspglib::spin_group::CollinearSpinGroup,
) -> Vec<PairOperation> {
    group
        .operations()
        .iter()
        .map(|operation| PairOperation {
            spatial: CrystalSymmetryOperation {
                rotation: operation.rotation,
                translation: operation.translation,
                time_reversal: operation.kind.is_antiunitary(),
            },
            channel: operation.action,
        })
        .collect()
}
```

3. `pair_operations_from_unitary_spin_group` 不再有独立类型；从完整群过滤：

```rust
fn pair_operations_from_unitary_spin_group(
    group: &cryspglib::spin_group::CollinearSpinGroup,
) -> Vec<PairOperation> {
    group
        .operations()
        .iter()
        .filter(|operation| operation.kind.is_unitary())
        .map(|operation| PairOperation {
            spatial: CrystalSymmetryOperation {
                rotation: operation.rotation,
                translation: operation.translation,
                time_reversal: false,
            },
            channel: operation.action,
        })
        .collect()
}
```

4. 两个公开方法的签名：

```rust
pub fn symmetrize_pair_with_collinear_spin_group(
    a: &mut Self,
    b: &mut Self,
    group: &cryspglib::spin_group::CollinearSpinGroup,
    basis: &(impl BasisSymmetryRepresentation<false, R> + Sync),
    parameters: &HamiltonianSymmetrizationParameters,
) -> Result<PairSymmetrizationReport>;

pub fn symmetrize_pair_with_unitary_collinear_spin_group(
    a: &mut Self,
    b: &mut Self,
    group: &cryspglib::spin_group::CollinearSpinGroup,
    basis: &(impl BasisSymmetryRepresentation<false, R> + Sync),
    parameters: &HamiltonianSymmetrizationParameters,
) -> Result<PairSymmetrizationReport>;
```

第一个使用完整群（含 `A0` 和所有 `gA0`），第二个内部过滤 `kind.is_unitary()`。**不要**再把
旧 `with_spin_only_time_reversal()` 后的群当作“完整群”传入并再次扩展；新检测已经是完整群。
注意完整群的 `len()` 是旧幺正子群的 2 倍（FM/无 A0 破缺时）。

5. 旧 `symmetrize_pair`（吃 `MagneticCrystalSymmetry`，把 `time_reversal` 当作 channel swap）
   保持原语义，不要用它表达无 SOC 共线群；它仍可服务真正的 MSG 输入。

## 4. 检测/验证调用

旧：

```rust
let group = cryspglib::detect_collinear_spin_group(
    &crystal, CollinearMoments::Signed(&moments), CollinearSpinTolerances::default(),
)?;
let full = group.with_spin_only_time_reversal();
```

新：

```rust
let group = cryspglib::spin_group::detect_collinear_spin_group(
    &crystal,
    cryspglib::spin_group::CollinearMoments::Signed(&moments),
    cryspglib::spin_group::SpinTolerances::default(),
)?;
// group 就是完整群：既有幺正操作也有 kind == Antiunitary 的 gA0 伙伴
```

显式验证的 target 操作类型改为 `SpinOperation<CollinearAction>`：

```rust
let target: Vec<cryspglib::spin_group::SpinOperation<cryspglib::spin_group::CollinearAction>> = ...;
let group = cryspglib::spin_group::validate_collinear_spin_group(
    &crystal, moments, &target, tolerances,
)?;
```

每个 target 用 `rotation`、`translation`、`action`、`kind` 四个公开字段；不再有
`channel_action` 或 `operation.operation` 包装。

## 5. 测试迁移

`tests/collinear_spin_pair.rs`：

- 导入改为 `use cryspglib::spin_group::{CollinearMoments, SpinTolerances, detect_collinear_spin_group};`
  （若 Rustb 在 `lib.rs` re-export，也可从 `Rustb::{...}` 导入）。
- `CollinearSpinTolerances::default()` → `SpinTolerances::default()`。
- 旧的“full group”测试不再调用 `.with_spin_only_time_reversal()`；检测结果直接是完整群。
- 旧的“unitary subgroup”测试：把完整群传给新的 `symmetrize_pair_with_unitary_collinear_spin_group`，
  该方法内部过滤幺正操作；不要手工从完整群删掉反幺正操作再另造类型。
- 同步更新 Rustb 文档注释里仍引用旧类型名（`ValidatedCollinearSpinGroup`、
  `CollinearSpinTimeReversalGroup`、`SpinChannelAction`）的地方，避免 rustdoc 断链。
- 断言操作数量时注意完整群是 2n；若测试只想验证幺正半群，改为
  `group.operations().iter().filter(|op| op.kind.is_unitary()).count()`。

`src/lib.rs` 的 re-export 改为：

```rust
pub use cryspglib::spin_group::{
    CollinearAction, CollinearMoments, CollinearSpinGroup, SpinTolerances,
    detect_collinear_spin_group,
};
```

## 6. 可选：共线 SSG 编号

若 Rustb 需要诊断/记录磁性结构的自旋空间群编号，可调用：

```rust
let identification = cryspglib::spin_group::identify_collinear_spin_space_group(&group)?;
println!("{}", identification.identifier()); // 例如 221.221.1.1.L
```

这是**标准胞精确**入口：输入晶体必须已经用数据库记录的基（`G0` 或 `L0`）和常规胞描述；
非标准胞返回 `NoExactMatch`。它不影响 pair 投影，pair 层不需要它。

## 7. 边界提醒

- 新 API 的完整群描述“磁结构的自旋对称性”，不证明实际 Hamiltonian 没有 SOC/外磁场/轨道 TR 破缺；
  Rustb 仍必须自行声明/验证这些前提。
- 共线 SSG 识别不需要有向磁轴；全局自旋旋转不改变编号。
- 共面/非共面、一般 setting 搜索、SU(2) 中心相位、轨道 sewing、H 投影均不在当前 cryspglib 实现内。
