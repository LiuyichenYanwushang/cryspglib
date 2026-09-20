# ISOTROPY 数据语义（isotropy subgroup 几何与分导）

本文记录 `isotropy_subgroup/iso.zip`（pinned，sha256 见
`scripts/generate_irrep_data.py`）中 isotropy 数据的字段语义、坐标系约定，以及
把它们钉死在官方程序输出上的验证方法。结论全部由随包 ISOTROPY 二进制
（`iso`，Version 9.6.1, Jan 2022，x86-64 静态链接）做 oracle 得到，可复现：

```bash
python3 scripts/verify_isotropy_oracle.py     # 42 行，6 种 centering（B 面心在标准 ITA setting 中不出现）
```

## 1. 表结构与索引

| 段 | 条数 | 含义 |
|---|---|---|
| `isotropy_subgroup` | 15239 | 非磁 isotropy 子群（SG 1–230） |
| `isotropy_parent` | 15239 | 母群 SG |
| `isotropy_irrep` | 15239 | 母群 irrep 序号（1–4777，ISO `data_irreps.txt` 顺序） |
| `isotropy_irrep_pointer` | 4778 | 每个 irrep 的记录区间起点 |
| `isotropy_basis` | 15239×9 | 子群格的 primitive 基（见 §2） |
| `isotropy_origin` | 15239×4 | 原点平移 `(x, y, z, d)`（见 §3） |
| `isotropy_direction` / `orderparam_dim` / `freeparam` / `orderparam_label` | 15239 | 方向 code、维数、自由参数数、ISO 方向标签 |
| `isotropy_domains` / `arms` | 15239 | domain 数、星臂数 |
| `isotropy_subduce_*` | 94271 | 分导频率表（见 §4） |
| `mag_iso_*` | 16721 | 磁子群表，**按同一批 4777 个非磁母群 irrep 索引**，子群为 UNI 1–1651 |

同一 irrep 内，方向描述串与 ISO 方向标签都是唯一的（已全表校验），因此方向可以
直接作为选择键。

## 2. Basis 的坐标系：母群 primitive 胞

`isotropy_basis` 的每一行是**子群格的一个 primitive 基矢**，表达在**母群
primitive 胞**基底下（ISOTROPY 程序内部使用的坐标系），系数为整数。

- 母群为 P 格时它与 conventional 基一致；C/A/B/I/F/R 格时不一致。
- 官方程序打印的是**子群的 ITA conventional 胞**（在母群 conventional 基下），
  因此存在 centering 的子群时两者数值不同但张成同一格：
  #16 P222 → #22 F222，机器数据 `(0,1,1),(1,0,1),(1,1,0)`（det 2），
  官方打印 `(2,0,0),(0,2,0),(0,0,2)`（det 8 = Z_F 4 × size 2）。
- 体积关系：`|det Basis_官方| = Z(子群) · Size / Z(母群)`，其中
  `Size = |det basis_机器|`（§5 的 oracle 逐行校验）。
- 想换到母群 conventional 基：`basis_in_parent_conventional(sg, basis)`
  （结果仍是子群格的 primitive 基，不是子群的 ITA conventional 胞）。

各 centering 的母群 primitive 基（官方 `P1` 行反推，Rust 里由
`parent_primitive_basis` 提供）：

```
P  (1,0,0),(0,1,0),(0,0,1)
I  (-1/2,1/2,1/2),(1/2,-1/2,1/2),(1/2,1/2,-1/2)
F  (0,1/2,1/2),(1/2,0,1/2),(1/2,1/2,0)
A  (0,1/2,1/2),(0,-1/2,1/2),(1,0,0)
B  (1/2,0,1/2),(-1/2,0,1/2),(0,1,0)
C  (1/2,1/2,0),(-1/2,1/2,0),(0,0,1)
R  (2/3,1/3,1/3),(-1/3,1/3,1/3),(-1/3,-2/3,1/3)   # hexagonal axes
```

## 3. Origin 的编码、坐标系，以及尚未钉死的打印约定

`isotropy_origin` 每条 4 个整数 `(x, y, z, d)`，即
`(x/d, y/d, z/d)`，同样是**母群 primitive 胞**下的坐标；分母集合为
`{1,2,3,4,6,8,12,16,24}`（生成器对此做了 fail-closed 校验）。

官方程序打印的是母群 **conventional** 基下的 origin（书中 "Origin" 列），例如
#167 R-3c 的 Γ3+ `(a,0)`：机器数据 `(0, 1/2, 0)`，官方打印
`(-1/6, 1/6, 1/6)`，正是 `1/2 · p2`（`p2 = (-1/3,1/3,1/3)`）。

Rust API：

- `IsotropyRecord::origin_rational()` / `origin_shift()`：原始（primitive）帧，
  精确有理数 / f64。
- `origin_shift_in_parent_conventional(sg, origin)`：纯坐标帧换算
  `w_conv = w_prim · P`。

**重要更正（对抗性审查后，2026-09-20）**：`w_prim · P` **不普遍等于**官方程序
打印的 Origin 列。对全部 4777 个 (SG, irrep) 做扫描（15044 条可比对记录）后，
约 **18.6%** 的记录两者之差既不是母群格矢量、也不是"每 SG 常数"：同一空间群内
不同记录可以有不同偏移（SG 141、227、230 内部都出现多种偏移，多数记录偏移为
0）。反例：父群 SG 139（I4/mmm）irrep `M1-` 方向 `P1` 在 primitive 帧存入
`(2,2,2)`（等价于母群原点），官方打印 `(1/4,1/4,1/4)`，差向量不属于 I 格。
原先记为"唯一例外"的 `SG 230 GM5+ P1` 只是这一大类中的一个实例。

因此：该函数的文档已改为"仅坐标帧换算"并显式声明该差异；
`scripts/verify_isotropy_oracle.py` 只是 **21 组 / 42 行抽样**，其通过不构成
origin 全表一致性的证据；官方打印 Origin 与存储值之间的确切约定**仍未钉死**
（见 §6）。

## 4. 分导（subduction）：`isotropy_subduce_*`

每条 isotropy 记录附带若干分导条目，四列并行：

- `isotropy_subduce_irrep`：母群 irrep 序号（1–4777，ISO 顺序）；
- `isotropy_subduce_frequency`：`i(G)`，该 irrep 的分导表示中包含子群单位表示
  （恒等表示）的次数；
- `isotropy_subduce_domain`：domain 序号；
- `isotropy_subduce_subgroup`：**已解开**——它是 1-based 的 isotropy 记录序号，
  它指向的那条记录的 `direction` 标签正是官方打印的 Dir 列（94271/94271 落在被
  分导 irrep 自己的记录区间内；30/30 抽样逐字符复现官方输出）。API 暴露为
  `IdentitySubduction::direction_label`。

这正是官方 `SHOW FREQ`（加 `DIR` 时附 `Dir(domain)`）的输出，例如 SG 221 `GM4+`
方向 `P1` → 子群 `#83 P4/m`：

```
官方:  83 P4/m  1 GM1+ P1(1), 1 GM3+ P1(3), 1 GM4+ P1(1)
API :  [("GM1+", 1, "P1", 1), ("GM3+", 1, "P1", 3), ("GM4+", 1, "P1", 1)]
```

**双值（spinor）分导条目**同样属于这一行输出：`isotropy_w_subduce_*` 共 5756 条，
覆盖 1006/15239 条记录，此前被整族丢弃；现在由
`double_valued_subduction(ordinal)` / `IsotropySubgroup::double_valued_subduction()`
提供。回归样例：SG 225 `W5` 方向 `S60` → 9 条标量 + `3 DT5, 3 SM3, 3 SM4`。

**边界**：本数据集只给“哪些母群 irrep 包含子群的恒等表示、重数多少”（Landau /
铁性分类所需），**不给**某个母群 irrep 分解成子群全部 irrep 的完整分导表示
（例如 Γ3+ ↓ P4/m = ?）。完整分解需要 ISODISTORT，或在本仓库自行实现字符表
分导引擎；`data_little.txt` 中的 `little_subduce_*` 是可能的线索，但其索引语义
未文档化，本轮未采用（未做猜测性实现）。

磁子群表没有对应的 `mag_iso_subduce_*`，因此磁子群只能给出几何 + UNI/BNS，
给不出分导。

## 5. 验证 gate

`scripts/verify_isotropy_oracle.py` 会：

1. 用 pinned 数据复现每条记录（`Size = |det W|`、方向标签、子群号）；
2. 对 21 组 (SG, irrep)（覆盖 6 种 centering：A/C/F/I/P/R；B 面心不出现在标准
   ITA setting）运行
   官方 `iso`，逐行比对：
   - 子群号、方向标签一致；
   - `Size == |det W|`；
   - `|det Basis_官方| == Z(子群)·Size/Z(母群)`；
   - `w_机器 · P − w_官方` 是母群格矢量（1 行已知例外，见 §3）。

当前结果：`oracle rows checked: 42`，全部通过。其中 1 行的 origin 比较被显式
allowlist 跳过（`SG 230 GM5+ P1`，见 §3），脚本会在结论行分别报告子群/Size/Basis/
标签的比较行数与 origin 的比较行数，不会把豁免行混进"全部通过"。规范表述是
"6 种 centering"（A/C/F/I/P/R）；脚本里定义的 B 面心在标准 ITA setting 中不出现，
因此没有对应用例。

## 6. 已知未决项（对抗性审查发现，未修复）

1. **官方 Origin 与存储 origin 的换算约定未钉死**（§3）：约 18.6% 记录不满足
   `w_prim·P ≡ w_printed (mod L_parent)`，偏移逐记录变化而非每 SG 常数。钉死它
   需要逐记录类测定 affine offset（对存储 origin 施加已知扰动、观察打印值），或
   改用 ISODISTORT 的等价定义。在此之前不得把 `w_prim·P` 说成"书中 Origin 列"。
2. **compound CIR irrep 无 oracle 覆盖**：官方对 SG199 `P1P1/P2P2/P3P3` 之类
   compound 标签打印空表（112 对 / 195 条记录），当前 21 组用例中没有 compound。
3. **抽样规模与分页**：gate 仅 42 行（0.28%）；程序分页上限 `PAGE ≤ 1000`，同一
   进程连续查询会被分页提示吞掉输入，全表验收必须按 (SG, irrep) 逐进程调用。

## 7. Rust API 位置

| 功能 | API |
|---|---|
| 按 (SG, irrep, 方向/标签/序号) 取子群 | `irrep::isotropy::isotropy_subgroup_for_direction` |
| 某 irrep 的全部子群 | `irrep::isotropy::isotropy_subgroups`（可用 `_at_k` 校验 k） |
| 磁子群 | `irrep::isotropy::magnetic_isotropy_subgroups(_for_direction)` |
| 几何换算 | `subgroup_size`、`parent_primitive_basis`、`basis_in_parent_conventional`、`origin_shift_in_parent_conventional` |
| 分导（恒等表示） | `identity_subduction`、`double_valued_subduction`、`format_identity_subduction` |
| 表格输出 | `format_isotropy_subgroups`、`format_magnetic_isotropy_subgroups` |
