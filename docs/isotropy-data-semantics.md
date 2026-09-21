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

## 3. Origin 的编码、坐标系，以及官方打印约定（已钉死）

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

**结论（对抗性审查第二轮，2026-09-20，已由主线程独立复现）**：官方 Origin 列不是
"另一个物理位置"，而是**同一存储 origin 在程序当前 ITA setting 下的表达式**。
`SET I` 切换母群/子群的 origin choice；出厂默认是"所有空间群 origin choice 2"，
而 pinned 数据表记录在混合 setting 下——230 个母群中 189 个是 choice 1，SG 227/228
一类是 choice 2。这个差值正好解释了整表的 18.6%。

复现（主线程已跑）：

```bash
cd isotropy_subgroup
printf 'PAGE 1000\nSC 250\nSET I ALL OR 1\nVALUE PARENT 139\nVALUE IRREP M1-\nSHOW SUBGROUP\nSHOW BASIS\nSHOW ORIGIN\nSHOW SIZE\nSHOW DIRECTION\nDISPLAY ISOTROPY\nQUIT\n' | ISODATA=$PWD/ ./iso
# -> 126 P4/nnc 2  P1  (1,0,0),(0,1,0),(0,0,1) (1,1,1)  == 存储 (2,2,2)·P 逐位相同
# 去掉 SET I ALL OR 1（出厂默认）同一行打印 (1/4,1/4,1/4)
```

- `SET I ALL OR 2` 的输出与默认逐行相同（15035/15035），即默认就是 choice 2。
- `SET I ALL OR 1` 下全表 13978/15035 逐位相同（92.97%）、14040/15035 模母群格相同
  （93.38%）；189/230 个母群完全一致，且不存在"默认一致但 OR1 不一致"的母群。
- `SHOW ELEMENTS`（需先用 `VALUE DIRECTION <lab>` 选方向）在默认与 OR1 下打印
  **模母群格相同的同一组母群操作**，而 Origin 列移动 `(3/4,3/4,3/4) ∉ L_parent`：
  该列是 setting 标签而非位置，所以任何纯坐标帧换算都不可能普遍复现默认输出。
- 已排除的假设：`w · B_printed`（SG 139 `P1` 的 `(a,0)` 与 `(a,a)` 存储
  `(basis,origin)` 完全相同却打印不同 origin）、domain/arm 差异、`*_old` 字段。
- 残余：约 995 条 / 41 个母群来自**其它 ITA setting**（单斜与三方晶系的 cell/axis
  choice）。其中 SG 227/228 的 545 条可用 `SET I 227/228 OR 2` 关闭；单斜/三方约
  450 条用 `SET I <sg> CELL k` / `AX RH|HEX` 未关闭，是唯一仍未钉死的一类。

契约：`origin_shift_in_parent_conventional` 给出**记录 setting 下**的值，等于官方在
`SET I ALL OR 1` 下的打印；官方出厂默认会打印另一个代表元（相差一个 ITA setting
变换）。`scripts/verify_isotropy_oracle.py` 现已显式运行 `SET I ALL OR 1`：42 行
抽样全部逐位相同、**无豁免**（此前依赖 mod L 谓词 + 1 行 allowlist）。

## 4. 分导（subduction）：`isotropy_subduce_*`

每条 isotropy 记录附带若干分导条目，四列并行：

- `isotropy_subduce_irrep`：母群 irrep 序号（1–4777，ISO 顺序）；
- `isotropy_subduce_frequency`：`i(G)`，该 irrep 的分导表示中包含子群单位表示
  （恒等表示）的次数；
- `isotropy_subduce_domain`：domain 序号；
- `isotropy_subduce_subgroup`：**已解开**——它是 1-based 的 isotropy 记录序号，
  它指向的那条记录的 `direction` 标签正是官方打印的 Dir 列（94271/94271 锚点与
  被分导 irrep 属于**同一母群 SG**；其中 79033 条落在同一 irrep 内，15238 条指向
  同一 SG 的另一个 irrep，这是合法的 compound/高对称共享情形；30/30 抽样逐字符
  复现官方输出）。API 暴露为 `IdentitySubduction::direction_label`。

这正是官方 `SHOW FREQ`（加 `DIR` 时附 `Dir(domain)`）的输出，例如 SG 221 `GM4+`
方向 `P1` → 子群 `#83 P4/m`：

```
官方:  83 P4/m  1 GM1+ P1(1), 1 GM3+ P1(3), 1 GM4+ P1(1)
API :  [("GM1+", 1, "P1", 1), ("GM3+", 1, "P1", 3), ("GM4+", 1, "P1", 1)]
```

**其它波矢（other-wave-vector）条目**同样属于这一行输出：`isotropy_w_subduce_*`
共 5756 条，覆盖 1006/15239 条记录，此前被整族丢弃；现在由
`other_wave_vector_subduction(ordinal)` /
`IsotropySubgroup::other_wave_vector_subduction()` 提供。回归样例：SG 225 `W5`
方向 `S60` → 9 条同 k 条目 + `3 DT5, 3 SM3, 3 SM4`。

**"双值/spinor"是错误命名（已纠正）**：`DT`、`SM` 是 SG 225 k 列表里的波矢标签
（Δ、Σ 线），这些条目是**同一母群 SG 在别的波矢上的单值 irrep**；用户复核取出
它们的纯二重旋转矩阵，全部满足 `D(C₂)² = +I`，与 spinor 语义不符。旧名
`double_valued_subduction` / `DoubleValuedSubduction` 已重命名为
`other_wave_vector_subduction` / `OtherWaveVectorSubduction`。

**边界**：本数据集只给“哪些母群 irrep 包含子群的恒等表示、重数多少”（Landau /
铁性分类所需），**不给**某个母群 irrep 分解成子群全部 irrep 的完整分导表示
（例如 Γ3+ ↓ P4/m = ?）。官方 `iso` 也没有这个能力：`SHOW FREQUENCY` 配
`DISPLAY IRREP` 给的是 **Wyckoff 位置**的诱导点群 irrep（手册 §SHOW FREQUENCY），
`SHOW COMPATIBILITY` 是 k 点兼容关系，`VALUE SUBGROUP` / `VALUE FREQUENCY` 只是
`DISPLAY ISOTROPY` 的过滤器。

`data_little.txt` 的 `little_subduce_*` 结构已经解开但**仍不可用**：

- 索引空间是 `little_irr_full_label` / `little_irr_old_map` 的 **10294 个紧凑
  little irrep**（不是本文档的 4777 个母群 irrep；两者交集为空 0/5517）；
- `little_subduce_irr_pointer[i]`（1-based，0 = 无数据，5517/10294 非零）给出该
  irrep 块的起始行，块大小是**每个空间群的常数**（SG221→14、SG225→12、SG230→8、
  SG1→1 …），与 k 和 irrep 无关；
- 每个块的最后一行恒为 `[(frequency, dim), pg_irrep=1)]`（5517/5517，dim =
  `little_irr_full_dim`），即分解到平凡点群 C1 的退化情形；
- 载荷只有 `(frequency, pg_irrep)` 对，`pg_irrep ∈ [1,12]`，**没有任何
  SG/basis/origin/direction 键**，因此无法与具体 isotropy 子群关联；
- 官方二进制没有任何命令打印该载荷（`SHOW COMPATIBILITY/STAR/MODES/KDEGREE` 都不
  打印），所以它也拿不到独立的 oracle 校验。

完整分解只能自行计算：用随包 PIR/CIR 矩阵限制到子群操作上（帧由
`isotropy_basis` + `isotropy_origin` 给出），并用 §4 的 94271 条恒等分导作为
trivial 列的回归 oracle。

磁子群表没有对应的 `mag_iso_subduce_*`，因此磁子群只能给出几何 + UNI/BNS，
给不出分导。

## 5. 验证 gate

`scripts/verify_isotropy_oracle.py` 会：

1. 用 pinned 数据复现每条记录（`Size = |det W|`、方向标签、子群号）；
2. 对 21 组 (SG, irrep)（覆盖 6 种 centering：A/C/F/I/P/R；B 面心不出现在标准
   ITA setting）运行官方 `iso`（显式 `SET I ALL OR 1`，见 §3），逐行比对：
   - 子群号、方向标签一致；
   - `Size == |det W|`；
   - `|det Basis_官方| == Z(子群)·Size/Z(母群)`；
   - `W·P_母群 == P_子群·B_官方` 作为**格**相等（`P_子群·B_官方` 给出打印胞的
     primitive 格；体积/行列式检查看不出基取向错误，这一条能看出来）；
   - `w_机器 · P − w_官方` 逐位相同（不再需要 mod L 或豁免）。

当前结果：`oracle rows checked: 42`，全部通过，且 origin 比较 **42/42 逐位相同**、
无任何 allowlist 豁免（脚本显式运行 `SET I ALL OR 1`，见 §3）。规范表述是
"6 种 centering"（A/C/F/I/P/R）；脚本里定义的 B 面心在标准 ITA setting 中不出现，
因此没有对应用例。

## 6. 已知未决项（对抗性审查发现，未修复）

1. **其它 ITA setting 的残余**（§3）：origin choice 已钉死，但约 995 条 / 41 个
   母群来自单斜与三方晶系的 **cell/axis choice**。SG 227/228 的 545 条可用
   `SET I 227/228 OR 2` 关闭；其余约 450 条（SG 3–41、63–68、151/152、178/179）
   用 `SET I <sg> CELL k` / `AX RH|HEX` 未能关闭。在这些记录上
   `origin_shift_in_parent_conventional` 仍是记录 setting 的值，与官方默认输出
   相差一个 setting 变换。
2. **compound CIR irrep 无 oracle 覆盖**：官方对 SG199 `P1P1/P2P2/P3P3` 之类
   compound 标签打印空表（112 对 / 195 条记录），当前 21 组用例中没有 compound。
3. **抽样规模与分页**：gate 仅 42 行（0.28%）；程序分页上限 `PAGE ≤ 1000`，同一
   进程连续查询会被分页提示吞掉输入，全表验收必须按 (SG, irrep) 逐进程调用。
4. **完整分导分解不在数据中**（§4）：`isotropy_subduce_*` 只给"包含子群恒等表示"
   的母群 irrep 与频率 i(G)；`little_subduce_*` 虽然结构已解开，但没有 subgroup
   键、与 4777 个母群 irrep 交集为空，任何官方命令都不打印它；官方 `iso` 也没有
   打印完整分解的命令（`SHOW FREQUENCY` 配 `DISPLAY IRREP` 给的是 **Wyckoff
   位置**的诱导点群 irrep）。要得到"母群 irrep → 子群全部 irrep + 重数"必须自行
   计算。
5. **磁 isotropy 表的空洞**（§4）：16721 条只覆盖 1421 个 UNI，**230 个 Type-II
   （grey）磁群完全没有记录**，`magnetic_isotropy_subgroups` 对它们只能返回空。
   另有 1482/16721 条磁记录的方向标签在同母群 irrep 的常表里不存在（例如 UNI 3
   `M1` 的 `C1` vs 常表 `P1, P3, S1`），所以**不能**按方向标签把磁表 join 到常表；
   可 join 的 15239 条中只有 8887 条 basis 矩阵逐项相同，其余 6352 条是**同一格的
   unimodular 换基**（`|det|` 15239/15239 相同），按矩阵相等 join 会产生假不匹配。
6. **112 个母群 irrep（195 条记录）的 ML 标签官方不接受**（SG23 `W1W1` vs 紧凑
   `W1WA1`，SG82 `P1P1` vs `P1PA1`）：`VALUE IRREP W1W1` 打印空表，紧凑拼写才有
   4 行且与存储数据一致。`IrrepRecord::ml` 因此在 oracle gate 之外还要看
   `little_irr_full_label` 的紧凑拼写。

## 7. Rust API 位置

| 功能 | API |
|---|---|
| 按 (SG, irrep, 方向/标签/序号) 取子群 | `irrep::isotropy::isotropy_subgroup_for_direction` |
| 某 irrep 的全部子群 | `irrep::isotropy::isotropy_subgroups`（可用 `_at_k` 校验 k） |
| 磁子群 | `irrep::isotropy::magnetic_isotropy_subgroups(_for_direction)` |
| 几何换算 | `subgroup_size`、`parent_primitive_basis`、`basis_in_parent_conventional`、`origin_shift_in_parent_conventional` |
| 分导（恒等表示） | `identity_subduction`、`other_wave_vector_subduction`、`format_identity_subduction` |
| 表格输出 | `format_isotropy_subgroups`、`format_magnetic_isotropy_subgroups` |
