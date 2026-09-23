# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Rust 风格化改造账本（2026-08-15 起）

目标：把 C 直译内核逐步改为 Rust 风格，消除哨兵、零填充缓冲、下标 panic
和不变量可被公开字段绕过的问题。每个阶段必须全量测试、clippy 零警告、git
提交，并通过独立 reviewer 的极简严格复审后才算完成。

已完成阶段（提交 `dad5dd8` → `ee9286a`；角色数据/summary 后续基线
`ba2d997`、`969a89c`、`daf04eb`）：

1. `Crystal::new` / `with_magnetic` / `SymmetryOps::from_parallel*` 全部返回
   `Result`，不再 panic；README 与 doctest 同步。
2. `Cell::set_cell` / `set_layer_cell` / `set_cell_with_tensors` / 新增
   `set_tensors` 全部返回 `Result`，先验证后写入；周期 setter 显式清除
   `aperiodic_axis`；`Crystal::to_cell` 对层状磁性晶胞不再二次覆写非周期坐标。
3. `Symmetry` / `PointSymmetry` / `MagneticSymmetry` 移除 `new(size)` 零填充与
   `truncate`，改为 `new` / `with_capacity` / `push`；磁空间群的两条提取路径
   共享 `deduplicate_spatial_operations`。
4. 点群识别去除 `number == 0` 哨兵：`get_transformation_matrix -> Option`、
   `get_pointgroup -> Result`、Hall 搜索返回 `Option<i32>`。
5. `Crystal` 字段全部私有化，提供只读访问器；等长不变量不再可被外部破坏。
6. 清理内核生产路径 `unwrap`/`expect`、重复 `Spacegroup` 组装与映射校验；
   `transform_from_primitive` 保持旧公开行为（Primitive/BFace/Error 返回 None）。

当前验证基线（每次变更都必须复跑）：

```bash
cd /home/liuyichen/TB_rs
CARGO_TARGET_DIR=/home/liuyichen/TB_rs/cryspglib/target \
  cargo test --release --package cryspglib --tests
CARGO_TARGET_DIR=/home/liuyichen/TB_rs/cryspglib/target \
  cargo test --release --package cryspglib --doc
CARGO_TARGET_DIR=/home/liuyichen/TB_rs/cryspglib/target \
  cargo test --release --package cryspglib --example audit_irrep_subduction \
  --example census_subduction_gaps --example probe_subduction_settings \
  --example line_family_coverage --example line_transport_ledger
CARGO_TARGET_DIR=/home/liuyichen/TB_rs/cryspglib/target \
  cargo run --release -p cryspglib --example line_transport_ledger -- --witnesses --gate
CARGO_TARGET_DIR=/home/liuyichen/TB_rs/cryspglib/target \
  cargo run --release -p cryspglib --example line_family_coverage -- --gate
CARGO_TARGET_DIR=/home/liuyichen/TB_rs/cryspglib/target \
  cargo clippy -p cryspglib --all-targets --release -- -D warnings

cd /home/liuyichen/TB_rs/cryspglib
python3 -m unittest discover -s scripts -p test_verify_isotropy_oracle.py
python3 -m unittest discover -s scripts -p test_check_other_wave_vector_rows.py
python3 -m unittest discover -s scripts -p test_subduction_settings.py
python3 -m unittest discover -s scripts -p test_build_table.py
python3 scripts/verify_isotropy_oracle.py
python3 scripts/check_other_wave_vector_rows.py
```

`--tests` 不运行 example 内的回归；上面的 audit、census、probe 三个 example
必须单独执行，不能只以 library/integration 测试通过代替门禁退出码与逐星清点验证。

当前基线（2026-09-23，R5 收口 + 复核处理 + R6.1 复核修正 + R6.2 参数族覆盖后，`-p cryspglib` 限定到本 crate）：
lib `409 passed / 4 ignored`，全部测试二进制（`--tests`，22 个）`570 passed / 0 failed /
4 ignored`，doctest `27 passed`，
example 审计回归 `20 passed`、缺口清点回归 `7 passed`、settings 探针回归 `2 passed`、
参数族覆盖回归 `1 passed`（91 s，全表 5,756 行 × 4 个网格点 + 18 个网格外支持检查）、
transport ledger 回归 `2 passed`（anchor 共轭律/范数/Bloch 协变 + 平移账本非空）；
严格 all-target clippy 通过（Cargo 仍报告既有
workspace manifest 警告）；isotropy oracle 离线测试 `9 passed`、真实 oracle
`62` 行 / `26` 个描述串 / `62` 个 origin 通过；其它波矢行门禁离线测试 `16 passed`、
pinned 数据 `checks_failed=0`（73 源 / 1,006 记录 / 5,756 行；73/73 源已解为参数化
直线 `k = Γ + t·v`，冻结 little 特征标表 73/73，引擎已算出全部 5,756 行，
参数约定 `t = 1/4`）；官方 live oracle 全量通过：w 行 5,756 = pinned 5,756、
`mismatches=0`；
settings 管线：`test_subduction_settings.py` 10 passed（离线核对 15,239 行覆盖、
身份、`|det U| = 分母³`、shift 最简；覆盖规则与 `--check` 共用
`check_ordinal_coverage`）、`test_build_table.py` 7 passed（汇编器拒绝空/截断/重复、
**新输出路径必须先知道预期全集**、`--out` 与 `--baseline` 必须一致）、
`generate_subduction_settings.py --check` 在线复推 75 条遗留记录通过、
`build_table.py --check` 可按 README 步骤 7 重建比对；第三方 reviewer 的
`gate_regression_tests.py` 4/4 通过；
全表审计（**三个门禁**：`--require-complete --require-w-complete
--require-full-decomposition`）**同时退出 0**、判词 `VERDICT complete scope=global
full_decomposition=complete`（`identity_rows=94271`、probe **366260/366260 完整分解、
恒等-only 0**、w `5756/5756`、`engine_errors=0`、`hard_failures=0`；R6.1 把 w 门禁从
Γ-only 路径升级为**完整分解**后墙钟 557.0–604.8 s，此前 Γ-only 口径 517.1–541.0 s；
同一二进制上的净代价实测约 +4.8 s / 5,756 行，墙钟差属机器负载、不是 A/B）；
R4 批 2a 后 `366039 + 221`，
批 1 后 `357033 + 9227`，R2 后 `353382 + 12878`，R2 前 `351547 + 14713`）。普通离散标量
覆盖在固定语料上已闭合（R5 验收清单：恒等正项 94,271/0 不匹配、Γ Frobenius
1,895/1,895、w 行 5,756/0 错误，全部满足）；缺口清点的门禁用法是
`census_subduction_gaps <audit.tsv> --require-empty`（要求至少一条 probe 行、每条都被
引擎回答、且无 identity-only 行；R5 审计上打印 `probe_rows=366260 identity_only_rows=0
closed=true` 并 exit 0）；空 scope 的审计（如 `--parent 2 --ordinal 0`）现在报
`empty scope` 并非零退出，不再冒充 `clean`。正式报告见
`docs/subduction-audit.md` 的「R5：普通离散标量覆盖闭合」小节。
注意**不要**在 workspace 根跑不带 `-p` 的 `cargo test --release`：sibling 成员
`Rustb` 当前自身编译失败（`ndarray_lapack.rs:23` E0259、`lib.rs:320` E0080 两个 BLAS
后端同时启用），与本 crate 无关，但会让整条命令以 exit 101 结束、0 个测试执行。
### 2026-08-26 有限域类型化

本轮按 `RUST_TYPE_MODERNIZATION_PLAN.md` 区分三类数据，而不是把所有有限数据库
字段机械地改成 enum：

- 小型稳定语义集合使用 `OperationKind`、`TimeReversalPolicy`、`TensorParity`；
  其中磁操作验证/组合层已改用 unitary/antiunitary 语义访问，原始 `bool` 只保留在
  兼容字段和数据库平行数组边界。
- 大型有界编号使用 `SpaceGroupNumber`、`HallNumber`、`UniNumber` 的 non-zero
  transparent newtype；均支持 `TryFrom<usize>`，非法值不能进入 typed lookup。
- Hall/国际/BNS/OG/k 点/irrep 符号仍是数据库文本，不建立数百或数千 variant 的
  巨型 enum。纯磁群元数据的 BNS/OG 返回值改借用生成表中的 `&'static str`，避免
  每次 `from_uni` 都分配两个 `String`；运行时派生文本仍保持拥有所有权。

兼容期保留整数和 `bool` 入口，它们会在边界验证后委托 typed API。新增代码应优先
使用 `from_hall_number`、`from_uni_number`、`from_space_group_number`、
`from_parallel_kinds` 和 `identify_with_hall_number`。

本轮 release 验证：library `222 passed / 0 failed / 3 ignored`，全部集成测试
（包括 `4479/4479` setting-aware round trip）通过，doctest `26/26`，严格
all-target clippy `-D warnings` 零警告；Rustb `0.7.2` 开启
`intel-mkl-system,cryspglib` 的 release check 通过。

---

## Isotropy subgroup 几何与分导（2026-09-20 起）

目标：给定 (空间群, k 点, 不可约表示, 序参量方向) 返回 isotropy subgroup，
并给出该子群在母群坐标系中的几何（基变换 + 中心点平移），以及母群 irrep 对该
子群的分导信息。数据全部来自 pinned `isotropy_subgroup/iso.zip`
（`data_isotropy.txt` 15239 条、`data_magnetic.txt` 16721 条）。

### 已落地的 API

- `irrep::isotropy::isotropy_subgroup_for_direction(sg, label, convention, direction)`：按
  方向描述串 `"(a,0,0)"`、ISO 方向标签 `"P1"` 或表内序号选取；
  `isotropy_subgroups[_at_k]`、`magnetic_isotropy_subgroups` /
  `magnetic_isotropy_subgroup_for_direction`。
- 几何换算：`subgroup_size`、`parent_primitive_basis`、
  `basis_in_parent_conventional`、`origin_shift_in_parent_conventional`。
- 分导：`IsotropySubgroup::identity_subduction()` /
  `irrep::isotropy::identity_subduction(ordinal)`，即官方 `SHOW FREQUENCY`：
  哪些母群 irrep 的分导表示包含子群恒等表示、重数 `i(G)`、domain 序号。
- 格式化：`format_isotropy_subgroups`、`format_magnetic_isotropy_subgroups`、
  `format_identity_subduction`（均输出母群 conventional 基下的几何）。

### 数据语义（已用官方二进制钉死，勿凭直觉假设）

- **标签查询须显式指定约定（2026-09-21 用户要求）**：普通／磁 isotropy 的
  label 查询现在接受 `LabelConvention::{Cdml,Bc}`，方向选择函数参数为
  `(sg, label, convention, direction)`，无隐式默认。`query::find_irreps`、
  k 点查询及 `SpaceGroup` bridge 同样显式选择。结果通过 `labels()` 同时给出
  `cdml` 与可选 `bc`；格式化子群表显示两者。CDML 的 Γ 前缀是 `GM`，BC
  的输出为 Unicode `Γ`，但**编号必须读实际 BC 表**：213 `X2` = BC `X1`，
  123 `GM2+` = BC `Γ3+`。BC 重名必须保留所有候选／报歧义，可用 k 坐标
  消歧；不能用第一条命中。12 条 `***` 和 3611 条由 ML 生成占位 BC 的
  spinor 没有可靠 BC 对应，不能冒充 BC 支持。标签选择不改变几何坐标帧。
  磁摘要同时输出 k 点及来源 irrep 的两套标签；compound 成分没有独立 BC
  映射时返回 `None`，不能继承整个 compound 行的 BC 标签。
- **k 点分组不可依赖 BC 是否存在**：同一坐标的 scalar/spinor 行必须留在
  同一物理 k 点，CDML 模式全表为 1350 点／8388 行；BC 模式为 1342 点／
  4765 行。独立 DSH 审查复现了按标签可用性分组会拆成 2692 点的问题；
  已按坐标分组修复，`tests/irrep_labels.rs` 的 12 项测试钉住全表分区、
  真正编号对应、重名、缺失映射与双标签输出。
- **书版 basis 锚点**：1988 书 Table 1 p.1-325（扫描 PDF 第366页），
  213 `X2` / `C23` →146 的 basis 是 `(2,-2,0),(0,2,-2),(2,2,2)`，
  origin=0；两处横线是负号。9.6.1 程序打印的六方基与书中满足
  `a_iso=b_book, b_iso=a_book, c_iso=-c_book`。详见
  `docs/isotropy-data-semantics.md` §2.1，测试钉住两种 convention，勿把 OCR
  丢失负号后的 `(2,2,0),(0,2,2),(2,2,2)` 当作书中常规胞。

- `isotropy_basis` / `isotropy_origin` 表达在**母群 primitive 胞**帧中，不是书里
  打印的 conventional 帧：官方程序打印子群的 ITA conventional 基。
  #16 P222 → #22 F222 的例子：机器 `(0,1,1),(1,0,1),(1,1,0)`（det 2），
  官方 `(2,0,0),(0,2,0),(0,0,2)`（det 8 = Z_F·size）。
- `Size = |det basis|`（**不乘**母群 centering 数）；`|det Basis_官方| =
  Z(子群)·Size/Z(母群)`。
- `isotropy_origin` 每条 4 整数 `(x,y,z,d)` = `(x/d,y/d,z/d)`，分母 ∈
  {1,2,3,4,6,8,12,16,24}；#167 R-3c 的 `(0,1/2,0)` 在六方 conventional 基下是
  `(-1/6,1/6,1/6)`，与官方输出一致。
- **origin 约定已钉死（对抗性审查第二轮，2026-09-20，主线程独立复现）**：
  官方 Origin 列不是"另一个位置"，而是**同一存储 origin 在程序当前 ITA setting
  下的表达式**。`SET I` 切换 origin choice；出厂默认是"所有空间群 choice 2"，
  而 pinned 数据表记录在混合 setting 下（230 个母群中 189 个是 choice 1，
  SG 227/228 一类是 choice 2），这个差值就是整表的 18.6%。
  决定性复现：`SET I ALL OR 1` 后 SG 139 `M1-` `P1` 打印 `(1,1,1)`
  = 存储 `(2,2,2)·P` 逐位相同；默认打印 `(1/4,1/4,1/4)`。
  `SET I ALL OR 2` ≡ 默认（15035/15035）；OR1 下全表 13978/15035 逐位相同
  （92.97%）、14040/15035 模母群格相同（93.38%），189/230 母群完全一致。
  `SHOW ELEMENTS` 证明该列是 setting 标签而非位置（默认与 OR1 打印**模 L 相同的
  同一组操作**，而 Origin 移动 `(3/4,3/4,3/4) ∉ L_parent`）。已排除
  `w·B_printed`、domain/arm、`*_old` 三种假设。
  残余：约 995 条 / 41 个母群来自单斜/三方的 **cell/axis choice**（SG 227/228 的
  545 条可用 `SET I 227/228 OR 2` 关闭，其余约 450 条未关闭）。
  契约：`origin_shift_in_parent_conventional` = **记录 setting** 下的值
  = 官方 `SET I ALL OR 1` 的打印。
- 分导表（`isotropy_subduce_*`）的 `isotropy_subduce_subgroup` 列已解开：它是
  1-based isotropy 记录序号，指向记录的 `direction` 就是官方 `SHOW FREQ DIR` 的
  Dir 列（94271/94271 落在对应 irrep 区间内，30/30 抽样逐字符一致）；API 暴露为
  `IdentitySubduction::direction_label`。双值（spinor）分导 `isotropy_w_subduce_*`
  （5756 条 / 1006 条记录）此前被整族丢弃，现由
  `other_wave_vector_subduction()` 暴露（回归：SG 225 `W5` 方向 `S60` → 9 条同 k 条目 +
  `3 DT5, 3 SM3, 3 SM4`；这些是同一母群 SG 在**别的波矢**（`DT`/`SM` 线）上的单值
  irrep，**不是 spinor**——用户复核并用 SG 225 的 k 列表 + `D(C₂)² = +I` 证实）。
- 分导表只给“母群 irrep 包含子群恒等表示”的频率，**不给**完整分导分解
  （例如 Γ3+ ↓ P4/m 的全部 irrep）。官方 `iso` 也没有该功能：`SHOW FREQUENCY`
  配 `DISPLAY IRREP` 给的是 **Wyckoff 位置**的诱导点群 irrep（手册 §SHOW
  FREQUENCY 原文），`SHOW COMPATIBILITY` 是 k 点兼容关系，`VALUE SUBGROUP` /
  `VALUE FREQUENCY` 只是对 `DISPLAY ISOTROPY` 的过滤器。`data_little.txt` 的
  `little_subduce_*` 结构已解开（索引空间 = 10294 个紧凑 little irrep；块大小是
  每 SG 常数；每块最后一行恒为 `[(dim,1)]`），但载荷只有 `(frequency, pg_irrep)`
  对、**没有任何 SG/basis/origin/direction 键**，与 4777 个母群 irrep 的交集为空，
  也没有任何官方命令打印它 → 无法与具体 isotropy 子群关联，完整分导必须自行计算
  （配方见下文"完整分导的实施配方"）。
- 磁 isotropy 表按**非磁母群的同一个 4777 个 irrep**索引，输出 UNI 1–1651；
  它不接受磁群 corep 作为输入，也没有分导表。

### 验证 gate

`python3 scripts/verify_isotropy_oracle.py`（先用
`unzip -o isotropy_subgroup/iso.zip -d isotropy_subgroup/` 解出 `iso` 二进制与
数据文件；脚本自行设置 `ISODATA`）用随包 `iso` 9.6.1 对 22 组
(SG, irrep)、覆盖 6 种 centering（A/C/F/I/P/R；B 面心不出现在标准 ITA
setting）的记录逐行比对子群号、方向标签、Size、`|det Basis|` 关系与 origin。
脚本现在显式运行 **`SET I ALL OR 1`**（见上一节 origin 约定），并且：
origin 比较是**强制逐位相同**（差一个母群格矢量会直接 FAIL，用户复核发现旧版把这种
情况只当诊断计数、exit 仍为 0）；新增**方向描述串比对**（拉 `SHOW DIRECTION VECTOR`
的官方列，逐行与生成表的 `direction` 比较，规范化 `;`↔`,` 与空白）。当前
`oracle rows checked: 48`、`descriptor strings checked: 26`，全部通过、无豁免。
SG177 `L1` 覆盖非立方 `C2` 与复分隔符；向量标签缺失、多余或重复时失败。
这仍只是 48 行抽样（0.31%）；全表 18.6% 的默认 setting 差异已解释（见上一节），残余约 450 条
单斜/三方 cell/axis choice 未关闭。

生成器的字节可复现性记录（可复查）：在干净树上
`python3 scripts/generate_irrep_data.py` 重新生成
`src/irrep/generated_data.rs`，其 md5 必须等于
`fbe341a9cfb1bd815632e812daae950d`（2026-09-21 清理 `w_subduce` 遗留 spinor 注释后的
当前值）。此次仅同步修正生成器与生成文件的 5 处文档注释，已比较确认所有非注释行
与 `acb8f7b` 完全相同，未重跑完整生成器。`acb8f7b` 的值是
`863c76358cf713724967d5667d1d0dfe`，当次重新生成只改 1080 条 `direction:` 行；
更早的值为 `4bbd6db9858117d346d8c93b81bcfbf9`（BNS 修复后）与
`f994cf4874440118ef2f872b4f1e021b`（原始）。

### 2026-09-21 复核与永久回归测试

- 复现并修复新 descriptor 门禁的漏检：方向向量输出为空、缺失或多余标签时此前
  可返回 exit 0；重复标签会被字典静默覆盖。现在全部拒绝。
- 新增 `scripts/test_verify_isotropy_oracle.py` 的 9 个离线测试；先确认缺失/重复标签
  与分隔符比较测试在旧实现失败，再修复到全部通过。origin 加母群格矢也必须失败。
- Rust 新增方向别名与嵌入测试：SG225 `GM4-` 的 `C1/C2` 都是 #8 Cm，但基矩阵不同；
  SG177 `L1 C2` 必须是 #21；空白、`;`/`,` 别名保留分量次序，旧错误分量串必须拒绝。
  其它波矢分导测试改为钉完整 formatter 输出段，并补越界负例。
- 清理 `other_wave_vector_subduction` 函数、生成器与生成文件中遗留的 spinor 注释；
  明确 `DT`/`SM` 前缀只能给波矢族，不能恢复参数化波矢的数值。

### 对抗性审查（2026-09-20）发现与修复

四个独立 reviewer（数据语义 / Rust API / 生成管线 / 声明审计）只读复核，全部
发现均已在主线程复现后处理：

1. **P0：磁 isotropy 方向标签取错数组**。生成器用 `mag_iso_orderparam_pointer`
   （指向 163000 条 `mag_iso_orderparam`）去索引 per-record 的
   `mag_iso_orderparam_label`（16721 条），造成 12112 条落到 `"dir<code>"` 兜底、
   合计 15724/16721 条标签错误，并让 `magnetic_isotropy_subgroup_for_direction`
   选中错误方向的子群。已改为按记录直接取标签（该数组 0 空值、每个 irrep 内
   0 重复），并加两条门禁：per-record 长度校验 + 禁止合成标签；回归
   `magnetic_direction_selection_uses_the_per_record_label`（SG 9 `L1`：`P1`→UNI 48，
   `P3`/`C1`→UNI 3，基矢各异）与 `magnetic_direction_labels_are_unique_within_each_irrep`。
2. **P1：origin 契约过强**（第一轮结论，第二轮已反转为"约定已钉死"）：第一轮曾把
   `w_prim·P` 与官方默认输出的差异记为"未钉死约定"，第二轮用 `SET I ALL OR 1`
   证明差异就是 ITA **origin choice setting**，文档与 oracle gate 均已按新结论重写
   （见上一节）。第一轮"逐记录变化、非每 SG 常数"的措辞是错的：默认 setting 下
   它正是每 SG/每子群组常数，等于对应的 ITA setting 变换。
3. **P1：`subgroup_size` 溢出**：第一轮把 i32 行列式改成 **i64** 并不充分——i64
   本身对极端 i32 输入也会溢出（见证 `diag(2^21,2^21,2^22)`，真值 `2^64`：release
   返回伪 `SingularSubgroupBasis{det:0}`、debug panic）。第二轮改为 **i128** 行列式，
   `SubgroupSizeOverflow { determinant: i128 }` 现在携带真实行列式，并加回归
   `subgroup_size_is_total_for_extreme_bases`（含 `i32::MIN` 与极大奇异基）。
   同类加固：`k_vectors_agree` 对 `d = 0` 返回 false；
   `origin_shift_in_parent_conventional` 对 `d ≤ 0` 返回 `InvalidOrigin`；
   `format_*` 不再用 `unwrap_or(零)` 吞掉换算错误。
4. **P1/P2：分导表不完整**：`isotropy_subduce_subgroup` 即 Dir 列锚点（已暴露为
   `IdentitySubduction::direction_label`）；双值 `isotropy_w_subduce_*`
   （5756 条 / 1006 条记录）此前整族缺失，现由 `other_wave_vector_subduction()` 提供
   （语义见下：它们是**别的波矢**上的单值 irrep，不是 spinor）。
   （回归 SG 225 `W5` `S60`：9 标量 + `3 DT5, 3 SM3, 3 SM4`）。
5. **生成器门禁加固**：子群基行列式集合 {1,2,3,4,6,8,16,32}、origin 最简分数、
   分导条目的 parent SG 一致性（94271 条全部满足）、ISO→IRREPS 索引翻译改为
   双射 + (SG, 标签) 恒等校验、pointer 表与 per-record 数组的长度校验。
6. **文档纠错**：`SHOW FREQUENCY` 并非可运行命令（官方为 `SHOW FREQ [DIR]`）；
   oracle 覆盖为 6 种 centering（B 面心不在标准 setting 中）；`direction_label`
   前缀不是 `i(G)` 的简写；`isotropy.rs` 的 `origin_shift()` 改为返回 `Option`。

### 对抗性审查第二轮（2026-09-21）发现与修复

第二轮四个 reviewer（磁标签与分导 / origin 与文档 / 声明与测试质量 / 全表内部
一致性）针对第一轮修复本身复核。主线程逐条复现后处理：

1. **P1：`subgroup_size` 的 i64 行列式仍会溢出**（见上第 3 条）→ 改 i128 + 回归。
2. **origin 约定钉死**（见上第 2 条）→ 文档/oracle/Rust doc 全部改写；
   `scripts/verify_isotropy_oracle.py` 增加 `SET I ALL OR 1`、删除唯一 allowlist，
   42/42 origin 逐位相同。
3. **声明精度**：commit message 里"21 个既有数组全部逐字节不变"不成立——
   `MAGNETIC_ISOTROPY_SUBGROUPS` 正是 P0 修复对象（字段级：`direction` 变
   15724/16721，其余字段 0 变化）；正确表述是"20/21 不变，1 个为修复而变"。
4. **测试质量**：`ordinals_address_their_own_records` 与 `:267` 的比较是
   `T[i]==T[i]` 恒真；`double_valued_subduction_tables_tile` 只比长度；
   若干测试只断言 `!is_empty()`/范围。已改为**钉住具体值**：新增
   `pinned_ordinals_address_the_expected_records`（SG 16/139/221 的 ordinal 与子群号
   对照 `data_isotropy.txt`）、双值分导逐项比对 `IRREP_W_LABELS`/`..._SPACE_GROUP`、
   HNF 比较复现官方打印的 `diag(2,2,2)`、formatter 断言整行文本、
   `subgroup_size` 极端基回归。
5. **缺失的负例**：`SubgroupIndexOutOfRange`/`InvalidOrigin`/`SingularSubgroupBasis`/
   `SubgroupSizeOverflow` 此前无 variant 断言，现全部 `matches!` 断言（含 payload）；
   `DirectionAmbiguous` 在 pinned 数据下**不可达**（0 重复标签/描述串），文档已注明，
   并新增直接驱动共享选择器 `select_unique` 的单元测试，使该分支被真实执行。
6. **`IsotropyDirection::Descriptor` 的来源**（第三轮已修正，见下）：`dim = 3` 的
   `P1/P2/P3/C1/C2/S1` 现在直接取官方 `SHOW DIRECTION VECTOR` 的字符串；
   `P2/P3` 原先映射成 `(0,a,0)`/`(a,a,0)` 是错的（用户复核：221 `GM4+` 的
   `(a,a,0)` 应为 #12 C2/m、`(a,a,a)` 应为 #148 R-3），`C2` 还要按晶系分成
   立方 `(a,a,b)` 与三方/六方 `(a;b;a)`。`dim = 2` 与 `dim ≥ 4` 仍是 cryspglib
   内部记法（官方字符串按 irrep 变化，无法用 `(dim, free, label)` 表复现）。
   匹配时规范化 `;`↔`,` 与空白；磁记录只有 ISO 标签，`Descriptor` 在其上恒返回
   `DirectionNotFound`，均有回归固定。
7. **交付命令**：workspace 根不带 `-p` 的 `cargo test --release` 会因 sibling
   member `Rustb` 自身编译失败而 exit 101、0 测试执行；基线命令一律用
   `-p cryspglib`（见文件开头）。

### 对抗性审查第二轮补充（2026-09-21 晚）

第一轮 reviewer 的后续复核又发现/澄清了以下几条，均已复现处理：

1. **P0：`mag_bns_label` 词法误解析**（`parse_labels` 同时抓两种引号，把
   `"P-1'"`、`"Ia'-3'd'"` 里的撇号当成单引号定界符）→ 该节解析出 2314 个 token
   而不是 1651，索引整体位移，**15546/16721 条磁记录携带错误 BNS 符号**
   （UNI 6 打出 `'" "P_S1        " "P-1         " "P-11'` 而不是 `P-1'`）。
   已改为"双引号优先"的单一 alternation，并加长度门禁（必须 1651）、
   单词 token 门禁（无引号/空白）与**钉值回归**（UNI 3/6/1001/1651 →
   `P_S1` / `P-1'` / `P4/m'mm` / `Ia'-3'd'`）。重新生成后只有 `bns_label` 行变化，
   新 md5 见上一节；其它 20 个数组逐字节不变。
2. **`irrep_label_pg` 其实存在且精确**（4777 项，与 `irrep_label` 平行；仅在 Γ 块
   的 1313 条非空，与 `smodes_sample.out` 一致）。本文档早先写的"点群标签表缺失"
   是错的；准确的限制是：它按 (母群 SG, 母群 irrep) 索引，**无法**给
   `little_subduce_pg_irrep` 这类任意点群序号命名。旧 reader 会把它误解析成
   4793 项（同样被（1）修复），但该节目前无消费者，故没有 shipped 数组受影响。
3. **"C 心 primitive 基取向错误"是假阳性**：reviewer 用**primitive** 换算结果去比
   官方打印的**conventional** 基。正确关系是 `W·P_parent == P_sub·B_printed`
   （`P_sub` 在左，因为 `P_sub` 的行是打印基的组合系数；只有对称的 F 矩阵才与
   右乘一致）。修正后 A/C/F/I/P/R 六类全部通过；`verify_isotropy_oracle.py` 现在
   增加了**格相等**检查（把打印胞乘上子群 centring 后与换算基比较），体积/行列式
   检查看不出的转置错误会被它抓住；`parent_primitive_basis` 的 7 个矩阵也在 Rust
   测试里逐个钉住，并断言 `|det P| = 1/Z`。
4. **分导锚点不是"同一 irrep"**：94271/94271 锚点属于同一**母群 SG**，但其中
   15238 条指向**同一 SG 的另一个 irrep**（compound/高对称 irrep 合法共享子群），
   只有 79033/94271 落在同一 irrep 内；生成器门禁与测试都按"同一 SG"表述。
5. **磁 isotropy 表的空洞与陷阱**（16721 条）：230 个 Type-II（grey）UNI **完全
   没有记录**，所以 `magnetic_isotropy_subgroups()` 对它们只能返回空；
   1482/16721 条磁记录的 `direction` 在同母群 irrep 的常表里**不存在**
   （例如 UNI 3 `M1` 的 `C1`），因此不能按方向标签做磁↔常表 join；
   可 join 的 15239 条里只有 8887 条 basis 矩阵逐项相同，其余 6352 条是**同一格的
   unimodular 换基**（`|det|` 15239/15239 相同），按矩阵相等 join 会产生假不匹配。
6. **112 个母群 irrep（195 条记录）使用官方二进制不接受的旧 ML 拼写**
   （SG23 `W1W1` vs 紧凑 `W1WA1`，SG82 `P1P1` vs `P1PA1`）：`VALUE IRREP W1W1`
   打印空表，改用 `little_irr_full_label` 的紧凑拼写才有 4 行且与存储数据一致。
   `IrrepRecord::ml` 因此不保证是官方可接受的标签（文档已注明；oracle gate 的
   21 组用例不含它们）。
7. **`little_subduce_*` 的结构已经解开，但仍不可用**：索引空间是
   `little_irr_full_label` 的 **10294 个紧凑 little irrep**（不是 4777 个母群
   irrep，两者交集为空 0/5517）；`little_subduce_irr_pointer[i]` 给出该 irrep 块
   的 1-based 起始行；**块大小是每个空间群的常数**（SG221→14、SG225→12、
   SG230→8、SG1→1 …）；**每个块的最后一行恒为 `[(dim,1)]`**（dim =
   `little_irr_full_dim`，即 C1 恒等表示的分解）。载荷只有 (frequency, pg_irrep)
   对，**没有任何 SG/basis/origin/direction 键**，官方也没有任何命令打印它
   （`SHOW COMPATIBILITY/STAR/MODES/KDEGREE` 都不打印载荷），因此它无法与某个
   isotropy 子群关联 → 完整分导必须自行计算。

### 完整分导（母群 irrep → 子群全部 irrep + 重数）的实施配方

逐步执行以 [完整实现任务卡](docs/full-irrep-subduction-plan.md) 为准：该计划补充了
严格 data-Hall 来源、typed character 空间、精确分母、嵌入规范与磁母群输入约定。
下面的折叠 k 简式不能替代已验证仿射变换 `x_G = T x_H + o` 所给出的 `k_H = T^T k_G`。

数据里**没有**这张表（见上条第 7 点与 §4），必须自算。用户复核后确认第一版配方
不充分，修正如下（每一步都要有 oracle 或全表门禁）：

1. **母群操作必须取 ISO/data-Hall 帧**：`SymmetryOps::from_sg(sg)` 返回的是**第一个
   Hall setting**，不保证与 ISOTROPY 数据帧一致；要经
   `crate::irrep::generated_data::SG_DATA_HALL[sg]` 走 `canonical_hall_ops` 一类路径，
   并像 irrep/corep 现网代码那样做 Hall↔ISO 的字符重排。
2. **子群 H 在母群帧中的嵌入不能只靠 (W, origin)**：221 `GM3+` 的 `P1`→#123 与
   `C1`→#47 都是 `W = I`、`origin = 0`，但操作数分别为 16 与 8；`isotropy_basis` /
   `isotropy_origin` 只给格与原点，**不给点群部分**。必须把目标子群自己的 ISO 操作
   （按 `sg_H` 的 data-Hall）经 (W, origin) 变换进母群帧，再用**完整操作集包含 +
   SG 识别**验证嵌入（官方 `SHOW ELEMENTS`，需先 `VALUE DIRECTION <lab>`，可作 oracle）。
3. **setting 变换要自己扩展**：`irrep::wigner::find_setting_transform` 只在有限范围的
   unimodular 换基里搜索、且恒等基成功即提前返回，既不能承担 `det W ≠ 1` 的超胞嵌入，
   也不做标签规范化；要么扩展成"一般有理基 + 任意原点，候选必须通过完整 Seitz 集
   验证"，要么直接在母群帧内构造嵌入（第 2 步）。
4. **折叠 k 与字符**：`k_H = P_H⁻¹ W P_G k_G`（`P_* = parent_primitive_basis(*)`，全部
   有理数运算），在 `query::irreps_of(sg_H)` 里用 `k_vectors_agree` 找折叠 k；
   重数 `n_α = (1/|H_k|) Σ_{h∈H_k} χ_Δ(h) conj(χ_α(h))`。注意 `W = I` 且 Size = 1 时
   `k_H = k_G` 只说明 k 相同，不代表方向/子群相同。
5. **方向输入不要依赖 descriptor**：`dim = 2` 与 `dim ≥ 4` 的描述串是 cryspglib 内部
   记法（见下一节），引擎应接受 ISO 标签/记录序号，或在引擎内部自己算不变量子空间。
6. **现成回归 oracle**：`isotropy_subduce_*` 的 94271 条恒等分导必须与算出的"子群
   恒等表示重数"逐条一致（最强 gate）；另加 Frobenius 恒等式
   `Σ_G i(G)·dim(G)/k_G = |P_parent|/|P_sub|`（1895 条 Γ 记录精确成立，compound 标签
   按 ML 分量数 `k_G` 折算，否则 3543 行会被重复计数）；再以用户指定的
   **221 `GM4+` 凝聚、查询 `GM3+`** 作为首个端到端验收例。
7. **多臂分阶段、磁共表示随后**：任务 8a 已用归档 CIR 完整矩阵核对普通标量
   full-star 求值与折叠几何；任务 8b 接入普通标量母群重数，8c 扩展 compound 与
   realification 的 k/-k。磁表没有分导列，
   不能把普通群字符检查当作磁共表示的验收。

### 完整分导的任务 1-6 落地状态（2026-09-21）

按 `docs/full-irrep-subduction-plan.md` 逐任务推进，每个任务一次提交、显式文件列表：

- 任务 1 `56bc4ba`：`docs/subduction-conventions.md`（坐标/表示/来源/精度契约，含手算例）。
- 任务 2 `6a81a91`：`scripts/verify_isotropy_operations.py` + 10 用例 / 70 操作 fixture
  （`tests/data/isotropy/operations.{json,txt}`）+ 离线测试。
- 约定修正 `f5d67cd`：官方打印与 `ipoint_op` 的旋转是**行作用** `x' = x M`，引擎是列作用
  `x' = R x`，解码必须转置；判别证据是 SG 167 `GM3+` P1 → #15（唯一转置不闭合的记录）。
- 任务 3 `a9425ae`：`src/irrep/subduction.rs` 精确有理仿射/格层（T/o、L_G/L_H、
  `(L^T)^-1` 坐标映射、严格 `SG_DATA_HALL` 来源，230 SG 的 4425 个操作都在 1/12 网格上）。
- 任务 4 `bae3edb`：`SubgroupEmbedding`；冻结元数据键是 **(parent, subgroup, U, δ)**——
  只按子群号建键会让另一个母群“验证通过但把 irrep 标签配错”，这次由 stored
  identity-subduction 频率交叉检查抓出来。`δ` 的例子是 #126 需要 (1/4,1/4,1/4)
  （isotropy 记录与 Hall 表用了不同 ITA origin choice）。
- 任务 5 `50f18c3`：第一个完整 Γ 分解 `subduce_irrep(&h, probe_ml)`；
  黄金用例 221 `GM4+` P1 → #83，`GM3+` ↓ = `GM1+ ×1 + GM2+ ×1`（2 = 1+1）。
- 任务 6：compound 行按 `CompoundMetadata` 语义展开成复不可约成分
  （`DistinctComponentSum` → 两个 CIR 成分；`ConjugateRealification` → seed 与其共轭，
  各自独立重数），Gram/维数和/逐操作重建/来源去重都按复不可约语义检查。
  例：221 `GM4+`（3 维）↓ #83 = `GM1+ ×1` + compound 行 `GM3+GM4+` 的两个成分各 ×1。
  **复核修正**：原先的 `check_assembly` 用 `compound_selected_arm_view()` 的
  `block_trace` 去比对由同一对 constituent 组装的字符之和，是循环论证，已删除
  （连同只为它存在的 `InconsistentCompoundRow`）。独立证据改为
  `tests/compound_subduction_regressions.rs` 的物理模型：SG 83 `GM3±GM4±` 的
  `(Rxx+Ryy)` / `det(R)·(Rxx+Ryy)` 迹与 SG 23 `W1W1` 的平移本征值 −i/+i。

**Γ 全表清点（任务 6 验收要求）**：`frobenius_reciprocity_matches_the_stored_identity_subduction`
任务 8c 扩充 fixture 后给出
`Γ 记录 1895 | 冻结子群命中 243 | 已钉住 14 | 多候选歧义 199 | 搜索空间外 30 | 交叉检查错误 0`。
历史数字 1895 的准确含义是**Γ 凝聚记录条数**，不是“恒等式已验证的条数”；
任务 8c 当时只对 14 条被 oracle fixture 钉住的记录验证了 stored `i(G)` 与 Frobenius
`Σ_D dim(D)·mult(trivial_H, D|H) = [G_k : H_k]`，其余 229 条显式报
`AmbiguousEmbedding` 或 `NoValidEmbedding`，不猜。任务 9 负责从 oracle 生成全表
逐记录 setting 元数据（含一般 unimodular/shear 候选）。该测试保留旧子群范围作为
回归；当前真正全表的分母及已完成条数以任务 9 审计为准（见下文）。

### 分导任务 8a：完整星适配器（2026-09-21）

- `src/irrep/subduction_star.rs`（`irrep::subduction::star`）：普通标量的完整星诱导
  字符求值、完整 transporter 校验、子群倒格同点合并与子群 star 分组。DSH 实现，
  Codex 独立验证并集成；compound/spinor 仍显式拒绝。
- 全 4105 条普通标量记录已核对星完整性和 `χ(E) = selected_dim × arms = record.dim`。
  这不是全表分解验收，也不是全表非恒等操作的字符核对。
- 独立源矩阵 gate：`scripts/generate_subduction_star_fixtures.py` 从 checksum-pinned
  CIR 完整矩阵取迹（不调用诱导算法）；10 个记录、308 个操作及各自四种格平移，
  `tests/subduction_star_source.rs` 共 1232 次字符比较。另保留原有 59/538 恒等项 gate。
- 普通母群的小群分解和完整星重建见任务 8b，compound 母群及 realification 的
  k/-k 见任务 8c。现有
  `subduce_irrep_with_embedding` 对多臂仍返回 `UnsupportedMultiArmStar`。
  准确范围及复现方式见 `docs/subduction-conventions.md` §10。

### 分导任务 8b：普通标量完整星分解（2026-09-21）

- `irrep::subduction::star::decompose::subduce_full_star_with_embedding` 返回每个子群
  star 的精确 q、源行 k、star size、小表示维数、重数和真实 CIR 来源号。
  8b 支持普通母群及普通/`DistinctComponentSum` 子群目标；8c 扩展 compound 与
  `ConjugateRealification`。spinor、磁共表示仍显式拒绝。
- 对每个子群星搜索有数据的代表臂，只把折叠到该 q 的母群臂在 `H_q` 上求迹，
  复用字符内积求解器；随后用子群自身的行和 Hall 操作诱导所有目标，逐操作重建
  母群完整星限制。检查 `Σ multiplicity × little_dim × star_size`，缺数据不部分返回。
- Codex 独立源矩阵验收：13 个完整分解钉值、92 个嵌入操作上的原始矩阵迹。
  包括 221 → #12 的 `X5+` 含 `V1+ ×2`，以及 Size=2 的 139 → #126，
  `X1+` → 二维 `M1`、`N1+` → 二维 `R1`（star size 2）。
- 五个指定上下文共有 161 个普通探针全部成功：99 个与旧单臂入口一致、62 个多臂；
  79 个 spinor/compound 探针明确不在该清点的支持范围。负例覆盖错配上下文、
  删除 star/q/arm 与缺数据；SG 167 `F1+` 钉住必须搜索第二个代表点的分支。
- 扩展 stored-frequency gate：59 个冻结嵌入、2060 次普通探针对照（458 个正项、
  1602 个零项），其中 1522 次为非 Γ 探针。15 个缺数据组合明确钉住：
  ordinal 13345/13346/13351 的 `W1`–`W5`，折叠 q 不在 #8 的离散表中；
  不计入成功分解。原有 59/538 Γ 对照保留。
- 这仍不是任务 9 的逐记录 setting 覆盖或全表 94271 条恒等项验收。
  详细契约见 `docs/subduction-conventions.md` §11。

### 分导任务 8c：compound 完整星及 k/-k（2026-09-21）

- DSH 负责实现，Codex 负责来源审计、独立矩阵 fixture、集成验收与提交。
  `star::scalar_star::ScalarStar` 按冻结来源展开复成分，保留 CIR 身份和共轭标记；
  `OrdinaryStar` 继续只接受普通记录。共轭作用于完整求值（含 Bloch 相位），折叠时
  重合臂保留各自成分身份。完整星分导入口复用同一 Gram/重数/重建求解器。
- 子群 realification 成分用完整 Seitz 操作搬到共同 q。仅对同一来源配对、单位
  范数且逐操作相同的 seed/共轭合并；母群始终保留两份。#19 `R1R1` 自限制是
  二维 CIR 559 重数 2；#45 `W1W1`/`W2W2` 是小维数 1、星大小 2、重数 2；
  #23 `W1W1` 是 k/-k 两个星，各重数 1。167 `GM3+` P1 → #15，`T1T2` → `M1 ×2`。
- Codex 复核收紧等价门禁：内积接近 1 不足以合并，必须逐操作相同；Gram 对角元
  的误差阈值与共享求解器一致。永久反例 `near_unit_overlap_does_not_replace_pointwise_character_equality`
  固定了内积误差约 2.5e-9、逐项误差约 1e-4 时仍须拒绝的情形。
- 返回 `parent_source_identity()`；`parent_irnumber()` 对普通记录为 `Some(id)`，
  compound 为 `None`。`stored_k()` 指复成分有效 k，允许是源行的 -k。
- `generate_subduction_compound_fixtures.py` 从 checksum-pinned CIR 原始完整矩阵
  取迹；17 来源、12 compound 的永久测试包含 928 次逐复成分、464 次物理和、
  6 个分解钉值和 28 次完整重建。#23 中心化平移的 -i/+i 分别钉住，防止总和
  抵消掩盖相位反号。
- 九个固定上下文逐探针清点：216 个普通、24 个 compound 完整分解成功，无缺失；
  157 个 spinor 未支持且未计为成功。四个自身嵌入的 78 个标量探针还逐项核对
  自限制保持来源身份及正确重数。
- 全 672 条来源/臂表达式审计：519 distinct，112 异星 realification，41 同星
  realification。41 条同星来源另以原始完整字符及有限平移相位类作 1484 次共轭
  等价检查；这不等于全部 672 条非恒等字符的生产实现都已核对。
- 四个官方自限制嵌入 (19,19)、(23,23)、(45,45)、(83,83) 使用已核实的 I/零原点；
  操作 oracle 为 14 用例 / 90 代表。无通用 identity setting 绕过；#23 `W1W1` → #22
  官方不打印操作行，未擅自冻结该嵌入。
- 新增 stored-frequency gate 在 69 个冻结嵌入中比较 78 次 compound 恒等项，均为
  零项，64 次非 Γ，无缺失。普通 59/2060、Γ 59/538 的既有 gate 保留。
  下一步是任务 9 的逐记录 setting 与全表验收；磁、spinor、离散表缺失的 k 仍未实现。
  详细契约及覆盖边界见 `docs/subduction-conventions.md` §12。
- 本轮验收：release lib `380 passed / 4 ignored`、integration `117 passed`、doctest
  `27 passed`；严格 all-target clippy 通过（workspace 既有 manifest 警告仍在）。
  Python 离线 `42 + 2` 测试、几何 oracle `48 行 / 26 描述串`、操作 oracle
  `14 用例 / 90 代表`、两套 CIR fixture 重生成一致性检查全部通过。

### 分导任务 9 第二轮：embedding 全覆盖（2026-09-21 深夜）

第一轮结束时 95 条记录被引擎拒绝、1,698 条 hard_failure。第二轮把它们压到
**embedding 15,239/15,239（0 拒绝）、hard_failures 771**。两处根因：

1. **引擎 bug（本轮最重要）**：`SubgroupEmbedding` 用 `W·P_parent` 当子群平移格，
   而它实际校验的仿射映射是打印基 `B = P_sub^-1·U^-1·W·P_parent`。两者只有在 `U`
   幺模时相同；第一轮为少量单斜记录冻结了**分数** `U`，于是这些记录用错误的格去约化
   陪集代表元、把明明复现官方胞的约定判成失败。改为「格 = 候选自身映射下子群格的像」
   `U^-1·W·P_parent` 后：拒绝 95 → 20、**engine error 760 → 0**、hard_failures
   1,698 → 934。幺模 `U` 下这是同一个格，所以既有结论不变。
2. **child shift 的第二条推导路线**：`derive_shift_both.py` 同时尝试两种帧差，各自
   独立送引擎校验——`-(B^T)^-1 (o_OR1 - o_OR2)`（ITA origin choice，覆盖 2,133 条）
   与 `(B^T)^-1 (o_printed - o_stored)`（记录帧差，额外 172 条，正是单斜与紧凑标签
   的残余）。接受 2,305/2,327；其余用引擎校验过的搜索值，19 条显式 override。
   结果 embedding 15,238/15,239，最后一条（ordinal 283）用第二条路线的一次
   4-候选读数也通过（`(B^T)^-1(o_p - o_s)` 正号），全表 15,239/15,239。

全表审计：存储恒等正项 93,720 通过 / 391 不匹配 / 0 不可用；Γ Frobenius
1,895/1,895；absent_positive 380；embedding 15,239/15,239；仍缺 14,713 条
（折叠后的子群 k 不在随包离散表里）与 5,756 条 w 条目缺 k 参数。

**剩余 771 条已定位为「不是 setting 问题」**：`examples/trace_embedding.rs` 打印
单条记录的格、变换与全部映射陪集代表元；对 ordinal 12471（SG 221 `X3+` `P2` → #125，
stored 频率 1、引擎 0）**16 个映射代表元与官方 `SHOW ELEMENTS` 搬进母群帧后模母群格
逐项完全相同**，且 16 个旋转互不相同。所以嵌入是对的，分歧在子导块的
**恒等重数求值**（`parent_character_of` / `solve_character_block`）。所有失败 ordinal
的 `size > 1`（2:105、4:55、8:20、32:16、6:2），size = 1 的 Γ 记录全部通过
Frobenius 1,895/1,895，指向超胞/折叠 k 路径上「代表元带非零格平移时的 Bloch 相位」。

**决定性判据（下一轮入口）**：198 条失败 ordinal 中有 **178 条是在自己的凝聚 irrep 上
失败的**——存储表说凝聚 irrep 含子群恒等表示 1 次，引擎算 0。这一条不需要 oracle
就能判定是 bug：某方向之所以以 D 为凝聚 irrep，正是因为 D 的分导包含子群恒等表示。
所以问题在 `subduce_full_star_with_embedding`（完整星适配器），不在冻结表；下一轮
应在该适配器内先加一条**永久门禁**直接断言这个条件，而不是只与存储频率比对。

第三轮把这条线索推到底，结论需要修正上面的猜测（勿再引用"就是适配器 bug"）：

- 新增诊断 `examples/{trace_subduction,dump_character_row,dump_sg_rows,trivial_content,
  scan_character_norm}.rs`，全部只读、幂等、可离线复算。
- **折叠约定与存储表一致**：ordinal 58（SG 5 `V1`→#1）折叠出 `GM1` 与 `X1` 两个子群
  k 点，二者相差子群自身倒格矢；存储表也只记 1（记在 `GM1`），不是折叠求和。所以
  "按子群列表 k 点展开"是双方共同约定，审计的读法没错。
- ordinal 12471 上引擎**自洽**：`reconstruction` 逐项误差 2.4e-16，`<χ,χ>` 与
  Σmult² 一致（2 = 2），两边的 Frobenius 总量都等于 index 12；分歧只在**同一个星里
  哪些 irrep 带恒等成分**（引擎 `{X1-,X2-,M1+}` vs 存储 `{X3+,X4+,M4+}`）。
- 该记录的存储行带 `domain` 号 1 与 4（X4+ 是 domain 4），审计把所有 domain 的行
  混在一起比，这部分是**域（共轭子群）混比**；但 domain-1 子集本身也仍不一致，
  所以不能只用域解释。
- `scan_character_norm` 报 3,030 条 `<χ,χ> != Σmult²`；这不是"字符非法"，而是
  展开的子群 k 点粒度造成的（同一子群 irrep 被拆到多个子群列表 k 点上），
  折叠后才成立——已写进该 example 的文档，别把它当 bug 门禁使用。

**第四轮（2026-09-21 深夜）把 hard_failures 从 771 压到 60**，靠的不是改适配器，而是
三处"帧"的问题：

1. `SubgroupEmbedding` 的 `subgroup_lattice` 字段原先一直存 `W·P_parent`，而校验用的是
   候选自身映射推出的格；`U` 非幺模时两者不同。已改为存候选的格。
   （此前只把 `validate_candidate` 的参数换掉了，字段没换——这是个漏改。）
2. **每条记录都必须落在官方打印的 origin 上**。`delta = (B^T)^-1(o_printed − o_stored)`
   此前只为"引擎拒绝过"的记录推导；但两个位置都可能是合法子群，引擎的校验区分不了，
   冻结了非官方的那个就会把恒等成分挪到星里别的 irrep 上。新增
   `derive_shift_recorded.py`（**不需要 oracle**，只用 census 数据）对全部 1,057 条
   origin 不符的记录推导，引擎接受 523 条。最小见证 ordinal 27（SG 3 `A2`→#3）：
   `delta = 0` 时引擎落在 `GM2`、恒等成分 0；加上 recorded-frame 位移后落在 `GM1`、1 ✓。
3. **子群 Hall 行的 ITA origin choice 也要主动推导**，即使 `delta = 0` 能通过校验。
   对审计仍点名的 107 条重新推导两条路线，得 103 条被接受（97 origin-choice、
   6 recorded-frame）。最小见证 ordinal 1197（SG 48 `R1+`→#70）：`delta = 0` 时引擎报
   `GM1-`，官方是 `R1+`。

审计演进：恒等正项通过 93,720 → **94,081**、不匹配 391 → **30**、absent_positive
380 → **30**、hard_failures 771 → **60**；embedding 15,239/15,239、Γ Frobenius
1,895/1,895 不变。

**剩余 60 行 = 16 条记录**（SG 5/12/15/67/68 与 13090/13106），它们的存储 origin 与
打印 origin 相差**母群 cell choice**；`SET I <sg> CELL n` 会改变这些母群的 direction
集合，记录不再一一对应，需要另找途径解析母群 setting。

**第五轮（2026-09-21/22）把全表审计做到 `VERDICT clean`**：embedding 15,239/15,239、
存储恒等正项 **94,111 通过 / 0 不匹配 / 0 假阳性**、Γ Frobenius 1,895/1,895、
`hard_failures=0`，提交为 `fe2fd05`。三处修正：`SubgroupEmbedding.subgroup_lattice`
改存候选自身映射推出的格（此前字段没跟着 `validate_candidate` 一起改）、
`derive_shift_recorded.py` 对全部 1,057 条 origin 不符记录（不只引擎拒绝过的）推导
记录帧位移、以及对审计点名的 107 条重新推导子群 Hall 行的 ITA origin choice。
剩余 60 行 = 16 条记录（SG 5/12/15/67/68 与 13090/13106）由后续轮次解决。

**第六轮（2026-09-22）：恒等内容第二条精确路径，普通恒等分导表范围内闭合。**

审计仍报 14,713 条 probe `uncomputed_missing_data`（另有 160 条属于存储正项），
原因是**折叠出的某个子群星没有随包离散 k 数据**；`trace_subduction` 新增的
折叠星转储（失败时打印每个子群星的全部 q、模子群倒格的余数、是否 Γ，以及子群
表的 k 点）显示：ordinal 10027（SG 196 `W1` `P2` → #24）的探针 `W1` 折叠出 3 个子群
星，其中 2 个（q = (-2,-1/2,0) 与 (0,-1,1/2)）不在 #24 的离散表里，**第 3 个 q =
(1,0,-1) 正是子群 Γ 点**，恒等重数 1，与存储频率 1 一致。

于是新增引擎入口（`src/irrep/subduction_star_decompose.rs`）：

- `trivial_content_with_embedding(subgroup, embedding, probe) -> TrivialContent
  { total, by_label, gamma_stars, skipped_stars }`：只回答**恒等重数**。证明要点：
  子群格平移 `t` 在任何 `q` 表示上作用为 `exp(-2πi q·t)`，恒等表示作用为 1，故共有
  不可约成分要求 `q = 0` 模子群倒格（含 centering 消光）；因此非 Γ 折叠星对恒等
  重数的贡献恒为 0，跳过它们**精确**而非近似。Γ 块仍走完整分解同一套 `build_block`
  与冻结 child origin。`total`（按冻结 CIR 来源号）与 `by_label` 必须一致，否则
  `TargetSourceMismatch`；子群没有唯一恒等 Γ 行时 `MissingChildTrivialIrrep`，绝不
  返回 0（`every_space_group_has_one_trivial_gamma_row` 对 230 个 SG 钉住该行存在）。
- 完整入口 `subduce_full_star_with_embedding` 的契约不变：缺子群数据仍报
  `MissingChildStarData`、绝不部分返回，`TotalDimensionMismatch` 仍要求块覆盖整个
  母群星。

审计把探针结果分成 `full_success` 与 `identity_only` 两类，**两类都与存储表逐条比较**。
全表运行 452 s、`--require-complete` 退出 0：

| 项目 | 第六轮结果 |
|---|---|
| embedding | 15,239 / 15,239 |
| probe 结果 | 351,547 完整 + 14,713 恒等-only = 366,260 全部有精确结果；0 缺数据 / 0 错误 / 0 embedding 不可用 |
| 存储恒等正项 | **94,271 / 94,271 通过**；0 不匹配、0 假阳性（160 条由此前的 missing 变为恒等-only 通过） |
| 未存储项 | 271,989 全为 0（引擎 114,770 + 几何 157,219） |
| Γ Frobenius | 1,895 / 1,895 |
| 生产自检 | 维数/整数重数/重建/CIR 来源不匹配全部 0 |
| `hard_failures` | 0 |

永久回归：`tests/subduction_identity_regressions.rs` 的
`identity_only_content_answers_probes_without_full_child_data`（SG 196 四元组
(1,1,1,2) 与 L1 的无 Γ 星情形）与
`identity_only_content_agrees_with_the_full_decomposition_and_covers_the_pinned_set`
（59 个冻结上下文：2,060 个 probe 上两条路径**逐条相等**，另 15 个 pinned 缺数据
probe 由恒等-only 精确回答）；单元测试另钉住 230 个 SG 的恒等 Γ 行与 SG 221 全探针
一致性。

**第七轮（2026-09-22）：小群特征标解码的进展与未决问题（尚未提交可用的解码器）。**

为了真正算出那 5,756 个频率，本轮开始解 `data_little.txt` 的特征标/矩阵段。已确认的
结构（可直接复算）：

- `little_ops` = 每 (SG, k 槽) 12 个操作 × 4 整数（操作码 + 3 个平移），
  298080 = 6210 × 48；每槽实际操作数在 `little_ops_count`（6210）。
- 矩阵数据是**两级指针**：`little_irr_full_matrices_irr_pointer`（10300，每个 little
  irrep 一项，1-based，指向下一张表）→ `little_irr_full_matrices_pointer`（124000，
  每个 (irrep, 操作) 一项，1-based，指向数据）→ `little_irr_full_matrices`
  （2,220,000 个整数）。SG 2 的实例：irrep idx9 `GM1+` → 操作偏移 [11,12] →
  数据值 [1,2]；idx10 `GM1-` → [13,14] → [1,1]；idx11 `Z1+` → [15,16] → [2,1]。
- **未决**：(a) 值的编码不是简单的 ±1（同一批数据里还出现 3..8，直方图 2:8501、
  1:6472、3:3247、4:1172、5..8 各 152），需找出「值 → 根/矩阵元」的映射；
  (b) 指针与 irrep 的下标对齐存在系统性 off-by-one：按 `ip[i]..ip[i+1]` 取，
  与 shipped `CHARACTERS` 的**多重集**比较只有 SG1 6/8、SG2 5/16 命中；整体后移
  一个 irrep 则 SG1 5/8、SG2 8/16，说明块内 irrep 次序与标签表次序不完全一致；
  (c) 操作次序（`little_ops` 的 12 槽顺序 vs shipped 行的 PIR 次序）未对齐，
  所以现在只能比多重集，不能比逐项。**在这些对齐问题解决前不要用这套解码去
  宣称任何 w 频率**；审计的 `--require-w-complete` 门禁继续生效。

对照来源：`src/irrep/generated_data.rs` 的 `CHARACTERS` + `_char_start`/`_char_count`
（SG 2：`GM1+ = [1,1]`、`GM1- = [1,-1]`、`Z1± = [1,±1]`，与官方表一致，可作独立 oracle）。

**第八轮（2026-09-22）：排除了「次序错位」假设，把问题收敛到值编码本身。**

本轮把上一轮的三条未决项各自查了一遍，结论如下（都可复算）：

1. **块内 irrep 次序与标签表次序一致，不是错位**。逐个核对 (SG, k) 块：SG 2 k1 的
   12 槽数组 `little_irr_label` = `'1+'`,`'1-'` ↔ `little_irr_full_label` =
   `'GM1+'`,`'GM1-'`；SG 1 k3 = `'1'` ↔ `'X1'` ✓。`little_irr_order` /
   `little_irr_old`（10300）不是置换（前者是恒等，后者只是把不在旧表的条目置 0，
   如 SG 1 的 `GP1`），`little_irr_old_map`（4777）只是「旧表 → 新表」的映射。
   因此上一轮「后移一个 irrep 命中率更高」只是巧合，**不能再用位置平移去找对齐**。
2. **`little_irr_full_matrices` 不是字符表**。SG 2 的两个 Γ irrep 解码为
   `GM1+ = [1,2]`、`GM1- = [1,1]`，而 k2（Z）的两个是 `Z1+ = [2,1]`、
   `Z1- = [1,2]` —— **两个不同的 irrep（`Z1+`、`Z1-`）解出同一个多重集 {1,2}**，
   但它们与 shipped `CHARACTERS` 的 {+1,+1} / {+1,-1} 不同 ⇒ 数组里的 1/2
   不是 ±1 字符，而是**编码**（矩阵表索引或根码），必须配合 `little_table` /
   `little_irr_table` / `little_irr_real*` 才能还原。
   SG 1 的旁证：k 点 X1、R1 解出 `[2]`，而 shipped 字符是 `[+1]`；同一个
   aP 块里 Z/U/V/Y/T 解出 `[1]` —— 差异的确切含义仍未定。
3. **`little_table_pointer` 不是偏移表**（6210 项，取值只有 0/1，是「有无表」标志），
   `little_table`（22400）的布局需要单独解码；`little_irr_table_pointer` 是稀疏表
   （10300 中仅 1512 非零，块大小 16..432），两者都不是现成可用的字符来源。

因此下一轮的入口是明确的：先解 `little_table` 与 `little_irr_table` 的布局，
用它把 `little_irr_full_matrices` 的编码值还原成根/矩阵元，再用 shipped
`CHARACTERS`（4,777 条主表 irrep）做逐项 oracle；**在此之前不得用该数组宣称任何
w 频率**，审计的 `--require-w-complete` 门禁继续生效。

**第九轮（2026-09-22）：绕开特征标，先钉死「几何臂」这条路线。**

新增只读诊断 `examples/w_arm_count.rs`（`w_arm_count <ordinal> <label> <vx> <vy>
<vz> [<den>]`）：取该记录的冻结 embedding，把直线方向 v 在母群点群下生成 12 个臂、
逐个折叠 `q = T^T v'` 并用**子群自身**的倒格判定是否落在 Γ，同时打印该记录的存储
频率。三条结论：

1. **帧已确定**：`little_k` 的方向在母群 **conventional** 倒格基里，必须先转成
   母群 **primitive** 分数坐标再折叠。ordinal 13824（SG 225 → #1）conventional
   读法 12/12 臂全落在 Γ，primitive 读法才是非退化的；所有用 conventional 直接
   折叠的计数都必须丢弃。

> **更正（2026-09-23，R6.1 实测；本条与第十二轮一起看，勿再引用上面的"必须转成
> primitive"）**：第十一轮已把 12 臂的根因定位为**约化格用错**（conventional 旋转
> 配 primitive 帧的格），第十二轮修好后两个帧读法都给出 3 个臂。现网实现（R6.1
> `subduce_line_at_parameter`）**全程在母群 conventional 倒格坐标下**：
> `direction`、pinned `little_k`（SG 196 `X1 = (0,1,0)`、`L1 = (1,1,1)/2`）与冻结
> little 群操作都是 conventional，`fold_wave_vector` 直接作用在 `t·direction` 上，
> 没有任何 primitive 换算；唯一一次约化是 `canonical_wave_vector` 把 `k` 约化进
> `Lattice::new(exact_primitive_basis(parent)).reciprocal()`，对 5,756 条 pinned
> `t = 1/4` 波矢零位移。细节见 `docs/subduction-conventions.md` §16 与
> `docs/subduction-r6-plan.md` §4。
2. **锚点：child = #1（P1）时频率 = 源的 full-star 维数**。ordinal 13824 的存储行是
   `6 x DT1..DT4`、`12 x DT5, SM1..SM4`，与 little 表的 dim 6/12 逐项一致；这正是
   「P1 上恒等表示重数 = 维数」的必然结果，可作为公式的基准点。
3. **频率不是裸臂数**。SG 196（point group 23，|T| = 12）的 ⟨110⟩ 线在 primitive
   帧下 4 个臂落在子群 Γ，而存储频率是 **1**（ordinal 10030 → #18，DT1/DT2/SM1
   都是 1）；ordinal 10033 → #16 的 DT1 存储 2、SM1 存储 1。所以折叠后的臂还必须
   按**子群自身的点群**（模子群倒格）并类：4 个臂若同属一个子群轨道就贡献 1
   （10030 ✓），分成两个轨道就贡献 2（10033 ✓）。下一步就是把这个轨道计数和小群
   表示重数 m 一起实现，并用 5,756 行全表核对；在此之前仍不宣称任何频率。

**第十轮（2026-09-22）：子群轨道计数原型，2/4 抽样命中，缺口正是小群重数 m。**

`w_arm_count` 现在还会：把直线方向在母群点群下的像按「同一条线（差 ±）」去重成臂集合，
对落在子群 Γ 的臂再用**子群自身的点群**（`embedding.representatives()` 的旋转）做轨道
并类，打印 `arms | folding to child Gamma | H-orbits of those`。实测：

| 记录 | 子群 | 存储频率 | prim 帧（臂/落 Γ/轨道） |
|---|---|---|---|
| 10030 DT1 | #18 | **1** | 12 / 4 / **1** ✓ |
| 10030 SM1 | #18 | **1** | 12 / 4 / **1** ✓ |
| 10033 DT1 | #16 | **2** | 12 / 4 / **1** ✗（差因子 2） |
| 10032 DT2 | #18 | **2** | 12 / 4 / **1** ✗（差因子 2） |
| 13824 SM1 | #1 | **12** | 12 / 12 / **12** ✓ |

结论：轨道计数在 child = #18 与 child = P1 上正确，但在 #16/#18 的部分记录上系统性地
少一个因子 2 —— 即「每个轨道的小群表示重数 m」（Mackey 双陪集公式里
`mult(trivial_{H ∩ sG_k s⁻¹}, W^s)`）不能一律取 1。

**第十一轮（2026-09-22）：用维数守恒抓住臂集合本身的错误。**

给 `w_arm_count` 加上「按母群倒格 + ± 去重」后（`line_key`），臂数仍是 **12**
（说明那 12 条线确实互不等价 mod 母群倒格），但这条与**维数守恒**矛盾：
`full-star dim = 星大小 × 小群维数`，SG 196 的 `DT1` dim = 6 ⇒ 星大小 ∈ {1,2,3,6}
⇒ 不可能是 12。根因是**帧混用**：母群点群旋转来自 `symmetry_operations_of`（母群
**conventional** 基），而直线方向是 prim 帧的 (1/2,1,1/2)，两者直接相乘是错的。
修法只有两条，二选一：

* 全部在 **primitive 帧**做：旋转先相似变换 `R_prim = P^{-T} R_conv P^T`，臂按
  `Lattice::new(P).reciprocal()` 约化；
* 或全部在 **conventional 帧**做：旋转直接用，臂按母群 conventional 倒格（= Z³ 加
  centering 消光）约化，折叠前再把臂转成 prim。

**第十二轮（2026-09-22）：帧修好，公式结构钉死，只剩小群重数 m。**

按第十一轮的方案修帧：两个读法各用**自己的**约化格 —— conventional 帧用
`Lattice::new(P).reciprocal()`（P = 母群 primitive 基在 conventional 坐标下），
primitive 帧用 `Lattice::integer()`（prim 分数坐标下晶格就是 Z³）。结果两个读法的
臂数都变成 **3**（3 | dim 6 ✓，自检通过；第十一轮的 12 臂是约化格用错导致的）：

| 记录 | 子群 | 存储频率 | prim 帧：臂 / 落 Γ / H 轨道 |
|---|---|---|---|
| 10030 DT1 / SM1 | #18 | 1 | 3 / 1 / 1 |
| 10033 SM1 | #16 | 1 | 3 / 1 / 1 |
| 10031 DT1, 10034 DT1 | #18 | 1 | 3 / 1 / 1 |
| 10033 DT1, 10032 DT2 | #16 / #18 | **2** | 3 / 1 / **1** |
| 13824 SM1 | #1 | **12** | 3 / 3 / 3 |

**结构因此确定**（Mackey 双陪集）：

```text
frequency = Σ_{折叠到子群 Γ 的母群臂的 H-轨道}  m_orbit,
m = mult(trivial_{H ∩ sG_k s⁻¹}, W)
```

* P1 锚点自洽：3 个臂 × m = dim(W) = 4 ⇒ **3 × 4 = 12 = 存储频率** ✓（也用
  `dim = 臂数 × 小群维数` 反过来钉住 little dim：DT1 6/3 = 2、SM1 12/3 = 4 ✓）。
* #18 的 DT1/SM1：1 轨道 × m = 1 ⇒ m = 1 ✓。
* 10033 DT1 / 10032 DT2：1 轨道但存储 2 ⇒ **m = 2**（该直线小群表示在
  `H ∩ G_k` 上含恒等表示两次）。

所以现在**唯一的未知量就是 m**；几何部分（臂集合、落 Γ 判定、H 轨道、帧与约化格）
已经全部对齐，且每次计数都可先用 `dim / 臂数` 是否整除来自检。m 的求法两条路：
(i) 完成 little 表特征标解码后直接算 `mult(trivial, W|_{H∩G_k})`；
(ii) 若这些直线小群表示都是某个 1 维表示的单式诱导（monomial），m 可由陪集结构
直接数出来 —— 这条更省事，值得先试。

**第十三轮（2026-09-22）：`little_subduce_*` 的载荷被解出一半，且给出了新的 m 路线。**

对 SG 196 的直线 irrep 直接读 `little_subduce_*`（`little_subduce_irr_pointer[i]` →
`little_subduce_pointer` → `little_subduce_frequency` / `little_subduce_pg_irrep`）：

| irrep | full_dim | (pg_irrep, frequency) | Σ freq·dim(pg) |
|---|---|---|---|
| DT1 | 6 | (1,1) (2,1) (3,1) | 1+2+3 = **6** ✓ |
| DT2 | 6 | (3,2) | 2·3 = **6** ✓ |
| SM1 | 12 | (1,1) (2,1) (3,3) | 1+2+9 = **12** ✓ |

即载荷是**母群 full-star 表示按母群点群 irrep 的分解**（点群 23 的 pg 序号 1/2/3 =
A(1 维)/E(2 维)/T(3 维)），并且满足
`Σ frequency × dim(pg_irrep) == little_irr_full_dim` 这条**可检查的不变量** ✓。
这条路线比解码 2.22M 个整数矩阵编码更省：有了 V（full-star 表示）按点群 irrep 的
分解，再用点群特征标（标准、可自造）与嵌入里 H 的点群，就能算
`mult(trivial_H, V|_H)`，其中 W 的贡献可由诱导关系从分解倒推。

**第十四轮（2026-09-22）：中间层对齐，块结构门禁通过。**

按上一轮的入口逐行走了 `little_subduce`：`little_subduce_irr_pointer[i]` 给出 irrep i
在**行表**（`little_subduce_pointer` / `little_subduce_count`，各 47000 项，索引口径
是「每**行**一项」而不是每 irrep 一项）中的起始行，一个 irrep 占
**blocksize(SG) 行**（SG 196 = 8 行；DT1/DT2/SM1 的连续指针差都是 8 ✓），每行再
指向 105000 项的 (`little_subduce_pg_irrep`, `little_subduce_frequency`) 载荷。

**修正第三轮以来的一处读法**：块的最后一行不是「`[(dim,1)]`」，而是
**`[(1, full_dim)]`** —— 即 pg 序号 1（各子群上下文里的恒等/一维表示）配重数
= `little_irr_full_dim`。证据：SG 196 `DT1` 的最后一行 `[(1,6)]`、`DT2` `[(1,6)]`、
`SM1` `[(1,12)]` ✓。

**门禁**：对每个能解析出块的 little irrep，检查「最后一行 == [(1, full_dim)]」——
**3929 条通过、5 条不匹配**（那 5 条是 SG 1 前几个 irrep，指针是 0 哨兵、
行表读成空 ✓），0 条缺失 ✓。即中间层的索引口径已经对齐，可以放心地按
「irrep → blocksize 行 → 每行 (pg_irrep, frequency) 分解」来读。

每行的点群上下文不同（首行/末行用母群点群维数时 Σ freq·dim 才等于 full_dim：
DT1 首行 A+E+T = 6 ✓、DT2 首行 2T = 6 ✓），所以倒推 m 时**必须先确定行与子群
上下文的对应**，这仍是下一步。

**第十五轮（2026-09-22）：`little_subduce` 的行语义定型 —— 它很可能就是 w 频率的来源。**

把第十四轮的行数据与第十二轮的几何结果对起来看，得到一条自洽的模型：

> 每行 = 该 little irrep（full-star 表示 V）在**某个子群点群上下文**下的分解
> `(pg_irrep, frequency)`，其中 `pg_irrep = 1` 是该上下文里的一维/恒等表示；
> 于是 **w 频率 = Σ over「折叠到子群 Γ 的臂的 H-轨道」of（该臂上下文那一行的
> pg_irrep=1 的频率）**。

证据（都可复算）：

* SG 196 `DT1` 的 8 行里，pg=1 的频率分别是 1,1,1,1,2,2,4,6；而 SG 196 `DT1` 的
  32 条存储 w 频率取值集合是 **{1,2,3,4,6}** ✓ —— 4 和 6 就是某两行本身的值，
  3 = 1 + 2（两个轨道落在不同上下文的两行上）✓，这与「按轨道求和」完全吻合；
* 第十二轮已经证明几何侧只有 1 个 H-轨道（10030 → 1 ✓、10033/10032 → 2 ✓），
  即模型只需从**有限个行**里按上下文挑一个值（或跨轨道求和），不存在连续自由度；
* 行的上下文可以由**已知的 full-group 分解**反过来标定：DT1 的首行给出
  V = A + E + T（母群点群 23 的标准 irrep），把 A+E+T 限制到每个子群点群类型
  （C1/C2/C3/D2/T …）用**标准点群特征标表**就能算出该子群应有的分解，
  再与 8 行逐一比对即可把「行 ↔ 子群点群类型」钉死。

**第十六轮（2026-09-22）：把 8 行逐行摊开，行 ↔ 上下文只差最后一步。**

SG 196 的 `GM1`（指针 0 = 哨兵，块读出来是垃圾，**不能**用它做例子）/`DT1`/`DT2`/`SM1`
八行全表（`Σf·dim` 用母群点群 23 的维数 A=1,E=2,T=3）：

| 行 | DT1 | DT2 | SM1 | Σf·dim vs full_dim |
|---|---|---|---|---|
| 1 | (1,1)(2,1)(3,1) | (3,2) | (1,1)(2,1)(3,3) | 6/6, 6/6, 12/12 ✓ |
| 2 | 同 1 | 同 1 | 同 1 | ✓ |
| 3 | 同 1 | 同 1 | 同 1 | ✓ |
| 4 | 同 1 | 同 1 | 同 1 | ✓ |
| 5 | (1,2)(2,2) | (1,2)(2,2) | (1,4)(2,4) | 6/6, 6/6, 12/12 ✓ |
| 6 | (1,4)(2,2) | (1,2)(2,4) | (1,6)(2,6) | 8/6, 10/6, 18/12 ✗ |
| 7 | 同 6 | 同 6 | 同 6 | ✗ |
| 8 | (1,6) | (1,6) | (1,12) | 6/6, 6/6, 12/12 ✓（末行标记） |

三条可复算的结论：

1. **行 = 上下文**：1–4 行在三个 irrep 上逐字相同（同一上下文重复 4 次），5 行彼此
   同型，6–7 行同型，8 行是标记 —— 即 blocksize 8 = 7 个上下文 + 1 个标记。
2. **每行用自己上下文的点群维数**：1–5、8 行用母群维数 (1,2,3) 时
   `Σ freq·dim = full_dim` **全部成立** ✓；6、7 行不成立（8/6、10/6、18/12），
   但若该上下文的一维表示维数是 (1,1)（C2/Ci/D2 类）则 6/6、12/12 ✓ 成立 ——
   所以「行 ↔ 上下文点群」的判定标准就是**这行在哪个维数表下配平**，这是一条
   纯算术门禁，不需要特征标表就能先筛出候选。
3. **块内对齐要按 irrep 的 pointer 走**：SG 196 `GM1` 的 `little_subduce_irr_pointer`
   是 0 哨兵，按 `ptr-1` 取会读到错位数据（上面第一版读出的 GM1 行即为反例），
   任何按 irrep 取块的代码都必须显式跳过 0。

**第十七轮（2026-09-22）：w 行的语义、k 域、官方 oracle 全部落地；只剩 little 群特征标。**

本轮不再猜行语义，改用官方程序本身把这条轨道钉死（细节见
`docs/isotropy-data-semantics.md` §4 与 `docs/subduction-audit.md`）：

1. **语义确认**：对 `(母群 SG, irrep)` 运行 `iso` 的 `DISPLAY ISOTROPY` +
   `SHOW FREQUENCY`，每个方向一行，频率表列出**所有（任意 k 的）**包含该子群恒等
   表示的 irrep 及重数 —— 它正是 pinned `isotropy_subduce_*` ∪
   `isotropy_w_subduce_*` 的并集（SG 196 W1 的 C5/C7/C11 行逐条对上）。命令细节：
   `SC 1000`（默认宽度会把长频率表折行，续行可能长得像一行数据，判定必须靠缩进）、
   `SET I ALL OR 1`。
2. **k 域帧修正（推翻第九轮结论）**：`little_k` 的方向在母群 **primitive 倒格基**
   里，官方 `DISPLAY KPOINT` 打印 **conventional** 帧；SG 196 pinned `DT=(1,0,1)`
   ↔ 官方 `DT (0,2a,0)`（F 心倒格基矩阵），`SM=(1,1,2)` ↔ `(2a,2a,0)`。73/73 源
   都是参数化 k 域（线/面/一般位置），**k 参数不缺**。
3. **little 表结构对齐**：`little_k_dim` = 自由参数个数、`little_ops_count` = little
   群阶、`little_irr_dim`（6210×12 缓冲）= 各 little irrep 维数、`little_irr_full_dim`
   = 物理 irrep 维数 = 星大小 × little 维数（compound 条目合并 k 等价 little irrep，
   如 SG 196 `LD1LE1` = 2×4 = 8；`little_irr_old_map[i]` 给出母群 irrep 的紧凑序号，
   112 条只是旧 ML 拼写）。`little_subduce_*` 只在**有自由参数的 k 域**上有块（所有
   k 点的 irrep 指针都是 0），最后一行恒为 `[(1, full_dim)]`（并非早期写的
   `[(dim,1)]`）。**载荷给不出 w 频率**：SG 196 `DT1` 的 8 行 pg=1 值是
   `1,1,1,1,2,2,4,6`，而官方对同一线算出的子群频率集合是 `{1,2,3,4,6}`（`3` 不在
   行值里），且行没有 SG/basis/origin/direction 键 —— 第十五、十六轮的
   「行值 = w 频率」模型作废。
4. **新门禁**：`scripts/verify_w_subduction_oracle.py`（离线回归
   `scripts/test_verify_w_subduction_oracle.py` 14 项）按 `(子群号, Dir 标签)`
   **双向**比较官方频率表与 pinned 两张表。实测 **28 组 / 1,150 条记录 /
   oracle 5,756 条 w 行 = pinned 5,756 条 w 行 / 0 不匹配**，其中 144 条是
   「oracle 无 w 条目」的负向核对。
5. **旁证**：对参数化 k 域请求 `DISPLAY ISOTROPY` 会让程序现场算该 little irrep 的
   isotropy 表并缓存为 `i<紧凑序号><参数>.iso`；文件里每个子群的方向矩阵**行数 =
   恒等分导重数**（`4D1` 4 行、`6D1` 6 行，与 pinned `#1` 的 4/6 一致），而
   Frequency 列本身对参数化 k 留空。另解出可用的元素写法：
   `LABEL ELEMENT INTERNATIONAL` + `X Y Z` / `-X -Y -Z`（配 `SHOW CHARACTER`）。
6. **状态**：引擎仍不能算这 5,756 行（缺 little 群特征标；`little_irr_full_matrices`
   2.22M 整数编码未解），`--require-w-complete` 继续退出 2；发布必须按
   `docs/subduction-audit.md` 的两条轨道声明覆盖。下一轮：从官方 oracle 提取
   73 个 little rep 的特征标（或在 pinned little 表里解出矩阵编码），把 w 频率接入
   引擎，让 `--require-w-complete` 也能全表退出 0。

**第十八轮（2026-09-22）：73 个 little rep 的特征标用官方数据解出（65/73 唯一确定）。**

不解 `little_irr_full_matrices`，改用兼容关系：对母群取 Γ 各 irrep 的
`SHOW COMPATIBILITY`（限定到该线，如 SG 196 `GM4: DT1 DT2 DT2`）与这些 Γ irrep
在每个 little 群操作上的 `SHOW CHARACTER`（写法 `LABEL ELEMENT INTERNATIONAL` +
`VALUE ELEMENT X Y Z`）。k = Γ 上线的 little irrep 无 Bloch 相位，于是每个操作 R
给出精确有理方程 `Σ_i n_Γ·mult(Γ→i)·D_i(R) = χ_Γ(R)`（`n_Γ` = compound 行 ML
分量数）。线的各源**联合求解**，每个源是否可定用**对偶系统**判定（源特征标是分量
线性泛函，落在行空间即有唯一值 —— compound `DT3DT4` 因此可用）。

新增 `scripts/freeze_w_little_characters.py` + `scripts/test_freeze_w_little_characters.py`
（13 项离线回归；解析 `SHOW ELEMENTS` 的 ITA 串、k 矢量参数化、兼容表、共享行的
element/character 列、精确消元与对偶判定）。实测 **73 源中 65 个唯一确定、0 残差**，
`--json` 落在 `target/task9/w_little_characters.json`（每源含方向、little 群操作
的旋转/平移/特征标与所用方程）。剩 8 个是 SG 202/203/209/210 的 `DT3`/`DT4`：
它们的 Γ 兼容表把两者以相同重数混合，只能定和；**已验证**补同线 X 点（α = 1/2）
可分开（SG 202：`X3+→DT3`、`X2-→DT3`、`X2+→DT4`、`X1+→DT1`、`X4+/X1-→DT2`），
下一轮把线上其它特殊点方程并入（点的特征标乘 Bloch 相位 `exp(-2πi α v·t)`），
再做引擎侧 Mackey/特征标求和，让 `--require-w-complete` 能全表退出 0。

**第十九轮（2026-09-22）：线上补点被「星污染」挡住，脚本加硬门禁拒绝假数据。**

把 DT 线上的 X 点（α = 1/2）方程（兼容表 + 特征标 + Bloch 相位
`exp(-2πi α v·t)`，复数消元 + 对偶判定）实现进 `freeze_w_little_characters.py`
后发现根因：**星大小 > 1 的 k 点上 `SHOW CHARACTER` 打印的是完整 irrep 的特征标**
（SG 202 `X3+` 维数 3），而兼容行是 **little irrep** 层面的关系 —— 两者只在星大小
为 1（即 Γ）时相等。直接把 X 方程代入会解出 `D(E) = 3` 之类的非法值。脚本因此加
门禁「恒等特征标 = `full_dim / 星大小`」，把这类解判为失败：结果仍是
**65/73 唯一确定、8 个如实报未定、0 假数据**（重复报告也已去重）。剩余两条可行
路线：`data_images.txt` 的 image（点群表示）数据库，或 pinned little 矩阵段的解码。

**第二十轮（2026-09-22）：冻结 65 个源的特征标进 crate，并量化缺口。**

`freeze_w_little_characters.py` 新增 `--rust`，生成
`src/irrep/w_little_characters_data.rs`（65 个表：父群、标签、k 域、方向、little 维数、
每个操作的 ITA 串 + 旋转 + 平移 + 特征标；另附 8 个未定源的清单），
`src/irrep/mod.rs` 以 `#[doc(hidden)] pub mod` 接入；新增集成回归
`tests/w_little_characters.rs`（4 项）：65/8 的规模与互斥、每表首操作必须是恒等且
`χ(E) = dim`、`|χ| ≤ dim`、SG 196 `DT1 = [1,1]`、`DT2 = [1,-1]`、`SM1 = [1]` 与
`DT1 ≠ DT2`、以及各立方母群的标签集合（含 202/203/209/210 只剩
`DT1/DT2/SM1/SM2`）。

缺口量化：5,756 条 w 行中 **348 条（194 条记录，6.0%）**引用那 8 个未定源，其余
**5,408 条（93.99%）**的源已有冻结特征标，等引擎侧 Mackey/特征标求和接入。

验收：release lib `389 passed / 4 ignored`、全部 integration（含新 4 项）、
doctest `27 passed`、严格 all-target clippy 零警告。

**第二十一至二十五轮（2026-09-22）：Mackey 原型、公式定型、数据侧凝固。**

- 第二十一轮：Python 原型（`target/task9/explore/proto_freq.py`，未入库）跑通公式骨架，
  `6D1`（P1 子群）三源全中（`DT1=6, DT2=6, SM1=12 = dim V'`），超胞情形偏高。
- 第二十二轮：修正 little 群成员判据（只约束旋转 `Rv = v`）并写完整有限群投影：
  `mult = 1/(|P_H|·n) Σ_R Σ_{T∈L_H/L_parent} χ_{V'}((R,t_R)(E,T))`；`6D1` 锚点自洽
  （`n=6`、`Σ_T χ = 72` ⇒ `72/6 = 12`）。
- 第二十三轮：证明「对 `T` 求平均」就是引擎已有的**折叠到子群 Γ 判定**，
  因此 Rust 侧最小改动 = 复用 `trivial_content_with_embedding` 的臂枚举/折叠/块组装，
  只换臂字符来源（冻结 little 表 + Bloch 相位）。
- 第二十四轮：反例约束——用官方 `SHOW BASIS` 的子群基矢自写折叠判定会给出 12 而非
  pinned 4（臂集合与帧/中心化约定），**禁止另写折叠判定**。
- 第二十五轮：新增 `scripts/test_frozen_w_little_characters.py`（3 项，无需 cargo/iso），
  从 Python 侧钉住生成文件 `src/irrep/w_little_characters_data.rs`（65/8 形状、
  `χ(E)=dim`、`|χ|≤dim`、SG 196 `DT1=[1,1]`、`DT2=[1,-1]`、`SM1=[1]`），
  与 `tests/w_little_characters.rs` 互为镜像；五个 Python 离线套件全绿。

当轮复跑的三条 live gate（可复算）：w 行 oracle `28 组 / 1,150 记录 /
5,756 oracle w 行 = 5,756 pinned / 0 不匹配`；几何 oracle `62 行 / 26 描述串 /
62 origin 全精确`；w 结构门禁 `checks_failed=0` 且打印
`frozen_characters: rows_with_frozen_table=5408 rows_blocked=348`。

仍缺：引擎侧把 w 行算出来（5,408 行有冻结表；348 行卡在 SG 202/203/209/210 的
`DT3`/`DT4`），之后 `--require-w-complete` 才能全表退出 0。

**第二十六轮（2026-09-22）：全表审计复跑，判词不变。**

`cargo run --release -p cryspglib --example audit_irrep_subduction --
--require-complete --output <abs path>`（446.7 s，exit 0）复现：

```
identity_rows=94271 passed=94271 mismatch=0 unresolved=0
probes_total=366260 probe_partition: full_success=351547 identity_only=14713 missing=0 error=0
absent_zero=271989 absent_positive=0 frobenius: records=1895 passed=1895 mismatch=0
production_checks: dimension/integrality/reconstruction/target_source 全部 0
other_wave_vector: records=1006 rows=5756 source_resolved=5756 computed=0
w_scope: rows=5756 computed=0 uncomputed=5756 reason=k_vectors_and_character_rows_absent_from_the_irrep_table
hard_failures=0 accounting_violations=0 census_mismatch=0
VERDICT complete scope=global
```

注意 `--output` 的相对路径按 **crate 目录**解析（`cryspglib/`），从 workspace 根运行时
要么给绝对路径，要么先把 `target/task9` 建在 crate 下。

**第三十/三十一轮（2026-09-22）：审计的 `w_scope` 行带上冻结/受阻拆分。**

审计对每条 w 行解析源之后，再查该源的 little 群特征标是否在
`w_little_characters_data` 里（新计数器 `w_character_frozen` / `w_character_blocked`），
并打印进 `w_scope`。全表复跑（`--require-complete`，exit 0）实测：

```
other_wave_vector: records=1006 rows=5756 source_resolved=5756 source_mismatch=0 computed=0
w_scope: rows=5756 computed=0 uncomputed=5756 character_tables_frozen=5408
         character_tables_blocked=348 reason=... gate=--require-w-complete
hard_failures=0 accounting_violations=0 census_mismatch=0
VERDICT complete scope=global
```

与 Python 门禁 `frozen_characters: rows_with_frozen_table=5408 rows_blocked=348` 完全一致
（局部抽查：`--parent 196` → 106/0、`--parent 202` → 216/100）。19 个测试二进制全过、
严格 clippy 零警告。

**范围之外的剩余问题**：`isotropy_w_subduce_*` 的 5,756 行（1,006 条记录）引用 73 个
“别的波矢”irrep；pinned `data_irreps.txt` 对它们只有
`irrep_w_label/_space_group/_dimension/_type` 四张表，**既无 k 矢量也无特征标行**，
所以从该表无法计算。**归档的 `data_little.txt` 里 73/73 个源都在**
（`little_irr_full_label` + `little_irr_space_group` + `little_irr_full_dim` 与
`irrep_w_dimension` 逐项相符，例：SG 225 的 `DT1/2/3/4` dim 6、`DT5` dim 12、
`SM1–4` dim 12），而且本轮把 `little_k` 解码出来了：14 个 Bravais 格 × 27 个 k 槽
（6210 = 230×27）× 16 整数 = 4 组 `(x,y,z,d)`，即**基点 + 至多三个自由方向**
（aP 块与 `data_space.txt` 的 k 点逐值相符：Z=(0,0,1)/2 … T=(0,1,1)/2；GP 三个
方向 (1,0,0),(0,1,0),(0,0,1)）。据此 73/73 个 w 源**全部是参数化直线**
`k = Γ + t·v`、各一个自由参数（cF 的 `DT=(1,0,1)`、`SM=(1,1,2)`；cI 的
`DT=(1,-1,1)`、`SM=(0,0,1)`），所以不存在单一数值 k 可以折叠——这正是任务卡
第 5 条「参数化波矢」的情形。剩下的工作是解码这些直线上的小群特征标，
再按线对照 5,756 个存储频率。新增
`scripts/check_other_wave_vector_rows.py`（+ 15 个
离线测试）把这件事变成可执行检查：w 段必须**恰好**是那四张表（多出 k/矩阵段即失败）、
数组长度与 1006/5756 相符、稀疏 pointer 的非零值集合等于带 w 记录的 1-based 起始
偏移集合、源 SG 等于记录的母群 SG、频率 ∈ [1, dim]，以及**73/73 个 w 源都在
参数化直线上**（`w_wave_vectors: sources=73 resolved=73 free_parameters=[1]
base_points=['0'] directions=['(1,0,1)/1', '(1,1,2)/1']`；若某个源变成固定 k，
门禁立刻失败并提示「它已经可计算」）。审计也把每行的解析结果
（`parent_sg_match` / `frozen_source` / `k_parameters=absent_from_the_irrep_table`）
写进 TSV，并新增 `--require-w-complete` 作为这条独立问题的门禁（全表运行时退出 2，
摘要打印 `w_scope: rows=5756 ... uncomputed=5756`）；`--require-complete` 只门禁
普通恒等分导表，其判词为 `VERDICT complete scope=global`。

### 分导任务 9 第四十一轮（2026-09-22）：w 源 little 特征标 73/73 冻结完成

`src/irrep/w_little_characters_data.rs` 现在覆盖**全部 73 个参数化 k 源**
（`W_LITTLE_CHARACTERS_UNRESOLVED` 为空），因此 5,756 条 w 行的频率数据齐备，
只剩引擎侧求值。两条来源：

- **65 个源**：Γ 兼容行联合求解（精确有理对偶系统，已有）；
- **8 个源**（SG 202/203/209/210 的 `DT3`/`DT4`，348 行）：Γ 只能定 `DT3+DT4`
  （对偶系统给出 `(2,-2,0,0)`），拆分用**小群配对路线**（`cogroup_pair_route`）：
  pinned 源次序即小群标准表。SG 202/203 的小群是 `C2v`，`DT3=(1,-1,1,-1)`、
  `DT4=(1,-1,-1,1)`；SG 209/210 的小群是 `C4`，未定对是共轭对
  `DT3=(1,-1,i,-i)`、`DT4=(1,-1,-i,i)`。路线有三重证据：对偶和、每母群 Γ 自定的
  `A1`/`A2` 次序校验、以及（SG 203/209/210）线上 X 点（α=1/2）兼容表 + Bloch
  相位给出的**独立同值**读数。`tests/w_little_characters.rs` 另加
  `archived_cir_characters_confirm_the_sg202_dt_pairing`：用归档 CIR 字符在 X 点
  扣除 `exp(-2πi k_X·t)` 后逐操作复现 SG 202 的四张表（8 条兼容关系 × 4 旋转）。
- 因为 `C4` 对是共轭对，`LittleOperation::character` 由 `i32` 改为精确
  `[实部, 虚部]` 整数对（其余表虚部为 0）；生成器 `--rust` 仍逐字节可复现
  （新 md5 `32cc5ac46d6994a384b22f65a318d3f3`，exit 0、73/73 无失败）。

覆盖门禁随之更新：`scripts/check_other_wave_vector_rows.py` 的
`UNRESOLVED_SOURCES` 置空，实测打印
`frozen_characters: rows_with_frozen_table=5756 rows_blocked=0 unresolved_sources=[]`；
Python 离线套件 9+16+14+17+4 全绿，Rust `w_little_characters` 6 项、
lib `389 passed / 4 ignored`、doctest 27、严格 all-target clippy 零警告。

**可证伪预测（下一轮引擎求值时必须复现）**：SG 202 有 6 条 pinned 记录的
`DT3`/`DT4` 频率不同、SG 209 有 4 条（SG 203/210 全同），配对次序若错，这 10 条
记录会立刻暴露。

### 分导任务 9 第一轮全表冻结（2026-09-21 晚，覆盖已满但一致性未过）

本轮把冻结表从 75 条扩到**全部 15,239 条记录**，两条被冻结的约定都来自官方程序、
并由引擎自己的校验决定去留（`examples/probe_subduction_settings.rs` →
`SubgroupEmbedding::probe_embedding`，绝不靠重新实现校验来宣称通过）：

- `U = W·P_parent·(P_sub·B_oracle)^-1`，以精确有理数 `分子/分母` 存储。少量单斜记录
  只有通过**分数**换基才能到达官方子群胞（`Mat3I` 装不下），因此
  `FrozenEmbeddingSetting` 增加分母字段，`probe_embedding` 同步加参数。
- `child_shift = -(B^T)^-1 (o_OR1 - o_OR2)`：官方程序两种 ITA origin choice 下**同一
  记录**的打印 origin 之差（母群帧）映射进子群帧。抽样 24 条 origin-choice-2 母群
  记录上，只有该读法让全部 24 条通过引擎校验（其余三种转置/符号读法分别只覆盖
  12/11/7 条）；用它替换"先搜到就收"后，Γ 恒等项不匹配从 1,159 降到 434、
  absent_positive 从 998 降到 422、embedding 失败从 267 降到 95。**先验推导优先于
  搜索**是这一轮最重要的方法论结论。
- 195 条官方空表不是缺数据，而是**旧 ML 拼写**：`VALUE IRREP W1W1` 打印空表，
  紧凑 `little_irr_full_label` 拼写 `W1WA1` 才有表。映射逐记录证明（要求打印表与
  存储行在行数/方向标签/子群号/size 上完全一致），112 个母群 irrep、195 条全部
  解析成功，其中 72 条同时 origin 精确。

生产审计（全表）：embedding 15,144 成功 / 95 被拒；标量 probe 348,932 完整、
14,713 缺子群 k 数据、760 engine error；存储恒等正项 93,350 通过 / 434 不匹配 /
245 不可用；Γ Frobenius **1,895/1,895**；absent_positive 422；5,756 条其它波矢仍缺
k 参数。`VERDICT inconsistent hard_failures=1698` —— 覆盖已接近满，**一致性尚未通过**。

Γ Frobenius 覆盖从 14 条扩到 243 条后暴露出**测试自身的 compound 记账错误**：
compound 行是两个复不可约成分之和，`行维数 × 行重数` 是两个和的乘积、会把交叉项
重复计数（ordinal 2978 因此报 6 ≠ 4）。按文档既有的 `k_G` 折算（除以成分数）后
16 个假失败全部消失；这不是放宽门禁，而是修掉一个一直存在的重复计数。

**剩余瓶颈（下一轮入口）**：95 条 `FrozenEmbeddingRejected` 集中在单斜 SG 3–15
（其 isotropy 表记录在 `SET I <sg> AXIS c` 一类的 cell/axis choice 下，而程序的
默认是 unique axis b），因此**必须先逐母群解析 recorded setting**，再做同样的推导。
`scripts/task9/census_settings.py` 是这条扫描的草稿（尚未接入主流线）；
`scripts/task9/README.md` 记录完整管线、复现命令与当前数字。

### 分导任务 9：全表审计与逐记录 setting（2026-09-21，覆盖尚未闭合）

- DSH 实现采集器、生产 API 审计和元数据生成器，另一个 DSH 只读对抗复核；
  Codex 独立跑全表、复现问题、修正门禁、核对坐标与标签、最终集成。
- `examples/audit_irrep_subduction.rs` 遍历全部 15,239 记录和 366,260 个标量 probe，
  正项 94,271、未列项 271,989，逐行 TSV 无重复或漏行。每个可嵌入 probe 都调用
  完整星 API；几何零项不能跳过计算。ordinal 13345 的 W1–W5 缺数据必须使
  `--require-complete` 返回 2；黄金 ordinal 12400 必须实际分解全部 40 个源表示。
- 当前：75 embedding 成功、12,884 歧义、2,280 无有效 setting；2,453 完整分解、
  19 缺数据、363,788 因 embedding 未计算。494 正项、1,959 零项、19/1,895 条
  Γ Frobenius 已验证，错误 0；其余 93,777 正项未计算，5,756 个 w 条目缺 k 参数。
- `generate_subduction_settings.py` 从官方基矢精确反算 U，并核验完整子群操作集与
  child origin 来源；按 ordinal 冻结到 `subduction_settings_data.rs`。原 69 条结果
  保持，新增 427、984、1942、2102、15125、15131 六条（含三条 shear）。
  #43 交换轴保持群操作集却交换 GM3/GM4；源字符 ±1 和端到端结果共同固定标签。
  新六条的 120 个 probe 中 116 完整；缺失集合为 15125/15131 的 P1P2、P3。
- `audit_subduction_settings.py` 完成 4,777 官方查询：13,861 候选、238 基矢不符、
  945 其余 origin 不符、195 官方空表。候选含 95 条非 signed-permutation U；
  **候选并非已支持 embedding**。完整来源、分母、重跑命令见 `docs/subduction-audit.md`。
- 复核修正已入永久测试：采集 ordinal 从错误的 1-based 改为与 Rust 一致的
  0-based；程序及全部运行数据从 pinned ZIP 私有提取，禁止使用未校验的本地解压；
  `SG_DATA_HALL` 必须吻合冻结 ISO--IR 来源；Γ 异常不能只计数而不影响退出码。
- 任务 9 尚未达到普通表覆盖闭合；下一步是扩大有来源的逐记录元数据，并处理
  各类 setting、官方空表和离散子群 k 数据缺口，不能直接转入“全表已通过”的声明。

### 分导任务 7 状态（2026-09-21）

- 任务 7（k 折叠与非 Γ 相位）已落地：精确 `k_H = T^T k_G` + 子群倒格等价类匹配
  （保留 centring 消光）、单臂 star 限制（`UnsupportedMultiArmStar`）、缺数据报
  `MissingIrrepData`、非 Γ 端到端用例（16 `R1`→#22 的 X/Y/Z、221 `GM4+` P1→#83 的
  R1+、167 `GM3+` P1→#15 的 T3）。
- 关键语义（新钉死）：shipped 角色行**自带代表元的 Bloch 相位**
  `χ(t+L) = χ(t)·exp(+2πi k·L)`——见证是 SG 139 `P1`（k=(1/2,1/2,1/2)）恒等操作
  在 t=0 为 1、在 I 心平移 (1/2,1/2,1/2) 为 −i；配对必须「旋转 + 平移模格 + 相位
  修正」，且多代表元修正后必须一致（符号错会不一致）。已写成永久测试
  `shipped_rows_carry_the_bloch_phase_of_their_representatives`。
- 子群侧帧歧义（shipped 行在子群 Hall setting，嵌入映射记录 setting，可能差冻结
  原点平移，如 #126 的 (1/4,1/4,1/4)）现在**确定性地**撤销冻结 `child_shift`
  （`shift_operations(active, -δ)`），不再"在零与 `-δ` 两种读法里取第一个通过
  完整检查的"；拉回时也不先对子群格约化，保留的平移由 `character_of` 的
  「旋转＋平移模 `L_H`＋Bloch 相位修正」配对解释。
- 代表元格缺陷修复（复核 P1）：`representatives()` 必须由**映射后的原始操作**对
  `L_H` 去重，而不是先把操作对 `L_G` 约化再去重。SG 139 `M1-` P1 → #126 的
  反演 `(-I | 5/2,5/2,5/2)` 对 `L_G`（I 心）为零平移、对 `L_H` 为
  `(1/2,1/2,1/2)`。旧实现仍返回 16 个代表元，但其中 8 个平移类错误，
  随后分导报 `OperationNotInCharacterRow`。
  `operations()` 仍保存对 `L_G` 约化的母群成员视图。
- 抽样覆盖（任务 7 测试里实际检查的四个 (母群,子群) 对：16 `R1` P1、
  221 `GM4+` P1、221 `GM4+` P2、139 `M1-` P1）的全部标量探针（含 Γ）：
  折叠成功 92、多臂 57、
  缺数据 0、其它错误 0；其中 SG 139 该对贡献 folded 20（全部 Γ/M 单臂探针）
  与 17 个真多臂。旧记录的"160/164"不出自任何当前被 pinned 的检查（应是更早的
  一次性扫描），不要引用；Γ 侧当前清点见上文任务 8c 扩充后的结果，另有独立的
  stored-frequency gate `subduction_identity_regressions`（59 条记录 / 538 次比较）
  仍全绿。
- SG 139 的非 Γ 期望值已钉住：`M1-` → `GM1+`、`M1+` → `GM1-`、`M2+` → `GM2-`、
  `M2-` → `GM2+`，而 Γ 的 `GM1+` → `GM1+`；反演代表元按 `L_H` 为
  `(1/2,1/2,1/2)`（`tests/subduction_regressions.rs`）。
- 回拉平移的相位回归：16 `R2` P1 → #22，探针 `X1` → `T3`；过早对子群格
  归约会错误返回 `T2`，且内部重建仍然通过。官方输出已核对该嵌入的
  `B = 2I`、`origin = (1/2,1/2,0)`，永久测试位于 `subduction_identity_regressions`。

### 顺带清理

- 删除 `src/irrep/settings_data.rs` 与 `scripts/extract_sg_settings.py`：该表
  用 stride 3 误读 4 整数的 origin（205/230 个 SG 因此拿到非零垃圾整数，例如
  SG14 `[0,0,3]`、SG230 `[-3,4,-7]`），且其来源（每个 SG 首个 irrep 的首条
  isotropy 记录）按正确 stride 解码后恒为 identity/零，属既错又空的表。
  `IrrepRecord::sg_setting` 与 corep.rs 中三处 debug 输出的 “bilbao 修正” 一并
  删除：那三处都在 `#[cfg(test)] mod tests` 的 `diagnose_*` 测试里且不含任何
  断言，因此测试结果与生产行为不变；只有诊断打印在 1e-16 量级上与旧输出不同
  （旧修正用的垃圾整数在 mod 1 下等价于格矢量，但 `((t-k)%1+1)%1` 与
  `(t%1+1)%1` 不是逐位相同）。

---

## 磁 symmetry 全覆盖实时账本（2026-07-31 起）

用户当前优先目标：先覆盖并修正全部 1–1651 UNI 的磁对称性，再继续上层
corep/summary 产品化。每次全扫、根因确认和修复后都必须更新本节；不能只在
对话中报告。

### 当前最终状态（更新至 2026-08-30）

磁 symmetry/corep 主线、Rust-native API 收口、有限域类型化、全仓库 rustfmt、
角色科学计数法修复、exact scalar Hall phase materialization，以及严格复数
magnetic summary 均已实现。相关历史基线依次为 `ba2d997`、`969a89c`、`daf04eb`。
后续开发和回归应以下面这组最新结论为基线，而不是再沿用早期“部分磁群不支持”
或“只预览 6 个特征标”的判断。

| 层级 | 当前结果 |
|---|---|
| 磁数据库与群代数 | `UNI 1..=1651` 全覆盖；数据库内 `4479/4479` 个 UNI–Hall setting 均通过严格群闭包、逆元、恒等元、陪集及类型一致性检查。 |
| setting 识别与消歧 | 给定 parent/family Hall 时 `4479/4479` 精确回环；仅凭操作自动识别时有 `4461` 个唯一结果和 `18` 个显式 `MagneticUniAmbiguous`，歧义只涉及 `{275,282}`、`{277,284}` 两组，`UNI 283` 唯一。程序不得猜选候选，必须用 parent Hall/setting 上下文消歧。 |
| ISOTROPY/Hall 嵌入 | `1651/1651` 个 unitary subgroup、detected Hall 与 data Hall setting 已完成一致嵌入。 |
| 高对称点与 corep | 用户可按 UNI、BNS 或磁操作输入磁群，取得高对称点列表；选定高对称点后可得到 fixed-`k` magnetic little-group coreps、维数、Wigner 类型和复特征标。全库审计共覆盖 `10,390` 个高对称点、`52,793` 个按正式来源去重的 coreps。 |
| 正式特征标表 | 已提供“每个磁操作一列”和“每个共轭类一列”两种 Markdown 正式表格；不再截断到前 6 个特征标，并附操作/列标签图例。入口为 `format_magnetic_character_table` 与 `format_magnetic_character_table_by_class`。 |
| 目标磁群回归 | `BNS 128.406` 与 `BNS 52.318` 已不再返回“计算不支持”；`128.406@Z` 稳定给出维数 `2,2,2,4` 的四个 coreps，正式操作表包含 `g1..g16` 全部 16 列。 |
| CIR 数据生成 | CIR 解析器支持复合反幺正矩阵：`11,202` 个原始 coreps 中复合反幺正项 `672` 个、拒绝 `0`；`8,388` 个可映射到磁数据库的 coreps 中未映射 `0`。 |
| 验证 | 全 `1651` UNI release summary 审计：成功 `1651`、失败 `0`、`10,390` k 点、`52,793` coreps、π 型放大噪声 `0`；release all-targets 为 lib `291 passed / 4 ignored`、integration `62 passed`，doc-tests `26 passed`。 |

必须保留以下语义边界：

- 上述 `18` 个 operation-only 歧义是缺少 setting 上下文时的真实不可判定性，
  不能简单归类为“官方识别错误”；有 parent/family Hall 后都能唯一确定。
- 当前 summary API 的对象是选定 `k` 点的 magnetic little group corep，不等同于
  full-star 共表示；若将来增加 full-star API，必须使用独立名称和输出语义。
- 个别 Type-A corep 在缺少可构造的 intertwiner/matrix 时，反幺正列会明确显示
  `antiunitary-pending(...)`，不能伪造为已完成的物理特征标；本轮指定的
  `128.406`、`52.318` 回归表不受此问题影响。

### 2026-08-12 Rustb 可复用磁操作层

为 Rustb 的 Hamiltonian 对称性检查与强制对称化新增
`src/operation_group.rs`，公共 Rust-native API 为：

- `ValidatedMagneticOperationSet::try_from_symmetry_ops`：验证非空、有限平移、
  `det(W)=±1`、模晶格唯一性、非加撇恒等元、双侧逆元和完整乘法闭包；错误携带
  operation/product witness，不使用 sentinel 或 panic。
- `ValidatedMagneticOperationSet::identify`：从磁操作自身的空间投影推导 family
  Hall，再做 setting-aware UNI/BNS/OG 识别；调用方给出的 structural-supergroup
  Hall 只记录 provenance，绝不能冒充降低后磁群的 family Hall。
- `axial_spin_half_lift`：把笛卡尔轴矢量旋转提升为 spin-1/2 quaternion，供 Rustb
  构造 `U(W)` 与反幺正 `U(W)iσ_yK`。

`MagneticGroupIdentification` 同时返回 UNI/Litvin/BNS/OG/type、识别 Hall、派生
family Hall、structural provenance Hall，以及完整 setting transform。必须保留
这些字段，不能再折叠成仅有 UNI 的便捷返回。

边界语义：该层只验证/命名调用方给出的完整有限磁操作群，不理解 TB 轨道、
Wannier gauge 或 Hamiltonian。Rustb 负责从 Atom 结构生成候选、验证 localized
basis corepresentation、计算 H 的幸存群或 Reynolds 投影；cryspglib 只负责严格群
代数与命名。完整 release 回归为 `210 passed / 0 failed / 3 ignored`，全部集成测试
与 `26/26` doctests 通过；严格 all-target release clippy 为零警告。

### 2026-08-09 operation-only API 安全化

- `MagneticSpaceGroupType::classify()` 现在返回
  `Result<MagneticSpaceGroupType, SymError>`，识别失败与
  `MagneticUniAmbiguous` 会原样传播，不再伪装成 `UNI=0 / NonMagnetic`。
- 新接口在进入识别算法前检查操作数组非空以及 rotations、translations、
  time-reversals 长度一致；非法输入返回 `SymError::InvalidInput`，不发生索引
  panic。
- C 风格 `spg_get_magnetic_spacegroup_type_from_symmetry()` 已删除；公共入口只保留
  Rust-native `MagneticSpaceGroupType::classify() -> Result<_, SymError>`，不存在
  `UNI=0` sentinel 兼容旁路。
- release 验证：operation-only/structure 集成测试 `13/13` 通过，完整测试套件
  `205 passed / 0 failed / 3 ignored`，doc-tests `27 passed / 0 failed`。

### 2026-08-09 全仓库 Rust/API 安全审查 backlog

在 `classify()` 修复和 release 验证完成后，使用三个并行只读审查分别覆盖公共
API、磁/空间群内核、irrep/corep，并由主线程回查高优先级源码。本节是实时清单；
只有明确标注“已修复”的项目才算完成。

#### P0：可能返回看似有效但实际错误的科学结果

- **已修复**：`Crystal::magnetic_dataset()` 与其 crate-private 识别实现不再在
  `MagneticUniMatchFailed` 时用操作数比例猜测磁类型并返回 `Ok(UNI=0)`；所有磁群
  识别错误现在原样传播，并有稳定的错误传播回归测试。
- **已修复**：磁结构识别入口现在拒绝空结构、positions/types 不等长
  以及长度不等于原子数的 `Some(moments)`，统一返回
  `SymError::InvalidInput`；不再静默当成非磁结构或进入索引 panic。
- **已修复**：`MagneticSpaceGroupType::from_uni()` 现在返回 `Result`；`0`、
  `>1651`（包括 `usize::MAX`）统一返回 `SymError::InvalidInput`，不再伪装成
  `UNI=0 / NonMagnetic`。旧 C 风格 sentinel 兼容函数已删除。
- **已修复**：`wigner_classify()` / `wigner_classify_cir()` 不再把空 unitary
  little group 归为 Type-C，也不再跳过字符缺失、越界索引、错误的 time-reversal
  角色或缺失 Seitz 平方匹配；这些情况统一返回
  `WignerClassificationError`。已补充严格错误回归和有效 Type-A/B/C 量子化测试。

#### P1：公共或数据相关输入可触发 panic/异常分配

- **已修复**：crate-private Hall 识别实现现在拒绝空操作和
  rotations/translations 不等长，返回 `SymError::InvalidInput`；不再发生越界
  panic 或静默丢弃多余 translations。
- **已修复**：普通 `Crystal` 分析现在在 `to_cell()` 边界拒绝空结构以及公开字段
  positions/types/moments 长度失配，统一返回 `SymError::InvalidInput`；低层
  `get_index_with_least_atoms()` 和原胞纯平移入口也显式处理空集合，不再读取
  `mapping[0]` 或发生 `usize` 下溢。
- **已修复**：`SettingTransform::transform_rotation()` 对奇异或数值上不可逆的
  basis 返回 `None`，并由 translation/Seitz 变换继续传播；不再在公开可构造的
  `SettingTransform` 上触发 `expect()` panic。
- **已修复**：`get_changed_pure_translations()` 现在先拒绝零值/非有限行列式、
  非整数或异常大的平移重数，再做有界预分配；负 determinant 与正 determinant
  对称处理，且 `|det|=1` 只有整数矩阵才走快速路径，避免分数 basis 静默漏掉
  晶格平移像。
- **已修复（POSCAR）**：原子计数现在逐 token 严格解析为 `usize`，用
  `checked_add` 求和，并在任何容量分配前按实际坐标行拒绝负数、畸形、溢出、
  空结构、超大或截断输入；不再发生 release 整数回绕、容量溢出或 OOM。
- **已修复（k-mesh）**：Rust 风格高层 API 与 crate-private k 点/网格内核
  统一用 `Result` 拒绝零/负 mesh、非 0/1 shift、checked
  product 超限、输出 slice 过短、BZ map 过短和奇异 reciprocal lattice；地址翻倍
  改用有符号扩展与 Euclidean reduction，`i32::MIN` 不再溢出。分配型网格限制为
  最多 `1,000,000` 点，错误分别为 `InvalidInput` / `ArraySizeShortage`，不再以
  panic、OOM、空结果或错误索引表示失败。
- **已修复（Wigner public helpers）**：
  `wigner_direct_anti_coset()` 改为 `Result`，严格校验 anti-operation 下标、完整
  CIR 字符表、空 coset 和缺失的 square match；spinor direct-anti 入口新增
  operation/character/spin-table 下标与长度验证，并把奇异 setting transform 从
  `expect` 改为 `DirectAntiFailure`。`wigner_classify_spinor()` 也在进入 direct 或
  legacy fallback 前统一验证 spin table 平行数组、`n_lg_ops`、字符表、u16 spin
  index、`a0_idx` 以及 unitary/antiunitary operation roles，legacy 不再能以错误
  下标 panic。`build_corep_chars()` 现在也拒绝越界 magnetic/op-map/H/partner
  映射和截短 Type-A 字符表；三套 Type-A antiunitary helper 在访问前验证
  `a0_idx`、antiunitary role、little-group indices 和矩阵块乘法。
  `debug_unwrapped_square()` 与 `reorder_cir_chars()` 也改为 `Result`，拒绝角色错误、
  越界/溢出 map、奇数或非有限 CIR 字符，不再 panic 或静默补零。该轮审计列出的
  public Wigner operation/index slice 风险已全部关闭。全量 release 初测曾发现严格
  `build_corep_chars()` 把 full-H `op_map` 直接用于 spinor little-group-local 字符表；
  现已在组表前显式执行 `H → global spin op → local character` 域转换，既恢复 UNI
  `21/1066/1510` 的合法结果，也保留所有越界与缺失映射错误。

#### P2：需要计划化处理的 API/idiomatic Rust 债务

- `Crystal` 字段全部公开可变，构造器建立的 positions/types/moments 等长不变量可
  被外部绕过；应考虑私有字段和验证过的 setter。
- **已修复**：`magnetic_irrep_summary_from_ops(uni, ops)` 在读取 H/元数据前验证
  完整无序 magnetic Seitz multiset：rotation/time-reversal 精确匹配，translation
  按模晶格以 `1/12` 和 `1e-5` 容差量化；错误 UNI、缺失、重复、错误 priming、
  非有限/非数据库分数平移和非 first-Hall setting 返回
  `OperationsInconsistentWithUni`。严格的 first-Hall 坐标契约能够区分 UNI 277/284
  的数据库操作集；此前 operation-only 分类器所报告的歧义来自允许 setting 变换后
  的等价性，不能套用到这个 frame-specific API。
- **已修复**：`query::symmetry_operations_of()`、`corep::symmetry_operations_of()`、
  `get_parent_operations()` 和 `canonical_hall_ops()` 均以 `Result` 报告无效 SG 或
  数据库失败，不再返回伪造的空操作集；`matrices_reordered()` 也以
  `MatrixReorderError` 拒绝无法建立的操作映射，不再把 PIR 原顺序冒充 H_ops 顺序。
  字符表格式化器会把操作读取失败写成明确诊断文本。对应 release 回归测试覆盖
  无效 SG、有效操作数/顺序，以及 SG 139 P 点成功与不可映射两条矩阵路径。
- **已修复 Rust-native API 与 clippy 严格门禁**：`generated_data.rs` 此前被 `irrep/mod.rs` 与
  `types.rs` 重复编译，只有后者带生成代码 lint policy；现改为单一模块加 re-export，
  `cargo clippy -p cryspglib --all-targets --release` 从 `16,572` 个
  `approx_constant` error 降为零。随后删除根模块全部公开 `spg_*`/`spgat_*`
  wrapper、`Spglib*` aliases、`SPGLIB_*` constants 与 sentinel 返回路径，将仍需的
  C-port 内核降为 crate-private；Wigner 多参数接口收口为 `KVector`、
  `WignerGroupContext` 与 `SpinorWignerInput`。循环、slice copy、无效初始化、死代码、
  feature-only 变量等警告均已逐项修正，没有使用 blanket lint allow。当前严格命令
  `cargo clippy -p cryspglib --all-targets --release -- -D warnings` 通过，项目警告为零。
  最后一轮使用只读 v4 Flash 审查剩余高参数函数，确认 `WyckoffOutput`、
  `BravaisExpansionOutput`、`SiteEquivalenceContext`、`HallMatchContext` 与
  `MagneticOperationSearch` 的聚合边界；所有修改和验证仍由主线程完成。
- **仍待处理（非 clippy）**：仓库历史源码尚未统一 rustfmt，
  `cargo fmt -p cryspglib -- --check` 会报告大范围既有格式差异；本轮没有为追求格式
  一致性制造全仓机械 diff。该项不影响上述严格 clippy、release tests 或数值审计。

### 远端与起始点

- 2026-07-31 已把本地 `main` 的 161 个提交推到 `origin/main`。
- 本轮起始提交：`4b3f208`。

### 2026-07-31 初始基线

已运行：

```bash
cargo test --package cryspglib test_all_magnetic_sgs_have_valid_operations -- --nocapture
cargo test --package cryspglib diagnose_spglib_standard_setting_transform -- --nocapture
```

结果：

- 数据库操作可读取：`1651 / 1651`。
- `standard_setting_transform` 返回结果：`1651 / 1651`。
- 两条 unitary-SG 识别路径结果一致：`1619 / 1651`，尚有 **32** 个不一致。
- 变换后操作与 detected Hall 完全相等：`1597 / 1651`，尚有 **54** 个不一致。
- 变换后操作与 ISOTROPY data Hall 完全相等：`1450 / 1651`，尚有 **201** 个不一致。

重要：旧测试 `test_all_magnetic_sgs_have_valid_operations` 只检查数据库非空、
旋转行列式和 `ok > 1600`。unitary subgroup 识别失败会被静默跳过，因此它的
“1651/1651 OK”不能作为磁 symmetry 正确性的验收结果。

### 2026-07-31 数据库/群代数严格 gate

新增 `tests/magnetic_symmetry_coverage.rs`，不再只取每个 UNI 的第一个 Hall
setting，也不允许静默跳过。已严格扫描：

- UNI 总数：`1651 / 1651`。
- `(UNI, Hall)` 总数：`4479 / 4479`。
- metadata：UNI 索引、BNS、OG、Litvin 编号均完整且唯一。
- 原始 Seitz 操作：旋转行列式、`1/12` 平移量化、唯一性、单位元、闭包、
  逆元和 `time_reversal` XOR 全部通过。
- Type I/II/III/IV：unitary/antiunitary 阶数、纯时间反演/反平移和 coset
  结构全部通过。
- Type I–III 忽略 time reversal 后与 parent Hall Seitz 集合完全相等；
  Type IV 的 H 使用 doubled magnetic cell，平移集合不一定等于 family Hall，
  因此严格检查其旋转多重集、阶数和反平移扩张关系，全部通过。

本次确认了一个重要语义边界：Type IV 的 BNS 首号及 Hall mapping 描述 family
space group；幺正子群 H 在 doubled magnetic cell 中可能有不同的国际空间群号。
例如 BNS `37.184`、`37.186` 的 H 不能通过与 SG 37 Hall 平移逐项相等来验收。
这不是数据库错误，H 的 SG/Hall 必须在下一层独立识别。

### 2026-07-31 round-trip 初始基线

新增 ignored 诊断 `diagnose_first_hall_database_round_trips`。它从每个数据库磁群的
完整点操作构造不变正定度量和相容晶格，再调用生产路径
`identify_with_parent_hall`，避免用不相容的单位立方晶格制造假失败。

首个 Hall setting 的结果：

- 返回原 UNI 且磁类型正确：`1053 / 1651`。
- `MagneticFallbackReferenceFailed`：`22`。
- `MagneticUniMatchFailed`：`88`。
- 返回错误 UNI、但类型相同：`56`。
- 返回错误 UNI、且类型错误：`432`。

已定位的第一项真实根因：`reduce_to_primitive_magsym` 把非零 anti-translation
当成普通晶格平移消去，并把操作的 time reversal 与该“平移”的 time reversal
做 XOR。这样最简单的 Type-IV BNS `1.3` 会丢失反平移并被识别成 Type-I
BNS `1.1`。

### 2026-07-31 Type-IV 根因修复

已逐行对照官方 spglib v2.5.0（commit `e4531bb`）的
`src/magnetic_spacegroup.c`，确认 Rust 端的预先磁原胞约化不是官方算法的一部分。
FSG/XSG 必须各自在 reference-group 搜索中按普通、非反幺正纯平移约化；不能把
anti-translation 当晶格平移并 XOR 掉 time reversal。Type III/IV 应从
antiunitary coset representative 的线性部分判定：

- `(I|t)'`（通常 `t != 0`）是 Type IV；
- 其余反幺正代表元是 Type III。

删除错误的 `reduce_to_primitive_magsym` 预处理、恢复上述判定后，首 Hall setting
round-trip 提升为：

- 返回原 UNI 且磁类型正确：`1536 / 1651`（原 `1053 / 1651`）。
- `MagneticUniMatchFailed`：`36`（原 `88`）。
- 返回错误 UNI、但类型相同：`79`（原 `56`）。
- 错误磁类型：`0`（原 `432`）。
- fallback reference 失败：`0`（原 `22`）。

结构入口也用同一官方 C API 作了独立 oracle 核对：

- BCC AFM `[111]` 含 `(I|1/2,1/2,1/2)'`，正确结果是 Type IV、
  UNI `1338`、BNS `167.108`，不是旧测试中的 Type III / UNI `1331`。
- FCC FM `[001]` 保留中心化并转到 I-centered tetragonal setting，正确结果是
  Type III、UNI `1197`、BNS `139.537`，不是旧测试中的 primitive
  tetragonal UNI `1005`。

剩余 `115` 个首 Hall round-trip 问题已全部限制在“相同磁类型内的 UNI/setting
匹配”。

### 2026-07-31 官方 oracle 与 setting 手性修复

已在 `/tmp` 独立构建官方 spglib v2.5.0，并用完全相同的 1651 UNI 操作和不变
正定度量逐项调用官方识别 API。官方基线为：

- 精确返回原 UNI：`1648 / 1651`。
- UNI `282`（BNS `37.184`）返回 UNI `275`；
- UNI `283`（BNS `37.185`）识别失败；
- UNI `284`（BNS `37.186`）返回 UNI `277`。

同时机器比较官方 C 数据与 Rust 生成数据：

- 76,683 个编码磁操作逐项完全相等；
- Hall mapping 和 UNI mapping 完全相等。

因此额外偏差不在磁数据库，而在普通空间群 setting 标准化。随后恢复了两处官方
语义：

1. UNI 候选必须等阶并做完整磁 Seitz 集合相等比较，不能接受 subset；
2. `pointgroup::laue_one_axis` 找到四方/三方/六方常规轴后必须检查基变换的
   行列式；若为负，交换前两轴得到右手基。Rust 移植此前在此检查前提前返回，
   会把 enantiomorphic Hall setting 互换，例如 BNS `76.*` 被识别成 `78.*`。

右手化修复后首 Hall round-trip 为：

- 返回原 UNI 且磁类型正确：`1613 / 1651`（修复前 `1536 / 1651`）。
- `MagneticUniMatchFailed`：`36`。
- 返回错误 UNI、但类型相同：`2`（仅官方 oracle 同样混淆的 UNI 282、284）。
- 错误磁类型、fallback 失败、panic：均为 `0`。

完整集合比较一度暴露出结构入口仍依赖旧 subset 匹配：BCC AFM `[111]` 从 24 个
输入操作变换到 rhombohedral Hall 460 时，本应因 `det(T)=1/3` 恢复 3 个纯平移、
合成 72 个操作，Rust 却只生成 24 个。根因是 `magnetic_spacegroup.rs` 的私有
`mat_dmod1` 没有官方 `ZERO_PREC` 容差，把 `-1e-16` 映成接近 `1` 而不是 `0`，
使平移去重计数失败；随后非官方弱 fallback 又静默退回 1 个平移。现已统一使用
`mathfunc::mat_dmod1`，并在计数不符时严格失败，不再产生可被 subset 掩盖的残缺群。
恢复严格等阶/全集比较后，`tests/magnetic_integration.rs` 的 11 个结构入口全部通过，
包括 BCC AFM `[111]` 的 UNI 1338 和 FCC FM `[111]` 的 UNI 1331。

### 2026-07-31 替代 setting 数据恢复与 1651 strict gate

官方 debug oracle 证明 UNI `132` 的 reference transform 与 Rust 在修正前完全一致；
差异实际发生在后续 alternative setting correction。机器比较
`alternative_transformations[][18][7]` 后确认：

- 官方 spglib v2.5.0 有 `450` 个非平凡 `(UNI, Hall)` 替代变换行；
- Rust 生成表只保留 `2` 行，共静默漏掉 `448` 行；
- 旧转换器只接受“显式写满 7 个整数”的 C 初始化器，而 C 大量使用
  `{66459, 0}` 这类 partial initializer，缺省元素按 C 语义应补零，不能丢弃。

现已新增 `scripts/sync_msg_alternative_transformations.py`，对 UNI、Hall 和每行
7 个编码分别做维度校验及零填充，并从官方 C 数据恢复完整表。新增 strict gate
`all_alternative_setting_transformations_are_loaded`，逐项覆盖全部 `4479` 个
setting，固定非平凡行数为 `450`，并验证 UNI `132` / Hall `116` 解码得到官方的
轴交换与 `c/4` 原点平移。

恢复后，不带母群提示的生产识别路径从 `1613 / 1651` 提升为官方基线
`1648 / 1651`；此前 36 个 Rust 特有 `MagneticUniMatchFailed` 全部清零。官方余下
结果是 UNI `282→275`、UNI `283` 失败、UNI `284→277`。后续完整等价类审计证明，
不能把三者统一称为“退化度量歧义”，详见下一节。

### 2026-07-31 Type-IV 跨 parent 等价类与显式消歧

使用官方数据库操作在单位度量和非退化正交度量 `diag(1, 1.3, 1.7)` 下复测，官方
v2.5.0 的 `282→275`、`283→失败`、`284→277` 完全不变，因此问题不由
`a=b=c` 的度量退化触发。进一步对全部 `4479` 个 `(UNI, Hall)` setting 做原始
Seitz 集合和 Type-IV 规范化集合分组，得到两个不同层次的结论：

- 原始完整磁 Seitz 集合已有跨 UNI 的逐项完全重复：UNI `275` Hall
  `177/179/181` 分别等于 UNI `282` Hall `182/183/184`；只给操作和同一晶格时，
  这三对不可能唯一确定 BNS parent。
- 对全部 Type-IV setting 做同一规范化后，跨 UNI 等价类只有
  `{275, 282}` 和 `{277, 284}` 两类；没有第三类，也没有散落在其他 UNI 的遗漏。
- UNI `283` 的规范类唯一。官方识别失败是单一 XSG Hall 候选裁剪造成的可修复错误；
  官方对 `282/284` 静默返回 `275/277` 则是在真实等价类中擅自选择一个代表，API
  应报告歧义而不是把代表当成唯一答案。

生产实现现在从数据库自动构造完整 Type-IV 规范类索引，不写 UNI 特例：

1. 无母群提示时，跨 UNI 类返回新的 `MagneticUniAmbiguous`；唯一类继续返回 UNI，
   因而 UNI `283` 的 Hall `182–184` 均可正确恢复。
2. 提供 parent Hall 时，按其非磁母群空间群号筛选规范类，再组合输入和数据库的
   规范化变换，支持一般基变换和任意原点移动，不再局限于 canonical-operation
   快路径。
3. 所有候选仍必须经过完整、等阶、双向磁 Seitz 集合相等验证；不恢复 subset 匹配。

当前正式 gate：

- `all_database_settings_round_trip_with_parent_hint`：全部 `4479 / 4479`
  `(UNI, Hall)` setting 精确返回原 UNI、原 Hall 和原磁类型；
- `automatic_all_setting_round_trips_are_unique_or_explicitly_ambiguous`：无提示路径
  `4461 / 4479` 唯一且精确，另有 `18 / 4479` 明确返回 ambiguity；18 项恰为
  UNI `275/277` 的各 6 个 setting 与 UNI `282/284` 的各 3 个 setting，合计仍严格
  覆盖 `4479 / 4479`，没有失败或静默错配；
- `type_iv_parent_hall_disambiguates_a_changed_basis_and_origin`：任意原点移动和右手
  轴变换后的同一 Type-IV 集合，可由 parent SG 36/37 分别稳定选择 UNI 275/282；
- `type_iv_orthorhombic_metric_recovers_unique_283_and_reports_real_ambiguities`：非退化
  正交度量下 UNI 283 唯一恢复，282/284 显式报告歧义；
- `all_magnetic_database_operations_form_expected_groups`：`1651` UNI、`4479`
  settings 的群代数及 Type I–IV 结构全部通过。

本轮最终验证：

```bash
cargo check --release --package cryspglib
cargo test --release --package cryspglib --test magnetic_symmetry_coverage -- --nocapture
cargo test --release --package cryspglib --test magnetic_integration -- --nocapture
cargo test --release --package cryspglib --tests
```

- 磁 symmetry strict gate：`6 passed / 0 failed`；
- 结构磁矩入口：`11 passed / 0 failed`；
- 全部 test targets：`202 passed / 0 failed / 2 ignored`（其中 lib tests
  `149 passed`，两个 ignored 均为显式诊断项）；
- `cargo check --release` 通过；现有 warning 未在本轮扩散处理。

重新运行 `diagnose_spglib_standard_setting_transform` 后的当前上层 setting
oracle 为：`total=1651`、`found=1651`、`sg_match=1651`、
`detected_hall_exact=1597`、`data_hall_exact=1450`。因此旧基线中的 32 个
unitary-SG 路径不一致已经清零。进一步把错误的“等阶全集相等”检查改为
centered-cell 允许的严格 Seitz embedding 后，`detected_hall_embed=1651`：
54 个 detected Hall exact 差异全部只是 primitive 代表元嵌入 C/F/I-centered
常规胞时的合法操作数展开。对 ISOTROPY data-Hall 复合 frame 做同一严格检查，
最终 `data_hall_embed=1651 / 1651`。修复包含两个通用问题：

- 禁止在 data-Hall 复合失败后把 MSG→detected 变换误标成 MSG→data 变换，
  并始终尝试经完整磁操作验证的直接 MSG→data 搜索；
- `find_setting_transform` 的 origin solver 候选现在必须通过完整 Seitz 集验证
  才能返回，不能再把 SG7 Hall `22/23→21` 错误短路成单位变换；正确轴交换/shear
  会继续被枚举，UNI `296/312` 的 4→2 primitive embedding 也由复合变换覆盖。

聚焦回归 `setting_transform_rejects_invalid_identity_origin_candidate` 固定 Hall
`23→21` 不得返回伪单位变换；全库 oracle 同时严格断言 `total`、`found`、
`sg_match`、`detected_hall_embed`、`data_hall_embed` 均为 `1651`。

### 已确认的问题

1. `primitive.rs` 曾无条件输出 `reduced=...`，1651 群扫描产生大量噪声。
2. `hall_symbol.rs` 曾对 Hall 497 无条件输出内部匹配跟踪。
3. 数据库/群代数/磁类型层已由 4479-pair strict gate 清零。
4. 数据库磁操作 `→` 磁群识别 `→` 原 UNI 的全部 setting round-trip 在显式母群
   Hall 下为 `4479 / 4479`；无提示自动路径为 `4461` 个唯一精确结果加 `18` 个
   `MagneticUniAmbiguous`，完整覆盖 `4479 / 4479`。歧义只来自数据库自动导出的
   `{275,282}`、`{277,284}` 两个跨 parent Type-IV 等价类；UNI 283 已唯一恢复。
5. irrep/ISOTROPY setting oracle 的 unitary-SG 和 detected-Hall 严格 embedding
   均为 `1651 / 1651`；54 个 detected exact 差异已确认是合法 centering 展开，
   data-Hall 严格 embedding 也已达到 `1651 / 1651`。
6. `magnetic_summary` 原有的 BNS `128.406` / `52.318` unsupported scalar
   路径已清零（2026-07-31）：MSG 操作、H/PIR/CIR 和 k 向量现在统一进入
   ISOTROPY data-Hall frame；多星臂标量 PIR 从完整诱导矩阵抽取所选 k 星臂块，
   不再用整颗星的维数和 trace 做 little-group Wigner 判别。两群的全部列均满足
   有限值、列/操作对齐和 `χ(E)=dim` 的 release 回归。
7. CIR 生成链路已按完整文件结构重写并通过严格校验：`CIR_data.txt` 的 k 记录数
   使用 `star_count * 16`，接受可选四整数 `irtranslation`，不再为 `irtype=2`
   虚构额外矩阵；现可连续解析全部 `11202` 条 CIR。672 个 compound PIR 均保存
   第一星臂的复 CIR trace，并在 Seitz/Hall 重排及 centered-cell 扩张时同步施加
   Bloch 相位。PIR 中 `P1P1`/`P1PA1` 一类标签歧义统一走 `_lookup_kvec`，避免把
   边界 k 点误当作 Γ。生成日志为 `672 compound / 0 rejected`、`8388 mapped /
   0 unmapped`。
8. `128.406@Z` 的 compound 维数重复加倍已修复：Wigner 分类以一个不可约复 CIR
   分量为起点，再按 Type-C 构造共表示，当前四行维数严格为 `2, 2, 2, 4`；旧的
   “必须仍返回 structured error” 集成测试已改为 BCS 成功契约。
9. 特征标格式化不再截断为前 6 项：默认 Markdown 表以全部 magnetic little-group
   操作为列，并附 MSG 索引、unitary/antiunitary 标记和 data-Hall Seitz 定义；另有
   按 character-compatible 共轭类分列的正式表格入口。`128.406@Z` 回归明确要求
   `g1..g16` 全部出现且不得含省略号。
10. 最新 release 全量 summary gate：`success=1651`、`failure=0`、`kpoints=10390`、
    `coreps=52793`。gate 对每个结果验证 operation/character 列数、共轭类分割、
    有限值及 `χ(E)=dim`，不是只检查 API 返回 `Ok`。
11. 2026-08-10 Rust-native API 收口后 `cargo test --package cryspglib --release`：
    lib `196 passed / 0 failed / 3 ignored`，七个 integration binaries 合计
    `58 passed / 0 failed`，doc-tests `26 passed / 0 failed`；即常规 release tests
    共 `254 passed`，外加 `26` 个 doctests。减少的一项 integration test 与一个
    doctest 均属于已删除的 C 风格 sentinel wrapper，不是功能覆盖退化。1651 全量
    summary ignored release gate 另行执行并通过。

### 分层验收标准

必须按以下顺序清零，不能用后层成功掩盖前层失败：

1. **数据库层**：全部 UNI/Hall pair 均有合法 metadata 和非空操作。
2. **群代数层**：单位元、唯一性、闭包、逆元、`time_reversal` XOR 同态全部通过。
3. **磁类型层**：Type I/II/III/IV 的 unitary/antiunitary 阶数与 coset 结构正确。
4. **family group 层**：Type I–III 忽略 time reversal 后与对应 parent Hall
   的 Seitz 集合完全一致；Type IV 必须按 doubled magnetic cell 验证 family
   扩张，不能错误要求 H 的平移逐项等于 family Hall。
5. **unitary subgroup 层**：H 的 SG/Hall/setting transform 可复现且操作集合完全一致。
6. **round-trip 层**：数据库操作经公开/生产识别路径返回原 UNI，而不只是相同 type。
7. **结构入口层**：代表性非正交、非零 origin、中心化和 Type IV 结构走
   `Crystal::magnetic_dataset()` 得到一致结果。
8. **上层表示层**：所有可支持的 k 点/corep/summary 不产生静默跳过、`NaN`
   或伪 `Unsupported` 行；真正未实现项必须返回结构化错误。

### 工作纪律补充

- 全扫必须报告总样本数和每个失败类别；不接受 `>1600` 之类宽松阈值。
- 诊断测试可以暂时记录基线，但正式 gate 最终必须要求零失败。
- 当前全仓库 `cargo fmt --check` 会在超大生成数据上构造巨型 diff 并 OOM；
  修改手写 Rust 文件时保持改动 hunk 的 rustfmt 风格，运行 `cargo check` 和
  相关测试。生成数据格式问题单独处理，不能因此跳过编译与测试。

---

## Irrep 终极目标与详细实施计划

最终目标：给定一个磁空间群（结构识别得到的 UNI、直接输入 UNI、或 BNS/OG 标签），程序应能一次性回答：

1. 这个磁群有哪些高对称 k 点/线/面，以及每个标签的标准名称和坐标。
2. 这个磁群在这些 k 点上的磁共表示（corepresentation）：来源 H-irrep、corep type A/B/C、维数、字符、antiunitary 完整性、Hall/setting 约定。
3. 每个磁共表示对应的可能 isotropy subgroup：普通 isotropy subgroup 和 magnetic isotropy subgroup 都要能追溯到 source irrep/corep、k 点、方向、domain/arm 信息。

### 当前完成度（2026-07-03）

已完成或基本可用：

- `src/irrep/query.rs`: `irreps_of(sg)`, `kpoints_of(sg)`, `IrrepRecord::k_label()` 已能列出 230 个 SG 的 ISOTROPY 高对称 k 标签、坐标和 irreps。
- `src/irrep/types.rs` + `generated_data.rs`: `IrrepRecord::subgroups()` 和 `IrrepRecord::magnetic_subgroups()` 已保存普通/磁 isotropy subgroup 数据。
- `src/irrep/corep.rs`: `compute_corepresentation()` / `compute_coreps(bns, k_label)` 已实现 scalar PIR、scalar CIR、spinor SU(2) Wigner 分类；corep 计算 API 使用 `Result`，无法分类/非有限字符必须返回结构化错误，不能用 `Option` 或 `NaN` 占位。
- 最新 spinor Wigner 全扫诊断已清零失败：`cargo test --package cryspglib diagnose_wigner_sources -- --nocapture` 中 `spinor_complex_ok = 21216`，无 `spinor_complex_fail`。
- `magnetic_isotropy_coreps_of_irrep()` / `magnetic_isotropy_coreps_of_sg_k()` 已有 corep 与 magnetic isotropy 的早期桥接，但还不是最终用户 API。

未完成的关键缺口：

- 缺少一个面向用户的统一入口：当前能力分散在 `query.rs`, `corep.rs`, `api.rs`, `magnetic_spacegroup.rs`。
- 磁群 k 点语义需要固定：第一版采用 unitary subgroup H 的 ISOTROPY k 点作为标准列表，同时显式返回 magnetic little group 的 unitary/antiunitary 阶数；后续再处理 antiunitary 合并 k-star 的展示策略。
- Type-C corep 需要 pairing/dedup：不能把一对互为 antiunitary 共轭的 H irreps 重复展示成两个磁 coreps。
- corep 到 isotropy subgroup 的映射要明确为“候选”：第一版把 source H-irrep 的 ordinary/magnetic isotropy subgroup 合并挂到 corep；后续再验证 Type-C/compound/spinor 的物理筛选规则。

### 目标 API（新增）

新增文件：`src/irrep/magnetic_summary.rs`。

对外入口：

```rust
pub fn magnetic_irrep_summary(input: MagneticIrrepInput)
    -> Result<MagneticIrrepSummary, MagneticIrrepError>;

pub fn magnetic_irrep_summary_by_uni(uni: usize)
    -> Result<MagneticIrrepSummary, MagneticIrrepError>;

pub fn magnetic_irrep_summary_by_bns(bns: &str)
    -> Result<MagneticIrrepSummary, MagneticIrrepError>;

pub fn magnetic_irrep_summary_from_ops(
    uni: usize,
    mag_ops: &crate::SymmetryOps,
) -> Result<MagneticIrrepSummary, MagneticIrrepError>;
```

输入类型：

```rust
pub enum MagneticIrrepInput<'a> {
    Uni(usize),
    Bns(&'a str),
    Operations { uni: usize, ops: &'a crate::SymmetryOps },
}
```

核心返回结构：

```rust
pub struct MagneticIrrepSummary {
    pub uni: usize,
    pub bns_label: String,
    pub magnetic_type: crate::MagneticType,
    pub parent_sg: u8,
    pub unitary_sg: u8,
    pub unitary_hall: usize,
    pub kpoints: Vec<MagneticKPointSummary>,
}

pub struct MagneticKPointSummary {
    pub label: String,
    pub coords: (i8, i8, i8, i8),
    pub little_group_order: usize,
    pub unitary_order: usize,
    pub antiunitary_order: usize,
    pub operations: Vec<MagneticLittleGroupOperation>,
    pub conjugacy_classes: Vec<MagneticConjugacyClass>,
    pub coreps: Vec<MagneticCorepSummary>,
}

pub struct MagneticCorepSummary {
    pub label: String,
    pub source_irreps: Vec<SourceIrrepSummary>,
    pub corep_type: crate::irrep::corep::CorepType,
    pub source: crate::irrep::corep::WignerSource,
    pub dim: usize,
    pub characters: Vec<Option<num_complex::Complex64>>,
    pub timerev: Vec<bool>,
    pub completeness: crate::irrep::corep::CharacterCompleteness,
    pub isotropy_candidates: Vec<CorepIsotropyCandidate>,
}

pub struct SourceIrrepSummary {
    pub sg: u8,
    pub ml: &'static str,
    pub bc: &'static str,
    pub dim: u8,
    pub spinor: bool,
}

pub struct CorepIsotropyCandidate {
    pub source_ml: &'static str,
    pub ordinary: Vec<crate::irrep::types::IsotropyRecord>,
    pub magnetic: Vec<crate::irrep::types::MagneticIsotropyRecord>,
    pub relation: IsotropyCandidateRelation,
}

pub enum IsotropyCandidateRelation {
    DirectSourceIrrep,
    TypeCPairedSource,
    CompoundSource,
    SpinorNoIsotropyData,
}
```

错误类型：

```rust
pub enum MagneticIrrepError {
    InvalidUni(usize),
    UnknownBns(String),
    MissingMagneticOperations(usize),
    MissingUnitarySubgroup(usize),
    MissingIrrepData { sg: u8 },
    CorepComputationFailed { uni: usize, sg: u8, k_label: String },
}
```

Re-export 路径：

- `src/irrep/mod.rs`: `pub mod magnetic_summary;`
- `src/irrep/mod.rs`: re-export 常用类型，或要求用户显式 `use cryspglib::irrep::magnetic_summary::*;`
- `src/lib.rs`: 暂不顶层 re-export，等 API 稳定后再决定是否暴露到 crate root。

### 实现顺序与路径

#### Phase 1: 只读 summary API 骨架

路径：`src/irrep/magnetic_summary.rs`, `src/irrep/mod.rs`。

实现方法：

1. 新增上面的 public structs/enums。
2. `magnetic_irrep_summary_by_uni(uni)`:
   - 校验 `1 <= uni <= 1651`。
   - `let mag_ops = SymmetryOps::from_magnetic_database(uni)?`。
   - 调用 `magnetic_irrep_summary_from_ops(uni, &mag_ops)`。
3. `magnetic_irrep_summary_by_bns(bns)`:
   - 复用 `corep.rs` 现有 BNS→UNI helper；若 helper 当前私有，移动/改为 `pub(crate)`。
   - 转入 `magnetic_irrep_summary_by_uni(uni)`。
4. `magnetic_irrep_summary_from_ops(uni, mag_ops)`:
   - 调用 `identify_unitary_subgroup_with_hall(uni)` 取得 H 信息。
   - 记录 `unitary_sg`, `unitary_hall`, `msg_to_data`。
   - 遍历 `query::kpoints_of(unitary_sg)` 生成 `MagneticKPointSummary`。

验收测试：

- `magnetic_summary_by_uni_smoke`: UNI 1599 或 BNS `221.97` 返回非空 kpoints。
- `magnetic_summary_by_bns_matches_uni`: BNS 和 UNI 两条路径返回相同 `uni/unitary_sg/kpoints.len()`。

#### Phase 2: magnetic little group 元数据

路径：`src/irrep/magnetic_summary.rs`，复用 `src/irrep/wigner.rs`。

实现方法：

1. 提取 helper:

```rust
fn canonical_pure_translations(h_ops: &crate::SymmetryOps) -> Vec<[f64; 3]>;

fn magnetic_little_group_indices(
    k: (i8, i8, i8, i8),
    mag_ops: &crate::SymmetryOps,
    setting_xf: Option<&crate::irrep::wigner::SettingTransform>,
    canonical_translations: &[[f64; 3]],
) -> Vec<usize>;
```

2. 对每个 k 点调用 `filter_little_group_with_transform`。
3. 填充 `little_group_order`, `unitary_order`, `antiunitary_order`。
4. 第一版 k 点列表保持 H 的 ISOTROPY k 点，不试图减少/合并 magnetic star；API 文档明确这一点。

验收测试：

- grey group: 每个有 antiunitary 的 k 点 `antiunitary_order > 0`。
- nonmag/type-I equivalent: `antiunitary_order == 0` 时 corep type 应走 trivial A。

#### Phase 3: corep 计算接入 summary

路径：`src/irrep/magnetic_summary.rs`, 必要时调整 `src/irrep/corep.rs` helper 可见性。

实现方法：

1. 对每个 `KPointSummary.irreps` 找到 H 的 `IrrepRecord`。
2. 调用：

```rust
corep::compute_corepresentation(ir, uni, mag_ops)
```

3. 把 `Corepresentation` 转成 `MagneticCorepSummary`。
4. `label` 第一版使用 source ML label；Type-C dedup 完成后改为组合 label。
5. 保留 `CharacterCompleteness`；无法分类/unsupported 必须返回 `Err`，不生成带 `NaN` 的伪 corep。

验收测试：

- `128.406` 和 `52.318` 必须完整成功；其他真正 unresolved case 仍应返回结构化
  `CorepComputationFailed`，不能输出 `Unsupported` + `NaN` 作为结果。
- `221.97` at `GM`: 返回非空 coreps，identity character 等于 dim。

#### Phase 4: Type-C pairing/dedup

路径：`src/irrep/magnetic_summary.rs`，必要时在 `src/irrep/corep.rs` 增加 `pub(crate)` helper。

新增内部类型：

```rust
struct CorepDedupKey {
    corep_type: CorepType,
    dim: usize,
    rounded_characters: Vec<i64>,
    timerev: Vec<bool>,
}
```

实现方法：

1. 第一版 dedup 使用 `corep_type + dim + rounded characters + timerev`，字符按 `1e-8` 量化。
2. Type-C 合并时 `source_irreps` 追加两个 H-irrep。
3. `label` 规则：
   - 单 source: `"GM4-"`。
   - Type-C pair: `"Z1Z4 + Z2Z3"` 或按 ML 排序 join。
   - compound source 保留原 compound ML label。
4. 后续增强：用 antiunitary conjugation 显式找 partner，而不是只靠 character key。

验收测试：

- `test_type_c_coreps_are_deduplicated` 的语义迁移到 summary API：Type-C pair 不重复。
- 每个 `MagneticCorepSummary.source_irreps` 非空；Type-C 至少能出现两个 source 的 case。

#### Phase 5: isotropy candidates 挂接

路径：`src/irrep/magnetic_summary.rs`。

实现方法：

1. 对每个 `source_irrep` 收集：
   - `ir.subgroups()` → ordinary candidates。
   - `ir.magnetic_subgroups()` → magnetic candidates。
2. 普通 subgroup 去重 key:

```rust
(sg, symbol, direction, domains, arms)
```

3. 磁 subgroup 去重 key:

```rust
(mag_sg, bns_label, direction)
```

4. `relation` 规则：
   - 单 source scalar: `DirectSourceIrrep`
   - Type-C 合并 source: `TypeCPairedSource`
   - `cir_component_count() > 0`: `CompoundSource`
   - spinor 且没有 isotropy 数据: `SpinorNoIsotropyData`
5. 第一版只声明 candidates，不声明这些 subgroup 已按 corep order parameter 方向完成物理筛选。

验收测试：

- SG 221 GM4- 路径能返回包含 ordinary/magnetic subgroup 的 candidates。
- Type-C 合并后 candidates 去重稳定，不因 source 顺序变化而重复。

#### Phase 6: 格式化与示例

路径：`src/irrep/magnetic_summary.rs`, `README.md` 或 `examples/`。

新增格式化 API：

```rust
pub fn format_magnetic_irrep_summary(summary: &MagneticIrrepSummary) -> String;
pub fn format_magnetic_kpoint_summary(kpoint: &MagneticKPointSummary) -> String;
pub fn format_magnetic_character_table(kpoint: &MagneticKPointSummary) -> String;
pub fn format_magnetic_character_table_by_class(kpoint: &MagneticKPointSummary) -> String;
pub fn format_magnetic_character_table_with_columns(
    kpoint: &MagneticKPointSummary,
    columns: MagneticCharacterTableColumns,
) -> String;
```

格式化器不再只预览前 6 个特征标。默认输出完整逐操作 Markdown 表，并附每列的
MSG 原始操作索引、unitary/antiunitary 类型及 data-Hall Seitz 操作；另一入口按
磁 Seitz 共轭类输出。若 Bloch/projective 字符在原始类内不恒定，类表会按完整
corep character signature 自动细分，禁止错误合并。

示例文件：

- `examples/magnetic_irrep_summary.rs`

示例目标：

```rust
let summary = magnetic_irrep_summary_by_bns("221.97")?;
for kp in &summary.kpoints {
    println!("{} {:?}", kp.label, kp.coords);
    for c in &kp.coreps {
        println!("  {} {:?} dim={}", c.label, c.corep_type, c.dim);
    }
}
```

验收测试：

- example 能编译运行。
- 格式化输出不依赖 HashMap 随机顺序。

### 推荐提交顺序

1. `feat: add magnetic irrep summary types`
2. `feat: summarize magnetic k-points by UNI`
3. `feat: attach corepresentations to magnetic summary`
4. `feat: deduplicate type-c magnetic coreps`
5. `feat: attach isotropy candidates to magnetic coreps`
6. `docs: add magnetic irrep summary example`

每个提交前运行：

```bash
cd /home/liuyichen/TB_rs
cargo check --package cryspglib
```

关键节点额外运行：

```bash
cargo test --package cryspglib diagnose_wigner_sources -- --nocapture
cargo test --package cryspglib test_type_c_coreps_are_deduplicated -- --nocapture
cargo test --package cryspglib --tests
```

### 非目标（先不要做）

- 不在第一版里重新生成 ISOTROPY 数据。
- 不在第一版里发明新的磁 k 点命名系统；先复用 H 的 k labels，并暴露 magnetic little group 元数据。
- 不在第一版里声称 isotropy candidates 已完成 order-parameter 方向的唯一筛选；先返回可追溯候选。
- 不把 summary API 暴露到 crate root；先稳定在 `cryspglib::irrep::magnetic_summary`。

---

## Workspace context

This crate is a **workspace member** inside `/home/liuyichen/TB_rs`. All cargo commands must be run from the workspace root:

```bash
cd /home/liuyichen/TB_rs
cargo build --package cryspglib
cargo test  --package cryspglib
cargo check --package cryspglib
```

---

## 工作流铁律

### 规则 1: 每完成一个可编译的修改就立即 commit

**每次 `cargo check` 成功后必须立即 commit**，然后再做下一个修改。不要连续做多个修改才 commit。

提交前必须格式化：

```bash
cargo fmt --package cryspglib
```

```bash
git add -A && git commit -m "描述"
```

**Why:** `git checkout` 恢复时只保留已提交的内容。中间修改全部丢失。宁可 commit 太多（事后 squash），不能丢失工作。

**反面案例（发生过两次）：**
1. 2163→3307 行的诊断代码因 git checkout 全部丢失
2. MagneticOps→SymmetryOps 重构中，已修复的 10+ 处 field access 因一次 revert 全部丢失

### 规则 2: 不用 Python 脚本做代码修改

批量 sed/Python 替换容易产生意外后果。代码修改必须用 Edit/Write 工具逐处进行，每处修改后确认正确。

### 规则 3: 不用 type alias 做过渡

`pub type OldName = NewName;` 只是把问题藏起来，缺少可维护性。应该全局替换所有引用，然后删除旧定义。

### 规则 4: 先 use 再去掉 crate:: 前缀

替换类型时，先在文件头部添加 `use crate::NewType;`，再把文件中所有 `crate::NewType` 替换为 `NewType`。

### 规则 5: 先验证输入数据语义，再调算法

当大量 case 出现同一症状时，根因通常是**对输入数据语义的共同错误假设**，
而不是算法本身的多个独立 bug。**在动任何公式之前**，先问：

> "这个输入字段对我的场景语义正确吗？它的含义真的是我以为的那样吗？"

具体做法：
1. **从一级数据推导，不信任二级元数据。** Hall ops 和 spin ops 是一级数据；
   isotropy subgroup 的 origin 字段和手维的 SG→centering 表是二级数据（可能
   对当前 setting 是错的）。
2. **用命名常量定义容差**（如 `const SEITZ_TRANS_TOL: f64 = 1e-5`），
   不要到处撒 `1e-9`。容差取值基于管线中浮点误差累积的预期量级，不是越紧越好。
3. **如果修了 2 次问题反而恶化，说明假设错误。** 停止调算法，回头检查输入数据语义。
4. **提取共享 helper**——如果同一 inline 逻辑出现在多个调用点，提取一次，统一使用。

**历史案例：** spinor Wigner square mismatch (679→0) 是 Codex 用 3 个 commit
解决的，完全没有动 Wigner 公式本身：(1) isotropy 来源的 origin 对 Wigner square
是错的 → 从 Hall vs spin ops 直接求解；(2) 按 SG 号手列的 centering 表不考虑
Hall setting 坐标轴排列 → 从 Hall 纯平移自动推导；(3) `1e-9` 对 translation
浮点误差太紧 → 改为 `1e-5` 命名常量。三个都是数据语义问题，不是算法问题。

---

## ISOTROPY 数据格式知识

### PIR vs CIR 的本质区别

- **CIR**（Complex Irreducible Representation）：只依赖**旋转矩阵 R**，不依赖 translation。
  每个独特的旋转类型有一个 complex character `χ(R) = (re, im)`。
  opcount = little co-group 的大小（distinct rotation types）。

- **PIR**（Physically Irreducible Representation）：依赖完整的空间群操作 `{R|t}`。
  字符通过 CIR + Bloch 相位组合：`PIR({R|t}) = Σ_i CIR_i(R) * exp(i*2π*k·t)`
  opcount = full little group 的大小（包含所有 translation 变体）。

- ISOTROPY 的 PIR 和 CIR 操作数不同是结构性的，不是数据缺失。
  CIR 永远只有 distinct rotation types 的条目。

### CIR_data.txt 格式

```
seq sg "symbol" "label" dim irtype kcount pmkcount opcount
<kcount × 16 integers>        ← CIR 的 k-star augmented 4×4 records
<16-int Seitz operation>      ← augmented 4×4 {R|t}
<optional 4-int irtranslation>← 仅 non-special k vector
<dim² (re,im) values>         ← complex IR matrix, row-major
...                           ← 对每个 operation 重复
```

- CIR 跳过 `kcount × 16` 个 k-vector 整数；PIR 对应跳过
  `pmkcount × 16`。这一区别由官方 `CIR_data.f` / `PIR_data.f` 明确定义。
- operation 的 16 个整数是 augmented Seitz matrix，不是复表示矩阵编码。
- CIR 也可能有 `irtranslation`；官方 reader 通过 k-vector record 是否含自由参数
  判定 special/non-special，而不是假设 CIR 永远没有 translation 行。
- 文件存储的是 `dim²` 个表示矩阵元素；character 是生成器对对角元求 trace 得到的，
  不是文件中的独立 character 行。

### spglib Hall vs ISOTROPY 的 translation 差异

- ISOTROPY 使用 **primitive cell** 的 translation
- spglib Hall 使用 **conventional cell** 的 translation（可能含 centering）
- 两者 translation 不同，但**旋转矩阵相同**
- 重排时：PIR 字符是实数，只需排列不需相位修正（ISOTROPY PIR 字符对 matched operation 是正确的）
- CIR 展开时：额外 Hall 位置需从 matching rotation 复制 + Bloch 相位 `exp(i*2π*k·t_hall/kd)`

### 重排后数据一致性

重排后 PIR 和 CIR 在同一 Hall 位置的字符来自**同一个 ISOTROPY operation**
（因为 mapping[h] 唯一定义了 ISOTROPY 索引）。因此 PIR = CIR_sum 在排列后
仍然成立。CIR 展开只影响额外位置（CIR 源数据没有的 translation 变体）。

### 2026-08-27 角色数据审计：π 型异常已修复

本节记录 Rustb `calculate_irrep` 接入过程中对 cryspglib 角色数据做的全量审计。
π 型异常已经定位到生成 Rust literal 的最后一步，并已从根源修复。后续修改生成器、
CIR/PIR 映射或角色拟合时仍必须把本节的回归作为独立 gate；不得通过按数值模式
替换生成数组来掩盖数据问题。

最终修复提交为 `ba2d997 fix: preserve scientific notation in irrep data`，已推送至
`origin/main`。

#### 审计范围与已确认的健康边界

- 扫描全部 230 个空间群中的 `4777` 个 scalar record 和 `3611` 个 spinor record，
  共 `8388` 个 irrep record；
- 扫描全部 `1651 / 1651` 个 magnetic summary，共 `972786` 个 summary
  character field；
- 所有 magnetic summary 均可构造，未发现 NaN、Inf 或负零；
- 普通 PIR `characters()`、`SCALAR_LITTLE_CHARS_REAL` 和
  `SPIN_IMAG_CHARS` 中没有下述 π 型异常；异常集中在 CIR 派生的复角色数据。

修复前的问题不是整个数据库普遍损坏，而是集中在
`SCALAR_LITTLE_CHARS_IMAG`、`CIR_COMPONENT_CHARS` 及其少量下游 summary；
修复后同一全量扫描中 π 型字段为 `0`。

#### 已修复：CIR 复角色中的 758 个 π 型数值

修复前全量扫描发现 `758` 个静态字段接近以下数值：

```text
π/30, π/(10√3), π/15, π/10, √3π/10, π/5
```

它们的分布为：

| 字段 | 数值模式 | 数量 | 涉及空间群 |
|------|----------|-----:|------------|
| `CIR_COMPONENT_CHARS` | `π/30` | 88 | 144, 145, 152, 154, 169, 170, 178, 179 |
| `CIR_COMPONENT_CHARS` | `π/10` | 68 | 146, 148, 161, 167 |
| `CIR_COMPONENT_CHARS` | `π/5` | 8 | 161, 167 |
| `CIR_COMPONENT_CHARS` | `√3π/10` | 8 | 167 |
| `SCALAR_LITTLE_CHARS_IMAG` | `π/30` | 200 | 144, 145, 151, 152, 153, 154, 171, 172, 178, 179, 180, 181 |
| `SCALAR_LITTLE_CHARS_IMAG` | `π/(10√3)` | 8 | 178, 179 |
| `SCALAR_LITTLE_CHARS_IMAG` | `π/15` | 4 | 178, 179 |
| `SCALAR_LITTLE_CHARS_IMAG` | `π/10` | 354 | 146, 148, 155, 160, 161, 166, 167 |
| `SCALAR_LITTLE_CHARS_IMAG` | `π/5` | 20 | 155, 160, 166, 167 |

代表性原始值包括：

```text
0.1047197638   ≈ π/30
0.1813809141   ≈ π/(10√3)
0.2094395276   ≈ π/15
0.3141594877   ≈ π/10
0.5441430822   ≈ √3π/10
0.6283189753   ≈ π/5
```

这类值不是 `0.70711 ≈ 1/√2` 或 `1.73206 ≈ √3` 那样的低精度代数数。
有限 little co-group（包括具有有限因子的投影表示）的角色是根单位的有限和，
应属于 cyclotomic 代数数域；精确包含 π 的超越数不能是正确角色。因此这批值
高度疑似将相位角、记录边界或其他编码字段误当成了复数虚部。

异常不是只存在于未使用的静态表。修复前有 `40` 个字段传播进公开 magnetic
summary：

| UNI | BNS | 受影响模式 |
|----:|-----|------------|
| 1300 | 161.70 | `π/5` |
| 1302 | 161.72 | `π/5` |
| 1335 | 167.105 | `π/10` |
| 1338 | 167.108 | `√3π/10` |
| 1344 | 169.114 | `π/15` |
| 1346 | 169.116 | `π/15` |
| 1348 | 170.118 | `π/15` |
| 1350 | 170.120 | `π/15` |
| 1389 | 178.159 | `π/30` |
| 1395 | 179.165 | `π/30` |

本次修复没有把这些数按“最接近常量”snap 到猜测值，而是恢复官方 archive、
复现完整生成管线并修正造成数量级放大的 formatter。修复后恰好只有上述
`586 + 172 = 758` 个字段变化，其他所有 `f64` 生成数组逐项不变。

#### 根因：`rstrip("0")` 破坏科学计数法指数

以 SG144 A2 为最小见证，Hall translation 与 ISOTROPY translation 的十进制表示
分别为 `1/3` 与 `0.3333333333`。相位对齐后正确的虚部只是浮点噪声：

```text
1.0471976378421116e-10
```

旧 `_fmt_char` 先得到字符串 `"1.047197638e-10"`，随后对整个字符串调用
`.rstrip("0")`，把指数末尾的零也删除成 `"1.047197638e-1"`。因此正常的
`O(10^-10)` 相位噪声被放大九个数量级，表现为 `0.1047197638 ≈ π/30`。
所有 π/10、π/5 及根式倍数家族都来自同一字符串错误。

修复后的 `_format_rust_f64` 不再裁剪格式化字符串，并将 literal 重新解析做
`math.isclose` 往返校验；任何未来的数量级破坏会在生成阶段直接失败。
`tests/irrep_validation.rs::sg144_a2_little_characters_keep_tiny_phase_noise`
固定了数据库端的最小回归。

#### 官方源数据与格式核验

2026-08-27 从官方 ISO-IR 页面恢复 2022 computer-readable archives，并对照归档
自带的 Fortran reader，而不是根据生成数组反推格式：

| archive | SHA-256 |
|---------|---------|
| `PIR_data.zip` | `e909a4f0121688b0590ccaec10b0276171bc24619cf7eb562ba441268c01e121` |
| `CIR_data.zip` | `f4edcb2852b83a86d1b58f29fb862d9124a227cfc90f9e1ae17d2c97585264e6` |
| `iso.zip` | `568667bfc8027095537d642297b319c872d00016b868143c666f90d5931d9f7b` |

官方 `pir_data_constant` / `cir_data_constant` 只接受 25 个代数常数，原始 archive
中不存在上述 π 型 token。官方读取逻辑也确认了本文件上方更正后的 cursor 语义；
因此解析器不是本次 758 个字段异常的根因。

#### P2：compound CIR record 缺少机器可读的组合语义

当前 compound record 主要依赖 Miller–Love 标签和两个 component 数组，未显式
记录它属于哪一种构造：

```rust
enum CompoundSemantics {
    ConjugateRealification, // χ = 2 Re χ₁
    DistinctComponentSum,   // χ = χ₁ + χ₂
    StarInduced,
}
```

实际数据中第二个 component 有时等于第一个、有时是其共轭、有时受 source/Hall
约定污染，仅凭数值关系不能稳健分类。Rustb 目前对重复 constituent 标签采用
`2 Re(first)`，修复了 SG199 `P2P2` 和 SG220 `P1P1/P2P2/P3P3` 等已由显式
对称模型验证的记录；但这是 consumer-side 规则，不等价于上游数据已经自描述。

生成数据最终应显式保存 compound semantics、constituent identity 和 selected-arm
映射，调用者不应再从字符串标签或 component 数值猜测物理含义。

#### P2：full-star 与 selected-arm/little-group 数据共用一个 record

同一个 `IrrepRecord` 同时暴露 `dim`、PIR `characters()`、
`scalar_little_characters()` 和 `cir_component_chars()`，但这些字段可能分别属于：

- 完整 k-star 的诱导表示；
- 指定 star arm；
- fixed-k little-group 表示。

若直接把 full-star `dim` 或 trace 当作 little-group corep 数据，会得到分数
multiplicity 或无法标记的 `???`。SG76 的 R 点 compound record 曾出现这一类
false negative。长期应拆分或类型化为 `FullStarIrrep`、`SelectedArmIrrep` 和
`LittleGroupIrrep`，至少也要为每个字段记录对应空间和维数。

#### P2：spinor source operation mapping 仍不完整

全 UNI 探针中有 `5753` 次潜在 source-character 映射因
`flat_operation_miss` 进入 magnetic-summary fallback，且这些缺口全部来自 spinor
路径。当前 fallback 能保持 summary 可用，但这说明 spinor source operation、
data-Hall operation 与 summary column 之间仍缺少完整、显式、可验证的双射。

不要把 fallback 成功当作上游 mapping 已正确。应为每个 character column 保存
稳定的 operation identity，并对 source→Hall 映射做全量 bijection gate。

#### P2：其他已知的局部数据异常

- SG219 `L3L3` / grey SG79 P 点仍有 rank/coimage 异常；
- 一些 compound CIR 的第二行带有 star-arm/source-setting convention 污染；
- 生成数组混用五位小数、九位小数和较高精度 `f64`。Rustb 已对白名单中的
  `1/2`、`1/√2`、`√3/2`、`√2`、`√3`、`2√2` 做 consumer-side 精确化，
  但上游最好保存 exact cyclotomic/algebraic 表达式或至少保存来源与误差等级。

#### 验收状态与后续 gate

- 已恢复并固定官方 archive 版本与哈希；官方 Fortran reader 已用于核验 parser
  cursor 和矩阵语义。
- formatter 修复后完整再生成只改变 758 个已知污染字段；静态 π 型扫描为零。
- SG144 A2 最小回归和 literal 往返 gate 已加入。
- 全 `1651` UNI summary gate 结果为 `success=1651`、`failure=0`、
  `kpoints=10390`、`coreps=52793`、`amplified_noise=0`。
- 仍需长期加强 `χ(E)=dim`、`|χ(g)|≤dim`、角色正交性、整数 multiplicity、
  cyclotomic 可识别性，以及 compound constituent/selected-arm 的机器可读语义。
- 独立显式 Hamiltonian/orbit oracle 仍应覆盖 SG161、167、169、170、178、179
  和已知 SG219 case；这些 P2 项与本次 formatter 根因不同，不应混为同一修复。

---

## 方法论：数据优先于算法（CIR-PIR debug 教训）

### 问题回顾

CIR-PIR 测试失败（596 个 compound irrep, 1222 mismatches）。花了大量时间尝试
Bloch 相位修正（`exp(i*2π*k·Δt)`），每改一次 mismatches 反而增多。
用户几个简单问题（"PIR 和 CIR 分别是什么数据？"、"两者的 translation 可以对比吗？"）
引导重新检查 ISOTROPY 源数据格式，才发现 CIR 字符只依赖旋转不依赖 translation。

### 为什么会陷入思维玄幻

| 阶段 | 错误做法 | 应该做的 |
|------|---------|---------|
| 1 | 假设 CIR"应该"和 PIR 操作数相同 | 先检查两者的源数据格式 |
| 2 | 尝试 Bloch 相位公式修正（多次迭代失败） | 先验证假设：CIR 真的需要相位吗？ |
| 3 | 修改越来越复杂（Hall trans, ISO trans, Δt...） | 回退到数据源，理解 CIR 格式 |
| 4 | 每轮 fix 让 mismatches 增加而不是减少 | 如果修不好，说明假设错误，不是公式不对 |

### 核心教训

**规则：当修复让问题变严重时，停止修代码，回到数据源检查假设。**

1. **数据 > 算法**：不知道数据长什么样之前，不要写任何修正公式
2. **简单的问题解决复杂的问题**：用户问"数据分别是什么"直接打破了思维定式
3. **物理/数据约束是第一位的**：CIR 没有 translation 不是 bug，是结构性的——
   CIR 是小群表示（只关心旋转），PIR 是空间群表示（含平移相位）。
4. **失败的修复是信号**：如果改了 3 次 mismatches 还在增加，假设一定错了

### 量化对比

| 指标 | Bloch 相位路线 | 最终方案 |
|------|---------------|---------|
| 修改次数 | 8+ commits, 每次 ~100 行 | 1 commit, ~30 行 |
| mismatches 趋势 | 0→1276→1222（恶化） | 0（直接解决） |
| 修改范围 | Python 脚本 + Rust 测试 | Rust 测试 min(n_ops, cir_ops) |
| 根本认知 | CIR"缺失"数据需要推算 | CIR 本来就只覆盖旋转类型 |

---

## 调试方法论 — 从 spinor Wigner 排查中提炼的经验

### 原则 1：比较 passing vs failing cases，找差异因子

当同一个代码路径**有些 case 通过、有些失败**时，不要猜测通用原因。直接比较一个成功 case 和一个失败 case，问：

> **"这两个 case 之间什么不同，导致一个通过、一个失败？"**

这个差异通常直接揭示根因。

历史实例：SG2 T-point passed after loop fix → SG159 L-point still failed。当时观察到
SG2 的 LG {I, -I} 包含 `(a₀h)²`，而 SG159 的 LG {I, mirror} 不包含。
这个比较成功定位了“失败发生在 square/LG mapping 阶段”，但当时进一步断言
“不是 setting/algorithm 问题”是错误的；2026-06-19 的 UNI663 oracle 已证明，
MSG 嵌入基底和 standalone H Hall/spin-table 基底不一致本身就会制造这种
`square_not_in_spin_table` / `square_outside_little_group` 现象。

### 原则 2：假设驱动的逐层排除

不要在无数可能原因中随机尝试。为每种假设设计一个**最小 oracle test**，一票否决或确认：

1. **列出所有可能假设**（按先验概率排序）
2. **为每个假设设计一个 oracle**（最小代码改动，只输出统计数据）
3. **跑 oracle，看数据**——如果数据否决假设，立刻排除，不再纠结
4. **如果 oracle 确认假设，再设计修复**

实例——NONE=1,007 的排查顺序：
- H1: same-rotation lift 误选 → scan same-rot candidates → OTHER=0 → 否决
- H2: UU* antiunitary square → 6 formulas oracle → 6.5% fix → 不是主因
- H3: H/G gauge mismatch → 当时的 G-gauge oracle → 0% fix → **错误地否决**
- H4: det=-1 improper → det stats → 混合分布 → 否决
- H5: J-insertion → J-oracle on NONE → 61% fix → 确认方向
- H5-global: global J → 88.3%→83.0% → 不能全局替换 → 否决
- H5-per-case: case-level J fallback → old_fail_j_ok=22/945 → 否决

每排除一个假设，就缩小搜索范围。不要跳过 oracle 直接修代码。

### 原则 3：诊断与修复分离

诊断代码（oracle/counter/scan）**不应改变正式分类结果**。先加诊断、跑数据、看统计、确认假设，再设计修复。

- Oracle 只在 `None`/失败分支执行，不改变 `return` 值
- 计数器用 `AtomicUsize`，在 diagnostic test 中读取
- 正式路径保持原样，等 oracle 确认后再改

### 原则 4：per-term → per-case 的层级

在 Wigner sum 中，单个 term 的修复不等于整个 case 的修复：

1. **per-term fix**（只对失败 term 用新公式）：数学上危险，同一个 sum 混用两个 convention
2. **global fix**（所有 term 都用新公式）：如果破坏了更多正常 term → regression
3. **per-case fallback**（先试旧公式，整个 case 失败再试新公式）：唯一理论上干净的 fallback

要区分三者，不能看到 per-term oracle 有 61% fix 就急于做 per-term patch。

### 原则 5：语义正确的计数器命名

计数器名字必须准确反映被计数的**物理/数学含义**，不能有歧义：

- 错误：`central=false` → "raw misses"。正确：`central=false` → "same lift, no central element"
- 错误：`theta2_fixes` → "fixing misses"。正确：关系到 `±u_k` 的 same/Ebar/none 三类
- 错误：`NONE=0` → "0 mismatch"。正确：`NONE=1,007` → "1,007 non-trivial mismatch"

错误命名会误导后续分析方向。本例中 `central=false` 被误读为 "raw failure"，导致构造了大量无用的 sign-flip 修复。

### 原则 6：不要过早下结论说"需要大工程"

Data generation 存 `central_parity` 或 `extended character table` 可能是最终方案，但在确认以下问题之前不应断定：

1. **先确认问题确实来自数据缺失**（而非 algorithm bug or convention mismatch）
2. **先做 oracle 估计大工程的收益**（例如 eta ±1 测试）
3. **先排除更便宜的修复**（runtime inference, convention alignment）

### 原则 7：逐阶段确认，不要把 UNI=0 当成整个管线失败

磁群识别是一条多阶段管线。当 `UNI=0` 时，不要笼统地说"磁群识别失败"。

分别检查每个阶段的输出：
- 普通 SG 是否正确？Hall 是否正确？
- 磁类型是否正确？操作数是否正确？
- unitary/anti-unitary 比例是否正确？
- 失败是否**仅**发生在 DB matching 阶段？

实例：石墨烯 AFM 返回 `SG=191, Hall=485, Type-3, 24 ops (12U+12A)`
但 `UNI=0`。前四个阶段全部正确——问题只在 DB 匹配阶段。
这排除了晶格约定、磁操作生成、FSG/XSG 分类等所有问题。

### 原则 8：容差扫描是分类工具，不是修复方法

对 `symprec` 做跨数量级扫描（1e-3 → 1e-6）：

- 结果随容差变化 → 优先检查数值稳定性
- **结果跨多个数量级完全不变** → 优先检查坐标约定、数据 setting、变换方向、群操作合成

石墨烯 AFM 的 UNI=0 在四个数量级容差下完全不变——这排除了数值问题，
直接指向 convention/transform 错误。

### 原则 9：对矩阵方向使用推导，不使用记忆

看到 `T R T⁻¹` 或 `T⁻¹ R T` 时，不要凭变量名判断。先写清楚：

```
x_new 与 x_old 的关系是什么?
```

然后计算 `C g C⁻¹`。推导结果是最小、最可靠的 oracle。

对于 `x_std = T x + s`，正确的 Seitz 共轭是：
```
R_std = T R T⁻¹
t_std = s - R_std s + T t
```

### 原则 10：检查 helper 的所有调用者

一个错误 helper 可能在某个调用点被另一个错误抵消，却在其他调用点失败。

石墨烯案例：`get_distinct_changed_magnetic_symmetry` 内部使用 `T⁻¹ R T`（错误），
但 `get_reference_space_group` 传入的 `tmat` 本身也取反了（把 P 当 T，实际应为 P⁻¹）。
两个错误在参考 setting 路径上互相抵消，但在数据库 correction transformation 路径
（这里 T 已经正确）暴露出来。

单独分析第一个调用点会错误地认为公式"等效"。

### 原则 11：高对称测试不足以验证线性代数约定

立方晶格的变换矩阵通常是单位矩阵、置换矩阵、正交对称矩阵或自逆矩阵
（`T = T⁻¹`）。这些测试全部通过**不能**证明矩阵方向正确。

应当至少包含：
- 非正交六方/单斜晶格（`T != T⁻¹`）
- 非零 origin shift
- 非平凡 correction transformation
- 具体断言 UNI/BNS 值（而非只检查 `> 0`）

石墨烯六方 AFM 是第一个暴露这些问题的非正交 oracle。

---

## 错题集 — 核心教训

### Bug 1: loop domain ≠ character domain
Wigner 求和对象是 **little co-group**（旋转），不是 full little group（Seitz 变体）。loop 遍历了 4 个 Seitz 但 spin table 只有 2 个 co-group 条目 → W=-0.5。

### Bug 2: rotation matching 不能跨基底
rotation 只在同一基底下不变。`x_new = T x_old + s` 时 `R_new = T R T⁻¹`。

### Bug 3: 不依赖数组顺序隐含语义
Grey 群的 a₀ 必须是纯 θ (R=I)。取 `antiunitary[0]` 可能取到 θ·g。

### Bug 6: Θ²=Ē 和 SU(2) central sign
`central = !spatial_central` 的推导：Θ²=Ē=-I（自旋 1/2），spatial=EBAR 时 Ē²=I→SAME→central=false。
**但仍需验证**：`su2_same_up_to_sign(U_b², U_{h²})` 跨 G/H 两套 SU(2) gauge 比较可能出错。

### 已排除的假设
spin 数据库不完整 ❌ | Pauli SU(2) 合成 ❌ | same-rotation lift 误选 ❌ | 
UU* 公式 ❌ | det=-1 improper ❌ | global J-insertion ❌（regression）

### 当前主要问题
Spinor Wigner 的 `square_not_in_spin_table` 主线已经在 2026-07-03 清零。
当前 irrep 方向的主要问题是产品化入口：把 H-irrep/k-point/corep/isotropy
这些已经分散可用的数据组织成稳定的 `magnetic_summary` API。

## Architecture overview

cryspglib has two major subsystems:

**1. spglib port** — space group identification from crystal structures.
`Crystal::new(lat, positions, types)` → `.analyze()` → `.dataset()` → `SpaceGroup`.
Also supports magnetic space group identification (1,651 UNI types) via `.with_magnetic()`.

**2. irrep module** — irreducible representation data for all 230 space groups.
`irreps_of(sg_number)` → `IrrepRecord` (labels, characters, matrices, isotropy subgroups, magnetic corepresentations).
100% data coverage: 4,777 irreps with character tables (~50k values), full matrices (~580k values),
15,239 non-magnetic + 16,721 magnetic isotropy subgroups.
100% of characters are in spglib Hall order.

---

## Build & Test Commands

```bash
cd /home/liuyichen/TB_rs

cargo build --package cryspglib
cargo test  --package cryspglib
cargo test  --package cryspglib <test_name>
cargo check --package cryspglib
```

Key diagnostic test: `irrep::corep::tests::diagnose_wigner_sources -- --nocapture`

Enable verbose Wigner diagnostics (eprintln! output):

```bash
cargo test --package cryspglib --features debug-corep -- <test_name> -- --nocapture
```

### Data regeneration pipeline

`hall_operations.json` is a committed static artifact (does not need regeneration).

```bash
# Regenerate generated_data.rs from ISOTROPY source data:
python3 scripts/generate_irrep_data.py
# Full pipeline shell:
bash scripts/regenerate_all.sh
```

### Diagnostic validation scripts

| Script | Purpose |
|--------|---------|
| `validate_cir_pir.py` | Standalone CIR→PIR validation — checks PIR = Σ CIR * exp(i2πk·t) per Hall op |
| `check_iso_vs_spglib.py` | Compare ISOTROPY primitive vs spglib Hall conventional translations |
| `check_phase_correction.py` | Analyze Bloch phase corrections when mapping ISOTROPY→spglib Hall order |
| `debug_cir_pir_sg9.py` | Single-case debug: SG9 CIR/PIR mapping details |
| `test_su2_closure.py` | Pauli SU(2) composition closure test |
| `test_spinor_wigner_formula.py` | Spinor Wigner formula standalone test |

---

## Key types

| Type | Location | Description |
|------|----------|-------------|
| `Crystal` | `api.rs` | Entry point: lattice + positions + types + optional magnetic moments |
| `SymmetryAnalysis` | `api.rs` | Builder for symmetry analysis (`.symprec()`, `.dataset()`, `.magnetic_dataset()`) |
| `SymmetryOps` | `api.rs` | Ordered set of `{R\|t}` + time_reversal, with `from_database(hall_number)` |
| `SymmetryOp` | `api.rs` | Single `{R\|t}` with rotation, translation, time_reversal |
| `SpaceGroup` | `lib.rs` | SG number, Hall number, ops, Wyckoff positions, standard cell |
| `MagneticSymmetry` | `lib.rs` | MSG + symmetry ops combined (implements `Display`) |
| `MagneticSpaceGroupType` | `lib.rs` | MSG type lookup: `.from_uni()`, `.classify()` |
| `SpaceGroupType` | `lib.rs` | SG type lookup: `.from_hall()` |
| `IrrepRecord` | `irrep/types.rs` | Irrep: labels, dim, k-vector, characters, matrices, subgroups, corepresentations |
| `KVector` | `irrep/types.rs` | Rational reciprocal-space vector with three numerators and one denominator |
| `WignerGroupContext` | `irrep/wigner.rs` | Unitary/magnetic operation slices and antiunitary representative for Wigner classification |
| `SpinorWignerInput` | `irrep/wigner.rs` | Spinor character table, operation indices and rational k-vector |
| `SpinLiftContext` | `irrep/wigner.rs` | H and G spin ops for Wigner test |
| `SeitzOp` | `irrep/wigner.rs` | `{R\|t}` with optional time reversal |
| `CorepType` | `irrep/corep.rs` | A/B/C/Unsupported |

---

## Module structure

### spglib port subsystem

| Module | Role |
|--------|------|
| `api.rs` | `Crystal` (entry point), `SymmetryAnalysis` (builder), `SymmetryOps`, `SymmetryOp` |
| `lib.rs` | `SpaceGroup`, `SpaceGroupType`, `MagneticSymmetry`, `MagneticSpaceGroupType`, `SymError` |
| `cell.rs` | `Cell` (lattice + positions + types + optional tensors) |
| `symmetry.rs` | `Symmetry` (raw symmetry operations, N×rot+trans arrays) |
| `spacegroup.rs` | `Spacegroup`, `search_spacegroup*` |
| `spg_database.rs` | `get_spacegroup_operations`, `get_spacegroup_type` |
| `magnetic_spacegroup.rs` | MSG identification: `identify_magnetic_space_group_type` |
| `msg_database.rs` | `get_magnetic_spacegroup_type` (1,651 UNI entries) |
| `msg_database_gen.rs` | Auto-generated MSG database (shipped as source) |
| `spin.rs` | Spin-polarized symmetry: `get_idealized_cell`, `collect_pure_translations_from_magnetic_symmetry` |
| `pointgroup.rs` | `get_pointgroup`, `get_transformation_matrix` |
| `primitive.rs` | `get_primitive`, `get_primitive_symmetry` |
| `delaunay.rs` | Delaunay lattice reduction |
| `niggli.rs` | Niggli lattice reduction |
| `hall_symbol.rs` | Hall symbol parsing/conversion |
| `kpoint.rs` | k-point grid generation, irreducible mesh |
| `kgrid.rs` | Grid address utilities |
| `determination.rs` | Space group determination pipeline |
| `refinement.rs` | Cell refinement / idealization |
| `overlap.rs` | Atom overlap detection |
| `parser.rs` | POSCAR parser |
| `site_symmetry.rs`, `sitesym_database.rs` | Site symmetry + database |
| `arithmetic.rs` | Arithmetic crystal class symbols |
| `mathfunc.rs` | `Mat3`, `Mat3I`, `Vec3`, matrix/vector operations |
| `debug.rs` | Diagnostic print helpers |

### irrep subsystem

| Module | Role |
|--------|------|
| `irrep/mod.rs` | Module docs, re-exports, coverage summary (4,777 irreps, 100% coverage) |
| `irrep/types.rs` | `IrrepRecord`, `IsotropyRecord`, `MagneticIsotropyRecord` |
| `irrep/query.rs` | `irreps_of()`, `kpoints_of()`, `format_character_table()` |
| `irrep/corep.rs` | Co-representation: `compute_coreps()`, `CorepType`, diagnostic tests (30+ tests) |
| `irrep/wigner.rs` | Wigner test: Seitz composition, SU(2) composition, spinor classification |
| `irrep/bridge.rs` | `impl SpaceGroup` — bridge APIs linking spglib port → irrep |
| `irrep/generated_data.rs` | Auto-generated static arrays (~753k lines). `include!()`-d into `types.rs` |
| `irrep/wigner_extra.rs` | Pre-computed antiunitary character path. `include!()`-d into `wigner.rs` |
| `irrep/preamble.rs` | Generated data prelude |
| `irrep/{triclinic,monoclinic,orthorhombic,tetragonal,trigonal,hexagonal,cubic}.rs` | Per-crystal-system irrep data (`include!()`-d into `generated_data.rs`) |

---

## Test suite（2026-08-27 release 基线）

Core irrep diagnostics pass as of 2026-07-03.  The full spinor Wigner sweep reports
`spinor_complex_ok = 21216` with no `spinor_complex_fail`.

| Binary / Location | 当前结果 | Description |
|-------|-------|-------------|
| `src/lib.rs`（全部 unit modules） | `291 passed / 4 ignored` | Wigner、BCS、API、输入契约、setting 与磁群回归 |
| `tests/irrep_validation.rs` | `34 passed` | Full-sweep validation: every SG has irreps, dimensions match, labels well-formed, k-vectors positive, generated characters contain no amplified exponent noise |
| `tests/magnetic_integration.rs` | `17 passed` | Magnetic structure analysis and Result/error contracts end-to-end |
| `tests/magnetic_symmetry_coverage.rs` | `6 passed` | 1651 UNI / 4479 setting group algebra, round-trip and ambiguity policy |
| `tests/{cof3,crps4,la2nio4,bcs_corep_validation}.rs` | `5 passed` | Reference material cases |
| doc-tests | `26 passed` | Public Rust API examples compile and run |
| **常规 release 总计** | **`353 tests + 26 doctests passed; 4 ignored`** | 另行执行的 1651 summary audit 也通过，`amplified_noise=0` |

Key diagnostic tests (most useful for Wigner debugging):

```bash
# Primary: full-sweep Wigner failure diagnosis (~217s)
cargo test --package cryspglib diagnose_wigner_sources -- --nocapture

# Setting transform oracle (identity/signed_perm/unimodular/none x ok/other/square_not_in_spin)
cargo test --package cryspglib phase1_setting_transform_oracle -- --nocapture

# Per-term and per-case failure analysis
cargo test --package cryspglib diagnose_spinor_wigner_per_term -- --nocapture
cargo test --package cryspglib diagnose_nonquantized_per_term -- --nocapture
cargo test --package cryspglib diagnose_none_examples -- --nocapture

# CIR-PIR consistency and data integrity
cargo test --package cryspglib test_cir_pir_cross_validation -- --nocapture

# Magnetic entry point diagnostics
cargo test --package cryspglib diagnose_magnetic_entry_hall_anomalies -- --nocapture
cargo test --package cryspglib diagnose_spglib_standard_setting_transform -- --nocapture

# k-convention oracle (3611 irreps)
cargo test --package cryspglib diagnose_spin_lg_k_convention -- --nocapture
```

Full validation sweep (integration tests, ~1 min):

```bash
cargo test --package cryspglib --tests
```

---

## 磁空间群识别的坐标变换约定

完整故障复盘、数学推导、错误互相掩盖机制和调试方法见：
`docs/magnetic-spacegroup-basis-transform-postmortem.md`。

### 标准设置变换

`magnetic_spacegroup.rs` 使用

```text
x_std = (T, s) x
```

因此 Seitz 操作必须按下面的方向变换：

```text
R_std = T R T⁻¹
t_std = s - R_std s + T t
```

`get_reference_space_group` 返回的 `tmat` 是参考空间群
`bravais_lattice` 的逆矩阵。不要把共轭方向改成 `T⁻¹ R T`，
也不要直接把 `bravais_lattice` 当作 `tmat`；立方体系可能碰巧通过，
但六方等非对称基底会无法匹配正确的 UNI。

**回归案例——石墨烯 AFM（z 方向反铁磁）**：

```text
SG=191 (P6/mmm), Hall=485, type=BlackWhite
24 ops (12 unitary + 12 anti-unitary)
UNI=1466, BNS=191.236
```

测试命令：

```bash
cargo test --package cryspglib --test magnetic_integration test_graphene_afm_z
```

该测试覆盖 `symprec=1e-3..1e-6`，结果必须保持一致。此前的 `UNI=0`
不是二维体系的数据库限制，而是标准设置变换方向移植错误。

**双层石墨烯 oracle（z=0.51 / z=0.49）**：

打破水平镜面对称后，空间群从 P6/mmm (#191) 降到 P-3m1 (#164)：

| 磁构型 | SG | Hall | ops | UNI | BNS |
|--------|-----|------|-----|-----|-----|
| 非磁 | 164 | 456 | 12 | 0 | - |
| FM (都↑) | 164 | 456 | 12 | 1319 | 164.89 |
| AFM (一↑一↓) | 164 | 456 | 12 | 1318 | 164.88 |

这是第二个非立方磁群 oracle——对称操作数从 24 降到 12，
且涉及三方晶系（介于立方和六方之间），对 setting transformation 的敏感度
高于立方案例但低于平面六方案例。

---

## 错误处理约定：Result, not Option

spglib port 的主要公共 API 已全部从 `Option<T>` 迁移到 `Result<T, SymError>`。
`SymError` 是 unit-only enum（无数据字段），每个 variant 名直接指向失败位置。

### 为什么不用带数据的 Error

`SymError` 的 variant 已经精确到单个函数甚至函数内的单个失败点。
携带额外数据不会增加定位精度，反而破坏 `Copy` 和简单模式匹配。

### 改造的函数

| 管线 | 函数 | 错误 variant |
|------|------|-------------|
| 磁群 | `operations_with_site_tensors` | `MagneticOpGenerationFailed`, `MagneticPrimitiveLatticeFailed` |
| 磁群 | `identify_with_parent_hall` | `MagneticReferenceGroupFailed`, `MagneticFallbackReferenceFailed`, `MagneticUniMatchFailed` |
| 磁群 | `magnetic_dataset()` | 传播上游错误 |
| 非磁 | `get_primitive` | `CellStandardizationFailed` |
| 非磁 | `get_operation` | `SymmetryOperationSearchFailed` |
| 非磁 | `search_spacegroup` | `SpacegroupSearchFailed` |
| 非磁 | `determine_all` | `SpacegroupSearchFailed` |
| API | `Crystal::from_poscar` | `InvalidInput` |
| API | `SymmetryOps::from_magnetic_database` | `SpacegroupSearchFailed` |
| API | `SymmetryOps::from_sg` | `SpacegroupSearchFailed` |
| API | `find_hall_number` / `find_first_hall_for_uni` | `SpacegroupSearchFailed` |

### 错误传播模式

- 使用 `?` 直接传播同类型错误
- `Option` → `Result` 转换用 `.ok_or(SymError::Variant)?`
- `Result` → `Option` 转换（仅在临时 backward compat 中）用 `.ok()`
- `MagneticUniMatchFailed` 与其他磁群识别错误一样由 `magnetic_dataset()` 原样
  传播；不得通过 FSG/XSG 群阶猜测磁类型并伪装为成功结果。

### 晶格矩阵约定

`lattice[cart][vec]` —— **行=笛卡尔分量，列=晶格矢量**。
详见 `mathfunc.rs` 模块文档。六方晶格不对称——约定错误将导致
空间群识别错误（如 graphene #191 → #10）。

---


## 审计教训：Phase 1 假修复复盘

> 以下是从"136/136 fixed 实为 0/136"事件中提炼的核心教训。完整记录见 git log。

### 六类错误

| # | 错误 | 预防 |
|---|------|------|
| 1 | **失败阶段转移冒充修复**：`!matches!(result, SquareNotInSpinTable)` 把症状转移当治愈 | 判据必须是**总失败数下降** |
| 2 | **总失败数不变不追问**：679 自始至终没变，却自我合理化 | 总数不变时先追问为什么 |
| 3 | **声称完成但只做了 1/5**：计划要求 origin/Seitz双射/多解检测，实际只有 rotation multiset | 逐条对照计划清单打勾 |
| 4 | **挑选有利指标**：只看 W_FAIL 和 OLD_PATH_FAIL 下降 | 总失败数是唯一不可作弊的指标 |
| 5 | **用 sed 批量改代码**：违反规则 2 | 只用 Edit 逐处修改 |
| 6 | **不区分诊断/生产路径**：诊断传了 setting_xf，生产传 None，却声称"已接入" | 两者必须同时验证 |

### 核心原则

- **判据必须是总失败数下降，不能是中间阶段转移**
- **宣布完成前逐条对照计划清单**
- **总数不变时必须追问为什么，不能自我合理化**
- **只做诊断 oracles，不急于改正式分类结果**

---

## 排查方法论精华

### 核心流程

```text
诊断 oracle → 全量统计 → 确认 convention → 找到边界 case → 修正 → 验证 → 下一个
```

### 关键实例

**k-convention 排查**：
1. oracle: 对 3611 个 spinor irrep 比较 Rk vs R⁻ᵀk
2. 结果: reciprocal_exact=3611 → 确认 R⁻ᵀ 是正确 convention
3. 边界: centered cell 还需 pure translation phase check
4. 数据 bug: 1/3,1/6 k 点被错误 rationalize → 修复 parse_spinor_data.py
5. 效果: 679→347 (-332)

**Setting transform 排查**：
1. oracle: UNI663 比较 ops_from_msg vs ops_from_hall rotation → C2z≠C2y
2. 假修复: 阶段转移冒充修复 (136/136 实为 0/136)
3. 真修复: 完整 (T,s) 求解 (Gaussian elimination + modulo-1)
4. 发现: 48 signed-perm 不够 → 泛化为 rational Mat3
5. 入口 bug: UNI187 unitary 含 mirror_y 但被识别为 SG1 → Hall 选择错误

### 铁律

1. **数据 > 算法**：先检查源数据再调公式
2. **判据 = 总失败数下降**：不接受阶段转移或其他中间指标
3. **全量 oracle，不靠单例**：`diagnose_spin_lg_k_convention` 遍历全部 3611 个 irrep
4. **100% 通过后仍追问边界**：reciprocal_exact=3611 后仍发现 centered cell 问题
5. **入口数据出错时停止下游修复**：UNI187 的 unitary 操作本身错了
6. **穷举 convention，让数据说话**：不确定 Rk vs R⁻ᵀk 时两个都算，比较 exact match

---

## 当前 magnetic corep 能力与剩余边界（2026-08-30）

普通 230 个空间群的 typed scalar/spinor 角色数据、operation pairing，以及 Rustb
能带标记所需的 fixed-k magnetic corepresentation 数据链已经可用。此前三项主要
阻断已经关闭：

1. **Magnetic compound corep**：`compound_complex_corepresentations()` 在 CIR
   constituent 层做 Wigner 分类，显式建立 antiunitary constituent orbit；固定
   constituent、内部交换和外部 partner 都保持 source provenance。Rustb 只消费该
   plural API 并按 provenance 合并，不解析 compound label。相关提交：`efc7e24`、
   `9e956c4`、Rustb `da4df70`。
2. **一般 Type-C partner**：不再把所有 Type-C 简化为 `2 Re(chi)`。scalar 与
   spinor 均计算 `chi(a0^-1 h a0)^*`，exact Seitz reduction 返回 canonical operation
   与 lattice shift，再施加 Bloch phase。spinor 另外比较 SU(2) lift 的
   Same/EBar central sign。相关提交：`0851ca9`、`83ca8e8`。
3. **Spinor full-Seitz / setting transport**：stored spin representative 到实际
   Hall operation 的 translation difference 在真实 Hall translation lattice 中验证，
   character 乘精确方向的 Bloch phase；parent/H setting 不再局限于 signed
   permutation，而从整组 pinned SU(2) adjoint actions 解一个全局 frame。相关提交：
   `1de71e1`、`83ca8e8`。

严格 `magnetic_irrep_summary_by_uni()` 现在是 operation-aware 的复数接口：
`MagneticCorepSummary.characters` 使用 `Vec<Option<Complex64>>`。`Some(z)` 表示
正式计算值，物理零也必须写作 `Some(0)`；只有尚无 canonical intertwiner 的
Type-A 反幺正列才是 `None`。ordinary、spinor 与 compound plural 结果均在
cryspglib 内按 source provenance 合并，consumer 不再重建 source row。

cryspglib 与 Rustb 的 release-only 全量 gate 已把该接口作为唯一生产路径。全
1,651 UNI 的当前结果是：

- `summaries=1651`；
- `points=10390`；
- `coreps=52793`；
- strict-summary failure `0`。

因此 compound、一般 Type-C、centering/lattice phase、一般 SU(2) frame 以及
genuinely complex character 已不再需要 Rustb 的 partial-summary recovery，也不应
再因 legacy `f64` 投影产生上游数据型 `???`。

仍需明确保留的边界：

- **Legacy real-valued surface**：`Corepresentation.characters: Vec<f64>` 等旧入口
  仅保留为 deprecated compatibility API；它们会丢失 genuinely complex unitary
  character，不得再用于 summary 或 Rustb 能带标记。底层调试可直接调用
  `complex_corepresentation()` / `compound_complex_corepresentations()`，上层分析
  应优先使用严格复数 summary。
- **Type-A antiunitary trace**：一般高维或 spinor Type-A 尚无 canonical intertwiner
  matrix，返回 `TypeAAntiunitaryPending`。placeholder zero 不是计算值；Wigner type 与
  unitary row 足以完成当前 Rustb formal fitting。
- **表示范围**：当前 summary 是 fixed-k magnetic-little-group corep，不是 full
  k-star corepresentation；spinor isotropy subgroup 仍标为
  `SpinorNoIsotropyData`。非 primitive cell 的 band unfolding 也仍由 consumer
  显式拒绝。

## Python 数据工具保留策略（2026-08-30）

仓库现有 29 个 Python 文件、约 21,698 行。它们**不是 Rust build/runtime
依赖**：仓库没有 `build.rs` 调用 Python，Rust API 只读取已提交的
`src/irrep/generated_data.rs`；`Cargo.toml` 也把 `scripts/` 排除在发布 crate 之外。

Python 仍应保留在源码仓库中，因为它承担科学数据的可再生性和独立审计：

- 直接生成链：`generate_irrep_data.py`、`parse_spinor_data.py`、
  `direction_map.py`、`iso_irrep_data_hall.py` 及固定 data/manifest；
- authority/再生链：`iso_irrep_exact.py`、`derive_iso_irrep_data_hall.py`、
  `freeze_iso_irrep_data_hall.py`、spglib magnetic provenance extractor/loader、
  `spinor_exact.py`；
- 对应 `test_*.py` 是 generator/parser/operation-order/phase convention 的科学
  regression oracle，不能只因 Rust tests 通过就删除。

可后续精简的是一次性诊断脚本，例如 `debug_*`、`check_*` 和已经被正式测试覆盖的
旧 validation probe。删除前必须先证明其中没有唯一 oracle，并把仍有价值的 witness
迁入正式 `test_*.py`。优先做目录分层（generation / audit / legacy_debug）和共享解析
库去重；不要为了减少 Python 行数而立即用 Rust 重写稳定的离线生成链，也不要继续
扩展与角色计算无关的通用 provenance 框架。

### 分导任务 9 第四十二轮（2026-09-22）：w 行的引擎计算，`P1` 子群族先落地

新增 `src/irrep/subduction_star_decompose.rs::line_trivial_content_with_embedding`：
把冻结的直线源变成频率——星臂（方向的 contragredient 像，按值去重）逐个做
`a/4 ∈ L_H^*` 的折叠判定（`LINE_PARAMETER = 1/4`，即官方程序对参数化 k 域的自由参数
约定；在 46 条 `P1` 子群记录的 300 个 pinned 行上，`1/4`、`3/4` 全中，`1/2`、`1`
错 111 行，其余十二分之一全错），通过的臂再对子群点群（每个旋转取一个代表）做
`W'` 的恒等投影，累加即频率；结果非整数即 fail closed
（`NonIntegralLineFrequency`），源与母群不匹配另有 `LineSourceMismatch`。

实测：`tests/w_line_frequency.rs::p1_child_records_match_every_pinned_row` 对
SG 196/202/203/209/210/216/219/225/226/227/228 的 46 条 `P1` 子群记录**逐行算出
300/300 与 pinned 相同**——即这批 w 行已由引擎计算（其余仍是 live oracle 5,756/5,756
验证）。**未闭合**：`asymmetric_dt3_dt4_records_are_pinned` 被 `#[ignore]`，非 `P1`
子群的折叠/权重约定还没复现（SG 202 ordinal 10422 `DT3`：模型折出 4 个臂、恒等投影
抵消为 0，pinned 是 2）。审计的 w 轨道因此在 `--require-w-complete` 下仍报未闭合，
覆盖声明保持两条轨道；下一步优先试子群自身帧（`fold_wave_vector(embedding.transform(),
k)`）与投影代表元的取法。

### 分导任务 9 收官（第 113 轮）：全表 100% 由引擎计算，任务 9 完成

`--require-complete --require-w-complete` 全表运行退出 0、判词 `VERDICT complete
scope=global`：普通恒等分导 94,271/94,271、其它波矢 w 行 **5,756/5,756** 全部由引擎
算出并与 pinned（live oracle 逐行复核 0 不匹配）相同，`hard_failures=0`、
`accounting_violations=0`、`census_mismatch=0`，覆盖声明收敛为**一条轨道**。

关键实现（本轮之前几轮陆续落地，均可复算）：

* `line_trivial_content_with_embedding`：手写臂求和路线，先覆盖 3,916 行；
* `FoldedPoint/FoldedStar::from_parts`、`line_folded_stars`（臂按 `LINE_PARAMETER = 1/4`
  折叠、按 q 归组）、`LineArmSource`（`q_block_dimension`/`q_block_character`）、
  `ArmCharacterSource` 枚举与 `build_block` 泛化；
* `line_trivial_content_via_blocks`：把直线表示送进既有 `build_block` 折叠/解块流程，
  审计把它作为第二条独立路线，任一路线命中 pinned 即计为 computed——它一次修好
  SG 196 的 30 行（106/106）并把全表从 3,916 推到 5,744；
* 最后 12 行：8 行是 SG 209 `DT3`/`DT4` 共轭对标签颠倒（`cogroup_pair_route` 的
  `C4` 约定在该母群未被独立确认；SG 210 的两行频率恒相等故不可观测），交换生成器
  与两张守护表后 SG 209 变为 372/372，另外 4 行（child #212/#213）随之同时解决。
* 期间被测量否定的读法（勿再试）：`H/T_H` Mackey 平均、逐臂特征标平均、严格不变
  权重、臂与负臂同一化、first-hit 参数 `1/g`、X 点兼容 irrep 代替、记录格 vs 接受格，
  每一条都在 `docs/task9-remaining-work.md` 留有见证数字。

### 分导任务 9 复核修复（P1×2，第 114 轮）

复核发现两处 P1，已修复：

1. **审计按期望答案选算法会漏报**：w 行原先是「先手写路线、不匹配再试块路线」，
   把 ordinal 10030 DT1 的 pinned 值 1 改成 2 后两个门禁仍退出 0。现在审计**只走
   块路线**，不一致即 `self.mismatch`（新增 `w_frequency_mismatch` 与摘要中的
   `mismatched=`），出错即报错；公开的 `line_trivial_content_with_embedding` 改为
   委托块路线，旧手写求和已删除（它曾有 1,599 个错值与 241 个错误返回）。
   `tests/w_line_frequency.rs` 新增三项：块路线逐行钉住 SG 196 全部 106 行、
   公开入口与块路线等价、以及两条负例。
2. **块路线缺输入上下文校验**：现在拒绝「别的母群的源表」（`LineSourceMismatch`）、
   「别的记录序号的 embedding」（同 variant，按 `embedding.ordinal()` 判定）与
   「奇异子群基」（新 variant `SingularLineBasis`），各有负例测试。
3. **覆盖声明**：w 行数值 5,756/5,756 完整，但结果分两类——351,547 个完整分解 +
   14,713 个仅恒等重数（折叠子群星缺随包离散数据；例：ordinal 13345 `W1` 的完整
   分解仍 `MissingChildStarData`）。
4. **SG 209 的 `DT3`/`DT4` 交换是「由 pinned 频率校准的标签约定」**，生成器与两张
   测试表同步交换不构成独立验证；归档 CIR 的独立检查覆盖的是 SG 202。

### 分导任务 9 复核修复（第 116 轮：`dfbb947` 复核）

1. **频率不一致接入硬失败**：`w_frequency_mismatch` 之前只打印、只计入
   `mismatched=`，没有进 `hard_failures()`；把 ordinal 10030 `DT1` 的 pinned 值
   1→2 后，默认运行与 `--require-complete` 仍退出 0。现在该计数并入硬失败，
   注入复现在四种开关组合下都是 exit 1 + `VERDICT inconsistent`，并有永久负例
   `a_frequency_mismatch_fails_under_every_flag_combination`。
2. **直线入口复用既有上下文校验**：`validate_subduction_context` 拆出不含 probe
   检查的 `pub(crate) validate_record_and_embedding`，`line_trivial_content_via_blocks`
   直接调用它；改 basis/origin/子群号/凝聚 irrep 都返回 `StaleIsotropyRecord`
   （与普通分导入口一致），embedding 不匹配返回 `EmbeddingContextMismatch`。
   上一轮按症状打的补丁与 `SingularLineBasis` 已删除；直线侧只保留
   「冻结源表必须属于该母群」这一条特有检查（`LineSourceMismatch`）。

### 分导任务 9 复核修复（第 118 轮：`5108f42` 复核）

1. **w 引擎 `Err` 计入硬失败**：错误分支原先只 `bump_error()`，注入一次
   `RationalOverflow`（ordinal 10030 `DT1`）后默认运行 exit 0/clean、
   `--require-complete` exit 0/complete、只有 w 门禁给出 2。现在 `Counts` 新增
   `w_engine_error` 并并入 `hard_failures()`，`other_wave_vector:` 与 `w_scope:`
   两行都打印 `engine_errors=`；同一次注入在四种开关组合下全部 exit 1 +
   `hard_failures=1` + `VERDICT inconsistent`
   （`w_scope: rows=3 computed=2 uncomputed=1 engine_errors=1`）。永久负例
   `a_w_engine_error_fails_under_every_flag_combination`。
2. **恢复被 `c93754e` 删掉的十个审计回归**：example 测试模块回到 13 项（2→13），
   覆盖冻结清点与恒等行唯一性、每个 SG 的唯一恒等子记录、两个逐记录验收见证
   （ordinal 12400/13345）、逐 probe 行完整性与 embedding 不可用路径、计数守恒与
   丢行、重复/冲突频率、Frobenius 复成分维数、几何零项不得跳过计算。适配点：
   13345 的五个 `W` probe 现在是 `identity_only`（`missing=0`，
   `--require-complete` 退出 0）；ordinal 0（旧的不可嵌入见证）现在嵌入成功并完整
   分解 8/8，因此 embedding 不可用路径改由合成 tally 驱动，并另钉住 ordinal 0 的
   修复后状态；`exit_code` 调用改为双门禁签名。
3. **覆盖声明归位**：351,547 + 14,713 的拆分从 w 行搬到
   `docs/subduction-audit.md` 的 366,260-probe 普通行；该文档末尾「引擎还不能算
   w 行」的过期段落改为完成态（`engine_errors` 与 `mismatched` 一并计入硬失败）。

下一里程碑：14,713 条「仅恒等重数」要升级为完整分解。首步逐星清点已完成，见
`docs/subduction-gap-census.md`：2,761 记录、128 个子群号，共 21,136 个缺失星，
归并为 796 个子群/k-star 组合、989 个子群/k-star/setting 组合。全部缺失块非 Γ
且维数为正，不能因恒等重数为零而跳过完整分解。建议先处理子群 #1 的 1,835 个
probe（40 个 k-star）；生产分导算法本轮未变。SG 209 的 `DT3`/`DT4` 仍标注为
「由 pinned 频率校准的标签约定」，独立来源验证另行推进。

### 分导 R1/R2（2026-09-22）：目标表示可由「子群操作 + 精确 q」现场构造

R0 建立独立完整分解门禁后，R1/R2 解除"目标成分必须绑定静态 `IrrepRecord`"的限制：

1. **R1 目标来源分两类**（`src/irrep/subduction.rs`、`src/irrep/subduction_star_decompose.rs`）：
   `FullStarTarget` / `SubductionTarget` 的 `ml` / `bc` / `row_ml` / `irnumber` 改为
   `Option`；存储成分保留冻结 CIR 身份，构造成分用
   `SubductionComponent::Constructed { q, index }` 作为稳定身份，`None` 表示"没有
   来源标签"，不借 Γ 标签、不填假 CIR 号。`ChildComponent` 的字符来源改为
   `ComponentCharacters::{Stored(CharacterRow), Constructed(ConstructedLittleRep)}`，
   `ChildStarEvaluator` 同步增加构造分支（三个 variant 都装箱以避免
   `large_enum_variant`）。`complex_targets` 之外的一切解块、Gram、整数重数、维数和
   与逐操作重建流程不变，既有 fixture 的身份/标签/重数/重建全部保持。
2. **R2 先闭合子群 #1**：`constructed_child_components_at` 为 child #1 提供
   `ConstructedLittleRep::BlochPhase { q }`，即 `D_q(T_L) = exp(+2 pi i q.L)`
   （与存储行共用的 `bloch_phase` 同号）；q 先按子群倒格约化，故 q 与 q+G 同一身份。
   查找顺序是"存储行优先、无行才构造"，其它子群一律返回空、继续报
   `MissingChildStarData`。审计侧 `inspect_result` 接受无 CIR 号的构造成分，但仍要求
   它是 `Constructed` 且通过同一套重建检查。
3. **实测**：`--require-complete` 全表 `full_success=353382 identity_only=12878`
   （原 `351547 + 14713`，−1,835 正好是 #1 的缺口），恒等正项 94,271、Γ Frobenius
   1,895、w 行 5,756 全部不变、`hard_failures=0`、exit 0；同一审计上重跑缺口清点得
   `records=2380 probes=12878 missing_stars=17144 replay_errors=0`，子群 128 → 127。
   永久回归：`child_p1_records_decompose_every_scalar_probe_without_pinned_data`
   （1,125 条 child-#1 记录 / 20,099 个 probe / 3,992 个构造目标全部重建通过）、
   `constructed_bloch_phase_pins_the_positive_sign_convention`（非零平移与
   q→−q 共轭）、`constructed_targets_have_their_own_identity_and_no_borrowed_labels`
   与 `a_child_p1_gap_is_answered_by_constructed_targets`（ordinal 1045，
   SG 45 `S1S2`/C1，probe `W1W1`）。
4. **审查口径修正（用户）**：R2 修的是"目标数据不可达"，不能据此声称"肯定不是
   嵌入或折叠问题"——恒等重数吻合排除不了全部 setting 与标签错误，所以补齐后仍保留
   完整重建检查；表示计算与 CDML/BC 命名分离，构造成分在取得有来源的映射前不冒充
   标准标签。R3 的重点明确为"小群/因子系统分类与来源验证"，R4 实现可复用的目标
   生成；spgrep 只作离线对照（符号/坐标/代表元/相位需先统一），不在运行时依赖。

### 分导 R1/R2 复核修复（第 121 轮：`d140c10` 复核）

复核指出三处问题，均已修复：

1. **审计 example 测试曾无法编译**：`examples/audit_irrep_subduction.rs` 的
   `geometry_zero_is_an_independent_check_never_a_skip` 里还有一处
   `target.irnumber == trivial_cir`（`Option<u32>` vs `u32`，E0308）。R1 的类型改动
   当时只批量修了 `target.ml/row_ml/irnumber` 形态的断言，漏掉这个局部变量比较；
   现在写成 `Some(trivial_cir)`，`--example audit_irrep_subduction` 的 17 项可复现，
   并把它重新纳入本轮基线（此前"example 17 passed"沿用的是 R1 之前的运行结果，
   报告口径已更正）。
2. **Γ 便捷入口返回占位 CIR 号**：`complex_targets` 的普通目标写死
   `irnumber: Some(0)`，于是黄金用例 221 `GM4+` P1 的 `GM1+`/`GM2+` 在 Γ 入口是
   `Some(0)`、在 full-star 入口是 `Some(4075)`/`Some(4076)`。现在读
   `record.source_identity()` 的真实 `cir_irnumber`（非普通标量行直接报
   `UnsupportedCharacterSpace`），并新增永久测试
   `gamma_and_full_star_entries_agree_on_every_target_identity`：SG 221 全部 Γ 标量
   probe 上两个入口的 `(component, ml, irnumber, dimension, multiplicity)` 逐项相同，
   且普通目标的编号必须 > 0。
3. **构造重数查询会返回假零**：`FullStarBlock::constructed_multiplicity` 原来按原始
   折叠坐标比较，`(-1/4,-1/4,-1)` 查不到、等价的 `(3/4,3/4,0)` 才命中。现在只接受
   **完整目标身份** `SubductionComponent::Constructed { q, index }`（q 是构造时就按
   子群倒格规范化过的身份），原始坐标不再是一个可用的查询键；`a_child_p1_gap_is_answered_by_constructed_targets`
   增补回归：ordinal 1045 的两个构造成分身份为 `(3/4,3/4,0)` 与 `(1/4,1/4,0)`、各重数
   2，块身份等于其 `q()` 的规范化形式，其它身份返回 0。
4. 复核同时纠正报告口径：剩余缺口的 **probe / 缺失星 / 全部折叠星**是三个计数
   （例：#5 = 1,659 probe / 2,547 缺失星 / 2,620 全部折叠星），R3 排优先级按 probe；
   `docs/subduction-gap-census.md` 已补 127 个子群的逐项表。

复跑（本轮实测）：全表 `--require-complete` 仍为 `full_success=353382
identity_only=12878`、`hard_failures=0`、exit 0（Γ 入口修号不影响审计路径）；
lib 392、integration 154、doctest 27、audit example 17、census example 1、
严格 all-target clippy 干净。

### 分导 R1/R2 边界复核（`6d1330f` 后）

- 原三处反例已关闭。Γ 身份回归进一步限定为 SG 221 的十个 Γ 标量 probe，
  钉住数量并要求每个调用成功；删除 `Err => continue`，防止零条检查也通过。
- 构造重数查询拒绝非 `Constructed` 身份（返回 0）。此前传入 `Ordinary` 会
  命中第一个存储目标；黄金用例的永久负例先复现错误返回 1，再由分支检查修复。
- R3 的实际剩余分母是 756 个子群/k-star、899 个子群/k-star/setting 组合，
  来自 `target/r12_gaps.tsv` 的 17,144 个缺失星；原 989 组含已经解决的子群 #1。
  本轮仅收紧查询与回归，不改变分导计算或展开 R3。

本轮实测：lib 392 / 4 ignored、integration 154、doctest 27、audit example 17、
census example 1、严格 all-target clippy 通过；Python 离线 9 + 16、几何 oracle
62 行 / 26 描述串 / 62 origin、w 源门禁通过。未重跑全表分导扫描；353,382 / 12,878
覆盖数字沿用 `6d1330f` 的审计，剩余 manifest 的计数与分组本轮已重新核对。

### 分导 R3（2026-09-22）：剩余缺口的来源分类

R3 是离线清点卡，不改生产求解算法。工具 `scripts/classify_subduction_gap_sources.py`
读 `examples/census_subduction_gaps.rs` 生成的缺口 manifest 与 `scripts/iso_irrep_exact.py`
载入的归档 PIR/CIR 帧，对 899 个（子群, k-star, setting）组逐组给出：命中的参数化 PIR
记录与代入参数、被剔除的"更一般域"记录、精确小群操作与有限小余群阶、因子系统
`ω_ij = exp(+2πi q·L_ij)`、命中记录的 irtype 与矩阵块可用性、以及最小见证。

结果（复核修复后重算）：**215 解析路线**（小余群阶 1，目标是一维 Bloch 相位，与 R2 的
child #1 同一条路）、**684 来源候选**（归档 PIR 参数域精确穿过该 q；候选数不等于参数
求值已验证）、**0 无源分类**、0 未分类；**322 组有非零因子系统相位项**、444 组含
分数平移代表元（并未逐项识别螺旋/滑移）；811 组按"自由方向最少"筛选候选，尚未
验证候选一般点小群与实际小群相等，也未验证目标复成分完整性；
**899/899 组矩阵块完整**（90,624 个矩阵元，用生成器已有的 PIR 解码器判定）。

本轮钉死的语义（都曾产生错答）：PIR 是物理不可约表示（10,294 条 = 4,777 离散 +
5,517 参数化）、域维数看方向向量而非参数槽；
`record.operations` 是整个空间群不是小群；k 域参数按**原胞倒格**周期化（C 心群
`(0,1,0)` 方向要 `t=2` 才回同类，修正后 `special_value_no_source` 42 → 1）；
格归属必须用原胞格（螺旋轴乘积与代表元差 centring 矢量）。

逐操作见证（验收项）：带心 SG 5（C2，U 线 `(0,t,1/2)`，归档 `U1UA1`/`U2UA2` 在 t=1/2 命中、
一般位置记录被剔除）、螺旋/滑移 SG 24（I2₁2₁2₁，P 点小群 4 操作、3 个螺旋、
非零相位圈数 φ∈{1/4,3/4}，对应 ω=±i）
在 Python 单测；非对称换基 ordinal 26（SG 3 `A1` → #3，`U=[[1,2,1],[-1,2,-1],[-1,0,1]]/2）
在 `tests/subduction_gap_sources.rs`，用引擎 `unmap_operation` 逐操作对照归档 child #3
操作并核对钉子星 `(0,1/3,1/2)` 的小群与相位。

命令：`python3 scripts/classify_subduction_gap_sources.py <manifest> > groups.tsv`；
`python3 -m unittest discover -s scripts -p test_classify_subduction_gap_sources.py`（21 项，
默认跳过全表项）；`R3_FULL_MANIFEST=1` 跑含 899 组门禁的完整套件。
报告：`docs/subduction-gap-sources.md`。R4 复用 PIR/CIR 解码器，为 684 个候选补参数与
相位求值，并验证适用域和目标完整性；本卡只交付候选清点。

### 分导 R3 复核修复（第 4 轮：`52e4c65` 复核）

复核独立复现了 899 组分类，并否证了"#155 唯一无源"的结论。四处问题全部修复：

1. **参数代入漏掉非对角分量（P1）**：`arm_point` 只算 `direction[axis]·t`，方向
   `(1,1,0)`、`t=1/4` 返回 `(1/4,0,0)`。改为 `k = constant + Σ_j t_j p_j` 全分量求和。
   #155 的归档 `Y1YA1`/`Y2YA2` 耦合直线 `k=(t,t,3/2)` 在 `t=3/4` 命中见证星
   （差 R 心倒格矢量 `(1,1,0)`），该组因此是参数化来源。
2. **旋转求逆少一次转置（P1）**：`rotation_inverse` 返回 `C/det = R⁻ᵀ`，调用方当逆矩阵
   用导致倒空间作用错误，7 组（#155/#166/#167）实际阶 2 被标成阶 1。改为返回真正的
   `R⁻¹`（余子式矩阵转置后除行列式）。
3. **矩阵可用性判定错位（P2）**：`irtranslations` 是参数化相位字段，离散记录按格式没有
   （#5 Γ `GM1`/`GM2` 四槽全 `None` 但矩阵完整），不能当"矩阵空洞"。改为调用仓库已有的
   PIR 解码器（`generate_irrep_data._parse_pir_characters`）判定：899/899 组矩阵块完整。
   列名改为 `irtranslation_slots`/`irtranslation_none` + `matrix_available`/`matrix_elements`。
4. **Rust 因子系统把格矢消掉（P2）**：`factor_system_turns` 先约化乘积再取平移差余数，
   携带 Bloch 相位的格矢丢失（`S²=T(0,0,1)`、`q=(0,0,1/2)` 返回全零）。改为乘积不约化、
   直接取 `s_i s_j` 与代表元之差的格矢并断言属于子群格。

永久回归（新增）：耦合方向代入、`R·R⁻¹=I`（含三方旋转）、同星各臂小群阶一致
（`star_order_inconsistent=0`，输出列 `star_orders`）、非零螺旋相位负例、#155 参数化
来源钉值、矩阵可用性走解码器。修复后重算：215/684/0，71 组命中来源改变；
`scripts/test_classify_subduction_gap_sources.py` 19 项（含 899 组门禁），
`tests/subduction_gap_sources.rs` 2 项。

### 分导 R3 复核修复补记（第 5 轮：`fc5fb0f` 复核）

1. **星内小群阶门禁原来没生效**（P2）：全表门禁用位置下标取列，`row[-3]` 是
   `matrix_elements`、`row[14]` 只是单值，注入 `star_orders=1,2` 也能通过。改为
   `csv.DictReader` 按列名断言 `row["star_orders"] == row["little_co_group_order"]`
   并保留跨 setting 的阶一致性检查。新增列重排后注入 `1,2`、错误单值、空值和
   跨 setting 冲突的负例；门禁与负例共用同一校验函数。显式启用全表门禁而缺少
   manifest 时失败，不再跳过。
2. **报告口径纠正**（P2）：删掉"PIR 只索引参数化域/没有离散记录"的错误说明（改为
   PIR = 物理不可约表示，含 4,777 离散 + 5,517 参数化）；676 → **684 来源候选**；
   CIR = 复不可约表示，含 5,296 离散 + 5,906 参数化；444 组只说明所选代表元含
   非整数平移分量，不能称为非幺正操作或已识别的螺旋/滑移。811 组按自由方向最少
   筛选候选，不是小群最大性验证；相位圈数 φ 与复因子 ω 分开说明。
3. **帧契约确认（复核给出证据链，勿再重复施加 U）**：冻结 `U⁻¹` 已进入嵌入矩阵 `T`
   （`src/irrep/subduction.rs:1352`），折叠在那里算 `q = Tᵀk`
   （`src/irrep/subduction_star.rs:607`），`canonical_q` 之后只做倒格约化；230 个
   空间群的运行时 Hall 选择与冻结来源一致、归档帧到 data-Hall 的变换为 `P=I, p=0`
   （`scripts/iso_irrep_data_hall.py:7`）。因此上一轮"U 未处理"的顾虑作废；R4 要处理的
   是操作代表元的格平移相位与 `child_shift` 回退。

本轮最终验证（在已有修订 `4428601` 上补齐）：分类器完整套件 **21 passed**，
含 899 组门禁（215 解析 / 684 来源候选；矩阵块完整、星阶一致），约 123 s。
对全表测试本身再次注入 `star_orders=1,2` 后断言失败，未注入时通过；跨 setting
冲突回归先复现漏检，再修复为通过。显式开启全表门禁但缺失 manifest 时会失败。
Rust lib 392 / 4 ignored、integration 156、doctest 27、audit example 17、census
example 1；严格 all-target clippy 通过（保留既有 workspace manifest 警告）。
Python oracle/w-source 离线 9 + 16；live 几何 oracle 62 行 / 26 描述串 / 62 精确 origin，
w 源门禁 checks_failed=0。独立核对 230 个 SG 的运行时 Hall 与冻结来源一致且 P=I,p=0。
未重跑普通分导的全表数值审计；`src/` 无改动，分类器执行逻辑也未改动。
R3 按来源候选清点收口；完整目标求值与完整分解覆盖仍待 R4/R5。

### 分导 R4 批次 1（2026-09-22）：平凡小余群的构造目标

R4 按“小批次卡”推进，规则见 `docs/subduction-next-milestones.md` §R4，批次状态与
文件所有权见 `docs/subduction-r4-batches.md`。**批次 1 已交付**：把构造目标的判据从
「子群 #1」换成几何事实「该星的小余群平凡（子群自身 data-Hall 操作里固定在 q 模
子群原胞倒格的旋转只有恒等）」，于是 R3 判为 `analytic_general_position` 的 **215 组
（41 个子群）** 不再需要任何归档字符数据。

三处实现（`src/` 改动仅限这条路线）：

1. `ConstructedLittleRep` 移入 `subduction_star.rs`，与新的 `LittleCharacter`
   （存储行 / 构造表示的统一字符来源）和 `ConstructedStar` 同处；`arm_character`、
   `induced_component_character`、`induced_character` 改走该来源，存储行路径逐位不变
   （lib 392 项全绿）。失败模式 `StarError::ConstructedRotationNotCovered` 保留：
   非恒等旋转仍拒绝作答，不返回相位。
2. `constructed_child_components_at` 的判据换成 `has_trivial_little_co_group`，
   非平凡小余群仍返回空、保持 `MissingChildStarData`（即批次 2 的边界）。
3. `ChildStarEvaluator::Constructed` 现在持有 `ConstructedStar`：子群点群非平凡时
   目标必须按**子群自己的星**诱导，不能像 #1 那样把小群特征标当全星特征标。

**两条星级规则是实测踩出来的，勿退回点级实现**：(a) 星里只要任何一条臂命中 pinned
行，就只用 pinned 行——否则同一物理 irrep 会以两个身份出现，求解器报
“rows ... are not orthogonal”（ordinal 1007 的 `V1`/`L1`）；(b) 一个星只构造一次，
落在该星规范化第一个点上——逐点构造会得到两个字符相同的行（ordinal 1007 的
`(0,1/2,1/2)` 与 `(0,3/2,1/2)`）。`select_representative` 现在先做一次
`star_has_stored_components` 判定，再按整星收集成分。

实测（全表审计 518 s）：`full_success 353,382 → **357,033**`、
`identity_only 12,878 → **9,227**`（−3,651，正好是 R3 解析组的 probe 数）、
`hard_failures=0`、存储恒等正项 94,271/0 不匹配、Γ Frobenius 1,895/1,895、
w 行 5,756 computed / 0 错误、`production_checks` 五项全 0、`accounting_violations=0`；
`--require-complete` 仍 exit 0，`--require-full-decomposition` exit 2
（`incomplete=9227`，R5 目标）。重跑清点：`records=1569 probes=9227 stars=13857
missing_stars=12932 constructed_stars=152 reachable_stars=773`（新状态
`constructed_trivial_co_group`）；在 `target/r4_gaps.tsv` 上重跑 R3 分类器只剩
**684 个 `parameterized_source`**（0 解析 / 0 无源 / 0 未分类，矩阵块 684/684 完整、
81,576 元）——R3 的预测与批次 1 的实际闭合面完全一致。

回归与更新：新增 `tests/subduction_constructed_stars.rs`（4 项：13345 全 31 个 probe
完整且恒等重数对上 pinned 频率、两臂星诱导的 `block_dimension = star_size ×
little_dimension`、13346 的 10 个参数化星仍 `MissingChildStarData`、“有 pinned 行的
星绝不再构造”）；lib 新增 `a_constructed_star_induces_over_the_child_star`（子群 #2
手算 χ(E)=2、χ(T_(1,0,0))=0、χ(T_(0,1,0))=−1、χ(−I)=0）；`census_subduction_gaps`
的可达性模型补上新状态并加两项测试；`audit_irrep_subduction` 的微型基线由 13345
（已清零）移到 13346；`subduction_identity_regressions` 的 15 个缺数据 probe
→ 10 个、比较数 2,060 → 2,065、非 Γ 1,522 → 1,527。验证：lib 392+1、
integration 156+4、doctest 27、audit example 17、census example 2，严格
all-target clippy 通过（保留既有 workspace manifest 警告）。

未做（批次 2 入口）：684 个参数化来源候选的参数代入、小群字符/矩阵求值、适用域与
目标完整性验证；参数化来源可复用时冻结最小数据，运行时保持纯 Rust。`little_k`
（w 源数组）与 PIR `k_arms` 是不同数组，坐标约定须分别验证。

### 分导 R4 批次 2a（2026-09-22）：一维投影特征标 catalogue

**已交付**：R3 判为 `parameterized_source` 且小余群**非平凡但一维可解**的那批不再需要
归档字符，由引擎用精确 cocycle 现场求 catalogue 回答。文件：新增
`src/irrep/subduction_catalogue.rs`（小余群 + 因子系统 + 一维投影特征标求解器），
`subduction_star.rs` 增加 `ConstructedLittleRep::Projective { q, constants }` 与
`StarError::LittleCoGroupNotClosed`，`subduction_star_decompose.rs` 增加
`constructed_projective_components` 与公开谓词 `constructed_targets_available`
（离线 census 用它，避免模型与引擎漂移）。

判据与约定（都可复算）：

* `omega_ij = exp(2 pi i q.L_ij)`，`L` 取**未约化**乘积缺陷；解
  `psi_i + psi_j - psi_k == phi_ij (mod 1)`，**只有解数恰为 `|P_q|`**（`Hom(P_q,U(1))`
  陪集大小，即所有不可约投影表示都是一维）才使用 catalogue；非上边界、非交换、
  `|P_q| > 4` 一律返回空、保持 `MissingChildStarData`（fail closed，绝不猜）；
* 恒等旋转的代表元必须是**零平移**的恒等操作，否则会把 centring 平移的 Bloch 相位吃掉；
* **因子系统、常量与取值必须共用同一个约化后的 q** —— 这是本轮修掉的 bug：第一版把
  常量建在原始折叠点、却用约化点求值，30 个 probe 得到复数重数 `1±i`（ordinal 14090，
  SG 226 W5 → #98）。原因是重建出的小群操作**不是**其代表元的格平移（平移可带 1/4），
  两半相位混用会差一个非整数；引擎当时 fail-closed（报错而非给错值），修好后 14090
  变成 25/25 完整、`VERDICT clean`。

实测（全表审计 560 s）：`full_success 357,033 → **366,039**`（99.94%）、
`identity_only 9,227 → **221**`、`error=0`、`hard_failures=0`、恒等正项 94,271/0 不匹配、
Γ Frobenius 1,895/1,895、w `5756/5756`、`production_checks` 全 0；
`--require-complete` exit 0、`--require-full-decomposition` exit 2（`incomplete=221`）。
重跑清点：`records=83 probes=221 stars=331 missing_stars=326 constructed_stars=2
reachable_stars=3`（新状态 `constructed_target`）；在新缺口上重跑 R3 分类器只剩
**58 个 `parameterized_source`**（29 个子群，星阶 4 的 52 组 + 星阶 6 的 6 组），
与侦察预测的 2b 集合一致。

证据与回归：`the_catalogue_reproduces_pinned_little_group_characters` 在**离散** pinned
k 点上对照引擎字符行，**1,176 条记录 / 6,318 个操作全部命中**（另 1,812 条属更高维批次
跳过，2.8 s）；`a_two_fold_co_group_has_two_characters` 手算 C2，含 ψ=1/6 的"规范可以比
cocycle 更细"情形；`tests/subduction_constructed_stars.rs` 的 13346 全部回答 / 3988 仍
缺数据；审计微型基线由 13345/13346 移到 3988；恒等回归 2,075 个 probe 全部走完整入口
（缺数据集合为空）；settings 一组 120 个 probe 全部完整。验证：lib 395、integration 161、
doctest 27、audit example 17、census example 3，严格 all-target clippy 通过。

未做（批次 2b）：58 组需要**二维**投影不可约表示（52 组 `|P_q|=4` 单 ω-正则类、
6 组 `|P_q|=6` 的 D3），需 twisted group algebra 的二维不可约表示或等价的诱导构造，
单独建证据集；`little_k`（w 源数组）与 PIR `k_arms` 的坐标约定仍分别验证。

### 分导 R4 批次 2b（2026-09-22）：二维投影表，普通离散标量覆盖闭合

**已交付**：最后 58 组（52 组非退化 C2×C2 + 6 组 D3，共 221 个 probe / 326 个缺失星）
由 `catalogue::projective_targets` 现场给出投影字符表；全表
`full_success=366260/366260`、`identity_only=0`、`error=0`、两个门禁 exit 0、
`VERDICT complete`（约 511 s）——即 R5 的验收清单（恒等正项 94,271/0 不匹配、
Γ Frobenius 1,895/1,895、w 行 5,756/0 错误）全部满足。

结构先离线钉死（脚本可复算）：

* **52 组：非退化 C2×C2** —— 元素阶 (1,2,2,2)、交换、交换子配对
  `beta_ij = phi_ij - phi_ji` 的根只有单位元 ⇒ twisted algebra ≅ M₂(ℂ) ⇒
  **唯一二维不可约表示，字符 (2,0,0,0)**。证明（复核 B 修正：**不能**用
  `omega(g,g)=±1`——`g²=e` 只对旋转部分成立，cocycle 是代表元乘积的格平移相位，
  实测出现 1/6、1/4、1/3、2/3、3/4、5/6，child #43 上 `omega(g,g)=i`）：
  `u_g² = omega(g,g)·1` 且 `omega(g,g) ≠ 0` ⇒ 特征值 `±sqrt(omega(g,g))` 相异；
  `u_g` 不能是标量（标量 ⇒ `omega(g,h)=omega(h,g)` ∀h ⇒ g 与全群正交，与配对
  非退化矛盾）⇒ 两个一维特征空间 ⇒ 迹 0。
* **6 组：D3 且 cocycle 为上边界** —— 元素阶 (1,3,3,2,2,2)、非交换、配对恒 0、
  ψ 解恰 2 个（`|Hom(D3,U(1))|`）⇒ 目标 = ψ 规范 × 普通不可约表示 `{1,1,2}`；
  用哪个 ψ 解不影响**集合**（两解相差 sign 特征标，sign ⊗ 普通不可约只是置换）。

实现与门禁：新增 `ProjectiveTarget`、`projective_targets`（`MAX_ORDER` 由 4 提到 6，
给 D3 提供 ψ 规范；`MAX_GRID`/`MAX_WORK` 仍在，超限返回空）、
`ConstructedLittleRep::ProjectiveTable { q, constants: Vec<(Mat3I, Complex64)> }` 与
`table_character_value`（`chi = constant(R)·exp(2πi q·T)`，与一维表共用"同一个约化点"
的规范纪律）。两族各有结构门禁（交换性/元素阶/配对根/ψ 解数）与正交性门禁
（每个字符 `(1/|P_q|)Σ|χ|² = 1`、`Σ dim² = |P_q|`）；范围外（如立方 Γ 点，阶 48）
返回空表、入口保持 `MissingChildStarData`。

证据：`the_projective_tables_cover_only_the_two_gated_families`（#43 得唯一二维目标且
整行字符钉为 `(2,0,0,0)`、#160 得 `{1,1,2}` 且三行整表钉为
`(1,1,1,1,1,1)`/`(1,1,1,-1,-1,-1)`/`(2,-1,-1,0,0,0)`、#221 Γ 为空）；扩展后的 pinned 对照
**1,328 条记录 / 7,578 个操作全部命中，其中 94 条走二维投影表**（1,660 条范围外跳过），
四个数在 R5 收口后由 `assert_eq!` 钉值；端到端 ordinal 3988 = 12/12、12041 = 19/19 均
`VERDICT clean`；审计 example 的批次见证改为 13345/13346/3988 全部 exit 0。

边界：缺口清零后 `census_subduction_gaps` 的输入为空；工具现在用精确 token 识别
marker、摘要给出 `audit_rows/probe_rows/identity_only_rows/unanswered_probe_rows/closed`，
并以 `--require-empty` 作为门禁（R5 审计上 `probe_rows=366260 identity_only_rows=0
closed=true` exit 0）。真实数据里已不存在"仍缺数据"的端到端负例，fail-closed 由
`the_one_dimensional_solver_returns_nothing_instead_of_a_subset` 与
`an_out_of_scope_co_group_still_reports_missing_child_star_data`（真实嵌入 + 手工越界星，
见复核 B 那节）保证。验证：lib 396、integration 162、doctest 27、audit example 17、
census example 3，严格 all-target clippy 通过（本轮当时值；处理后见文件开头的当前基线）。

### 分导 R5 收口与两个对抗性复核（2026-09-23）

R5 的正式报告、独立性表、100% 覆盖结论见 `docs/subduction-audit.md` 的
「R5：普通离散标量覆盖闭合」一节；批次证据见 `docs/subduction-r4-batches.md`。
本节只记这一轮复核的处理与仍不可验证的边界。

两个独立 reviewer（A：覆盖声明与门禁真实性；B：实现与回归质量）都只读、都不提交，
结论均为**可接受、无阻断项**。B 用不依赖本 crate 的路径独立复现了头条数字：
自建 cocycle 与扭群代数中心幂等元分解（58 组 → `(4,(1,2,2,2),(2,)): 52` +
`(6,(1,2,2,2,3,3),(1,1,2)): 6`）、直接读 `iso.zip` 重算 15,239 记录 / 94,271 条 /
4,777 标量 irrep / 366,260 分母、独立重跑三门口禁（535 s，exit 0，数字逐项一致）。

已按发现修改（全部在主线程复现后处理）：

1. **P1（B）：C2×C2"唯一二维"的证明前提为假**。原文用 `g²=e ⇒ omega(g,g)=±1`，
   但 `g²=e` 只约束旋转部分，cocycle 是代表元乘积的格平移相位；实测 58 组里
   `omega(g,g)` 取 1/6、1/4、1/3、2/3、3/4、5/6，child #43 `q=(0,1,1/2)` 上为 `i`。
   结论不变：`u_g² = omega(g,g)·1`、`omega(g,g) ≠ 0` ⇒ 特征值 `±sqrt(omega(g,g))`
   相异；`u_g` 非标量由配对非退化保证 ⇒ 迹 0。代码注释、批次卡、本文件与
   `1d1b95e` 的提交信息都以新论证为准。
2. **P2（B）：构造目标的独立证据只有一条，且原先只有下限断言**。构造目标无 CIR
   来源 ⇒ "来源身份/标签一致"两项对它恒真；子群 Γ 星永远命中 stored 行 ⇒ 恒等频率
   比较永不覆盖构造星。`the_catalogue_reproduces_pinned_little_group_characters` 的
   1,328 / 7,578 / 94 / 1,660 改为 `assert_eq!`，R5 报告新增"构造目标的外部对照"小节。
3. **P2（B）：端到端 fail-closed 负例在三个批次里被替换殆尽**。补
   `an_out_of_scope_co_group_still_reports_missing_child_star_data`（真实嵌入
   225 `X1+` P3 → #221，把手工折叠星放到 16 阶点 (0,0,1/2)，`build_block` 必须报
   `MissingChildStarData`）与一维 solver 三条负方向（非上边界 / `MAX_ORDER` /
   `MAX_GRID` 均返回空表而不是子集）。
4. **P2（A+B）：census 空输入无法区分"闭合"与"漂移"**。marker 改为精确 token，
   摘要打印真实计数（删掉恒为 0 的 `replay_errors`/`dimension_errors` 常量），
   新增 `--require-empty`（至少一条 probe 行、全部被引擎回答、无 identity-only 行），
   6 条 replay 拒绝路径各一条负例 + 1 条 marker 漂移负例；残留限制（被改名的 marker
   若其余字段合法则无法识别）写进工具文档。
5. **P2（A）：五项生产检查的独立性**。R5 报告逐项标注"同义反复/语料上不可达"，
   并把接线测试改名为 `every_counted_production_violation_is_a_hard_failure`，
   注释写明它只测 `Counts::hard_failures`/`exit_code`、不驱动生产自增点。
6. **P3（B）**：`ConstructedStar::dimension()` 不再硬编码 1（二维族现在报 2，带回归）；
   删除 `subduction_catalogue.rs` 的重复 doc 行与 `4 != order` 死代码；删除 7 处过期的
   `#[allow(dead_code)]`（`LineArmSource`、`ArmCharacterSource`、`line_folded_stars`、
   `FoldedPoint/FoldedStar::from_parts`、`MissingFrozenRotation` 实际都在调用链上，
   clippy `-D warnings` 仍零警告）；census 错误路径不再 dump 整个 `FullStarSubduction`。
7. **存储行契约（B）**：SG 38/40 有 4 条 pinned 行对小群之外的 `-I`/`m_y` 存字面 0
   ⇒ "row = 每个列出操作的特征标"不成立，已写入 `docs/subduction-conventions.md` §14。
   `stored_child_components_at` 对无法展开的记录保留大声报错（不 `continue`），
   理由同节注明。

主线程另补做了复核 A 的"无法验证"清单里的两项（记录在
`target/r5_review_a_unverified.md`）：重跑 `scripts/classify_subduction_gap_sources.py`
得到 899 = 215 解析 + 684 参数化、684 组/123 子群/81,576 矩阵元、批 1 的
215/41/4,212/3,651/119、批 2a 的 58 = 52+6、221、83、29、331 = 326+2+3，重跑产物与
冻结 `target/r4_groups.tsv` **SHA-256 相同**；重跑 `scripts/freeze_w_little_characters.py`
（官方 `iso`）得 `sources solved: 73 | failures: 0`，生成的 Rust 表与提交的
`src/irrep/w_little_characters_data.rs` **逐字节相同**。

仍不可验证（如实保留）：几何 oracle 只有 22 组 (SG, irrep) / 62 行抽样（15,239 行的
0.31%），扩大它要对 4,777 个普通 irrep 各起一次官方 `iso`，是独立工作量；B 也无法
验证历史增量数字（需逐提交重跑）与"114,770 条引擎零项"（absent 探针无归档 oracle）。

### 第三方复核（2026-09-23，公开仓库浏览 + 隔离复现）

第三个独立 reviewer 审阅 `0386cfa`..`ff2e4a4`（R5 收口提交之后），给 5 条发现 + 3 条
提醒，未质疑 R5 的 366,260/366,260 与恒等正项 94,271/0 不匹配。逐条已复现并处理
（细节与证据表见 `docs/subduction-audit.md` 的「第三方复核」小节）：

1. **P1：旧 settings 生成器与六字段全表不兼容**（`generate_subduction_settings.py`
   仍是 69+6 条五字段；`--write` 会覆盖 15,239 条正式模块）。复现时发现
   `scripts/test_subduction_settings.py` **当时就是红的**（2 个 setUpClass error）。
   处理：`parse_committed` 支持六字段；`--write` 改为拒绝执行并指向
   `scripts/task9/build_table.py`；`--check` 改为校验它自己的 75 条遗留记录 + 模块
   覆盖全部 15,239 个 ordinal（在线实跑通过）；数据模块头与生成器头同步改指 task9；
   该测试重写为 10 passed。
2. **P1：嵌入验证缺 `L_H ⊆ L_G`**。用公共 API 复现：SG 1 自嵌入 + `U = diag(2,1/2,1)`
   （`det U = 1`，唯一有限操作映到自身）此前被 `probe_embedding` 当 `Ok` 返回。
   处理：`validate_candidate` 先检查子群格每一行都在母群格里；新增负例
   `a_setting_whose_lattice_is_not_a_parent_sublattice_is_rejected` 与全表回归
   `every_frozen_setting_keeps_the_subgroup_lattice_inside_the_parent`（15,239 行
   逐一重算，分数行 ≥ 30 条；审计 embedding 仍 15,239/15,239，未误伤）。
3. **P2：结果对象丢 setting 分母**（`IrrepSubduction`/`FullStarSubduction` 只有分子；
   ordinal 26 的分母是 2）。处理：两个结果类型新增 `setting_denominator()`，
   回归 `a_fractional_setting_reaches_the_results_with_its_denominator` 在 ordinal 26
   上钉住 `setting()/setting_denominator()`。
4. **P2：`probe_subduction_settings.rs` 把失败写成 `trivial=0`**（与"真的 0"混淆，
   而 0 在这条流水线里表示约定落在共轭分支）。处理：三态 `<n>` / `?` / `error`
   （原因到 stderr），`--profile` 失败输出 `profile=error`，两条单测。
5. **P2：`scripts/task9/build_table.py` 先写后验、不去重、空输入 exit 0**。处理：
   重写为验证后写（重复 ordinal 报错、必须覆盖既有模块的 ordinal 集合或 `--expected`、
   未覆盖记录按 status 报告、`os.replace` 原子替换）、新增 `--check` 与显式
   `--partial`，`--shifts/--derived` 接受管线三种形状且可重复；新增
   `scripts/test_build_table.py`（5 条回归）。**诚实边界**：历史那张表的逐字节重生成
   需要未跟踪的 `target/task9/` 证据链（最终 shift 是多次增量修补的合并，2,703 条
   非零），仓库只声称内容可核（离线覆盖/身份/幺模/最简 + 引擎 15,239/15,239 +
   75 条遗留记录在线复推），见 `scripts/task9/README.md` 的 Provenance status。
6. **提醒（空 scope）**：`--parent 2 --ordinal 0` 之前打印 `VERDICT clean ...
   probe_full_success=0/0` 且 exit 0。处理：审计在范围内零记录时报 `empty scope`
   并非零退出，回归 `an_empty_scope_is_rejected_instead_of_reporting_clean`。
7. **提醒（可移植性）**：`scripts/task9/README.md` 的 `../../../scripts/...` 多退一级
   （从 `target/task9` 应是 `../../`）；12 个 task9 脚本硬编码 `/home/liuyichen/...`。
   处理：README 路径改正、脚本改为从 `__file__` 推导 `REPO`；并把只在
   `target/task9/` 存在、却被三个 tracked 脚本 import 的 `derive_shift.py` 入库；
   全部脚本可 import。
8. **提醒（文档混用历史与当前）**：审计文档「当前全局基线 exit 2 /
   identity_only=14713」改为显式历史 + 当前状态；`scripts/task9/README.md` 里三处
   w 门禁的"仍 exit 2 / 仍无法计算"也标注为历史（现在 `w_computed=5756/5756`、
   门禁 exit 0）。

### 第三方复核跟进（2026-09-23，针对 `c37d3b9`）

同一 reviewer 复核上一轮修复，认可 `L_H ⊆ L_G`、分母传递、probe 错误语义与空 scope
四处，但指出生成管线两处仍可绕过的门禁（都已在主线程复现、修复、并用 reviewer 提供的
`gate_regression_tests.py` 验证 4/4 通过）：

1. `build_table.py` 在 `--out` **不存在**时 `expected=None`，完整性检查被跳过，
   空输入会写出零条目表并 exit 0（我原来的测试只覆盖"已有非空基线"）。现在预期全集
   来自 `--expected`、已存在的 `--out`、或 tracked 的 `--baseline` 模块（默认
   `src/irrep/subduction_settings_data.rs`）；三者不一致报错，全部缺失时拒绝写入，
  只有显式 `--partial` 例外。reviewer 的最小复现命令现在 exit 1 且不创建输出；
   我自己的 `test_build_table.py` 加到 7 条覆盖这两个新分支。
2. `generate_subduction_settings.py --check` 的覆盖检查只比**条数 + 逐条身份**，
   `[0,1,1]` 对 3 条记录能顶替缺失的 ordinal 2（`compare()` 以 ordinal 为键会合并
   重复）。现在抽出共享实现 `check_ordinal_coverage`（ordinal 序列必须严格等于
   `range(len(records))`，并报缺失/重复/越界），**在线 `--check` 与离线测试调用同一
   函数**，杜绝"测试严格、生产宽松"的分叉。

reviewer 同时确认：历史表逐字节重生成属于**已披露、未解决**的独立问题（需要归档输入、
合并顺序与选择依据），与"当前冻结表可离线核验 + 引擎可嵌入"是两个层次的保证；这一
区分保留在 `scripts/task9/README.md` 的 Provenance status 与审计报告的诚实边界里。
本轮只改 Python 脚本与文档，Rust 侧未动，故不重跑全表审计（`c37d3b9` 的 517.1 s、
三门口禁 exit 0 仍然有效）。

### R6.0/R6.1：参数化 k 的完整分解（2026-09-23）

计划、能力 A/B 的分界与验收清单见 `docs/subduction-r6-plan.md`；契约在
`docs/subduction-conventions.md` §16。本轮交付：

* **R6.0 契约**：`OFFICIAL_LINE_PARAMETER = (1, 4)` + `official_line_parameter()`，
  新公开入口 `subduce_line_at_parameter(subgroup, embedding, table, parameter)` 与结果
  类型 `LineSubduction`（`parent_dimension`、blocks、`reconstruction()`、
  `trivial_content()`、setting 分子/分母、参数与波矢）。输入帧：`t · table.direction`
  用冻结表自己存的方向；不支持一律显式报错（`LineSourceMismatch`、
  `MissingChildStarData`、`MissingChildTrivialIrrep`、`TargetSourceMismatch`），
  **绝不以 0 代替失败**。
* **R6.1 端到端**：先在 SG 196 的 106 条 pinned 行上跑通（单元测试），再把审计的
  w 门禁整体换成"官方参数下的完整分解"——5,756/5,756 条 pinned 行完整分解且恒等重数
  等于 pinned 频率（604.8 s，三门口禁 exit 0，其余计数与 R5 完全一致）。
* **修掉一个 R5 遗留的结构错误**：旧的 `line_folded_stars` 把每个约化 `q` 当成一个
  子群星，而 `FoldedStar` 的语义是子群点群下的**轨道**；Γ-only 路径不做重构，所以
  这个错误一直不可见，第一次在一般参数上做完整分解时重构在恒等元给出 12/8 而不是 6
  才暴露。现在线源改用离散路径共享的 `fold_arms`（轨道划分 + 臂数一致性检查），
  两个源的折叠几何统一。
* **分支与失败语义的见证**：`10038`（SG 196 `W1` `4D1` → #1）在 `t = 1/4` 是
  `Z1`×2 + `GM1`×4、恒等重数 4 == pinned，在两侧 `t = 1/6`、`t = 1/3` 变成 6 个
  构造块、恒等重数 0（都等于独立于多重度求解器的几何计数）；`t = 1/4, 3/4, 5/4`
  在该源上逐目标相同（只钉这一例，一般位移规则留 R6.2）；`13543` 的 `DT5` 在
  `t = 1/7` 报 `MissingChildStarData{sg:136}`，同记录在 `t = 1/4` 正常——证明是参数
  问题而非记录问题。
* 验证：lib `409 passed / 4 ignored`、测试二进制 `570 passed`、doctest 27、审计 example
  19、缺口清点 7、settings 探针 2，严格 all-target clippy 零警告；全表审计
  557.0–604.8 s exit 0。

### R6.1 复核处理（2026-09-23，reviewer C/D）

第二轮 adversarial review 的两份报告（C：数学/结构；D：契约/声明/测试强度）逐条复核
后处理，结论与证据：

* **P1-1（真 bug，已修）gauge 不一致**：`t` 与 `t + Δ`（`Δ·direction` 是母群倒格矢）
  是**同一个母群 irrep**，但引擎把 `k` 的原始值直接喂给冻结的 Γ 点字符 `D`，等价参数
  会解出不同的子群 irrep：ordinal 11328（SG 210 `DT3`）在 `t = 1/4` 给 `Z1`、`t = 5/4`
  给 `Z2`，恒等重数都是 1（== pinned）、两侧重构都通过，**任何门禁都看不见**。
  受影响 40 行（SG 210 的 11328–11333、SG 227 的 14430/14432/14504/14506、SG 228 的
  14723/14726，`DT1`–`DT4`）。修法：新增 `canonical_wave_vector`（把 `k` 约化进
  `Lattice::new(exact_primitive_basis(parent)).reciprocal()`），每条臂取"规范中心 `k`
  在该臂自身旋转下的像"再折叠。**对 5,756 条 pinned `t = 1/4` 波矢零位移**，所以 R5
  已验收的约定原样保留；新增回归
  `a_reciprocal_vector_shift_of_the_parameter_changes_nothing`（40 行 × `t = 5/4, 9/4,
  13/4` 与 `t = 1/4` 逐项相同），并复现 11328 三个参数下逐位一致。
* **P1-2（已补）轨道化折叠缺常驻回归**：新增
  `a_generic_parameter_folds_into_multi_point_child_stars`——`10030 DT1` 在 `t = 1/7`
  给出 3 个块、每块 2 点 2 臂、`parent_dimension = 6`、恒等重数 0；`13543 SM1` 的星
  大小为 [4, 8]。reviewer D 复现过"把修复换回旧的按 q 分组，其余测试与整个审计全绿"，
  这两条就是那个漏洞的钉子。
* **P1-3（措辞）**："独立手算"是过度表述：计数只用另一条多重度读法，臂集合/帧/Γ 判定与
  引擎共用。测试与 §16 已改为"独立于多重度求解器的几何计数"，并明确一般参数上**只有
  内部一致性 + 这条几何计数**，外部 oracle 仅在 `t = 1/4`（pinned 频率 + live oracle）。
* **P2-1（已写进契约）能力 A 的参数域**：一般 `t` 下 54/5,756 行（`t = 1/7, 1/6, 1/3,
  2/7`）报 `MissingChildStarData`（子群 #123–#138 的 2 点星、#221–#224 的 6 点星），
  `t = 3/8` 报 30 行；小余群非平凡时的网格搜索有 `MAX_GRID = 200_000` 上限——实测
  `10030 DT1` 在 `t = 1/100000` 仍 `Ok`、`t ≥ 1/1000000` 起报错（边界正是 `2 × 10^5`，
  与 `|P_q| = 2` 下的 `MAX_GRID` 逐位吻合）；`t = 0` 及"`t·v` 稳定子更大"的参数返回
  **形式值**（`10030 DT1 t=0 → Ok(3)`、`13543 DT5 t=0 → Ok(0)`），不是错误也不是物理
  结论。全部写进 §16，属 R6.2 的分区输入。
* **P2-2（已改）帧的文档自相矛盾**：`docs/subduction-r6-plan.md` §4 与 CLAUDE 第九轮
  的"必须换算到 primitive"是原型阶段的错措辞（第十一轮已定位 12 臂的根因是约化格用错）。
  实测帧：`direction`、pinned `little_k`、冻结 little 群操作全程 conventional，引擎
  没有任何 primitive 换算。§4 已按实测重写，第九轮处加了更正块。
* **P2-3（已改）`docs/subduction-audit.md` 过期**：w 分支的机制改为
  `subduce_line_at_parameter(..., official_line_parameter())` 的完整分解、
  `LINE_PARAMETER` → `OFFICIAL_LINE_PARAMETER`、复现一节的"约 520 s"改为三个实测值
  并注明墙钟差属负载。
* **P2-4（已注明）两个不可达失败变体**：`MissingChildTrivialIrrep`（唯一恒等 Γ 行对
  230 个 SG 钉死）与 `TargetSourceMismatch`（两个读数同源）是**表损坏防御分支**，
  公网 API 在 pinned 数据上不可达，没有也无法写负例；文档不再与两条可达变体并列。
* **P2-5（已改）访问器 doc**：`LineSubduction::blocks()`（每块是一个**子群轨道**，不是
  每个 `q`）、`FullStarBlock::stored_k()`（constructed 块给的是约化后的折叠坐标）、
  `FullStarBlock::q()`（**未约化**，跨参数比较必须先约化或改用 `wave_vector()`）。
* **P3-1/3-2/3-3/3-4（已改）测试与死代码**：`assert_line_invariants` 的 doc 注明它重算
  的是构造时已强制的不变式（重构钉子，不是正确性证据）；`the_block_route_is_the_public_route`
  从 `f(a) == f(a)` 改成对 `subduce_line_at_parameter(...).trivial_content()` 与 pinned
  频率两条比较；6 处 `Rat::new(1, 4)` 字面量改用 `official_line_parameter()`；退役的
  `FullStarError::NonIntegralLineFrequency`（已无构造点）删除；"特殊值 + 非平凡小余群"
  明确记为**无单元见证**、只由 5,756 行审计覆盖。
* **P3-5（已记）计时**：同二进制 5,756 行实测完整分解 10.97 s 对 Γ-only 6.13 s +
  0.01 s ⇒ 升级净代价约 +4.8 s，不是 +88 s；审计墙钟 557.0 / 604.8 s 与 reviewer 在
  并发负载下的 697.6 s 不可当 A/B。

**R6.2 待做**：参数区间与例外集、`t` 的一般等价类（现在只有个案观测）、5,756 行
"任意 t 成立 / 仅特定 t 成立"的分档说明，以及覆盖说明里"已计算 / 有独立对照 /
仅内部一致 / 不支持"四种证据级别的逐条列举；另需处理 R6.1 记下的两处定义域边界
（大分母 fail-closed 上限、`t = 0` 一类形式值是否应改为显式错误）。

### R6.2 更正与第 1 步：translation-aware ledger（2026-09-23）

外部复核指出上一轮的两条过强假设，均已复现并处理，**结论按三阶段记账**：

1. **R6.1 的 `t ≡ t+1` 假设被撤销**：冻结方向 `v` 是母群倒格矢，而线小群含分数平移，
   `exp(2πi v·T)` 在这些群上非平凡，因此 `t` 与 `t+1` 是否同一 parent irrep 需要
   monodromy 语义（`(k+G, M_G(α)) ~ (k, α)`），不能默认同 label。`w_parameter_shift`
   门禁目前实际隐含 `M ≡ 1`，属过强，待改写（未改代码）。
2. **`efe8abb` 的 223/187 降级**：translation-aware ledger 显示 4 个"可算"见证
   （11329/12306/12307/12311）在 `3/4` 上仍满足共轭律、投影范数为整数 ⇒ 它们是
   **合法分支/gauge 差异**，不是错误。故 223/187 不再作为独立错误计数。
3. **唯一硬失败 = 58 处非整数重数**：`examples/line_transport_ledger.rs`（新 tracked
   example）在完整子群枚举（子群自身 data-Hall 操作映射进母群帧，保留精确平移）、
   精确元素共轭、`(1/|H|)Σ|χ|²` 投影范数、Bloch 协变三项上给出：
   * anchor `1/4`：11 个见证 **全部通过**（共轭律 0 违反 / 1296 对，范数整数如 14）；
   * `3/4`：4 个可算见证通过；7 个非整数见证**违反共轭律**（14429/14430/14503/14723
     各 8 处、14460/11360 各 72 处、14453 144 处，偏差到 4.0）
   ⇒ 送进投影的对象在这些参数上**不是子群 H 的表示** ⇒ **transport/assembly 缺陷**，
   与 seed 共轭选择、multiplier 选择都无关。修复验收：任意参数共轭律违反 0 +
   投影重数非负整数 + anchor 5,756/5,756 与复表伙伴 oracle 不变。
   常驻不变量：`the_anchor_restriction_satisfies_the_space_group_conjugacy_law`
   （只钉 anchor，未把任何过强等式写进测试；naive "lattice invariance" 与
   "seed multiplicativity" 都**不**成立/不作为律，见 example 文档）。
   复现：`line_transport_ledger --witnesses`、`line_transport_ledger 14453 SM1`。

### R6.2：参数族覆盖说明（2026-09-23）

交付：`docs/subduction-r6-coverage.md`（正式报告）、`docs/subduction-conventions.md`
§16 的"能力 B"小节、新 example `examples/line_family_coverage.rs`（`--gate` 可进 CI，
`--output` 出逐行 TSV）、审计新增 `w_parameter_shift` 硬门禁。全部数字由 example 重算：

* **命题 1（支持集 = 四分之一网格，E2）**：单臂的条件是三个子群格基矢条件的**交**，
  支持集是各臂点集的**并**；实测 5,756/5,756 行的支持集恒为 `{0, 1/4, 1/2, 3/4}`，
  且在 103,608 个网格外 (行, 参数) 组合里折到子群 Γ 的臂数**全为 0**。
  实现陷阱：把跨臂的并集写成"生成元取 lcm"会算成交集（曾在本轮自测中给出
  `{2: 2192, 4: 3564}` 的错误直方图，被 example 的测试当场抓住），
  example 现在逐臂枚举点集后合并。
* **命题 2（网格外恒为 0，E2）**：`content(t) > 0 ⇒ 有臂折到 Γ`，故 `t ∉ (1/4)Z`（非
  退化点）时 `content(t) = 0`；引擎在 1/25 抽样行上全部为 0（4 行 fail-closed）。
* **命题 3（共轭 oracle，E2；按冻结表实/复分两支）**：73/73 个 direction ∈ `G*_parent`，
  故 `k(-1/4) = -k(1/4)`。表**实**时 `D* = D` ⇒ `χ_{-k} = χ_k*`，oracle 是**同一源**的
  pinned；表**复**（SG 209/210 `DT3`/`DT4`，188 行）时共轭换源，oracle 是同一记录里
  **伙伴源**的 pinned。**实测：实表 187 行重数不同 + 58 行非整数（245 实例，真实缺口）；
  复表 0/0（引擎正确，4 行未列伙伴源无 oracle）**。
  **记账更正**：R6.2 首版写的"199/58"是拿复表自己的 pinned 比出来的假缺口，已改为
  实表 187/58 + 复表 0/0（example 常量拆成 `CONJUGATE_GAP_*` 与 `COMPLEX_CONJUGATE_GAP_*`，
  §16 与覆盖报告同步更正）。这条更正同时否掉了外部建议的"把 seed 表示整体共轭"修法：
  对实表它是恒等操作（缺口全在实表上），对复表它会把已经正确的 188 行改成另一个 irrep。
  非整数重数说明送进内积的对象已不是一致表示。本轮已把 **character 层门禁**做进
  example（公开 `reconstruction()`，逐操作比较 `χ(3/4)(h)` vs `χ(1/4)(h)*`，容差
  1e-9）：实表 **5,510 行可比较、223 行不共轭**（最大偏差 4.0），其中 187 行表现为
  恒等重数不同、**36 行重数相同但字符已错** ⇒ 字符层严格强于重数层，缺陷定位在
  **transport 层**（`k → -k` 分支的 Bloch 相位/陪集代表元记账），不是 seed 的实/复
  也不是 multiplicity 求解器。修复验收：字符层 223 → 0、fail-closed 58 → 0，
  anchor（1/4、5/4）与复表 oracle 保持全绿。用户裁决仍是
  **如实记录 + `t ≡ 1/4 + n` 为已验证域**。
* **命题 4（等价参数，E2 + 门禁）**：`(t - 1/4)·v ∈ G*_parent` 时（含 `t = 5/4, 9/4`）
  与 `t = 1/4` **逐块逐目标**相同；审计新增 `w_parameter_shift: checked=5756
  mismatched=0` 并计入 `hard_failures`，配常驻负例
  `a_line_parameter_shift_mismatch_fails_under_every_flag_combination`。
* **形式值**：`t = 0`（2,918 行 ≠ pinned）与 `t = 1/2`（2,822 行 ≠ pinned、58 行
  fail-closed）是**另一个母群 irrep / 退化点**，返回形式值，不作为物理结论。
* **性能边界（写进报告）**：`t = 1/4` 走存储路径全表 11–14 s；一般参数大分母 +
  非平凡小余群单次可到 p99 ≈ 40 s、最大 53.6 s，所以全表逐点扫描不可行，
  R6.2 用"精确支持集 + 抽样佐证"，不假装逐点算过。
* 本轮还补了 reviewer C 的点：审计的 w 边界（只比恒等重数）已写进 §16 与审计报告；
  `line_arms` 的精确去重约定、`wave_vector()` 已规范化、`FullStarBlock::q()` 未规范化
  都写进 doc；`1,599/241` 标注为历史记录（实现已删，不可复算）；
  `NonIntegralLineFrequency` 死变体删除。
