# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## 磁 symmetry 全覆盖实时账本（2026-07-31 起）

用户当前优先目标：先覆盖并修正全部 1–1651 UNI 的磁对称性，再继续上层
corep/summary 产品化。每次全扫、根因确认和修复后都必须更新本节；不能只在
对话中报告。

### 当前最终状态（更新至 2026-08-10）

磁 symmetry/corep 主线与第一轮 API 安全化已提交为
`e8610d3 fix: harden magnetic symmetry APIs`，本地 `main` 当前领先远端一个提交；
2026-08-10 的 Rust-native API 收口与 clippy 清零仍在工作树中，尚未提交或推送。
后续开发和回归应以下面这组最新结论为基线，而不是再沿用早期“部分磁群不支持”
或“只预览 6 个特征标”的判断。

| 层级 | 当前结果 |
|---|---|
| 磁数据库与群代数 | `UNI 1..=1651` 全覆盖；数据库内 `4479/4479` 个 UNI–Hall setting 均通过严格群闭包、逆元、恒等元、陪集及类型一致性检查。 |
| setting 识别与消歧 | 给定 parent/family Hall 时 `4479/4479` 精确回环；仅凭操作自动识别时有 `4461` 个唯一结果和 `18` 个显式 `MagneticUniAmbiguous`，歧义只涉及 `{275,282}`、`{277,284}` 两组，`UNI 283` 唯一。程序不得猜选候选，必须用 parent Hall/setting 上下文消歧。 |
| ISOTROPY/Hall 嵌入 | `1651/1651` 个 unitary subgroup、detected Hall 与 data Hall setting 已完成一致嵌入。 |
| 高对称点与 corep | 用户可按 UNI、BNS 或磁操作输入磁群，取得高对称点列表；选定高对称点后可得到 fixed-`k` magnetic little-group coreps、维数、Wigner 类型和特征标。全库审计共覆盖 `10,390` 个高对称点、`54,717` 个 coreps。 |
| 正式特征标表 | 已提供“每个磁操作一列”和“每个共轭类一列”两种 Markdown 正式表格；不再截断到前 6 个特征标，并附操作/列标签图例。入口为 `format_magnetic_character_table` 与 `format_magnetic_character_table_by_class`。 |
| 目标磁群回归 | `BNS 128.406` 与 `BNS 52.318` 已不再返回“计算不支持”；`128.406@Z` 稳定给出维数 `2,2,2,4` 的四个 coreps，正式操作表包含 `g1..g16` 全部 16 列。 |
| CIR 数据生成 | CIR 解析器支持复合反幺正矩阵：`11,202` 个原始 coreps 中复合反幺正项 `672` 个、拒绝 `0`；`8,388` 个可映射到磁数据库的 coreps 中未映射 `0`。 |
| 验证 | 全 `1651` UNI release summary 审计：成功 `1651`、失败 `0`、`10,390` k 点、`54,717` coreps；常规 release 测试：通过 `254`、失败 `0`、忽略 `3`，另有 doc-tests `26/26` 通过；`cargo clippy -p cryspglib --all-targets --release -- -D warnings` 为零警告。 |

必须保留以下语义边界：

- 上述 `18` 个 operation-only 歧义是缺少 setting 上下文时的真实不可判定性，
  不能简单归类为“官方识别错误”；有 parent/family Hall 后都能唯一确定。
- 当前 summary API 的对象是选定 `k` 点的 magnetic little group corep，不等同于
  full-star 共表示；若将来增加 full-star API，必须使用独立名称和输出语义。
- 个别 Type-A corep 在缺少可构造的 intertwiner/matrix 时，反幺正列会明确显示
  `antiunitary-pending(...)`，不能伪造为已完成的物理特征标；本轮指定的
  `128.406`、`52.318` 回归表不受此问题影响。

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
- **已修复（k-mesh）**：Rust 风格高层 API 与 crate-private `kpt_*`/`kgd_*` 内核
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
`msg_identify_with_parent_hall`，避免用不相容的单位立方晶格制造假失败。

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
    `coreps=54717`。gate 对每个结果验证 operation/character 列数、共轭类分割、
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
    pub characters: Vec<f64>,
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
seq sg "name" "label" dim kx ky kz kd
<16-int complex matrix line>  ← 每个 operation 一行
<16-int complex matrix line>  ← 第二行（和第一行一起编码复矩阵）
(re, im)                      ← complex character
...                            ← 重复（如果有更多 operation）
```

- 第一行 16 个整数是复矩阵的编码（ISOTROPY 特有格式）
- 第二行 16 个整数也是复矩阵的一部分
- 每个 operation 由 2 行矩阵 + 1 行字符组成
- CIR 无 irtranslation（字符不依赖 translation）

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
| `spacegroup.rs` | `Spacegroup`, `spa_search_spacegroup*` |
| `spg_database.rs` | `spgdb_get_spacegroup_operations`, `spgdb_get_spacegroup_type` |
| `magnetic_spacegroup.rs` | MSG identification: `msg_identify_magnetic_space_group_type` |
| `msg_database.rs` | `msgdb_get_magnetic_spacegroup_type` (1,651 UNI entries) |
| `msg_database_gen.rs` | Auto-generated MSG database (shipped as source) |
| `spin.rs` | Spin-polarized symmetry: `spn_get_operations_with_site_tensors` |
| `pointgroup.rs` | `ptg_get_pointgroup`, `ptg_get_transformation_matrix` |
| `primitive.rs` | `prm_get_primitive`, `prm_get_primitive_symmetry` |
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
| `irrep/settings_data.rs` | Hall→setting mappings. `include!()`-d into `generated_data.rs` |
| `irrep/wigner_extra.rs` | Pre-computed antiunitary character path. `include!()`-d into `wigner.rs` |
| `irrep/preamble.rs` | Generated data prelude |
| `irrep/{triclinic,monoclinic,orthorhombic,tetragonal,trigonal,hexagonal,cubic}.rs` | Per-crystal-system irrep data (`include!()`-d into `generated_data.rs`) |

---

## Test suite（2026-08-10 release 基线）

Core irrep diagnostics pass as of 2026-07-03.  The full spinor Wigner sweep reports
`spinor_complex_ok = 21216` with no `spinor_complex_fail`.

| Binary / Location | 当前结果 | Description |
|-------|-------|-------------|
| `src/lib.rs`（全部 unit modules） | `196 passed / 3 ignored` | Wigner、BCS、API、输入契约、setting 与磁群回归 |
| `tests/irrep_validation.rs` | `31 passed` | Full-sweep validation: every SG has irreps, dimensions match, labels well-formed, k-vectors positive |
| `tests/magnetic_integration.rs` | `16 passed` | Magnetic structure analysis and Result/error contracts end-to-end |
| `tests/magnetic_symmetry_coverage.rs` | `6 passed` | 1651 UNI / 4479 setting group algebra, round-trip and ambiguity policy |
| `tests/{cof3,crps4,la2nio4,bcs_corep_validation}.rs` | `5 passed` | Reference material cases |
| doc-tests | `26 passed` | Public Rust API examples compile and run |
| **常规 release 总计** | **`254 tests + 26 doctests passed; 3 ignored`** | 另行执行的 1651 summary audit 也通过 |

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
| 磁群 | `spn_get_operations_with_site_tensors` | `MagneticOpGenerationFailed`, `MagneticPrimitiveLatticeFailed` |
| 磁群 | `msg_identify_with_parent_hall` | `MagneticReferenceGroupFailed`, `MagneticFallbackReferenceFailed`, `MagneticUniMatchFailed` |
| 磁群 | `magnetic_dataset()` | 传播上游错误 |
| 非磁 | `prm_get_primitive` | `CellStandardizationFailed` |
| 非磁 | `sym_get_operation` | `SymmetryOperationSearchFailed` |
| 非磁 | `spa_search_spacegroup` | `SpacegroupSearchFailed` |
| 非磁 | `det_determine_all` | `SpacegroupSearchFailed` |
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
