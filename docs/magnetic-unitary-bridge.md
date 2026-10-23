# 精确磁嵌入的几何桥接与记录格障碍（R36–R42）

[公开入口](../src/irrep/magnetic_embedding.rs)是
`embed_via_unitary_subgroup(geometry, magnetic_set, ordinary_embedding)`。
它复用现有普通分导的幺正子群嵌入，不改变 R24 的五种候选 map 族。

## 经幺正子群构造（R36）

设普通嵌入的线性部分为 T，已有磁群识别结果的 M → H data-Hall 线性部分为 Q。
只提出一个线性候选 **A = T Q**，不枚举新的小矩阵族：

1. 检查磁 UNI、母群身份、磁操作的数据库 provenance 与群公理；
2. 普通嵌入的子群号必须等于独立识别的幺正子群 H；
3. Q 的浮点数只提供候选，必须能落在既有有理网格上；
4. 精确、双向检查 **A L_M = L_record**，包括整格平移与心矢；
5. 复用完整联立同余求解器，针对全部幺正及反幺正类重算共同原点；
6. 对每个离散原点分支的代表元重新执行完整磁嵌入验证。

不会直接沿用普通记录的原点或浮点 setting 的 shift。返回的
`ParentEmbedding` 标记为 `EmbeddingKind::UnitarySubgroup`，保存实际有理 map/shift。
一个分支有连续核时仍只返回一个见证，与既有 shift API 的语义一致。

空列表仅说明**这一个**派生 map 未取得 exact 见证；不证明别的 map 不存在。
身份/数据/帧不适用通过 `GeometryError::UnitaryBridgeUnavailable` 报告；
精确算术或 shift 枚举上限仍是错误/未决，不能转成空解集。
新增的两个 enum variant 意味着下游穷尽 match 需要增加分支。

## 使用方式与见证

调用者先用已有 isotropy 查询取得磁记录与普通 H 记录，再构造
`MagneticGeometry`、磁操作集与 `SubgroupEmbedding`，传给新入口。
常见选择是同一凝聚 irrep、同一方向标签的普通记录；没有同名记录或 H 身份不一致时，
不能假定该桥接适用。入口本身也允许调用者显式提供其他已验证的 H 嵌入。

SG 3 的 E1/P1 磁记录（UNI 24，H = SG 5）得到一个确切见证：

```text
A = [[0,0,1], [0,2,0], [-2,0,1]], shift = (0,0,0)
```

它在旧 R24 候选族之外：旧的“族内无成员”认证仍然为真，但这个 map 完整验证成功。
[常驻公开测试](../tests/magnetic_unitary_bridge.rs)还覆盖 SG48 的四分之一原点、
SG227 含心设置的八分之一原点、错误身份/原始操作顺序、错误格与算术上限。
对 SG3 的两个 Γ 源、两种 η、两种标签约定，共 8 个高层调用逐类检查非平凡字符
及复向量 intertwiner 重建；没有把凝聚 E1 误当成待分导的 Γ 源。

单独的全表桥接回归覆盖 16,721 条磁记录：1,482 条无同名普通方向，16 条 H 身份不同；
其余 15,223 条中，12,906 条取得完整 exact 见证，2,317 条派生 map 未通过。
这是**桥接路线自身**的几何统计，不是新增记录数或全磁分导成功率。

## Primitive 格归一化（R37）

另一个公开入口 `embed_via_primitive_lattices(geometry, magnetic_set)` 不需要普通记录。
令 B 的列为磁记录的 primitive 基矢（母群 conventional 坐标），C 的列为磁数据库参考胞的
primitive 基矢，提出 **A = B C⁻¹**。不能只用 B：磁操作使用 conventional 坐标，完整源格
可能含心。这里不枚举基矢置换或其它小矩阵族，也不把磁参考胞误当成 H 的标准 data-Hall 胞。

首先用原始操作重验 **L_M = C Z³** 的两个包含关系，再使用两条桥接共用的 exact 格、
联立原点与完整磁操作验收路径。参考胞约定与全部 **1,651 个 UNI** 的源平移格相符；
这只是源格约定的全库核对，**不是 1,651 个磁群都能嵌入参考母群**的声明。

[公开回归](../tests/magnetic_primitive_bridge.rs)钉住 SG21 T1/P1 的含心格归一化见证、
T3/P1 必需的四分之一原点、SG19 Y1/P3 的“该 map 无解但另一 map 可嵌入”反例、
身份/原始数据/元数据拒绝、周期缩放与算术上限。SG21 的四个 Γ 源 × 两 η × 两标签，
共 16 个高层调用逐类重建非平凡作用。单独全表回归有 **11,766/16,721** 条取得几何见证，
其余返回空列表；这个数不是在旧路线之外新增的记录数。

新见证使用 `EmbeddingKind::PrimitiveLattice`；源数据或参考格不可验证时返回
`GeometryError::PrimitiveBridgeUnavailable`，算术/容量错误仍保留。穷尽 match 需补分支。
空结果只针对这一种基对应；没有一般 setting 或 OP stabilizer 的否证。

## 字面记录格的无解证书（R38）

公开入口 `record_lattice_obstruction(geometry, magnetic_set)` 不构造候选 map。
它复用 R37 的 provenance 与完整源格校验，检查是否有某个源旋转，在所有保持记录格的
母群旋转中找不到相同的 trace、determinant 与模二秩。

证明：令 C、B 分别为完整源格和给定记录格的 primitive 列矩阵。若存在 exact map A，
则 U = B⁻¹ A C ∈ GL(3,Z)。对应的旋转必须满足

**B⁻¹ S B = U (C⁻¹ R C) U⁻¹**。

因此在 F₂ 上，R−I 的秩相同。母群的完整旋转集合是有限的；若某个源旋转的秩在所有
同 trace/determinant 且保持记录格的母群旋转中都不出现，就排除了**对此字面格的所有
exact map**，而不只是某个候选族。原点平移不改变这项障碍。

- 返回 `Some(RecordLatticeObstruction)`：证书保存源旋转、源模二秩和全部兼容母群秩。
- 返回 `None`：只有“这个必要条件没有发现障碍”；不是存在性证明，也没有解联立旋转、
  平移或原点条件。计算/来源错误仍是错误，不转成 `None` 或证书。
- SG3 A1/P1、UNI11 的源秩为 0，母群可用秩只有 1；D1/P1、UNI24 则反过来。
  A1 仍能以 A=2I、零 shift 作非 exact 包含；两 Γ 源、两 η 的高层分导均通过。
  因此这**不否证磁子群本身**，也不拒绝不同周期格的 `ContainmentOnly` 调用。

[公开回归](../tests/magnetic_lattice_obstruction.rs)遍历全部 16,721 条记录，得到 **56 条**
字面格障碍（母群均在 SG3–15）。另有完整 GL(3,F₂) 的 168 个元素直接检验两个反例的
共轭方程、168 种整数 unimodular 重选基控制、所有 512 个二进制矩阵的独立行空间计数，
以及非保格旋转、来源与算术错误控制。没有枚举新的物理 map 候选族。

这 56 条应继续核对原始表与官方 cell/axis convention，不能靠更大的 map 搜索消除。
R37 的 1,134 条缺口据此细分为 56 条字面格障碍与 **1,078 条仍未解决**；没有增加正见证，
也没有把这些字面格结论扩写成其它 convention 或任意周期格的不存在性证明。

## 反幺正平移与成对辅助群（R39）

`embed_via_antitranslations(geometry, set, ordinary_embedding)` 仅以普通 H 嵌入的
**旋转子群**为候选，不再要求其校准周期与字面记录格相同。仍先检验母群/H/UNI 身份、
数据库 provenance、群公理与完整源格 C；不适用的 Type I/II/III 返回空搜索。

在源与目标 primitive 坐标中，分别以 H 的旋转构造**去掉 H 平移的辅助群**，保留二阶
反幺正平移。目标遍历 (½L)/L 的七个非零类，只保留落入母群格且被 H 模 L 固定的类；
这不是枚举小矩阵。现成识别器将两辅助群放入同一无 Hall 限制的 canonical frame，
其线性变换 D_s、D_t 提议 **A = B D_t⁻¹ D_s C⁻¹**。不能用实际 M 的 Hall 作为辅助群
的强制身份：去掉螺旋/滑移平移后可能变成另一空间群。

辅助群只用于提议线性 map。验收使用**原始 M 的全部平移与操作**，重算完整格双向相等、
全部共同原点分支和包含，因此非 symmorphic H 也可使用，不会把螺旋/滑移条件丢掉。
不采用浮点原点；浮点标签和矩阵都不是接受证据，也不声称穷尽 setting/normalizer。

返回 `AntitranslationSearch`：`embeddings` 均为 exact 见证，kind 为
`Antitranslation`；`unresolved` 逐个保存识别失败的半格平移和原始 `SymError`。
即使其它候选成功，也不清空这些未决项。准备/精确算术失败仍返回错误，不能变成无解。
新增 `AntitranslationBridgeUnavailable`；调用者的穷尽 match 需更新。

[公开回归](../tests/magnetic_antitranslation_bridge.rs)包含 SG3 B1/C1 与 SG4 B1、
完整原点分支、168 种重选基、原始 screw 不能换成辅助群、身份/算术/范围负例、
保留与成功并存的识别失败，以及 48 个 source/η/标签/嵌入的高层逐操作复向量重建。
全 16,721 条记录各取**首个 H 号匹配普通行**，R42 选轴修复后为 11,539 条有见证、
5,170 条空结果、12 条无对应 H 行；准备失败 0、候选未决 114。**11,539 不是新增覆盖数**，
也不是所有 H 方向的穷尽结果。R41 对照为 11,525 /5,184 /130，旧有成功记录全部保留；
新增 14 条首个 H 行见证，所有返回的 map/origin 条目均回验原始 M。
普查只在全部旧路线失败后使用，新增数与最新冻结验收见账本。

## 分级点扩张（R40）

[`embed_via_point_extensions(geometry, set, ordinary_embedding)`](<../src/irrep/magnetic_embedding/point_extension.rs>)
处理 Type III，仍验证 UNI/母群/H 身份、原始数据和完整源格 C。
Type I/II/IV 不在此路线内；非母格内周期或 H 不保记录格也返回空搜索。

在记录 primitive 基 B 中枚举**所有保格母群旋转** a。对于给定普通 H 的点子群，
只接纳 a∉H、a²∈H 且 aH=Ha 的情形：这正是 H∪aH 成为分级点群的条件。
按完整 aH 去重，不是枚举任意小矩阵；四阶反幺正生成元也允许，平方只须属于 H。
源与目标辅助群都去掉平移但保留分级；识别度量使用**整个分级点群**的旋转，
而非仅 H。复用现成群验证、无 Hall 限制的识别、有限有理帧与
**A=B D_t⁻¹ D_s C⁻¹**。随后重验**原始 M**的全格、全部共同原点分支和所有操作，
不以 symmorphic 辅助群替代实际螺旋/滑移。

返回 `PointExtensionSearch`：exact 见证的 kind 为 `PointExtension`，
`unresolved` 保存每个未识别陪集的代表旋转（记录 primitive 坐标）与 `SymError`；
其它候选成功不会抹掉这些项。准备错误为 `PointExtensionUnavailable`，精确算术错误
仍传播；新增 enum variant 需更新穷尽 match。没有穷尽识别器的 normalizer，空结果不证无解。

[公开回归](<../tests/magnetic_point_extension_bridge.rs>)覆盖 SG79/87/88 含心记录、
SG88 必需非零原点、SG4 的原始反幺正 screw、身份/算术/范围控制，以及 16 个
source/η/标签高层调用的逐操作复向量重建。模块另以独立的有色乘法表直接检查陪集条件。
全 16,721 条记录各取**首个 H 号匹配普通行**，得到 2,872 条有见证、13,837 条空结果、
12 条无 H 行，准备错误/候选未决均 0。**2,872 不是新增覆盖数**。

**R40 识别器边界的历史实证**：SG79 GM3GM4 的一个已验证 map 在
全部 168 种二进制 unimodular 重选基后仍逐一通过 exact 检查；重新运行此提议器则
132 种成功、36 种保留 `MagneticUniMatchFailed`，没有静默空结果。比如
U=[[1,1,1],[0,1,0],[1,0,0]] 下反幺正旋转为 [[0,1,0],[1,0,0],[-2,-2,-1]]。
不能把这些识别失败误报成不存在，也不以“每个重选基都重新找到见证”为本路线承诺。
R41 在共享定心函数修复了这一具体根因，见下节；原 36 个失败现全部恢复，168 种都重新得到 exact 见证且无未决。这不证明一般识别器/normalizer 已完备。

## 指数二的精确定心（R41）

旧普通/磁识别共享的 [`get_base_center`](<../src/spacegroup.rs>) 以单位坐标行或行绝对值和
猜测心型。这不保持 primitive 重选基不变性：对
T=[[0,−1,1],[0,−1,−1],[1,2,0]]，det T=2，旧函数退回 Primitive，
导致合法 invariant 度量下的 C2 被普通接口错认成 SG3，磁接口返回未决；同一群经
unimodular 变换后却可识别为 SG5 /UNI22。不是扩大候选族就能修复的问题。

T 的列是 conventional 基在 primitive 坐标中的表示，所以 conventional 坐标内完整格为
T⁻¹ Z³。det T=±2 时，格商 T⁻¹ Z³/Z³ 恰有两个元素，非零类必为 v/2，且
**T v = 0 mod 2**。遍历七个非零二进制 v 得到唯一心矢，分别确定 A/B/C/I；
左乘任意整数 unimodular U 不改变这个核。单轴半平移意味着轴未 primitive 化，
返回 Error，不再静默丢掉半格平移。逐元素先取奇偶，求和不受 i32 极值影响。
原有 A/B/I→C 校正与指数 1/3/4 分支保持不变，无新增 API、依赖或设置搜索。

常驻回归以**有理矩阵逆的格商**独立核对全部 3,792 个 {−1,0,1}、det=±2 的矩阵，
对四种心型各施加 168 个二进制 unimodular 左变换并核验原有校正，还直接比较
普通/磁公开识别在两种等价输入坐标下的 SG5 /UNI22 结果。更上层的原始 M、完整格与
共同原点验收没有删减；两条辅助群路线的实测覆盖、残余未决与旧输出比较见最新账本。

## 选轴成功必须有完整基（R42）

共享 [pointgroup](<../src/pointgroup.rs>) 的单斜选轴曾把 norm²≥8 的已有轴丢弃。
对 R=[[3,−1,1],[4,−2,1],[−4,1,−2]]，目录中的负特征轴 norm² 分别为 2、11、11，
第二个独立轴因此找不到；调用层又丢弃失败标志，以默认坐标轴补出了错误矩阵。
修复后搜索完整的**原有 73 轴目录**，确认主轴存在，传播单斜/层群/单轴/多轴分支失败，
多轴分支要求收齐三个轴，拒绝不合法层平面和奇异矩阵；不扩目录、不猜 fallback。
原实例普通识别得到 SG5，带半格反平移的磁辅助群得到 UNI23。

五项常驻回归包含上述公开识别、530 Hall 的整非奇异点群帧，以及 C2/D2/C3/C4/立方和
层群的拒绝控制。168 个二进制 unimodular 重选基以独立整数特征空间检查目录可用性：
**138 个有效帧、30 个目录不足而返回 None**。这是有界目录的准确结果，不是一般选轴完备性。
超出目录的合法群仍可未决；尤其不能把已知点群号加奇异 frame 当成识别成功。

R42 的 Γ 普查新增 10 条 fixed 记录和 144 个成功源/η 行；还有一条已成功的
SG20 T1 /UNI23 记录选中了新出现、按既有全序更小的嵌入。其 8 行 map/origin 改变，
其中 4 行目标 GM1/GM2 交换。**这不是同一嵌入的 gauge 换名**：两种嵌入可把 H 放在母群
不同位置，限制的作用因而不同。旧、新两组明确 map 都通过原始 M/全格与逐操作重建，
高层 API 仍接受旧 map；比较固定嵌入的结果应显式传入它。此处不声称已确定该 OP 的稳定子。

## 普查接入与边界

[受限 Γ 普查](../examples/magnetic_gamma_subduction_census.rs)按以下顺序采用首个成功路线，
不因后续路线另有见证而替换它；同一路线增加候选后，全序最小见证可能改变：R24 exact → 同方向的 H 桥接 → primitive 格归一化 →
同母群/凝聚 irrep 的其它普通方向 → R39 反幺正平移 → R40 分级点扩张。R37 的其它方向桥接要求 H 号相同、primitive 记录格双向相等，
再调用公开 H 桥接逐一验证；按既有精确全序选一个代表，不能用同体积代替同格。
R39/R40 则只取同 irrep、H 号匹配的旋转子群（按旋转集合去重），不复用该普通行的格；
所有候选仍按同一全序选择。普通/磁方向标签不是同一编号系统，也不因找到其它方向的 H
就宣称两个 OP 向量等价。

选择只读取记录与几何，不读取待分导源、η 或分导结果。两种 H 来源合计报告为
`records_fixed_via_unitary_bridge`，primitive 来源单列 `records_fixed_via_primitive_lattice`。
R36 路线的准备/算术错误仍是 `unitary_bridge`，新增路线错误是 `additional_bridge`，
均保留诊断，而不归入无见证。R39 的准备/精确算术错误为 `antitranslation_bridge`；
新增见证计入 `records_fixed_via_antitranslation`。候选识别未决单列为
`antitranslation_unresolved_proposals`，逐项打印，不假装候选不存在；旧 R24 族内证书
不因此扩张。已有成功不会因为其它浮点提议失败而被撤销。R40 的对应计数为
`records_fixed_via_point_extension` 与 `point_extension_unresolved_proposals`，准备错误为
`point_extension_bridge`；两条路线共用普通 H 点子群去重与原有的确定性选择序。
`certified_absence` 仍只指 R24 的有限族，绝不扩写为所有设置的无解证明。
R38 的独立格审计在选取见证之后运行，不改变选择或 TSV。证书计数单列为
`records_literal_lattice_obstructed`；`lattice_obstruction_errors` 包含计算失败和证书/正见证
冲突，非零会使 gate 失败。最新冻结统计见[账本](../AGENTS.md)，旧 R28 数字保留为历史。

这些桥接入口证明 grey 母群包含和记录格相等，**不证明特定 order-parameter 向量的 stabilizer**，
不保证唯一 basis/origin，也不提供完整 setting 分类。它没有改变任何表示域或进行非 Γ
相位输运；高维、一般 B/C、非 Γ/full-star 与 spinor 分导仍须各自的作用契约与算法。
