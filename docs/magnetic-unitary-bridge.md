# 精确磁嵌入的几何桥接（R36–R37）

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

## 普查接入与边界

[受限 Γ 普查](../examples/magnetic_gamma_subduction_census.rs)按以下顺序选择，并保留
每一旧路线已经选出的见证：R24 exact → 同方向的 H 桥接 → primitive 格归一化 →
同母群/凝聚 irrep 的其它普通方向。最后一步要求 H 号相同、primitive 记录格双向相等，
再调用公开 H 桥接逐一验证；按既有精确全序选一个代表，不能用同体积代替同格。
普通/磁方向标签不是同一编号系统，也不因找到其它方向的 H 就宣称两个 OP 向量等价。

选择只读取记录与几何，不读取待分导源、η 或分导结果。两种 H 来源合计报告为
`records_fixed_via_unitary_bridge`，primitive 来源单列 `records_fixed_via_primitive_lattice`。
R36 路线的准备/算术错误仍是 `unitary_bridge`，新增路线错误是 `additional_bridge`，
均保留诊断，而不归入无见证。
`certified_absence` 仍只指 R24 的有限族，绝不扩写为所有设置的无解证明。
最新冻结统计见[账本](../AGENTS.md)，旧 R28 数字保留为历史。

这些入口证明 grey 母群包含和记录格相等，**不证明特定 order-parameter 向量的 stabilizer**，
不保证唯一 basis/origin，也不提供完整 setting 分类。它没有改变任何表示域或进行非 Γ
相位输运；高维、一般 B/C、非 Γ/full-star 与 spinor 分导仍须各自的作用契约与算法。
