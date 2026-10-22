# 经幺正子群构造精确磁嵌入（R36）

[公开入口](../src/irrep/magnetic_embedding.rs)是
`embed_via_unitary_subgroup(geometry, magnetic_set, ordinary_embedding)`。
它复用现有普通分导的幺正子群嵌入，不改变 R24 的五种候选 map 族。

## 构造与验证

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

## 普查接入与边界

[受限 Γ 普查](../examples/magnetic_gamma_subduction_census.rs)先保留 R24 原来选出的
exact 见证；只有原先没有 exact 见证时才尝试此路线，并按原有精确全序选一个代表。
选择只读取记录与几何，不读取待分导源、η 或分导结果；新增数量单独报告为
`records_fixed_via_unitary_bridge`。准备或算术失败作为 `unitary_bridge` 错误保留。
`certified_absence` 仍只指 R24 的有限族，绝不扩写为所有设置的无解证明。
最新冻结统计见[账本](../AGENTS.md)，旧 R28 数字保留为历史。

本入口证明 grey 母群包含和记录格相等，**不证明特定 order-parameter 向量的 stabilizer**，
不保证唯一 basis/origin，也不提供完整 setting 分类。它没有改变任何表示域或进行非 Γ
相位输运；高维、一般 B/C、非 Γ/full-star 与 spinor 分导仍须各自的作用契约与算法。
