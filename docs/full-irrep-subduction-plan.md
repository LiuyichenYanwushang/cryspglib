# 完整 irrep 分导：供独立 DeepSeek CLI 顺序执行的任务卡

> 2026-09-22 更新：原任务 1–8 已落地，任务 9 已闭合恒等内容审计；R4 的三个批次
> （构造目标 / 一维投影 catalogue / 二维投影表）随后把普通离散标量的**完整分解**
> 做到 366,260/366,260 = 100%（三门口禁 exit 0，见
> [subduction-audit.md](subduction-audit.md) 的 R5 报告）。当前状态和后续 R6–R12
> 执行卡见 [subduction-next-milestones.md](subduction-next-milestones.md)。本文保留
> 历史设计，其中的旧基线（如"14,713 个 probe 缺完整分解"）和示例签名不代表当前
> 实现；不要从任务 0 重新执行。

本计划针对能直接读写仓库的独立 DeepSeek CLI。按任务卡顺序执行，每次只完成一张；
本文件是实施计划，不代表下面的引擎已经实现。基准日期：2026-09-21。

最终目标：给定凝聚 irrep 及其 isotropy 方向，确定嵌入母群的子群，再计算另一母群
表示在该子群中的**全部不可约成分与重数**；支持普通子群，并逐步扩展到磁子群。
运行时全部在 cryspglib 内完成，不调用官方 iso，不做数值序参量稳定子求解。

首个验收链：

```text
母群：221 Pm-3m
凝聚：GM4+，P1 / (a,0,0)
子群：83 P4/m，必须保留这个方向的具体嵌入
被分导的母群表示：GM3+
预期：83 的 GM1+ × 1 + GM2+ × 1
维数：2 = 1 + 1
```

此预期已经用官方打印的 P1 基与现有 typed character rows 做过独立临时探针；
正式实现必须加入永久端到端测试，不能把该答案硬编码进生产函数。

**执行纪律与完成标准**

- 每步先读指定入口和调用者，再修改；不要按文件名猜 API 是否存在。
- 原有用户改动必须保留。提交用显式文件清单，不使用 `git add -A`。
- 错误、数据缺失、歧义必须返回具体错误；空分解不能冒充计算成功。
- 原子提交按可验证的行为拆分；输出提交、变更文件、验证结果与剩余限制。
- 测试先证明能抓住错误，再用于验收。来源相同的数组互相比对不构成独立物理验证。
- 未完成某类输入时，公开返回 Unsupported/MissingData，并在覆盖报告中计数。
  不能把跳过行计为通过，也不能通过放宽误差或增加豁免来消除失败。
- 不为后续阶段提前添加大框架、泛型表示引擎或新的大型依赖。

建议里程碑：任务 0–5 得到首个完整分解；任务 6 完成 Γ 分解；任务 7–9 扩展普通
非 Γ / 多臂；任务 10–11 完成磁分导；任务 12 补齐正式方向输入与公开交付。

**现有可复用入口与限制**

| 入口 | 用途 | 必须保留的限制 |
|---|---|---|
| `src/irrep/isotropy.rs` | 选子群、几何、恒等列、其它波矢条目 | 方向标签和记录携带凝聚上下文；不是完整分解引擎 |
| `src/irrep/types.rs::CharacterRow` | 特征标与完整 Seitz 操作成对存储 | 使用其声明的 representation space 和 dimension |
| `ordinary_scalar_selected_arm_block_trace()` | 普通标量 selected-arm 行 | 非 Γ 时只在相应小群操作上解释；不是 full-star 行 |
| `compound_selected_arm_view()` | compound constituent / realification 元数据 | 不按 ML 字符串长度猜重复计数 |
| `spinor_selected_arm_view()` | 双群操作与 spinor 特征标 | 必须保留 SU(2) lift；放到专门扩展阶段 |
| `src/irrep/bridge.rs::canonical_hall_ops()` | 使用 `SG_DATA_HALL` 的操作入口 | 当前有 first-Hall fallback；新引擎要求严格来源，不能静默回退 |
| `src/irrep/wigner.rs` | Seitz 代数、相位、setting 候选、磁共表示帮助函数 | 分母、代表元、变换搜索范围各有边界 |
| `scripts/iso_irrep_exact.py` | 精确 source operations、k arms、irtranslations 与来源编号 | 当前解析模型不等于已可用的任意操作 full-star matrix evaluator |
| `scripts/iso_irrep_data_hall.py` 及 `scripts/data/iso_irrep_data_hall_v1.*` | 已冻结的 source / data-Hall 来源 | 优先使用正式来源映射，不另造猜测表 |
| `scripts/verify_isotropy_oracle.py` | 官方几何与方向 oracle | 目前只有 48 行；没有验证完整操作嵌入 |
| `src/operation_group.rs` | 群闭包、逆元、操作集验证 | 当前操作归约使用的晶胞必须与待验群一致 |

禁止把 `IrrepRecord::characters()` 或 `matrices()` 与 `pir_rotations()`、Hall 操作
直接 zip。它们的文档明确说明 legacy 存储没有这样的帧与相位保证。

**固定数学约定：任务 1 必须落实成文档和测试**

坐标一律为列向量。令真正已验证的、子群 conventional 帧到母群 conventional 帧
的变换为：

```text
x_G = T x_H + o
R_G = T R_H T^-1
t_G = T t_H + o - R_G o
k_H = T^T k_G
```

如果子群 conventional 基矢按行存为 B，则 T = B^T。最后一个式子由 Bloch 相位
配对不变性推导，不能混用实空间和倒空间变换。origin 不改变 k 的坐标变换，
但会改变操作平移代表元。

原表 W 给出母群 primitive 坐标中的子群 primitive 基。必须先证明该基与子群
规范 ISO/Hall 基的**有向对应**，才能把它用于 T；格相等或体积相等不足以证明
标签对应。需要时把额外 unimodular 变换 U 显式保存，再组合 source/Hall 变换。
不要未经验证便把 `P_H^-1 W P_G` 当作规范 B。

区分 L_G（母群平移格）和 L_H（子群平移格）。操作集是否属于 G 可以用母群格
判断，但 H 的闭包、去重必须使用 H 的平移格。归约丢掉的平移必须返回给调用者，
以便非 Γ 阶段恢复 Bloch 相位。

**任务 0：保存当前基线**

输入：当前工作区、`AGENTS.md`、`docs/isotropy-data-semantics.md`。

工作：

1. 检查 `git status` 和 diff；本计划编写时，上一轮永久测试与修复仍未提交。
2. 将这些既有改动单独保存为基线提交；不要和新引擎混在一个提交里。
3. 跑本文件末尾的基线命令，记录实际结果与当前生成文件哈希。

验收：lib 311 passed / 4 ignored、integration 84 passed、doctest 27 passed；
离线 oracle 测试 9 passed；真实 oracle 48 行 / 26 个描述串通过；严格 clippy
通过。已有 Cargo workspace manifest 警告单独记录。

当前生成文件 MD5：`fbe341a9cfb1bd815632e812daae950d`。它相对 `acb8f7b` 仅改
5 处文档注释，非注释行相同；这不是完整生成器已经重跑的声明。

停止条件：基线不符时先定位差异。不能通过更新预计计数来掩盖测试丢失。

**任务 1：确定表示语义与数据契约**

负责：新建 `docs/subduction-conventions.md`；暂不改数值算法。

工作：

1. 写明“凝聚 irrep”和“被分导 irrep”是两个独立输入。
2. 固定上述坐标约定、L_G / L_H、Γ / selected-arm / full-star 的含义。
3. 首版选择：输入为现有物理母群记录；结果按复不可约表示分解。对实物理表示，
   分解的是其复化；维数校验仍使用原实空间的维数。物理实表示的重新分组另行展示。
4. 区分稳定的 source identity 与显示标签。compound constituent 不能用一个拼接
   ML 标签冒充一个复 irrep；目标条目必须可以追溯到正式 CIR 来源。
5. 核查已有精确类型：`ExactSeitzOp` 的固定分母是 12；原表 origin 包括分母 16，
   变换过程中也可能产生新分母。决定局部 checked rational 表示，不改全库代数。
6. 列明不支持边界：首版只接收 Γ 标量请求；未知输入报错，不返回部分项。

验收：文档给出 T、o、k、一个非对称换基的手算例；列出所用 API 的实际签名、
帧来源、维数来源和精度。不得把 `canonical_hall_ops` 当作无 fallback 的严格 API。

停止条件：同一个字段被解释为两种维数或两种坐标帧时，先解决契约。

**任务 2：建立操作与嵌入的官方证据集**

负责：新增 `scripts/verify_isotropy_operations.py` 和对应小型离线解析测试；
精简 fixture 放在 `tests/data/isotropy/`，沿用已有 pinned archive。

工作：

1. 在固定 setting 下查询 `SHOW BASIS`、`SHOW ORIGIN`、`SHOW ELEMENTS`，先明确
   `VALUE DIRECTION <label>`，并固定操作标签格式。保存命令与原始输出来源。
2. 固定母群、子群 setting 各自的来源；不要把所有记录都假定为 OR1。
3. 将输出解析为精确 R、t、B、o，明确程序是在母群帧还是子群帧打印操作。
4. 每个查询独立进程；如并行，使用独立工作目录，避免 `iso.log` 互相覆盖。
5. 解析失败、进程异常、分页截断、缺行、重复操作均失败；旧 ML 拼写的修正必须
   来自正式 source identity 映射，不能模糊匹配标签。

首批用例：221 GM4+ P1/P2/P3；221 GM3+ P1/C1；225 GM4- C1/C2；16 R1；
167 GM3+ P1；139 M1- P1。包含 P/C/F/I/R、相同子群号的不同嵌入和非零 origin。

验收：fixture 可以精确重放；221 GM3+ 的两个相同 W/o 记录分别得到 16 与 8 个
点群商代表；221 GM4+ P1 的 H 有 8 个。普通 conventional 操作数与点群商阶数
不能混为一谈，尤其是中心化胞。

停止条件：无法解释一条操作的帧或平移时，保存最小反例，不猜补。

**任务 3：实现局部精确仿射变换与晶格归约**

负责：先放在新模块 `src/irrep/subduction.rs` 的私有实现；必要时再拆私有子模块。

工作：

1. 实现所需的最小有理数/3×3 运算，使用约分和 checked 算术；整数运算溢出报错。
   Python 证据工具用标准库 `Fraction` 作为独立对照。
2. 实现 T/o 的正向、逆向、组合与 Seitz 变换；验证 R_G 是否真为整数，不做强转。
3. 实现针对显式晶格的 membership / reduction，归约返回代表元和被移除的平移。
4. 新引擎严格读取 `SG_DATA_HALL`；缺失或加载失败直接报错。尽量复用操作读取底层，
   不为兼容方便静默转到 first Hall。

验收：非对称换基、负行列式、超胞、非零 origin、含 1/16 的有理输入均能往返；
奇异矩阵和溢出被拒绝；模 L_G 相同、模 L_H 不同的两个操作不能被误合并。

停止条件：必须把分母硬截到 12 或把平移取模后丢弃才能通过测试时，设计不成立。

**任务 4：构造并验证具体子群嵌入**

负责：`subduction.rs`、`tests/irrep_subduction.rs`；读取任务 2 的 fixture。

工作：

1. 设计私有构造的 `SubgroupEmbedding`，保留 parent/subgroup source setting、
   isotropy ordinal、T、o、L_H、操作集及代表元映射。
2. 用目标子群的规范操作，经经过验证的 T/o 搬入母群帧；同时验证格、完整 Seitz
   集包含、闭包、逆元、群号和官方 `SHOW ELEMENTS`。
3. 不能只凭 SG 识别或第一组成功的 setting 候选选标签；不同自同构可能给出同一
   操作群但置换 irrep 标签。规范选择必须以已记录的 setting/基方向为依据。
4. 先检验 W 是否足以恢复规范基。如果额外 U 或 setting 无法由现有数据唯一确定，
   生成并冻结精简的逐记录嵌入元数据。优先存 T/o/setting 或相对修正，避免复制
   整套大操作表；这是离线生成数据，运行时仍为纯 Rust。
5. 校验公开 `IsotropySubgroup` 的 ordinal 与母群上下文，不能信任可修改公开字段
   组成的任意伪造记录。

验收：任务 2 全部用例精确匹配；故意交换 C1/C2、转置 B、改 origin 或删除操作
必须失败；对无法规范消歧的记录返回明确错误。

停止条件：只验证“同一个格/同一个 SG”还不算完成；禁止直接采用候选列表第一项。

**任务 5：交付第一个完整 Γ 分解**

负责：`subduction.rs`、`src/irrep/mod.rs`、`tests/irrep_subduction.rs`。

建议先只提供一个公开入口，暂不同时开发三个包装 API：

```rust,ignore
let h = isotropy_subgroup_for_direction(
    221, "GM4+", IsotropyDirection::Label("P1")
)?;
let decomposition = subduce_irrep(&h, "GM3+")?;
```

工作：

1. 复用 typed character rows，按完整 Seitz 操作配对，不能只按数组顺序或旋转匹配。
2. 在 Γ、普通单位群情形，对 H/T_H 的完整代表集 K 做复特征标内积：
   `m_a = sum_h chi_parent(h) * conj(chi_a(h)) / |K|`。
3. 先证明候选目标行是完整、互异、复不可约的集合，或明确限定已证明的输入集合。
   `compound_selected_arm_view()` 的 constituent 可用于提供复目标行。
4. 检查虚部、非负整数重数、维数和，以及重建后的每个操作特征标。只有在数值误差
   界内才可取整；错误中保留原始内积、操作/表示标识与残差。
5. 输出至少含：凝聚子群/嵌入上下文、被分导母群 source identity、其表示空间与
   维数、每个目标 source identity/ML/BC/k/维数/重数。没有来源的标签不要合成。

验收：黄金用例给出 GM1+ 与 GM2+ 各一次；2=1+1；其余项为零；逐操作重建成立；
更改被分导母群表示后，结果由字符决定，不能依赖黄金用例硬编码。

停止条件：结果只给恒等项、漏了其它项或只校验维数，均不是“完整分解”。

**任务 6：扩展 Γ 与 compound 的正确计数**

负责：表示行适配与 Γ 扫描测试，不改磁算法。

工作：

1. 按 `CompoundMetadata` 区分 `ConjugateRealification` 与 `DistinctComponentSum`。
2. 使用稳定 CIR 来源身份展开输入与目标，保留复化语义。目标 Gram 矩阵应符合复
   不可约特征标正交性；不能把 norm=2 的 compound 行套进普通复 irrep 内积。
3. 扫描可验证的 Γ 嵌入，比较所有目标项；按规范来源做去重，不按显示标签去重。
4. 对完整复 Γ irrep 集验证 Frobenius reciprocity：
   `sum_D dim(D) * mult(trivial_H, D|H) = [P_G : P_H]`。
   该式使用 Γ 点群范围；不把 ML 拼接分量数作为未经证明的通用除数。

验收：非零项、零项、Gram 矩阵、维数和、逐操作重建均通过；对交接所述的
1895 条 Γ Frobenius 记录重新清点适用范围，清点结果与历史数字不同必须解释。

停止条件：通过除以 2、除以标签分量数或删 compound 行才能得到整数时，先修语义。

**任务 7：实现 k 折叠与非 Γ 的平移相位**

负责：`subduction.rs` 的 k/相位适配；必要的 source 数据生成扩展独立提交。

工作：

1. 从已验证 T 计算 `k_H = T^T k_G`；按子群倒格等价类匹配，不能任意逐坐标 mod 1
   后就当作 conventional 中心化胞的完整等价判据。
2. 初期限制母群 star 只有一臂的非 Γ 情形；若做 selected-arm 子程序，名称和
   输出语义必须明确，不能把它作为 full-star 公开结果。
3. 完整保存母群操作代表元与映射平移，用与已有 typed 数据一致的 Bloch 相位符号。
4. 在同一折叠 k、同一 projective factor system 与同一操作代表集中做内积。
5. 子群折叠 k 可能不在 `query::irreps_of(sg_H)` 的离散高对称点记录中；先枚举真实
   缺口。缺数据时报告 MissingIrrepData，不能选最近 k 点或把它当零分解。
6. 若覆盖需要参数化 k 的目标表示，复用 exact source 的 k 参数模型并实现经验证的
   参数代入与操作求值；无法由当前 source 获得时，作为明确的数据/算法子任务处理。

验收：Bloch 配对在变换前后相同；带滑移/螺旋或中心化平移的用例对丢失相位敏感；
一个 Γ 结果保持不变；超出现有窄整数 KVector 范围时不能截断。

停止条件：浮点“最近匹配”、任意取整或忽略平移才能得到分解时，不进入下一步。

**任务 8：实现 full-star 限制与子群 k-star 分组**

负责：full-star 适配器及其测试；缺失数据通过生成器明确补齐。

工作：

1. 从有来源的母群 star arms 与 coset representatives 构建 full-star 特征标求值。
   优先复用小表示与诱导表示公式，避免生成每个群元的大稠密矩阵。
2. 对 h，只有被 h 固定的臂向 trace 贡献；计算各共轭小群元素时保存完整平移相位。
3. 将所有母群臂折叠到子群倒空间，按 H 的作用分成子群 stars。同一点的多条母群臂
   必须合并保留；不能只选一条，也不能重复计算整个 full-star。
4. 对每个子群 star 的代表 q，在其小群上分解该 q 的表示块，再恢复子群 full-star
   维数。维数和统一检查在 full-star 空间中进行。
5. 引入一个多臂折叠到同一点、一个折叠成多个子群 stars 的实际样例；候选由扫描
   选出后记录固定 ordinal/source identity，不凭标签猜样例。

验收：chi(E) 等于正确 full-star 维数；各输出项的
`multiplicity × child_little_dim × child_star_size` 之和等于母群 full-star 维数；
字符重建成立；臂枚举顺序改变不改变结果；去掉一条臂的负例失败。

停止条件：需要将 full-star trace 除以臂数来假造 selected-arm 行时，算法错误。

**任务 9：普通表全量验收与覆盖闭合**

负责：新增 `scripts/verify_irrep_subduction.py` 或一个 Rust 审计 example，复用生产 API。

工作：

1. 分别清点嵌入、母群 source 表示、目标 k/irrep 和完整计算结果；报告 total、
   passed、failed、unsupported、missing_data、ambiguous，禁止只报通过百分比。
2. 94271 条既有恒等分导逐条比较固定子空间重数。先审计相同源表示的多行、domain
   与 compound 语义，避免把不同记录当成可无条件相加的同一项。
3. 对适用母群表示中未列出的恒等项检查零值；只检查已列正项会漏掉假阳性。
4. 每次成功的完整分解都检查维数、非负整数、逐操作字符重建。使用独立数学检查
   补足“恒等列相同但非恒等标签错误”的盲区。
5. `other_wave_vector_subduction` 的 5756 条是附加单值条目，不是 spinor gate。
   只有其实际 k 参数已确定时才参与对应回归，否则单独报告参数缺失。
6. 约 450 条 setting 残余、旧 ML 拼写及官方空表必须分别处理。缺官方输出的记录
   可以用其它有来源的证据验证，但不能伪称已有官方逐行对照。

验收：已定义覆盖范围内全部正项、零项和独立恒等式通过；如果仍有未支持行，
发布范围明确写成部分覆盖。声称覆盖普通 pinned 表前，必须清零范围内未支持项。

停止条件：新引擎结果只有已知 94271 条正项正确，不能单独作为完整引擎验收结论。

**任务 10：单独建立磁子群嵌入**

负责：磁 embedding 分支与官方磁操作 fixture，暂不计算 corep 重数。

工作：

1. 从磁 isotropy 记录的 UNI 与 BNS setting 取得磁操作，保留每个操作的幺正/
   反幺正标志、T、o、母群帧和原始来源。
2. 不按普通表方向标签 join 磁表；不因普通表中存在同号空间群就借用其 basis。
3. 校验磁群闭包、逆元、unitary subgroup、antiunitary coset、UNI 与 setting。
4. 加入 Type-I、Type-III、Type-IV 的小型固定样例；Type-II 在该 isotropy 表中
   没有记录，明确报告这一数据范围，不伪造 230 条记录补覆盖。

验收：完整磁操作集与来源对应，反幺正标志变动的负例必须失败；适用时与既有
`ValidatedMagneticOperationSet` 和磁识别流程交叉核对。

停止条件：只能匹配空间旋转、不能匹配 time reversal 与平移时，不算磁嵌入完成。

**任务 11：实现磁共表示分导**

负责：独立的磁分导入口，复用 `corep.rs`、`wigner.rs` 中有来源的能力。

先明确输入：磁 H 不是普通单位母群 G 的普通子群，单个普通 G irrep 不足以定义
反幺正操作的作用。必须规定母群磁表示/grey-group 扩展、time-reversal action
及其平方；Landau 实序参量还应说明 time-even/time-odd 约定。

工作：

1. 先限制到磁 H 的 unitary subgroup，复用普通分导引擎得到完整成分。
2. 利用反幺正作用与 Wigner A/B/C 结构配对成分、检查维数、确定 corep 重数。
3. 反幺正矩阵的普通 trace 不是可直接套用的普通群字符，不能对所有磁操作照搬
   单位群内积。等价性和约束使用正确的共表示变换规则。
4. 缺少 intertwiner 或反幺正作用证据时返回明确错误，不能用零矩阵补齐。
5. 若同时要求 spinor 被分导表示，另加 SU(2) lift、中心元与 projective phase
   专门任务；不得将 `_w_subduce` 当作其数据来源。

验收：先完成一个 Γ、time-even/odd 约定明确的标量磁样例，再覆盖 A/B/C；
检查 unitary 限制、反幺正配对及 full dimension；只有通过相应证据的类别才开放。

停止条件：普通母群 irrep 的反幺正扩展未定义时，不给出看似唯一的磁分导答案。

**任务 12：补齐方向输入、公开 API 与交付文档**

负责：公开 wrapper、方向生成数据、formatter、README 和最终覆盖说明。

工作：

1. 官方 descriptor 按 `(母群 source irrep, ordinary/magnetic, direction label)`
   逐记录取得，覆盖 dim=2、dim>=4 与磁表。缓存按 pinned archive 与命令 setting
   校验；普通表与磁表分开取证，不能凭同标签互相推导。
2. 将官方 descriptor 与 legacy 内部记法显式区分。若保留旧别名，必须检查冲突；
   同一字符串若指向不同记录，应明确报歧义，不能按遍历顺序选一个。
3. 名称不存在或来源无法给出 descriptor 时，Label/Index 仍可用；不能用看似物理
   分量的合成串冒充官方输入。
4. 核心入口稳定后再加方向与整表 wrapper。方向 wrapper 必须区分两个 irrep：

```rust,ignore
subduce_irrep(&subgroup, probe_ml)
subduction_for_direction(parent_sg, condensing_ml, direction, probe_ml)
subduction_table(&subgroup)
```

5. 磁入口单独接收磁母群表示上下文，不能把普通函数名称不变地复用到不同数学对象。
6. 表格注明母群/子群 setting、复/实/共表示语义、k 与 star 维数；输出稳定排序，
   不把调试实现细节放进普通用户结果。
7. 冻结版本、来源哈希和覆盖报告。实际 profiling 后再缓存 embedding / character
   对齐；先验证正确性，再处理测得的热点。

验收：用户可以用官方 descriptor 走完整 golden path；CLI 示例、doctest、
错误示例、普通/磁覆盖说明与真实行为一致；没有意外的运行时文件或网络依赖。

**每张任务卡通用的验证与报告格式**

局部开发先跑能抓住该变化的测试；本步逻辑通过后，在提交前跑一次相关全套。
不要在同一份代码完全没变化时反复全扫。不要用 sibling Rustb 的失败替代本 crate 验证。

```bash
cd /home/liuyichen/TB_rs/cryspglib
export CARGO_TARGET_DIR=/home/liuyichen/TB_rs/cryspglib/target
cargo test --release -p cryspglib --tests
cargo test --release -p cryspglib --doc
cargo clippy -p cryspglib --all-targets --release -- -D warnings
python3 -m unittest discover -s scripts -p test_verify_isotropy_oracle.py
python3 scripts/verify_isotropy_oracle.py
git diff --check
```

新增脚本必须附可运行测试与具体运行命令。小型已固定的端到端样例进入普通测试；
昂贵全表扫描保留单独命令，覆盖扩展或数据生成变更时必须执行并记录覆盖分母。

每步交付：

```text
任务编号 / 提交：
本步已实现行为：
修改文件：
输入/输出/坐标/表示语义：
永久测试与负例：
实际运行命令及结果：
覆盖：total / passed / failed / unsupported / missing_data / ambiguous
剩余问题与最小反例：
下一步是否满足前置条件：
```

**可直接复制到 DeepSeek CLI 的单步提示词**

```text
你在 /home/liuyichen/TB_rs/cryspglib 中工作。
请读取 docs/full-irrep-subduction-plan.md，只执行任务 N。
先检查工作区和任务的前置条件，保留所有既有改动。
按任务卡的负责范围实现；遇到已存在的帮助函数，先核对语义再复用。
不能把 legacy characters()/matrices() 与 Hall 操作直接配对，不能静默回退 setting，
不能丢弃平移相位，不能把 selected-arm 结果冒充 full-star 结果。
为本步新增实际可失败的永久测试，并展示至少一个相关错误输入被拒绝的结果。
完成验收后按显式文件清单单独提交，并按计划末尾模板报告。
若前置证据不足，交付最小反例及缺失的数据，不猜测、放宽误差或跳过失败。
完成任务 N 后停止，不自动进入下一张任务卡。
```

第一轮填 N=0，随后填 N=1。不要一次要求模型完成全部任务；尤其不要在嵌入规范、
复/实表示和 full-star 语义尚未分别验收时并行修改这些相互依赖的模块。
