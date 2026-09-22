# R3：剩余缺口的来源候选清点

状态：**来源候选清点完成，复核门禁通过**。本卡交付 899 组的来源候选、矩阵块可用性与
几何见证，不改生产求解算法。参数求值、特殊参数点的适用性与目标表示完整性由 R4
验证；候选清点完成不等于新增完整分解覆盖。

## 输入分母（R2 后，2026-09-22）

| 项目 | 数量 |
|---|---:|
| probe | 12,878 |
| 缺失星 | 17,144 |
| (子群, k-star) 组合 | 756 |
| (子群, k-star, setting) 组合 | 899 |

逐子群的 probe / 缺失星 / 全部折叠星三个计数见
[subduction-gap-census.md](subduction-gap-census.md) 末节；行级 manifest 由
`examples/census_subduction_gaps.rs` 从全表审计 TSV 生成（列含 `ordinal`、
`probe_cdml`、`child_sg`、`canonical_q`、`status`、`setting_numerator`、
`setting_denominator`、`child_shift`）。

## 工具

- `scripts/iso_irrep_exact.py`：归档 ISO-IR PIR/CIR 源帧的严格加载器
  （校验 SHA-256 与内部一致性，不做 Hall 选择或相位修正）。**PIR 是物理不可约
  表示记录**（10,294 条：4,777 条离散 k + 5,517 条参数化 k 域），**CIR 是复不可约
  表示记录**（11,202 条：5,296 条离散 k + 5,906 条参数化 k 域）。两者都含空间群
  表示数据，不是按"小群 vs 全群"区分。
- `scripts/generate_irrep_data.py`：pinned 生成表读取/生成。
- `examples/census_subduction_gaps.rs`：缺口清单（R2 后重跑）。
- Rust 侧 `strict_sg_hall_ops` / `Lattice` / `ExactSeitz`：与引擎同帧的精确操作与
  倒格运算（分类工具应与这些约定一致；不得另立一套帧）。

## 已复现的侦察结果（可复算）

1. **PIR 归档同时含离散与参数化记录**：10,294 条 PIR 记录里 4,777 条是离散 k
   （参数方向全为零），5,517 条带非零参数方向。此前的 7,739/5,268/8,482 计数
   未正确区分零方向，已撤回；域维数应由方向向量的秩决定。
2. **899 个缺失组均找到 PIR 来源候选**。按归档 k 臂与原胞倒格精确匹配
   `constant + Σ t_j·p_j = q + G`，不使用最近点近似。以下子群的分组见证为：

   | 子群 | PIR 记录 | 缺口组（含 setting） | 命中组 |
   |---|---:|---:|---:|
   | #5 | 16 | 72 | **72** |
   | #2 | 17 | 40 | **40** |
   | #12 | 27 | 25 | **25** |
   | #8 | 15 | 41 | **41** |
   | #6 | 25 | 14 | **14** |

   这是归档域的坐标可达性证据，还不能代替表示求值。
3. **候选按非零参数方向数最少筛选**：`record.operations` 是整个空间群的代表元，
   不能用它直接比较各参数域的小群。分类器没有验证候选的一般点小群是否等于
   当前 q 的实际小群，也没有验证候选复成分是否完整；这些检查属于 R4 的输入验收。

## 分类结果（2026-09-22，899 组全部有分类）

工具：`scripts/classify_subduction_gap_sources.py`（离线，只读归档，不改生产求解路径），
输入是 `examples/census_subduction_gaps.rs` 的全表缺口 manifest：

```bash
CARGO_TARGET_DIR=$PWD/target cargo run --release -p cryspglib \
  --example audit_irrep_subduction -- --require-complete --output target/audit.tsv
CARGO_TARGET_DIR=$PWD/target cargo run --release -p cryspglib \
  --example census_subduction_gaps -- target/audit.tsv > target/r12_gaps.tsv
python3 scripts/classify_subduction_gap_sources.py target/r12_gaps.tsv > target/r3_groups.tsv
python3 -m unittest discover -s scripts -p test_classify_subduction_gap_sources.py     # 21 项，默认跳过全表项
R3_FULL_MANIFEST=1 python3 -m unittest discover -s scripts -p test_classify_subduction_gap_sources.py  # 899 组门禁
```

| 分类 | 组数 | 含义与后续路线 |
|---|---:|---|
| `analytic_general_position` | **215** | 小余群阶为 1（只有平移）：目标就是一维 Bloch 相位，与 R2 的子群 #1 同一条解析路线，不需要任何字符表 |
| `parameterized_source` | **684** | 归档 PIR 参数域精确命中且矩阵块可读；是来源候选，R4 尚须验证参数求值、适用性和目标完整性 |
| `special_value_no_source` | **0** | 无 |

其余统计：**322 组的当前代表元因子系统有非零相位项**；444 组含平移分量非整数的
所选代表元。输出保留旧列名 `nonsymmorphic_ops`，但该计数没有逐项识别螺旋/滑移，
也不是非幺正或反幺正操作的计数。811 组存在自由方向数更多的匹配，被候选筛选剔除；
这不是小群最大性验证。非零因子系统项也不自动证明其上同调类非平凡。
899/899 组的星内各臂给出同一个**小余群阶**（`star_orders == little_co_group_order`）。

覆盖率数字是复核修复后的重算值：修复前（含两处计算错误）是 222/676/1，71 组的命中
来源因此改变。修复详情见下节。

### 关键语义（本轮钉死，勿再重新踩）

1. **PIR 含离散与参数化记录**：离散记录是"三个零方向"的退化参数域（如 `GM1`
   常数 Γ、零方向）。域的维数由**方向向量的秩**决定，不能数非空参数槽。
2. **`record.operations` 是整个空间群**（约化 centring 后的代表元，conventional 基），
   不是小群；小群必须自己按 q 从这些操作里筛。
3. **参数周期不等于 1**：k 域参数按**原胞倒格**周期化。C 心群的 `(0,1,0)` 方向要
   `t = 2` 才回到同一类——按 `[0,1]` 采样会把 `(0,3/2,0)` 这类点误判成无源。
   历史上修正这一条后无源计数从 42 降到 1；修正耦合方向代入后剩余一组也已找到来源。
4. **格归属必须用原胞格**：螺旋轴乘积与所选代表元之间会差一个 centring 矢量
   （SG 24 的 `(1/2,-1/2,1/2)`），用 `Z³` 判定会误报"乘积离开小群"。
5. **帧契约（复核确认，勿再重复施加 U）**：冻结 setting 的 `U⁻¹` 已经进入引擎的
   嵌入矩阵 `T` 构造（`src/irrep/subduction.rs:1352`），折叠过程在那里算
   `q = Tᵀk`（`src/irrep/subduction_star.rs:607`），`canonical_q` 之后只做倒格约化；
   230 个空间群的运行时 Hall 选择都与冻结来源一致，归档帧到 data-Hall 帧的坐标变换
   为 `P=I, p=0`（来源契约 `scripts/iso_irrep_data_hall.py:7`）。因此
   `canonical_q` **已经**使用归档的坐标基，再对它应用一次 U 是重复变换。R4 仍须
   沿用操作代表元的格平移相位修正与 `child_shift` 回退，不能把坐标基一致理解为
   所有操作代表元的平移都逐位相同。
6. **因子系统**按任务卡公式用同一条乘法表算：SG 24 的 P 点 `(1/2,1/2,1/2)` 小余群阶 4、
   3 个带非格平移的螺旋操作；非零相位圈数 `φ=q·L mod 1` 为 1/4、3/4，
   对应 `ω=exp(2πiφ)=±i`（I 心 + 三个 2₁ 螺旋），可作螺旋/滑移见证；
   SG 5 的 U 线 `(0,t,1/2)` 阶 2、ω 全 1（C 心、对称操作），可作带心见证。

### 字符/矩阵可用性（逐组）

矩阵可用性用仓库**已有的 PIR 解码器**（`scripts/generate_irrep_data.py` 的
`_parse_pir_characters`）逐组核对：**899/899 组的命中记录矩阵块完整**，合计
90,624 个矩阵元；判定条件是解码器给出该 `(SG, label)` 且展平块长度等于
`dim² × 操作数`。输出列 `matrix_available` / `matrix_elements`。

**更正**：此前用 `irtranslations` 是否为 `None` 当作"矩阵空洞"是错的——它是参数化
相位字段，离散记录按格式本来就没有（例如 #5 Γ 的 `GM1`/`GM2` 四个槽全 `None`，
矩阵却完整）。现在这两列改名为 `irtranslation_slots` / `irtranslation_none`，
只描述该参数化相位字段本身。

限制：`scripts/iso_irrep_exact.py` 按设计只校验不物化；分类器改为直接调用仓库已有的
PIR 解码器判定完整性。因此"可用"的准确含义是：*数据在归档里且解码器能读出完整
矩阵块*；把某个参数值下的字符/矩阵求值接成分导目标，是 R4 的工作。

### 复核修复（第 4 轮，2026-09-22）

复核复现了两处计算错误与两处验证问题，均已修复并加永久回归：

1. **参数代入漏掉方向向量的非对角分量**（P1）：`arm_point` 原来只算
   `direction[axis]·t`，方向 `(1,1,0)`、`t=1/4` 会返回 `(1/4,0,0)`。现在按
   `k = constant + Σ_j t_j·p_j` 全分量求和。这直接推翻了"#155 唯一无源"的结论：
   归档 `Y1YA1`(7219)/`Y2YA2`(7220) 的耦合直线 `k=(t,t,3/2)` 在 `t=3/4` 给出
   `(3/4,3/4,3/2)`，与见证 `(-1/4,-1/4,3/2)` 相差 R 心倒格矢量 `(1,1,0)`。
   回归：`test_coupled_direction_contributes_to_every_component`、
   `test_the_reviewed_witness_is_a_parameterized_source`。
2. **旋转求逆少一次转置**（P1）：`rotation_inverse` 返回的是余子式矩阵除行列式，
   即 `R⁻ᵀ`；调用方再当作逆矩阵用，三方旋转上的倒空间作用因此错误。现在返回真正的
   `R⁻¹`（余子式矩阵转置后除行列式），`preserves_q` 的 `(R⁻¹)ᵀ q` 随之正确。
   7 个实际阶为 2 的组曾被标成阶 1 的解析目标（涉及 #155、#166、#167）。
   回归：`test_rotation_inverse_is_the_matrix_inverse`（对 6 个含三方旋转的空间群
   逐元素验证 `R·R⁻¹ = I`）。
3. **矩阵可用性判定错位**（P2）：见上节更正。
4. **Rust 因子系统测试把格矢消掉了**（P2）：`factor_system_turns` 先把乘积操作按格
   约化、再对平移差取余，真正携带 Bloch 相位的格矢因此丢失。现在乘积**不约化**，
   直接从 `s_i s_j` 与代表元之差取格矢并断言它属于子群格。
   回归：`a_screw_relation_produces_a_non_trivial_phase`（`S²=T(0,0,1)`、
   `q=(0,0,1/2)` 必须给出半圈相位，修复前返回全零）与既有的 SG 3 见证。

修复后 899 组重算：215 解析 / **684 来源候选** / 0 无源；811 组需要按自由方向数
筛选候选（原 787）；444 组含分数平移代表元（原 442）；星内各臂小余群阶全部一致
（`star_orders == little_co_group_order`）。

### 门禁修复（第 5 轮，`fc5fb0f` 复核后）

全表测试原来把 `matrix_elements` 当成 `star_orders`，注入 `star_orders=1,2` 仍通过。
现在小清单和全表测试都用 `csv.DictReader` 按列名读取，共用星阶校验函数：每行的
`star_orders` 必须等于 `little_co_group_order`，相同子群和星在不同 setting 下的阶
也必须相同。永久负例覆盖列重排后的 `1,2`、错误单值和空值，以及跨 setting 阶冲突。
显式启用 `R3_FULL_MANIFEST=1` 时若 manifest 缺失，测试失败，不再跳过后报告成功。

### 逐操作对照（验收项）

- **带心**：SG 5（C2，C 心）U 线：小群 2 个操作（E、C2），归档记录 `U1UA1`/`U2UA2`
  在 t=1/2 精确命中，一般位置 `GP1GQ1` 被剔除（见单测）。
- **螺旋/滑移**：SG 24（I2₁2₁2₁，I 心 + 三条 2₁）：P 点小群 4 个操作、其中 3 个带
  非格（螺旋）平移，非零相位圈数为 1/4、3/4（见单测 `test_screw_little_group_is_projective`）。
- **非对称换基**：**完成**（`tests/subduction_gap_sources.rs`）。冻结表里有 392 条
  shear（非 signed permutation）setting，都不落在 899 个缺口组里（缺口组的 346 条
  冻结 setting 全是 signed permutation），因此见证取冻结表中的第一条 shear：
  ordinal 26（SG 3 `A1` → #3，`U = [[1,2,1],[-1,2,-1],[-1,0,1]]/2`）。单测用引擎自己的
  路径逐个操作对照：`SubgroupEmbedding::transform().unmap_operation` 把该记录的每个
  母群代表元映到子群帧、模子群格约化后，**逐项等于归档里 child #3 的两个操作**
  （E 与绕 b 的 C2），且两个方向的包含都成立（不是子集）；在钉子星
  `(0,1/3,1/2)` 上小群 2 个操作、因子系统 4 个有序对的相位全为 0（SG 3 对称）。
  与离线分类器在该星上的结果（阶 2、ω 全 1）一致。
  复现：`cargo test --release -p cryspglib --test subduction_gap_sources`；
  引擎 side 的 setting 可用 `cargo run --release -p cryspglib --example trace_embedding -- 26` 打印。

## 结论与交付

R3 的输出（全部已核对）：

* 工具：`scripts/classify_subduction_gap_sources.py`（离线、只读归档、不改生产求解路径）。
* 报告：本文件；逐组清单 `target/r3_groups.tsv`（25 列，899 行；由 manifest 生成，
  manifest 由 `examples/census_subduction_gaps.rs` 从全表审计生成）。
* 测试：`scripts/test_classify_subduction_gap_sources.py` 21 项（含 `R3_FULL_MANIFEST=1`
  的 899 组门禁：分类分布、矩阵完整、星内阶一致，证明无静默漏项）；
  `tests/subduction_gap_sources.rs` 的非对称换基逐操作见证与螺旋相位负例。
* 三类结果：215 解析（小余群阶 1，Bloch 相位路线）/ 684 来源候选（归档 PIR 精确
  命中，含代入参数与完整矩阵块）/ **0 无源分类**；0 未分类。
  因子系统与代表元平移单列：322 组有非零相位项、444 组含分数平移代表元。

**R3 边界之外（下一张卡 R4 的入口，本卡不实现）**：把来源候选在代入参数后的
字符/矩阵真正求值成可用的分导目标，并在有离散表可对照的 q 上交叉验证；同时处理
操作代表元的格平移相位与既有的 `child_shift` 回退。**口径：候选数只是来源候选数，
不代表参数求值与目标表示完整性已验证**。

**R4 批次 1 与 2a 之后（2026-09-22）**：R3 判为 `analytic_general_position` 的 215 组
（小余群阶 1）由批次 1 现场构造 Bloch 相位关掉；R3 判为 `parameterized_source`、
小余群非平凡但**一维投影特征标可解**的那批由批次 2a 用精确 cocycle 关掉（不读归档
字符）。两者合计 `full_success 353,382 → 366,039`、`identity_only 12,878 → 221`、
0 错误。在新缺口上重跑本分类器：**58 组全部是 `parameterized_source`**（29 个子群、
星阶 4 的 52 组 + 星阶 6 的 6 组，矩阵块 58/58 完整、8,448 个矩阵元），它们需要
**二维**投影不可约表示，是批次 2b 的输入。批次状态见
[subduction-r4-batches.md](subduction-r4-batches.md)，新分母见
[subduction-gap-census.md](subduction-gap-census.md) 的“R4 批 2a 后”小节。

运行成本随矩阵解码和机器负载变化；本轮实际验证结果另记于 `CLAUDE.md`。
