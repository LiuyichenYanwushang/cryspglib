# R3：剩余缺口的来源分类（工作文档）

状态：**进行中**。本文记录 R3 的输入分母、工具、已复现的侦察结果与尚未钉死的
语义问题；分类表与逐操作对照在后续轮次补齐后才算交付。本卡不改生产求解算法。

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
  （校验 SHA-256 与内部一致性，不做 Hall 选择或相位修正）。PIR = 小群
  （little-group）不可约表示记录，CIR = 全表示记录。
- `scripts/generate_irrep_data.py`：pinned 生成表读取/生成。
- `examples/census_subduction_gaps.rs`：缺口清单（R2 后重跑）。
- Rust 侧 `strict_sg_hall_ops` / `Lattice` / `ExactSeitz`：与引擎同帧的精确操作与
  倒格运算（分类工具应与这些约定一致；不得另立一套帧）。

## 已复现的侦察结果（可复算）

1. **PIR 归档只索引参数化 k 域**：10,294 条 PIR 记录的 `k_arms` 全部带自由参数
   （自由参数个数 1/2/3 分别为 7,739 / 5,268 / 8,482；**没有** 0 参数的离散 k 臂）。
   离散 k 点上的小群表示走 CIR（11,202 条）。这与早期结论一致：`little_subduce_*`
   的块只在有自由参数的 k 域上出现。
2. **每个缺失星的点都能被某个 PIR 参数域精确覆盖**。抽样（精确有理线性求解
   `constant + Σ t_j·p_j = q`，未做任何"最近点"匹配）：

   | 子群 | PIR 记录 | 缺口组 | 命中组 |
   |---|---:|---:|---:|
   | #5 | 16 | 35 | **35** |
   | #2 | 17 | 24 | **24** |
   | #12 | 27 | 25 | **25** |
   | #8 | 15 | 40 | **40** |
   | #6 | 25 | 15 | **15** |

   即"缺表"不等于"缺数据"：归档里有参数化小群源穿过这些 q。
3. **命中记录的"小群操作集"一致**：把每条命中记录的 `operations` 作为签名分组，
   每组只有一个最大签名（一般位置记录的小群更小，会一并穿过同一个 q，必须按
   最大小群筛选，否则会把一般位置的表示当成目标——正是任务卡警告的陷阱）。

## 分类结果（2026-09-22，899 组全部有分类）

工具：`scripts/classify_subduction_gap_sources.py`（离线，只读归档，不改生产求解路径），
输入是 `examples/census_subduction_gaps.rs` 的全表缺口 manifest：

```bash
CARGO_TARGET_DIR=$PWD/target cargo run --release -p cryspglib \
  --example audit_irrep_subduction -- --require-complete --output target/audit.tsv
CARGO_TARGET_DIR=$PWD/target cargo run --release -p cryspglib \
  --example census_subduction_gaps -- target/audit.tsv > target/r12_gaps.tsv
python3 scripts/classify_subduction_gap_sources.py target/r12_gaps.tsv > target/r3_groups.tsv
python3 -m unittest discover -s scripts -p test_classify_subduction_gap_sources.py     # 15 项
R3_FULL_MANIFEST=1 python3 -m unittest discover -s scripts -p test_classify_subduction_gap_sources.py  # 899 组门禁
```

| 分类 | 组数 | 含义与后续路线 |
|---|---:|---|
| `analytic_general_position` | **222** | 小余群阶为 1（只有平移）：目标就是一维 Bloch 相位，与 R2 的子群 #1 同一条解析路线，不需要任何字符表 |
| `parameterized_source` | **676** | 归档 PIR 的参数化小群记录精确穿过该 q（含参数值 t）：数据存在，R4 只需把该参数下的记录materialize/解码，不需要新算法 |
| `special_value_no_source` | **1** | 小余群阶 > 1，但没有任何归档参数域穿过该 q：真缺口，必须走构造路线 |

其余统计：小余群阶分布 1/2/3/4/6 = 222/445/4/222/6；**322 组因子系统非平凡**、
442 组含非幺正（螺旋/滑移）操作；787 组存在"一般位置记录也穿过同一 q"的情形，已被
最大小群过滤剔除（不筛就会把一般位置的表示当成目标）。

唯一无源组：child **#155**（R32），见证 ordinal 11067、probe `W1`、
星 `(-1/4,-1/4,3/2; -1/4,1/2,3/2; 1/2,-1/4,3/2)`，实际小余群阶 2，只匹配到一般位置
记录 `GP1GQ1`（dim 12）——不得用它冒充目标。

### 关键语义（本轮钉死，勿再重新踩）

1. **PIR 只索引参数化 k 域**：10,294 条记录里离散点是"三个零方向"的退化参数域
   （如 `GM1` 常数 Γ、零方向）。域的维数由**方向向量**决定，不能数非空参数槽。
2. **`record.operations` 是整个空间群**（约化 centring 后的代表元，conventional 基），
   不是小群；小群必须自己按 q 从这些操作里筛。
3. **参数周期不等于 1**：k 域参数按**原胞倒格**周期化。C 心群的 `(0,1,0)` 方向要
   `t = 2` 才回到同一类——按 `[0,1]` 采样会把 `(0,3/2,0)` 这类点误判成无源。
   修正这一条后 `special_value_no_source` 由 42 降到 1。
4. **格归属必须用原胞格**：螺旋轴乘积与所选代表元之间会差一个 centring 矢量
   （SG 24 的 `(1/2,-1/2,1/2)`），用 `Z³` 判定会误报"乘积离开小群"。
5. **因子系统**按任务卡公式用同一条乘法表算：SG 24 的 P 点 `(1/2,1/2,1/2)` 小余群阶 4、
   3 个非幺正操作、ω ∈ {1/4, 3/4}（I 心 + 三个 2₁ 螺旋），可作螺旋/滑移见证；
   SG 5 的 U 线 `(0,t,1/2)` 阶 2、ω 全 1（C 心、对称操作），可作带心见证。

### 字符/矩阵可用性（逐组）

每个组的命中记录都带**完整的逐操作 token 槽**：`token_slots` 全部非空、
`token_missing = 0`（899/899 组）。归档 PIR 记录本身有 10,294 条，其中 5,517 条
token 槽全满，其余有空洞；但**缺口组的命中记录全部全满**，所以 676 个
`parameterized_source` 组的字据确实在 pinned 归档里。

限制：`scripts/iso_irrep_exact.py` 按设计**只校验、不物化**矩阵/字符 token，
仓库目前没有 PIR 物化器。所以"可用"的准确含义是：*数据在归档里，代入参数后的
物化（解码）是 R4 的数据工程任务，不是新的数学*。工具输出的
`matched_irtypes`（PIR 记录类型 1/2/3）与 `token_slots`/`token_missing` 列给出逐组依据。

### 逐操作对照（验收项）

- **带心**：SG 5（C2，C 心）U 线：小群 2 个操作（E、C2），归档记录 `U1UA1`/`U2UA2`
  在 t=1/2 精确命中，一般位置 `GP1GQ1` 被剔除（见单测）。
- **螺旋/滑移**：SG 24（I2₁2₁2₁，I 心 + 三条 2₁）：P 点小群 4 个操作、3 个非幺正、
  ω = 1/4 与 3/4（见单测 `test_screw_little_group_is_projective`）。
- **非对称换基**：**部分完成**。冻结表里有 392 条 shear（非 signed permutation）
  setting，但都不落在 899 个缺口组里（缺口组的 346 条冻结 setting 全是 signed
  permutation）。已钉住见证 ordinal 26（SG 3 → #3，`U = [[1,2,1],[-1,2,-1],[-1,0,1]]/2`）：
  冻结表的值与引擎 `trace_embedding 26` 打印的 `setting` 一致（建议用
  `cargo run --release -p cryspglib --example trace_embedding -- 26` 复现）。
  **尚未完成**的是把归档子群帧的小群操作经该 setting 变换到引擎子群帧的逐操作对照：
  实测用 U 直接共轭会得到分母为 4 的非整旋转，说明 U 不是两个子群帧之间的直接
  setting 变换（引擎的 `transform` 是整矩阵 `[[1,0,1],[0,1,0],[-1,0,1]]`，det 2，
  与 U 不是同一个对象）。这一条按"不猜约定"的原则留作 R3 收口前的最后一项，
  下一步应直接用引擎的 `SubgroupEmbedding::transform()` 与归档操作做逐操作映射，
  而不是继续试矩阵方向。

## 剩余工作

1. 非对称换基的逐操作对照（上节第三项，约定未钉死前不猜）。
2. PIR 物化器：把 676 个 `parameterized_source` 组在代入参数后的字符/矩阵真正解码出来
   （归档 token 已确认完整），并在有离散表可对照的 q 上交叉验证；这是 R4 的数据路线入口。
3. 单测默认套件约 41 s（归档加载 39 s）；全 899 组门禁由 `R3_FULL_MANIFEST=1` 显式开启。
