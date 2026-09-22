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

## 尚未钉死、下一轮必须先解决

1. **维数语义**：PIR 记录的 `dimension` 是物理（full-star）维数还是小群维数？
   抽样中最大签名的小群只有 2–4 个操作，却出现 `dimension` 4–6 的记录，
   Σdim² 远大于小群阶，说明不能直接用 `dimension` 当小群维数。需要：
   `k_arms` 的星大小（star size）与 `little_irr_full_dim` 的关系逐条核对，
   并用一个手算见证（例如 #5 的某条 2 操作小群）钉死。
2. **完整性与特殊参数点**：对每个 q，(a) 实际小群是否严格大于命中记录的
   *generic* 小群（特殊参数点，小群增大）；(b) 最大小群记录的
   Σ(小群维数)² 是否等于小群（有限模型）阶。两者都通过才可标
   `parameterized_source_complete`，否则标 `special_parameter`。
3. **帧核对**：PIR k 臂的坐标基与缺口 `canonical_q`（子群 primitive 倒格基）必须
   逐条对齐；抽样是精确命中，但要以"同一 q 的两条独立来源"（归档 k 臂 vs 引擎
   折叠 q）显式记录，不能只靠数值巧合。
4. **因子系统**：按任务卡公式 `s_i s_j = T_{L_ij} s_k ⇒ D(s_i)D(s_j) =
   exp(+2πi q·L_ij) D(s_k)`，用子群自身操作建有限小余群乘法表并输出 ω_ij；
   `ω ≠ 1` 的组单独标 nonsymmorphic/projective（螺旋/滑移），不能与 symmorphic
   组共用结论。

## 计划中的分类输出（每组的必需字段）

`child_sg`、`canonical_q`、star size / arm count / block dimension / parent dimension、
实际小群操作与代表元、有限小余群结构（阶、abelian 与否）、因子系统（ω 值集合、
是否平凡、nonsymmorphic 标志）、命中的 PIR irnumber 与维数、完整性判定
（Σdim² 与小群阶）、特殊参数点标志、可代入参数 t（精确有理）、字符/矩阵可用性、
分类标签（`analytic_constructible` / `parameterized_source_available` /
`needs_new_source_or_algorithm`）、最小见证（ordinal + probe + q）。

验收附加项（任务卡）：至少一个带心、一个非对称换基、一个螺旋/滑移样例完成
逐操作对照。为此分类工具必须复用引擎同帧的 `strict_sg_hall_ops`/`Lattice`，
并保留可回查的见证命令。
