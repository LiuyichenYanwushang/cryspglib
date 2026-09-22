# R4 批次卡（按可验证的目标族扩展）

R4 的规则来自 [subduction-next-milestones.md](subduction-next-milestones.md) §R4：一次只做
一个明确目标族，每批先离线证明特征标、代表元与相位，再接生产求解，并给出该批
**原缺口 probe 的完整入口成功数 / 仍缺数据数 / 新错误数**；既有成功集合不得退化。
批次顺序以 R3 的可行性为准（[subduction-gap-sources.md](subduction-gap-sources.md)）。

全表基线（R2/R3 后）：`full_success = 353,382`、`identity_only = 12,878`
（`probe = 366,260`）、`hard_failures = 0`。

## 批次 1：平凡小余群的构造目标（215 组解析路线）

**状态：已交付（2026-09-22）。** 实测 `full_success = 357,033`
（+3,651）、`identity_only = 9,227`（−3,651）、`hard_failures = 0`；
存储恒等正项 94,271/0 不匹配、Γ Frobenius 1,895/1,895、w 行 5,756/0 错误不变。
重跑清点后剩余 9,227 个 probe、12,932 个缺失星（原 12,878 / 17,144），
R3 分类器在新缺口上只剩 **684 个 `parameterized_source`**（0 个解析、0 个无源），
即批次 1 恰好关掉了 R3 判定为解析的那一批。

### 范围（R3 输入，可复算）

R3 的 899 个缺口组里 `analytic_general_position = 215`（41 个子群），它们的**每个**
星臂小余群阶都是 1：小群只剩平移，因此 q 上的目标就是一维 Bloch 相位
`D(T_L) = exp(+2πi q·L)`，不需要任何归档字符数据。

| 项目 | 数量 |
|---|---:|
| 组 | 215（41 个子群） |
| manifest 行（缺失星） | 4,212 |
| 只含解析星、本批可完全闭合的 probe | 3,651 |
| 同时含参数化星的混合 probe（本批只闭合其解析星） | 119 |

命令：

```bash
python3 scripts/classify_subduction_gap_sources.py target/r12_gaps.tsv > target/r3_groups.tsv
```

### 机制

R2 只对子群 #1 构造目标（`child_sg != 1` 直接返回空，保持 `MissingChildStarData`）。
本批把判据从「子群号」换成**几何事实**：子群自身 data-Hall 操作里，固定在 q（模子群
原胞倒格）的旋转只有恒等 → 构造一维 Bloch 相位目标；否则仍返回空，保持
`MissingChildStarData`，由批次 2 处理。

两条**星级**规则是本批实现中实测踩出来的，写在这里以免重犯：

1. **构造是星级回退，不是点级回退**：一个物理子群 irrep 活在**整个星**上。星里只要
   任何一条臂命中了 pinned 行，就必须只用 pinned 行（否则同一个 irrep 会以两个身份
   出现，求解器的正交性检查直接报错：ordinal 1007 的 `V1`/`L1` 就是见证）；
2. **一个星只构造一次**，落在该星规范化排序的第一个点上：平凡小余群在整个星上只有
   一个小群 irrep，逐点各构造一个会得到两个相同字符的行（ordinal 1007 的
   `(0,1/2,1/2)` 与 `(0,3/2,1/2)` 就是见证）。

两处必须同时成立才能算「目标完整」：

1. **小群块（`H_q`）**：`transported_values` 用构造成的相位表示在子群操作上取值（已有）；
2. **子群全星**：重建阶段把目标按**子群自己的星**诱导后逐子群陪集代表元比较。子群点群
   非平凡时星大小 > 1，因此构造目标也必须有诱导星（本批新增 `ConstructedStar`），
   不能像 #1 那样把「小群特征标 = 全星特征标」。

### 文件所有权

- `src/irrep/subduction_star.rs`：`ConstructedLittleRep` 移入此处；新增
  `LittleCharacter`（存储行 / 构造表示的统一字符来源）与 `ConstructedStar`；
  `arm_character`、`induced_component_character`、`induced_character` 改为走该来源，
  存储行路径行为不变。
- `src/irrep/subduction_star_decompose.rs`：`constructed_child_components_at` 的判据、
  `ChildStarEvaluator::Constructed` 持有构造星、单元测试。
- `tests/subduction_constructed_stars.rs`（新）：端到端与负例。
- `docs/subduction-r4-batches.md`、`docs/subduction-gap-sources.md`、`CLAUDE.md`。

### 验收（实测）

- `cargo test --release -p cryspglib --tests --doc` 全绿（lib 392、integration 160、
  doctest 27；audit example 17、census example 2）；严格 all-target clippy 通过。
- 全表审计（`--require-complete --require-full-decomposition`，518 s）：
  `full_success=357033 identity_only=9227`、`hard_failures=0`、
  恒等正项 `94271/94271 passed, mismatch=0`、Γ Frobenius `1895/1895`、
  w 行 `5756 computed / 0 mismatched / 0 engine_errors`、
  `production_checks` 五项全 0、`accounting_violations=0`。
  `--require-complete` 仍 exit 0；`--require-full-decomposition` exit 2
  （`incomplete=9227`，R5 的目标）。
- 新增回归：
  - `tests/subduction_constructed_stars.rs`（4 项）：ordinal 13345（SG 225 → #8）
    31/31 完整且恒等重数与 pinned 频率一致、两臂星的诱导（`block_dimension =
    star_size × little_dimension`）、ordinal 13346 的 10 个参数化星仍是
    `MissingChildStarData`、以及「有 pinned 行的星绝不再构造」；
  - `src/irrep/subduction_star_decompose.rs`：
    `a_constructed_star_induces_over_the_child_star`（子群 #2，手算
    `χ(E)=2`、`χ(T_(1,0,0))=0`、`χ(T_(0,1,0))=−1`、`χ(−I)=0`）；
  - `examples/census_subduction_gaps.rs`：清点的可达性模型改为「无 pinned 行的星，
    平凡小余群 → `constructed_trivial_co_group`，否则 → `missing_discrete_scalar_data`」，
    两条测试分别钉住 13345（0/3/0）与 13346（0/1/4）；
  - `examples/audit_irrep_subduction.rs`：13345 由「五个 W 恒等-only」改为
    「31/31 完整、全门禁 clean」，缺口的微型基线移到 13346；
  - `tests/subduction_identity_regressions.rs`：冻结上下文里 15 个缺数据 probe → 10 个，
    比较数 2,060 → 2,065、非 Γ 1,522 → 1,527。
- 负例：非平凡小余群且无 pinned 行的折叠星**仍然**是 `MissingChildStarData`
  （ordinal 13346 的 `B` 线四星，R3 的 `parameterized_source`），批次边界未被越过。

### 尚未做（批次 2 入口）

- 684 个 `parameterized_source` 组（123 个子群、9,227 个 probe）：需要把归档参数域
  代入求解参数、取出该 q 上的小群字符/矩阵并验证适用域与目标完整性；
  R3 已确认 684/684 的 PIR 矩阵块可读（81,576 个矩阵元）。
- 参数化来源可复用则冻结最小数据，运行时保持纯 Rust（R4 规则）。
- `little_k`（w 源数组）与 PIR `k_arms` 是不同数组，各自的坐标约定要分别验证，
  不能互相套用。
