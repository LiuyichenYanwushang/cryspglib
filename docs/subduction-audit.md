# 完整分导的全表审计（任务 9）

任务 9 的工具分别检查几何来源和生产 API 的计算覆盖。**收集到官方基矢、恒等项
吻合、完整分解成功是三个不同的结论**；任何一个都不能替代另外两个。

## 生产 API 审计

在工作区根目录运行：

```bash
CARGO_TARGET_DIR=$PWD/cryspglib/target cargo run --release -p cryspglib \
  --example audit_irrep_subduction -- --output /tmp/subduction.tsv
```

`--parent 221` 或 `--ordinal 12400` 限定诊断范围。`--require-complete` 要求所选范围
内全部标量源表示都有完整结果，并且附加波矢记录的参数已经解决；退出码 2 表示
覆盖不完整，退出码 1 表示发现不一致。默认运行允许明确报告的覆盖缺口，退出 0
只表示未发现不一致，不能据此声称全表实现完成。

每条子群记录遍历该母群的全部标量源表示，调用实际的
`subduce_full_star_with_embedding`。恒等项为零的几何证明作为交叉检查保留，但
**不能跳过完整分解调用**。ordinal 13345 的 W1–W5 是永久反例：它们没有恒等项，
却缺少完整分解所需的子群 k 数据；把它们仅记成零项会错误宣告完整覆盖。

报告分别列出：

- embedding 成功、歧义、无有效 setting 及其它错误；
- 每个源表示的完整结果、数据缺失、计算错误和 embedding 不可用；
- 94,271 条存储正项的比较状态，以及未存储的零项；
- Γ Frobenius 检查实际完成的条数、未计算条数；
- 5,756 条 `other_wave_vector_subduction` 的 k 参数缺失。

每个总数都必须由互斥类别之和复原。重复源表示不能直接累加频率；同一源的不同
频率必须报错。当前 pinned 原始表独立清点得到 94,271 个不同的
`(isotropy ordinal, source irrep)`，没有重复行。

成功的完整分解经过维数、整数重数和逐操作特征标重建检查。这些计算检查及
恒等频率比较仍不足以固定非恒等 irrep 标签；标签依赖下面的官方坐标约定证据。

## setting 的来源与冻结

```bash
python3 scripts/audit_subduction_settings.py \
  --output /tmp/subduction-settings.jsonl --workers 8
python3 scripts/generate_subduction_settings.py --check
```

采集器从校验 SHA-256 的 `iso.zip` 使用官方程序及数据，每个 irrep 查询放在独立
临时目录，固定 `SET I ALL OR 1`。输出保留实际 ordinal、源 irrep 身份、方向标签、
官方基矢及 origin、查询失败原因和统计分母。官方空表、格式错误、origin 不一致
与不能得到整数 unimodular 换基的情况分别记录。产物是候选证据，不直接启用引擎。

冻结项按 ordinal 选择。换基由官方 conventional 基矢 `B` 精确反算：

```text
U = W P_parent (P_sub B)^-1
T = B^T
```

不能用“从母群里找到了一个同构的子群”替代这条有向基矢关系；也不能从操作集
匹配中任意挑选 child origin。冻结的 child origin 约定须保留来源，并逐操作验证。
冻结项若验证失败必须报错，不能回退到另一个候选。

`tests/subduction_settings.rs` 固定了 #43 的坐标轴陷阱：交换两个子群轴保持完整
操作集，却交换 GM3/GM4。源特征标在 x 法向镜面上的迹分别为 +1/−1，因此这种
错误可以通过所有群闭包检查，却必须被非恒等标签回归拦住。另有真正的 shear
矩阵见证，保证冻结机制不会被无意限制成 signed permutations。

## 本轮全表结果

扫描遍历所有 230 个母群和 15,239 条普通 isotropy 记录。运行时冻结项从旧的 69 条
扩充至 75 条，完整结果从 2,337 增至 2,453；这仍是部分覆盖，不是全表通过声明。

| 项目 | 实际结果 |
|---|---|
| embedding | 75 成功；12,884 歧义；2,280 无有效 setting |
| 所有标量 probe 请求 | 366,260 |
| probe 结果 | 2,453 完整；19 缺数据；363,788 因 embedding 未计算 |
| 恒等正项 | 494 比较通过；93,777 未计算；合计 94,271 |
| 已完整计算的零项 | 1,959 通过 |
| Γ Frobenius | 19 通过；1,876 未计算；合计 1,895 |
| 其它波矢记录 | 5,756 条参数缺失，0 条完成分解 |
| 不一致及计数错误 | 0 |

新增六条记录的全部 120 个标量 probe 已逐项清点：116 个完整结果，其中 18 个正项、
98 个零项；4 个缺数据请求精确为 ordinal 15125/15131 的 `P1P2`、`P3`。
缺数据总数由 15 增至 19，原因是计算覆盖扩大；永久测试固定了完整缺失集合。

官方采集器完成 4,777 个查询，输出 **0-based ordinal 0..15238** 的全部记录。
13,861 条候选通过 origin 精确相等与整数 unimodular 换基初检（其中 95 条 U
不是 signed permutation）；238 条基矢关系不符、945 条其余 origin 不符、195 条
官方空表，共 15,239 条。基矢不符项中有 112 条同时 origin 不符，因此不区分优先级
的 origin 不符总数是 1,057。这些候选数不代表引擎覆盖；未验证完整操作及坐标
约定的候选不会自动冻结。

后续扩充必须重新运行审计并更新实际覆盖。磁群、spinor 和离散子群 irrep 数据
未提供的 k 不会因这些工具而自动获得支持；剩余类别没有清零前，任务 9 的普通表
覆盖闭合仍未完成。
