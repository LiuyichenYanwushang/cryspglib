# 普通完整分导的剩余缺口：逐星清点

2026-09-22，在 `f7941e5` 的引擎上重放全部 14,713 个 `identity_only` probe。
这份清点不扩大生产 API 的覆盖声明，也不包含参数化 w 行或磁共表示。

| 项目 | 数量 |
|---|---:|
| probe（ordinal + 母群 CDML 标签） | 14,713 |
| isotropy 记录 | 2,761 |
| 子群号 | 128 |
| 全部折叠子群星 | 22,410 |
| 无随包离散标量 k 数据的星 | 21,136 |
| 能匹配随包标量 k 数据的星 | 1,274 |
| 缺失的（子群号，规范化 k-star）组合 | 796 |
| 再细分 `(U, denominator, child_shift)` 后的组合 | 989 |
| 有不止一个缺失星的 probe | 5,784 |
| 单个 probe 最多缺失星数 | 6 |

全部 14,713 次生产入口调用仍返回 `MissingChildStarData`，其首个缺失星与清点一致。
每个 probe 的全部块维数之和等于母群 full-star 维数。缺失块全部非 Γ，最小维数为
1，累计维数 51,690；有存储 k 可达的块累计维数为 3,952。**恒等重数为零不能用于
删除这些正维数块**，它只说明它们不含子群平凡表示。

`stored_k_reachable` 只说明至少一个标量复成分的有效 k 落在该星上，不保证该星的
完整分解已经通过 Gram、整数重数、维数和与逐操作重建。特别是引擎在首个缺失星
返回后，后续星可能尚未进入求解。这 1,274 个块不能直接算作新增完整分解。

## 方法与数据帧

工具 `examples/census_subduction_gaps.rs` 从审计 TSV 选取 `identity_only` 行，
重建当前 embedding 并重放生产入口，随后复用 `ScalarStar::folded_stars` 展开
**所有**子群星。若 probe 已修复、变成别的错误、上下文改变或维数不守恒，工具退出
非零，不能把过期审计当作当前覆盖情况。它报告输入所选范围，不自行宣称范围全局；
复现下面的全局清点须使用无 `--parent` / `--ordinal` 过滤的全局审计。

匹配遍历一个星的所有 q，并按子群自己的 primitive reciprocal lattice 判等价，
不按 conventional 坐标逐分量模 1。spinor 排除；`DistinctComponentSum` 的成分用
存储 k，`ConjugateRealification` 同时提供 k 与 -k。这一步只读取 k 可达性，
不合成特征标。其规则对应引擎的 `child_components_at` / `select_representative`。

`q` 列保留原始精确有理坐标；`canonical_q` 列将各点按子群倒格约化后排序，保存
完整星，不能只拿报错中的 first q 作分组键。约化后的 conventional 分量可以为负
或大于 1（例如带心格），这不表示约化失败。setting 单独保存精确 U、分母与
`child_shift`，ordinal 可回查完整的 T/origin 和凝聚记录。

永久回归钉住 ordinal 13345、母群 #225 `W1` → 子群 #8：引擎首条错误只指出一个
两点星，但实际有 **3 个**缺失星，各维数 2，母群维数 6。这样可防止以后退回只数
首条错误的清点方式。

## 按影响 probe 数排序的前十个子群

| 子群号 | probe | 记录 | 缺失星 | k-star 组合 | 含 setting 的组合 |
|---|---:|---:|---:|---:|---:|
| 1 | 1,835 | 381 | 3,992 | 40 | 90 |
| 5 | 1,659 | 360 | 2,547 | 33 | 72 |
| 2 | 998 | 213 | 1,164 | 22 | 40 |
| 12 | 825 | 143 | 901 | 23 | 25 |
| 8 | 770 | 125 | 1,405 | 38 | 41 |
| 6 | 622 | 63 | 1,238 | 14 | 14 |
| 9 | 536 | 116 | 896 | 38 | 50 |
| 15 | 509 | 119 | 551 | 23 | 33 |
| 3 | 467 | 90 | 508 | 14 | 25 |
| 4 | 446 | 92 | 489 | 14 | 25 |

## 下一步实施范围

1. **先处理子群 #1。** 目标是升级这 1,835 个 probe；它们有 40 个缺失 q，坐标分母
   只出现 1、2、3、4、6。#1 的标量小群只有平移，其任意 q 的不可约表示是一维
   Bloch 相位。可直接按已钉死的相位约定提供目标字符，不必采集另一整张数据库。
   首例：ordinal 1045（#45 `S1S2` / C1），probe `W1W1`；第一个缺失 q 是
   `(-1/4,-1/4,-1)`，按子群倒格约化为 `(3/4,3/4,0)`，块维数 2。
2. **把新目标送入既有解块与重建检查。** 每个 q/star 的重数、维数和及逐操作字符
   必须通过；平移相位需有非零平移的独立钉值测试。没有可靠 CDML/BC 名称的数据，
   用精确 k 和目标身份标识，并明确标签缺失；不能借用 Γ 的标签伪装成非 Γ 标签。
3. **按当前 1,835 个 probe 清单逐条验收，然后重跑全局审计。** 当前的
   `--require-complete` 接受恒等-only，不能单独作为这个里程碑的验收。实际升级数
   必须由完整入口成功、维数和与重建检查来确认；其余恒等频率不得回归。
4. **再处理 #5/#2/#12 等。** 先按缺失 q 求实际 little group、检查 nonsymmorphic
   因子系统和 setting，再选择参数化目标字符或冻结离散数据。不要把 #1 的一维
   平移结论直接推广到其它子群，也不要用 pinned 恒等重数反推完整分解。

SG 209 `DT3`/`DT4` 标签的独立来源验证保持为另一工作项，本清点不增强其证据等级。

## 复现

本次输入为 `target/task9/audit_r2.txt`，SHA-256：
`ec8d21048fb830fc6de320709d53152f7696241c4b437f811bb633ff911fc07c`。
以下从 crate 目录运行；也可先用 `audit_irrep_subduction` 的两个门禁重建全局 TSV。
输出保留在被忽略的 `target/task9` 中，只有退出 0 的清点输出才可使用。

```bash
CARGO_TARGET_DIR=$PWD/target cargo test --release -p cryspglib \
  --example census_subduction_gaps
CARGO_TARGET_DIR=$PWD/target cargo run --release -p cryspglib \
  --example census_subduction_gaps -- target/task9/audit_r2.txt \
  > target/task9/gaps_all_stars.tsv
```

摘要应为：

```text
records=2761 probes=14713 stars=22410 missing_stars=21136 reachable_stars=1274 replay_errors=0 dimension_errors=0
```

从逐星输出生成可直接分配实现任务的 989 组清单，每组带复现 ordinal/probe：

```bash
python3 - <<'PY'
import collections, csv
groups = collections.defaultdict(list)
keys = ['child_sg', 'canonical_q', 'setting_numerator',
        'setting_denominator', 'child_shift']
with open('target/task9/gaps_all_stars.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['status'] == 'missing_discrete_scalar_data':
            groups[tuple(row[k] for k in keys)].append(row)
with open('target/task9/gaps_by_setting.tsv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(keys + ['probes', 'records', 'blocks', 'dimension_sum',
                            'witness_ordinal', 'witness_probe'])
    for key, rows in sorted(groups.items(), key=lambda kv: (int(kv[0][0]), kv[0][1:])):
        writer.writerow([*key,
            len({(r['ordinal'], r['probe_cdml']) for r in rows}),
            len({r['ordinal'] for r in rows}), len(rows),
            sum(int(r['block_dimension']) for r in rows),
            rows[0]['ordinal'], rows[0]['probe_cdml']])
print('missing k/setting groups:', len(groups))
PY
```

## R2 后（2026-09-22）：子群 #1 清零，余 12,878

`--require-complete` 全表审计（`hard_failures=0`、exit 0）实测
`full_success=353,382 identity_only=12,878`；在该审计上重跑同一清点工具：

```
records=2380 probes=12878 stars=18387 missing_stars=17144 reachable_stars=1243 replay_errors=0 dimension_errors=0
```

即 gap probe 14,713 → **12,878**（−1,835，正好是子群 #1 的全部缺口）、记录
2,761 → 2,380、缺失星 21,136 → 17,144、子群号 128 → 127（#1 消失）。子群 #1 的
3,992 个缺失星由现场构造的一维 Bloch 相位回答，构造来源、身份/标签契约与永久回归
见 `docs/subduction-audit.md` 的 R1/R2 小节。清点仍然只是"当前引擎的缺失星清单"，
不是覆盖声明。

按影响排序的剩余缺口（**probe、缺失星、全部折叠星是三个不同的计数**，排优先级看
probe 列；"全部折叠星"含该 probe 集合里已经能由存储 k 到达的星）：

| 子群号 | probe | 记录 | 缺失星 | 全部折叠星 |
|---|---:|---:|---:|---:|
| 5 | 1659 | 360 | 2547 | 2620 |
| 2 | 998 | 213 | 1164 | 1220 |
| 12 | 825 | 143 | 901 | 959 |
| 8 | 770 | 125 | 1405 | 1419 |
| 6 | 622 | 63 | 1238 | 1246 |
| 9 | 536 | 116 | 896 | 914 |
| 15 | 509 | 119 | 551 | 581 |
| 3 | 467 | 90 | 508 | 531 |
| 4 | 446 | 92 | 489 | 514 |
| 38 | 312 | 32 | 472 | 480 |
| 10 | 292 | 33 | 292 | 312 |
| 21 | 281 | 58 | 287 | 311 |
| … 其余 115 个子群 | 5161 | 936 | 6394 | 7280 |
| **合计** | **12878** | **2380** | **17144** | **18387** |

全部 127 个子群的 manifest 在 `target/r12_gaps.tsv`（列 `ordinal`/`probe_cdml`/
`child_sg`/`canonical_q`/`status`/`setting_*`），R3 的分类按行逐组进行，不要用
"全部折叠星"当 probe 数。
