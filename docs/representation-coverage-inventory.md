# P1 / R52：源、目标及数据边界清点

基线0929424，2026-10-27。这里计数的是固定数据目录和已验证的数据前提，**不是新增分导成功数**。
生产接口及其拒绝边界不变；两项独立只读审计分别负责目标目录和参数/spinor档案，主实现者
亲自复核汇总，新增仅测试编译的Γ矩阵/帧检查和固定见证。最终验收状态见[当前账本](<../AGENTS.md>)。

## 后续验收增量（R60，2026-10-29）

以下主体保留P1/R52原始清点语义，其中“未接生产/下一见证”属于当时状态。
P2a/P3a/P2b/P3b随后已闭合四条显式Γ路线：980+177+69+87=1,313源，固定R47几何上240,496双η调用。
R57–R60复跑四条路线，全部结果仍保持；不能据此称一般复源、spinor或非Γ磁分导完成。

R57在W之外新增[SG3 LD1LE1逐点物理线](<physical-line-subduction.md>)，一条CIR29与其派生共轭，
不伪造LE1存储行。16固定上下文×720参数采样4/8逐字一致，原始2×2矩阵重建与独立复核通过。
原5,444未注册族中此1个已有显式逐点入口，另5,443待接；原73W统计及高层行为不变。
R58补齐此物理对的精确几何/复标量等价分区：半开周期36点代表/36开区间，计各记录右端点为52点出现；
切点由精确同余证明，不由采样推断。动量碰撞与完整字符等价分开，数值算术范围、标签/gauge插值未扩大。
R59显式接PIR45 U1UA1/CIR45 U1的k₀=(1/2,0,1/2)：同样的精确切点、不同的真实目标标签，
4/8各11,520调用/223,200原始矩阵通过。seed_table()改Option：U1为None，避免不完整零偏置表；
旧数值/目标保持，Debug新增offset另说明。R59时原5,444缺口有2个独立物理入口/分区、余5,442，均不混入73W。
R60加入SG3 LD2/U2/V1/V2/W1/W2，闭合该单空间臂、标量β线八来源子域；其余F1/G1参数面及GP1GQ1参数体仍未接。
四组横向offset与±C₂点值经过原始矩阵认证；LD2同样可给完整零offset种子表，六个非零offset来源None。
全416边界CIR/owner/component与288开区间有独立夹具，4/8各92,160调用/1,785,600矩阵，旧29/45完整载荷保持。
当前原5,444缺口中8个显式入口/分区已验收、余5,436；不扩旧W73、高层路由或整个参数族能力声明。
本批按用户要求收尾后停止，不启动下一阶段。若后续恢复，先认证SG4螺旋线的Seitz/irtranslation、不同PIR/CIR编号
和完整字符周期，再继续其它线/面/体与磁SG1 Z1；
不能把P2整数平移/模1证明外推到一般affine来源。
最新v2完整电池1,052/0；夹具换列/三处测试索引修正另经完整字段等价复核与全部门禁重跑。证据与边界以[账本](<../AGENTS.md>)和[计划](<completion-plan.md>)为准。

## 1. 编译数据与磁Γ源：没有笼统的“高维都缺数据”

公共IRREPS为8,388条：4,105普通+672 compound+3,611 spinor；不是4,777条总记录。
4,777是非spinor物理记录数，所有这些记录的ScalarStar构造成功，完整维数等于物理dim。
其中1,326条selected-arm维数不同（普通1,188、compound138）；不能把dim直接当little维数。
所有8,388条对应的typed字符视图均可构造；这不是分导或矩阵群律的全域证明。

[常驻Γ源审计](<../src/irrep/magnetic_matrix_subduction/source_inventory.rs>)逐项检查非spinor Γ的
raw矩阵形状、有限值、正交性、迹与typed行、单位元、完整有序对群律、实commutant，随后复用
原ParentFrame验证完整data-Hall/心格与共同原点。结果：

| Γ源家族 | 源数 | 零/非零共同原点 | 实commutant维数 | 同母群记录×源对数 |
|---|---:|---:|---:|---:|
| 实一维普通 | 980 | 828 / 152 | 1 | 75,534 |
| 实二维非compound | 177 | 149 / 28 | 1 | 21,746 |
| 实三维普通 | 87 | 67 / 20 | 1 | 19,190 |
| 二维compound物理载体 | 69 | 59 / 10 | 2 | 3,778 |
| **合计** | **1,313** | **1,103 / 210** | 分开计数 | **120,248** |

69个compound全部是两条不同且共轭的复一维CIR的DistinctComponentSum；共138个具名CIR身份。
不是缺失的第二条伪造共轭行，也不是复域不可约二维表示。来源标签、component/irnumber与原点
必须保留。当前实一维及非compound实二维域已有全量分导验收；后三维/compound数据检查通过
**不解除**SourceDimensionNotTwo/SourceIsCompound。未来候选分别38,380和7,556个双η调用，
尚未作为生产成功计数。也未检查所有一般复输入的反幺正扩张；单独CIR不必能承载所请求的Θ。

固定三维见证SG195 GM4：12个操作、108个矩阵元素，数据与完整帧通过；交换两条同迹、不同矩阵
的行仍保持正交性与字符，却破坏群律，常驻负例拒绝它。不能以字符匹配代替矩阵证明。

## 2. 全部1,651 UNI的Γ目标目录

[清点器](<../examples/magnetic_target_inventory.rs>)及[独立计数审计](<../scripts/audit_magnetic_target_inventory.py>)：
1,651=230 Type-I+230 grey+674 Type-III+517 Type-IV；每个恰有一个Γ摘要。
返回52,793个目标摘要=Γ12,627+其它固定k40,166，fatal/未决Γ/未决其它k均0。
“没有未决摘要”不等于其全部源×记录分导已经支持。

| 目标家族 | 非spinor | spinor |
|---|---:|---:|
| Type-I，单H源：维数1 / 2 / 3 / 4 | 980 / 177 / 87 / 0 | 261 / 451 / 0 / 34 |
| Type-I，两个复H源的配对载体 | 69 | 0 |
| 有anti的A：维数1 / 2 / 3 / 4 | 5,819 / 646 / 233 / 0 | 991 / 1,933 / 0 / 80 |
| B：维数2，H一维重复两次 | 20 | 25 |
| C：维数2 / 4 | 346 / 0 | 436 / 39 |
| **合计** | **8,377** | **4,250** |

69个纯幺正packet的raw enum虽为A，但H字符范数2、无anti；不是磁C，也不是复不可约二维。
例如UNI1231 GM2+GM3，与grey UNI1232同名的磁C不同。B的H范数为4；A单源为1、C为2。
现有目标家族**形状**8,037 +非spinor三维A320 +非spinor B20 +spinor4,250=12,627；
“形状”不替代逐目标矩阵认证、完整重建及数据帧验证。

3,905摘要标Complete；8,722标TypeAAntiunitaryPending（非spinor5,718、spinor3,004）。
后者不是遗漏corep错误；反幺正字符未列出与矩阵作用尚未构造是不同概念。

### 2.1 Type-B有真实非spinor见证，不能沿用旧110-UNI抽样的零计数

20个非spinor B来自14个Type-III UNI：663,669,675,681,687,691,695,701,705,713,721,729,737,743。
H源全部是实±1一维行。独立探针检查86个anti代表均有χ_H(a²)=-1；不能以一维λK承载这个平方。

常驻最小见证：**SG75源GM3GM4 → GM2记录0/P1 → UNI663 (75.3)，目标GM2，H3 GM2**。
map=I、shift=0且exact，源矩阵和完整群/格均检查。独立构造
T(h)=χ(h)I₂、T(a₀h)=χ(h)J K，J=[[0,-1],[1,0]]、J²=-I；显式正交U对两η逐个有色操作
验证S(g)U=U T(g)（U实），而公开输入仍返回SourceIsCompound。**这是下一步的可行性证据，不是B已交付。**

14个UNI共17条isotropy记录；14条在SG75–88，需先支持compound Γ源。
另3条UNI701记录在SG220/230，仅有5个旧实二维record/source对；旧10个η调用全部是2×GM1(A)，
不是B的漏解。因此“20个B目录项”不能写成“旧43,492域新增20成功”。先做compound源，再接此窄B。

### 2.2 H数据与spinor不是同一种缺口

Γ目标共13,517个H源出现次数，2,128个唯一H身份；每个都有唯一具名owner与typed行，缺失/歧义0。
唯一身份分区：980实一维+177实二维+87实三维+138复一维CIR+261 spinor一维+451 spinor二维+
34 spinor四维。物理owner共2,059=1,313非spinor+746 spinor。

非spinor raw矩阵都有形状/有限性/typed迹证据，Γ源的进一步操作认证见§1。
全部3,611个公共spinor记录的raw irrep matrices为空（Γ746）；SU(2)操作lift仍存在，不能当作
一般不可约矩阵。Γ中485个高维spinor H身份缺存储矩阵；一维字符可成为标量矩阵，但仍须证明
双群中心符号、projective群律和反幺正扩张，不能套ηK冒充Θ²=-1。
SG16 GM5/UNI99是二维缺矩阵见证；SG3 A3的χ(C2y)=-i是窄一维见证。

621个使用CIR标签的目标附带isotropy候选为空，是组分没有独立物理记录，不是缺字符；
需保留owner/component来源，不伪造独立IrrepRecord。

## 3. 普通参数域与磁非Γ：分开推进

[公共源清点器](<../examples/irrep_inventory.rs>)及[精确档案审计](<../scripts/audit_irrep_archives.py>)
复用已有严格exact loader与已钉住的官方压缩档案，不重新生成数据。

| 原始PIR参数秩 | 记录数 | 说明 |
|---|---:|---|
| 0 | 4,777 | Γ1,313+固定非Γ3,464 |
| 1 | 3,980 | 其中2,049条有非零k₀ |
| 2 | 1,307 | 其中383个面有非零k₀ |
| 3 | 230 | 体参数族 |
| **总数** | **10,294** | little_full同数；额外6个数组槽是空槽，不是缺记录 |

CIR另为11,202=5,296固定+4,331线+1,345面+230体。PIR物理载体数不能与CIR不可约组分数混用。
编译IRREPS的全部非spinor CIR来源参数秩为0；112个物理标签与原PIR拼写不同，不是112条缺记录。

当前W层只注册73条rank1、k₀=0源，11个F立方母群，1,006上下文/5,756行；这些源不在公共IRREPS。
P1清点共有5,517个连续PIR族，其中当时**5,444个未具名注册**（R57已另接其中1个逐点入口）；73条本身也不是所有参数取值下都已变成
具名增强little-group不可约表示。现有schema为k=t·v，不能表达一般k₀+Σtᵢvᵢ。
低层允许调用者提供表，并非73条指针白名单；但表、相位、帧、little群和完整星仍需认证。

P1下一普通最小见证曾是SG3 PIR#29 LD1LE1，k=(0,b,0)；R57已认证实际CIR/物理配对/irtranslation并接逐点入口；
随后#45 U1UA1，k=(1/2,b,1/2)。面#51 F1=(a,0,c)、体SG1 #9 GP1GQ1=(a,b,c)另立任务。
这些raw-only标签当前公开查询无记录，高层报CondensingNotFound，不是缺少原始档案。

普通固定非Γ/full-star已经工作：SG213 X2/C23记录11600→SG146，selected2×3臂=6；
SG46 W1W2/C12记录1085→SG9，完整4维；SG23 W1W1/C1记录500→SG22，完整2维。
磁非Γ仍拒绝SourceNotGamma（如SG1 Z1、SG5 V1），不是普通full-star也未实现。

SG196 DT1的参数路线须固定记录10038或10062；10038,t=1/7成功，给出6个未具名构造目标；
t=0仍为Formal，不自动宣称增强点irrep。省略参数则使用固定W1上下文，不是推断未知参数。
保留这个接口差异，不从标签猜连续参数。

## 4. 元数据缺口不是数学无解

普通方向串15,044/15,239；195个例外全是当前采集下的TableEmpty，分布于112个物理记录。
原官方紧凑拼写可能有输出，这是后续采集与身份关联问题，不能说字符串原则上不存在。
磁方向串选择器仍NoOfficialDescriptor。16,721条磁记录覆盖全部1,421个非grey UNI；
230个grey都有群操作/共表示目录但没有该isotropy记录，不能伪造230条补数，也不能说磁群数据库缺它们。

## 5. 下一实施顺序与复现

1. P2a：69个二维compound Γ载体，复用认证/共同原点/复Hom；显式保存两CIR来源，保持旧入口契约。
2. P3a：与P2a配套的实一维H、χ(a²)=-1窄Type-B；用UNI663见证、全部14 UNI和固定记录测量，
   重数不能直接取实Hom维数。随后处理三维源/三维A目标，不预造无限维框架。
3. P4拆成普通连续参数注册与磁非Γ两条支路；各先做上述最小见证，再面/体与全星。
4. P5先区分一维spinor字符可构造与高维矩阵缺数据，另证anti sewing；P6/P7的原契约不取消。

从仓库根目录运行（输出目录先创建；保留不同前缀做比较）：

```sh
mkdir -p target/coverage-inventory
CARGO_TARGET_DIR=target cargo test --offline --release --lib source_inventory -- --nocapture --test-threads=1
CARGO_TARGET_DIR=target cargo build --offline --release --example irrep_inventory --example magnetic_target_inventory
target/release/examples/irrep_inventory > target/coverage-inventory/irrep-inventory.tsv
python3 scripts/audit_irrep_archives.py target/coverage-inventory/irrep
RAYON_NUM_THREADS=4 target/release/examples/magnetic_target_inventory target/coverage-inventory/magnetic-
python3 scripts/audit_magnetic_target_inventory.py target/coverage-inventory/magnetic-
python3 -m unittest discover -s scripts -p test_audit_irrep_archives.py
```

清点器不把typed视图成功当成分导成功；审计脚本任一不变量失败都退出非零。
本轮验收另比较独立2线程探针与正式4/8线程目录输出，并做丢行、错误B范数、错分纯幺正packet、
缺H来源等破坏控制。最终证据统一保存于target/logs/r52-accepted，数字以账本冻结结果为准。
