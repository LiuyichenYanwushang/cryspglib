# 显式物理参数线：SG3 LD1LE1

## 范围与来源

[独立入口](<../src/irrep/subduction_physical_line.rs>)位于
`irrep::subduction::star::decompose::physical_line`。
当前只注册 SG3 PIR#29 LD1LE1，k=(0,b,0)，接受一个明确的有理参数 b、固定 isotropy 记录及其已认证 embedding。
记录定义限制到哪个子群，不是 LD1LE1 自己的凝聚记录。原73条 W 表、查询、参数分区及高层路由不变。
此入口返回普通字符分导；不返回生产矩阵 intertwiner，不声称磁非Γ、全部连续族或整个项目已完成。

原始PIR#29为 dim2、irtype3、kcount=pmkcount=1；唯一存储CIR#29 LD1为dim1。
第二分量是这个CIR的数学共轭，不是另一条LE1记录。报告用既有 RealificationSeed/RealificationConjugate
标记同一个CIR编号；PIR编号与CIR编号属于不同档案命名空间。查询不伪造 IrrepRecord 或 BC 别名。
不可变来源表明确保留源身份，seed_table()只表示一个复分量，不能冒充完整物理载体。

[精确源核验器](<../scripts/verify_physical_line_sources.py>)复用严格档案loader、矩阵reader及现有Q(ζ₂₄)算术：
- PIR ZIP SHA256 e909a4f0121688b0590ccaec10b0276171bc24619cf7eb562ba441268c01e121；
- CIR ZIP SHA256 f4edcb2852b83a86d1b58f29fb862d9124a227cfc90f9e1ae17d2c97585264e6。

两档案的原始操作均为 E、C₂y，Seitz平移及irtranslation均严格为零，P格。
PIR点矩阵为I₂，CIR点矩阵为1；参数只有b槽，系数(0,1,0)，k₀=0。
这些精确字段在Python与[Rust公开回归](<../tests/physical_line_subduction.rs>)两侧分别钉住。
原始参考算法先检查 t−τ 是格矢，再对非特殊模板用 u=t−τ+irtranslation，
以平移矩阵左乘点矩阵。这里才可约成u=t；没有把irtranslation一般等同于τ。
核验器还检查未接入的#45 U1UA1作为非零offset对照，不把该源的数学见证计作生产支持。

## 数学与端点

令θ=2πb n_y。CIR D(g)=exp(iθ)，原始PIR D(g)=[[cosθ,sinθ],[-sinθ,cosθ]]。
这与常见实二维化约定的正弦符号相反。以U=(1/√2)[[1,−i],[1,i]]，
U D_PIR U†=diag(D,conj(D))，字符为2cosθ，不是2D。
C₂y保持n_y，故此公式满足整个P2空间群的乘法，对任意b成立；有限精确测试另验解码矩阵与Seitz乘积。

两分量的动量为+k、−k。它们在母群倒格模意义下重合当且仅当2b为整数。
母群little co-group始终为C₂，**这不是little群增强**。特殊值仍保留二维载体及两个分量：
b=0或1为2×GM1，b=1/2为2×Z1；不把延拓载体误叫作一个新的二维不可约表示。
conjugate_stars_coincide()只报告母群动量轨道的几何重合，子群折叠可以更早重合。

例如固定Z2记录39限制到SG4：b=1/4时两动量折叠到同一点，但screw映到母群(C₂y|0,1,0)，
特征值为+i和−i。结果保留真实CIR70/Z1及71/Z2，各一份，物理owner均为Z1Z2。
错误的“一个+b表维数改2”会给2×CIR70及字符2i，而真实字符是0。

## 一次共同分导与验证

两个分量的arm索引分别保留，进入同一次既有fold_arms、build_block及reconstruct_on。
重复或重合的动量不删除分量，不拼接两个已经完成的报告。目标仍复用既有存储来源检查、
projective catalogue、整数重数、维数守恒及各臂重建；未命名构造目标保持未命名。
新增内部适配器只将种子**整个**字符（含Bloch相位及点字符）取共轭；旧单复线与固定源行为不改。

仅比较点操作不足以识别错误动量。新入口额外比较每个代表乘以三条L_H基矢的0/1/2线性组合，
保留全部验证操作及源/重建字符。包含平移平方和混合平移：b=1/4时错误的0、1/2动量对在
第一平移矩上与±1/4一样给0，却在平方上给+2而非−2。两类错误folding都由模块控制实际拒绝。
旧入口的验证操作列表及错误分类不变；ReconstructionMismatch的index在新入口对应新验证列表。

公开回归使用全部16个固定SG3记录，包括P1、倍胞、SG4 screw、SG5心格及非零原点。
对正负、重合、远平移参数另取生产验证立方体外的负/混合格矢，独立消费每个目标copy，
以固定复化基从目标字符重建完整原始2×2实矩阵。它是测试oracle，不是一个额外生产intertwiner API。
来源编号、旧73表、目标CIR身份、重复重数、错误记录/母群与算术溢出均有控制；没有放宽容差。

## 尚未接入

- 这是**逐个参数值**的接口，不是整个物理配对族的精确参数分区。旧单CIR分区不能描述跨分量碰撞。
- #45 U1UA1的k₀=(1/2,0,1/2)使x/z原胞平移取−1；必须显式支持affine offset，不能藏在方向里。
- 其它来源、非零原始irtranslation、面/体、多维little矩阵、一般setting及磁非Γ仍独立推进。
- 高层以凝聚标签查W行的旧请求仍不接受LD1LE1；使用本显式来源+固定几何入口，不伪造来源关系。

复现：

```text
CARGO_TARGET_DIR=target cargo test -p cryspglib --test physical_line_subduction
CARGO_TARGET_DIR=target cargo test -p cryspglib --lib physical_line::tests
python3 scripts/verify_physical_line_sources.py
```
