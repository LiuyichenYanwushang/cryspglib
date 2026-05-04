# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

# v0.2.0: Current API State

> **Status**: DONE. C-style `spg_*` replaced with Rust-idiomatic `Crystal` + `SymmetryAnalysis`.

## Public API

### Core entry point

```rust
use cryspglib::Crystal;

let si = Crystal::new(
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    vec![[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
    vec![14, 14],
);
let ds = si.analyze().symprec(1e-5).dataset()?;
```

### Crystal methods

| Method | Returns | Description |
|--------|---------|-------------|
| `Crystal::new(lat, pos, types)` | `Self` | Non-magnetic 3D crystal |
| `.with_magnetic(moments)` | `Self` | Add magnetic moments `[mx,my,mz]` per atom |
| `.with_layer(axis)` | `Self` | Mark as 2D slab |
| `.analyze()` | `SymmetryAnalysis` | Begin symmetry analysis |
| `.delaunay_reduce(prec)` | `Result<Mat3>` | Delaunay reduction |
| `.niggli_reduce(prec)` | `Result<Mat3>` | Niggli reduction |
| `.natom()` | `usize` | Number of atoms |
| `Crystal::from_poscar(data)` | `Option<Self>` | Parse POSCAR format |

### SymmetryAnalysis terminal methods

| Method | Returns | Description |
|--------|---------|-------------|
| `.dataset()` | `Result<SpaceGroup>` | Full space group dataset |
| `.symmetry()` | `Result<SymmetryOps>` | Symmetry operations only |
| `.primitive_cell()` | `Result<Crystal>` | Primitive cell |
| `.standardize(to_prim, no_ideal)` | `Result<Crystal>` | Standardized cell |
| `.hall_number()` | `Result<usize>` | Hall number (1-530) |
| `.international()` | `Result<(usize, String)>` | Space group number + symbol |
| `.irreducible_mesh(mesh, shift, tr)` | `Result<IrMesh>` | Irreducible k-point grid |
| `.magnetic_dataset()` | `Option<MagneticSymmetry>` | Magnetic space group |

### Key types

| Type | Description |
|------|-------------|
| `SpaceGroup` | Space group data (number, symbol, Wyckoff, symmetry ops) |
| `SpaceGroupType` | Space group info from database. Has `.from_hall(n)` constructor |
| `MagneticDataset` | Magnetic space group dataset |
| `MagneticSymmetry` | Magnetic symmetry analysis result. Implements `Display` |
| `MagneticSpaceGroupType` | Magnetic SG from database. Has `.from_uni(n)` constructor |
| `SymmetryOp` | Single operation `{R\|t}` + `time_reversal: bool` |
| `SymmetryOps` | Collection of `SymmetryOp` |
| `IrMesh` | Irreducible k-point mesh |
| `SymError` | Error enum (SpacegroupSearchFailed, AtomsTooClose, ...) |
| `MagneticType` | NonMagnetic/Ordinary/Grey/BlackWhite/AntiTranslation |

### Old functions still present (to be deprecated)

| Function | Replacement |
|----------|-------------|
| `spg_get_magnetic_dataset(...)` | `crystal.analyze().magnetic_dataset()` |
| `spg_format_magnetic_symmetry(&r)` | `r.to_string()` (Display impl) |
| `spg_get_pointgroup(&rots)` | Keep as pure math utility |
| `spg_get_symmetry_from_database(n)` | To be added to `SymmetryOps` |
| `spg_delaunay_reduce` / `spg_niggli_reduce` | `crystal.delaunay_reduce()` / `crystal.niggli_reduce()` |
| k-point grid functions (6) | `crystal.analyze().irreducible_mesh()` |

## API design rules

1. **No `spg_` prefix** on new types/functions
2. **No `spgat_` variants** — use `.angle_tolerance()` builder
3. **No pre-allocated output buffers** — functions return owned types
4. **Crystal is the single entry point** — not naked (lattice, pos, types) tuples
5. **Builder pattern for optional parameters** — `SymmetryAnalysis` with `.symprec()`, `.angle_tolerance()`
6. **Deprecated aliases** — old names (`SpglibDataset`, `SpglibError`) kept as type aliases

## Integer type convention

- **All integer types use `usize` / `isize`** except for enums which may use `u8`.
- Never use `i16`, `i32`, `i64`, `u16`, `u32`, `u64` in new code — use `usize` or `isize`.
- The `IrrepRecord` fields (`kx`, `ky`, `kz`, `kd`, `dim`, `sg`) use `i8`/`u8` only because they are packed into a large static array (8,388 entries) where memory is a concern.
- Internal pointer fields (`_char_start`, `_mat_start`, etc.) use `u32`/`u16` for the same memory-compaction reason.

## Module structure

| Module | Role |
|--------|------|
| `api.rs` | Crystal, SymmetryAnalysis, SymmetryOp — new public API |
| `parser.rs` | POSCAR parser + 118-element periodic table |
| `lib.rs` | Re-exports, type definitions, remaining old functions |
| `cell.rs` | Internal `Cell` struct (lattice + positions + tensors) |
| `mathfunc.rs` | `Mat3`, `Mat3I`, `Vec3` aliases + matrix ops |
| `symmetry.rs` | Symmetry operation detection |
| `spacegroup.rs` | Space group search pipeline |
| `primitive.rs` | Primitive cell detection |
| `determination.rs` | Top-level retry loop |
| `refinement.rs` | Wyckoff + exact structure |
| `spg_database.rs` | Hall symbol database |
| `hall_symbol.rs` | Hall symbol matching |
| `pointgroup.rs` | Point group classification |
| `delaunay.rs`, `niggli.rs` | Lattice reduction |
| `kgrid.rs`, `kpoint.rs` | k-point grid generation |
| `magnetic_spacegroup.rs`, `spin.rs` | Magnetic symmetry |
| `overlap.rs`, `sitesym_database.rs`, `site_symmetry.rs` | Site symmetry |

## Matrix convention

`lattice[cart][vec]`: rows = Cartesian (x,y,z), columns = lattice vectors (a,b,c).

```
lattice = [[a_x, b_x, c_x],
           [a_y, b_y, c_y],
           [a_z, b_z, c_z]]
```

## Build & Test Commands

```bash
cd /home/liuyichen/TB_rs
cargo build --package cryspglib
cargo test --package cryspglib
cargo check --package cryspglib
```

## Deprecation Status (v0.2.0)

All `spg_*` C-style functions now have `#[deprecated]` annotations pointing to Rust-idiomatic equivalents. Internal callers use `#[allow(deprecated)]` where the old function serves as the implementation.

### Deprecation summary (17 functions total)

| New API | Old function (deprecated) |
|--------|---------------------------|
| `Crystal::new(...)` | entry point (no old equivalent) |
| `cry.analyze().dataset()` | `spg_get_dataset` |
| `cry.analyze().symmetry()` | `spg_get_symmetry` |
| `cry.analyze().primitive_cell()` | `spgat_find_primitive` |
| `cry.analyze().standardize(to_prim, no_ideal)` | `spgat_standardize_cell`, `spgat_refine_cell` |
| `cry.analyze().hall_number()` | (new) |
| `cry.analyze().international()` | (new) |
| `cry.analyze().magnetic_dataset()` | `spg_get_magnetic_dataset` |
| `cry.analyze().irreducible_mesh(...)` | `spg_get_ir_reciprocal_mesh`, `spg_get_dense_ir_reciprocal_mesh` |
| `cry.delaunay_reduce(prec)` | `spg_delaunay_reduce` |
| `cry.niggli_reduce(prec)` | `spg_niggli_reduce` |
| `Crystal::from_poscar(data)` | `spg_read_structure` |
| `SpaceGroupType::from_hall(n)` | `spg_get_spacegroup_type` |
| `MagneticSpaceGroupType::from_uni(n)` | `spg_get_magnetic_spacegroup_type` |
| `MagneticSpaceGroupType::classify(...)` | `spg_get_magnetic_spacegroup_type_from_symmetry` |
| `SymmetryOps::from_database(n)` | `spg_get_symmetry_from_database` |
| `result.to_string()` (Display) | `spg_format_magnetic_symmetry` |
| `grid_point_from_address(addr, mesh)` | `spg_get_grid_point_from_address` |
| `stabilized_reciprocal_mesh(...)` | `spg_get_stabilized_reciprocal_mesh` |
| `dense_grid_points_by_rotations(...)` | `spg_get_dense_grid_points_by_rotations` |
| `dense_bz_grid_points_by_rotations(...)` | `spg_get_dense_BZ_grid_points_by_rotations` |
| `relocate_bz_grid_address(...)` | `spg_relocate_BZ_grid_address` |

### Kept as utility (no deprecation)

- `spg_get_pointgroup(&rotations)` — pure math utility, no Crystal dependency
- `spg_get_hall_number_from_symmetry(...)` — low-level symmetry classifier

### New types added this session

| Type | Description |
|------|-------------|
| `StabilizedMesh` | Stabilized reciprocal mesh with q-points |
| `BzMesh` | BZ-relocated grid addresses + mapping |

## v0.3.0 需求: 完整交互式查询 API

### 用户工作流

```
输入 SG# → 查看所有高对称 k 点 → 选 k 点 → 查看特征标表 → 选 irrep → 查看 isotropy subgroup
```

### Rust API 需求

```rust
// 1. 输入 SG → 所有对称操作 + k 点列表
let sg = SpaceGroupView::new(221);
sg.symmetry_operations()  // &[SymmetryOp] (from spg_database)
sg.kpoints()              // &[KPointSummary] (name + coords, deduplicated)

// 2. 选 k 点 → 特征标表 (两种模式)
let kp = sg.kpoint("GM");
kp.irreps()               // &[IrrepAtK] (label + dim + characters + full_matrix)
kp.character_table()      // 特征标表 (trace only)
kp.full_matrix_table()    // 完整矩阵表 (dim×dim per operator)

// 3. 选 irrep → isotropy subgroups
let ir = kp.irrep("GM4-");
ir.characters()           // &[f64] — trace per operator
ir.matrices()             // &[f64] — dim² values per operator, flattened
ir.subgroups()            // &[IsotropyRecord] (direction, SG, size, basis, origin)
```

## v0.3.0: Irrep (不可约表示) 模块

### 概述

基于 ISOTROPY Suite (Stokes & Campbell, 2022) 数据，为全部 230 个空间群提供：
- 所有高对称 k 点的不可约表示（含 CDML/BC/Kovalev 三种标签）
- 每个 irrep 的特征标表 χ(g) = Tr(D(g)) 和完整矩阵 D(g)
- 每个 irrep 的 isotropy subgroup（方向、SG编号、domain 数量）
- `format_character_table()` 以 {R|t} 操作标签为列头的特征标表

### 数据流

### zip 档案及其内容

#### 核心数据 (3D 空间群, 有理 k 矢量)

| 文件 | 大小 | 内容 |
|------|------|------|
| `iso.zip` | ~90 MB | 27 个文件: 元数据 + Fortran 工具 |
| ├ `data_irreps.txt` | 1 MB | 4,777 个 irrep 的 ML/BC/Kovalev 标签, SG编号, image码, Lifshitz |
| ├ `data_isotropy.txt` | 6 MB | 15,239 个 isotropy subgroup (方向, domain数, basis, origin) |
| ├ `data_little.txt` | 10 MB | 每个 SG 的 k 点小群: k_count, k_label, k_kov |
| ├ `data_images.txt` | 144 KB | 141 个 image 符号 (A1a, C24c, ...) |
| ├ `data_space.txt` | 570 KB | 点操作标签 (E, C2x, C3+, ...) + 空间群符号 |
| ├ `data_magnetic.txt` | 2.5 MB | 磁空间群数据 (1,651 个) |
| ├ `data_wyckoff.txt` | 712 KB | Wyckoff 位置数据 |
| ├ `data_ssg.txt` | 25 MB | 超空间群数据 (16,697 个) |
| ├ `data_ssgmag.txt` | 30 MB | 磁超空间群数据 (325,127 个) |
| ├ `data_diperiodic.txt` | 11 KB | 双周期 (层) 群数据 |
| ├ `const.dat` | 3 KB | pir_data_constant 表 (25 个浮点数, 矩阵元素解码) |
| ├ `*.f` | — | Fortran 源码 (未打包在 zip 中, 仅 PIR/CIR zip 有) |
| ├ `comsubs`, `findsym`, `iso`, `smodes` | ~8 MB | Fortran 可执行文件 (Linux) |
| └ `*sample*`, `*.txt` | ~30 KB | 示例输入/输出 + 文档 |
| `PIR_data.zip` | 29 MB | **实不可约表示** (Physical IR) |
| ├ `PIR_data.txt` | 29 MB | 10,294 条 irrep 矩阵: k-vector + 操作矩阵 + dim² 浮点矩阵 |
| └ `PIR_data.f` | 30 KB | Fortran 读取源码; 定义 `pir_data_constant` 编解码表 |
| `CIR_data.zip` | 54 MB | **复不可约表示** (Complex IR, PIR 的回退数据源) |
| ├ `CIR_data.txt` | 54 MB | 11,202 条 irrep 矩阵: 同样格式, 矩阵元素为 `(实部,虚部)` |
| └ `CIR_data.f` | 28 KB | Fortran 读取源码 |

**注意**: PIR 存实表示 (H2H3 为 4D 实矩阵), CIR 存复表示 (H2, H3 为独立 2D 复矩阵)。CIR 是 PIR 的超集: CIR→PIR 可转换 (复共轭对→实块矩阵), 反之不可。当前生成器先用 PIR, 缺失时回退到 CIR 并做复→实转换。

#### 超空间群数据 (3+d 维, 无理 k 矢量)

| 文件 | 大小 | 内容 |
|------|------|------|
| `PIR_SSG_data.zip` | 37 MB | 超空间群 PIR |
| ├ `PIR_SSG_data.txt` | 37 MB | 超空间群实不可约表示矩阵 |
| └ `PIR_SSG_data.f` | 28 KB | Fortran 读取源码 |
| `CIR_SSG_data.zip` | 61 MB | 超空间群 CIR |
| ├ `CIR_SSG_data.txt` | 61 MB | 超空间群复不可约表示矩阵 |
| └ `CIR_SSG_data.f` | 25 KB | Fortran 读取源码 |

**注意**: SSG (Superspace Group) 数据目前**未使用**, 仅保留以备将来扩展。

### 数据流

```
iso.zip/data_irreps.txt     → 元数据 (标签, SG, image, isotropy)
iso.zip/data_isotropy.txt   → isotropy subgroups
PIR_data.zip/PIR_data.txt   → 实表示矩阵 (主数据源)
CIR_data.zip/CIR_data.txt   → 复表示矩阵 (回退数据源)

                    ↓ scripts/generate_irrep_data.py ↓

src/irrep/generated_data.rs   → 5 个静态数组, ~8.7 MB, 编译进二进制
```

### Rust 源文件

| 文件 | 功能 |
|------|------|
| `src/irrep/mod.rs` | 模块入口, 230 SG 索引表 + rustdoc |
| `src/irrep/types.rs` | `IrrepRecord`, `IsotropyRecord`, `MagneticIsotropyRecord` |
| `src/irrep/query.rs` | `irreps_of()`, `kpoints_of()`, `format_character_table()` |
| `src/irrep/wigner.rs` | Wigner 测试: Seitz 组合, Bloch 相位, 分类函数, SU(2) 辅助 |
| `src/irrep/corep.rs` | Co-representation API: `compute_coreps()`, `MagneticOps`, `CorepType` |
| `src/irrep/generated_data.rs` | 自动生成: 15+ 个静态数组, ~18 MB |
| `src/irrep/preamble.rs` | 230 SG 的 HM 符号 + Schoenflies 符号 |
| `src/irrep/{triclinic,monoclinic,...,cubic}.rs` | 7 个晶系文件的 rustdoc 页面 |

### 静态数组 (全部在 generated_data.rs, ~18 MB 总计)

| 数组 | 大小 | 说明 |
|------|------|------|
| `SG_IRREP_INDEX` | 231 | SG# → (start, count) into IRREPS |
| `IRREPS` | 8,388 | 所有 irreps (4,777 scalar + 3,611 spinor) |
| `ISOTROPY_SUBGROUPS` | 15,239 | 所有 isotropy subgroups |
| `MAGNETIC_ISOTROPY_SUBGROUPS` | 16,721 | 磁 isotropy subgroups |
| `CHARACTERS` | ~106K f64 | 所有特征标值 (含 spinor) |
| `MATRICES` | ~634K f64 | 所有矩阵元素 |
| `CIR_COMPONENT_CHARS` | ~23K f64 | Compound irrep CIR 复字符 |
| `CIR_ROTS` | ~104K i32 | CIR 操作旋转矩阵 |
| `PIR_ROTS` | ~500K i32 | PIR 操作旋转矩阵 |
| `SPIN_OP_ROTS/TRANS/SU2` | ~73K total | Spinor 操作数据 |
| `SPIN_EXTRA_CHARS` | ~17K f64 | Spinor Wigner extra 值 |
| `SPIN_OP_SG_INDEX` | 231 | SG# → spin op 索引 |

### IrrepRecord 字段

```rust
pub struct IrrepRecord {
    pub sg: u8,        // SG 编号 1-230
    pub ml: &str,      // CDML/Miller-Love 标签 "GM4+"
    pub bc: &str,      // Bradley-Cracknell 标签 (LaTeX)
    pub kov: &str,     // Kovalev 标签 (LaTeX)
    pub dim: u8,       // 维度 1-24
    pub image: &str,   // Stokes-Hatch image "C24c"
    pub lifshitz: bool,
    pub kx, ky, kz, kd: i8,  // k 矢量坐标
    // 内部指针:
    _char_start, _char_count,  // → CHARACTERS
    _mat_start, _mat_count,    // → MATRICES
    _iso_start, _iso_count,    // → ISOTROPY_SUBGROUPS
}
// 方法: .characters(), .matrices(), .subgroups(), .k_label(), .is_point()
```

### Python 生成脚本

| 文件 | 功能 |
|------|------|
| `scripts/generate_irrep_data.py` | **主生成器**: 解析 iso/PIR/CIR zip → generated_data.rs |
| `scripts/direction_map.py` | 657 个方向码 → 人类可读标签 |
| `scripts/parse_characters.py` | 单独的 PIR 特征标解析器 (调试用) |
| `scripts/export_irrep_json.py` | 导出 JSON 用于交互式 HTML 浏览器 |
| `scripts/generate_irrep_docs.py` | 生成各晶系文件的 rustdoc 表格 |
| `scripts/extract_*.py` | PDF 数据提取 (实验性) |

### 匹配策略 (ISO 标签 → PIR/CIR 数据)

1. 精确匹配
2. 复合标签模糊匹配 (去掉尾部数字/符号)
3. 前缀匹配
4. k-vector 坐标 + 数字编号匹配 (跨命名惯例: X/Y/Z ↔ A/B/C)
5. k-vector 坐标 + 字母前缀 + 符号匹配 (编号偏移)
6. CIR 回退: 精确匹配 + 复合标签拆解 (H2H3 → H2 + H3)
7. CIR 复→实矩阵转换: [[A,-B],[B,A]] 或块对角复制

### 重新生成数据

```bash
cd cryspglib
python3 scripts/generate_irrep_data.py
cargo test --package cryspglib --test irrep_validation
```

### 注意事项

- `isotropy_subgroup/*.zip` 不在 git 中 (.gitignore 排除), 也不在 cargo publish 中
- PIR_data 存储**实表示** (PIR), CIR_data 存储**复表示** (CIR)
- CIR 是 PIR 的超集: CIR→PIR 可转换 (取实部), PIR→CIR 不可逆
- 生成器从 zip 直接读取, 不需要解压
- `domains` 和 `arms` 是 `usize`, 因为值可达 384

### 磁群编号系统

**ISOTROPY `mag_nlabel` 索引 = spglib UNI 编号** (1-1651)。两者使用同一套编号系统。
`MagneticIsotropyRecord.mag_sg` 存储的就是 UNI 号，可以直接传给 `corepresentation()`。

| 编号类型 | 128.406 示例 | 说明 |
|---------|-------------|------|
| UNI | 1066 | spglib 内部编号，`MagneticSpaceGroupType.uni_number` |
| BNS | `"128.406"` | Belov-Neronova-Smirnova 标签，`MagneticSpaceGroupType.bns_number` |
| OG | `"128.8.1073"` | Opechowski-Guccione 标签，末尾 1073 是 **Litvin 编号**（非 UNI！） |
| Litvin | 1073 | Litvin 编号，出现于 OG 标签末尾，易与 UNI 混淆 |

### corep API (`src/irrep/corep.rs`)

| 函数 | 签名 | 说明 |
|------|------|------|
| `identify_unitary_subgroup` | `(uni: usize) -> Option<usize>` | 从磁群 UNI 识别其 unitary subgroup 的 SG 号 |
| `uni_from_bns` | `(bns: &str) -> Option<usize>` | BNS 标签 → UNI 编号 |
| `uni_from_og` | `(og: &str) -> Option<usize>` | OG 标签 → UNI 编号 |
| `get_magnetic_operations` | `(uni: usize) -> Option<MagneticOps>` | 获取磁群对称操作 |

`IrrepRecord::corepresentation(uni)` 对非磁 irrep 现场计算磁共表示，无需预存表格。

**已验证**: 128.406 (UNI 1066) 的 unitary subgroup = SG 118 (P-4n2)，与 BCS 一致。

---

## v0.3.0: Co-representation (corep) / Wigner 测试实现

### 概述

基于 Wigner (1959) 和 Bradley & Cracknell (1972) 理论，从磁空间群的 BNS 标签自动计算磁共表示 (corepresentation)。核心流程：

```
BNS label ("128.406") + k-point label ("Z")
  → uni_from_bns()              // BNS → UNI number
  → identify_unitary_subgroup()  // UNI → H (unitary subgroup) SG number
  → get_magnetic_operations()    // 获取磁群对称操作 {R|t,θ}
  → irreps_of(H) at k-point      // H 的不可约表示
  → compute_corepresentation()   // Wigner 测试 → corep type + character table
```

### BCS 128.406 Z-point 验证结果 (全部通过)

| Irrep | Type | W | Dim | BCS | 路径 |
|-------|------|---|-----|-----|------|
| Z1Z4 | C | W=(0,0) | 4 | C | CIR 路径 |
| Z2Z3 | C | W=(0,0) | 4 | C | CIR 路径 |
| Z5 | A | W=1.0000 | 2 | A | PIR 路径+重排 |
| Z6 | C | extra_sum=0 | 4 | C | Spinor extra chars |
| Z7 | C | extra_sum=0 | 4 | C | Spinor extra chars |

### 核心问题与解决方案

#### 问题 1: 操作顺序不匹配 (Operation Ordering Mismatch)

**现象**: Wigner 测试给出非量子化 W=0.5，而非预期的 0, +1, -1。

**根因**: ISOTROPY 数据 (PIR/CIR) 和 spglib 数据库使用**不同的操作排列顺序**。
`find_seitz()` 返回的 `m.op_index` 是 spglib 顺序的索引，但 `h_chars[m.op_index]`
读取的是 ISOTROPY 顺序的特征标，两者不对应。

例如：SG 118 的 H[2] = {2_001|0,0,0} 在 spglib 中位于索引 2，但在 CIR 数据中
2_001 操作位于索引 1。`cir_chars[2]` 读到了错误操作的特征标。

**修复**: 通过匹配旋转矩阵建立 H_ops(spglib) → CIR/PIR(ISOTROPY) 索引映射。
- 生成阶段：从 PIR_data.txt 和 CIR_data.txt 提取操作旋转矩阵 (9 i32/op)
- 存储为 `PIR_ROTS` 和 `CIR_ROTS` 静态数组
- 运行时：`build_h_to_cir_map()` 通过旋转矩阵匹配建立映射
- `reorder_cir_chars()` 将字符表重排到 spglib 顺序
- Wigner 测试使用重排后的字符表，`m.op_index` 直接对应正确字符

#### 问题 2: CIR 解析器的 conjugate 矩阵跳过

**现象**: CIR 字符表包含错误值（如 dim=12 的 irrep 给出 χ(id)=1.0）。

**根因**: CIR_data.txt 中复不可约表示在 ordinary 复矩阵之后还有 **conjugate 复矩阵**。
旧解析器未跳过 conjugate 部分，导致后续操作读取了错误的数据行。

CIR_data.txt 格式 (每个操作):
```
[operator 4×4 matrix: 16 ints]
[ordinary complex matrix: dim² (re,im) pairs]
[conjugate complex matrix: dim² (re,im) pairs]  ← 仅 irtype==2 时存在
```

**修复**: 解析 header 中的 `irtype` 字段 (1=real, 2=complex)。当 irtype==2 时，
读取 ordinary 矩阵后额外跳过 dim² 个 token (conjugate 部分)。

#### 问题 3: 复合标签分解过于激进

**现象**: 251 个 compound irrep 的 CIR 字符和与 PIR 字符不匹配。

**根因**: `_decompose_compound_label("W1W2")` 纯基于字符串规则拆分为 ["W1", "W2"]，
但某些标签 (如 SG226 W1W2) 可能不是真的复合标签，或需要特殊处理 (同标签取共轭)。

**修复 (三层校验)**:
1. **Identity 字符校验**: `sum(CIR_χ_id) == PIR_χ_id` — 不满足则拒绝 (69 个被拒)
2. **同标签共轭**: `P3P3 = P3 ⊕ P3*` (第二份取共轭使虚部抵消)
3. **全字符表校验**: straight sum 和 conjugate sum 双候选，必须全操作匹配 PIR

结果: 251 → 41 → **0 mismatches**

#### 问题 4: Spinor Wigner 测试

**现象**: Z6/Z7 给出 W=0.5 (非量子化)，无法分类。

**根因**: Spinor (双值) 不可约表示需要 double-group 处理。每个空间操作 {R|t} 在
双群中有两个 lift: g 和 Ēg (Ē = 2π 旋转 = -1)。SeitzOp 只包含空间部分，
无法区分中心元，导致部分 Wigner 贡献项少乘了 -1。

**发现**: Bilbao spin.dat 文件中 BZ 边界 k 点的 irrep 行包含 **额外字符值** (extra chars)。
例如 Z6 的 16 个字符值中，后 8 个 `[0,0,0.5,-0.5,0,0,0,0]` 求和为 0，
恰好等于 Wigner indicator。分析 2141 个有 extra 的 spinor irrep 发现：
- 31% extra_sum = +1 (Type A)
- 18% extra_sum = 0 (Type C)
- 8.5% extra_sum = -1 (Type B)
- 其余为未量子化值

这表明 extra chars 是 Bilbao 预计算的 Wigner test contributions。

**修复**: 
1. 拆分 spin.dat 字符为标准字符 (前 n_lg) + extra (后 n_extra)
2. 存储 `SPIN_EXTRA_CHARS` 数组
3. `wigner_classify_spinor_extra()` 直接求和 → 分类
4. 优先使用 extra 通道，回退到 SU(2) 双群路径 (框架已搭建)

### Wigner 测试三条路径

| 路径 | 适用 | 关键函数 | 数据源 |
|------|------|---------|--------|
| CIR 路径 | Compound irreps (Z1Z4) | `wigner_classify_cir` | CIR_COMPONENT_CHARS + CIR_ROTS |
| PIR 路径 | Non-compound scalar (Z5) | `wigner_classify` | CHARACTERS + PIR_ROTS |
| Spinor extra | Spinor with extra (Z6/Z7) | `wigner_classify_spinor_extra` | SPIN_EXTRA_CHARS |
| Spinor SU(2) | Spinor w/o extra (fallback) | `wigner_classify_spinor` | SPIN_OP_ROTS/TRANS/SU2 |

### 严格 Wigner 分类

所有 Wigner 测试现在使用**严格量子化检查**:
- `|W-1| < 1e-6` → Type A
- `|W+1| < 1e-6` → Type B
- `|W| < 1e-6` → Type C
- 否则 → `CorepType::Unsupported`

### 新增/修改的源文件

| 文件 | 功能 |
|------|------|
| `src/irrep/wigner.rs` | Wigner 测试核心: Seitz 操作组合、Bloch 相位、分类函数 |
| `src/irrep/corep.rs` | Co-representation API: `compute_coreps()`, `CorepType`, `MagneticOps` |
| `src/irrep/types.rs` | `IrrepRecord` 新增: `_cir_start/count/ops`, `_pir_rot_start`, `_spin_*` 字段 |
| `scripts/parse_spinor_data.py` | 解析 Bilbao spin.dat: SU(2) lifts, op_indices, characters |

### 新增静态数组 (generated_data.rs)

| 数组 | 大小 | 说明 |
|------|------|------|
| `CIR_COMPONENT_CHARS` | ~27K f64 | Compound irrep 的 CIR 分量复字符 (re,im 对) |
| `CIR_ROTS` | ~104K i32 | CIR 操作旋转矩阵 (9/op)，与 CIR_COMPONENT_CHARS 同序 |
| `PIR_ROTS` | ~500K i32 | PIR 操作旋转矩阵 (9/op)，与 CHARACTERS 同序 |
| `SPIN_OP_ROTS` | ~40K i32 | Spinor 操作旋转矩阵 (9/op)，按 SG 索引 |
| `SPIN_OP_TRANS` | ~14K f64 | Spinor 操作平移矢量 (3/op) |
| `SPIN_OP_SU2` | ~19K f64 | Spinor SU(2) lift (4/op) |
| `SPIN_EXTRA_CHARS` | ~17K f64 | Pre-computed Wigner contributions for spinor irreps |
| `SPIN_OP_SG_INDEX` | 231 | SG# → (start, count) into SPIN_OP_* arrays |

### 数据生成关键修复 (generate_irrep_data.py)

1. **CIR 解析器**: 跳过 conjugate 复矩阵 (irtype==2 时)
2. **CIR_ROTS 对齐**: 当某 irrep 缺少旋转数据时零填充，保持与 CIR_COMPONENT_CHARS 同步
3. **PIR 旋转提取**: 从 PIR_data.txt 操作矩阵行提取旋转矩阵 → `PIR_ROTS`
4. **同标签共轭**: 自动检测并应用第二份取共轭
5. **全字符校验**: straight sum + conjugate sum 双候选，全操作匹配 PIR
6. **Spinor 旋转继承**: spinor irrep 从同 SG 同 k 点 scalar irrep 继承 PIR 旋转
7. **Spinor extra 拆分**: 按 op_indices 长度切分标准字符和 extra

### 重要教训

1. **永远不要假设数据顺序一致**: 不同数据源 (ISOTROPY vs spglib vs Bilbao) 的操作排列顺序不同，必须通过旋转矩阵匹配来建立映射
2. **字符串规则拆分不可靠**: `_decompose_compound_label()` 必须用字符表校验 (identity + full)，不能纯靠字符串匹配
3. **非量子化 W 必须报错**: 不能把 W=0.5 当做 Type A (W>0)，这是数据不完整的标志
4. **Bilbao extra chars 是预计算值**: 不需要自己实现 SU(2) 组合即可获得正确的 spinor Wigner 分类
5. **CIR 解析器边界处理**: 高维矩阵跨多行时，token-count 推进必须精确，否则逐操作累积偏移

