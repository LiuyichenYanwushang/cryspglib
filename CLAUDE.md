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

```
isotropy_subgroup/
├── iso.zip/data_irreps.txt    → 元数据 (标签, SG, image, isotropy)
├── PIR_data.zip/PIR_data.txt  → 实表示矩阵 (PIR)
└── CIR_data.zip/CIR_data.txt  → 复表示矩阵 (CIR, 回退用)

                    ↓ scripts/generate_irrep_data.py ↓

src/irrep/generated_data.rs   → 5 个静态数组, ~8.7 MB, 编译进二进制
```

### Rust 源文件

| 文件 | 功能 |
|------|------|
| `src/irrep/mod.rs` | 模块入口, 230 SG 索引表 + rustdoc |
| `src/irrep/types.rs` | `IrrepRecord`, `IsotropyRecord`, `IrrepData`, `KPointData` |
| `src/irrep/query.rs` | `irreps_of()`, `kpoints_of()`, `format_character_table()`, `symmetry_operations_of()` |
| `src/irrep/generated_data.rs` | 自动生成: `IRREPS[4777]`, `ISOTROPY_SUBGROUPS[15239]`, `CHARACTERS[f64]`, `MATRICES[f64]`, `SG_IRREP_INDEX[231]` |
| `src/irrep/preamble.rs` | 230 SG 的 HM 符号 + Schoenflies 符号 |
| `src/irrep/{triclinic,monoclinic,...,cubic}.rs` | 7 个晶系文件的 rustdoc 页面 |

### 静态数组 (全部在 generated_data.rs)

| 数组 | 大小 | 说明 |
|------|------|------|
| `SG_IRREP_INDEX` | 231 | SG# → (start, count) into IRREPS |
| `IRREPS` | 4,777 | 所有 irreps (sg, ml, bc, kov, dim, image, kvec, pointers) |
| `ISOTROPY_SUBGROUPS` | 15,239 | 所有 isotropy subgroups (sg, symbol, direction, domains, arms) |
| `CHARACTERS` | ~50K f64 | 所有特征标值, 100% 覆盖 |
| `MATRICES` | ~580K f64 | 所有矩阵元素, 100% 覆盖 |

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
- `domains` 和 `arms` 是 `u16` (不是 u8), 因为值可达 384

