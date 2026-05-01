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

## Remaining TODOs (in-progress, may be broken mid-refactor)

### 1. Fix test imports (CURRENTLY BROKEN)
Tests (`tests/cof3.rs`, `tests/crps4.rs`, `tests/la2nio4.rs`) have broken imports
after migration to new API. Need to restore `use cryspglib::{Crystal, SpaceGroupType,
spg_get_pointgroup, ...}`.

### 2. Migrate magnetic tests to new API
- `spg_get_magnetic_dataset(lattice, pos, types, moments, prec)` →
  `Crystal::new(lat, pos, types).with_magnetic(moments).analyze().symprec(prec).magnetic_dataset()`
- `spg_format_magnetic_symmetry(&result)` → `result.to_string()` (Display impl added)
- `magnetic_integration.rs` still uses old API exclusively

### 3. Migrate `spg_get_magnetic_spacegroup_type_from_symmetry`
Add as method on `MagneticSymmetry` or `SymmetryAnalysis`.

### 4. Add database lookup methods
- `SymmetryOps::from_database(hall_number)` for `spg_get_symmetry_from_database`
- Keep `spg_get_pointgroup(&rotations)` as standalone (pure math utility)

### 5. Clean up old k-point functions
`spg_get_ir_reciprocal_mesh` etc. still public. Already wrapped by
`SymmetryAnalysis::irreducible_mesh()`, but old functions not yet deleted.

### 6. Clean up old `spg_*` k-point BZ functions
`spg_get_dense_BZ_grid_points_by_rotations`, `spg_relocate_BZ_grid_address` etc.

