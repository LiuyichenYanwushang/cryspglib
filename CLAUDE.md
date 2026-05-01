# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

# REFACTORING: Rust-idiomatic API (v0.2.0)

> **Status**: IN PROGRESS. Replace C-style `spg_*` free functions with a clean, struct-based Rust API.

## 1. Current Problems

| Problem | Example |
|---------|---------|
| C-style prefix | `spg_get_dataset()`, `spg_get_symmetry()` — 40+ functions |
| Parameter soup | Every function repeats `(lattice, position, types, symprec)` |
| Pre-allocated output | `spg_get_ir_reciprocal_mesh(grid_address: &mut [..])` |
| Duplicate API surface | `spg_*` + `spgat_*` double every function for angle_tolerance |
| Mixed concerns | `spg_read_structure`, `element_to_number` in lib.rs |
| Version getter spam | 6 functions for version constants |

## 2. New API Design

### 2.1 Core struct: `Crystal`

Replaces bare `(lattice, position, types, tensors?)` tuples:

```rust
pub struct Crystal {
    lattice: Mat3,              // [cart][vec] convention unchanged
    positions: Vec<[f64; 3]>,   // fractional coordinates
    types: Vec<i32>,            // atomic numbers
    tensors: Option<Vec<f64>>,  // magnetic moments (None = non-magnetic)
    aperiodic_axis: Option<AperiodicAxis>, // None = 3D bulk
}

impl Crystal {
    pub fn new(lattice: Mat3, positions: Vec<[f64; 3]>, types: Vec<i32>) -> Self;

    // Builder-style setters
    pub fn with_magnetic(mut self, moments: Vec<[f64; 3]>) -> Self;
    pub fn with_layer(mut self, axis: AperiodicAxis) -> Self;
}
```

### 2.2 Analysis builder: `Crystal::analyze()`

```rust
impl Crystal {
    pub fn analyze(&self) -> SymmetryAnalysis<'_> { ... }
}

pub struct SymmetryAnalysis<'a> {
    crystal: &'a Crystal,
    symprec: f64,           // default 1e-5
    angle_tolerance: f64,   // default -1.0 (auto)
    hall_number_hint: Option<i32>,
}
impl<'a> SymmetryAnalysis<'a> {
    pub fn symprec(mut self, val: f64) -> Self;
    pub fn angle_tolerance(mut self, val: f64) -> Self;
    pub fn hall_number(mut self, val: i32) -> Self;

    // Terminal methods — consume the builder, return result
    pub fn dataset(&self) -> Result<SpglibDataset, SpglibError>;
    pub fn symmetry(&self) -> Result<SymmetryOps, SpglibError>;
    pub fn primitive_cell(&self) -> Result<Crystal, SpglibError>;
    pub fn standardize(&self, to_primitive: bool, no_idealize: bool) -> Result<Crystal, SpglibError>;
    pub fn irreducible_mesh(&self, mesh: [i32; 3], is_shift: [i32; 3], time_reversal: bool) -> Result<IrMesh, SpglibError>;
}
```

### 2.3 Proper types instead of parallel arrays

```rust
/// A single symmetry operation {R|t}
pub struct SymmetryOp {
    pub rotation: [[i32; 3]; 3],  // integer rotation matrix
    pub translation: [f64; 3],    // fractional translation
    pub time_reversal: bool,      // for magnetic symmetry
}

/// Ordered set of symmetry operations
pub struct SymmetryOps {
    pub operations: Vec<SymmetryOp>,
}

/// Irreducible k-point mesh
pub struct IrMesh {
    pub grid_addresses: Vec<[i32; 3]>,
    pub mapping_table: Vec<usize>,   // full grid → irreducible map
    pub num_ir: usize,
}
```

### 2.4 What to REMOVE from public API

- All `spg_get_*` free functions → replaced by `Crystal::analyze().xxx()`
- All `spgat_*` variants → replaced by `.angle_tolerance()` builder
- `spg_get_version()`, `spg_get_major_version()`, etc. → single `pub const VERSION: &str`
- `spg_get_error_message()` → `thiserror::Error` already provides `Display`

### 2.5 What STAYS (re-exported from lib)

- `SpglibError` enum (already Rust-idiomatic)
- `SpglibDataset` struct — but rename to `SpaceGroup`?
- `SpglibMagneticDataset` → rename to `MagneticSpaceGroup`?
- `MagneticType` enum
- `AperiodicAxis` enum
- `TensorRank` enum
- `spg_read_structure` → move to `parser` module, rename to `Crystal::from_poscar(data)`
- `spg_format_magnetic_symmetry` → `impl Display for MagneticSymmetry`

### 2.6 Usage comparison

**Before (C-style):**
```rust
let ds = spg_get_dataset(&lattice, &positions, &types, 1e-5)?;
let prim = spg_find_primitive(&lattice, &positions, &types, 1e-5)?;
let sym = spg_get_symmetry(&lattice, &positions, &types, 1e-5)?;
```

**After (Rust-idiomatic):**
```rust
let crystal = Crystal::new(lattice, positions, types);
let ds = crystal.analyze().symprec(1e-5).dataset()?;
let prim = crystal.analyze().symprec(1e-5).primitive_cell()?;
let sym = crystal.analyze().symprec(1e-5).symmetry()?;
```

## 3. Migration Steps (ordered)

### Step 1: `Crystal` struct + `SymmetryAnalysis` builder
- Create `Crystal` in a new top-level module (or in `cell.rs`)
- Implement `SymmetryAnalysis` with all builder methods
- Wire up `dataset()`, `symmetry()` to call existing internal functions
- Keep old `spg_*` functions working (delegate to new API internally)

### Step 2: New output types
- `SymmetryOp`, `SymmetryOps` — replace parallel `rotations`/`translations` arrays
- `IrMesh` — owned output for k-point grid
- Update internal pipeline to use new types

### Step 3: Move parser
- `spg_read_structure` → `crate::parser::from_poscar(data) -> Option<Crystal>`
- `element_to_number` → `crate::parser` (private)

### Step 4: Clean up public API
- Remove old `spg_*`/`spgat_*` free functions
- Remove version getter functions, use `const VERSION`
- Remove `spg_get_error_message` (thiserror already provides Display)

### Step 5: Rename for clarity
- `SpglibDataset` → `SpaceGroup` (the `Spglib` prefix is redundant)
- `SpglibMagneticDataset` → `MagneticSpaceGroup`
- All names drop the `Spglib`/`spg` prefix

### Step 6: Update all tests and docs

## 4. Risk Points

- **Type alias `Mat3`** — is `[[f64; 3]; 3]`, fine for 3D but cannot express 1D/2D. This is a pre-existing limitation, not introduced by this refactor.
- **Internal code uses `Cell` struct** heavily — we keep `Cell` as an internal type, `Crystal` is the public-facing wrapper that converts to/from `Cell`.
- **Backward compatibility** — old `spg_*` functions can be kept as `#[deprecated]` shims during migration, then removed.
- **`TensorRank`** — currently has `NoSpin = -1`, `Collinear = 0`, `NonCollinear = 1`. This maps to C enums and may need cleanup (negative discriminant is unusual in Rust).

## 5. Naming Convention

| Old | New |
|-----|-----|
| `spg_get_dataset` | `crystal.analyze().dataset()` |
| `spg_get_symmetry` | `crystal.analyze().symmetry()` |
| `spg_find_primitive` | `crystal.analyze().primitive_cell()` |
| `spg_standardize_cell` | `crystal.analyze().standardize(to_prim, no_ideal)` |
| `spg_refine_cell` | `crystal.analyze().standardize(false, false)` |
| `spg_get_international` | `crystal.analyze().dataset()?.spacegroup_number` |
| `spg_get_schoenflies` | `crystal.analyze().dataset()?.schoenflies` |
| `spg_get_ir_reciprocal_mesh` | `crystal.analyze().irreducible_mesh(mesh, shift, tr)` |
| `spg_delaunay_reduce` | `crystal.delaunay_reduce(symprec)` |
| `spg_niggli_reduce` | `crystal.niggli_reduce(symprec)` |

---

## Build & Test Commands

This crate is part of a workspace at `/home/liuyichen/TB_rs/`. Run all commands from the workspace root with `--package cryspglib`:

```bash
cd /home/liuyichen/TB_rs
cargo build --package cryspglib          # Build the library
cargo test --package cryspglib           # Run all tests (57 unit + 13 integration + 5 doctests)
cargo test --package cryspglib -- --nocapture  # Run tests with debug output
cargo test --package cryspglib <test_name>     # Run a single test by name
cargo test --package cryspglib --test cof3     # Run a specific integration test file
cargo check --package cryspglib          # Check compilation without building
cargo clippy --package cryspglib         # Run lints
```

The project uses `edition = "2024"` (stable since Rust 1.85). No nightly required.

## Critical Convention: Matrix Layout `lattice[cart][vec]`

All 3x3 matrices (`Mat3 = [[f64; 3]; 3]`) use the convention **rows = Cartesian components, columns = lattice vectors**:

```
lattice = [[a_x, b_x, c_x],
           [a_y, b_y, c_y],
           [a_z, b_z, c_z]]
```

- `lattice[i][j]` — the i-th Cartesian component (0=x, 1=y, 2=z) of the j-th lattice vector (0=a, 1=b, 2=c)
- Each **column** is a lattice vector; each **row** is a Cartesian direction

**Why this matters:**

- `mat_get_metric(lattice) = L^T * L` computes the metric tensor G where G[i][j] = dot(vec_i, vec_j). This only works when columns are vectors.
- `mat_multiply_matrix_vector_d3(lattice, frac)` transforms fractional coordinates to Cartesian: `cart = L * frac`
- Cubic matrices are symmetric (a=b=c, all angles 90°), so both conventions look identical — this is why cubic tests passed before the convention was corrected
- Hexagonal matrices are asymmetric (`a != b` direction), so swapping conventions silently produces wrong results (e.g., graphene identified as space group #10 instead of #191)

**When constructing test lattices**, always use `lattice[cart][vec]`:
```rust
// Hexagonal: a along x, b at 120° in xy-plane, c along z
fn hexagonal_lattice(a: f64, c: f64) -> Mat3 {
    let sqrt3 = 3.0_f64.sqrt();
    [
        [a, -a / 2.0, 0.0],           // x components of a, b, c
        [0.0, a * sqrt3 / 2.0, 0.0],  // y components
        [0.0, 0.0, c],                // z components
    ]
}
```

## Architecture: Space Group Identification Pipeline

This is a pure Rust port of the C library `spglib` with a single dependency (`thiserror`).

The pipeline for identifying a crystal's space group:

1. **Cell** (`cell.rs`) — crystal structure input: lattice vectors, atomic positions, types, overlap checks
2. **Symmetry detection** (`symmetry.rs`) — find symmetry operations (rotations + translations) of the cell
3. **Point group classification** (`pointgroup.rs`) — classify the rotation part into a crystallographic point group, determine Laue class and holohedry
4. **Lattice reduction** (`delaunay.rs`, `niggli.rs`) — reduce lattice to canonical form (Delaunay for general, Niggli for triclinic/monoclinic)
5. **Primitive cell** (`primitive.rs`) — find primitive cell and its symmetry operations
6. **Hall symbol matching** (`hall_symbol.rs`) — match primitive cell symmetry against the Hall symbol database (`spg_database.rs`) to identify the Hall number (1–530)
7. **Space group search** (`spacegroup.rs`) — iterate Hall number candidates, perform origin shift and conventional cell transformation, select best match; returns the space group number (1–230)
8. **Structure refinement** (`refinement.rs`) — rotate Bravais lattice to standard orientation, determine exact Wyckoff positions via site symmetry database, recover symmetry operations in the original cell. Produces `ExactStructure` with standardized positions, Wyckoff labels, equivalent atom mappings.
9. **Determination** (`determination.rs`) — top-level retry loop that combines primitive + spacegroup + refinement steps with decreasing tolerance
10. **Public API** (`lib.rs`) — wraps the pipeline behind a clean Rust API: `SpglibDataset`, `spg_get_dataset`, `spg_standardize_cell`, etc.

Key supporting modules:
- `mathfunc.rs` — 3x3 matrix/vector operations with type aliases `Mat3`, `Mat3I`, `Vec3`, `Vec3I` (no external BLAS/LAPACK)
- `arithmetic.rs` — integer arithmetic for symmetry checks
- `overlap.rs` — overlap detection for atomic positions
- `site_symmetry.rs` — Wyckoff position exact location and equivalent atom assignment using site symmetry database
- `sitesym_database.rs` — site symmetry database (Wyckoff positions, site symmetry symbols)
- `spin.rs` — magnetic symmetry types
- `kgrid.rs`, `kpoint.rs` — k-point grid generation (Monkhorst-Pack)
- `debug.rs` — debug/warning print macros

## Error Handling: `SpglibError` + `thiserror`

```rust
#[derive(thiserror::Error, Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpglibError {
    Success,                      // no error (used internally)
    SpacegroupSearchFailed,       // no matching space group found
    CellStandardizationFailed,    // could not standardize cell
    SymmetryOperationSearchFailed,// symmetry detection failed
    AtomsTooClose,                // overlapping atoms beyond tolerance
    PointgroupNotFound,           // rotation part doesn't match any point group
    NiggliFailed,                 // Niggli reduction did not converge
    DelaunayFailed,               // Delaunay reduction did not converge
    ArraySizeShortage,            // internal buffer too small
    InvalidInput,                 // malformed input (e.g., POSCAR parse error)
    MathFailed,                   // singular matrix, non-invertible, etc.
}
```

All fallible public APIs return `Result<T, SpglibError>`. Internal functions that can fail in multiple ways (e.g., `mat_inverse_matrix_d3`, `mat_get_similar_matrix_d3`) also return `Result`.

## API Naming Conventions

| Pattern | Meaning |
|---------|---------|
| `spg_<name>` | Standard API with default angle tolerance (`-1.0` = auto) |
| `spgat_<name>` | Same API with explicit `angle_tolerance` parameter |
| `spg_<name>_with_hall_number` | Returns dataset + Hall number hint for magnetic matching |
| `spg_get_layer_<name>` | 2D slab variant (one aperiodic axis, periodicity in the other two) |

## `TensorRank` and `AperiodicAxis`

`TensorRank` controls how magnetic moments (or site tensors) are treated:

- `TensorRank::Scalar` — non-magnetic (no spin)
- `TensorRank::Collinear` — collinear magnetic (spin up/down along one axis)
- `TensorRank::NonCollinear` — non-collinear magnetic (arbitrary 3D spin vectors)

`AperiodicAxis` marks which axis lacks periodicity in 2D slab systems (`X`, `Y`, or `Z`). Layer functions accept `Option<AperiodicAxis>` where `None` means full 3D periodicity.

## Test Coverage

- **Unit tests** — `#[cfg(test)] mod tests` inside source files, can access private APIs
- **Integration tests** — `tests/*.rs` files, use only the public API (real material tests)

**Unit tests in source files:**

| Module | What's tested |
|--------|--------------|
| `spacegroup.rs` | Cubic (#221/229/225/227), hexagonal (#191/194/164), supercells |
| `magnetic_spacegroup.rs` | DB matching Type-1/2/3, edge cases (empty sym, missing identity) |
| `mathfunc.rs` | Matrix inverse, multiply, det, trace, similar, nint, dabs |
| `kgrid.rs`, `kpoint.rs` | Grid generation, irreducible mesh, mesh symmetry |
| `cell.rs` | Cell creation, overlap detection, trim |
| `delaunay.rs` | Delaunay reduction, extended basis |
| `niggli.rs` | Niggli reduction identity/swap |
| `overlap.rs` | Overlap checker, layer overlap |
| `symmetry.rs` | Pure translation detection, supercell |
| `debug.rs`, `sitesym_database.rs`, `arithmetic.rs`, `determination.rs` | Smoke/correctness |

**Integration tests in `tests/`:**

| File | System | Space group | Notes |
|------|--------|------------|-------|
| `tests/cof3.rs` | CoF₃ | R-3c (#167) | Non-magnetic + AFM [001] magnetic (UNI=1333) |
| `tests/crps4.rs` | CrPS₄ | C2 (#5) | Monoclinic, 48 atoms, C-centering |
| `tests/la2nio4.rs` | La₂NiO₄ | P4₂/ncm (#138) | Tetragonal, 28 atoms, D₄ₕ |
| `tests/magnetic_integration.rs` | Fe SC/BCC/FCC | Various | Type-1/2/3 magnetic, FM/AFM, public API end-to-end |

## Magnetic Space Group: Complete Pipeline

End-to-end: **POSCAR-like input → non-magnetic space group + magnetic space group (BNS) + symmetry operations**.

### Pipeline Flow

```
POSCAR → Cell (TensorRank::NonCollinear)
  → prm_get_primitive + spa_search_spacegroup → 非磁空间群
  → sym_get_operation(&cell)                  → 非磁对称操作（常规晶胞基）
  → spn_get_operations_with_site_tensors      → 磁对称操作（含 timerev 标记）
  → reduce_to_primitive_magsym                → 纯平移约化到原胞表示
  → msg_identify_magnetic_space_group_type    → 磁空间群 UNI + BNS
```

### Critical Design Rules

1. **Must use `sym_get_operation(&cell)` for non-magnetic symmetry**, not `sym_get_operation` on the primitive cell. Using primitive cell basis to map conventional cell positions gives wrong rotation mapping (FCC drops from 64 to 4 ops).

2. **Magnetic matching uses subset, not equality.** Magnetic moments may break some symmetries, so `changed_symmetry` may have fewer operations than the database. Use `is_subset` instead of `is_equal` when matching against DB.

3. **`reduce_to_primitive_magsym` is mandatory** before calling `msg_identify_magnetic_space_group_type`. The non-magnetic symmetry operations come in the conventional cell representation (e.g., FCC 4-atom cell has 64 ops with 4 centering translations), but the magnetic DB uses primitive cell representation (16 ops). This function detects pure translations (identity rotation) and reduces operations to shortest translations with deduplication.

Key rules for lattice translation collection:
- timerev=0 identity ops → always included as reduction translations
- timerev=1 identity ops with **non-zero** translation (e.g., BCC (0.5,0.5,0.5) with spin flip in AFM) → included, since they represent centering translations that carry spin reversal
- timerev=1 identity ops with **zero** translation → excluded, as these are part of Type-2 (Grey) group structure, not lattice translations

When applying such a translation `(lt, ltr)`: `(R, t, tr) → (R, t - lt, tr XOR ltr)`.

4. **Always call `match_hall_symbol_db`** (not `hal_match_hall_symbol_db` directly) in `search_hall_number`. The intermediate layer does: (a) point-group filter to skip unmatching Hall numbers, (b) centering expansion for RCenter (expands 12 primitive ops → 36 with centering translations). Without expansion, all 7 rhombohedral Hall numbers fail matching because DB ops (36) can't match symmetry ops (12). Other centerings are pre-expanded by `get_initial_conventional_symmetry`.

5. **UNI matching must prioritize parent Hall number.** After primitive reduction, `changed_symmetry` has the same operations regardless of centering. Since `is_subset` allows matching both, the result is non-unique. The parent (non-magnetic) Hall number breaks this tie—the magnetic space group must be a subgroup of the parent. `msg_identify_with_parent_hall` tries `parent_hall_number` first for UNI candidates.

## Notes

- `debug_test.rs` and similar ad-hoc files in the repo root are gitignored/untracked — use them for quick experiments but don't commit.
- Integration tests in `tests/` use only the public API. When adding new material tests, prefer `tests/` over `src/`.
- Internal DB-level or edge-case tests go in the corresponding source file's `#[cfg(test)] mod tests`.
