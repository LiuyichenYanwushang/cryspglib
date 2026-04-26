# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test Commands

```bash
rustup default nightly    # Required: edition 2024 needs nightly Rust
cargo build              # Build the library
cargo test               # Run all 70 tests + 1 doctest
cargo test -- --nocapture  # Run tests with debug output (eprintln visible)
cargo test <test_name>   # Run a single test (e.g., `cargo test test_graphene_p6mmm`)
cargo check              # Check compilation without building
cargo clippy             # Run lints
```

The project uses `edition = "2024"` (requires nightly Rust). All commands must use nightly toolchain.

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
        [a, -a / 2.0, 0.0],      // x components of a, b, c
        [0.0, a * sqrt3 / 2.0, 0.0],  // y components
        [0.0, 0.0, c],            // z components
    ]
}
```

**When reading C reference code**, note that C uses the same memory layout (`double lattice[3][3]` with first index as row/cartesian), but the C code's array initialization may be transposed relative to Rust's `[[f64; 3]; 3]` depending on how it's written.

## Architecture: Space Group Identification Pipeline

This is a pure Rust port of the C library `spglib` with zero external dependencies.

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

## Critical Porting Pattern: C Output Parameters → Rust Return Values

In the C source (`src_c/*.txt`), many functions write results through pointer parameters:
```c
void mat_cast_matrix_3i_to_3d(double out[3][3], int const a[3][3]);
```

The Rust port converts these to return values:
```rust
pub fn mat_cast_matrix_3i_to_3d(a: &Mat3I) -> Mat3;
```

**When comparing Rust against C code, always check that the return value is actually captured.** The most severe bug found in this port was ignoring the return value of `mat_cast_matrix_3i_to_3d` in `get_conventional_symmetry`, which caused all rotation matrices to be zero and silently broke space group identification.

## Rust API Design Conventions (C → Rust mapping)

The `lib.rs` public API follows these conventions when porting from C:

| C Pattern | Rust Pattern |
|-----------|-------------|
| `malloc`/`free` | `Vec` (ownership, auto-drop) |
| `NULL` return on failure | `Option<T>` (None = failure) |
| `thread_local` error code + return 0 | `Option<T>` return |
| Output parameter `double out[3][3]` | Return `Mat3` |
| Output parameter `int *types` + `int *num_atom` | Return `Cell` (owns lattice, positions, types) |
| `spg_free_dataset(dataset)` | Not needed — Rust drops `SpglibDataset` automatically |
| `char[7]` / `char[11]` buffers | `String` |
| `int *mapping_table` pre-allocated by caller | `Vec<i32>` returned |
| K-point pre-allocated grid arrays | `&mut [i32]` / `&mut [usize]` slices |

## Reference C Code

`src_c/` contains the original spglib C source (`*.txt` files, renamed from `*.c` for tooling compatibility). Most `.rs` files correspond to a `.txt` file of the same name. The magnetic modules (`msg_database`, `magnetic_spacegroup`) are partially ported.

| C file | Rust module | Purpose |
|--------|------------|---------|
| `spglib.txt` | `lib.rs` | Public API, dataset, standardization |
| `determination.txt` | `determination.rs` | Top-level retry loop |
| `refinement.txt` | `refinement.rs` | Bravais lattice refinement, Wyckoff positions |
| `site_symmetry.txt` | `site_symmetry.rs` | Site symmetry exact positions |
| `spacegroup.txt` | `spacegroup.rs` | Space group search and identification |
| `primitive.txt` | `primitive.rs` | Primitive cell finding |
| `hall_symbol.txt` | `hall_symbol.rs` | Hall symbol matching |
| `pointgroup.txt` | `pointgroup.rs` | Point group classification |
| `spg_database.txt` | `spg_database.rs` | Space group database |
| `sitesym_database.txt` | `sitesym_database.rs` | Site symmetry database |
| `mathfunc.txt` | `mathfunc.rs` | Matrix/vector math |
| `delaunay.txt` | `delaunay.rs` | Delaunay lattice reduction |
| `niggli.txt` | `niggli.rs` | Niggli lattice reduction |
| `kgrid.txt` | `kgrid.rs` | k-point grid generation |
| `kpoint.txt` | `kpoint.rs` | k-point address mapping |
| `symmetry.txt` | `symmetry.rs` | Symmetry operation search |
| `cell.txt` | `cell.rs` | Crystal cell and overlap checks |
| `msg_database.txt` | `msg_database.rs` + `msg_database_gen.rs` | Magnetic space group database (auto-converted from C) |
| `magnetic_spacegroup.txt` | `magnetic_spacegroup.rs` | Magnetic space group identification |
| `spin.txt` | `spin.rs` | Magnetic symmetry operations with site tensors |

When debugging, compare Rust implementation against the C original line-by-line, paying special attention to:
- Output parameter → return value conversions
- Array indexing differences (C is row-major, Rust uses same layout)
- Floating point tolerance constants

## Test Coverage

Tests are inline (`#[cfg(test)] mod tests` in most `.rs` files) and in dedicated test files.

**Non-magnetic space group tests** (in `spacegroup.rs`):
- Cubic: simple (#221), bcc (#229), fcc (#225), diamond (#227), rocksalt (#225)
- Hexagonal: graphene (#191), hcp (#194), AB-buckled silicene (#164)
- Supercells: 2x2x2 simple cubic (#221), 2x2x2 CsCl (#221), 2x2x1 graphene (#191)

**Magnetic space group tests** (in `magnetic_spacegroup_test.rs`):
- Type-1/2/3 identification via Pm-3m (#221) symmetry operations from database
- FCC FM [001] and [111] end-to-end tests
- Public API `spg_get_magnetic_spacegroup_type_from_symmetry`
- Edge cases: empty symmetry, missing identity

**CoF3 tests** (in `cof3_test.rs`):
- Non-magnetic: R-3c (#167) with 24 atoms (18 F + 6 Cr)
- Magnetic: Cr AFM [001] alternating ±z → Type-1, UNI=1333, BNS=167.103
- POSCAR test files in `test/CoF3/` (POSCAR, BPOSCAR, PPOSCAR)

**CrPS₄ tests** (in `crps4_test.rs`):
- Monoclinic C2 (#5) with 48 atoms (8 Cr + 8 P + 32 S)
- Primitive cell validation against phonopy PPOSCAR (12 atoms)

**La₂NiO₄ tests** (in `la2nio4_test.rs`):
- Tetragonal P4₂/ncm (#138) with 28 atoms (8 La + 4 Ni + 16 O)

**Unit tests** for individual modules: `kgrid.rs`, `kpoint.rs`, `cell.rs`, `delaunay.rs`, `debug.rs`, `overlap.rs`, `symmetry.rs`.

## Magnetic Space Group: Complete Pipeline

End-to-end: **POSCAR-like input → non-magnetic space group + magnetic space group (BNS) + symmetry operations**.

### API

```rust
pub fn spg_get_magnetic_dataset(
    lattice, positions, types, magnetic_moments, symprec
) -> Option<SpglibMagneticSymmetry>;
pub fn spg_read_structure(data) -> Option<(Mat3, Vec<Vec3>, Vec<i32>, Option<Vec<[f64;3]>>)>;
pub fn spg_format_magnetic_symmetry(result) -> String;
```

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

4. **Always call `match_hall_symbol_db`** (not `hal_match_hall_symbol_db` directly) in `search_hall_number`. The intermediate layer does: (a) point-group filter to skip unmatching Hall numbers, (b) centering expansion for RCenter (expands 12 primitive ops → 36 with centering translations `(2/3,1/3,1/3)` and `(1/3,2/3,2/3)`). Without expansion, all 7 rhombohedral Hall numbers fail matching because DB ops (36) can't match symmetry ops (12). Other centerings are pre-expanded by `get_initial_conventional_symmetry`.

5. **UNI matching must prioritize parent Hall number.** After primitive reduction, `changed_symmetry` has the same operations regardless of centering (e.g., P-3c1 reduced=12 ≡ R-3c reduced=12). Since `is_subset` allows matching both, the result is non-unique (e.g., CoF3 AFM [001] matches both 165.91 Hall=457 and 167.103 Hall=460). The parent (non-magnetic) Hall number breaks this tie—the magnetic space group must be a subgroup of the parent. `msg_identify_with_parent_hall` now tries `parent_hall_number` first for UNI candidates, only falling back to the FSG-discovered Hall number if the parent gives no match.

## Next Task: C-Style `&mut` Refactoring

Many Rust functions still use C-style patterns that don't need `&mut`:

1. **`&mut` output parameter pattern:** `fn foo(&mut out: &mut Mat3, input...)` — these should return `Mat3` instead. When a function takes `&mut` solely to write a result (not to modify an existing value in place), convert it to `-> ReturnType`.

2. **`&mut` for temporary mutation:** Some functions take `&mut` for a local variable that gets mutated and never read after. These can use `let mut` locally and return by value.

3. **`Vec` vs pointer pre-allocation:** The C code pre-allocates arrays passed by pointer. Rust ports should return `Vec<T>` instead of taking `&mut Vec<T>` and pushing into it.

### Systematic approach

For each file in `src/`, audit every `fn.*&mut` signature and classify:
- **Needs `&mut`**: genuinely mutates external state (e.g., `Cell`, `Spacegroup` fields)
- **Should return by value**: output-only `&mut` parameters → convert to return
- **Borderline**: flag for discussion

Priority files (heavy C-style porting):
- `mathfunc.rs` — many matrix functions with `&mut` out params
- `primitive.rs` — lattice vector search with pre-allocated output
- `spacegroup.rs` — contains `spa_copy_spacegroup` etc.
- `symmetry.rs`, `refinement.rs`

### Example transformation

```rust
// C-style (current):
pub fn foo(out: &mut Mat3, input: &Mat3, symprec: f64) { ... }

// Rust-style (target):
pub fn foo(input: &Mat3, symprec: f64) -> Option<Mat3> { ... }
```

## Notes

- `debug_test.rs` and similar ad-hoc files in the repo root are gitignored/untracked — use them for quick experiments but don't commit.
- `periodic_axes` in `cell.rs` layer functions uses `[usize; 2]` through the `AperiodicAxis` enum (not raw `i32`).
