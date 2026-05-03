# cryspglib

A pure-Rust port of [spglib](https://github.com/spglib/spglib) (2.x) with zero external dependencies, **ported entirely with the assistance of AI (DeepSeek-V4-Pro, driven by Claude Code)**.

> **Disclaimer:** cryspglib is **not** intended to replace or compete with spglib. It exists to give Rust projects a zero-dependency, pure-Rust alternative that compiles without a C toolchain. It will continue to track upstream spglib as both evolve.

## Background

spglib is a widely-used C library in crystallography for identifying space groups, symmetry operations, and standardizing crystal cells. cryspglib is a line-by-line Rust port that preserves the same algorithms and data pipeline as the original C implementation.

Starting from the C source, the AI-assisted porting effort delivered:

- **30+ modules** ported from C to Rust (~15,000 lines of code)
- Automatic conversion of the **530 Hall symbol database** and **1,651 magnetic space group (UNI) entries**
- Complete identification pipelines for both **non-magnetic space groups** (230) and **magnetic space groups** (1,651 UNI numbers)
- Space group search, primitive cell reduction, Bravais lattice classification, Wyckoff position refinement, and more
- 70 unit tests + 1 doc test, covering cubic, hexagonal, tetragonal, monoclinic systems and multiple magnetic structures

## Features

### Non-magnetic space group identification

- Lattice + positions + types → space group (1–230), Hall symbol, international symbol
- Primitive cell finding (`spg_find_primitive`)
- Cell standardization (`spg_standardize_cell`)
- Wyckoff position refinement (`spg_refine_cell`)
- Delaunay / Niggli lattice reduction
- k-point grid generation (Monkhorst-Pack)

### Magnetic space group identification

- POSCAR-style input (lattice + positions + types + magnetic moments) → non-magnetic space group + magnetic space group (BNS) + symmetry operations
- Supports Type-1 (ordinary), Type-2 (grey), Type-3 (black-white / anti-rotation), Type-4 (anti-translation)
- Atomic magnetic moments treated as axial vectors (`is_axial=true`), with correct transformation under spatial inversion

### Verified physical systems

| System | Non-magnetic SG | Magnetic structure | Magnetic SG |
|--------|----------------|-------------------|-------------|
| Fe BCC AFM [111] | Im-3m (#229) | AFM along [111] | UNI=1331, BNS=166.101 |
| Fe SC [001] | Pm-3m (#221) | FM along [001] | UNI=1005, BNS=123.345 |
| FCC FM [001] | Fm-3m (#225) | FM along [001] | UNI=1005, BNS=123.345 |
| FCC FM [111] | Fm-3m (#225) | FM along [111] | UNI=1331, BNS=166.101 |
| CoF₃ AFM [001] | R-3c (#167) | Cr AFM alternating ±z | UNI=1333, BNS=167.103 |
| Graphene | P6/mmm (#191) | — | — |
| HCP | P6₃/mmc (#194) | — | — |
| Diamond | Fd-3m (#227) | — | — |
| Rocksalt | Fm-3m (#225) | — | — |
| CrPS₄ | C2 (#5) | — | — |
| La₂NiO₄ | P4₂/ncm (#138) | — | — |

## Usage

### Non-magnetic

```rust
use cryspglib::Crystal;

// FCC Al (space group Fm-3m, #225)
let al = Crystal::new(
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]],
    vec![13, 13, 13, 13],
);
let ds = al.analyze().symprec(1e-5).dataset()?;
// → SpaceGroup { spacegroup_number: 225, international_symbol: "Fm-3m", ... }
# Ok::<(), cryspglib::SymError>(())
```

### Magnetic

```rust
use cryspglib::Crystal;

// BCC AFM [111]: Fe at [0,0,0] and [0.5,0.5,0.5], opposite spins along [111]
let n = (3.0_f64).sqrt();
let fe = Crystal::new(
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    vec![26, 26],
).with_magnetic(vec![
    [1.0/n, 1.0/n, 1.0/n],
    [-1.0/n, -1.0/n, -1.0/n],
]);

let result = fe.analyze().symprec(1e-5).magnetic_dataset().unwrap();
// → MagneticSymmetry {
//     spacegroup_number: 229 (Im-3m),
//     uni_number: 1331,
//     bns_number: "166.101",
//     magnetic_type: BlackWhite (Type-3),
//     num_operations: 24,
//     ...
//   }
```

### Key API

| Type / Method | Description |
|---------------|-------------|
| `Crystal::new(lat, pos, types)` | Create a non-magnetic 3D crystal |
| `.with_magnetic(moments)` | Add magnetic moments `[mx,my,mz]` per atom |
| `.with_layer(axis)` | Mark as 2D slab |
| `.analyze()` | Begin symmetry analysis |
| `.delaunay_reduce(prec)` | Delaunay lattice reduction |
| `.niggli_reduce(prec)` | Niggli lattice reduction |
| `Crystal::from_poscar(data)` | Parse POSCAR-format input |
| `SymmetryAnalysis::symprec(val)` | Set symmetry tolerance |
| `.dataset()` | Full space group dataset |
| `.symmetry()` | Symmetry operations only |
| `.primitive_cell()` | Primitive cell |
| `.standardize(to_prim, no_ideal)` | Standardized cell |
| `.hall_number()` | Hall number (1-530) |
| `.international()` | Space group number + symbol |
| `.magnetic_dataset()` | Magnetic space group |
| `.irreducible_mesh(mesh, shift, tr)` | Irreducible k-point grid |
| `SpaceGroupType::from_hall(n)` | Look up space group type by Hall number |
| `MagneticSpaceGroupType::from_uni(n)` | Look up magnetic SG type by UNI number |

## Build

```bash
# This crate is part of a workspace; run from workspace root:
cd /home/liuyichen/TB_rs
cargo build --package cryspglib
cargo test --package cryspglib   # 57 unit + 13 integration + 5 doc tests
```

## Relationship to upstream spglib

- cryspglib is **not a fork or a replacement** — it is a companion port that aims to make spglib's functionality available to the Rust ecosystem without requiring a C compiler or FFI bindings
- Algorithms and databases are kept identical to spglib 2.x and will track upstream updates
- Naming conventions preserve the `spg_` prefix with snake_case variants for easy migration
- Rust type safety: `Option<T>` replaces NULL returns, `Vec` replaces malloc/free, automatic memory management

## Acknowledgments

~95% of the code in this project was generated by **DeepSeek-V4-Pro** via Claude Code. The porting process involved approximately 50+ conversations and hundreds of iterative corrections, with the AI handling everything from line-by-line translation to logic debugging to test authoring.

## License

BSD-3-Clause (inherited from spglib)
