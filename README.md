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

```rust
use cryspglib::spg_get_magnetic_dataset;

// BCC AFM [111]: Fe at [0,0,0] and [0.5,0.5,0.5], opposite spins along [111]
let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
let positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];
let types = [26, 26];
let n = (3.0_f64).sqrt();
let moments = [[1.0/n, 1.0/n, 1.0/n], [-1.0/n, -1.0/n, -1.0/n]];

let result = spg_get_magnetic_dataset(&lattice, &positions, &types, Some(&moments), 1e-5);
// → SpglibMagneticSymmetry {
//     spacegroup_number: 229 (Im-3m),
//     uni_number: 1331,
//     bns_number: "166.101",
//     magnetic_type: BlackWhite (Type-3),
//     num_operations: 24,
//     ...
//   }
```

### Key API

| Function | Description |
|----------|-------------|
| `spg_get_dataset` | Non-magnetic space group identification |
| `spg_standardize_cell` | Standardize crystal cell |
| `spg_find_primitive` | Find primitive cell |
| `spg_refine_cell` | Wyckoff position refinement |
| `spg_get_symmetry` | Get symmetry operations |
| `spg_get_magnetic_dataset` | Full magnetic space group identification |
| `spg_get_magnetic_spacegroup_type` | Look up magnetic space group type by UNI number |
| `spg_get_magnetic_spacegroup_type_from_symmetry` | Identify magnetic space group type from symmetry operations |
| `spg_delaunay_reduce` | Delaunay lattice reduction |
| `spg_niggli_reduce` | Niggli lattice reduction |
| `spg_read_structure` | Parse POSCAR-format input |

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
