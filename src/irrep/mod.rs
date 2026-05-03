//! # Irreducible representations of the 230 crystallographic space groups
//!
//! This module provides irrep data at high-symmetry **k**-points for all 230
//! space groups, based on **Stokes & Hatch (1988)**, *Isotropy Subgroups of
//! the 230 Crystallographic Space Groups* (World Scientific).
//!
//! ## Three labeling conventions
//!
//! Each irrep is labeled in three widely-used conventions (cross-referenced
//! from Stokes & Hatch Table 7):
//!
//! | Convention | Reference | Example |
//! |-----------|-----------|---------|
//! | **Miller & Love** (ML) | Miller & Love (1967); Cracknell, Davies, Miller & Love (1979) | `GM1+`, `X3-` |
//! | **Kovalev** | Kovalev (1965, 1980, 1986) | `τ1`, `k6τ3` |
//! | **Bradley & Cracknell** (B&C) | Bradley & Cracknell (1972) | `Γ1+`, `X1` |
//!
//! ## Quick lookup
//!
//! Use the table below to find which crystal-system page contains your space
//! group.  After navigating to the page, use **Ctrl+F** to jump to the
//! space group number.
//!
//! | SG# | Symbol | Schoenflies | Crystal system | Page |
//! |-----|--------|------------|----------------|------|
//! | 1–2 | P1 … P-1 | C₁¹ … Cᵢ¹ | Triclinic | [`triclinic`] |
//! | 3–15 | P2 … C2/c | C₂¹ … C₂ₕ⁶ | Monoclinic | [`monoclinic`] |
//! | 16–74 | P222 … Imma | D₂¹ … D₂ₕ²⁸ | Orthorhombic | [`orthorhombic`] |
//! | 75–142 | P4 … I4₁/acd | C₄¹ … D₄ₕ²⁰ | Tetragonal | [`tetragonal`] |
//! | 143–167 | P3 … R-3c | C₃¹ … D₃ₔ⁶ | Trigonal | [`trigonal`] |
//! | 168–194 | P6 … P6₃/mmc | C₆¹ … D₆ₕ⁴ | Hexagonal | [`hexagonal`] |
//! | 195–230 | P23 … Ia-3d | T¹ … Oₕ¹⁰ | Cubic | [`cubic`] |
//!
//! ## Data format
//!
//! For each space group we list:
//!
//! 1. **k-points** — high-symmetry points in the Brillouin zone with their
//!    little co-groups.
//! 2. **Irreps** — each k-point's irreducible representations in all three
//!    label conventions, plus dimension, Stokes-Hatch image, and basis
//!    functions.
//! 3. **Isotropy subgroups** — for each irrep: the lower-symmetry space
//!    groups that can result from a distortion transforming as that irrep,
//!    together with the order-parameter direction, number of domains, cell
//!    basis, and origin shift.
//!
//! ## Programmatic access
//!
//! The [`types`] module provides [`IrrepData`], [`KPointData`], and
//! [`IsotropySubgroup`] structs that can be constructed from the
//! tabulated data for use in your own symmetry analysis code.
//!
//! [`types`]: crate::irrep::types
//! [`triclinic`]: crate::irrep::triclinic
//! [`monoclinic`]: crate::irrep::monoclinic
//! [`orthorhombic`]: crate::irrep::orthorhombic
//! [`tetragonal`]: crate::irrep::tetragonal
//! [`trigonal`]: crate::irrep::trigonal
//! [`hexagonal`]: crate::irrep::hexagonal
//! [`cubic`]: crate::irrep::cubic

pub mod types;
pub mod generated_data;
pub mod preamble;
pub mod query;

pub mod triclinic;
pub mod monoclinic;
pub mod orthorhombic;
pub mod tetragonal;
pub mod trigonal;
pub mod hexagonal;
pub mod cubic;
