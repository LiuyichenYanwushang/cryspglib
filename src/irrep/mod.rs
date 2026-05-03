//! # Irreducible representations of the 230 crystallographic space groups
//!
//! This module provides irrep data at high-symmetry **k**-points for all 230
//! space groups, based on the **ISOTROPY Suite** (Stokes, Campbell & Hatch,
//! 2022 version).  Data is stored as static arrays and compiled directly into
//! the binary — no runtime I/O required.
//!
//! ## Coverage
//!
//! | Data | Entries | Coverage |
//! |------|---------|----------|
//! | Irreps | 4,777 | 100% |
//! | Character tables χ(g) | ~50,000 f64 values | 100% |
//! | Full matrices D(g) | ~580,000 f64 values | 100% |
//! | Non-magnetic isotropy subgroups | 15,239 | 100% |
//! | Magnetic isotropy subgroups | 16,721 | 100% |
//!
//! Data sources: `PIR_data.txt` (physical / real irreps), `CIR_data.txt`
//! (complex irreps, used as fallback), `data_irreps.txt` (metadata),
//! `data_isotropy.txt` (non-mag subgroups), `data_magnetic.txt` (mag subgroups).
//!
//! ## Three labeling conventions
//!
//! Each irrep is labeled in three widely-used conventions:
//!
//! | Convention | Reference | Example |
//! |-----------|-----------|---------|
//! | **Miller & Love** (ML) | Miller & Love (1967) | `GM1+`, `X3-` |
//! | **Kovalev** | Kovalev (1986) | `τ1`, `k6τ3` |
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
//! ## Programmatic access
//!
//! The [`query`] module provides the main entry point:
//!
//! ```rust,no_run
//! use cryspglib::irrep::query::*;
//!
//! // SG 221 (Pm-3m) at Γ point
//! let irreps = irreps_of(221);
//! let gm4m = irreps.iter().find(|r| r.ml == "GM4-").unwrap();
//!
//! // Character table
//! for (i, &chi) in gm4m.characters().iter().enumerate() {
//!     println!("op {}: χ = {}", i, chi);
//! }
//!
//! // Non-magnetic isotropy subgroups
//! for sub in gm4m.subgroups() {
//!     println!("#{} {} — dir={}", sub.sg, sub.symbol, sub.direction);
//! }
//!
//! // Magnetic isotropy subgroups
//! for sub in gm4m.magnetic_subgroups() {
//!     println!("{} {} — dir={}", sub.bns_label, sub.iso_label, sub.direction);
//! }
//! ```
//!
//! ## Key types
//!
//! | Type | Description |
//! |------|-------------|
//! | [`IrrepRecord`] | Single irrep with labels, k-vector, character table, matrices, subgroups |
//! | [`IsotropyRecord`] | Non-magnetic isotropy subgroup (SG#, symbol, direction, domains) |
//! | [`MagneticIsotropyRecord`] | Magnetic isotropy subgroup (BNS label, direction) |
//! | [`KPointSummary`] | k-point with label, coordinates, and irrep indices |
//!
//! [`IrrepRecord`]: crate::irrep::types::IrrepRecord
//! [`IsotropyRecord`]: crate::irrep::types::IsotropyRecord
//! [`MagneticIsotropyRecord`]: crate::irrep::types::MagneticIsotropyRecord
//! [`KPointSummary`]: crate::irrep::query::KPointSummary
//! [`query`]: crate::irrep::query
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
pub mod corep;

pub mod triclinic;
pub mod monoclinic;
pub mod orthorhombic;
pub mod tetragonal;
pub mod trigonal;
pub mod hexagonal;
pub mod cubic;
