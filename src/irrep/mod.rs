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
pub use types::generated_data;
pub mod bridge;
pub mod corep;
pub mod isotropy;
pub mod labels;
pub use labels::{IrrepLabels, LabelConvention};
pub mod magnetic_summary;
pub mod preamble;
pub mod query;
pub mod wigner;

// Exact rational affine/lattice layer, the concrete star decomposition and the
// parametric-k line engine (docs/full-irrep-subduction-plan.md, tasks 3-12).
// Promoted by milestone R8: the supported entry point is `subduction_api`, and
// the types its signatures use are public here.
pub mod subduction;

// Staged entry point: `subduce_irrep(&subgroup, probe_ml)` decomposes a parent
// irrep on the subgroup an isotropy record selects.  Gamma scalar irreps only,
// for now; the compound, k-folding and magnetic stages are tasks 6-11.
#[doc(hidden)]
pub use subduction::{
    IrrepSubduction, SubductionComponent, SubductionError, SubductionTarget, SubgroupEmbedding,
    subduce_irrep, subduce_irrep_with_embedding,
};

// Frozen little-group character tables of the parametric-k irreps behind the
// `other_wave_vector_subduction` rows (docs/subduction-audit.md).  Generated by
// `scripts/freeze_w_little_characters.py`; `tests/w_little_characters.rs` pins
// the anchors and the deliberately unresolved set.
#[doc(hidden)]
pub mod w_little_characters_data;

// The ISOTROPY program's **own** order-parameter component strings, per
// (parent space group, Miller-Love label, direction label).  Generated by
// `scripts/freeze_direction_descriptors.py` from the collected oracle output
// (`scripts/data/direction_descriptors_v1.json`); the public query is
// `isotropy::official_direction_descriptor` and the selector is
// `IsotropyDirection::OfficialDescriptor`.  Milestone R7.
pub mod direction_descriptors_data;

// R9 stage 1: the magnetic input contract.  One magnetic space group (UNI)
// turned into exact operations with the time-reversal flag preserved, and a
// proof that the set really is a group whose unitary part is a subgroup and
// whose antiunitary part is one coset of it.  See
// `docs/subduction-next-milestones.md` R9.
pub mod magnetic_embedding;

// R25 card 1: the input contract of a scalar Gamma magnetic subduction.  One
// verified embedding plus one explicit antiunitary sign are turned into the
// parent irrep's real +-1 corepresentation, with the group law checked on the
// operation set's classes and the parent's tabulated unitary character.  It is a
// separate entry point from `subduction_api` on purpose: no subduction happens
// here, and the ordinary path is untouched.  See
// `docs/magnetic-gamma-input.md`.
pub mod magnetic_input;

// R26 card 1: the first restricted scalar Gamma magnetic subduction.  It consumes
// the verified R25 `ScalarGammaModel`, takes its target identity independently
// from the magnetic corep summary at Gamma, reconstructs the one-dimensional
// Wigner type-A target's operators (including the antiunitary ones, which the
// summary leaves pending) from the target's own unitary characters and the actual
// multiplication of the magnetic group, and exhibits the complex intertwiner
// between the model's restriction and that target.  One entry point, private
// result fields, no matrix framework.  See `docs/magnetic-gamma-subduction.md`.
pub mod magnetic_subduction;

// R29 card 1: the two-dimensional real Gamma magnetic **matrix** subduction.
// It certifies the source's stored operation-aware matrices (orthogonality,
// character agreement with two storage surfaces, the complete group law and
// real irreducibility), binds every magnetic class to its unique stored
// operation, verifies the parent grey corepresentation law for either eta,
// decomposes the restriction into the summary's type-A real Gamma coreps with
// independently constructed antiunitary extensions, and exhibits normalized
// matrix intertwiners that cover the whole source space with a complex-vector
// reconstruction.  See `docs/magnetic-matrix-gamma-subduction.md`.
pub mod magnetic_matrix_subduction;

// R30: the supported high-level entry point for Gamma magnetic subduction.  It
// resolves a source irrep and a condensing magnetic record (by an existing
// direction or a UNI number, each with its own label convention), dispatches on
// the source dimension to the verified R25/R26 scalar or R29 matrix entry point
// unchanged, and reports the source, condensing, magnetic-record and target
// identities together.  The magnetic corep display label and the target's
// ordinary H-irrep CDML/BC labels are kept separate; no magnetic BC name is
// fabricated.  The ordinary `subduction_api` is not modified.
pub mod magnetic_subduction_api;

// Monodromy of those frozen line sources under a parent reciprocal shift: the
// label map `M_K` computed from character fingerprints (never from label
// names), the contract `(k + K, alpha) ~ (k, M_K(alpha))` it encodes, and the
// conjugation map it commutes with.  See `docs/subduction-conventions.md` §16.
// Promoted by milestone R8 (the parametric-k entry point reports its kind).
pub mod line_monodromy;

// The supported ordinary-subduction entry point (milestone R8): one request
// type, one computation entry, and a report carrying its own provenance.
pub mod subduction_api;
pub use subduction_api::{
    BlockView, CondensingIdentity, DirectionIdentity, SubductionApiError, SubductionReport,
    SubductionRequest, TargetView, subduce, subduction_for_direction,
};

pub mod cubic;
pub mod hexagonal;
pub mod monoclinic;
pub mod orthorhombic;
pub mod tetragonal;
pub mod triclinic;
pub mod trigonal;
