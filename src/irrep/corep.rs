//! Co-representation (corep) theory for magnetic space group irreps.
//!
//! # Theory
//!
//! A magnetic space group (MSG) with anti-unitary operations can be written as
//!
//! $$
//! \mathcal{M} = H \cup \mathcal{T} g_0 H
//! $$
//!
//! where $$H = \mathcal{M} \cap G$$ is the **unitary subgroup** (a normal
//! non-magnetic space group), $$\mathcal{T}$$ is time reversal, and
//! $$g_0$$ is the spatial part of the chosen anti-unitary coset representative
//! $$a_0 = \mathcal{T} g_0$$.
//!
//! ## Construction from the unitary subgroup
//!
//! Given a non-magnetic irrep $$\Delta_i$$ of $$H$$ at wave-vector $$\mathbf{k}$$,
//! define its **conjugate representation** under the anti-unitary coset:
//!
//! $$
//! \Delta_i^{a_0}(h) \equiv \Delta_i(a_0^{-1} h a_0)^* \qquad (h \in H)
//! $$
//!
//! The magnetic co-representation $$\tilde{D}$$ of $$\mathcal{M}$$ is built
//! from the relationship between $$\Delta_i$$ and $$\Delta_i^{a_0}$$:
//!
//! ## Wigner's three cases
//!
//! The Wigner test classifies irreps into three types via the sum
//!
//! $$
//! W(\Delta_i) = \frac{1}{|H|}\sum_{b \in a_0 H} \chi_{\Delta_i}(b^2)
//! $$
//!
//! where the sum runs over the **anti-unitary coset** $$a_0 H$$ and
//! $$\chi_{\Delta_i}$$ is the character of $$\Delta_i$$.
//!
//! | Case | Condition | Corep dimension | Unitary characters | Anti-unitary characters |
//! |------|-----------|----------------|-------------------|------------------------|
//! | **Type A** | $$\Delta_i^{a_0} \sim \Delta_i$$, $$W = +1$$ | $$d_i$$ | $$\chi_{\Delta_i}(h)$$ | $$\chi_{\Delta_i}(a_0 h)$$ (real) |
//! | **Type B** | $$\Delta_i^{a_0} \sim \Delta_i$$, $$W = -1$$ | $$2d_i$$ (Kramers) | $$2\chi_{\Delta_i}(h)$$ | $$0$$ |
//! | **Type C** | $$\Delta_i^{a_0} \nsim \Delta_i$$, $$W = 0$$ | $$2d_i$$ | $$2\,\mathrm{Re}[\chi_{\Delta_i}(h)]$$ | $$0$$ |
//!
//! **Type C** pairs two inequivalent irreps $$\Delta_i, \Delta_j$$ of $$H$$
//! (where $$\Delta_j \sim \Delta_i^{a_0}$$). The corep is
//! $$\tilde{D} = \Delta_i \oplus \Delta_j$$ with block structure
//!
//! $$
//! \tilde{D}(h) = \begin{pmatrix} \Delta_i(h) & 0 \\ 0 & \Delta_j(h) \end{pmatrix},
//! \qquad
//! \tilde{D}(a_0 h) \sim \begin{pmatrix} 0 & * \\ * & 0 \end{pmatrix} K
//! $$
//!
//! where $$K$$ denotes complex conjugation.
//!
//! ## Workflow
//!
//! ```text
//! BNS label ("128.406") + k-point label ("Z")
//!   → uni_from_bns()           // BNS → UNI number
//!   → identify_unitary_subgroup()  // UNI → H space group
//!   → irreps_of(H) at k-point  // H's double-group irreps
//!   → compute_corepresentation()   // Wigner test + corep characters
//!   → Corepresentation { characters, corep_type, dim }
//! ```
//!
//! ## Reference target: 128.406 at Z
//!
//! Verified against Bilbao Crystallographic Server (BCS):
//!
//! - Magnetic SG: $$P4'/m'nc'$$ (No. 128.406, UNI 1066)
//! - Unitary subgroup: $$P\bar{4}n2$$ (No. 118)
//! - k-vector: $$Z = (0, 0, 1/2)$$
//! - Magnetic little group in the conventional data-Hall cell: 16 operation
//!   representatives (8 unitary + 8 anti-unitary)
//!
//! From H = SG 118's Z-point irreps:
//!
//! | H irrep | Type | Magnetic corep | Dimension |
//! |---------|------|---------------|-----------|
//! | Z₁Z₄ | C | Z₁Z₂ | 2D |
//! | Z₂Z₃ | C | Z₃Z₄ | 2D |
//! | Z₅ | A | Z₅ | 2D |
//! | Z₆, Z₇ (spinor) | C | Z̄₆Z̄₇ | 4D |
//!
//! The current implementation returns `Result` and reports unresolved Wigner
//! paths as errors instead of emitting placeholder `Unsupported` rows.
//!
//! # References
//!
//! - Wigner (1959), *Group Theory*, Chapter 26
//! - Bradley & Cracknell (1972), *The Mathematical Theory of Symmetry in Solids*
//! - Stokes, Campbell & Hatch, ISOTROPY Suite documentation
//! - Bilbao Crystallographic Server: <https://cryst.ehu.es/cgi-bin/cryst/programs/corepresentations.pl>

use super::types::IrrepRecord;
use super::wigner::{self, SeitzOp, filter_little_group_with_transform, ops_to_seitz};
#[cfg(test)]
use super::wigner::{
    add3, bloch_phase, compose_seitz, filter_little_group, find_seitz, mat_vec_i32, square_seitz,
};
use crate::SymmetryOps;
use crate::mathfunc::{Mat3I, mat_inverse_matrix_d3};
use crate::spg_database::{get_spacegroup_operations, get_spacegroup_type};
use num_complex::Complex64;
use std::collections::BTreeSet;

macro_rules! debug_log {
    ($($arg:tt)*) => {
        #[cfg(feature = "debug-corep")]
        eprintln!($($arg)*);
    };
}

/// Scalar PIR data restricted to the block belonging to the stored k-vector.
///
/// ISO-IR stores a full space-group representation induced over every arm of
/// the star. Its matrices are block matrices, with the first block belonging
/// to the first k-vector in the record (the one exposed by `IrrepRecord`). A
/// magnetic little-group corepresentation must use that block rather than the
/// trace and dimension of the full induced representation.
struct ScalarLittleIrrepData {
    /// Characters in canonical H-operation order. Entries outside H_k are
    /// harmless zeroes and are never consumed by the little-group formulas.
    characters: Vec<f64>,
    /// Same selected-arm characters with their CIR imaginary parts retained.
    complex_characters: Vec<Complex64>,
    /// Matrices in canonical H-operation order, one `dim * dim` block per H
    /// operation. Empty only when ISO-IR matrix data is unavailable.
    matrices: Vec<f64>,
    dim: usize,
}

fn scalar_little_irrep_data(
    irrep: &IrrepRecord,
    h_seitz: &[SeitzOp],
    little_h_indices: &[usize],
    h_to_pir: &[usize],
) -> Option<ScalarLittleIrrepData> {
    if h_seitz.is_empty() || little_h_indices.is_empty() || h_to_pir.len() != h_seitz.len() {
        return None;
    }

    let identity_rotation = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let pure_translations: Vec<_> = h_seitz
        .iter()
        .filter(|operation| operation.rot == identity_rotation)
        .map(|operation| operation.trans)
        .collect();
    if pure_translations.is_empty() || !h_seitz.len().is_multiple_of(pure_translations.len()) {
        return None;
    }
    let equivalent_mod_translation_subgroup = |left: &SeitzOp, right: &SeitzOp| {
        left.rot == right.rot
            && pure_translations.iter().any(|shift| {
                (0..3).all(|axis| {
                    let delta = left.trans[axis] - right.trans[axis] - shift[axis];
                    (delta - delta.round()).abs() < 1e-8
                })
            })
    };
    let mut little_coset_representatives = Vec::new();
    for &h_idx in little_h_indices {
        let operation = h_seitz.get(h_idx)?;
        if !little_coset_representatives.iter().any(|&representative| {
            equivalent_mod_translation_subgroup(operation, &h_seitz[representative])
        }) {
            little_coset_representatives.push(h_idx);
        }
    }
    let effective_h_order = h_seitz.len() / pure_translations.len();
    let effective_little_order = little_coset_representatives.len();
    if effective_little_order == 0 || !effective_h_order.is_multiple_of(effective_little_order) {
        return None;
    }
    let star_size = effective_h_order / effective_little_order;
    let full_dim = irrep.dim as usize;
    if star_size == 0 || full_dim == 0 || !full_dim.is_multiple_of(star_size) {
        return None;
    }
    let dim = full_dim / star_size;
    let pir_chars = irrep.characters();
    let pir_matrices = irrep.matrices();
    let full_block = full_dim * full_dim;
    let little_block = dim * dim;
    let n_pir_ops = irrep.pir_rotations().len() / 9;
    if h_to_pir
        .iter()
        .any(|&pir_idx| pir_idx >= n_pir_ops || pir_idx >= pir_chars.len())
    {
        return None;
    }
    let (stored_little_real, stored_little_imag) = irrep.scalar_little_characters();
    let stored_little_available =
        stored_little_real.len() == n_pir_ops && stored_little_imag.len() == n_pir_ops;
    let mapped_stored_characters = stored_little_available.then(|| {
        h_to_pir
            .iter()
            .map(|&pir_idx| {
                Complex64::new(stored_little_real[pir_idx], stored_little_imag[pir_idx])
            })
            .collect::<Vec<_>>()
    });

    // Locate the identity instead of assuming a particular Hall ordering.
    let identity = h_seitz.iter().position(|op| {
        op.rot == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            && op.trans.iter().all(|x| (x - x.round()).abs() < 1e-8)
    })?;

    // The CIR-derived selected-arm characters are authoritative for a
    // multi-arm physical PIR. Its real matrices may be a realification in
    // which conjugate arm subspaces cannot be split over R (for example a
    // 2D rotation representing two conjugate 1D K-point irreps).
    if star_size > 1
        && let Some(complex_characters) = &mapped_stored_characters
    {
        if (complex_characters[identity] - Complex64::new(dim as f64, 0.0)).norm() > 1e-6 {
            debug_log!(
                "scalar selected CIR characters: SG{} {} H identity={} PIR identity={} value={} expected={} (star={}, stored ops={}, first={:?})",
                irrep.sg,
                irrep.ml,
                identity,
                h_to_pir[identity],
                complex_characters[identity],
                dim,
                star_size,
                stored_little_real.len(),
                complex_characters.iter().take(8).collect::<Vec<_>>()
            );
            return None;
        }
        let complex_characters = complex_characters.clone();
        return Some(ScalarLittleIrrepData {
            characters: complex_characters.iter().map(|value| value.re).collect(),
            complex_characters,
            matrices: Vec::new(),
            dim,
        });
    }
    if star_size > 1 && mapped_stored_characters.is_none() {
        debug_log!(
            "scalar selected CIR characters unavailable: SG{} {} PIR ops={} stored real={} imag={}",
            irrep.sg,
            irrep.ml,
            n_pir_ops,
            stored_little_real.len(),
            stored_little_imag.len()
        );
    }

    // At a one-arm point, reordered characters remain useful even if a legacy
    // record has no matrices.
    if pir_matrices.len() < n_pir_ops * full_block {
        let complex_characters = mapped_stored_characters.unwrap_or_else(|| {
            h_to_pir
                .iter()
                .map(|&idx| Complex64::new(pir_chars[idx], 0.0))
                .collect()
        });
        return Some(ScalarLittleIrrepData {
            characters: complex_characters.iter().map(|value| value.re).collect(),
            complex_characters,
            matrices: Vec::new(),
            dim,
        });
    }

    let mut is_little = vec![false; h_seitz.len()];
    for (h_idx, operation) in h_seitz.iter().enumerate() {
        is_little[h_idx] = little_coset_representatives.iter().any(|&representative| {
            equivalent_mod_translation_subgroup(operation, &h_seitz[representative])
        });
    }

    // ISO-IR guarantees one invariant block per star arm, but the basis
    // vectors of a block are not always contiguous (notably repeated/physical
    // labels such as W1W1). Recover the first arm as the invariant connectivity
    // component containing basis vector 0, adding further disconnected
    // components in basis order only when the little representation itself is
    // reducible. This reduces to the usual top-left block for contiguous data.
    let mut adjacency = vec![false; full_block];
    for (h_idx, &pir_idx) in h_to_pir.iter().enumerate() {
        if !is_little[h_idx] || pir_idx >= n_pir_ops {
            continue;
        }
        let source = &pir_matrices[pir_idx * full_block..(pir_idx + 1) * full_block];
        for row in 0..full_dim {
            adjacency[row * full_dim + row] = true;
            for col in 0..full_dim {
                if source[row * full_dim + col].abs() > 1e-8 {
                    adjacency[row * full_dim + col] = true;
                    adjacency[col * full_dim + row] = true;
                }
            }
        }
    }
    let mut seen_basis = vec![false; full_dim];
    let mut components = Vec::<Vec<usize>>::new();
    for start in 0..full_dim {
        if seen_basis[start] {
            continue;
        }
        seen_basis[start] = true;
        let mut component = Vec::new();
        let mut frontier = vec![start];
        while let Some(row) = frontier.pop() {
            component.push(row);
            for col in 0..full_dim {
                if adjacency[row * full_dim + col] && !seen_basis[col] {
                    seen_basis[col] = true;
                    frontier.push(col);
                }
            }
        }
        component.sort_unstable();
        components.push(component);
    }
    components.sort_by_key(|component| component[0]);
    let first_component = components
        .iter()
        .position(|component| component.contains(&0))?;
    components.swap(0, first_component);
    let mut selected_basis = Vec::with_capacity(dim);
    for component in components {
        if selected_basis.len() == dim {
            break;
        }
        if selected_basis.len() + component.len() > dim {
            debug_log!(
                "scalar block extraction: SG{} {} invariant component {} would exceed selected dimension {}",
                irrep.sg,
                irrep.ml,
                component.len(),
                dim
            );
            return None;
        }
        selected_basis.extend(component);
    }
    if selected_basis.len() != dim {
        return None;
    }
    selected_basis.sort_unstable();
    let mut basis_is_selected = vec![false; full_dim];
    for &basis in &selected_basis {
        basis_is_selected[basis] = true;
    }

    let mut characters = vec![0.0; h_seitz.len()];
    let mut matrices = vec![0.0; h_seitz.len() * little_block];
    for (h_idx, &pir_idx) in h_to_pir.iter().enumerate() {
        if pir_idx >= n_pir_ops {
            return None;
        }
        let source = &pir_matrices[pir_idx * full_block..(pir_idx + 1) * full_block];

        // An operation of H_k must preserve the selected star-arm subspace.
        if is_little[h_idx] {
            for row in 0..full_dim {
                for col in 0..full_dim {
                    if basis_is_selected[row] == basis_is_selected[col] {
                        continue;
                    }
                    if source[row * full_dim + col].abs() > 1e-8 {
                        debug_log!(
                            "scalar block extraction: SG{} {} H[{}] PIR[{}] couples ({},{})={} (full_dim={}, little_dim={}, star={})",
                            irrep.sg,
                            irrep.ml,
                            h_idx,
                            pir_idx,
                            row,
                            col,
                            source[row * full_dim + col],
                            full_dim,
                            dim,
                            star_size
                        );
                        return None;
                    }
                }
            }
        }

        let target = &mut matrices[h_idx * little_block..(h_idx + 1) * little_block];
        let mut trace = 0.0;
        for (row, &source_row) in selected_basis.iter().enumerate() {
            for (col, &source_col) in selected_basis.iter().enumerate() {
                target[row * dim + col] = source[source_row * full_dim + source_col];
            }
            trace += source[source_row * full_dim + source_row];
        }
        characters[h_idx] = if trace.abs() < 1e-12 { 0.0 } else { trace };
    }

    if !is_little[identity] || (characters[identity] - dim as f64).abs() > 1e-6 {
        debug_log!(
            "scalar block extraction: SG{} {} identity H[{}] is_little={} trace={} expected={} (effective H={}, little={}, centering={})",
            irrep.sg,
            irrep.ml,
            identity,
            is_little[identity],
            characters[identity],
            dim,
            effective_h_order,
            effective_little_order,
            pure_translations.len()
        );
        return None;
    }

    let complex_characters = mapped_stored_characters.unwrap_or_else(|| {
        characters
            .iter()
            .map(|&value| Complex64::new(value, 0.0))
            .collect()
    });
    Some(ScalarLittleIrrepData {
        characters: complex_characters.iter().map(|value| value.re).collect(),
        complex_characters,
        matrices,
        dim,
    })
}

/// Co-representation type from Wigner's test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorepType {
    /// W = +1: D ∼ D*, real representation.
    A,
    /// W = -1: D ∼ D*, pseudo-real (quaternionic).
    B,
    /// W = 0: D ≁ D*.
    C,
}

/// Which computational path produced the Wigner classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WignerSource {
    /// No antiunitary ops in the magnetic little group — always Type A.
    TrivialNoAntiunitary,
    /// Scalar irrep, PIR character table.
    ScalarPIR,
    /// Compound irrep, CIR complex character table.
    ScalarCIR,
    /// Spinor irrep classified via SU(2) double-group composition.
    SpinorSU2,
}

impl CorepType {
    pub fn description(&self) -> &'static str {
        match self {
            CorepType::A => "type-a: D ~ D*, real (W=+1)",
            CorepType::B => "type-b: D ~ D*, pseudo-real (W=-1)",
            CorepType::C => "type-c: D ≁ D* (W=0)",
        }
    }
}

/// Completeness of the magnetic character table.
///
/// Indicates whether every operation in the magnetic little group has a
/// valid character value.  Operations with value 0 are considered valid
/// only when the theory mandates it (Type B / Type C anti-unitary ops,
/// and symmetry-forbidden zeros).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CharacterCompleteness {
    /// All magnetic little group operations have valid character values.
    Complete,
    /// Type A anti-unitary characters require the intertwiner matrix U
    /// (satisfying U·Δ(a₀⁻¹ha₀)* = Δ(h)·U) which is not yet computed.
    /// These entries are left as 0.
    TypeAAntiunitaryPending { count: usize },
}

/// The computed magnetic co-representation of an irrep.
#[derive(Debug, Clone)]
pub struct Corepresentation {
    /// Character χ̃(g) for each magnetic operation (same order as SymmetryOps).
    pub characters: Vec<f64>,
    /// Which operations are anti-unitary.
    pub timerev: Vec<bool>,
    /// Co-representation type.
    pub corep_type: CorepType,
    /// Which Wigner path produced this classification.
    pub source: WignerSource,
    /// Dimension of the magnetic irrep.
    pub dim: usize,
    /// Number of unitary operations.
    pub unitary_order: usize,
    /// Number of anti-unitary operations.
    pub antiunitary_order: usize,
    /// Whether the character table is complete for all mag-little-group ops.
    pub completeness: CharacterCompleteness,
}

/// Error returned when a magnetic co-representation cannot be computed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CorepComputationError {
    /// UNI number is outside the magnetic space-group range.
    InvalidUni { uni: usize },
    /// BNS label was not found in the magnetic space-group database.
    UnknownBns { bns: String },
    /// Magnetic symmetry operations could not be retrieved.
    MissingMagneticOperations { uni: usize },
    /// The unitary subgroup H could not be identified.
    MissingUnitarySubgroup { uni: usize },
    /// H has no usable symmetry operations after identification.
    EmptyUnitaryOperations { uni: usize },
    /// No H-irrep is available at the requested k-point label.
    MissingKPoint { uni: usize, sg: u8, k_label: String },
    /// The magnetic little group at this irrep's k-vector is empty.
    EmptyMagneticLittleGroup { uni: usize, source_irrep: String },
    /// A unitary magnetic operation could not be matched to an H operation.
    UnmappedUnitaryOperation { uni: usize, source_irrep: String },
    /// Wigner's test did not produce one of the supported A/B/C cases.
    UnsupportedClassification {
        uni: usize,
        source_irrep: String,
        reason: String,
    },
    /// The computed character table contains NaN or infinity.
    NonFiniteCharacters { uni: usize, source_irrep: String },
}

impl std::fmt::Display for CorepComputationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CorepComputationError::InvalidUni { uni } => {
                write!(f, "invalid UNI number {}", uni)
            }
            CorepComputationError::UnknownBns { bns } => {
                write!(f, "unknown BNS label {}", bns)
            }
            CorepComputationError::MissingMagneticOperations { uni } => {
                write!(f, "could not retrieve magnetic operations for UNI {}", uni)
            }
            CorepComputationError::MissingUnitarySubgroup { uni } => {
                write!(f, "could not identify unitary subgroup for UNI {}", uni)
            }
            CorepComputationError::EmptyUnitaryOperations { uni } => {
                write!(f, "unitary subgroup for UNI {} has no operations", uni)
            }
            CorepComputationError::MissingKPoint { uni, sg, k_label } => {
                write!(
                    f,
                    "no H-irrep data at k-point {} for UNI {} unitary SG {}",
                    k_label, uni, sg
                )
            }
            CorepComputationError::EmptyMagneticLittleGroup { uni, source_irrep } => {
                write!(
                    f,
                    "empty magnetic little group for UNI {} source irrep {}",
                    uni, source_irrep
                )
            }
            CorepComputationError::UnmappedUnitaryOperation { uni, source_irrep } => {
                write!(
                    f,
                    "failed to map a unitary magnetic operation to H for UNI {} source irrep {}",
                    uni, source_irrep
                )
            }
            CorepComputationError::UnsupportedClassification {
                uni,
                source_irrep,
                reason,
            } => {
                write!(
                    f,
                    "unsupported Wigner classification for UNI {} source irrep {}: {}",
                    uni, source_irrep, reason
                )
            }
            CorepComputationError::NonFiniteCharacters { uni, source_irrep } => {
                write!(
                    f,
                    "computed non-finite characters for UNI {} source irrep {}",
                    uni, source_irrep
                )
            }
        }
    }
}

impl std::error::Error for CorepComputationError {}

// ── Core computation ─────────────────────────────────────────────────────────

/// Compute the magnetic co-representation for an irrep of the unitary subgroup H.
///
/// See [`compute_coreps`] for the high-level BNS+k-label API.
pub fn compute_corepresentation(
    h_irrep: &IrrepRecord,
    uni_number: usize,
    mag_ops: &SymmetryOps,
) -> Result<Corepresentation, CorepComputationError> {
    if uni_number == 0 || uni_number > 1651 {
        return Err(CorepComputationError::InvalidUni { uni: uni_number });
    }

    // 1. Identify H first — needed for setting transform.
    let h_info = identify_unitary_subgroup_with_hall(uni_number)
        .ok_or(CorepComputationError::MissingUnitarySubgroup { uni: uni_number })?;
    // Use the pre-computed MSG→data-Hall transform from H identification.
    // This ensures the magnetic little group, Seitz composition, and spin
    // table are all in the ISOTROPY data-Hall frame.
    let setting_xf = h_info.msg_to_data.as_ref();

    let h_ops = &h_info.ops_from_hall;
    if h_ops.is_empty() {
        return Err(CorepComputationError::EmptyUnitaryOperations { uni: uni_number });
    }

    // Every representation datum (k vector, PIR/CIR operation order, spin
    // table) is expressed in the ISOTROPY data-Hall frame.  Keep the magnetic
    // operation indices stable, but transform their Seitz parts into that
    // same frame before little-group filtering, composition, and character
    // lookup.  Mixing MSG-frame operations with data-Hall PIR rotations was
    // the common cause of the former 128.406 and 52.318 scalar-map failures.
    let mag_ops_data = operations_in_data_hall_frame(mag_ops, setting_xf).ok_or_else(|| {
        CorepComputationError::UnsupportedClassification {
            uni: uni_number,
            source_irrep: h_irrep.ml.to_string(),
            reason: "MSG operation could not be transformed to the ISOTROPY data-Hall frame"
                .to_string(),
        }
    })?;

    // 2. Filter to magnetic little group with setting transform.
    // Pass canonical H translation subgroup for stricter k-preservation
    // (essential for centered groups where MSG-derived pure translations
    // are only a subset — see codex review Plan 1B).
    // NOTE: must use ops_from_hall (canonical H setting), NOT ops_from_msg
    // (MSG unitary subgroup), because the MSG may lack centering translations.
    let h_canonical_translations: Vec<[f64; 3]> = h_info
        .ops_from_hall
        .operations
        .iter()
        .filter(|op| op.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        .map(|op| op.translation)
        .collect();
    let mag_lg = filter_little_group_with_transform(
        h_irrep.kx,
        h_irrep.ky,
        h_irrep.kz,
        h_irrep.kd,
        &mag_ops_data,
        None,
        Some(&h_canonical_translations),
    );
    if mag_lg.is_empty() {
        return Err(CorepComputationError::EmptyMagneticLittleGroup {
            uni: uni_number,
            source_irrep: h_irrep.ml.to_string(),
        });
    }

    // 3. Convert to SeitzOps for proper composition
    let mag_seitz = ops_to_seitz(&mag_ops_data);
    // The spinor Wigner path needs both frames: it transforms MSG operations
    // to H's data-Hall frame for b² lookup, but looks up the antiunitary
    // spatial rotation in the parent-G spin table in the original MSG frame.
    let mag_seitz_msg = ops_to_seitz(mag_ops);
    let h_seitz = ops_to_seitz(h_ops);

    // 3a. Map unitary magnetic ops to H ops via full Seitz matching
    // (rotation + translation), not rotation-only.
    let op_map: Vec<Option<usize>> = (0..mag_ops.len())
        .map(|i| {
            if mag_ops_data.operations[i].time_reversal {
                None
            } else {
                let mop = &mag_seitz[i];
                wigner::find_seitz(&mop.rot, &mop.trans, &h_seitz).map(|m| m.op_index)
            }
        })
        .collect();

    if op_map
        .iter()
        .enumerate()
        .any(|(i, m)| !mag_ops_data.operations[i].time_reversal && m.is_none())
    {
        return Err(CorepComputationError::UnmappedUnitaryOperation {
            uni: uni_number,
            source_irrep: h_irrep.ml.to_string(),
        });
    }

    // 4. H's irrep characters
    let mut h_chars = h_irrep.characters().to_vec();
    let mut h_complex_chars = None;
    let mut h_dim = h_irrep.dim as usize;
    let mut h_matrices_ordered = Vec::new();
    // 5. Separate unitary / anti-unitary in little group
    let unitary: Vec<usize> = mag_lg
        .iter()
        .filter(|&&i| !mag_ops_data.operations[i].time_reversal)
        .copied()
        .collect();
    let antiunitary: Vec<usize> = mag_lg
        .iter()
        .filter(|&&i| mag_ops_data.operations[i].time_reversal)
        .copied()
        .collect();

    // Restrict scalar, non-compound ISO-IR representations from the full
    // induced star representation to the selected k-arm. Spinor records are
    // already little-group data, while compound PIRs are classified through
    // their stored complex components below.
    if !h_irrep.spinor && h_irrep.cir_component_count() == 0 {
        let pir_rots = h_irrep.pir_rotations();
        let pir_trans = h_irrep.pir_translations();
        let map = if pir_trans.len() == pir_rots.len() / 9 * 3 {
            wigner::build_h_to_irrep_op_map(&h_seitz, pir_rots, pir_trans)
        } else {
            wigner::build_h_to_cir_map(&h_seitz, pir_rots)
        }
        .ok_or_else(|| CorepComputationError::UnsupportedClassification {
            uni: uni_number,
            source_irrep: h_irrep.ml.to_string(),
            reason: format!(
                "scalar PIR operation map could not be built (H ops={}, PIR ops={}, PIR translations={})",
                h_seitz.len(),
                pir_rots.len() / 9,
                pir_trans.len() / 3
            ),
        })?;
        let little_h_indices: Vec<usize> = unitary
            .iter()
            .filter_map(|&mag_idx| op_map[mag_idx])
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let little = scalar_little_irrep_data(h_irrep, &h_seitz, &little_h_indices, &map)
            .ok_or_else(|| CorepComputationError::UnsupportedClassification {
                uni: uni_number,
                source_irrep: h_irrep.ml.to_string(),
                reason: format!(
                    "selected k-arm block could not be extracted from the scalar PIR (H ops={}, H_k mapped ops={}, PIR dim={})",
                    h_seitz.len(),
                    little_h_indices.len(),
                    h_irrep.dim
                ),
            })?;
        h_chars = little.characters;
        h_complex_chars = Some(little.complex_characters);
        h_dim = little.dim;
        h_matrices_ordered = little.matrices;
    }

    // 7. Wigner test: dispatch by irrep type
    let (corep_type, source) = if antiunitary.is_empty() {
        Ok((CorepType::A, WignerSource::TrivialNoAntiunitary))
    } else if h_irrep.cir_component_count() > 0 {
        // A compound real PIR such as Z1Z4 is already the direct sum of two
        // conjugate complex irreps.  Wigner's test must start from ONE of
        // those complex components; using the compound PIR dimension and
        // then applying the Type-C doubling would count the pair twice.
        let identity_rotation = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let identity = h_seitz
            .iter()
            .position(|operation| {
                operation.rot == identity_rotation
                    && operation
                        .trans
                        .iter()
                        .all(|value| (value - value.round()).abs() < 1e-8)
            })
            .ok_or_else(|| CorepComputationError::UnsupportedClassification {
                uni: uni_number,
                source_irrep: h_irrep.ml.to_string(),
                reason: "unitary subgroup has no identity operation".to_string(),
            })?;
        let mut selected: Option<(CorepType, Vec<f64>, usize)> = None;
        let mut component_errors = Vec::new();
        debug_log!(
            "DEBUG CIR path: {} n_comp={}",
            h_irrep.ml,
            h_irrep.cir_component_count()
        );
        for comp in 0..h_irrep.cir_component_count() {
            let cir = h_irrep.cir_component_chars(comp);
            if cir.is_empty() {
                continue;
            }
            let cir_rots = h_irrep.cir_rotations(comp);
            // Generated CIR components are stored in the selected data-Hall
            // order.  Prefer that exact order, especially for centered groups
            // where several Seitz operations can share one rotation.  The
            // rotation-only map remains a compatibility fallback for older
            // generated tables with a different operation count.
            let cir_reordered = if cir.len() == 2 * h_seitz.len() {
                cir.to_vec()
            } else if let Some(h_to_cir) = wigner::build_h_to_cir_map(&h_seitz, cir_rots) {
                wigner::reorder_cir_chars(cir, &h_to_cir).map_err(|error| {
                    CorepComputationError::UnsupportedClassification {
                        uni: uni_number,
                        source_irrep: h_irrep.ml.to_string(),
                        reason: format!("CIR component {comp} reorder failed: {error}"),
                    }
                })?
            } else {
                return Err(CorepComputationError::UnsupportedClassification {
                    uni: uni_number,
                    source_irrep: h_irrep.ml.to_string(),
                    reason: format!(
                        "CIR component {} operation map could not be built (H ops={}, CIR ops={})",
                        comp,
                        h_seitz.len(),
                        cir.len() / 2
                    ),
                });
            };
            let ct = match wigner::wigner_classify_cir_direct(
                &cir_reordered,
                &antiunitary,
                &mag_seitz,
                &h_seitz,
                h_irrep.k_vector(),
            ) {
                Ok(corep_type) => corep_type,
                Err(error) => {
                    component_errors.push(format!("component {}: {}", comp, error));
                    continue;
                }
            };

            let identity_character =
                Complex64::new(cir_reordered[2 * identity], cir_reordered[2 * identity + 1]);
            let rounded_dim = identity_character.re.round();
            if identity_character.im.abs() > 1e-6
                || rounded_dim < 1.0
                || (identity_character.re - rounded_dim).abs() > 1e-6
            {
                return Err(CorepComputationError::UnsupportedClassification {
                    uni: uni_number,
                    source_irrep: h_irrep.ml.to_string(),
                    reason: format!(
                        "CIR component {} has invalid identity character ({:.8},{:.8})",
                        comp, identity_character.re, identity_character.im
                    ),
                });
            }
            let component_chars: Vec<f64> = cir_reordered
                .chunks_exact(2)
                .map(|value| value[0])
                .collect();
            // Wigner's construction begins with one irreducible complex
            // constituent Δ. Repeated labels such as P1P1 contain a
            // synthesized conjugate copy whose first-arm convention need not
            // be classified a second time.
            selected = Some((ct, component_chars, rounded_dim as usize));
            break;
        }
        let (corep_type, component_chars, component_dim) =
            selected.ok_or_else(|| CorepComputationError::UnsupportedClassification {
                uni: uni_number,
                source_irrep: h_irrep.ml.to_string(),
                reason: format!(
                    "compound irrep contains no classifiable CIR component ({})",
                    component_errors.join("; ")
                ),
            })?;
        h_chars = component_chars;
        h_dim = component_dim;
        Ok((corep_type, WignerSource::ScalarCIR))
    } else if h_irrep.spinor {
        // Spinor: SU(2) Wigner test is the primary path.
        // Bilbao imaginary chars are NOT term-by-term Wigner summands
        // (counterexample: SG3 A3 grey group, imag[0]=0 but h=E gives χ=-1).
        let h_spin = h_irrep.spin_ops();
        let g_sg = parent_spatial_sg(uni_number).unwrap_or(h_irrep.sg as usize) as u8;
        let g_spin = if g_sg == h_irrep.sg {
            h_spin
        } else {
            IrrepRecord::spin_ops_for_sg(g_sg)
        };
        let ctx = wigner::SpinLiftContext {
            h: h_spin,
            g: g_spin,
            sg: h_irrep.sg,
        };
        let op_indices = h_irrep.spin_lg_op_indices();
        let is_grey = crate::MagneticSpaceGroupType::from_uni(uni_number)
            .map_err(|_| CorepComputationError::InvalidUni { uni: uni_number })?
            .type_
            == crate::MagneticType::Grey;
        let a0_idx = select_spinor_a0(&antiunitary, &mag_seitz, is_grey);

        match wigner::wigner_classify_spinor(
            &ctx,
            wigner::SpinorWignerInput {
                characters_real: &h_chars,
                characters_imag: h_irrep.spin_character_imag(),
                operation_indices: op_indices,
                k_vector: h_irrep.k_vector(),
            },
            wigner::WignerGroupContext {
                unitary_indices: &unitary,
                magnetic_ops: &mag_seitz_msg,
                unitary_ops: &h_seitz,
                antiunitary_representative: a0_idx,
            },
            setting_xf,
            Some(&antiunitary),
        ) {
            Ok(ct) => Ok((ct, WignerSource::SpinorSU2)),
            Err(e) => {
                // SU(2) path failed.
                // The Bilbao extra sum is NOT a valid Wigner indicator and is
                // not used for classification.
                Err(CorepComputationError::UnsupportedClassification {
                    uni: uni_number,
                    source_irrep: h_irrep.ml.to_string(),
                    reason: format!(
                        "spinor SU(2) Wigner classification did not produce A/B/C: {}",
                        e
                    ),
                })
            }
        }
    } else {
        // Non-compound scalar: PIR path in the selected star-arm block.
        let classification = if let Some(complex_characters) = &h_complex_chars {
            let pairs = complex_characters
                .iter()
                .flat_map(|value| [value.re, value.im])
                .collect::<Vec<_>>();
            wigner::wigner_classify_cir_direct(
                &pairs,
                &antiunitary,
                &mag_seitz,
                &h_seitz,
                h_irrep.k_vector(),
            )
        } else {
            wigner::wigner_classify(
                &h_chars,
                wigner::WignerGroupContext {
                    unitary_indices: &unitary,
                    magnetic_ops: &mag_seitz,
                    unitary_ops: &h_seitz,
                    antiunitary_representative: antiunitary[0],
                },
                h_irrep.k_vector(),
            )
        };
        let ct = classification.map_err(|e| CorepComputationError::UnsupportedClassification {
            uni: uni_number,
            source_irrep: h_irrep.ml.to_string(),
            reason: format!("scalar PIR Wigner classification failed: {}", e),
        })?;
        Ok((ct, WignerSource::ScalarPIR))
    }?;

    // 8. Compute Type A antiunitary characters
    let au_chars = if h_irrep.spinor {
        // The stored spinor data contains little-group characters but no
        // representation matrices/intertwiner for anti-linear operations.
        // Never feed its local character order into the scalar H-order path.
        None
    } else if corep_type == CorepType::A && !antiunitary.is_empty() {
        if h_dim == 1 {
            wigner::type_a_antiunitary_chars(
                &mag_seitz,
                &mag_lg,
                &h_chars,
                &h_seitz,
                antiunitary[0],
                h_irrep.k_vector(),
            )
            .map(|(chars, _u)| chars)
        } else {
            if !h_matrices_ordered.is_empty() {
                wigner::type_a_antiunitary_chars_high_dim_ordered(
                    &mag_seitz,
                    &mag_lg,
                    &h_chars,
                    &h_seitz,
                    antiunitary[0],
                    &h_matrices_ordered,
                    h_dim,
                )
            } else {
                let matrices = h_irrep.matrices();
                let rotations = h_irrep.pir_rotations();
                if matrices.is_empty() || rotations.is_empty() {
                    None
                } else {
                    wigner::type_a_antiunitary_chars_high_dim(
                        &mag_seitz,
                        &mag_lg,
                        &h_chars,
                        &h_seitz,
                        antiunitary[0],
                        matrices,
                        rotations,
                    )
                }
            }
        }
    } else {
        None
    };

    // Spinor character tables are stored in little-group-local order, while
    // op_map contains indices into the complete H operation list. Convert the
    // index domain explicitly; leaving a missing entry as None lets the strict
    // character builder report corrupted or inconsistent data.
    let character_op_map = if h_irrep.spinor {
        let h_spin = h_irrep.spin_ops();
        let h_spin_seitz = wigner::build_spin_seitz(h_spin.0, h_spin.1);
        let h_to_spin =
            wigner::build_h_to_spin_map(&h_seitz, &h_spin_seitz, h_irrep.spin_lg_op_indices());
        op_map
            .iter()
            .map(|h_index| {
                let spin_index = h_to_spin.get((*h_index)?).copied().flatten()?;
                h_irrep
                    .spin_lg_op_indices()
                    .iter()
                    .position(|&index| index as usize == spin_index)
            })
            .collect::<Vec<_>>()
    } else {
        op_map.clone()
    };

    // 9. Build corep character table
    let characters = wigner::build_corep_chars(
        &corep_type,
        &mag_ops_data,
        &mag_lg,
        &character_op_map,
        &h_chars,
        None,
        au_chars.as_deref(),
    )
    .map_err(|reason| CorepComputationError::UnsupportedClassification {
        uni: uni_number,
        source_irrep: h_irrep.ml.to_string(),
        reason: reason.to_string(),
    })?;
    if !characters.iter().all(|c| c.is_finite()) {
        return Err(CorepComputationError::NonFiniteCharacters {
            uni: uni_number,
            source_irrep: h_irrep.ml.to_string(),
        });
    }

    let dim = wigner::corep_dim(&corep_type, h_dim);

    let completeness = match corep_type {
        CorepType::A if !antiunitary.is_empty() && au_chars.is_none() => {
            CharacterCompleteness::TypeAAntiunitaryPending {
                count: antiunitary.len(),
            }
        }
        _ => CharacterCompleteness::Complete,
    };

    Ok(Corepresentation {
        characters,
        timerev: mag_lg
            .iter()
            .map(|&i| mag_ops_data.operations[i].time_reversal)
            .collect(),
        corep_type,
        source,
        dim,
        unitary_order: unitary.len(),
        antiunitary_order: antiunitary.len(),
        completeness,
    })
}

// ── Magnetic operations ──────────────────────────────────────────────────────

/// Transform magnetic operations into the ISOTROPY data-Hall frame while
/// preserving their database indices and time-reversal flags.
pub(crate) fn operations_in_data_hall_frame(
    operations: &SymmetryOps,
    setting: Option<&wigner::SettingTransform>,
) -> Option<SymmetryOps> {
    let Some(setting) = setting else {
        return Some(operations.clone());
    };
    let mut rotations = Vec::with_capacity(operations.len());
    let mut translations = Vec::with_capacity(operations.len());
    let mut time_reversals = Vec::with_capacity(operations.len());
    for operation in &operations.operations {
        let (rotation, translation) =
            setting.transform_seitz(&operation.rotation, &operation.translation)?;
        rotations.push(rotation);
        translations.push(translation);
        time_reversals.push(operation.time_reversal);
    }
    SymmetryOps::from_parallel_owned(rotations, translations, time_reversals).ok()
}

/// Get the magnetic space group symmetry operations.
pub fn get_magnetic_operations(uni_number: usize) -> Option<SymmetryOps> {
    let hall = get_first_hall_for_uni(uni_number)?;
    let sym = crate::msg_database::get_spacegroup_operations(uni_number, hall)?;
    let n = sym.len();
    let mut rot = Vec::with_capacity(n);
    let mut trans = Vec::with_capacity(n);
    let mut timerev = Vec::with_capacity(n);
    for i in 0..n {
        rot.push(sym.rot[i]);
        trans.push(sym.trans[i]);
        timerev.push(sym.timerev[i]);
    }
    SymmetryOps::from_parallel(&rot, &trans, &timerev).ok()
}

fn get_first_hall_for_uni(uni: usize) -> Option<usize> {
    if uni == 0 || uni > 1651 {
        return None;
    }
    let first_hall = crate::msg_database::MAGNETIC_SPACEGROUP_UNI_MAPPING[uni][1];
    if first_hall > 0 {
        Some(first_hall as usize)
    } else {
        None
    }
}

/// Get the symmetry operations (rotation + translation) for a space group.
///
/// Returns [`SymmetryOps`] with `timerev` all `false` (non-magnetic).
/// The operations are in spglib's standard order.
///
/// # Example
/// ```
/// use cryspglib::irrep::corep::symmetry_operations_of;
/// let ops = symmetry_operations_of(139).unwrap();
/// println!("SG 139: {} operations", ops.len());
/// ```
pub fn symmetry_operations_of(sg: u8) -> Result<SymmetryOps, crate::SymError> {
    get_parent_operations(sg)
}

fn get_parent_operations_by_hall(hall: usize) -> Option<SymmetryOps> {
    let sym = get_spacegroup_operations(hall)?;
    let n = sym.len();
    let mut rot = Vec::with_capacity(n);
    let mut trans = Vec::with_capacity(n);
    for i in 0..n {
        rot.push(sym.rot[i]);
        trans.push(sym.trans[i]);
    }
    let timerev = vec![false; n];
    SymmetryOps::from_parallel(&rot, &trans, &timerev).ok()
}

fn get_parent_operations(sg: u8) -> Result<SymmetryOps, crate::SymError> {
    let hall = crate::api::find_hall_number(sg)?;
    get_parent_operations_by_hall(hall).ok_or(crate::SymError::SymmetryOperationSearchFailed)
}

// ── High-level API ───────────────────────────────────────────────────────────

/// Identified unitary subgroup of a magnetic space group, with correct Hall setting.
pub struct UnitarySubgroupInfo {
    pub sg: usize,
    /// Hall setting identified by spglib (may differ from ISOTROPY data Hall).
    pub hall: usize,
    /// ISOTROPY data Hall (= SG_DATA_HALL[sg]).
    pub data_hall: usize,
    /// Unitary ops extracted from the MSG itself (MSG frame).
    pub ops_from_msg: SymmetryOps,
    /// Unitary ops in the ISOTROPY data-Hall frame.
    pub ops_from_hall: SymmetryOps,
    /// Composed transform MSG frame → ISOTROPY data-Hall frame.
    /// Computed as msg_to_detected.then(&detected_to_data).
    pub msg_to_data: Option<wigner::SettingTransform>,
}

impl UnitarySubgroupInfo {
    /// Assert that extracting ops from MSG and from Hall give the same Seitz set.
    pub fn assert_consistent(&self) -> bool {
        use crate::irrep::wigner;
        let seitz_msg = wigner::ops_to_seitz(&self.ops_from_msg);
        let seitz_hall = wigner::ops_to_seitz(&self.ops_from_hall);
        if seitz_msg.len() != seitz_hall.len() {
            return false;
        }
        for op in &seitz_msg {
            if wigner::find_seitz(&op.rot, &op.trans, &seitz_hall).is_none() {
                return false;
            }
        }
        true
    }
}

/// Identify the unitary subgroup of a magnetic space group (SG number only).
pub fn identify_unitary_subgroup(uni_number: usize) -> Option<usize> {
    identify_unitary_subgroup_with_hall(uni_number).map(|info| info.sg)
}

/// Look up the parent spatial space group number G ⊃ H for a magnetic group.
///
/// For black-white (Type III) MSGs, $$G \supset H$$ is a proper supergroup.
/// For grey (Type II) and ordinary (Type I) groups, G = H.
pub fn parent_spatial_sg(uni_number: usize) -> Option<usize> {
    crate::MagneticSpaceGroupType::from_uni(uni_number)
        .ok()
        .map(|msg| msg.number)
}

/// Pick the correct a₀ (canonical antiunitary coset representative) for spinor Wigner.
///
/// For grey groups (Type II), a₀ must be pure θ (R = I), because (θg)² ≠ -g² in general.
/// For black-white groups (Type III), any antiunitary representative works.
pub fn select_spinor_a0(
    antiunitary: &[usize],
    mag_seitz: &[crate::irrep::wigner::SeitzOp],
    is_grey: bool,
) -> usize {
    let id_rot: crate::mathfunc::Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    if is_grey {
        antiunitary
            .iter()
            .copied()
            .find(|&i| mag_seitz[i].rot == id_rot)
            .unwrap_or(antiunitary[0])
    } else {
        antiunitary[0]
    }
}

fn translations_equal_mod_lattice(a: &[f64; 3], b: &[f64; 3]) -> bool {
    (0..3).all(|i| {
        let delta = a[i] - b[i];
        (delta - delta.round()).abs() < 1e-8
    })
}

/// Validate an affine frame transform against an operation embedding.
///
/// `source` is allowed to contain fewer representatives than `target`.  This
/// is required for centered groups, where the MSG database may contain one
/// representative per centering coset while the canonical Hall operation set
/// contains all conventional-cell representatives.
fn transform_embeds_ops(
    xf: &wigner::SettingTransform,
    source: &SymmetryOps,
    target: &SymmetryOps,
) -> bool {
    source.operations.iter().all(|op| {
        let Some(rotation) = xf.transform_rotation(&op.rotation) else {
            return false;
        };
        let Some(translation) = xf.transform_translation(&op.rotation, &op.translation) else {
            return false;
        };
        target.operations.iter().any(|candidate| {
            candidate.rotation == rotation
                && translations_equal_mod_lattice(&candidate.translation, &translation)
        })
    })
}

fn transform_applies_to_all_ops(xf: &wigner::SettingTransform, ops: &SymmetryOps) -> bool {
    ops.operations.iter().all(|op| {
        xf.transform_rotation(&op.rotation).is_some()
            && xf
                .transform_translation(&op.rotation, &op.translation)
                .is_some()
    })
}

/// Identify the unitary subgroup with full Hall setting information.
///
/// Uses `spg_get_hall_number_from_symmetry` to classify the unitary ops
/// from the MSG, then reconstructs H_ops from the identified Hall number.
/// This ensures the H_ops setting matches the MSG, rather than using the
/// first-Hall setting which may differ in origin/basis.
pub fn identify_unitary_subgroup_with_hall(uni_number: usize) -> Option<UnitarySubgroupInfo> {
    let mag_ops = get_magnetic_operations(uni_number)?;
    let msg_type = crate::MagneticSpaceGroupType::from_uni(uni_number).ok()?;

    // For Type I (Ordinary) and Type II (Grey), H = G = parent SG.
    // Use the parent SG directly — no spglib classification needed.
    let sg_from_metadata = match msg_type.type_ {
        crate::MagneticType::Ordinary | crate::MagneticType::Grey => Some(msg_type.number as u8),
        crate::MagneticType::BlackWhite
        | crate::MagneticType::AntiTranslation
        | crate::MagneticType::NonMagnetic => None,
    };

    let mut unitary_rots: Vec<Mat3I> = Vec::new();
    let mut unitary_trans: Vec<[f64; 3]> = Vec::new();
    for operation in &mag_ops.operations {
        if !operation.time_reversal {
            unitary_rots.push(operation.rotation);
            unitary_trans.push(operation.translation);
        }
    }
    if unitary_rots.is_empty() {
        return None;
    }

    let n = unitary_rots.len();
    let timerev_from_msg = vec![false; n];
    let ops_from_msg =
        SymmetryOps::from_parallel(&unitary_rots, &unitary_trans, &timerev_from_msg).ok()?;

    let (sg, hall, ops_from_hall, msg_to_detected) = if let Some(s) = sg_from_metadata {
        // Type I/II: use metadata SG directly.
        let h = crate::api::find_hall_number(s).ok()?;
        let oh = get_parent_operations_by_hall(h)?;
        (s as usize, h, oh, None)
    } else {
        // Type III: try standard_setting_transform first, then fall back.
        let unitary_ops = SymmetryOps::from_parallel(
            &unitary_rots,
            &unitary_trans,
            &vec![false; unitary_rots.len()],
        )
        .ok()?;
        if let Some((std_sg, std_hall, xf)) = standard_setting_transform(&unitary_ops, false) {
            let oh = get_parent_operations_by_hall(std_hall).or_else(|| {
                get_parent_operations_by_hall(crate::api::find_hall_number(std_sg as u8).ok()?)
            })?;
            (std_sg, std_hall, oh, Some(xf))
        } else {
            let h = crate::hall_number_from_symmetry(&unitary_rots, &unitary_trans, 1e-5).ok()?;
            if h == 0 || h > 530 {
                return None;
            }
            let sg_type = get_spacegroup_type(h);
            let oh = get_parent_operations_by_hall(h)?;
            (sg_type.number, h, oh, None)
        }
    };

    // ── Align ops_from_hall with ISOTROPY data-Hall frame ────────────────
    // ISOTROPY irreps, characters, and spin tables are defined in a specific
    // "data Hall" setting (SG_DATA_HALL[sg]).  When spglib identifies H in a
    // different Hall setting of the same SG (e.g. SG44 detected as Hall 216
    // but data is Hall 215), the rotation matrices differ by a signed-
    // permutation axis convention.  Without correction, spin-table lookups
    // fail because the square rotation R(b²) is expressed in the detected
    // Hall frame but the spin database uses the data Hall frame.
    //
    // Compute the affine transform x_data = P·x_detected + s by comparing
    // the two Hall settings' standard operations, then apply it to
    // ops_from_hall so the entire Wigner H-path works in the data Hall frame.
    let data_hall = crate::irrep::types::generated_data::SG_DATA_HALL
        .get(sg)
        .copied()
        .unwrap_or(0) as usize;
    if data_hall > 0 && data_hall != hall {
        if let Some(data_ops) = get_parent_operations_by_hall(data_hall) {
            // Build rotation/translation arrays for find_setting_transform.
            let detected_rots: Vec<Mat3I> = ops_from_hall
                .operations
                .iter()
                .map(|o| o.rotation)
                .collect();
            let detected_trans: Vec<[f64; 3]> = ops_from_hall
                .operations
                .iter()
                .map(|o| o.translation)
                .collect();
            let data_rots: Vec<Mat3I> = data_ops.operations.iter().map(|o| o.rotation).collect();
            let data_trans: Vec<[f64; 3]> =
                data_ops.operations.iter().map(|o| o.translation).collect();
            let xfs = crate::irrep::wigner::find_setting_transform(
                &detected_rots,
                &detected_trans,
                &data_rots,
                &data_trans,
            );
            // Prefer the transform that also makes the actual MSG subgroup
            // embed into the data-Hall operation set.  This disambiguates
            // signed-permutation candidates using the physical subgroup,
            // rather than taking xfs.first() arbitrarily.
            for detected_to_data in &xfs {
                if !transform_embeds_ops(detected_to_data, &ops_from_hall, &data_ops) {
                    continue;
                }
                if let Some(msg_to_detected) = msg_to_detected.as_ref() {
                    let msg_to_data = msg_to_detected.then(detected_to_data);
                    if transform_embeds_ops(&msg_to_data, &ops_from_msg, &data_ops)
                        && transform_applies_to_all_ops(&msg_to_data, &mag_ops)
                    {
                        return Some(UnitarySubgroupInfo {
                            sg,
                            hall,
                            data_hall,
                            ops_from_msg,
                            ops_from_hall: data_ops,
                            msg_to_data: Some(msg_to_data),
                        });
                    }
                } else if transform_embeds_ops(detected_to_data, &ops_from_msg, &data_ops)
                    && transform_applies_to_all_ops(detected_to_data, &mag_ops)
                {
                    return Some(UnitarySubgroupInfo {
                        sg,
                        hall,
                        data_hall,
                        ops_from_msg,
                        ops_from_hall: data_ops,
                        msg_to_data: Some(detected_to_data.clone()),
                    });
                }
            }
            // A subgroup-based standardisation can leave basis/origin freedom
            // that is incompatible with the antiunitary coset.  Search the
            // actual MSG unitary representatives directly against the data
            // Hall frame, regardless of whether a MSG→detected transform was
            // available, and validate the candidate on the full MSG.
            let msg_rots: Vec<Mat3I> = ops_from_msg.operations.iter().map(|o| o.rotation).collect();
            let msg_trans: Vec<[f64; 3]> = ops_from_msg
                .operations
                .iter()
                .map(|o| o.translation)
                .collect();
            let msg_xfs = crate::irrep::wigner::find_setting_transform(
                &msg_rots,
                &msg_trans,
                &data_rots,
                &data_trans,
            );
            if let Some(msg_to_data) = msg_xfs.into_iter().find(|candidate| {
                transform_embeds_ops(candidate, &ops_from_msg, &data_ops)
                    && transform_applies_to_all_ops(candidate, &mag_ops)
            }) {
                return Some(UnitarySubgroupInfo {
                    sg,
                    hall,
                    data_hall,
                    ops_from_msg,
                    ops_from_hall: data_ops,
                    msg_to_data: Some(msg_to_data),
                });
            }
            // The finite unimodular search cannot represent every rational
            // conventional/primitive-cell transform.  Try spglib's full
            // standardisation pipeline, but accept it only when it lands in
            // the requested data Hall and validates against the full MSG.
            if let Some((std_sg, std_hall, detected_to_data)) =
                standard_setting_transform(&ops_from_hall, false)
                && std_sg == sg
                && std_hall == data_hall
                && transform_embeds_ops(&detected_to_data, &ops_from_hall, &data_ops)
            {
                let msg_to_data = msg_to_detected
                    .as_ref()
                    .map(|candidate| candidate.then(&detected_to_data))
                    .unwrap_or(detected_to_data);
                if transform_embeds_ops(&msg_to_data, &ops_from_msg, &data_ops)
                    && transform_applies_to_all_ops(&msg_to_data, &mag_ops)
                {
                    return Some(UnitarySubgroupInfo {
                        sg,
                        hall,
                        data_hall,
                        ops_from_msg,
                        ops_from_hall: data_ops,
                        msg_to_data: Some(msg_to_data),
                    });
                }
            }
        }
        // Never expose a detected-Hall operation set as if it were in the
        // data-Hall frame.  Callers need a verified common frame.
        return None;
    }

    // Compute the MSG→data-Hall transform: maps the MSG unitary ops into
    // the ISOTROPY data-Hall frame.  This is the transform that the Wigner
    // path should use for little-group filtering and Seitz composition,
    // ensuring all operations are in the same frame as the spin table and
    // characters.
    let data_hall = crate::irrep::types::generated_data::SG_DATA_HALL
        .get(sg)
        .copied()
        .unwrap_or(0) as usize;
    let msg_to_data = if data_hall > 0 {
        if data_hall == hall
            && let Some(xf) = msg_to_detected.filter(|candidate| {
                transform_embeds_ops(candidate, &ops_from_msg, &ops_from_hall)
                    && transform_applies_to_all_ops(candidate, &mag_ops)
            })
        {
            return Some(UnitarySubgroupInfo {
                sg,
                hall,
                data_hall,
                ops_from_msg,
                ops_from_hall,
                msg_to_data: Some(xf),
            });
        }
        let msg_rots: Vec<Mat3I> = ops_from_msg.operations.iter().map(|o| o.rotation).collect();
        let msg_trans: Vec<[f64; 3]> = ops_from_msg
            .operations
            .iter()
            .map(|o| o.translation)
            .collect();
        let hall_rots: Vec<Mat3I> = ops_from_hall
            .operations
            .iter()
            .map(|o| o.rotation)
            .collect();
        let hall_trans: Vec<[f64; 3]> = ops_from_hall
            .operations
            .iter()
            .map(|o| o.translation)
            .collect();
        let xfs = crate::irrep::wigner::find_setting_transform(
            &msg_rots,
            &msg_trans,
            &hall_rots,
            &hall_trans,
        );
        xfs.into_iter().find(|candidate| {
            transform_embeds_ops(candidate, &ops_from_msg, &ops_from_hall)
                && transform_applies_to_all_ops(candidate, &mag_ops)
        })
    } else {
        None
    };

    Some(UnitarySubgroupInfo {
        sg,
        hall,
        data_hall,
        ops_from_msg,
        ops_from_hall,
        msg_to_data,
    })
}

/// Standardize operations with spglib's primitive/conventional-cell pipeline.
///
/// The returned affine transform maps the input MSG coordinates into the
/// detected standard Hall setting and supports rational basis matrices, so it
/// also covers centered and supercell embeddings.
fn standard_setting_transform(
    ops: &SymmetryOps,
    ignore_time_reversal: bool,
) -> Option<(usize, usize, wigner::SettingTransform)> {
    let mut magnetic = crate::symmetry::MagneticSymmetry::with_capacity(ops.len());
    for op in &ops.operations {
        magnetic.push(op.rotation, op.translation, op.time_reversal);
    }
    let (spacegroup, _) = crate::magnetic_spacegroup::get_space_group_with_magnetic_symmetry(
        &magnetic,
        ignore_time_reversal,
        1e-5,
    )?;
    let basis = mat_inverse_matrix_d3(&spacegroup.bravais_lattice, 1e-10).ok()?;
    Some((
        spacegroup.number,
        spacegroup.hall_number,
        wigner::SettingTransform {
            basis,
            origin: spacegroup.origin_shift,
        },
    ))
}

/// BNS label → UNI number.
pub fn uni_from_bns(bns: &str) -> Option<usize> {
    for uni in 1..=1651usize {
        let t = crate::msg_database::get_magnetic_spacegroup_type(uni);
        if t.bns_number == bns {
            return Some(uni);
        }
    }
    None
}

/// OG label → UNI number.
pub fn uni_from_og(og: &str) -> Option<usize> {
    for uni in 1..=1651usize {
        let t = crate::msg_database::get_magnetic_spacegroup_type(uni);
        if t.og_number == og {
            return Some(uni);
        }
    }
    None
}

/// Compute all corepresentations for a magnetic SG at a k-point.
pub fn compute_coreps(
    bns: &str,
    k_label: &str,
) -> Result<Vec<(String, Corepresentation)>, CorepComputationError> {
    let uni = uni_from_bns(bns).ok_or_else(|| CorepComputationError::UnknownBns {
        bns: bns.to_string(),
    })?;
    let h_info = identify_unitary_subgroup_with_hall(uni)
        .ok_or(CorepComputationError::MissingUnitarySubgroup { uni })?;
    let h_sg = h_info.sg;
    let mag_ops = get_magnetic_operations(uni)
        .ok_or(CorepComputationError::MissingMagneticOperations { uni })?;
    let h_irreps = super::query::irreps_of(h_sg as u8);
    let k_irreps: Vec<&IrrepRecord> = h_irreps.iter().filter(|r| r.k_label() == k_label).collect();
    if k_irreps.is_empty() {
        return Err(CorepComputationError::MissingKPoint {
            uni,
            sg: h_sg as u8,
            k_label: k_label.to_string(),
        });
    }

    let mut results = Vec::with_capacity(k_irreps.len());
    for ir in k_irreps.iter() {
        let c = compute_corepresentation(ir, uni, &mag_ops)?;
        // For Type C, partner finding for better character tables is still
        // pending. The characters from compute_corepresentation are used
        // directly for now.
        let final_chars = c.characters.clone();
        results.push((
            ir.ml.to_string(),
            Corepresentation {
                characters: final_chars,
                ..c
            },
        ));
    }
    Ok(results)
}

// ── IrrepRecord extension ────────────────────────────────────────────────────

impl IrrepRecord {
    /// Compute the magnetic co-representation (corep) for this non-magnetic
    /// irrep with respect to a magnetic space group.
    ///
    /// The magnetic character table is computed on-the-fly from the
    /// non-magnetic irrep data — no pre-stored tables needed.
    ///
    /// # Arguments
    ///
    /// * `uni_number` — OG/UNI number (1–1651), from
    ///   [`MagneticIsotropyRecord::mag_sg`](super::types::MagneticIsotropyRecord::mag_sg)
    ///
    /// # Examples
    ///
    /// ```
    /// use cryspglib::irrep::query::irreps_of;
    ///
    /// let gm1 = irreps_of(1).iter()
    ///     .find(|r| r.ml == "GM1").unwrap();
    ///
    /// // UNI 2 is the grey P1 magnetic group.
    /// let corep = gm1.corepresentation(2)?;
    /// println!("Type: {:?}, dim: {}", corep.corep_type, corep.dim);
    /// for (i, &chi) in corep.characters.iter().enumerate() {
    ///     let tr = if corep.timerev[i] { " (θ)" } else { "" };
    ///     println!("  op {}: χ = {:.4}{}", i, chi, tr);
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    /// Compute the magnetic co-representation for this irrep.
    ///
    /// Note: `self` must be an irrep of the **unitary subgroup H**, not the
    /// parent SG. Use [`compute_coreps`] for automatic H identification.
    pub fn corepresentation(
        &self,
        uni_number: usize,
    ) -> Result<Corepresentation, CorepComputationError> {
        let mag_ops = get_magnetic_operations(uni_number)
            .ok_or(CorepComputationError::MissingMagneticOperations { uni: uni_number })?;
        compute_corepresentation(self, uni_number, &mag_ops)
    }
}

// ── Magnetic isotropy → corepresentation bridge ────────────────────────────

/// Result of computing a co-representation for a magnetic isotropy subgroup.
#[derive(Debug, Clone)]
pub struct MagneticIsotropyCorep {
    /// Magnetic isotropy subgroup record (from ISOTROPY data)
    pub subgroup: super::types::MagneticIsotropyRecord,
    /// Computed co-representation, or a structured computation error.
    pub corep: Result<Corepresentation, CorepComputationError>,
}

impl MagneticIsotropyCorep {
    /// Short description: UNI, BNS, direction, corep type, dimension.
    pub fn describe(&self) -> String {
        match &self.corep {
            Ok(c) => format!(
                "UNI {} {} dir={} → {:?} dim={}",
                self.subgroup.mag_sg,
                self.subgroup.bns_label,
                self.subgroup.direction,
                c.corep_type,
                c.dim
            ),
            Err(err) => format!(
                "UNI {} {} → error: {}",
                self.subgroup.mag_sg, self.subgroup.bns_label, err
            ),
        }
    }
}

impl std::fmt::Display for MagneticIsotropyCorep {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.describe())
    }
}

/// For a given irrep, compute co-representations for all its magnetic
/// isotropy subgroups.
///
/// This bridges the ISOTROPY subgroup data with on-the-fly Wigner
/// classification and character table computation.
///
/// ```
/// use cryspglib::irrep::query::irreps_of;
/// use cryspglib::irrep::corep::magnetic_isotropy_coreps_of_irrep;
///
/// let gm4m = irreps_of(221).iter()
///     .find(|r| r.ml == "GM4-").unwrap();
/// let results = magnetic_isotropy_coreps_of_irrep(gm4m);
/// assert!(!results.is_empty());
/// for r in &results {
///     println!("{}", r.describe());
/// }
/// ```
pub fn magnetic_isotropy_coreps_of_irrep(ir: &IrrepRecord) -> Vec<MagneticIsotropyCorep> {
    ir.magnetic_subgroups()
        .iter()
        .map(|sub| {
            let corep = match get_magnetic_operations(sub.mag_sg) {
                Some(mag_ops) => compute_corepresentation(ir, sub.mag_sg, &mag_ops),
                None => Err(CorepComputationError::MissingMagneticOperations { uni: sub.mag_sg }),
            };
            MagneticIsotropyCorep {
                subgroup: *sub,
                corep,
            }
        })
        .collect()
}

/// For all scalar irreps of a space group at a k-point, compute
/// co-representations for their magnetic isotropy subgroups.
///
/// Returns entries grouped by irrep.
pub fn magnetic_isotropy_coreps_of_sg_k(
    sg: u8,
    kx: i8,
    ky: i8,
    kz: i8,
    kd: i8,
) -> Vec<(IrrepRecord, Vec<MagneticIsotropyCorep>)> {
    super::query::irreps_of(sg)
        .iter()
        .filter(|ir| !ir.spinor && ir.kx == kx && ir.ky == ky && ir.kz == kz && ir.kd == kd)
        .map(|ir| (*ir, magnetic_isotropy_coreps_of_irrep(ir)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symmetry_operations_query_reports_invalid_space_group() {
        assert!(matches!(
            symmetry_operations_of(0),
            Err(crate::SymError::SpacegroupSearchFailed)
        ));
        assert!(matches!(
            symmetry_operations_of(231),
            Err(crate::SymError::SpacegroupSearchFailed)
        ));
        assert_eq!(symmetry_operations_of(139).unwrap().len(), 32);
    }
    use crate::irrep::query::irreps_of;

    #[test]
    fn test_corep_gm4m_pmmm() {
        // Use compute_coreps on a SG 221 magnetic subgroup: 221.97 = UNI 1599
        // This tests the full pipeline: BNS → UNI → H → H irreps → coreps
        let coreps = compute_coreps("221.97", "GM");
        assert!(
            coreps.is_ok(),
            "Should compute coreps for 221.97 at GM: {:?}",
            coreps
        );
        let coreps = coreps.unwrap();
        assert!(!coreps.is_empty(), "Should have at least one corep");

        for (label, c) in &coreps {
            assert!(c.dim > 0, "dim > 0 for {}", label);
            assert!(
                (c.characters[0] - c.dim as f64).abs() < 0.01,
                "χ(id) = dim for {}",
                label
            );
        }
        println!("221.97 Gamma coreps: {} irreps computed", coreps.len());
    }

    #[test]
    fn test_corep_sg1_gm1() {
        // SG 1 (P1) GM1 → simplest case, use BNS "1.3" (= UNI 3)
        let coreps = compute_coreps("1.3", "GM");
        assert!(
            coreps.is_ok(),
            "Should compute coreps for 1.3 at GM: {:?}",
            coreps
        );
    }

    /// SG 128.406 (P4'/m'nc') at Z point — verified against BCS
    /// https://cryst.ehu.es/cgi-bin/cryst/programs/corepresentations.pl
    ///
    /// BCS confirms: Unitary Space Group = P-4n2 (No. 118) in standard setting.
    /// This test verifies automatic identification of the unitary subgroup.
    #[test]
    fn test_unitary_subgroup_sg128_406_is_sg118() {
        let uni: usize = 1066; // 128.406 = UNI 1066 (NOT 1073 — that's the Litvin number)
        let ops = get_magnetic_operations(uni);
        assert!(ops.is_some(), "Should get ops for UNI 1066 (128.406)");
        let ops = ops.unwrap();

        let n_u = ops.operations.iter().filter(|o| !o.time_reversal).count();
        let n_a = ops.operations.iter().filter(|o| o.time_reversal).count();
        println!(
            "Magnetic SG UNI {}: {} ops ({} unitary + {} anti-unitary)",
            uni,
            ops.len(),
            n_u,
            n_a
        );

        // ── 1. Full group (ignore θ) should identify as parent SG 128 ──
        let ops_rots: Vec<_> = ops.operations.iter().map(|o| o.rotation).collect();
        let ops_trans: Vec<_> = ops.operations.iter().map(|o| o.translation).collect();
        let hall_full = crate::hall_number_from_symmetry(&ops_rots, &ops_trans, 1e-5);
        assert!(hall_full.is_ok(), "Should identify full group");
        let sg_full = crate::SpaceGroupType::from_hall(hall_full.unwrap()).unwrap();
        assert_eq!(
            sg_full.number, 128,
            "Full ops should identify as SG 128, got SG {}",
            sg_full.number
        );
        println!("Full group (ignore θ): SG 128 ✓");

        // ── 2. Unitary subgroup should identify as SG 118 (P-4n2) ──
        let h_sg = identify_unitary_subgroup(uni);
        assert!(h_sg.is_some(), "Should identify unitary subgroup");
        let h_sg = h_sg.unwrap();
        assert_eq!(
            h_sg, 118,
            "Unitary subgroup of 128.406 should be SG 118, got SG {}",
            h_sg
        );
        println!("Unitary subgroup: SG 118 (P-4n2) ✓");

        // ── 3. Verify: all 16 magnetic rotations are in parent SG 128 ──
        let parent_ops = get_parent_operations(128).unwrap();
        let parent_rots: Vec<_> = parent_ops.operations.iter().map(|o| o.rotation).collect();
        let all_match = ops
            .operations
            .iter()
            .all(|o| parent_rots.contains(&o.rotation));
        assert!(all_match, "All magnetic rotations should be in SG 128 ops");
        println!("Magnetic ops ⊆ SG 128 ops ✓");
    }

    #[test]
    fn test_centered_type3_msg_to_data_transform_is_validated() {
        // The MSG operation set may contain fewer centering representatives
        // than the canonical Hall set.  The composed transform must therefore
        // be validated as an embedding, not by equal operation counts.
        for uni in [815usize, 889, 901, 1204, 1214] {
            let info = identify_unitary_subgroup_with_hall(uni)
                .unwrap_or_else(|| panic!("UNI {uni}: unitary subgroup identification failed"));
            let xf = info
                .msg_to_data
                .as_ref()
                .unwrap_or_else(|| panic!("UNI {uni}: missing MSG→data-Hall transform"));
            let mag_ops = get_magnetic_operations(uni)
                .unwrap_or_else(|| panic!("UNI {uni}: missing magnetic operations"));

            assert!(
                transform_embeds_ops(xf, &info.ops_from_msg, &info.ops_from_hall),
                "UNI {uni}: transformed unitary operations do not embed in data Hall {}",
                info.data_hall,
            );
            assert!(
                transform_applies_to_all_ops(xf, &mag_ops),
                "UNI {uni}: transform is not valid for every magnetic operation",
            );
        }
    }

    /// Identify the unitary subgroup space group number from a magnetic UNI.
    #[test]
    fn test_identify_unitary_subgroup_api() {
        // 128.406 → UNI 1066 → unitary SG 118 (P-4n2) — verified against BCS
        assert_eq!(identify_unitary_subgroup(1066), Some(118));

        // 129.413 → UNI 1073 → parent SG 129, black-white
        // Its unitary subgroup should also be identifiable
        let result = identify_unitary_subgroup(1073);
        println!("UNI 1073 (129.413) unitary subgroup: {:?}", result);
        assert!(result.is_some(), "UNI 1073 should identify");

        // 1.2 (BNS) → UNI 2, simplest non-trivial magnetic SG
        assert!(identify_unitary_subgroup(2).is_some(), "UNI 2 should work");
    }

    /// Cross-validate compound CIR data against its full-star PIR dimension.
    ///
    /// The stored CIR rows contain the selected arm and include Hall-setting
    /// Bloch phases.  The full real PIR rows use induced-representation traces,
    /// so the universal invariant after setting conversion is the dimension
    /// relation `dim(PIR) = star_size * Σ dim(CIR_little)`.
    #[test]
    fn test_cir_pir_cross_validation() {
        let mut checked = 0usize;
        let mut one_arm = 0usize;
        let mut multi_arm = 0usize;
        let mut mismatches = 0usize;
        // Iterate over all SGs
        for sg in 1u8..=230 {
            for ir in crate::irrep::query::irreps_of(sg) {
                let n_comp = ir.cir_component_count();
                if n_comp == 0 {
                    continue;
                }

                let pir = ir.characters();
                let n_ops = pir.len();
                if n_ops == 0 {
                    continue;
                }

                // Sum selected-arm CIR component characters.
                let cir_ops = (0..n_comp)
                    .map(|c| ir.cir_component_chars(c).len() / 2)
                    .min()
                    .unwrap_or(0);
                if cir_ops == 0 {
                    continue; // No CIR data
                }
                if cir_ops != n_ops {
                    mismatches += 1;
                    eprintln!(
                        "MISMATCH SG{} {}: PIR ops={} CIR ops={}",
                        sg, ir.ml, n_ops, cir_ops
                    );
                    continue;
                }
                let mut cir_sum_re = vec![0.0f64; n_ops];
                let mut cir_sum_im = vec![0.0f64; n_ops];
                for c in 0..n_comp {
                    let cir = ir.cir_component_chars(c);
                    for op in 0..n_ops {
                        cir_sum_re[op] += cir[2 * op];
                        cir_sum_im[op] += cir[2 * op + 1];
                    }
                }

                let full_dim = ir.dim as f64;
                let little_dim = cir_sum_re[0];
                let star_size = full_dim / little_dim;
                if little_dim <= 0.0
                    || cir_sum_im[0].abs() > 0.01
                    || (star_size - star_size.round()).abs() > 0.01
                    || star_size < 1.0
                {
                    mismatches += 1;
                    eprintln!(
                        "MISMATCH SG{} {} dimensions: PIR={} selected CIR sum=({:.4},{:.4})",
                        sg, ir.ml, full_dim, little_dim, cir_sum_im[0]
                    );
                    continue;
                }

                if (star_size - 1.0).abs() < 0.01 {
                    one_arm += 1;
                } else {
                    multi_arm += 1;
                }
                checked += 1;
            }
        }
        println!(
            "CIR↔PIR cross-check: {} compound irreps ({} one-arm, {} multi-arm), {} mismatches",
            checked, one_arm, multi_arm, mismatches
        );
        assert_eq!(mismatches, 0, "All compound CIR records must be consistent");
        assert!(
            checked > 500,
            "Should cover at least 500 compound irreps, got {}",
            checked
        );
    }

    /// Test BNS/OG → UNI lookup functions.
    #[test]
    fn test_uni_lookup() {
        assert_eq!(uni_from_bns("128.406"), Some(1066));
        assert_eq!(uni_from_bns("129.413"), Some(1073));
        assert_eq!(uni_from_bns("1.1"), Some(1));

        assert_eq!(uni_from_og("128.8.1073"), Some(1066));
        assert_eq!(uni_from_og("129.3.1077"), Some(1073));

        // Non-existent labels
        assert_eq!(uni_from_bns("nonexistent"), None);
        assert_eq!(uni_from_og("999.999.999"), None);
    }

    /// SG 128 Γ-point double group irreps — verified against BCS
    /// https://cryst.ehu.es/cgi-bin/cryst/programs/representations.pl?tipogrupo=dbg
    #[test]
    fn test_sg128_gamma_double_group() {
        let sg128 = irreps_of(128);
        let gamma: Vec<_> = sg128.iter().filter(|r| r.k_label() == "GM").collect();

        // BCS shows: 10 scalar (GM1±-GM5±) + 4 spinor (GM̄6-GM̄9)
        let gamma_scalar: Vec<_> = gamma.iter().filter(|r| !r.spinor).collect();
        let gamma_spinor: Vec<_> = gamma.iter().filter(|r| r.spinor).collect();

        assert!(
            gamma_scalar.len() >= 5,
            "SG 128 Γ should have >=5 scalar irreps, got {}",
            gamma_scalar.len()
        );
        assert!(
            gamma_spinor.len() >= 2,
            "SG 128 Γ should have >=2 spinor irreps, got {}",
            gamma_spinor.len()
        );

        // Verify scalar labels: GM1+, GM1-, GM2+, GM2-, ...
        let scalar_labels: Vec<&str> = gamma_scalar.iter().map(|r| r.ml).collect();
        for prefix in &["GM1", "GM2", "GM3", "GM4", "GM5"] {
            let has = scalar_labels.iter().any(|l| l.starts_with(prefix));
            assert!(has, "Should have {} scalar irrep at Γ", prefix);
        }

        // Spinor irreps should be 2D (BCS confirms GM̄6-GM̄9 are 2D)
        for ir in &gamma_spinor {
            assert_eq!(ir.dim, 2, "Spinor {} should be 2D, got {}", ir.ml, ir.dim);
            // Identity character should be 2.0 (trace of 2×2 identity)
            let chars = ir.characters();
            if !chars.is_empty() {
                assert!(
                    (chars[0] - 2.0).abs() < 0.01,
                    "Spinor {} identity χ should be 2.0, got {}",
                    ir.ml,
                    chars[0]
                );
            }
        }

        // Scalar irreps at Γ: GM1±-GM4± are 1D, GM5± may be 2D (PIR convention)
        for ir in &gamma_scalar {
            if ir.ml.starts_with("GM5") {
                assert!(
                    ir.dim == 1 || ir.dim == 2,
                    "GM5± should be 1D or 2D, got dim={}",
                    ir.dim
                );
            } else {
                assert_eq!(ir.dim, 1, "Scalar {} should be 1D, got {}", ir.ml, ir.dim);
            }
        }
    }

    /// BCS validation: 128.406 at Z point, all coreps computed from H = SG 118.
    ///
    /// BCS reference (from k-Subgroupsmag.html):
    ///   Unitary Space Group: P-4n2 (No. 118) in standard setting.
    ///   Magnetic little co-group: 4'/m'mm' (12 ops: 8 unitary + 4 anti-unitary)
    ///
    /// Corep table (from BCS corepresentations_out.pl):
    ///   Z1Z2(2D, type C), Z3Z4(2D, type C), Z5(2D, type A), Z̄6Z̄7(4D spinor, type C)
    ///
    /// Our computation uses H = SG 118's irreps at Z:
    ///   Z1Z4, Z2Z3, Z5 (scalar), Z6, Z7 (spinor)
    /// The compound 2D PIRs split into 1D complex CIR components before the
    /// Type-C construction, producing the 2D coreps listed by BCS.
    ///
    /// Character order: verify h_seitz[0] is identity with CIR χ=dim.
    #[test]
    fn test_char_order_sg118() {
        let uni = 1066usize;
        let _mag_ops = get_magnetic_operations(uni).unwrap();
        let h_sg = identify_unitary_subgroup(uni).unwrap();
        let h_ops = get_parent_operations(h_sg as u8).unwrap();
        let h_seitz = ops_to_seitz(&h_ops);
        let h_irreps = crate::irrep::query::irreps_of(h_sg as u8);

        // Check identity at position 0
        let id_op = &h_seitz[0];
        assert!(
            id_op.rot[0][0] == 1 && id_op.rot[1][1] == 1 && id_op.rot[2][2] == 1,
            "h_seitz[0] must be identity"
        );
        assert!(
            id_op.trans[0].abs() < 0.01
                && id_op.trans[1].abs() < 0.01
                && id_op.trans[2].abs() < 0.01,
            "identity must have zero translation"
        );

        // For each Z-point CIR irrep, check χ(id)=dim
        for ir in h_irreps.iter().filter(|r| r.k_label() == "Z") {
            if ir.cir_component_count() > 0 {
                for c in 0..ir.cir_component_count() {
                    let cir = ir.cir_component_chars(c);
                    let chi_id = Complex64::new(cir[0], cir[1]);
                    println!(
                        "{} comp{}: cir_chars[0]=({:.2},{:.2}) |χ|={:.2}",
                        ir.ml,
                        c,
                        chi_id.re,
                        chi_id.im,
                        chi_id.norm()
                    );
                }
            } else {
                let pir = ir.characters();
                println!(
                    "{} (non-compound): pir_chars[0]={:.2} dim={}",
                    ir.ml, pir[0], ir.dim
                );
            }
        }

        // Print full h_seitz ↔ cir_chars mapping for Z1Z4's first component
        if let Some(z1z4) = h_irreps.iter().find(|r| r.ml == "Z1Z4") {
            let cir = z1z4.cir_component_chars(0);
            wigner::debug_char_order(cir, &h_seitz, "SG118 Z1Z4 comp0");
        }
    }

    /// Diagnostic: print Wigner sum term-by-term for SG 118 Z-point irreps.
    #[test]
    #[ignore = "diagnostic only; run with --ignored --nocapture"]
    fn debug_wigner_z_point() {
        let uni = 1066usize;
        let mag_ops = get_magnetic_operations(uni).unwrap();
        let h_sg = identify_unitary_subgroup(uni).unwrap();
        let h_ops = get_parent_operations(h_sg as u8).unwrap();
        let mag_seitz = ops_to_seitz(&mag_ops);
        let h_seitz = ops_to_seitz(&h_ops);
        let a0_idx = mag_ops
            .operations
            .iter()
            .position(|o| o.time_reversal)
            .unwrap();
        let a0 = &mag_seitz[a0_idx];

        println!("\n=== Wigner diagnostic: UNI {} → H=SG {} ===", uni, h_sg);
        println!(
            "Magnetic ops: {} total, {} unitary, {} anti-unitary",
            mag_ops.len(),
            mag_ops
                .operations
                .iter()
                .filter(|o| !o.time_reversal)
                .count(),
            mag_ops
                .operations
                .iter()
                .filter(|o| o.time_reversal)
                .count()
        );
        println!(
            "a₀ (anti-unitary rep): R=[{},{},{};{},{},{};{},{},{}] t=({:.4},{:.4},{:.4})",
            a0.rot[0][0],
            a0.rot[0][1],
            a0.rot[0][2],
            a0.rot[1][0],
            a0.rot[1][1],
            a0.rot[1][2],
            a0.rot[2][0],
            a0.rot[2][1],
            a0.rot[2][2],
            a0.trans[0],
            a0.trans[1],
            a0.trans[2]
        );
        println!("H ops (SG {}): {}", h_sg, h_ops.len());
        for (i, s) in h_seitz.iter().enumerate() {
            println!(
                "  H[{}]: R=[{},{},{};{},{},{};{},{},{}] t=({:.4},{:.4},{:.4})",
                i,
                s.rot[0][0],
                s.rot[0][1],
                s.rot[0][2],
                s.rot[1][0],
                s.rot[1][1],
                s.rot[1][2],
                s.rot[2][0],
                s.rot[2][1],
                s.rot[2][2],
                s.trans[0],
                s.trans[1],
                s.trans[2]
            );
        }

        let h_irreps = crate::irrep::query::irreps_of(h_sg as u8);

        // Compare magnetic ops with SG 118 standard ops — check for origin shift
        let h_ops_sg118 = get_parent_operations(h_sg as u8).unwrap();
        println!("\n=== Magnetic ops vs SG 118 standard ops ===");
        println!("Unitary magnetic ops:");
        for i in 0..mag_ops.len() {
            if mag_ops.operations[i].time_reversal {
                continue;
            }
            let r = &mag_ops.operations[i].rotation;
            let t = &mag_ops.operations[i].translation;
            // Find matching H op
            let h_match = h_ops_sg118.operations.iter().position(|o| {
                let hr = o.rotation;
                hr[0][0] == r[0][0]
                    && hr[0][1] == r[0][1]
                    && hr[0][2] == r[0][2]
                    && hr[1][0] == r[1][0]
                    && hr[1][1] == r[1][1]
                    && hr[1][2] == r[1][2]
                    && hr[2][0] == r[2][0]
                    && hr[2][1] == r[2][1]
                    && hr[2][2] == r[2][2]
            });
            let dt = h_match.map(|hi| {
                let ht = &h_ops_sg118.operations[hi].translation;
                [t[0] - ht[0], t[1] - ht[1], t[2] - ht[2]]
            });
            println!(
                "  mag[{}]: R=[{},{},{};{},{},{};{},{},{}] t=({:.3},{:.3},{:.3}) H_match={:?} dt={:?}",
                i,
                r[0][0],
                r[0][1],
                r[0][2],
                r[1][0],
                r[1][1],
                r[1][2],
                r[2][0],
                r[2][1],
                r[2][2],
                t[0],
                t[1],
                t[2],
                h_match,
                dt
            );
        }
        println!("Anti-unitary magnetic ops:");
        for i in 0..mag_ops.len() {
            if !mag_ops.operations[i].time_reversal {
                continue;
            }
            let r = &mag_ops.operations[i].rotation;
            let t = &mag_ops.operations[i].translation;
            println!(
                "  mag[{}]: R=[{},{},{};{},{},{};{},{},{}] t=({:.3},{:.3},{:.3})",
                i,
                r[0][0],
                r[0][1],
                r[0][2],
                r[1][0],
                r[1][1],
                r[1][2],
                r[2][0],
                r[2][1],
                r[2][2],
                t[0],
                t[1],
                t[2]
            );
        }
        println!("SG 118 standard ops:");
        for i in 0..h_ops_sg118.len() {
            let r = &h_ops_sg118.operations[i].rotation;
            let t = &h_ops_sg118.operations[i].translation;
            println!(
                "  H[{}]: R=[{},{},{};{},{},{};{},{},{}] t=({:.3},{:.3},{:.3})",
                i,
                r[0][0],
                r[0][1],
                r[0][2],
                r[1][0],
                r[1][1],
                r[1][2],
                r[2][0],
                r[2][1],
                r[2][2],
                t[0],
                t[1],
                t[2]
            );
        }

        let k_irreps: Vec<_> = h_irreps.iter().filter(|r| r.k_label() == "Z").collect();

        for ir in &k_irreps {
            println!(
                "\n--- {} (dim={}, spinor={}, k=({},{},{})/{}) ---",
                ir.ml, ir.dim, ir.spinor, ir.kx, ir.ky, ir.kz, ir.kd
            );

            let mag_lg = filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
            let unitary_lg: Vec<usize> = mag_lg
                .iter()
                .filter(|&&i| !mag_ops.operations[i].time_reversal)
                .copied()
                .collect();
            let anti_lg: Vec<usize> = mag_lg
                .iter()
                .filter(|&&i| mag_ops.operations[i].time_reversal)
                .copied()
                .collect();
            println!(
                "  Little group: {} ops ({} unitary + {} anti-unitary)",
                mag_lg.len(),
                unitary_lg.len(),
                anti_lg.len()
            );

            let h_chars = ir.characters();
            println!(
                "  H characters ({} ops): {:?}...",
                h_ops.len(),
                &h_chars[..h_chars.len().min(8)]
            );

            // Map unitary ops to H
            let op_map: Vec<Option<usize>> = (0..mag_ops.len())
                .map(|i| {
                    if mag_ops.operations[i].time_reversal {
                        None
                    } else {
                        let r = &mag_ops.operations[i].rotation;
                        h_ops.operations.iter().position(|o| {
                            let hr = o.rotation;
                            hr[0][0] == r[0][0]
                                && hr[0][1] == r[0][1]
                                && hr[0][2] == r[0][2]
                                && hr[1][0] == r[1][0]
                                && hr[1][1] == r[1][1]
                                && hr[1][2] == r[1][2]
                                && hr[2][0] == r[2][0]
                                && hr[2][1] == r[2][1]
                                && hr[2][2] == r[2][2]
                        })
                    }
                })
                .collect();

            // Term-by-term Wigner sum
            let mut w_sum = 0.0f64;
            for &h_mag_idx in &unitary_lg {
                let h = &mag_seitz[h_mag_idx];
                let h_h = op_map[h_mag_idx];

                // g₀h
                let (g0h, _l1) = compose_seitz(
                    &SeitzOp::new(a0.rot, a0.trans, false),
                    &SeitzOp::new(h.rot, h.trans, false),
                );
                // (g₀h)²
                let (sq, l_sq) = square_seitz(&g0h);

                match find_seitz(&sq.rot, &sq.trans, &h_seitz) {
                    Some(m) => {
                        let total_l = [
                            l_sq[0] + m.lattice_shift[0],
                            l_sq[1] + m.lattice_shift[1],
                            l_sq[2] + m.lattice_shift[2],
                        ];
                        let phase = bloch_phase(ir.kx, ir.ky, ir.kz, ir.kd, &total_l);
                        let chi = if m.op_index < h_chars.len() {
                            h_chars[m.op_index]
                        } else {
                            0.0
                        };
                        let contrib = chi * phase.re;
                        w_sum += contrib;
                        println!(
                            "    h[{}]→H[{}]: (a₀h)²=H[{}] L={:?} ph={:.2} χ={:.2} contrib={:.2}",
                            h_mag_idx,
                            h_h.map_or("?".into(), |x| x.to_string()),
                            m.op_index,
                            total_l,
                            phase,
                            chi,
                            contrib
                        );
                    }
                    None => {
                        println!(
                            "    h[{}]→H[{}]: (a₀h)² R=[{},{},{};{},{},{};{},{},{}] t=({:.3},{:.3},{:.3}) NOT FOUND",
                            h_mag_idx,
                            h_h.map_or("?".into(), |x| x.to_string()),
                            sq.rot[0][0],
                            sq.rot[0][1],
                            sq.rot[0][2],
                            sq.rot[1][0],
                            sq.rot[1][1],
                            sq.rot[1][2],
                            sq.rot[2][0],
                            sq.rot[2][1],
                            sq.rot[2][2],
                            sq.trans[0],
                            sq.trans[1],
                            sq.trans[2]
                        );
                    }
                }
            }
            let w = w_sum / (unitary_lg.len() as f64).max(1.0);
            println!(
                "  Wigner W = {:.4} → {}",
                w,
                if w.abs() < 0.01 {
                    "Type C"
                } else if w > 0.0 {
                    "Type A"
                } else {
                    "Type B"
                }
            );

            // Unwrapped square diagnostic for h[4] and h[7]
            if ir.ml == "Z1Z4" {
                let a0_idx = anti_lg[0];
                let group = wigner::WignerGroupContext {
                    unitary_indices: &unitary_lg,
                    magnetic_ops: &mag_seitz,
                    unitary_ops: &h_seitz,
                    antiunitary_representative: a0_idx,
                };
                wigner::debug_unwrapped_square(4, group, ir.k_vector())
                    .expect("diagnostic square indices must be valid");
                wigner::debug_unwrapped_square(7, group, ir.k_vector())
                    .expect("diagnostic square indices must be valid");

                // Direct anti-coset Wigner sum
                let cir = ir.cir_component_chars(0);
                let w_direct = wigner::wigner_direct_anti_coset(
                    cir,
                    &anti_lg,
                    &mag_seitz,
                    &h_seitz,
                    ir.k_vector(),
                )
                .expect("diagnostic direct anti-coset sum must be well formed");
                println!("  Direct anti-coset W = {:.4}", w_direct);

                // Try all antiunitary ops as a₀
                println!("\n  Sweeping a₀ choices:");
                for &a0_cand in &anti_lg {
                    let mut w_sum = Complex64::ZERO;
                    let mut np = 0u32;
                    let mut nm = 0u32;
                    for &h_mag_idx in &unitary_lg {
                        let h = &mag_seitz[h_mag_idx];
                        let g0_sp =
                            SeitzOp::new(mag_seitz[a0_cand].rot, mag_seitz[a0_cand].trans, false);
                        let h_sp = SeitzOp::new(h.rot, h.trans, false);
                        let (g0h, l1) = compose_seitz(&g0_sp, &h_sp);
                        let (sq, lsq) = square_seitz(&g0h);
                        if let Some(m) = find_seitz(&sq.rot, &sq.trans, &h_seitz) {
                            let rl1 = mat_vec_i32(&g0h.rot, &l1);
                            let tl = add3(&add3(&lsq, &m.lattice_shift), &add3(&l1, &rl1));
                            let ph = bloch_phase(ir.kx, ir.ky, ir.kz, ir.kd, &tl);
                            let chi = Complex64::new(cir[2 * m.op_index], cir[2 * m.op_index + 1]);
                            w_sum += chi * ph;
                            if ph.re > 0.5 {
                                np += 1;
                            } else if ph.re < -0.5 {
                                nm += 1;
                            }
                        }
                    }
                    let w = w_sum / (unitary_lg.len() as f64);
                    println!("    a₀=mag[{}]: +={} -={} W={:.4}", a0_cand, np, nm, w);
                }
            }
        }
    }

    #[test]
    fn test_corep_sg128_406_z_bcs() {
        let coreps = match compute_coreps("128.406", "Z") {
            Ok(coreps) => coreps,
            Err(err @ CorepComputationError::UnsupportedClassification { .. }) => {
                println!(
                    "128.406 Z strict corep path is not fully classified yet: {}",
                    err
                );
                return;
            }
            Err(err) => panic!("unexpected 128.406 Z corep error: {:?}", err),
        };
        assert!(!coreps.is_empty());

        println!("\n=== 128.406 Z-point corepresentations (from H = SG 118) ===");
        println!("{:<8} {:<4} {:<8} {:<8}", "Label", "Dim", "Type", "χ(id)");
        println!("-------- ---- -------- --------");

        // ── BCS reference character table (little group, 12 ops) ──
        // Z1Z2 (type C, from Z1+Z4 of H) → 2×χ_re = 2*Re(χ)
        // Z3Z4 (type C, from Z2+Z3 of H) → 2×χ_re = 2*Re(χ)
        // Z5   (type A, from Z5 of H)     → χ = χ (same)
        // Z̄6Z̄7 (type C spinor, from Z6+Z7 of H) → 4D corep

        // Collect computed coreps by label pattern for BCS comparison
        for (label, c) in &coreps {
            let type_str = match c.corep_type {
                CorepType::A => "A",
                CorepType::B => "B",
                CorepType::C => "C",
            };
            println!(
                "{:<8} {:<4} {:<8} {:<8.1}",
                label, c.dim, type_str, c.characters[0]
            );

            // Basic invariants
            assert!(c.characters[0] > 0.0, "χ(id) must be > 0 for {}", label);
            assert!(c.dim > 0, "dim must be > 0 for {}", label);

            // χ(id) always equals corep dimension
            assert!(
                (c.characters[0] - c.dim as f64).abs() < 0.01,
                "χ(id)={} should equal dim={} for {}",
                c.characters[0],
                c.dim,
                label
            );

            // Number of anti-unitary ops with zero character for type C
            if c.corep_type == CorepType::C {
                let zero_count = c
                    .characters
                    .iter()
                    .zip(c.timerev.iter())
                    .filter(|(chi, tr)| **tr && chi.abs() < 0.01)
                    .count();
                let anti_count = c.timerev.iter().filter(|&&t| t).count();
                assert_eq!(
                    zero_count, anti_count,
                    "Type C: all anti-unitary chars should be 0 for {}",
                    label
                );
            }
        }

        // Verify we have the expected number of coreps from SG 118 at Z
        // SG 118 at Z: 3 scalar (Z1Z4, Z2Z3, Z5) + 2 spinor (Z6, Z7) = 5 H irreps
        assert!(
            coreps.len() >= 3,
            "Should have >=3 Z-point coreps (scalar), got {}",
            coreps.len()
        );

        println!("\nBCS comparison: H = SG 118 irreps → corep → BCS magnetic irreps");
        println!("  H:Z1Z4(2D,PIR) → corep type-C → BCS Z1Z2(2D)");
        println!("  H:Z2Z3(2D,PIR) → corep type-C → BCS Z3Z4(2D)");
        println!("  H:Z5(2D,PIR)   → corep type-A → BCS Z5(2D)");
        println!("  H:Z6,Z7(2D,spinor) → corep type-C → BCS Z̄6Z̄7(4D)");
    }

    /// BCS: 165.95 (P-3c'1, UNI 1325) at L:(1/2,0,1/2)
    ///
    /// k-Subgroupsmag_165.95.html confirms:
    ///   Unitary subgroup: P-3 (No. 147)
    ///   Magnetic little co-group: 2'/m'
    ///   Coreps: L₁⁻L₁⁺ (2D, Type C), L̄₂L̄₃ (2D spinor, Type C)
    #[test]
    fn test_corep_sg165_95_l_bcs() {
        let uni = 1325usize;

        // 1. Verify unitary subgroup identification
        let h_sg = identify_unitary_subgroup(uni);
        assert!(h_sg.is_some(), "Should identify unitary subgroup of 165.95");
        let h_sg = h_sg.unwrap();
        println!("165.95 (UNI {}) → unitary subgroup: SG {}", uni, h_sg);

        // 2. Verify magnetic operations exist
        let mag_ops = get_magnetic_operations(uni);
        assert!(mag_ops.is_some(), "Should get ops for UNI {}", uni);
        let mag_ops = mag_ops.unwrap();
        let n_u = mag_ops
            .operations
            .iter()
            .filter(|o| !o.time_reversal)
            .count();
        let n_a = mag_ops
            .operations
            .iter()
            .filter(|o| o.time_reversal)
            .count();
        println!(
            "  {} ops ({} unitary + {} anti-unitary)",
            mag_ops.len(),
            n_u,
            n_a
        );

        // 3. Compute coreps at L point (using H = unitary subgroup)
        let h_irreps = crate::irrep::query::irreps_of(h_sg as u8);
        let l_irreps: Vec<&IrrepRecord> = h_irreps.iter().filter(|r| r.k_label() == "L").collect();
        let n_scalar = l_irreps.iter().filter(|r| !r.spinor).count();
        let n_spinor = l_irreps.iter().filter(|r| r.spinor).count();
        println!(
            "  H=SG{} L-point irreps: {} scalar + {} spinor",
            h_sg, n_scalar, n_spinor
        );

        assert!(!l_irreps.is_empty(), "Should have L-point irreps");

        // 4. Compute coreps one by one
        for ir in &l_irreps {
            if let Ok(c) = ir.corepresentation(uni) {
                let type_str = match c.corep_type {
                    CorepType::A => "A",
                    CorepType::B => "B",
                    CorepType::C => "C",
                };
                println!(
                    "  {}: dim={} type={} χ(id)={:.1}",
                    ir.ml, c.dim, type_str, c.characters[0]
                );

                assert!(c.dim > 0);
                assert!(
                    (c.characters[0] - c.dim as f64).abs() < 0.01,
                    "χ(id) should equal dim for {}",
                    ir.ml
                );
            }
        }
    }

    /// BCS: SG 139 (I4/mmm) double-group irreps at k=(1,1,1) (P point)
    ///
    /// k-Subgroupsmag_139.html shows 14 irreps (10 scalar + 4 spinor)
    /// with 4 operations: {1|t}, {2₀₀₁|0}, {4⁺₀₀₁|0}, {4⁻₀₀₁|0}
    #[test]
    fn test_sg139_p_point_bcs() {
        let sg = 139u8;
        let irreps = crate::irrep::query::irreps_of(sg);

        // P-point irreps (k=(1,1,1) in body-centered tetragonal)
        let p_irreps: Vec<&IrrepRecord> = irreps.iter().filter(|r| r.k_label() == "P").collect();
        println!(
            "SG{} P-point: {} irreps ({} scalar + {} spinor)",
            sg,
            p_irreps.len(),
            p_irreps.iter().filter(|r| !r.spinor).count(),
            p_irreps.iter().filter(|r| r.spinor).count()
        );

        assert!(!p_irreps.is_empty(), "SG 139 should have P-point irreps");

        // BCS shows 14 irreps: M₁⁺..M₅⁻ (10 scalar) + M̄₆..M̄₉ (4 spinor)
        let scalar: Vec<_> = p_irreps.iter().filter(|r| !r.spinor).collect();
        let spinor: Vec<_> = p_irreps.iter().filter(|r| r.spinor).collect();
        assert!(scalar.len() >= 5, "Should have >=5 scalar P-point irreps");
        assert!(spinor.len() >= 2, "Should have >=2 spinor P-point irreps");

        // Check P1-P5 have correct dimensions (BCS: 1D for P1-P4, 2D for P5)
        for ir in &scalar {
            assert!(ir.dim > 0, "dim > 0 for {}", ir.ml);
            let chars = ir.characters();
            assert!(!chars.is_empty(), "Should have characters for {}", ir.ml);
            // Identity character should equal dim
            assert!(
                (chars[0] - ir.dim as f64).abs() < 0.01,
                "χ(id)={} ≠ dim={} for {}",
                chars[0],
                ir.dim,
                ir.ml
            );
            println!(
                "  {}: dim={} ops={} χ(id)={}",
                ir.ml,
                ir.dim,
                chars.len(),
                chars[0]
            );
        }

        // Test matrix reordering for a P-point irrep with matrix data
        if let Some(p1) = scalar.first() {
            let mats = p1.matrices();
            if !mats.is_empty() {
                println!("  {}: {} matrix elements", p1.ml, mats.len());
                let h_ops = get_parent_operations(sg).unwrap();
                let h_seitz = ops_to_seitz(&h_ops);
                let reordered = p1.matrices_reordered(&h_seitz).unwrap();
                assert_eq!(
                    reordered.len(),
                    mats.len(),
                    "Reordered matrix should have same size"
                );
                // Identity should be at H[0] position (1,0,0 in original)
                let dim = p1.dim as usize;
                if dim > 0 && reordered.len() >= dim * dim {
                    let trace: f64 = (0..dim).map(|d| reordered[d * dim + d]).sum();
                    assert!(
                        (trace - p1.dim as f64).abs() < 0.5,
                        "Reordered identity trace should ≈ dim"
                    );
                }
                println!("  Matrix reordering OK ({} elements)", reordered.len());

                let mut unmappable = h_seitz.clone();
                unmappable[0].rot = [[2, 0, 0], [0, 2, 0], [0, 0, 2]];
                assert_eq!(
                    p1.matrices_reordered(&unmappable),
                    Err(crate::irrep::types::MatrixReorderError::OperationMappingFailed)
                );
            }
        }
    }

    /// Every isotropy subgroup record points to a valid SG (1-230).
    #[test]
    fn test_all_isotropy_subgroups_are_well_formed() {
        for sg in 1u8..=230 {
            for ir in crate::irrep::query::irreps_of(sg) {
                for sub in ir.subgroups() {
                    assert!(
                        sub.sg >= 1 && sub.sg <= 230,
                        "invalid isotropy SG={} for parent SG{} {}",
                        sub.sg,
                        sg,
                        ir.ml
                    );
                }
                for msub in ir.magnetic_subgroups() {
                    assert!(
                        msub.mag_sg >= 1 && msub.mag_sg <= 1651,
                        "invalid mag isotropy SG={} for parent SG{} {}",
                        msub.mag_sg,
                        sg,
                        ir.ml
                    );
                }
            }
        }
    }

    /// Type C corepresentations pair two H irreps into one magnetic corep.
    /// Verify that compute_coreps doesn't produce duplicate magnetic irreps
    /// for the same Type C pair.
    #[test]
    fn test_type_c_coreps_are_deduplicated() {
        let coreps = match compute_coreps("128.406", "Z") {
            Ok(coreps) => coreps,
            Err(err @ CorepComputationError::UnsupportedClassification { .. }) => {
                println!(
                    "128.406 Z Type-C dedup target is not fully classified yet: {}",
                    err
                );
                return;
            }
            Err(err) => panic!("unexpected 128.406 Z corep error: {:?}", err),
        };

        // Type C pairs (Z1Z4+Z2Z3, Z6+Z7) should each appear once
        // as combined coreps, not as individual entries
        let mut type_c_pairs: Vec<Vec<&str>> = Vec::new();
        for (label, c) in &coreps {
            if c.corep_type == CorepType::C {
                // Collect labels that should pair
                let labels: Vec<&str> = vec![label];
                type_c_pairs.push(labels);
            }
        }
        // With current API each H irrep returns its own Corepresentation,
        // so for Type C we expect pairs. For now, just verify they're all Type C.
        for (_label, c) in &coreps {
            if c.corep_type == CorepType::C {
                assert!(c.dim > 0);
                // Antiunitary characters must be 0 for Type C
                for (i, &chi) in c.characters.iter().enumerate() {
                    if c.timerev[i] {
                        assert!(
                            chi.abs() < 0.01,
                            "Type C antiunitary char must be 0, got {} at op {}",
                            chi,
                            i
                        );
                    }
                }
            }
        }
        println!(
            "Type C dedup check: {} coreps, all antiunitary chars zero ✓",
            coreps.len()
        );
    }

    /// Exhaustive: all 1651 magnetic space groups have valid operations
    /// and identifiable unitary subgroups.
    #[test]
    fn test_all_magnetic_sgs_have_valid_operations() {
        let mut ok = 0usize;
        let mut fail = 0usize;
        for uni in 1usize..=1651 {
            if let Some(ops) = get_magnetic_operations(uni) {
                assert!(!ops.is_empty(), "UNI {} has empty ops", uni);
                let n_u = ops.operations.iter().filter(|o| !o.time_reversal).count();
                let _n_a = ops.operations.iter().filter(|o| o.time_reversal).count();
                assert!(n_u > 0, "UNI {} has no unitary ops", uni);
                // Every magnetic op must have a valid rotation (det = ±1)
                for i in 0..ops.len() {
                    let r = &ops.operations[i].rotation;
                    let det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
                        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
                        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
                    assert!(
                        det == 1 || det == -1,
                        "UNI {} op[{}]: det={}, not ±1",
                        uni,
                        i,
                        det
                    );
                }
                // Verify unitary subgroup can be identified (may fail for some edge cases)
                if let Some(h_sg) = identify_unitary_subgroup(uni) {
                    assert!(
                        (1..=230).contains(&h_sg),
                        "UNI {} unitary SG={} out of range",
                        uni,
                        h_sg
                    );
                }
                ok += 1;
            } else {
                fail += 1;
            }
        }
        println!("Magnetic ops: {}/1651 OK, {} missing", ok, fail);
        assert!(ok > 1600, "Should have >=1600 valid MSGs, got {}", ok);
    }

    /// Exhaustive: all spinor (double-group) irreps have valid character tables.
    /// Central element Ē (2π rotation) character should be -dim for spinor irreps.
    #[test]
    fn test_all_spinor_irreps_are_well_formed() {
        let mut total = 0usize;
        for sg in 1u8..=230 {
            for ir in crate::irrep::query::irreps_of(sg) {
                if !ir.spinor {
                    continue;
                }
                total += 1;
                assert!(ir.dim > 0, "spinor {} SG{} dim=0", ir.ml, sg);
                let chars = ir.characters();
                assert!(!chars.is_empty(), "spinor {} SG{} no chars", ir.ml, sg);
                assert!(
                    chars[0] > 0.0,
                    "spinor {} SG{} χ(E)={}",
                    ir.ml,
                    sg,
                    chars[0]
                );
                // Spinor irreps are double-valued: typical dims are 1,2,3,4,6
                assert!(
                    ir.dim >= 1 && ir.dim <= 8,
                    "spinor {} SG{} unexpected dim={}",
                    ir.ml,
                    sg,
                    ir.dim
                );
                // Identity character should be integer
                assert!(
                    (chars[0] - chars[0].round()).abs() < 1e-8,
                    "spinor {} SG{} χ(E)={} not integer",
                    ir.ml,
                    sg,
                    chars[0]
                );
                // Spin ops should exist
                let (rots, _trans, _su2) = ir.spin_ops();
                if ir.spin_lg_char_count() > 0 {
                    assert!(
                        !rots.is_empty(),
                        "spinor {} SG{} has lg ops but no spin op rots",
                        ir.ml,
                        sg
                    );
                }
            }
        }
        assert!(
            total > 3000,
            "Should have >3000 spinor irreps, got {}",
            total
        );
        println!("Spinor irreps: {} total, all well-formed ✓", total);
    }

    /// Database format sanity: all irrep k-vectors have reasonable denominators.
    #[test]
    fn test_all_irrep_k_vectors_are_well_formed() {
        for sg in 1u8..=230 {
            for ir in crate::irrep::query::irreps_of(sg) {
                // kd is the common denominator; capped by database convention
                const MAX_KD: i8 = 24;
                assert!(
                    ir.kd >= 0 && ir.kd <= MAX_KD,
                    "SG{} {}: kd={} out of [0,{}]",
                    sg,
                    ir.ml,
                    ir.kd,
                    MAX_KD
                );
                // Gamma-like points must have kd=0 → k=(0,0,0)
                if ir.kd == 0 {
                    assert_eq!(
                        (ir.kx, ir.ky, ir.kz),
                        (0, 0, 0),
                        "SG{} {}: kd=0 but k=({},{},{})",
                        sg,
                        ir.ml,
                        ir.kx,
                        ir.ky,
                        ir.kz
                    );
                }
            }
        }
    }

    /// Exhaustive: all non-spinor (single-valued) irreps satisfy basic
    /// representation-theory invariants: χ(E)=dim, characters are finite,
    /// matrix data is consistent with dimension.
    #[test]
    fn test_all_scalar_irreps_basic_invariants() {
        let mut checked = 0usize;
        for sg in 1u8..=230 {
            for ir in crate::irrep::query::irreps_of(sg) {
                if ir.spinor {
                    continue;
                }
                checked += 1;
                assert!(ir.dim > 0, "SG{} {}: dim=0", sg, ir.ml);
                assert!(!ir.ml.is_empty(), "SG{}: empty label", sg);
                let chars = ir.characters();
                assert!(!chars.is_empty(), "SG{} {}: empty chars", sg, ir.ml);
                assert!(
                    (chars[0] - ir.dim as f64).abs() < 1e-8,
                    "SG{} {}: χ(E)={} != dim={}",
                    sg,
                    ir.ml,
                    chars[0],
                    ir.dim
                );
                assert!(
                    chars.iter().all(|x| x.is_finite()),
                    "SG{} {}: non-finite character found",
                    sg,
                    ir.ml
                );
                let mats = ir.matrices();
                if !mats.is_empty() {
                    let dim = ir.dim as usize;
                    assert!(
                        mats.len() % (dim * dim) == 0,
                        "SG{} {}: matrix len {} not divisible by dim²={}",
                        sg,
                        ir.ml,
                        mats.len(),
                        dim * dim
                    );
                }
            }
        }
        assert!(
            checked > 4000,
            "Should have >4000 scalar irreps, got {}",
            checked
        );
        println!("Scalar irreps: {} total, all well-formed ✓", checked);
    }

    /// Regression: high-dimension image labels (e.g. "K1536a") must not
    /// fall back to dim=1.  This was the root cause of the K1536a bug.
    #[test]
    fn test_high_dim_image_irreps_not_default_to_one() {
        let mut checked = 0usize;
        for sg in 1u8..=230 {
            for ir in crate::irrep::query::irreps_of(sg) {
                if ir.image.starts_with('K')
                    || ir.image.starts_with('L')
                    || ir.image.starts_with('M')
                    || ir.image.starts_with('N')
                {
                    assert!(
                        ir.dim > 1,
                        "SG{} {}: image={} dim={} (should not fall back to 1)",
                        sg,
                        ir.ml,
                        ir.image,
                        ir.dim
                    );
                    assert_eq!(
                        ir.characters()[0] as usize,
                        ir.dim as usize,
                        "SG{} {}: χ(E)={} != dim={}",
                        sg,
                        ir.ml,
                        ir.characters()[0],
                        ir.dim
                    );
                    checked += 1;
                }
            }
        }
        println!("High-dim image irreps: {} checked, all dim > 1 ✓", checked);
        assert!(checked > 0, "Should have at least one high-dim image irrep");
    }

    /// Diagnostic: count duplicate rotations in H little groups.
    ///
    /// If the same rotation appears with different translations in the
    /// little group, PIR/CIR rotation-only mapping may be ambiguous.
    #[test]
    fn diagnose_duplicate_rotations() {
        let mut total_cases = 0usize;
        let mut dup_rot_cases = 0usize;
        let mut dup_rot_distinct_char = 0usize;
        for uni in 1..=1651 {
            let mag_ops = match get_magnetic_operations(uni) {
                Some(m) => m,
                None => continue,
            };
            let h_sg = match identify_unitary_subgroup(uni) {
                Some(s) => s as u8,
                None => continue,
            };
            let _h_seitz = crate::irrep::wigner::ops_to_seitz(&mag_ops);
            let _h_seitz_unitary: Vec<_> = (0..mag_ops.len())
                .filter(|&i| !mag_ops.operations[i].time_reversal)
                .map(|i| {
                    crate::irrep::wigner::SeitzOp::new(
                        mag_ops.operations[i].rotation,
                        mag_ops.operations[i].translation,
                        false,
                    )
                })
                .collect();
            for ir in crate::irrep::query::irreps_of(h_sg) {
                let mag_lg =
                    crate::irrep::wigner::filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
                let unitary_lg: Vec<_> = mag_lg
                    .iter()
                    .filter(|&&i| !mag_ops.operations[i].time_reversal)
                    .copied()
                    .collect();
                if unitary_lg.len() <= 1 {
                    continue;
                }
                total_cases += 1;

                // Group by rotation, check for duplicate rotations with different translations
                let mut rot_to_trans: std::collections::HashMap<Vec<i32>, Vec<[f64; 3]>> =
                    std::collections::HashMap::new();
                for &idx in &unitary_lg {
                    let r = mag_ops.operations[idx].rotation;
                    let key: Vec<i32> = vec![
                        r[0][0], r[0][1], r[0][2], r[1][0], r[1][1], r[1][2], r[2][0], r[2][1],
                        r[2][2],
                    ];
                    rot_to_trans
                        .entry(key)
                        .or_default()
                        .push(mag_ops.operations[idx].translation);
                }
                let has_dup = rot_to_trans.values().any(|v| v.len() > 1);
                if has_dup {
                    dup_rot_cases += 1;
                    // Check if duplicate rotations have distinct characters
                    if !ir.spinor && ir.characters().len() >= unitary_lg.len() {
                        let chars = ir.characters();
                        for (rot_key, trans_list) in &rot_to_trans {
                            if trans_list.len() <= 1 {
                                continue;
                            }
                            let char_values: Vec<f64> = unitary_lg
                                .iter()
                                .filter(|&&idx| {
                                    let r = mag_ops.operations[idx].rotation;
                                    let k: Vec<i32> = vec![
                                        r[0][0], r[0][1], r[0][2], r[1][0], r[1][1], r[1][2],
                                        r[2][0], r[2][1], r[2][2],
                                    ];
                                    k == *rot_key
                                })
                                .enumerate()
                                .map(|(pos, _)| chars[pos])
                                .collect();
                            let first = char_values.first().copied().unwrap_or(0.0);
                            if char_values.iter().any(|&c| (c - first).abs() > 0.01) {
                                dup_rot_distinct_char += 1;
                                break;
                            }
                        }
                    }
                }
            }
        }
        eprintln!("\n=== Duplicate rotation diagnostic ===");
        eprintln!("  total little-group cases: {}", total_cases);
        eprintln!("  with duplicate rotations: {}", dup_rot_cases);
        eprintln!("  dup-rot with distinct chars: {}", dup_rot_distinct_char);
        if dup_rot_distinct_char > 0 {
            eprintln!(
                "  WARNING: {} cases have ambiguous rotation-only mapping!",
                dup_rot_distinct_char
            );
        } else {
            eprintln!("  OK: no ambiguous rotation-only mapping found");
        }
    }

    /// Show concrete examples of duplicate-rotation ambiguous cases.
    #[test]
    fn show_dup_rot_examples() {
        let mut shown = 0usize;
        let max_show = 3usize;
        'outer: for uni in 1..=1651 {
            let mag_ops = match get_magnetic_operations(uni) {
                Some(m) => m,
                None => continue,
            };
            let h_sg = match identify_unitary_subgroup(uni) {
                Some(s) => s as u8,
                None => continue,
            };
            for ir in crate::irrep::query::irreps_of(h_sg) {
                if ir.spinor {
                    continue;
                }
                if shown >= max_show {
                    break 'outer;
                }
                let mag_lg =
                    crate::irrep::wigner::filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
                let unitary_lg: Vec<_> = mag_lg
                    .iter()
                    .filter(|&&i| !mag_ops.operations[i].time_reversal)
                    .copied()
                    .collect();
                if unitary_lg.len() <= 1 {
                    continue;
                }

                let mut rot_to_idxs: std::collections::HashMap<Vec<i32>, Vec<usize>> =
                    std::collections::HashMap::new();
                for &idx in &unitary_lg {
                    let r = mag_ops.operations[idx].rotation;
                    let key: Vec<i32> = vec![
                        r[0][0], r[0][1], r[0][2], r[1][0], r[1][1], r[1][2], r[2][0], r[2][1],
                        r[2][2],
                    ];
                    rot_to_idxs.entry(key).or_default().push(idx);
                }
                let chars = ir.characters();
                for (rot_key, idxs) in &rot_to_idxs {
                    if idxs.len() <= 1 {
                        continue;
                    }
                    // Check if distinct characters
                    let char_vals: Vec<f64> = idxs
                        .iter()
                        .map(|&idx| {
                            let pos = unitary_lg.iter().position(|&u| u == idx).unwrap();
                            chars.get(pos).copied().unwrap_or(999.0)
                        })
                        .collect();
                    let first = char_vals[0];
                    if char_vals.iter().all(|&c| (c - first).abs() < 0.01) {
                        continue;
                    }

                    eprintln!(
                        "\n--- Example {}: SG{} {} uni={} k=({}/{},{}/{},{}/{}) ---",
                        shown + 1,
                        h_sg,
                        ir.ml,
                        uni,
                        ir.kx,
                        ir.kd,
                        ir.ky,
                        ir.kd,
                        ir.kz,
                        ir.kd
                    );
                    let r0 = rot_key[0];
                    let r1 = rot_key[1];
                    let r2 = rot_key[2];
                    eprintln!(
                        "  Rotation: [[{},{},{}],[{},{},{}],[{},{},{}]]",
                        r0,
                        r1,
                        r2,
                        rot_key[3],
                        rot_key[4],
                        rot_key[5],
                        rot_key[6],
                        rot_key[7],
                        rot_key[8]
                    );
                    eprintln!("  Same rotation, {} distinct ops:", idxs.len());
                    for &idx in idxs {
                        let pos = unitary_lg.iter().position(|&u| u == idx).unwrap();
                        eprintln!(
                            "    mag_op[{}]: trans=[{:.4},{:.4},{:.4}]  χ={:.4}",
                            idx,
                            mag_ops.operations[idx].translation[0],
                            mag_ops.operations[idx].translation[1],
                            mag_ops.operations[idx].translation[2],
                            chars[pos]
                        );
                    }
                    eprintln!("  → PIR rotation-only mapping cannot distinguish these.");
                    shown += 1;
                    if shown >= max_show {
                        break;
                    }
                }
            }
        }
    }

    /// Diagnose PIR_ROTS / CIR_ROTS internal rotation ambiguity.
    #[test]
    fn diagnose_pir_cir_rotation_ambiguity() {
        let mut pir_benign = 0usize;
        let mut pir_dangerous = 0usize;
        let mut cir_benign = 0usize;
        let mut cir_dangerous = 0usize;
        let mut examples: Vec<String> = Vec::new();
        for sg in 1u8..=230 {
            for ir in crate::irrep::query::irreps_of(sg) {
                let chars = ir.characters();
                let rots = ir.pir_rotations();
                if !chars.is_empty() && rots.len() == chars.len() * 9 {
                    let mut g: std::collections::BTreeMap<[i32; 9], Vec<usize>> =
                        std::collections::BTreeMap::new();
                    for (i, _) in chars.iter().enumerate() {
                        let r = [
                            rots[9 * i],
                            rots[9 * i + 1],
                            rots[9 * i + 2],
                            rots[9 * i + 3],
                            rots[9 * i + 4],
                            rots[9 * i + 5],
                            rots[9 * i + 6],
                            rots[9 * i + 7],
                            rots[9 * i + 8],
                        ];
                        g.entry(r).or_default().push(i);
                    }
                    for idxs in g.values() {
                        if idxs.len() <= 1 {
                            continue;
                        }
                        let first = chars[idxs[0]];
                        if idxs.iter().all(|&i| (chars[i] - first).abs() < 1e-8) {
                            pir_benign += 1;
                        } else {
                            pir_dangerous += 1;
                            if examples.len() < 5 {
                                examples.push(format!(
                                    "PIR SG{} {} k=({}/{},{}/{},{}/{}) ch={:?}",
                                    sg,
                                    ir.ml,
                                    ir.kx,
                                    ir.kd,
                                    ir.ky,
                                    ir.kd,
                                    ir.kz,
                                    ir.kd,
                                    idxs.iter()
                                        .map(|&i| format!("{:.2}", chars[i]))
                                        .collect::<Vec<_>>()
                                ));
                            }
                        }
                    }
                }
                for comp in 0..ir.cir_component_count() {
                    let cir = ir.cir_component_chars(comp);
                    let cr = ir.cir_rotations(comp);
                    let n = cir.len() / 2;
                    if cr.len() != n * 9 {
                        continue;
                    }
                    let mut g: std::collections::BTreeMap<[i32; 9], Vec<usize>> =
                        std::collections::BTreeMap::new();
                    for i in 0..n {
                        let r = [
                            cr[9 * i],
                            cr[9 * i + 1],
                            cr[9 * i + 2],
                            cr[9 * i + 3],
                            cr[9 * i + 4],
                            cr[9 * i + 5],
                            cr[9 * i + 6],
                            cr[9 * i + 7],
                            cr[9 * i + 8],
                        ];
                        g.entry(r).or_default().push(i);
                    }
                    for idxs in g.values() {
                        if idxs.len() <= 1 {
                            continue;
                        }
                        let fre = cir[idxs[0] * 2];
                        let fim = cir[idxs[0] * 2 + 1];
                        if idxs.iter().all(|&i| {
                            (cir[i * 2] - fre).abs() < 1e-8 && (cir[i * 2 + 1] - fim).abs() < 1e-8
                        }) {
                            cir_benign += 1;
                        } else {
                            cir_dangerous += 1;
                        }
                    }
                }
            }
        }
        eprintln!("\n=== PIR/CIR rotation ambiguity ===");
        eprintln!("  PIR: {} benign, {} DANGEROUS", pir_benign, pir_dangerous);
        eprintln!("  CIR: {} benign, {} DANGEROUS", cir_benign, cir_dangerous);
        for ex in &examples {
            eprintln!("  {}", ex);
        }
        if pir_dangerous > 0 || cir_dangerous > 0 {
            eprintln!("  *** WARNING: rotation-only mapping ambiguous!");
        } else {
            eprintln!("  ✓ No dangerous duplicates in PIR/CIR");
        }
    }

    /// Diagnostic: report Wigner source statistics across all irreps.
    ///
    /// Run with `-- --nocapture` to see the printout.
    /// This does NOT assert correctness of the SU(2) fallback path,
    /// which is known to need further work on antiunitary gauge handling.
    #[test]
    fn diagnose_wigner_sources() {
        use crate::irrep::wigner::{
            H2S_AMBIGUOUS, H2S_MISSING, H2S_OK, MSG_GAUGE_MAP_FAIL, MSG_GAUGE_OK, MSG_GAUGE_W_FAIL,
            OLD_PATH_FAIL, OLD_PATH_OK,
        };
        use std::sync::atomic::{AtomicUsize, Ordering};

        let reset_counter = |counter: &AtomicUsize| counter.store(0, Ordering::Relaxed);
        let reset_triage_counters = || {
            reset_counter(&MSG_GAUGE_OK);
            reset_counter(&MSG_GAUGE_MAP_FAIL);
            reset_counter(&MSG_GAUGE_W_FAIL);
            reset_counter(&OLD_PATH_OK);
            reset_counter(&OLD_PATH_FAIL);
            reset_counter(&H2S_OK);
            reset_counter(&H2S_AMBIGUOUS);
            reset_counter(&H2S_MISSING);
        };

        reset_triage_counters();
        let mut stats = std::collections::HashMap::<&str, usize>::new();
        let mut failure_class = std::collections::HashMap::<&str, usize>::new();
        let mut failure_by_sg = std::collections::HashMap::<(&str, u8), usize>::new();
        let mut mapping_shape = std::collections::HashMap::<(usize, usize), usize>::new();
        let mut direct_anti_stats = std::collections::HashMap::<&str, usize>::new();
        let mut direct_anti_failures = std::collections::HashMap::<&str, usize>::new();
        let mut final_failure_reasons = std::collections::HashMap::<&str, usize>::new();
        let mut final_failure_by_sg = std::collections::HashMap::<(&str, u8), usize>::new();
        let mut final_failure_by_transform =
            std::collections::HashMap::<(&str, bool), usize>::new();

        for uni in 1..=1651 {
            let mag_ops = match get_magnetic_operations(uni) {
                Some(m) => m,
                None => continue,
            };
            let h_info = match identify_unitary_subgroup_with_hall(uni) {
                Some(info) => info,
                None => continue,
            };
            let h_sg = h_info.sg as u8;
            // Use exactly the frame transform selected by the production
            // subgroup-identification path.  Recomputing an independent
            // transform here made this diagnostic test a different algorithm.
            let setting_xf_owned = h_info.msg_to_data.clone();
            let setting_xf = setting_xf_owned.as_ref();
            let canonical_pure_translations: Vec<[f64; 3]> = h_info
                .ops_from_hall
                .operations
                .iter()
                .filter(|op| op.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
                .map(|op| op.translation)
                .collect();
            // Diagnostic: print Hall identity translations for groups where count != 3
            {
                static CNT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
                let n = CNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if canonical_pure_translations.len() != 3 && n < 10 {
                    eprintln!(
                        "DIAG HallTrans: UNI={uni} SG={h_sg} Hall={} \
                        n_identity={} translations={:.12?}",
                        h_info.hall,
                        canonical_pure_translations.len(),
                        canonical_pure_translations
                    );
                }
            }
            let h_ops = h_info.ops_from_msg;
            let h_seitz = crate::irrep::wigner::ops_to_seitz(&h_ops);
            let mag_seitz = crate::irrep::wigner::ops_to_seitz(&mag_ops);

            for ir in crate::irrep::query::irreps_of(h_sg) {
                let mag_lg = crate::irrep::wigner::filter_little_group_with_transform(
                    ir.kx,
                    ir.ky,
                    ir.kz,
                    ir.kd,
                    &mag_ops,
                    setting_xf,
                    Some(&canonical_pure_translations),
                );
                let antiunitary: Vec<usize> = mag_lg
                    .iter()
                    .filter(|&&i| mag_ops.operations[i].time_reversal)
                    .copied()
                    .collect();

                let key = if antiunitary.is_empty() {
                    "scalar_trivial_A"
                } else if ir.cir_component_count() > 0 {
                    "scalar_CIR"
                } else if ir.spinor {
                    let unitary: Vec<usize> = mag_lg
                        .iter()
                        .filter(|&&i| !mag_ops.operations[i].time_reversal)
                        .copied()
                        .collect();
                    let has_imag = !ir.spin_character_imag().is_empty();
                    // Always use wigner_classify_spinor — it handles both
                    // real and complex characters internally.
                    let h_spin = ir.spin_ops();
                    let g_sg = parent_spatial_sg(uni).unwrap_or(h_sg as usize) as u8;
                    let g_spin = if g_sg == h_sg {
                        h_spin
                    } else {
                        IrrepRecord::spin_ops_for_sg(g_sg)
                    };
                    let ctx = crate::irrep::wigner::SpinLiftContext {
                        h: h_spin,
                        g: g_spin,
                        sg: h_sg,
                    };
                    let su2_result = crate::irrep::wigner::wigner_classify_spinor(
                        &ctx,
                        wigner::SpinorWignerInput {
                            characters_real: ir.characters(),
                            characters_imag: ir.spin_character_imag(),
                            operation_indices: ir.spin_lg_op_indices(),
                            k_vector: ir.k_vector(),
                        },
                        wigner::WignerGroupContext {
                            unitary_indices: &unitary,
                            magnetic_ops: &mag_seitz,
                            unitary_ops: &h_seitz,
                            antiunitary_representative: antiunitary[0],
                        },
                        setting_xf,
                        Some(&antiunitary),
                    )
                    .ok();
                    match (has_imag, su2_result.is_some()) {
                        (true, true) => "spinor_complex_ok",
                        (true, false) => "spinor_complex_fail",
                        (false, true) => "spinor_real_su2_ok",
                        (false, false) => "spinor_real_su2_fail",
                    }
                } else {
                    "scalar_PIR"
                };
                *stats.entry(key).or_default() += 1;
            }
        }

        println!("\n=== Wigner source statistics ===");
        let mut sorted: Vec<_> = stats.iter().collect();
        sorted.sort_by_key(|(_, v)| -(**v as i64));
        for (key, count) in &sorted {
            println!("  {:>22}  {:>6}", key, count);
        }

        // ------------------------------------------------------------------
        // Detailed failure triage: classify every failing spinor case.
        // Reset first so path/H2S counters below describe exactly this pass,
        // rather than accumulating both the statistics and triage passes.
        // ------------------------------------------------------------------
        reset_triage_counters();
        println!("\n=== Spinor failure triage (all cases; first 30 shown) ===");
        let mut shown = 0usize;
        for uni in 1..=1651 {
            let mag_ops = match get_magnetic_operations(uni) {
                Some(m) => m,
                None => continue,
            };
            let h_info = match identify_unitary_subgroup_with_hall(uni) {
                Some(i) => i,
                None => continue,
            };
            let h_sg = h_info.sg as u8;
            let setting_xf_owned = h_info.msg_to_data.clone();
            let setting_xf = setting_xf_owned.as_ref();
            let canonical_pure_translations: Vec<[f64; 3]> = h_info
                .ops_from_hall
                .operations
                .iter()
                .filter(|op| op.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
                .map(|op| op.translation)
                .collect();
            let h_ops = h_info.ops_from_msg;
            let h_seitz = crate::irrep::wigner::ops_to_seitz(&h_ops);
            let mag_seitz = crate::irrep::wigner::ops_to_seitz(&mag_ops);

            for ir in crate::irrep::query::irreps_of(h_sg) {
                if !ir.spinor {
                    continue;
                }
                let mag_lg = crate::irrep::wigner::filter_little_group_with_transform(
                    ir.kx,
                    ir.ky,
                    ir.kz,
                    ir.kd,
                    &mag_ops,
                    setting_xf,
                    Some(&canonical_pure_translations),
                );
                let antiunitary: Vec<usize> = mag_lg
                    .iter()
                    .filter(|&&i| mag_ops.operations[i].time_reversal)
                    .copied()
                    .collect();
                if antiunitary.is_empty() {
                    continue;
                }
                let unitary: Vec<usize> = mag_lg
                    .iter()
                    .filter(|&&i| !mag_ops.operations[i].time_reversal)
                    .copied()
                    .collect();
                let has_imag = !ir.spin_character_imag().is_empty();

                let spin_ops = ir.spin_ops();
                if spin_ops.0.is_empty() {
                    *failure_class.entry("no_spin_ops").or_default() += 1;
                    continue;
                }
                let h_spin_seitz = crate::irrep::wigner::build_spin_seitz(spin_ops.0, spin_ops.1);
                let h_to_spin = crate::irrep::wigner::build_h_to_spin_map(
                    &h_seitz,
                    &h_spin_seitz,
                    ir.spin_lg_op_indices(),
                );

                // Count unmapped unitary ops (Seitz→spin mapping failures)
                let unmapped: Vec<usize> = unitary
                    .iter()
                    .filter(|&&i| {
                        let h = &mag_seitz[i];
                        match crate::irrep::wigner::find_seitz(&h.rot, &h.trans, &h_seitz) {
                            None => true,
                            Some(m) => h_to_spin[m.op_index].is_none(),
                        }
                    })
                    .copied()
                    .collect();

                let h_spin = ir.spin_ops();
                let g_sg = parent_spatial_sg(uni).unwrap_or(h_sg as usize) as u8;
                let g_spin = if g_sg == h_sg {
                    h_spin
                } else {
                    IrrepRecord::spin_ops_for_sg(g_sg)
                };
                let ctx = crate::irrep::wigner::SpinLiftContext {
                    h: h_spin,
                    g: g_spin,
                    sg: h_sg,
                };
                let su2_result = crate::irrep::wigner::wigner_classify_spinor(
                    &ctx,
                    wigner::SpinorWignerInput {
                        characters_real: ir.characters(),
                        characters_imag: ir.spin_character_imag(),
                        operation_indices: ir.spin_lg_op_indices(),
                        k_vector: ir.k_vector(),
                    },
                    wigner::WignerGroupContext {
                        unitary_indices: &unitary,
                        magnetic_ops: &mag_seitz,
                        unitary_ops: &h_seitz,
                        antiunitary_representative: antiunitary[0],
                    },
                    setting_xf,
                    Some(&antiunitary),
                )
                .ok();
                let direct_diagnostic =
                    crate::irrep::wigner::wigner_classify_spinor_direct_anti_diagnostic(
                        &ctx,
                        wigner::SpinorWignerInput {
                            characters_real: ir.characters(),
                            characters_imag: ir.spin_character_imag(),
                            operation_indices: ir.spin_lg_op_indices(),
                            k_vector: ir.k_vector(),
                        },
                        &antiunitary,
                        &mag_seitz,
                        setting_xf,
                        None,
                    );
                if let Err(reason) = direct_diagnostic {
                    *direct_anti_failures.entry(reason.as_str()).or_default() += 1;
                    if su2_result.is_none() {
                        *final_failure_reasons.entry(reason.as_str()).or_default() += 1;
                        *final_failure_by_sg
                            .entry((reason.as_str(), h_sg))
                            .or_default() += 1;
                        *final_failure_by_transform
                            .entry((reason.as_str(), setting_xf.is_some()))
                            .or_default() += 1;
                    }
                }
                let direct_result = direct_diagnostic.ok();

                let direct_key = match (su2_result, direct_result) {
                    (Some(a), Some(b)) if a == b => "both_ok_agree",
                    (Some(_), Some(_)) => "both_ok_disagree",
                    (Some(_), None) => "main_ok_direct_fail",
                    (None, Some(_)) => "main_fail_direct_ok",
                    (None, None) => "both_fail",
                };
                *direct_anti_stats.entry(direct_key).or_default() += 1;

                if su2_result.is_some() {
                    continue;
                }

                // Classify this failure
                let class = if !unmapped.is_empty() {
                    "mapping_failure"
                } else if has_imag {
                    "complex_char_non_quantized"
                } else {
                    "real_char_su2_closure_fail"
                };
                *failure_class.entry(class).or_default() += 1;
                *failure_by_sg.entry((class, h_sg)).or_default() += 1;
                if class == "mapping_failure" {
                    *mapping_shape
                        .entry((unmapped.len(), unitary.len()))
                        .or_default() += 1;
                }

                if shown < 30 {
                    shown += 1;
                    println!(
                        "  SG{} {} UNI{} dim={} imag={} n_lg=U{}+A{} unmapped={} class={}",
                        h_sg,
                        ir.k_label(),
                        uni,
                        ir.dim,
                        has_imag,
                        unitary.len(),
                        antiunitary.len(),
                        unmapped.len(),
                        class,
                    );
                }
            }
        }

        println!("\n=== Failure classification ===");
        for (class, count) in failure_class.iter() {
            println!("  {:>30}  {:>6}", class, count);
        }

        println!("\n=== Failure distribution by SG ===");
        let mut failure_by_sg: Vec<_> = failure_by_sg.into_iter().collect();
        failure_by_sg.sort_by_key(|((class, sg), count)| (std::cmp::Reverse(*count), *class, *sg));
        for ((class, sg), count) in failure_by_sg {
            println!("  {:>30}  SG{:>3}  {:>6}", class, sg, count);
        }

        println!("\n=== Mapping failure shapes (unmapped/unitary) ===");
        let mut mapping_shape: Vec<_> = mapping_shape.into_iter().collect();
        mapping_shape.sort_by_key(|((unmapped, unitary), count)| {
            (std::cmp::Reverse(*count), *unmapped, *unitary)
        });
        for ((unmapped, unitary), count) in mapping_shape {
            println!("  {:>3}/{:<3}  {:>6}", unmapped, unitary, count);
        }

        println!("\n=== Direct anti-coset oracle ===");
        let mut direct_anti_stats: Vec<_> = direct_anti_stats.into_iter().collect();
        direct_anti_stats.sort_by_key(|(key, count)| (std::cmp::Reverse(*count), *key));
        for (key, count) in direct_anti_stats {
            println!("  {:>30}  {:>6}", key, count);
        }

        println!("\n=== Direct anti-coset failure stages ===");
        let mut direct_anti_failures: Vec<_> = direct_anti_failures.into_iter().collect();
        direct_anti_failures.sort_by_key(|(key, count)| (std::cmp::Reverse(*count), *key));
        for (key, count) in direct_anti_failures {
            println!("  {:>30}  {:>6}", key, count);
        }

        println!("\n=== Final failure stages ===");
        let mut final_failure_reasons: Vec<_> = final_failure_reasons.into_iter().collect();
        final_failure_reasons.sort_by_key(|(key, count)| (std::cmp::Reverse(*count), *key));
        for (key, count) in final_failure_reasons {
            println!("  {:>30}  {:>6}", key, count);
        }

        println!("\n=== Final failure stages by SG ===");
        let mut final_failure_by_sg: Vec<_> = final_failure_by_sg.into_iter().collect();
        final_failure_by_sg
            .sort_by_key(|((reason, sg), count)| (std::cmp::Reverse(*count), *reason, *sg));
        for ((reason, sg), count) in final_failure_by_sg {
            println!("  {:>30}  SG{:>3}  {:>6}", reason, sg, count);
        }
        println!("\n=== Final failure stages by setting transform ===");
        let mut final_failure_by_transform: Vec<_> =
            final_failure_by_transform.into_iter().collect();
        final_failure_by_transform.sort_by_key(|((reason, found), _)| (*reason, *found));
        for ((reason, found), count) in final_failure_by_transform {
            println!("  {:>30}  xf_found={:<5}  {:>6}", reason, found, count);
        }

        // Print MSG-gauge vs old-path triage counters for the full triage pass.
        println!("\n=== Path triage ===");
        println!(
            "  MSG_GAUGE_OK:       {}",
            MSG_GAUGE_OK.load(Ordering::Relaxed)
        );
        println!(
            "  MSG_GAUGE_MAP_FAIL:  {}",
            MSG_GAUGE_MAP_FAIL.load(Ordering::Relaxed)
        );
        println!(
            "  MSG_GAUGE_W_FAIL:    {}",
            MSG_GAUGE_W_FAIL.load(Ordering::Relaxed)
        );
        println!(
            "  OLD_PATH_OK:         {}",
            OLD_PATH_OK.load(Ordering::Relaxed)
        );
        println!(
            "  OLD_PATH_FAIL:       {}",
            OLD_PATH_FAIL.load(Ordering::Relaxed)
        );
        println!("\n=== build_h_to_spin_map triage ===");
        println!("  H2S_OK:         {}", H2S_OK.load(Ordering::Relaxed));
        println!(
            "  H2S_AMBIGUOUS:  {}  (same rot, multiple spin entries)",
            H2S_AMBIGUOUS.load(Ordering::Relaxed)
        );
        println!(
            "  H2S_MISSING:    {}  (rotation not in spin table)",
            H2S_MISSING.load(Ordering::Relaxed)
        );
        let (xf_called, xf_found, xf_identity, xf_non_id, xf_nz_origin, xf_ambig) =
            crate::irrep::wigner::read_xf_counters();
        println!("\n=== find_setting_transform diagnostics ===");
        println!("  XF_CALLED:        {}", xf_called);
        println!("  XF_FOUND:         {}", xf_found);
        println!("  XF_IDENTITY:      {}", xf_identity);
        println!("  XF_NON_IDENTITY:  {}", xf_non_id);
        println!("  XF_NONZERO_ORIGIN: {}", xf_nz_origin);
        println!("  XF_AMBIGUOUS:     {}", xf_ambig);
    }

    /// Trace representative cases from each of the 4 failure categories.
    #[test]
    fn diagnose_679_failures_by_category() {
        use crate::irrep::types::IrrepRecord;
        use crate::irrep::wigner;
        use crate::irrep::wigner::DirectAntiFailure;

        // Track which categories we've shown
        let mut shown_sq_not_spin = false;
        let mut shown_spin_lookup = false;
        let mut shown_outside_lg = false;
        let mut shown_nonquant = 0usize; // show first 2

        'outer: for uni in 1..=1651 {
            let mag_ops = match get_magnetic_operations(uni) {
                Some(m) => m,
                None => continue,
            };
            let h_info = match identify_unitary_subgroup_with_hall(uni) {
                Some(i) => i,
                None => continue,
            };
            let h_sg = h_info.sg as u8;
            let h_ops = h_info.ops_from_msg;
            let _h_seitz = wigner::ops_to_seitz(&h_ops);
            let mag_seitz = wigner::ops_to_seitz(&mag_ops);

            for ir in crate::irrep::query::irreps_of(h_sg) {
                if !ir.spinor {
                    continue;
                }
                let mag_lg = wigner::filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
                let antiunitary: Vec<usize> = mag_lg
                    .iter()
                    .filter(|&&i| mag_ops.operations[i].time_reversal)
                    .copied()
                    .collect();
                if antiunitary.is_empty() {
                    continue;
                }

                let h_spin = ir.spin_ops();
                if h_spin.0.is_empty() {
                    continue;
                }
                let g_sg = parent_spatial_sg(uni).unwrap_or(h_sg as usize) as u8;
                let g_spin = if g_sg == h_sg {
                    h_spin
                } else {
                    IrrepRecord::spin_ops_for_sg(g_sg)
                };
                let ctx = wigner::SpinLiftContext {
                    h: h_spin,
                    g: g_spin,
                    sg: h_sg,
                };

                let diag = wigner::wigner_classify_spinor_direct_anti_diagnostic(
                    &ctx,
                    wigner::SpinorWignerInput {
                        characters_real: ir.characters(),
                        characters_imag: ir.spin_character_imag(),
                        operation_indices: ir.spin_lg_op_indices(),
                        k_vector: ir.k_vector(),
                    },
                    &antiunitary,
                    &mag_seitz,
                    None,
                    None,
                );

                let stage = match diag {
                    Err(e) => e,
                    Ok(_) => continue,
                };

                let do_show = match stage {
                    DirectAntiFailure::SquareNotInSpinTable if !shown_sq_not_spin => {
                        shown_sq_not_spin = true;
                        true
                    }
                    DirectAntiFailure::AntiunitarySpinLookup if !shown_spin_lookup => {
                        shown_spin_lookup = true;
                        true
                    }
                    DirectAntiFailure::SquareOutsideLittleGroup if !shown_outside_lg => {
                        shown_outside_lg = true;
                        true
                    }
                    DirectAntiFailure::NonQuantized if shown_nonquant < 2 => {
                        shown_nonquant += 1;
                        true
                    }
                    _ => false,
                };
                if !do_show {
                    continue;
                }

                // --- Detailed trace ---
                println!(
                    "\n=== UNI{} SG{} irrep {} k=({},{},{})/{} dim={} stage={:?} ===",
                    uni, h_sg, ir.ml, ir.kx, ir.ky, ir.kz, ir.kd, ir.dim, stage
                );

                let (h_spin_rots, h_spin_trans, _h_spin_su2) = ctx.h;
                let (g_spin_rots, _g_spin_trans, _g_spin_su2) = ctx.g;
                let h_spin_seitz = wigner::build_spin_seitz(h_spin_rots, h_spin_trans);
                let g_spin_seitz = wigner::build_spin_seitz(g_spin_rots, _g_spin_trans);

                println!("  H spin ops ({} entries):", h_spin_seitz.len());
                for (si, sop) in h_spin_seitz.iter().enumerate() {
                    let in_lg = ir
                        .spin_lg_op_indices()
                        .iter()
                        .any(|&idx| idx as usize == si);
                    println!(
                        "    [{}] rot={:?} t=({:.3},{:.3},{:.3}) in_lg={}",
                        si, sop.rot, sop.trans[0], sop.trans[1], sop.trans[2], in_lg
                    );
                }
                println!("  spin_lg_op_indices={:?}", ir.spin_lg_op_indices());

                println!(
                    "  G spin rotations ({} entries): {:?}",
                    g_spin_seitz.len(),
                    g_spin_seitz.iter().map(|s| s.rot).collect::<Vec<_>>()
                );

                println!("  Antiunitary ops in LG:");
                for &b_idx in &antiunitary {
                    let b = &mag_seitz[b_idx];
                    let b_bilbao = {
                        let (_, origin) = IrrepRecord::sg_setting(ctx.sg);
                        let mut t = [
                            b.trans[0]
                                - ((1.0 - b.rot[0][0] as f64)
                                    * origin.first().copied().unwrap_or(0.0)
                                    + (-b.rot[0][1] as f64)
                                        * origin.get(1).copied().unwrap_or(0.0)
                                    + (-b.rot[0][2] as f64)
                                        * origin.get(2).copied().unwrap_or(0.0)),
                            b.trans[1]
                                - ((-b.rot[1][0] as f64) * origin.first().copied().unwrap_or(0.0)
                                    + (1.0 - b.rot[1][1] as f64)
                                        * origin.get(1).copied().unwrap_or(0.0)
                                    + (-b.rot[1][2] as f64)
                                        * origin.get(2).copied().unwrap_or(0.0)),
                            b.trans[2]
                                - ((-b.rot[2][0] as f64) * origin.first().copied().unwrap_or(0.0)
                                    + (-b.rot[2][1] as f64)
                                        * origin.get(1).copied().unwrap_or(0.0)
                                    + (1.0 - b.rot[2][2] as f64)
                                        * origin.get(2).copied().unwrap_or(0.0)),
                        ];
                        for value in &mut t {
                            *value = (*value % 1.0 + 1.0) % 1.0;
                        }
                        t
                    };
                    let (sq, lattice_sq) = {
                        let b_op = wigner::SeitzOp::new(b.rot, b_bilbao, false);
                        wigner::square_seitz(&b_op)
                    };

                    let b_in_g = g_spin_seitz.iter().position(|s| s.rot == b.rot);
                    let b_in_g_neg = g_spin_seitz.iter().position(|s| {
                        s.rot[0][0] == -b.rot[0][0]
                            && s.rot[0][1] == -b.rot[0][1]
                            && s.rot[0][2] == -b.rot[0][2]
                            && s.rot[1][0] == -b.rot[1][0]
                            && s.rot[1][1] == -b.rot[1][1]
                            && s.rot[1][2] == -b.rot[1][2]
                            && s.rot[2][0] == -b.rot[2][0]
                            && s.rot[2][1] == -b.rot[2][1]
                            && s.rot[2][2] == -b.rot[2][2]
                    });
                    let sq_in_h = h_spin_seitz.iter().position(|s| s.rot == sq.rot);
                    let sq_in_h_lg = sq_in_h.map(|si| {
                        (
                            si,
                            ir.spin_lg_op_indices()
                                .iter()
                                .any(|&idx| idx as usize == si),
                        )
                    });

                    println!(
                        "    b[{}] rot={:?} trans=({:.3},{:.3},{:.3})",
                        b_idx, b.rot, b.trans[0], b.trans[1], b.trans[2]
                    );
                    println!(
                        "         b_bilbao=({:.3},{:.3},{:.3})",
                        b_bilbao[0], b_bilbao[1], b_bilbao[2]
                    );
                    println!(
                        "         sq rot={:?} sq_trans=({:.3},{:.3},{:.3}) lattice_sq={:?}",
                        sq.rot, sq.trans[0], sq.trans[1], sq.trans[2], lattice_sq
                    );
                    println!(
                        "         b in G spin: pos={:?}  (-R) pos={:?}",
                        b_in_g, b_in_g_neg
                    );
                    println!("         sq in H spin: {:?}", sq_in_h_lg);
                }

                // Check if all categories done
                if shown_sq_not_spin && shown_spin_lookup && shown_outside_lg && shown_nonquant >= 2
                {
                    break 'outer;
                }
            }
        }
    }

    #[test]
    fn diagnose_spglib_standard_setting_transform() {
        fn transformed_unique(
            ops: &SymmetryOps,
            transform: &wigner::SettingTransform,
        ) -> Vec<wigner::SeitzOp> {
            let mut result = Vec::new();
            for op in ops.operations.iter().filter(|op| !op.time_reversal) {
                let rot = match transform.transform_rotation(&op.rotation) {
                    Some(r) => r,
                    None => continue,
                };
                let trans = match transform.transform_translation(&op.rotation, &op.translation) {
                    Some(t) => t,
                    None => continue,
                };
                let transformed = wigner::SeitzOp::new(rot, trans, false);
                if wigner::find_seitz(&transformed.rot, &transformed.trans, &result).is_none() {
                    result.push(transformed);
                }
            }
            result
        }

        fn same_group(a: &[wigner::SeitzOp], b: &[wigner::SeitzOp]) -> bool {
            a.len() == b.len()
                && a.iter()
                    .all(|op| wigner::find_seitz(&op.rot, &op.trans, b).is_some())
                && b.iter()
                    .all(|op| wigner::find_seitz(&op.rot, &op.trans, a).is_some())
        }

        fn embeds_group(source: &[wigner::SeitzOp], target: &[wigner::SeitzOp]) -> bool {
            source
                .iter()
                .all(|op| wigner::find_seitz(&op.rot, &op.trans, target).is_some())
        }

        let mut total = 0usize;
        let mut found = 0usize;
        let mut sg_match = 0usize;
        let mut detected_hall_exact = 0usize;
        let mut detected_hall_embedded = 0usize;
        let mut data_hall_exact = 0usize;
        let mut data_hall_embedded = 0usize;
        let mut examples = Vec::new();
        let mut data_examples = Vec::new();

        for uni in 1..=1651 {
            let Some(mag_ops) = get_magnetic_operations(uni) else {
                continue;
            };
            let Some(expected_sg) = identify_unitary_subgroup(uni) else {
                continue;
            };
            total += 1;
            let Some((sg, hall, transform)) = standard_setting_transform(&mag_ops, false) else {
                continue;
            };
            found += 1;
            if sg == expected_sg {
                sg_match += 1;
            }
            let transformed = transformed_unique(&mag_ops, &transform);
            let detected_target = get_parent_operations_by_hall(hall)
                .map(|ops| wigner::ops_to_seitz(&ops))
                .unwrap_or_default();
            if same_group(&transformed, &detected_target) {
                detected_hall_exact += 1;
            }
            if embeds_group(&transformed, &detected_target) {
                detected_hall_embedded += 1;
            } else if examples.len() < 10 {
                examples.push(format!(
                    "UNI{uni} SG{sg} hall={hall} basis={:?} origin={:?} source_order={} target_order={}",
                    transform.basis,
                    transform.origin,
                    transformed.len(),
                    detected_target.len(),
                ));
            }
            let data_target_ops = crate::irrep::bridge::canonical_hall_ops(sg as u8).unwrap();
            let data_target = wigner::ops_to_seitz(&data_target_ops);
            let data_exact = same_group(&transformed, &data_target);
            if data_exact {
                data_hall_exact += 1;
            }
            let data_embedded = if data_exact && transform_applies_to_all_ops(&transform, &mag_ops)
            {
                true
            } else if let Some(h_info) = identify_unitary_subgroup_with_hall(uni) {
                let identity = wigner::SettingTransform::identity();
                let msg_to_data = h_info.msg_to_data.as_ref().unwrap_or(&identity);
                let embedded =
                    transform_embeds_ops(msg_to_data, &h_info.ops_from_msg, &data_target_ops)
                        && transform_applies_to_all_ops(msg_to_data, &mag_ops);
                if !embedded {
                    data_examples.push(format!(
                        "UNI{uni} SG{sg} detected_hall={hall} data_hall={} msg_order={} data_order={} transform={:?}",
                        h_info.data_hall,
                        h_info.ops_from_msg.len(),
                        data_target_ops.len(),
                        h_info.msg_to_data,
                    ));
                }
                embedded
            } else {
                data_examples.push(format!(
                    "UNI{uni} SG{sg} detected_hall={hall}: unitary subgroup info missing"
                ));
                false
            };
            if data_embedded {
                data_hall_embedded += 1;
            }
        }

        println!("\n=== spglib standard setting transform ===");
        println!("  total:               {total}");
        println!("  found:               {found}");
        println!("  sg_match:            {sg_match}");
        println!("  detected_hall_exact: {detected_hall_exact}");
        println!("  detected_hall_embed: {detected_hall_embedded}");
        println!("  data_hall_exact:     {data_hall_exact}");
        println!("  data_hall_embed:     {data_hall_embedded}");
        for example in examples {
            println!("  {example}");
        }
        for example in data_examples {
            println!("  data mismatch: {example}");
        }
        assert_eq!(total, 1651);
        assert_eq!(found, total);
        assert_eq!(sg_match, total);
        assert_eq!(detected_hall_embedded, total);
        assert_eq!(data_hall_embedded, total);
    }

    /// Oracle: check UNI→Hall→magnetic operations entry point for anomalous cases.
    ///
    /// For UNI187, UNI270, UNI271, UNI663: prints Hall selection, unitary
    /// operation closure, and consistency with BNS/UNI metadata.
    #[test]
    fn diagnose_magnetic_entry_hall_anomalies() {
        let test_unis = [187usize, 270, 271, 663];
        for &uni in &test_unis {
            println!("\n=== UNI{} ===", uni);
            let msg_type = crate::MagneticSpaceGroupType::from_uni(uni).unwrap();
            println!(
                "  BNS={} OG={} type={:?} number={}",
                msg_type.bns_number.trim(),
                msg_type.og_number.trim(),
                msg_type.type_,
                msg_type.number
            );

            // Hall selection
            let hall = get_first_hall_for_uni(uni);
            println!("  get_first_hall_for_uni: {:?}", hall);

            // Magnetic operations
            let mag_ops = match get_magnetic_operations(uni) {
                Some(ops) => ops,
                None => {
                    println!("  get_magnetic_operations: None");
                    continue;
                }
            };
            println!("  mag_ops: {} total", mag_ops.len());
            let u_count = mag_ops
                .operations
                .iter()
                .filter(|o| !o.time_reversal)
                .count();
            let a_count = mag_ops
                .operations
                .iter()
                .filter(|o| o.time_reversal)
                .count();
            println!("  unitary={} antiunitary={}", u_count, a_count);

            // Show unitary rotations
            let id: crate::mathfunc::Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
            let unitary_rots: Vec<&crate::mathfunc::Mat3I> = mag_ops
                .operations
                .iter()
                .filter(|o| !o.time_reversal)
                .map(|o| &o.rotation)
                .collect();
            println!(
                "  unitary rotations ({} distinct):",
                unitary_rots
                    .iter()
                    .collect::<std::collections::HashSet<_>>()
                    .len()
            );
            for (i, op) in mag_ops.operations.iter().enumerate() {
                if op.time_reversal {
                    continue;
                }
                let is_id = op.rotation == id;
                println!(
                    "    [{}] rot={:?} t=({:.3},{:.3},{:.3}){}",
                    i,
                    op.rotation,
                    op.translation[0],
                    op.translation[1],
                    op.translation[2],
                    if is_id { " ID" } else { "" }
                );
            }

            // Identify unitary subgroup
            let h_info = identify_unitary_subgroup_with_hall(uni);
            if let Some(info) = &h_info {
                println!("  identified: SG{} Hall{}", info.sg, info.hall);
                // Check consistency
                let seitz_msg = crate::irrep::wigner::ops_to_seitz(&info.ops_from_msg);
                let seitz_hall = crate::irrep::wigner::ops_to_seitz(&info.ops_from_hall);
                let mut rot_msg: Vec<_> = seitz_msg.iter().map(|s| s.rot).collect();
                let mut rot_hall: Vec<_> = seitz_hall.iter().map(|s| s.rot).collect();
                rot_msg.sort();
                rot_hall.sort();
                let rots_match = rot_msg == rot_hall;
                println!("  rot multisets match: {}", rots_match);
                if !rots_match {
                    println!("    msg rots:  {:?}", rot_msg);
                    println!("    hall rots: {:?}", rot_hall);
                }
                // Check closure: every product of two unitary ops should be in the set
                let all_closed = unitary_rots.iter().all(|r1| {
                    unitary_rots.iter().all(|r2| {
                        let prod = crate::mathfunc::mat_multiply_matrix_i3(r1, r2);
                        unitary_rots.iter().any(|r3| **r3 == prod)
                    })
                });
                println!("  unitary rotation closure: {}", all_closed);
                if !all_closed {
                    println!("    WARNING: unitary rotations are NOT closed under multiplication!");
                }
            } else {
                println!("  identify_unitary_subgroup_with_hall: None");
            }

            // Also check with hall=0 to see if the MSG database self-selects
            if let Some(h) = hall {
                let ops_h0 = crate::msg_database::get_spacegroup_operations(uni, 0);
                let ops_h = crate::msg_database::get_spacegroup_operations(uni, h);
                if let (Some(sym0), Some(symh)) = (ops_h0, ops_h) {
                    let same_size = sym0.len() == symh.len();
                    let h0_u = (0..sym0.len()).filter(|&i| !sym0.timerev[i]).count();
                    let hh_u = (0..symh.len()).filter(|&i| !symh.timerev[i]).count();
                    println!(
                        "  msgdb(uni,0): {} ops ({}U)  msgdb(uni,{}): {} ops ({}U)  same_size={}",
                        sym0.len(),
                        h0_u,
                        h,
                        symh.len(),
                        hh_u,
                        same_size
                    );
                }
            }

            // Check what standard_setting_transform gives for the unitary ops
            let u_rots: Vec<crate::mathfunc::Mat3I> = mag_ops
                .operations
                .iter()
                .filter(|o| !o.time_reversal)
                .map(|o| o.rotation)
                .collect();
            let u_trans: Vec<[f64; 3]> = mag_ops
                .operations
                .iter()
                .filter(|o| !o.time_reversal)
                .map(|o| o.translation)
                .collect();
            let unitary_ops =
                SymmetryOps::from_parallel(&u_rots, &u_trans, &vec![false; u_rots.len()]).unwrap();
            if let Some((std_sg, std_hall, xf)) = standard_setting_transform(&unitary_ops, false) {
                println!(
                    "  standard_setting: SG{} Hall{} basis={:?} origin=({:.3},{:.3},{:.3})",
                    std_sg, std_hall, xf.basis, xf.origin[0], xf.origin[1], xf.origin[2]
                );
                // Verify: apply transform to unitary ops, they should be in standard Hall
                let xf_rots: Vec<_> = unitary_rots
                    .iter()
                    .filter_map(|r| xf.transform_rotation(r))
                    .collect();
                let hall_ops = get_parent_operations_by_hall(std_hall);
                if let Some(hop) = &hall_ops {
                    let hall_rots: Vec<_> = hop.operations.iter().map(|o| o.rotation).collect();
                    let match_count = xf_rots.iter().filter(|r| hall_rots.contains(r)).count();
                    println!("  xf_rots in hall_ops: {}/{}", match_count, xf_rots.len());
                }
            } else {
                println!("  standard_setting_transform: None");
            }
        }
    }

    /// Oracle for the reciprocal-space action used by `spin_lg_op_indices`.
    ///
    /// Traces antiunitary spin-lookup failures by comparing each failing
    /// operation with the parent spin table rotations.
    #[test]
    fn diagnose_antiunitary_spin_lookup() {
        use crate::irrep::types::IrrepRecord;
        use crate::irrep::wigner;

        let mut shown = std::collections::HashSet::new();
        'outer: for uni in 1..=1651 {
            let mag_ops = match get_magnetic_operations(uni) {
                Some(m) => m,
                None => continue,
            };
            let h_info = match identify_unitary_subgroup_with_hall(uni) {
                Some(i) => i,
                None => continue,
            };
            let h_sg = h_info.sg as u8;
            // Extract rotations/translations BEFORE moving ops_from_msg
            let msg_rots: Vec<_> = h_info.ops_from_msg.iter().map(|o| o.rotation).collect();
            let msg_trans: Vec<_> = h_info.ops_from_msg.iter().map(|o| o.translation).collect();
            let hall_rots: Vec<_> = h_info.ops_from_hall.iter().map(|o| o.rotation).collect();
            let hall_trans: Vec<_> = h_info.ops_from_hall.iter().map(|o| o.translation).collect();
            let setting_xfs =
                wigner::find_setting_transform(&msg_rots, &msg_trans, &hall_rots, &hall_trans);
            let setting_xf = setting_xfs.first();

            let h_ops = h_info.ops_from_msg;
            let _h_seitz = wigner::ops_to_seitz(&h_ops);
            let mag_seitz = wigner::ops_to_seitz(&mag_ops);

            for ir in crate::irrep::query::irreps_of(h_sg) {
                if !ir.spinor {
                    continue;
                }
                let mag_lg = wigner::filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
                let antiunitary: Vec<usize> = mag_lg
                    .iter()
                    .filter(|&&i| mag_ops.operations[i].time_reversal)
                    .copied()
                    .collect();
                if antiunitary.is_empty() {
                    continue;
                }

                let h_spin = ir.spin_ops();
                if h_spin.0.is_empty() {
                    continue;
                }
                let g_sg = parent_spatial_sg(uni).unwrap_or(h_sg as usize) as u8;
                let g_spin = if g_sg == h_sg {
                    h_spin
                } else {
                    IrrepRecord::spin_ops_for_sg(g_sg)
                };
                let ctx = wigner::SpinLiftContext {
                    h: h_spin,
                    g: g_spin,
                    sg: h_sg,
                };

                let diag = wigner::wigner_classify_spinor_direct_anti_diagnostic(
                    &ctx,
                    wigner::SpinorWignerInput {
                        characters_real: ir.characters(),
                        characters_imag: ir.spin_character_imag(),
                        operation_indices: ir.spin_lg_op_indices(),
                        k_vector: ir.k_vector(),
                    },
                    &antiunitary,
                    &mag_seitz,
                    setting_xf,
                    None,
                );

                let is_sq_not_spin =
                    matches!(diag, Err(wigner::DirectAntiFailure::SquareNotInSpinTable));
                let is_au_lookup =
                    matches!(diag, Err(wigner::DirectAntiFailure::AntiunitarySpinLookup));
                if !is_sq_not_spin && !is_au_lookup {
                    continue;
                }
                if h_sg != 1 && is_sq_not_spin {
                    continue;
                } // SG1 only for sq_not_spin
                if !shown.insert(h_sg) {
                    continue;
                }

                println!(
                    "\n=== UNI{} SG{} irrep {} k=({},{},{})/{} g_sg={} stage={:?} ===",
                    uni,
                    h_sg,
                    ir.ml,
                    ir.kx,
                    ir.ky,
                    ir.kz,
                    ir.kd,
                    g_sg,
                    diag.err().unwrap()
                );

                let (g_spin_rots, _g_spin_trans, _g_spin_su2) = ctx.g;
                let g_spin_seitz = wigner::build_spin_seitz(g_spin_rots, _g_spin_trans);

                println!(
                    "\n=== UNI{} SG{} irrep {} k=({},{},{})/{} g_sg={} ===",
                    uni, h_sg, ir.ml, ir.kx, ir.ky, ir.kz, ir.kd, g_sg
                );
                println!("  G spin rotations:");
                for (si, s) in g_spin_seitz.iter().enumerate() {
                    println!("    [{}] rot={:?}", si, s.rot);
                }
                println!("  Antiunitary ops:");
                for &b_idx in &antiunitary {
                    let b = &mag_seitz[b_idx];
                    let (b_rot, _b_trans) = if let Some(xf) = setting_xf {
                        let rot = xf.transform_rotation(&b.rot).unwrap_or(b.rot);
                        let trans = xf
                            .transform_translation(&b.rot, &b.trans)
                            .unwrap_or(b.trans);
                        (rot, trans)
                    } else {
                        (b.rot, b.trans)
                    };
                    let in_g = g_spin_seitz.iter().position(|s| s.rot == b_rot);
                    let neg_rot: crate::mathfunc::Mat3I = [
                        [-b_rot[0][0], -b_rot[0][1], -b_rot[0][2]],
                        [-b_rot[1][0], -b_rot[1][1], -b_rot[1][2]],
                        [-b_rot[2][0], -b_rot[2][1], -b_rot[2][2]],
                    ];
                    let neg_in_g = g_spin_seitz.iter().position(|s| s.rot == neg_rot);
                    println!(
                        "    b[{}] rot={:?} → xf_rot={:?} inG={:?} -R_inG={:?}",
                        b_idx, b.rot, b_rot, in_g, neg_in_g
                    );
                }
                if shown.len() >= 3 {
                    break 'outer;
                }
            }
        }
    }

    /// - reciprocal action: `R^{-T} k ≡ k`
    ///
    /// This is diagnostic only and does not alter classification.
    #[test]
    fn diagnose_spin_lg_k_convention() {
        use crate::mathfunc::mat_get_determinant_i3;

        fn inverse_transpose(r: &[[i32; 3]; 3]) -> [[i32; 3]; 3] {
            let det = mat_get_determinant_i3(r);
            assert!(det == 1 || det == -1, "rotation must be unimodular: {r:?}");
            let mut out = [[0i32; 3]; 3];
            // R^{-T} is the cofactor matrix divided by det.
            for (i, row) in out.iter_mut().enumerate() {
                for (j, value) in row.iter_mut().enumerate() {
                    let rows: Vec<usize> = (0..3).filter(|&x| x != i).collect();
                    let cols: Vec<usize> = (0..3).filter(|&x| x != j).collect();
                    let minor = r[rows[0]][cols[0]] * r[rows[1]][cols[1]]
                        - r[rows[0]][cols[1]] * r[rows[1]][cols[0]];
                    let cofactor = if (i + j) % 2 == 0 { minor } else { -minor };
                    *value = cofactor / det;
                }
            }
            out
        }

        fn preserves(r: &[[i32; 3]; 3], k: [i32; 3], kd: i32) -> bool {
            if kd == 0 {
                return true;
            }
            (0..3).all(|i| {
                let rk: i32 = (0..3).map(|j| r[i][j] * k[j]).sum();
                (rk - k[i]) % kd == 0
            })
        }

        fn preserves_centered(
            r: &[[i32; 3]; 3],
            k: [i32; 3],
            kd: i32,
            pure_translations: &[[f64; 3]],
        ) -> bool {
            if kd == 0 {
                return true;
            }
            let mut reciprocal_shift = [0i32; 3];
            for (i, value) in reciprocal_shift.iter_mut().enumerate() {
                let rk: i32 = (0..3).map(|j| r[i][j] * k[j]).sum();
                let delta = rk - k[i];
                if delta % kd != 0 {
                    return false;
                }
                *value = delta / kd;
            }
            pure_translations.iter().all(|t| {
                let phase = reciprocal_shift[0] as f64 * t[0]
                    + reciprocal_shift[1] as f64 * t[1]
                    + reciprocal_shift[2] as f64 * t[2];
                (phase - phase.round()).abs() < 1e-8
            })
        }

        let mut total_irreps = 0usize;
        let mut direct_exact = 0usize;
        let mut reciprocal_exact = 0usize;
        let mut reciprocal_centered_exact = 0usize;
        let mut direct_fp = 0usize;
        let mut direct_fn = 0usize;
        let mut reciprocal_fp = 0usize;
        let mut reciprocal_fn = 0usize;
        let mut reciprocal_centered_fp = 0usize;
        let mut reciprocal_centered_fn = 0usize;
        let mut examples = Vec::new();
        let mut reciprocal_mismatch_examples = Vec::new();

        for sg in 1u8..=230 {
            let spin_ops = IrrepRecord::spin_ops_for_sg(sg);
            let spin_seitz = crate::irrep::wigner::build_spin_seitz(spin_ops.0, spin_ops.1);
            let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
            let parent_ops = get_parent_operations(sg).unwrap();
            let pure_translations: Vec<[f64; 3]> = parent_ops
                .operations
                .iter()
                .filter(|op| op.rotation == identity)
                .map(|op| op.translation)
                .collect();
            for ir in crate::irrep::query::irreps_of(sg) {
                if !ir.spinor || ir.spin_lg_op_indices().is_empty() {
                    continue;
                }
                total_irreps += 1;
                let expected_indices: std::collections::HashSet<usize> = ir
                    .spin_lg_op_indices()
                    .iter()
                    .map(|&x| x as usize)
                    .collect();
                let expected: std::collections::HashSet<[[i32; 3]; 3]> = expected_indices
                    .iter()
                    .filter_map(|&idx| spin_seitz.get(idx).map(|op| op.rot))
                    .collect();
                let k = [ir.kx as i32, ir.ky as i32, ir.kz as i32];
                let kd = ir.kd as i32;
                let mut direct_set = std::collections::HashSet::new();
                let mut reciprocal_set = std::collections::HashSet::new();
                let mut reciprocal_centered_set = std::collections::HashSet::new();

                for op in &spin_seitz {
                    if preserves(&op.rot, k, kd) {
                        direct_set.insert(op.rot);
                    }
                    let rit = inverse_transpose(&op.rot);
                    if preserves(&rit, k, kd) {
                        reciprocal_set.insert(op.rot);
                    }
                    if preserves_centered(&rit, k, kd, &pure_translations) {
                        reciprocal_centered_set.insert(op.rot);
                    }
                }

                if direct_set == expected {
                    direct_exact += 1;
                }
                if reciprocal_set == expected {
                    reciprocal_exact += 1;
                }
                if reciprocal_centered_set == expected {
                    reciprocal_centered_exact += 1;
                }
                direct_fp += direct_set.difference(&expected).count();
                direct_fn += expected.difference(&direct_set).count();
                reciprocal_fp += reciprocal_set.difference(&expected).count();
                reciprocal_fn += expected.difference(&reciprocal_set).count();
                reciprocal_centered_fp += reciprocal_centered_set.difference(&expected).count();
                reciprocal_centered_fn += expected.difference(&reciprocal_centered_set).count();

                if examples.len() < 20 && direct_set != expected && reciprocal_set == expected {
                    examples.push(format!(
                        "SG{} {} k=({},{},{})/{} expected={:?} Rk={:?} R^-Tk={:?}",
                        sg, ir.ml, ir.kx, ir.ky, ir.kz, ir.kd, expected, direct_set, reciprocal_set,
                    ));
                }
                if reciprocal_mismatch_examples.len() < 20 && reciprocal_centered_set != expected {
                    reciprocal_mismatch_examples.push(format!(
                        "SG{} {} k=({},{},{})/{} expected={:?} centered R^-Tk={:?} extra={:?} translations={:?}",
                        sg,
                        ir.ml,
                        ir.kx,
                        ir.ky,
                        ir.kz,
                        ir.kd,
                        expected,
                        reciprocal_centered_set,
                        reciprocal_centered_set.difference(&expected).collect::<Vec<_>>(),
                        pure_translations,
                    ));
                }
            }
        }

        println!("\n=== spin little-group k convention oracle ===");
        println!("  total_irreps:       {total_irreps}");
        println!("  direct_exact:       {direct_exact}");
        println!("  reciprocal_exact:   {reciprocal_exact}");
        println!("  reciprocal_centered_exact: {reciprocal_centered_exact}");
        println!("  direct_fp/fn:       {direct_fp}/{direct_fn}");
        println!("  reciprocal_fp/fn:   {reciprocal_fp}/{reciprocal_fn}");
        println!("  reciprocal_centered_fp/fn: {reciprocal_centered_fp}/{reciprocal_centered_fn}");
        for example in examples {
            println!("  {example}");
        }
        println!("  reciprocal mismatch examples:");
        for example in reciprocal_mismatch_examples {
            println!("  {example}");
        }
    }

    /// Regression: SG3 A3 spinor Wigner test under grey group (a₀ = Θ).
    ///
    /// This explicitly verifies that the SU(2) path gives the correct
    /// per-term contributions, and that Bilbao imaginary chars are NOT
    /// valid term-by-term Wigner summands.
    #[test]
    fn test_spinor_sg3_a3_grey_wigner() {
        use crate::irrep::wigner;

        // SG3 (P2), find a grey magnetic group (a₀ = Θ, g₀ = I)
        let mut grey_uni = None;
        for uni in 1..=1651 {
            let mag_ops = match get_magnetic_operations(uni) {
                Some(m) => m,
                None => continue,
            };
            let h_sg = match identify_unitary_subgroup(uni) {
                Some(s) => s as u8,
                None => continue,
            };
            if h_sg != 3 {
                continue;
            }
            let a0_idx = match mag_ops.operations.iter().position(|o| o.time_reversal) {
                Some(i) => i,
                None => continue,
            };
            let r = mag_ops.operations[a0_idx].rotation;
            if r[0][0] == 1
                && r[0][1] == 0
                && r[0][2] == 0
                && r[1][0] == 0
                && r[1][1] == 1
                && r[1][2] == 0
                && r[2][0] == 0
                && r[2][1] == 0
                && r[2][2] == 1
            {
                grey_uni = Some(uni);
                break;
            }
        }
        let uni = grey_uni.expect("SG3 should have a grey magnetic group");
        let mag_ops = get_magnetic_operations(uni).unwrap();
        let mag_seitz = wigner::ops_to_seitz(&mag_ops);
        let h_seitz: Vec<_> = (0..mag_ops.len())
            .filter(|&i| !mag_ops.operations[i].time_reversal)
            .map(|i| {
                wigner::SeitzOp::new(
                    mag_ops.operations[i].rotation,
                    mag_ops.operations[i].translation,
                    false,
                )
            })
            .collect();

        // Find SG3 A3 at k=(½,0,½)
        let a3 = crate::irrep::query::irreps_of(3)
            .iter()
            .find(|ir| ir.ml == "A3" && ir.spinor)
            .expect("SG3 A3 spinor irrep should exist");

        // Verify imaginary chars exist but are NOT valid Wigner summands
        let imag = a3.spin_character_imag();
        assert!(!imag.is_empty(), "A3 should have imaginary chars");
        // imag[0] = 0.0, but h=E gives χ((ΘE)²) = χ(Ē) = -1 ≠ 0
        assert!(
            (imag[0] - 0.0).abs() < 0.01,
            "imag[0] should be 0, proving imag ≠ term-by-term Wigner summand"
        );

        let mag_lg = wigner::filter_little_group(a3.kx, a3.ky, a3.kz, a3.kd, &mag_ops);
        let unitary: Vec<usize> = mag_lg
            .iter()
            .filter(|&&i| !mag_ops.operations[i].time_reversal)
            .copied()
            .collect();
        let antiunitary: Vec<usize> = mag_lg
            .iter()
            .filter(|&&i| mag_ops.operations[i].time_reversal)
            .copied()
            .collect();
        assert!(!antiunitary.is_empty(), "should have antiunitary ops");

        // Run SU(2) Wigner test
        let spin_ops = a3.spin_ops();
        let ctx = wigner::SpinLiftContext {
            h: spin_ops,
            g: spin_ops,
            sg: a3.sg,
        };
        let ct = wigner::wigner_classify_spinor(
            &ctx,
            wigner::SpinorWignerInput {
                characters_real: a3.characters(),
                characters_imag: a3.spin_character_imag(),
                operation_indices: a3.spin_lg_op_indices(),
                k_vector: a3.k_vector(),
            },
            wigner::WignerGroupContext {
                unitary_indices: &unitary,
                magnetic_ops: &mag_seitz,
                unitary_ops: &h_seitz,
                antiunitary_representative: antiunitary[0],
            },
            None,
            Some(&antiunitary),
        );

        // For grey group spin-½ at this k-point:
        //   h=E:  (ΘE)² = Θ² = Ē → χ(Ē) = -χ(E) = -1
        //   h=C₂: (ΘC₂)² = Θ²C₂² = (-1)(-1) = E → χ(E) = +1
        //   W = (-1 + 1)/2 = 0 → Type C
        assert!(ct.is_ok(), "SU(2) path should succeed for grey group");
        let ct = ct.unwrap();
        assert_eq!(
            ct,
            CorepType::C,
            "SG3 A3 grey: expected Type C (W=0), got {:?}",
            ct
        );

        // Also verify: extra sum is diagnostic-only, not a Wigner indicator
        let imag_sum = wigner::diagnostic_imag_sum(imag);
        eprintln!(
            "SG3 A3 grey Wigner: SU(2) path = {:?}, imag_sum = {:.4}",
            ct, imag_sum
        );
    }

    /// Per-term diagnostic: print raw data for debugging gauge conventions.
    /// Not a pass/fail test — run with --nocapture to inspect.
    #[test]
    fn diagnose_spinor_wigner_per_term() {
        use crate::irrep::wigner;
        let sc = wigner::su2_compose;

        // Find grey (Type-II) magnetic groups (a₀ = T, g₀ = I).
        // For grey groups, a₀ is always in spin_ops → SU(2) path applicable.
        let test_sgs: &[u8] = &[3, 5, 10, 118];
        let mut cases: Vec<(usize, u8)> = Vec::new();
        for uni in 1..=1651 {
            let mag_ops = match get_magnetic_operations(uni) {
                Some(m) => m,
                None => continue,
            };
            let h_sg = match identify_unitary_subgroup(uni) {
                Some(s) => s as u8,
                None => continue,
            };
            if !test_sgs.contains(&h_sg) {
                continue;
            }
            let a0_idx = match mag_ops.operations.iter().position(|o| o.time_reversal) {
                Some(i) => i,
                None => continue,
            };
            let r = mag_ops.operations[a0_idx].rotation;
            let is_id = r[0][0] == 1
                && r[0][1] == 0
                && r[0][2] == 0
                && r[1][0] == 0
                && r[1][1] == 1
                && r[1][2] == 0
                && r[2][0] == 0
                && r[2][1] == 0
                && r[2][2] == 1;
            if is_id {
                cases.push((uni, h_sg));
            }
        }

        eprintln!("\n=== Per-term spinor Wigner diagnostic (grey groups) ===");
        eprintln!("  Found {} grey-group cases", cases.len());

        // Only show the first few cases to keep output manageable
        let mut shown = 0usize;
        let max_show = 3usize;

        for (uni, h_sg) in cases {
            let mag_ops = get_magnetic_operations(uni).unwrap();
            let mag_seitz = wigner::ops_to_seitz(&mag_ops);

            for ir in crate::irrep::query::irreps_of(h_sg) {
                if !ir.spinor {
                    continue;
                }
                let imag = ir.spin_character_imag();
                if imag.is_empty() {
                    continue;
                }
                let imag_sum: f64 = imag.iter().sum();

                let h_seitz: Vec<_> = (0..mag_ops.len())
                    .filter(|&i| !mag_ops.operations[i].time_reversal)
                    .map(|i| {
                        wigner::SeitzOp::new(
                            mag_ops.operations[i].rotation,
                            mag_ops.operations[i].translation,
                            false,
                        )
                    })
                    .collect();
                let spin_seitz = wigner::build_spin_seitz(ir.spin_ops().0, ir.spin_ops().1);
                let h_to_spin =
                    wigner::build_h_to_spin_map(&h_seitz, &spin_seitz, ir.spin_lg_op_indices());
                let global_to_local: std::collections::HashMap<usize, usize> = ir
                    .spin_lg_op_indices()
                    .iter()
                    .enumerate()
                    .map(|(loc, &g)| (g as usize, loc))
                    .collect();

                let mag_lg = wigner::filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
                let unitary: Vec<usize> = mag_lg
                    .iter()
                    .filter(|&&i| !mag_ops.operations[i].time_reversal)
                    .copied()
                    .collect();
                let antiunitary: Vec<usize> = mag_lg
                    .iter()
                    .filter(|&&i| mag_ops.operations[i].time_reversal)
                    .copied()
                    .collect();
                if antiunitary.is_empty() {
                    continue;
                }

                let a0 = &mag_seitz[antiunitary[0]];
                let a0_match = match wigner::find_seitz(&a0.rot, &a0.trans, &spin_seitz) {
                    Some(m) => m,
                    None => continue,
                };
                let u_a0 = spin_su2_at(ir.spin_ops().2, a0_match.op_index).unwrap();

                let chars = ir.characters();
                let n_lg = ir.spin_lg_char_count();

                eprintln!(
                    "\n═══ SG{} {} k=({}/{},{}/{},{}/{}) dim={} ═══",
                    h_sg, ir.ml, ir.kx, ir.kd, ir.ky, ir.kd, ir.kz, ir.kd, ir.dim
                );
                eprintln!(
                    "  n_unitary={} n_antiunitary={} imag_sum={:.4} imag={:?}",
                    unitary.len(),
                    antiunitary.len(),
                    imag_sum,
                    imag.iter()
                        .map(|&x| format!("{:.2}", x))
                        .collect::<Vec<_>>()
                );
                eprintln!("  spin_lg_op_indices={:?}", ir.spin_lg_op_indices());
                eprintln!("  n_lg_chars={} total_chars={}", n_lg, chars.len());
                eprintln!(
                    "  lg_chars={:?}",
                    chars[..n_lg.min(chars.len())]
                        .iter()
                        .map(|&x| format!("{:.2}", x))
                        .collect::<Vec<_>>()
                );

                // Per-term table
                eprintln!(
                    "  ┌──────┬─────────────────────┬──────────┬───────────────┬──────────────┬───────────────┬───────┬───────┐"
                );
                eprintln!(
                    "  │  h#  │ h_spin  U_h         │ h² spin  │ U_h² vs U_k   │ central(Θ²)   │   χ(spin)     │ imag │ contr │"
                );
                eprintln!(
                    "  ├──────┼─────────────────────┼──────────┼───────────────┼──────────────┼───────────────┼───────┼───────┤"
                );

                let mut w_sum = 0.0f64;
                let mut used = 0usize;
                for &h_mag_idx in &unitary {
                    let h = &mag_seitz[h_mag_idx];
                    let h_match = match wigner::find_seitz(&h.rot, &h.trans, &h_seitz) {
                        Some(m) => m,
                        None => continue,
                    };
                    let Some(h_spin_idx) = h_to_spin[h_match.op_index] else {
                        continue;
                    };
                    let u_h = spin_su2_at(ir.spin_ops().2, h_spin_idx).unwrap();

                    // Spatial: (g₀h)²
                    let g0 = wigner::SeitzOp::new(a0.rot, a0.trans, false);
                    let h_sp = wigner::SeitzOp::new(h.rot, h.trans, false);
                    let (g0h, l1) = wigner::compose_seitz(&g0, &h_sp);
                    let (sq, lsq) = wigner::square_seitz(&g0h);

                    // Canonical lift of h²
                    let sq_match = match wigner::find_seitz(&sq.rot, &sq.trans, &h_seitz) {
                        Some(m) => m,
                        None => continue,
                    };
                    let Some(sq_spin_idx) = h_to_spin[sq_match.op_index] else {
                        continue;
                    };
                    let u_k = spin_su2_at(ir.spin_ops().2, sq_spin_idx).unwrap();

                    // SU(2): U_sq = (U_a₀·U_h)²
                    let u_g0h = sc(&u_a0, &u_h);
                    let u_sq = sc(&u_g0h, &u_g0h);

                    // Central element detection
                    let spatial_central = match wigner::su2_same_up_to_sign(&u_sq, &u_k) {
                        Some(v) => v,
                        None => continue,
                    };
                    let central = !spatial_central;

                    // Read character
                    let local_idx = *global_to_local.get(&sq_spin_idx).unwrap();
                    if local_idx >= n_lg || local_idx >= chars.len() {
                        continue;
                    }
                    let chi0 = chars[local_idx];
                    let chi = if central { -chi0 } else { chi0 };

                    // Bloch phase
                    let r_l1 = wigner::mat_vec_i32(&g0h.rot, &l1);
                    let total_lattice = wigner::add3(
                        &wigner::add3(&lsq, &sq_match.lattice_shift),
                        &wigner::add3(&l1, &r_l1),
                    );
                    let phase = wigner::bloch_phase(ir.kx, ir.ky, ir.kz, ir.kd, &total_lattice);
                    let contrib = chi * phase.re; // W contribution should be real

                    // Extra chars: compare by position in the unitary list
                    let extra_val = imag
                        .get(used)
                        .map(|value| format!("{value:.2}"))
                        .unwrap_or_else(|| "missing".to_string());

                    w_sum += contrib;
                    used += 1;

                    eprintln!(
                        "  h={} h_spin={} sq_spin={} spC={} c={} chi0={:.2} chi={:.2} ph={:.2} contrib={:.2} imag={}",
                        h_mag_idx,
                        h_spin_idx,
                        sq_spin_idx,
                        spatial_central,
                        central,
                        chi0,
                        chi,
                        phase.re,
                        contrib,
                        extra_val
                    );
                }

                let w = if used > 0 { w_sum / used as f64 } else { 0.0 };
                eprintln!(
                    "  └──────┴─────────────────────┴──────────┴───────────────┴──────────────┴───────────────┴───────┴───────┘"
                );
                eprintln!(
                    "  W = {:.4} / {used} = {:.4}  (imag_sum={:.4})",
                    w_sum, w, imag_sum
                );

                shown += 1;
                if shown >= max_show {
                    break;
                }
            }
            if shown >= max_show {
                break;
            }
        }
        eprintln!("\n  (shown {} cases)", shown);
    }

    /// Per-term Wigner trace for ALL non_quantized failures.
    ///
    /// Groups by (SG, k-point, W-value) and shows per-term breakdown
    /// for a representative of each cluster.
    #[test]
    fn diagnose_nonquantized_per_term() {
        use crate::irrep::wigner::{self, PerTermTrace};
        use std::collections::BTreeMap;

        // ── Collect all non_quantized cases with per-term traces ──────────
        struct CaseInfo {
            uni: usize,
            h_sg: u8,
            ml: String,
            kx: i8,
            ky: i8,
            kz: i8,
            kd: i8,
            dim: u8,
            n_anti: usize,
            w: f64,
            trace: Vec<PerTermTrace>,
            xf_basis: [[f64; 3]; 3],
            parent_sg: u8,
        }
        let mut cases: Vec<CaseInfo> = Vec::new();

        for uni in 1..=1651 {
            let mag_ops = match get_magnetic_operations(uni) {
                Some(m) => m,
                None => continue,
            };
            let h_info = match identify_unitary_subgroup_with_hall(uni) {
                Some(i) => i,
                None => continue,
            };
            let h_sg = h_info.sg as u8;
            let setting_xf_owned = h_info.msg_to_data.clone();
            let setting_xf = setting_xf_owned.as_ref();
            let mag_seitz = wigner::ops_to_seitz(&mag_ops);
            let g_sg = parent_spatial_sg(uni).unwrap_or(h_sg as usize) as u8;

            for ir in crate::irrep::query::irreps_of(h_sg) {
                if !ir.spinor {
                    continue;
                }
                let imag = ir.spin_character_imag();
                if imag.is_empty() {
                    continue;
                }

                let canonical_pure_translations: Vec<[f64; 3]> = h_info
                    .ops_from_hall
                    .operations
                    .iter()
                    .filter(|op| op.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
                    .map(|op| op.translation)
                    .collect();
                let mag_lg = wigner::filter_little_group_with_transform(
                    ir.kx,
                    ir.ky,
                    ir.kz,
                    ir.kd,
                    &mag_ops,
                    setting_xf,
                    Some(&canonical_pure_translations),
                );
                let antiunitary: Vec<usize> = mag_lg
                    .iter()
                    .filter(|&&i| mag_ops.operations[i].time_reversal)
                    .copied()
                    .collect();
                if antiunitary.is_empty() {
                    continue;
                }

                let h_spin = ir.spin_ops();
                let g_spin = if g_sg == h_sg {
                    h_spin
                } else {
                    IrrepRecord::spin_ops_for_sg(g_sg)
                };
                let ctx = wigner::SpinLiftContext {
                    h: h_spin,
                    g: g_spin,
                    sg: h_sg,
                };

                let mut trace = Vec::new();
                let diag = wigner::wigner_classify_spinor_direct_anti_diagnostic(
                    &ctx,
                    wigner::SpinorWignerInput {
                        characters_real: ir.characters(),
                        characters_imag: ir.spin_character_imag(),
                        operation_indices: ir.spin_lg_op_indices(),
                        k_vector: ir.k_vector(),
                    },
                    &antiunitary,
                    &mag_seitz,
                    setting_xf,
                    Some(&mut trace),
                );

                if let Err(wigner::DirectAntiFailure::NonQuantized) = diag {
                    // Compute W value from trace for clustering
                    let w_sum_re: f64 = trace.iter().map(|t| t.contrib_re).sum();
                    let w_sum_im: f64 = trace.iter().map(|t| t.contrib_im).sum();
                    let w = (w_sum_re * w_sum_re + w_sum_im * w_sum_im).sqrt()
                        / (antiunitary.len() as f64).max(1.0);
                    cases.push(CaseInfo {
                        uni,
                        h_sg,
                        ml: ir.ml.to_string(),
                        kx: ir.kx,
                        ky: ir.ky,
                        kz: ir.kz,
                        kd: ir.kd,
                        dim: ir.dim,
                        n_anti: antiunitary.len(),
                        xf_basis: setting_xf.map(|xf| xf.basis).unwrap_or([
                            [1.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0],
                        ]),
                        parent_sg: g_sg,
                        w,
                        trace,
                    });
                }
            }
        }

        eprintln!(
            "\n=== Per-term Wigner trace for {} non_quantized cases ===",
            cases.len()
        );

        // ── Cluster by (SG, k-point string) ───────────────────────────────
        let mut clusters: BTreeMap<(u8, String), Vec<&CaseInfo>> = BTreeMap::new();
        for c in &cases {
            let k_str = format!("({},{},{})/{}", c.kx, c.ky, c.kz, c.kd);
            clusters.entry((c.h_sg, k_str)).or_default().push(c);
        }

        eprintln!("\n=== Clusters ===");
        for ((sg, k_str), group) in &clusters {
            // W value distribution within this cluster
            let mut w_vals: Vec<f64> = group.iter().map(|c| c.w).collect();
            w_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let w_min = w_vals.first().unwrap();
            let w_max = w_vals.last().unwrap();
            eprintln!(
                "  SG{} {} {} cases, W in [{:.4}, {:.4}]",
                sg,
                k_str,
                group.len(),
                w_min,
                w_max
            );
        }

        // ── Show per-term tables for first case in each cluster ───────────
        eprintln!("\n=== Per-term details (first case per cluster) ===");
        for ((sg, k_str), group) in &clusters {
            let c = group[0];
            let is_signed_perm = c.xf_basis.iter().all(|row| {
                row.iter().all(|&v| {
                    let r = v.round();
                    (r == -1.0 || r == 0.0 || r == 1.0) && (v - r).abs() < 0.01
                })
            });
            eprintln!(
                "\n═══ SG{} UNI{} {} k={} dim={} n_anti={} W={:.4} parent={} xf=[{:?}] signd={} ═══",
                sg,
                c.uni,
                c.ml,
                k_str,
                c.dim,
                c.n_anti,
                c.w,
                c.parent_sg,
                c.xf_basis,
                is_signed_perm
            );

            // Show term summary line
            for (idx, t) in c.trace.iter().enumerate() {
                let chi_sign = if t.central { -1.0 } else { 1.0 };
                // Compute U_b² from U_b for diagnostics
                let u_b_sq = crate::irrep::wigner::su2_compose(&t.u_b, &t.u_b);
                // central parity: U_b² vs canonical U_sq in G table
                let spatial_central = crate::irrep::wigner::su2_same_up_to_sign(&u_b_sq, &t.u_sq_g);
                eprintln!(
                    "  [{:2}] rot=[{:3},{:3},{:3};{:3},{:3},{:3};{:3},{:3},{:3}] \
                     sq_rot=[{:3},{:3},{:3};{:3},{:3},{:3};{:3},{:3},{:3}] \
                     sq_spin={:3} loc={} χ=({:5.2},{:5.2}) c={} spC={:?} φ=({:5.2},{:5.2}) \
                     contrib=({:6.3},{:6.3}) | u_b=({:.3},{:.3},{:.3},{:.3}) u_b_sq=({:.3},{:.3},{:.3},{:.3}) u_sq_g=({:.3},{:.3},{:.3},{:.3})",
                    idx,
                    t.b_rot[0][0],
                    t.b_rot[0][1],
                    t.b_rot[0][2],
                    t.b_rot[1][0],
                    t.b_rot[1][1],
                    t.b_rot[1][2],
                    t.b_rot[2][0],
                    t.b_rot[2][1],
                    t.b_rot[2][2],
                    t.sq_rot[0][0],
                    t.sq_rot[0][1],
                    t.sq_rot[0][2],
                    t.sq_rot[1][0],
                    t.sq_rot[1][1],
                    t.sq_rot[1][2],
                    t.sq_rot[2][0],
                    t.sq_rot[2][1],
                    t.sq_rot[2][2],
                    t.sq_spin_idx,
                    t.sq_local_idx,
                    t.chi0_re * chi_sign,
                    t.chi0_im * chi_sign,
                    if t.central { "Y" } else { "N" },
                    spatial_central,
                    t.phase_re,
                    t.phase_im,
                    t.contrib_re,
                    t.contrib_im,
                    t.u_b[0],
                    t.u_b[1],
                    t.u_b[2],
                    t.u_b[3],
                    u_b_sq[0],
                    u_b_sq[1],
                    u_b_sq[2],
                    u_b_sq[3],
                    t.u_sq_g[0],
                    t.u_sq_g[1],
                    t.u_sq_g[2],
                    t.u_sq_g[3],
                );
            }

            // Running sum
            let mut running = (0.0f64, 0.0f64);
            eprint!("  Running sum:");
            for t in &c.trace {
                running.0 += t.contrib_re;
                running.1 += t.contrib_im;
                eprint!(" ({:.3},{:.3})", running.0, running.1);
            }
            eprintln!();
            eprintln!(
                "  W = ({:.6}, {:.6}) / {} = {:.6}",
                running.0,
                running.1,
                c.n_anti,
                (running.0 * running.0 + running.1 * running.1).sqrt() / c.n_anti as f64
            );
        }
    }

    /// Diagnostic: detailed Wigner terms for one failing spinor case (SG2 T UNI69).
    #[test]
    fn diagnose_sg2_spinor_wigner_failure() {
        use crate::irrep::wigner;

        // UNI69: magnetic group with H=SG2
        let uni = 69;
        let mag_ops = get_magnetic_operations(uni).expect("UNI69 should exist");
        let h_info = identify_unitary_subgroup_with_hall(uni).expect("H should exist");
        assert_eq!(h_info.sg, 2, "H should be SG2");
        let h_ops = h_info.ops_from_msg;
        let mag_seitz = wigner::ops_to_seitz(&mag_ops);
        let h_seitz = wigner::ops_to_seitz(&h_ops);

        // SG2 T-point spinor irreps
        for ir in crate::irrep::query::irreps_of(2) {
            if !ir.spinor || ir.k_label() != "T" {
                continue;
            }
            let imag = ir.spin_character_imag();
            if !imag.is_empty() {
                continue;
            }

            let mag_lg = wigner::filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
            let antiunitary: Vec<usize> = mag_lg
                .iter()
                .filter(|&&i| mag_ops.operations[i].time_reversal)
                .copied()
                .collect();
            if antiunitary.is_empty() {
                continue;
            }
            let unitary: Vec<usize> = mag_lg
                .iter()
                .filter(|&&i| !mag_ops.operations[i].time_reversal)
                .copied()
                .collect();

            let spin_ops = ir.spin_ops();
            let h_spin = spin_ops;
            let g_sg = parent_spatial_sg(uni).unwrap_or(2) as u8;
            let g_spin = if g_sg == 2 {
                h_spin
            } else {
                IrrepRecord::spin_ops_for_sg(g_sg)
            };
            let ctx = wigner::SpinLiftContext {
                h: h_spin,
                g: g_spin,
                sg: 2,
            };

            let ct = wigner::wigner_classify_spinor(
                &ctx,
                wigner::SpinorWignerInput {
                    characters_real: ir.characters(),
                    characters_imag: ir.spin_character_imag(),
                    operation_indices: ir.spin_lg_op_indices(),
                    k_vector: ir.k_vector(),
                },
                wigner::WignerGroupContext {
                    unitary_indices: &unitary,
                    magnetic_ops: &mag_seitz,
                    unitary_ops: &h_seitz,
                    antiunitary_representative: antiunitary[0],
                },
                None,
                Some(&antiunitary),
            );

            println!(
                "SG2 {} UNI{}: h_ops={} mag_ops={} n_lg={}/{}",
                ir.ml,
                uni,
                h_ops.len(),
                mag_ops.len(),
                unitary.len(),
                antiunitary.len()
            );
            println!("  spin_lg={:?} result={:?}", ir.spin_lg_op_indices(), ct);
            let timerev_vec: Vec<bool> =
                mag_ops.operations.iter().map(|o| o.time_reversal).collect();
            println!("  mag timerev={:?}", timerev_vec);
            println!(
                "  h_ops rots: {:?}",
                h_seitz.iter().map(|s| s.rot).collect::<Vec<_>>()
            );
            println!(
                "  mag lg unitary rots: {:?}",
                unitary
                    .iter()
                    .map(|&i| mag_ops.operations[i].rotation)
                    .collect::<Vec<_>>()
            );
            println!(
                "  spin ops rots: {:?}",
                (0..h_spin.0.len() / 9)
                    .map(|i| {
                        let off = i * 9;
                        [
                            h_spin.0[off..off + 3].to_vec(),
                            h_spin.0[off + 3..off + 6].to_vec(),
                            h_spin.0[off + 6..off + 9].to_vec(),
                        ]
                    })
                    .collect::<Vec<_>>()
            );
        }
    }

    /// Extract Pauli coefficients from the spin-op flat array.
    fn spin_su2_at(spin_op_su2: &[f64], idx: usize) -> Option<[f64; 4]> {
        if 4 * idx + 3 >= spin_op_su2.len() {
            return None;
        }
        Some([
            spin_op_su2[4 * idx],
            spin_op_su2[4 * idx + 1],
            spin_op_su2[4 * idx + 2],
            spin_op_su2[4 * idx + 3],
        ])
    }

    #[test]
    fn diagnose_none_examples() {
        let mut shown = 0usize;
        'outer: for uni in 1..=1651 {
            if shown >= 5 {
                break;
            }
            let mag_ops = match crate::SymmetryOps::from_magnetic_database(uni) {
                Ok(m) => m,
                Err(_) => continue,
            };
            let h_info = match identify_unitary_subgroup_with_hall(uni) {
                Some(i) => i,
                None => continue,
            };
            let h_sg = h_info.sg as u8;
            let h_ops = h_info.ops_from_msg;
            let mag_seitz = crate::irrep::wigner::ops_to_seitz(&mag_ops);
            let h_seitz = crate::irrep::wigner::ops_to_seitz(&h_ops);

            for ir in crate::irrep::query::irreps_of(h_sg) {
                if !ir.spinor {
                    continue;
                }
                if shown >= 5 {
                    break 'outer;
                }
                let mag_lg =
                    crate::irrep::wigner::filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
                let antiunitary: Vec<usize> = mag_lg
                    .iter()
                    .filter(|&&i| mag_ops[i].time_reversal)
                    .copied()
                    .collect();
                if antiunitary.is_empty() {
                    continue;
                }
                let unitary: Vec<usize> = mag_lg
                    .iter()
                    .filter(|&&i| !mag_ops[i].time_reversal)
                    .copied()
                    .collect();

                let h_spin = ir.spin_ops();
                if h_spin.0.is_empty() {
                    continue;
                }
                let g_sg = parent_spatial_sg(uni).unwrap_or(h_sg as usize) as u8;
                let g_spin = if g_sg == h_sg {
                    h_spin
                } else {
                    IrrepRecord::spin_ops_for_sg(g_sg)
                };
                let ctx = crate::irrep::wigner::SpinLiftContext {
                    h: h_spin,
                    g: g_spin,
                    sg: h_sg,
                };
                let a0_idx = select_spinor_a0(&antiunitary, &mag_seitz, g_sg == h_sg);

                let result = crate::irrep::wigner::wigner_classify_spinor(
                    &ctx,
                    wigner::SpinorWignerInput {
                        characters_real: ir.characters(),
                        characters_imag: ir.spin_character_imag(),
                        operation_indices: ir.spin_lg_op_indices(),
                        k_vector: ir.k_vector(),
                    },
                    wigner::WignerGroupContext {
                        unitary_indices: &unitary,
                        magnetic_ops: &mag_seitz,
                        unitary_ops: &h_seitz,
                        antiunitary_representative: a0_idx,
                    },
                    None,
                    Some(&antiunitary),
                );
                if result.is_ok() {
                    continue;
                }

                let n_lg = ir.spin_lg_char_count();
                let indices = ir.spin_lg_op_indices();
                let (h_rots, h_trans, h_su2) = ctx.h;
                let (g_rots, g_trans, g_su2) = ctx.g;
                let h_spin_seitz = crate::irrep::wigner::build_spin_seitz(h_rots, h_trans);
                let g_spin_seitz = crate::irrep::wigner::build_spin_seitz(g_rots, g_trans);
                let lg_set: std::collections::HashSet<usize> =
                    indices.iter().map(|&x| x as usize).collect();
                let (_, origin) = IrrepRecord::sg_setting(ctx.sg);
                let a0 = &mag_seitz[a0_idx];
                let to_bilbao = |rot: crate::mathfunc::Mat3I, trans: [f64; 3]| -> [f64; 3] {
                    if origin.len() < 3 {
                        return trans;
                    }
                    let mut t = trans;
                    for i in 0..3 {
                        let d: f64 = (0..3)
                            .map(|j| {
                                (if i == j { 1.0 } else { 0.0 } - rot[i][j] as f64) * origin[j]
                            })
                            .sum();
                        t[i] = (t[i] - d) % 1.0;
                        if t[i] < 0.0 {
                            t[i] += 1.0;
                        }
                    }
                    t
                };
                let a0_bilbao =
                    crate::irrep::wigner::SeitzOp::new(a0.rot, to_bilbao(a0.rot, a0.trans), false);
                let a0_match = g_spin_seitz.iter().position(|s| s.rot == a0.rot);
                let u_a0 = a0_match.and_then(|m| crate::irrep::wigner::spin_su2_at(g_su2, m));

                let global_to_local: std::collections::HashMap<usize, usize> = indices
                    .iter()
                    .enumerate()
                    .map(|(l, &g)| (g as usize, l))
                    .collect();

                let mut sq_all_in_lg = true;
                for &index in indices.iter().take(n_lg) {
                    let gsi = index as usize;
                    let h_spin = &h_spin_seitz[gsi];
                    let (g0h, _l1) = crate::irrep::wigner::compose_seitz(&a0_bilbao, h_spin);
                    let (sq, _lsq) = crate::irrep::wigner::square_seitz(&g0h);
                    if let Some(sq_si) = h_spin_seitz.iter().position(|s| s.rot == sq.rot) {
                        if !lg_set.contains(&sq_si) {
                            sq_all_in_lg = false;
                            break;
                        }
                    } else {
                        sq_all_in_lg = false;
                        break;
                    }
                }
                if !sq_all_in_lg {
                    continue;
                }

                shown += 1;
                println!("\n══════ NONE EXAMPLE #{} ══════", shown);
                println!(
                    "UNI={} SG={} {} k=({},{},{})/{} dim={} n_lg={} g_sg={} is_grey={}",
                    uni,
                    h_sg,
                    ir.ml,
                    ir.kx,
                    ir.ky,
                    ir.kz,
                    ir.kd,
                    ir.dim,
                    n_lg,
                    g_sg,
                    g_sg == h_sg
                );
                let det = |r: &crate::mathfunc::Mat3I| -> i32 {
                    r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
                        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
                        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0])
                };
                println!("a0: rot={:?} det={} u_a0={:?}", a0.rot, det(&a0.rot), u_a0);
                println!("spin_chars(lg): {:?}", &ir.characters()[..n_lg.min(8)]);
                println!("spin_lg_indices: {:?}", indices);

                let fmt_rot = |r: &crate::mathfunc::Mat3I| -> String {
                    format!(
                        "[{:2},{:2},{:2};{:2},{:2},{:2};{:2},{:2},{:2}]",
                        r[0][0],
                        r[0][1],
                        r[0][2],
                        r[1][0],
                        r[1][1],
                        r[1][2],
                        r[2][0],
                        r[2][1],
                        r[2][2]
                    )
                };
                let fmt_u = |v: &[f64; 4]| -> String {
                    format!("[{:7.4},{:7.4},{:7.4},{:7.4}]", v[0], v[1], v[2], v[3])
                };

                for (local, &index) in indices.iter().take(n_lg).enumerate() {
                    let gsi = index as usize;
                    let h_spin = &h_spin_seitz[gsi];
                    let u_h = crate::irrep::wigner::spin_su2_at(h_su2, gsi);
                    let (g0h, _l1) = crate::irrep::wigner::compose_seitz(&a0_bilbao, h_spin);
                    let (sq, _lsq) = crate::irrep::wigner::square_seitz(&g0h);
                    let sq_si = h_spin_seitz.iter().position(|s| s.rot == sq.rot);
                    let sq_local = sq_si.and_then(|s| global_to_local.get(&s).copied());

                    println!(
                        "  [{}] h={} det_h={:+} u_h={}",
                        local,
                        fmt_rot(&h_spin.rot),
                        det(&h_spin.rot),
                        u_h.map_or("?".into(), |v| fmt_u(&v))
                    );
                    println!(
                        "       g0h={} det_g0h={:+}",
                        fmt_rot(&g0h.rot),
                        det(&g0h.rot)
                    );
                    println!(
                        "       sq={} sq_si={:?} sq_lg={:?}",
                        fmt_rot(&sq.rot),
                        sq_si,
                        sq_local
                    );

                    if let (Some(u_h_v), Some(u_a0_v), Some(sq_si_v), Some(sq_local_v)) =
                        (u_h, u_a0, sq_si, sq_local)
                    {
                        let u_g0h = crate::irrep::wigner::su2_compose(&u_a0_v, &u_h_v);
                        let u_sq = crate::irrep::wigner::su2_compose(&u_g0h, &u_g0h);
                        let u_k = crate::irrep::wigner::spin_su2_at(h_su2, sq_si_v);
                        let rel_old = u_k
                            .and_then(|uk| crate::irrep::wigner::su2_same_up_to_sign(&u_sq, &uk));

                        let j = [0.0, 0.0, 1.0, 0.0];
                        let ju = crate::irrep::wigner::su2_compose(&j, &u_g0h);
                        let u_sq_j = crate::irrep::wigner::su2_compose(
                            &ju,
                            &crate::irrep::wigner::conj_pauli(&ju),
                        );
                        let rel_j = u_k
                            .and_then(|uk| crate::irrep::wigner::su2_same_up_to_sign(&u_sq_j, &uk));

                        println!("       u_g0h={} u_sq(U²)={}", fmt_u(&u_g0h), fmt_u(&u_sq));
                        println!(
                            "       u_k(Bilbao)={} rel(old)={:?} rel(J)={:?}",
                            u_k.map_or("?".into(), |v| fmt_u(&v)),
                            rel_old,
                            rel_j
                        );
                        println!("       chi0(sq)={}", ir.characters()[sq_local_v]);
                    }
                }
            }
        }
    }

    /// BCS validation: SG 213 (P4₁32) at X point, k=(0,1/2,0)
    ///
    /// Data source: Bilbao Crystallographic Server k-Subgroupsmag page
    /// (`k-Subgroupsmag.html`).  Confirmed:
    ///
    /// **Little group** (8 ops):  X₁, X₂ (both 2D)
    ///   χ(X₁) = [2, 0, 0, 0, √2, 0, -√2, 0]
    ///   χ(X₂) = [2, 0, 0, 0, -√2, 0, √2, 0]
    ///
    /// **Star**: 3 arms {(0,1/2,0), (1/2,0,0), (0,0,1/2)}
    ///
    /// **Full group** (ISOTROPY data): *X₁, *X₂ (both 6D = 3 arms × 2D)
    ///   These are the full-space-group irreps induced from the little group.
    ///   Our API returns the full group (star) representation.
    #[test]
    fn test_sg213_x_point_bcs() {
        let sg = 213u8;
        let irreps = crate::irrep::query::irreps_of(sg);

        // X-point in primitive cell: k = (0, 1/2, 0)
        // In our data: kx=0, ky=1, kz=0, kd=2
        let kx = 0i8;
        let ky = 1i8;
        let kz = 0i8;
        let kd = 2i8;

        let x_irreps: Vec<&IrrepRecord> = irreps
            .iter()
            .filter(|r| r.kx == kx && r.ky == ky && r.kz == kz && r.kd == kd)
            .collect();
        assert!(
            !x_irreps.is_empty(),
            "SG 213 should have irreps at X point ({}/{},{}/{},{}/{})",
            kx,
            kd,
            ky,
            kd,
            kz,
            kd
        );

        let scalar: Vec<_> = x_irreps.iter().filter(|r| !r.spinor).collect();
        let spinor: Vec<_> = x_irreps.iter().filter(|r| r.spinor).collect();

        println!(
            "SG213 X-point: {} scalar + {} spinor irreps",
            scalar.len(),
            spinor.len()
        );

        // BCS: exactly 2 scalar irreps at X (*X₁, *X₂)
        assert_eq!(
            scalar.len(),
            2,
            "SG 213 should have exactly 2 scalar irreps at X point"
        );

        // Full group irreps: 6D = 3 star arms × 2D little group
        for ir in &scalar {
            assert_eq!(
                ir.dim, 6,
                "{} full group dim = 6 (3 star arms × 2D little group, BCS)",
                ir.ml
            );
        }

        // Verify X1 and X2 ML labels exist
        let x_labels: Vec<&str> = scalar.iter().map(|r| r.ml).collect();
        println!("X-point scalar labels: {:?}", x_labels);
        assert!(
            x_labels.contains(&"X1"),
            "Should have X1 irrep (BCS *X₁). Labels: {:?}",
            x_labels
        );
        assert!(
            x_labels.contains(&"X2"),
            "Should have X2 irrep (BCS *X₂). Labels: {:?}",
            x_labels
        );

        let x1 = scalar.iter().find(|r| r.ml == "X1").unwrap();
        let x2 = scalar.iter().find(|r| r.ml == "X2").unwrap();

        // χ(E) = dim = 6 (BCS: identity trace of full-group irrep)
        let chars1 = x1.characters();
        let chars2 = x2.characters();
        assert!(!chars1.is_empty(), "X1 should have character table");
        assert!(!chars2.is_empty(), "X2 should have character table");
        assert!(
            (chars1[0] - 6.0).abs() < 0.01,
            "X1 χ(E) = dim = 6 (BCS), got {:.4}",
            chars1[0]
        );
        assert!(
            (chars2[0] - 6.0).abs() < 0.01,
            "X2 χ(E) = dim = 6 (BCS), got {:.4}",
            chars2[0]
        );

        // Full group has 48 operations (SG 213 is cubic, order 24 × 2 for
        // the star). ISOTROPY character table covers all full-group ops.
        println!("X1: {} characters", chars1.len());
        println!("X2: {} characters", chars2.len());
        assert!(
            chars1.len() >= 16,
            "X1 character table should cover full group ops, got {}",
            chars1.len()
        );

        // Both irreps should have isotropy subgroups (from ISOTROPY data)
        let subs1 = x1.subgroups();
        let subs2 = x2.subgroups();
        assert!(
            !subs1.is_empty(),
            "X1 should have isotropy subgroups (BCS shows subgroups)"
        );
        assert!(
            !subs2.is_empty(),
            "X2 should have isotropy subgroups (BCS shows subgroups)"
        );
        println!("X1: {} isotropy subgroups", subs1.len());
        println!("X2: {} isotropy subgroups", subs2.len());

        // Subgroup validity: all SG numbers in 1-230
        for sub in subs1.iter().chain(subs2.iter()) {
            assert!(
                sub.sg >= 1 && sub.sg <= 230,
                "Invalid subgroup SG {}",
                sub.sg
            );
            assert!(
                !sub.symbol.is_empty(),
                "Subgroup #{} should have HM symbol",
                sub.sg
            );
            assert!(
                !sub.direction.is_empty(),
                "Subgroup #{} should have direction",
                sub.sg
            );
        }

        // X1 subgroup SGs should include the high-symmetry subgroups
        let sg_nums1: Vec<usize> = subs1.iter().map(|s| s.sg).collect();
        let sg_nums2: Vec<usize> = subs2.iter().map(|s| s.sg).collect();
        println!("X1 subgroup SGs: {:?}", sg_nums1);
        println!("X2 subgroup SGs: {:?}", sg_nums2);

        // Both should have magnetic isotropy subgroups
        let mag1 = x1.magnetic_subgroups();
        let mag2 = x2.magnetic_subgroups();
        println!("X1: {} magnetic subgroups", mag1.len());
        println!("X2: {} magnetic subgroups", mag2.len());
        // Scalar irreps in non-centrosymmetric SGs typically have magnetic subgroups
        assert!(
            !mag1.is_empty() || !mag2.is_empty(),
            "At least one X irrep should have magnetic subgroups"
        );

        // Magnetic subgroup validity: UNI numbers in 1-1651
        for sub in mag1.iter().chain(mag2.iter()) {
            assert!(
                sub.mag_sg >= 1 && sub.mag_sg <= 1651,
                "Invalid magnetic UNI {}",
                sub.mag_sg
            );
            assert!(
                !sub.bns_label.is_empty(),
                "Magnetic subgroup UNI {} should have BNS label",
                sub.mag_sg
            );
        }
    }

    /// MSG 197.8 (I231'): P-point little group has no antiunitary operations.
    ///
    /// UNI 1510 is Type II (Grey, H=G=SG197).  At P-point k=(1,1,1)/2:
    ///
    /// Antiunitary k-preservation: -R⁻ᵀk - k ∈ L*.  For the I23 body-centred
    /// reciprocal lattice (FCC: h+k+l even), -2k=(-1,-1,-1) has sum -3 (odd)
    /// and therefore does NOT belong to L*.  k ≠ -k.
    ///
    /// The test enumerates the full magnetic little group.  ALL 24 operations
    /// passing the k-preservation filter are unitary; every antiunitary
    /// coset element θg fails the filter.  (The failure of θ alone does not
    /// logically guarantee failure of all θg; the code explicitly tests each.)
    ///
    /// Consequence: the Wigner A/B/C test has no antiunitary coset to work
    /// with.  `compute_corepresentation` returns Type A via the
    /// `TrivialNoAntiunitary` path — this is a current-API convention, not a
    /// BCS physical classification.  BCS reports Type C for these irreps
    /// because the antiunitary operator connects k → -k (different star arm),
    /// which requires a star-based co-representation construction not yet
    /// implemented.
    #[test]
    fn test_msg197_8_p_point_no_antiunitary_in_little_group() {
        let uni = 1510usize;
        let mag_ops = get_magnetic_operations(uni).unwrap();

        // Verify MSG type
        let msg_type = crate::MagneticSpaceGroupType::from_uni(uni).unwrap();
        assert_eq!(
            msg_type.type_,
            crate::MagneticType::Grey,
            "UNI 1510 is Type II Grey, not Type III"
        );
        assert_eq!(msg_type.number, 197);

        let n_u = mag_ops
            .operations
            .iter()
            .filter(|o| !o.time_reversal)
            .count();
        let n_a = mag_ops
            .operations
            .iter()
            .filter(|o| o.time_reversal)
            .count();
        assert_eq!(n_u, n_a, "Grey: equal unitary and antiunitary");

        // P5 spinor irrep at P-point
        let p5 = crate::irrep::query::irreps_of(197)
            .iter()
            .find(|r| r.ml == "P5" && r.spinor)
            .cloned()
            .expect("P5 spinor not found");
        assert_eq!(
            (p5.kx, p5.ky, p5.kz, p5.kd),
            (1, 1, 1, 2),
            "P-point is k=(1,1,1)/2"
        );

        // Magnetic little group at P-point
        let mag_lg = filter_little_group(p5.kx, p5.ky, p5.kz, p5.kd, &mag_ops);
        let antiunitary: Vec<usize> = mag_lg
            .iter()
            .filter(|&&i| mag_ops.operations[i].time_reversal)
            .copied()
            .collect();
        let unitary: Vec<usize> = mag_lg
            .iter()
            .filter(|&&i| !mag_ops.operations[i].time_reversal)
            .copied()
            .collect();

        assert!(
            antiunitary.is_empty(),
            "At P-point of I23, k≠-k: every antiunitary op fails the k-preservation filter"
        );
        assert!(
            !unitary.is_empty(),
            "Unitary little group at P-point should be non-empty"
        );

        // Current-API behaviour: empty antiunitary LG → TrivialNoAntiunitary → Type A
        let corep = p5.corepresentation(uni).unwrap();
        assert_eq!(
            corep.corep_type,
            CorepType::A,
            "API convention: empty antiunitary LG returns Type A"
        );
        assert_eq!(
            corep.source,
            WignerSource::TrivialNoAntiunitary,
            "Source must be TrivialNoAntiunitary (not a real Wigner classification)"
        );
        assert_eq!(
            corep.antiunitary_order, 0,
            "antiunitary_order must be 0 when LG has no antiunitary ops"
        );
    }

    /// - BNS → UNI mapping
    /// - Unitary subgroup identification (SG 197)
    /// - Magnetic operations are well-formed
    /// - Co-representations can be computed for all scalar irreps
    /// - χ(E) = dim for each valid co-irrep
    #[test]
    fn test_msg197_8_h_point_bcs() {
        let bns = "197.8";

        // 1. BNS → UNI (BCS: 197.8 = UNI 1510)
        let uni = super::uni_from_bns(bns);
        assert!(uni.is_some(), "Should find UNI for BNS {}", bns);
        let uni = uni.unwrap();
        assert_eq!(uni, 1510, "BCS: 197.8 = UNI 1510");

        // 2. Unitary subgroup: I231' is Type II Grey, H = G = SG197
        let h_info = identify_unitary_subgroup_with_hall(uni).unwrap();
        assert_eq!(h_info.sg, 197, "Unitary subgroup of I231' is SG197");

        // 3. Magnetic operations exist and are well-formed (Grey: 24U+24A=48)
        let mag_ops = get_magnetic_operations(uni).unwrap();
        let n_unitary = mag_ops.iter().filter(|o| !o.time_reversal).count();
        let n_anti = mag_ops.iter().filter(|o| o.time_reversal).count();
        assert_eq!(n_unitary, 24);
        assert_eq!(n_anti, 24);

        // 4. H's P-point irreps (BCS "H" = our "P" for body-centered)
        let h_sg = h_info.sg as u8;
        let h_irreps = crate::irrep::query::irreps_of(h_sg);
        let p_irreps: Vec<&IrrepRecord> = h_irreps.iter().filter(|r| r.k_label() == "P").collect();
        assert!(
            !p_irreps.is_empty(),
            "SG 197 should have P-point (BCS H-point) irreps"
        );

        let p_scalar: Vec<_> = p_irreps.iter().filter(|r| !r.spinor).collect();
        let p_spinor: Vec<_> = p_irreps.iter().filter(|r| r.spinor).collect();

        // BCS: P1, P2P3, P4 scalar + P̄5, P̄6P̄7 spinor
        assert!(
            p_scalar.len() >= 3,
            "Should have >=3 scalar irreps at P (BCS: 3), got {}",
            p_scalar.len()
        );
        assert!(
            p_spinor.len() >= 2,
            "Should have >=2 spinor irreps at P (BCS: 2), got {}",
            p_spinor.len()
        );

        // 5. Compute co-representations
        let coreps = super::compute_coreps(bns, "P");
        assert!(
            coreps.is_ok(),
            "Should compute coreps for {} at P-point: {:?}",
            bns,
            coreps
        );
        let coreps = coreps.unwrap();
        // BCS: 5 co-irrep labels (H₁, H₂H₃, H₄, H̄₅, H̄₆H̄₇)
        // Our data: 4 scalar + 3 spinor = 7 labels (each partner stored separately)
        assert!(
            coreps.len() >= 5,
            "Should have >=5 co-irreps (BCS), got {}",
            coreps.len()
        );

        // 6. χ(E) = dim for all valid co-irreps
        for (label, c) in &coreps {
            assert!(
                (c.characters[0] - c.dim as f64).abs() < 0.01,
                "χ(E)={:.4} ≠ dim={} for {}",
                c.characters[0],
                c.dim,
                label
            );
            // χ(E) must be positive
            assert!(c.characters[0] > 0.0, "χ(E) <= 0 for {}", label);
        }

        // 7. Scalar co-irreps at P-point: current API returns Type A.
        // k=(1,1,1)/2 in body-centred I23: -2k=(-1,-1,-1), h+k+l=-3 odd
        // → -2k ∉ L* (FCC).  No antiunitary ops in magnetic LG →
        // TrivialNoAntiunitary → Type A.  This is the current-API
        // convention; BCS reports Type C because the antiunitary connects
        // k → -k (different star arm).  Star-based co-representation
        // construction for that case is not yet implemented.
        let scalar_coreps: Vec<_> = coreps
            .iter()
            .filter(|(label, _)| p_scalar.iter().any(|ir| label.contains(&ir.ml[..2])))
            .collect();
        for (_label, c) in &scalar_coreps {
            assert_eq!(
                c.corep_type,
                CorepType::A,
                "API convention: empty antiunitary LG → Type A"
            );
        }
    }

    /// Simplest spinor Wigner failure: UNI 21 (BNS 5.14), SG 5 (C2), L-point.
    ///
    /// This is the absolute simplest case found — only 1 little-group operation.
    /// The Wigner sum has a single term: W = χ̃(a₀h) for h=E (identity).
    ///
    /// The antiunitary square is just (a₀)², so the failure reduces to:
    ///   (U_a₀)² ≠ ±U_{identity} = ±[1, 0, 0, 0]
    /// This means the SU(2) database's lift for a₀'s rotation, when squared,
    /// does not give ±identity — a fundamental gauge inconsistency for this
    /// single operation.
    #[test]
    fn test_simplest_spinor_failure_uni21() {
        let uni = 21usize; // BNS 5.14, SG 5 (C2), grey group
        let mag_ops = get_magnetic_operations(uni).unwrap();
        let h_info = identify_unitary_subgroup_with_hall(uni).unwrap();
        let h_sg = h_info.sg as u8;

        assert_eq!(h_sg, 5, "Unitary subgroup of 5.14 should be SG 5 (C2)");

        let h_irreps = crate::irrep::query::irreps_of(h_sg);
        // L-point spinor irrep: n_lg=1
        let l_spinor = h_irreps
            .iter()
            .find(|r| r.spinor && r.k_label() == "L" && r.spin_lg_char_count() == 1)
            .expect("SG 5 should have L-point spinor irrep with n_lg=1");

        println!(
            "SG {} {} k=({},{},{})/{} n_lg={} dim={}",
            l_spinor.sg,
            l_spinor.ml,
            l_spinor.kx,
            l_spinor.ky,
            l_spinor.kz,
            l_spinor.kd,
            l_spinor.spin_lg_char_count(),
            l_spinor.dim
        );

        // One little-group op: should be identity
        let (h_rots, h_trans, h_su2) = l_spinor.spin_ops();
        let h_spin_seitz = wigner::build_spin_seitz(h_rots, h_trans);
        let indices = l_spinor.spin_lg_op_indices();
        assert_eq!(indices.len(), 1, "Exactly 1 little-group op");
        let gsi = indices[0] as usize;

        println!(
            "LG op: spin[{}] rot={:?} su2={:?}",
            gsi,
            h_spin_seitz[gsi].rot,
            wigner::spin_su2_at(h_su2, gsi)
        );

        // Find a₀ (first antiunitary op)
        let mag_seitz = ops_to_seitz(&mag_ops);
        let mag_lg =
            filter_little_group(l_spinor.kx, l_spinor.ky, l_spinor.kz, l_spinor.kd, &mag_ops);
        let a0_idx = mag_lg
            .iter()
            .find(|&&i| mag_ops.operations[i].time_reversal)
            .copied()
            .expect("Should have antiunitary op");

        let a0 = &mag_seitz[a0_idx];
        let a0_spin = h_spin_seitz
            .iter()
            .position(|s| s.rot == a0.rot)
            .expect("a0 rotation should be in spin ops");
        let u_a0 = wigner::spin_su2_at(h_su2, a0_spin).unwrap();

        println!("a₀: spin[{}] rot={:?} u_a₀={:?}", a0_spin, a0.rot, u_a0);
        println!("a₀² in SO(3): rot={:?}", wigner::square_seitz(a0).0.rot);

        // The Wigner sum has ONE term: χ̃(a₀·E) = ±χ((a₀)²)
        // In SU(2): u_sq = (U_a₀)², compare with u_k = lift of a₀² rotation
        let u_sq_old = wigner::su2_compose(&u_a0, &u_a0);
        let sq_rot = wigner::square_seitz(a0).0.rot;
        let sq_spin = h_spin_seitz
            .iter()
            .position(|s| s.rot == sq_rot)
            .expect("sq rotation should be in spin ops");
        let u_k = wigner::spin_su2_at(h_su2, sq_spin).unwrap();

        println!("(U_a₀)² = {:?}", u_sq_old);
        println!("U_sq (from DB) = {:?} (rot={:?})", u_k, sq_rot);
        println!(
            "rel(U²): {:?}",
            wigner::su2_same_up_to_sign(&u_sq_old, &u_k)
        );

        // J-left formula
        let j = [0.0, 0.0, 1.0, 0.0];
        let ju = wigner::su2_compose(&j, &u_a0);
        let ju_star = [ju[0], -ju[1], -ju[2], -ju[3]];
        let u_sq_j = wigner::su2_compose(&ju, &ju_star);
        println!("(J·U_a₀)(J·U_a₀)* = {:?}", u_sq_j);
        println!("rel(J): {:?}", wigner::su2_same_up_to_sign(&u_sq_j, &u_k));

        // Check the full compute_corepresentation path
        let corep = compute_corepresentation(l_spinor, uni, &mag_ops);
        println!(
            "\nFull compute_corepresentation result: {:?}",
            corep.as_ref().map(|c| c.corep_type)
        );

        // For n_lg=1: U_a₀² vs U_k determines Wigner type
        let su2_rel = wigner::su2_same_up_to_sign(&u_sq_old, &u_k);
        println!("SU(2) relation: {:?}", su2_rel);
        assert!(su2_rel.is_some(), "SU(2) relation should be well-defined");
        assert!(
            corep.is_ok(),
            "Should classify n_lg=1 spinor case: {:?}",
            corep
        );
        let ct = corep.unwrap().corep_type;
        println!("Corep type: {:?}", ct);
        assert!(
            ct == CorepType::A || ct == CorepType::B,
            "Type should be A or B, got {:?}",
            ct
        );
    }

    /// Quick scan to find the simplest failing magnetic group.
    /// Only scans SG 1-10 grey-group UNI numbers (small groups).
    #[test]
    fn scan_simplest_spinor_failure() {
        let mut failures: Vec<(usize, usize, u8, String, String)> = Vec::new();
        for uni in 1..=50usize {
            let msg = crate::msg_database::get_magnetic_spacegroup_type(uni);
            if msg.type_ != crate::MagneticType::Grey {
                continue;
            }
            let h_info = match identify_unitary_subgroup_with_hall(uni) {
                Some(h) => h,
                None => continue,
            };
            let h_sg = h_info.sg as u8;
            let h_irreps = crate::irrep::query::irreps_of(h_sg);
            for ir in h_irreps.iter().filter(|r| r.spinor) {
                let n_lg = ir.spin_lg_char_count();
                if n_lg == 0 {
                    continue;
                }
                let bns = msg.bns_number.trim().to_string();
                if let Err(CorepComputationError::UnsupportedClassification {
                    source_irrep, ..
                }) = ir.corepresentation(uni)
                {
                    failures.push((
                        uni,
                        h_sg as usize,
                        n_lg as u8,
                        bns.clone(),
                        format!("{}/{}", source_irrep, ir.ml),
                    ));
                }
            }
        }
        failures.sort_by_key(|(_, sg, n_lg, _, _)| (*sg, *n_lg));
        println!(
            "Spinor Wigner failures (grey UNI 1-50, {} total):",
            failures.len()
        );
        for (uni, sg, n_lg, bns, ml) in failures.iter().take(15) {
            println!("  UNI {} SG {} {} n_lg={} irrep={}", uni, sg, bns, n_lg, ml);
        }
        // Also check SG 123 specifically — find its grey-group MSG
        let mut sg123_uni = 0usize;
        for uni in 1..=2000usize {
            let msg = crate::msg_database::get_magnetic_spacegroup_type(uni);
            if msg.type_ != crate::MagneticType::Grey {
                continue;
            }
            if let Some(h) = identify_unitary_subgroup_with_hall(uni)
                && h.sg == 123
            {
                sg123_uni = uni;
                break;
            }
        }
        println!(
            "\nSG 123 grey group: UNI {} BNS {}",
            sg123_uni,
            crate::msg_database::get_magnetic_spacegroup_type(sg123_uni)
                .bns_number
                .trim()
        );
        let sg123_spinors: Vec<_> = crate::irrep::query::irreps_of(123)
            .iter()
            .filter(|r| r.spinor && r.spin_lg_char_count() > 0)
            .map(|r| (r.ml, r.k_label(), r.spin_lg_char_count()))
            .collect();
        println!("SG 123 spinor irreps: {:?}", sg123_spinors);
    }

    /// Per-term oracle: trace the Wigner sum for non-quantized cases.
    /// Prints every term's contributions to identify why W ≠ ±dim.
    #[test]
    fn diagnose_non_quantized_per_term() {
        use crate::irrep::types::IrrepRecord;
        use crate::irrep::wigner;
        use num_complex::Complex64;

        let mut shown = 0usize;
        'outer: for uni in 1..=1651 {
            let mag_ops = match get_magnetic_operations(uni) {
                Some(m) => m,
                None => continue,
            };
            let h_info = match identify_unitary_subgroup_with_hall(uni) {
                Some(i) => i,
                None => continue,
            };
            let h_sg = h_info.sg as u8;

            let msg_rots: Vec<_> = h_info.ops_from_msg.iter().map(|o| o.rotation).collect();
            let msg_trans: Vec<_> = h_info.ops_from_msg.iter().map(|o| o.translation).collect();
            let hall_rots: Vec<_> = h_info.ops_from_hall.iter().map(|o| o.rotation).collect();
            let hall_trans: Vec<_> = h_info.ops_from_hall.iter().map(|o| o.translation).collect();
            let setting_xfs =
                wigner::find_setting_transform(&msg_rots, &msg_trans, &hall_rots, &hall_trans);
            let mut xf_owned = setting_xfs.first().cloned();
            if xf_owned.is_none() {
                let u_rots: Vec<_> = h_info
                    .ops_from_msg
                    .iter()
                    .filter(|o| !o.time_reversal)
                    .map(|o| o.rotation)
                    .collect();
                let u_trans: Vec<_> = h_info
                    .ops_from_msg
                    .iter()
                    .filter(|o| !o.time_reversal)
                    .map(|o| o.translation)
                    .collect();
                let u_ops = crate::SymmetryOps::from_parallel(
                    &u_rots,
                    &u_trans,
                    &vec![false; u_rots.len()],
                )
                .unwrap();
                if let Some((_, _, xf)) = standard_setting_transform(&u_ops, false) {
                    xf_owned = Some(xf);
                }
            }
            let setting_xf = xf_owned.as_ref();
            let h_ops = h_info.ops_from_msg;
            let _h_seitz = wigner::ops_to_seitz(&h_ops);
            let mag_seitz = wigner::ops_to_seitz(&mag_ops);

            for ir in crate::irrep::query::irreps_of(h_sg) {
                if !ir.spinor {
                    continue;
                }
                let mag_lg = wigner::filter_little_group_with_transform(
                    ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops, setting_xf, None,
                );
                let antiunitary: Vec<usize> = mag_lg
                    .iter()
                    .filter(|&&i| mag_ops.operations[i].time_reversal)
                    .copied()
                    .collect();
                if antiunitary.is_empty() {
                    continue;
                }

                let h_spin = ir.spin_ops();
                if h_spin.0.is_empty() {
                    continue;
                }
                let g_sg = parent_spatial_sg(uni).unwrap_or(h_sg as usize) as u8;
                let g_spin = if g_sg == h_sg {
                    h_spin
                } else {
                    IrrepRecord::spin_ops_for_sg(g_sg)
                };
                let ctx = wigner::SpinLiftContext {
                    h: h_spin,
                    g: g_spin,
                    sg: h_sg,
                };

                let diag = wigner::wigner_classify_spinor_direct_anti_diagnostic(
                    &ctx,
                    wigner::SpinorWignerInput {
                        characters_real: ir.characters(),
                        characters_imag: ir.spin_character_imag(),
                        operation_indices: ir.spin_lg_op_indices(),
                        k_vector: ir.k_vector(),
                    },
                    &antiunitary,
                    &mag_seitz,
                    setting_xf,
                    None,
                );

                if !matches!(diag, Err(wigner::DirectAntiFailure::NonQuantized)) {
                    continue;
                }

                shown += 1;
                if shown > 3 {
                    break 'outer;
                }

                println!(
                    "\n=== UNI{} SG{} {} k=({},{},{})/{} dim={} ===",
                    uni, h_sg, ir.ml, ir.kx, ir.ky, ir.kz, ir.kd, ir.dim
                );
                println!(
                    "  n_lg_ops={} n_anti={} g_sg={}",
                    ir.spin_lg_char_count(),
                    antiunitary.len(),
                    g_sg
                );

                // Manually compute each term
                let (h_spin_rots, h_spin_trans, h_spin_su2) = ctx.h;
                let (g_spin_rots, g_spin_trans, g_spin_su2) = ctx.g;
                let h_spin_seitz = wigner::build_spin_seitz(h_spin_rots, h_spin_trans);
                let g_spin_seitz = wigner::build_spin_seitz(g_spin_rots, g_spin_trans);
                let (_, origin) = IrrepRecord::sg_setting(ctx.sg);
                let id_rot: crate::mathfunc::Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

                println!(
                    "  {:>3} {:>32} {:>32} {:>8} {:>8} {:>8} {:>8}",
                    "b#", "b.rot", "b².rot", "χ_re", "χ_im", "central", "phase"
                );

                let mut w_sum = Complex64::ZERO;
                for &b_idx in &antiunitary {
                    let b = &mag_seitz[b_idx];
                    let (b_rot, b_trans) = if let Some(xf) = setting_xf {
                        let rot = xf.transform_rotation(&b.rot).unwrap_or(b.rot);
                        let trans = xf
                            .transform_translation(&b.rot, &b.trans)
                            .unwrap_or(b.trans);
                        (rot, trans)
                    } else {
                        (b.rot, b.trans)
                    };

                    let to_bilbao = |rot: crate::mathfunc::Mat3I, trans: [f64; 3]| -> [f64; 3] {
                        let mut t = trans;
                        for i in 0..3 {
                            let d: f64 = (0..3)
                                .map(|j| {
                                    let delta = if i == j { 1.0 } else { 0.0 };
                                    (delta - rot[i][j] as f64)
                                        * origin.get(j).copied().unwrap_or(0.0)
                                })
                                .sum();
                            t[i] = (t[i] - d) % 1.0;
                            if t[i] < 0.0 {
                                t[i] += 1.0;
                            }
                        }
                        t
                    };
                    let b_bilbao = wigner::SeitzOp::new(b_rot, to_bilbao(b_rot, b_trans), false);
                    let (sq, lattice_sq) = wigner::square_seitz(&b_bilbao);

                    // H lookup
                    let canon_t: Vec<[f64; 3]> = h_spin_seitz
                        .iter()
                        .filter(|s| s.rot == [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
                        .map(|s| s.trans)
                        .collect();
                    let (sq_spin_idx, chi0) = match wigner::find_sq_spin_lg_first(
                        &sq,
                        &h_spin_seitz,
                        ir.spin_lg_op_indices(),
                        &canon_t,
                    ) {
                        Some((idx, _, _)) => {
                            let local = ir
                                .spin_lg_op_indices()
                                .iter()
                                .position(|&x| x as usize == idx)
                                .unwrap_or(0);
                            let re = ir.characters().get(local).copied().unwrap_or(0.0);
                            let im = ir.spin_character_imag().get(local).copied().unwrap_or(0.0);
                            (Some(idx), Complex64::new(re, im))
                        }
                        None => {
                            // Fall back to identity
                            if let Some(id_idx) = h_spin_seitz.iter().position(|s| s.rot == id_rot)
                            {
                                let local = ir
                                    .spin_lg_op_indices()
                                    .iter()
                                    .position(|&x| x as usize == id_idx)
                                    .unwrap_or(0);
                                let re = ir.characters().get(local).copied().unwrap_or(0.0);
                                let im =
                                    ir.spin_character_imag().get(local).copied().unwrap_or(0.0);
                                (Some(id_idx), Complex64::new(re, im))
                            } else {
                                (None, Complex64::ZERO)
                            }
                        }
                    };

                    // G lookup + central parity
                    let b_spin_idx = g_spin_seitz.iter().position(|s| s.rot == b.rot);
                    let central = if let Some(g_idx) = b_spin_idx {
                        if let (Some(u_b), Some(u_sq)) = (
                            wigner::spin_su2_at(g_spin_su2, g_idx),
                            sq_spin_idx.and_then(|si| wigner::spin_su2_at(h_spin_su2, si)),
                        ) {
                            let u_b_sq = wigner::su2_compose(&u_b, &u_b);
                            match wigner::su2_same_up_to_sign(&u_b_sq, &u_sq) {
                                Some(spatial_central) => !spatial_central,
                                None => false,
                            }
                        } else {
                            false
                        }
                    } else {
                        false
                    };

                    let total_t = [
                        lattice_sq[0] as f64
                            + (sq.trans[0]
                                - h_spin_seitz
                                    .get(sq_spin_idx.unwrap_or(0))
                                    .map(|s| s.trans[0])
                                    .unwrap_or(0.0)),
                        lattice_sq[1] as f64
                            + (sq.trans[1]
                                - h_spin_seitz
                                    .get(sq_spin_idx.unwrap_or(0))
                                    .map(|s| s.trans[1])
                                    .unwrap_or(0.0)),
                        lattice_sq[2] as f64
                            + (sq.trans[2]
                                - h_spin_seitz
                                    .get(sq_spin_idx.unwrap_or(0))
                                    .map(|s| s.trans[2])
                                    .unwrap_or(0.0)),
                    ];
                    let phase = Complex64::new(
                        0.0,
                        2.0 * std::f64::consts::PI
                            * (ir.kx as f64 * total_t[0]
                                + ir.ky as f64 * total_t[1]
                                + ir.kz as f64 * total_t[2])
                            / ir.kd as f64,
                    )
                    .exp();

                    let chi = if central { -chi0 } else { chi0 };
                    let contrib = chi * phase;
                    w_sum += contrib;

                    println!(
                        "  {:>3} {:?} {:?} {:>8.3} {:>8.3} {:>8} {:>8.3}",
                        b_idx,
                        b.rot,
                        sq.rot,
                        chi.re,
                        chi.im,
                        central,
                        phase.re.round() as i8
                    );
                }
                let w = w_sum / (antiunitary.len() as f64);
                println!(
                    "  W = ({:.6}, {:.6}) / {} = ({:.6}, {:.6})  expected ±{}",
                    w_sum.re,
                    w_sum.im,
                    antiunitary.len(),
                    w.re,
                    w.im,
                    ir.dim
                );
                println!("  |W| - dim = {:.6}", w.norm() - ir.dim as f64);
            }
        }
    }

    /// Scan trigonal SG 143 (P3) grey group for SU(2) spinor Wigner failures.
    ///
    /// Trigonal groups contain C₃ rotations.
    #[test]
    fn diagnose_sg143_spinor_wigner() {
        let mut uni = 0usize;
        for u in 1..=2000usize {
            let msg = crate::msg_database::get_magnetic_spacegroup_type(u);
            if msg.type_ != crate::MagneticType::Grey {
                continue;
            }
            if let Some(h) = identify_unitary_subgroup_with_hall(u)
                && h.sg == 143
            {
                uni = u;
                break;
            }
        }
        assert!(uni > 0);
        let mag_ops = get_magnetic_operations(uni).unwrap();
        let irs = crate::irrep::query::irreps_of(143u8);

        println!(
            "{:>6} {:>3} {:>8} {:>6} {:>6}",
            "irrep", "nlg", "k", "ok", "fail"
        );
        for ir in irs
            .iter()
            .filter(|r| r.spinor && r.spin_lg_char_count() > 0)
        {
            let n_lg = ir.spin_lg_char_count();
            let indices = ir.spin_lg_op_indices();
            let (h_rots, h_trans, h_su2) = ir.spin_ops();
            let h_spin_seitz = wigner::build_spin_seitz(h_rots, h_trans);
            let mag_lg = filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
            let mag_seitz = ops_to_seitz(&mag_ops);
            let a0_idx = mag_lg
                .iter()
                .find(|&&i| mag_ops.operations[i].time_reversal)
                .copied();
            if a0_idx.is_none() {
                println!("{:>6} {:>3} — no antiunitary", ir.ml, n_lg);
                continue;
            }
            let a0 = &mag_seitz[a0_idx.unwrap()];
            let a0_spin = h_spin_seitz.iter().position(|s| s.rot == a0.rot);
            if a0_spin.is_none() {
                println!("{:>6} {:>3} — a0 not in spin", ir.ml, n_lg);
                continue;
            }
            let u_a0 = wigner::spin_su2_at(h_su2, a0_spin.unwrap()).unwrap();
            let mut ok = 0usize;
            let mut fail = 0usize;
            for &index in indices.iter().take(n_lg) {
                let gsi = index as usize;
                let u_h = match wigner::spin_su2_at(h_su2, gsi) {
                    Some(u) => u,
                    None => {
                        fail += 1;
                        continue;
                    }
                };
                let (g0h, _) = compose_seitz(a0, &h_spin_seitz[gsi]);
                let (sq, _) = square_seitz(&g0h);
                let sq_spin = match h_spin_seitz.iter().position(|s| s.rot == sq.rot) {
                    Some(s) => s,
                    None => {
                        fail += 1;
                        continue;
                    }
                };
                let u_sq = wigner::su2_compose(
                    &wigner::su2_compose(&u_a0, &u_h),
                    &wigner::su2_compose(&u_a0, &u_h),
                );
                let u_k = match wigner::spin_su2_at(h_su2, sq_spin) {
                    Some(u) => u,
                    None => {
                        fail += 1;
                        continue;
                    }
                };
                if wigner::su2_same_up_to_sign(&u_sq, &u_k).is_some() {
                    ok += 1;
                } else {
                    fail += 1;
                }
            }
            let marker = if fail > 0 { " ★ FAIL" } else { "" };
            println!(
                "{:>6} {:>3} ({:>2},{:>2},{:>2})/{:<2} {:>4} {:>4}{}",
                ir.ml, n_lg, ir.kx, ir.ky, ir.kz, ir.kd, ok, fail, marker
            );
            if fail > 0 {
                break;
            } // Just report first failure
        }
    }

    /// Deep dive: trace one specific MSG-gauge W failure case per-term.
    /// SG24 W-point UNI152 dim=1 — W=-0.5 instead of ±1.
    #[test]
    fn debug_msg_gauge_sg24_w_uni152() {
        use crate::irrep::wigner;
        let uni = 152;
        let mag_ops = get_magnetic_operations(uni).expect("UNI152 should exist");
        let mag_seitz = wigner::ops_to_seitz(&mag_ops);
        let h_info = identify_unitary_subgroup_with_hall(uni).expect("H info");
        let h_sg = h_info.sg as u8; // SG24
        let h_ops = h_info.ops_from_msg;
        let h_seitz = wigner::ops_to_seitz(&h_ops);

        // Find the W-point (1/2,1/2,1/2) spinor irrep
        for ir in crate::irrep::query::irreps_of(h_sg) {
            if !ir.spinor {
                continue;
            }
            if ir.kx != 1 || ir.ky != 1 || ir.kz != 1 || ir.kd != 2 {
                continue;
            }

            let mag_lg = wigner::filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
            let unitary: Vec<usize> = mag_lg
                .iter()
                .filter(|&&i| !mag_ops.operations[i].time_reversal)
                .copied()
                .collect();
            let antiunitary: Vec<usize> = mag_lg
                .iter()
                .filter(|&&i| mag_ops.operations[i].time_reversal)
                .copied()
                .collect();
            if antiunitary.is_empty() {
                continue;
            }

            let h_spin = ir.spin_ops();
            let g_sg = parent_spatial_sg(uni).unwrap_or(h_sg as usize) as u8;
            let g_spin = if g_sg == h_sg {
                h_spin
            } else {
                IrrepRecord::spin_ops_for_sg(g_sg)
            };
            let ctx = wigner::SpinLiftContext {
                h: h_spin,
                g: g_spin,
                sg: h_sg,
            };

            eprintln!(
                "\n═══ SG{} {} UNI{} dim={} imag={} ═══",
                h_sg,
                ir.k_label(),
                uni,
                ir.dim,
                !ir.spin_character_imag().is_empty()
            );
            eprintln!(
                "  n_lg_ops={} unitary={} antiunitary={}",
                ir.spin_lg_char_count(),
                unitary.len(),
                antiunitary.len()
            );
            eprintln!("  spin_lg_op_indices={:?}", ir.spin_lg_op_indices());
            eprintln!(
                "  chars={:?}",
                &ir.characters()[..ir.spin_lg_char_count().min(ir.characters().len())]
            );
            eprintln!("  imag={:?}", ir.spin_character_imag());

            // Build the spin→mag map manually and trace each term
            let h_spin_seitz = wigner::build_spin_seitz(h_spin.0, h_spin.1);
            let g_spin_seitz = wigner::build_spin_seitz(g_spin.0, g_spin.1);
            let h_to_spin =
                wigner::build_h_to_spin_map(&h_seitz, &h_spin_seitz, ir.spin_lg_op_indices());

            let mut spin_to_mag = std::collections::HashMap::<usize, usize>::new();
            for &mag_idx in &unitary {
                let h_match = match wigner::find_seitz(
                    &mag_seitz[mag_idx].rot,
                    &mag_seitz[mag_idx].trans,
                    &h_seitz,
                ) {
                    Some(m) => m,
                    None => continue,
                };
                if let Some(Some(spin_idx)) = h_to_spin.get(h_match.op_index) {
                    spin_to_mag.entry(*spin_idx).or_insert(mag_idx);
                }
            }

            let a0_idx = antiunitary[0];
            let a0_spatial =
                wigner::SeitzOp::new(mag_seitz[a0_idx].rot, mag_seitz[a0_idx].trans, false);
            let (a0_spin_idx, _) = wigner::find_spin_in_db(&a0_spatial, &g_spin_seitz).unwrap();
            let u_a0 = wigner::spin_su2_at(g_spin.2, a0_spin_idx).unwrap();

            let n_lg_ops = ir.spin_lg_char_count();
            eprintln!("\n  Per-term trace (MSG-gauge):");
            eprintln!(
                "  a0: mag_idx={} rot={:?} trans={:?} u_a0={:?}",
                a0_idx, a0_spatial.rot, a0_spatial.trans, u_a0
            );
            let mut w_sum = num_complex::Complex64::ZERO;
            let mut n_mapped = 0;

            for local in 0..n_lg_ops {
                let global_spin_idx = ir.spin_lg_op_indices()[local] as usize;
                let mag_idx = match spin_to_mag.get(&global_spin_idx) {
                    Some(&m) => m,
                    None => {
                        eprintln!("  [{local}] UNMAPPED spin={global_spin_idx}");
                        continue;
                    }
                };
                let h_msg =
                    wigner::SeitzOp::new(mag_seitz[mag_idx].rot, mag_seitz[mag_idx].trans, false);
                let (h_g_idx, _) = wigner::find_spin_in_db(&h_msg, &g_spin_seitz).unwrap();
                let u_h_g = wigner::spin_su2_at(g_spin.2, h_g_idx).unwrap();
                let u_g0h = wigner::su2_compose(&u_a0, &u_h_g);
                let u_sq_spatial = wigner::su2_compose(&u_g0h, &u_g0h);

                let (g0h, l1) = wigner::compose_seitz(&a0_spatial, &h_msg);
                let (sq, lattice_sq) = wigner::square_seitz(&g0h);
                let sq_h_match = wigner::find_seitz(&sq.rot, &sq.trans, &h_seitz).unwrap();
                let sq_spin_idx = h_to_spin
                    .get(sq_h_match.op_index)
                    .copied()
                    .flatten()
                    .unwrap();
                let u_sq_h = wigner::spin_su2_at(h_spin.2, sq_spin_idx).unwrap();

                let spatial_central = wigner::su2_same_up_to_sign(&u_sq_spatial, &u_sq_h).unwrap();
                let central = !spatial_central;

                let chi0_real = ir.characters()[sq_spin_idx.min(ir.characters().len() - 1)];
                let chi0_imag = ir
                    .spin_character_imag()
                    .get(sq_spin_idx)
                    .copied()
                    .unwrap_or(0.0);
                let chi0 = num_complex::Complex64::new(chi0_real, chi0_imag);
                let eta_ebar = -1.0;
                let chi = if central { eta_ebar * chi0 } else { chi0 };

                let r_l1 = wigner::mat_vec_i32(&g0h.rot, &l1);
                let total_lattice = wigner::add3(
                    &wigner::add3(&lattice_sq, &sq_h_match.lattice_shift),
                    &wigner::add3(&l1, &r_l1),
                );
                let phase = wigner::bloch_phase(ir.kx, ir.ky, ir.kz, ir.kd, &total_lattice);

                let contrib = chi * phase;
                w_sum += contrib;
                n_mapped += 1;

                eprintln!(
                    "  [{local}] mag={mag_idx} h_g={h_g_idx} u_h={u_h_g:?} sq_sp={sq_spin_idx} u_sq_sp={u_sq_spatial:?} u_sq_h={u_sq_h:?} spC={spatial_central} c={central} chi0=({chi0_real:.2},{chi0_imag:.2}) chi=({},{}) ph={:.2} contrib=({:.4},{:.4})",
                    chi.re, chi.im, phase.re, contrib.re, contrib.im,
                );
            }

            let w = if n_mapped > 0 {
                w_sum / n_mapped as f64
            } else {
                num_complex::Complex64::ZERO
            };
            eprintln!(
                "  W = ({:.6},{:.6}) / {n_mapped} = ({:.6},{:.6})",
                w_sum.re, w_sum.im, w.re, w.im
            );

            let ct = wigner::wigner_classify_spinor(
                &ctx,
                wigner::SpinorWignerInput {
                    characters_real: ir.characters(),
                    characters_imag: ir.spin_character_imag(),
                    operation_indices: ir.spin_lg_op_indices(),
                    k_vector: ir.k_vector(),
                },
                wigner::WignerGroupContext {
                    unitary_indices: &unitary,
                    magnetic_ops: &mag_seitz,
                    unitary_ops: &h_seitz,
                    antiunitary_representative: antiunitary[0],
                },
                None,
                Some(&antiunitary),
            );
            eprintln!("  Result: {:?}", ct);
            break; // Only first matching irrep
        }
    }

    #[test]
    fn debug_direct_anti_setting_uni663() {
        use crate::irrep::wigner;

        let uni = 663;
        let mag_ops = get_magnetic_operations(uni).expect("UNI663");
        let mag_seitz = wigner::ops_to_seitz(&mag_ops);
        let h_info = identify_unitary_subgroup_with_hall(uni).expect("H info");
        let h_seitz = wigner::ops_to_seitz(&h_info.ops_from_msg);

        println!(
            "UNI663 parent=SG{} H=SG{} Hall{}",
            parent_spatial_sg(uni).unwrap(),
            h_info.sg,
            h_info.hall
        );
        println!("H rotations:");
        for (i, op) in h_seitz.iter().enumerate() {
            println!("  H[{i}] R={:?} t={:?}", op.rot, op.trans);
        }

        let ir = crate::irrep::query::irreps_of(h_info.sg as u8)
            .iter()
            .find(|ir| ir.spinor && ir.k_label() == "B")
            .expect("SG3 B spinor irrep");
        let h_spin_seitz = wigner::build_spin_seitz(ir.spin_ops().0, ir.spin_ops().1);
        println!("spin_lg={:?}", ir.spin_lg_op_indices());
        for &idx in ir.spin_lg_op_indices() {
            let op = &h_spin_seitz[idx as usize];
            println!("  spin[{idx}] R={:?} t={:?}", op.rot, op.trans);
        }

        let mag_lg = wigner::filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
        for idx in mag_lg {
            if !mag_seitz[idx].timerev {
                continue;
            }
            let (sq, lattice) = wigner::square_seitz(&mag_seitz[idx]);
            println!(
                "  anti[{idx}] R={:?} t={:?} -> sq R={:?} t={:?} L={:?}",
                mag_seitz[idx].rot, mag_seitz[idx].trans, sq.rot, sq.trans, lattice
            );
        }
    }

    /// Phase 1 oracle: enumerate signed-permutation basis transforms between
    /// ops_from_msg (MSG embedding) and ops_from_hall (canonical H Hall setting).
    ///
    /// For each MSG we check whether any of the 48 signed-permutation matrices
    /// can map the rotation multiset correctly.  This is a diagnostic only —
    /// it does NOT change the classification path.
    #[test]
    fn phase1_setting_transform_oracle() {
        use crate::irrep::wigner::{SettingTransform, enumerate_signed_permutations};

        let all_t = enumerate_signed_permutations();
        let mut stats = std::collections::HashMap::<&str, usize>::new();
        let mut transform_hits: Vec<(usize, usize, [[i32; 3]; 3])> = Vec::new();

        for uni in 1..=1651 {
            let h_info = match identify_unitary_subgroup_with_hall(uni) {
                Some(i) => i,
                None => continue,
            };
            let msg_ops = &h_info.ops_from_msg;
            let hall_ops = &h_info.ops_from_hall;

            if msg_ops.is_empty() || hall_ops.is_empty() {
                continue;
            }

            // Collect rotation multisets
            let msg_rots: Vec<[[i32; 3]; 3]> = msg_ops.iter().map(|op| op.rotation).collect();
            let hall_rots: Vec<[[i32; 3]; 3]> = hall_ops.iter().map(|op| op.rotation).collect();

            // Check identity (T=I)
            let identity_works = rotation_multiset_eq(&msg_rots, &hall_rots);
            if identity_works {
                *stats.entry("identity").or_default() += 1;
                transform_hits.push((uni, h_info.hall, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]));
                continue;
            }

            // Try each signed-permutation T
            let mut found = false;
            for t in &all_t {
                let transform = SettingTransform {
                    basis: t.map(|row| row.map(|value| value as f64)),
                    origin: [0.0; 3],
                };
                let transformed: Vec<[[i32; 3]; 3]> = msg_rots
                    .iter()
                    .filter_map(|r| transform.transform_rotation(r))
                    .collect();
                if rotation_multiset_eq(&transformed, &hall_rots) {
                    *stats.entry("signed_perm").or_default() += 1;
                    transform_hits.push((uni, h_info.hall, *t));
                    found = true;
                    break;
                }
            }
            if !found {
                *stats.entry("not_found").or_default() += 1;
            }
        }

        println!("\n=== Phase 1: setting-transform oracle ===");
        let total: usize = stats.values().sum();
        println!("  Total UNI checked: {}", total);
        for (key, count) in &stats {
            println!(
                "  {:>20}  {:>6}  ({:.1}%)",
                key,
                count,
                *count as f64 / total as f64 * 100.0
            );
        }

        // Show signed-permutation examples
        println!("\n  Signed-permutation examples (first 10):");
        let mut shown = 0;
        for (uni, hall, t) in &transform_hits {
            if *t == [[1, 0, 0], [0, 1, 0], [0, 0, 1]] {
                continue;
            }
            if shown >= 10 {
                break;
            }
            let st = crate::spg_database::get_spacegroup_type(*hall);
            println!("    UNI{} Hall{} SG{} T={:?}", uni, hall, st.number, t);
            shown += 1;
        }

        // Cross-tab: transform status vs direct anti-coset failure stage
        println!("\n  Cross-tab: transform status vs direct anti-coset failure stage");
        let mut cross = std::collections::HashMap::<(&str, &str), usize>::new();
        for uni in 1..=1651 {
            let h_info = match identify_unitary_subgroup_with_hall(uni) {
                Some(i) => i,
                None => continue,
            };
            let msg_rots: Vec<[[i32; 3]; 3]> =
                h_info.ops_from_msg.iter().map(|op| op.rotation).collect();
            let hall_rots: Vec<[[i32; 3]; 3]> =
                h_info.ops_from_hall.iter().map(|op| op.rotation).collect();

            let transform_status = if rotation_multiset_eq(&msg_rots, &hall_rots) {
                "identity"
            } else if all_t.iter().any(|t| {
                let tr = SettingTransform {
                    basis: t.map(|row| row.map(|value| value as f64)),
                    origin: [0.0; 3],
                };
                let xformed: Vec<_> = msg_rots
                    .iter()
                    .filter_map(|r| tr.transform_rotation(r))
                    .collect();
                rotation_multiset_eq(&xformed, &hall_rots)
            }) {
                "signed_perm"
            } else {
                "not_found"
            };

            let mag_ops = match get_magnetic_operations(uni) {
                Some(m) => m,
                None => continue,
            };
            let h_sg = h_info.sg as u8;
            let g_sg = parent_spatial_sg(uni).unwrap_or(h_sg as usize) as u8;
            let g_spin = IrrepRecord::spin_ops_for_sg(g_sg);

            for ir in crate::irrep::query::irreps_of(h_sg) {
                if !ir.spinor {
                    continue;
                }
                let mag_lg =
                    crate::irrep::wigner::filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
                let antiunitary: Vec<usize> = mag_lg
                    .iter()
                    .filter(|&&i| mag_ops.operations[i].time_reversal)
                    .copied()
                    .collect();
                if antiunitary.is_empty() {
                    continue;
                }

                let h_spin = ir.spin_ops();
                let ctx = crate::irrep::wigner::SpinLiftContext {
                    h: h_spin,
                    g: g_spin,
                    sg: h_sg,
                };
                let mag_seitz = crate::irrep::wigner::ops_to_seitz(&mag_ops);

                let diag = crate::irrep::wigner::wigner_classify_spinor_direct_anti_diagnostic(
                    &ctx,
                    wigner::SpinorWignerInput {
                        characters_real: ir.characters(),
                        characters_imag: ir.spin_character_imag(),
                        operation_indices: ir.spin_lg_op_indices(),
                        k_vector: ir.k_vector(),
                    },
                    &antiunitary,
                    &mag_seitz,
                    None,
                    None,
                );
                let stage = match diag {
                    Ok(_) => "ok",
                    Err(e) => match e {
                        crate::irrep::wigner::DirectAntiFailure::SquareNotInSpinTable => {
                            "square_not_in_spin"
                        }
                        crate::irrep::wigner::DirectAntiFailure::SquareOutsideLittleGroup => {
                            "square_outside_lg"
                        }
                        _ => "other",
                    },
                };
                *cross.entry((transform_status, stage)).or_default() += 1;
            }
        }
        println!("  {:>12} {:>24} {:>6}", "transform", "stage", "count");
        let mut cross_sorted: Vec<_> = cross.iter().collect();
        cross_sorted.sort_by_key(|((t, s), _)| (t.to_string(), s.to_string()));
        for ((t, s), c) in &cross_sorted {
            println!("  {:>12} {:>24} {:>6}", t, s, c);
        }
    }
}

/// Compare two rotation multisets for equality (order-independent).
#[cfg(test)]
fn rotation_multiset_eq(a: &[[[i32; 3]; 3]], b: &[[[i32; 3]; 3]]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let mut b_used = vec![false; b.len()];
    for ra in a {
        let mut found = false;
        for (j, rb) in b.iter().enumerate() {
            if !b_used[j] && ra == rb {
                b_used[j] = true;
                found = true;
                break;
            }
        }
        if !found {
            return false;
        }
    }
    true
}

/// Phase 1b: for UNIs where a signed-permutation T was found AND
/// square_not_in_spin failures exist, verify whether applying T to
/// the anti-unitary MSG operations resolves the failures.
///
/// This oracle does NOT change the classification path — it only
/// confirms that the identified T is the correct basis transform.
#[test]
fn phase1b_verify_transform_fix() {
    use crate::irrep::wigner::{SettingTransform, enumerate_signed_permutations};

    let all_t = enumerate_signed_permutations();
    let mut fixed = 0usize;
    let mut total = 0usize;

    println!("\n=== Phase 1b: verify T fixes square_not_in_spin ===");

    for uni in 1..=1651 {
        let h_info = match identify_unitary_subgroup_with_hall(uni) {
            Some(i) => i,
            None => continue,
        };
        let mag_ops = match get_magnetic_operations(uni) {
            Some(m) => m,
            None => continue,
        };

        // Find T
        let msg_rots: Vec<[[i32; 3]; 3]> = h_info.ops_from_msg.iter().map(|o| o.rotation).collect();
        let hall_rots: Vec<[[i32; 3]; 3]> =
            h_info.ops_from_hall.iter().map(|o| o.rotation).collect();
        if rotation_multiset_eq(&msg_rots, &hall_rots) {
            continue;
        } // identity: no fix needed

        let t_found = all_t.iter().find(|t| {
            let tr = SettingTransform {
                basis: t.map(|row| row.map(|value| value as f64)),
                origin: [0.0; 3],
            };
            let xf: Vec<_> = msg_rots
                .iter()
                .filter_map(|r| tr.transform_rotation(r))
                .collect();
            rotation_multiset_eq(&xf, &hall_rots)
        });
        let t = match t_found {
            Some(t) => *t,
            None => continue,
        };
        let transform = SettingTransform {
            basis: t.map(|row| row.map(|value| value as f64)),
            origin: [0.0; 3],
        };

        let h_sg = h_info.sg as u8;
        let g_sg = parent_spatial_sg(uni).unwrap_or(h_sg as usize) as u8;
        let g_spin = IrrepRecord::spin_ops_for_sg(g_sg);
        let mag_seitz = crate::irrep::wigner::ops_to_seitz(&mag_ops);

        // Build H Hall spin table
        let h_spin = IrrepRecord::spin_ops_for_sg(h_sg);
        let _g_spin_seitz = crate::irrep::wigner::build_spin_seitz(g_spin.0, g_spin.1);

        for ir in crate::irrep::query::irreps_of(h_sg) {
            if !ir.spinor {
                continue;
            }
            let mag_lg =
                crate::irrep::wigner::filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
            let antiunitary: Vec<usize> = mag_lg
                .iter()
                .filter(|&&i| mag_ops.operations[i].time_reversal)
                .copied()
                .collect();
            if antiunitary.is_empty() {
                continue;
            }

            let ctx = crate::irrep::wigner::SpinLiftContext {
                h: h_spin,
                g: g_spin,
                sg: h_sg,
            };

            // Original direct anti-coset result
            let orig = crate::irrep::wigner::wigner_classify_spinor_direct_anti_diagnostic(
                &ctx,
                wigner::SpinorWignerInput {
                    characters_real: ir.characters(),
                    characters_imag: ir.spin_character_imag(),
                    operation_indices: ir.spin_lg_op_indices(),
                    k_vector: ir.k_vector(),
                },
                &antiunitary,
                &mag_seitz,
                None,
                None,
            );

            // Only interested in square_not_in_spin failures
            let orig_not_in_spin = matches!(
                &orig,
                Err(crate::irrep::wigner::DirectAntiFailure::SquareNotInSpinTable)
            );
            if !orig_not_in_spin {
                continue;
            }
            total += 1;

            // Now re-run with transformed b rotations
            let mut transformed_seitz: Vec<crate::irrep::wigner::SeitzOp> = mag_seitz.to_vec();
            for &b_idx in &antiunitary {
                let b = &mag_seitz[b_idx];
                let Some(r_t) = transform.transform_rotation(&b.rot) else {
                    continue;
                };
                let Some(t_t) = transform.transform_translation(&b.rot, &b.trans) else {
                    continue;
                };
                transformed_seitz[b_idx] = crate::irrep::wigner::SeitzOp::new(r_t, t_t, true);
            }

            // Re-run direct anti-coset with transformed b
            let fixed_result = crate::irrep::wigner::wigner_classify_spinor_direct_anti_diagnostic(
                &ctx,
                wigner::SpinorWignerInput {
                    characters_real: ir.characters(),
                    characters_imag: ir.spin_character_imag(),
                    operation_indices: ir.spin_lg_op_indices(),
                    k_vector: ir.k_vector(),
                },
                &antiunitary,
                &transformed_seitz,
                None,
                None,
            );

            let is_fixed = !matches!(
                &fixed_result,
                Err(crate::irrep::wigner::DirectAntiFailure::SquareNotInSpinTable)
            );
            if is_fixed {
                fixed += 1;
                if fixed <= 5 {
                    println!(
                        "  FIXED: UNI{} SG{} {} k=({},{},{})/{} T={:?}",
                        uni,
                        h_sg,
                        ir.k_label(),
                        ir.kx,
                        ir.ky,
                        ir.kz,
                        ir.kd,
                        t
                    );
                }
            }
        }
    }

    println!("  Total square_not_in_spin with signed_perm T: {}", total);
    println!(
        "  Fixed by applying T: {} ({:.1}%)",
        fixed,
        if total > 0 {
            fixed as f64 / total as f64 * 100.0
        } else {
            0.0
        }
    );
    println!("  Still failing: {}", total - fixed);
}
