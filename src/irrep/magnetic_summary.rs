//! Magnetic irrep summary — unified entry point for magnetic space group irreps.
//!
//! Given a magnetic space group (by UNI number, BNS label, or explicit operations),
//! this module returns a little-group corep summary of:
//!
//! 1. High-symmetry k-points with labels and fractional coordinates
//! 2. Magnetic co-representations (coreps) at the little-group level,
//!    classified by Wigner's test (star-based co-representations not yet implemented)
//! 3. Source H-irreps with Miller-Love / Bradley-Cracknell labels
//! 4. Isotropy subgroup candidates (ordinary and magnetic)
//!
//! # Example
//!
//! Query UNI 2 / BNS 1.2 (grey P1) and list its high-symmetry k-points
//! with their co-representations:
//!
//! ```
//! use cryspglib::irrep::magnetic_summary::*;
//!
//! let s = magnetic_irrep_summary_by_uni(2).unwrap();
//! println!("BNS {}  UNI={}  type={:?}  H=SG{}",
//!     s.bns_label, s.uni, s.magnetic_type, s.unitary_sg);
//!
//! for kp in &s.kpoints {
//!     let (kx, ky, kz, kd) = kp.coords;
//!     println!();
//!     println!("k-point {}  ({}/{}, {}/{}, {}/{})  |LG|={} ({}U+{}A)  coreps={}",
//!         kp.label, kx, kd, ky, kd, kz, kd,
//!         kp.little_group_order, kp.unitary_order, kp.antiunitary_order,
//!         kp.coreps.len());
//!
//!     for c in &kp.coreps {
//!         let srcs: Vec<&str> = c.source_irreps.iter().map(|s| s.ml).collect();
//!         let chi0 = c.characters.first().and_then(|value| *value)
//!             .map_or("N/A".to_string(), |value| format!("{}", value));
//!         println!("  {:20}  type={:?}  dim={}  χ(E)={}  src=[{}]",
//!             c.label, c.corep_type, c.dim, chi0, srcs.join(", "));
//!     }
//! }
//! ```

use std::collections::{BTreeSet, HashMap};

use num_complex::Complex64;

use crate::SymmetryOps;

// ── Error type ─────────────────────────────────────────────────────────────────

/// Errors that can occur during magnetic irrep summary computation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MagneticIrrepError {
    /// UNI number out of valid range (1–1651).
    InvalidUni(usize),
    /// BNS label not found in the magnetic space group database.
    UnknownBns(String),
    /// Could not retrieve magnetic symmetry operations for this UNI.
    MissingMagneticOperations(usize),
    /// Explicit operations do not match this UNI's standard BNS setting.
    OperationsInconsistentWithUni { uni: usize, reason: String },
    /// Could not identify the unitary subgroup H for this UNI.
    MissingUnitarySubgroup(usize),
    /// Magnetic operations could not be transformed to H's data-Hall frame.
    OperationSettingTransformFailed(usize),
    /// No irrep data available for this space group.
    MissingIrrepData { sg: u8 },
    /// Corep computation failed for a specific source H-irrep.
    CorepComputationFailed {
        uni: usize,
        sg: u8,
        k_label: String,
        source_irrep: String,
        reason: String,
    },
}

// ── Input type ─────────────────────────────────────────────────────────────────

/// Input specification for magnetic irrep summary.
pub enum MagneticIrrepInput<'a> {
    /// Look up by UNI number (1–1651).
    Uni(usize),
    /// Look up by BNS label (e.g. `"221.97"`).
    Bns(&'a str),
    /// Provide explicit magnetic symmetry operations with their UNI number.
    Operations { uni: usize, ops: &'a SymmetryOps },
}

// ── Summary output types ───────────────────────────────────────────────────────

/// Little-group corep summary for a magnetic space group.
///
/// This covers fixed-k little-group co-representations.  Star-based
/// (full Brillouin zone) co-representations are not yet implemented.
#[derive(Debug, Clone)]
pub struct MagneticIrrepSummary {
    /// UNI number (1–1651).
    pub uni: usize,
    /// BNS label (e.g. `"221.97"`).
    pub bns_label: String,
    /// Magnetic type: Grey, BlackWhite, Ordinary, or AntiTranslation.
    pub magnetic_type: crate::MagneticType,
    /// Parent spatial space group G ⊇ H.
    pub parent_sg: u8,
    /// Unitary subgroup H = M ∩ G.
    pub unitary_sg: u8,
    /// Hall number of the unitary subgroup.
    pub unitary_hall: usize,
    /// High-symmetry k-points with their coreps.
    pub kpoints: Vec<MagneticKPointSummary>,
}

/// A magnetic-irrep summary together with source corepresentations that could
/// not be constructed safely.
///
/// The ordinary [`magnetic_irrep_summary_by_uni`] API remains strict and
/// returns the first such failure.  This partial form is intended for band
/// analysis programs that can conservatively report an unresolved label for
/// only the affected band dimensions while retaining independently verified
/// labels elsewhere.
#[derive(Debug, Clone)]
pub struct PartialMagneticIrrepSummary {
    pub summary: MagneticIrrepSummary,
    pub unresolved_coreps: Vec<UnresolvedMagneticCorep>,
}

/// One source H-irrep whose magnetic corepresentation is unavailable.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnresolvedMagneticCorep {
    pub uni: usize,
    pub sg: u8,
    pub k_label: String,
    pub source_irrep: String,
    pub spinor: bool,
    /// Authoritative selected-arm source dimension when it remains available.
    /// Every corepresentation induced from this source has at least this
    /// dimension.
    pub minimum_dimension: Option<usize>,
    /// Wigner type when classification completed before a legacy output
    /// surface rejected the character row.
    pub classified_type: Option<crate::irrep::corep::CorepType>,
    /// Classification backend paired with [`Self::classified_type`].
    pub wigner_source: Option<crate::irrep::corep::WignerSource>,
    /// Full corepresentation dimension paired with a completed classification.
    pub classified_dimension: Option<usize>,
    pub reason: String,
}

/// Summary of magnetic corepresentations at a single k-point.
#[derive(Debug, Clone)]
pub struct MagneticKPointSummary {
    /// k-point label (e.g. `"GM"`, `"X"`, `"Z"`).
    pub label: String,
    /// Fractional reciprocal coordinates `(kx, ky, kz, denom)`.
    pub coords: (i8, i8, i8, i8),
    /// Total number of operations in the magnetic little group.
    pub little_group_order: usize,
    /// Number of unitary operations in the magnetic little group.
    pub unitary_order: usize,
    /// Number of anti-unitary operations in the magnetic little group.
    pub antiunitary_order: usize,
    /// Ordered magnetic little-group operations. Character entry `i` in every
    /// corep is the character of `operations[i]`.
    pub operations: Vec<MagneticLittleGroupOperation>,
    /// Ordinary Seitz conjugacy classes modulo lattice translations. The
    /// class formatter may refine these when projective/Bloch characters are
    /// not constant on a raw class.
    pub conjugacy_classes: Vec<MagneticConjugacyClass>,
    /// Magnetic co-representations at this k-point.
    pub coreps: Vec<MagneticCorepSummary>,
}

/// One column of a magnetic little-group character table.
#[derive(Debug, Clone, PartialEq)]
pub struct MagneticLittleGroupOperation {
    /// Zero-based character-table column.
    pub column: usize,
    /// Index in the full magnetic operation list supplied to the summary API.
    pub magnetic_operation_index: usize,
    /// Rotation in the ISOTROPY data-Hall frame used by the calculation.
    pub rotation: [[i32; 3]; 3],
    /// Fractional translation in the same data-Hall frame.
    pub translation: [f64; 3],
    /// Whether the operation is anti-unitary (contains time reversal).
    pub time_reversal: bool,
}

/// A conjugacy class of magnetic little-group operation columns.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MagneticConjugacyClass {
    /// Zero-based class number.
    pub class: usize,
    /// Representative operation column.
    pub representative: usize,
    /// Operation columns belonging to this class.
    pub members: Vec<usize>,
    /// Whether all members are anti-unitary (conjugation preserves this flag).
    pub time_reversal: bool,
}

/// Column layout for the formal character-table formatter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MagneticCharacterTableColumns {
    /// One column for every magnetic little-group operation.
    Operations,
    /// One column per character-compatible conjugacy class.
    ConjugacyClasses,
}

/// A single magnetic co-representation (corep).
#[derive(Debug, Clone)]
pub struct MagneticCorepSummary {
    /// Label for this corep (e.g. `"GM4-"` for single source, `"Z1Z4 + Z2Z3"` for Type-C pair).
    pub label: String,
    /// Source H-irreps that compose this corep.
    pub source_irreps: Vec<SourceIrrepSummary>,
    /// Wigner classification for successful coreps: A, B, or C.
    pub corep_type: crate::irrep::corep::CorepType,
    /// Which computational path produced the classification.
    pub source: crate::irrep::corep::WignerSource,
    /// Dimension of the magnetic irrep.
    pub dim: usize,
    /// Character χ̃(g) for each magnetic operation.
    ///
    /// `None` is used only for anti-unitary Type-A columns whose intertwiner
    /// has not been constructed. Defined zeroes remain `Some(0)` and
    /// unitary characters retain their full complex value.
    pub characters: Vec<Option<Complex64>>,
    /// Which operations are anti-unitary.
    pub timerev: Vec<bool>,
    /// Whether the character table is complete.
    pub completeness: crate::irrep::corep::CharacterCompleteness,
    /// Isotropy subgroup candidates (ordinary and magnetic).
    pub isotropy_candidates: Vec<CorepIsotropyCandidate>,
}

/// Summary of a source H-irrep (non-magnetic irrep of the unitary subgroup).
#[derive(Debug, Clone)]
pub struct SourceIrrepSummary {
    /// Space group number of H.
    pub sg: u8,
    /// Miller-Love label (e.g. `"GM4-"`, `"Z1Z4"`).
    pub ml: &'static str,
    /// Bradley-Cracknell label (e.g. `"\\Gamma_4^-"`).
    pub bc: &'static str,
    /// Selected-arm/little representation dimension, not the full-star image
    /// dimension stored in [`crate::irrep::types::IrrepRecord::dim`].
    pub dim: u8,
    /// Whether this is a spinor (double-valued) irrep.
    pub spinor: bool,
}

/// An isotropy subgroup candidate for a magnetic corep.
#[derive(Debug, Clone)]
pub struct CorepIsotropyCandidate {
    /// Source Miller-Love label from which this candidate originates.
    pub source_ml: &'static str,
    /// Ordinary (non-magnetic) isotropy subgroups.
    pub ordinary: Vec<crate::irrep::types::IsotropyRecord>,
    /// Magnetic isotropy subgroups.
    pub magnetic: Vec<crate::irrep::types::MagneticIsotropyRecord>,
    /// How this candidate relates to the source irrep(s).
    pub relation: IsotropyCandidateRelation,
}

/// How an isotropy candidate relates to the source irrep(s).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum IsotropyCandidateRelation {
    /// Directly from a single source irrep.
    DirectSourceIrrep,
    /// From a Type-C paired source (two H-irreps).
    TypeCPairedSource,
    /// From a compound irrep (CIR components).
    CompoundSource,
    /// Spinor irrep with no isotropy data available.
    SpinorNoIsotropyData,
}

fn summary_characters(
    corep: &crate::irrep::corep::ComplexCorepresentation,
    uni: usize,
    source_irrep: &str,
) -> Result<Vec<Option<Complex64>>, crate::irrep::corep::CorepComputationError> {
    if corep.characters.len() != corep.timerev.len() {
        return Err(
            crate::irrep::corep::CorepComputationError::UnsupportedClassification {
                uni,
                source_irrep: source_irrep.to_string(),
                reason: format!(
                    "complex character/timerev length mismatch: {} vs {}",
                    corep.characters.len(),
                    corep.timerev.len()
                ),
            },
        );
    }

    let pending_antiunitary = match corep.completeness {
        crate::irrep::corep::CharacterCompleteness::Complete => false,
        crate::irrep::corep::CharacterCompleteness::TypeAAntiunitaryPending { count } => {
            let actual = corep.timerev.iter().filter(|&&value| value).count();
            if count != actual {
                return Err(
                    crate::irrep::corep::CorepComputationError::UnsupportedClassification {
                        uni,
                        source_irrep: source_irrep.to_string(),
                        reason: format!(
                            "pending antiunitary count {count} disagrees with operation columns {actual}"
                        ),
                    },
                );
            }
            true
        }
    };

    corep
        .characters
        .iter()
        .copied()
        .zip(corep.timerev.iter().copied())
        .enumerate()
        .map(|(column, (value, time_reversal))| {
            if !value.re.is_finite() || !value.im.is_finite() {
                return Err(
                    crate::irrep::corep::CorepComputationError::UnsupportedClassification {
                        uni,
                        source_irrep: source_irrep.to_string(),
                        reason: format!("complex character column {column} is not finite: {value}"),
                    },
                );
            }
            Ok((!pending_antiunitary || !time_reversal).then_some(value))
        })
        .collect()
}

fn compound_coreps_for_summary(
    irrep: &crate::irrep::types::IrrepRecord,
    uni: usize,
    mag_ops: &SymmetryOps,
    operations: &[MagneticLittleGroupOperation],
) -> Result<Vec<MagneticCorepSummary>, crate::irrep::corep::CorepComputationError> {
    let branches =
        crate::irrep::corep::compute_compound_corepresentations_complex(irrep, uni, mag_ops)?;
    let mut summaries = Vec::with_capacity(branches.len());
    for branch in branches {
        let corep = branch.corep;
        let aligned = corep.characters.len() == operations.len()
            && corep.timerev.len() == operations.len()
            && corep.magnetic_operation_indices.len() == operations.len()
            && corep.operations.len() == operations.len()
            && corep
                .timerev
                .iter()
                .zip(operations)
                .all(|(time_reversal, operation)| *time_reversal == operation.time_reversal)
            && corep
                .magnetic_operation_indices
                .iter()
                .zip(operations)
                .all(|(index, operation)| *index == operation.magnetic_operation_index)
            && corep.operations.iter().zip(operations).all(
                |(corep_operation, summary_operation)| {
                    corep_operation.rotation == summary_operation.rotation
                        && corep_operation.translation == summary_operation.translation
                        && corep_operation.time_reversal == summary_operation.time_reversal
                },
            );
        if !aligned {
            return Err(
                crate::irrep::corep::CorepComputationError::UnsupportedClassification {
                    uni,
                    source_irrep: irrep.ml.to_string(),
                    reason: "compound plural character columns are not paired with the magnetic little-group operations"
                        .to_string(),
                },
            );
        }
        let characters = summary_characters(&corep, uni, irrep.ml)?;
        let source_irreps = branch
            .sources
            .iter()
            .map(|source| SourceIrrepSummary {
                sg: irrep.sg,
                ml: source.label,
                bc: irrep.bc,
                dim: u8::try_from(source.dimension).unwrap_or(u8::MAX),
                spinor: false,
            })
            .collect::<Vec<_>>();
        if source_irreps.iter().any(|source| source.dim == u8::MAX) {
            return Err(
                crate::irrep::corep::CorepComputationError::UnsupportedClassification {
                    uni,
                    source_irrep: irrep.ml.to_string(),
                    reason: "compound constituent dimension exceeds u8".to_string(),
                },
            );
        }
        let label = branch
            .sources
            .iter()
            .map(|source| match source.kind {
                crate::irrep::corep::CompoundCorepSourceKind::AuthoritativeCir => {
                    source.label.to_string()
                }
                crate::irrep::corep::CompoundCorepSourceKind::ConjugateRealification => {
                    format!("conj({})", source.label)
                }
                crate::irrep::corep::CompoundCorepSourceKind::DerivedAntiunitaryPartner => {
                    format!("a0({})", source.label)
                }
            })
            .collect::<Vec<_>>()
            .join(" + ");
        summaries.push(MagneticCorepSummary {
            label,
            source_irreps,
            corep_type: corep.corep_type,
            source: corep.source,
            dim: corep.dim,
            characters,
            timerev: corep.timerev,
            completeness: corep.completeness,
            isotropy_candidates: Vec::new(),
        });
    }
    Ok(summaries)
}

// ── Type-C dedup ───────────────────────────────────────────────────────────────

/// Deduplicate coreps by
/// `(corep_type, spinor-family, dim, rounded characters, timerev)`.
///
/// Only Type-C coreps with finite character values are eligible for merging.
/// Other valid coreps pass through without deduplication. Unsupported and
/// non-finite coreps should be rejected before this function is called.
fn dedup_coreps(coreps: Vec<MagneticCorepSummary>) -> Vec<MagneticCorepSummary> {
    let mut groups: Vec<Vec<MagneticCorepSummary>> = Vec::new();

    for c in coreps {
        // Only Type-C with finite characters can be deduplicated.
        let spinor_family = c.source_irreps.first().map(|source| source.spinor);
        let can_dedup = c.corep_type == crate::irrep::corep::CorepType::C
            && c.dim > 0
            && spinor_family.is_some()
            && c.source_irreps
                .iter()
                .all(|source| Some(source.spinor) == spinor_family)
            && c.characters.iter().all(|character| {
                character.is_none_or(|value| value.re.is_finite() && value.im.is_finite())
            });
        if !can_dedup {
            // Pass through — never merge non-C entries.
            groups.push(vec![c]);
            continue;
        }

        let key = (
            spinor_family.expect("dedup requires a uniform source family"),
            c.dim,
            round_chars(&c.characters),
            c.timerev.clone(),
        );
        let found = groups.iter_mut().find(|g| {
            let first = &g[0];
            first.corep_type == crate::irrep::corep::CorepType::C
                && first
                    .source_irreps
                    .iter()
                    .all(|source| source.spinor == key.0)
                && first.dim == key.1
                && round_chars(&first.characters) == key.2
                && first.timerev == key.3
        });
        match found {
            Some(group) => group.push(c),
            None => groups.push(vec![c]),
        }
    }

    groups
        .into_iter()
        .map(|mut group| {
            if group.len() == 1 {
                return group.remove(0);
            }
            // Type-C pair: merge source_irreps and update label.
            let mut merged = group.remove(0);
            let mut extra_sources: Vec<_> =
                group.into_iter().flat_map(|c| c.source_irreps).collect();
            merged.source_irreps.append(&mut extra_sources);
            // Several parent PIR records can expose the same reciprocal
            // constituent pair.  They are alternate provenance for one
            // Type-C corepresentation, not additional copies of its source
            // irreps.  Keep each typed source identity exactly once so that
            // the reported dimension remains the sum over the pair rather
            // than over duplicate parent records.
            merged
                .source_irreps
                .sort_by_key(|source| (source.sg, source.ml, source.bc, source.dim, source.spinor));
            merged.source_irreps.dedup_by_key(|source| {
                (source.sg, source.ml, source.bc, source.dim, source.spinor)
            });
            // Build combined label: sort source ML labels and join with " + ".
            let mut labels: Vec<&str> = merged.source_irreps.iter().map(|s| s.ml).collect();
            labels.sort();
            merged.label = labels.join(" + ");
            merged
        })
        .collect()
}

/// Round character values for dedup comparison while preserving undefined
/// Type-A anti-unitary columns and the complex plane.
fn round_chars(chars: &[Option<Complex64>]) -> Vec<Option<(i64, i64)>> {
    chars
        .iter()
        .map(|character| {
            character.map(|value| {
                (
                    (value.re * 1e8).round() as i64,
                    (value.im * 1e8).round() as i64,
                )
            })
        })
        .collect()
}

/// Return the authoritative dimension of an irrep on the selected k arm.
///
/// [`crate::irrep::types::IrrepRecord::dim`] is the dimension of the full-star/induced image and is
/// therefore not suitable for source dimensions in a little-group summary.
/// The typed character views carry the selected-arm dimension and validate the
/// generated data while constructing it.
fn selected_arm_dimension(
    irrep: &crate::irrep::types::IrrepRecord,
) -> Result<usize, crate::irrep::types::CharacterViewError> {
    if irrep.spinor {
        return Ok(irrep.spinor_selected_arm_view()?.dimension());
    }

    match irrep.ordinary_scalar_selected_arm_block_trace() {
        Ok(row) => Ok(row.dimension()),
        Err(crate::irrep::types::CharacterViewError::NotApplicable) => Ok(irrep
            .compound_selected_arm_view()?
            .block_trace()
            .dimension()),
        Err(error) => Err(error),
    }
}

// ── Isotropy candidates ────────────────────────────────────────────────────────

/// Attach isotropy subgroup candidates to each corep.
fn attach_isotropy_candidates(
    coreps: Vec<MagneticCorepSummary>,
    h_irreps: &[crate::irrep::types::IrrepRecord],
) -> Vec<MagneticCorepSummary> {
    coreps
        .into_iter()
        .map(|mut c| {
            let multi_source = c.source_irreps.len() > 1;
            let mut candidates: Vec<CorepIsotropyCandidate> = Vec::new();

            for src in &c.source_irreps {
                // Find the original IrrepRecord.
                let ir = h_irreps.iter().find(|r| r.sg == src.sg && r.ml == src.ml);

                let relation = if src.spinor {
                    IsotropyCandidateRelation::SpinorNoIsotropyData
                } else if multi_source {
                    IsotropyCandidateRelation::TypeCPairedSource
                } else if ir.is_some_and(|r| r.cir_component_count() > 0) {
                    IsotropyCandidateRelation::CompoundSource
                } else {
                    IsotropyCandidateRelation::DirectSourceIrrep
                };

                let (ordinary, magnetic) = match ir {
                    Some(rec) if !src.spinor => {
                        (rec.subgroups().to_vec(), rec.magnetic_subgroups().to_vec())
                    }
                    _ => (Vec::new(), Vec::new()),
                };

                // Always push a candidate so spinor sources are explicitly
                // marked SpinorNoIsotropyData rather than silently absent.
                candidates.push(CorepIsotropyCandidate {
                    source_ml: src.ml,
                    ordinary,
                    magnetic,
                    relation,
                });
            }

            // Dedup candidates by key.
            candidates = dedup_isotropy_candidates(candidates);
            c.isotropy_candidates = candidates;
            c
        })
        .collect()
}

/// Deduplicate isotropy candidates: keep the one with the richest data per source_ml.
fn dedup_isotropy_candidates(
    candidates: Vec<CorepIsotropyCandidate>,
) -> Vec<CorepIsotropyCandidate> {
    // Group by (source_ml, relation) and merge ordinary/magnetic sets.
    let mut seen: BTreeSet<(String, IsotropyCandidateRelation)> = BTreeSet::new();
    let mut result: Vec<CorepIsotropyCandidate> = Vec::new();

    for cand in candidates {
        let key = (cand.source_ml.to_string(), cand.relation);
        if seen.contains(&key) {
            // Merge into existing entry.
            if let Some(existing) = result
                .iter_mut()
                .find(|e| e.source_ml == cand.source_ml && e.relation == cand.relation)
            {
                existing.ordinary.extend(cand.ordinary);
                existing.magnetic.extend(cand.magnetic);
            }
        } else {
            seen.insert(key);
            result.push(cand);
        }
    }

    // Dedup ordinary subgroups by (sg, symbol, direction, domains, arms).
    for cand in &mut result {
        let mut ord_seen: BTreeSet<(usize, &str, &str, usize, usize)> = BTreeSet::new();
        cand.ordinary
            .retain(|s| ord_seen.insert((s.sg, s.symbol, s.direction, s.domains, s.arms)));

        // Dedup magnetic subgroups by (mag_sg, bns_label, direction).
        let mut mag_seen: BTreeSet<(usize, &str, &str)> = BTreeSet::new();
        cand.magnetic
            .retain(|s| mag_seen.insert((s.mag_sg, s.bns_label, s.direction)));
    }

    result
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct OperationKey {
    rotation: [[i32; 3]; 3],
    translation: [i64; 3],
    time_reversal: bool,
}

fn periodic_translation_key(value: f64) -> i64 {
    const SCALE: f64 = 1_000_000_000.0;
    let mut key = (value.rem_euclid(1.0) * SCALE).round() as i64;
    if key == SCALE as i64 {
        key = 0;
    }
    key
}

fn operation_key(operation: &crate::irrep::wigner::SeitzOp) -> OperationKey {
    OperationKey {
        rotation: operation.rot,
        translation: operation.trans.map(periodic_translation_key),
        time_reversal: operation.timerev,
    }
}

fn singleton_conjugacy_classes(
    operations: &[MagneticLittleGroupOperation],
) -> Vec<MagneticConjugacyClass> {
    operations
        .iter()
        .map(|operation| MagneticConjugacyClass {
            class: operation.column,
            representative: operation.column,
            members: vec![operation.column],
            time_reversal: operation.time_reversal,
        })
        .collect()
}

fn magnetic_conjugacy_classes(
    operations: &[MagneticLittleGroupOperation],
) -> Vec<MagneticConjugacyClass> {
    use crate::irrep::wigner::{SeitzOp, compose_seitz};

    let n = operations.len();
    if n == 0 {
        return Vec::new();
    }
    let seitz: Vec<_> = operations
        .iter()
        .map(|operation| {
            SeitzOp::new(
                operation.rotation,
                operation.translation,
                operation.time_reversal,
            )
        })
        .collect();
    let lookup: HashMap<_, _> = seitz
        .iter()
        .enumerate()
        .map(|(index, operation)| (operation_key(operation), index))
        .collect();
    if lookup.len() != n {
        return singleton_conjugacy_classes(operations);
    }

    let mut multiplication = vec![0usize; n * n];
    for left in 0..n {
        for right in 0..n {
            let (product, _) = compose_seitz(&seitz[left], &seitz[right]);
            let Some(&index) = lookup.get(&operation_key(&product)) else {
                // A numerically non-closed operation set must never be merged
                // into purported conjugacy classes.
                return singleton_conjugacy_classes(operations);
            };
            multiplication[left * n + right] = index;
        }
    }

    let identity_rotation = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let Some(identity) = seitz.iter().position(|operation| {
        operation.rot == identity_rotation
            && !operation.timerev
            && operation
                .trans
                .iter()
                .all(|value| periodic_translation_key(*value) == 0)
    }) else {
        return singleton_conjugacy_classes(operations);
    };
    let mut inverses = vec![usize::MAX; n];
    for element in 0..n {
        if let Some(inverse) =
            (0..n).find(|&candidate| multiplication[element * n + candidate] == identity)
        {
            inverses[element] = inverse;
        } else {
            return singleton_conjugacy_classes(operations);
        }
    }

    let mut visited = vec![false; n];
    let mut classes = Vec::new();
    for representative in 0..n {
        if visited[representative] {
            continue;
        }
        let mut members = BTreeSet::new();
        for conjugator in 0..n {
            let left = multiplication[conjugator * n + representative];
            let conjugate = multiplication[left * n + inverses[conjugator]];
            members.insert(conjugate);
        }
        let members: Vec<_> = members.into_iter().collect();
        for &member in &members {
            visited[member] = true;
        }
        classes.push(MagneticConjugacyClass {
            class: classes.len(),
            representative: members[0],
            time_reversal: operations[members[0]].time_reversal,
            members,
        });
    }
    classes
}

// ── Entry points ───────────────────────────────────────────────────────────────

/// Compute magnetic irrep summary from any input type.
pub fn magnetic_irrep_summary(
    input: MagneticIrrepInput,
) -> Result<MagneticIrrepSummary, MagneticIrrepError> {
    match input {
        MagneticIrrepInput::Uni(uni) => magnetic_irrep_summary_by_uni(uni),
        MagneticIrrepInput::Bns(bns) => magnetic_irrep_summary_by_bns(bns),
        MagneticIrrepInput::Operations { uni, ops } => magnetic_irrep_summary_from_ops(uni, ops),
    }
}

/// Compute magnetic irrep summary from a UNI number (1–1651).
pub fn magnetic_irrep_summary_by_uni(
    uni: usize,
) -> Result<MagneticIrrepSummary, MagneticIrrepError> {
    if uni == 0 || uni > 1651 {
        return Err(MagneticIrrepError::InvalidUni(uni));
    }
    let mag_ops = SymmetryOps::from_magnetic_database(uni)
        .map_err(|_| MagneticIrrepError::MissingMagneticOperations(uni))?;
    magnetic_irrep_summary_from_ops(uni, &mag_ops)
}

/// Compute every safely available magnetic corepresentation for a UNI number.
///
/// Unlike [`magnetic_irrep_summary_by_uni`], a source-specific
/// [`MagneticIrrepError::CorepComputationFailed`] does not discard unrelated
/// k-points and corepresentations.  Each omitted source is returned in
/// [`PartialMagneticIrrepSummary::unresolved_coreps`]; structural/database
/// errors remain fatal.
pub fn magnetic_irrep_summary_by_uni_partial(
    uni: usize,
) -> Result<PartialMagneticIrrepSummary, MagneticIrrepError> {
    if uni == 0 || uni > 1651 {
        return Err(MagneticIrrepError::InvalidUni(uni));
    }
    let mag_ops = SymmetryOps::from_magnetic_database(uni)
        .map_err(|_| MagneticIrrepError::MissingMagneticOperations(uni))?;
    magnetic_irrep_summary_from_ops_partial(uni, &mag_ops)
}

/// Compute magnetic irrep summary from a BNS label (e.g. `"221.97"`).
pub fn magnetic_irrep_summary_by_bns(
    bns: &str,
) -> Result<MagneticIrrepSummary, MagneticIrrepError> {
    let uni = crate::irrep::corep::uni_from_bns(bns)
        .ok_or_else(|| MagneticIrrepError::UnknownBns(bns.to_string()))?;
    magnetic_irrep_summary_by_uni(uni)
}

/// Compute magnetic irrep summary using explicit magnetic operations.
///
/// `ops` must be the complete operation set in the UNI database's first-Hall
/// BNS setting. Operation order and integer lattice shifts of translations are
/// ignored, while rotations and time-reversal flags must match exactly.
/// Setting-equivalent operation sets in another Hall frame are rejected because
/// the corepresentation pipeline consumes first-Hall coordinates. This strict
/// frame-specific contract also distinguishes database first-Hall sets that an
/// operation-only classifier may find equivalent after a setting transform.
pub fn magnetic_irrep_summary_from_ops(
    uni: usize,
    mag_ops: &SymmetryOps,
) -> Result<MagneticIrrepSummary, MagneticIrrepError> {
    Ok(magnetic_irrep_summary_from_ops_impl(uni, mag_ops, false)?.summary)
}

/// Partial counterpart of [`magnetic_irrep_summary_from_ops`].
pub fn magnetic_irrep_summary_from_ops_partial(
    uni: usize,
    mag_ops: &SymmetryOps,
) -> Result<PartialMagneticIrrepSummary, MagneticIrrepError> {
    magnetic_irrep_summary_from_ops_impl(uni, mag_ops, true)
}

fn magnetic_irrep_summary_from_ops_impl(
    uni: usize,
    mag_ops: &SymmetryOps,
    retain_partial: bool,
) -> Result<PartialMagneticIrrepSummary, MagneticIrrepError> {
    validate_operations_for_uni(uni, mag_ops)?;

    // 1. Identify H (unitary subgroup) with Hall setting information.
    let h_info = crate::irrep::corep::identify_unitary_subgroup_with_hall(uni)
        .ok_or(MagneticIrrepError::MissingUnitarySubgroup(uni))?;

    // 2. Get MSG metadata.
    let msg = crate::MagneticSpaceGroupType::from_uni(uni)
        .map_err(|_| MagneticIrrepError::InvalidUni(uni))?;

    // 3. Get k-points from H's irrep data.
    let h_kpoints = crate::irrep::query::kpoints_of(h_info.sg as u8);
    if h_kpoints.is_empty() {
        return Err(MagneticIrrepError::MissingIrrepData {
            sg: h_info.sg as u8,
        });
    }

    // 4. Canonical pure translations from H Hall setting.
    //    Needed for centered groups (F/I/C/A) where MSG-derived translations
    //    are only a subset of the full centering translation subgroup.
    let canonical_translations: Vec<[f64; 3]> = h_info
        .ops_from_hall
        .operations
        .iter()
        .filter(|op| op.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        .map(|op| op.translation)
        .collect();
    let setting_xf = h_info.msg_to_data.as_ref();
    let mag_ops_data = crate::irrep::corep::operations_in_data_hall_frame(mag_ops, setting_xf)
        .ok_or(MagneticIrrepError::OperationSettingTransformFailed(uni))?;

    // 5. Get H's irreps for corep computation.
    let h_irreps = crate::irrep::query::irreps_of(h_info.sg as u8);

    // 6. Build k-point summaries with little group metadata and coreps.
    let mut unresolved_coreps = Vec::new();
    let kpoints: Result<Vec<MagneticKPointSummary>, MagneticIrrepError> = h_kpoints
        .into_iter()
        .map(|kp| {
            let (kx, ky, kz, kd) = kp.coords;
            let mag_lg = crate::irrep::wigner::filter_little_group_with_transform(
                kx,
                ky,
                kz,
                kd,
                &mag_ops_data,
                None,
                Some(&canonical_translations),
            );
            let unitary_order = mag_lg
                .iter()
                .filter(|&&i| !mag_ops_data.operations[i].time_reversal)
                .count();
            let antiunitary_order = mag_lg
                .iter()
                .filter(|&&i| mag_ops_data.operations[i].time_reversal)
                .count();
            let operations: Vec<_> = mag_lg
                .iter()
                .enumerate()
                .map(|(column, &magnetic_operation_index)| {
                    let operation = &mag_ops_data.operations[magnetic_operation_index];
                    MagneticLittleGroupOperation {
                        column,
                        magnetic_operation_index,
                        rotation: operation.rotation,
                        translation: operation.translation,
                        time_reversal: operation.time_reversal,
                    }
                })
                .collect();
            let conjugacy_classes = magnetic_conjugacy_classes(&operations);

            // Compute coreps for each irrep at this k-point.
            let mut raw_coreps: Vec<MagneticCorepSummary> = Vec::new();
            for &idx in &kp.irreps {
                let ir = &h_irreps[idx];
                if ir.cir_component_count() > 0 {
                    match compound_coreps_for_summary(ir, uni, mag_ops, &operations) {
                        Ok(coreps) => raw_coreps.extend(coreps),
                        Err(err) => {
                            let failure = UnresolvedMagneticCorep {
                                uni,
                                sg: ir.sg,
                                k_label: kp.label.clone(),
                                source_irrep: ir.ml.to_string(),
                                spinor: false,
                                minimum_dimension: selected_arm_dimension(ir).ok(),
                                classified_type: None,
                                wigner_source: Some(
                                    crate::irrep::corep::WignerSource::ScalarCIR,
                                ),
                                classified_dimension: None,
                                reason: err.to_string(),
                            };
                            if retain_partial {
                                unresolved_coreps.push(failure);
                            } else {
                                return Err(MagneticIrrepError::CorepComputationFailed {
                                    uni: failure.uni,
                                    sg: failure.sg,
                                    k_label: failure.k_label,
                                    source_irrep: failure.source_irrep,
                                    reason: failure.reason,
                                });
                            }
                        }
                    }
                    continue;
                }
                match crate::irrep::corep::compute_corepresentation_complex(ir, uni, mag_ops) {
                    Ok(c) => {
                        let aligned = c.characters.len() == operations.len()
                            && c.timerev.len() == operations.len()
                            && c.magnetic_operation_indices.len() == operations.len()
                            && c.operations.len() == operations.len()
                            && c.timerev.iter().zip(&operations).all(
                                |(time_reversal, operation)| {
                                    *time_reversal == operation.time_reversal
                                },
                            )
                            && c
                                .magnetic_operation_indices
                                .iter()
                                .zip(&operations)
                                .all(|(index, operation)| {
                                    *index == operation.magnetic_operation_index
                                })
                            && c.operations.iter().zip(&operations).all(
                                |(corep_operation, summary_operation)| {
                                    corep_operation.rotation == summary_operation.rotation
                                        && corep_operation.translation
                                            == summary_operation.translation
                                        && corep_operation.time_reversal
                                            == summary_operation.time_reversal
                                },
                            );
                        if !aligned {
                            return Err(MagneticIrrepError::CorepComputationFailed {
                                uni,
                                sg: h_info.sg as u8,
                                k_label: kp.label.clone(),
                                source_irrep: ir.ml.to_string(),
                                reason: "corepresentation character columns are not paired with the constructed magnetic little-group operations"
                                    .to_string(),
                            });
                        }
                        let selected_dim = selected_arm_dimension(ir).map_err(|error| {
                            MagneticIrrepError::CorepComputationFailed {
                                uni,
                                sg: h_info.sg as u8,
                                k_label: kp.label.clone(),
                                source_irrep: ir.ml.to_string(),
                                reason: format!(
                                    "selected-arm source dimension lookup failed: {error}"
                                ),
                            }
                        })?;
                        let selected_dim_u8 = u8::try_from(selected_dim).map_err(|_| {
                            MagneticIrrepError::CorepComputationFailed {
                                uni,
                                sg: h_info.sg as u8,
                                k_label: kp.label.clone(),
                                source_irrep: ir.ml.to_string(),
                                reason: format!(
                                    "selected-arm source dimension {selected_dim} exceeds u8"
                                ),
                            }
                        })?;
                        let expected_dim = match c.corep_type {
                            crate::irrep::corep::CorepType::A => selected_dim,
                            crate::irrep::corep::CorepType::B
                            | crate::irrep::corep::CorepType::C => {
                                selected_dim.checked_mul(2).ok_or_else(|| {
                                    MagneticIrrepError::CorepComputationFailed {
                                        uni,
                                        sg: h_info.sg as u8,
                                        k_label: kp.label.clone(),
                                        source_irrep: ir.ml.to_string(),
                                        reason: format!(
                                            "raw {:?} corepresentation dimension overflow for selected-arm source dimension {selected_dim}",
                                            c.corep_type
                                        ),
                                    }
                                })?
                            }
                        };
                        if c.dim != expected_dim {
                            return Err(MagneticIrrepError::CorepComputationFailed {
                                uni,
                                sg: h_info.sg as u8,
                                k_label: kp.label.clone(),
                                source_irrep: ir.ml.to_string(),
                                reason: format!(
                                    "raw {:?} corepresentation dimension {} disagrees with selected-arm source dimension {} (expected {})",
                                    c.corep_type, c.dim, selected_dim, expected_dim
                                ),
                            });
                        }
                        let characters = summary_characters(&c, uni, ir.ml).map_err(|error| {
                            MagneticIrrepError::CorepComputationFailed {
                                uni,
                                sg: h_info.sg as u8,
                                k_label: kp.label.clone(),
                                source_irrep: ir.ml.to_string(),
                                reason: error.to_string(),
                            }
                        })?;
                        raw_coreps.push(MagneticCorepSummary {
                            label: ir.ml.to_string(),
                            source_irreps: vec![SourceIrrepSummary {
                                sg: ir.sg,
                                ml: ir.ml,
                                bc: ir.bc,
                                dim: selected_dim_u8,
                                spinor: ir.spinor,
                            }],
                            corep_type: c.corep_type,
                            source: c.source,
                            dim: c.dim,
                            characters,
                            timerev: c.timerev,
                            completeness: c.completeness,
                            isotropy_candidates: Vec::new(),
                        });
                    }
                    Err(err) => {
                        let (classified_type, wigner_source, classified_dimension) = match &err {
                            crate::irrep::corep::CorepComputationError::ComplexUnitaryCharacters {
                                corep_type,
                                source,
                                dimension,
                                ..
                            } => (Some(*corep_type), Some(*source), Some(*dimension)),
                            _ => (None, None, None),
                        };
                        let failure = UnresolvedMagneticCorep {
                            uni,
                            sg: ir.sg,
                            k_label: kp.label.clone(),
                            source_irrep: ir.ml.to_string(),
                            spinor: ir.spinor,
                            minimum_dimension: selected_arm_dimension(ir).ok(),
                            classified_type,
                            wigner_source,
                            classified_dimension,
                            reason: err.to_string(),
                        };
                        if retain_partial {
                            unresolved_coreps.push(failure);
                        } else {
                            return Err(MagneticIrrepError::CorepComputationFailed {
                                uni: failure.uni,
                                sg: failure.sg,
                                k_label: failure.k_label,
                                source_irrep: failure.source_irrep,
                                reason: failure.reason,
                            });
                        }
                    }
                }
            }
            let coreps = dedup_coreps(raw_coreps);
            let coreps = attach_isotropy_candidates(coreps, h_irreps);

            Ok(MagneticKPointSummary {
                label: kp.label,
                coords: kp.coords,
                little_group_order: mag_lg.len(),
                unitary_order,
                antiunitary_order,
                operations,
                conjugacy_classes,
                coreps,
            })
        })
        .collect();
    let kpoints = kpoints?;

    Ok(PartialMagneticIrrepSummary {
        summary: MagneticIrrepSummary {
            uni,
            bns_label: msg.bns_number.trim().to_string(),
            magnetic_type: msg.type_,
            parent_sg: msg.number as u8,
            unitary_sg: h_info.sg as u8,
            unitary_hall: h_info.hall,
            kpoints,
        },
        unresolved_coreps,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct MagneticOperationKey {
    rotation: [i32; 9],
    translation_twelfths: [i32; 3],
    time_reversal: bool,
}

fn magnetic_operation_keys(
    uni: usize,
    ops: &SymmetryOps,
) -> Result<Vec<MagneticOperationKey>, MagneticIrrepError> {
    const TRANSLATION_TOLERANCE: f64 = 1e-5;
    let mut keys = Vec::with_capacity(ops.len());
    for (index, operation) in ops.operations.iter().enumerate() {
        let mut translation_twelfths = [0; 3];
        for (axis, value) in operation.translation.iter().copied().enumerate() {
            if !value.is_finite() {
                return Err(MagneticIrrepError::OperationsInconsistentWithUni {
                    uni,
                    reason: format!("operation {index} has a non-finite translation"),
                });
            }
            let reduced = value.rem_euclid(1.0);
            let scaled = reduced * 12.0;
            let rounded = scaled.round();
            if (scaled - rounded).abs() / 12.0 > TRANSLATION_TOLERANCE {
                return Err(MagneticIrrepError::OperationsInconsistentWithUni {
                    uni,
                    reason: format!("operation {index} translation is not quantized in twelfths"),
                });
            }
            translation_twelfths[axis] = (rounded as i32).rem_euclid(12);
        }
        keys.push(MagneticOperationKey {
            rotation: [
                operation.rotation[0][0],
                operation.rotation[0][1],
                operation.rotation[0][2],
                operation.rotation[1][0],
                operation.rotation[1][1],
                operation.rotation[1][2],
                operation.rotation[2][0],
                operation.rotation[2][1],
                operation.rotation[2][2],
            ],
            translation_twelfths,
            time_reversal: operation.time_reversal,
        });
    }
    keys.sort_unstable();
    Ok(keys)
}

fn validate_operations_for_uni(uni: usize, ops: &SymmetryOps) -> Result<(), MagneticIrrepError> {
    if uni == 0 || uni > 1651 {
        return Err(MagneticIrrepError::InvalidUni(uni));
    }
    let database = SymmetryOps::from_magnetic_database(uni)
        .map_err(|_| MagneticIrrepError::MissingMagneticOperations(uni))?;
    if ops.len() != database.len() {
        return Err(MagneticIrrepError::OperationsInconsistentWithUni {
            uni,
            reason: format!(
                "operation count {} does not match database count {}",
                ops.len(),
                database.len()
            ),
        });
    }
    if magnetic_operation_keys(uni, ops)? != magnetic_operation_keys(uni, &database)? {
        return Err(MagneticIrrepError::OperationsInconsistentWithUni {
            uni,
            reason: "full magnetic Seitz multiset does not match the UNI first-Hall setting"
                .to_string(),
        });
    }
    Ok(())
}

// ── Formatting ─────────────────────────────────────────────────────────────────

/// Format a full magnetic irrep summary as human-readable text.
pub fn format_magnetic_irrep_summary(summary: &MagneticIrrepSummary) -> String {
    let mut lines = Vec::new();
    lines.push(format!(
        "Magnetic space group: UNI={}  BNS={}  type={:?}",
        summary.uni, summary.bns_label, summary.magnetic_type
    ));
    lines.push(format!(
        "  parent SG: {}  unitary subgroup H: {}  Hall: {}",
        summary.parent_sg, summary.unitary_sg, summary.unitary_hall
    ));
    lines.push(format!("  {} k-points", summary.kpoints.len()));

    for kp in &summary.kpoints {
        lines.push(String::new());
        lines.push(format_magnetic_kpoint_summary(kp));
    }

    lines.join("\n")
}

fn format_fraction(value: f64) -> String {
    let value = value.rem_euclid(1.0);
    let value = if (1.0 - value).abs() < 1e-9 {
        0.0
    } else {
        value
    };
    for denominator in 1..=48i64 {
        let numerator = (value * denominator as f64).round() as i64;
        if (value - numerator as f64 / denominator as f64).abs() < 1e-9 {
            return if denominator == 1 {
                numerator.to_string()
            } else {
                format!("{numerator}/{denominator}")
            };
        }
    }
    let rendered = format!("{value:.8}");
    rendered
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string()
}

fn format_operation(operation: &MagneticLittleGroupOperation) -> String {
    let rotation = operation.rotation;
    let translation = operation.translation.map(format_fraction);
    let prefix = if operation.time_reversal { "θ·" } else { "" };
    format!(
        "{prefix}{{R=[[{},{},{}],[{},{},{}],[{},{},{}]]; t=({},{},{})}}",
        rotation[0][0],
        rotation[0][1],
        rotation[0][2],
        rotation[1][0],
        rotation[1][1],
        rotation[1][2],
        rotation[2][0],
        rotation[2][1],
        rotation[2][2],
        translation[0],
        translation[1],
        translation[2]
    )
}

fn format_real_character(value: f64) -> String {
    if value.abs() < 1e-10 {
        return "0".to_string();
    }
    let integer = value.round();
    if (value - integer).abs() < 1e-9 {
        return format!("{integer:.0}");
    }
    let rendered = format!("{value:.8}");
    rendered
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string()
}

fn format_character(value: Complex64) -> String {
    if value.im.abs() < 1e-10 {
        return format_real_character(value.re);
    }
    if value.re.abs() < 1e-10 {
        let imaginary = format_real_character(value.im);
        return match imaginary.as_str() {
            "1" => "i".to_string(),
            "-1" => "-i".to_string(),
            _ => format!("{imaginary}i"),
        };
    }

    let real = format_real_character(value.re);
    let imaginary = format_real_character(value.im.abs());
    let sign = if value.im.is_sign_negative() {
        "-"
    } else {
        "+"
    };
    let imaginary = if imaginary == "1" {
        "i".to_string()
    } else {
        format!("{imaginary}i")
    };
    format!("{real}{sign}{imaginary}")
}

fn completeness_label(completeness: &crate::irrep::corep::CharacterCompleteness) -> String {
    match completeness {
        crate::irrep::corep::CharacterCompleteness::Complete => "complete".to_string(),
        crate::irrep::corep::CharacterCompleteness::TypeAAntiunitaryPending { count } => {
            format!("antiunitary-pending({count})")
        }
    }
}

fn same_character_signature(kpoint: &MagneticKPointSummary, left: usize, right: usize) -> bool {
    kpoint.coreps.iter().all(|corep| {
        corep
            .characters
            .get(left)
            .zip(corep.characters.get(right))
            .is_some_and(|(left, right)| match (left, right) {
                (Some(left), Some(right)) => (*left - *right).norm() < 1e-8,
                (None, None) => true,
                _ => false,
            })
    })
}

/// Return conjugacy classes refined only when the computed (possibly
/// projective) character rows are not constant on a raw Seitz class.
fn character_compatible_classes(kpoint: &MagneticKPointSummary) -> Vec<(String, Vec<usize>)> {
    let raw_classes = if kpoint.conjugacy_classes.is_empty() {
        singleton_conjugacy_classes(&kpoint.operations)
    } else {
        kpoint.conjugacy_classes.clone()
    };
    let mut result = Vec::new();
    for raw_class in raw_classes {
        let mut buckets: Vec<Vec<usize>> = Vec::new();
        for member in raw_class.members {
            if let Some(bucket) = buckets
                .iter_mut()
                .find(|bucket| same_character_signature(kpoint, bucket[0], member))
            {
                bucket.push(member);
            } else {
                buckets.push(vec![member]);
            }
        }
        let split = buckets.len() > 1;
        for (part, members) in buckets.into_iter().enumerate() {
            let label = if split {
                format!("C{}.{}", raw_class.class + 1, part + 1)
            } else {
                format!("C{}", raw_class.class + 1)
            };
            result.push((label, members));
        }
    }
    result
}

/// Format a complete magnetic character table with a selectable column layout.
///
/// No character values are truncated. Operation columns follow
/// [`MagneticKPointSummary::operations`] exactly. Conjugacy-class columns use
/// raw magnetic Seitz classes, refined when Bloch/projective characters are not
/// constant on a raw class.
pub fn format_magnetic_character_table_with_columns(
    kpoint: &MagneticKPointSummary,
    columns: MagneticCharacterTableColumns,
) -> String {
    if kpoint.operations.is_empty() {
        return "(no magnetic little-group operations)".to_string();
    }

    let display_columns: Vec<(String, Vec<usize>)> = match columns {
        MagneticCharacterTableColumns::Operations => kpoint
            .operations
            .iter()
            .map(|operation| (format!("g{}", operation.column + 1), vec![operation.column]))
            .collect(),
        MagneticCharacterTableColumns::ConjugacyClasses => character_compatible_classes(kpoint),
    };

    let mut lines = Vec::new();
    let mut header = vec![
        "corep".to_string(),
        "type".to_string(),
        "dim".to_string(),
        "status".to_string(),
    ];
    header.extend(display_columns.iter().map(|(label, members)| {
        if members.len() == 1 {
            label.clone()
        } else {
            format!("{label} (×{})", members.len())
        }
    }));
    lines.push(format!("| {} |", header.join(" | ")));
    lines.push(format!(
        "| {} |",
        std::iter::repeat_n("---", header.len())
            .collect::<Vec<_>>()
            .join(" | ")
    ));
    for corep in &kpoint.coreps {
        let mut row = vec![
            corep.label.clone(),
            format!("{:?}", corep.corep_type),
            corep.dim.to_string(),
            completeness_label(&corep.completeness),
        ];
        row.extend(display_columns.iter().map(
            |(_, members)| match corep.characters.get(members[0]) {
                Some(Some(value)) => format_character(*value),
                Some(None) => "N/A".to_string(),
                None => "?".to_string(),
            },
        ));
        lines.push(format!("| {} |", row.join(" | ")));
    }

    lines.push(String::new());
    lines.push("Column definitions:".to_string());
    match columns {
        MagneticCharacterTableColumns::Operations => {
            lines.push(
                "| column | MSG op index | kind | Seitz operation (data-Hall frame) |".to_string(),
            );
            lines.push("| --- | ---: | --- | --- |".to_string());
            for operation in &kpoint.operations {
                lines.push(format!(
                    "| g{} | {} | {} | {} |",
                    operation.column + 1,
                    operation.magnetic_operation_index,
                    if operation.time_reversal {
                        "antiunitary"
                    } else {
                        "unitary"
                    },
                    format_operation(operation)
                ));
            }
        }
        MagneticCharacterTableColumns::ConjugacyClasses => {
            lines.push("| class | size | member operation columns | representative |".to_string());
            lines.push("| --- | ---: | --- | --- |".to_string());
            for (label, members) in &display_columns {
                let member_labels = members
                    .iter()
                    .map(|member| format!("g{}", member + 1))
                    .collect::<Vec<_>>()
                    .join(", ");
                lines.push(format!(
                    "| {label} | {} | {member_labels} | {} |",
                    members.len(),
                    format_operation(&kpoint.operations[members[0]])
                ));
            }
        }
    }
    lines.join("\n")
}

/// Format the complete table with one column per magnetic little-group
/// operation.
pub fn format_magnetic_character_table(kpoint: &MagneticKPointSummary) -> String {
    format_magnetic_character_table_with_columns(kpoint, MagneticCharacterTableColumns::Operations)
}

/// Format the complete table with character-compatible conjugacy classes as
/// columns.
pub fn format_magnetic_character_table_by_class(kpoint: &MagneticKPointSummary) -> String {
    format_magnetic_character_table_with_columns(
        kpoint,
        MagneticCharacterTableColumns::ConjugacyClasses,
    )
}

/// Format a single k-point summary as human-readable text.
pub fn format_magnetic_kpoint_summary(kpoint: &MagneticKPointSummary) -> String {
    let mut lines = Vec::new();
    let (kx, ky, kz, kd) = kpoint.coords;
    lines.push(format!(
        "k-point {}  ({}/{}, {}/{}, {}/{})  |LG|= {}  ({}U + {}A)",
        kpoint.label,
        kx,
        kd,
        ky,
        kd,
        kz,
        kd,
        kpoint.little_group_order,
        kpoint.unitary_order,
        kpoint.antiunitary_order
    ));

    if kpoint.coreps.is_empty() {
        lines.push("  (no coreps)".to_string());
        return lines.join("\n");
    }

    lines.push(String::new());
    lines.push(format_magnetic_character_table(kpoint));

    for c in &kpoint.coreps {
        if !c.isotropy_candidates.is_empty() {
            for ic in &c.isotropy_candidates {
                let n_ord = ic.ordinary.len();
                let n_mag = ic.magnetic.len();
                // Always show spinor entries (with 0 subgroups) so users can
                // distinguish "no data" from "data not checked".
                if n_ord > 0
                    || n_mag > 0
                    || ic.relation == IsotropyCandidateRelation::SpinorNoIsotropyData
                {
                    lines.push(format!(
                        "isotropy (corep {}, source {} {:?}): {} ordinary + {} magnetic subgroups",
                        c.label, ic.source_ml, ic.relation, n_ord, n_mag
                    ));
                }
            }
        }
    }

    lines.join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn magnetic_summary_by_uni_smoke() {
        // UNI 2 = BNS 1.2 (grey P1)
        let s = magnetic_irrep_summary_by_uni(2).unwrap();
        assert_eq!(s.uni, 2);
        assert_eq!(s.bns_label, "1.2");
        assert!(!s.kpoints.is_empty(), "should have at least one k-point");
    }

    #[test]
    fn magnetic_summary_by_bns_matches_uni() {
        let by_uni = magnetic_irrep_summary_by_uni(2).unwrap();
        let by_bns = magnetic_irrep_summary_by_bns("1.2").unwrap();
        assert_eq!(by_uni.uni, by_bns.uni);
        assert_eq!(by_uni.unitary_sg, by_bns.unitary_sg);
        assert_eq!(by_uni.kpoints.len(), by_bns.kpoints.len());
    }

    #[test]
    fn summary_corep_columns_remain_aligned_with_operation_columns() {
        let summary = magnetic_irrep_summary_by_bns("1.2").unwrap();
        for kpoint in &summary.kpoints {
            for corep in &kpoint.coreps {
                assert_eq!(corep.characters.len(), kpoint.operations.len());
                assert_eq!(corep.timerev.len(), kpoint.operations.len());
                for (character, time_reversal, operation) in corep
                    .characters
                    .iter()
                    .zip(&corep.timerev)
                    .zip(&kpoint.operations)
                    .map(|((character, time_reversal), operation)| {
                        (character, time_reversal, operation)
                    })
                {
                    match character {
                        Some(value) => {
                            assert!(value.re.is_finite());
                            assert!(value.im.is_finite());
                        }
                        None => {
                            assert!(operation.time_reversal);
                            assert!(matches!(
                                corep.completeness,
                                crate::irrep::corep::CharacterCompleteness::TypeAAntiunitaryPending { .. }
                            ));
                        }
                    }
                    assert_eq!(*time_reversal, operation.time_reversal);
                }
            }
        }
    }

    #[test]
    fn invalid_uni_returns_error() {
        assert!(magnetic_irrep_summary_by_uni(0).is_err());
        assert!(magnetic_irrep_summary_by_uni(1652).is_err());
    }

    #[test]
    fn explicit_ops_validate_uni_before_summary_metadata() {
        let empty = SymmetryOps::default();
        assert!(matches!(
            magnetic_irrep_summary_from_ops(0, &empty),
            Err(MagneticIrrepError::InvalidUni(0))
        ));
        assert!(matches!(
            magnetic_irrep_summary_from_ops(1652, &empty),
            Err(MagneticIrrepError::InvalidUni(1652))
        ));
        assert!(matches!(
            magnetic_irrep_summary_from_ops(2, &empty),
            Err(MagneticIrrepError::OperationsInconsistentWithUni { uni: 2, .. })
        ));
    }

    #[test]
    fn explicit_ops_accept_reordering_and_integer_lattice_shifts() {
        let database = SymmetryOps::from_magnetic_database(2).unwrap();
        let mut reordered = database.clone();
        reordered.operations.reverse();
        validate_operations_for_uni(2, &reordered).unwrap();

        let summary = magnetic_irrep_summary_from_ops(2, &reordered).unwrap();
        assert_eq!(summary.uni, 2);
        assert_eq!(summary.bns_label, "1.2");

        let mut shifted = database;
        for (index, operation) in shifted.operations.iter_mut().enumerate() {
            operation.translation[0] += 1.0;
            operation.translation[1] -= (index % 2) as f64;
            operation.translation[2] += 2.0;
        }
        validate_operations_for_uni(2, &shifted).unwrap();
    }

    #[test]
    fn explicit_ops_reject_wrong_uni_flags_missing_and_duplicate_operations() {
        let ops_282 = SymmetryOps::from_magnetic_database(282).unwrap();
        assert!(matches!(
            validate_operations_for_uni(283, &ops_282),
            Err(MagneticIrrepError::OperationsInconsistentWithUni { uni: 283, .. })
        ));

        let mut wrong_flag = ops_282.clone();
        wrong_flag.operations[0].time_reversal = !wrong_flag.operations[0].time_reversal;
        assert!(validate_operations_for_uni(282, &wrong_flag).is_err());

        let mut incomplete = ops_282.clone();
        incomplete.operations.pop();
        assert!(validate_operations_for_uni(282, &incomplete).is_err());

        let mut duplicate = ops_282.clone();
        duplicate.operations[0] = duplicate.operations[1];
        assert!(validate_operations_for_uni(282, &duplicate).is_err());
    }

    #[test]
    fn explicit_ops_distinguish_uni_277_and_284_in_first_hall_frame() {
        let ops_277 = SymmetryOps::from_magnetic_database(277).unwrap();
        let ops_284 = SymmetryOps::from_magnetic_database(284).unwrap();
        assert_ne!(
            magnetic_operation_keys(277, &ops_277).unwrap(),
            magnetic_operation_keys(284, &ops_284).unwrap()
        );
        assert!(validate_operations_for_uni(277, &ops_284).is_err());
        assert!(validate_operations_for_uni(284, &ops_277).is_err());

        assert!(validate_operations_for_uni(275, &ops_277).is_err());
        let ops_283 = SymmetryOps::from_magnetic_database(283).unwrap();
        validate_operations_for_uni(283, &ops_283).unwrap();
        assert!(validate_operations_for_uni(283, &ops_284).is_err());
    }

    #[test]
    fn unknown_bns_returns_error() {
        assert!(magnetic_irrep_summary_by_bns("999.999").is_err());
    }

    #[test]
    fn grey_group_has_antiunitary_ops() {
        // UNI 2 = BNS 1.2 (grey P1, Type II): every k-point should
        // have antiunitary operations because time reversal is a symmetry.
        let s = magnetic_irrep_summary_by_uni(2).unwrap();
        assert_eq!(s.magnetic_type, crate::MagneticType::Grey);
        for kp in &s.kpoints {
            assert!(
                kp.antiunitary_order > 0,
                "grey group: k-point {} should have antiunitary ops, got 0",
                kp.label
            );
            assert!(kp.little_group_order > 0);
            assert_eq!(
                kp.little_group_order,
                kp.unitary_order + kp.antiunitary_order
            );
        }
    }

    #[test]
    fn ordinary_group_has_no_antiunitary_ops() {
        // UNI 1 = BNS 1.1 (ordinary P1, Type I): all ops are unitary.
        let s = magnetic_irrep_summary_by_uni(1).unwrap();
        assert_eq!(s.magnetic_type, crate::MagneticType::Ordinary);
        for kp in &s.kpoints {
            assert_eq!(
                kp.antiunitary_order, 0,
                "ordinary group: k-point {} should have no antiunitary ops",
                kp.label
            );
        }
    }

    #[test]
    fn isotropy_candidates_attached_to_coreps() {
        // UNI 2 = BNS 1.2 (grey P1): fully classified and includes spinor
        // sources, whose isotropy data should be explicitly marked unavailable.
        let s = magnetic_irrep_summary_by_uni(2).unwrap();
        let spinor_corep = s
            .kpoints
            .iter()
            .flat_map(|kp| kp.coreps.iter())
            .find(|c| c.source_irreps.iter().any(|src| src.spinor))
            .expect("UNI 2 should have at least one spinor corep");
        assert!(
            spinor_corep
                .isotropy_candidates
                .iter()
                .any(|ic| ic.relation == IsotropyCandidateRelation::SpinorNoIsotropyData),
            "spinor corep should explicitly report missing isotropy data"
        );
    }

    fn assert_well_formed_summary(summary: &MagneticIrrepSummary) {
        for kpoint in &summary.kpoints {
            assert_eq!(kpoint.operations.len(), kpoint.little_group_order);
            for (column, operation) in kpoint.operations.iter().enumerate() {
                assert_eq!(operation.column, column);
            }

            let mut class_members = kpoint
                .conjugacy_classes
                .iter()
                .flat_map(|class| class.members.iter().copied())
                .collect::<Vec<_>>();
            class_members.sort_unstable();
            assert_eq!(
                class_members,
                (0..kpoint.operations.len()).collect::<Vec<_>>()
            );

            let identity = kpoint
                .operations
                .iter()
                .position(|operation| {
                    !operation.time_reversal
                        && operation.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                        && operation
                            .translation
                            .iter()
                            .all(|value| periodic_translation_key(*value) == 0)
                })
                .expect("little group must contain identity");
            for corep in &kpoint.coreps {
                assert_eq!(corep.characters.len(), kpoint.operations.len());
                assert_eq!(corep.timerev.len(), kpoint.operations.len());
                assert!(corep.characters.iter().all(|value| {
                    value.is_none_or(|value| value.re.is_finite() && value.im.is_finite())
                }));
                let identity_character =
                    corep.characters[identity].expect("unitary identity character must be defined");
                assert!(
                    (identity_character - Complex64::new(corep.dim as f64, 0.0)).norm() < 1e-6,
                    "{} {}: χ(E)={} != dim={}",
                    kpoint.label,
                    corep.label,
                    identity_character,
                    corep.dim
                );
            }
        }
    }

    fn assert_source_dimension_relations(summary: &MagneticIrrepSummary) {
        for kpoint in &summary.kpoints {
            for corep in &kpoint.coreps {
                assert!(!corep.source_irreps.is_empty());
                let source_dim_sum = corep
                    .source_irreps
                    .iter()
                    .map(|source| source.dim as usize)
                    .sum::<usize>();
                match corep.corep_type {
                    crate::irrep::corep::CorepType::A => {
                        assert_eq!(
                            corep.dim, source_dim_sum,
                            "{} {}: Type A dimension must equal the sum of its selected-arm source dimensions",
                            kpoint.label, corep.label
                        );
                    }
                    crate::irrep::corep::CorepType::B => {
                        assert_eq!(
                            corep.dim,
                            2 * source_dim_sum,
                            "{} {}: Type B dimension must double the sum of its selected-arm source dimensions",
                            kpoint.label,
                            corep.label
                        );
                    }
                    crate::irrep::corep::CorepType::C => {
                        let source_dim = corep.source_irreps[0].dim as usize;
                        assert!(
                            corep
                                .source_irreps
                                .iter()
                                .all(|source| source.dim as usize == source_dim)
                        );
                        assert_eq!(
                            corep.dim,
                            2 * source_dim,
                            "{} {}: Type C must double the common partner-source dimension",
                            kpoint.label,
                            corep.label
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn bns_128_406_preserves_complex_compound_characters() {
        let summary = magnetic_irrep_summary_by_bns("128.406")
            .expect("the strict summary must preserve compound complex rows");
        let z = summary
            .kpoints
            .iter()
            .find(|point| point.label == "Z")
            .expect("missing Z point");
        let compound = z
            .coreps
            .iter()
            .find(|corep| {
                corep.source == crate::irrep::corep::WignerSource::ScalarCIR
                    && corep.corep_type == crate::irrep::corep::CorepType::C
            })
            .expect("missing compound Type-C branch");
        assert!(
            compound
                .characters
                .iter()
                .flatten()
                .any(|value| value.im.abs() > 1.0),
            "the genuinely complex compound columns must not be projected to f64"
        );
    }

    #[test]
    fn no_duplicate_coreps_at_any_kpoint() {
        // For a fully classified group, each k-point should have no duplicate coreps
        // (duplicated = same type + dim + rounded characters + timerev).
        let s = magnetic_irrep_summary_by_uni(2).unwrap();
        for kp in &s.kpoints {
            for c in &kp.coreps {
                assert!(
                    !c.source_irreps.is_empty(),
                    "k-point {}: corep has no source irrep",
                    kp.label
                );
            }
            // Check no duplicates among computed coreps.
            for i in 0..kp.coreps.len() {
                for j in (i + 1)..kp.coreps.len() {
                    let ci = &kp.coreps[i];
                    let cj = &kp.coreps[j];
                    let same_type = ci.corep_type == cj.corep_type;
                    let same_dim = ci.dim == cj.dim;
                    let same_chars = round_chars(&ci.characters) == round_chars(&cj.characters);
                    let same_tr = ci.timerev == cj.timerev;
                    assert!(
                        !(same_type && same_dim && same_chars && same_tr),
                        "k-point {}: duplicate coreps {} and {} at indices {} and {}",
                        kp.label,
                        ci.label,
                        cj.label,
                        i,
                        j
                    );
                }
            }
        }
    }

    #[test]
    fn coreps_at_gm_for_grey_p1() {
        // UNI 2 = BNS 1.2 (grey P1): time reversal is a symmetry.
        let s = magnetic_irrep_summary_by_uni(2).unwrap();
        let gm_kp = s
            .kpoints
            .iter()
            .find(|kp| kp.label == "GM")
            .expect("should have GM k-point");
        assert!(!gm_kp.coreps.is_empty(), "should have coreps at GM");
        for c in &gm_kp.coreps {
            assert!(c.dim > 0, "corep {} has zero dimension", c.label);
            // Identity character (first) should be close to dim.
            if !c.characters.is_empty() {
                let chi_e = c.characters[0].expect("identity character must be defined");
                assert!(
                    (chi_e - Complex64::new(c.dim as f64, 0.0)).norm() < 1e-6,
                    "corep {}: χ(E)={} != dim={}",
                    c.label,
                    chi_e,
                    c.dim
                );
            }
        }
    }

    #[test]
    fn pending_type_a_columns_are_distinct_from_defined_zeroes() {
        let pending_summary = magnetic_irrep_summary_by_uni(3).unwrap();
        let pending_gm = pending_summary
            .kpoints
            .iter()
            .find(|point| point.label == "GM")
            .expect("missing anti-translation GM point");
        let pending_antiunitary = pending_gm
            .operations
            .iter()
            .position(|operation| operation.time_reversal)
            .expect("UNI 3 must contain an antiunitary operation");
        let type_a = pending_gm
            .coreps
            .iter()
            .find(|corep| corep.corep_type == crate::irrep::corep::CorepType::A)
            .expect("missing Type-A corepresentation");
        assert_eq!(type_a.characters[pending_antiunitary], None);
        assert!(format_magnetic_character_table(pending_gm).contains("N/A"));

        let defined_summary = magnetic_irrep_summary_by_uni(2).unwrap();
        let defined_gm = defined_summary
            .kpoints
            .iter()
            .find(|point| point.label == "GM")
            .expect("missing grey-P1 GM point");
        let defined_antiunitary = defined_gm
            .operations
            .iter()
            .position(|operation| operation.time_reversal)
            .expect("grey P1 must contain time reversal");
        let type_b = defined_gm
            .coreps
            .iter()
            .find(|corep| corep.corep_type == crate::irrep::corep::CorepType::B)
            .expect("missing Type-B corepresentation");
        assert_eq!(
            type_b.characters[defined_antiunitary],
            Some(Complex64::new(0.0, 0.0))
        );
    }

    #[test]
    fn bns_182_183_includes_compound_branches() {
        let summary = magnetic_irrep_summary_by_bns("182.183")
            .expect("compound branches must be part of the strict summary");
        let gm = summary
            .kpoints
            .iter()
            .find(|point| point.label == "GM")
            .expect("missing GM point");
        assert!(gm.coreps.iter().any(|corep| {
            corep.source == crate::irrep::corep::WignerSource::ScalarCIR
                && corep
                    .source_irreps
                    .iter()
                    .any(|source| source.ml == "GM3" || source.ml == "GM5")
        }));
    }

    #[test]
    fn partial_summary_has_no_compound_gap_after_complex_migration() {
        let partial = magnetic_irrep_summary_by_uni_partial(1413)
            .expect("partial compatibility surface must share the strict complex path");
        assert_eq!(partial.summary.uni, 1413);
        assert!(
            partial
                .summary
                .kpoints
                .iter()
                .flat_map(|point| &point.coreps)
                .next()
                .is_some(),
            "the partial summary must not discard every safe corep"
        );
        assert!(partial.unresolved_coreps.is_empty());
    }

    #[test]
    fn strict_summary_preserves_completed_complex_classification() {
        let partial = magnetic_irrep_summary_by_uni_partial(1440)
            .expect("complex source rows must be available to typed consumers");
        assert!(partial.unresolved_coreps.is_empty());
        let corep = partial
            .summary
            .kpoints
            .iter()
            .find(|point| point.label == "K")
            .expect("missing K point")
            .coreps
            .iter()
            .find(|corep| {
                corep
                    .source_irreps
                    .iter()
                    .any(|source| source.ml == "K3" && !source.spinor)
            })
            .expect("missing SG187 K3 corepresentation");
        assert_eq!(corep.corep_type, crate::irrep::corep::CorepType::A);
        assert_eq!(corep.dim, 1);
        assert!(
            corep
                .characters
                .iter()
                .flatten()
                .any(|value| value.im.abs() > 1.0e-6)
        );
    }

    #[test]
    fn formal_formatters_emit_every_operation_and_class_column() {
        let summary = magnetic_irrep_summary_by_bns("1.2").unwrap();
        let z = summary
            .kpoints
            .iter()
            .find(|kpoint| kpoint.label == "GM")
            .expect("missing GM point");
        assert_eq!(z.operations.len(), 2);

        let operations = format_magnetic_character_table(z);
        assert!(operations.contains("| corep | type | dim | status | g1 |"));
        assert!(operations.contains("| g2 |"));
        assert!(operations.contains("Seitz operation (data-Hall frame)"));
        assert!(
            !operations.contains("..."),
            "formatter must not truncate rows"
        );

        let classes = format_magnetic_character_table_by_class(z);
        assert!(classes.contains("member operation columns"));
        assert!(classes.contains("| C1"));
    }

    /// Exhaustive release-only audit for the strict summary boundary.
    ///
    /// Run with:
    /// `cargo test --release --package cryspglib audit_all_1651_magnetic_summaries -- --ignored --nocapture`
    #[test]
    #[ignore = "exhaustive 1651-UNI magnetic summary audit"]
    fn audit_all_1651_magnetic_summaries() {
        let selected_uni = std::env::var("CRYSPGLIB_AUDIT_UNI")
            .ok()
            .and_then(|value| value.parse::<usize>().ok());
        let unis: Vec<_> = selected_uni.map_or_else(|| (1..=1651).collect(), |uni| vec![uni]);
        let mut failures = Vec::new();
        let mut kpoint_count = 0usize;
        let mut corep_count = 0usize;
        let amplified_noise_targets = [
            std::f64::consts::PI / 30.0,
            std::f64::consts::PI / (10.0 * 3.0_f64.sqrt()),
            std::f64::consts::PI / 15.0,
            std::f64::consts::PI / 10.0,
            3.0_f64.sqrt() * std::f64::consts::PI / 10.0,
            std::f64::consts::PI / 5.0,
        ];
        let mut amplified_noise_count = 0usize;
        for &uni in &unis {
            match magnetic_irrep_summary_by_uni(uni) {
                Ok(summary) => {
                    assert_well_formed_summary(&summary);
                    assert_source_dimension_relations(&summary);
                    amplified_noise_count += summary
                        .kpoints
                        .iter()
                        .flat_map(|kpoint| &kpoint.coreps)
                        .flat_map(|corep| &corep.characters)
                        .flatten()
                        .filter(|value| {
                            amplified_noise_targets.iter().any(|target| {
                                (value.re.abs() - target).abs() < 2e-6
                                    || (value.im.abs() - target).abs() < 2e-6
                            })
                        })
                        .count();
                    kpoint_count += summary.kpoints.len();
                    corep_count += summary
                        .kpoints
                        .iter()
                        .map(|kpoint| kpoint.coreps.len())
                        .sum::<usize>();
                }
                Err(error) => failures.push((uni, error)),
            }
        }
        eprintln!(
            "magnetic summary audit: success={} failure={} kpoints={} coreps={} amplified_noise={}",
            unis.len() - failures.len(),
            failures.len(),
            kpoint_count,
            corep_count,
            amplified_noise_count
        );
        let mut categories = std::collections::BTreeMap::<&str, usize>::new();
        for (_, error) in &failures {
            let category = match error {
                MagneticIrrepError::CorepComputationFailed { reason, .. }
                    if reason.contains("scalar PIR operation map") =>
                {
                    "scalar PIR operation map"
                }
                MagneticIrrepError::CorepComputationFailed { reason, .. }
                    if reason.contains("selected k-arm block") =>
                {
                    "selected k-arm block"
                }
                MagneticIrrepError::CorepComputationFailed { reason, .. }
                    if reason.contains("AntiunitarySpinLookup") =>
                {
                    "spinor AntiunitarySpinLookup"
                }
                MagneticIrrepError::CorepComputationFailed { reason, .. }
                    if reason.contains("spinor SU(2)") =>
                {
                    "other spinor SU(2)"
                }
                MagneticIrrepError::CorepComputationFailed { reason, .. }
                    if reason.contains("scalar PIR Wigner") =>
                {
                    "scalar PIR Wigner"
                }
                MagneticIrrepError::CorepComputationFailed { .. } => "other corep",
                _ => "summary setup",
            };
            *categories.entry(category).or_default() += 1;
        }
        for (category, count) in categories {
            eprintln!("  category {category}: {count}");
        }
        let failure_limit = if std::env::var_os("CRYSPGLIB_AUDIT_VERBOSE").is_some() {
            usize::MAX
        } else {
            50
        };
        for (uni, error) in failures.iter().take(failure_limit) {
            eprintln!("  UNI {uni}: {error:?}");
        }
        assert_eq!(
            amplified_noise_count, 0,
            "amplified exponent noise propagated into magnetic summaries"
        );
    }

    /// Regression: spinor coreps must carry SpinorNoIsotropyData candidates
    /// when no isotropy data is available for spinor source irreps.
    #[test]
    fn spinor_coreps_have_spinor_no_isotropy_data() {
        let s = magnetic_irrep_summary_by_uni(2).unwrap();
        let c = s
            .kpoints
            .iter()
            .flat_map(|kp| kp.coreps.iter())
            .find(|c| c.source_irreps.iter().any(|src| src.spinor))
            .expect("missing spinor corep");
        let n_spinor_sources = c.source_irreps.iter().filter(|src| src.spinor).count();
        let spinor_candidates: Vec<_> = c
            .isotropy_candidates
            .iter()
            .filter(|ic| ic.relation == IsotropyCandidateRelation::SpinorNoIsotropyData)
            .collect();
        assert_eq!(
            spinor_candidates.len(),
            n_spinor_sources,
            "spinor corep should have one SpinorNoIsotropyData candidate per spinor source"
        );
    }

    /// Interactive demo: print the full UNI 2 summary.
    ///
    /// Run with `-- --nocapture` to see the output.
    #[test]
    fn demo_uni2_summary() {
        let s = magnetic_irrep_summary_by_uni(2).unwrap();
        println!(
            "BNS {}  UNI={}  type={:?}  H=SG{}",
            s.bns_label, s.uni, s.magnetic_type, s.unitary_sg
        );
        for kp in &s.kpoints {
            let (kx, ky, kz, kd) = kp.coords;
            println!();
            println!(
                "k-point {}  ({}/{}, {}/{}, {}/{})  |LG|={} ({}U+{}A)  coreps={}",
                kp.label,
                kx,
                kd,
                ky,
                kd,
                kz,
                kd,
                kp.little_group_order,
                kp.unitary_order,
                kp.antiunitary_order,
                kp.coreps.len()
            );
            for c in &kp.coreps {
                let srcs: Vec<&str> = c.source_irreps.iter().map(|s| s.ml).collect();
                let chi0 = c
                    .characters
                    .first()
                    .and_then(|value| *value)
                    .map_or("N/A".to_string(), |value| format!("{}", value));
                println!(
                    "  {:20}  type={:?}  dim={}  χ(E)={}  src=[{}]",
                    c.label,
                    c.corep_type,
                    c.dim,
                    chi0,
                    srcs.join(", ")
                );
            }
        }
    }

    #[test]
    fn sg221_x1plus_summary_reports_selected_arm_dimension() {
        let irrep = crate::irrep::query::irreps_of(221)
            .iter()
            .find(|irrep| !irrep.spinor && irrep.k_label() == "X" && irrep.ml == "X1+")
            .expect("SG221 X1+ scalar irrep");
        assert_eq!(
            irrep.dim, 3,
            "X1+ must retain its full-star image dimension"
        );
        assert_eq!(selected_arm_dimension(irrep).unwrap(), 1);

        // UNI 1594 (BNS 221.92) is a currently successful SG221 summary at X.
        let summary = magnetic_irrep_summary_by_uni(1594).expect("UNI 1594 summary");
        let x = summary
            .kpoints
            .iter()
            .find(|kpoint| kpoint.label == "X")
            .expect("SG221 X k-point");
        let corep = x
            .coreps
            .iter()
            .find(|corep| corep.source_irreps.iter().any(|source| source.ml == "X1+"))
            .expect("X1+ source corep");
        assert_eq!(corep.corep_type, crate::irrep::corep::CorepType::A);
        assert_eq!(corep.source_irreps.len(), 1);
        assert_eq!(corep.source_irreps[0].dim, 1);
        assert_eq!(corep.dim, corep.source_irreps[0].dim as usize);
    }

    #[test]
    fn type_b_summary_doubles_selected_arm_source_dimension() {
        let summary = magnetic_irrep_summary_by_uni(2).expect("UNI 2 summary");
        let corep = summary
            .kpoints
            .iter()
            .flat_map(|kpoint| &kpoint.coreps)
            .find(|corep| corep.corep_type == crate::irrep::corep::CorepType::B)
            .expect("UNI 2 should contain a Type B corep");
        assert_eq!(corep.source_irreps.len(), 1);
        assert_eq!(
            corep.dim,
            2 * corep.source_irreps[0].dim as usize,
            "Type B must double the selected-arm source dimension"
        );
    }

    #[test]
    fn type_c_summary_merges_two_selected_arm_source_dimensions() {
        let summary = magnetic_irrep_summary_by_uni(9).expect("UNI 9 summary");
        let corep = summary
            .kpoints
            .iter()
            .flat_map(|kpoint| &kpoint.coreps)
            .find(|corep| corep.corep_type == crate::irrep::corep::CorepType::C)
            .expect("UNI 9 should contain a Type C corep");
        assert_eq!(corep.source_irreps.len(), 2);
        let source_dim_sum = corep
            .source_irreps
            .iter()
            .map(|source| source.dim as usize)
            .sum::<usize>();
        assert_eq!(
            corep.dim, source_dim_sum,
            "deduplicated Type C dimension must equal source dimension sum"
        );
    }

    #[test]
    fn type_c_dedup_keeps_scalar_and_spinor_families_separate() {
        let summary = magnetic_irrep_summary_by_uni(7).expect("UNI 7 summary");
        let z = summary
            .kpoints
            .iter()
            .find(|kpoint| kpoint.label == "Z")
            .expect("UNI 7 Z point");
        let type_c = z
            .coreps
            .iter()
            .filter(|corep| corep.corep_type == crate::irrep::corep::CorepType::C)
            .collect::<Vec<_>>();
        assert_eq!(
            type_c.len(),
            2,
            "Z has one scalar and one spin Type-C corep"
        );
        let scalar = type_c
            .iter()
            .find(|corep| corep.source_irreps.iter().all(|source| !source.spinor))
            .expect("scalar Z Type-C pair");
        let spinor = type_c
            .iter()
            .find(|corep| corep.source_irreps.iter().all(|source| source.spinor))
            .expect("spinor Z Type-C pair");
        let sorted_labels = |corep: &MagneticCorepSummary| {
            let mut labels = corep
                .source_irreps
                .iter()
                .map(|source| source.ml)
                .collect::<Vec<_>>();
            labels.sort_unstable();
            labels
        };
        assert_eq!(sorted_labels(scalar), vec!["Z1+", "Z1-"]);
        assert_eq!(sorted_labels(spinor), vec!["Z2", "Z3"]);
    }
}
