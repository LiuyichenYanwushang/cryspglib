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
//!         let chi0 = c.characters.first().map_or("N/A".to_string(), |v| format!("{:.0}", v));
//!         println!("  {:20}  type={:?}  dim={}  χ(E)={}  src=[{}]",
//!             c.label, c.corep_type, c.dim, chi0, srcs.join(", "));
//!     }
//! }
//! ```

use std::collections::{BTreeSet, HashMap};

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
    pub characters: Vec<f64>,
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
    /// Irrep dimension.
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

// ── Type-C dedup ───────────────────────────────────────────────────────────────

/// Deduplicate coreps by `(corep_type, dim, rounded characters, timerev)`.
///
/// Only Type-C coreps with finite character values are eligible for merging.
/// Other valid coreps pass through without deduplication. Unsupported and
/// non-finite coreps should be rejected before this function is called.
fn dedup_coreps(coreps: Vec<MagneticCorepSummary>) -> Vec<MagneticCorepSummary> {
    let mut groups: Vec<Vec<MagneticCorepSummary>> = Vec::new();

    for c in coreps {
        // Only Type-C with finite characters can be deduplicated.
        let can_dedup = c.corep_type == crate::irrep::corep::CorepType::C
            && c.dim > 0
            && c.characters.iter().all(|&ch| ch.is_finite());
        if !can_dedup {
            // Pass through — never merge non-C entries.
            groups.push(vec![c]);
            continue;
        }

        let key = (c.dim, round_chars(&c.characters), c.timerev.clone());
        let found = groups.iter_mut().find(|g| {
            let first = &g[0];
            first.corep_type == crate::irrep::corep::CorepType::C
                && first.dim == key.0
                && round_chars(&first.characters) == key.1
                && first.timerev == key.2
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
            // Build combined label: sort source ML labels and join with " + ".
            let mut labels: Vec<&str> = merged.source_irreps.iter().map(|s| s.ml).collect();
            labels.sort();
            merged.label = labels.join(" + ");
            merged
        })
        .collect()
}

/// Round character values to integers for dedup comparison.
fn round_chars(chars: &[f64]) -> Vec<i64> {
    chars.iter().map(|&c| (c * 1e8).round() as i64).collect()
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
                match crate::irrep::corep::compute_corepresentation(ir, uni, mag_ops) {
                    Ok(c) => raw_coreps.push(MagneticCorepSummary {
                        label: ir.ml.to_string(),
                        source_irreps: vec![SourceIrrepSummary {
                            sg: ir.sg,
                            ml: ir.ml,
                            bc: ir.bc,
                            dim: ir.dim,
                            spinor: ir.spinor,
                        }],
                        corep_type: c.corep_type,
                        source: c.source,
                        dim: c.dim,
                        characters: c.characters,
                        timerev: c.timerev,
                        completeness: c.completeness,
                        isotropy_candidates: Vec::new(),
                    }),
                    Err(err) => {
                        return Err(MagneticIrrepError::CorepComputationFailed {
                            uni,
                            sg: ir.sg,
                            k_label: kp.label.clone(),
                            source_irrep: ir.ml.to_string(),
                            reason: err.to_string(),
                        });
                    }
                }
            }
            for corep in &raw_coreps {
                let aligned = corep.characters.len() == operations.len()
                    && corep.timerev.len() == operations.len()
                    && corep
                        .timerev
                        .iter()
                        .zip(&operations)
                        .all(|(time_reversal, operation)| {
                            *time_reversal == operation.time_reversal
                        });
                if !aligned {
                    return Err(MagneticIrrepError::CorepComputationFailed {
                        uni,
                        sg: h_info.sg as u8,
                        k_label: kp.label.clone(),
                        source_irrep: corep.label.clone(),
                        reason: "character columns are not aligned with magnetic little-group operations"
                            .to_string(),
                    });
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

    Ok(MagneticIrrepSummary {
        uni,
        bns_label: msg.bns_number.trim().to_string(),
        magnetic_type: msg.type_,
        parent_sg: msg.number as u8,
        unitary_sg: h_info.sg as u8,
        unitary_hall: h_info.hall,
        kpoints,
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

fn format_character(value: f64) -> String {
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
            .is_some_and(|(left, right)| (left - right).abs() < 1e-8)
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
        row.extend(display_columns.iter().map(|(_, members)| {
            corep
                .characters
                .get(members[0])
                .map_or_else(|| "?".to_string(), |value| format_character(*value))
        }));
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
                assert!(corep.characters.iter().all(|value| value.is_finite()));
                assert!(
                    (corep.characters[identity] - corep.dim as f64).abs() < 1e-6,
                    "{} {}: χ(E)={} != dim={}",
                    kpoint.label,
                    corep.label,
                    corep.characters[identity],
                    corep.dim
                );
            }
        }
    }

    #[test]
    fn bns_128_406_has_official_dimensions_and_pending_status() {
        let error = magnetic_irrep_summary_by_bns("128.406")
            .expect_err("compound corepresentation must fail closed in summaries");
        match error {
            MagneticIrrepError::CorepComputationFailed {
                uni,
                sg,
                k_label,
                source_irrep,
                reason,
            } => {
                assert_eq!(uni, 1066);
                assert_eq!(sg, 118);
                assert_eq!(k_label, "Z");
                assert_eq!(source_irrep, "Z1Z4");
                assert!(reason.contains("compound corepresentations"));
                assert!(reason.contains("complex operation-aware"));
            }
            other => panic!("unexpected 128.406 summary error: {other:?}"),
        }
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
                let chi_e = c.characters[0];
                assert!(
                    (chi_e - c.dim as f64).abs() < 1e-6,
                    "corep {}: χ(E)={} != dim={}",
                    c.label,
                    chi_e,
                    c.dim
                );
            }
        }
    }

    #[test]
    fn bns_52_318_reports_compound_corep_error() {
        let error = magnetic_irrep_summary_by_bns("52.318")
            .expect_err("compound corepresentation must fail closed in summaries");
        match error {
            MagneticIrrepError::CorepComputationFailed {
                uni,
                sg,
                k_label,
                source_irrep,
                reason,
            } => {
                assert_eq!(uni, 416);
                assert_eq!(sg, 52);
                assert_eq!(k_label, "S");
                assert_eq!(source_irrep, "S1S2");
                assert!(reason.contains("compound corepresentations"));
            }
            other => panic!("unexpected 52.318 summary error: {other:?}"),
        }
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
                    amplified_noise_count += summary
                        .kpoints
                        .iter()
                        .flat_map(|kpoint| &kpoint.coreps)
                        .flat_map(|corep| &corep.characters)
                        .filter(|value| {
                            amplified_noise_targets
                                .iter()
                                .any(|target| (value.abs() - target).abs() < 2e-6)
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
        assert!(
            failures.is_empty(),
            "{} UNI summaries failed",
            failures.len()
        );
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
                    .map_or("N/A".to_string(), |v| format!("{:.0}", v));
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
}
