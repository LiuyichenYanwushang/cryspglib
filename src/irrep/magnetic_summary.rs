//! Magnetic irrep summary — unified entry point for magnetic space group irreps.
//!
//! Given a magnetic space group (by UNI number, BNS label, or explicit operations),
//! this module returns a complete summary of:
//!
//! 1. High-symmetry k-points with labels and fractional coordinates
//! 2. Magnetic co-representations (coreps) classified by Wigner's test
//! 3. Source H-irreps with Miller-Love / Bradley-Cracknell labels
//! 4. Isotropy subgroup candidates (ordinary and magnetic)

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
    /// Could not identify the unitary subgroup H for this UNI.
    MissingUnitarySubgroup(usize),
    /// No irrep data available for this space group.
    MissingIrrepData { sg: u8 },
    /// Corep computation failed for a specific case.
    CorepComputationFailed {
        uni: usize,
        sg: u8,
        k_label: String,
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

/// Complete magnetic irrep summary for a magnetic space group.
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
    /// Magnetic co-representations at this k-point.
    pub coreps: Vec<MagneticCorepSummary>,
}

/// A single magnetic co-representation (corep).
#[derive(Debug, Clone)]
pub struct MagneticCorepSummary {
    /// Label for this corep (e.g. `"GM4-"` for single source, `"Z1Z4 + Z2Z3"` for Type-C pair).
    pub label: String,
    /// Source H-irreps that compose this corep.
    pub source_irreps: Vec<SourceIrrepSummary>,
    /// Wigner classification: A, B, C, or Unsupported.
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

// ── Entry points ───────────────────────────────────────────────────────────────

/// Compute magnetic irrep summary from any input type.
pub fn magnetic_irrep_summary(
    input: MagneticIrrepInput,
) -> Result<MagneticIrrepSummary, MagneticIrrepError> {
    match input {
        MagneticIrrepInput::Uni(uni) => magnetic_irrep_summary_by_uni(uni),
        MagneticIrrepInput::Bns(bns) => magnetic_irrep_summary_by_bns(bns),
        MagneticIrrepInput::Operations { uni, ops } => {
            magnetic_irrep_summary_from_ops(uni, ops)
        }
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

/// Compute magnetic irrep summary from explicit symmetry operations and UNI number.
pub fn magnetic_irrep_summary_from_ops(
    uni: usize,
    mag_ops: &SymmetryOps,
) -> Result<MagneticIrrepSummary, MagneticIrrepError> {
    // 1. Identify H (unitary subgroup) with Hall setting information.
    let h_info = crate::irrep::corep::identify_unitary_subgroup_with_hall(uni)
        .ok_or(MagneticIrrepError::MissingUnitarySubgroup(uni))?;

    // 2. Get MSG metadata.
    let msg = crate::MagneticSpaceGroupType::from_uni(uni);

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

    // 5. Build k-point summaries with magnetic little group metadata.
    let kpoints: Vec<MagneticKPointSummary> = h_kpoints
        .into_iter()
        .map(|kp| {
            let (kx, ky, kz, kd) = kp.coords;
            let mag_lg = crate::irrep::wigner::filter_little_group_with_transform(
                kx,
                ky,
                kz,
                kd,
                mag_ops,
                setting_xf,
                Some(&canonical_translations),
            );
            let unitary_order = mag_lg
                .iter()
                .filter(|&&i| !mag_ops.operations[i].time_reversal)
                .count();
            let antiunitary_order = mag_lg
                .iter()
                .filter(|&&i| mag_ops.operations[i].time_reversal)
                .count();
            MagneticKPointSummary {
                label: kp.label,
                coords: kp.coords,
                little_group_order: mag_lg.len(),
                unitary_order,
                antiunitary_order,
                coreps: Vec::new(),
            }
        })
        .collect();

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
}
