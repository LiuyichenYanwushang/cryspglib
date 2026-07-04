//! Magnetic irrep summary — unified entry point for magnetic space group irreps.
//!
//! Given a magnetic space group (by UNI number, BNS label, or explicit operations),
//! this module returns a complete summary of:
//!
//! 1. High-symmetry k-points with labels and fractional coordinates
//! 2. Magnetic co-representations (coreps) classified by Wigner's test
//! 3. Source H-irreps with Miller-Love / Bradley-Cracknell labels
//! 4. Isotropy subgroup candidates (ordinary and magnetic)

use std::collections::BTreeSet;

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
/// Type-C coreps from two antiunitary-conjugate H-irreps produce identical
/// character tables.  We merge them into a single entry with combined labels
/// and source_irreps, rather than showing the same magnetic corep twice.
fn dedup_coreps(coreps: Vec<MagneticCorepSummary>) -> Vec<MagneticCorepSummary> {
    let mut groups: Vec<Vec<MagneticCorepSummary>> = Vec::new();

    for c in coreps {
        let key = (
            c.corep_type,
            c.dim,
            round_chars(&c.characters),
            c.timerev.clone(),
        );
        // Find existing group with matching key.
        let found = groups.iter_mut().find(|g| {
            let first = &g[0];
            key.0 == first.corep_type
                && key.1 == first.dim
                && key.2 == round_chars(&first.characters)
                && key.3 == first.timerev
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
            let mut extra_sources: Vec<_> = group
                .into_iter()
                .flat_map(|c| c.source_irreps)
                .collect();
            merged.source_irreps.append(&mut extra_sources);
            // Build combined label: sort source ML labels and join with " + ".
            let mut labels: Vec<&str> = merged
                .source_irreps
                .iter()
                .map(|s| s.ml)
                .collect();
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
                let ir = h_irreps
                    .iter()
                    .find(|r| r.sg == src.sg && r.ml == src.ml);

                let relation = if src.spinor {
                    IsotropyCandidateRelation::SpinorNoIsotropyData
                } else if multi_source {
                    IsotropyCandidateRelation::TypeCPairedSource
                } else if ir.map_or(false, |r| r.cir_component_count() > 0) {
                    IsotropyCandidateRelation::CompoundSource
                } else {
                    IsotropyCandidateRelation::DirectSourceIrrep
                };

                let (ordinary, magnetic) = match ir {
                    Some(rec) if !src.spinor => (
                        rec.subgroups().to_vec(),
                        rec.magnetic_subgroups().to_vec(),
                    ),
                    _ => (Vec::new(), Vec::new()),
                };

                if !ordinary.is_empty() || !magnetic.is_empty() || !src.spinor {
                    candidates.push(CorepIsotropyCandidate {
                        source_ml: src.ml,
                        ordinary,
                        magnetic,
                        relation,
                    });
                }
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
            if let Some(existing) = result.iter_mut().find(|e| {
                e.source_ml == cand.source_ml && e.relation == cand.relation
            }) {
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
        cand.ordinary.retain(|s| {
            ord_seen.insert((s.sg, s.symbol, s.direction, s.domains, s.arms))
        });

        // Dedup magnetic subgroups by (mag_sg, bns_label, direction).
        let mut mag_seen: BTreeSet<(usize, &str, &str)> = BTreeSet::new();
        cand.magnetic.retain(|s| {
            mag_seen.insert((s.mag_sg, s.bns_label, s.direction))
        });
    }

    result
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

    // 5. Get H's irreps for corep computation.
    let h_irreps = crate::irrep::query::irreps_of(h_info.sg as u8);

    // 6. Build k-point summaries with little group metadata and coreps.
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

            // Compute coreps for each irrep at this k-point.
            let raw_coreps: Vec<MagneticCorepSummary> = kp
                .irreps
                .iter()
                .filter_map(|&idx| {
                    let ir = &h_irreps[idx];
                    let c = crate::irrep::corep::compute_corepresentation(ir, uni, mag_ops)?;
                    Some(MagneticCorepSummary {
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
                    })
                })
                .collect();
            let coreps = dedup_coreps(raw_coreps);
            let coreps = attach_isotropy_candidates(coreps, h_irreps);

            MagneticKPointSummary {
                label: kp.label,
                coords: kp.coords,
                little_group_order: mag_lg.len(),
                unitary_order,
                antiunitary_order,
                coreps,
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

    for c in &kpoint.coreps {
        let src_labels: Vec<&str> = c.source_irreps.iter().map(|s| s.ml).collect();
        lines.push(format!(
            "  {}  type={:?}  source={:?}  dim={}  src=[{}]",
            c.label,
            c.corep_type,
            c.source,
            c.dim,
            src_labels.join(", ")
        ));
        // Show first few characters.
        let char_preview: Vec<String> = c
            .characters
            .iter()
            .take(6)
            .map(|&ch| {
                if ch.abs() < 1e-12 {
                    "0".to_string()
                } else {
                    format!("{:.2}", ch)
                }
            })
            .collect();
        let char_str = if c.characters.len() > 6 {
            format!("[{}...]", char_preview.join(", "))
        } else {
            format!("[{}]", char_preview.join(", "))
        };
        lines.push(format!("    chars: {}", char_str));

        // Isotropy candidates summary.
        if !c.isotropy_candidates.is_empty() {
            for ic in &c.isotropy_candidates {
                let n_ord = ic.ordinary.len();
                let n_mag = ic.magnetic.len();
                if n_ord > 0 || n_mag > 0 {
                    lines.push(format!(
                        "    isotropy ({} {:?}): {} ordinary + {} magnetic subgroups",
                        ic.source_ml, ic.relation, n_ord, n_mag
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
        // BNS 128.406 (UNI 1066): verify isotropy candidates are populated.
        let s = magnetic_irrep_summary_by_bns("128.406").unwrap();
        // Z1Z4 is a compound irrep with CIR components → CompoundSource.
        let z_kp = s.kpoints.iter().find(|k| k.label == "Z").unwrap();
        for c in &z_kp.coreps {
            if c.label == "Z1Z4" {
                assert!(
                    !c.isotropy_candidates.is_empty(),
                    "Z1Z4 should have isotropy candidates"
                );
                let has_compound = c
                    .isotropy_candidates
                    .iter()
                    .any(|ic| ic.relation == IsotropyCandidateRelation::CompoundSource);
                assert!(has_compound, "Z1Z4 should have CompoundSource relation");
            }
        }
    }

    #[test]
    fn coreps_at_z_for_128_406() {
        // BNS 128.406 (UNI 1066) at Z: verified against Bilbao BCS.
        // Expected: Z1Z4 (C), Z2Z3 (C), Z5 (A, currently Unsupported),
        // Z6 (C spinor), Z7 (C spinor).
        let s = magnetic_irrep_summary_by_bns("128.406").unwrap();
        let z_kp = s
            .kpoints
            .iter()
            .find(|kp| kp.label == "Z")
            .expect("should have Z k-point");
        assert!(!z_kp.coreps.is_empty(), "should have coreps at Z");

        // Type-C coreps come from the scalar CIR path (Z1Z4, Z2Z3).
        let has_type_c = z_kp
            .coreps
            .iter()
            .any(|c| c.corep_type == crate::irrep::corep::CorepType::C);
        assert!(has_type_c, "should have at least one Type-C corep at Z");

        // Every corep should have non-empty source_irreps.
        for c in &z_kp.coreps {
            assert!(!c.source_irreps.is_empty(), "corep {} has no source irrep", c.label);
        }
    }

    #[test]
    fn no_duplicate_coreps_at_any_kpoint() {
        // For 128.406, each k-point should have no duplicate coreps
        // (duplicated = same type + dim + rounded characters + timerev).
        let s = magnetic_irrep_summary_by_bns("128.406").unwrap();
        for kp in &s.kpoints {
            for c in &kp.coreps {
                assert!(
                    !c.source_irreps.is_empty(),
                    "k-point {}: corep has no source irrep",
                    kp.label
                );
            }
            // Check no duplicates.
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
                        kp.label, ci.label, cj.label, i, j
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
}
