//! Geometry-derived Type-IV proposals, not a general setting classifier.

use std::collections::BTreeSet;

use super::{
    EmbeddingKind, GeometryError, MagneticGeometry, MagneticOperation, MagneticOperationSet, Mat3I,
    Mat3R, ParentEmbedding, Rat, Vec3R, basis_in_parent_conventional,
    canonical_translation_checked, exact_embeddings_for_map, exact_fraction,
    verified_source_primitive_basis, verify_group,
};
use crate::SymError;
use crate::irrep::subduction::SubgroupEmbedding;
use crate::mathfunc::mat_get_invariant_lattice as invariant_lattice;
use crate::symmetry::MagneticSymmetry;

const IDENTITY: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// A valid auxiliary candidate which the floating identifier could not classify.
/// This is neither a rejection of that candidate nor an absence certificate.
#[derive(Debug, Clone)]
pub struct AntitranslationProposalFailure {
    /// Half-lattice representative in the record's primitive coordinates.
    pub translation: [Rat; 3],
    pub error: SymError,
}

/// Verified witnesses and unresolved *proposal generation*, kept separate.
#[derive(Debug, Clone, Default)]
pub struct AntitranslationSearch {
    /// Every returned entry passes full original-source, lattice and origin checks.
    pub embeddings: Vec<ParentEmbedding>,
    pub unresolved: Vec<AntitranslationProposalFailure>,
}

/// Derive Type-IV maps from an ordinary H point subgroup and antitranslations.
///
/// B and C are the record and verified full source primitive column bases. In
/// each primitive frame form an AUXILIARY symmorphic group from the unitary
/// rotations, retaining the antitranslation but discarding H's translations.
/// On the target side enumerate the seven nonzero elements of (1/2 L)/L,
/// retaining those in the parent lattice and fixed modulo L by H's rotations.
/// The ordinary H embedding supplies only these rotations, NOT its period lattice.
///
/// Identify both auxiliary groups in the same unrestricted canonical frame.
/// Matching identifications propose A = B D_target^-1 D_source C^-1. Neither
/// auxiliary group is substituted for the actual magnetic group: every map is
/// checked against its ORIGINAL translations, full lattice, and all common origin
/// components. Thus nonsymmorphic H is allowed, but its screw/glide constraints
/// are not lost. Floating identifications only propose maps, never certify them.
///
/// Only an actual Type-IV antitranslation is used; Type I/II/III return an empty
/// search. An empty search never disproves other maps or H point subgroups.
/// Unclassified target proposals remain in `unresolved` even if a different
/// proposal succeeds. Preparation, provenance and exact-arithmetic errors remain
/// errors. No R24 family, OP-stabilizer claim or representation scope is changed.
pub fn embed_via_antitranslations(
    geometry: &MagneticGeometry,
    set: &MagneticOperationSet,
    ordinary: &SubgroupEmbedding,
) -> Result<AntitranslationSearch, GeometryError> {
    let unavailable = |reason: String| GeometryError::AntitranslationBridgeUnavailable {
        uni: geometry.uni,
        reason,
    };
    if geometry.uni != set.uni || geometry.parent_sg != ordinary.parent_sg() {
        return Err(unavailable(
            "magnetic UNI or parent identity differs".into(),
        ));
    }
    let c = verified_source_primitive_basis(geometry, set)?;
    let h = crate::irrep::corep::cached_unitary_subgroup_with_hall(set.uni)
        .ok_or_else(|| unavailable("unitary subgroup identification unavailable".into()))?;
    if h.sg != usize::from(ordinary.subgroup_sg()) {
        return Err(unavailable(
            "ordinary subgroup differs from identified H".into(),
        ));
    }
    let mut out = AntitranslationSearch::default();
    let Some(a) = set
        .classes
        .iter()
        .find(|o| o.time_reversal && o.rotation == IDENTITY)
    else {
        return Ok(out);
    };
    let ci = c.inverse().map_err(|_| GeometryError::Arithmetic)?;
    let alpha = ci
        .checked_mul_vector(&Vec3R::new(a.translation))
        .map_err(|_| GeometryError::Arithmetic)?;
    let alpha = canonical_translation_checked(*alpha.as_array(), &[])
        .map_err(|e| unavailable(e.to_string()))?;
    if alpha == [Rat::ZERO; 3] || geometry.basis.iter().flatten().any(|r| !r.is_integer()) {
        return Ok(out);
    }
    let b = Mat3R::new(basis_in_parent_conventional(geometry)?).transpose();
    let bi = b.inverse().map_err(|_| GeometryError::Arithmetic)?;
    let mut target_rotations = BTreeSet::new();
    for op in ordinary.operations() {
        let r = bi
            .checked_mul(&Mat3R::from_ints(op.rotation()))
            .and_then(|r| r.checked_mul(&b))
            .map_err(|_| GeometryError::Arithmetic)?;
        if !r.is_integral() {
            return Ok(out);
        }
        target_rotations.insert(r.to_int_matrix().map_err(|_| GeometryError::Arithmetic)?);
    }
    let target_rotations: Vec<_> = target_rotations.into_iter().collect();
    let mut source_rotations = BTreeSet::new();
    for op in set.classes.iter().filter(|o| !o.time_reversal) {
        let r = ci
            .checked_mul(&Mat3R::from_ints(op.rotation))
            .and_then(|r| r.checked_mul(&c))
            .map_err(|_| GeometryError::Arithmetic)?;
        source_rotations.insert(r.to_int_matrix().map_err(|_| GeometryError::Arithmetic)?);
    }
    let source_rotations: Vec<_> = source_rotations.into_iter().collect();
    let source = virtual_group(&source_rotations, alpha).map_err(&unavailable)?;
    let source_metric = invariant_lattice(&source_rotations)
        .ok_or_else(|| unavailable("source invariant metric unavailable".into()))?;
    let source =
        crate::magnetic_spacegroup::identify_with_parent_hall(&source_metric, &source, None, 1e-5)
            .map_err(|e| unavailable(format!("source auxiliary identification: {e:?}")))?;
    let ds = rational_frame(source.transformation_matrix)
        .ok_or_else(|| unavailable("source auxiliary frame off the rational grid".into()))?;
    let target_metric = invariant_lattice(&target_rotations)
        .ok_or_else(|| unavailable("target invariant metric unavailable".into()))?;
    for code in 1..8 {
        let alpha =
            std::array::from_fn(|axis| Rat::new((code >> axis) & 1, 2).expect("denominator two"));
        let parent = Mat3R::new(geometry.basis)
            .transpose()
            .checked_mul_vector(&Vec3R::new(alpha))
            .map_err(|_| GeometryError::Arithmetic)?;
        if parent.as_array().iter().any(|r| !r.is_integer()) {
            continue;
        }
        let mut fixed = true;
        for &r in &target_rotations {
            let image = Mat3R::from_ints(r)
                .checked_mul_vector(&Vec3R::new(alpha))
                .map_err(|_| GeometryError::Arithmetic)?;
            for (image, original) in image.as_array().iter().zip(alpha) {
                fixed &= image
                    .checked_sub(original)
                    .map_err(|_| GeometryError::Arithmetic)?
                    .is_integer();
            }
        }
        if !fixed {
            continue;
        }
        let target = virtual_group(&target_rotations, alpha).map_err(&unavailable)?;
        let target = match crate::magnetic_spacegroup::identify_with_parent_hall(
            &target_metric,
            &target,
            None,
            1e-5,
        ) {
            Ok(target) => target,
            Err(error) => {
                out.unresolved.push(AntitranslationProposalFailure {
                    translation: alpha,
                    error,
                });
                continue;
            }
        };
        if (target.uni_number, target.hall_number) != (source.uni_number, source.hall_number) {
            continue;
        }
        let dt = rational_frame(target.transformation_matrix)
            .ok_or_else(|| unavailable("target auxiliary frame off the rational grid".into()))?;
        let map = dt
            .inverse()
            .and_then(|dt| b.checked_mul(&dt))
            .and_then(|m| m.checked_mul(&ds))
            .and_then(|m| m.checked_mul(&ci))
            .map_err(|_| GeometryError::Arithmetic)?;
        for entry in exact_embeddings_for_map(geometry, set, map, EmbeddingKind::Antitranslation)? {
            if !out.embeddings.contains(&entry) {
                out.embeddings.push(entry);
            }
        }
    }
    Ok(out)
}

pub(super) fn rational_frame(matrix: [[f64; 3]; 3]) -> Option<Mat3R> {
    let mut out = [[Rat::ZERO; 3]; 3];
    for (row, values) in out.iter_mut().enumerate() {
        for (column, value) in values.iter_mut().enumerate() {
            *value = exact_fraction(matrix[row][column])?;
        }
    }
    Some(Mat3R::new(out))
}

fn virtual_group(rotations: &[Mat3I], alpha: [Rat; 3]) -> Result<MagneticSymmetry, String> {
    let mut ops = Vec::new();
    for time_reversal in [false, true] {
        for &rotation in rotations {
            ops.push(MagneticOperation {
                rotation,
                time_reversal,
                translation: if time_reversal { alpha } else { [Rat::ZERO; 3] },
            });
        }
    }
    verified_auxiliary_group(ops)
}

// Shared by the translation and point-extension proposals. Preserve operation
// order, and never treat a failed group/arithmetic check as a rejected setting.
pub(super) fn verified_auxiliary_group(
    ops: Vec<MagneticOperation>,
) -> Result<MagneticSymmetry, String> {
    let n = ops.iter().filter(|o| !o.time_reversal).count();
    let a = ops.len() - n;
    verify_group(&MagneticOperationSet {
        uni: 0,
        operations: ops.clone(),
        classes: ops.clone(),
        unitary_count: n,
        antiunitary_count: a,
        unitary_classes: n,
        antiunitary_classes: a,
    })
    .map_err(|e| e.to_string())?;
    let mut symmetry = MagneticSymmetry::with_capacity(ops.len());
    for op in ops {
        symmetry.push(
            op.rotation,
            op.translation.map(|r| r.to_f64()),
            op.time_reversal,
        );
    }
    Ok(symmetry)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_parent_metrics_are_positive_and_rotation_invariant() {
        use crate::mathfunc::{
            mat_multiply_matrix_d3 as mul, mat_transpose_matrix_d3 as transpose,
        };
        for sg in 1..=230 {
            let rotations = super::super::parent_frame(sg).unwrap().rotations;
            let lattice = invariant_lattice(&rotations).unwrap();
            let gram = mul(&transpose(&lattice), &lattice);
            for r in rotations {
                let r = r.map(|row| row.map(f64::from));
                let transformed = mul(&mul(&transpose(&r), &gram), &r);
                for i in 0..3 {
                    for j in 0..3 {
                        assert!((gram[i][j] - transformed[i][j]).abs() < 1e-10);
                    }
                }
            }
        }
        assert!(invariant_lattice(&[]).is_none());
    }

    #[test]
    fn frame_conversion_rejects_nonfinite_and_off_grid_values() {
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.123456789] {
            let mut frame = [[0.0; 3]; 3];
            frame[0][0] = invalid;
            assert!(rational_frame(frame).is_none());
        }
    }
}
