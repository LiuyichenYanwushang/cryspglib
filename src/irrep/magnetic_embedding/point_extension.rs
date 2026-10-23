//! Type-III proposals from finite, lattice-preserving parent point extensions.

use super::antitranslation::{rational_frame, verified_auxiliary_group};
use super::{
    EmbeddingKind, GeometryError, MagneticGeometry, MagneticOperation, MagneticOperationSet, Mat3I,
    Mat3R, ParentEmbedding, Rat, basis_in_parent_conventional, exact_embeddings_for_map,
    parent_frame, verified_source_primitive_basis,
};
use crate::SymError;
use crate::irrep::subduction::SubgroupEmbedding;
use crate::mathfunc::mat_get_invariant_lattice as invariant_lattice;
use std::collections::BTreeSet;

const IDENTITY: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// A valid graded point extension which the floating identifier did not resolve.
#[derive(Debug, Clone)]
pub struct PointExtensionProposalFailure {
    /// Representative of its antiunitary point coset, in record primitive coordinates.
    pub anti_rotation: Mat3I,
    pub error: SymError,
}

/// Proven embeddings and unresolved proposal generation, never a nonexistence certificate.
#[derive(Debug, Clone, Default)]
pub struct PointExtensionSearch {
    pub embeddings: Vec<ParentEmbedding>,
    pub unresolved: Vec<PointExtensionProposalFailure>,
}

/// Propose Type-III embeddings from the supplied ordinary H's point subgroup.
///
/// In record primitive coordinates, enumerate ALL parent rotations preserving
/// that lattice. For a outside H, H union aH is a graded point group exactly
/// when a^2 belongs to H and aH=Ha. Enumerate distinct such cosets, not a guessed
/// matrix family. The ordinary embedding supplies rotations, not its period.
///
/// Set the translations of both the source and target auxiliary groups to zero,
/// retaining their time-reversal flags. Matching unrestricted identifications
/// propose A=B D_target^-1 D_source C^-1, where C is the verified full source
/// primitive basis and B is the record basis. Every proposal is then checked
/// with the ORIGINAL magnetic translations, exact full lattice equality, and
/// ALL common-origin components. Auxiliary groups and floating frames are never
/// accepted in place of the actual nonsymmorphic source.
///
/// Only Type III applies; Type I/II/IV return an empty search after identity/data
/// validation. Identifier failures remain in `unresolved` even if another coset
/// succeeds. Preparation and exact-arithmetic failures remain errors. Neither
/// an empty search nor a failed auxiliary frame proves nonexistence: different
/// normalizer choices or H point subgroups may work. This does not certify an OP
/// stabilizer or extend the representation domain.
pub fn embed_via_point_extensions(
    geometry: &MagneticGeometry,
    set: &MagneticOperationSet,
    ordinary: &SubgroupEmbedding,
) -> Result<PointExtensionSearch, GeometryError> {
    let unavailable = |reason: String| GeometryError::PointExtensionUnavailable {
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
    let mut out = PointExtensionSearch::default();
    if set.antiunitary_classes == 0
        || set
            .classes
            .iter()
            .any(|o| o.time_reversal && o.rotation == IDENTITY)
        || geometry.basis.iter().flatten().any(|r| !r.is_integer())
    {
        return Ok(out);
    }
    let b = Mat3R::new(basis_in_parent_conventional(geometry)?).transpose();
    let bi = b.inverse().map_err(|_| GeometryError::Arithmetic)?;
    let ci = c.inverse().map_err(|_| GeometryError::Arithmetic)?;
    let mut target_h = BTreeSet::new();
    for op in ordinary.operations() {
        let r = in_basis(op.rotation(), &b, &bi)?;
        if !r.is_integral() {
            return Ok(out);
        }
        target_h.insert(r.to_int_matrix().map_err(|_| GeometryError::Arithmetic)?);
    }
    let source_ops = set
        .classes
        .iter()
        .map(|op| {
            Ok(MagneticOperation {
                rotation: in_basis(op.rotation, &c, &ci)?
                    .to_int_matrix()
                    .map_err(|_| GeometryError::Arithmetic)?,
                translation: [Rat::ZERO; 3],
                time_reversal: op.time_reversal,
            })
        })
        .collect::<Result<Vec<_>, GeometryError>>()?;
    // ponytail: one canonical frame, not a normalizer search. Keep failures
    // visible; broader setting coverage needs exact normalizer witnesses.
    let source = verified_auxiliary_group(source_ops).map_err(&unavailable)?;
    let lattice = invariant_lattice(&source.rot)
        .ok_or_else(|| unavailable("source invariant metric unavailable".into()))?;
    let source =
        crate::magnetic_spacegroup::identify_with_parent_hall(&lattice, &source, None, 1e-5)
            .map_err(|e| unavailable(format!("source auxiliary identification: {e:?}")))?;
    let ds = rational_frame(source.transformation_matrix)
        .ok_or_else(|| unavailable("source auxiliary frame off the rational grid".into()))?;
    let mut parent = BTreeSet::new();
    for r in parent_frame(geometry.parent_sg)?.rotations {
        let r = in_basis(r, &b, &bi)?;
        if r.is_integral() {
            parent.insert(r.to_int_matrix().map_err(|_| GeometryError::Arithmetic)?);
        }
    }
    let mut seen = BTreeSet::new();
    for a in parent {
        let Some(coset) = extension_coset(&target_h, a)? else {
            continue;
        };
        if !seen.insert(coset.clone()) {
            continue;
        }
        let mut operations: Vec<_> = target_h
            .iter()
            .map(|&r| (r, false))
            .chain(coset.into_iter().map(|r| (r, true)))
            .collect();
        operations.sort_unstable();
        let operations = operations
            .into_iter()
            .map(|(rotation, time_reversal)| MagneticOperation {
                rotation,
                translation: [Rat::ZERO; 3],
                time_reversal,
            })
            .collect();
        let target = verified_auxiliary_group(operations).map_err(&unavailable)?;
        let lattice = invariant_lattice(&target.rot)
            .ok_or_else(|| unavailable("target invariant metric unavailable".into()))?;
        let target = match crate::magnetic_spacegroup::identify_with_parent_hall(
            &lattice, &target, None, 1e-5,
        ) {
            Ok(target) => target,
            Err(error) => {
                out.unresolved.push(PointExtensionProposalFailure {
                    anti_rotation: a,
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
        for entry in exact_embeddings_for_map(geometry, set, map, EmbeddingKind::PointExtension)? {
            if !out.embeddings.contains(&entry) {
                out.embeddings.push(entry);
            }
        }
    }
    Ok(out)
}

fn in_basis(r: Mat3I, b: &Mat3R, bi: &Mat3R) -> Result<Mat3R, GeometryError> {
    bi.checked_mul(&Mat3R::from_ints(r))
        .and_then(|r| r.checked_mul(b))
        .map_err(|_| GeometryError::Arithmetic)
}
fn multiply(a: Mat3I, b: Mat3I) -> Result<Mat3I, GeometryError> {
    Mat3R::from_ints(a)
        .checked_mul(&Mat3R::from_ints(b))
        .and_then(|r| r.to_int_matrix())
        .map_err(|_| GeometryError::Arithmetic)
}

// H is the ordinary embedding's full finite point subgroup. Algebraic
// inadmissibility is checked explicitly, never by swallowing verify_group errors.
fn extension_coset(h: &BTreeSet<Mat3I>, a: Mat3I) -> Result<Option<Vec<Mat3I>>, GeometryError> {
    if h.contains(&a) || !h.contains(&multiply(a, a)?) {
        return Ok(None);
    }
    let left = h
        .iter()
        .map(|&r| multiply(a, r))
        .collect::<Result<BTreeSet<_>, _>>()?;
    let right = h
        .iter()
        .map(|&r| multiply(r, a))
        .collect::<Result<BTreeSet<_>, _>>()?;
    Ok((left == right).then(|| left.into_iter().collect()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mathfunc::mat_multiply_matrix_i3 as product;
    use std::collections::BTreeMap;

    #[test]
    fn coset_conditions_agree_with_independent_coloured_product_tables() {
        let parent = parent_frame(221).unwrap().rotations;
        let mut subgroups = BTreeSet::new();
        for &generator in &parent {
            let mut h = BTreeSet::from([IDENTITY]);
            let mut r = generator;
            while h.insert(r) {
                r = product(&r, &generator);
            }
            subgroups.insert(h);
        }
        subgroups.insert(parent.iter().copied().collect());
        subgroups.insert(
            parent
                .iter()
                .copied()
                .filter(|&r| Mat3R::from_ints(r).determinant().unwrap() == Rat::ONE)
                .collect(),
        );
        let (mut accepted, mut rejected) = (0, 0);
        for h in &subgroups {
            for &a in &parent {
                let members: BTreeMap<_, _> = h
                    .iter()
                    .map(|&r| (r, false))
                    .chain(h.iter().map(|&r| (product(&a, &r), true)))
                    .collect();
                let expected = members.len() == 2 * h.len()
                    && members.iter().all(|(r, t)| {
                        members
                            .iter()
                            .all(|(s, u)| members.get(&product(r, s)) == Some(&(t ^ u)))
                    });
                assert_eq!(extension_coset(h, a).unwrap().is_some(), expected);
                if expected {
                    accepted += 1;
                } else {
                    rejected += 1;
                }
            }
        }
        assert!(accepted > 0 && rejected > 0);
        // Order-four antiunitary generators are allowed: their square need
        // only be in H, not the identity.
        let a = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        let h = BTreeSet::from([IDENTITY, product(&a, &a)]);
        assert_ne!(product(&a, &a), IDENTITY);
        assert!(extension_coset(&h, a).unwrap().is_some());
    }
}
