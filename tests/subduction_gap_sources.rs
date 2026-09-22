//! R3 per-operation frame witnesses for the gap-source classification.
//!
//! The offline classifier (`scripts/classify_subduction_gap_sources.py`) works
//! in the archived ISO-IR frame, while the engine works in the subgroup's
//! data-Hall frame.  These tests pin the mapping for the witnesses the R3 card
//! requires, operation by operation, using the engine's own embedding.

use cryspglib::irrep::LabelConvention;
use cryspglib::irrep::isotropy::{self, IsotropySubgroup};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::{
    ExactSeitz, Lattice, Rat, SubductionError, SubgroupEmbedding, Vec3R,
};

/// The archived PIR operations of child #3 (P2, unique axis b) in its own
/// primitive frame: the identity and the twofold rotation about b.
fn archive_child_3() -> Vec<ExactSeitz> {
    let zero = Vec3R::new([Rat::from_integer(0); 3]);
    vec![
        ExactSeitz::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], zero),
        ExactSeitz::new([[-1, 0, 0], [0, 1, 0], [0, 0, -1]], zero),
    ]
}

fn subgroup_by_ordinal(parent: u8, ordinal: usize) -> IsotropySubgroup {
    for record in query::irreps_of(parent) {
        if record.spinor {
            continue;
        }
        for subgroup in isotropy::isotropy_subgroups(parent, record.ml, LabelConvention::Cdml)
            .expect("isotropy table")
        {
            if subgroup.ordinal == ordinal {
                return subgroup;
            }
        }
    }
    panic!("ordinal {ordinal} of space group {parent} not found");
}

fn q_vector(numerators: [i128; 3], denominator: i128) -> Vec3R {
    Vec3R::new(numerators.map(|value| Rat::new(value, denominator).expect("rational")))
}

/// Projective factor system of a finite little group at `q`, as rational turns.
///
/// One representative per rotation is kept; for every ordered pair the lattice
/// translation of the product relative to the representative is reduced modulo
/// the child lattice and paired with `q` exactly as the offline classifier
/// does (`omega = exp(+2 pi i q.L)`).
fn factor_system_turns(
    little: &[ExactSeitz],
    q: &Vec3R,
    lattice: &Lattice,
) -> Result<Vec<Rat>, SubductionError> {
    let mut representatives: Vec<ExactSeitz> = Vec::new();
    for operation in little {
        if !representatives
            .iter()
            .any(|kept| kept.rotation() == operation.rotation())
        {
            representatives.push(*operation);
        }
    }
    let mut turns = Vec::new();
    for first in &representatives {
        for second in &representatives {
            // The product is taken *unreduced*: reducing it first would subtract
            // the very lattice vector whose Bloch phase is the factor we are
            // measuring, and a screw relation S^2 = T(0,0,1) at q = (0,0,1/2)
            // would come out as a trivial phase.
            let product = first.compose(second)?;
            let target = representatives
                .iter()
                .find(|candidate| candidate.rotation() == product.rotation())
                .expect("little group is closed under multiplication");
            let difference = product.translation().checked_sub(target.translation())?;
            assert!(
                lattice.contains(&difference)?,
                "the operation product must differ from its representative by a lattice vector"
            );
            let value = q
                .get(0)
                .checked_mul(difference.get(0))?
                .checked_add(q.get(1).checked_mul(difference.get(1))?)?
                .checked_add(q.get(2).checked_mul(difference.get(2))?)?;
            // Into [0, 1): a phase in turns.  `Rat` is not `Ord`, so the
            // normalisation compares through the exact value converted once.
            let mut value = value;
            while value.to_f64() < 0.0 {
                value = value.checked_add(Rat::from_integer(1))?;
            }
            while value.to_f64() >= 1.0 {
                value = value.checked_sub(Rat::from_integer(1))?;
            }
            turns.push(value);
        }
    }
    turns.sort_by_key(|value| (value.numerator(), value.denominator()));
    Ok(turns)
}

fn little_group_at(
    operations: &[ExactSeitz],
    q: &Vec3R,
    reciprocal: &Lattice,
) -> Result<Vec<ExactSeitz>, SubductionError> {
    let mut little = Vec::new();
    for operation in operations {
        if reciprocal.preserves(operation.rotation(), q)? {
            little.push(*operation);
        }
    }
    Ok(little)
}

/// Ordinal 26 (SG 3 `A1` -> child #3) is the first frozen record whose setting
/// is a non-symmetric change of basis, `U = [[1,2,1],[-1,2,-1],[-1,0,1]] / 2`.
/// Every mapped child operation must still be one of the archived child #3
/// operations, and the little group and factor system at a pinned star must
/// agree with what the offline classifier computes in the archive frame.
#[test]
fn a_non_symmetric_setting_maps_every_child_operation_into_the_engine_frame() {
    let subgroup = subgroup_by_ordinal(3, 26);
    assert_eq!(subgroup.parent_sg, 3);
    assert_eq!(subgroup.record.sg, 3);
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).expect("embedding");
    let lattice = Lattice::integer();
    let archive = archive_child_3();
    let mut mapped = Vec::new();
    for operation in embedding.representatives() {
        let child = embedding
            .transform()
            .unmap_operation(operation)
            .expect("unmap")
            .reduce(&lattice)
            .expect("reduce");
        mapped.push(child);
    }
    assert_eq!(mapped.len(), 2, "child #3 has two operations in this record");
    for operation in &mapped {
        assert!(
            archive.iter().any(|candidate| candidate == operation),
            "mapped operation {operation:?} is not an archived child #3 operation"
        );
    }
    // The mapped set is the archived set, not a subset.
    for operation in &archive {
        assert!(
            mapped.iter().any(|candidate| candidate == operation),
            "archived operation {operation:?} has no mapped counterpart"
        );
    }

    // Pinned star from the R3 manifest for child #3: order two, trivial factor
    // system (SG 3 is symmorphic).
    let q = q_vector([0, 2, 3], 6); // (0, 1/3, 1/2)
    assert_eq!(q.get(1), Rat::new(1, 3).unwrap());
    assert_eq!(q.get(2), Rat::new(1, 2).unwrap());
    let reciprocal = lattice.reciprocal().expect("reciprocal");
    let little = little_group_at(&mapped, &q, &reciprocal).expect("little group");
    assert_eq!(little.len(), 2);
    let turns = factor_system_turns(&little, &q, &lattice).expect("factor system");
    assert_eq!(turns.len(), 4, "two representatives give four ordered pairs");
    assert!(
        turns.iter().all(|value| *value == Rat::ZERO),
        "SG 3 is symmorphic: every phase must be zero, found {turns:?}"
    );
}

/// The factor system must see the lattice vector of a product relation, not its
/// remainder: a screw-like little group with ``S^2 = T(0,0,1)`` at
/// ``q = (0, 0, 1/2)`` has the half-turn phase ``exp(2 pi i * 1/2) = -1``.
#[test]
fn a_screw_relation_produces_a_non_trivial_phase() {
    let lattice = Lattice::integer();
    let zero = Vec3R::new([Rat::from_integer(0); 3]);
    let half = Rat::new(1, 2).unwrap();
    let identity = ExactSeitz::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], zero);
    let screw = ExactSeitz::new(
        [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
        Vec3R::new([Rat::ZERO, Rat::ZERO, half]),
    );
    let q = Vec3R::new([Rat::ZERO, Rat::ZERO, half]);
    let little = vec![identity, screw];
    let turns = factor_system_turns(&little, &q, &lattice).expect("factor system");
    assert_eq!(turns.len(), 4);
    assert_eq!(
        turns.iter().filter(|value| **value == half).count(),
        1,
        "the screw pair must carry the half-turn phase, found {turns:?}"
    );
    assert_eq!(
        turns.iter().filter(|value| value.is_zero()).count(),
        3,
        "the other three ordered pairs are trivial, found {turns:?}"
    );
}
