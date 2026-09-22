//! R4 batch 1: constructed targets for every trivial-little-co-group child star.
//!
//! R3 classified 215 of the 899 remaining gap groups as
//! `analytic_general_position`: every arm of their star has a trivial little
//! co-group, so the child irrep at that point is the one-dimensional Bloch phase
//! `D(T_L) = exp(+2 pi i q.L)` and no archived character data is needed.  R2 had
//! wired that route for child #1 only; this batch replaces the child-number
//! guard with the geometric one, so a child with a **non-trivial point group**
//! is answered too — as long as the star's little co-group is trivial.
//!
//! The two contexts below are the batch's own boundary witnesses, both SG 225 ->
//! #8 (Cm, C-centred, point group `m`):
//!
//! * ordinal 13345 (`C2` direction): all 31 ordinary probes are answered
//!   completely, including folded stars that are nowhere in the discrete #8
//!   table;
//! * ordinal 13346 (`C1` direction): probes whose stars include points on the
//!   parametric `B` line are answered by R4 batch 2a, the exact one-dimensional
//!   projective catalogue of their two-fold little co-group;
//! * ordinal 3988 (SG 109 -> #43): a four-element little co-group with a single
//!   omega-regular class, answered by R4 batch 2b's gated projective table.
//!
//! Every reported block is validated by the production entry point itself
//! (dimension conservation, integral multiplicities, per-operation
//! reconstruction), and the identity multiplicities are additionally compared
//! with the pinned `isotropy_subduce_*` frequency, which is the only quantity the
//! frozen table constrains.

use cryspglib::irrep::isotropy::{IsotropySubgroup, isotropy_subgroups};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::star::decompose::{
    subduce_full_star_with_embedding, trivial_content_with_embedding,
};
use cryspglib::irrep::subduction::{SubgroupEmbedding, SubductionComponent};
use cryspglib::irrep::types::{IrrepRecord, IrrepSourceIdentity};
use cryspglib::irrep::LabelConvention;

/// The isotropy record with this ordinal, with its validated embedding.
fn ordinal_context(parent: u8, ordinal: usize) -> (IsotropySubgroup, SubgroupEmbedding) {
    for record in query::irreps_of(parent).iter().filter(|record| !record.spinor) {
        for subgroup in
            isotropy_subgroups(parent, record.ml, LabelConvention::Cdml).unwrap_or_else(|error| {
                panic!("SG {parent} {}: isotropy subgroups: {error}", record.ml)
            })
        {
            if subgroup.ordinal == ordinal {
                let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup)
                    .unwrap_or_else(|error| panic!("ordinal {ordinal}: embedding: {error}"));
                return (subgroup, embedding);
            }
        }
    }
    panic!("space group {parent} has no isotropy record {ordinal}");
}

fn probe(sg: u8, ml: &str) -> &'static IrrepRecord {
    query::irreps_of(sg)
        .iter()
        .find(|record| record.ml == ml && !record.spinor)
        .unwrap_or_else(|| panic!("SG {sg} has no scalar irrep {ml}"))
}

fn ordinary_probes(sg: u8) -> impl Iterator<Item = &'static IrrepRecord> {
    query::irreps_of(sg).iter().filter(|record| {
        matches!(
            record.source_identity(),
            IrrepSourceIdentity::OrdinaryScalar { .. }
        )
    })
}

fn stored_frequency(subgroup: &IsotropySubgroup, ml: &str) -> u32 {
    subgroup
        .identity_subduction()
        .unwrap_or_else(|error| panic!("ordinal {}: frequency: {error}", subgroup.ordinal))
        .iter()
        .filter(|entry| entry.parent_ml == ml)
        .map(|entry| entry.frequency)
        .next()
        .unwrap_or(0)
        .into()
}

/// Ordinal 13345 (SG 225 -> #8, `C2` direction): the constructed source answers
/// every folded star, so all ordinary probes return a complete decomposition
/// whose trivial content still matches the pinned frequency.
#[test]
fn a_trivial_co_group_context_is_answered_completely() {
    let (subgroup, embedding) = ordinal_context(225, 13345);
    assert_eq!(embedding.subgroup_sg(), 8);
    let trivial = query::irreps_of(8)
        .iter()
        .find(|record| {
            !record.spinor
                && record.k_vector().numerators == [0, 0, 0]
                && record.dim == 1
        })
        .expect("#8 has a trivial Gamma row");

    let mut answered = 0usize;
    let mut constructed_blocks = 0usize;
    for probe in ordinary_probes(225) {
        let result = subduce_full_star_with_embedding(&subgroup, &embedding, probe)
            .unwrap_or_else(|error| panic!("ordinal 13345, probe {}: {error}", probe.ml));
        // The only quantity the pinned table constrains.
        let content = trivial_content_with_embedding(&subgroup, &embedding, probe)
            .unwrap_or_else(|error| panic!("ordinal 13345, probe {}: {error}", probe.ml));
        assert_eq!(content.total, content.by_label);
        let full_total: u32 = result
            .blocks()
            .iter()
            .map(|block| block.multiplicity(trivial.ml))
            .sum();
        assert_eq!(
            full_total,
            stored_frequency(&subgroup, probe.ml),
            "ordinal 13345, probe {}: stored frequency",
            probe.ml
        );
        for block in result.blocks() {
            for target in block.targets() {
                if matches!(
                    target.component,
                    SubductionComponent::Constructed { .. }
                ) {
                    assert!(
                        target.ml.is_none() && target.row_ml.is_none() && target.irnumber.is_none(),
                        "a constructed target must not borrow a stored identity"
                    );
                    assert_eq!(target.dimension, 1, "a trivial co-group target is one-dimensional");
                    assert_eq!(
                        block.constructed_multiplicity(target.component),
                        target.multiplicity
                    );
                    constructed_blocks += 1;
                }
            }
        }
        answered += 1;
    }
    assert_eq!(answered, 31, "every ordinary SG 225 probe of this context");
    assert!(
        constructed_blocks > 0,
        "the batch must actually construct targets here"
    );
}

/// The decisive multi-arm witness of this batch: SG 225 `W1` -> #8 (ordinal
/// 13345) folds onto the two points `(0, 1/2, 1/2)` and `(0, 3/2, 1/2)`, which
/// differ by `(0, 1, 0)` — *not* a #8 reciprocal lattice vector (the C-centring
/// extinction), so they are two arms of one child star.  The star carries one
/// constructed target of dimension one, and its induced character is zero on a
/// rotation that moves the arm while the identity is the star size.
#[test]
fn a_constructed_target_is_induced_over_the_whole_child_star() {
    let (subgroup, embedding) = ordinal_context(225, 13345);
    let result = subduce_full_star_with_embedding(&subgroup, &embedding, probe(225, "W1"))
        .expect("ordinal 13345 W1 is answered by constructed targets");
    let two_arm = result
        .blocks()
        .iter()
        .find(|block| block.star_size() == 2)
        .expect("W1 folds onto a two-arm child star");
    assert_eq!(two_arm.little_dimension(), 1);
    assert_eq!(two_arm.block_dimension(), 2);
    assert_eq!(two_arm.targets().len(), 1);
    let target = &two_arm.targets()[0];
    assert!(matches!(
        target.component,
        SubductionComponent::Constructed { .. }
    ));
    assert_eq!(target.dimension, 1);
    assert_eq!(target.multiplicity, 1);
    // The identity character of the induced star is the star size; the batch's
    // reconstruction check inside the entry point already compared every
    // subgroup operation against the parent full-star character.
    let (parent_characters, reconstructed) = result.reconstruction();
    assert_eq!(parent_characters.len(), reconstructed.len());
    assert!(
        parent_characters
            .iter()
            .zip(reconstructed)
            .all(|(expected, found)| (expected - found).norm() < 1e-9),
        "the reported child irreps must rebuild the subduced parent character"
    );
    // One component induced over two arms: the block dimension is exactly
    // `star size x little dimension`, which a per-arm duplicate could not give.
    assert_eq!(
        two_arm.little_dimension() * u32::try_from(two_arm.star_size()).unwrap(),
        two_arm.block_dimension()
    );
}

/// The parametric `B` line of the same child: its two-fold little co-group has
/// a cocycle that is a coboundary, so R4 batch 2a answers every folded star with
/// the exact one-dimensional projective catalogue while the pinned identity
/// frequency still comes out of the same run.
#[test]
fn a_parametric_co_group_is_answered_by_the_one_dimensional_catalogue() {
    let (subgroup, embedding) = ordinal_context(225, 13346);
    assert_eq!(embedding.subgroup_sg(), 8);
    let trivial = query::irreps_of(8)
        .iter()
        .find(|record| {
            !record.spinor && record.k_vector().numerators == [0, 0, 0] && record.dim == 1
        })
        .expect("#8 has a trivial Gamma row");
    for ml in ["W1", "W2", "W3", "W4", "W5"] {
        let probe = probe(225, ml);
        let result = subduce_full_star_with_embedding(&subgroup, &embedding, probe)
            .unwrap_or_else(|error| panic!("probe {ml}: {error}"));
        let content = trivial_content_with_embedding(&subgroup, &embedding, probe)
            .unwrap_or_else(|error| panic!("probe {ml}: {error}"));
        assert_eq!(content.total, content.by_label);
        assert_eq!(content.total, stored_frequency(&subgroup, ml));
        let full_total: u32 = result
            .blocks()
            .iter()
            .map(|block| block.multiplicity(trivial.ml))
            .sum();
        assert_eq!(full_total, content.total, "probe {ml}");
    }
}

/// R4 batch 2b: ordinal 3988 (SG 109 -> #43) folds onto a four-element little
/// co-group with a single omega-regular class, whose projective table has one
/// two-dimensional irrep with character `(2, 0, 0, 0)`.  The probe is complete
/// and its identity content still matches the pinned frequency.
#[test]
fn a_two_dimensional_co_group_is_answered_by_the_projective_table() {
    let (subgroup, embedding) = ordinal_context(109, 3988);
    assert_eq!(embedding.subgroup_sg(), 43);
    let probe = probe(109, "P1");
    let result = subduce_full_star_with_embedding(&subgroup, &embedding, probe)
        .expect("the projective table answers this context");
    let trivial = query::irreps_of(43)
        .iter()
        .find(|record| {
            !record.spinor && record.k_vector().numerators == [0, 0, 0] && record.dim == 1
        })
        .expect("#43 has a trivial Gamma row");
    let content = trivial_content_with_embedding(&subgroup, &embedding, probe)
        .expect("identity-only content");
    assert_eq!(content.total, content.by_label);
    assert_eq!(content.total, stored_frequency(&subgroup, "P1"));
    let full_total: u32 = result
        .blocks()
        .iter()
        .map(|block| block.multiplicity(trivial.ml))
        .sum();
    assert_eq!(full_total, content.total);
}

/// A constructed star is not a per-point fallback: a star that reaches a pinned
/// row anywhere keeps answering from pinned rows only.  Ordinal 1007 (SG 44 ->
/// #8, `P1`) folds onto stars that contain the stored `V1`/`L1` rows on one arm
/// and no stored row on another; presenting both a stored row and a constructed
/// component for one physical irrep is rejected by the solver, so the star-level
/// choice is what keeps this context consistent.
#[test]
fn a_star_with_a_stored_row_never_also_constructs() {
    let (subgroup, embedding) = ordinal_context(44, 1007);
    assert_eq!(embedding.subgroup_sg(), 8);
    let mut stored_targets = 0usize;
    let mut constructed_targets = 0usize;
    for probe in ordinary_probes(44) {
        let result = subduce_full_star_with_embedding(&subgroup, &embedding, probe)
            .unwrap_or_else(|error| panic!("ordinal 1007, probe {}: {error}", probe.ml));
        for block in result.blocks() {
            let constructed = block
                .targets()
                .iter()
                .filter(|target| matches!(target.component, SubductionComponent::Constructed { .. }))
                .count();
            if constructed > 0 {
                assert_eq!(
                    constructed,
                    block.targets().len(),
                    "probe {}: one star is either stored or constructed, never mixed",
                    probe.ml
                );
                // One constructed component per star, not one per arm.
                assert_eq!(constructed, 1, "probe {}", probe.ml);
                constructed_targets += constructed;
            } else {
                stored_targets += block.targets().len();
            }
        }
    }
    assert!(stored_targets > 0, "stored rows must still answer their stars");
    assert!(
        constructed_targets > 0,
        "the all-analytic stars of this context are constructed"
    );
}
