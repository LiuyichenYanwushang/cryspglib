//! Oracle basis conventions must fix labels even when the operation set does not.
use cryspglib::irrep::LabelConvention;
use cryspglib::irrep::isotropy::{IsotropyDirection, isotropy_subgroup_for_direction};
use cryspglib::irrep::query;
use cryspglib::irrep::subduce_irrep;
use cryspglib::irrep::subduction::star::decompose::{
    FullStarError, subduce_full_star_with_embedding,
};
use cryspglib::irrep::subduction::{Mat3R, SeitzTransform, SubgroupEmbedding, Vec3R};

#[test]
fn swapping_child_axes_preserves_the_group_but_changes_irrep_labels() {
    let subgroup =
        isotropy_subgroup_for_direction(43, "GM1", LabelConvention::Cdml, IsotropyDirection::Label("P1")).unwrap();
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
    let swap = SeitzTransform::new(
        Mat3R::from_ints([[0, 1, 0], [1, 0, 0], [0, 0, 1]]),
        Vec3R::zero(),
    );
    let operations = embedding.representatives();
    for operation in operations {
        let exchanged = swap
            .map_operation(operation)
            .unwrap()
            .reduce(embedding.subgroup_lattice())
            .unwrap();
        assert!(operations.contains(&exchanged));
    }
    // Raw Gamma traces on the child's x-normal mirror distinguish the labels.
    // An operation-set-only tie breaker cannot detect this permutation.
    let mirror = [-1, 0, 0, 0, 1, 0, 0, 0, 1];
    for (ml, expected) in [("GM3", 1.0), ("GM4", -1.0)] {
        let record = query::irreps_of(43).iter().find(|r| r.ml == ml).unwrap();
        let row = record.ordinary_scalar_selected_arm_block_trace().unwrap();
        let mut checked = 0;
        for (operation, value) in row.operations().iter().zip(row.values()) {
            if operation.rotation == mirror {
                assert!((*value - expected).norm() < 1e-9);
                checked += 1;
            }
        }
        assert!(checked > 0);
    }
    let result = subduce_irrep(&subgroup, "GM3").unwrap();
    assert_eq!(result.multiplicity("GM3"), 1);
    assert_eq!(result.multiplicity("GM4"), 0);
}

#[test]
fn official_axes_pin_nontrivial_labels_from_two_cubic_condensates() {
    // SHOW BASIS fixes each child's x axis to a coordinate mirror of the
    // parent, where the raw GM2+ trace is +1.  The diagonal mirror has trace
    // -1, selecting GM3 (and excluding the swapped-axis GM4 convention).
    for condensing in ["GM4-", "GM5-"] {
        let subgroup =
            isotropy_subgroup_for_direction(230, condensing, LabelConvention::Cdml, IsotropyDirection::Label("P2"))
                .unwrap();
        let result = subduce_irrep(&subgroup, "GM2+").unwrap();
        assert_eq!(result.multiplicity("GM3"), 1, "{condensing}");
        assert_eq!(result.multiplicity("GM4"), 0, "{condensing}");
        assert_eq!(result.parent_dimension(), 1);
    }
}

#[test]
fn oracle_settings_include_genuine_shears() {
    for (sg, ml, child, setting) in [
        (21, "R1", 22, [[1, -1, 0], [1, 0, -1], [1, 0, 0]]),
        (64, "GM2+", 14, [[2, 0, 1], [-1, 0, -1], [0, 1, 0]]),
        (67, "GM2+", 13, [[2, 0, 1], [-1, 0, -1], [0, 1, 0]]),
    ] {
        let subgroup =
            isotropy_subgroup_for_direction(sg, ml, LabelConvention::Cdml, IsotropyDirection::Label("P1")).unwrap();
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
        assert_eq!(embedding.subgroup_sg(), child);
        assert_eq!(embedding.setting(), setting);
        assert_eq!(embedding.representatives().len(), 4);
    }
}

#[test]
fn every_new_record_probe_is_computed_or_matches_the_exact_missing_set() {
    let mut missing = Vec::new();
    let (mut complete, mut positive) = (0, 0);
    for (sg, ml, direction) in [
        (21, "R1", "P1"),
        (43, "GM1", "P1"),
        (64, "GM2+", "P1"),
        (67, "GM2+", "P1"),
        (230, "GM4-", "P2"),
        (230, "GM5-", "P2"),
    ] {
        let subgroup =
            isotropy_subgroup_for_direction(sg, ml, LabelConvention::Cdml, IsotropyDirection::Label(direction)).unwrap();
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
        let trivial = if matches!(embedding.subgroup_sg(), 13 | 14) {
            "GM1+"
        } else {
            "GM1"
        };
        let stored = subgroup.identity_subduction().unwrap();
        for probe in query::irreps_of(sg).iter().filter(|record| !record.spinor) {
            let result = match subduce_full_star_with_embedding(&subgroup, &embedding, probe) {
                Ok(result) => result,
                Err(FullStarError::MissingChildStarData { .. }) => {
                    missing.push((subgroup.ordinal, probe.ml));
                    continue;
                }
                Err(error) => panic!("ordinal {} probe {}: {error}", subgroup.ordinal, probe.ml),
            };
            let actual: u32 = result
                .blocks()
                .iter()
                .map(|block| block.multiplicity(trivial))
                .sum();
            let expected = stored
                .iter()
                .find(|entry| entry.parent_ml == probe.ml)
                .map_or(0, |entry| u32::from(entry.frequency));
            assert_eq!(actual, expected, "SG {sg} {ml}, probe {}", probe.ml);
            assert_eq!(result.covered_dimension(), result.parent_dimension());
            complete += 1;
            positive += usize::from(expected != 0);
        }
    }
    assert_eq!((complete, positive), (116, 18));
    assert_eq!(
        missing,
        [
            (15125, "P1P2"),
            (15125, "P3"),
            (15131, "P1P2"),
            (15131, "P3")
        ]
    );
}
