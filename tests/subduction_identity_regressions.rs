//! Stored identity frequencies check the embeddings independently of reconstruction.

use cryspglib::irrep::isotropy::{
    IsotropyDirection, isotropy_subgroup_for_direction, isotropy_subgroups,
};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::star::decompose::{
    FullStarError, subduce_full_star_with_embedding,
};
use cryspglib::irrep::subduction::{SubgroupEmbedding, subduce_irrep_with_embedding};
use cryspglib::irrep::types::{IrrepRecord, IrrepSourceIdentity};

fn gamma(record: &IrrepRecord) -> bool {
    !record.spinor && record.kx == 0 && record.ky == 0 && record.kz == 0
}

#[test]
fn non_gamma_pullback_keeps_the_phase_of_removed_child_translations() {
    let subgroup =
        isotropy_subgroup_for_direction(16, "R2", IsotropyDirection::Label("P1")).unwrap();
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
    let probe = query::irreps_of(16)
        .iter()
        .find(|record| record.ml == "X1")
        .unwrap();
    let result = subduce_irrep_with_embedding(&subgroup, &embedding, probe).unwrap();
    // SHOW ELEMENTS pins B=2I, origin=(1/2,1/2,0), and the C2z translation
    // (1,1,0). With the correct subgroup representatives but prematurely reduced
    // child-frame operations, this silently becomes T2 and still reconstructs
    // the incorrectly paired characters. The unreduced Hall mapping gives T3.
    assert_eq!(result.subgroup_sg(), 22);
    assert_eq!(result.targets().len(), 1);
    assert_eq!(result.multiplicity("T3"), 1);
    assert_eq!(result.multiplicity("T2"), 0);
    assert_eq!(result.parent_dimension(), 1);
}

#[test]
fn frozen_pairs_preserve_stored_gamma_frequencies_for_all_condensing_k() {
    let mut records = 0;
    let mut comparisons = 0;
    // Include non-Gamma condensates: they can remove parent translations and
    // expose errors hidden by the existing Gamma-condensate Frobenius gate.
    for sg in [16, 139, 167, 221, 225] {
        for condensing in query::irreps_of(sg).iter().filter(|record| !record.spinor) {
            for subgroup in isotropy_subgroups(sg, condensing.ml).unwrap() {
                if !matches!(
                    (sg, subgroup.record.sg),
                    (221, 83 | 12 | 148 | 123 | 47) | (225, 8) | (16, 22) | (167, 15) | (139, 126)
                ) {
                    continue;
                }
                let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup)
                    .unwrap_or_else(|error| panic!("ordinal {}: {error}", subgroup.ordinal));
                let trivial = query::irreps_of(embedding.subgroup_sg())
                    .iter()
                    .find(|record| {
                        gamma(record)
                            && record
                                .ordinary_scalar_selected_arm_block_trace()
                                .is_ok_and(|row| {
                                    row.values().iter().all(|value| (value - 1.0).norm() < 1e-9)
                                })
                    })
                    .expect("the child's trivial Gamma irrep");
                let stored = subgroup.identity_subduction().unwrap();
                for probe in query::irreps_of(sg).iter().filter(|record| gamma(record)) {
                    let frequencies: Vec<_> = stored
                        .iter()
                        .filter(|entry| {
                            entry.parent_ml == probe.ml && entry.parent_k.numerators == [0, 0, 0]
                        })
                        .map(|entry| entry.frequency)
                        .collect();
                    let expected = frequencies.first().copied().unwrap_or(0);
                    assert!(frequencies.iter().all(|frequency| *frequency == expected));
                    let result = subduce_irrep_with_embedding(&subgroup, &embedding, probe)
                        .unwrap_or_else(|error| {
                            panic!("ordinal {}, probe {}: {error}", subgroup.ordinal, probe.ml)
                        });
                    assert_eq!(
                        result.multiplicity(trivial.ml),
                        u32::from(expected),
                        "ordinal {}, probe {}: stored identity frequency",
                        subgroup.ordinal,
                        probe.ml
                    );
                    comparisons += 1;
                }
                records += 1;
            }
        }
    }
    assert_eq!((records, comparisons), (59, 538));
}

#[test]
fn full_stars_preserve_stored_frequencies_including_non_gamma_probes() {
    let mut records = 0;
    let mut comparisons = 0;
    let mut positive = 0;
    let mut non_gamma = 0;
    let mut missing = Vec::new();
    for sg in [16, 139, 167, 221, 225] {
        for condensing in query::irreps_of(sg).iter().filter(|record| !record.spinor) {
            for subgroup in isotropy_subgroups(sg, condensing.ml).unwrap() {
                if !matches!(
                    (sg, subgroup.record.sg),
                    (221, 83 | 12 | 148 | 123 | 47) | (225, 8) | (16, 22) | (167, 15) | (139, 126)
                ) {
                    continue;
                }
                let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
                let trivial = query::irreps_of(embedding.subgroup_sg())
                    .iter()
                    .find(|record| {
                        gamma(record)
                            && record
                                .ordinary_scalar_selected_arm_block_trace()
                                .is_ok_and(|row| {
                                    row.values().iter().all(|value| (value - 1.0).norm() < 1e-9)
                                })
                    })
                    .unwrap();
                let stored = subgroup.identity_subduction().unwrap();
                for probe in query::irreps_of(sg).iter().filter(|record| {
                    matches!(
                        record.source_identity(),
                        IrrepSourceIdentity::OrdinaryScalar { .. }
                    )
                }) {
                    let result =
                        match subduce_full_star_with_embedding(&subgroup, &embedding, probe) {
                            Ok(result) => result,
                            Err(FullStarError::MissingChildStarData { .. }) => {
                                missing.push((subgroup.ordinal, probe.ml));
                                continue;
                            }
                            Err(error) => {
                                panic!("ordinal {}, probe {}: {error}", subgroup.ordinal, probe.ml)
                            }
                        };
                    let frequencies: Vec<_> = stored
                        .iter()
                        .filter(|entry| entry.parent_ml == probe.ml)
                        .map(|entry| entry.frequency)
                        .collect();
                    let expected = frequencies.first().copied().unwrap_or(0);
                    assert!(frequencies.iter().all(|frequency| *frequency == expected));
                    let actual: u32 = result
                        .blocks()
                        .iter()
                        .map(|block| block.multiplicity(trivial.ml))
                        .sum();
                    assert_eq!(
                        actual,
                        u32::from(expected),
                        "ordinal {}, probe {}: full-star stored frequency",
                        subgroup.ordinal,
                        probe.ml
                    );
                    comparisons += 1;
                    positive += usize::from(expected != 0);
                    non_gamma += usize::from(!gamma(probe));
                }
                records += 1;
            }
        }
    }
    assert_eq!(
        (records, comparisons, positive, non_gamma),
        (59, 2060, 458, 1522)
    );
    // These W stars fold to wave vectors absent from the discrete #8 table.
    // Pin the complete missing set: newly unsupported probes must fail this
    // gate rather than silently shrinking the tested coverage.
    let expected_missing: Vec<_> = [13345, 13346, 13351]
        .into_iter()
        .flat_map(|ordinal| {
            ["W1", "W2", "W3", "W4", "W5"]
                .into_iter()
                .map(move |ml| (ordinal, ml))
        })
        .collect();
    assert_eq!(missing, expected_missing);
}

#[test]
fn compound_full_stars_preserve_stored_identity_frequencies() {
    let (mut records, mut comparisons, mut positive, mut non_gamma) = (0, 0, 0, 0);
    let mut missing = Vec::new();
    for sg in [16, 19, 23, 45, 83, 139, 167, 221, 225] {
        for condensing in query::irreps_of(sg).iter().filter(|record| !record.spinor) {
            for subgroup in isotropy_subgroups(sg, condensing.ml).unwrap() {
                if !matches!(
                    (sg, subgroup.record.sg),
                    (221, 83 | 12 | 148 | 123 | 47)
                        | (225, 8)
                        | (16, 22)
                        | (167, 15)
                        | (139, 126)
                        | (19, 19)
                        | (23, 23)
                        | (45, 45)
                        | (83, 83)
                ) {
                    continue;
                }
                let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
                let trivial = query::irreps_of(embedding.subgroup_sg())
                    .iter()
                    .find(|record| {
                        gamma(record)
                            && record
                                .ordinary_scalar_selected_arm_block_trace()
                                .is_ok_and(|row| {
                                    row.values().iter().all(|value| (value - 1.0).norm() < 1e-9)
                                })
                    })
                    .unwrap();
                let stored = subgroup.identity_subduction().unwrap();
                for probe in query::irreps_of(sg).iter().filter(|record| {
                    matches!(
                        record.source_identity(),
                        IrrepSourceIdentity::Compound { .. }
                    )
                }) {
                    let result =
                        match subduce_full_star_with_embedding(&subgroup, &embedding, probe) {
                            Ok(result) => result,
                            Err(FullStarError::MissingChildStarData { .. }) => {
                                missing.push((subgroup.ordinal, probe.ml));
                                continue;
                            }
                            Err(error) => panic!(
                                "ordinal {}, compound {}: {error}",
                                subgroup.ordinal, probe.ml
                            ),
                        };
                    let frequencies: Vec<_> = stored
                        .iter()
                        .filter(|entry| entry.parent_ml == probe.ml)
                        .map(|entry| entry.frequency)
                        .collect();
                    let expected = frequencies.first().copied().unwrap_or(0);
                    assert!(frequencies.iter().all(|frequency| *frequency == expected));
                    let actual: u32 = result
                        .blocks()
                        .iter()
                        .map(|block| block.multiplicity(trivial.ml))
                        .sum();
                    assert_eq!(
                        actual,
                        u32::from(expected),
                        "ordinal {}, compound {}",
                        subgroup.ordinal,
                        probe.ml
                    );
                    comparisons += 1;
                    positive += usize::from(expected != 0);
                    non_gamma += usize::from(!gamma(probe));
                }
                records += 1;
            }
        }
    }
    assert_eq!(
        (records, comparisons, positive, non_gamma),
        (69, 78, 0, 64),
        "missing: {missing:?}"
    );
    // These frozen embeddings have no positive compound identity frequency:
    // all 78 comparisons are zero terms, 64 at non-Gamma probe wave vectors.
    // Positive complex multiplicities are pinned by the raw-CIR decomposition
    // witnesses in subduction_compound_stars; this is not full-table coverage.
    assert!(missing.is_empty(), "missing: {missing:?}");
}
