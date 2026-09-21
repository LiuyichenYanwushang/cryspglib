//! Stored identity frequencies check the embeddings independently of reconstruction.

use cryspglib::irrep::LabelConvention;
use cryspglib::irrep::isotropy::{
    IsotropyDirection, isotropy_subgroup_for_direction, isotropy_subgroups,
};
use cryspglib::irrep::query;
use cryspglib::irrep::isotropy::IsotropySubgroup;
use cryspglib::irrep::subduction::star::decompose::{
    FullStarError, subduce_full_star_with_embedding, trivial_content_with_embedding,
};
use cryspglib::irrep::subduction::star::scalar_star::ScalarStar;
use cryspglib::irrep::subduction::{SubgroupEmbedding, subduce_irrep_with_embedding};
use cryspglib::irrep::types::{IrrepRecord, IrrepSourceIdentity};

fn gamma(record: &IrrepRecord) -> bool {
    !record.spinor && record.kx == 0 && record.ky == 0 && record.kz == 0
}

#[test]
fn non_gamma_pullback_keeps_the_phase_of_removed_child_translations() {
    let subgroup =
        isotropy_subgroup_for_direction(16, "R2", LabelConvention::Cdml, IsotropyDirection::Label("P1")).unwrap();
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
            for subgroup in isotropy_subgroups(sg, condensing.ml, LabelConvention::Cdml).unwrap() {
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
            for subgroup in isotropy_subgroups(sg, condensing.ml, LabelConvention::Cdml).unwrap() {
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
            for subgroup in isotropy_subgroups(sg, condensing.ml, LabelConvention::Cdml).unwrap() {
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

/// One isotropy record of SG 196 by its global ordinal.
fn sg196_record(ordinal: usize) -> (IsotropySubgroup, SubgroupEmbedding) {
    let subgroup = isotropy_subgroups(196, "W1", LabelConvention::Cdml)
        .unwrap()
        .into_iter()
        .find(|subgroup| subgroup.ordinal == ordinal)
        .unwrap_or_else(|| panic!("SG 196 W1 has no isotropy record {ordinal}"));
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
    (subgroup, embedding)
}

/// The pinned frequency of one probe in a record's stored identity subduction.
fn stored_frequency(subgroup: &IsotropySubgroup, ml: &str) -> u32 {
    let entries = subgroup.identity_subduction().unwrap();
    let frequencies: Vec<u32> = entries
        .iter()
        .filter(|entry| entry.parent_ml == ml)
        .map(|entry| u32::from(entry.frequency))
        .collect();
    assert!(
        frequencies.iter().all(|value| *value == frequencies[0]),
        "probe {ml} appears with two frequencies"
    );
    frequencies.first().copied().unwrap_or(0)
}

/// A folded child star whose wave vector is not the child's Gamma point cannot
/// carry the subgroup's trivial representation, so the pinned identity
/// frequency stays exactly computable for probes whose *full* decomposition has
/// no child data at all.
#[test]
fn identity_only_content_answers_probes_without_full_child_data() {
    // SG 196 W1 -> #24 (I2_12_12_1), direction P2, Size 4.
    let (subgroup, embedding) = sg196_record(10027);
    assert_eq!((subgroup.record.sg, subgroup.record.direction_label), (24, "P2"));
    let probe = |ml: &str| {
        query::irreps_of(196)
            .iter()
            .find(|record| record.ml == ml)
            .unwrap()
    };

    // W1 folds onto three child stars; two of them have no stored child k.
    let w1 = probe("W1");
    assert!(matches!(
        subduce_full_star_with_embedding(&subgroup, &embedding, w1),
        Err(FullStarError::MissingChildStarData { sg: 24, .. })
    ));
    let content = trivial_content_with_embedding(&subgroup, &embedding, w1).unwrap();
    assert_eq!(
        (content.total, content.by_label, content.gamma_stars, content.skipped_stars),
        (1, 1, 1, 2)
    );
    assert_eq!(content.total, stored_frequency(&subgroup, "W1"));
    assert_eq!(stored_frequency(&subgroup, "W1"), 1);

    // The same geometry with a probe the pinned table does not list: the Gamma
    // block really is decomposed, and it carries no trivial term.
    let w2 = probe("W2");
    assert!(matches!(
        subduce_full_star_with_embedding(&subgroup, &embedding, w2),
        Err(FullStarError::MissingChildStarData { sg: 24, .. })
    ));
    let content = trivial_content_with_embedding(&subgroup, &embedding, w2).unwrap();
    assert_eq!((content.total, content.gamma_stars, content.skipped_stars), (0, 1, 2));
    assert_eq!(stored_frequency(&subgroup, "W2"), 0);

    // L1 folds onto no Gamma star at all, and its blocks do have child data, so
    // the two entry points must agree on the zero.
    let l1 = probe("L1");
    let full = subduce_full_star_with_embedding(&subgroup, &embedding, l1).unwrap();
    let trivial = trivial_record_of(embedding.subgroup_sg());
    let full_total: u32 = full
        .blocks()
        .iter()
        .map(|block| block.multiplicity(trivial.ml))
        .sum();
    let content = trivial_content_with_embedding(&subgroup, &embedding, l1).unwrap();
    assert_eq!((content.total, content.gamma_stars, full_total), (0, 0, 0));
    let folded = ScalarStar::new(l1).unwrap().folded_stars(&embedding).unwrap();
    assert_eq!(
        content.gamma_stars + content.skipped_stars,
        folded.len(),
        "every folded child star is either evaluated or skipped"
    );
    assert_eq!(stored_frequency(&subgroup, "L1"), 0);
}

/// The child's trivial Gamma irrep, found the way the audit finds it: a
/// one-dimensional Gamma row that is `+1` on every operation.
fn trivial_record_of(sg: u8) -> &'static IrrepRecord {
    query::irreps_of(sg)
        .iter()
        .find(|record| {
            gamma(record)
                && record.dim == 1
                && record
                    .ordinary_scalar_selected_arm_block_trace()
                    .is_ok_and(|row| row.values().iter().all(|value| (value - 1.0).norm() < 1e-9))
        })
        .expect("the child's trivial Gamma irrep")
}

/// The identity-only entry point must agree with the fully validated full-star
/// decomposition wherever that one has data, and it must answer every probe of
/// the pinned record set that the full path cannot.
#[test]
fn identity_only_content_agrees_with_the_full_decomposition_and_covers_the_pinned_set() {
    let (mut agreements, mut covered_missing, mut positives) = (0, 0, 0);
    for sg in [16, 139, 167, 221, 225] {
        for condensing in query::irreps_of(sg).iter().filter(|record| !record.spinor) {
            for subgroup in isotropy_subgroups(sg, condensing.ml, LabelConvention::Cdml).unwrap() {
                if !matches!(
                    (sg, subgroup.record.sg),
                    (221, 83 | 12 | 148 | 123 | 47) | (225, 8) | (16, 22) | (167, 15) | (139, 126)
                ) {
                    continue;
                }
                let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
                let trivial = trivial_record_of(embedding.subgroup_sg());
                for probe in query::irreps_of(sg).iter().filter(|record| {
                    matches!(
                        record.source_identity(),
                        IrrepSourceIdentity::OrdinaryScalar { .. }
                    )
                }) {
                    let content = trivial_content_with_embedding(&subgroup, &embedding, probe)
                        .unwrap_or_else(|error| {
                            panic!("ordinal {}, probe {}: {error}", subgroup.ordinal, probe.ml)
                        });
                    assert_eq!(content.total, content.by_label);
                    positives += usize::from(content.total != 0);
                    match subduce_full_star_with_embedding(&subgroup, &embedding, probe) {
                        Ok(result) => {
                            let full_total: u32 = result
                                .blocks()
                                .iter()
                                .map(|block| block.multiplicity(trivial.ml))
                                .sum();
                            assert_eq!(
                                content.total, full_total,
                                "ordinal {}, probe {}: identity-only content",
                                subgroup.ordinal, probe.ml
                            );
                            agreements += 1;
                        }
                        Err(FullStarError::MissingChildStarData { .. }) => {
                            // The pinned table is the oracle for exactly these.
                            assert_eq!(
                                content.total,
                                stored_frequency(&subgroup, probe.ml),
                                "ordinal {}, probe {}: covered pinned frequency",
                                subgroup.ordinal,
                                probe.ml
                            );
                            covered_missing += 1;
                        }
                        Err(error) => {
                            panic!("ordinal {}, probe {}: {error}", subgroup.ordinal, probe.ml)
                        }
                    }
                }
            }
        }
    }
    // Pinned from the round-6 census: 2075 probes in total, of which 15
    // (ordinals 13345/13346/13351, probes W1-W5) have no child data for their
    // full star and are now answered exactly by the identity-only entry point;
    // the other 2060 agree with the fully validated full-star decomposition, and
    // 458 of all of them carry a positive trivial content.
    assert_eq!((covered_missing, agreements, positives), (15, 2060, 458));
    assert_eq!(covered_missing + agreements, 2075);
}
