//! Stored identity frequencies check the embeddings independently of reconstruction.

use cryspglib::irrep::isotropy::{
    IsotropyDirection, isotropy_subgroup_for_direction, isotropy_subgroups,
};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::{SubgroupEmbedding, subduce_irrep_with_embedding};
use cryspglib::irrep::types::IrrepRecord;

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
