//! The `other_wave_vector_subduction` frequencies computed by the engine.
//!
//! The 73 sources behind these rows are little irreps of lines `k = t*v`
//! through Gamma, so they have no discrete wave vector to fold and no character
//! rows in the pinned irrep table.  `src/irrep/w_little_characters_data.rs`
//! carries their frozen little-group character tables, and
//! `line_trivial_content_with_embedding` turns one of them into the pinned
//! frequency.  This test pins the measured behaviour of the `P1`-child family
//! (all 300 pinned rows of the 46 such records) and the ten records whose
//! `DT3`/`DT4` frequencies differ, which is what pins the frozen ordering of
//! that pair.

use cryspglib::irrep::isotropy::{self, IsotropySubgroup};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::SubgroupEmbedding;
use cryspglib::irrep::subduction::star::decompose::{
    line_trivial_content_via_blocks, line_trivial_content_with_embedding,
};
use cryspglib::irrep::w_little_characters_data::{LittleCharacterTable, W_LITTLE_CHARACTERS};
use cryspglib::irrep::LabelConvention;

/// Every isotropy record of one parent space group, as the audit enumerates it.
fn subgroups_of(sg: u8) -> Vec<IsotropySubgroup> {
    let mut out = Vec::new();
    for record in query::irreps_of(sg) {
        if record.spinor || record.subgroups().is_empty() {
            continue;
        }
        let subgroups = isotropy::isotropy_subgroups(sg, record.ml, LabelConvention::Cdml)
            .unwrap_or_else(|error| panic!("SG {sg} {}: {error}", record.ml));
        out.extend(subgroups);
    }
    out
}

/// The engine's frequency for one w row of one record.
fn frequency(subgroup: &IsotropySubgroup, label: &str) -> Result<u32, String> {
    let table = W_LITTLE_CHARACTERS
        .iter()
        .find(|table| {
            usize::from(table.space_group) == usize::from(subgroup.parent_sg)
                && table.label == label
        })
        .ok_or_else(|| format!("no frozen table for SG {} {label}", subgroup.parent_sg))?;
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup)
        .map_err(|error| format!("ordinal {}: {error}", subgroup.ordinal))?;
    line_trivial_content_with_embedding(subgroup, &embedding, table)
        .map_err(|error| format!("ordinal {} {label}: {error}", subgroup.ordinal))
}

#[test]
fn p1_child_records_match_every_pinned_row() {
    let mut checked = 0usize;
    let mut mismatched: Vec<String> = Vec::new();
    for sg in [196u8, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228] {
        for subgroup in subgroups_of(sg) {
            // A `P1` child: the record's own space group number is 1.
            if subgroup.record.sg != 1 {
                continue;
            }
            let Ok(rows) = subgroup.other_wave_vector_subduction() else {
                continue;
            };
            for entry in rows {
                checked += 1;
                match frequency(&subgroup, entry.parent_ml) {
                    Ok(value) if value == u32::from(entry.frequency) => {}
                    other => mismatched.push(format!(
                        "ordinal {} SG {sg} {}: computed {other:?}, pinned {}",
                        subgroup.ordinal, entry.parent_ml, entry.frequency
                    )),
                }
            }
        }
    }
    assert_eq!(
        checked, 300,
        "the P1-child family is 300 pinned rows; got {checked}"
    );
    assert!(
        mismatched.is_empty(),
        "{} mismatches, first few: {:?}",
        mismatched.len(),
        &mismatched[..mismatched.len().min(8)]
    );
}

#[test]
fn asymmetric_dt3_dt4_records_are_pinned() {
    // These ten records are the ones whose `DT3` and `DT4` frequencies differ,
    // so they are the falsifiable prediction of the frozen pair ordering.
    let expected: &[(usize, u8, u16, u16)] = &[
        (10422, 202, 2, 1),
        (10423, 202, 1, 2),
        (10446, 202, 1, 2),
        (10447, 202, 2, 1),
        (10485, 202, 1, 5),
        (10509, 202, 5, 1),
        (11171, 209, 2, 1),
        (11173, 209, 1, 2),
        (11202, 209, 2, 1),
        (11204, 209, 1, 2),
    ];
    for (ordinal, sg, dt3, dt4) in expected {
        let subgroup = subgroups_of(*sg)
            .into_iter()
            .find(|subgroup| subgroup.ordinal == *ordinal)
            .unwrap_or_else(|| panic!("ordinal {ordinal} of SG {sg} is missing"));
        assert_eq!(
            frequency(&subgroup, "DT3"),
            Ok(u32::from(*dt3)),
            "ordinal {ordinal} SG {sg} DT3"
        );
        assert_eq!(
            frequency(&subgroup, "DT4"),
            Ok(u32::from(*dt4)),
            "ordinal {ordinal} SG {sg} DT4"
        );
    }
}

#[test]
fn a_source_of_another_parent_is_rejected() {
    let subgroup = subgroups_of(196)
        .into_iter()
        .find(|subgroup| !subgroup.other_wave_vector_subduction().unwrap().is_empty())
        .expect("SG 196 has w records");
    let foreign = W_LITTLE_CHARACTERS
        .iter()
        .find(|table| table.space_group != subgroup.parent_sg)
        .expect("another parent's table");
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
    assert!(line_trivial_content_with_embedding(&subgroup, &embedding, foreign).is_err());
}

fn find(sg: u8, ordinal: usize) -> IsotropySubgroup {
    subgroups_of(sg)
        .into_iter()
        .find(|subgroup| subgroup.ordinal == ordinal)
        .unwrap_or_else(|| panic!("ordinal {ordinal} of SG {sg} is missing"))
}

fn table_for(subgroup: &IsotropySubgroup, label: &str) -> &'static LittleCharacterTable {
    W_LITTLE_CHARACTERS
        .iter()
        .find(|table| {
            usize::from(table.space_group) == usize::from(subgroup.parent_sg)
                && table.label == label
        })
        .unwrap_or_else(|| panic!("no frozen table for SG {} {label}", subgroup.parent_sg))
}

#[test]
fn an_embedding_of_another_record_is_rejected() {
    // The same source table and the same parent, but the embedding of a
    // different isotropy record: the answer would be computed from another
    // record's geometry, so it must be refused rather than returned.
    let subgroup = find(196, 10032);
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&find(196, 10030)).unwrap();
    let table = table_for(&subgroup, "DT1");
    assert!(line_trivial_content_via_blocks(&subgroup, &embedding, table).is_err());
}

#[test]
fn a_singular_record_basis_is_rejected() {
    let mut broken = find(196, 10030);
    broken.record.basis = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&find(196, 10030)).unwrap();
    let table = table_for(&broken, "DT1");
    assert!(line_trivial_content_via_blocks(&broken, &embedding, table).is_err());
}

#[test]
fn the_block_route_is_the_public_route() {
    // The retired hand-written sum disagreed with the pinned rows; the public
    // entry point delegates, so both must answer identically everywhere.
    for ordinal in [10030usize, 10032, 10033] {
        let subgroup = find(196, ordinal);
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
        for entry in subgroup.other_wave_vector_subduction().unwrap() {
            let table = table_for(&subgroup, entry.parent_ml);
            assert_eq!(
                line_trivial_content_with_embedding(&subgroup, &embedding, table),
                line_trivial_content_via_blocks(&subgroup, &embedding, table),
                "ordinal {ordinal} {}",
                entry.parent_ml
            );
        }
    }
}

#[test]
fn the_block_route_matches_every_pinned_row_of_sg196() {
    // The retired hand-written sum returned a wrong value for 30 of SG 196's w
    // rows, so the audit used to reach the pinned value through a second route
    // and never noticed.  The block route has to get all of them right on its
    // own; this is the assertion that makes an expected-answer-selected audit
    // unnecessary.
    let mut checked = 0usize;
    for subgroup in subgroups_of(196) {
        let Ok(rows) = subgroup.other_wave_vector_subduction() else {
            continue;
        };
        if rows.is_empty() {
            continue;
        }
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
        for entry in rows {
            let table = table_for(&subgroup, entry.parent_ml);
            assert_eq!(
                line_trivial_content_via_blocks(&subgroup, &embedding, table),
                Ok(u32::from(entry.frequency)),
                "ordinal {} {}",
                subgroup.ordinal,
                entry.parent_ml
            );
            checked += 1;
        }
    }
    assert_eq!(checked, 106, "SG 196 has 106 w rows");
}
