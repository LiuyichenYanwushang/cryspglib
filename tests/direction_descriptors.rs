//! R7: the ISOTROPY program's **own** direction descriptors.
//!
//! The strings users copy out of the ISOTROPY suite come from `SHOW DIRECTION
//! VECTOR`; the bundled tables carry a different notation for `dim = 2`, for
//! `dim >= 4` and in their separators.  These tests pin the collected table
//! (`scripts/data/direction_descriptors_v1.json` ->
//! `irrep::direction_descriptors_data`), the official-string selector, and the
//! honest gap reporting for the records and labels the program does not print.

use cryspglib::irrep::direction_descriptors_data::{
    OFFICIAL_RECORDS, OFFICIAL_TABLE_EMPTY, official_of,
};
use cryspglib::irrep::isotropy::{
    self, IsotropyDirection, IsotropyError, OfficialDescriptor, magnetic_isotropy_subgroup_for_direction,
    magnetic_isotropy_subgroups, official_direction_descriptor,
};
use cryspglib::irrep::{LabelConvention, query};

/// The milestone's acceptance records, pinned from the collected oracle output.
const ACCEPTANCE: &[(u8, &str, &str, &str)] = &[
    (221, "GM4+", "P1", "(a,0,0)"),
    (221, "GM4+", "P2", "(a,a,0)"),
    (221, "GM4+", "P3", "(a,a,a)"),
    (5, "L1", "P1", "(a;a)"),
    (91, "A1", "P1", "(a,0)"),
    (194, "GM6+", "P1", "(a,0)"),
    (177, "L1", "P1", "(a;0;0)"),
    (177, "L1", "S1", "(a;b;c)"),
];

fn subgroup(sg: u8, ml: &str, direction: IsotropyDirection<'_>) -> isotropy::IsotropySubgroup {
    isotropy::isotropy_subgroup_for_direction(sg, ml, LabelConvention::Cdml, direction)
        .unwrap_or_else(|error| panic!("SG {sg} {ml}: {error}"))
}

#[test]
fn the_acceptance_records_carry_their_official_strings() {
    for (sg, ml, label, text) in ACCEPTANCE {
        assert_eq!(
            official_direction_descriptor(*sg, ml, label),
            OfficialDescriptor::Official(text),
            "sg {sg} {ml} {label}"
        );
        // The official string selects the row it came from...
        let selected = subgroup(
            *sg,
            ml,
            IsotropyDirection::OfficialDescriptor(text),
        );
        assert_eq!(selected.record.direction_label, *label);
        // ... and so does the label.
        assert_eq!(
            subgroup(*sg, ml, IsotropyDirection::Label(label))
                .record
                .direction_label,
            *label
        );
    }
    // The complex-separator case: components of one vector use `;`.
    let official = official_of(177, "L1", "S1").expect("SG 177 L1 S1 is collected");
    assert_eq!(official, "(a;b;c)");
    // Two independent vectors separate them with `,`.
    let two_vectors = official_of(100, "R1", "4D1").expect("SG 100 R1 4D1 is collected");
    assert_eq!(two_vectors, "(a,b;c,d)");
}

/// Every collected official string selects its own row through the official
/// selector -- the whole table, not a sample.
#[test]
fn every_official_string_selects_its_own_direction() {
    let mut checked = 0usize;
    let mut wrong: Vec<String> = Vec::new();
    for record in OFFICIAL_RECORDS {
        for (label, text) in record.directions {
            let selected = isotropy::isotropy_subgroup_for_direction(
                record.sg,
                record.ml,
                LabelConvention::Cdml,
                IsotropyDirection::OfficialDescriptor(text),
            );
            checked += 1;
            match selected {
                Ok(subgroup) if subgroup.record.direction_label == *label => {}
                Ok(subgroup) => wrong.push(format!(
                    "{} {} {text} -> {}",
                    record.sg, record.ml, subgroup.record.direction_label
                )),
                Err(error) => wrong.push(format!("{} {} {text}: {error}", record.sg, record.ml)),
            }
        }
    }
    assert_eq!(wrong, Vec::<String>::new(), "{} wrong selections", wrong.len());
    assert_eq!(checked, 15044, "the collected table has 15,044 directions");
    assert_eq!(OFFICIAL_RECORDS.len(), 4665);
    assert_eq!(OFFICIAL_TABLE_EMPTY.len(), 112);
}

/// The official selector is a **separate** selector, not another spelling of the
/// internal one: the program's strings rarely normalize to the stored
/// descriptor of the row they belong to.
#[test]
fn the_official_and_internal_selectors_are_different() {
    let normalize_internal = |text: &str| -> String {
        text.chars()
            .filter(|character| !character.is_whitespace())
            .map(|character| if character == ';' { ',' } else { character })
            .collect()
    };
    let mut official_rows = 0usize;
    let mut same_as_internal = 0usize;
    let mut legacy_refuses = 0usize;
    let mut legacy_other_row = 0usize;
    for record in OFFICIAL_RECORDS {
        for (label, text) in record.directions {
            official_rows += 1;
            let rows = isotropy::isotropy_subgroups(record.sg, record.ml, LabelConvention::Cdml)
                .unwrap_or_else(|error| panic!("{} {}: {error}", record.sg, record.ml));
            let own = rows
                .iter()
                .find(|row| row.record.direction_label == *label)
                .unwrap_or_else(|| panic!("{} {} {label} missing", record.sg, record.ml));
            if normalize_internal(text) == normalize_internal(own.record.direction) {
                same_as_internal += 1;
            }
            match isotropy::isotropy_subgroup_for_direction(
                record.sg,
                record.ml,
                LabelConvention::Cdml,
                IsotropyDirection::Descriptor(text),
            ) {
                Err(_) => legacy_refuses += 1,
                Ok(selected) if selected.record.direction_label == *label => {}
                Ok(_) => legacy_other_row += 1,
            }
        }
    }
    assert_eq!(official_rows, 15044);
    // Measured on the collected corpus (scripts/data/direction_descriptors_v1.json):
    assert_eq!(same_as_internal, 7214, "rows whose two notations coincide");
    assert_eq!(legacy_refuses, 7158, "official strings the internal selector refuses");
    assert_eq!(legacy_other_row, 672, "official strings the internal selector sends elsewhere");

    // The witness: SG 101 `X1` swaps `P1` and `P3`.  Their official strings are
    // each other's internal descriptor, so the legacy selector answers with the
    // wrong row -- the official selector gets both right.
    assert_eq!(official_of(101, "X1", "P1"), Some("(a;a)"));
    assert_eq!(official_of(101, "X1", "P3"), Some("(a;0)"));
    let legacy = isotropy::isotropy_subgroup_for_direction(
        101,
        "X1",
        LabelConvention::Cdml,
        IsotropyDirection::Descriptor("(a;a)"),
    )
    .expect("the internal notation still answers");
    assert_eq!(legacy.record.direction_label, "P3", "the legacy row is the other one");
    assert_eq!(
        subgroup(
            101,
            "X1",
            IsotropyDirection::OfficialDescriptor("(a;a)")
        )
        .record
        .direction_label,
        "P1"
    );
}

/// Records whose program table is empty are reported as such, and the label
/// selector keeps working -- the internal notation is not passed off as the
/// program's string.
#[test]
fn empty_program_tables_are_reported_and_labels_still_work() {
    assert_eq!(OFFICIAL_TABLE_EMPTY.len(), 112);
    let (sg, ml, label) = (121u8, "P1P1", "C1");
    assert_eq!(
        official_direction_descriptor(sg, ml, label),
        OfficialDescriptor::TableEmpty
    );
    let error = isotropy::isotropy_subgroup_for_direction(
        sg,
        ml,
        LabelConvention::Cdml,
        IsotropyDirection::OfficialDescriptor("(a,b)"),
    )
    .expect_err("no official string exists for this record");
    assert!(matches!(error, IsotropyError::NoOfficialDescriptor { .. }));
    // The label selector answers, with the internal notation.
    let selected = subgroup(sg, ml, IsotropyDirection::Label(label));
    assert_eq!(selected.record.direction_label, label);
    assert!(
        !selected.record.direction.is_empty(),
        "the internal notation exists and is deliberately not labelled official"
    );
}

/// The 195 labels that have no program row all belong to the 112 records whose
/// program table is **empty**; they are reported as `TableEmpty` (a collected
/// fact), never silently replaced by the internal string.  The `LabelNotInTable`
/// outcome is defensive: 0 of the 15,044 official rows sit in a non-empty
/// record that lacks the label, and `isotropy`'s unit tests cover it.
#[test]
fn empty_program_tables_cover_every_label_without_a_program_row() {
    let mut empty_labels = 0usize;
    let mut not_in_table = 0usize;
    let mut witness = None;
    assert_eq!(OFFICIAL_TABLE_EMPTY.len(), 112);
    for (sg, ml) in OFFICIAL_TABLE_EMPTY {
        let rows = isotropy::isotropy_subgroups(*sg, ml, LabelConvention::Cdml)
            .unwrap_or_else(|error| panic!("SG {sg} {ml}: {error}"));
        for row in rows {
            empty_labels += 1;
            witness.get_or_insert((*sg, *ml, row.record.direction_label));
            assert_eq!(
                official_direction_descriptor(*sg, ml, row.record.direction_label),
                OfficialDescriptor::TableEmpty
            );
        }
    }
    for record in OFFICIAL_RECORDS {
        let rows = isotropy::isotropy_subgroups(record.sg, record.ml, LabelConvention::Cdml)
            .unwrap_or_else(|error| panic!("{} {}: {error}", record.sg, record.ml));
        for row in rows {
            if official_of(record.sg, record.ml, row.record.direction_label).is_none() {
                not_in_table += 1;
            }
        }
    }
    assert_eq!(empty_labels, 195, "labels inside empty program tables");
    assert_eq!(not_in_table, 0, "labels missing from a non-empty program table");

    let (sg, ml, label) = witness.expect("some empty program table exists");
    let selected = subgroup(sg, ml, IsotropyDirection::Label(label));
    assert_eq!(selected.record.direction_label, label);
    assert!(
        !selected.record.direction.is_empty(),
        "the internal notation exists and is deliberately not labelled official"
    );
}

/// Magnetic records print no direction-vector table in the program, so the
/// magnetic selector refuses an official string while the label keeps working;
/// the ordinary query judges the *ordinary* record of the same `(sg, ml)`.
#[test]
fn magnetic_records_have_no_official_direction_strings() {
    let mut witness = None;
    'scan: for sg in 1..=230u8 {
        for irrep in query::irreps_of(sg) {
            if irrep.spinor {
                continue;
            }
            let Ok(rows) = magnetic_isotropy_subgroups(sg, irrep.ml, LabelConvention::Cdml) else {
                continue;
            };
            if let Some(first) = rows.first() {
                witness = Some((sg, irrep.ml, first.record.direction));
                break 'scan;
            }
        }
    }
    let (sg, ml, label) = witness.expect("some magnetic record exists");
    let error = magnetic_isotropy_subgroup_for_direction(
        sg,
        ml,
        LabelConvention::Cdml,
        IsotropyDirection::OfficialDescriptor("(a)"),
    )
    .expect_err("magnetic records have no collected official strings");
    assert!(matches!(error, IsotropyError::NoOfficialDescriptor { .. }));
    // The magnetic label selector answers.
    let selected = magnetic_isotropy_subgroup_for_direction(
        sg,
        ml,
        LabelConvention::Cdml,
        IsotropyDirection::Label(label),
    )
    .unwrap_or_else(|error| panic!("SG {sg} {ml} magnetic {label}: {error}"));
    assert_eq!(selected.record.direction, label);
}
