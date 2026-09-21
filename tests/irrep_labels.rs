//! Explicit CDML vs Bradley–Cracknell label queries.
//!
//! These tests pin the scientific facts the API must not paper over:
//!
//! * CDML labels are stored literally and available for every record, and a
//!   CDML query only matches that literal spelling (`"GM"`, never `"Γ"`).
//! * Stored BC numbering is not ML numbering (SG 213 `X1` ↔ `X2`,
//!   SG 123 `GM2+` ↔ `Γ3+`), so BC must be read from `IrrepRecord::bc`.
//! * BC misses are honest: the 12 `"***"` placeholder rows and all 3611
//!   spinor rows (whose `.bc` is synthesized from `.ml`) have no BC label.
//! * Duplicate BC spellings return every match (SG 1 `Γ1` = GM1 and Z1,
//!   SG 2 `Γ1±` = GM1± and Z1±), and reverse searches never dedup them.
//! * The generator LaTeX-escapes only the first component of a compound BC
//!   label (`\Gamma_{2}GM3`); BC normalizes every component to `Γ2Γ3`.
//! * k-point summaries carry both spellings and partition their rows exactly,
//!   under the requested convention.

use std::collections::BTreeSet;

use cryspglib::irrep::LabelConvention;
use cryspglib::irrep::generated_data::IRREPS;
use cryspglib::irrep::query::{
    find_irreps, irreps_at_k_label, irreps_of, irreps_with_magnetic_subgroup, irreps_with_subgroup,
    kpoints_of,
};
use cryspglib::irrep::types::IrrepRecord;

fn ids(records: &[&IrrepRecord]) -> Vec<u16> {
    records.iter().map(|record| record.id().index()).collect()
}

fn mls(records: &[&IrrepRecord]) -> Vec<&'static str> {
    records.iter().map(|record| record.ml).collect()
}

fn record(sg: u8, ml: &str) -> &'static IrrepRecord {
    irreps_of(sg)
        .iter()
        .find(|record| record.ml == ml)
        .unwrap_or_else(|| panic!("SG {sg} has no record {ml}"))
}

#[test]
fn aggregate_and_magnetic_outputs_keep_both_labels_and_honest_gaps() {
    use cryspglib::irrep::magnetic_summary::{
        SourceIrrepSummary, format_magnetic_kpoint_summary, magnetic_irrep_summary_by_uni,
    };
    use cryspglib::irrep::query::{isotropy_subgroups_of, magnetic_isotropy_subgroups_of};

    let ordinary = isotropy_subgroups_of(1);
    let magnetic = magnetic_isotropy_subgroups_of(1);
    for (ml, bc) in [("GM1", Some("Γ1")), ("U1", None)] {
        let ordinary_labels = ordinary
            .iter()
            .find(|row| row.ml_label == ml)
            .unwrap()
            .labels();
        let magnetic_labels = magnetic
            .iter()
            .find(|row| row.ml_label == ml)
            .unwrap()
            .labels();
        assert_eq!(ordinary_labels.cdml, ml);
        assert_eq!(ordinary_labels.bc.as_deref(), bc);
        assert_eq!(magnetic_labels, ordinary_labels);
    }

    let summary = magnetic_irrep_summary_by_uni(2).unwrap();
    assert_eq!(summary.kpoints.len(), 8);
    let gamma = summary
        .kpoints
        .iter()
        .find(|kp| kp.coords == (0, 0, 0, 1))
        .unwrap();
    assert_eq!(gamma.label, "GM");
    assert_eq!(gamma.bc_label.as_deref(), Some("Γ"));
    let sources: Vec<_> = gamma
        .coreps
        .iter()
        .flat_map(|corep| &corep.source_irreps)
        .collect();
    assert_eq!(
        sources
            .iter()
            .find(|source| source.ml == "GM1")
            .unwrap()
            .labels()
            .bc
            .as_deref(),
        Some("Γ1")
    );
    assert_eq!(
        sources
            .iter()
            .find(|source| source.ml == "GM2")
            .unwrap()
            .labels()
            .bc,
        None
    );
    let table = format_magnetic_kpoint_summary(gamma);
    assert!(table.contains("k-point CDML=GM / BC=Γ"), "{table}");
    assert!(table.contains("| source CDML | source BC |"), "{table}");
    assert!(table.contains("| GM1 | Γ1 |"), "{table}");
    assert!(table.contains("| GM2 | unavailable |"), "{table}");

    // A compound component must not inherit the whole compound's BC label.
    let compound = record(83, "GM3+GM4+");
    let component = SourceIrrepSummary {
        sg: 83,
        ml: "GM3+",
        bc: compound.bc,
        dim: 1,
        spinor: false,
    };
    assert_eq!(component.labels().bc, None);
}

/// Count of raw `.bc` fields that still carry the generator's CDML `GM` token
/// in a later compound component.
const RAW_BC_WITH_GM: usize = 69;

#[test]
fn cdml_labels_are_literal_and_available_for_every_record() {
    assert_eq!(IRREPS.len(), 8388);
    for ir in IRREPS.iter() {
        let label = ir
            .label(LabelConvention::Cdml)
            .expect("CDML is available for every record");
        assert_eq!(label, ir.ml);
        assert_eq!(ir.labels().cdml, ir.ml);
        assert_eq!(ir.k_labels().cdml, ir.k_label());
        let found = find_irreps(ir.sg, &label, LabelConvention::Cdml);
        assert!(
            found.iter().any(|hit| hit.id() == ir.id()),
            "SG {} {} was not found by its own CDML label",
            ir.sg,
            ir.ml
        );
    }

    // CDML is matched literally: BC spellings, markup and whitespace are not
    // silently normalized into it.
    assert!(find_irreps(123, "Γ3+", LabelConvention::Cdml).is_empty());
    assert!(find_irreps(123, "\\Gamma_{3}^+", LabelConvention::Cdml).is_empty());
    assert!(find_irreps(1, " GM1", LabelConvention::Cdml).is_empty());
    assert!(find_irreps(1, "gm1", LabelConvention::Cdml).is_empty());
    assert!(irreps_at_k_label(213, "Γ", LabelConvention::Cdml).is_empty());
    assert!(irreps_at_k_label(213, "X ", LabelConvention::Cdml).is_empty());
    assert_eq!(irreps_at_k_label(213, "GM", LabelConvention::Cdml).len(), 8);
}

#[test]
fn bc_availability_is_exactly_scalar_non_placeholder_rows() {
    let mut scalar_available = 0usize;
    let mut missing = 0usize;
    let mut spinor = 0usize;
    let mut raw_with_gm = 0usize;

    for ir in IRREPS.iter() {
        if ir.bc.contains("GM") {
            raw_with_gm += 1;
        }
        if ir.spinor {
            spinor += 1;
            assert!(
                ir.label(LabelConvention::Bc).is_none(),
                "SG {} {} exposed a synthesized BC label",
                ir.sg,
                ir.ml
            );
            assert!(ir.labels().bc.is_none());
            assert!(ir.k_labels().bc.is_none());
            continue;
        }
        if ir.bc == "***" {
            missing += 1;
            assert!(
                ir.label(LabelConvention::Bc).is_none(),
                "SG {} {} invented a BC label",
                ir.sg,
                ir.ml
            );
            assert!(ir.labels().bc.is_none());
            assert!(ir.k_labels().bc.is_none());
            continue;
        }

        scalar_available += 1;
        let label = ir.label(LabelConvention::Bc).unwrap();
        assert_eq!(ir.labels().bc.as_deref(), Some(label.as_str()));
        // Every available BC label is fully spelled: the generator's repeated
        // CDML `GM` token never survives into a BC label.
        assert!(
            !label.contains("GM"),
            "SG {} {} left a raw GM token in its BC label {label}",
            ir.sg,
            ir.ml
        );
        // The stored field itself is never rewritten; the raw LaTeX round-trip
        // below is the witness that it still carries the generator spelling.

        // Stored raw LaTeX, normalized Unicode and the label itself are one
        // query, and every roundtrip finds this exact record.
        let by_label = find_irreps(ir.sg, &label, LabelConvention::Bc);
        let by_raw = find_irreps(ir.sg, ir.bc, LabelConvention::Bc);
        assert_eq!(ids(&by_raw), ids(&by_label));
        assert!(
            by_label.iter().any(|hit| hit.id() == ir.id()),
            "SG {} {} was not found by its own BC label {label}",
            ir.sg,
            ir.ml
        );
    }

    assert_eq!(scalar_available, 4765);
    assert_eq!(missing, 12);
    assert_eq!(spinor, 3611);
    assert_eq!(scalar_available + missing + spinor, IRREPS.len());
    assert_eq!(raw_with_gm, RAW_BC_WITH_GM);

    // Pinned witness: the raw field is untouched while the BC label is honest.
    let gamma1 = record(1, "GM1");
    assert_eq!(gamma1.bc, "\\Gamma_{1}");
    assert_eq!(gamma1.label(LabelConvention::Bc).as_deref(), Some("Γ1"));
    let spinor = record(1, "GM2");
    assert!(spinor.spinor);
    assert_eq!(spinor.bc, "\\Gamma_{2}");
    assert!(spinor.label(LabelConvention::Bc).is_none());
    assert!(find_irreps(1, "Γ2", LabelConvention::Bc).is_empty());
}

#[test]
fn bc_duplicate_labels_return_every_match() {
    let gamma1 = find_irreps(1, "\\Gamma_{1}", LabelConvention::Bc);
    assert_eq!(mls(&gamma1), vec!["GM1", "Z1"]);
    assert_eq!(ids(&gamma1).len(), 2);
    assert_ne!(
        (gamma1[0].kx, gamma1[0].ky, gamma1[0].kz, gamma1[0].kd),
        (gamma1[1].kx, gamma1[1].ky, gamma1[1].kz, gamma1[1].kd),
        "the two SG 1 Γ1 records are different k-points"
    );

    assert_eq!(
        mls(&find_irreps(2, "Γ1+", LabelConvention::Bc)),
        vec!["GM1+", "Z1+"]
    );
    assert_eq!(
        mls(&find_irreps(2, "Γ1-", LabelConvention::Bc)),
        vec!["GM1-", "Z1-"]
    );
    assert!(find_irreps(2, "Γ1", LabelConvention::Bc).is_empty());

    // The corresponding CDML labels are unique, so only BC needs all matches.
    assert_eq!(find_irreps(1, "GM1", LabelConvention::Cdml).len(), 1);
    assert_eq!(find_irreps(1, "Z1", LabelConvention::Cdml).len(), 1);
}

#[test]
fn stored_bc_numbering_is_not_ml_numbering() {
    // SG 213 X: CDML X1 is stored BC X_{2}; CDML X2 is stored BC X_{1}.
    let cdml_x1 = find_irreps(213, "X1", LabelConvention::Cdml);
    assert_eq!(mls(&cdml_x1), vec!["X1"]);
    assert_eq!(cdml_x1[0].bc, "X_{2}");
    let bc_x1 = find_irreps(213, "X1", LabelConvention::Bc);
    assert_eq!(mls(&bc_x1), vec!["X2"]);
    assert_eq!(bc_x1[0].bc, "X_{1}");
    assert_ne!(cdml_x1[0].id(), bc_x1[0].id());

    // SG 123 Γ: CDML GM2+ is stored BC Γ3+; CDML GM3+ is stored BC Γ2+.
    assert_eq!(
        mls(&find_irreps(123, "GM2+", LabelConvention::Cdml)),
        vec!["GM2+"]
    );
    assert_eq!(
        mls(&find_irreps(123, "Γ3+", LabelConvention::Bc)),
        vec!["GM2+"]
    );
    assert_eq!(
        mls(&find_irreps(123, "Γ2+", LabelConvention::Bc)),
        vec!["GM3+"]
    );
    assert!(find_irreps(123, "GM2+", LabelConvention::Bc).is_empty());
    assert!(find_irreps(123, "Γ3+", LabelConvention::Cdml).is_empty());
}

#[test]
fn bc_search_accepts_unicode_and_latex_without_collapsing_labels() {
    let latex = find_irreps(123, "\\Gamma_{3}^+", LabelConvention::Bc);
    assert_eq!(mls(&latex), vec!["GM2+"]);
    for spelling in [
        "Γ3+",
        "Γ₃⁺",
        "$\\Gamma_{3}^{+}$",
        " \\Gamma _ { 3 } ^{ + } ",
    ] {
        assert_eq!(
            ids(&find_irreps(123, spelling, LabelConvention::Bc)),
            ids(&latex),
            "spelling {spelling:?} did not normalize to the stored BC label"
        );
    }

    // Signs stay distinct.
    assert_eq!(
        mls(&find_irreps(123, "Γ3-", LabelConvention::Bc)),
        vec!["GM2-"]
    );
    assert_ne!(
        ids(&find_irreps(123, "Γ3+", LabelConvention::Bc)),
        ids(&find_irreps(123, "Γ3-", LabelConvention::Bc))
    );

    // Compound components stay compound.
    assert_eq!(
        mls(&find_irreps(213, "M1M2", LabelConvention::Bc)),
        vec!["M1M4"]
    );
    assert!(find_irreps(213, "M1", LabelConvention::Bc).is_empty());
    assert_eq!(
        mls(&find_irreps(182, "H1H2", LabelConvention::Bc)),
        vec!["H1H2"]
    );

    // Unparsed text is preserved, never dropped into a match.
    assert!(find_irreps(123, "\\Unknown_{3}^+", LabelConvention::Bc).is_empty());
    assert!(find_irreps(123, "\\UnknownGM_{1}", LabelConvention::Bc).is_empty());
    assert!(find_irreps(1, "Γ1x", LabelConvention::Bc).is_empty());

    // A leading GM is the CDML spelling of Γ and never crosses into BC.
    assert!(find_irreps(1, "GM1", LabelConvention::Bc).is_empty());
    assert!(irreps_at_k_label(213, "GM", LabelConvention::Bc).is_empty());
    assert_eq!(
        ids(&irreps_at_k_label(213, "\\Gamma", LabelConvention::Bc)),
        ids(&irreps_at_k_label(213, "Γ", LabelConvention::Bc))
    );
    assert!(irreps_at_k_label(213, "Γ", LabelConvention::Cdml).is_empty());
}

#[test]
fn gamma_compound_bc_labels_normalize_every_component() {
    // SG 143 GM2GM3: the generator LaTeX-escaped only the first component.
    let compound = record(143, "GM2GM3");
    assert_eq!(compound.bc, "\\Gamma_{2}GM3");
    assert_eq!(
        compound.label(LabelConvention::Bc).as_deref(),
        Some("Γ2Γ3"),
        "the stored raw LaTeX must normalize to a fully Γ-spelled BC label"
    );
    assert_eq!(compound.k_labels().bc.as_deref(), Some("Γ"));
    assert_eq!(compound.k_labels().cdml, "GM");

    let by_raw = find_irreps(143, "\\Gamma_{2}GM3", LabelConvention::Bc);
    let by_unicode = find_irreps(143, "Γ2Γ3", LabelConvention::Bc);
    assert_eq!(mls(&by_raw), vec!["GM2GM3"]);
    assert_eq!(ids(&by_raw), ids(&by_unicode));
    assert_eq!(
        ids(&by_raw),
        ids(&find_irreps(143, "Γ₂Γ₃", LabelConvention::Bc))
    );
    // The CDML compound spelling is not a BC spelling.
    assert!(find_irreps(143, "GM2GM3", LabelConvention::Bc).is_empty());

    // Signed compound witnesses keep their signs in every component.
    let plus = record(147, "GM2+GM3+");
    let minus = record(147, "GM2-GM3-");
    assert_eq!(plus.bc, "\\Gamma_{2}^+GM3+");
    assert_eq!(minus.bc, "\\Gamma_{2}^-GM3-");
    assert_eq!(plus.label(LabelConvention::Bc).as_deref(), Some("Γ2+Γ3+"));
    assert_eq!(minus.label(LabelConvention::Bc).as_deref(), Some("Γ2-Γ3-"));
    assert_eq!(
        mls(&find_irreps(147, "Γ2+Γ3+", LabelConvention::Bc)),
        vec!["GM2+GM3+"]
    );
    assert_eq!(
        mls(&find_irreps(147, "Γ2-Γ3-", LabelConvention::Bc)),
        vec!["GM2-GM3-"]
    );
    assert_ne!(plus.id(), minus.id());
}

#[test]
fn kpoints_preserve_physical_coordinates_and_partition_every_record() {
    let mut cdml_total = 0usize;
    let mut bc_total = 0usize;
    let mut point_total = 0usize;

    for sg in 1u8..=230 {
        let irreps = irreps_of(sg);

        // One physical k-point must include scalar AND spinor rows. Missing
        // spinor BC labels cannot split the magnetic-summary input in two.
        let mut covered = vec![false; irreps.len()];
        let mut summaries = 0usize;
        let mut coordinates = BTreeSet::new();
        for kpoint in kpoints_of(sg, LabelConvention::Cdml) {
            summaries += 1;
            point_total += 1;
            assert!(
                coordinates.insert(kpoint.coords),
                "SG {sg}: duplicate physical k-point"
            );
            assert!(!kpoint.irreps.is_empty());
            for &index in &kpoint.irreps {
                assert!(!covered[index], "SG {sg} row {index} repeated in CDML");
                covered[index] = true;
                let ir = &irreps[index];
                assert_eq!(kpoint.coords, (ir.kx, ir.ky, ir.kz, ir.kd));
                assert_eq!(kpoint.labels.cdml, ir.k_label());
                if let Some(bc) = ir.k_labels().bc {
                    assert_eq!(kpoint.labels.bc.as_deref(), Some(bc.as_str()));
                }
            }
            cdml_total += kpoint.irreps.len();
        }
        assert!(summaries > 0, "SG {sg} has no CDML k-point");
        assert!(
            covered.iter().all(|&seen| seen),
            "SG {sg} CDML k-points lost a record"
        );

        // BC: exactly the rows with a genuine BC label, each once.
        let mut bc_covered = vec![false; irreps.len()];
        for kpoint in kpoints_of(sg, LabelConvention::Bc) {
            assert!(kpoint.labels.bc.is_some());
            for &index in &kpoint.irreps {
                assert!(!bc_covered[index], "SG {sg} row {index} repeated in BC");
                bc_covered[index] = true;
                let ir = &irreps[index];
                assert!(!ir.spinor && ir.bc != "***");
                assert_eq!(
                    kpoint.labels.bc.as_deref(),
                    ir.k_label_with_convention(LabelConvention::Bc).as_deref()
                );
                assert_eq!(kpoint.labels.cdml, ir.k_label());
                assert_eq!(kpoint.coords, (ir.kx, ir.ky, ir.kz, ir.kd));
            }
            bc_total += kpoint.irreps.len();
        }
        for (index, ir) in irreps.iter().enumerate() {
            let available = !ir.spinor && ir.bc != "***";
            assert_eq!(
                bc_covered[index], available,
                "SG {sg} {} BC k-point coverage",
                ir.ml
            );
        }
    }

    assert_eq!(cdml_total, IRREPS.len());
    assert_eq!(bc_total, 4765);
    assert_eq!(point_total, 1350); // Independent census of frozen (SG,k) tuples.
}

#[test]
fn tables_show_cdml_and_unicode_bc_together() {
    use cryspglib::irrep::query::{
        format_character_table, format_isotropy_table, format_magnetic_isotropy_table,
    };
    for table in [
        format_character_table(123, 0, 0, 0, 1),
        format_isotropy_table(123, 0, 0, 0, 1),
        format_magnetic_isotropy_table(123, 0, 0, 0, 1),
    ] {
        assert!(table.contains("| CDML | BC |"), "{table}");
        assert!(table.contains("| GM2+ | Γ3+ |"), "{table}");
        assert!(!table.contains(r"\Gamma"), "{table}");
    }
    let spinor_table = format_character_table(1, 0, 0, 0, 1);
    assert!(
        spinor_table.contains("| GM2 | unavailable |"),
        "{spinor_table}"
    );
}

#[test]
fn bc_kpoints_use_stored_prefixes_and_skip_unavailable_rows() {
    let cdml = kpoints_of(1, LabelConvention::Cdml);
    let cdml_labels: BTreeSet<_> = cdml.iter().map(|kpoint| kpoint.labels.cdml).collect();
    assert!(cdml_labels.contains("GM"));
    assert!(cdml_labels.contains("U"));
    assert_eq!(
        cdml.iter().map(|kpoint| kpoint.irreps.len()).sum::<usize>(),
        16
    );

    let bc = kpoints_of(1, LabelConvention::Bc);
    assert_eq!(
        bc.iter()
            .map(|kpoint| kpoint.labels.bc.clone().unwrap())
            .collect::<BTreeSet<_>>(),
        BTreeSet::from(["Γ".to_string(), "B".to_string(), "F".to_string()])
    );
    // BC names both GM and Z Γ: same symbol, different coordinates.
    let gammas: Vec<_> = bc
        .iter()
        .filter(|kpoint| kpoint.labels.bc.as_deref() == Some("Γ"))
        .collect();
    assert_eq!(gammas.len(), 2);
    assert_ne!(gammas[0].coords, gammas[1].coords);
    assert_eq!(
        gammas
            .iter()
            .map(|kpoint| kpoint.labels.cdml)
            .collect::<Vec<_>>(),
        vec!["GM", "Z"]
    );
    assert_eq!(
        bc.iter().map(|kpoint| kpoint.irreps.len()).sum::<usize>(),
        4
    );
    for kpoint in &bc {
        for &index in &kpoint.irreps {
            let ir = &irreps_of(1)[index];
            assert!(!ir.spinor && ir.bc != "***");
        }
    }

    // A known scalar BC k-prefix labels the physical point, even when the
    // spinor irrep at that point has no verified BC irrep label.
    let gamma_cdml = cdml
        .iter()
        .find(|kpoint| kpoint.labels.cdml == "GM" && kpoint.labels.bc.is_some())
        .expect("SG 1 Γ summary");
    assert_eq!(gamma_cdml.labels.bc.as_deref(), Some("Γ"));
    assert_eq!(gamma_cdml.irreps.len(), 2);
    assert_eq!(irreps_of(1)[gamma_cdml.irreps[0]].ml, "GM1");
    let unavailable = cdml
        .iter()
        .find(|kpoint| kpoint.labels.cdml == "T")
        .expect("SG 1 T summary");
    assert!(unavailable.labels.bc.is_none());
    assert_eq!(unavailable.irreps.len(), 2, "T1 placeholder + T spinors");

    // SG 213: the BC list keeps only the 12 scalar rows; spinors stay out.
    let bc_213 = kpoints_of(213, LabelConvention::Bc);
    assert_eq!(
        bc_213
            .iter()
            .map(|kpoint| kpoint.labels.bc.clone().unwrap())
            .collect::<Vec<_>>(),
        vec!["Γ", "X", "M", "R"]
    );
    assert_eq!(
        bc_213
            .iter()
            .map(|kpoint| kpoint.irreps.len())
            .sum::<usize>(),
        12
    );
    assert_eq!(
        kpoints_of(213, LabelConvention::Cdml)
            .iter()
            .map(|kpoint| kpoint.irreps.len())
            .sum::<usize>(),
        27
    );
}

#[test]
fn reverse_searches_return_all_source_identities() {
    let cdml = irreps_with_subgroup(1, LabelConvention::Cdml);
    assert_eq!(cdml.len(), 843);
    let mut seen = BTreeSet::new();
    for ir in &cdml {
        assert!(seen.insert(ir.id()), "duplicate source identity returned");
        assert!(!ir.spinor);
        assert!(ir.subgroups().iter().any(|sub| sub.sg == 1));
        assert_eq!(ir.labels().cdml, ir.ml);
    }
    // SG 1 GM1 and Z1 share the BC label Γ1; both must survive as separate
    // source identities instead of being deduplicated by label.
    let shared_bc: Vec<_> = cdml
        .iter()
        .filter(|ir| ir.labels().bc.as_deref() == Some("Γ1"))
        .map(|ir| (ir.sg, ir.ml))
        .collect();
    assert!(shared_bc.contains(&(1, "GM1")));
    assert!(shared_bc.contains(&(1, "Z1")));

    let bc = irreps_with_subgroup(1, LabelConvention::Bc);
    assert_eq!(bc.len(), 839);
    assert_eq!(cdml.len() - bc.len(), 4, "SG 1's four *** rows lose BC");
    for ir in &bc {
        assert!(!ir.spinor);
        let label = ir.labels().bc.expect("BC rows carry a BC label");
        assert!(!label.contains("GM"));
        assert!(ir.subgroups().iter().any(|sub| sub.sg == 1));
        let hits = find_irreps(ir.sg, &label, LabelConvention::Bc);
        assert!(hits.iter().any(|hit| hit.id() == ir.id()));
        // Both conventions remain readable on the returned record.
        assert_eq!(ir.labels().cdml, ir.ml);
        assert_eq!(ir.k_labels().cdml, ir.k_label());
    }
    let bc_ids: BTreeSet<_> = bc.iter().map(|ir| ir.id()).collect();
    let expected: BTreeSet<_> = cdml
        .iter()
        .filter(|ir| ir.labels().bc.is_some())
        .map(|ir| ir.id())
        .collect();
    assert_eq!(bc_ids, expected);

    let magnetic_cdml = irreps_with_magnetic_subgroup(3, LabelConvention::Cdml);
    let magnetic_bc = irreps_with_magnetic_subgroup(3, LabelConvention::Bc);
    assert_eq!(magnetic_cdml.len(), 593);
    assert_eq!(magnetic_bc.len(), 589);
    assert_eq!(magnetic_cdml.len() - magnetic_bc.len(), 4);
    for ir in &magnetic_bc {
        assert!(ir.magnetic_subgroups().iter().any(|sub| sub.mag_sg == 3));
        assert!(ir.labels().bc.is_some());
    }
}

#[test]
fn bridge_methods_require_an_explicit_convention() {
    let crystal = cryspglib::Crystal::new(
        [[1.0, 0.0, 0.0], [0.3, 1.1, 0.0], [0.2, 0.4, 0.9]],
        vec![[0.1, 0.2, 0.3], [0.4, 0.15, 0.25]],
        vec![1, 2],
    )
    .expect("triclinic P1 crystal");
    let dataset = crystal.analyze().symprec(1e-5).dataset().expect("dataset");
    assert_eq!(dataset.spacegroup_number, 1);

    let cdml: Vec<_> = dataset
        .kpoints(LabelConvention::Cdml)
        .into_iter()
        .map(|kpoint| (kpoint.labels.cdml, kpoint.coords, kpoint.irreps))
        .collect();
    assert_eq!(cdml.len(), 8); // Eight physical k-points, shared by scalar/spinor rows.
    assert_eq!(
        cdml.iter()
            .map(|(_, _, irreps)| irreps.len())
            .sum::<usize>(),
        16
    );
    let bc: Vec<_> = dataset
        .kpoints(LabelConvention::Bc)
        .into_iter()
        .map(|kpoint| (kpoint.labels.cdml, kpoint.coords, kpoint.irreps))
        .collect();
    assert_eq!(bc.len(), 4);

    // CDML is literal; BC reads the stored prefix and skips the spinors.
    assert_eq!(
        mls(&dataset.irreps_at_k("GM", LabelConvention::Cdml)),
        vec!["GM1", "GM2"]
    );
    assert!(dataset.irreps_at_k("Γ", LabelConvention::Cdml).is_empty());
    assert_eq!(
        mls(&dataset.irreps_at_k("Z", LabelConvention::Cdml)),
        vec!["Z1", "Z2"]
    );
    assert!(dataset.irreps_at_k("Z", LabelConvention::Bc).is_empty());
    assert_eq!(
        mls(&dataset.irreps_at_k("Γ", LabelConvention::Bc)),
        vec!["GM1", "Z1"]
    );

    // The isotropy bridge keeps the same subgroup description per record; BC
    // returns both physical k-points that share the Γ symbol.
    let cdml_iso = dataset.irreps_with_isotropy_at_k("GM", LabelConvention::Cdml);
    assert_eq!(cdml_iso.len(), 1);
    assert_eq!(cdml_iso[0].0.ml, "GM1");
    assert!(cdml_iso[0].1.starts_with("#1 "));
    let bc_iso = dataset.irreps_with_isotropy_at_k("Γ", LabelConvention::Bc);
    assert_eq!(bc_iso.len(), 2);
    assert_eq!(bc_iso[0].0.ml, "GM1");
    assert_eq!(bc_iso[0].1, cdml_iso[0].1);
    assert_eq!(bc_iso[1].0.ml, "Z1");
}
