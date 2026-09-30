//! R8 acceptance: the supported ordinary-subduction entry point.
//!
//! These tests pin the promises of `cryspglib::irrep::subduction_api`:
//! explicit label conventions, both label spellings in the report, provenance
//! for every input, an explicit parameter semantics, and errors instead of a
//! silent identity-only answer.

use cryspglib::irrep::subduction::SubductionError;
use cryspglib::irrep::subduction::star::decompose::ParameterKind;
use cryspglib::irrep::subduction_api::{
    SubductionApiError, SubductionRequest, subduce, subduce_table, subduction_for_direction,
};
use cryspglib::irrep::{
    LabelConvention, isotropy, isotropy::IsotropyDirection, line_monodromy, query, subduce_irrep,
};

/// The milestone's golden case: SG 221 `GM4+` along `P1`.
#[test]
fn the_gm4_plus_p1_golden_case_reports_its_provenance() {
    let report = subduction_for_direction(
        221,
        "GM4+",
        LabelConvention::Cdml,
        IsotropyDirection::Label("P1"),
    )
    .expect("SG 221 GM4+ has a P1 isotropy subgroup");

    assert_eq!(report.condensing().parent_sg, 221);
    assert_eq!(report.condensing().cdml, "GM4+");
    assert_eq!(report.condensing().dimension, 3);
    assert!(
        report.condensing().bc.is_some(),
        "SG 221 GM4+ carries a BC spelling"
    );
    assert_eq!(report.direction().label, "P1");
    assert!(report.direction().subgroup_sg >= 1 && report.direction().subgroup_sg <= 230);
    assert!(report.parameter().is_none());
    assert!(report.parameter_kind().is_none());
    assert!(report.covered_dimension() > 0);
    assert!(!report.blocks().is_empty());

    // The reported provenance must re-find the same record through a different
    // selector: the table index and the ISOTROPY label agree.
    let by_index = subduce(
        &SubductionRequest::new(221, "GM4+", LabelConvention::Cdml)
            .direction(IsotropyDirection::Index(report.direction().index_in_irrep)),
    )
    .expect("the reported index selects the same record");
    assert_eq!(by_index.direction(), report.direction());

    // Every target either carries a pinned CDML label or is explicitly
    // constructed; the report never invents a label.
    let mut targets = 0usize;
    for block in report.blocks() {
        for target in &block.targets {
            targets += 1;
            assert!(
                target.cdml.is_some() || target.row_ml.is_none(),
                "a labelled row must report its CDML label"
            );
        }
    }
    assert!(targets > 0, "the golden case has at least one target");
}

/// The dimension-2 probe named by the milestone.
#[test]
fn the_probe_gm3_plus_case_works() {
    let report = subduction_for_direction(
        221,
        "GM3+",
        LabelConvention::Cdml,
        IsotropyDirection::Label("P1"),
    )
    .expect("SG 221 GM3+ has a P1 isotropy subgroup");
    assert_eq!(report.condensing().cdml, "GM3+");
    assert_eq!(report.condensing().dimension, 2);
    assert!(report.covered_dimension() > 0);
}

/// The milestone's cross-convention case: SG 213 CDML `X2` and BC `X1` are the
/// same record and must select the same context **and the `#146` geometry**
/// (child space group 146, which is the `C23` direction of this irrep).
#[test]
fn cdml_x2_and_bc_x1_select_the_same_context() {
    let cdml = subduce(
        &SubductionRequest::new(213, "X2", LabelConvention::Cdml)
            .direction(IsotropyDirection::Label("C23")),
    )
    .expect("CDML X2 resolves");
    let bc = subduce(
        &SubductionRequest::new(213, "X1", LabelConvention::Bc)
            .direction(IsotropyDirection::Label("C23")),
    )
    .expect("BC X1 resolves");
    assert_eq!(cdml.direction().subgroup_sg, 146, "the #146 geometry");
    assert_eq!(bc.direction().subgroup_sg, 146);
    assert_eq!(cdml.condensing().cdml, bc.condensing().cdml);
    assert_eq!(cdml.direction(), bc.direction());
    assert_eq!(cdml.covered_dimension(), bc.covered_dimension());
    assert_eq!(cdml.blocks(), bc.blocks());
}

/// A BC spelling that names more than one record is an input error, never a
/// first-hit guess.  The scan finds a real duplicate in the generated tables.
#[test]
fn a_duplicated_bc_label_is_refused_as_ambiguous() {
    let duplicate = (1..=230u8).find_map(|sg| {
        let records = query::irreps_of(sg);
        let mut groups: Vec<(&str, Vec<&'static str>)> = Vec::new();
        for record in records {
            let Some(bc) = record.label(LabelConvention::Bc) else {
                continue;
            };
            match groups.iter_mut().find(|(label, _)| *label == bc) {
                Some((_, members)) => members.push(record.ml),
                None => groups.push((
                    Box::leak(bc.into_boxed_str()),
                    vec![record.ml],
                )),
            }
        }
        groups
            .into_iter()
            .find(|(_, members)| members.len() > 1)
            .map(|(label, members)| (sg, label.to_string(), members))
    });
    let (sg, label, members) =
        duplicate.expect("the generated tables contain a duplicated BC label");
    let error = subduce(
        &SubductionRequest::new(sg, &label, LabelConvention::Bc)
            .direction(IsotropyDirection::Label("P1")),
    )
    .expect_err("an ambiguous BC label must be refused");
    match error {
        SubductionApiError::AmbiguousCondensing {
            sg: reported,
            candidates,
            ..
        } => {
            assert_eq!(reported, sg);
            assert!(
                candidates.len() > 1,
                "the error must list the candidates, got {candidates:?}"
            );
        }
        other => panic!("expected AmbiguousCondensing, got {other}"),
    }
    // Selecting one of them by its CDML label still works.
    let record = subduce(
        &SubductionRequest::new(sg, members[0], LabelConvention::Cdml)
            .direction(IsotropyDirection::Label("P1")),
    )
    .expect("the CDML label of a candidate resolves");
    assert_eq!(record.condensing().cdml, members[0]);
}

/// Unknown names and directions are errors that name the input.
#[test]
fn unknown_labels_and_directions_are_errors() {
    let unknown_label = subduce(
        &SubductionRequest::new(221, "NOPE1", LabelConvention::Cdml)
            .direction(IsotropyDirection::Label("P1")),
    )
    .expect_err("an unknown CDML label is an error");
    assert!(matches!(
        unknown_label,
        SubductionApiError::CondensingNotFound { .. }
    ));

    let unknown_direction = subduce(
        &SubductionRequest::new(221, "GM4+", LabelConvention::Cdml)
            .direction(IsotropyDirection::Label("NOPE1")),
    )
    .expect_err("an unknown direction is an error");
    assert!(matches!(
        unknown_direction,
        SubductionApiError::Isotropy(_)
    ));

    let unknown_descriptor = subduce(
        &SubductionRequest::new(221, "GM4+", LabelConvention::Cdml)
            .direction(IsotropyDirection::Descriptor("(nope)")),
    )
    .expect_err("an unknown descriptor is an error");
    assert!(matches!(
        unknown_descriptor,
        SubductionApiError::Isotropy(_)
    ));
}

/// A parameter is only meaningful on the frozen parametric-k lines, and when it
/// is accepted the report says which kind of parameter it was.
#[test]
fn a_parameter_is_only_accepted_on_the_frozen_lines() {
    let official = cryspglib::irrep::subduction::star::decompose::official_line_parameter()
        .expect("the official parameter parses");
    // The frozen source label and the condensing irrep label are two different
    // spellings of the same request: `DT1` is the source, `W1` the irrep whose
    // subgroup carries it.
    // The source label alone is ambiguous: SG 196 `DT1` runs through two
    // subgroups at `4D1` (W1 and W2), so the request is refused with the
    // candidate ordinals and one of them must be pinned.
    let ambiguous = subduce(
        &SubductionRequest::new(196, "DT1", LabelConvention::Cdml)
            .direction(IsotropyDirection::Label("4D1"))
            .parameter(official),
    )
    .expect_err("a source through two subgroups must be refused");
    let candidates = match ambiguous {
        SubductionApiError::AmbiguousLineContext {
            sg: 196,
            source,
            candidates,
        } => {
            assert_eq!(source, "DT1");
            candidates
        }
        other => panic!("expected AmbiguousLineContext, got {other}"),
    };
    assert_eq!(candidates, vec![10038, 10062]);

    let by_source = subduce(
        &SubductionRequest::new(196, "DT1", LabelConvention::Cdml)
            .direction(IsotropyDirection::Label("4D1"))
            .parameter(official)
            .ordinal(10038),
    )
    .expect("SG 196 DT1 is a frozen line source");
    assert_eq!(by_source.parameter(), Some(official));
    assert_eq!(by_source.parameter_kind(), Some(ParameterKind::LineIrrep));
    assert_eq!(by_source.condensing().cdml, "W1");
    assert_eq!(by_source.direction().subgroup_sg, 1);
    assert_eq!(by_source.direction().ordinal, 10038);

    // Naming the irrep is *not* enough here: three frozen lines (DT1, DT2,
    // SM1) pass through SG 196 W1 4D1, so the request is refused as ambiguous
    // rather than answering for an arbitrary one of them.
    let by_irrep = subduce(
        &SubductionRequest::new(196, "W1", LabelConvention::Cdml)
            .direction(IsotropyDirection::Label("4D1"))
            .parameter(official),
    )
    .expect_err("an irrep label that carries three sources is ambiguous");
    match by_irrep {
        SubductionApiError::AmbiguousLineSource {
            sg: 196,
            candidates,
            ..
        } => assert_eq!(candidates, vec!["DT1", "DT2", "SM1"]),
        other => panic!("expected AmbiguousLineSource, got {other}"),
    }

    let not_a_line = subduce(
        &SubductionRequest::new(221, "GM4+", LabelConvention::Cdml)
            .direction(IsotropyDirection::Label("P1"))
            .parameter(official),
    )
    .expect_err("a non-line irrep takes no parameter");
    assert!(matches!(
        not_a_line,
        SubductionApiError::NotAParameterizedLine { sg: 221, .. }
    ));
}

/// A condensing irrep the engine cannot decompose is an error, not an
/// identity-only report.
#[test]
fn a_spinor_condensing_irrep_is_refused_not_downgraded() {
    let spinor = (1..=230u8).find_map(|sg| {
        query::irreps_of(sg)
            .iter()
            .find(|record| record.spinor)
            .map(|record| (sg, record.ml))
    });
    let (sg, ml) = spinor.expect("the generated tables contain spinor records");
    let error = subduce(
        &SubductionRequest::new(sg, ml, LabelConvention::Cdml)
            .direction(IsotropyDirection::Label("P1")),
    )
    .expect_err("a spinor condensing irrep is refused");
    // The refusal may come from any stage (the isotropy tables do not define
    // subgroups for double-valued irreps), but it must be an error and it must
    // name the reason -- never an identity-only report.
    let message = error.to_string();
    assert!(
        message.contains("double-valued"),
        "the refusal must name the reason, got {message:?}"
    );
}

/// The whole-table wrapper answers every direction in table order, and a
/// parameterized sweep reports per entry which directions carry that line.
#[test]
fn the_whole_table_wrapper_covers_every_direction() {
    let entries = subduce_table(&SubductionRequest::new(221, "GM4+", LabelConvention::Cdml))
        .expect("SG 221 GM4+ has an isotropy table");
    assert_eq!(entries.len(), 4, "GM4+ has four isotropy directions");
    for (position, entry) in entries.iter().enumerate() {
        assert_eq!(entry.index, position, "entries are in table order");
        assert!(
            !entry.label.is_empty(),
            "every entry carries its direction label"
        );
        let report = entry
            .result
            .as_ref()
            .unwrap_or_else(|error| panic!("direction {} ({}) failed: {error}", entry.index, entry.label));
        assert_eq!(report.direction().index_in_irrep, entry.index);
        assert_eq!(report.direction().label, entry.label);
    }

    // A parameterized sweep: only the directions that carry the frozen line
    // answer, the others say why, and nothing is dropped silently.
    let official = cryspglib::irrep::subduction::star::decompose::official_line_parameter()
        .expect("the official parameter parses");
    let sweep = subduce_table(
        &SubductionRequest::new(196, "DT1", LabelConvention::Cdml)
            .parameter(official)
            .ordinal(10038),
    )
    .expect("the W1 table resolves through the source label and a pinned context");
    assert!(!sweep.is_empty());
    let answered = sweep.iter().filter(|entry| entry.result.is_ok()).count();
    assert!(
        answered > 0,
        "at least one direction of W1 carries the DT1 line"
    );
    assert!(
        answered < sweep.len(),
        "the DT1 line does not run along every direction of W1"
    );
    for entry in &sweep {
        match &entry.result {
            Ok(report) => {
                // The sweep was pinned to ordinal 10038, i.e. to SG 196 `W1`:
                // every answered entry must be that record's own subgroup at
                // that direction.  Round-15e found the opposite (entries
                // reporting `W2` with a foreign `index_in_irrep`) when the
                // sweep re-resolved the source label per direction.
                assert_eq!(
                    report.condensing().cdml,
                    "W1",
                    "entry {} answered with {}",
                    entry.index,
                    report.condensing().cdml
                );
                assert_eq!(report.direction().index_in_irrep, entry.index);
                assert_eq!(report.direction().label, entry.label);
            }
            Err(error) => {
                let message = error.to_string();
                assert!(
                    message.contains("frozen line source")
                        || message.contains("frozen line")
                        || message.contains("ambiguous"),
                    "an unanswered direction must say why, got {message:?}"
                );
            }
        }
    }
}

/// The R8 review's F1 witness: one source along one direction can run through
/// two subgroups with **different** decompositions.  The source label is then
/// ambiguous (never a first hit), and both contexts are reachable by ordinal --
/// including the one the first version could not reach at all.
#[test]
fn a_source_through_two_subgroups_is_ambiguous_and_both_are_reachable() {
    let official = cryspglib::irrep::subduction::star::decompose::official_line_parameter()
        .expect("the official parameter parses");
    let request = SubductionRequest::new(196, "DT1", LabelConvention::Cdml)
        .direction(IsotropyDirection::Label("4D2"))
        .parameter(official);

    let candidates = match subduce(&request).expect_err("two contexts must be refused") {
        SubductionApiError::AmbiguousLineContext {
            sg: 196, candidates, ..
        } => candidates,
        other => panic!("expected AmbiguousLineContext, got {other}"),
    };
    assert_eq!(candidates, vec![10049, 10072]);

    let first = subduce(&request.ordinal(10049)).expect("the W1 context answers");
    let second = subduce(&request.ordinal(10072)).expect("the W2 context answers");
    assert_eq!(first.condensing().cdml, "W1");
    assert_eq!(second.condensing().cdml, "W2");
    assert_eq!(first.direction().ordinal, 10049);
    assert_eq!(second.direction().ordinal, 10072);
    assert_ne!(
        first.direction().subgroup_sg, second.direction().subgroup_sg,
        "the two contexts are different subgroups"
    );
    assert_ne!(
        first.blocks(),
        second.blocks(),
        "the two contexts decompose differently, which is why the first hit was wrong"
    );

    // An ordinal that belongs to neither context is an error, not a fallback.
    match subduce(&request.ordinal(1)).expect_err("an unrelated ordinal is an error") {
        SubductionApiError::UnknownContext { sg: 196, ordinal: 1 } => {}
        other => panic!("expected UnknownContext, got {other}"),
    }
}

/// F2: a frozen source that exists in the parent but not along the selected
/// direction says exactly that, instead of claiming the label does not exist.
#[test]
fn a_source_missing_from_the_direction_says_so() {
    let official = cryspglib::irrep::subduction::star::decompose::official_line_parameter()
        .expect("the official parameter parses");
    // `P1` is a direction of SG 196's irreps along which no subgroup carries
    // the DT1 source (the source exists in this parent, just not here).
    let error = subduce(
        &SubductionRequest::new(196, "DT1", LabelConvention::Cdml)
            .direction(IsotropyDirection::Label("P1"))
            .parameter(official),
    )
    .expect_err("no SG 196 subgroup carries DT1 along P1");
    let message = error.to_string();
    assert!(
        message.contains("does not run along"),
        "the refusal must name the real reason, got {message:?}"
    );
}

/// The milestone's "newly added gap case": a request the **old** public entry
/// could not answer at all.
///
/// `subduce_irrep` is the Γ-only convenience surface: it refuses any condensing
/// irrep away from Γ before attempting a decomposition.  SG 213 `X2` sits at
/// k = (0,1,0)/2, so the old entry refuses the very context the new one
/// answers -- with the full star, both label spellings and the `#146` geometry.
#[test]
fn the_old_gamma_only_entry_refused_what_the_new_entry_answers() {
    let subgroups = isotropy::isotropy_subgroups(213, "X2", LabelConvention::Cdml)
        .expect("SG 213 X2 has an isotropy table");
    let subgroup = subgroups
        .iter()
        .find(|subgroup| subgroup.record.direction_label == "C23")
        .expect("X2 has the C23 direction");
    assert_eq!(subgroup.record.sg, 146);

    // The old surface: refused, because the irrep is not at Γ.
    let old = subduce_irrep(subgroup, "X2");
    assert!(
        matches!(old, Err(SubductionError::ProbeNotAtGamma { .. })),
        "the Γ-only entry must refuse a non-Γ irrep, got {old:?}"
    );

    // The supported entry: answered, with the full decomposition.
    let report = subduce(
        &SubductionRequest::new(213, "X2", LabelConvention::Cdml)
            .direction(IsotropyDirection::Label("C23")),
    )
    .expect("the supported entry answers a non-Γ irrep");
    assert_eq!(report.condensing().cdml, "X2");
    assert_eq!(report.condensing().k, [[0, 1, 0], [2, 2, 2]]);
    assert_eq!(report.covered_dimension(), u32::from(report.condensing().dimension));
    assert!(!report.blocks().is_empty());

    // The parameterized facet had no public entry at all before R8: the frozen
    // line is answered with its parameter kind.
    let official = cryspglib::irrep::subduction::star::decompose::official_line_parameter()
        .expect("the official parameter parses");
    let line = subduce(
        &SubductionRequest::new(196, "DT1", LabelConvention::Cdml)
            .direction(IsotropyDirection::Label("4D1"))
            .parameter(official)
            .ordinal(10038),
    )
    .expect("the frozen line answers at its official parameter");
    assert_eq!(line.parameter_kind(), Some(ParameterKind::LineIrrep));
}

/// Round-15e: a sweep pinned by a frozen **source** label must answer only the
/// resolved record's own subgroups.  Without the per-entry ordinal, a direction
/// whose own subgroup does not carry the source re-resolved the name and
/// answered with **another irrep's** context (SG 196 `DT1` `.ordinal(10038)`
/// reported `W2` for five directions), which is the worst kind of wrong answer:
/// a complete decomposition of the wrong subgroup.
#[test]
fn a_source_pinned_sweep_never_answers_with_another_record() {
    let mut checked = 0usize;
    let mut wrong: Vec<String> = Vec::new();
    let mut sources = 0usize;
    for parent in line_monodromy::line_parents() {
        for source in line_monodromy::line_sources(parent) {
            // The first context (record + row) whose rows carry this source.
            let mut found = None;
            'scan: for record in query::irreps_of(parent) {
                if record.spinor || record.subgroups().is_empty() {
                    continue;
                }
                let Ok(list) =
                    isotropy::isotropy_subgroups(parent, record.ml, LabelConvention::Cdml)
                else {
                    continue;
                };
                for subgroup in list {
                    let Ok(rows) = subgroup.other_wave_vector_subduction() else {
                        continue;
                    };
                    if rows.iter().any(|row| row.parent_ml == source.label) {
                        found = Some((record.ml, subgroup.ordinal));
                        break 'scan;
                    }
                }
            }
            let Some((owning, ordinal)) = found else {
                continue;
            };
            sources += 1;
            let entries = subduce_table(
                &SubductionRequest::new(parent, source.label, LabelConvention::Cdml)
                    .ordinal(ordinal),
            )
            .unwrap_or_else(|error| {
                panic!("sg {parent} {} ordinal {ordinal}: {error}", source.label)
            });
            for entry in &entries {
                let Ok(report) = &entry.result else {
                    continue;
                };
                checked += 1;
                if report.condensing().cdml != owning
                    || report.direction().index_in_irrep != entry.index
                {
                    wrong.push(format!(
                        "sg {parent} {} ordinal {ordinal}: entry {} ({}) answered {} index {}",
                        source.label,
                        entry.index,
                        entry.label,
                        report.condensing().cdml,
                        report.direction().index_in_irrep
                    ));
                }
            }
        }
    }
    assert!(sources > 0, "the frozen sources must be reachable");
    assert_eq!(wrong, Vec::<String>::new(), "{} wrong entries", wrong.len());
    assert!(checked > 0, "the sweep must answer something");
}
