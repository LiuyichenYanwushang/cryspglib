//! R8 acceptance: the supported ordinary-subduction entry point.
//!
//! These tests pin the promises of `cryspglib::irrep::subduction_api`:
//! explicit label conventions, both label spellings in the report, provenance
//! for every input, an explicit parameter semantics, and errors instead of a
//! silent identity-only answer.

use cryspglib::irrep::subduction::star::decompose::ParameterKind;
use cryspglib::irrep::subduction_api::{
    SubductionApiError, SubductionRequest, subduce, subduction_for_direction,
};
use cryspglib::irrep::{LabelConvention, isotropy::IsotropyDirection, query};

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
/// same record, so the two requests must select the same context.
#[test]
fn cdml_x2_and_bc_x1_select_the_same_context() {
    let cdml = subduce(
        &SubductionRequest::new(213, "X2", LabelConvention::Cdml)
            .direction(IsotropyDirection::Label("P1")),
    )
    .expect("CDML X2 resolves");
    let bc = subduce(
        &SubductionRequest::new(213, "X1", LabelConvention::Bc)
            .direction(IsotropyDirection::Label("P1")),
    )
    .expect("BC X1 resolves");
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
    let by_source = subduce(
        &SubductionRequest::new(196, "DT1", LabelConvention::Cdml)
            .direction(IsotropyDirection::Label("4D1"))
            .parameter(official),
    )
    .expect("SG 196 DT1 is a frozen line source");
    assert_eq!(by_source.parameter(), Some(official));
    assert_eq!(by_source.parameter_kind(), Some(ParameterKind::LineIrrep));
    assert_eq!(by_source.condensing().cdml, "W1");
    assert_eq!(by_source.direction().subgroup_sg, 1);

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
