//! R25 card 1: the scalar Γ magnetic **input contract** is real.
//!
//! The model under test is the one the round prescribes: the ordinary parent `P`
//! and its grey group `P × {1, Θ}`, whose action is `ρ(p, Θ^t) = χ(p)·η^t·K^t` on
//! a fixed real one-dimensional Γ basis (`Θ = η·K`, so `Θ² = +1`).  A magnetic
//! subgroup `M` (Type I, III or IV) embeds into that grey parent by the verified
//! `φ` and need not contain `Θ`; the model is the restriction of `ρ` to `φ(M)`,
//! `A(g) = η^{t(g)}·χ(π(g))` for the mapped class `π(g)`.
//!
//! The tests pin: a Type I input accepted with and without `η` (the parent datum
//! is preserved either way); Type III inputs
//! accepted for **both** `η = ±1` with every class's value checked against the
//! parent's own character row read a second time in this file; a non-identity map
//! carried through to the model; BC and CDML label spellings of the same record
//! giving the same data; and the refusals (unknown parent/label, spinor,
//! non-scalar, non-Γ source, a missing `η` on a set that has a coset, an
//! unsatisfied lattice policy, a tampered map, an out-of-range row).  The corpus
//! sweep at the end walks every scalar Γ magnetic isotropy row and requires each
//! to end in a model or in a refusal that names a property of the *input*.

use cryspglib::irrep::magnetic_embedding::{
    self, MagneticOperationSet, ParentEmbedding, embed_in_parent_conventional, verify_embedding,
};
use cryspglib::irrep::magnetic_input::{
    ClassAction, EmbeddingEvidence, LatticePolicy, MagneticInputError, RealSign, RecordIdentity,
    ScalarGammaModel, ScalarGammaRequest, SourceIdentity, TimeParity, verify_scalar_gamma_input,
};
use cryspglib::irrep::subduction::{Mat3R, Rat};
use cryspglib::irrep::{LabelConvention, isotropy, query};

/// One source irrep plus one magnetic isotropy row of it, with the derived set.
struct Witness {
    sg: u8,
    source: &'static str,
    direction: usize,
    uni: usize,
    set: MagneticOperationSet,
    candidates: Vec<ParentEmbedding>,
}

impl Witness {
    fn request(
        &self,
        eta: Option<RealSign>,
        embedding: EmbeddingEvidence,
    ) -> ScalarGammaRequest<'static> {
        ScalarGammaRequest::new(
            SourceIdentity {
                parent_sg: self.sg,
                label: self.source,
                convention: LabelConvention::Cdml,
            },
            RecordIdentity {
                parent_sg: self.sg,
                condensing_label: self.source,
                convention: LabelConvention::Cdml,
                direction_index: self.direction,
            },
            eta,
        )
        .with_embedding(embedding)
    }

    fn geometry(&self) -> magnetic_embedding::MagneticGeometry {
        let row =
            isotropy::magnetic_isotropy_subgroups(self.sg, self.source, LabelConvention::Cdml)
                .expect("the label resolves")[self.direction];
        magnetic_embedding::geometry_of(self.sg, &row.record).expect("geometry")
    }

    /// The identity selection with the containment policy.
    fn identity(&self) -> EmbeddingEvidence {
        EmbeddingEvidence::identity()
    }

    /// Every candidate of the search, as evidence with the exact policy that
    /// matches the candidate's own lattice verdict.
    fn candidate_evidence(&self, candidate: &ParentEmbedding) -> EmbeddingEvidence {
        EmbeddingEvidence {
            map: candidate.map,
            shift: candidate.shift,
            lattice_policy: if candidate.realises_record_lattice {
                LatticePolicy::RequireExact
            } else {
                LatticePolicy::ContainmentOnly
            },
        }
    }
}

/// Collect scalar Γ magnetic isotropy rows up to `limit`, in
/// (space group, table order).  The embedding search is only run when
/// `with_candidates` is set: it is by far the most expensive part and most tests
/// only need the operation set.
fn witnesses(limit: usize, with_candidates: bool) -> Vec<Witness> {
    let mut out = Vec::new();
    for sg in 1..=230u8 {
        for irrep in query::irreps_of(sg) {
            if irrep.spinor || irrep.dim != 1 || irrep.k_vector().numerators != [0, 0, 0] {
                continue;
            }
            let rows = isotropy::magnetic_isotropy_subgroups(sg, irrep.ml, LabelConvention::Cdml)
                .unwrap_or_default();
            for (direction, row) in rows.iter().enumerate() {
                let Ok(set) = magnetic_embedding::magnetic_operations(row.record.mag_sg) else {
                    continue;
                };
                let Ok(geometry) = magnetic_embedding::geometry_of(sg, &row.record) else {
                    continue;
                };
                let candidates = if with_candidates {
                    embed_in_parent_conventional(&geometry, &set).unwrap_or_default()
                } else {
                    Vec::new()
                };
                out.push(Witness {
                    sg,
                    source: irrep.ml,
                    direction,
                    uni: set.uni,
                    set,
                    candidates,
                });
                if out.len() >= limit {
                    return out;
                }
            }
        }
    }
    out
}

/// A witness whose set is purely unitary (Type I) with the largest order found.
fn a_type_one_witness() -> Witness {
    let mut best: Option<Witness> = None;
    for witness in witnesses(usize::MAX, false) {
        if witness.set.unitary_index() != 1 {
            continue;
        }
        if best
            .as_ref()
            .map(|current| witness.set.unitary_classes > current.set.unitary_classes)
            .unwrap_or(true)
        {
            best = Some(witness);
        }
    }
    best.expect("the database has a purely unitary scalar Gamma magnetic row")
}

/// A Type III witness: two classes, one of them antiunitary.
fn a_type_three_witness() -> Witness {
    witnesses(usize::MAX, false)
        .into_iter()
        .find(|witness| witness.set.unitary_index() == 2 && witness.set.classes.len() == 2)
        .expect("the database has a two-class Type III scalar Gamma row")
}

/// The parent's character row of a witness's source, for independent checking.
fn parent_row(witness: &Witness) -> cryspglib::irrep::types::CharacterRow {
    query::irreps_of(witness.sg)
        .iter()
        .find(|irrep| irrep.ml == witness.source)
        .expect("the source resolves")
        .ordinary_scalar_selected_arm_block_trace()
        .expect("the source row is readable")
}

/// The source row's `±1` for a classification's rotation, or `None` if the row
/// does not cover that rotation.
fn row_sign(row: &cryspglib::irrep::types::CharacterRow, rotation: [[i32; 3]; 3]) -> Option<i32> {
    let flat = [
        rotation[0][0],
        rotation[0][1],
        rotation[0][2],
        rotation[1][0],
        rotation[1][1],
        rotation[1][2],
        rotation[2][0],
        rotation[2][1],
        rotation[2][2],
    ];
    row.operations()
        .iter()
        .zip(row.values())
        .find(|(operation, _)| operation.rotation == flat)
        .map(|(_, value)| if value.re > 0.0 { 1 } else { -1 })
}

/// The mapped `(rotation, translation)` of a class under an embedding.
fn mapped_operation(
    evidence: &EmbeddingEvidence,
    class: &cryspglib::irrep::magnetic_embedding::MagneticOperation,
) -> ([[i32; 3]; 3], Rat3) {
    let inverse = evidence.map.inverse().expect("a change of basis");
    let mapped_rotation = evidence
        .map
        .checked_mul(&Mat3R::from_ints(class.rotation))
        .and_then(|product| product.checked_mul(&inverse))
        .expect("a conjugated rotation")
        .to_int_matrix()
        .expect("integral");
    let moved = evidence
        .map
        .checked_mul_vector(&cryspglib::irrep::subduction::Vec3R::new(class.translation))
        .expect("a mapped translation");
    let mut translation = [Rat::ZERO; 3];
    for row in 0..3 {
        let mut sum = Rat::ZERO;
        for (column, shift) in evidence.shift.iter().enumerate() {
            let entry = i128::from(if row == column { 1 } else { 0 })
                - i128::from(mapped_rotation[row][column]);
            sum = sum
                .checked_add(
                    shift
                        .checked_mul(Rat::from_integer(entry))
                        .expect("small product"),
                )
                .expect("small sum");
        }
        translation[row] = moved.get(row).checked_add(sum).expect("small sum");
    }
    (mapped_rotation, translation)
}

type Rat3 = [Rat; 3];

/// Check every class of a model against the *literal* rule
/// `A(g) = η^{t(g)}·χ(π(g))`, with `χ` read independently from the source row.
fn assert_actions_match_the_rule(
    witness: &Witness,
    model: &ScalarGammaModel,
    evidence: &EmbeddingEvidence,
    eta: Option<RealSign>,
) {
    let row = parent_row(witness);
    for (index, action) in model.actions().iter().enumerate() {
        let class = &witness.set.classes[index];
        let (rotation, _) = mapped_operation(evidence, class);
        let chi = row_sign(&row, rotation).unwrap_or_else(|| {
            panic!(
                "SG {} {}: the source row does not cover the mapped rotation {rotation:?}",
                witness.sg, witness.source
            )
        });
        let expected = if class.time_reversal {
            eta.expect("an antiunitary class implies eta").value() * chi
        } else {
            chi
        };
        assert_eq!(
            action.value.value(),
            expected,
            "SG {} {}: class {index} must satisfy A(g) = eta^t * chi(mapped g)",
            witness.sg,
            witness.source
        );
        assert!(matches!(action.value, RealSign::Plus | RealSign::Minus));
        assert_eq!(action.class.rotation, class.rotation);
        assert_eq!(action.class.time_reversal, class.time_reversal);
        let _: &ClassAction = action;
    }
}

#[test]
fn a_type_one_input_is_accepted_without_claiming_a_time_reversal() {
    let witness = a_type_one_witness();
    assert_eq!(witness.set.unitary_index(), 1);
    assert!(witness.set.unitary_classes > 2);
    verify_embedding(
        &witness.geometry(),
        &witness.set,
        witness.identity().map,
        witness.identity().shift,
    )
    .expect("the canonical selection embeds the Type I row");

    // A Type I subgroup has no antiunitary class, so nothing is squared and no
    // square is claimed; the parent `eta` is still a datum and is preserved when
    // the caller states it.
    let plain = verify_scalar_gamma_input(&witness.request(None, witness.identity()))
        .unwrap_or_else(|error| {
            panic!("SG {} {} without eta: {error}", witness.sg, witness.source)
        });
    assert_eq!(plain.time_reversal(), None);
    assert_eq!(plain.time_parity(), None);
    assert_eq!(plain.time_reversal_square(), None);
    assert_eq!(plain.actions().len(), witness.set.classes.len());
    assert_actions_match_the_rule(&witness, &plain, &witness.identity(), None);

    for eta in [RealSign::Plus, RealSign::Minus] {
        let model = verify_scalar_gamma_input(&witness.request(Some(eta), witness.identity()))
            .unwrap_or_else(|error| {
                panic!("SG {} {} eta {eta:?}: {error}", witness.sg, witness.source)
            });
        assert_eq!(
            model.time_reversal(),
            Some(eta),
            "the parent eta is preserved"
        );
        assert_eq!(
            model.time_parity(),
            Some(if eta == RealSign::Plus {
                TimeParity::Even
            } else {
                TimeParity::Odd
            })
        );
        // The square is a property of the parent action's `Theta`, so a stated
        // `eta` reports `+1` even for a Type I subgroup; only the parity of the
        // restriction is uninformed here.
        assert_eq!(model.time_reversal_square(), Some(RealSign::Plus));
        assert_actions_match_the_rule(&witness, &model, &witness.identity(), Some(eta));
    }
}

#[test]
fn a_type_three_input_succeeds_for_both_signs_and_never_confuses_parity_with_the_square() {
    let witness = a_type_three_witness();
    assert_eq!(witness.set.antiunitary_classes, 1);
    let evidence = witness.identity();

    let even = verify_scalar_gamma_input(&witness.request(Some(RealSign::Plus), evidence))
        .unwrap_or_else(|error| panic!("eta = +1: {error}"));
    let odd = verify_scalar_gamma_input(&witness.request(Some(RealSign::Minus), evidence))
        .unwrap_or_else(|error| panic!("eta = -1: {error}"));

    assert_eq!(even.time_reversal(), Some(RealSign::Plus));
    assert_eq!(odd.time_reversal(), Some(RealSign::Minus));
    assert_eq!(even.time_parity(), Some(TimeParity::Even));
    assert_eq!(odd.time_parity(), Some(TimeParity::Odd));
    // `Theta^2 = +1` in this model, checked on the class the representative
    // squares to: the odd model is time-odd with a `+1` square, which is the
    // case that must not be reported as `T^2 = -1`.
    assert_eq!(even.time_reversal_square(), Some(RealSign::Plus));
    assert_eq!(odd.time_reversal_square(), Some(RealSign::Plus));

    assert_actions_match_the_rule(&witness, &even, &evidence, Some(RealSign::Plus));
    assert_actions_match_the_rule(&witness, &odd, &evidence, Some(RealSign::Minus));

    // Flipping `eta` is legal: it is the parent action's datum, and it must change
    // exactly the antiunitary classes of the restriction.
    let mut flipped = 0usize;
    for (index, class) in witness.set.classes.iter().enumerate() {
        let even_value = even.actions()[index].value.value();
        let odd_value = odd.actions()[index].value.value();
        if class.time_reversal {
            assert_eq!(odd_value, -even_value, "eta must flip class {index}");
            flipped += 1;
        } else {
            assert_eq!(odd_value, even_value, "eta must not touch class {index}");
        }
    }
    assert_eq!(flipped, witness.set.antiunitary_classes);
    assert!(flipped > 0);
}

#[test]
fn every_reachable_row_answers_final_the_contract_or_a_named_input_refusal() {
    // The sweep at the input layer: not a coverage claim, but a statement that
    // nothing panics and that no row is refused for an internal reason.  The
    // canonical (identity map, zero shift, containment) selection is a genuine
    // selection the caller may make; a row it does not embed is refused by the
    // embedding check, which is a property of that selection, not a bug.
    let mut rows = 0usize;
    let mut built = 0usize;
    let mut refusals: std::collections::BTreeMap<&'static str, usize> =
        std::collections::BTreeMap::new();
    for witness in witnesses(usize::MAX, false) {
        rows += 1;
        let eta = if witness.set.antiunitary_classes > 0 {
            Some(RealSign::Minus)
        } else {
            None
        };
        match verify_scalar_gamma_input(&witness.request(eta, witness.identity())) {
            Ok(model) => {
                built += 1;
                assert_eq!(model.actions().len(), witness.set.classes.len());
                assert_eq!(model.uni(), witness.uni);
                assert_actions_match_the_rule(&witness, &model, &witness.identity(), eta);
            }
            Err(error) => {
                let name = refusal_name(&error);
                // The allowed refusals are all properties of the *input*: a
                // selection that does not embed, a source outside the family, a
                // record index outside the table, or a value off the exact grid.
                // Anything else (`NotACorepresentation`, a group proof, an
                // arithmetic failure) would be a bug in the checker.
                assert!(
                    matches!(
                        name,
                        "embedding rejected"
                            | "missing time reversal"
                            | "lattice policy unsatisfied"
                            | "label not found"
                            | "ambiguous label"
                            | "direction out of range"
                            | "parent mismatch"
                            | "unknown parent"
                            | "source not Gamma"
                            | "not a scalar source"
                            | "spinor source"
                            | "off grid"
                            | "parent translation off grid"
                            | "character row setting mismatch"
                    ),
                    "SG {} {} row {}: refused for an internal reason: {error}",
                    witness.sg,
                    witness.source,
                    witness.direction
                );
                *refusals.entry(name).or_insert(0) += 1;
            }
        }
    }
    assert!(rows > 0, "the sweep visited no row");
    assert!(built > 0, "the sweep built no model at all");
}

/// The error's own name, for the sweep's histogram.
fn refusal_name(error: &MagneticInputError) -> &'static str {
    match error {
        MagneticInputError::UnknownParentSpaceGroup { .. } => "unknown parent",
        MagneticInputError::LabelNotFound { .. } => "label not found",
        MagneticInputError::AmbiguousLabel { .. } => "ambiguous label",
        MagneticInputError::SpinorSource { .. } => "spinor source",
        MagneticInputError::NotScalarSource { .. } => "not a scalar source",
        MagneticInputError::SourceNotGamma { .. } => "source not Gamma",
        MagneticInputError::ParentMismatch { .. } => "parent mismatch",
        MagneticInputError::DirectionOutOfRange { .. } => "direction out of range",
        MagneticInputError::DatabaseMismatch(_) => "database mismatch",
        MagneticInputError::GroupProof(_) => "group proof failed",
        MagneticInputError::ClassListInconsistent { .. } => "class list inconsistent",
        MagneticInputError::OffGrid { .. } => "off grid",
        MagneticInputError::ParentTranslationOffGrid { .. } => "parent translation off grid",
        MagneticInputError::NonRealCharacter { .. } => "non real character",
        MagneticInputError::IdentityCharacterNotOne { .. } => "identity character not one",
        MagneticInputError::CharacterRowSettingMismatch { .. } => "character row setting mismatch",
        MagneticInputError::EmbeddingRejected(_) => "embedding rejected",
        MagneticInputError::LatticePolicyUnsatisfied { .. } => "lattice policy unsatisfied",
        MagneticInputError::MissingTimeReversal { .. } => "missing time reversal",
        MagneticInputError::NotACorepresentation { .. } => "not a corepresentation",
        MagneticInputError::Arithmetic(_) => "arithmetic",
    }
}

#[test]
fn a_non_identity_verified_map_is_carried_through_to_the_model() {
    let witness = witnesses(usize::MAX, true)
        .into_iter()
        .find(|witness| {
            witness
                .candidates
                .iter()
                .any(|candidate| candidate.map != Mat3R::identity())
        })
        .expect("a witness with a non-identity verified map");
    let candidate = witness
        .candidates
        .iter()
        .find(|candidate| candidate.map != Mat3R::identity())
        .expect("the non-identity candidate");
    let evidence = witness.candidate_evidence(candidate);
    let model = verify_scalar_gamma_input(&witness.request(Some(RealSign::Minus), evidence))
        .unwrap_or_else(|error| panic!("SG {} {}: {error}", witness.sg, witness.source));

    // The evidence is kept, not summarized away.
    assert_eq!(model.map(), &candidate.map);
    assert_eq!(model.shift(), &candidate.shift);
    assert_eq!(model.uni(), witness.uni);
    assert_eq!(model.direction_index(), witness.direction);
    assert_eq!(model.parent_sg(), witness.sg);
    assert_eq!(model.source_label(), witness.source);
    assert_eq!(model.source_convention(), LabelConvention::Cdml);
    assert_eq!(model.record_label(), witness.source);
    assert_eq!(model.lattice_policy(), evidence.lattice_policy);
    assert_ne!(model.map(), &Mat3R::identity(), "a real non-identity map");
    assert_actions_match_the_rule(&witness, &model, &evidence, Some(RealSign::Minus));
}

#[test]
fn the_bc_and_cdml_spellings_of_one_record_give_the_same_data() {
    // Both spellings must resolve to a single irrep each, or the request is
    // ambiguous by design and not a spelling test.  The candidates are found by
    // trying the spelling on the source and on the record row and keeping only
    // what the resolver accepts.
    let unique = |sg: u8, label: &str, convention: LabelConvention| {
        query::irreps_of(sg)
            .iter()
            .filter(|candidate| candidate.label(convention).as_deref() == Some(label))
            .count()
            == 1
    };
    let mut found = None;
    'search: for sg in 1..=230u8 {
        for irrep in query::irreps_of(sg) {
            if irrep.spinor || irrep.dim != 1 || irrep.k_vector().numerators != [0, 0, 0] {
                continue;
            }
            let (Some(bc), Some(cdml)) = (
                irrep.label(LabelConvention::Bc),
                irrep.label(LabelConvention::Cdml),
            ) else {
                continue;
            };
            if bc == cdml {
                continue;
            }
            if !unique(sg, &bc, LabelConvention::Bc) || !unique(sg, &cdml, LabelConvention::Cdml) {
                continue;
            }
            // The isotropy query keys on the irrep's own CDML label; the request
            // may instead name the BC spelling, which the resolver must map to
            // the same irrep.
            let Ok(trial) = isotropy::magnetic_isotropy_subgroups(sg, &cdml, LabelConvention::Cdml)
            else {
                continue;
            };
            if trial.is_empty() {
                continue;
            }
            found = Some((sg, cdml.to_string(), bc, trial[0].record.mag_sg));
            break 'search;
        }
    }
    let (sg, cdml, bc, uni) = found.expect("a row with two usable label spellings");
    let eta = match magnetic_embedding::magnetic_operations(uni)
        .expect("operations")
        .unitary_index()
    {
        2 => Some(RealSign::Minus),
        _ => None,
    };
    let build = |label: &str, convention: LabelConvention| {
        verify_scalar_gamma_input(
            &ScalarGammaRequest::new(
                SourceIdentity {
                    parent_sg: sg,
                    label,
                    convention,
                },
                RecordIdentity {
                    parent_sg: sg,
                    condensing_label: label,
                    convention,
                    direction_index: 0,
                },
                eta,
            )
            .with_embedding(EmbeddingEvidence::identity()),
        )
    };
    let by_cdml = build(&cdml, LabelConvention::Cdml)
        .unwrap_or_else(|error| panic!("CDML spelling {cdml:?}: {error}"));
    let by_bc = build(&bc, LabelConvention::Bc)
        .unwrap_or_else(|error| panic!("BC spelling {bc:?}: {error}"));
    assert_eq!(by_cdml.uni(), by_bc.uni());
    assert_eq!(by_cdml.parent_sg(), by_bc.parent_sg());
    // Both report canonical labels under `Cdml`, so each pair resolves again;
    // what the caller typed is recorded separately as provenance.
    assert_eq!(by_cdml.source_label(), by_bc.source_label());
    assert_eq!(by_cdml.record_label(), by_bc.record_label());
    assert_eq!(by_cdml.direction_index(), by_bc.direction_index());
    assert_eq!(by_cdml.source_convention(), LabelConvention::Cdml);
    assert_eq!(by_bc.source_convention(), LabelConvention::Cdml);
    assert_eq!(by_cdml.record_convention(), LabelConvention::Cdml);
    assert_eq!(by_bc.record_convention(), LabelConvention::Cdml);
    assert_eq!(by_cdml.input_convention(), LabelConvention::Cdml);
    assert_eq!(by_bc.input_convention(), LabelConvention::Bc);
    for (label, convention) in [
        (by_bc.source_label(), by_bc.source_convention()),
        (by_bc.record_label(), by_bc.record_convention()),
    ] {
        let hits: Vec<_> = query::irreps_of(sg)
            .iter()
            .filter(|irrep| irrep.label(convention).as_deref() == Some(label))
            .collect();
        assert_eq!(hits.len(), 1, "the reported identity must resolve: {label}");
    }
    assert_eq!(by_cdml.actions().len(), by_bc.actions().len());
    for (left, right) in by_cdml.actions().iter().zip(by_bc.actions()) {
        assert_eq!(left.value, right.value);
        assert_eq!(left.class.rotation, right.class.rotation);
        assert_eq!(left.class.time_reversal, right.class.time_reversal);
    }
}

#[test]
fn input_identity_errors_are_specific() {
    let (sg, label, direction) = {
        let witness = witnesses(1, false).into_iter().next().expect("a witness");
        (witness.sg, witness.source, witness.direction)
    };
    let source = |sg, label: &'static str, convention| SourceIdentity {
        parent_sg: sg,
        label,
        convention,
    };
    let record = |sg, label: &'static str, convention, direction| RecordIdentity {
        parent_sg: sg,
        condensing_label: label,
        convention,
        direction_index: direction,
    };
    let cdml = LabelConvention::Cdml;

    let error = verify_scalar_gamma_input(&ScalarGammaRequest::new(
        source(0, "GM1", cdml),
        record(0, "GM1", cdml, 0),
        None,
    ))
    .expect_err("SG 0");
    assert!(matches!(
        error,
        MagneticInputError::UnknownParentSpaceGroup { sg: 0 }
    ));

    let error = verify_scalar_gamma_input(&ScalarGammaRequest::new(
        source(sg, "NO-SUCH-LABEL", cdml),
        record(sg, label, cdml, direction),
        None,
    ))
    .expect_err("no such label");
    assert!(
        matches!(error, MagneticInputError::LabelNotFound { .. }),
        "{error}"
    );

    let (spinor_sg, spinor_label) = first_irrep(|irrep| irrep.spinor);
    let error = verify_scalar_gamma_input(&ScalarGammaRequest::new(
        source(spinor_sg, spinor_label, cdml),
        record(spinor_sg, spinor_label, cdml, 0),
        None,
    ))
    .expect_err("spinor source");
    assert!(
        matches!(error, MagneticInputError::SpinorSource { .. }),
        "{error}"
    );

    let (big_sg, big_label) = first_irrep(|irrep| !irrep.spinor && irrep.dim > 1);
    let error = verify_scalar_gamma_input(&ScalarGammaRequest::new(
        source(big_sg, big_label, cdml),
        record(big_sg, big_label, cdml, 0),
        None,
    ))
    .expect_err("non-scalar source");
    assert!(
        matches!(error, MagneticInputError::NotScalarSource { .. }),
        "{error}"
    );

    let (k_sg, k_label) = first_irrep(|irrep| {
        !irrep.spinor && irrep.dim == 1 && irrep.k_vector().numerators != [0, 0, 0]
    });
    let error = verify_scalar_gamma_input(&ScalarGammaRequest::new(
        source(k_sg, k_label, cdml),
        record(k_sg, k_label, cdml, 0),
        None,
    ))
    .expect_err("source not at Gamma");
    assert!(
        matches!(error, MagneticInputError::SourceNotGamma { .. }),
        "{error}"
    );

    // A record index one past the table.
    let witness = a_type_three_witness();
    let rows =
        isotropy::magnetic_isotropy_subgroups(witness.sg, witness.source, cdml).unwrap_or_default();
    let error = verify_scalar_gamma_input(&ScalarGammaRequest::new(
        source(witness.sg, witness.source, cdml),
        record(witness.sg, witness.source, cdml, rows.len()),
        Some(RealSign::Plus),
    ))
    .expect_err("row index out of range");
    match error {
        MagneticInputError::DirectionOutOfRange { index, len, .. } => {
            assert_eq!(index, rows.len());
            assert_eq!(len, rows.len());
        }
        other => panic!("unexpected error {other}"),
    }

    // Different parents on the two identities.
    let error = verify_scalar_gamma_input(&ScalarGammaRequest::new(
        source(witness.sg, witness.source, cdml),
        record(1, "GM1", cdml, 0),
        Some(RealSign::Plus),
    ))
    .expect_err("parent mismatch");
    assert!(
        matches!(error, MagneticInputError::ParentMismatch { .. }),
        "{error}"
    );
}

fn first_irrep(
    predicate: impl Fn(&cryspglib::irrep::types::IrrepRecord) -> bool,
) -> (u8, &'static str) {
    for sg in 1..=230u8 {
        for irrep in query::irreps_of(sg) {
            if predicate(irrep) {
                return (sg, irrep.ml);
            }
        }
    }
    panic!("no irrep matches the predicate");
}

#[test]
fn a_coset_without_eta_is_refused_and_a_type_one_does_not_need_one() {
    let witness = a_type_three_witness();
    let error = verify_scalar_gamma_input(&witness.request(None, witness.identity()))
        .expect_err("an antiunitary coset needs eta");
    assert!(
        matches!(error, MagneticInputError::MissingTimeReversal { .. }),
        "{error}"
    );

    let type_one = a_type_one_witness();
    verify_scalar_gamma_input(&type_one.request(None, type_one.identity()))
        .expect("a Type I input needs no eta");
}

#[test]
fn a_tampered_selection_or_an_unsatisfied_lattice_policy_is_refused() {
    // The witness must have a non-identity class, or a shift would be invisible
    // (`(I - R) delta` vanishes for every rotation when the set is trivial) and
    // the mutation would not be a mutation.
    let witness = witnesses(usize::MAX, false)
        .into_iter()
        .find(|witness| {
            witness.set.unitary_index() == 2
                && witness
                    .set
                    .classes
                    .iter()
                    .any(|class| class.rotation != [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        })
        .expect("a Type III witness with a non-identity rotation");
    let good = witness.identity();
    verify_scalar_gamma_input(&witness.request(Some(RealSign::Plus), good))
        .expect("the canonical selection is accepted");

    // A shift moved off the solution set.
    let mut moved = good;
    moved.shift[0] = moved.shift[0]
        .checked_add(Rat::new(1, 8).expect("1/8"))
        .expect("small sum");
    let error = verify_scalar_gamma_input(&witness.request(Some(RealSign::Plus), moved))
        .expect_err("a moved shift is not the same embedding");
    assert!(
        matches!(
            error,
            MagneticInputError::EmbeddingRejected(_)
                | MagneticInputError::LatticePolicyUnsatisfied { .. }
        ),
        "unexpected error {error}"
    );

    // A map that doubles the cell on one axis is not a change of basis of the
    // record lattice, yet for this witness it still places every class on a parent
    // operation: the two policies separate exactly there.  `ContainmentOnly` is
    // "the proven embedding is enough" and accepts it, `RequireExact` adds the
    // record-lattice equality and refuses it — with the recomputed verdict carried
    // by the accepted model.
    let geometry = witness.geometry();
    let scaled = EmbeddingEvidence {
        map: Mat3R::diagonal([2, 1, 1]),
        shift: [Rat::ZERO; 3],
        lattice_policy: LatticePolicy::ContainmentOnly,
    };
    let exact =
        magnetic_embedding::realises_record_lattice_checked(&geometry, &witness.set, &scaled.map)
            .expect("the exact predicate runs");
    let accepted = verify_scalar_gamma_input(&witness.request(Some(RealSign::Plus), scaled))
        .unwrap_or_else(|error| panic!("the containment policy: {error}"));
    assert_eq!(accepted.lattice_policy(), LatticePolicy::ContainmentOnly);
    assert_eq!(
        accepted.realises_record_lattice(),
        exact,
        "the model carries the recomputed verdict whatever the policy"
    );
    if !exact {
        let exact_probe = EmbeddingEvidence {
            lattice_policy: LatticePolicy::RequireExact,
            ..scaled
        };
        let error = verify_scalar_gamma_input(&witness.request(Some(RealSign::Plus), exact_probe))
            .expect_err("the exact policy refuses the same map");
        assert!(
            matches!(error, MagneticInputError::LatticePolicyUnsatisfied { .. }),
            "unexpected error {error}"
        );
    }
}

#[test]
fn a_different_origin_is_aligned_without_discarding_seitz_translations() {
    // R35 replaces the old SG48 refusal with a joint, verified origin witness.
    let witness = witnesses(usize::MAX, false)
        .into_iter()
        .find(|witness| witness.sg == 48 && witness.source == "GM1+")
        .unwrap();
    let model =
        verify_scalar_gamma_input(&witness.request(Some(RealSign::Plus), witness.identity()))
            .unwrap();
    assert_eq!(*model.source_origin_shift(), [Rat::new(1, 4).unwrap(); 3]);
    assert_eq!(*model.shift(), [Rat::ZERO; 3]);
    verify_embedding(
        &witness.geometry(),
        &witness.set,
        witness.identity().map,
        witness.identity().shift,
    )
    .unwrap();

    // The rotation sets agree while the operation sets do not: that is the
    // measured shape of the boundary, restated here on one row.
    let row = parent_row(&witness);
    let row_rotations: std::collections::BTreeSet<[i32; 9]> =
        row.operations().iter().map(|op| op.rotation).collect();
    let parent = query::symmetry_operations_of(48).expect("SG 48 operations");
    let parent_rotations: std::collections::BTreeSet<[i32; 9]> = parent
        .operations
        .iter()
        .map(|op| {
            let r = op.rotation;
            [
                r[0][0], r[0][1], r[0][2], r[1][0], r[1][1], r[1][2], r[2][0], r[2][1], r[2][2],
            ]
        })
        .collect();
    assert_eq!(row_rotations, parent_rotations);
    let row_operations: std::collections::BTreeSet<([i32; 9], [i64; 3])> = row
        .operations()
        .iter()
        .map(|op| {
            (
                op.rotation,
                [
                    (op.translation[0] * 12.0).round() as i64,
                    (op.translation[1] * 12.0).round() as i64,
                    (op.translation[2] * 12.0).round() as i64,
                ],
            )
        })
        .collect();
    let parent_operations: std::collections::BTreeSet<([i32; 9], [i64; 3])> = parent
        .operations
        .iter()
        .map(|op| {
            let r = op.rotation;
            (
                [
                    r[0][0], r[0][1], r[0][2], r[1][0], r[1][1], r[1][2], r[2][0], r[2][1], r[2][2],
                ],
                [
                    (op.translation[0] * 12.0).round() as i64,
                    (op.translation[1] * 12.0).round() as i64,
                    (op.translation[2] * 12.0).round() as i64,
                ],
            )
        })
        .collect();
    assert_ne!(
        row_operations, parent_operations,
        "this witness is pinned because the two operation sets differ"
    );
}

#[test]
fn no_scalar_gamma_row_of_this_database_has_a_complex_character() {
    // The family is the real +-1 one and the module refuses a non-real character
    // by name instead of snapping it.  Measured across the whole database: no
    // one-dimensional scalar Γ source row is complex, so that refusal is a
    // contract on future data; this scan is what makes the statement falsifiable.
    let mut scanned = 0usize;
    for sg in 1..=230u8 {
        for irrep in query::irreps_of(sg) {
            if irrep.spinor || irrep.dim != 1 || irrep.k_vector().numerators != [0, 0, 0] {
                continue;
            }
            let Ok(row) = irrep.ordinary_scalar_selected_arm_block_trace() else {
                continue;
            };
            scanned += 1;
            for value in row.values() {
                assert!(
                    value.im.abs() < 1e-9,
                    "SG {sg} {}: the source row is complex at {value:?}, which this family refuses",
                    irrep.ml
                );
            }
        }
    }
    assert_eq!(
        scanned, 980,
        "scalar Gamma source rows in the shipped tables"
    );
}

#[test]
fn the_number_of_magnetic_gamma_rows_with_a_coset_is_pinned() {
    // A measurement, not a coverage claim: it fixes what the sweep above walks,
    // so a change in the tables or in the filter shows up here.
    let mut rows = 0usize;
    let mut with_coset = 0usize;
    for witness in witnesses(usize::MAX, false) {
        rows += 1;
        if witness.set.unitary_index() == 2 {
            with_coset += 1;
        }
    }
    assert_eq!(rows, 980, "scalar Gamma magnetic isotropy rows");
    assert_eq!(
        with_coset, 750,
        "of those, rows whose set has an antiunitary coset"
    );
}
