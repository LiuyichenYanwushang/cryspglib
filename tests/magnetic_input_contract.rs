//! Independent public-contract regressions: identities and translation-only magnetism.
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, ScalarGammaRequest, SourceIdentity,
    verify_scalar_gamma_input,
};
use cryspglib::irrep::subduction::{Mat3R, Rat};
use cryspglib::irrep::{LabelConvention, isotropy, query};

#[test]
fn reported_label_and_convention_resolve_to_the_selected_irrep() {
    let convention = LabelConvention::Bc;
    let (sg, irrep, model) = (1..=230)
        .find_map(|sg| {
            query::irreps_of(sg).iter().find_map(|irrep| {
                if irrep.spinor || irrep.dim != 1 || irrep.k_vector().numerators != [0; 3] {
                    return None;
                }
                let label = irrep.label(convention)?;
                if label == irrep.ml {
                    return None;
                }
                let request = ScalarGammaRequest::new(
                    SourceIdentity {
                        parent_sg: sg,
                        label: &label,
                        convention,
                    },
                    RecordIdentity {
                        parent_sg: sg,
                        condensing_label: &label,
                        convention,
                        direction_index: 0,
                    },
                    Some(RealSign::Minus),
                );
                verify_scalar_gamma_input(&request)
                    .ok()
                    .map(|model| (sg, irrep, model))
            })
        })
        .expect("an unambiguous accepted BC request");
    for (reported, spelling) in [
        (model.source_label(), model.source_convention()),
        (model.record_label(), model.record_convention()),
    ] {
        let selected: Vec<_> = query::irreps_of(sg)
            .iter()
            .filter(|ir| ir.label(spelling).as_deref() == Some(reported))
            .collect();
        assert_eq!(
            selected.len(),
            1,
            "identity must be resolvable: {reported} {spelling:?}"
        );
        assert_eq!(selected[0].ml, irrep.ml);
    }
}

#[test]
fn gamma_source_on_type_four_record_preserves_the_time_reversal_coset() {
    // SG1 X1 condenses into UNI3 with an antiunitary half translation.
    // The representation being restricted is instead the trivial Gamma source.
    let rows = isotropy::magnetic_isotropy_subgroups(1, "X1", LabelConvention::Cdml).unwrap();
    let index = rows.iter().position(|row| row.record.mag_sg == 3).unwrap();
    let source = SourceIdentity {
        parent_sg: 1,
        label: "GM1",
        convention: LabelConvention::Cdml,
    };
    let record = RecordIdentity {
        parent_sg: 1,
        condensing_label: "X1",
        convention: LabelConvention::Cdml,
        direction_index: index,
    };
    let two = Mat3R::from_ints([[2, 0, 0], [0, 2, 0], [0, 0, 2]]);
    let make = |eta| {
        ScalarGammaRequest::new(source, record, Some(eta)).with_embedding(EmbeddingEvidence {
            map: two,
            shift: [Rat::ZERO; 3],
            lattice_policy: LatticePolicy::ContainmentOnly,
        })
    };
    for eta in [RealSign::Plus, RealSign::Minus] {
        let model = verify_scalar_gamma_input(&make(eta)).unwrap();
        assert_eq!(model.uni(), 3);
        assert!(!model.realises_record_lattice());
        assert_eq!(model.time_reversal_square(), Some(RealSign::Plus));
        assert!(model.actions().iter().any(|a| a.class.time_reversal));
        for action in model.actions() {
            assert_eq!(
                action.value,
                if action.class.time_reversal {
                    eta
                } else {
                    RealSign::Plus
                }
            );
        }
    }
    let exact = make(RealSign::Minus).with_embedding(EmbeddingEvidence {
        map: two,
        shift: [Rat::ZERO; 3],
        lattice_policy: LatticePolicy::RequireExact,
    });
    assert!(matches!(
        verify_scalar_gamma_input(&exact),
        Err(cryspglib::irrep::magnetic_input::MagneticInputError::LatticePolicyUnsatisfied { .. })
    ));
}
