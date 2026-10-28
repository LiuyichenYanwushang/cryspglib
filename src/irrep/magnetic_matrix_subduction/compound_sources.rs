//! Explicit physical compound-source admission; no new target family or fallback.
use super::*;
use crate::irrep::types::CompoundSelectedArmCharacter;

/// Source/real-A boundary for the explicit compound input, not a complex irreducible source.
pub const COMPOUND_SCOPE: &str = "Certified real non-spinor two-dimensional Gamma physical compound carrier of two distinct authoritative nonreal conjugate scalar CIRs; explicit real type-A target entry. Both CIR laws and physical matrix/frame evidence retained; no type-B, spinor, non-Gamma, higher-dimensional source or setting search.";
/// Scalar-A target boundary for an explicitly admitted compound physical source.
pub const COMPOUND_SCALAR_TARGET_SCOPE: &str = "Certified two-dimensional non-spinor Gamma compound physical source with both named scalar CIRs; complete decomposition into real/complex scalar type-A targets. No type-B/C target, non-Gamma/spinor source or setting search.";
/// Conjugate-pair target boundary for an explicitly admitted compound physical source.
pub const COMPOUND_CONJUGATE_PAIR_SCOPE: &str = "Certified two-dimensional non-spinor Gamma compound physical source with both named scalar CIRs; one complete nonreal conjugate scalar H pair: magnetic type C with anti, reducible complex unitary packet without anti. No type-B, real-H pair, non-Gamma/spinor source or setting search.";
pub(super) const COMPOUND_TYPE_C_SCOPE: &str = "Certified two-dimensional non-spinor Gamma compound physical source with both named scalar CIRs; restricted type-C target of two inequivalent real scalar H rows only. No type-B, non-Gamma/spinor source or setting search.";

/// Admit only a real 2D physical carrier of two distinct nonreal conjugate scalar
/// CIRs. The old verifier still refuses compounds. The returned model can be
/// passed explicitly to the existing target engines; this does not implement B.
pub fn verify_compound_matrix_gamma_input(
    request: &MatrixGammaRequest<'_>,
) -> Result<MatrixGammaModel, MatrixSubductionError> {
    verify_matrix_gamma_input_for(request, true)
}

pub(super) fn certify(
    irrep: &IrrepRecord,
) -> Result<(CertifiedTable, [CompoundConstituentCharacter; 2]), MatrixSubductionError> {
    let invalid = |reason: &str| MatrixSubductionError::CompoundSourceInvalid {
        sg: irrep.sg,
        label: irrep.ml,
        reason: reason.into(),
    };
    // Only the immutable authoritative constructor supplies named rows. Algebra
    // alone cannot distinguish a global conjugation/relabeling with the same trace.
    let view = irrep
        .compound_selected_arm_view()
        .map_err(|_| invalid("authoritative compound view unavailable"))?;
    let CompoundSelectedArmCharacter::DistinctComponentSum {
        first,
        second,
        block_trace,
    } = view
    else {
        return Err(invalid(
            "requires two distinct stored CIRs, not a generated conjugate",
        ));
    };
    let parts = [first, second];
    let metadata = irrep
        .compound_metadata()
        .ok_or_else(|| invalid("compound metadata missing"))?;
    if metadata.sg != irrep.sg
        || metadata.record_label != irrep.ml
        || metadata.cir_dimensions != [1, 1]
    {
        return Err(invalid(
            "owner metadata does not describe these scalar CIRs",
        ));
    }
    if parts[0].label == parts[1].label || parts[0].irnumber == parts[1].irnumber {
        return Err(invalid("CIR identities are not distinct"));
    }
    let table = certify_stored_matrices(irrep, SOURCE_DIMENSION, Some(&block_trace)).map_err(
        |problem| MatrixSubductionError::SourceMatrixUncertified {
            sg: irrep.sg,
            label: irrep.ml,
            problem,
        },
    )?;
    if real_commutant_dimension(&table.matrices, SOURCE_DIMENSION) != 2 {
        return Err(invalid(
            "physical source must have real commutant dimension two",
        ));
    }
    for (i, part) in parts.iter().enumerate() {
        if part.component != i
            || part.dimension != 1
            || metadata.cir_labels[i] != part.label
            || metadata.cir_irnumbers[i] != part.irnumber
        {
            return Err(invalid(
                "scalar CIR identity differs from frozen owner metadata",
            ));
        }
        let scalar = scalar_targets::certify_scalar_characters(
            part.row.dimension(),
            part.row.operations(),
            part.row.values(),
        )
        .map_err(|problem| MatrixSubductionError::SourceMatrixUncertified {
            sg: irrep.sg,
            label: part.label,
            problem,
        })?;
        if scalar.operations != table.operations {
            return Err(invalid(
                "complete CIR/PIR operation order or translation differs",
            ));
        }
    }
    if !parts[0].row.values().iter().any(|z| !close(z.im, 0.0)) {
        return Err(invalid("real scalar rows are not this compound family"));
    }
    for (i, (a, b)) in parts[0]
        .row
        .values()
        .iter()
        .zip(parts[1].row.values())
        .enumerate()
    {
        if !close(a.re, b.re)
            || !close(a.im, -b.im)
            || !close(a.re + b.re, matrix_trace(&table.matrices[i], 2))
        {
            return Err(invalid(
                "CIRs are not conjugate or do not give the physical matrix trace",
            ));
        }
    }
    Ok((table, parts))
}

#[cfg(test)]
mod tests {
    use super::scalar_targets::certify_scalar_characters;
    use super::*;

    #[test]
    fn all_69_compounds_certify_both_cirs_without_relaxing_the_ordinary_gate() {
        let mut count = 0;
        for sg in 1..=230 {
            for ir in query::irreps_of(sg).iter().filter(|r| {
                !r.spinor
                    && r.dim == 2
                    && r.compound_metadata().is_some()
                    && r.k_vector().numerators == [0; 3]
                    && r.k_vector().denominator == 1
            }) {
                let (table, parts) = certify(ir).unwrap();
                assert_eq!(parts.len(), 2);
                assert_eq!(real_commutant_dimension(&table.matrices, 2), 2);
                assert!(matches!(
                    certify_irrep_matrices(ir, 2),
                    Err(MatrixCertificationProblem::Reducible {
                        commutant_dimension: 2
                    })
                ));
                assert!(matches!(
                    resolve_matrix_source(SourceIdentity {
                        parent_sg: sg,
                        label: ir.ml,
                        convention: LabelConvention::Cdml
                    }),
                    Err(MatrixSubductionError::SourceIsCompound { .. })
                ));
                count += 1;
            }
        }
        assert_eq!(count, 69);
    }

    #[test]
    fn unchanged_total_trace_does_not_certify_either_constituent_law() {
        let ir = query::irreps_of(75)
            .iter()
            .find(|r| r.ml == "GM3GM4")
            .unwrap();
        let (_, parts) = certify(ir).unwrap();
        let i = parts[0]
            .row
            .values()
            .iter()
            .position(|z| z.im.abs() > 0.1)
            .unwrap();
        let mut changed = [
            parts[0].row.values().to_vec(),
            parts[1].row.values().to_vec(),
        ];
        for p in 0..2 {
            changed[p][i] = changed[p][i].conj();
            let row = &parts[p].row;
            assert!(matches!(
                certify_scalar_characters(1, row.operations(), &changed[p]),
                Err(MatrixCertificationProblem::LawViolated { .. })
            ));
            // Global conjugation is a valid representation: labels come from
            // authoritative metadata, not from this algebraic check.
            assert!(
                certify_scalar_characters(
                    1,
                    row.operations(),
                    &row.values().iter().map(|z| z.conj()).collect::<Vec<_>>()
                )
                .is_ok()
            );
        }
        for (n, _) in parts[0].row.values().iter().enumerate() {
            assert_eq!(
                changed[0][n] + changed[1][n],
                parts[0].row.values()[n] + parts[1].row.values()[n]
            );
        }
    }

    #[test]
    fn scalar_certificate_keeps_invalid_shapes_values_and_arithmetic_errors() {
        let ir = query::irreps_of(75)
            .iter()
            .find(|r| r.ml == "GM3GM4")
            .unwrap();
        let (_, parts) = certify(ir).unwrap();
        let row = &parts[0].row;
        assert!(matches!(
            certify_scalar_characters(2, row.operations(), row.values()),
            Err(MatrixCertificationProblem::TypedRowShapeMismatch { .. })
        ));
        assert!(certify_scalar_characters(1, row.operations(), &row.values()[1..]).is_err());
        for value in [
            Complex64::new(0.5, 0.0),
            Complex64::new(f64::NAN, 0.0),
            Complex64::new(f64::INFINITY, 0.0),
        ] {
            let mut values = row.values().to_vec();
            values[1] = value;
            assert!(certify_scalar_characters(1, row.operations(), &values).is_err());
        }
        let mut ops = row.operations().to_vec();
        ops[0].translation[0] = 1.0 / 101.0;
        assert!(matches!(
            certify_scalar_characters(1, &ops, row.values()),
            Err(MatrixCertificationProblem::TranslationOffGrid { .. })
        ));
        ops[0].translation[0] = i128::MAX as f64;
        assert!(matches!(
            certify_scalar_characters(1, &ops, row.values()),
            Err(MatrixCertificationProblem::Arithmetic { .. })
        ));
    }
}
