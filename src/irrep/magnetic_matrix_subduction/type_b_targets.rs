//! Narrow type B: a repeated real scalar H character and a quaternionic anti action.
use super::conjugate_pairs::realified_scalar_table;
use super::scalar_targets::{bind_h_row, resolve_h_row};
use super::*;
use crate::irrep::magnetic_summary::SourceIrrepSummary;

/// Explicit B target boundary for the original noncompound matrix source.
pub const TYPE_B_SCOPE: &str = "Certified real non-spinor non-compound two-dimensional Gamma source, completely reconstructed as a two-dimensional type-B corep from one authoritative real scalar H character. T(h)=chi(h) I, T(a0 h)=chi(h) J K, J^2=-I. No complex or higher-dimensional H, spinor, non-Gamma, full-star or setting search.";
/// The same B target boundary with explicitly certified compound source provenance.
pub const COMPOUND_TYPE_B_SCOPE: &str = "Explicitly certified non-spinor two-dimensional compound Gamma physical source with two named conjugate scalar CIRs, completely reconstructed as a two-dimensional type-B corep from one authoritative real scalar H character. T(h)=chi(h) I, T(a0 h)=chi(h) J K, J^2=-I. No complex or higher-dimensional H, spinor, non-Gamma, full-star or setting search.";

/// One complete type-B target, not two independent scalar target copies.
#[derive(Debug)]
pub struct MatrixTypeBSubduction<'a> {
    model: &'a MatrixGammaModel,
    piece: MatrixSubductionPiece,
    h_source: SourceIrrepSummary,
    storage: HCharacterStorage,
}
impl MatrixTypeBSubduction<'_> {
    /// The original certified source, geometry, eta and source-frame evidence.
    pub fn model(&self) -> &MatrixGammaModel {
        self.model
    }
    /// This explicit target boundary, distinct from the legacy model scope.
    pub fn scope(&self) -> &'static str {
        self.model.source_scope(TYPE_B_SCOPE, COMPOUND_TYPE_B_SCOPE)
    }
    /// Reuses the real matrix target payload; its dimension is two, not H's one.
    pub fn piece(&self) -> &MatrixSubductionPiece {
        &self.piece
    }
    /// The single genuine H identity, repeated twice in the restriction to H.
    pub fn h_source(&self) -> &SourceIrrepSummary {
        &self.h_source
    }
    /// The authoritative H row, never a fabricated second constituent.
    pub fn storage(&self) -> &HCharacterStorage {
        &self.storage
    }
    /// Source dimension, two.
    pub fn source_dimension(&self) -> u32 {
        self.model.source_dimension()
    }
    /// Actual exhibited complete-carrier coverage.
    pub fn covered_dimension(&self) -> u32 {
        self.piece.covered_dimension()
    }
}
fn failure(context: &TargetContext, label: &str, reason: &str) -> MatrixSubductionError {
    MatrixSubductionError::TargetCorepInconsistent {
        uni: context.uni,
        label: label.into(),
        reason: reason.into(),
    }
}

/// Values have already passed full authoritative H-row certification and binding.
/// Check all summary data before excluding nonreal H from this narrow family.
fn type_b_table(
    corep: &MagneticCorepSummary,
    context: &TargetContext,
    values: &[Complex64],
    anti: usize,
) -> Result<Option<TargetTable>, MatrixSubductionError> {
    let fail = |s| failure(context, &corep.label, s);
    if values.len() != context.classes.len() {
        return Err(fail("incomplete type-B H values"));
    }
    let mut real = true;
    for (i, _) in context
        .classes
        .iter()
        .enumerate()
        .filter(|(_, g)| !g.time_reversal)
    {
        let value = values[i];
        if !value.re.is_finite() || !value.im.is_finite() || !close(value.norm_sqr(), 1.0) {
            return Err(fail("invalid type-B scalar H value"));
        }
        if context.class_columns[i].is_empty() {
            return Err(fail("uncovered type-B unitary class"));
        }
        for &column in &context.class_columns[i] {
            let summary = corep
                .characters
                .get(column)
                .copied()
                .flatten()
                .ok_or_else(|| fail("missing type-B summary character"))?;
            if !close(summary.re, 2.0 * value.re) || !close(summary.im, 2.0 * value.im) {
                return Err(fail(
                    "type-B summary does not equal twice the authoritative H row",
                ));
            }
        }
        real &= close(value.im, 0.0) && close(value.re.abs(), 1.0);
    }
    if !real {
        return Ok(None);
    }
    // Normalize only the vanishing imaginary part, including -0, to fix sqrt(-1)=+i.
    let real_values: Vec<_> = values.iter().map(|z| Complex64::new(z.re, 0.0)).collect();
    let table = realified_scalar_table(corep, context, &real_values, Some(anti))?;
    if !table.operators[anti]
        .entries
        .iter()
        .zip([0.0, -1.0, 1.0, 0.0])
        .all(|(&a, b)| close(a, b))
    {
        return Err(fail(
            "real scalar type B requires chi(a0 squared)=-1 and T(a0)=J K",
        ));
    }
    Ok(Some(table))
}

/// Explicitly reconstruct the complete real 2D source as one narrow type-B corep.
/// Existing A/C/scalar/pair routes do not dispatch here. H norm four and the Hom
/// real kernel dimension are not multiplicities; copies must cover the source.
pub fn subduce_matrix_gamma_type_b(
    model: &MatrixGammaModel,
) -> Result<MatrixTypeBSubduction<'_>, MatrixSubductionError> {
    let (gamma, context) = prepare_matrix_targets(model)?;
    let anti = context.classes.iter().position(|g| g.time_reversal);
    let mut found = None;
    let mut examined = Vec::new();
    for corep in &gamma.coreps {
        if corep.corep_type != CorepType::B
            || corep.dim != 2
            || corep.source_irreps.iter().any(|s| s.spinor)
        {
            continue;
        }
        if !has_unitary_overlap(model, corep, &context)? {
            continue;
        }
        let h_source = &corep.source_irreps[0];
        let h = resolve_h_row(h_source, model.uni)?;
        let values = bind_h_row(&h, h_source, &corep.label, &context)?;
        let a = anti.ok_or_else(|| {
            failure(
                &context,
                &corep.label,
                "type B requires an antiunitary coset",
            )
        })?;
        let Some(table) = type_b_table(corep, &context, &values, a)? else {
            continue;
        };
        examined.push(corep.label.clone());
        let (_, copies) = matching_copies(&RealCorep::from_source(&model.actions), &table);
        if copies.is_empty() {
            continue;
        }
        let piece = MatrixSubductionPiece {
            label: corep.label.clone(),
            source_label: h_source.ml,
            unitary_sg: h_source.sg,
            dimension: 2,
            multiplicity: copies.len() as u32,
            copies,
            operators: table.operators,
        };
        if !verify_full_isometry(std::slice::from_ref(&piece)) {
            return Err(MatrixSubductionError::SourceSpaceNotIsometric {
                covered_dimension: piece.covered_dimension(),
            });
        }
        verify_reconstruction_data(model, &corep.label, 2, piece.copies(), piece.operators())?;
        if found.is_some() {
            return Err(failure(
                &context,
                &corep.label,
                "more than one complete type-B target",
            ));
        }
        found = Some(MatrixTypeBSubduction {
            model,
            piece,
            h_source: h_source.clone(),
            storage: h.storage,
        });
    }
    found.ok_or(MatrixSubductionError::NoMatchingTarget {
        parent_sg: model.parent_sg,
        source: model.source_label,
        uni: model.uni,
        examined,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn model(sg: u8, label: &str, eta: RealSign) -> MatrixGammaModel {
        verify_compound_matrix_gamma_input(
            &MatrixGammaRequest::new(
                SourceIdentity {
                    parent_sg: sg,
                    label,
                    convention: LabelConvention::Cdml,
                },
                RecordIdentity {
                    parent_sg: sg,
                    condensing_label: if sg < 83 { "GM2" } else { "GM2+" },
                    convention: LabelConvention::Cdml,
                    direction_index: 0,
                },
                Some(eta),
            )
            .with_embedding(EmbeddingEvidence {
                map: crate::irrep::subduction::Mat3R::identity(),
                shift: [Rat::ZERO; 3],
                lattice_policy: LatticePolicy::RequireExact,
            }),
        )
        .unwrap()
    }
    #[test]
    fn all_twenty_targets_and_every_anti_representative_have_one_copy_not_four() {
        let (mut models, mut representatives) = (0, 0);
        for sg in 75..=88 {
            let labels: &[&str] = if sg < 83 {
                &["GM3GM4"]
            } else {
                &["GM3+GM4+", "GM3-GM4-"]
            };
            for label in labels {
                for eta in [RealSign::Plus, RealSign::Minus] {
                    let m = model(sg, label, eta);
                    let (gamma, ctx) = prepare_matrix_targets(&m).unwrap();
                    let c = gamma
                        .coreps
                        .iter()
                        .find(|c| {
                            c.corep_type == CorepType::B
                                && has_unitary_overlap(&m, c, &ctx).unwrap()
                        })
                        .unwrap();
                    let h = resolve_h_row(&c.source_irreps[0], m.uni).unwrap();
                    let v = bind_h_row(&h, &c.source_irreps[0], &c.label, &ctx).unwrap();
                    for (a, _) in ctx
                        .classes
                        .iter()
                        .enumerate()
                        .filter(|(_, g)| g.time_reversal)
                    {
                        let table = type_b_table(c, &ctx, &v, a).unwrap().unwrap();
                        assert_eq!(table.operators[a].entries, vec![0., -1., 1., 0.]);
                        let (kernel, copies) =
                            matching_copies(&RealCorep::from_source(&m.actions), &table);
                        assert_eq!(kernel.len(), 4);
                        assert_eq!(copies.len(), 1);
                        let u = &copies[0];
                        assert!(verify_column_isometry(&[
                            vec![u[0], u[2]],
                            vec![u[1], u[3]]
                        ]));
                        verify_reconstruction_data(&m, &c.label, 2, &copies, &table.operators)
                            .unwrap();
                        representatives += 1;
                    }
                    let out = subduce_matrix_gamma_type_b(&m).unwrap();
                    assert_eq!(out.piece().multiplicity(), 1);
                    assert_eq!(out.covered_dimension(), 2);
                    assert_eq!(out.piece().label(), c.label);
                    models += 1;
                }
            }
        }
        assert_eq!((models, representatives), (40, 128));
    }
    #[test]
    fn negative_squares_do_not_replace_twisted_or_coset_laws_and_complex_h_stays_outside() {
        let m = model(75, "GM3GM4", RealSign::Plus);
        let (g, mut ctx) = prepare_matrix_targets(&m).unwrap();
        let mut c = g
            .coreps
            .iter()
            .find(|c| c.corep_type == CorepType::B)
            .unwrap()
            .clone();
        let h = resolve_h_row(&c.source_irreps[0], m.uni).unwrap();
        let v = bind_h_row(&h, &c.source_irreps[0], &c.label, &ctx).unwrap();
        let a = ctx.classes.iter().position(|g| g.time_reversal).unwrap();
        let mut table = type_b_table(&c, &ctx, &v, a).unwrap().unwrap();
        for value in &mut table.operators[a].entries {
            *value = -*value;
        }
        // Every anti still squares to -I, but one local sign violates the coset law.
        for op in table.operators.iter().filter(|op| op.antiunitary) {
            let x = &op.entries;
            let square = [
                x[0] * x[0] + x[1] * x[2],
                x[0] * x[1] + x[1] * x[3],
                x[2] * x[0] + x[3] * x[2],
                x[2] * x[1] + x[3] * x[3],
            ];
            assert!(
                square
                    .iter()
                    .zip([-1., 0., 0., -1.])
                    .all(|(a, b)| close(*a, b))
            );
        }
        assert!(check_target_law(&table, &ctx).is_err());
        // D8: H=<r²,s>, a0=r, chi(r²)=-1, chi(s)=1. H law and a0² pass,
        // but r s r^-1 = r²s changes chi. This is an algebraic control, not H provenance.
        let id = MagneticOperation {
            rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            translation: [Rat::ZERO; 3],
            time_reversal: false,
        };
        let r = MagneticOperation {
            rotation: [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
            translation: [Rat::ZERO; 3],
            time_reversal: true,
        };
        let s = MagneticOperation {
            rotation: [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
            translation: [Rat::ZERO; 3],
            time_reversal: false,
        };
        let mut power = id;
        let mut classes = Vec::new();
        let mut values = Vec::new();
        for i in 0..4 {
            classes.push(compose(&power, &id).unwrap());
            classes.push(compose(&power, &s).unwrap());
            let value = if i % 2 == 1 {
                0.
            } else if i == 0 {
                1.
            } else {
                -1.
            };
            values.extend([Complex64::new(value, 0.); 2]);
            power = compose(&power, &r).unwrap();
        }
        ctx.classes = classes;
        ctx.class_columns = (0..8).map(|i| vec![i]).collect();
        c.characters = ctx
            .classes
            .iter()
            .zip(&values)
            .map(|(g, v)| if g.time_reversal { None } else { Some(2. * v) })
            .collect();
        let index = |g: &MagneticOperation| ctx.classes.iter().position(|x| x == g).unwrap();
        for (i, g) in ctx
            .classes
            .iter()
            .enumerate()
            .filter(|(_, g)| !g.time_reversal)
        {
            for (j, h) in ctx
                .classes
                .iter()
                .enumerate()
                .filter(|(_, h)| !h.time_reversal)
            {
                assert!(
                    (values[i] * values[j] - values[index(&compose(g, h).unwrap())]).norm() < 1e-10
                );
            }
        }
        assert_eq!(
            values[index(&compose(&r, &r).unwrap())],
            Complex64::new(-1., 0.)
        );
        assert!(matches!(
            type_b_table(&c, &ctx, &values, 2),
            Err(MatrixSubductionError::TargetCorepInconsistent { .. })
        ));
        // A genuinely certified nonreal CIR H row is outside this explicit real-H family.
        let m = verify_matrix_gamma_input(
            &MatrixGammaRequest::new(
                SourceIdentity {
                    parent_sg: 149,
                    label: "GM3",
                    convention: LabelConvention::Cdml,
                },
                RecordIdentity {
                    parent_sg: 149,
                    condensing_label: "H1",
                    convention: LabelConvention::Cdml,
                    direction_index: 2,
                },
                Some(RealSign::Plus),
            )
            .with_embedding(EmbeddingEvidence {
                map: crate::irrep::subduction::Mat3R::from_ints([
                    [-2, 1, 0],
                    [-1, 2, 0],
                    [0, 0, -2],
                ]),
                shift: [Rat::ZERO; 3],
                lattice_policy: LatticePolicy::RequireExact,
            }),
        )
        .unwrap();
        let (g, ctx) = prepare_matrix_targets(&m).unwrap();
        let mut c = g
            .coreps
            .iter()
            .find(|c| {
                c.corep_type == CorepType::C
                    && c.dim == 2
                    && c.source_irreps.len() == 2
                    && c.source_irreps.iter().all(|s| !s.spinor && s.dim == 1)
            })
            .unwrap()
            .clone();
        let h = resolve_h_row(&c.source_irreps[0], m.uni).unwrap();
        let values = bind_h_row(&h, &c.source_irreps[0], &c.label, &ctx).unwrap();
        assert!(values.iter().any(|v| v.im.abs() > 0.1));
        c.corep_type = CorepType::B;
        c.source_irreps.truncate(1);
        for (i, columns) in ctx
            .class_columns
            .iter()
            .enumerate()
            .filter(|(i, _)| !ctx.classes[*i].time_reversal)
        {
            for &col in columns {
                c.characters[col] = Some(2. * values[i]);
            }
        }
        let a = ctx.classes.iter().position(|g| g.time_reversal).unwrap();
        assert!(type_b_table(&c, &ctx, &values, a).unwrap().is_none());
        let unit = ctx.classes.iter().position(|g| !g.time_reversal).unwrap();
        c.characters[ctx.class_columns[unit][0]] = None;
        assert!(matches!(
            type_b_table(&c, &ctx, &values, a),
            Err(MatrixSubductionError::TargetCorepInconsistent { .. })
        ));
    }
    #[test]
    fn b_summary_norm_identity_missing_and_nonfinite_data_are_not_zero_overlap() {
        let m = model(75, "GM3GM4", RealSign::Plus);
        let (g, mut ctx) = prepare_matrix_targets(&m).unwrap();
        let c = g
            .coreps
            .iter()
            .find(|c| c.corep_type == CorepType::B)
            .unwrap()
            .clone();
        let h = resolve_h_row(&c.source_irreps[0], m.uni).unwrap();
        let v = bind_h_row(&h, &c.source_irreps[0], &c.label, &ctx).unwrap();
        let a = ctx.classes.iter().position(|g| g.time_reversal).unwrap();
        assert!(has_unitary_overlap(&m, &c, &ctx).unwrap());
        let i = ctx.classes.iter().position(|g| !g.time_reversal).unwrap();
        let col = ctx.class_columns[i][0];
        let is_data_error = |r: Result<bool, MatrixSubductionError>| {
            assert!(
                matches!(r,Err(MatrixSubductionError::TargetCorepInconsistent{label,..}) if label==c.label)
            )
        };
        for bad in [
            None,
            Some(Complex64::new(f64::NAN, 0.)),
            Some(Complex64::new(0., f64::INFINITY)),
            Some(Complex64::new(1., 0.)),
        ] {
            let mut bad_c = c.clone();
            bad_c.characters[col] = bad;
            is_data_error(has_unitary_overlap(&m, &bad_c, &ctx));
            assert!(
                matches!(type_b_table(&bad_c,&ctx,&v,a),Err(MatrixSubductionError::TargetCorepInconsistent{label,..}) if label==c.label)
            );
        }
        let mut bad_c = c.clone();
        bad_c.source_irreps.clear();
        is_data_error(has_unitary_overlap(&m, &bad_c, &ctx));
        let mut bad_c = c.clone();
        bad_c.source_irreps[0].sg = 1;
        is_data_error(has_unitary_overlap(&m, &bad_c, &ctx));
        let mut bad_c = c.clone();
        let identity = find_identity_class(&ctx.classes).unwrap();
        for (i, columns) in ctx.class_columns.iter().enumerate() {
            if !ctx.classes[i].time_reversal && i != identity {
                for &col in columns {
                    bad_c.characters[col] = Some(Complex64::new(0., 0.));
                }
            }
        }
        is_data_error(has_unitary_overlap(&m, &bad_c, &ctx));
        for z in [
            Complex64::new(0.5, 0.),
            Complex64::new(f64::NAN, 0.),
            Complex64::new(1., f64::INFINITY),
            Complex64::new(0., 1.),
        ] {
            let mut bad = v.clone();
            bad[i] = z;
            assert!(
                matches!(type_b_table(&c,&ctx,&bad,a),Err(MatrixSubductionError::TargetCorepInconsistent{label,..}) if label==c.label)
            );
        }
        assert!(type_b_table(&c, &ctx, &v[..1], a).is_err());
        ctx.class_columns[i].clear();
        is_data_error(has_unitary_overlap(&m, &c, &ctx));
        assert!(type_b_table(&c, &ctx, &v, a).is_err());
    }
    #[test]
    fn positive_square_is_not_b_and_bad_representatives_coefficients_and_arithmetic_fail() {
        let m = model(75, "GM3GM4", RealSign::Minus);
        let (g, mut ctx) = prepare_matrix_targets(&m).unwrap();
        let c = g
            .coreps
            .iter()
            .find(|c| c.corep_type == CorepType::B)
            .unwrap()
            .clone();
        let h = resolve_h_row(&c.source_irreps[0], m.uni).unwrap();
        let v = bind_h_row(&h, &c.source_irreps[0], &c.label, &ctx).unwrap();
        let a = ctx.classes.iter().position(|g| g.time_reversal).unwrap();
        let unit = ctx.classes.iter().position(|g| !g.time_reversal).unwrap();
        assert!(type_b_table(&c, &ctx, &v, unit).is_err());
        assert!(type_b_table(&c, &ctx, &v, ctx.classes.len()).is_err());
        let mut positive = c.clone();
        let mut trivial = v.clone();
        for (i, op) in ctx
            .classes
            .iter()
            .enumerate()
            .filter(|(_, op)| !op.time_reversal)
        {
            let _ = op;
            trivial[i] = Complex64::new(1., 0.);
            for &col in &ctx.class_columns[i] {
                positive.characters[col] = Some(Complex64::new(2., 0.));
            }
        }
        // Two A copies obey every product law; that does not turn them into B.
        let two_a = realified_scalar_table(&positive, &ctx, &trivial, Some(a)).unwrap();
        assert_eq!(two_a.operators[a].entries, vec![1., -0., 0., 1.]);
        assert!(
            matches!(type_b_table(&positive,&ctx,&trivial,a),Err(MatrixSubductionError::TargetCorepInconsistent{label,..}) if label==c.label)
        );
        let mut signed_zero = v.clone();
        for z in &mut signed_zero {
            z.im = -0.0;
        }
        let mut table = type_b_table(&c, &ctx, &signed_zero, a).unwrap().unwrap();
        assert_eq!(table.operators[a].entries, vec![0., -1., 1., 0.]);
        table.operators[a].entries[0] += 0.25;
        assert!(check_target_law(&table, &ctx).is_err());
        // An arbitrary scalar anti phase cannot fix the square: (i K)^2=+1.
        let z = Complex64::new(0.3, -0.7);
        let ik = |z: Complex64| Complex64::new(0., 1.) * z.conj();
        assert!((ik(ik(z)) - z).norm() < 1e-10);
        assert!((ik(ik(z)) + z).norm() > 1.);
        ctx.classes[a].translation[2] = Rat::new(i128::MAX, 1).unwrap();
        assert!(matches!(
            type_b_table(&c, &ctx, &v, a),
            Err(MatrixSubductionError::Arithmetic { .. })
        ));
    }
}
