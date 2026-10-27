//! Realification targets for authoritative conjugate scalar H pairs.
//! The purely unitary result is a packet, not a complex irreducible corep.
use super::scalar_targets::{bind_h_row, resolve_h_row};
use super::*;
use crate::irrep::magnetic_summary::SourceIrrepSummary;

/// Exact scope of the conjugate-pair companion.
pub const CONJUGATE_PAIR_SCOPE: &str = "Certified real non-spinor non-compound two-dimensional Gamma source, with one complete authoritative nonreal conjugate scalar H pair. Magnetic type-C corep when an antiunitary coset exists; an explicitly reducible complex ordinary pair otherwise. No type-B, real H pairs, higher-dimensional H sources, spinor, non-Gamma, full-star or setting search.";

/// Whether the two-dimensional carrier is a magnetic corep or an ordinary packet.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConjugatePairKind {
    /// Irreducible magnetic type C, with a genuinely antiunitary coset.
    MagneticTypeC,
    /// Two inequivalent complex ordinary irreps; NOT a complex irreducible 2D irrep.
    UnitaryPair,
}

/// A fully reconstructed target in the canonical realification basis.
#[derive(Debug, Clone)]
pub struct ConjugatePairPiece {
    label: String,
    kind: ConjugatePairKind,
    h_sources: [SourceIrrepSummary; 2],
    storage: [HCharacterStorage; 2],
    copies: Vec<Vec<Complex64>>,
    operators: Vec<TargetMatrixOperator>,
}
impl ConjugatePairPiece {
    /// Original summary label for the pair, never a fabricated standalone label.
    pub fn label(&self) -> &str {
        &self.label
    }
    /// Distinguishes magnetic irreducibility from a purely unitary packet.
    pub fn kind(&self) -> ConjugatePairKind {
        self.kind
    }
    /// Both H constituents; use labels() for genuine labels, not legacy raw bc.
    pub fn h_sources(&self) -> &[SourceIrrepSummary; 2] {
        &self.h_sources
    }
    /// Authoritative storage identities in the same order as h_sources().
    pub fn storage(&self) -> &[HCharacterStorage; 2] {
        &self.storage
    }
    /// The full pair's carrier dimension, not a constituent's dimension.
    pub fn dimension(&self) -> u32 {
        2
    }
    /// Number of complete pair copies (one for the accepted 2D source).
    pub fn multiplicity(&self) -> u32 {
        self.copies.len() as u32
    }
    /// Row-major 2 by 2 complex intertwiners into the original source basis.
    pub fn copies(&self) -> &[Vec<Complex64>] {
        &self.copies
    }
    /// Real 2 by 2 coefficients with the original antiunitary flags.
    pub fn operators(&self) -> &[TargetMatrixOperator] {
        &self.operators
    }
}

/// Complete decomposition of the fixed model into one conjugate pair.
#[derive(Debug)]
pub struct ConjugatePairSubduction<'a> {
    model: &'a MatrixGammaModel,
    piece: ConjugatePairPiece,
}
impl ConjugatePairSubduction<'_> {
    /// The original verified source and embedding.
    pub fn model(&self) -> &MatrixGammaModel {
        self.model
    }
    /// This entry's explicit target boundary, not model.scope().
    pub fn scope(&self) -> &'static str {
        CONJUGATE_PAIR_SCOPE
    }
    /// The full pair with both authoritative H identities.
    pub fn piece(&self) -> &ConjugatePairPiece {
        &self.piece
    }
    /// Certified source dimension.
    pub fn source_dimension(&self) -> u32 {
        self.model.source_dimension()
    }
    /// Actual exhibited carrier coverage.
    pub fn covered_dimension(&self) -> u32 {
        self.piece.dimension() * self.piece.multiplicity()
    }
}

fn failure(context: &TargetContext, label: &str, reason: &str) -> MatrixSubductionError {
    MatrixSubductionError::TargetCorepInconsistent {
        uni: context.uni,
        label: label.into(),
        reason: reason.into(),
    }
}
fn realification(z: Complex64) -> Vec<f64> {
    vec![z.re, -z.im, z.im, z.re]
}

/// A nonreal character and its conjugate have the real carrier R(χ). In the
/// magnetic C case χ is invariant under a0 conjugation, so R(sqrt(χ(a0²))) K
/// supplies a0. The complete law below proves this condition instead of assuming it.
fn pair_table(
    corep: &MagneticCorepSummary,
    context: &TargetContext,
    values: &[Complex64],
    anti: Option<usize>,
) -> Result<TargetTable, MatrixSubductionError> {
    let fail = |s| failure(context, &corep.label, s);
    if values.len() != context.classes.len() || anti.is_some_and(|a| a >= context.classes.len()) {
        return Err(fail("invalid pair table dimensions or representative"));
    }
    let product = |left: usize, right: usize| -> Result<usize, MatrixSubductionError> {
        let op = compose(&context.classes[left], &context.classes[right]).map_err(|e| {
            MatrixSubductionError::Arithmetic {
                reason: e.to_string(),
            }
        })?;
        find_model_class(&op, &context.classes, &context.lattice)
            .ok_or_else(|| fail("a product is not a magnetic class"))
    };
    let mut coefficients = values.to_vec();
    if let Some(a) = anti {
        if !context.classes[a].time_reversal {
            return Err(fail("the coset representative is not antiunitary"));
        }
        let identity =
            find_identity_class(&context.classes).ok_or_else(|| fail("no identity class"))?;
        let mut inverse = None;
        for j in 0..context.classes.len() {
            if product(j, a)? == identity {
                inverse = Some(j);
                break;
            }
        }
        let inverse =
            inverse.ok_or_else(|| fail("no inverse of the antiunitary representative"))?;
        let square = product(a, a)?;
        if context.classes[square].time_reversal {
            return Err(fail("the antiunitary square is not unitary"));
        }
        let lambda = values[square].sqrt();
        if !lambda.re.is_finite()
            || !lambda.im.is_finite()
            || !close(lambda.norm_sqr(), 1.0)
            || (lambda * lambda - values[square]).norm() > GRID_TOLERANCE
        {
            return Err(fail("no unit-modulus square root"));
        }
        for (i, g) in context.classes.iter().enumerate() {
            if g.time_reversal {
                let h = product(inverse, i)?;
                if context.classes[h].time_reversal {
                    return Err(fail("antiunitary coset factor is not unitary"));
                }
                // R(χ(h)) is real: K does NOT conjugate its scalar parameter χ(h).
                coefficients[i] = lambda * values[h];
            }
        }
    }
    if coefficients.len() != context.classes.len()
        || coefficients
            .iter()
            .any(|z| !z.re.is_finite() || !z.im.is_finite() || !close(z.norm_sqr(), 1.0))
    {
        return Err(fail("incomplete or nonunitary pair coefficients"));
    }
    let table = TargetTable {
        dimension: 2,
        operators: coefficients
            .into_iter()
            .zip(&context.classes)
            .map(|(z, g)| TargetMatrixOperator {
                dimension: 2,
                entries: realification(z),
                antiunitary: g.time_reversal,
            })
            .collect(),
    };
    check_target_law(&table, context).map_err(|e| label_error(e, context.uni, &corep.label))?;
    Ok(table)
}

fn is_conjugate_pair(
    corep: &MagneticCorepSummary,
    context: &TargetContext,
    values: &[Vec<Complex64>; 2],
) -> Result<bool, MatrixSubductionError> {
    if values.iter().any(|v| v.len() != context.classes.len()) {
        return Err(failure(context, &corep.label, "incomplete paired H values"));
    }
    let unitary: Vec<_> = (0..context.classes.len())
        .filter(|&i| !context.classes[i].time_reversal)
        .collect();
    for &i in &unitary {
        if context.class_columns[i].is_empty() {
            return Err(failure(context, &corep.label, "uncovered unitary class"));
        }
        for &column in &context.class_columns[i] {
            let summary = corep
                .characters
                .get(column)
                .copied()
                .flatten()
                .ok_or_else(|| failure(context, &corep.label, "missing pair summary character"))?;
            let value = values[0][i] + values[1][i];
            if !close(summary.re, value.re) || !close(summary.im, value.im) {
                return Err(failure(
                    context,
                    &corep.label,
                    "typed H pair and summary disagree",
                ));
            }
        }
    }
    // Check data first: inconsistent H rows are never a scope exclusion.
    Ok(unitary
        .iter()
        .any(|&i| values[0][i].im.abs() > GRID_TOLERANCE)
        && unitary
            .iter()
            .all(|&i| (values[0][i].conj() - values[1][i]).norm() <= GRID_TOLERANCE))
}

/// Reconstruct the complete source as one authoritative nonreal conjugate H
/// pair. The unitary case is explicitly a reducible complex ordinary packet;
/// this does not widen any existing type-A or restricted real-H type-C entry.
pub fn subduce_matrix_gamma_conjugate_pair(
    model: &MatrixGammaModel,
) -> Result<ConjugatePairSubduction<'_>, MatrixSubductionError> {
    let (gamma, context) = prepare_matrix_targets(model)?;
    let anti = context.classes.iter().position(|g| g.time_reversal);
    let kind = if anti.is_some() {
        ConjugatePairKind::MagneticTypeC
    } else {
        ConjugatePairKind::UnitaryPair
    };
    let mut found = None;
    let mut examined = Vec::new();
    for corep in &gamma.coreps {
        let expected = if anti.is_some() {
            CorepType::C
        } else {
            CorepType::A
        };
        if corep.corep_type != expected
            || corep.dim != 2
            || corep.source_irreps.len() != 2
            || corep.source_irreps.iter().any(|s| s.spinor || s.dim != 1)
        {
            continue;
        }
        if !has_unitary_overlap(model, corep, &context)? {
            continue;
        }
        let sources = [
            corep.source_irreps[0].clone(),
            corep.source_irreps[1].clone(),
        ];
        let rows = [
            resolve_h_row(&sources[0], model.uni)?,
            resolve_h_row(&sources[1], model.uni)?,
        ];
        let values = [
            bind_h_row(&rows[0], &sources[0], &corep.label, &context)?,
            bind_h_row(&rows[1], &sources[1], &corep.label, &context)?,
        ];
        if !is_conjugate_pair(corep, &context, &values)? {
            continue;
        }
        examined.push(corep.label.clone());
        let table = pair_table(corep, &context, &values[0], anti)?;
        let (_, copies) = matching_copies(&RealCorep::from_source(&model.actions), &table);
        if copies.is_empty() {
            continue;
        }
        let columns: Vec<_> = copies
            .iter()
            .flat_map(|u| (0..2).map(move |j| vec![u[j], u[2 + j]]))
            .collect();
        if !verify_column_isometry(&columns) {
            return Err(MatrixSubductionError::SourceSpaceNotIsometric {
                covered_dimension: 2 * copies.len() as u32,
            });
        }
        verify_reconstruction_data(model, &corep.label, 2, &copies, &table.operators)?;
        if found.is_some() {
            return Err(failure(
                &context,
                &corep.label,
                "more than one complete conjugate-pair target",
            ));
        }
        found = Some(ConjugatePairPiece {
            label: corep.label.clone(),
            kind,
            h_sources: sources,
            storage: [rows[0].storage.clone(), rows[1].storage.clone()],
            copies,
            operators: table.operators,
        });
    }
    let piece = found.ok_or(MatrixSubductionError::NoMatchingTarget {
        parent_sg: model.parent_sg,
        source: model.source_label,
        uni: model.uni,
        examined,
    })?;
    Ok(ConjugatePairSubduction { model, piece })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::irrep::subduction::Mat3R;
    fn model(unitary: bool) -> MatrixGammaModel {
        verify_matrix_gamma_input(
            &MatrixGammaRequest::new(
                SourceIdentity {
                    parent_sg: 149,
                    label: "GM3",
                    convention: LabelConvention::Cdml,
                },
                RecordIdentity {
                    parent_sg: 149,
                    condensing_label: if unitary { "K1" } else { "H1" },
                    convention: LabelConvention::Cdml,
                    direction_index: 2,
                },
                Some(RealSign::Plus),
            )
            .with_embedding(EmbeddingEvidence {
                map: Mat3R::new(
                    [
                        [-2, 1, 0],
                        [-1, 2, 0],
                        [0, 0, if unitary { -1 } else { -2 }],
                    ]
                    .map(|r| r.map(|x| Rat::new(x, 1).unwrap())),
                ),
                shift: [Rat::ZERO; 3],
                lattice_policy: LatticePolicy::RequireExact,
            }),
        )
        .unwrap()
    }
    fn data(m: &MatrixGammaModel) -> (MagneticCorepSummary, TargetContext, Vec<Complex64>) {
        let (gamma, context) = prepare_matrix_targets(m).unwrap();
        let c = gamma
            .coreps
            .into_iter()
            .find(|c| c.dim == 2 && c.source_irreps.len() == 2 && !c.source_irreps[0].spinor)
            .unwrap();
        let h = resolve_h_row(&c.source_irreps[0], m.uni).unwrap();
        let v = bind_h_row(&h, &c.source_irreps[0], &c.label, &context).unwrap();
        (c, context, v)
    }
    #[test]
    fn unitary_packet_overlap_rejects_invalid_norm_missing_and_nonfinite_summary() {
        let m = model(true);
        let (c, mut context, values) = data(&m);
        assert!(matches!(c.corep_type, CorepType::A));
        assert_eq!(c.source_irreps.len(), 2);
        assert!(has_unitary_overlap(&m, &c, &context).unwrap());
        let reject = |bad: &MagneticCorepSummary, context: &TargetContext, expected: &str| {
            match has_unitary_overlap(&m, bad, context).unwrap_err() {
                MatrixSubductionError::TargetCorepInconsistent { uni, label, reason } => {
                    assert_eq!(uni, m.uni);
                    assert_eq!(label, c.label);
                    assert!(reason.contains(expected), "{reason}");
                }
                other => panic!("wrong error: {other:?}"),
            }
        };
        let mut single = c.clone();
        single.dim = 1;
        single.source_irreps.truncate(1);
        for (class, columns) in context.class_columns.iter().enumerate() {
            for &column in columns {
                single.characters[column] = Some(values[class]);
            }
        }
        assert!(
            has_unitary_overlap(&m, &single, &context).unwrap(),
            "single-source A keeps norm one"
        );
        let mut bad = c.clone();
        bad.characters.fill(Some(Complex64::new(2.0, 0.0)));
        reject(&bad, &context, "invalid unitary character norm");
        let identity = find_identity_class(&context.classes).unwrap();
        let other = (0..context.classes.len()).find(|&i| i != identity).unwrap();
        let column = context.class_columns[other][0];
        for value in [
            None,
            Some(Complex64::new(f64::NAN, 0.0)),
            Some(Complex64::new(0.0, f64::INFINITY)),
        ] {
            let mut bad = c.clone();
            bad.characters[column] = value;
            reject(
                &bad,
                &context,
                if value.is_none() {
                    "no summary character"
                } else {
                    "non-finite"
                },
            );
        }
        let mut bad = c.clone();
        bad.characters.clear();
        reject(&bad, &context, "no summary character");
        context.class_columns[other].clear();
        reject(&c, &context, "no columns");
    }

    #[test]
    fn every_anti_gauge_and_unitary_packet_reconstruct_without_counting_kernel_as_multiplicity() {
        let mut nonreal_squares = 0;
        let mut negative_squares = 0;
        let fourfold = verify_matrix_gamma_input(
            &MatrixGammaRequest::new(
                SourceIdentity {
                    parent_sg: 100,
                    label: "GM5",
                    convention: LabelConvention::Cdml,
                },
                RecordIdentity {
                    parent_sg: 100,
                    condensing_label: "A5",
                    convention: LabelConvention::Cdml,
                    direction_index: 0,
                },
                Some(RealSign::Plus),
            )
            .with_embedding(EmbeddingEvidence {
                map: Mat3R::new(
                    [[-1, -1, 0], [-1, 1, 0], [0, 0, -2]]
                        .map(|r| r.map(|x| Rat::new(x, 1).unwrap())),
                ),
                shift: [Rat::ZERO; 3],
                lattice_policy: LatticePolicy::RequireExact,
            }),
        )
        .unwrap();
        for (m, unitary) in [
            (model(false), false),
            (model(true), true),
            (fourfold, false),
        ] {
            let (c, context, v) = data(&m);
            let gauges: Vec<_> = if unitary {
                vec![None]
            } else {
                context
                    .classes
                    .iter()
                    .enumerate()
                    .filter(|(_, g)| g.time_reversal)
                    .map(|(i, _)| Some(i))
                    .collect()
            };
            for a in gauges {
                if let Some(a) = a {
                    let square = product_class(
                        &context.classes[a],
                        &context.classes[a],
                        &context.classes,
                        &context.lattice,
                    )
                    .unwrap();
                    if v[square].im.abs() > GRID_TOLERANCE {
                        nonreal_squares += 1;
                    }
                    if (v[square] + Complex64::new(1.0, 0.0)).norm() < GRID_TOLERANCE {
                        negative_squares += 1;
                    }
                }
                let table = pair_table(&c, &context, &v, a).unwrap();
                let (kernel, copies) = matching_copies(&RealCorep::from_source(&m.actions), &table);
                assert_eq!(kernel.len(), if unitary { 4 } else { 2 });
                assert_eq!(copies.len(), 1);
                verify_reconstruction_data(&m, &c.label, 2, &copies, &table.operators).unwrap();
                let u = &copies[0];
                assert!(verify_column_isometry(&[
                    vec![u[0], u[2]],
                    vec![u[1], u[3]]
                ]));
            }
        }
        assert_eq!(nonreal_squares, 2);
        assert_eq!(negative_squares, 2);
    }
    #[test]
    fn a_scalar_type_a_twist_cannot_be_relabelled_as_a_conjugate_c_pair() {
        let m = verify_matrix_gamma_input(
            &MatrixGammaRequest::new(
                SourceIdentity {
                    parent_sg: 89,
                    label: "GM5",
                    convention: LabelConvention::Cdml,
                },
                RecordIdentity {
                    parent_sg: 89,
                    condensing_label: "GM3",
                    convention: LabelConvention::Cdml,
                    direction_index: 0,
                },
                Some(RealSign::Plus),
            )
            .with_embedding(EmbeddingEvidence {
                map: Mat3R::identity(),
                shift: [Rat::ZERO; 3],
                lattice_policy: LatticePolicy::RequireExact,
            }),
        )
        .unwrap();
        let (gamma, context) = prepare_matrix_targets(&m).unwrap();
        let c = gamma
            .coreps
            .iter()
            .find(|c| c.dim == 1 && c.source_irreps[0].ml == "GM3")
            .unwrap();
        let h = resolve_h_row(&c.source_irreps[0], m.uni).unwrap();
        let v = bind_h_row(&h, &c.source_irreps[0], &c.label, &context).unwrap();
        let a = context
            .classes
            .iter()
            .position(|g| g.time_reversal)
            .unwrap();
        assert!(
            matches!(pair_table(c,&context,&v,Some(a)),Err(MatrixSubductionError::TargetCorepInconsistent{label,..}) if label==c.label)
        );
    }
    #[test]
    fn malformed_pair_data_and_target_product_overflow_remain_errors() {
        let m = model(false);
        let (c, mut context, v) = data(&m);
        let corrupt_real_rows = [
            vec![Complex64::new(1.0, 0.0); v.len()],
            vec![Complex64::new(1.0, 0.0); v.len()],
        ];
        assert!(
            is_conjugate_pair(&c, &context, &corrupt_real_rows).is_err(),
            "inconsistent H rows are not a scope exclusion"
        );
        let a = context
            .classes
            .iter()
            .position(|g| g.time_reversal)
            .unwrap();
        assert!(pair_table(&c, &context, &v[..1], Some(a)).is_err());
        assert!(pair_table(&c, &context, &v, Some(context.classes.len())).is_err());
        for bad in [0.5, f64::NAN, f64::INFINITY] {
            let mut v = v.clone();
            let i = context
                .classes
                .iter()
                .position(|g| !g.time_reversal)
                .unwrap();
            v[i] = Complex64::new(bad, 0.);
            assert!(
                matches!(pair_table(&c,&context,&v,Some(a)),Err(MatrixSubductionError::TargetCorepInconsistent{label,..}) if label==c.label)
            );
        }
        let table = pair_table(&c, &context, &v, Some(a)).unwrap();
        let mut tampered = table.clone();
        tampered.operators[a].entries[0] += 0.25;
        assert!(check_target_law(&tampered, &context).is_err());
        context.classes[0].rotation = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        context.classes[0].translation = [Rat::new(i128::MAX, 1).unwrap(), Rat::ZERO, Rat::ZERO];
        assert!(
            matches!(
                check_target_law(&table, &context),
                Err(MatrixSubductionError::Arithmetic { .. })
            ),
            "overflow is not evidence against the target law"
        );
    }
}
