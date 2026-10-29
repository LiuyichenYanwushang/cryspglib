//! Real/complex scalar type-A targets of the existing certified real 2D source.
//! Characters determine scalar matrices, not arbitrary higher-dimensional ones.
use super::*;
use crate::irrep::labels::IrrepLabels;
use crate::irrep::magnetic_subduction::TargetOperator;
use crate::irrep::magnetic_summary::SourceIrrepSummary;
use crate::irrep::types::CompoundSelectedArmCharacter;

/// Scope of the explicit scalar-target companion (not the legacy real-matrix scope).
pub const SCALAR_TARGET_SCOPE: &str = "Certified real non-spinor non-compound two-dimensional Gamma source; non-spinor one-dimensional type-A Gamma targets with real or complex H characters from standalone scalar rows or authoritative named typed CIR constituents. No type-B/C or higher-dimensional targets, unstored conjugates, non-Gamma/spinor/compound sources, or setting search.";

/// Authoritative storage of an H character; a CIR component is not a standalone irrep row.
#[derive(Debug, Clone)]
pub enum HCharacterStorage {
    /// An ordinary scalar record.
    Standalone {
        /// Genuine labels of the scalar storage record.
        record: IrrepLabels,
    },
    /// A named component exposed by the owning record's typed compound view.
    Constituent {
        /// Labels of the compound owner, not labels of the component.
        record: IrrepLabels,
        /// Zero-based constituent position in that owner's frozen metadata.
        component: usize,
        /// Authoritative CIR source identifier from the typed view.
        irnumber: u32,
    },
}

/// One scalar type-A target and its normalized independent 2×1 copies.
#[derive(Debug, Clone)]
pub struct ScalarTargetPiece {
    label: String,
    source: SourceIrrepSummary,
    storage: HCharacterStorage,
    copies: Vec<Vec<Complex64>>,
    operators: Vec<TargetOperator>,
}
impl ScalarTargetPiece {
    /// Magnetic summary label (not an invented BC label).
    pub fn label(&self) -> &str {
        &self.label
    }
    /// H identity; labels() preserves absent component BC mappings.
    pub fn source(&self) -> &SourceIrrepSummary {
        &self.source
    }
    /// The independently resolved owner and, when applicable, CIR identity.
    pub fn storage(&self) -> &HCharacterStorage {
        &self.storage
    }
    /// Target dimension, always one.
    pub fn dimension(&self) -> u32 {
        1
    }
    /// Number of independent copies, not the real dimension of Hom.
    pub fn multiplicity(&self) -> u32 {
        self.copies.len() as u32
    }
    /// Row-major 2×1 complex isometries.
    pub fn copies(&self) -> &[Vec<Complex64>] {
        &self.copies
    }
    /// Scalar unitary/antiunitary actions in model class order.
    pub fn operators(&self) -> &[TargetOperator] {
        &self.operators
    }
}

/// Complete scalar-target decomposition; the source/embedding remains the verified model's.
#[derive(Debug)]
pub struct ScalarTargetSubduction<'a> {
    model: &'a MatrixGammaModel,
    pieces: Vec<ScalarTargetPiece>,
}
impl ScalarTargetSubduction<'_> {
    /// The explicit target-family boundary; model.scope() describes the legacy entry.
    pub fn scope(&self) -> &'static str {
        self.model
            .source_scope(SCALAR_TARGET_SCOPE, COMPOUND_SCALAR_TARGET_SCOPE)
    }
    /// Original certified real 2D model.
    pub fn model(&self) -> &MatrixGammaModel {
        self.model
    }
    /// All exhibited scalar targets.
    pub fn pieces(&self) -> &[ScalarTargetPiece] {
        &self.pieces
    }
    /// Certified source dimension.
    pub fn source_dimension(&self) -> u32 {
        SOURCE_DIMENSION as u32
    }
    /// Sum of the exhibited copy counts, verified to be two.
    pub fn covered_dimension(&self) -> u32 {
        self.pieces
            .iter()
            .map(ScalarTargetPiece::multiplicity)
            .sum()
    }
}

fn data_error(
    context: &TargetContext,
    label: &str,
    reason: impl Into<String>,
) -> MatrixSubductionError {
    MatrixSubductionError::TargetCorepInconsistent {
        uni: context.uni,
        label: label.into(),
        reason: reason.into(),
    }
}

pub(super) struct HRow {
    row: CharacterRow,
    pub(super) storage: HCharacterStorage,
}

/// Resolve only metadata-named owners. Never suppress a malformed typed view,
/// choose the first ambiguous owner, or expose a generated conjugate as stored data.
pub(super) fn resolve_h_row(
    source: &SourceIrrepSummary,
    uni: usize,
) -> Result<HRow, MatrixSubductionError> {
    let error = |reason: String| MatrixSubductionError::TargetSourceUnresolved {
        uni,
        sg: source.sg,
        label: source.ml.into(),
        reason,
    };
    let mut found = Vec::new();
    for owner in query::irreps_of(source.sg) {
        let metadata = owner.compound_metadata();
        if owner.ml != source.ml && !metadata.is_some_and(|m| m.cir_labels.contains(&source.ml)) {
            continue;
        }
        if owner.spinor
            || source.spinor
            || source.dim != 1
            || owner.k_vector().numerators != [0; 3]
            || owner.k_vector().denominator != 1
        {
            return Err(error(
                "the named H owner is not a non-spinor Gamma scalar source".into(),
            ));
        }
        if metadata.is_some() {
            let view = owner
                .compound_selected_arm_view()
                .map_err(|e| error(format!("typed compound row: {e:?}")))?;
            let parts = match view {
                CompoundSelectedArmCharacter::ConjugateRealification { seed, .. } => vec![seed],
                CompoundSelectedArmCharacter::DistinctComponentSum { first, second, .. } => {
                    vec![first, second]
                }
            };
            for part in parts {
                if part.label == source.ml {
                    if part.dimension != 1 || part.row.dimension() != 1 {
                        return Err(error("the CIR dimension is not one".into()));
                    }
                    found.push(HRow {
                        row: part.row,
                        storage: HCharacterStorage::Constituent {
                            record: owner.labels(),
                            component: part.component,
                            irnumber: part.irnumber,
                        },
                    });
                }
            }
        } else {
            if owner.dim != 1 {
                return Err(error("the standalone owner dimension is not one".into()));
            }
            let row = owner
                .ordinary_scalar_selected_arm_block_trace()
                .map_err(|e| error(format!("typed scalar row: {e:?}")))?;
            found.push(HRow {
                row,
                storage: HCharacterStorage::Standalone {
                    record: owner.labels(),
                },
            });
        }
    }
    if found.len() != 1 {
        return Err(error(format!(
            "{} authoritative scalar H rows carry this label",
            found.len()
        )));
    }
    Ok(found.remove(0))
}

/// Certify the entire H row (including redundant centered translations), then
/// bind through the same strict H data-Hall Seitz map used by the real targets.
pub(super) fn unitary_values(
    h: &HRow,
    corep: &MagneticCorepSummary,
    context: &TargetContext,
) -> Result<Vec<Complex64>, MatrixSubductionError> {
    let values = bind_h_row(h, &corep.source_irreps[0], &corep.label, context)?;
    let fail = |s| data_error(context, &corep.label, s);
    for (class, value) in values.iter().enumerate() {
        if context.classes[class].time_reversal {
            continue;
        }
        if context.class_columns[class].is_empty() {
            return Err(fail("uncovered unitary class".into()));
        }
        for &column in &context.class_columns[class] {
            let summary = corep
                .characters
                .get(column)
                .copied()
                .flatten()
                .ok_or_else(|| fail("missing summary character".into()))?;
            if !close(summary.re, value.re) || !close(summary.im, value.im) {
                return Err(fail(format!(
                    "typed H scalar and summary disagree at class {class}"
                )));
            }
        }
    }
    Ok(values)
}

/// Faithful scalar realification, shared by H rows and explicitly admitted CIR sources.
/// The caller establishes the row's authoritative identity and full-frame binding.
pub(super) fn certify_scalar_characters(
    dimension: usize,
    source_operations: &[crate::irrep::types::SeitzOperation],
    values: &[Complex64],
) -> Result<CertifiedTable, MatrixCertificationProblem> {
    if dimension != 1 || values.is_empty() || values.len() != source_operations.len() {
        return Err(MatrixCertificationProblem::TypedRowShapeMismatch {
            row: values.len(),
            operators: source_operations.len(),
            dimension,
        });
    }
    let mut operations = Vec::new();
    let mut lattice = Vec::new();
    let mut matrices = Vec::new();
    for (i, (op, value)) in source_operations.iter().zip(values).enumerate() {
        if !value.re.is_finite() || !value.im.is_finite() {
            return Err(MatrixCertificationProblem::NonFinite { operation: i });
        }
        if !close(value.norm_sqr(), 1.0) {
            return Err(MatrixCertificationProblem::NotOrthogonal { operation: i });
        }
        let rotation = std::array::from_fn(|r| std::array::from_fn(|c| op.rotation[3 * r + c]));
        let mut translation = [Rat::ZERO; 3];
        for (axis, x) in translation.iter_mut().enumerate() {
            *x = exact_fraction(op.translation[axis])
                .ok_or(MatrixCertificationProblem::TranslationOffGrid { operation: i })?;
        }
        if rotation == identity_rotation()
            && translation != [Rat::ZERO; 3]
            && !lattice.contains(&translation)
        {
            lattice.push(translation);
        }
        operations.push((rotation, translation));
        matrices.push(vec![value.re, -value.im, value.im, value.re]);
    }
    check_stored_matrix_law(&operations, &matrices, 2, &lattice)?;
    Ok(CertifiedTable {
        dimension: 2,
        operations,
        matrices,
    })
}

/// Certify a complete authoritative scalar H row and bind it to M's unitary
/// classes. The caller must independently check its target's summary character.
pub(super) fn bind_h_row(
    h: &HRow,
    source: &SourceIrrepSummary,
    label: &str,
    context: &TargetContext,
) -> Result<Vec<Complex64>, MatrixSubductionError> {
    let fail = |s| data_error(context, label, s);
    certify_scalar_characters(h.row.dimension(), h.row.operations(), h.row.values()).map_err(
        |problem| match problem {
            MatrixCertificationProblem::TypedRowShapeMismatch { .. } => {
                fail("invalid scalar H row shape".into())
            }
            MatrixCertificationProblem::NonFinite { operation }
            | MatrixCertificationProblem::NotOrthogonal { operation } => fail(format!(
                "H row {operation} is not a finite unit-modulus scalar"
            )),
            MatrixCertificationProblem::TranslationOffGrid { operation } => {
                fail(format!("H translation {operation} is off the exact grid"))
            }
            problem => MatrixSubductionError::TargetMatrixUncertified {
                uni: context.uni,
                sg: context.h_sg,
                label: source.ml,
                problem,
            },
        },
    )?;
    let rotations: Vec<_> = h.row.operations().iter().flat_map(|o| o.rotation).collect();
    let translations: Vec<_> = h
        .row
        .operations()
        .iter()
        .flat_map(|o| o.translation)
        .collect();
    let map = wigner::build_h_to_irrep_op_map(&context.h_seitz, &rotations, &translations)
        .ok_or_else(|| {
            fail("the complete typed H row differs from the H data-Hall Seitz list".into())
        })?;
    let mut values = vec![Complex64::new(0.0, 0.0); context.classes.len()];
    for (class, slot) in values.iter_mut().enumerate() {
        if context.classes[class].time_reversal {
            continue;
        }
        let h_index = context
            .unitary_operation_index(class)
            .map_err(|error| label_error(error, context.uni, label))?;
        *slot = h.row.values()[map[h_index]];
    }
    Ok(values)
}

pub(super) fn scalar_table(
    corep: &MagneticCorepSummary,
    context: &TargetContext,
    values: &[Complex64],
    anti: Option<usize>,
) -> Result<Vec<TargetOperator>, MatrixSubductionError> {
    let operators = magnetic_subduction::build_scalar_target(
        corep,
        &context.classes,
        &context.lattice,
        values,
        anti,
        context.uni,
    )
    .map_err(|e| match e {
        magnetic_subduction::MagneticSubductionError::Arithmetic { reason } => {
            MatrixSubductionError::Arithmetic { reason }
        }
        other => data_error(context, &corep.label, other.to_string()),
    })?;
    // Checks the square and twisted conjugation too; no phase is guessed for a
    // type-B square. Composition conjugates the right coefficient when needed.
    for (left, g) in context.classes.iter().enumerate() {
        for (right, h) in context.classes.iter().enumerate() {
            let product = compose(g, h).map_err(|e| MatrixSubductionError::Arithmetic {
                reason: e.to_string(),
            })?;
            let index = find_model_class(&product, &context.classes, &context.lattice)
                .ok_or_else(|| data_error(context, &corep.label, "missing product class"))?;
            let actual = operators[left].compose(&operators[right]);
            let expected = operators[index];
            if actual.is_antiunitary() != expected.is_antiunitary()
                || !close(actual.value().re, expected.value().re)
                || !close(actual.value().im, expected.value().im)
            {
                return Err(data_error(
                    context,
                    &corep.label,
                    format!("scalar target law fails at {left},{right}"),
                ));
            }
        }
    }
    Ok(operators)
}

/// Decompose the existing real 2D Γ model into non-spinor one-dimensional
/// type-A targets with real or complex H characters. Higher target dimensions,
/// type B/C and derived/unavailable CIR rows are not silently substituted.
/// The legacy real-matrix and paired-C entry points are unchanged.
pub fn subduce_matrix_gamma_scalar_targets(
    model: &MatrixGammaModel,
) -> Result<ScalarTargetSubduction<'_>, MatrixSubductionError> {
    let (gamma, context) = prepare_matrix_targets(model)?;
    let source = RealCorep::from_source(&model.actions);
    let mut pieces = Vec::new();
    let mut examined = Vec::new();
    for corep in &gamma.coreps {
        if corep.corep_type != CorepType::A
            || corep.dim != 1
            || corep.source_irreps.len() != 1
            || corep.source_irreps[0].spinor
            || corep.source_irreps[0].dim != 1
        {
            continue;
        }
        examined.push(corep.label.clone());
        if !has_unitary_overlap(model, corep, &context)? {
            continue;
        }
        let h = resolve_h_row(&corep.source_irreps[0], context.uni)?;
        let values = unitary_values(&h, corep, &context)?;
        let operators = scalar_table(
            corep,
            &context,
            &values,
            context.classes.iter().position(|g| g.time_reversal),
        )?;
        let target = RealCorep {
            dimension: 1,
            matrices: operators.iter().map(|o| vec![o.value().re]).collect(),
            antiunitary: operators
                .iter()
                .map(TargetOperator::is_antiunitary)
                .collect(),
        };
        let imaginary: Vec<_> = operators.iter().map(|o| vec![o.value().im]).collect();
        let (rows, width) = hom_rows_with_imaginary_target(&source, &target, Some(&imaginary));
        let copies = copies_from_kernel(&real_kernel(&rows, width), 1);
        if copies.is_empty() {
            continue;
        }
        verify_reconstruction_with(model, &corep.label, 1, &copies, |class, v| {
            vec![operators[class].apply(v[0])]
        })?;
        pieces.push(ScalarTargetPiece {
            label: corep.label.clone(),
            source: corep.source_irreps[0].clone(),
            storage: h.storage,
            copies,
            operators,
        });
    }
    if pieces.is_empty() {
        return Err(MatrixSubductionError::NoMatchingTarget {
            parent_sg: model.parent_sg,
            source: model.source_label,
            uni: model.uni,
            examined,
        });
    }
    let columns: Vec<_> = pieces
        .iter()
        .flat_map(|p| p.copies.iter().cloned())
        .collect();
    if columns.len() != SOURCE_DIMENSION {
        return Err(MatrixSubductionError::IncompleteDecomposition {
            source_dimension: SOURCE_DIMENSION as u32,
            covered_dimension: columns.len() as u32,
        });
    }
    if !verify_column_isometry(&columns) {
        return Err(MatrixSubductionError::SourceSpaceNotIsometric {
            covered_dimension: columns.len() as u32,
        });
    }
    Ok(ScalarTargetSubduction { model, pieces })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::irrep::subduction::Mat3R;

    fn witness() -> (MatrixGammaModel, MagneticCorepSummary, TargetContext) {
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
        let corep = gamma
            .coreps
            .into_iter()
            .find(|c| {
                c.corep_type == CorepType::A
                    && c.dim == 1
                    && c.source_irreps.len() == 1
                    && c.source_irreps[0].ml == "GM3"
                    && !c.source_irreps[0].spinor
            })
            .unwrap();
        (m, corep, context)
    }

    #[test]
    fn missing_and_wrong_dimensional_h_sources_are_errors() {
        let (_, corep, context) = witness();
        let good = &corep.source_irreps[0];
        assert!(matches!(
            resolve_h_row(good, context.uni).unwrap().storage,
            HCharacterStorage::Constituent { .. }
        ));
        let mut missing = good.clone();
        missing.ml = "not-a-source";
        assert!(matches!(
            resolve_h_row(&missing, context.uni),
            Err(MatrixSubductionError::TargetSourceUnresolved { .. })
        ));
        let mut bad = good.clone();
        bad.dim = 2;
        assert!(resolve_h_row(&bad, context.uni).is_err());
        let mut bad = good.clone();
        bad.spinor = true;
        assert!(resolve_h_row(&bad, context.uni).is_err());
    }

    #[test]
    fn typed_binding_preserves_missing_nonfinite_and_frame_errors() {
        let (_, corep, context) = witness();
        let h = resolve_h_row(&corep.source_irreps[0], context.uni).unwrap();
        assert!(unitary_values(&h, &corep, &context).is_ok());
        let index = context
            .classes
            .iter()
            .position(|g| !g.time_reversal)
            .unwrap();
        let column = context.class_columns[index][0];
        for value in [
            None,
            Some(Complex64::new(f64::NAN, 0.0)),
            Some(Complex64::new(f64::INFINITY, 0.0)),
            Some(Complex64::new(0.5, 0.0)),
        ] {
            let mut bad = corep.clone();
            bad.characters[column] = value;
            assert!(unitary_values(&h, &bad, &context).is_err());
        }
        let (_, corep, mut context) = witness();
        context.h_seitz.clear();
        assert!(unitary_values(&h, &corep, &context).is_err());
        let (_, corep, mut context) = witness();
        context.classes[index].translation[0] = Rat::new(1, 7).unwrap();
        let error = unitary_values(&h, &corep, &context).expect_err("the class is not in H");
        assert!(
            matches!(error, MatrixSubductionError::TargetCorepInconsistent { label, .. }
            if label == corep.label),
            "operation errors must retain the candidate identity"
        );
    }

    #[test]
    fn scalar_coset_gauges_reconstruct_and_negative_square_is_rejected() {
        let (m, corep, context) = witness();
        let h = resolve_h_row(&corep.source_irreps[0], context.uni).unwrap();
        let values = unitary_values(&h, &corep, &context).unwrap();
        let source = RealCorep::from_source(&m.actions);
        for anti in (0..context.classes.len()).filter(|&i| context.classes[i].time_reversal) {
            let table = scalar_table(&corep, &context, &values, Some(anti)).unwrap();
            let target = RealCorep {
                dimension: 1,
                matrices: table.iter().map(|o| vec![o.value().re]).collect(),
                antiunitary: table.iter().map(TargetOperator::is_antiunitary).collect(),
            };
            let imaginary: Vec<_> = table.iter().map(|o| vec![o.value().im]).collect();
            let (rows, width) = hom_rows_with_imaginary_target(&source, &target, Some(&imaginary));
            let copies = copies_from_kernel(&real_kernel(&rows, width), 1);
            assert_eq!(copies.len(), 1);
            verify_reconstruction_with(&m, &corep.label, 1, &copies, |i, v| {
                vec![table[i].apply(v[0])]
            })
            .unwrap();
        }
        // Abstract 4' with H=2 and χ(C2)=-1 cannot have a scalar type-A extension.
        let rotations = [
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
            [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
            [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
        ];
        let context = TargetContext {
            uni: 750,
            h_sg: 3,
            h_seitz: vec![],
            setting: None,
            lattice: vec![],
            class_columns: vec![],
            classes: rotations
                .into_iter()
                .enumerate()
                .map(|(i, rotation)| MagneticOperation {
                    rotation,
                    translation: [Rat::ZERO; 3],
                    time_reversal: i % 2 == 1,
                })
                .collect(),
        };
        let mut values = vec![Complex64::new(1.0, 0.0); 4];
        assert!(scalar_table(&corep, &context, &values, Some(1)).is_ok());
        values[2] = Complex64::new(-1.0, 0.0);
        assert!(scalar_table(&corep, &context, &values, Some(1)).is_err());
    }
}
