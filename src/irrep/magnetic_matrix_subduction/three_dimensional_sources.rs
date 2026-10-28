//! Explicit real, nonspinor, noncompound three-dimensional Gamma source certification.
//! This module certifies the entire parent table and its restricted semilinear
//! action. It does not yet decompose that action into magnetic target blocks.
use super::*;

/// Scope of the explicit three-dimensional source verifier (not a subduction report).
pub const THREE_DIMENSIONAL_SOURCE_SCOPE: &str =
    "certified real nonspinor noncompound three-dimensional Gamma source; no target decomposition";

/// A certified real 3x3 coefficient matrix, optionally followed by conjugation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThreeDimensionalSourceOperator {
    entries: [f64; 9],
    antiunitary: bool,
}
impl ThreeDimensionalSourceOperator {
    /// A row-major entry, or None outside the 3x3 block.
    pub fn entry(&self, row: usize, column: usize) -> Option<f64> {
        (row < 3 && column < 3).then(|| self.entries[row * 3 + column])
    }
    /// All nine row-major real coefficients.
    pub fn entries(&self) -> [f64; 9] {
        self.entries
    }
    /// Whether complex conjugation acts before the coefficient matrix.
    pub fn is_antiunitary(&self) -> bool {
        self.antiunitary
    }
    /// Apply the operator to a three-component complex vector.
    pub fn apply(&self, vector: [Complex64; 3]) -> [Complex64; 3] {
        let vector = if self.antiunitary {
            vector.map(|z| z.conj())
        } else {
            vector
        };
        std::array::from_fn(|r| (0..3).map(|j| self.entries[r * 3 + j] * vector[j]).sum())
    }
}

/// Errors of the separate three-dimensional source path; old errors retain their causes.
#[derive(Debug)]
pub enum ThreeDimensionalGammaError {
    /// A source of any dimension other than exactly three is outside this entry point.
    SourceDimensionNotThree {
        sg: u8,
        label: &'static str,
        dimension: u8,
    },
    /// An existing identity, data, frame, arithmetic or geometry certificate failed.
    Validation(Box<MatrixSubductionError>),
    /// The source's semilinear operators fail one ordered group product.
    SourceLawViolated {
        left: usize,
        right: usize,
        expected_class: usize,
    },
}
impl From<MatrixSubductionError> for ThreeDimensionalGammaError {
    fn from(value: MatrixSubductionError) -> Self {
        Self::Validation(Box::new(value))
    }
}
impl std::fmt::Display for ThreeDimensionalGammaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SourceDimensionNotThree {
                sg,
                label,
                dimension,
            } => write!(
                f,
                "SG{sg} source {label} has dimension {dimension}, not three"
            ),
            Self::Validation(error) => error.fmt(f),
            Self::SourceLawViolated {
                left,
                right,
                expected_class,
            } => write!(
                f,
                "three-dimensional source law fails at {left} * {right} -> {expected_class}"
            ),
        }
    }
}
impl std::error::Error for ThreeDimensionalGammaError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Validation(error) => Some(error.as_ref()),
            _ => None,
        }
    }
}

/// A verified three-dimensional source with independent source and embedding origins.
/// Fields are private; only the explicit verifier below constructs this model.
#[derive(Debug)]
pub struct ThreeDimensionalGammaModel {
    parent_sg: u8,
    source_label: &'static str,
    source_convention: LabelConvention,
    input_convention: LabelConvention,
    record_label: &'static str,
    record_convention: LabelConvention,
    direction_index: usize,
    ordinal: usize,
    uni: usize,
    map: crate::irrep::subduction::Mat3R,
    shift: [Rat; 3],
    source_origin_shift: [Rat; 3],
    lattice_policy: LatticePolicy,
    realises_record_lattice: bool,
    time_reversal: Option<RealSign>,
    actions: Vec<ThreeDimensionalSourceOperator>,
    classes: Vec<MagneticOperation>,
    lattice: Vec<[Rat; 3]>,
}

impl ThreeDimensionalGammaModel {
    /// Parent space group of the source and the record.
    pub fn parent_sg(&self) -> u8 {
        self.parent_sg
    }

    /// The source parent irrep label, in its canonical CDML spelling.
    pub fn source_label(&self) -> &'static str {
        self.source_label
    }

    /// The convention [`Self::source_label`] is written in (`Cdml`).
    pub fn source_convention(&self) -> LabelConvention {
        self.source_convention
    }

    /// The convention the caller used to name the source label.
    pub fn input_convention(&self) -> LabelConvention {
        self.input_convention
    }

    /// The condensing irrep label of the record, in its canonical CDML spelling.
    pub fn record_label(&self) -> &'static str {
        self.record_label
    }

    /// The convention [`Self::record_label`] is written in (`Cdml`).
    pub fn record_convention(&self) -> LabelConvention {
        self.record_convention
    }

    /// The record's row index in the condensing irrep's magnetic isotropy table.
    pub fn direction_index(&self) -> usize {
        self.direction_index
    }

    /// The record's ordinal in the generated magnetic isotropy table.
    pub fn ordinal(&self) -> usize {
        self.ordinal
    }

    /// The record's magnetic space group (UNI) number.
    pub fn uni(&self) -> usize {
        self.uni
    }

    /// The verified map of the embedding.
    pub fn map(&self) -> &crate::irrep::subduction::Mat3R {
        &self.map
    }

    /// The verified origin shift of the embedding.
    pub fn shift(&self) -> &[Rat; 3] {
        &self.shift
    }

    /// A common origin witness from the certified source operation frame to the
    /// parent frame: t_parent = t_source + (I - R) delta, modulo the parent lattice.
    /// This is distinct from the magnetic embedding's [`Self::shift`]. Gamma
    /// matrices remain unchanged; this does not certify a raw PIR basis/gauge.
    pub fn source_origin_shift(&self) -> &[Rat; 3] {
        &self.source_origin_shift
    }

    /// The lattice policy the model was checked against.
    pub fn lattice_policy(&self) -> LatticePolicy {
        self.lattice_policy
    }

    /// Whether the selected map realises the record's lattice exactly.
    pub fn realises_record_lattice(&self) -> bool {
        self.realises_record_lattice
    }

    /// The parent action's `η`, or `None` when the caller stated none.
    pub fn time_reversal(&self) -> Option<RealSign> {
        self.time_reversal
    }

    /// The source dimension, always three for this entry point.
    pub fn source_dimension(&self) -> u32 {
        3
    }

    /// The certified source action of every class, in the set's class order.
    pub fn source_operators(&self) -> &[ThreeDimensionalSourceOperator] {
        &self.actions
    }

    /// The certified source action of one class, by class index.
    pub fn source_operator_at(&self, index: usize) -> Option<&ThreeDimensionalSourceOperator> {
        self.actions.get(index)
    }

    /// The classes of the magnetic group, in the set's class order.
    pub fn classes(&self) -> &[MagneticOperation] {
        &self.classes
    }

    /// The group's translation lattice, as the class comparisons reduced it.
    pub fn translation_lattice(&self) -> &[[Rat; 3]] {
        &self.lattice
    }

    /// One line stating what this entry point covers.
    pub fn scope(&self) -> &'static str {
        THREE_DIMENSIONAL_SOURCE_SCOPE
    }
}

/// Certify exactly a real nonspinor noncompound 3D Gamma source.
/// Uses the existing request/evidence container, never the old 2D admission gate.
/// No magnetic target decomposition is returned by this function.
pub fn verify_three_dimensional_gamma_input(
    request: &MatrixGammaRequest<'_>,
) -> Result<ThreeDimensionalGammaModel, ThreeDimensionalGammaError> {
    let source_identity = request.source;
    let record_identity = request.record;
    check_matrix_parent(source_identity, record_identity)?;
    let parent_sg = source_identity.parent_sg;
    let source_irrep = resolve_three_dimensional_source(source_identity)?;

    let PreparedMatrixRecord {
        record_label,
        ordinal,
        set,
        evidence,
        realises_record_lattice,
        time_reversal,
    } = prepare_matrix_record(request)?;

    let mut certified = certify_irrep_matrices(source_irrep, 3).map_err(|problem| {
        MatrixSubductionError::SourceMatrixUncertified {
            sg: parent_sg,
            label: source_irrep.ml,
            problem,
        }
    })?;

    let parent = ParentFrame::of(parent_sg).map_err(|error| {
        MatrixSubductionError::ParentFrameUnavailable {
            sg: parent_sg,
            reason: format!("{error:?}"),
        }
    })?;
    // Keep the matrix/operation row correspondence. The shared scalar proof
    // anchors the whole source to its data Hall, including centring and every
    // rotation, before solving one common origin for the complete group.
    let source_origin_shift = parent
        .align_source_operations(parent_sg, set.uni, &mut certified.operations)
        .map_err(|error| MatrixSubductionError::SourceOriginAlignment {
            sg: parent_sg,
            label: source_irrep.ml,
            reason: Box::new(error),
        })?;
    let mapped = magnetic_input::mapped_classes(&set, &evidence, &parent).map_err(|error| {
        MatrixSubductionError::ClassMapping {
            reason: format!("{error:?}"),
        }
    })?;

    let eta = time_reversal.map(|sign| sign.value() as f64).unwrap_or(1.0);
    let mut actions = Vec::with_capacity(mapped.len());
    for mapped_class in &mapped {
        let index = bind_source_operation(
            &certified,
            source_irrep.ml,
            parent_sg,
            mapped_class,
            &parent,
        )?;
        let mut entries: [f64; 9] = certified.matrices[index]
            .as_slice()
            .try_into()
            .expect("certified 3x3 source");
        if mapped_class.time_reversal {
            for value in entries.iter_mut() {
                *value *= eta;
            }
        }
        actions.push(ThreeDimensionalSourceOperator {
            entries,
            antiunitary: mapped_class.time_reversal,
        });
    }

    let classes = set.classes.clone();
    let lattice = magnetic_embedding::translation_lattice(&set);
    check_three_dimensional_source_law(&classes, &actions, &lattice)?;

    Ok(ThreeDimensionalGammaModel {
        parent_sg,
        source_label: source_irrep.ml,
        source_convention: LabelConvention::Cdml,
        input_convention: source_identity.convention,
        record_label,
        record_convention: LabelConvention::Cdml,
        direction_index: record_identity.direction_index,
        ordinal,
        uni: set.uni,
        map: evidence.map,
        shift: evidence.shift,
        source_origin_shift,
        lattice_policy: evidence.lattice_policy,
        realises_record_lattice,
        time_reversal,
        actions,
        classes,
        lattice,
    })
}

fn resolve_three_dimensional_source(
    source: SourceIdentity<'_>,
) -> Result<&'static IrrepRecord, ThreeDimensionalGammaError> {
    let irrep = resolve_convention_irrep(source.parent_sg, source.label, source.convention)?;
    if irrep.spinor {
        return Err(MatrixSubductionError::SourceIsSpinor {
            sg: source.parent_sg,
            label: irrep.ml,
        }
        .into());
    }
    if irrep.compound_metadata().is_some() {
        return Err(MatrixSubductionError::SourceIsCompound {
            sg: source.parent_sg,
            label: irrep.ml,
        }
        .into());
    }
    if irrep.dim != 3 {
        return Err(ThreeDimensionalGammaError::SourceDimensionNotThree {
            sg: source.parent_sg,
            label: irrep.ml,
            dimension: irrep.dim,
        });
    }
    let k = irrep.k_vector();
    if k.numerators != [0, 0, 0] || k.denominator != 1 {
        return Err(MatrixSubductionError::SourceNotGamma {
            sg: source.parent_sg,
            label: irrep.ml,
            numerators: k.numerators,
            denominator: k.denominator,
        }
        .into());
    }
    Ok(irrep)
}

fn check_three_dimensional_source_law(
    classes: &[MagneticOperation],
    actions: &[ThreeDimensionalSourceOperator],
    lattice: &[[Rat; 3]],
) -> Result<(), ThreeDimensionalGammaError> {
    for left in 0..classes.len() {
        for right in 0..classes.len() {
            let product = compose(&classes[left], &classes[right]).map_err(|error| {
                MatrixSubductionError::Arithmetic {
                    reason: error.to_string(),
                }
            })?;
            let expected_class = find_model_class(&product, classes, lattice)
                .ok_or(MatrixSubductionError::SourceLawUnmatched { left, right })?;
            let entries = matrix_product(&actions[left].entries, &actions[right].entries, 3);
            if actions[left].antiunitary ^ actions[right].antiunitary
                != actions[expected_class].antiunitary
                || !matrices_close(&entries, &actions[expected_class].entries)
            {
                return Err(ThreeDimensionalGammaError::SourceLawViolated {
                    left,
                    right,
                    expected_class,
                });
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::irrep::subduction::{Mat3R, SeitzTransform, Vec3R};

    #[test]
    fn three_dimensional_storage_keeps_dimension_and_full_source_frame_gates() {
        let ir = query::irreps_of(195)
            .iter()
            .find(|r| r.ml == "GM4")
            .unwrap();
        for dimension in [0, 4] {
            let mut fake = *ir;
            fake.dim = dimension as u8;
            assert!(
                matches!(certify_irrep_matrices(&fake,dimension),Err(MatrixCertificationProblem::DimensionOutOfScope{dimension:d}) if d==dimension)
            );
        }
        assert!(matches!(
            certify_irrep_matrices(ir, 2),
            Err(MatrixCertificationProblem::DeclaredDimensionMismatch {
                declared: 3,
                expected: 2
            })
        ));
        let mut table = certify_irrep_matrices(ir, 3).unwrap();
        let parent = ParentFrame::of(195).unwrap();
        let change = SeitzTransform::new(
            Mat3R::identity(),
            Vec3R::new([Rat::new(1, 7).unwrap(), Rat::ZERO, Rat::ZERO]),
        );
        for (r, t) in &mut table.operations {
            *t = *change
                .map_translation(*r, &Vec3R::new(*t))
                .unwrap()
                .as_array();
        }
        // Algebra survives an origin change, but that does not authorize a
        // different stored data-Hall frame for this named source.
        check_stored_matrix_law(&table.operations, &table.matrices, 3, &parent.lattice).unwrap();
        assert!(matches!(
            parent.align_source_operations(195, 0, &mut table.operations),
            Err(magnetic_input::MagneticInputError::CharacterRowSettingMismatch { .. })
        ));
    }

    #[test]
    fn three_dimensional_source_law_rejects_wrong_squares_and_preserves_arithmetic() {
        let id = MagneticOperation {
            rotation: identity_rotation(),
            translation: [Rat::ZERO; 3],
            time_reversal: false,
        };
        let mut anti = id;
        anti.time_reversal = true;
        let identity = ThreeDimensionalSourceOperator {
            entries: [1., 0., 0., 0., 1., 0., 0., 0., 1.],
            antiunitary: false,
        };
        let a = ThreeDimensionalSourceOperator {
            entries: [1., 0., 0., 0., -1., 0., 0., 0., -1.],
            antiunitary: true,
        };
        check_three_dimensional_source_law(&[id, anti], &[identity, a], &[]).unwrap();
        let mut wrong = a;
        wrong.entries[0] = 2.;
        assert!(matches!(
            check_three_dimensional_source_law(&[id, anti], &[identity, wrong], &[]),
            Err(ThreeDimensionalGammaError::SourceLawViolated {
                left: 1,
                right: 1,
                expected_class: 0
            })
        ));
        let mut overflow = id;
        overflow.translation[0] = Rat::from_integer(i128::MAX);
        assert!(
            matches!(check_three_dimensional_source_law(&[overflow],&[identity],&[]),Err(ThreeDimensionalGammaError::Validation(e)) if matches!(*e,MatrixSubductionError::Arithmetic{..}))
        );
        let mut not_closed = id;
        not_closed.rotation = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        assert!(
            matches!(check_three_dimensional_source_law(&[not_closed],&[identity],&[]),Err(ThreeDimensionalGammaError::Validation(e)) if matches!(*e,MatrixSubductionError::SourceLawUnmatched{left:0,right:0}))
        );
    }
}
