//! Bounded mixed-target decomposition of the independently certified real 3D source.
use super::scalar_targets::{bind_h_row, resolve_h_row, scalar_table, unitary_values};
use super::*;
use crate::irrep::labels::IrrepLabels;
use crate::irrep::magnetic_subduction::TargetOperator;
use crate::irrep::magnetic_summary::SourceIrrepSummary;

/// Explicit scope; the old two-dimensional entry points are not broadened.
pub const THREE_DIMENSIONAL_TARGET_SCOPE: &str = "Certified real nonspinor noncompound 3D Gamma source; mixed scalar A, real matrix A2/A3, real-scalar B2, real or conjugate-scalar C2 and pure-unitary conjugate packets; full-source Hom, Gram, projector and complex reconstruction. No general complex source, spinor, non-Gamma or setting search.";

/// A pure-unitary packet is not a complex irreducible representation or magnetic C.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThreeDimensionalTargetKind {
    TypeA,
    TypeB,
    TypeC,
    UnitaryPair,
}

/// The actual stored authority, not an invented standalone CIR record.
#[derive(Debug, Clone)]
pub enum ThreeDimensionalHStorage {
    /// A complete certified real matrix table.
    Matrix { record: IrrepLabels },
    /// A certified scalar row, with its standalone or named CIR storage.
    Scalar(HCharacterStorage),
}

/// Reuse the established scalar/real operators rather than a new matrix framework.
#[derive(Debug, Clone)]
pub enum ThreeDimensionalTargetOperator {
    Scalar(TargetOperator),
    Real(TargetMatrixOperator),
}
impl ThreeDimensionalTargetOperator {
    /// Target dimension, one through three in a successful report.
    pub fn dimension(&self) -> usize {
        match self {
            Self::Scalar(_) => 1,
            Self::Real(o) => o.dimension(),
        }
    }
    /// Whether conjugation precedes the coefficient matrix.
    pub fn is_antiunitary(&self) -> bool {
        match self {
            Self::Scalar(o) => o.is_antiunitary(),
            Self::Real(o) => o.is_antiunitary(),
        }
    }
    /// A complex coefficient, or None outside the block.
    pub fn entry(&self, row: usize, column: usize) -> Option<Complex64> {
        match self {
            Self::Scalar(o) => (row == 0 && column == 0).then(|| o.value()),
            Self::Real(o) => o.entry(row, column).map(|x| Complex64::new(x, 0.0)),
        }
    }
    /// Apply to a vector of exactly the target dimension; wrong lengths are refused.
    pub fn apply(&self, vector: &[Complex64]) -> Option<Vec<Complex64>> {
        if vector.len() != self.dimension() {
            return None;
        }
        Some(match self {
            Self::Scalar(o) => vec![o.apply(vector[0])],
            Self::Real(o) => o.apply(vector),
        })
    }
}

/// One target, its authoritative H identities and normalized row-major 3xd copies.
#[derive(Debug, Clone)]
pub struct ThreeDimensionalSubductionPiece {
    label: String,
    kind: ThreeDimensionalTargetKind,
    h_sources: Vec<SourceIrrepSummary>,
    storage: Vec<ThreeDimensionalHStorage>,
    dimension: usize,
    copies: Vec<Vec<Complex64>>,
    operators: Vec<ThreeDimensionalTargetOperator>,
}
impl ThreeDimensionalSubductionPiece {
    pub fn label(&self) -> &str {
        &self.label
    }
    pub fn kind(&self) -> ThreeDimensionalTargetKind {
        self.kind
    }
    pub fn h_sources(&self) -> &[SourceIrrepSummary] {
        &self.h_sources
    }
    pub fn storage(&self) -> &[ThreeDimensionalHStorage] {
        &self.storage
    }
    pub fn dimension(&self) -> u32 {
        self.dimension as u32
    }
    /// Actual orthonormal copies, not the real Hom kernel dimension.
    pub fn multiplicity(&self) -> u32 {
        self.copies.len() as u32
    }
    pub fn copies(&self) -> &[Vec<Complex64>] {
        &self.copies
    }
    pub fn operators(&self) -> &[ThreeDimensionalTargetOperator] {
        &self.operators
    }
}

/// A single complete three-dimensional decomposition, never concatenated old reports.
#[derive(Debug)]
pub struct ThreeDimensionalSubduction<'a> {
    model: &'a ThreeDimensionalGammaModel,
    pieces: Vec<ThreeDimensionalSubductionPiece>,
}
impl ThreeDimensionalSubduction<'_> {
    pub fn model(&self) -> &ThreeDimensionalGammaModel {
        self.model
    }
    pub fn pieces(&self) -> &[ThreeDimensionalSubductionPiece] {
        &self.pieces
    }
    pub fn source_dimension(&self) -> u32 {
        3
    }
    pub fn covered_dimension(&self) -> u32 {
        self.pieces
            .iter()
            .map(|p| p.dimension() * p.multiplicity())
            .sum()
    }
    pub fn scope(&self) -> &'static str {
        THREE_DIMENSIONAL_TARGET_SCOPE
    }
}

/// Old target data/arithmetic failures retain their causes; new reconstruction is separate.
#[derive(Debug)]
pub enum ThreeDimensionalSubductionError {
    Target(Box<MatrixSubductionError>),
    RealTypeC(Box<MatrixTypeCError>),
    Reconstruction {
        uni: usize,
        label: String,
        reason: String,
    },
}
impl From<MatrixSubductionError> for ThreeDimensionalSubductionError {
    fn from(e: MatrixSubductionError) -> Self {
        Self::Target(Box::new(e))
    }
}
impl From<MatrixTypeCError> for ThreeDimensionalSubductionError {
    fn from(e: MatrixTypeCError) -> Self {
        Self::RealTypeC(Box::new(e))
    }
}
impl std::fmt::Display for ThreeDimensionalSubductionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Target(e) => e.fmt(f),
            Self::RealTypeC(e) => e.fmt(f),
            Self::Reconstruction { uni, label, reason } => {
                write!(f, "UNI {uni} target {label}: {reason}")
            }
        }
    }
}
impl std::error::Error for ThreeDimensionalSubductionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Target(e) => Some(e.as_ref()),
            Self::RealTypeC(e) => Some(e.as_ref()),
            _ => None,
        }
    }
}
fn failure(uni: usize, label: &str, reason: &str) -> ThreeDimensionalSubductionError {
    ThreeDimensionalSubductionError::Reconstruction {
        uni,
        label: label.into(),
        reason: reason.into(),
    }
}

// Recheck the complex matrix equation independently of the real nullspace rows.
fn intertwines(
    source: &RealCorep,
    target: &RealCorep,
    imaginary: &[Vec<f64>],
    u: &[Complex64],
) -> bool {
    let d = target.dimension;
    if source.dimension != 3 || u.len() != 3 * d || source.antiunitary != target.antiunitary {
        return false;
    }
    for (g, source_matrix) in source.matrices.iter().enumerate() {
        for r in 0..3 {
            for c in 0..d {
                let left: Complex64 = (0..3)
                    .map(|k| {
                        source_matrix[r * 3 + k]
                            * if source.antiunitary[g] {
                                u[k * d + c].conj()
                            } else {
                                u[k * d + c]
                            }
                    })
                    .sum();
                let right: Complex64 = (0..d)
                    .map(|k| {
                        u[r * d + k]
                            * Complex64::new(target.matrices[g][k * d + c], imaginary[g][k * d + c])
                    })
                    .sum();
                if !(left - right).norm().is_finite() || (left - right).norm() > 1e-8 {
                    return false;
                }
            }
        }
    }
    true
}

fn copies_from_basis(
    source: &RealCorep,
    target: &RealCorep,
    imaginary: &[Vec<f64>],
    kernel: &[Vec<f64>],
    unitary_pair: bool,
) -> Vec<Vec<Complex64>> {
    let d = target.dimension;
    let mut copies: Vec<Vec<Complex64>> = Vec::new();
    for vector in kernel {
        let u = intertwiner_from_vector(vector, d);
        // I+iJ can be rank one on a complex-reducible unitary packet. Its real
        // components are real intertwiners of the irreducible real H carrier.
        let candidates = if unitary_pair {
            vec![
                u.iter().map(|z| Complex64::new(z.re, 0.0)).collect(),
                u.iter().map(|z| Complex64::new(z.im, 0.0)).collect(),
            ]
        } else {
            vec![u]
        };
        for mut u in candidates {
            if d == 1 {
                // A repeated scalar may have a nonorthogonal RREF basis. Remove
                // the already represented image instead of discarding that vector.
                for previous in &copies {
                    let overlap: Complex64 =
                        previous.iter().zip(&u).map(|(a, b)| a.conj() * b).sum();
                    for (z, p) in u.iter_mut().zip(previous) {
                        *z -= p * overlap;
                    }
                }
            }
            if u.iter().map(|z| z.norm_sqr()).sum::<f64>() <= 1e-16 {
                continue;
            }
            let Some(u) = normalize_intertwiner_for_source(&u, d, 3) else {
                continue;
            };
            if !intertwines(source, target, imaginary, &u) {
                continue;
            }
            copies.push(u);
            // ponytail: in source3 a d>=2 block fits once; larger sources need
            // a real multiplicity-space algorithm, not truncating this loop.
            if d >= 2 {
                return copies;
            }
        }
    }
    copies
}

fn matching_three_dimensional_copies(
    source: &RealCorep,
    operators: &[ThreeDimensionalTargetOperator],
    unitary_pair: bool,
) -> Vec<Vec<Complex64>> {
    let d = operators[0].dimension();
    let entries: Vec<Vec<Complex64>> = operators
        .iter()
        .map(|o| {
            (0..d)
                .flat_map(|r| (0..d).map(move |c| o.entry(r, c).expect("inside target block")))
                .collect()
        })
        .collect();
    let target = RealCorep {
        dimension: d,
        matrices: entries
            .iter()
            .map(|v| v.iter().map(|z| z.re).collect())
            .collect(),
        antiunitary: operators
            .iter()
            .map(ThreeDimensionalTargetOperator::is_antiunitary)
            .collect(),
    };
    let imaginary: Vec<Vec<f64>> = entries
        .iter()
        .map(|v| v.iter().map(|z| z.im).collect())
        .collect();
    let (rows, width) = hom_rows_with_imaginary_target(source, &target, Some(&imaginary));
    copies_from_basis(
        source,
        &target,
        &imaginary,
        &real_kernel(&rows, width),
        unitary_pair,
    )
}

struct Candidate {
    kind: ThreeDimensionalTargetKind,
    storage: Vec<ThreeDimensionalHStorage>,
    gauges: Vec<Vec<ThreeDimensionalTargetOperator>>,
}
fn real_operators(table: TargetTable) -> Vec<ThreeDimensionalTargetOperator> {
    table
        .operators
        .into_iter()
        .map(ThreeDimensionalTargetOperator::Real)
        .collect()
}
fn candidate(
    corep: &MagneticCorepSummary,
    context: &TargetContext,
) -> Result<Candidate, ThreeDimensionalSubductionError> {
    let sources = &corep.source_irreps;
    let anti = context.classes.iter().position(|g| g.time_reversal);
    let fail = |s: &str| failure(context.uni, &corep.label, s);
    if corep.corep_type == CorepType::A
        && sources.len() == 1
        && usize::from(sources[0].dim) == corep.dim
    {
        if corep.dim == 1 {
            let h = resolve_h_row(&sources[0], context.uni)?;
            let values = unitary_values(&h, corep, context)?;
            let operators = scalar_table(corep, context, &values, anti)?;
            return Ok(Candidate {
                kind: ThreeDimensionalTargetKind::TypeA,
                storage: vec![ThreeDimensionalHStorage::Scalar(h.storage)],
                gauges: vec![
                    operators
                        .into_iter()
                        .map(ThreeDimensionalTargetOperator::Scalar)
                        .collect(),
                ],
            });
        }
        if matches!(corep.dim, 2 | 3) {
            let tables = build_target_tables(corep, context, context.uni)?;
            // Resolution and the whole matrix record have just been certified.
            let record = query::irreps_of(sources[0].sg)
                .iter()
                .find(|r| r.ml == sources[0].ml)
                .expect("certified H record");
            return Ok(Candidate {
                kind: ThreeDimensionalTargetKind::TypeA,
                storage: vec![ThreeDimensionalHStorage::Matrix {
                    record: record.labels(),
                }],
                gauges: tables.into_iter().map(real_operators).collect(),
            });
        }
    }
    if corep.corep_type == CorepType::B
        && corep.dim == 2
        && sources.len() == 1
        && sources[0].dim == 1
    {
        let h = resolve_h_row(&sources[0], context.uni)?;
        let values = bind_h_row(&h, &sources[0], &corep.label, context)?;
        let a = anti.ok_or_else(|| fail("type B has no antiunitary representative"))?;
        let table = type_b_targets::type_b_table(corep, context, &values, a)?.ok_or_else(|| {
            fail("positive-overlap type B requires nonreal H outside the scalar-real domain")
        })?;
        return Ok(Candidate {
            kind: ThreeDimensionalTargetKind::TypeB,
            storage: vec![ThreeDimensionalHStorage::Scalar(h.storage)],
            gauges: vec![real_operators(table)],
        });
    }
    if corep.dim == 2
        && sources.len() == 2
        && sources.iter().all(|s| s.dim == 1)
        && (corep.corep_type == CorepType::C
            || (corep.corep_type == CorepType::A && anti.is_none()))
    {
        let h = [
            resolve_h_row(&sources[0], context.uni)?,
            resolve_h_row(&sources[1], context.uni)?,
        ];
        let values = [
            bind_h_row(&h[0], &sources[0], &corep.label, context)?,
            bind_h_row(&h[1], &sources[1], &corep.label, context)?,
        ];
        let conjugate = conjugate_pairs::is_conjugate_pair(corep, context, &values)?;
        let table = if conjugate {
            conjugate_pairs::realified_scalar_table(corep, context, &values[0], anti)?
        } else {
            if corep.corep_type != CorepType::C
                || anti.is_none()
                || !values.iter().all(|v| {
                    v.iter().enumerate().all(|(i, z)| {
                        context.classes[i].time_reversal
                            || (close(z.im, 0.0) && close(z.re.abs(), 1.0))
                    })
                })
            {
                return Err(fail(
                    "positive-overlap pair is neither real scalar C nor a nonreal conjugate pair",
                ));
            }
            let pair =
                screen_type_c_pair(corep, context, context.uni)?.map_err(|reason| fail(&reason))?;
            build_type_c_table(
                &pair,
                context,
                context.uni,
                &corep.label,
                anti.expect("checked anti representative"),
            )?
        };
        let [first, second] = h;
        return Ok(Candidate {
            kind: if anti.is_some() {
                ThreeDimensionalTargetKind::TypeC
            } else {
                ThreeDimensionalTargetKind::UnitaryPair
            },
            storage: vec![
                ThreeDimensionalHStorage::Scalar(first.storage),
                ThreeDimensionalHStorage::Scalar(second.storage),
            ],
            gauges: vec![real_operators(table)],
        });
    }
    Err(fail(
        "positive-overlap target lies outside the bounded three-dimensional target families",
    ))
}

fn verify_complete(
    model: &ThreeDimensionalGammaModel,
    pieces: &[ThreeDimensionalSubductionPiece],
) -> Result<(), ThreeDimensionalSubductionError> {
    let fail = |s| failure(model.uni(), "<combined>", s);
    let columns: Vec<Vec<Complex64>> = pieces
        .iter()
        .flat_map(|p| {
            p.copies.iter().flat_map(move |u| {
                (0..p.dimension).map(move |c| (0..3).map(|r| u[r * p.dimension + c]).collect())
            })
        })
        .collect();
    if columns.len() != 3 || !verify_column_isometry_for_source(&columns, 3) {
        return Err(fail(
            "copies do not form one orthonormal complete 3D source (Gram and projector)",
        ));
    }
    let mut vectors = Vec::new();
    for j in 0..3 {
        for z in [Complex64::new(1.0, 0.0), Complex64::new(0.0, 1.0)] {
            let mut v = [Complex64::new(0.0, 0.0); 3];
            v[j] = z;
            vectors.push(v);
        }
    }
    vectors.push([
        Complex64::new(0.25, 0.5),
        Complex64::new(-0.75, 0.125),
        Complex64::new(0.375, -0.625),
    ]);
    for (g, source) in model.source_operators().iter().enumerate() {
        for &v in &vectors {
            let expected = source.apply(v);
            let mut actual = [Complex64::new(0.0, 0.0); 3];
            for p in pieces {
                for u in &p.copies {
                    let coefficients: Vec<_> = (0..p.dimension)
                        .map(|c| (0..3).map(|r| u[r * p.dimension + c].conj() * v[r]).sum())
                        .collect();
                    let image = p.operators[g]
                        .apply(&coefficients)
                        .expect("target-sized vector");
                    for r in 0..3 {
                        actual[r] += (0..p.dimension)
                            .map(|c| u[r * p.dimension + c] * image[c])
                            .sum::<Complex64>();
                    }
                }
            }
            if actual
                .iter()
                .zip(expected)
                .any(|(a, b)| !(a - b).norm().is_finite() || (a - b).norm() > 1e-8)
            {
                return Err(fail("full complex source action does not reconstruct"));
            }
        }
    }
    Ok(())
}

/// Decompose a certified 3D source in one mixed collector. Character overlaps only
/// exclude orthogonal candidates; every retained copy must solve Hom and the final
/// combined Gram/projector and complex reconstruction must cover all three dimensions.
pub fn subduce_three_dimensional_gamma(
    model: &ThreeDimensionalGammaModel,
) -> Result<ThreeDimensionalSubduction<'_>, ThreeDimensionalSubductionError> {
    let (gamma, context) = prepare_gamma_targets(
        model.uni(),
        model.parent_sg(),
        3,
        model.classes(),
        model.translation_lattice(),
    )?;
    let source = RealCorep {
        dimension: 3,
        matrices: model
            .source_operators()
            .iter()
            .map(|o| o.entries().to_vec())
            .collect(),
        antiunitary: model
            .source_operators()
            .iter()
            .map(ThreeDimensionalSourceOperator::is_antiunitary)
            .collect(),
    };
    let mut pieces = Vec::new();
    for corep in &gamma.coreps {
        if corep.dim > 3 || corep.source_irreps.iter().any(|s| s.spinor) {
            continue;
        }
        if !has_unitary_overlap_with_trace(
            |i| matrix_trace(&source.matrices[i], 3),
            corep,
            &context,
        )? {
            continue;
        }
        let candidate = candidate(corep, &context)?;
        let mut matched = None;
        for operators in candidate.gauges {
            let copies = matching_three_dimensional_copies(
                &source,
                &operators,
                candidate.kind == ThreeDimensionalTargetKind::UnitaryPair,
            );
            if !copies.is_empty() {
                matched = Some((copies, operators));
                break;
            }
        }
        let (copies, operators) = matched.ok_or_else(|| {
            failure(
                model.uni(),
                &corep.label,
                "positive-overlap target has no complete isometric Hom copy",
            )
        })?;
        pieces.push(ThreeDimensionalSubductionPiece {
            label: corep.label.clone(),
            kind: candidate.kind,
            h_sources: corep.source_irreps.clone(),
            storage: candidate.storage,
            dimension: corep.dim,
            copies,
            operators,
        });
    }
    verify_complete(model, &pieces)?;
    Ok(ThreeDimensionalSubduction { model, pieces })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn model(cond: &str) -> ThreeDimensionalGammaModel {
        verify_three_dimensional_gamma_input(
            &MatrixGammaRequest::new(
                SourceIdentity {
                    parent_sg: 195,
                    label: "GM4",
                    convention: LabelConvention::Cdml,
                },
                RecordIdentity {
                    parent_sg: 195,
                    condensing_label: cond,
                    convention: LabelConvention::Cdml,
                    direction_index: 0,
                },
                Some(RealSign::Plus),
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
    fn whole_source_gate_rejects_duplicate_missing_and_wrong_anti_phase_copies() {
        let m = model("GM4");
        assert!(m.classes().iter().any(|g| g.time_reversal));
        let r = subduce_three_dimensional_gamma(&m).unwrap();
        let mut duplicate = r.pieces.clone();
        let p = duplicate.iter_mut().find(|p| p.copies.len() == 2).unwrap();
        p.copies[1] = p.copies[0].clone();
        assert!(
            verify_complete(&m, &duplicate)
                .unwrap_err()
                .to_string()
                .contains("Gram and projector")
        );
        let mut missing = r.pieces.clone();
        missing.pop();
        assert!(verify_complete(&m, &missing).is_err());
        let mut wrong_phase = r.pieces.clone();
        for z in &mut wrong_phase[0].copies[0] {
            *z *= Complex64::new(0.0, 1.0);
        }
        // The phase preserves the Gram and projector but not antiunitary sewing.
        assert!(
            verify_complete(&m, &wrong_phase)
                .unwrap_err()
                .to_string()
                .contains("complex source action")
        );
    }

    #[test]
    fn three_dimensional_targets_keep_positive_overlap_data_failures() {
        use std::error::Error;
        let m = model("GM1");
        let (gamma, context) = prepare_gamma_targets(
            m.uni(),
            m.parent_sg(),
            3,
            m.classes(),
            m.translation_lattice(),
        )
        .unwrap();
        let corep = gamma
            .coreps
            .iter()
            .find(|c| c.dim == 3 && c.source_irreps[0].ml == "GM4")
            .unwrap();
        let trace = |i: usize| matrix_trace(&m.source_operators()[i].entries(), 3);
        assert!(has_unitary_overlap_with_trace(trace, corep, &context).unwrap());
        let mut missing = corep.clone();
        missing.source_irreps[0].ml = "MISSING";
        assert!(has_unitary_overlap_with_trace(trace, &missing, &context).unwrap());
        let error = match candidate(&missing, &context) {
            Err(e) => e,
            Ok(_) => panic!("missing H data accepted"),
        };
        assert!(
            matches!(error,ThreeDimensionalSubductionError::Target(ref e) if matches!(**e,MatrixSubductionError::TargetSourceUnresolved{..}))
        );
        assert!(
            error
                .source()
                .unwrap()
                .downcast_ref::<MatrixSubductionError>()
                .is_some()
        );
        let mut bad = corep.clone();
        bad.source_irreps[0].sg = 1;
        assert!(has_unitary_overlap_with_trace(trace, &bad, &context).is_err());
        let column = context
            .class_columns
            .iter()
            .find(|c| !c.is_empty())
            .unwrap()[0];
        for value in [
            None,
            Some(Complex64::new(f64::NAN, 0.0)),
            Some(Complex64::new(f64::INFINITY, 0.0)),
        ] {
            let mut bad = corep.clone();
            bad.characters[column] = value;
            assert!(has_unitary_overlap_with_trace(trace, &bad, &context).is_err());
        }
    }

    fn matrix_identity(n: usize) -> Vec<f64> {
        (0..n * n).map(|i| f64::from(i / n == i % n)).collect()
    }
    fn encode(u: &[Complex64]) -> Vec<f64> {
        u.iter().flat_map(|z| [z.re, z.im]).collect()
    }
    fn zero_imaginary(t: &RealCorep) -> Vec<Vec<f64>> {
        vec![vec![0.0; t.dimension * t.dimension]; t.matrices.len()]
    }
    fn embed(t: &[f64]) -> Vec<f64> {
        vec![1.0, 0.0, 0.0, 0.0, t[0], t[1], 0.0, t[2], t[3]]
    }

    #[test]
    fn repeated_scalar_projects_nonorthogonal_basis_and_preserves_anti_equations() {
        for anti in [false, true] {
            let s = RealCorep {
                dimension: 3,
                matrices: vec![
                    matrix_identity(3),
                    matrix_identity(3)
                        .iter()
                        .map(|x| if anti { -x } else { *x })
                        .collect(),
                ],
                antiunitary: vec![false, anti],
            };
            let t = RealCorep {
                dimension: 1,
                matrices: vec![vec![1.0]; 2],
                antiunitary: vec![false, anti],
            };
            let phase = if anti {
                Complex64::new(0.0, 1.0)
            } else {
                Complex64::new(1.0, 0.0)
            };
            let basis: Vec<_> = [[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 0.0]]
                .iter()
                .map(|v| encode(&v.map(|x| phase * x)))
                .collect();
            let copies = copies_from_basis(&s, &t, &zero_imaginary(&t), &basis, false);
            assert_eq!(copies.len(), 3);
            assert!(verify_column_isometry_for_source(&copies, 3));
            for u in &copies {
                assert!(intertwines(&s, &t, &zero_imaginary(&t), u));
            }
            let (rows, width) = hom_rows(&s, &t);
            let kernel = real_kernel(&rows, width);
            assert_eq!(kernel.len(), if anti { 3 } else { 6 });
            assert_eq!(
                copies_from_basis(&s, &t, &zero_imaginary(&t), &kernel, false).len(),
                3
            );
            if anti {
                let wrong = encode(&[Complex64::new(1.0, 0.0); 3]);
                assert!(copies_from_basis(&s, &t, &zero_imaginary(&t), &[wrong], false).is_empty());
            }
        }
    }

    #[test]
    fn unitary_packet_extracts_real_hom_from_rank_one_complex_basis() {
        let c = -0.5;
        let s = 3.0_f64.sqrt() / 2.0;
        let t = RealCorep {
            dimension: 2,
            matrices: vec![matrix_identity(2), vec![c, -s, s, c]],
            antiunitary: vec![false; 2],
        };
        let source = RealCorep {
            dimension: 3,
            matrices: t.matrices.iter().map(|m| embed(m)).collect(),
            antiunitary: vec![false; 2],
        };
        let z = Complex64::new(0.0, 0.0);
        let one = Complex64::new(1.0, 0.0);
        let i = Complex64::new(0.0, 1.0);
        let singular = vec![z, z, one, -i, i, one];
        assert!(intertwines(&source, &t, &zero_imaginary(&t), &singular));
        assert!(normalize_intertwiner_for_source(&singular, 2, 3).is_none());
        let basis = [
            encode(&singular),
            encode(&singular.iter().map(|z| z.conj()).collect::<Vec<_>>()),
        ];
        assert!(copies_from_basis(&source, &t, &zero_imaginary(&t), &basis, false).is_empty());
        let copies = copies_from_basis(&source, &t, &zero_imaginary(&t), &basis, true);
        assert_eq!(copies.len(), 1);
        assert!(intertwines(&source, &t, &zero_imaginary(&t), &copies[0]));
        assert!(normalize_intertwiner_for_source(&copies[0], 2, 3).is_some());
    }

    #[test]
    fn quaternionic_real_hom_dimension_four_is_one_two_dimensional_copy() {
        let j = vec![0.0, -1.0, 1.0, 0.0];
        let t = RealCorep {
            dimension: 2,
            matrices: vec![
                matrix_identity(2),
                j.clone(),
                vec![-1.0, 0.0, 0.0, -1.0],
                j.iter().map(|x| -x).collect(),
            ],
            antiunitary: vec![false, true, false, true],
        };
        let q = vec![
            6.0 / 7.0,
            -2.0 / 7.0,
            -3.0 / 7.0,
            -2.0 / 7.0,
            3.0 / 7.0,
            -6.0 / 7.0,
            -3.0 / 7.0,
            -6.0 / 7.0,
            -2.0 / 7.0,
        ];
        let source = RealCorep {
            dimension: 3,
            matrices: t
                .matrices
                .iter()
                .map(|m| matrix_product(&matrix_product(&q, &embed(m), 3), &q, 3))
                .collect(),
            antiunitary: t.antiunitary.clone(),
        };
        let (rows, width) = hom_rows(&source, &t);
        let kernel = real_kernel(&rows, width);
        assert_eq!(kernel.len(), 4);
        let copies = copies_from_basis(&source, &t, &zero_imaginary(&t), &kernel, false);
        assert_eq!(copies.len(), 1);
        let a = RealCorep {
            dimension: 1,
            matrices: vec![vec![1.0]; 4],
            antiunitary: t.antiunitary.clone(),
        };
        let (rows, width) = hom_rows(&source, &a);
        let scalar = copies_from_basis(
            &source,
            &a,
            &zero_imaginary(&a),
            &real_kernel(&rows, width),
            false,
        );
        assert_eq!(scalar.len(), 1);
        let mut columns = scalar;
        for c in 0..2 {
            columns.push((0..3).map(|r| copies[0][2 * r + c]).collect());
        }
        assert!(verify_column_isometry_for_source(&columns, 3));
        let wrong = RealCorep {
            dimension: 3,
            matrices: vec![matrix_identity(3); 2],
            antiunitary: vec![false, true],
        };
        let target = RealCorep {
            dimension: 2,
            matrices: vec![matrix_identity(2), j],
            antiunitary: vec![false, true],
        };
        let (rows, width) = hom_rows(&wrong, &target);
        assert!(real_kernel(&rows, width).is_empty());
    }
}
