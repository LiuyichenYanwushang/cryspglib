//! Test-only P1 source audit. This is not a wider production input contract.
use super::*;
use crate::irrep::types::CompoundSelectedArmCharacter;
use std::collections::BTreeMap;

fn gamma(ir: &IrrepRecord) -> bool {
    let k = ir.k_vector();
    k.numerators == [0; 3] && k.denominator == 1
}

fn row(ir: &IrrepRecord) -> CharacterRow {
    if ir.compound_metadata().is_some() {
        ir.compound_selected_arm_view()
            .unwrap()
            .block_trace()
            .clone()
    } else {
        ir.ordinary_scalar_selected_arm_block_trace().unwrap()
    }
}

// Historical P1 audit of all stored shapes, including compound carriers.
// Production verifiers separately enforce their explicit source admission gates.
fn audit_table(ir: &IrrepRecord) -> (CertifiedTable, usize, [Rat; 3]) {
    assert!(!ir.spinor && gamma(ir));
    let row = row(ir);
    let d = usize::from(ir.dim);
    assert_eq!(row.dimension(), d);
    assert_eq!(ir.pir_rotations().len(), row.len() * 9);
    assert_eq!(ir.pir_translations().len(), row.len() * 3);
    assert_eq!(ir.characters().len(), row.len());
    assert_eq!(ir.matrices().len(), row.len() * d * d);
    let mut operations = Vec::new();
    let mut lattice = Vec::new();
    let mut matrices = Vec::new();
    for (i, op) in row.operations().iter().enumerate() {
        assert_eq!(
            op.rotation.as_slice(),
            &ir.pir_rotations()[i * 9..(i + 1) * 9]
        );
        assert_eq!(
            op.translation.as_slice(),
            &ir.pir_translations()[i * 3..(i + 1) * 3]
        );
        let rotation = std::array::from_fn(|r| std::array::from_fn(|c| op.rotation[r * 3 + c]));
        let translation = op.translation.map(|v| exact_fraction(v).unwrap());
        if rotation == identity_rotation() && !lattice.contains(&translation) {
            lattice.push(translation);
        }
        let matrix = ir.matrices()[i * d * d..(i + 1) * d * d].to_vec();
        assert!(matrix.iter().all(|v| v.is_finite()));
        assert!(is_orthogonal(&matrix, d));
        let trace = matrix_trace(&matrix, d);
        assert!(close(trace, ir.characters()[i]));
        assert!((row.values()[i] - Complex64::new(trace, 0.0)).norm() < GRID_TOLERANCE);
        if rotation == identity_rotation() && translation == [Rat::ZERO; 3] {
            assert_eq!(
                matrix,
                (0..d * d)
                    .map(|j| f64::from(j % (d + 1) == 0))
                    .collect::<Vec<_>>()
            );
        }
        operations.push((rotation, translation));
        matrices.push(matrix);
    }
    assert!(operations.contains(&(identity_rotation(), [Rat::ZERO; 3])));
    check_stored_matrix_law(&operations, &matrices, d, &lattice).unwrap();
    let commutant = real_commutant_dimension(&matrices, d);
    // Source-only audit: 0 is error context, not a claimed magnetic UNI.
    let shift = ParentFrame::of(ir.sg)
        .unwrap()
        .align_source_operations(ir.sg, 0, &mut operations)
        .unwrap();
    (
        CertifiedTable {
            dimension: d,
            operations,
            matrices,
        },
        commutant,
        shift,
    )
}

#[test]
fn all_nonspinor_gamma_sources_have_complete_stored_action_and_frame_evidence() {
    let mut counts = BTreeMap::new();
    println!(
        "\nSOURCE\tsg\tlabel\tdimension\tfamily\trows\treal_commutant\torigin\tbc\tconstituents\tsource_identity"
    );
    for sg in 1..=230 {
        for ir in query::irreps_of(sg)
            .iter()
            .filter(|ir| !ir.spinor && gamma(ir))
        {
            let (table, commutant, shift) = audit_table(ir);
            let family = if ir.compound_metadata().is_some() {
                "compound"
            } else {
                "ordinary"
            };
            let constituents = if let Some(meta) = ir.compound_metadata() {
                assert_eq!(ir.dim, 2);
                assert_eq!(commutant, 2);
                let view = ir.compound_selected_arm_view().unwrap();
                let CompoundSelectedArmCharacter::DistinctComponentSum { first, second, .. } = view
                else {
                    panic!("Gamma inventory changed: {} {}", sg, ir.ml)
                };
                assert_eq!((first.row.dimension(), second.row.dimension()), (1, 1));
                assert_ne!(first.irnumber, second.irnumber);
                assert!(
                    first
                        .row
                        .values()
                        .iter()
                        .zip(second.row.values())
                        .all(|(a, b)| (*a - b.conj()).norm() < GRID_TOLERANCE)
                );
                assert!(
                    first
                        .row
                        .values()
                        .iter()
                        .any(|a| a.im.abs() > GRID_TOLERANCE)
                );
                format!(
                    "{}:{}+{}:{}",
                    meta.cir_labels[0],
                    meta.cir_irnumbers[0],
                    meta.cir_labels[1],
                    meta.cir_irnumbers[1]
                )
            } else {
                assert_eq!(commutant, 1);
                "-".into()
            };
            *counts
                .entry((family, ir.dim, shift != [Rat::ZERO; 3]))
                .or_insert(0usize) += 1;
            println!(
                "SOURCE\t{sg}\t{}\t{}\t{family}\t{}\t{commutant}\t{shift:?}\t{}\t{constituents}\t{:?}",
                ir.ml,
                ir.dim,
                table.operations.len(),
                ir.labels().bc.as_deref().unwrap_or("-"),
                ir.source_identity()
            );
        }
    }
    let count = |family, dim| {
        counts
            .iter()
            .filter(|((f, d, _), _)| *f == family && *d == dim)
            .map(|(_, n)| n)
            .sum::<usize>()
    };
    assert_eq!(count("ordinary", 1), 980);
    assert_eq!(count("ordinary", 2), 177);
    assert_eq!(count("ordinary", 3), 87);
    assert_eq!(count("compound", 2), 69);
    assert_eq!(counts.values().sum::<usize>(), 1313);
    println!("SOURCE_TOTAL {counts:?}");
}

#[test]
fn three_dimensional_source_stays_outside_old_input_and_trace_is_not_a_law_proof() {
    let ir = query::irreps_of(195)
        .iter()
        .find(|r| r.ml == "GM4")
        .unwrap();
    let (mut table, commutant, _) = audit_table(ir);
    assert_eq!((table.dimension, commutant), (3, 1));
    assert!(matches!(
        resolve_matrix_source(SourceIdentity {
            parent_sg: 195,
            label: ir.ml,
            convention: LabelConvention::Cdml
        }),
        Err(MatrixSubductionError::SourceDimensionNotTwo { dimension: 3, .. })
    ));
    let admitted_table = certify_irrep_matrices(ir, 3).unwrap();
    assert_eq!(admitted_table.operations, table.operations);
    assert_eq!(admitted_table.matrices, table.matrices);
    let (a, b) = (0..table.matrices.len())
        .flat_map(|a| (a + 1..table.matrices.len()).map(move |b| (a, b)))
        .find(|&(a, b)| {
            table.matrices[a] != table.matrices[b]
                && close(
                    matrix_trace(&table.matrices[a], 3),
                    matrix_trace(&table.matrices[b], 3),
                )
        })
        .unwrap();
    table.matrices.swap(a, b);
    assert!(
        table
            .matrices
            .iter()
            .enumerate()
            .all(|(i, m)| is_orthogonal(m, 3) && close(matrix_trace(m, 3), ir.characters()[i]))
    );
    assert!(matches!(
        check_stored_matrix_law(
            &table.operations,
            &table.matrices,
            3,
            &ParentFrame::of(195).unwrap().lattice
        ),
        Err(MatrixCertificationProblem::LawViolated { .. })
    ));
}

#[test]
fn compound_source_exposes_a_real_scalar_h_type_b_witness_without_widening_the_api() {
    let ir = query::irreps_of(75)
        .iter()
        .find(|r| r.ml == "GM3GM4")
        .unwrap();
    let (table, commutant, shift) = audit_table(ir);
    assert_eq!(commutant, 2);
    assert_eq!(shift, [Rat::ZERO; 3]);
    let record = query::irreps_of(75)
        .iter()
        .find(|r| r.ml == "GM2")
        .unwrap()
        .magnetic_subgroups()[0];
    assert_eq!(record.mag_sg, 663);
    let set = magnetic_operations(663).unwrap();
    magnetic_embedding::verify_group(&set).unwrap();
    magnetic_embedding::verify_against_database(&set).unwrap();
    let geometry = magnetic_embedding::geometry_of(75, &record).unwrap();
    magnetic_embedding::verify_embedding(
        &geometry,
        &set,
        crate::irrep::subduction::Mat3R::identity(),
        [Rat::ZERO; 3],
    )
    .unwrap();
    assert!(
        magnetic_embedding::realises_record_lattice_checked(
            &geometry,
            &set,
            &crate::irrep::subduction::Mat3R::identity()
        )
        .unwrap()
    );
    let index = |g: &MagneticOperation| {
        table
            .operations
            .iter()
            .position(|(r, t)| {
                *r == g.rotation
                    && differs_by_lattice(t, &g.translation, &ParentFrame::of(75).unwrap().lattice)
            })
            .unwrap()
    };
    let a = set.classes.iter().find(|g| g.time_reversal).unwrap();
    let a2 = compose(a, a).unwrap();
    assert_eq!(table.matrices[index(&a2)], vec![-1.0, 0.0, 0.0, -1.0]);
    // Canonical Type-B target: h -> chi(h) I, a0 h -> chi(h) J K, J²=-I.
    // These are independent of the source matrices. The real orthogonal U
    // below then checks every actual coloured operation for both eta values.
    let summary = magnetic_summary::magnetic_irrep_summary_by_uni_partial(663).unwrap();
    let gamma = summary
        .summary
        .kpoints
        .iter()
        .find(|k| k.coords == GAMMA)
        .unwrap();
    let b = gamma
        .coreps
        .iter()
        .find(|c| c.label == "GM2" && c.corep_type == CorepType::B)
        .unwrap();
    assert_eq!((b.dim, b.source_irreps.len()), (2, 1));
    let h_source = &b.source_irreps[0];
    assert_eq!(
        (h_source.sg, h_source.ml, h_source.dim, h_source.spinor),
        (3, "GM2", 1, false)
    );
    let owners: Vec<_> = query::irreps_of(h_source.sg)
        .iter()
        .filter(|r| r.ml == h_source.ml && !r.spinor)
        .collect();
    assert_eq!(owners.len(), 1);
    let h = owners[0]
        .ordinary_scalar_selected_arm_block_trace()
        .unwrap();
    for eta in [RealSign::Plus, RealSign::Minus] {
        let sign = f64::from(eta.value());
        let u = [1.0, 0.0, 0.0, sign * table.matrices[index(a)][2]];
        assert!(is_orthogonal(&u, 2));
        for col in gamma.operations.iter().filter(|g| !g.time_reversal) {
            let raw = &set.operations[col.magnetic_operation_index];
            let h_index = h
                .operations()
                .iter()
                .position(|op| {
                    op.rotation.as_slice() == col.rotation.as_flattened()
                        && op
                            .translation
                            .iter()
                            .zip(col.translation)
                            .all(|(x, y)| ((x - y) - (x - y).round()).abs() < GRID_TOLERANCE)
                })
                .unwrap();
            let chi = h.values()[h_index];
            assert!(close(chi.im, 0.0) && close(chi.re.abs(), 1.0));
            for (g, target) in [
                (*raw, vec![chi.re, 0.0, 0.0, chi.re]),
                (compose(a, raw).unwrap(), vec![0.0, -chi.re, chi.re, 0.0]),
            ] {
                let mut source = table.matrices[index(&g)].clone();
                if g.time_reversal {
                    source.iter_mut().for_each(|x| *x *= sign);
                }
                assert!(
                    matrix_product(&source, &u, 2)
                        .iter()
                        .zip(matrix_product(&u, &target, 2))
                        .all(|(x, y)| close(*x, y))
                );
            }
        }
        let request = MatrixGammaRequest::new(
            SourceIdentity {
                parent_sg: 75,
                label: ir.ml,
                convention: LabelConvention::Cdml,
            },
            RecordIdentity {
                parent_sg: 75,
                condensing_label: "GM2",
                convention: LabelConvention::Cdml,
                direction_index: 0,
            },
            Some(eta),
        );
        assert!(matches!(
            verify_matrix_gamma_input(&request),
            Err(MatrixSubductionError::SourceIsCompound { .. })
        ));
    }
}
