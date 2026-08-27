//! Regression tests against Bilbao Crystallographic Server corep metadata.
//!
//! Reference magnetic group: BNS 128.406 (`P4'/m'nc'`, UNI 1066).
//! Its unitary subgroup is SG 118 (`P-4n2`).  At Z, BCS lists four magnetic
//! corepresentations with dimensions 2, 2, 2, and 4.

use cryspglib::irrep::corep::{CharacterCompleteness, CorepType};
use cryspglib::irrep::magnetic_summary::{
    format_magnetic_character_table, format_magnetic_character_table_by_class,
    magnetic_irrep_summary_by_bns,
};

#[test]
fn bcs_sg128_406_z_has_official_dimensions() {
    let summary = magnetic_irrep_summary_by_bns("128.406")
        .expect("BNS 128.406 must have a valid magnetic-irrep summary");
    assert_eq!(summary.uni, 1066);
    assert_eq!(summary.parent_sg, 128);
    assert_eq!(summary.unitary_sg, 118);

    let z = summary
        .kpoints
        .iter()
        .find(|kpoint| kpoint.label == "Z")
        .expect("BNS 128.406 must contain the Z high-symmetry point");
    assert_eq!(z.coords, (0, 0, 1, 2));
    assert_eq!(z.little_group_order, 16);
    assert_eq!(z.unitary_order, 8);
    assert_eq!(z.antiunitary_order, 8);
    assert_eq!(z.operations.len(), 16);

    let identity = z
        .operations
        .iter()
        .position(|operation| {
            !operation.time_reversal
                && operation.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                && operation
                    .translation
                    .iter()
                    .all(|value| (value - value.round()).abs() < 1e-8)
        })
        .expect("magnetic little group must contain identity");

    let expected = [
        ("Z1Z4", CorepType::C, 2usize),
        ("Z2Z3", CorepType::C, 2usize),
        ("Z5", CorepType::A, 2usize),
        ("Z6 + Z7", CorepType::C, 4usize),
    ];
    assert_eq!(z.coreps.len(), expected.len());
    for (label, corep_type, dimension) in expected {
        let corep = z
            .coreps
            .iter()
            .find(|corep| corep.label == label)
            .unwrap_or_else(|| panic!("missing Z-point corep {label}"));
        assert_eq!(corep.corep_type, corep_type, "{label}");
        assert_eq!(corep.dim, dimension, "{label}");
        if corep_type == CorepType::A && dimension > 1 {
            assert_eq!(
                corep.completeness,
                CharacterCompleteness::TypeAAntiunitaryPending { count: 8 }
            );
        } else {
            assert_eq!(corep.completeness, CharacterCompleteness::Complete);
        }
        assert_eq!(corep.characters.len(), z.operations.len());
        assert_eq!(corep.timerev.len(), z.operations.len());
        assert!(corep.characters.iter().all(|value| value.is_finite()));
        assert!((corep.characters[identity] - dimension as f64).abs() < 1e-8);
    }

    // The public formatter is a formal table: every operation is a column,
    // followed by an explicit Seitz-operation legend.  No six-column preview
    // or ellipsis is allowed.
    let operation_table = format_magnetic_character_table(z);
    assert!(operation_table.contains("| g16 |"));
    assert!(operation_table.contains("Seitz operation (data-Hall frame)"));
    assert!(!operation_table.contains("..."));

    let class_table = format_magnetic_character_table_by_class(z);
    assert!(class_table.contains("member operation columns"));
    assert!(class_table.contains("| C1"));
}
