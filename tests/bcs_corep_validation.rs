//! Regression tests against Bilbao Crystallographic Server corep metadata.
//!
//! Reference magnetic group: BNS 128.406 (`P4'/m'nc'`, UNI 1066).
//! Its unitary subgroup is SG 118 (`P-4n2`).  At Z, BCS lists four magnetic
//! corepresentations with dimensions 2, 2, 2, and 4.

use cryspglib::irrep::magnetic_summary::{MagneticIrrepError, magnetic_irrep_summary_by_bns};

#[test]
fn bcs_sg128_406_z_has_official_dimensions() {
    let error = magnetic_irrep_summary_by_bns("128.406")
        .expect_err("compound corepresentation must fail closed in summaries");
    match error {
        MagneticIrrepError::CorepComputationFailed {
            uni,
            sg,
            k_label,
            source_irrep,
            reason,
        } => {
            assert_eq!(uni, 1066);
            assert_eq!(sg, 118);
            assert_eq!(k_label, "Z");
            assert_eq!(source_irrep, "Z1Z4");
            assert!(reason.contains("constituent-orbit Wigner analysis"));
            assert!(reason.contains("physical aggregate block trace"));
        }
        other => panic!("unexpected BNS 128.406 error: {other:?}"),
    }
}
