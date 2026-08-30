//! Regression tests against Bilbao Crystallographic Server corep metadata.
//!
//! Reference magnetic group: BNS 128.406 (`P4'/m'nc'`, UNI 1066).
//! Its unitary subgroup is SG 118 (`P-4n2`).  At Z, BCS lists four magnetic
//! corepresentations with dimensions 2, 2, 2, and 4.

use cryspglib::irrep::magnetic_summary::magnetic_irrep_summary_by_bns;

#[test]
fn bcs_sg128_406_z_has_official_dimensions() {
    let summary = magnetic_irrep_summary_by_bns("128.406")
        .expect("the strict complex summary must include the compound branches");
    let z = summary
        .kpoints
        .iter()
        .find(|point| point.label == "Z")
        .expect("missing Z point");
    let mut dimensions = z.coreps.iter().map(|corep| corep.dim).collect::<Vec<_>>();
    dimensions.sort_unstable();
    assert_eq!(dimensions, vec![2, 2, 2, 4]);
    assert!(
        z.coreps
            .iter()
            .flat_map(|corep| corep.characters.iter().flatten())
            .any(|character| character.im.abs() > 1.0),
        "the official compound table contains genuinely complex columns"
    );
}
