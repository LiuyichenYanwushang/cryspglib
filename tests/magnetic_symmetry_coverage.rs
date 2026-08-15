//! Exhaustive magnetic-symmetry database invariants.
//!
//! These tests deliberately cover every supported `(UNI, Hall)` pair rather
//! than only the first Hall setting of each UNI.  They validate the raw
//! magnetic Seitz groups independently of the higher-level identification and
//! corepresentation pipelines.

use std::collections::{BTreeMap, HashSet};

use cryspglib::mathfunc::{Mat3, Mat3I, mat_get_determinant_i3, mat_multiply_matrix_i3};
use cryspglib::msg_database::{
    ALTERNATIVE_TRANSFORMATIONS, MAGNETIC_SPACEGROUP_TYPES, MAGNETIC_SPACEGROUP_UNI_MAPPING,
    msgdb_get_spacegroup_operations, msgdb_get_std_transformations, msgdb_get_uni_candidates,
};
use cryspglib::{MagneticType, SymError, magnetic_spacegroup, msg_database, spg_database};

const TRANSLATION_DENOMINATOR: f64 = 12.0;
const TRANSLATION_TOLERANCE: f64 = 1e-8;
const IDENTITY_ROTATION: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

#[derive(Debug, Clone, Copy)]
struct TestOp {
    rotation: Mat3I,
    translation: [f64; 3],
    time_reversal: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct OpKey {
    rotation: [i32; 9],
    translation_twelfths: [i32; 3],
    time_reversal: bool,
}

#[derive(Default)]
struct Audit {
    counts: BTreeMap<&'static str, usize>,
    examples: BTreeMap<&'static str, Vec<String>>,
}

impl Audit {
    fn record(&mut self, category: &'static str, message: impl Into<String>) {
        *self.counts.entry(category).or_default() += 1;
        let examples = self.examples.entry(category).or_default();
        if examples.len() < 8 {
            examples.push(message.into());
        }
    }

    fn assert_clean(&self) {
        if self.counts.is_empty() {
            return;
        }

        self.report("magnetic symmetry audit failures");
        panic!(
            "magnetic symmetry audit found {} failure categories",
            self.counts.len()
        );
    }

    fn report(&self, heading: &str) {
        eprintln!("{heading}:");
        if self.counts.is_empty() {
            eprintln!("  none");
            return;
        }
        for (category, count) in &self.counts {
            eprintln!("  {category}: {count}");
            if let Some(examples) = self.examples.get(category) {
                for example in examples {
                    eprintln!("    {example}");
                }
            }
        }
    }
}

fn flatten_rotation(rotation: &Mat3I) -> [i32; 9] {
    [
        rotation[0][0],
        rotation[0][1],
        rotation[0][2],
        rotation[1][0],
        rotation[1][1],
        rotation[1][2],
        rotation[2][0],
        rotation[2][1],
        rotation[2][2],
    ]
}

fn quantize_translation(translation: &[f64; 3]) -> Option<[i32; 3]> {
    let mut result = [0; 3];
    for axis in 0..3 {
        if !translation[axis].is_finite() {
            return None;
        }
        let scaled = translation[axis] * TRANSLATION_DENOMINATOR;
        let rounded = scaled.round();
        if (scaled - rounded).abs() > TRANSLATION_TOLERANCE {
            return None;
        }
        result[axis] = (rounded as i32).rem_euclid(TRANSLATION_DENOMINATOR as i32);
    }
    Some(result)
}

fn op_key(op: &TestOp) -> Option<OpKey> {
    Some(OpKey {
        rotation: flatten_rotation(&op.rotation),
        translation_twelfths: quantize_translation(&op.translation)?,
        time_reversal: op.time_reversal,
    })
}

fn compose(left: &TestOp, right: &TestOp) -> TestOp {
    let rotation = mat_multiply_matrix_i3(&left.rotation, &right.rotation);
    let mut translation = left.translation;
    for (axis, coordinate) in translation.iter_mut().enumerate() {
        *coordinate += left.rotation[axis][0] as f64 * right.translation[0]
            + left.rotation[axis][1] as f64 * right.translation[1]
            + left.rotation[axis][2] as f64 * right.translation[2];
    }
    TestOp {
        rotation,
        translation,
        time_reversal: left.time_reversal ^ right.time_reversal,
    }
}

fn inverse_rotation(rotation: &Mat3I) -> Option<Mat3I> {
    let determinant = mat_get_determinant_i3(rotation);
    if determinant != 1 && determinant != -1 {
        return None;
    }

    let a = rotation[0][0];
    let b = rotation[0][1];
    let c = rotation[0][2];
    let d = rotation[1][0];
    let e = rotation[1][1];
    let f = rotation[1][2];
    let g = rotation[2][0];
    let h = rotation[2][1];
    let i = rotation[2][2];
    Some([
        [
            (e * i - f * h) / determinant,
            (c * h - b * i) / determinant,
            (b * f - c * e) / determinant,
        ],
        [
            (f * g - d * i) / determinant,
            (a * i - c * g) / determinant,
            (c * d - a * f) / determinant,
        ],
        [
            (d * h - e * g) / determinant,
            (b * g - a * h) / determinant,
            (a * e - b * d) / determinant,
        ],
    ])
}

fn inverse(op: &TestOp) -> Option<TestOp> {
    let rotation = inverse_rotation(&op.rotation)?;
    let mut translation = [0.0; 3];
    for axis in 0..3 {
        translation[axis] = -(rotation[axis][0] as f64 * op.translation[0]
            + rotation[axis][1] as f64 * op.translation[1]
            + rotation[axis][2] as f64 * op.translation[2]);
    }
    Some(TestOp {
        rotation,
        translation,
        time_reversal: op.time_reversal,
    })
}

fn is_zero_translation(translation: &[f64; 3]) -> bool {
    quantize_translation(translation) == Some([0, 0, 0])
}

fn rotation_multiset<'a>(rotations: impl Iterator<Item = &'a Mat3I>) -> BTreeMap<[i32; 9], usize> {
    let mut result = BTreeMap::new();
    for rotation in rotations {
        *result.entry(flatten_rotation(rotation)).or_default() += 1;
    }
    result
}

fn invariant_lattice(rotations: &[Mat3I]) -> Option<Mat3> {
    if rotations.is_empty() {
        return None;
    }

    // Average R^T R over the finite point group. The result is a positive
    // definite metric G satisfying R^T G R = G for every point operation.
    let mut metric = [[0.0; 3]; 3];
    for rotation in rotations {
        for row in 0..3 {
            for column in 0..3 {
                for rotation_row in rotation {
                    metric[row][column] +=
                        (rotation_row[row] * rotation_row[column]) as f64;
                }
            }
        }
    }
    let scale = rotations.len() as f64;
    for row in &mut metric {
        for value in row {
            *value /= scale;
        }
    }

    // Cholesky G = C C^T, then lattice A = C^T so A^T A = G.
    let c00 = metric[0][0].sqrt();
    if !c00.is_finite() || c00 <= 0.0 {
        return None;
    }
    let c10 = metric[1][0] / c00;
    let c20 = metric[2][0] / c00;
    let c11_squared = metric[1][1] - c10 * c10;
    if c11_squared <= 0.0 {
        return None;
    }
    let c11 = c11_squared.sqrt();
    let c21 = (metric[2][1] - c20 * c10) / c11;
    let c22_squared = metric[2][2] - c20 * c20 - c21 * c21;
    if c22_squared <= 0.0 {
        return None;
    }
    let c22 = c22_squared.sqrt();

    Some([[c00, c10, c20], [0.0, c11, c21], [0.0, 0.0, c22]])
}

fn identification_error_category(error: SymError) -> &'static str {
    match error {
        SymError::Success => "error_success",
        SymError::SpacegroupSearchFailed => "error_spacegroup_search",
        SymError::CellStandardizationFailed => "error_cell_standardization",
        SymError::SymmetryOperationSearchFailed => "error_symmetry_operation_search",
        SymError::AtomsTooClose => "error_atoms_too_close",
        SymError::PointgroupNotFound => "error_pointgroup_not_found",
        SymError::NiggliFailed => "error_niggli",
        SymError::DelaunayFailed => "error_delaunay",
        SymError::ArraySizeShortage => "error_array_size",
        SymError::InvalidInput => "error_invalid_input",
        SymError::MathFailed => "error_math",
        SymError::MagneticOpGenerationFailed => "error_magnetic_op_generation",
        SymError::MagneticReferenceGroupFailed => "error_magnetic_reference_group",
        SymError::MagneticFallbackReferenceFailed => "error_magnetic_fallback_reference",
        SymError::MagneticUniCandidatesNotFound => "error_magnetic_uni_candidates",
        SymError::MagneticUniMatchFailed => "error_magnetic_uni_match",
        SymError::MagneticPrimitiveLatticeFailed => "error_magnetic_primitive_lattice",
        SymError::MagneticUniAmbiguous => "error_magnetic_uni_ambiguous",
    }
}

#[test]
fn all_magnetic_database_metadata_is_complete_and_unique() {
    let mut bns_labels = HashSet::new();
    let mut og_labels = HashSet::new();
    let mut litvin_numbers = HashSet::new();

    for uni in 1usize..=1651 {
        let metadata = &MAGNETIC_SPACEGROUP_TYPES[uni];
        assert_eq!(
            metadata.uni_number, uni,
            "UNI {uni}: metadata index mismatch"
        );
        assert!(
            (1..=230).contains(&metadata.number),
            "UNI {uni}: invalid parent SG {}",
            metadata.number
        );
        assert!(
            metadata.type_ != MagneticType::NonMagnetic,
            "UNI {uni}: invalid non-magnetic type"
        );
        assert!(
            !metadata.bns_number.is_empty() && bns_labels.insert(metadata.bns_number),
            "UNI {uni}: empty or duplicate BNS label {}",
            metadata.bns_number
        );
        assert!(
            !metadata.og_number.is_empty() && og_labels.insert(metadata.og_number),
            "UNI {uni}: empty or duplicate OG label {}",
            metadata.og_number
        );
        assert!(
            (1..=1651).contains(&metadata.litvin_number)
                && litvin_numbers.insert(metadata.litvin_number),
            "UNI {uni}: invalid or duplicate Litvin number {}",
            metadata.litvin_number
        );

        let [num_halls, first_hall] = MAGNETIC_SPACEGROUP_UNI_MAPPING[uni];
        assert!(num_halls > 0, "UNI {uni}: no Hall settings");
        let last_hall = first_hall + num_halls - 1;
        assert!(
            first_hall >= 1 && last_hall <= 530,
            "UNI {uni}: invalid Hall range {first_hall}..={last_hall}"
        );
    }

    assert_eq!(bns_labels.len(), 1651);
    assert_eq!(og_labels.len(), 1651);
    assert_eq!(litvin_numbers.len(), 1651);
}

#[test]
fn all_magnetic_database_operations_form_expected_groups() {
    let mut audit = Audit::default();
    let mut hall_pair_count = 0usize;

    for (uni, &[num_halls, first_hall]) in MAGNETIC_SPACEGROUP_UNI_MAPPING
        .iter()
        .enumerate()
        .take(1652)
        .skip(1)
    {
        let metadata = msg_database::msgdb_get_magnetic_spacegroup_type(uni);

        for hall in first_hall as usize..(first_hall + num_halls) as usize {
            hall_pair_count += 1;
            let context = format!("UNI {uni} BNS {} Hall {hall}", metadata.bns_number);

            let Some([candidate_min, candidate_max]) = msgdb_get_uni_candidates(hall) else {
                audit.record("missing_hall_candidate_range", context);
                continue;
            };
            if !(candidate_min..=candidate_max).contains(&uni) {
                audit.record(
                    "uni_outside_hall_candidate_range",
                    format!("{context}: candidates={candidate_min}..={candidate_max}"),
                );
            }

            let hall_type = spg_database::spgdb_get_spacegroup_type(hall);
            if hall_type.number != metadata.number {
                audit.record(
                    "parent_sg_mismatch",
                    format!(
                        "{context}: metadata SG{} but Hall belongs to SG{}",
                        metadata.number, hall_type.number
                    ),
                );
            }

            let Some(magnetic) = msgdb_get_spacegroup_operations(uni, hall) else {
                audit.record("missing_magnetic_operations", context);
                continue;
            };
            if magnetic.is_empty() {
                audit.record("empty_magnetic_operations", context);
                continue;
            }

            let operations: Vec<TestOp> = (0..magnetic.len())
                .map(|index| TestOp {
                    rotation: magnetic.rot[index],
                    translation: magnetic.trans[index],
                    time_reversal: magnetic.timerev[index],
                })
                .collect();

            let mut keys = HashSet::with_capacity(operations.len());
            for (index, op) in operations.iter().enumerate() {
                if mat_get_determinant_i3(&op.rotation).abs() != 1 {
                    audit.record(
                        "invalid_rotation",
                        format!("{context} op {index}: rotation={:?}", op.rotation),
                    );
                }
                match op_key(op) {
                    Some(key) => {
                        if !keys.insert(key) {
                            audit.record(
                                "duplicate_operation",
                                format!("{context} op {index}: {op:?}"),
                            );
                        }
                    }
                    None => audit.record(
                        "invalid_translation",
                        format!("{context} op {index}: {:?}", op.translation),
                    ),
                }
            }

            let identity = OpKey {
                rotation: flatten_rotation(&IDENTITY_ROTATION),
                translation_twelfths: [0, 0, 0],
                time_reversal: false,
            };
            if !keys.contains(&identity) {
                audit.record("missing_identity", context.clone());
            }

            for (left_index, left) in operations.iter().enumerate() {
                let Some(inverse_key) = inverse(left).as_ref().and_then(op_key) else {
                    audit.record(
                        "missing_inverse",
                        format!("{context} op {left_index}: cannot construct inverse"),
                    );
                    continue;
                };
                if !keys.contains(&inverse_key) {
                    audit.record(
                        "missing_inverse",
                        format!("{context} op {left_index}: inverse={inverse_key:?}"),
                    );
                }

                for (right_index, right) in operations.iter().enumerate() {
                    let Some(product_key) = op_key(&compose(left, right)) else {
                        audit.record(
                            "invalid_composed_translation",
                            format!("{context}: product {left_index}*{right_index}"),
                        );
                        continue;
                    };
                    if !keys.contains(&product_key) {
                        audit.record(
                            "not_closed",
                            format!(
                                "{context}: product {left_index}*{right_index}={product_key:?}"
                            ),
                        );
                    }
                }
            }

            let n_unitary = operations
                .iter()
                .filter(|operation| !operation.time_reversal)
                .count();
            let n_antiunitary = operations.len() - n_unitary;
            let anti_identity_zero = operations.iter().any(|operation| {
                operation.time_reversal
                    && operation.rotation == IDENTITY_ROTATION
                    && is_zero_translation(&operation.translation)
            });
            let anti_identity_nonzero = operations.iter().any(|operation| {
                operation.time_reversal
                    && operation.rotation == IDENTITY_ROTATION
                    && !is_zero_translation(&operation.translation)
            });

            let Some(parent) = spg_database::spgdb_get_spacegroup_operations(hall) else {
                audit.record("missing_parent_operations", context);
                continue;
            };
            let parent_keys: HashSet<OpKey> = (0..parent.len())
                .filter_map(|index| {
                    op_key(&TestOp {
                        rotation: parent.rot[index],
                        translation: parent.trans[index],
                        time_reversal: false,
                    })
                })
                .collect();
            let family_keys: HashSet<OpKey> = operations
                .iter()
                .filter_map(|operation| {
                    op_key(&TestOp {
                        time_reversal: false,
                        ..*operation
                    })
                })
                .collect();
            let unitary_keys: HashSet<OpKey> = operations
                .iter()
                .filter(|operation| !operation.time_reversal)
                .filter_map(op_key)
                .collect();
            let reference_matches = if metadata.type_ == MagneticType::AntiTranslation {
                // For Type IV, the Hall mapping names the family space
                // group. H uses the doubled magnetic cell, so its
                // translations need not equal the parent Hall translations
                // and H can have a different international number.
                rotation_multiset(
                    operations
                        .iter()
                        .filter(|operation| !operation.time_reversal)
                        .map(|operation| &operation.rotation),
                ) == rotation_multiset((0..parent.len()).map(|index| &parent.rot[index]))
            } else {
                family_keys == parent_keys
            };
            if !reference_matches {
                audit.record(
                    "reference_group_mismatch",
                    format!(
                        "{context}: family={} unitary={} parent={}",
                        family_keys.len(),
                        unitary_keys.len(),
                        parent_keys.len()
                    ),
                );
            }

            let expected_family_order = parent_keys.len();
            let type_is_consistent = match metadata.type_ {
                MagneticType::Ordinary => {
                    operations.len() == expected_family_order
                        && n_unitary == expected_family_order
                        && n_antiunitary == 0
                        && !anti_identity_zero
                        && !anti_identity_nonzero
                }
                MagneticType::Grey => {
                    operations.len() == 2 * expected_family_order
                        && n_unitary == expected_family_order
                        && n_antiunitary == expected_family_order
                        && anti_identity_zero
                }
                MagneticType::BlackWhite => {
                    operations.len() == expected_family_order
                        && 2 * n_unitary == expected_family_order
                        && n_antiunitary == n_unitary
                        && !anti_identity_zero
                        && !anti_identity_nonzero
                }
                MagneticType::AntiTranslation => {
                    operations.len() == 2 * expected_family_order
                        && n_unitary == expected_family_order
                        && n_antiunitary == n_unitary
                        && !anti_identity_zero
                        && anti_identity_nonzero
                }
                MagneticType::NonMagnetic => false,
            };
            if !type_is_consistent {
                audit.record(
                    "magnetic_type_structure_mismatch",
                    format!(
                        "{context}: type={:?} parent={} total={} U={} A={} theta={} anti-translation={}",
                        metadata.type_,
                        expected_family_order,
                        operations.len(),
                        n_unitary,
                        n_antiunitary,
                        anti_identity_zero,
                        anti_identity_nonzero
                    ),
                );
            }
        }
    }

    println!("Audited all 1651 UNI groups across {hall_pair_count} UNI/Hall settings");
    audit.assert_clean();
}

#[test]
fn enantiomorphic_magnetic_settings_preserve_handedness() {
    for uni in [667usize, 679] {
        let metadata = msg_database::msgdb_get_magnetic_spacegroup_type(uni);
        let first_hall = MAGNETIC_SPACEGROUP_UNI_MAPPING[uni][1] as usize;
        let magnetic = msgdb_get_spacegroup_operations(uni, first_hall).unwrap();
        let lattice = invariant_lattice(&magnetic.rot[..magnetic.len()]).unwrap();
        let dataset = magnetic_spacegroup::msg_identify_with_parent_hall(
            &lattice,
            &magnetic,
            Some(first_hall),
            1e-5,
        )
        .unwrap();

        assert_eq!(
            dataset.uni_number, uni,
            "BNS {} crossed to its enantiomorphic setting",
            metadata.bns_number
        );
        assert_eq!(dataset.hall_number, first_hall);
    }
}

#[test]
fn all_alternative_setting_transformations_are_loaded() {
    let mut hall_pair_count = 0usize;
    let mut nontrivial_setting_count = 0usize;

    for (uni, &[num_halls, first_hall]) in MAGNETIC_SPACEGROUP_UNI_MAPPING
        .iter()
        .enumerate()
        .take(1652)
        .skip(1)
    {
        let num_halls = num_halls as usize;
        let first_hall = first_hall as usize;

        for (hall_offset, encoded) in ALTERNATIVE_TRANSFORMATIONS[uni]
            .iter()
            .take(num_halls)
            .enumerate()
        {
            hall_pair_count += 1;
            let hall = first_hall + hall_offset;
            let encoded_count = encoded.iter().take_while(|&&value| value != 0).count();

            assert!(
                encoded[encoded_count..].iter().all(|&value| value == 0),
                "UNI {uni} Hall {hall} has a nonzero transformation after its sentinel"
            );
            if encoded_count > 0 {
                nontrivial_setting_count += 1;
            }

            let transformations = msgdb_get_std_transformations(uni, hall)
                .unwrap_or_else(|| panic!("missing transformations for UNI {uni} Hall {hall}"));
            assert_eq!(
                transformations.len(),
                encoded_count + 1,
                "wrong transformation count for UNI {uni} Hall {hall}"
            );
            assert_eq!(transformations.rot[0], IDENTITY_ROTATION);
            assert_eq!(transformations.trans[0], [0.0; 3]);
        }
    }

    // Upstream spglib v2.5.0 has 450 nontrivial UNI/Hall rows. The old
    // converter retained only the two rows whose C initializers had all seven
    // integers and silently discarded the other 448 partial initializers.
    assert_eq!(hall_pair_count, 4479);
    assert_eq!(nontrivial_setting_count, 450);

    let transformations = msgdb_get_std_transformations(132, 116).unwrap();
    assert_eq!(transformations.len(), 2);
    assert_eq!(transformations.rot[1], [[0, -1, 0], [-1, 0, 0], [0, 0, -1]]);
    assert_eq!(transformations.trans[1], [0.0, 0.0, 0.25]);
}

#[test]
fn all_database_settings_round_trip_with_parent_hint() {
    let mut audit = Audit::default();
    let mut exact_matches = 0usize;

    for (uni, &[num_halls, first_hall]) in MAGNETIC_SPACEGROUP_UNI_MAPPING
        .iter()
        .enumerate()
        .take(1652)
        .skip(1)
    {
        let metadata = msg_database::msgdb_get_magnetic_spacegroup_type(uni);
        let num_halls = num_halls as usize;
        let first_hall = first_hall as usize;
        for hall_offset in 0..num_halls {
            let hall = first_hall + hall_offset;
            let context = format!("UNI {uni} BNS {} Hall {hall}", metadata.bns_number);
            let Some(magnetic) = msgdb_get_spacegroup_operations(uni, hall) else {
                audit.record("missing_magnetic_operations", context);
                continue;
            };
            let rotations = magnetic.rot[..magnetic.len()].to_vec();
            let Some(lattice) = invariant_lattice(&rotations) else {
                audit.record("invariant_lattice_failed", context);
                continue;
            };

            let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                magnetic_spacegroup::msg_identify_with_parent_hall(
                    &lattice,
                    &magnetic,
                    Some(hall),
                    1e-5,
                )
            }));
            match outcome {
                Ok(Ok(dataset))
                    if dataset.uni_number == uni
                        && dataset.msg_type == metadata.type_
                        && dataset.hall_number == hall =>
                {
                    exact_matches += 1;
                }
                Ok(Ok(dataset)) => {
                    let returned =
                        msg_database::msgdb_get_magnetic_spacegroup_type(dataset.uni_number);
                    let category = if dataset.msg_type == metadata.type_ {
                        "wrong_uni_same_type"
                    } else {
                        "wrong_uni_wrong_type"
                    };
                    audit.record(
                        category,
                        format!(
                            "{context}: returned UNI {} BNS {} type {:?} Hall {}",
                            dataset.uni_number,
                            returned.bns_number,
                            dataset.msg_type,
                            dataset.hall_number
                        ),
                    );
                }
                Ok(Err(error)) => {
                    audit.record(
                        identification_error_category(error),
                        format!("{context}: {error:?}"),
                    );
                }
                Err(_) => audit.record("identification_panicked", context),
            }
        }
    }

    println!("Exact all-setting round-trips: {exact_matches} / 4479");
    assert_eq!(exact_matches, 4479);
    audit.assert_clean();
}

#[test]
fn automatic_all_setting_round_trips_are_unique_or_explicitly_ambiguous() {
    let mut setting_count = 0usize;
    let mut exact_matches = 0usize;
    let mut explicit_ambiguities = 0usize;

    for (uni, &[num_halls, first_hall]) in MAGNETIC_SPACEGROUP_UNI_MAPPING
        .iter()
        .enumerate()
        .take(1652)
        .skip(1)
    {
        let metadata = msg_database::msgdb_get_magnetic_spacegroup_type(uni);
        let num_halls = num_halls as usize;
        let first_hall = first_hall as usize;
        for hall_offset in 0..num_halls {
            setting_count += 1;
            let hall = first_hall + hall_offset;
            let magnetic = msgdb_get_spacegroup_operations(uni, hall)
                .unwrap_or_else(|| panic!("missing operations for UNI {uni} Hall {hall}"));
            let lattice = invariant_lattice(&magnetic.rot[..magnetic.len()])
                .unwrap_or_else(|| panic!("invariant lattice failed for UNI {uni} Hall {hall}"));
            let result = magnetic_spacegroup::msg_identify_magnetic_space_group_type(
                &lattice, &magnetic, 1e-5,
            );

            match result {
                Ok(dataset) => {
                    assert_eq!(
                        dataset.uni_number, uni,
                        "automatic round-trip crossed UNI {uni} BNS {} Hall {hall}",
                        metadata.bns_number
                    );
                    assert_eq!(dataset.msg_type, metadata.type_);
                    exact_matches += 1;
                }
                Err(SymError::MagneticUniAmbiguous) => {
                    assert!(
                        matches!(uni, 275 | 277 | 282 | 284),
                        "unexpected ambiguity for UNI {uni} BNS {} Hall {hall}",
                        metadata.bns_number
                    );
                    explicit_ambiguities += 1;
                }
                Err(error) => {
                    panic!(
                        "automatic round-trip failed for UNI {uni} BNS {} Hall {hall}: {error:?}",
                        metadata.bns_number
                    );
                }
            }
        }
    }

    assert_eq!(setting_count, 4479);
    assert_eq!(exact_matches, 4461);
    assert_eq!(explicit_ambiguities, 18);
}
