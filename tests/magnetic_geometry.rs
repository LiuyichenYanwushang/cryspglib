//! R9 stage 2a: the magnetic isotropy record's geometry, exactly, and how far it
//! embeds into the parent.
//!
//! The record's `basis`/`origin` are integers in the parent's **primitive**
//! frame; this file pins their exact reading, the cross-table comparison against
//! the ordinary isotropy table (entry-level and setting-free), and the measured
//! parent-frame containment of the operation sets.  The containment counts are
//! a **measurement** on purpose: the naive `basis^T` map does not contain for
//! every record (224 parents have at least one unitary class outside the parent),
//! which is the stage-2b work -- aligning the magnetic group's tabulated setting
//! with the record's lattice.

use cryspglib::irrep::magnetic_embedding::{
    GeometryAgreement, GeometryError, compare_with_ordinary_geometry, embedded_translation_lattice,
    geometry_of, magnetic_operations, measure_parent_containment, search_parent_setting,
    translation_lattice,
};
use cryspglib::irrep::query;
use cryspglib::irrep::magnetic_embedding::canonical_translation;
use cryspglib::irrep::subduction::Rat;

fn first_record(sg: u8) -> cryspglib::irrep::types::MagneticIsotropyRecord {
    query::magnetic_isotropy_subgroups_of(sg)
        .first()
        .expect("a magnetic isotropy row")
        .subgroup
}

#[test]
fn every_magnetic_record_has_exact_readable_geometry() {
    let mut rows = 0usize;
    let mut indices: std::collections::BTreeMap<i128, usize> = std::collections::BTreeMap::new();
    let mut origin_denominators: std::collections::BTreeMap<i32, usize> =
        std::collections::BTreeMap::new();
    let mut agreements: std::collections::BTreeMap<&'static str, usize> =
        std::collections::BTreeMap::new();
    let mut labelled = 0usize;
    let mut same_lattice = 0usize;
    for sg in 1..=230u8 {
        for row in query::magnetic_isotropy_subgroups_of(sg) {
            rows += 1;
            let record = row.subgroup;
            let geometry = geometry_of(sg, &record)
                .unwrap_or_else(|error| panic!("sg {sg} UNI {}: {error}", record.mag_sg));
            assert!(geometry.index >= 1, "sg {sg}: a lattice index is positive");
            *indices.entry(geometry.index).or_insert(0) += 1;
            *origin_denominators.entry(record.origin[3]).or_insert(0) += 1;
            let agreement = compare_with_ordinary_geometry(sg, row.ml_label, record.direction, &geometry);
            *agreements
                .entry(match agreement {
                    GeometryAgreement::Agree => "agree",
                    GeometryAgreement::BasisDiffers => "basis_differs",
                    GeometryAgreement::OriginDiffers => "origin_differs",
                    GeometryAgreement::BothDiffer => "both_differ",
                    GeometryAgreement::NoOrdinaryRow => "no_ordinary_row",
                    GeometryAgreement::NoOrdinaryIrrep => "no_ordinary_irrep",
                })
                .or_insert(0) += 1;
            if let Some(irrep) = query::irreps_of(sg)
                .iter()
                .find(|irrep| irrep.ml == row.ml_label)
            {
                let labelled_rows: Vec<_> = irrep
                    .subgroups()
                    .iter()
                    .filter(|ordinary| ordinary.direction_label == record.direction)
                    .collect();
                if !labelled_rows.is_empty() {
                    labelled += 1;
                    if labelled_rows
                        .iter()
                        .any(|ordinary| geometry.spans_same_lattice_as(ordinary.basis))
                    {
                        same_lattice += 1;
                    }
                }
            }
        }
    }
    assert_eq!(rows, 16_721, "magnetic isotropy rows");
    assert_eq!(
        indices,
        std::collections::BTreeMap::from([
            (1, 2_064),
            (2, 4_859),
            (3, 667),
            (4, 5_483),
            (6, 496),
            (8, 2_146),
            (16, 413),
            (32, 593),
        ]),
        "lattice index |det(basis)| distribution"
    );
    assert_eq!(
        origin_denominators,
        std::collections::BTreeMap::from([
            (1, 7_849),
            (2, 4_691),
            (3, 297),
            (4, 2_969),
            (6, 133),
            (8, 687),
            (12, 95),
        ]),
        "origin denominators"
    );
    // Entry-level comparison with the ordinary isotropy table.
    assert_eq!(
        agreements,
        std::collections::BTreeMap::from([
            ("agree", 8_287),
            ("basis_differs", 3_852),
            ("origin_differs", 600),
            ("both_differ", 2_500),
            ("no_ordinary_row", 1_482),
        ])
    );
    // Setting-free comparison: the two tables keep the same lattice far more
    // often than they keep the same representative.
    assert_eq!(labelled, 15_239, "rows with a labelled ordinary row");
    // Round r9b, item 1: with the change of basis on the correct side every
    // labelled pair spans the same lattice; the earlier 12,893 was the column
    // -lattice count, i.e. 2,346 false negatives and 0 false positives.
    assert_eq!(same_lattice, 15_239, "rows whose ordinary row spans the same lattice");
}

#[test]
fn the_lattice_comparison_explains_the_entry_level_differences() {
    // SG 1 GM1: the ordinary row keeps the identity basis, the magnetic row
    // keeps a rotated unimodular one -- a setting change, not a different
    // subgroup.
    let record = first_record(1);
    let geometry = geometry_of(1, &record).expect("geometry");
    assert_eq!(record.mag_sg, 1);
    assert_eq!(record.basis, [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]);
    assert_eq!(geometry.index, 1);
    assert!(geometry.spans_same_lattice_as([[1, 0, 0], [0, 1, 0], [0, 0, 1]]));
    assert_eq!(
        compare_with_ordinary_geometry(1, "GM1", "P1", &geometry),
        GeometryAgreement::BasisDiffers
    );
    // A counterexample from the review: SG 1 X1 P1 keeps 2Z x Z x Z in two
    // different representatives, so the setting-free test must accept it.
    let repeated = query::magnetic_isotropy_subgroups_of(1)
        .into_iter()
        .find(|row| row.ml_label == "X1" && row.subgroup.direction == "P1")
        .expect("SG 1 X1 P1 row");
    let repeated_geometry = geometry_of(1, &repeated.subgroup).expect("geometry");
    assert_eq!(repeated.subgroup.mag_sg, 3);
    assert_eq!(repeated_geometry.index, 2);
    assert!(repeated_geometry.spans_same_lattice_as([[2, 0, 0], [0, 1, 0], [0, 0, 1]]));
    // Negative control: a genuine lattice change is refused.
    assert!(!repeated_geometry.spans_same_lattice_as([[1, 0, 0], [0, 1, 0], [0, 0, 1]]));
    assert!(!repeated_geometry.spans_same_lattice_as([[1, 0, 0], [0, 1, 0], [0, 0, 2]]));

    // A row where both tables agree outright: SG 14 GM1+ P1 (UNI 82).
    let agreed = query::magnetic_isotropy_subgroups_of(14)
        .into_iter()
        .find(|row| row.ml_label == "GM1+" && row.subgroup.direction == "P1")
        .expect("SG 14 GM1+ P1 row");
    let agreed_geometry = geometry_of(14, &agreed.subgroup).expect("geometry");
    assert_eq!(agreed.subgroup.mag_sg, 82);
    assert_eq!(
        compare_with_ordinary_geometry(14, "GM1+", "P1", &agreed_geometry),
        GeometryAgreement::Agree
    );
    assert!(agreed_geometry.spans_same_lattice_as([[1, 0, 0], [0, 1, 0], [0, 0, 1]]));
}

/// The setting-free test agrees with an independent membership computation:
/// same lattice iff every row of one basis lies in the other's lattice, judged
/// by `canonical_translation` (a different code path from the matrix inverse).
#[test]
fn the_setting_free_test_agrees_with_lattice_membership() {
    let mut checked = 0usize;
    for sg in 1..=230u8 {
        let irreps = query::irreps_of(sg);
        for row in query::magnetic_isotropy_subgroups_of(sg) {
            let record = row.subgroup;
            let geometry = geometry_of(sg, &record).expect("geometry");
            let Some(irrep) = irreps.iter().find(|irrep| irrep.ml == row.ml_label) else {
                continue;
            };
            for ordinary in irrep
                .subgroups()
                .iter()
                .filter(|ordinary| ordinary.direction_label == record.direction)
            {
                checked += 1;
                let magnetic_lattice: Vec<[Rat; 3]> = geometry.basis.to_vec();
                let ordinary_lattice: Vec<[Rat; 3]> = ordinary
                    .basis
                    .iter()
                    .map(|values| values.map(|value| Rat::from_integer(i128::from(value))))
                    .collect();
                let contains = |lattice: &[[Rat; 3]], row: [i32; 3]| {
                    let vector = row.map(|value| Rat::from_integer(i128::from(value)));
                    canonical_translation(vector, lattice)
                        .iter()
                        .all(|value| value.is_zero())
                };
                let to_ints = |row: &[Rat; 3]| {
                    let mut out = [0i32; 3];
                    for (axis, value) in row.iter().enumerate() {
                        out[axis] = value.to_i32().expect("integer basis");
                    }
                    out
                };
                // Row lattices are equal iff each basis lies in the other's.
                let forward = ordinary.basis.iter().all(|row| contains(&magnetic_lattice, *row));
                let backward = geometry
                    .basis
                    .iter()
                    .all(|row| contains(&ordinary_lattice, to_ints(row)));
                let expected = forward && backward;
                assert_eq!(
                    geometry.spans_same_lattice_as(ordinary.basis),
                    expected,
                    "sg {sg} UNI {} direction {}: membership says {expected}",
                    record.mag_sg,
                    record.direction
                );
            }
        }
    }
    assert_eq!(checked, 15_239, "labelled pairs compared");
}

/// R9 stage 2b measurement (round 15j): the candidate-setting search is an
/// instrument, not a decision procedure.  Pinning the histogram keeps the
/// negative result from being quietly upgraded into a containment claim.
#[test]
fn the_setting_search_family_is_not_a_decision_procedure() {
    let mut histogram: std::collections::BTreeMap<usize, usize> =
        std::collections::BTreeMap::new();
    let mut unique = 0usize;
    for sg in 1..=230u8 {
        for row in query::magnetic_isotropy_subgroups_of(sg) {
            let record = row.subgroup;
            let geometry = geometry_of(sg, &record).expect("geometry");
            let set = magnetic_operations(record.mag_sg).expect("operations");
            let report = search_parent_setting(&geometry, &set).expect("search");
            assert_eq!(report.candidates, 48);
            assert_eq!(report.naive.classes, set.class_count());
            *histogram.entry(report.survivors).or_insert(0) += 1;
            if report.survivors == 1 {
                unique += 1;
            }
        }
    }
    // Round 15k: this histogram was wrong while the candidate base multiplied
    // by the group's centring matrix; the identity map already contains
    // UNI 1221/1333, so "zero survivors" was an artefact, not containment
    // evidence.  The pinned numbers below are from the fixed base (the record's
    // own basis).
    assert_eq!(
        histogram,
        std::collections::BTreeMap::from([
            (0, 10_151),
            (2, 96),
            (3, 2),
            (4, 247),
            (6, 16),
            (8, 1_429),
            (16, 1_875),
            (24, 55),
            (32, 62),
            (48, 2_788),
        ]),
        "survivor histogram of the signed-permutation family"
    );
    assert_eq!(unique, 0, "no record has a unique survivor");

    // The invariant that would have caught the round-15k bug: whenever the
    // naive map already contains the record completely, the search must find at
    // least the identity candidate, so a zero survivor count is impossible.
    let mut contained_but_zero = 0usize;
    let mut permissive = 0usize;
    for sg in 1..=230u8 {
        for row in query::magnetic_isotropy_subgroups_of(sg) {
            let record = row.subgroup;
            let geometry = geometry_of(sg, &record).expect("geometry");
            let set = magnetic_operations(record.mag_sg).expect("operations");
            let report = search_parent_setting(&geometry, &set).expect("search");
            let full = report.naive.unitary_in_parent == report.naive.unitary_classes
                && report.naive.antiunitary_rotations_in_parent == set.antiunitary_classes;
            if full {
                permissive += 1;
                if report.survivors == 0 {
                    contained_but_zero += 1;
                }
            }
        }
    }
    // Fully contained = every unitary class AND every antiunitary rotation lands
    // in the parent; the earlier probe's 6,375 counted the unitary part alone.
    assert_eq!(permissive, 5_972, "records the naive map already contains");
    assert_eq!(
        contained_but_zero, 0,
        "a contained record must have at least one surviving candidate"
    );

    // The identity-map witnesses the reviewer used to refute the old claim.
    for (parent_sg, uni, survivors) in [(167u8, 1333usize, 4usize), (142, 1221, 8)] {
        let record = query::magnetic_isotropy_subgroups_of(parent_sg)
            .into_iter()
            .find(|row| row.subgroup.mag_sg == uni)
            .map(|row| row.subgroup)
            .expect("record");
        let geometry = geometry_of(parent_sg, &record).expect("geometry");
        let set = magnetic_operations(uni).expect("operations");
        let report = search_parent_setting(&geometry, &set).expect("search");
        assert_eq!(report.survivors, survivors, "SG {parent_sg} UNI {uni}");
        assert_eq!(
            report.naive.unitary_in_parent, report.naive.unitary_classes,
            "SG {parent_sg} UNI {uni} is fully contained by the naive map"
        );
    }

    // Two witnesses of the two failure modes.
    let failing = query::magnetic_isotropy_subgroups_of(3)
        .into_iter()
        .map(|row| row.subgroup)
        .find(|record| record.mag_sg == 24)
        .expect("SG 3 row for UNI 24");
    let geometry = geometry_of(3, &failing).expect("geometry");
    let set = magnetic_operations(24).expect("operations");
    let report = search_parent_setting(&geometry, &set).expect("search");
    assert_eq!(report.survivors, 0, "SG 3 UNI 24: no candidate contains");
    assert!(
        report.naive.unitary_in_parent < report.naive.unitary_classes,
        "and the naive map already misses a class"
    );

    let permissive = first_record(1);
    let geometry = geometry_of(1, &permissive).expect("geometry");
    let set = magnetic_operations(permissive.mag_sg).expect("operations");
    let report = search_parent_setting(&geometry, &set).expect("search");
    assert_eq!(report.survivors, 48, "P1 parent: every candidate contains");
}

#[test]
fn tampering_with_the_record_geometry_is_caught() {
    let record = first_record(14);
    let geometry = geometry_of(14, &record).expect("geometry");
    assert_eq!(
        compare_with_ordinary_geometry(14, "GM1+", "P1", &geometry),
        GeometryAgreement::Agree
    );

    // Move the origin by 1/4: the anchor sees it.
    let mut moved_origin = record;
    moved_origin.origin = [1, 0, 0, 4];
    let moved = geometry_of(14, &moved_origin).expect("geometry");
    assert_eq!(
        compare_with_ordinary_geometry(14, "GM1+", "P1", &moved),
        GeometryAgreement::OriginDiffers
    );
    assert_ne!(moved.origin, geometry.origin);

    // Replace the basis by a different lattice: the anchor sees that too.
    let mut moved_basis = record;
    moved_basis.basis = [[1, 0, 0], [0, 1, 0], [0, 0, 2]];
    let moved = geometry_of(14, &moved_basis).expect("geometry");
    assert_eq!(moved.index, 2);
    assert_eq!(
        compare_with_ordinary_geometry(14, "GM1+", "P1", &moved),
        GeometryAgreement::BasisDiffers
    );
    assert!(!moved.spans_same_lattice_as(record.basis));

    // A non-positive origin denominator and a singular basis are refused.
    let mut bad_origin = record;
    bad_origin.origin = [0, 0, 0, 0];
    assert_eq!(
        geometry_of(14, &bad_origin),
        Err(GeometryError::BadOriginEncoding {
            uni: record.mag_sg,
            denominator: 0,
        })
    );
    let mut singular = record;
    singular.basis = [[1, 0, 0], [1, 0, 0], [0, 0, 1]];
    assert_eq!(
        geometry_of(14, &singular),
        Err(GeometryError::DegenerateBasis { uni: record.mag_sg })
    );
    assert_eq!(
        geometry_of(0, &record),
        Err(GeometryError::UnknownParentSpaceGroup { sg: 0 })
    );
}

/// The containment counts are measured, not asserted: this test pins the
/// measurement so that the stage-2b gap cannot silently change.
#[test]
fn parent_containment_of_the_naive_map_is_measured() {
    let mut primitive = (0usize, 0usize, 0usize);
    let mut centred = (0usize, 0usize, 0usize);
    let mut failing_parents = 0usize;
    let mut per_parent: std::collections::BTreeMap<u8, (usize, usize, usize)> =
        std::collections::BTreeMap::new();
    for sg in 1..=230u8 {
        let parent = query::symmetry_operations_of(sg).expect("parent operations");
        let is_centred = parent.operations.iter().any(|operation| {
            operation.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                && operation
                    .translation
                    .iter()
                    .any(|value| (value - value.round()).abs() > 1e-9)
        });
        let rows = query::magnetic_isotropy_subgroups_of(sg);
        if rows.is_empty() {
            continue;
        }
        let mut totals = (0usize, 0usize, 0usize);
        for row in rows {
            let record = row.subgroup;
            let geometry = geometry_of(sg, &record).expect("geometry");
            let set = magnetic_operations(record.mag_sg).expect("operations");
            let containment =
                measure_parent_containment(&geometry, &set).expect("containment measurement");
            assert_eq!(containment.classes, set.class_count());
            assert_eq!(containment.unitary_classes, set.unitary_classes);
            assert!(containment.unitary_in_parent <= containment.unitary_classes);
            assert!(containment.antiunitary_rotations_in_parent <= set.antiunitary_classes);
            assert!(containment.integral_rotations <= set.class_count());
            totals.0 += 1;
            totals.1 += containment.unitary_in_parent;
            totals.2 += containment.unitary_classes;
        }
        if is_centred {
            centred.0 += totals.0;
            centred.1 += totals.1;
            centred.2 += totals.2;
        } else {
            primitive.0 += totals.0;
            primitive.1 += totals.1;
            primitive.2 += totals.2;
        }
        if totals.1 < totals.2 {
            failing_parents += 1;
            per_parent.insert(sg, totals);
        }
    }
    assert_eq!(primitive, (9_726, 33_894, 48_618), "primitive parents");
    assert_eq!(centred, (6_995, 13_918, 35_471), "centred parents");
    assert_eq!(failing_parents, 224, "parents with a unitary class outside the parent");
    assert_eq!(per_parent.get(&3), Some(&(16, 18, 31)));
}

#[test]
fn the_embedded_lattice_is_the_image_of_the_records_lattice() {
    // Pick a record whose source lattice is non-trivial: UNI 20 (centred).
    let mut chosen = None;
    'outer: for sg in 1..=230u8 {
        for row in query::magnetic_isotropy_subgroups_of(sg) {
            if row.subgroup.mag_sg == 20 {
                chosen = Some((sg, row.subgroup));
                break 'outer;
            }
        }
    }
    let (sg, record) = chosen.expect("a record referencing UNI 20");
    let geometry = geometry_of(sg, &record).expect("geometry");
    let set = magnetic_operations(20).expect("operations");
    let source = translation_lattice(&set);
    assert!(!source.is_empty(), "UNI 20 lists a centring translation");
    let image = embedded_translation_lattice(&geometry, &set);
    assert_eq!(image.len(), source.len(), "one image per generator");
    // Every image is a lattice vector of the image lattice, and each generator's
    // image is exactly the mapped generator.
    let map = geometry.map_matrix();
    for (generator, mapped) in source.iter().zip(image.iter()) {
        let moved = map
            .checked_mul_vector(&cryspglib::irrep::subduction::Vec3R::new(*generator))
            .expect("exact product");
        assert_eq!(moved.as_array(), mapped);
        assert!(
            canonical_translation(*mapped, &image)
                .iter()
                .all(|value| value.is_zero()),
            "an image must be a lattice vector of the image lattice"
        );
    }
    // An identity geometry maps the lattice to itself.
    let identity = geometry_of(
        14,
        &query::magnetic_isotropy_subgroups_of(14)
            .first()
            .expect("row")
            .subgroup,
    )
    .expect("geometry");
    let identity_set = magnetic_operations(
        query::magnetic_isotropy_subgroups_of(14)
            .first()
            .expect("row")
            .subgroup
            .mag_sg,
    )
    .expect("operations");
    assert_eq!(
        embedded_translation_lattice(&identity, &identity_set),
        translation_lattice(&identity_set)
    );
}
