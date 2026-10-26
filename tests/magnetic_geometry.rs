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

use cryspglib::irrep::magnetic_embedding::canonical_translation;
use cryspglib::irrep::magnetic_embedding::{
    EmbeddingKind, GeometryAgreement, GeometryError, MagneticGeometry, MagneticOperation,
    MagneticOperationSet, basis_in_parent_conventional, compare_with_ordinary_geometry,
    conventional_lattice_index, embed_in_parent_conventional, embedded_translation_lattice,
    embedding_shift_components, full_translation_lattice, geometry_of, magnetic_operations,
    measure_parent_containment, origin_in_parent_conventional, parent_embedding_verdict,
    search_parent_setting, translation_lattice, verify_embedding,
};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::{Lattice, Mat3R, Rat, Vec3R};

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
            let agreement =
                compare_with_ordinary_geometry(sg, row.ml_label, record.direction, &geometry);
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
    assert_eq!(
        same_lattice, 15_239,
        "rows whose ordinary row spans the same lattice"
    );
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
/// by [`Lattice::coordinates`] (an independent, separately audited code path).
///
/// Round 20: the previous version of this test asked `canonical_translation`
/// whether an **integer** row lies in the other lattice.  That function reduces
/// every component modulo one, so every integer vector is "in" every lattice
/// and the comparison was vacuously true -- `e_1` passed for `2Z x Z x Z`.  The
/// negative control below pins that the replacement rejects it.
#[test]
fn the_setting_free_test_agrees_with_lattice_membership() {
    // Negative control: the real membership test must reject a proper sublattice.
    let sublattice = Lattice::new(Mat3R::new([
        [Rat::from_integer(2), Rat::ZERO, Rat::ZERO],
        [Rat::ZERO, Rat::ONE, Rat::ZERO],
        [Rat::ZERO, Rat::ZERO, Rat::ONE],
    ]))
    .expect("non-singular");
    assert!(
        !sublattice
            .coordinates(&Vec3R::new([Rat::ONE, Rat::ZERO, Rat::ZERO]))
            .expect("coordinates")
            .as_array()
            .iter()
            .all(|value| value.is_integer()),
        "e_1 is not in 2Z x Z x Z"
    );

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
                let magnetic_lattice =
                    Lattice::new(Mat3R::new(geometry.basis)).expect("non-singular basis");
                let mut ordinary_rows = [[Rat::ZERO; 3]; 3];
                for (row_index, values) in ordinary.basis.iter().enumerate() {
                    for (column, value) in values.iter().enumerate() {
                        ordinary_rows[row_index][column] = Rat::from_integer(i128::from(*value));
                    }
                }
                let ordinary_lattice =
                    Lattice::new(Mat3R::new(ordinary_rows)).expect("non-singular basis");
                let contains = |lattice: &Lattice, row: [Rat; 3]| {
                    lattice
                        .coordinates(&Vec3R::new(row))
                        .expect("coordinates")
                        .as_array()
                        .iter()
                        .all(|value| value.is_integer())
                };
                // Row lattices are equal iff each basis lies in the other's.
                let forward = ordinary.basis.iter().all(|row| {
                    contains(
                        &magnetic_lattice,
                        row.map(|v| Rat::from_integer(i128::from(v))),
                    )
                });
                let backward = geometry
                    .basis
                    .iter()
                    .all(|row| contains(&ordinary_lattice, *row));
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
///
/// Round 20 note: this instrument keeps its **class-level** acceptance criterion
/// (and its raw, unconverted frames), so its histogram is deliberately unchanged
/// by the containment gate.  A survivor here is not a containment proof; use
/// `verify_embedding` for a verdict.
#[test]
fn the_setting_search_family_is_not_a_decision_procedure() {
    let mut histogram: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
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
    // Round 15o: the antiunitary classes now get the same translation equation
    // as the unitary ones (before, only their rotations were checked), so the
    // histogram moved.  The intermediate state -- search strict, naive criterion
    // still rotation-only -- reported 115 rows as contained-but-zero; the code
    // before that fix reported 0, because the identity candidate also satisfied
    // the weak check then (round 15p, L3).
    //
    // Round 15k: this histogram was wrong while the candidate base multiplied
    // by the group's centring matrix; the identity map already places
    // UNI 1221/1333, so "zero survivors" was an artefact, not containment
    // evidence.  The pinned numbers below are from the fixed base (the record's
    // own basis).
    assert_eq!(
        histogram,
        std::collections::BTreeMap::from([
            (0, 10_524),
            (2, 120),
            (3, 4),
            (4, 237),
            (6, 14),
            (8, 1_414),
            (16, 1_623),
            (24, 54),
            (32, 47),
            (48, 2_684),
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
                && report.naive.antiunitary_classes_in_parent == set.antiunitary_classes;
            if full {
                permissive += 1;
                if report.survivors == 0 {
                    contained_but_zero += 1;
                }
            }
        }
    }
    // Fully contained = every unitary class AND every antiunitary class lands in
    // the parent (round 15o extended the measurement to place the primed
    // classes, not just their rotations); the earlier probe's 6,375 counted the
    // unitary part alone.
    assert_eq!(permissive, 5_818, "records the naive map already contains");
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

/// The exact conversion to the parent's conventional frame must agree with the
/// isotropy helper that performs the same conversion in floating point, and the
/// conventional lattice index must be the pinned census.
#[test]
fn the_conventional_frame_conversion_matches_the_isotropy_helper() {
    let mut histogram: std::collections::BTreeMap<(i128, i128), usize> =
        std::collections::BTreeMap::new();
    let mut rows = 0usize;
    for sg in 1..=230u8 {
        for row in query::magnetic_isotropy_subgroups_of(sg) {
            rows += 1;
            let record = row.subgroup;
            let geometry = geometry_of(sg, &record).expect("geometry");
            let exact = basis_in_parent_conventional(&geometry).expect("conversion");
            // The ordinary helper has a canonical primitive-frame contract.
            // Re-express raw monoclinic coefficients independently before comparing.
            let mut canonical = record;
            if (3..=15).contains(&sg) {
                for v in &mut canonical.basis {
                    v.rotate_left(1);
                }
                canonical.origin[..3].rotate_left(1);
            }
            let helper =
                cryspglib::irrep::isotropy::basis_in_parent_conventional(sg, canonical.basis)
                    .expect("helper");
            for axis in 0..3 {
                for column in 0..3 {
                    let value = exact[axis][column].to_f64();
                    assert!(
                        (value - helper[axis][column]).abs() < 1e-12,
                        "sg {sg} UNI {}: basis[{axis}][{column}] {value} vs {}",
                        record.mag_sg,
                        helper[axis][column]
                    );
                }
            }
            let exact_origin = origin_in_parent_conventional(&geometry).expect("origin");
            let helper_origin = cryspglib::irrep::isotropy::origin_shift_in_parent_conventional(
                sg,
                canonical.origin,
            )
            .expect("helper origin");
            for axis in 0..3 {
                assert!(
                    (exact_origin[axis].to_f64() - helper_origin[axis]).abs() < 1e-12,
                    "sg {sg} UNI {}: origin[{axis}]",
                    record.mag_sg
                );
            }
            let index = conventional_lattice_index(&geometry).expect("index");
            assert!(index.numerator() > 0, "the index is positive");
            *histogram
                .entry((index.numerator(), index.denominator()))
                .or_insert(0) += 1;
        }
    }
    assert_eq!(rows, 16_721);
    assert_eq!(
        histogram,
        std::collections::BTreeMap::from([
            ((1, 1), 3_204),
            ((1, 2), 707),
            ((1, 3), 41),
            ((1, 4), 262),
            ((2, 1), 5_575),
            ((2, 3), 66),
            ((3, 1), 667),
            ((4, 1), 4_223),
            ((4, 3), 108),
            ((6, 1), 496),
            ((8, 1), 1_333),
            ((8, 3), 39),
        ]),
        "conventional-frame index |det(W * P)| histogram"
    );
}

/// The round-17 witnesses: records that a family built only from the record's
/// basis could not embed, but which the identity map embeds trivially (the
/// tabulated operations are already parent operations).  The search must find
/// the identity, and the embedding must also realise the record's lattice.
#[test]
fn the_parent_embedding_search_finds_the_identity_witnesses() {
    for (parent_sg, uni) in [(88u8, 741usize), (142, 1221), (227, 1630)] {
        let record = query::magnetic_isotropy_subgroups_of(parent_sg)
            .into_iter()
            .find(|row| row.subgroup.mag_sg == uni)
            .map(|row| row.subgroup)
            .expect("record");
        let geometry = geometry_of(parent_sg, &record).expect("geometry");
        let set = magnetic_operations(uni).expect("operations");
        let embeddings = embed_in_parent_conventional(&geometry, &set).expect("search");
        assert!(
            !embeddings.is_empty(),
            "SG {parent_sg} UNI {uni} must embed"
        );
        let identity = embeddings
            .iter()
            .find(|entry| entry.map == cryspglib::irrep::subduction::Mat3R::identity())
            .unwrap_or_else(|| panic!("SG {parent_sg} UNI {uni}: identity map missing"));
        assert_eq!(identity.kind, EmbeddingKind::Plain);
        assert!(
            identity.shift.iter().all(|value| value.is_zero()),
            "SG {parent_sg} UNI {uni}: the identity map needs no shift"
        );
        assert!(
            identity.realises_record_lattice,
            "SG {parent_sg} UNI {uni}: the identity map realises the record lattice"
        );
    }
}

/// SG 3 UNI 24: the family returns nothing, and round 23 certifies that as a
/// property of the family.  The parent (P 1 2 1) has only the identity and one
/// twofold axis, both with zero translation, while the record's source
/// translation lattice contains the centring shift (1/2,1/2,0) and its twofold
/// class carries it.  Every one of the family's 240 candidate maps is refuted by
/// an exact necessary condition, so the empty result is a proof about the
/// family.  It is **not** a proof that the record is not a subgroup: the
/// companion test below exhibits a map outside the family -- the search returns
/// nothing for this record, while the public verifier accepts that map -- and so
/// the record does embed.
#[test]
fn the_unresolved_witness_is_reported_as_unresolved() {
    let record = query::magnetic_isotropy_subgroups_of(3)
        .into_iter()
        .find(|row| row.subgroup.mag_sg == 24)
        .map(|row| row.subgroup)
        .expect("record");
    let geometry = geometry_of(3, &record).expect("geometry");
    let set = magnetic_operations(24).expect("operations");
    let verdict = parent_embedding_verdict(&geometry, &set).expect("search");
    assert_eq!(
        verdict.embeddings,
        Vec::new(),
        "the family does not resolve SG 3 UNI 24"
    );
    // Round 23: this is no longer "the walk found nothing" but a proof about the
    // candidate family.  Every one of its 240 maps is refuted by one of the three
    // exact necessary conditions of an embedding: the image of the full source
    // translation lattice must lie in the parent's lattice; a conjugated class
    // rotation must stay in the parent's rotation set; and every class must be
    // placeable -- the joint anchor congruence `(I - R_p) delta = tau_p - A tau`
    // must be solvable modulo the parent lattice for some parent operation with
    // the mapped rotation.  The counters are exclusive in the order the
    // production search screens them, so the split below is that order.  Round
    // 22's per-normal anchor test was withdrawn (round 22b): it never ran for
    // this record, and the version shipped was not the joint congruence.
    println!(
        "SG 3 UNI 24: source={} rotation={} anchor={} residual={}",
        verdict.excluded_by_source_lattice,
        verdict.excluded_by_rotation,
        verdict.excluded_by_anchor_congruence,
        verdict.candidates_without_witness,
    );
    assert!(
        verdict.absence_is_certified,
        "the absence of SG 3 UNI 24 must be certified, candidates_without_witness={}",
        verdict.candidates_without_witness
    );
    assert_eq!(verdict.candidates_without_witness, 0);
    assert_eq!(
        verdict.excluded_by_source_lattice
            + verdict.excluded_by_rotation
            + verdict.excluded_by_anchor_congruence,
        240,
        "every candidate map is excluded by an exact necessary condition"
    );
    assert_eq!(
        (
            verdict.excluded_by_rotation,
            verdict.excluded_by_source_lattice,
            verdict.excluded_by_anchor_congruence,
        ),
        (0, 240, 0),
        "the split depends on the order of the three screens: source lattice, then \
         rotation set, then per-class placeability.  This search tests the source \
         lattice first, so all 240 maps are attributed to it (the round-22 review, \
         which tested rotations first, measured 224 rotation + 16 source)"
    );
    // The parent really has no translation to absorb the centring shift.
    let parent = query::symmetry_operations_of(3).expect("parent");
    assert_eq!(parent.operations.len(), 2);
    for operation in &parent.operations {
        assert!(
            operation.translation.iter().all(|value| *value == 0.0),
            "the parent has a non-zero translation"
        );
    }
}

/// Round 24 pinned count: of the 399 rows in the stratified subset with no
/// candidate, this many are refuted exactly by one of the four conditions of an
/// embedding.  The complete shift solution set decides every candidate, so this
/// now equals the whole zero-row count and the undecided residual is 0 (round
/// 23 pinned 399 of 406, with 7 left open by the finite walk; round 22's
/// withdrawn per-normal anchor test pinned 276).
const CERTIFIED_ABSENCE_ON_THE_SUBSET: usize = 399;

/// The pinned census of the corrected embedding search on a deterministic
/// stratified subset: every 20th record of each parent plus the witnesses.
///
/// Two facts matter and both are pinned here: the family resolves the witnesses
/// a basis-only family missed (through the identity map), and no record in the
/// sample has a unique embedding that realises the record's lattice -- the
/// candidate family is not a decision procedure for subgrouphood, even though
/// (round 24) membership *in the family* is now decided exactly.
///
/// Round 15n then found that the consistency check placed only the **unitary**
/// classes (46% of the entries failed a checker that also places the
/// antiunitary ones); with that gap closed the counts were 402 / 416 and the
/// witness tuples were unchanged.
///
/// Round 15m corrected three things behind those numbers: the
/// `realises_record_lattice` flag tested the record's *column* lattice, the
/// magnetic lattice wrongly included the tabulated cell's unit vectors (so a
/// centred group could never realise), and the shift enumeration truncated a
/// coefficient range.  The replacement walked the affine set with the step set
/// `{+-1, +-1/2}`, which was a **finite sample** of a continuum whenever the
/// kernel is non-trivial -- not a closure proof; round 24 then replaced that
/// walk with the exact complete solution set.  With those fixes the witnesses
/// land on (24,24), 32 and 96 -- exactly what the reviewer predicted
/// independently.
///
/// Round 20 added the missing containment gate (the image of the full source
/// translation lattice) and turned the flag into lattice **equality**.  Numbers
/// moved to **406** rows with no embedding and **446** with no realising one (was
/// 402 / 416), the returned entries dropped from 150,003 to **137,875**, 40 rows
/// kept entries but none that realised the record, and the six witness tuples
/// were unchanged.  Every returned entry is re-checked through the public
/// verifier below, so the gate and the search cannot drift apart silently.
///
/// Round 23 made the empty result decidable through **three** exact necessary
/// conditions of an embedding: the image of the full source translation lattice
/// must lie in the parent's lattice; every class rotation must stay in the
/// parent's rotation set; and every class must be **placeable**, i.e. the joint
/// congruence `(I - R_p) delta = tau_p - A tau` must be solvable modulo the
/// parent lattice for some parent operation with the mapped rotation.  Any
/// candidate refuted by one of them can never embed, so a row is certified when
/// every one of its candidates is refuted.  Round 22's per-normal anchor test
/// was withdrawn in round 22b -- it never ran on its certified rows and was not
/// the joint congruence; the replacement (round 23a) counts antiunitary classes
/// too.
///
/// Round 24 replaced the finite `{+-1, +-1/2}` shift walk with the **complete**
/// solution set of each surviving candidate (a finite union of cosets of the
/// common continuous kernel), which adds the fourth exact condition -- the class
/// congruences must be jointly consistent -- and removes the residual
/// altogether: on this subset **399** of 399 zero rows are certified (8,116 of
/// 8,116 in the full corpus) and the undecided count is 0, so an empty result is
/// now a decision for the family.  The counting contract changed with it:
/// `entries` counts `(map, discrete component)` pairs (one representative per
/// coset of the kernel) instead of the walk's arbitrary sample points, so the
/// entry and pair totals moved even though the witness tuples did not.
///
/// The counts are pinned, and every row with embeddings asserts that the verdict
/// did **not** certify an absence.  SG 3 UNI 24 is the headline: all 240 members
/// of the family are refuted exactly, so the family's empty result is certified
/// -- and the companion test below still exhibits a non-family map that embeds
/// the record, so this is not a proof that the record is not a subgroup.
#[test]
fn the_embedding_census_on_the_stratified_subset_is_pinned() {
    let witnesses: [(u8, usize); 6] = [
        (167, 1333),
        (142, 1221),
        (88, 741),
        (227, 1630),
        (3, 24),
        (1, 1),
    ];
    let mut histogram: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
    let mut realising: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
    let mut unique_realising = 0usize;
    // Round 23: of the rows with no candidate, how many the three exact
    // necessary conditions settle and how many stay open.
    let mut certified_absence = 0usize;
    let mut unsettled_absence = 0usize;
    let mut rows = 0usize;
    let mut entries = 0usize;
    // Map rows as a hashable key; the type alias keeps clippy's complexity lint
    // quiet (the tuple itself is exactly the map's three rows).
    type MapRows = ([Rat; 3], [Rat; 3], [Rat; 3]);
    let mut distinct_pairs: std::collections::HashSet<(u8, usize, MapRows, [Rat; 3])> =
        std::collections::HashSet::new();
    let mut distinct_maps: std::collections::HashSet<(u8, usize, MapRows)> =
        std::collections::HashSet::new();
    let mut witness_counts: std::collections::BTreeMap<(u8, usize), (usize, usize)> =
        std::collections::BTreeMap::new();
    for sg in 1..=230u8 {
        for (index, row) in query::magnetic_isotropy_subgroups_of(sg)
            .into_iter()
            .enumerate()
        {
            let record = row.subgroup;
            let is_witness = witnesses
                .iter()
                .any(|(wsg, wuni)| *wsg == sg && *wuni == record.mag_sg);
            if !is_witness && (usize::from(sg) + index) % 20 != 0 {
                continue;
            }
            rows += 1;
            let geometry = geometry_of(sg, &record).expect("geometry");
            let set = magnetic_operations(record.mag_sg).expect("operations");
            let verdict = parent_embedding_verdict(&geometry, &set).expect("search");
            if verdict.embeddings.is_empty() {
                if verdict.absence_is_certified {
                    certified_absence += 1;
                } else {
                    unsettled_absence += 1;
                }
            } else if verdict.absence_is_certified {
                panic!(
                    "sg {sg} UNI {}: a row with embeddings cannot certify their absence",
                    record.mag_sg
                );
            }
            let embeddings = verdict.embeddings;
            entries += embeddings.len();
            *histogram.entry(embeddings.len()).or_insert(0) += 1;
            let good = embeddings
                .iter()
                .filter(|entry| entry.realises_record_lattice)
                .count();
            *realising.entry(good).or_insert(0) += 1;
            if good == 1 {
                unique_realising += 1;
            }
            // Re-derivation of the flag through the public API (round 20):
            // forward inclusion of the full source lattice into the record's row
            // lattice via the matrix inverse, plus the reverse inclusion of the
            // record's rows into the source lattice.  The reverse half uses
            // `canonical_translation`, the same primitive the production
            // `in_lattice` wraps -- it is an independent re-derivation of the
            // expression, not of the membership primitive.  This is what gives
            // the census teeth for the row/column convention of the production
            // call site and for a forward-only regression.
            let record_lattice = basis_in_parent_conventional(&geometry).expect("record lattice");
            let record_inverse_transposed = Mat3R::new(record_lattice)
                .inverse()
                .expect("non-singular record lattice")
                .transpose();
            let source_generators = full_translation_lattice(&set);
            let source_lattice = translation_lattice(&set);
            let row_membership = |vector: [Rat; 3]| {
                record_inverse_transposed
                    .checked_mul_vector(&Vec3R::new(vector))
                    .map(|image| image.as_array().iter().all(|value| value.is_integer()))
                    .unwrap_or(false)
            };
            // Round 20: the search and the public verifier must agree on every
            // returned pair -- the containment gate lives in both.
            for entry in &embeddings {
                let map_key = (*entry.map.row(0), *entry.map.row(1), *entry.map.row(2));
                distinct_pairs.insert((sg, index, map_key, entry.shift));
                distinct_maps.insert((sg, index, map_key));
                verify_embedding(&geometry, &set, entry.map, entry.shift).unwrap_or_else(|error| {
                    panic!(
                        "sg {sg} UNI {}: returned map {:?} fails the public verifier: {error}",
                        record.mag_sg, entry.map
                    )
                });
                let forward = source_generators.iter().all(|generator| {
                    entry
                        .map
                        .checked_mul_vector(&Vec3R::new(*generator))
                        .map(|image| row_membership(*image.as_array()))
                        .unwrap_or(false)
                });
                let reverse = entry
                    .map
                    .inverse()
                    .map(|inverse| {
                        record_lattice.iter().all(|row| {
                            inverse
                                .checked_mul_vector(&Vec3R::new(*row))
                                .map(|pullback| {
                                    canonical_translation(*pullback.as_array(), &source_lattice)
                                        .iter()
                                        .all(|value| value.is_zero())
                                })
                                .unwrap_or(false)
                        })
                    })
                    .unwrap_or(false);
                assert_eq!(
                    entry.realises_record_lattice,
                    forward && reverse,
                    "sg {sg} UNI {} map {:?}: the flag must be lattice equality",
                    record.mag_sg,
                    entry.map
                );
            }
            if is_witness {
                // Keep the FIRST occurrence: a witness may appear in several
                // rows (different direction labels) with different bases.
                witness_counts
                    .entry((sg, record.mag_sg))
                    .or_insert((embeddings.len(), good));
            }
        }
    }
    // Printed only when an assertion below fails (cargo captures stdout on
    // success): the numbers a re-pin needs, in one place.
    println!(
        "census rows={rows} entries={entries} pairs={} maps={} zero={:?} realising_zero={:?} \
         unique_realising={unique_realising} certified={certified_absence} \
         unsettled={unsettled_absence} witnesses={witness_counts:?}",
        distinct_pairs.len(),
        distinct_maps.len(),
        histogram.get(&0),
        realising.get(&0),
    );
    assert_eq!(rows, 851);
    assert_eq!(
        entries, 131_116,
        "returned entries: round 24 counts one per (map, discrete component)"
    );
    // The same map can be reached through several `kind`s, so the entry count is
    // not a pair count.  Both are pinned: 137,875 entries, 76,291 distinct
    // `(record, map, shift)`, 17,517 distinct `(record, map)`.  (The older
    // reviewer aggregation keyed on `(parent_sg, UNI, map, shift, kind)` and
    // merged distinct records, which is why its 3,696 duplicates were neither
    // the entry count nor the pair count.)
    assert_eq!(
        distinct_pairs.len(),
        73_432,
        "distinct (record, map, shift)"
    );
    assert_eq!(distinct_maps.len(), 18_684, "distinct (record, map)");
    assert_eq!(histogram.get(&0), Some(&399), "records with no embedding");
    assert_eq!(
        realising.get(&0),
        Some(&436),
        "records with no realising embedding"
    );
    assert_eq!(
        unique_realising, 0,
        "no record has a unique realising embedding"
    );
    // Round 23: the three-condition refutation must partition the zero rows, and
    // every row that has embeddings must have refused to certify (checked in the
    // loop).
    assert_eq!(
        certified_absence + unsettled_absence,
        histogram.get(&0).copied().unwrap_or(0),
        "certified and unsettled zero rows must partition the zero rows"
    );
    assert_eq!(
        certified_absence, CERTIFIED_ABSENCE_ON_THE_SUBSET,
        "zero rows settled by the exact refutation"
    );
    assert_eq!(
        unsettled_absence, 0,
        "round 24's complete solution set decides every candidate, so no row stays open"
    );
    assert_eq!(
        witness_counts,
        std::collections::BTreeMap::from([
            ((1u8, 1usize), (240usize, 240usize)),
            ((3, 24), (0, 0)),
            ((88, 741), (32, 32)),
            ((142, 1221), (96, 96)),
            ((167, 1333), (24, 24)),
            ((227, 1630), (288, 288)),
        ]),
        "witness embedding counts (total, realising the record lattice)"
    );
}

/// The round-15m witness: SG 3 UNI 24 has **no** member of the search family,
/// but a small explicit embedding does exist -- `A = [[0,0,-2],[0,-2,0],
/// [-2,0,-2]]` with a zero shift.  Both facts are pinned here, using the public
/// verifier for the map found by the independent sweep.
/// Round 15n's F1: the consistency check must place the **antiunitary** classes
/// too.  SG 6 UNI 29 with the identity map and a zero shift is the minimal
/// counterexample the reviewer gave -- its primed class `(m|(0,1/2,0), T)`
/// maps to `(m|(0,1/2,0))`, which is not a parent operation.
#[test]
fn the_antiunitary_classes_are_checked_too() {
    let record = query::magnetic_isotropy_subgroups_of(6)
        .into_iter()
        .find(|row| row.subgroup.mag_sg == 29)
        .map(|row| row.subgroup)
        .expect("record");
    let geometry = geometry_of(6, &record).expect("geometry");
    let set = magnetic_operations(29).expect("operations");
    assert!(
        set.classes.iter().any(|operation| operation.time_reversal),
        "the witness has an antiunitary class"
    );
    assert!(
        verify_embedding(
            &geometry,
            &set,
            Mat3R::identity(),
            [Rat::ZERO, Rat::ZERO, Rat::ZERO]
        )
        .is_err(),
        "the identity map must not place the antiunitary class"
    );
    let embeddings = embed_in_parent_conventional(&geometry, &set).expect("search");
    assert!(
        !embeddings.iter().any(|entry| {
            entry.map == Mat3R::identity() && entry.shift.iter().all(|value| value.is_zero())
        }),
        "a pair that fails the antiunitary check must not be returned"
    );
    // Every returned pair must pass the public verifier, antiunitary included.
    for entry in &embeddings {
        verify_embedding(&geometry, &set, entry.map, entry.shift)
            .expect("every returned pair is a full embedding");
    }
}

#[test]
fn the_unresolved_witness_does_have_a_small_embedding() {
    let record = query::magnetic_isotropy_subgroups_of(3)
        .into_iter()
        .find(|row| row.subgroup.mag_sg == 24)
        .map(|row| row.subgroup)
        .expect("record");
    let geometry = geometry_of(3, &record).expect("geometry");
    let set = magnetic_operations(24).expect("operations");
    assert!(
        embed_in_parent_conventional(&geometry, &set)
            .expect("search")
            .is_empty(),
        "the family does not contain this witness"
    );
    let map = Mat3R::from_ints([[0, 0, -2], [0, -2, 0], [-2, 0, -2]]);
    verify_embedding(&geometry, &set, map, [Rat::ZERO, Rat::ZERO, Rat::ZERO])
        .expect("the small explicit embedding is valid");
    // A map that is deliberately wrong must be rejected by the same verifier:
    // the identity does not work here, because the parent's twofold axis has a
    // zero translation while the record's twofold class carries (1/2,1/2,0).
    assert!(
        verify_embedding(
            &geometry,
            &set,
            Mat3R::identity(),
            [Rat::ZERO, Rat::ZERO, Rat::ZERO]
        )
        .is_err(),
        "the identity map must not verify for this record"
    );
}

#[test]
fn the_conversion_refuses_an_unknown_parent() {
    let record = first_record(14);
    let geometry = geometry_of(14, &record).expect("geometry");
    let mut unknown = geometry.clone();
    unknown.parent_sg = 0;
    assert_eq!(
        basis_in_parent_conventional(&unknown),
        Err(GeometryError::ParentPrimitiveBasisUnavailable { sg: 0 })
    );
    assert_eq!(
        origin_in_parent_conventional(&unknown),
        Err(GeometryError::ParentPrimitiveBasisUnavailable { sg: 0 })
    );
    assert_eq!(
        conventional_lattice_index(&unknown),
        Err(GeometryError::ParentPrimitiveBasisUnavailable { sg: 0 })
    );
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
    assert_eq!(
        failing_parents, 224,
        "parents with a unitary class outside the parent"
    );
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
    // Round 20: the full lattice is `Z^3` plus the listed pure translations, and
    // the image is taken over exactly those generators.
    let listed = translation_lattice(&set);
    let source = full_translation_lattice(&set);
    assert!(!listed.is_empty(), "UNI 20 lists a centring translation");
    assert_eq!(
        source.len(),
        3 + listed.len(),
        "Z^3 plus the listed generators"
    );
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
        full_translation_lattice(&identity_set)
    );
}

/// The full source lattice is `Z^3` plus the listed pure translations, and the
/// doc claim that an antiunitary pure translation cannot exist is checked
/// against the corpus rather than assumed (round 20).
#[test]
fn the_full_translation_lattice_contains_the_cell_translations() {
    let first = first_record(1);
    let set = magnetic_operations(first.mag_sg).expect("operations");
    assert!(
        translation_lattice(&set).is_empty(),
        "UNI 1 lists no pure translation"
    );
    let full = full_translation_lattice(&set);
    assert_eq!(
        full,
        vec![
            [Rat::ONE, Rat::ZERO, Rat::ZERO],
            [Rat::ZERO, Rat::ONE, Rat::ZERO],
            [Rat::ZERO, Rat::ZERO, Rat::ONE],
        ],
        "the cell translations are part of the lattice even when nothing is listed"
    );

    let centred = magnetic_operations(20).expect("UNI 20");
    let full = full_translation_lattice(&centred);
    assert_eq!(full.len(), 4, "Z^3 plus the one centring translation");
    assert!(
        full[3..].iter().all(|vector| vector
            .iter()
            .any(|value| !value.is_zero() && !value.is_integer())),
        "the fourth generator is the non-integral centring translation"
    );

    // Every referenced UNI is scanned: a set with an antiunitary pure
    // translation `(E | tau, T)` is the reason `translation_lattice` alone is
    // not the full lattice (its square is a unitary lattice translation).
    let mut referenced = std::collections::BTreeSet::new();
    for sg in 1..=230u8 {
        for row in query::magnetic_isotropy_subgroups_of(sg) {
            referenced.insert(row.subgroup.mag_sg);
        }
    }
    let mut antiunitary_translations = 0usize;
    for uni in referenced {
        let set = magnetic_operations(uni).expect("operations");
        antiunitary_translations += set
            .operations
            .iter()
            .filter(|operation| {
                operation.time_reversal && operation.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            })
            .count();
    }
    assert!(
        antiunitary_translations > 0,
        "the corpus does contain antiunitary pure translations"
    );
}

/// Round 20's containment counterexamples.  The per-class equations are
/// satisfied in both cases; the image of the full translation lattice is not.
#[test]
fn the_full_source_lattice_gate_rejects_the_round_20_counterexamples() {
    // (a) SG 1 UNI 1 has a single class, so every invertible map passes the
    //     class equation -- but `A = I/2` sends e_1 to (1/2, 0, 0), which is not
    //     a translation of the parent.
    let record = first_record(1);
    assert_eq!(record.mag_sg, 1);
    let geometry = geometry_of(1, &record).expect("geometry");
    let set = magnetic_operations(record.mag_sg).expect("operations");
    assert_eq!(set.classes.len(), 1, "the witness has a single class");
    let half = Rat::new(1, 2).expect("1/2");
    let half_identity = Mat3R::new([
        [half, Rat::ZERO, Rat::ZERO],
        [Rat::ZERO, half, Rat::ZERO],
        [Rat::ZERO, Rat::ZERO, half],
    ]);
    assert_eq!(
        verify_embedding(&geometry, &set, half_identity, [Rat::ZERO; 3]),
        Err(GeometryError::SourceLatticeNotInParent { uni: record.mag_sg })
    );
    // Positive control: the identity still verifies, so the gate is not a
    // blanket rejection.
    verify_embedding(&geometry, &set, Mat3R::identity(), [Rat::ZERO; 3])
        .expect("the identity embeds");

    // (b) SG 143 M1/S1 UNI 1: the family used to return the rational dual map
    //     `A = W^-1` with `realises_record_lattice = true`.
    let record = query::magnetic_isotropy_subgroups_of(143)
        .into_iter()
        .find(|row| {
            row.subgroup.mag_sg == 1 && row.subgroup.basis == [[0, 0, 1], [0, 2, 0], [-2, 0, 0]]
        })
        .map(|row| row.subgroup)
        .expect("the round-20 record");
    let geometry = geometry_of(143, &record).expect("geometry");
    let set = magnetic_operations(1).expect("UNI 1");
    let culprit = Mat3R::new(geometry.basis)
        .inverse()
        .expect("invertible basis");
    assert_eq!(
        verify_embedding(&geometry, &set, culprit, [Rat::ZERO; 3]),
        Err(GeometryError::SourceLatticeNotInParent { uni: 1 })
    );
    let image = culprit
        .checked_mul_vector(&Vec3R::new([Rat::ZERO, Rat::ONE, Rat::ZERO]))
        .expect("product");
    assert!(
        image.as_array().iter().any(|value| !value.is_integer()),
        "A e_2 leaves the parent's Z^3"
    );
    let embeddings = embed_in_parent_conventional(&geometry, &set).expect("search");
    assert!(
        !embeddings.iter().any(|entry| entry.map == culprit),
        "the search must not return a map whose source lattice leaves the parent"
    );
}

/// `realises_record_lattice` is lattice **equality** `A L_M = L_record`: the
/// forward inclusion alone was vacuous for a set with no listed pure
/// translation, and the reverse inclusion was never tested (round 20).
///
/// Teeth map (round-20 review, F3): this test pins the gate's error variant, the
/// `Z^3` + listed composition and the contained-but-not-equal corpus and
/// synthetic cases.  It does **not** discriminate the reverse inclusion or the
/// row/column convention -- those are pinned by
/// `irrep::magnetic_embedding::tests::the_lattice_realisation_needs_both_inclusions`
/// and by the census re-derivation in `the_embedding_census_on_the_stratified_subset_is_pinned`
/// respectively.
#[test]
fn the_record_lattice_flag_is_equality() {
    let zero = [Rat::ZERO; 3];
    let synthetic = |parent_sg: u8, uni: usize, basis: [[i32; 3]; 3], index: i128| {
        let mut exact = [[Rat::ZERO; 3]; 3];
        for (row, values) in basis.iter().enumerate() {
            for (column, value) in values.iter().enumerate() {
                exact[row][column] = Rat::from_integer(i128::from(*value));
            }
        }
        MagneticGeometry {
            parent_sg,
            uni,
            basis: exact,
            origin: zero,
            index,
        }
    };
    let entry_for = |geometry: &MagneticGeometry, map: &Mat3R| {
        let set = magnetic_operations(geometry.uni).expect("operations");
        let embeddings = embed_in_parent_conventional(geometry, &set).expect("search");
        let entry = embeddings
            .iter()
            .find(|entry| entry.map == *map && entry.shift == zero)
            .unwrap_or_else(|| panic!("map {map:?} is not returned"));
        (entry.realises_record_lattice, entry.kind)
    };

    // The identity on a primitive record realises `Z^3` exactly.
    let cell = synthetic(1, 1, [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 1);
    assert_eq!(
        entry_for(&cell, &Mat3R::identity()),
        (true, EmbeddingKind::Plain)
    );

    // The same group with a doubled record lattice: containment still holds
    // (`Z^3` is mapped into `Z^3`), but `A L_M = Z^3 != 2Z^3`, so the flag is
    // false -- the case where "contained" and "realises the record" differ.
    let doubled = synthetic(1, 1, [[2, 0, 0], [0, 2, 0], [0, 0, 2]], 8);
    assert_eq!(
        entry_for(&doubled, &Mat3R::identity()),
        (false, EmbeddingKind::Plain)
    );
    let set = magnetic_operations(1).expect("operations");
    verify_embedding(&doubled, &set, Mat3R::identity(), zero)
        .expect("the doubled record is still contained in the parent");

    // A centred source: the identity realises `L_C` when the record keeps the
    // parent's lattice, and fails when the record lattice is a proper sublattice
    // of it, although the whole group is contained either way.
    let centred_cell = synthetic(5, 20, [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 1);
    assert_eq!(
        entry_for(&centred_cell, &Mat3R::identity()),
        (true, EmbeddingKind::Plain)
    );
    let centred_doubled = synthetic(5, 20, [[2, 0, 0], [0, 2, 0], [0, 0, 2]], 8);
    assert_eq!(
        entry_for(&centred_doubled, &Mat3R::identity()),
        (false, EmbeddingKind::Plain)
    );
    let c2 = magnetic_operations(20).expect("UNI 20");
    verify_embedding(&centred_doubled, &c2, Mat3R::identity(), zero)
        .expect("the whole C2 group is contained despite the false flag");

    // The **listed** half of the source lattice is checked too, not just `Z^3`:
    // an integer shear keeps the cell translations inside the parent but sends
    // the C centring `(1/2,1/2,0)` to `(1,1/2,0)`, which is not a C lattice
    // vector.
    let shear = Mat3R::from_ints([[1, 1, 0], [0, 1, 0], [0, 0, 1]]);
    assert_eq!(
        verify_embedding(&centred_cell, &c2, shear, zero),
        Err(GeometryError::SourceLatticeNotInParent { uni: 20 })
    );

    // Corpus witness: SG 1 Z1/P1 UNI 3.  The integer map is contained (it maps
    // the classes onto P1 operations) but its image lattice `Z e_1 + 2Z e_2 +
    // Z e_3` is neither contained in nor containing the record lattice
    // `Z e_1 + Z e_2 + 2Z e_3` -- same index, different lattice, flag false.  The
    // old listed-translation-only test called it true.
    let record = query::magnetic_isotropy_subgroups_of(1)
        .into_iter()
        .find(|row| {
            row.subgroup.mag_sg == 3 && row.subgroup.basis == [[0, 1, 0], [-1, 0, 0], [0, 0, 2]]
        })
        .map(|row| row.subgroup)
        .expect("the Z1/P1 record");
    let geometry = geometry_of(1, &record).expect("geometry");
    let set = magnetic_operations(3).expect("UNI 3");
    let map = Mat3R::from_ints([[0, -1, 0], [0, 0, 2], [1, 0, 0]]);
    let entry = embed_in_parent_conventional(&geometry, &set)
        .expect("search")
        .into_iter()
        .find(|entry| entry.map == map && entry.shift == zero)
        .expect("the integer map is returned");
    assert!(
        !entry.realises_record_lattice,
        "A e_1 = e_3 is not in the record lattice, so the flag cannot be true"
    );
    verify_embedding(&geometry, &set, map, zero).expect("the map is still a containment");
}

/// Round 24's continuous boundary on real corpus data: for SG 1 UNI 1 (`P1`)
/// the identity map places the single identity class for **every** shift, so the
/// solution set is the whole torus -- one component whose kernel is three
/// dimensional.  The API must say exactly that instead of sampling the
/// continuum: one representative (the zero shift) plus a three-vector kernel.
///
/// This is the shape the withdrawn `{+-1, +-1/2}` walk reported as 1,920 sampled
/// points (the census pins 240 entries for this witness now: one per
/// `(map, component, kind)`).
#[test]
fn the_continuous_kernel_is_reported_not_sampled() {
    let record = first_record(1);
    assert_eq!(record.mag_sg, 1, "the SG 1 witness is UNI 1");
    let geometry = geometry_of(1, &record).expect("geometry");
    let set = magnetic_operations(1).expect("operations");
    let solution =
        embedding_shift_components(&geometry, &set, Mat3R::identity()).expect("components");
    assert_eq!(
        solution.representatives,
        vec![[Rat::ZERO, Rat::ZERO, Rat::ZERO]],
        "the whole torus is one component, represented once"
    );
    assert_eq!(solution.kernel.len(), 3, "the kernel is three dimensional");
    // The identity map is returned by the search, and every entry it returns for
    // it is the zero shift (one per kind that builds the identity map).
    let embeddings = embed_in_parent_conventional(&geometry, &set).expect("search");
    let identity: Vec<_> = embeddings
        .iter()
        .filter(|entry| entry.map == Mat3R::identity())
        .collect();
    assert!(!identity.is_empty(), "the identity map embeds SG 1 UNI 1");
    for entry in identity {
        assert_eq!(entry.shift, [Rat::ZERO; 3]);
        verify_embedding(&geometry, &set, entry.map, entry.shift).expect("verified");
    }
}

/// Round 23's review (F1) found two in-family `(map, shift)` pairs that the
/// public `verify_embedding` accepts while the shipped `{+-1, +-1/2}` shift walk
/// returned nothing: SG 56 UNI 119 with `map = diag(1,1,2)` and shift
/// `(-1/4,-1/4,-1/4)`, and SG 89 UNI 111 with
/// `map = [[1,0,1],[1,0,-1],[0,1,0]]` and shift `(-1/4,-1/4,0)`.  The walk could
/// not reach a lattice-shifted branch (the first) or a quarter-denominator
/// component (both).
///
/// Round 24's complete solver must return both maps, with the claimed shift in
/// one returned component, and every returned entry must pass the public
/// verifier.  The component counts are pinned: 8 and 4.  For these two records
/// the parent lattice is `Z^3` and the shift solution set is discrete, so "the
/// claimed shift is in the returned set" is exactly "it differs from a
/// representative by an integer vector" -- both facts are asserted, not assumed.
#[test]
fn the_public_api_captures_the_missed_quarter_shifts() {
    for (parent_sg, uni, rows, shift, components) in [
        (
            56u8,
            119usize,
            [[1, 0, 0], [0, 1, 0], [0, 0, 2]],
            [-1i128, -1, -1],
            8usize,
        ),
        (89, 111, [[1, 0, 1], [1, 0, -1], [0, 1, 0]], [-1, -1, 0], 4),
    ] {
        let record = query::magnetic_isotropy_subgroups_of(parent_sg)
            .into_iter()
            .find(|row| row.subgroup.mag_sg == uni)
            .map(|row| row.subgroup)
            .expect("record");
        let geometry = geometry_of(parent_sg, &record).expect("geometry");
        let set = magnetic_operations(uni).expect("operations");
        let map = Mat3R::from_ints(rows);
        let claimed = [
            Rat::new(shift[0], 4).expect("quarter"),
            Rat::new(shift[1], 4).expect("quarter"),
            Rat::new(shift[2], 4).expect("quarter"),
        ];
        // The claimed pair is a genuine embedding, independently of the search:
        // this is the evidence the round-23 review produced through the same
        // public verifier.
        verify_embedding(&geometry, &set, map, claimed)
            .unwrap_or_else(|error| panic!("SG {parent_sg} UNI {uni}: {error}"));

        // Round 24's complete solution set of the same map.
        let solution = embedding_shift_components(&geometry, &set, map).expect("components");
        assert_eq!(
            solution.representatives.len(),
            components,
            "SG {parent_sg} UNI {uni}: the discrete component count"
        );
        assert!(
            solution.kernel.is_empty(),
            "SG {parent_sg} UNI {uni}: this witness's shift set is discrete"
        );
        assert!(
            solution.representatives.iter().any(|representative| {
                (0..3).all(|axis| {
                    (claimed[axis]
                        .checked_sub(representative[axis])
                        .expect("exact"))
                    .is_integer()
                })
            }),
            "SG {parent_sg} UNI {uni}: the missed shift is not in any returned component"
        );

        // The search must return the map, one entry per component per `kind`
        // that builds the same map; both witnesses are reached through the
        // record's primitive and conventional bases, which coincide here.
        let embeddings = embed_in_parent_conventional(&geometry, &set).expect("search");
        let same_map: Vec<_> = embeddings.iter().filter(|entry| entry.map == map).collect();
        assert_eq!(
            same_map.len(),
            2 * components,
            "SG {parent_sg} UNI {uni}: entries with the missed map"
        );
        for entry in &same_map {
            verify_embedding(&geometry, &set, entry.map, entry.shift).unwrap_or_else(|error| {
                panic!("SG {parent_sg} UNI {uni}: a returned entry fails the verifier: {error}")
            });
        }
        let distinct: std::collections::BTreeSet<[(i128, i128); 3]> = same_map
            .iter()
            .map(|entry| {
                entry
                    .shift
                    .map(|value| (value.numerator(), value.denominator()))
            })
            .collect();
        assert_eq!(
            distinct.len(),
            components,
            "SG {parent_sg} UNI {uni}: the returned shifts are the components"
        );
    }
}

/// Round 24's arithmetic boundary at the public door, in both directions.
///
/// The door's empty solution set is documented as a **proof** that no shift
/// makes the map an embedding, so an exact-arithmetic limit must never be folded
/// into it: the huge-rational map below passes the map-only containment gate
/// (its entries are integers) and then overflows while inverting, and the answer
/// must be an explicit error.  A singular map, by contrast, is refuted for good
/// -- no shift can make it a change of basis -- and that empty set is genuine.
#[test]
fn the_public_door_distinguishes_limits_from_refutations() {
    let record = first_record(1);
    let geometry = geometry_of(1, &record).expect("geometry");
    let set = magnetic_operations(1).expect("operations");

    let huge = Rat::from_integer(i128::MAX / 2);
    let enormous = Mat3R::new([
        [huge, Rat::ZERO, Rat::ZERO],
        [Rat::ZERO, huge, Rat::ZERO],
        [Rat::ZERO, Rat::ZERO, Rat::ONE],
    ]);
    match embedding_shift_components(&geometry, &set, enormous) {
        Err(GeometryError::Arithmetic) | Err(GeometryError::ShiftEnumerationUnavailable { .. }) => {
        }
        other => panic!(
            "an exact-arithmetic limit must be unknown, not an empty solution set: {other:?}"
        ),
    }

    let singular = Mat3R::new([
        [Rat::ONE, Rat::ZERO, Rat::ZERO],
        [Rat::ZERO, Rat::ONE, Rat::ZERO],
        [Rat::ZERO, Rat::ZERO, Rat::ZERO],
    ]);
    let components = embedding_shift_components(&geometry, &set, singular).expect("decided");
    assert!(
        components.representatives.is_empty() && components.kernel.is_empty(),
        "a singular map cannot be a change of basis, so its solution set is empty"
    );
    assert!(
        verify_embedding(&geometry, &set, singular, [Rat::ZERO; 3]).is_err(),
        "the public verifier agrees that the singular map is not an embedding"
    );
}

/// The round-24 invariant at the **production** verdict
/// ([`parent_embedding_verdict`]): an exact-arithmetic limit leaves candidates
/// undecided, so a row can never be certified on a system that was not actually
/// decided.  The two controls below are genuine refutations -- a source-lattice
/// image that cannot land in the parent lattice, and a rotation that is not a
/// parent rotation -- and must still be counted as exclusions, not as undecided.
#[test]
fn the_verdict_never_certifies_a_row_with_an_arithmetic_limit() {
    let unit_rotation = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let identity_operation = MagneticOperation {
        rotation: unit_rotation,
        translation: [Rat::ZERO; 3],
        time_reversal: false,
    };
    let integer_geometry = MagneticGeometry {
        parent_sg: 1,
        uni: 999,
        basis: [
            [Rat::ONE, Rat::ZERO, Rat::ZERO],
            [Rat::ZERO, Rat::ONE, Rat::ZERO],
            [Rat::ZERO, Rat::ZERO, Rat::ONE],
        ],
        origin: [Rat::ZERO; 3],
        index: 1,
    };
    let real_set = magnetic_operations(1).expect("operations");

    // Synthetic geometry, `parent_sg = 1` so the frame conversion is the
    // identity and the only source of failure is arithmetic on the record's own
    // basis: `diag(huge, huge, 1)` has a determinant whose square overflows.
    let huge = Rat::from_integer(i128::MAX / 2);
    let overflowing = MagneticGeometry {
        parent_sg: 1,
        uni: 1,
        basis: [
            [huge, Rat::ZERO, Rat::ZERO],
            [Rat::ZERO, huge, Rat::ZERO],
            [Rat::ZERO, Rat::ZERO, Rat::ONE],
        ],
        origin: [Rat::ZERO; 3],
        index: 1,
    };
    match parent_embedding_verdict(&overflowing, &real_set) {
        // Either the prologue reports the limit directly, or the search carries
        // on and books the undecided candidates: both are "unknown".  What is
        // forbidden is a verdict that certifies the row.
        Err(GeometryError::Arithmetic) => {}
        Err(other) => panic!("an arithmetic limit must be reported as one: {other}"),
        Ok(verdict) => {
            assert!(
                !verdict.absence_is_certified,
                "an arithmetic limit must never certify a row"
            );
            assert!(
                verdict.candidates_without_witness > 0,
                "the limit has to show up as the undecided residual"
            );
        }
    }

    // Control 1: a genuine source-lattice refutation that also kills the `Plain`
    // (identity) base.  The synthetic set is C centred: its full source lattice
    // is `Z^3 + (1/2,1/2,0)Z`, and for the `P1` parent (`L_P = Z^3`) a signed
    // permutation sends that centring vector to `(+-1/2,+-1/2,0)`, never into
    // `Z^3`.  Every candidate is refuted by the containment gate, so the row is
    // certified with no undecided residual.
    let centred = MagneticOperationSet {
        uni: 998,
        operations: vec![
            identity_operation,
            MagneticOperation {
                rotation: unit_rotation,
                translation: [
                    Rat::new(1, 2).expect("half"),
                    Rat::new(1, 2).expect("half"),
                    Rat::ZERO,
                ],
                time_reversal: false,
            },
        ],
        classes: vec![identity_operation],
        unitary_count: 2,
        antiunitary_count: 0,
        unitary_classes: 1,
        antiunitary_classes: 0,
    };
    let verdict = parent_embedding_verdict(&integer_geometry, &centred).expect("decided");
    assert!(
        verdict.embeddings.is_empty(),
        "the centring vector cannot be pushed into Z^3 by a signed permutation"
    );
    assert_eq!(
        verdict.candidates_without_witness, 0,
        "no arithmetic limit here"
    );
    assert!(
        verdict.excluded_by_source_lattice > 0,
        "the centring is refuted by the containment gate"
    );
    assert!(verdict.absence_is_certified);

    // Control 2: a genuine rotation refutation.  `P1` has only the identity
    // rotation, and `A (-I) A^-1 = -I` for every map, so a set with an inversion
    // class is refuted by the rotation screen for every candidate.
    let inversion = MagneticOperation {
        rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
        translation: [Rat::ZERO; 3],
        time_reversal: false,
    };
    let inverting = MagneticOperationSet {
        uni: 997,
        operations: vec![identity_operation, inversion],
        classes: vec![identity_operation, inversion],
        unitary_count: 2,
        antiunitary_count: 0,
        unitary_classes: 2,
        antiunitary_classes: 0,
    };
    let verdict = parent_embedding_verdict(&integer_geometry, &inverting).expect("decided");
    assert_eq!(verdict.embeddings, Vec::new());
    assert_eq!(
        verdict.candidates_without_witness, 0,
        "no arithmetic limit here"
    );
    assert!(
        verdict.excluded_by_rotation > 0,
        "the inversion is refuted by the rotation screen"
    );
    assert!(verdict.absence_is_certified);
}

/// The per-class placeability **pre-screen** must not read an exact-arithmetic
/// failure in its right-hand side `tau_p - A tau` as "this class is
/// unplaceable": that control flow used to fall through every matching parent
/// operation, leave `any == false`, and book the map in
/// `excluded_by_anchor_congruence` -- a refutation manufactured from a
/// computation that never finished (round-24 review 696f4ec6).
///
/// The synthetic class below carries `tau = (0,0,-i128::MAX)` under the twofold
/// rotation of `P2_1` (SG 4), whose listed operation has translation
/// `(0,0,1/2)`: `1/2 - (-MAX)` overflows, while the map passes both map screens
/// (the source lattice is `Z^3` because the huge translation belongs to a
/// non-identity rotation, and `A R A^-1 = R`).  The verdict may report the limit
/// directly or leave a positive residual; it must never certify the row from
/// that failed subtraction.
#[test]
fn the_placeability_screen_never_reads_an_arithmetic_limit_as_a_refutation() {
    let huge = Rat::from_integer(-i128::MAX);
    let operations = vec![MagneticOperation {
        rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
        translation: [Rat::ZERO, Rat::ZERO, huge],
        time_reversal: false,
    }];
    let set = MagneticOperationSet {
        uni: 996,
        operations: operations.clone(),
        classes: operations,
        unitary_count: 1,
        antiunitary_count: 0,
        unitary_classes: 1,
        antiunitary_classes: 0,
    };
    let geometry = MagneticGeometry {
        parent_sg: 4,
        uni: 996,
        basis: [
            [Rat::ONE, Rat::ZERO, Rat::ZERO],
            [Rat::ZERO, Rat::ONE, Rat::ZERO],
            [Rat::ZERO, Rat::ZERO, Rat::ONE],
        ],
        origin: [Rat::ZERO; 3],
        index: 1,
    };
    match parent_embedding_verdict(&geometry, &set) {
        Err(GeometryError::Arithmetic) => {}
        Err(other) => panic!("an arithmetic limit must be reported as one: {other}"),
        Ok(verdict) => {
            assert!(
                !verdict.absence_is_certified,
                "a failed right-hand side must not certify the row"
            );
            assert!(
                verdict.candidates_without_witness > 0,
                "the limit has to show up as the undecided residual"
            );
            assert_eq!(
                verdict.excluded_by_anchor_congruence, 0,
                "no map may be refuted by a subtraction that overflowed"
            );
        }
    }
}
