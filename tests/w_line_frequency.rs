//! The `other_wave_vector_subduction` frequencies computed by the engine.
//!
//! The 73 sources behind these rows are little irreps of lines `k = t*v`
//! through Gamma, so they have no discrete wave vector to fold and no character
//! rows in the pinned irrep table.  `src/irrep/w_little_characters_data.rs`
//! carries their frozen little-group character tables, and
//! `line_trivial_content_with_embedding` turns one of them into the pinned
//! frequency.  This test pins the measured behaviour of the `P1`-child family
//! (all 300 pinned rows of the 46 such records) and the ten records whose
//! `DT3`/`DT4` frequencies differ, which is what pins the frozen ordering of
//! that pair.

use cryspglib::irrep::isotropy::{self, IsotropySubgroup};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::{Rat, SubgroupEmbedding, Vec3R};
use cryspglib::irrep::subduction::star::decompose::{
    line_trivial_content_via_blocks, line_trivial_content_with_embedding, official_line_parameter,
    subduce_line_at_parameter,
};
use cryspglib::irrep::line_monodromy::{line_direction, monodromy};
use cryspglib::irrep::w_little_characters_data::{LittleCharacterTable, W_LITTLE_CHARACTERS};
use cryspglib::irrep::LabelConvention;

/// Every isotropy record of one parent space group, as the audit enumerates it.
fn subgroups_of(sg: u8) -> Vec<IsotropySubgroup> {
    let mut out = Vec::new();
    for record in query::irreps_of(sg) {
        if record.spinor || record.subgroups().is_empty() {
            continue;
        }
        let subgroups = isotropy::isotropy_subgroups(sg, record.ml, LabelConvention::Cdml)
            .unwrap_or_else(|error| panic!("SG {sg} {}: {error}", record.ml));
        out.extend(subgroups);
    }
    out
}

/// The engine's frequency for one w row of one record.
fn frequency(subgroup: &IsotropySubgroup, label: &str) -> Result<u32, String> {
    let table = W_LITTLE_CHARACTERS
        .iter()
        .find(|table| {
            usize::from(table.space_group) == usize::from(subgroup.parent_sg)
                && table.label == label
        })
        .ok_or_else(|| format!("no frozen table for SG {} {label}", subgroup.parent_sg))?;
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup)
        .map_err(|error| format!("ordinal {}: {error}", subgroup.ordinal))?;
    line_trivial_content_with_embedding(subgroup, &embedding, table)
        .map_err(|error| format!("ordinal {} {label}: {error}", subgroup.ordinal))
}

#[test]
fn p1_child_records_match_every_pinned_row() {
    let mut checked = 0usize;
    let mut mismatched: Vec<String> = Vec::new();
    for sg in [196u8, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228] {
        for subgroup in subgroups_of(sg) {
            // A `P1` child: the record's own space group number is 1.
            if subgroup.record.sg != 1 {
                continue;
            }
            let Ok(rows) = subgroup.other_wave_vector_subduction() else {
                continue;
            };
            for entry in rows {
                checked += 1;
                match frequency(&subgroup, entry.parent_ml) {
                    Ok(value) if value == u32::from(entry.frequency) => {}
                    other => mismatched.push(format!(
                        "ordinal {} SG {sg} {}: computed {other:?}, pinned {}",
                        subgroup.ordinal, entry.parent_ml, entry.frequency
                    )),
                }
            }
        }
    }
    assert_eq!(
        checked, 300,
        "the P1-child family is 300 pinned rows; got {checked}"
    );
    assert!(
        mismatched.is_empty(),
        "{} mismatches, first few: {:?}",
        mismatched.len(),
        &mismatched[..mismatched.len().min(8)]
    );
}

#[test]
fn asymmetric_dt3_dt4_records_are_pinned() {
    // These ten records are the ones whose `DT3` and `DT4` frequencies differ,
    // so they are the falsifiable prediction of the frozen pair ordering.
    let expected: &[(usize, u8, u16, u16)] = &[
        (10422, 202, 2, 1),
        (10423, 202, 1, 2),
        (10446, 202, 1, 2),
        (10447, 202, 2, 1),
        (10485, 202, 1, 5),
        (10509, 202, 5, 1),
        (11171, 209, 2, 1),
        (11173, 209, 1, 2),
        (11202, 209, 2, 1),
        (11204, 209, 1, 2),
    ];
    for (ordinal, sg, dt3, dt4) in expected {
        let subgroup = subgroups_of(*sg)
            .into_iter()
            .find(|subgroup| subgroup.ordinal == *ordinal)
            .unwrap_or_else(|| panic!("ordinal {ordinal} of SG {sg} is missing"));
        assert_eq!(
            frequency(&subgroup, "DT3"),
            Ok(u32::from(*dt3)),
            "ordinal {ordinal} SG {sg} DT3"
        );
        assert_eq!(
            frequency(&subgroup, "DT4"),
            Ok(u32::from(*dt4)),
            "ordinal {ordinal} SG {sg} DT4"
        );
    }
}

#[test]
fn a_source_of_another_parent_is_rejected() {
    let subgroup = subgroups_of(196)
        .into_iter()
        .find(|subgroup| !subgroup.other_wave_vector_subduction().unwrap().is_empty())
        .expect("SG 196 has w records");
    let foreign = W_LITTLE_CHARACTERS
        .iter()
        .find(|table| table.space_group != subgroup.parent_sg)
        .expect("another parent's table");
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
    assert!(line_trivial_content_with_embedding(&subgroup, &embedding, foreign).is_err());
}

fn find(sg: u8, ordinal: usize) -> IsotropySubgroup {
    subgroups_of(sg)
        .into_iter()
        .find(|subgroup| subgroup.ordinal == ordinal)
        .unwrap_or_else(|| panic!("ordinal {ordinal} of SG {sg} is missing"))
}

fn table_for(subgroup: &IsotropySubgroup, label: &str) -> &'static LittleCharacterTable {
    W_LITTLE_CHARACTERS
        .iter()
        .find(|table| {
            usize::from(table.space_group) == usize::from(subgroup.parent_sg)
                && table.label == label
        })
        .unwrap_or_else(|| panic!("no frozen table for SG {} {label}", subgroup.parent_sg))
}

#[test]
fn an_embedding_of_another_record_is_rejected() {
    // The same source table and the same parent, but the embedding of a
    // different isotropy record: the answer would be computed from another
    // record's geometry, so it must be refused rather than returned.
    let subgroup = find(196, 10032);
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&find(196, 10030)).unwrap();
    let table = table_for(&subgroup, "DT1");
    assert!(line_trivial_content_via_blocks(&subgroup, &embedding, table).is_err());
}

#[test]
fn a_singular_record_basis_is_rejected() {
    let mut broken = find(196, 10030);
    broken.record.basis = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&find(196, 10030)).unwrap();
    let table = table_for(&broken, "DT1");
    assert!(line_trivial_content_via_blocks(&broken, &embedding, table).is_err());
}

#[test]
fn the_block_route_is_the_public_route() {
    // The legacy helper now delegates to the block route, so comparing the two
    // would be `f(a) == f(a)`.  Compare against the **R6.1 public route**
    // instead: the complete line decomposition at the official parameter, whose
    // trivial content must equal the legacy frequency on every pinned row.
    for ordinal in [10030usize, 10032, 10033] {
        let subgroup = find(196, ordinal);
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
        for entry in subgroup.other_wave_vector_subduction().unwrap() {
            let table = table_for(&subgroup, entry.parent_ml);
            let parameter = official_line_parameter().unwrap();
            let complete = subduce_line_at_parameter(&subgroup, &embedding, table, parameter)
                .unwrap_or_else(|error| panic!("ordinal {ordinal} {}: {error}", entry.parent_ml));
            let trivial = complete
                .trivial_content()
                .unwrap_or_else(|error| panic!("ordinal {ordinal} {}: {error}", entry.parent_ml));
            assert_eq!(
                trivial,
                u32::from(entry.frequency),
                "ordinal {ordinal} {} pinned frequency",
                entry.parent_ml
            );
            assert_eq!(
                line_trivial_content_with_embedding(&subgroup, &embedding, table),
                Ok(trivial),
                "ordinal {ordinal} {} legacy route",
                entry.parent_ml
            );
        }
    }
}

#[test]
fn the_block_route_matches_every_pinned_row_of_sg196() {
    // The retired hand-written sum returned a wrong value for 30 of SG 196's w
    // rows, so the audit used to reach the pinned value through a second route
    // and never noticed.  The block route has to get all of them right on its
    // own; this is the assertion that makes an expected-answer-selected audit
    // unnecessary.
    let mut checked = 0usize;
    for subgroup in subgroups_of(196) {
        let Ok(rows) = subgroup.other_wave_vector_subduction() else {
            continue;
        };
        if rows.is_empty() {
            continue;
        }
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
        for entry in rows {
            let table = table_for(&subgroup, entry.parent_ml);
            assert_eq!(
                line_trivial_content_via_blocks(&subgroup, &embedding, table),
                Ok(u32::from(entry.frequency)),
                "ordinal {} {}",
                subgroup.ordinal,
                entry.parent_ml
            );
            checked += 1;
        }
    }
    assert_eq!(checked, 106, "SG 196 has 106 w rows");
}

/// A parameter on the **conjugate** coset (`t = 3/4`) differs from the anchor by
/// a unit-modulus factor built from the frozen direction.
///
/// The frozen direction `v` is a parent reciprocal lattice vector and
/// `3/4 = -1/4 + 1`, so the band the row names at `t = 3/4` is the conjugate of
/// the anchor band, transported by one reciprocal step.  Writing the engine's
/// character convention out (`D` is the frozen little co-group matrix, real for
/// these two sources, and the parameter enters through the Bloch factor) gives
///
/// ```text
/// chi(3/4)(h) = conj(chi(1/4)(h)) * exp(2 pi i v . T_h),
/// ```
///
/// operation by operation, with a **unit-modulus** factor: a reciprocal-shift
/// twist, not a garbled character.  The factor is `-1` exactly on the operations
/// whose translation has a quarter along the direction, so `Phi_v` is a
/// **non-trivial** character of the line's little group -- which is the same
/// statement as `M_v != 1` on these SG 210/219 lines.  The trivial content at
/// `3/4` is therefore the content of the twisted band and is **not** required to
/// equal the pinned value.
///
/// History (three stages, deliberately kept):
///
/// * R6.1 measured this relation against a *canonicalized* wave vector, pinned
///   `12306 SM1: 1 -> 0`, and shipped `canonical_wave_vector` so that a parameter
///   step changed nothing;
/// * R6.2 found that the frozen characters are Gamma-point matrices, so the raw
///   `t v` is what the table describes, and revoked the canonicalization: the
///   step `t -> t + 1` is the monodromy image of the label
///   (`tests/line_monodromy.rs`, `docs/subduction-conventions.md` §16), and the
///   two numbers on this conjugate coset became `1 -> 1` and `1 -> 1`;
/// * the anchor `t = 1/4` stays the only parameter with an external oracle: this
///   test pins the *relation* between the cosets, not a second oracle.
///
/// The previous revision of this test closed with "landing such a change requires
/// deliberately rewriting this test and the contract in
/// `docs/subduction-conventions.md` §16"; that rewrite is this revision.
///
#[test]
fn a_conjugate_parameter_carries_the_reciprocal_gauge_factor() {
    // (parent SG, ordinal, source, pinned content, content at t = 3/4,
    //  monodromy image, negative factors among the enumerated representatives)
    //
    // Both regimes are pinned: SG 210 `DT1` is transported to `DT2` and sees the
    // twist as `-1` on half of the coset, while SG 219 `SM1` is its own image and
    // the twist is `+1` on every enumerated representative.  Pinning the counts
    // keeps the measurement from being silently replaced by an assumption in
    // either direction.
    let cases = [
        (210u8, 11329usize, "DT1", 1u32, 1u32, "DT2", 4usize),
        (219, 12306, "SM1", 1, 1, "SM1", 0),
    ];
    let mut total_negatives = 0usize;
    for (sg, ordinal, label, pinned, at_three_quarters, image, expected_negatives) in cases {
        let subgroup = find(sg, ordinal);
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
        let table = table_for(&subgroup, label);
        let quarter = subduce_line_at_parameter(
            &subgroup,
            &embedding,
            table,
            official_line_parameter().unwrap(),
        )
        .unwrap_or_else(|error| panic!("ordinal {ordinal} {label} t=1/4: {error}"));
        let three = subduce_line_at_parameter(
            &subgroup,
            &embedding,
            table,
            Rat::new(3, 4).unwrap(),
        )
        .unwrap_or_else(|error| panic!("ordinal {ordinal} {label} t=3/4: {error}"));

        assert_eq!(
            quarter.trivial_content().unwrap(),
            pinned,
            "ordinal {ordinal} {label}: the anchor content is the pinned frequency"
        );
        assert_eq!(
            three.trivial_content().unwrap(),
            at_three_quarters,
            "ordinal {ordinal} {label}: the conjugate parameter is the twisted multiplicity, \
             not the pinned one (making them equal changes what `t` means)"
        );

        // G = k_c(3/4) + k_pin, in the parent's own reciprocal coordinates.
        let g = sum(three.wave_vector(), quarter.wave_vector());
        let (quarter_characters, _) = quarter.reconstruction();
        let (three_characters, _) = three.reconstruction();
        let representatives = quarter.representatives();
        assert_eq!(representatives.len(), quarter_characters.len());
        assert_eq!(representatives.len(), three_characters.len());

        let mut negative_factors = 0usize;
        let mut non_unit = 0usize;
        for (index, h) in representatives.iter().enumerate() {
            let factor = gauge_phase(&g, h.translation());
            if (factor.0 * factor.0 + factor.1 * factor.1 - 1.0).abs() > 1e-12 {
                non_unit += 1;
            }
            if factor.0 < -0.5 {
                negative_factors += 1;
            }
            // conj(chi(1/4)) * factor
            let expected_re = quarter_characters[index].re * factor.0
                + quarter_characters[index].im * factor.1;
            let expected_im = quarter_characters[index].re * factor.1
                - quarter_characters[index].im * factor.0;
            let deviation = ((three_characters[index].re - expected_re).powi(2)
                + (three_characters[index].im - expected_im).powi(2))
            .sqrt();
            assert!(
                deviation <= 1e-9,
                "ordinal {ordinal} {label} h[{index}]: chi(3/4) = {:.6}{:+.6}i deviates from \
                 conj(chi(1/4)) * exp(2 pi i G.T) by {deviation:.3e}",
                three_characters[index].re,
                three_characters[index].im
            );
        }
        assert_eq!(non_unit, 0, "the gauge factor must be unit modulus");
        assert_eq!(
            negative_factors, expected_negatives,
            "ordinal {ordinal} {label}: the twisted operations are a measured property of the \
             frozen table, not a convention"
        );
        let map = monodromy(sg, &line_direction(table).expect("frozen direction"))
            .expect("frozen direction is a reciprocal-lattice vector");
        assert_eq!(
            map.unique_image(label),
            Some(image),
            "ordinal {ordinal} {label}: the label the transport names is the monodromy image"
        );
        total_negatives += negative_factors;
    }
    assert!(
        total_negatives > 0,
        "at least one witness must show the twist as a *non-trivial* character: otherwise \
         `t` and `t + 1` would name the same irrep and this test (and the contract) would have \
         to be rewritten again"
    );
}

/// Component-wise sum of two exact vectors.
fn sum(left: &Vec3R, right: &Vec3R) -> Vec3R {
    let mut values = [Rat::ZERO; 3];
    for (axis, value) in values.iter_mut().enumerate() {
        *value = left
            .get(axis)
            .checked_add(right.get(axis))
            .expect("exact sum");
    }
    Vec3R::new(values)
}

/// `exp(2 pi i g . t)` as a real/imaginary pair.
fn gauge_phase(g: &Vec3R, t: &Vec3R) -> (f64, f64) {
    let mut angle = 0.0f64;
    for axis in 0..3 {
        angle += g.get(axis).to_f64() * t.get(axis).to_f64();
    }
    let angle = std::f64::consts::TAU * angle;
    (angle.cos(), angle.sin())
}
