//! Probe: can a candidate-setting search recover the parent-frame convention?
//!
//! Candidates: `A_S = S * A` with `A = (W * C)^T` (W = record basis in the
//! parent's primitive frame, C = the magnetic group's own centring matrix read
//! off its pure translations) and `S` running over the 48 signed permutations.
//! A candidate is accepted when every class maps to an integral rotation, every
//! unitary class maps to a parent operation (mod the parent lattice) and every
//! antiunitary class maps to a parent rotation.

use cryspglib::irrep::magnetic_embedding::{geometry_of, magnetic_operations, translation_lattice};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::{Mat3I, Mat3R, Rat, Vec3R};

fn exact(value: f64) -> Option<Rat> {
    for denominator in 1..=96i128 {
        let scaled = value * denominator as f64;
        if (scaled - scaled.round()).abs() < 1e-9 {
            return Rat::new(scaled.round() as i128, denominator).ok();
        }
    }
    None
}

fn signed_permutations() -> Vec<Mat3I> {
    let mut out = Vec::new();
    let permutations = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    for permutation in permutations {
        for signs in 0..8 {
            let mut matrix = [[0i32; 3]; 3];
            for (row, column) in permutation.iter().enumerate() {
                let sign = if signs & (1 << row) == 0 { 1 } else { -1 };
                matrix[row][*column] = sign;
            }
            out.push(matrix);
        }
    }
    out
}

fn centring_matrix(set: &cryspglib::irrep::magnetic_embedding::MagneticOperationSet) -> Mat3R {
    // Read the magnetic group's centring off its own pure unitary translations.
    let identity = Mat3R::from_ints([[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    let lattice = translation_lattice(set);
    let has = |target: [Rat; 3]| {
        lattice.iter().any(|vector| {
            (0..3).all(|axis| {
                vector[axis].numerator() * target[axis].denominator()
                    == target[axis].numerator() * vector[axis].denominator()
            })
        })
    };
    let half = Rat::new(1, 2).expect("1/2");
    let zero = Rat::ZERO;
    let _ = has([zero, zero, zero]);
    // F: three face vectors; I: body vector; A/B/C: one face vector; R: two thirds.
    if has([half, half, zero]) && has([half, zero, half]) && has([zero, half, half]) {
        return Mat3R::from_ints([[0, 1, 1], [1, 0, 1], [1, 1, 0]]);
    }
    if has([half, half, half]) {
        return Mat3R::from_ints([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]);
    }
    if has([zero, half, half]) {
        return Mat3R::from_ints([[1, 0, 0], [0, 1, 1], [0, -1, 1]]);
    }
    if has([half, zero, half]) {
        return Mat3R::from_ints([[1, 0, 1], [0, 1, 0], [-1, 0, 1]]);
    }
    if has([half, half, zero]) {
        return Mat3R::from_ints([[1, 1, 0], [-1, 1, 0], [0, 0, 1]]);
    }
    let third = Rat::new(1, 3).expect("1/3");
    let two_thirds = Rat::new(2, 3).expect("2/3");
    if has([two_thirds, third, third]) {
        let negate = |value: Rat| value.checked_neg().expect("negate");
        return Mat3R::new([
            [two_thirds, negate(third), negate(third)],
            [third, two_thirds, negate(third)],
            [third, third, two_thirds],
        ]);
    }
    identity
}

fn main() {
    let permutations = signed_permutations();
    let mut reports = Vec::new();
    for sg in 1..=230u8 {
        let parent = query::symmetry_operations_of(sg).expect("parent");
        let parent_operations: Vec<(Mat3I, [Rat; 3])> = parent
            .operations
            .iter()
            .filter_map(|operation| {
                Some((
                    operation.rotation,
                    [
                        exact(operation.translation[0])?,
                        exact(operation.translation[1])?,
                        exact(operation.translation[2])?,
                    ],
                ))
            })
            .collect();
        let parent_rotations: Vec<Mat3I> = parent_operations
            .iter()
            .map(|(rotation, _)| *rotation)
            .collect();
        for row in query::magnetic_isotropy_subgroups_of(sg) {
            let record = row.subgroup;
            let geometry = geometry_of(sg, &record).expect("geometry");
            let set = magnetic_operations(record.mag_sg).expect("operations");
            let centring = centring_matrix(&set);
            // W_conv = W * C, then A = (W_conv)^T.
            let w = Mat3R::new(geometry.basis);
            let w_conv = w.checked_mul(&centring).expect("product");
            let mut a = [[Rat::ZERO; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    a[i][j] = w_conv.row(j)[i];
                }
            }
            let base = Mat3R::new(a);
            let origin = geometry.origin;
            let mut survivors = 0usize;
            for permutation in &permutations {
                let s = Mat3R::from_ints(*permutation);
                let candidate = match s.checked_mul(&base) {
                    Ok(matrix) => matrix,
                    Err(_) => continue,
                };
                let Ok(inverse) = candidate.inverse() else {
                    continue;
                };
                let mut ok = true;
                for operation in &set.classes {
                    let rotation = Mat3R::from_ints(operation.rotation);
                    let conjugated = candidate
                        .checked_mul(&rotation)
                        .and_then(|product| product.checked_mul(&inverse))
                        .expect("product");
                    let Ok(mapped_rotation) = conjugated.to_int_matrix() else {
                        ok = false;
                        break;
                    };
                    if operation.time_reversal {
                        if !parent_rotations.contains(&mapped_rotation) {
                            ok = false;
                            break;
                        }
                        continue;
                    }
                    let mut shift = [Rat::ZERO; 3];
                    for (index, cell) in shift.iter_mut().enumerate() {
                        let mut sum = Rat::ZERO;
                        for (column, value) in origin.iter().enumerate() {
                            let entry = i128::from(if index == column { 1 } else { 0 })
                                - i128::from(mapped_rotation[index][column]);
                            sum = sum
                                .checked_add(
                                    value
                                        .checked_mul(Rat::from_integer(entry))
                                        .expect("product"),
                                )
                                .expect("sum");
                        }
                        *cell = sum;
                    }
                    let moved = candidate
                        .checked_mul_vector(&Vec3R::new(operation.translation))
                        .expect("moved");
                    let mut mapped = [Rat::ZERO; 3];
                    for axis in 0..3 {
                        mapped[axis] =
                            moved.get(axis).checked_add(shift[axis]).expect("sum");
                    }
                    let found = parent_operations.iter().any(|(protor, protau)| {
                        *protor == mapped_rotation
                            && (0..3).all(|axis| {
                                let difference = mapped[axis]
                                    .checked_sub(protau[axis])
                                    .expect("difference");
                                difference.is_integer()
                            })
                    });
                    if !found {
                        ok = false;
                        break;
                    }
                }
                if ok {
                    survivors += 1;
                }
            }
            reports.push((sg, record.mag_sg, survivors));
        }
    }
    let mut histogram: std::collections::BTreeMap<usize, usize> =
        std::collections::BTreeMap::new();
    for (_, _, survivors) in &reports {
        *histogram.entry(*survivors).or_insert(0) += 1;
    }
    println!("records {}; survivor histogram {:?}", reports.len(), histogram);
    let zero: Vec<_> = reports.iter().filter(|(_, _, s)| *s == 0).take(8).collect();
    let many: Vec<_> = reports.iter().filter(|(_, _, s)| *s > 1).take(8).collect();
    println!("no survivor (first 8): {zero:?}");
    println!("several survivors (first 8): {many:?}");
}
