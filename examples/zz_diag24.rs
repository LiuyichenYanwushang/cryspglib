use cryspglib::irrep::magnetic_embedding::{
    basis_in_parent_conventional, geometry_of, magnetic_operations,
};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::{Mat3I, Mat3R, Rat};

fn main() {
    let (sg, uni) = (3u8, 24usize);
    let record = query::magnetic_isotropy_subgroups_of(sg)
        .into_iter()
        .find(|row| row.subgroup.mag_sg == uni)
        .map(|row| row.subgroup)
        .expect("record");
    let geometry = geometry_of(sg, &record).expect("geometry");
    let set = magnetic_operations(uni).expect("ops");
    println!("record basis {:?} origin {:?}", record.basis, record.origin);
    println!("conv basis {:?}", basis_in_parent_conventional(&geometry).expect("conv"));
    let parent = query::symmetry_operations_of(sg).expect("parent");
    println!("parent ops {}", parent.operations.len());
    for operation in &parent.operations {
        println!("  rot {:?} trans {:?}", operation.rotation, operation.translation);
    }
    println!("magnetic classes:");
    for operation in &set.classes {
        println!(
            "  rot {:?} trans ({}/{},{}/{},{}/{}) timerev {}",
            operation.rotation,
            operation.translation[0].numerator(), operation.translation[0].denominator(),
            operation.translation[1].numerator(), operation.translation[1].denominator(),
            operation.translation[2].numerator(), operation.translation[2].denominator(),
            operation.time_reversal
        );
    }
    // Which signed permutations conjugate the magnetic rotations into the parent?
    let parent_rotations: Vec<Mat3I> = parent.operations.iter().map(|op| op.rotation).collect();
    let permutations: Vec<Mat3I> = {
        let perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
        let mut out = Vec::new();
        for p in perms {
            for signs in 0..8 {
                let mut m = [[0i32; 3]; 3];
                for (row, column) in p.iter().enumerate() {
                    m[row][*column] = if signs & (1 << row) == 0 { 1 } else { -1 };
                }
                out.push(m);
            }
        }
        out
    };
    let mut rotation_ok = 0usize;
    for permutation in &permutations {
        let map = Mat3R::from_ints(*permutation);
        let Ok(inverse) = map.inverse() else { continue };
        let mut ok = true;
        for operation in &set.classes {
            let rotation = Mat3R::from_ints(operation.rotation);
            let Ok(conjugated) = map.checked_mul(&rotation).and_then(|p| p.checked_mul(&inverse)) else {
                ok = false;
                break;
            };
            let Ok(mapped) = conjugated.to_int_matrix() else {
                ok = false;
                break;
            };
            if !parent_rotations.contains(&mapped) {
                ok = false;
                break;
            }
        }
        if ok {
            rotation_ok += 1;
        }
    }
    println!("plain signed permutations passing the rotation test: {rotation_ok}/{}", permutations.len());
    // The magnetic group's own pure translations and their images.
    for generator in cryspglib::irrep::magnetic_embedding::translation_lattice(&set) {
        println!("magnetic pure translation ({}/{},{}/{},{}/{})",
            generator[0].numerator(), generator[0].denominator(),
            generator[1].numerator(), generator[1].denominator(),
            generator[2].numerator(), generator[2].denominator());
    }
    let identity = Mat3R::identity();
    for operation in &set.classes {
        let rotation = Mat3R::from_ints(operation.rotation);
        let _ = rotation;
        let moved = identity
            .checked_mul_vector(&cryspglib::irrep::subduction::Vec3R::new(operation.translation))
            .expect("moved");
        println!("identity image of translation = ({}/{},{}/{},{}/{})",
            moved.get(0).numerator(), moved.get(0).denominator(),
            moved.get(1).numerator(), moved.get(1).denominator(),
            moved.get(2).numerator(), moved.get(2).denominator());
    }
    let _ = Rat::ZERO;
}
