//! 磁性对称性支持。
//!
//! 提供磁性张量、共线/非共线磁结构的点群分析等功能。

use crate::cell::{Cell, TensorRank, cel_is_overlap_with_same_type};
use crate::debug;
use crate::mathfunc::{
    Mat3, Mat3I, Vec3, mat_check_identity_matrix_i3, mat_dabs,
    mat_get_determinant_d3, mat_inverse_matrix_d3, mat_multiply_matrix_d3, mat_multiply_matrix_id3,
    mat_multiply_matrix_vector_id3, mat_nint,
};
use crate::primitive::prm_get_primitive_lattice_vectors;
use crate::symmetry::{MagneticSymmetry, Symmetry};

static IDENTITY: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// 获取磁性对称操作
/// equivalent_atoms: 输出参数，等价原子映射表
/// permutations: 输出参数，对称操作的原子排列
/// prim_lattice: 输出参数，原胞晶格
pub fn spn_get_operations_with_site_tensors(
    equivalent_atoms: &mut Vec<i32>,
    permutations: &mut Vec<i32>,
    prim_lattice: &mut Mat3,
    sym_nonspin: &Symmetry,
    cell: &Cell,
    with_time_reversal: bool,
    is_axial: bool,
    symprec: f64,
    angle_tolerance: f64,
    mag_symprec_: f64,
) -> Option<MagneticSymmetry> {
    // TODO: More robust way to guess mag_symprec
    let mag_symprec = if mag_symprec_ < 0.0 {
        symprec
    } else {
        mag_symprec_
    };

    let magnetic_symmetry = get_operations(
        sym_nonspin,
        cell,
        with_time_reversal,
        is_axial,
        symprec,
        mag_symprec,
    )?;

    /* equivalent atoms */
    *permutations = get_symmetry_permutations(
        &magnetic_symmetry,
        cell,
        with_time_reversal,
        is_axial,
        symprec,
        mag_symprec,
    )?;

    *equivalent_atoms = get_orbits(permutations, magnetic_symmetry.size, cell.size)?;

    let pure_trans = spn_collect_pure_translations_from_magnetic_symmetry(&magnetic_symmetry)?;

    let Some((found_lattice, multi)) = prm_get_primitive_lattice_vectors(
        cell,
        &pure_trans,
        symprec,
        angle_tolerance,
    ) else {
        return None;
    };

    *prim_lattice = found_lattice;

    /* By definition, change of number of pure translations would */
    /* not be allowed. */
    if multi != pure_trans.len() {
        return None;
    }

    // pure_trans (VecDBL) is dropped automatically
    Some(magnetic_symmetry)
}

pub fn spn_collect_pure_translations_from_magnetic_symmetry(
    sym_msg: &MagneticSymmetry,
) -> Option<Vec<Vec3>> {
    let mut pure_trans = Vec::new();

    for i in 0..sym_msg.size {
        /* Take translation with rot=identity and timerev=false */
        /* time reversal should be considered for type-IV magnetic space group */
        if mat_check_identity_matrix_i3(&IDENTITY, &sym_msg.rot[i]) && !sym_msg.timerev[i] {
            pure_trans.push(sym_msg.trans[i]);
        }
    }

    if pure_trans.is_empty() {
        // In C code, it returns NULL if alloc fails, but here empty vec is valid but logic might expect something?
        // C code logic: if alloc fails return NULL. If num_pure_trans is 0, it returns a VecDBL with size 0.
        // But usually there is at least identity.
        // We return Some(vec) even if empty, or None if something critical failed (not applicable here).
        // To match C's "return NULL on failure" vs "return object", we return Some.
    }

    Some(pure_trans)
}

/// Apply special position operator to `cell`.
pub fn spn_get_idealized_cell(
    permutations: &[i32],
    cell: &Cell,
    magnetic_symmetry: &MagneticSymmetry,
    with_time_reversal: bool,
    is_axial: bool,
) -> Option<Cell> {
    let mut exact_cell = Cell::new(cell.size, cell.tensor_rank);
    exact_cell.lattice = cell.lattice;
    exact_cell.aperiodic_axis = cell.aperiodic_axis;
    exact_cell.types = cell.types.clone();

    let mut inv_perm = vec![0; cell.size];
    let mut rotations_cart = vec![[[0.0; 3]; 3]; magnetic_symmetry.size];
    set_rotations_in_cartesian(&mut rotations_cart, &cell.lattice, magnetic_symmetry);

    for i in 0..cell.size {
        let mut pos_res = [0.0; 3];
        let mut scalar_res = 0.0;
        let mut vector_res = [0.0; 3];

        for p in 0..magnetic_symmetry.size {
            // Inverse of the p-th permutation
            for j in 0..cell.size {
                inv_perm[permutations[p * cell.size + j] as usize] = j;
            }

            let j = inv_perm[i]; // p-th operation maps site-j to site-i

            let mut pos_tmp = [0.0; 3];
            pos_tmp = apply_symmetry_to_position(
                &cell.position[j],
                &magnetic_symmetry.rot[p],
                &magnetic_symmetry.trans[p],
            );

            for s in 0..3 {
                // To minimize rounding error, subtract by the original position[i][s].
                // `pos_tmp[s] - cell->position[i][s]` should be close to integer.
                let diff = pos_tmp[s] - cell.position[i][s];
                pos_res[s] += diff - mat_nint(diff) as f64;
            }

            if cell.tensor_rank == TensorRank::Collinear {
                let scalar_tmp = apply_symmetry_to_site_scalar(
                    cell.tensors[j],
                    &rotations_cart[p],
                    magnetic_symmetry.timerev[p],
                    with_time_reversal,
                    is_axial,
                );
                scalar_res += scalar_tmp - cell.tensors[i];
            } else if cell.tensor_rank == TensorRank::NonCollinear {
                let mut vector_tmp = [0.0; 3];
                apply_symmetry_to_site_vector(
                    &mut vector_tmp,
                    j,
                    &cell.tensors,
                    &rotations_cart[p],
                    magnetic_symmetry.timerev[p],
                    with_time_reversal,
                    is_axial,
                );
                for s in 0..3 {
                    vector_res[s] += vector_tmp[s] - cell.tensors[3 * i + s];
                }
            }
        }

        for s in 0..3 {
            exact_cell.position[i][s] =
                cell.position[i][s] + pos_res[s] / magnetic_symmetry.size as f64;
        }

        debug::debug_print(format_args!("Idealize position\n"));
        debug::debug_print_vector_d3(&cell.position[i]);
        debug::debug_print_vector_d3(&exact_cell.position[i]);

        debug::debug_print(format_args!("Idealize site tensor\n"));
        if cell.tensor_rank == TensorRank::Collinear {
            exact_cell.tensors[i] = cell.tensors[i] + scalar_res / magnetic_symmetry.size as f64;
            debug::debug_print(format_args!("{}\n", cell.tensors[i]));
            debug::debug_print(format_args!("{}\n", exact_cell.tensors[i]));
        } else if cell.tensor_rank == TensorRank::NonCollinear {
            for s in 0..3 {
                exact_cell.tensors[3 * i + s] =
                    cell.tensors[3 * i + s] + vector_res[s] / magnetic_symmetry.size as f64;
            }
        }
    }

    Some(exact_cell)
}

// In Rust, Vec handles allocation, so spn_alloc_site_tensors is not needed.

// --- Local functions ---

fn get_operations(
    sym_nonspin: &Symmetry,
    cell: &Cell,
    with_time_reversal: bool,
    is_axial: bool,
    symprec: f64,
    mag_symprec: f64,
) -> Option<MagneticSymmetry> {
    let mut rotations_cart = vec![[[0.0; 3]; 3]; sym_nonspin.size];

    // FIX: mat_inverse_matrix_d3 returns Option<Mat3>
    let inv_lat = match mat_inverse_matrix_d3(&cell.lattice, 0.0) {
        Some(m) => m,
        None => return None,
    };

    for i in 0..sym_nonspin.size {
        // rot_cart = lattice @ rot @ lattice^-1
        let temp = mat_multiply_matrix_id3(&sym_nonspin.rot[i], &inv_lat);
        rotations_cart[i] = mat_multiply_matrix_d3(&cell.lattice, &temp);
    }

    let mut rotations = Vec::new();
    let mut trans = Vec::new();
    let mut spin_flips = Vec::new();

    for i in 0..sym_nonspin.size {
        let mut found = true;
        let mut determined = false;
        let mut sign = 0;

        for j in 0..cell.size {
            let mut pos = [0.0; 3];
            pos = apply_symmetry_to_position(
                &cell.position[j],
                &sym_nonspin.rot[i],
                &sym_nonspin.trans[i],
            );

            let mut k = cell.size; // Default to not found
            for (idx, pos_k) in cell.position.iter().enumerate() {
                if cel_is_overlap_with_same_type(
                    pos_k,
                    &pos,
                    cell.types[idx],
                    cell.types[j],
                    &cell.lattice,
                    symprec,
                ) {
                    k = idx;
                    break;
                }
            }

            if k == cell.size {
                debug::debug_print(format_args!(
                    "Failed to overlap atom-{} by operation-{}\n",
                    j, i
                ));
                found = false;
                break;
            }

            // Skip if relevant tensors are zeros
            if cell.tensor_rank == TensorRank::Collinear {
                if is_zero(cell.tensors[j], 0.5 * mag_symprec)
                    && is_zero(cell.tensors[k], 0.5 * mag_symprec)
                {
                    continue;
                }
            }
            if cell.tensor_rank == TensorRank::NonCollinear {
                if is_zero(cell.tensors[j * 3], 0.5 * mag_symprec)
                    && is_zero(cell.tensors[j * 3 + 1], 0.5 * mag_symprec)
                    && is_zero(cell.tensors[j * 3 + 2], 0.5 * mag_symprec)
                    && is_zero(cell.tensors[k * 3], 0.5 * mag_symprec)
                    && is_zero(cell.tensors[k * 3 + 1], 0.5 * mag_symprec)
                    && is_zero(cell.tensors[k * 3 + 2], 0.5 * mag_symprec)
                {
                    continue;
                }
            }

            if !determined {
                if cell.tensor_rank == TensorRank::Collinear {
                    sign = get_operation_sign_on_scalar(
                        cell.tensors[j],
                        cell.tensors[k],
                        &rotations_cart[i],
                        with_time_reversal,
                        is_axial,
                        mag_symprec,
                    );
                }
                if cell.tensor_rank == TensorRank::NonCollinear {
                    sign = get_operation_sign_on_vector(
                        j,
                        k,
                        &cell.tensors,
                        &rotations_cart[i],
                        with_time_reversal,
                        is_axial,
                        mag_symprec,
                    );
                }
                determined = true;

                if sign == 0 || (!with_time_reversal && sign != 1) {
                    found = false;
                    break;
                }
            } else {
                // Check consistency
                if cell.tensor_rank == TensorRank::Collinear {
                    if get_operation_sign_on_scalar(
                        cell.tensors[j],
                        cell.tensors[k],
                        &rotations_cart[i],
                        with_time_reversal,
                        is_axial,
                        mag_symprec,
                    ) != sign
                    {
                        found = false;
                        break;
                    }
                }
                if cell.tensor_rank == TensorRank::NonCollinear {
                    if get_operation_sign_on_vector(
                        j,
                        k,
                        &cell.tensors,
                        &rotations_cart[i],
                        with_time_reversal,
                        is_axial,
                        mag_symprec,
                    ) != sign
                    {
                        found = false;
                        break;
                    }
                }
            }
        }

        if found {
            if determined {
                rotations.push(sym_nonspin.rot[i]);
                trans.push(sym_nonspin.trans[i]);
                if with_time_reversal {
                    spin_flips.push(sign);
                }
            } else if with_time_reversal {
                // sign=1
                rotations.push(sym_nonspin.rot[i]);
                trans.push(sym_nonspin.trans[i]);
                spin_flips.push(1);

                // sign=-1
                rotations.push(sym_nonspin.rot[i]);
                trans.push(sym_nonspin.trans[i]);
                spin_flips.push(-1);
            } else {
                rotations.push(sym_nonspin.rot[i]);
                trans.push(sym_nonspin.trans[i]);
            }
        }
    }

    let num_sym = rotations.len();
    let mut magnetic_symmetry = MagneticSymmetry::new(num_sym);
    for i in 0..num_sym {
        magnetic_symmetry.rot[i] = rotations[i];
        magnetic_symmetry.trans[i] = trans[i];
        if with_time_reversal {
            magnetic_symmetry.timerev[i] = spin_flips[i] == -1;
        } else {
            magnetic_symmetry.timerev[i] = false;
        }
        debug::debug_print(format_args!("-- {} -- \n", i));
        debug::debug_print_matrix_i3(&magnetic_symmetry.rot[i]);
        debug::debug_print_vector_d3(&magnetic_symmetry.trans[i]);
        debug::debug_print(format_args!("timerev {}\n", magnetic_symmetry.timerev[i]));
    }

    Some(magnetic_symmetry)
}

fn get_symmetry_permutations(
    magnetic_symmetry: &MagneticSymmetry,
    cell: &Cell,
    with_time_reversal: bool,
    is_axial: bool,
    symprec: f64,
    mag_symprec: f64,
) -> Option<Vec<i32>> {
    let mut permutations = vec![-1; magnetic_symmetry.size * cell.size];
    let mut rotations_cart = vec![[[0.0; 3]; 3]; magnetic_symmetry.size];

    set_rotations_in_cartesian(&mut rotations_cart, &cell.lattice, magnetic_symmetry);

    for p in 0..magnetic_symmetry.size {
        for i in 0..cell.size {
            let mut pos = [0.0; 3];
            pos = apply_symmetry_to_position(
                &cell.position[i],
                &magnetic_symmetry.rot[p],
                &magnetic_symmetry.trans[p],
            );

            let mut scalar = 0.0;
            let mut vector = [0.0; 3];

            if cell.tensor_rank == TensorRank::Collinear {
                scalar = apply_symmetry_to_site_scalar(
                    cell.tensors[i],
                    &rotations_cart[p],
                    magnetic_symmetry.timerev[p],
                    with_time_reversal,
                    is_axial,
                );
            } else if cell.tensor_rank == TensorRank::NonCollinear {
                apply_symmetry_to_site_vector(
                    &mut vector,
                    i,
                    &cell.tensors,
                    &rotations_cart[p],
                    magnetic_symmetry.timerev[p],
                    with_time_reversal,
                    is_axial,
                );
            }

            for j in 0..cell.size {
                if !cel_is_overlap_with_same_type(
                    &pos,
                    &cell.position[j],
                    cell.types[i],
                    cell.types[j],
                    &cell.lattice,
                    symprec,
                ) {
                    continue;
                }

                debug::debug_print(format_args!(
                    "Try to overlap site-{} ({}) with site-{} ({})\n",
                    i, scalar, j, cell.tensors[j]
                ));

                if cell.tensor_rank == TensorRank::Collinear
                    && !is_zero(cell.tensors[j] - scalar, mag_symprec)
                {
                    continue;
                }
                if cell.tensor_rank == TensorRank::NonCollinear {
                    let diff = [
                        cell.tensors[3 * j] - vector[0],
                        cell.tensors[3 * j + 1] - vector[1],
                        cell.tensors[3 * j + 2] - vector[2],
                    ];
                    if !is_zero_d3(&diff, mag_symprec) {
                        continue;
                    }
                }

                permutations[p * cell.size + i] = j as i32;
                break;
            }

            if permutations[p * cell.size + i] == -1 {
                debug::debug_print(format_args!(
                    "Failed to map site-{} by operation-{}\n",
                    i, p
                ));
                return None; // Unreachable in theory
            }
        }

        debug::debug_print(format_args!("Operation {}\n", p));
        for i in 0..cell.size {
            debug::debug_print(format_args!(" {}", permutations[p * cell.size + i]));
        }
        debug::debug_print(format_args!("\n"));
    }

    Some(permutations)
}

fn get_orbits(permutations: &[i32], num_sym: usize, num_atoms: usize) -> Option<Vec<i32>> {
    let mut equivalent_atoms = vec![-1; num_atoms];

    for i in 0..num_atoms {
        if equivalent_atoms[i] != -1 {
            continue;
        }

        equivalent_atoms[i] = i as i32;
        for s in 0..num_sym {
            let target = permutations[s * num_atoms + i] as usize;
            if target < num_atoms {
                equivalent_atoms[target] = i as i32;
            }
        }
    }

    Some(equivalent_atoms)
}

fn get_operation_sign_on_scalar(
    spin_j: f64,
    spin_k: f64,
    rot_cart: &Mat3,
    with_time_reversal: bool,
    is_axial: bool,
    mag_symprec: f64,
) -> i32 {
    for &timerev in &[false, true] {
        let spin_j_ops =
            apply_symmetry_to_site_scalar(spin_j, rot_cart, timerev, with_time_reversal, is_axial);
        if is_zero(spin_k - spin_j_ops, mag_symprec) {
            return if timerev { -1 } else { 1 }; // Spin-flip
        }
    }
    0
}

fn get_operation_sign_on_vector(
    j: usize,
    k: usize,
    vectors: &[f64],
    rot_cart: &Mat3,
    with_time_reversal: bool,
    is_axial: bool,
    mag_symprec: f64,
) -> i32 {
    for &timerev in &[false, true] {
        let mut vec_j_ops = [0.0; 3];
        apply_symmetry_to_site_vector(
            &mut vec_j_ops,
            j,
            vectors,
            rot_cart,
            timerev,
            with_time_reversal,
            is_axial,
        );
        let mut diff = [0.0; 3];
        for s in 0..3 {
            diff[s] = vectors[3 * k + s] - vec_j_ops[s];
        }
        if is_zero_d3(&diff, mag_symprec) {
            return if timerev { -1 } else { 1 }; // Spin-flip
        }
    }
    0
}

/// 将对称操作作用于分数坐标: `pos' = rot·pos + trans`。
fn apply_symmetry_to_position(pos_src: &Vec3, rot: &Mat3I, trans: &Vec3) -> Vec3 {
    let mut pos_dst = mat_multiply_matrix_vector_id3(rot, pos_src);
    for k in 0..3 {
        pos_dst[k] += trans[k];
    }
    pos_dst
}

fn apply_symmetry_to_site_scalar(
    src: f64,
    rot_cart: &Mat3,
    timerev: bool,
    with_time_reversal: bool,
    is_axial: bool,
) -> f64 {
    let mut dst = if with_time_reversal && timerev {
        -src
    } else {
        src
    };

    if is_axial {
        let det = mat_get_determinant_d3(rot_cart);
        dst *= det;
    }
    dst
}

fn apply_symmetry_to_site_vector(
    dst: &mut Vec3,
    idx: usize,
    tensors: &[f64],
    rot_cart: &Mat3,
    timerev: bool,
    with_time_reversal: bool,
    is_axial: bool,
) {
    let mut vec = [0.0; 3];
    for s in 0..3 {
        vec[s] = tensors[3 * idx + s];
    }

    let det = mat_get_determinant_d3(rot_cart);
    *dst = crate::mathfunc::mat_multiply_matrix_vector_d3(rot_cart, &vec);

    for s in 0..3 {
        if with_time_reversal && timerev {
            dst[s] *= -1.0;
        }
        if is_axial {
            dst[s] *= det;
        }
    }
}

fn set_rotations_in_cartesian(
    rotations_cart: &mut [Mat3],
    lattice: &Mat3,
    magnetic_symmetry: &MagneticSymmetry,
) {
    // FIX: mat_inverse_matrix_d3 returns Option<Mat3>
    let inv_lat = mat_inverse_matrix_d3(lattice, 0.0).unwrap_or([[0.0; 3]; 3]);

    for i in 0..magnetic_symmetry.size {
        // rot_cart = lattice @ rot @ lattice^-1
        let temp = mat_multiply_matrix_id3(&magnetic_symmetry.rot[i], &inv_lat);
        rotations_cart[i] = mat_multiply_matrix_d3(lattice, &temp);
    }
}

fn is_zero(a: f64, mag_symprec: f64) -> bool {
    mat_dabs(a) < mag_symprec
}

fn is_zero_d3(a: &[f64; 3], mag_symprec: f64) -> bool {
    for i in 0..3 {
        if !is_zero(a[i], mag_symprec) {
            return false;
        }
    }
    true
}
