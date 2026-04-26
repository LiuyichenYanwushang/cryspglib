//! Bravais 格子精细化和精确结构确定。
//!
//! 将空间群搜索找到的 Bravais 格子旋转到标准取向，
//! 利用位点对称性数据库确定原子的精确 Wyckoff 位置，
//! 并恢复原始晶胞中的对称操作。

use crate::cell::{cel_is_overlap, cel_is_overlap_with_same_type, cel_layer_is_overlap_with_same_type, AperiodicAxis, Cell};
use crate::debug;
use crate::mathfunc::{
    Mat3, Mat3I, Vec3, mat_cast_matrix_3d_to_3i, mat_cast_matrix_3i_to_3d,
    mat_check_identity_matrix_d3, mat_check_identity_matrix_i3,
    mat_cross_product_d3, mat_dmod1, mat_get_determinant_i3, mat_get_metric,
    mat_inverse_matrix_d3, mat_multiply_matrix_d3, mat_multiply_matrix_di3,
    mat_multiply_matrix_id3, mat_multiply_matrix_vector_d3,
    mat_multiply_matrix_vector_id3, mat_nint, mat_norm_squared_d3,
    mat_transpose_matrix_d3,
};
use crate::pointgroup::{ptg_get_pointgroup, Holohedry};
use crate::site_symmetry::ssm_get_exact_positions;
use crate::spacegroup::Spacegroup;
use crate::spg_database::spgdb_get_spacegroup_operations;
use crate::symmetry::Symmetry;

const IDENTITY: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// 精确结构信息。
///
/// 包含 Bravais 格子中的标准化原子位置、Wyckoff 标记、
/// 位点对称性符号以及旋转到标准取向的旋转矩阵。
pub struct ExactStructure {
    pub bravais: Cell,
    pub symmetry: Symmetry,
    pub wyckoffs: Vec<i32>,
    pub site_symmetry_symbols: Vec<String>,
    pub equivalent_atoms: Vec<i32>,
    pub crystallographic_orbits: Vec<i32>,
    pub std_mapping_to_primitive: Vec<i32>,
    pub rotation: Mat3,
}

/// 精细化结构并获取完整对称信息。
///
/// 重写 `spacegroup.bravais_lattice` 和 `spacegroup.origin_shift`。
/// 返回 `None` 表示失败。
pub fn ref_get_exact_structure_and_symmetry(
    spacegroup: &mut Spacegroup,
    primitive: &Cell,
    cell: &Cell,
    mapping_table: &[i32],
    symprec: f64,
) -> Option<ExactStructure> {
    if !ref_find_similar_bravais_lattice(spacegroup, symprec) {
        return None;
    }

    let symmetry = get_refined_symmetry_operations(cell, primitive, spacegroup, symprec)?;

    let n = cell.size;
    let mut wyckoffs = vec![0i32; n];
    let mut site_symmetry_symbols: Vec<String> = vec![String::new(); n];
    let mut equivalent_atoms = vec![0i32; n];
    let mut crystallographic_orbits = vec![0i32; n];
    let mut std_mapping_to_primitive = vec![0i32; primitive.size * 4];

    let bravais = get_wyckoff_positions(
        &mut wyckoffs,
        &mut site_symmetry_symbols,
        &mut equivalent_atoms,
        &mut crystallographic_orbits,
        &mut std_mapping_to_primitive,
        primitive,
        cell,
        spacegroup,
        &symmetry,
        mapping_table,
        symprec,
    )?;

    let mut rotation = [[0.0; 3]; 3];
    rotation = measure_rigid_rotation(&spacegroup.bravais_lattice, &bravais.lattice);

    Some(ExactStructure {
        bravais,
        symmetry,
        wyckoffs,
        site_symmetry_symbols,
        equivalent_atoms,
        crystallographic_orbits,
        std_mapping_to_primitive,
        rotation,
    })
}

/// 获取 Wyckoff 位置（内部函数）。
fn get_wyckoff_positions(
    wyckoffs: &mut [i32],
    site_symmetry_symbols: &mut [String],
    equiv_atoms: &mut [i32],
    crystallographic_orbits: &mut [i32],
    std_mapping_to_primitive: &mut [i32],
    primitive: &Cell,
    cell: &Cell,
    spacegroup: &Spacegroup,
    symmetry: &Symmetry,
    mapping_table: &[i32],
    symprec: f64,
) -> Option<Cell> {
    debug::debug_print(format_args!("get_Wyckoff_positions\n"));

    let prim_size_4 = primitive.size * 4;
    let mut wyckoffs_bravais = vec![0i32; prim_size_4];
    let mut site_sym_symbols_bravais: Vec<String> = vec![String::new(); prim_size_4];
    let mut equiv_atoms_bravais = vec![0i32; prim_size_4];

    let bravais = get_bravais_exact_positions_and_lattice(
        &mut wyckoffs_bravais,
        &mut site_sym_symbols_bravais,
        &mut equiv_atoms_bravais,
        std_mapping_to_primitive,
        spacegroup,
        primitive,
        symprec,
    )?;

    // Map bravais-level data to cell-level
    for i in 0..cell.size {
        wyckoffs[i] = wyckoffs_bravais[mapping_table[i] as usize];
        site_symmetry_symbols[i] = site_sym_symbols_bravais[mapping_table[i] as usize].clone();
    }

    // Set crystallographic orbits
    if !set_crystallographic_orbits(
        crystallographic_orbits,
        primitive,
        cell,
        &equiv_atoms_bravais,
        mapping_table,
    ) {
        debug::warning_print(format_args!(
            "spglib: set_crystallographic_orbits failed.\n"
        ));
        return None;
    }

    // Check symmetry breaking by unusual multiplicity of primitive cell
    let (op_count, _) = crate::spg_database::spgdb_get_operation_index(spacegroup.hall_number);
    let num_prim_sym = op_count as i32 / ((bravais.size / primitive.size) as i32);

    if cell.size * num_prim_sym as usize != symmetry.size * primitive.size {
        set_equivalent_atoms_broken_symmetry(
            equiv_atoms,
            cell,
            symmetry,
            mapping_table,
            symprec,
        );
    } else {
        for i in 0..cell.size {
            equiv_atoms[i] = crystallographic_orbits[i];
        }
    }

    Some(bravais)
}

/// 获取 Bravais 格子中的精确原子位置和格子。
fn get_bravais_exact_positions_and_lattice(
    wyckoffs: &mut [i32],
    site_symmetry_symbols: &mut [String],
    equiv_atoms: &mut [i32],
    std_mapping_to_primitive: &mut [i32],
    spacegroup: &Spacegroup,
    primitive: &Cell,
    symprec: f64,
) -> Option<Cell> {
    debug::debug_print(format_args!("get_bravais_exact_positions_and_lattice\n"));

    // Get conventional primitive cell
    let mut conv_prim = get_conventional_primitive(spacegroup, primitive)?;

    // Get symmetry operations from database
    let conv_sym = spgdb_get_spacegroup_operations(spacegroup.hall_number)?;
    let num_pure_trans = get_number_of_pure_translation(&conv_sym);

    // Set conventional lattice
    conv_prim.lattice = ref_get_conventional_lattice(spacegroup);

    // Set aperiodic axis
    conv_prim.aperiodic_axis = if spacegroup.hall_number > 0 { None } else { Some(AperiodicAxis::Z) };

    // Get exact positions via site symmetry
    let (exact_positions, wyckoffs_prim, equiv_atoms_prim, site_symmetry_symbols_prim) =
        ssm_get_exact_positions(
            &conv_prim,
            &conv_sym,
            num_pure_trans,
            spacegroup.hall_number,
            symprec,
        )?;

    // Copy exact positions back to conv_prim
    for i in 0..conv_prim.size {
        conv_prim.position[i] = exact_positions[i];
    }

    // Expand positions to Bravais cell
    let bravais = expand_positions_in_bravais(
        wyckoffs,
        site_symmetry_symbols,
        equiv_atoms,
        std_mapping_to_primitive,
        &conv_prim,
        &conv_sym,
        num_pure_trans,
        &wyckoffs_prim,
        &site_symmetry_symbols_prim,
        &equiv_atoms_prim,
    );

    Some(bravais)
}

/// 将原胞中的位置扩展到完整 Bravais 格子。
fn expand_positions_in_bravais(
    wyckoffs: &mut [i32],
    site_symmetry_symbols: &mut [String],
    equiv_atoms: &mut [i32],
    std_mapping_to_primitive: &mut [i32],
    conv_prim: &Cell,
    conv_sym: &Symmetry,
    num_pure_trans: i32,
    wyckoffs_prim: &[i32],
    site_symmetry_symbols_prim: &[String],
    equiv_atoms_prim: &[i32],
) -> Cell {
    let bravais_size = conv_prim.size * num_pure_trans as usize;
    let mut bravais = Cell::new(bravais_size, conv_prim.tensor_rank);

    let mut num_atom = 0usize;
    for i in 0..conv_sym.size {
        if mat_check_identity_matrix_i3(&IDENTITY, &conv_sym.rot[i]) {
            for j in 0..conv_prim.size {
                bravais.types[num_atom] = conv_prim.types[j];
                bravais.position[num_atom] = conv_prim.position[j];
                for k in 0..3 {
                    bravais.position[num_atom][k] += conv_sym.trans[i][k];
                }
                wyckoffs[num_atom] = wyckoffs_prim[j];
                site_symmetry_symbols[num_atom] = site_symmetry_symbols_prim[j].clone();
                equiv_atoms[num_atom] = equiv_atoms_prim[j];
                std_mapping_to_primitive[num_atom] = j as i32;
                num_atom += 1;
            }
        }
    }

    // Modulo 1 for periodic axes
    let lattice_rank = conv_prim.aperiodic_axis.map_or(3, |_| 2);
    for i in 0..num_atom {
        for k in 0..lattice_rank {
            bravais.position[i][k] = mat_dmod1(bravais.position[i][k]);
        }
    }

    bravais.lattice = conv_prim.lattice;
    bravais.aperiodic_axis = conv_prim.aperiodic_axis;

    bravais
}

/// 计算纯平移操作的数量。
fn get_number_of_pure_translation(conv_sym: &Symmetry) -> i32 {
    let mut count = 0;
    for i in 0..conv_sym.size {
        if mat_check_identity_matrix_i3(&IDENTITY, &conv_sym.rot[i]) {
            count += 1;
        }
    }
    count
}

/// 获取常规原胞：将原胞坐标变换到 Bravais 格子坐标系中。
fn get_conventional_primitive(spacegroup: &Spacegroup, primitive: &Cell) -> Option<Cell> {
    let mut conv_prim = Cell::new(primitive.size, primitive.tensor_rank);

    let inv_brv = mat_inverse_matrix_d3(&spacegroup.bravais_lattice, 0.0)?;
    let trans_mat = mat_multiply_matrix_d3(&inv_brv, &primitive.lattice);

    for i in 0..primitive.size {
        conv_prim.types[i] = primitive.types[i];
        conv_prim.position[i] = mat_multiply_matrix_vector_d3(&trans_mat, &primitive.position[i]);
        for j in 0..3 {
            conv_prim.position[i][j] += spacegroup.origin_shift[j];
            if primitive.aperiodic_axis.is_none() || j != 2 {
                conv_prim.position[i][j] = mat_dmod1(conv_prim.position[i][j]);
            }
        }
    }

    Some(conv_prim)
}

/// 计算给定空间群类型的标准化晶格。
///
/// 根据全形对称性 (holohedry) 从 metric tensor 计算出理想化的标准晶格矢量。
/// 根据空间群和 Bravais 晶格计算常规晶格矩阵。
///
/// 通过点群的 holohedry 类型确定晶格约束，
/// 将自由的 Bravais 晶格约化到标准常规形式。
///
/// # Returns
/// 满足空间群约束的标准常规晶格矩阵。
pub fn ref_get_conventional_lattice(spacegroup: &Spacegroup) -> Mat3 {
    let pointgroup = ptg_get_pointgroup(spacegroup.pointgroup_number);

    let mut lattice = [[0.0; 3]; 3];

    let metric = mat_get_metric(&spacegroup.bravais_lattice);

    debug::debug_print(format_args!("bravais lattice\n"));
    debug::debug_print_matrix_d3(&spacegroup.bravais_lattice);
    debug::debug_print(format_args!("{}\n", spacegroup.choice));

    match pointgroup.holohedry {
        Holohedry::Tricli => lattice = set_tricli(&metric),
        Holohedry::Monocli => {
            if spacegroup.hall_number > 0 {
                lattice = set_monocli(&metric, &spacegroup.choice);
            } else {
                lattice = set_layer_monocli(&metric, &spacegroup.choice);
            }
        }
        Holohedry::Ortho => lattice = set_ortho(&metric),
        Holohedry::Tetra => lattice = set_tetra(&metric),
        Holohedry::Trigo => {
            if spacegroup.choice.starts_with('R') {
                lattice = set_rhomb(&metric);
            } else {
                lattice = set_trigo(&metric);
            }
        }
        Holohedry::Hexa => lattice = set_trigo(&metric),
        Holohedry::Cubic => lattice = set_cubic(&metric),
        Holohedry::None => {}
    }
    lattice
}

// --- 标准晶格生成函数 ---

fn set_tricli(metric: &Mat3) -> Mat3 {
    let mut lattice = [[0.0; 3]; 3];
    let a = metric[0][0].sqrt();
    let b = metric[1][1].sqrt();
    let c = metric[2][2].sqrt();
    let alpha = (metric[1][2] / b / c).acos();
    let beta = (metric[0][2] / a / c).acos();
    let gamma = (metric[0][1] / a / b).acos();

    let cg = gamma.cos();
    let cb = beta.cos();
    let ca = alpha.cos();
    let sg = gamma.sin();

    lattice[0][0] = a;
    lattice[0][1] = b * cg;
    lattice[0][2] = c * cb;
    lattice[1][1] = b * sg;
    lattice[1][2] = c * (ca - cb * cg) / sg;
    lattice[2][2] = c * (1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg).sqrt() / sg;
    lattice
}

fn set_monocli(metric: &Mat3, choice: &str) -> Mat3 {
    let mut lattice = [[0.0; 3]; 3];
    debug::debug_print(format_args!("set_monocli:\n"));
    debug::debug_print_matrix_d3(metric);

    let pos_char: usize = if choice.starts_with('-') { 1 } else { 0 };

    let a = metric[0][0].sqrt();
    let b = metric[1][1].sqrt();
    let c = metric[2][2].sqrt();

    match choice.as_bytes().get(pos_char) {
        Some(b'a') => {
            let angle = (metric[1][2] / b / c).acos();
            lattice[0][2] = c;
            lattice[1][0] = a;
            lattice[0][1] = b * angle.cos();
            lattice[2][1] = b * angle.sin();
        }
        Some(b'b') => {
            let angle = (metric[0][2] / a / c).acos();
            lattice[0][0] = a;
            lattice[1][1] = b;
            lattice[0][2] = c * angle.cos();
            lattice[2][2] = c * angle.sin();
        }
        Some(b'c') => {
            let angle = (metric[0][1] / a / b).acos();
            lattice[0][1] = b;
            lattice[1][2] = c;
            lattice[0][0] = a * angle.cos();
            lattice[2][0] = a * angle.sin();
        }
        _ => {
            debug::warning_print(format_args!(
                "spglib: Monoclinic unique axis could not be found.\n"
            ));
        }
    }
    lattice
}

fn set_layer_monocli(metric: &Mat3, choice: &str) -> Mat3 {
    let mut lattice = [[0.0; 3]; 3];
    debug::debug_print(format_args!("set_layer_monocli:\n"));
    debug::debug_print_matrix_d3(metric);

    let a = metric[0][0].sqrt();
    let b = metric[1][1].sqrt();
    let c = metric[2][2].sqrt();

    match choice.as_bytes().first() {
        Some(b'a') => {
            let angle = (metric[1][2] / b / c).acos();
            lattice[0][0] = a;
            lattice[1][1] = b;
            lattice[1][2] = c * angle.cos();
            lattice[2][2] = c * angle.sin();
        }
        Some(b'b') => {
            let angle = (metric[0][2] / a / c).acos();
            lattice[0][0] = b;
            lattice[1][1] = a;
            lattice[0][2] = c * angle.cos();
            lattice[2][2] = c * angle.sin();
        }
        Some(b'c') => {
            let angle = (metric[0][1] / a / b).acos();
            lattice[0][0] = a;
            lattice[0][1] = b * angle.cos();
            lattice[1][1] = b * angle.sin();
            lattice[2][2] = c;
        }
        _ => {
            debug::warning_print(format_args!(
                "spglib: Monoclinic unique axis could not be found.\n"
            ));
        }
    }
    lattice
}

fn set_ortho(metric: &Mat3) -> Mat3 {
    let mut lattice = [[0.0; 3]; 3];
    lattice[0][0] = metric[0][0].sqrt();
    lattice[1][1] = metric[1][1].sqrt();
    lattice[2][2] = metric[2][2].sqrt();
    lattice
}

fn set_tetra(metric: &Mat3) -> Mat3 {
    let mut lattice = [[0.0; 3]; 3];
    let a = metric[0][0].sqrt();
    let b = metric[1][1].sqrt();
    let c = metric[2][2].sqrt();
    lattice[0][0] = (a + b) / 2.0;
    lattice[1][1] = (a + b) / 2.0;
    lattice[2][2] = c;
    lattice
}

fn set_rhomb(metric: &Mat3) -> Mat3 {
    let mut lattice = [[0.0; 3]; 3];
    let a = metric[0][0].sqrt();
    let b = metric[1][1].sqrt();
    let c = metric[2][2].sqrt();
    let angle = ((metric[0][1] / a / b + metric[0][2] / a / c + metric[1][2] / b / c) / 3.0).acos();

    // Reference: http://cst-www.nrl.navy.mil/lattice/struk/rgr.html
    let ahex = 2.0 * (a + b + c) / 3.0 * (angle / 2.0).sin();
    let chex = (a + b + c) / 3.0 * (3.0 * (1.0 + 2.0 * angle.cos())).sqrt();
    lattice[0][0] = ahex / 2.0;
    lattice[1][0] = ahex / (2.0 * (3.0f64).sqrt());
    lattice[2][0] = chex / 3.0;
    lattice[0][1] = -ahex / 2.0;
    lattice[1][1] = ahex / (2.0 * (3.0f64).sqrt());
    lattice[2][1] = chex / 3.0;
    lattice[0][2] = 0.0;
    lattice[1][2] = -ahex / (3.0f64).sqrt();
    lattice[2][2] = chex / 3.0;
    lattice
}

fn set_trigo(metric: &Mat3) -> Mat3 {
    let mut lattice = [[0.0; 3]; 3];
    let a = metric[0][0].sqrt();
    let b = metric[1][1].sqrt();
    let c = metric[2][2].sqrt();
    lattice[0][0] = (a + b) / 2.0;
    lattice[0][1] = -(a + b) / 4.0;
    lattice[1][1] = (a + b) / 4.0 * (3.0f64).sqrt();
    lattice[2][2] = c;
    lattice
}

fn set_cubic(metric: &Mat3) -> Mat3 {
    let mut lattice = [[0.0; 3]; 3];
    let a = metric[0][0].sqrt();
    let b = metric[1][1].sqrt();
    let c = metric[2][2].sqrt();
    let avg = (a + b + c) / 3.0;
    lattice[0][0] = avg;
    lattice[1][1] = avg;
    lattice[2][2] = avg;
    lattice
}

// --- 对称操作精细化 ---

/// 从数据库对称操作恢复原始晶胞中的对称操作。
fn get_refined_symmetry_operations(
    cell: &Cell,
    primitive: &Cell,
    spacegroup: &Spacegroup,
    symprec: f64,
) -> Option<Symmetry> {
    // Primitive symmetry from database
    let conv_sym = spgdb_get_spacegroup_operations(spacegroup.hall_number)?;

    let inv_prim_lat = mat_inverse_matrix_d3(&primitive.lattice, 0.0)?;
    let t_mat = mat_multiply_matrix_d3(&inv_prim_lat, &spacegroup.bravais_lattice);

    let mut conv_sym_mut = conv_sym;
    set_translation_with_origin_shift(&mut conv_sym_mut, &spacegroup.origin_shift);

    let prim_sym = get_primitive_db_symmetry(&t_mat, &conv_sym_mut)?;

    // Input cell symmetry from primitive symmetry
    let t_mat_cell = mat_multiply_matrix_d3(&inv_prim_lat, &cell.lattice);
    let t_mat_int = mat_cast_matrix_3d_to_3i(&t_mat_cell);

    let multiplicity = (cell.size / primitive.size) as i32;
    recover_symmetry_in_original_cell(
        &prim_sym,
        &t_mat_int,
        &cell.lattice,
        multiplicity,
        cell.aperiodic_axis,
        symprec,
    )
}

/// 设置晶体学轨道（等价原子映射）。
fn set_crystallographic_orbits(
    equiv_atoms_cell: &mut [i32],
    primitive: &Cell,
    cell: &Cell,
    equiv_atoms_prim: &[i32],
    mapping_table: &[i32],
) -> bool {
    let mut equiv_atoms = vec![0i32; primitive.size];

    for i in 0..primitive.size {
        for j in 0..cell.size {
            if mapping_table[j] == equiv_atoms_prim[i] {
                equiv_atoms[i] = j as i32;
                break;
            }
        }
    }
    for i in 0..cell.size {
        equiv_atoms_cell[i] = equiv_atoms[mapping_table[i] as usize];
    }

    true
}

/// 处理因超胞导致对称性破缺的等价原子分配。
fn set_equivalent_atoms_broken_symmetry(
    equiv_atoms_cell: &mut [i32],
    cell: &Cell,
    symmetry: &Symmetry,
    mapping_table: &[i32],
    symprec: f64,
) {
    if cell.aperiodic_axis.is_none() {
        for i in 0..cell.size {
            equiv_atoms_cell[i] = i as i32;
            for j in 0..cell.size {
                if mapping_table[i] == mapping_table[j] {
                    if i == j {
                        equiv_atoms_cell[i] =
                            equiv_atoms_cell[search_equivalent_atom(i, cell, symmetry, symprec) as usize];
                    } else {
                        equiv_atoms_cell[i] = equiv_atoms_cell[j];
                    }
                    break;
                }
            }
        }
    } else {
        let aperiodic = cell.aperiodic_axis.unwrap();
        for i in 0..cell.size {
            equiv_atoms_cell[i] = i as i32;
            for j in 0..cell.size {
                if mapping_table[i] == mapping_table[j] {
                    if i == j {
                        equiv_atoms_cell[i] = equiv_atoms_cell[search_layer_equivalent_atom(
                            i, cell, symmetry, aperiodic, symprec,
                        ) as usize];
                    } else {
                        equiv_atoms_cell[i] = equiv_atoms_cell[j];
                    }
                    break;
                }
            }
        }
    }
}

fn search_equivalent_atom(
    atom_index: usize,
    cell: &Cell,
    symmetry: &Symmetry,
    symprec: f64,
) -> usize {
    for i in 0..symmetry.size {
        let mut pos_rot = mat_multiply_matrix_vector_id3(&symmetry.rot[i], &cell.position[atom_index]);
        for j in 0..3 {
            pos_rot[j] += symmetry.trans[i][j];
        }
        for j in 0..atom_index {
            if cel_is_overlap_with_same_type(
                &cell.position[j],
                &pos_rot,
                cell.types[j],
                cell.types[atom_index],
                &cell.lattice,
                symprec,
            ) {
                return j;
            }
        }
    }
    atom_index
}

fn search_layer_equivalent_atom(
    atom_index: usize,
    cell: &Cell,
    symmetry: &Symmetry,
    aperiodic: AperiodicAxis,
    symprec: f64,
) -> usize {
    for i in 0..symmetry.size {
        let mut pos_rot = mat_multiply_matrix_vector_id3(&symmetry.rot[i], &cell.position[atom_index]);
        for j in 0..3 {
            pos_rot[j] += symmetry.trans[i][j];
        }
        for j in 0..atom_index {
            if cel_layer_is_overlap_with_same_type(
                &cell.position[j],
                &pos_rot,
                cell.types[j],
                cell.types[atom_index],
                &cell.lattice,
                aperiodic,
                symprec,
            ) {
                return j;
            }
        }
    }
    atom_index
}

/// 应用原点平移修改对称操作。
fn set_translation_with_origin_shift(conv_sym: &mut Symmetry, origin_shift: &Vec3) {
    // t' = t + (R - E) * w
    for i in 0..conv_sym.size {
        let mut tmp_mat = conv_sym.rot[i];
        tmp_mat[0][0] -= 1;
        tmp_mat[1][1] -= 1;
        tmp_mat[2][2] -= 1;
        let tmp_vec = mat_multiply_matrix_vector_id3(&tmp_mat, origin_shift);
        for j in 0..3 {
            conv_sym.trans[i][j] += tmp_vec[j];
        }
    }
}

/// 从常规对称操作变换到原胞对称操作。
fn get_primitive_db_symmetry(t_mat: &Mat3, conv_sym: &Symmetry) -> Option<Symmetry> {
    let inv_mat = mat_inverse_matrix_d3(t_mat, 0.0)?;

    let mut r_prim: Vec<Mat3I> = Vec::with_capacity(conv_sym.size);
    let mut t_prim: Vec<Vec3> = Vec::with_capacity(conv_sym.size);

    for i in 0..conv_sym.size {
        // Skip duplicate rotations
        let mut is_dup = false;
        for j in 0..i {
            if mat_check_identity_matrix_i3(&conv_sym.rot[i], &conv_sym.rot[j]) {
                is_dup = true;
                break;
            }
        }
        if is_dup {
            continue;
        }

        // R' = T * R * T^-1
        let tmp_mat = mat_multiply_matrix_di3(t_mat, &conv_sym.rot[i]);
        let tmp_mat = mat_multiply_matrix_d3(&tmp_mat, &inv_mat);
        r_prim.push(mat_cast_matrix_3d_to_3i(&tmp_mat));
        // t' = T * t
        t_prim.push(mat_multiply_matrix_vector_d3(t_mat, &conv_sym.trans[i]));
    }

    let num_op = r_prim.len();
    let mut prim_sym = Symmetry::new(num_op);
    for i in 0..num_op {
        prim_sym.rot[i] = r_prim[i];
        prim_sym.trans[i] = t_prim[i];
    }

    Some(prim_sym)
}

// --- 框架和角落计算 ---

fn get_surrounding_frame(t_mat: &Mat3I) -> [i32; 3] {
    let corners = get_corners(t_mat);

    let mut frame = [0; 3];
    for i_axis in 0..3 {
        let mut max = corners[i_axis][0];
        let mut min = corners[i_axis][0];
        for j in 1..8 {
            if max < corners[i_axis][j] {
                max = corners[i_axis][j];
            }
            if min > corners[i_axis][j] {
                min = corners[i_axis][j];
            }
        }
        frame[i_axis] = max - min;
    }
    frame
}

fn get_corners(t_mat: &Mat3I) -> [[i32; 8]; 3] {
    let mut corners = [[0i32; 8]; 3];
    // O
    for i in 0..3 { corners[i][0] = 0; }
    // a, b, c
    for i in 0..3 {
        for j in 0..3 { corners[j][i + 1] = t_mat[j][i]; }
    }
    // b+c, c+a, a+b
    for i in 0..3 {
        for j in 0..3 {
            corners[j][i + 4] = t_mat[j][(i + 1) % 3] + t_mat[j][(i + 2) % 3];
        }
    }
    // a+b+c
    for i in 0..3 { corners[i][7] = t_mat[i][0] + t_mat[i][1] + t_mat[i][2]; }
    corners
}

// --- 原始晶胞对称操作恢复 ---

fn recover_symmetry_in_original_cell(
    prim_sym: &Symmetry,
    t_mat: &Mat3I,
    lattice: &Mat3,
    multiplicity: i32,
    aperiodic_axis: Option<AperiodicAxis>,
    symprec: f64,
) -> Option<Symmetry> {
    let mut frame = [0i32; 3];
    frame = get_surrounding_frame(t_mat);

    let t_mat_d = mat_cast_matrix_3i_to_3d(t_mat);
    let inv_tmat = mat_inverse_matrix_d3(&t_mat_d, 0.0)?;

    let lattice_trans = get_lattice_translations(&frame, &inv_tmat);
    let pure_trans = remove_overlapping_lattice_points(lattice, &lattice_trans, symprec)?;

    let t_sym = get_symmetry_in_original_cell(t_mat, &inv_tmat, lattice, prim_sym, symprec)?;

    if pure_trans.len() == multiplicity as usize {
        copy_symmetry_upon_lattice_points(&pure_trans, &t_sym, aperiodic_axis)
    } else {
        None
    }
}

fn get_lattice_translations(frame: &[i32; 3], inv_tmat: &Mat3) -> Vec<Vec3> {
    let total = frame[0] as usize * frame[1] as usize * frame[2] as usize;
    let mut result = Vec::with_capacity(total);

    for i in 0..frame[0] {
        for j in 0..frame[1] {
            for k in 0..frame[2] {
                let v = [i as f64, j as f64, k as f64];
                let tv = mat_multiply_matrix_vector_d3(inv_tmat, &v);
                result.push([mat_dmod1(tv[0]), mat_dmod1(tv[1]), mat_dmod1(tv[2])]);
            }
        }
    }

    result
}

fn remove_overlapping_lattice_points(
    lattice: &Mat3,
    lattice_trans: &[Vec3],
    symprec: f64,
) -> Option<Vec<Vec3>> {
    let mut pure_trans: Vec<Vec3> = Vec::new();

    for i in 0..lattice_trans.len() {
        let mut is_found = false;
        for j in 0..pure_trans.len() {
            if cel_is_overlap(&lattice_trans[i], &pure_trans[j], lattice, symprec) {
                is_found = true;
                break;
            }
        }
        if !is_found {
            pure_trans.push(lattice_trans[i]);
        }
    }

    if pure_trans.is_empty() {
        None
    } else {
        Some(pure_trans)
    }
}

fn get_symmetry_in_original_cell(
    t_mat: &Mat3I,
    inv_tmat: &Mat3,
    lattice: &Mat3,
    prim_sym: &Symmetry,
    symprec: f64,
) -> Option<Symmetry> {
    let mut t_sym = Symmetry::new(prim_sym.size);
    let mut size_sym_orig = 0usize;

    for i in 0..prim_sym.size {
        // R' = T^-1 * R * T
        let tmp_mat = mat_multiply_matrix_di3(inv_tmat, &prim_sym.rot[i]);
        let tmp_rot_d = mat_multiply_matrix_di3(&tmp_mat, t_mat);

        let tmp_rot_i = mat_cast_matrix_3d_to_3i(&tmp_rot_d);
        let tmp_lat_i = mat_multiply_matrix_di3(lattice, &tmp_rot_i);
        let tmp_lat_d = mat_multiply_matrix_d3(lattice, &tmp_rot_d);

        if mat_check_identity_matrix_d3(&tmp_lat_i, &tmp_lat_d, symprec) {
            t_sym.rot[size_sym_orig] = tmp_rot_i;
            // t' = T^-1 * t
            t_sym.trans[size_sym_orig] =
                mat_multiply_matrix_vector_d3(inv_tmat, &prim_sym.trans[i]);
            size_sym_orig += 1;
        }
    }

    // Trim if symmetry was broken
    if size_sym_orig != prim_sym.size {
        let mut t_red_sym = Symmetry::new(size_sym_orig);
        for i in 0..size_sym_orig {
            t_red_sym.rot[i] = t_sym.rot[i];
            t_red_sym.trans[i] = t_sym.trans[i];
        }
        t_sym = t_red_sym;
    }

    // Update size
    t_sym.size = size_sym_orig;

    Some(t_sym)
}

fn copy_symmetry_upon_lattice_points(
    pure_trans: &[Vec3],
    t_sym: &Symmetry,
    aperiodic_axis: Option<AperiodicAxis>,
) -> Option<Symmetry> {
    let size_sym_orig = t_sym.size;
    let total = pure_trans.len() * size_sym_orig;
    let mut symmetry = Symmetry::new(total);

    for i in 0..pure_trans.len() {
        for j in 0..size_sym_orig {
            let idx = size_sym_orig * i + j;
            symmetry.rot[idx] = t_sym.rot[j];
            symmetry.trans[idx] = t_sym.trans[j];
            for k in 0..3 {
                symmetry.trans[idx][k] += pure_trans[i][k];
                if k as i32 != aperiodic_axis.map_or(-1, |ap| ap.axis_index() as i32) {
                    symmetry.trans[idx][k] = mat_dmod1(symmetry.trans[idx][k]);
                }
            }
        }
    }

    Some(symmetry)
}

// --- Bravais 格子相似性匹配 ---

/// 寻找最接近标准晶格的 Bravais 格子取向。
///
/// 重写 `spacegroup.bravais_lattice` 和 `spacegroup.origin_shift`。
/// 返回 false 表示失败。
pub fn ref_find_similar_bravais_lattice(spacegroup: &mut Spacegroup, symprec: f64) -> bool {
    let conv_sym = match spgdb_get_spacegroup_operations(spacegroup.hall_number) {
        Some(s) => s,
        None => return false,
    };

    let mut std_lattice = [[0.0; 3]; 3];
    std_lattice = ref_get_conventional_lattice(spacegroup);

    // Find best rotation
    let mut min_length2 = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            min_length2 += spacegroup.bravais_lattice[i][j] * spacegroup.bravais_lattice[i][j];
        }
    }
    let mut min_length = min_length2.sqrt();
    let mut rot_i: i32 = -1;
    let mut rot_lat = [[0.0; 3]; 3];

    for i in 0..conv_sym.size {
        if mat_get_determinant_i3(&conv_sym.rot[i]) < 0 {
            continue;
        }

        // (a_s', b_s', c_s') = (a_s, b_s, c_s) * Rs
        let tmp_mat = mat_multiply_matrix_di3(&spacegroup.bravais_lattice, &conv_sym.rot[i]);
        let mut length2 = 0.0;
        for j in 0..3 {
            for k in 0..3 {
                let diff = tmp_mat[j][k] - std_lattice[j][k];
                length2 += diff * diff;
            }
        }
        let length = length2.sqrt();
        if length < min_length - symprec {
            rot_lat = tmp_mat;
            rot_i = i as i32;
            min_length = length;
        }
    }

    // Find best origin shift
    let lattice_rank: i32 = if spacegroup.hall_number > 0 { 3 } else { 2 };
    if rot_i > -1 {
        let mut min_len = 2.0f64;
        let mut shortest_p = [0.0; 3];
        for i in 0..conv_sym.size {
            if !mat_check_identity_matrix_i3(&conv_sym.rot[i], &conv_sym.rot[rot_i as usize]) {
                continue;
            }
            let tmp_mat_d = mat_cast_matrix_3i_to_3d(&conv_sym.rot[i]);
            let tmp_inv = mat_inverse_matrix_d3(&tmp_mat_d, 0.0);
            let Some(tmp_inv) = tmp_inv else { continue; };

            let mut p = mat_multiply_matrix_vector_d3(&tmp_inv, &spacegroup.origin_shift);
            let tmp_vec = mat_multiply_matrix_vector_d3(&tmp_inv, &conv_sym.trans[i]);
            for j in 0..lattice_rank as usize {
                p[j] -= tmp_vec[j];
                p[j] -= (p[j]).round(); // mat_Nint equivalent for f64
            }
            for j in lattice_rank as usize..3 {
                p[j] -= tmp_vec[j];
            }
            let length = mat_norm_squared_d3(&p).sqrt();
            if length < min_len - symprec {
                min_len = length;
                for j in 0..lattice_rank as usize {
                    p[j] = mat_dmod1(p[j]);
                }
                shortest_p = p;
            }
        }
        spacegroup.origin_shift = shortest_p;
        spacegroup.bravais_lattice = rot_lat;
    }

    true
}

/// 计算 rigid rotation: std_lattice = rotation @ bravais_lattice
fn measure_rigid_rotation(bravais_lattice: &Mat3, std_lattice: &Mat3) -> Mat3 {
    let brv_basis = get_orthonormal_basis(bravais_lattice);
    let std_basis = get_orthonormal_basis(std_lattice);
    let inv_brv_basis = mat_inverse_matrix_d3(&brv_basis, 0.0);
    let mut rotation = [[0.0; 3]; 3];
    if let Some(inv) = inv_brv_basis {
        rotation = mat_multiply_matrix_d3(&std_basis, &inv);
    }
    rotation
}

fn get_orthonormal_basis(lattice: &Mat3) -> Mat3 {
    let lattice_t = mat_transpose_matrix_d3(lattice);
    let mut basis_t = [[0.0; 3]; 3];

    basis_t[0] = lattice_t[0];
    basis_t[2] = mat_cross_product_d3(&lattice_t[0], &lattice_t[1]);
    basis_t[1] = mat_cross_product_d3(&basis_t[2], &lattice_t[0]);

    for i in 0..3 {
        let length = mat_norm_squared_d3(&basis_t[i]).sqrt();
        for j in 0..3 {
            basis_t[i][j] /= length;
        }
    }

    mat_transpose_matrix_d3(&basis_t)
}
