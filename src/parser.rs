//! Parser for POSCAR-like structure files.

use crate::mathfunc::{Mat3, Vec3};

/// Parse a POSCAR-format string into `(lattice, positions, types, magnetic_moments)`.
///
/// Format:
/// ```text
/// comment line
/// scale_factor
/// a1x a1y a1z
/// a2x a2y a2z
/// a3x a3y a3z
/// atom_types  (e.g. "Fe O")
/// atom_counts (e.g. "2 1")
/// Direct|Cartesian
/// x y z [mx my mz]  # positions, optional 3 magnetic moment components
/// ```
///
/// Returns `None` on malformed input.
pub fn parse_poscar(data: &str) -> Option<(Mat3, Vec<Vec3>, Vec<i32>, Option<Vec<[f64; 3]>>)> {
    let lines: Vec<&str> = data.lines().collect();
    if lines.len() < 6 {
        return None;
    }

    let scale: f64 = lines.get(1)?.trim().parse().ok()?;

    let mut rows = [[0.0; 3]; 3];
    for i in 0..3 {
        let parts: Vec<f64> = lines[i + 2]
            .split_whitespace()
            .filter_map(|x| x.parse().ok())
            .collect();
        if parts.len() < 3 {
            return None;
        }
        rows[i] = [parts[0], parts[1], parts[2]];
    }
    let mut lattice = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            lattice[i][j] = rows[j][i];
        }
    }
    if scale != 1.0 {
        for i in 0..3 {
            for j in 0..3 {
                lattice[i][j] *= scale;
            }
        }
    }

    let type_line = lines.get(5)?;
    let counts: Vec<i32> = lines
        .get(6)?
        .split_whitespace()
        .filter_map(|x| x.parse().ok())
        .collect();
    if counts.is_empty() {
        return None;
    }
    let n_atoms: usize = counts.iter().map(|&c| c as usize).sum();

    let type_names: Vec<&str> = type_line.split_whitespace().collect();
    let mut types = Vec::with_capacity(n_atoms);
    for (ti, &cnt) in counts.iter().enumerate() {
        let type_num = if type_names.len() > ti && type_names[ti].parse::<i32>().is_err() {
            element_to_number(type_names[ti])
        } else {
            (ti + 1) as i32
        };
        for _ in 0..cnt {
            types.push(type_num);
        }
    }

    let mode_line = lines.get(7)?;
    let is_cartesian = mode_line.trim().to_uppercase().starts_with('C')
        || mode_line.trim().to_uppercase().starts_with('K');

    let mut positions = Vec::with_capacity(n_atoms);
    let mut moments = Vec::with_capacity(n_atoms);
    let mut has_moments = false;

    for i in 0..n_atoms {
        let line = lines.get(8 + i)?;
        let parts: Vec<f64> = line
            .split_whitespace()
            .filter_map(|x| x.parse().ok())
            .collect();
        if parts.len() < 3 {
            return None;
        }
        positions.push([parts[0], parts[1], parts[2]]);
        if parts.len() >= 6 {
            moments.push([parts[3], parts[4], parts[5]]);
            has_moments = true;
        } else {
            moments.push([0.0; 3]);
        }
    }

    if is_cartesian {
        let inv_lat = crate::mathfunc::mat_inverse_matrix_d3(&lattice, 1e-10).ok()?;
        for i in 0..n_atoms {
            let frac = crate::mathfunc::mat_multiply_matrix_vector_d3(&inv_lat, &positions[i]);
            positions[i] = frac;
        }
    }

    let mag_opt = if has_moments { Some(moments) } else { None };
    Some((lattice, positions, types, mag_opt))
}

fn element_to_number(symbol: &str) -> i32 {
    match symbol.trim() {
        "H" => 1, "He" => 2,
        "Li" => 3, "Be" => 4, "B" => 5, "C" => 6, "N" => 7, "O" => 8, "F" => 9, "Ne" => 10,
        "Na" => 11, "Mg" => 12, "Al" => 13, "Si" => 14, "P" => 15, "S" => 16, "Cl" => 17, "Ar" => 18,
        "K" => 19, "Ca" => 20, "Sc" => 21, "Ti" => 22, "V" => 23, "Cr" => 24, "Mn" => 25,
        "Fe" => 26, "Co" => 27, "Ni" => 28, "Cu" => 29, "Zn" => 30, "Ga" => 31, "Ge" => 32,
        "As" => 33, "Se" => 34, "Br" => 35, "Kr" => 36,
        "Rb" => 37, "Sr" => 38, "Y" => 39, "Zr" => 40, "Nb" => 41, "Mo" => 42,
        "Tc" => 43, "Ru" => 44, "Rh" => 45, "Pd" => 46, "Ag" => 47, "Cd" => 48,
        "In" => 49, "Sn" => 50, "Sb" => 51, "Te" => 52, "I" => 53, "Xe" => 54,
        "Cs" => 55, "Ba" => 56,
        "La" => 57, "Ce" => 58, "Pr" => 59, "Nd" => 60, "Pm" => 61, "Sm" => 62,
        "Eu" => 63, "Gd" => 64, "Tb" => 65, "Dy" => 66, "Ho" => 67, "Er" => 68,
        "Tm" => 69, "Yb" => 70, "Lu" => 71,
        "Hf" => 72, "Ta" => 73, "W" => 74, "Re" => 75, "Os" => 76, "Ir" => 77,
        "Pt" => 78, "Au" => 79, "Hg" => 80, "Tl" => 81, "Pb" => 82, "Bi" => 83,
        "Po" => 84, "At" => 85, "Rn" => 86,
        "Fr" => 87, "Ra" => 88,
        "Ac" => 89, "Th" => 90, "Pa" => 91, "U" => 92, "Np" => 93, "Pu" => 94,
        "Am" => 95, "Cm" => 96, "Bk" => 97, "Cf" => 98, "Es" => 99, "Fm" => 100,
        "Md" => 101, "No" => 102, "Lr" => 103,
        "Rf" => 104, "Db" => 105, "Sg" => 106, "Bh" => 107, "Hs" => 108, "Mt" => 109,
        "Ds" => 110, "Rg" => 111, "Cn" => 112, "Nh" => 113, "Fl" => 114, "Mc" => 115,
        "Lv" => 116, "Ts" => 117, "Og" => 118,
        _ => 1,
    }
}
