//! 空间群确定的主流程。
//!
//! 组合原胞查找、空间群搜索和结构精细化，
//! 通过递减容差的重试机制提高搜索成功率。

use crate::SymError;
use crate::cell::Cell;
use crate::debug;
use crate::primitive::{Primitive, get_primitive};
use crate::refinement::{ExactStructure, get_exact_structure_and_symmetry};
use crate::spacegroup::{Spacegroup, search_spacegroup};

const REDUCE_RATE_OUTER: f64 = 0.9;
const NUM_ATTEMPT_OUTER: i32 = 10;
const REDUCE_RATE: f64 = 0.95;
const ANGLE_REDUCE_RATE: f64 = 0.95;
const NUM_ATTEMPT: i32 = 20;

/// 空间群确定过程的中间结果容器。
pub struct DataContainer {
    pub spacegroup: Option<Spacegroup>,
    pub primitive: Option<Primitive>,
    pub exact_structure: Option<ExactStructure>,
}

/// 确定空间群的主入口。
///
/// 先通过 `get_spacegroup_and_primitive` 找到空间群和原胞，
/// 再通过 `get_exact_structure_and_symmetry` 精细化结构。
/// 如果失败，逐步降低容差重试。
pub fn determine_all(
    cell: &Cell,
    hall_number: i32,
    symprec: f64,
    angle_symprec: f64,
) -> Result<DataContainer, SymError> {
    if hall_number > 530 {
        return Err(SymError::SpacegroupSearchFailed);
    }

    let mut tolerance = symprec;
    for _ in 0..NUM_ATTEMPT_OUTER {
        if let Ok(mut container) =
            get_spacegroup_and_primitive(cell, hall_number, tolerance, angle_symprec)
        {
            let exstr = {
                let sg = container
                    .spacegroup
                    .as_mut()
                    .ok_or(SymError::SpacegroupSearchFailed)?;
                let prim = container
                    .primitive
                    .as_ref()
                    .ok_or(SymError::SpacegroupSearchFailed)?;
                let prim_cell = prim.cell.as_ref().ok_or(SymError::SpacegroupSearchFailed)?;
                get_exact_structure_and_symmetry(
                    sg,
                    prim_cell,
                    cell,
                    &prim.mapping_table,
                    prim.tolerance,
                )
            };
            if let Some(exstr) = exstr {
                container.exact_structure = Some(exstr);
                return Ok(container);
            }
            debug::debug_print(format_args!(
                "spglib: get_exact_structure_and_symmetry failed.\n"
            ));
        }
        tolerance *= REDUCE_RATE_OUTER;
    }

    Err(SymError::SpacegroupSearchFailed)
}

/// 在递减容差下寻找空间群和原胞。
fn get_spacegroup_and_primitive(
    cell: &Cell,
    hall_number: i32,
    symprec: f64,
    angle_symprec: f64,
) -> Result<DataContainer, SymError> {
    debug::debug_print(format_args!(
        "get_spacegroup_and_primitive (tolerance = {}):\n",
        symprec
    ));

    let mut tolerance = symprec;
    let mut angle_tolerance = angle_symprec;

    for attempt in 0..NUM_ATTEMPT {
        let primitive = get_primitive(cell, tolerance, angle_tolerance).ok();
        if let Some(primitive) = primitive {
            debug::debug_print(format_args!("primitive lattice\n"));

            let prim_tol = primitive.tolerance;
            let prim_angle_tol = primitive.angle_tolerance;

            let spacegroup =
                search_spacegroup(&primitive, hall_number, prim_tol, prim_angle_tol).ok();
            if let Some(spacegroup) = spacegroup {
                return Ok(DataContainer {
                    spacegroup: Some(spacegroup),
                    primitive: Some(primitive),
                    exact_structure: None,
                });
            }
        }

        debug::debug_print(format_args!(
            "spglib: Attempt {} tolerance = {} failed.\n",
            attempt, tolerance
        ));

        tolerance *= REDUCE_RATE;
        if angle_tolerance > 0.0 {
            angle_tolerance *= ANGLE_REDUCE_RATE;
        }
    }

    Err(SymError::SpacegroupSearchFailed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::Cell;

    #[test]
    fn test_det_hall_number_out_of_range() {
        let cell = Cell::new(1, crate::cell::TensorRank::NoSpin);
        let result = determine_all(&cell, 999, 1e-5, -1.0);
        assert!(result.is_err());
    }
}
