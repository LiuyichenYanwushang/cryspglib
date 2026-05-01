//! Rust-idiomatic public API for cryspglib.
//!
//! # Quick start
//!
//! ```no_run
//! use cryspglib::{Crystal, SymError};
//!
//! let cry = Crystal::new(
//!     [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
//!     vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
//!     vec![26, 26],
//! );
//! let ds = cry.analyze().symprec(1e-5).dataset()?;
//! println!("Space group #{}: {}", ds.spacegroup_number, ds.international_symbol);
//! # Ok::<(), SymError>(())
//! ```

use crate::cell::{AperiodicAxis, Cell, TensorRank};
use crate::debug;
use crate::delaunay::del_delaunay_reduce;
use crate::determination::det_determine_all;
use crate::mathfunc::{Mat3, Mat3I, Vec3};
use crate::niggli::niggli_reduce;
use crate::pointgroup::ptg_get_pointgroup;
use crate::primitive::Primitive;
use crate::spacegroup::Spacegroup;
use crate::spg_database::{Centering, spgdb_get_spacegroup_type};
use crate::{SpaceGroup, SymError};

// ── Crystal ──────────────────────────────────────────────────────────────────

/// Crystal structure: lattice + atomic positions + optional magnetic moments.
///
/// This is the primary entry point for all symmetry analysis.
///
/// # Lattice convention
/// `lattice[cart][vec]`: rows = Cartesian components (x, y, z), columns = lattice vectors (a, b, c).
///
/// # Examples
///
/// ```
/// use cryspglib::Crystal;
///
/// let si = Crystal::new(
///     [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
///     vec![[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
///     vec![14, 14],
/// );
/// ```
#[derive(Debug, Clone)]
pub struct Crystal {
    /// Lattice matrix, layout `[cart][vec]`
    pub lattice: Mat3,
    /// Atomic positions in fractional coordinates
    pub positions: Vec<Vec3>,
    /// Atomic numbers (e.g., 14 for Si, 26 for Fe)
    pub types: Vec<i32>,
    /// Magnetic moments per atom (`[mx, my, mz]`). `None` = non-magnetic.
    pub moments: Option<Vec<[f64; 3]>>,
    /// Aperiodic axis for 2D slabs. `None` = full 3D periodicity.
    pub aperiodic_axis: Option<AperiodicAxis>,
}

impl Crystal {
    /// Create a non-magnetic 3D crystal.
    ///
    /// # Panics
    /// Panics if `positions.len() != types.len()`.
    pub fn new(lattice: Mat3, positions: Vec<Vec3>, types: Vec<i32>) -> Self {
        assert_eq!(
            positions.len(),
            types.len(),
            "positions and types must have the same length"
        );
        Crystal {
            lattice,
            positions,
            types,
            moments: None,
            aperiodic_axis: None,
        }
    }

    /// Add collinear magnetic moments (one `[mx, my, mz]` per atom).
    pub fn with_magnetic(mut self, moments: Vec<[f64; 3]>) -> Self {
        assert_eq!(moments.len(), self.positions.len());
        self.moments = Some(moments);
        self
    }

    /// Mark as a 2D slab with the given aperiodic axis.
    pub fn with_layer(mut self, axis: AperiodicAxis) -> Self {
        self.aperiodic_axis = Some(axis);
        self
    }

    /// Number of atoms.
    pub fn natom(&self) -> usize {
        self.positions.len()
    }

    /// Begin symmetry analysis with default settings.
    ///
    /// Returns a [`SymmetryAnalysis`] builder that can be configured before
    /// executing any terminal operation.
    pub fn analyze(&self) -> SymmetryAnalysis<'_> {
        SymmetryAnalysis {
            crystal: self,
            symprec: 1e-5,
            angle_tolerance: -1.0,
        }
    }

    /// Delaunay lattice reduction.
    pub fn delaunay_reduce(&self, symprec: f64) -> Result<Mat3, SymError> {
        del_delaunay_reduce(&self.lattice, symprec).ok_or(SymError::DelaunayFailed)
    }

    /// Niggli lattice reduction.
    pub fn niggli_reduce(&self, symprec: f64) -> Result<Mat3, SymError> {
        let mut reduced = self.lattice;
        if niggli_reduce(&mut reduced, symprec, None) {
            Ok(reduced)
        } else {
            Err(SymError::NiggliFailed)
        }
    }

    // ── Internal: convert to Cell ──────────────────────────────────────────

    pub(crate) fn to_cell(&self) -> Cell {
        let tensor_rank = self.tensor_rank();
        let mut cell = Cell::new(self.positions.len(), tensor_rank);
        if self.aperiodic_axis().is_none() {
            cell.set_cell(&self.lattice, &self.positions, &self.types);
        } else {
            cell.set_layer_cell(
                &self.lattice,
                &self.positions,
                &self.types,
                self.aperiodic_axis,
            );
        }
        if let Some(ref moments) = self.moments {
            let tensors: Vec<f64> = moments.iter().flat_map(|m| [m[0], m[1], m[2]]).collect();
            cell.set_cell_with_tensors(&self.lattice, &self.positions, &self.types, &tensors);
        }
        cell
    }

    pub(crate) fn tensor_rank(&self) -> TensorRank {
        match self.moments {
            None => TensorRank::NoSpin,
            Some(_) => TensorRank::NonCollinear,
        }
    }

    pub(crate) fn aperiodic_axis(&self) -> Option<AperiodicAxis> {
        self.aperiodic_axis
    }

    // ── POSCAR parser ───────────────────────────────────────────────────────

    /// Parse a POSCAR-format string into a `Crystal`.
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
    pub fn from_poscar(data: &str) -> Option<Self> {
        crate::parser::parse_poscar(data).map(
            |(lattice, positions, types, moments)| Crystal {
                lattice,
                positions,
                types,
                moments,
                aperiodic_axis: None,
            },
        )
    }
}

// ── SymmetryAnalysis ─────────────────────────────────────────────────────────

/// Builder for symmetry analysis of a [`Crystal`].
///
/// Configure analysis parameters, then call a terminal method:
///
/// ```no_run
/// # use cryspglib::Crystal;
/// let cry = Crystal::new([[1.;3];3], vec![[0.;3]], vec![14]);
/// let ds = cry.analyze()
///     .symprec(1e-5)
///     .angle_tolerance(-0.1)
///     .dataset()?;
/// # Ok::<(), cryspglib::SymError>(())
/// ```
pub struct SymmetryAnalysis<'a> {
    crystal: &'a Crystal,
    symprec: f64,
    angle_tolerance: f64,
}

impl<'a> SymmetryAnalysis<'a> {
    /// Set symmetry tolerance (Cartesian distance, default `1e-5`).
    pub fn symprec(mut self, val: f64) -> Self {
        self.symprec = val;
        self
    }

    /// Set angle tolerance in radians (default `-1.0` = auto).
    pub fn angle_tolerance(mut self, val: f64) -> Self {
        self.angle_tolerance = val;
        self
    }

    // ── Terminal methods ────────────────────────────────────────────────────

    /// Full space group dataset.
    pub fn dataset(&self) -> Result<SpaceGroup, SymError> {
        let cell = self.crystal.to_cell();
        get_dataset_inner(&cell, self.symprec, self.angle_tolerance, 0)
    }

    /// Symmetry operations only (rotations + translations).
    pub fn symmetry(&self) -> Result<SymmetryOps, SymError> {
        let ds = self.dataset()?;
        let ops: Vec<SymmetryOp> = (0..ds.n_operations)
            .map(|i| SymmetryOp {
                rotation: ds.rotations[i],
                translation: ds.translations[i],
                time_reversal: false,
            })
            .collect();
        Ok(SymmetryOps { operations: ops })
    }

    /// Primitive cell.
    pub fn primitive_cell(&self) -> Result<Crystal, SymError> {
        let cell = self.crystal.to_cell();
        let prim_cell = standardize_primitive_inner(&cell, self.symprec, self.angle_tolerance)?;
        Ok(cell_to_crystal(&prim_cell))
    }

    /// Standardize cell: `to_primitive` returns primitive cell; `no_idealize` skips position idealization.
    pub fn standardize(&self, to_primitive: bool, no_idealize: bool) -> Result<Crystal, SymError> {
        let cell = self.crystal.to_cell();
        let cc = standardize_cell_inner(&cell, to_primitive, no_idealize, self.symprec, self.angle_tolerance)?;
        Ok(cell_to_crystal(&cc))
    }

    /// Space group hall number.
    pub fn hall_number(&self) -> Result<usize, SymError> {
        let ds = self.dataset()?;
        Ok(ds.hall_number)
    }

    /// Space group international symbol.
    pub fn international(&self) -> Result<(usize, String), SymError> {
        let ds = self.dataset()?;
        if ds.spacegroup_number > 0 {
            Ok((ds.spacegroup_number, ds.international_symbol))
        } else {
            Err(SymError::SpacegroupSearchFailed)
        }
    }

    /// Irreducible k-point mesh.
    pub fn irreducible_mesh(
        &self,
        mesh: [i32; 3],
        is_shift: [i32; 3],
        time_reversal: bool,
    ) -> Result<IrMesh, SymError> {
        let ds = self.dataset()?;
        use crate::mathfunc::MatINT;
        let mut rotations = MatINT::new(ds.n_operations);
        for i in 0..ds.n_operations {
            rotations.mat[i] = ds.rotations[i];
        }

        let rot_reciprocal = crate::kpoint::kpt_get_point_group_reciprocal(
            &rotations,
            if time_reversal { 1 } else { 0 },
        )
        .ok_or(SymError::SpacegroupSearchFailed)?;

        let total = (mesh[0] as usize) * (mesh[1] as usize) * (mesh[2] as usize);
        let mut grid_address = vec![[0i32; 3]; total];
        let mut mapping_table = vec![0usize; total];

        let num_ir = crate::kpoint::kpt_get_irreducible_reciprocal_mesh(
            &mut grid_address,
            &mut mapping_table,
            &mesh,
            &is_shift,
            &rot_reciprocal,
        );

        Ok(IrMesh {
            grid_addresses: grid_address,
            mapping_table,
            num_ir,
        })
    }
}

// ── New output types ─────────────────────────────────────────────────────────

/// A single symmetry operation {R|t}[θ] with optional time reversal.
#[derive(Debug, Clone, Copy)]
pub struct SymmetryOp {
    /// Integer rotation matrix (3×3)
    pub rotation: Mat3I,
    /// Fractional translation vector
    pub translation: Vec3,
    /// Time reversal: false = ordinary, true = primed (anti-unitary)
    pub time_reversal: bool,
}

/// Ordered set of symmetry operations.
#[derive(Debug, Clone)]
pub struct SymmetryOps {
    pub operations: Vec<SymmetryOp>,
}

/// Irreducible k-point mesh.
#[derive(Debug, Clone)]
pub struct IrMesh {
    /// Grid point addresses in fractional coordinates
    pub grid_addresses: Vec<[i32; 3]>,
    /// Full grid index → irreducible grid index mapping
    pub mapping_table: Vec<usize>,
    /// Number of irreducible points
    pub num_ir: usize,
}

// ── Internal helpers ─────────────────────────────────────────────────────────

fn cell_to_crystal(cell: &Cell) -> Crystal {
    Crystal {
        lattice: cell.lattice,
        positions: cell.position.clone(),
        types: cell.types.clone(),
        moments: None,
        aperiodic_axis: cell.aperiodic_axis,
    }
}

fn get_dataset_inner(
    cell: &Cell,
    symprec: f64,
    angle_tolerance: f64,
    hall_number: i32,
) -> Result<SpaceGroup, SymError> {
    let container = det_determine_all(cell, hall_number, symprec, angle_tolerance)
        .ok_or(SymError::SpacegroupSearchFailed)?;

    let spacegroup = container
        .spacegroup
        .as_ref()
        .ok_or(SymError::SpacegroupSearchFailed)?;
    let primitive = container
        .primitive
        .as_ref()
        .ok_or(SymError::SpacegroupSearchFailed)?;
    let exstr = container
        .exact_structure
        .as_ref()
        .ok_or(SymError::SpacegroupSearchFailed)?;

    let dataset = build_dataset(cell, primitive, spacegroup, exstr)
        .ok_or(SymError::SpacegroupSearchFailed)?;
    Ok(dataset)
}

fn build_dataset(
    cell: &Cell,
    primitive: &Primitive,
    spacegroup: &Spacegroup,
    exstr: &crate::refinement::ExactStructure,
) -> Option<SpaceGroup> {
    let n_atoms = cell.size;
    let n_operations = exstr.symmetry.size;

    let mut dataset = SpaceGroup {
        spacegroup_number: spacegroup.number,
        hall_number: spacegroup.hall_number,
        international_symbol: spacegroup.international_short.clone(),
        hall_symbol: spacegroup.hall_symbol.clone(),
        choice: spacegroup.choice.clone(),
        transformation_matrix: [[0.0; 3]; 3],
        origin_shift: spacegroup.origin_shift,
        n_operations,
        rotations: vec![[[0; 3]; 3]; n_operations],
        translations: vec![[0.0; 3]; n_operations],
        n_atoms,
        wyckoffs: vec![0i32; n_atoms],
        site_symmetry_symbols: vec![String::new(); n_atoms],
        equivalent_atoms: vec![0i32; n_atoms],
        crystallographic_orbits: vec![0i32; n_atoms],
        mapping_to_primitive: vec![0i32; n_atoms],
        n_std_atoms: exstr.bravais.size,
        std_lattice: exstr.bravais.lattice,
        std_positions: exstr.bravais.position.clone(),
        std_types: exstr.bravais.types.clone(),
        std_rotation_matrix: [[0.0; 3]; 3],
        std_mapping_to_primitive: vec![0i32; exstr.bravais.size],
        primitive_lattice: [[0.0; 3]; 3],
        pointgroup_symbol: String::new(),
    };

    let inv_lat =
        crate::mathfunc::mat_inverse_matrix_d3(&spacegroup.bravais_lattice, 0.0).ok()?;
    dataset.transformation_matrix =
        crate::mathfunc::mat_multiply_matrix_d3(&inv_lat, &cell.lattice);

    for i in 0..n_operations {
        dataset.rotations[i] = exstr.symmetry.rot[i];
        dataset.translations[i] = exstr.symmetry.trans[i];
    }

    for i in 0..n_atoms {
        dataset.wyckoffs[i] = exstr.wyckoffs[i];
        dataset.site_symmetry_symbols[i] = exstr.site_symmetry_symbols[i].clone();
        dataset.equivalent_atoms[i] = exstr.equivalent_atoms[i];
        dataset.crystallographic_orbits[i] = exstr.crystallographic_orbits[i];
    }

    if let Some(prim_cell) = &primitive.cell {
        dataset.primitive_lattice = prim_cell.lattice;
    }
    for i in 0..n_atoms {
        dataset.mapping_to_primitive[i] = primitive.mapping_table[i];
    }

    for i in 0..dataset.n_std_atoms {
        dataset.std_mapping_to_primitive[i] = exstr.std_mapping_to_primitive[i];
    }
    dataset.std_rotation_matrix = exstr.rotation;

    let pointgroup = ptg_get_pointgroup(spacegroup.pointgroup_number);
    dataset.pointgroup_symbol = pointgroup.symbol.to_string();

    Some(dataset)
}

fn standardize_primitive_inner(
    cell: &Cell,
    symprec: f64,
    angle_tolerance: f64,
) -> Result<Cell, SymError> {
    let dataset = get_dataset_inner(cell, symprec, angle_tolerance, 0)?;
    let centering = spgdb_get_spacegroup_type(dataset.hall_number).centering;

    let mut bravais = Cell::new(dataset.n_std_atoms, TensorRank::NoSpin);
    bravais.lattice = dataset.std_lattice;
    for i in 0..dataset.n_std_atoms {
        bravais.types[i] = dataset.std_types[i];
        bravais.position[i] = dataset.std_positions[i];
    }

    let mut mapping_table = vec![0usize; bravais.size];
    let identity: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let primitive = crate::spacegroup::spa_transform_to_primitive(
        &mut mapping_table,
        &bravais,
        &identity,
        centering,
        symprec,
    )
    .ok_or(SymError::CellStandardizationFailed)?;

    for i in 0..primitive.size {
        if mapping_table[i] != i {
            debug::warning_print(format_args!(
                "spglib: spa_transform_to_primitive failed ({} != {})\n",
                mapping_table[i], i
            ));
            return Err(SymError::CellStandardizationFailed);
        }
    }
    Ok(primitive)
}

fn standardize_cell_inner(
    cell: &Cell,
    to_primitive: bool,
    no_idealize: bool,
    symprec: f64,
    angle_tolerance: f64,
) -> Result<Cell, SymError> {
    if to_primitive && !no_idealize {
        return standardize_primitive_inner(cell, symprec, angle_tolerance);
    }

    let dataset = get_dataset_inner(cell, symprec, angle_tolerance, 0)?;

    if to_primitive && no_idealize {
        // Use existing standardize logic with dataset
        let centering = spgdb_get_spacegroup_type(dataset.hall_number).centering;
        let num_atom = cell.size;
        let mut work_cell = Cell::new(num_atom, TensorRank::NoSpin);
        work_cell.lattice = cell.lattice;
        for i in 0..num_atom {
            work_cell.types[i] = cell.types[i];
            work_cell.position[i] = cell.position[i];
        }
        let mut mapping_table = vec![0usize; num_atom];
        let primitive = crate::spacegroup::spa_transform_to_primitive(
            &mut mapping_table,
            &work_cell,
            &dataset.transformation_matrix,
            centering,
            symprec,
        )
        .ok_or(SymError::CellStandardizationFailed)?;

        for i in 0..num_atom {
            if mapping_table[i] != dataset.mapping_to_primitive[i] as usize {
                debug::warning_print(format_args!(
                    "spglib: spa_transform_to_primitive failed ({} != {})\n",
                    mapping_table[i], dataset.mapping_to_primitive[i]
                ));
                return Err(SymError::CellStandardizationFailed);
            }
        }
        Ok(primitive)
    } else if no_idealize {
        // no_idealize, not to_primitive
        let centering = spgdb_get_spacegroup_type(dataset.hall_number).centering;
        let num_atom = cell.size;
        let mut work_cell = Cell::new(num_atom, TensorRank::NoSpin);
        work_cell.lattice = cell.lattice;
        for i in 0..num_atom {
            work_cell.types[i] = cell.types[i];
            work_cell.position[i] = cell.position[i];
        }
        let mut mapping_table = vec![0usize; num_atom];
        let primitive = crate::spacegroup::spa_transform_to_primitive(
            &mut mapping_table,
            &work_cell,
            &dataset.transformation_matrix,
            centering,
            symprec,
        )
        .ok_or(SymError::CellStandardizationFailed)?;

        for i in 0..num_atom {
            if mapping_table[i] != dataset.mapping_to_primitive[i] as usize {
                debug::warning_print(format_args!(
                    "spglib: spa_transform_to_primitive failed ({} != {})\n",
                    mapping_table[i], dataset.mapping_to_primitive[i]
                ));
                return Err(SymError::CellStandardizationFailed);
            }
        }
        if matches!(centering, Centering::Primitive) {
            return Ok(primitive);
        }
        crate::spacegroup::spa_transform_from_primitive(&primitive, centering, symprec)
            .ok_or(SymError::CellStandardizationFailed)
    } else {
        // Standard refinement
        let n_std = dataset.n_std_atoms;
        let mut cc = Cell::new(n_std, TensorRank::NoSpin);
        cc.lattice = dataset.std_lattice;
        for i in 0..n_std {
            cc.types[i] = dataset.std_types[i];
            cc.position[i] = dataset.std_positions[i];
        }
        Ok(cc)
    }
}
