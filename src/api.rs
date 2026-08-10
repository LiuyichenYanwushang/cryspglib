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
use crate::{MagneticSymmetry, SpaceGroup, SymError};

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

    pub(crate) fn to_cell(&self) -> Result<Cell, SymError> {
        if self.positions.is_empty()
            || self.types.len() != self.positions.len()
            || self
                .moments
                .as_ref()
                .is_some_and(|moments| moments.len() != self.positions.len())
        {
            return Err(SymError::InvalidInput);
        }
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
        Ok(cell)
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
    ///
    /// # Examples
    ///
    /// ```
    /// use cryspglib::Crystal;
    ///
    /// let poscar = "\
    /// Si
    /// 1.0
    ///    5.4300000000    0.0000000000    0.0000000000
    ///    0.0000000000    5.4300000000    0.0000000000
    ///    0.0000000000    0.0000000000    5.4300000000
    /// Si
    /// 2
    /// Direct
    /// 0.00 0.00 0.00
    /// 0.25 0.25 0.25
    /// ";
    /// let cry = Crystal::from_poscar(poscar).unwrap();
    /// assert_eq!(cry.natom(), 2);
    /// assert_eq!(cry.types, vec![14, 14]);
    /// ```
    pub fn from_poscar(data: &str) -> Result<Self, crate::SymError> {
        crate::parser::parse_poscar(data)
            .ok_or(crate::SymError::InvalidInput)
            .map(|parsed| Crystal {
                lattice: parsed.lattice,
                positions: parsed.positions,
                types: parsed.types,
                moments: parsed.magnetic_moments,
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
        let cell = self.crystal.to_cell()?;
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
        let cell = self.crystal.to_cell()?;
        let prim_cell = standardize_primitive_inner(&cell, self.symprec, self.angle_tolerance)?;
        Ok(cell_to_crystal(&prim_cell))
    }

    /// Standardize cell: `to_primitive` returns primitive cell; `no_idealize` skips position idealization.
    pub fn standardize(&self, to_primitive: bool, no_idealize: bool) -> Result<Crystal, SymError> {
        let cell = self.crystal.to_cell()?;
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
        let total = crate::kgrid::validate_mesh(&mesh)?;
        crate::kgrid::validate_shift(&is_shift)?;
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

        let mut grid_address = vec![[0i32; 3]; total];
        let mut mapping_table = vec![0usize; total];

        let num_ir = crate::kpoint::kpt_get_irreducible_reciprocal_mesh(
            &mut grid_address,
            &mut mapping_table,
            &mesh,
            &is_shift,
            &rot_reciprocal,
        )?;

        Ok(IrMesh {
            grid_addresses: grid_address,
            mapping_table,
            num_ir,
        })
    }

    /// Magnetic space group dataset.
    ///
    /// With magnetic moments set via [`Crystal::with_magnetic`], identifies the
    /// magnetic group. Without moments, returns the corresponding non-magnetic
    /// symmetry result.
    ///
    /// # Errors
    ///
    /// Returns [`crate::SymError::InvalidInput`] for an empty structure or
    /// inconsistent positions/types/moments lengths. Other symmetry and
    /// magnetic-space-group identification errors, including
    /// [`crate::SymError::MagneticUniMatchFailed`], are propagated unchanged.
    pub fn magnetic_dataset(&self) -> Result<MagneticSymmetry, crate::SymError> {
        crate::magnetic_symmetry_from_crystal(
            &self.crystal.lattice,
            &self.crystal.positions,
            &self.crystal.types,
            self.crystal.moments.as_deref(),
            self.symprec,
        )
    }
}

// ── Display impls ────────────────────────────────────────────────────────────

use std::fmt;

impl fmt::Display for MagneticSymmetry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "--- Space group ---")?;
        writeln!(f, "  Number:          {}", self.spacegroup_number)?;
        writeln!(f, "  International:   {}", self.international_short)?;
        writeln!(f, "  Hall number:     {}", self.hall_number)?;
        writeln!(f, "  Hall symbol:     {}", self.hall_symbol)?;

        if self.magnetic_type != crate::MagneticType::NonMagnetic {
            let type_name = match self.magnetic_type {
                crate::MagneticType::Ordinary => "Type-1 (ordinary, no time reversal)",
                crate::MagneticType::Grey => "Type-2 (grey, with pure 1')",
                crate::MagneticType::BlackWhite => "Type-3 (black-white, anti-rotation)",
                crate::MagneticType::AntiTranslation => {
                    "Type-4 (black-white, anti-translation)"
                }
                crate::MagneticType::NonMagnetic => "none",
            };
            writeln!(f, "--- Magnetic space group ---")?;
            writeln!(f, "  UNI number:      {}", self.uni_number)?;
            writeln!(
                f,
                "  Magnetic type:   {} ({type_name})",
                self.magnetic_type as i32
            )?;
            writeln!(f, "  BNS symbol:      {}", self.bns_number)?;
            writeln!(f, "  OG number:       {}", self.og_number)?;
        } else {
            writeln!(f, "  (non-magnetic)")?;
        }

        writeln!(f, "--- Symmetry operations ({}) ---", self.num_operations)?;
        for (index, ((rotation, translation), &time_reversal)) in self
            .rotations
            .iter()
            .zip(&self.translations)
            .zip(&self.time_reversals)
            .enumerate()
        {
            let prime = if time_reversal { "'" } else { " " };
            writeln!(
                f,
                "  {}. rot=[{:2},{:2},{:2};{:2},{:2},{:2};{:2},{:2},{:2}] trans=[{:.3},{:.3},{:.3}]{}",
                index + 1,
                rotation[0][0],
                rotation[0][1],
                rotation[0][2],
                rotation[1][0],
                rotation[1][1],
                rotation[1][2],
                rotation[2][0],
                rotation[2][1],
                rotation[2][2],
                translation[0],
                translation[1],
                translation[2],
                prime,
            )?;
        }
        Ok(())
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
#[derive(Debug, Clone, Default)]
pub struct SymmetryOps {
    pub operations: Vec<SymmetryOp>,
}

impl SymmetryOps {
    /// Number of symmetry operations.
    pub fn len(&self) -> usize {
        self.operations.len()
    }

    /// Whether this is an empty set.
    pub fn is_empty(&self) -> bool {
        self.operations.is_empty()
    }

    /// Iterate over symmetry operations.
    pub fn iter(&self) -> impl Iterator<Item = &SymmetryOp> {
        self.operations.iter()
    }

    /// Build from parallel arrays (structure-of-arrays form).
    ///
    /// Panics if the three slices have different lengths.
    pub fn from_parallel(
        rot: &[Mat3I],
        trans: &[Vec3],
        timerev: &[bool],
    ) -> Self {
        assert_eq!(rot.len(), trans.len(), "rot and trans length mismatch");
        assert_eq!(rot.len(), timerev.len(), "rot and timerev length mismatch");
        let n = rot.len();
        let operations: Vec<SymmetryOp> = (0..n)
            .map(|i| SymmetryOp {
                rotation: rot[i],
                translation: trans[i],
                time_reversal: timerev[i],
            })
            .collect();
        SymmetryOps { operations }
    }

    /// Build from owned parallel vectors.
    pub fn from_parallel_owned(
        rot: Vec<Mat3I>,
        trans: Vec<Vec3>,
        timerev: Vec<bool>,
    ) -> Self {
        let n = rot.len();
        assert_eq!(trans.len(), n, "length mismatch");
        assert_eq!(timerev.len(), n, "length mismatch");
        let mut operations = Vec::with_capacity(n);
        for i in 0..n {
            operations.push(SymmetryOp {
                rotation: rot[i],
                translation: trans[i],
                time_reversal: timerev[i],
            });
        }
        SymmetryOps { operations }
    }

    /// Look up symmetry operations from the space group database by Hall number.
    ///
    /// Returns all symmetry operations for the given Hall number (1–530).
    ///
    /// # Examples
    ///
    /// ```
    /// use cryspglib::SymmetryOps;
    ///
    /// // Pm-3m (Hall number 517) has 48 symmetry operations
    /// let ops = SymmetryOps::from_database(517).unwrap();
    /// assert_eq!(ops.len(), 48);
    ///
    /// // Invalid Hall number returns an error
    /// assert!(SymmetryOps::from_database(999).is_err());
    /// ```
    pub fn from_database(hall_number: usize) -> Result<Self, crate::SymError> {
        let sym = crate::spg_database::spgdb_get_spacegroup_operations(hall_number)
            .ok_or(crate::SymError::SpacegroupSearchFailed)?;
        let ops: Vec<SymmetryOp> = (0..sym.size)
            .map(|i| SymmetryOp {
                rotation: sym.rot[i],
                translation: sym.trans[i],
                time_reversal: false,
            })
            .collect();
        Ok(SymmetryOps { operations: ops })
    }

    /// Look up magnetic symmetry operations from the MSG database by UNI number.
    ///
    /// Returns operations with `time_reversal` flags set from the magnetic
    /// space group database (1–1651).
    pub fn from_magnetic_database(uni_number: usize) -> Result<Self, crate::SymError> {
        let hall = find_first_hall_for_uni(uni_number)?;
        let sym = crate::msg_database::msgdb_get_spacegroup_operations(uni_number, hall)
            .ok_or(crate::SymError::SpacegroupSearchFailed)?;
        let n = sym.size;
        let ops: Vec<SymmetryOp> = (0..n)
            .map(|i| SymmetryOp {
                rotation: sym.rot[i],
                translation: sym.trans[i],
                time_reversal: sym.timerev[i],
            })
            .collect();
        Ok(SymmetryOps { operations: ops })
    }
}

impl std::ops::Index<usize> for SymmetryOps {
    type Output = SymmetryOp;

    fn index(&self, index: usize) -> &Self::Output {
        &self.operations[index]
    }
}

impl SymmetryOps {
    /// Convenience: get symmetry operations for a space group number.
    ///
    /// Looks up the first Hall number for the SG and returns its operations.
    pub fn from_sg(sg: u8) -> Result<Self, crate::SymError> {
        let hall = find_hall_number(sg)?;
        Self::from_database(hall)
    }
}

/// Find the first Hall number whose space group number matches `sg`.
pub fn find_hall_number(sg: u8) -> Result<usize, crate::SymError> {
    for hall in 1..=530 {
        let st = crate::spg_database::spgdb_get_spacegroup_type(hall);
        if st.number == sg as usize {
            return Ok(hall);
        }
    }
    Err(crate::SymError::SpacegroupSearchFailed)
}

/// Find the first Hall number for a magnetic UNI number.
pub fn find_first_hall_for_uni(uni: usize) -> Result<usize, crate::SymError> {
    if uni == 0 || uni > 1651 {
        return Err(crate::SymError::SpacegroupSearchFailed);
    }
    for hall in 1..=530 {
        if let Some([lo, hi]) = crate::msg_database::msgdb_get_uni_candidates(hall)
            && uni >= lo
            && uni <= hi
        {
            return Ok(hall);
        }
    }
    Err(crate::SymError::SpacegroupSearchFailed)
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

// ── K-point mesh types ────────────────────────────────────────────────────────

/// Result of stabilized reciprocal mesh generation (with q-points).
#[derive(Debug, Clone)]
pub struct StabilizedMesh {
    /// Grid point addresses in fractional coordinates
    pub grid_addresses: Vec<[i32; 3]>,
    /// Full grid index → irreducible grid index mapping
    pub mapping_table: Vec<usize>,
    /// Number of irreducible points
    pub num_ir: usize,
}

/// Result of Brillouin-zone grid address relocation.
#[derive(Debug, Clone)]
pub struct BzMesh {
    /// Grid addresses relocated into first BZ
    pub grid_addresses: Vec<[i32; 3]>,
    /// Mapping table (unmapped entries = `usize::MAX`)
    pub bz_map: Vec<usize>,
    /// Number of BZ grid points
    pub num_bz: usize,
}

// ── K-point free functions ────────────────────────────────────────────────────

/// Convert a 3D grid address to a linear grid index.
///
/// # Examples
///
/// ```
/// use cryspglib::grid_point_from_address;
///
/// // Γ point is always index 0
/// let idx = grid_point_from_address([0, 0, 0], [4, 4, 4]).unwrap();
/// assert_eq!(idx, 0);
/// ```
pub fn grid_point_from_address(
    grid_address: [i32; 3],
    mesh: [i32; 3],
) -> Result<usize, SymError> {
    let mut address_double = [0i32; 3];
    let is_shift = [0i32; 3];
    crate::kgrid::kgd_get_grid_address_double_mesh(
        &mut address_double,
        &grid_address,
        &mesh,
        &is_shift,
    )?;
    crate::kgrid::kgd_get_dense_grid_point_double_mesh(&address_double, &mesh)
}

/// Generate a stabilized irreducible reciprocal mesh for given q-points.
///
/// # Examples
///
/// ```
/// use cryspglib::{Crystal, stabilized_reciprocal_mesh};
///
/// // Get symmetry operations from a cubic crystal
/// let cry = Crystal::new(
///     [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
///     vec![[0.0, 0.0, 0.0]],
///     vec![14],
/// );
/// let ds = cry.analyze().symprec(1e-5).dataset().unwrap();
///
/// // 2x2x2 mesh, no q-point distortion
/// let sm = stabilized_reciprocal_mesh(
///     [2, 2, 2],
///     [0, 0, 0],
///     true,  // time-reversal symmetry
///     &ds.rotations,
///     &[[0.0, 0.0, 0.0]],  // q-points (Γ only)
/// ).unwrap();
/// assert!(sm.num_ir > 0);
/// assert!(sm.num_ir <= 8); // 8 total points in 2x2x2 mesh
/// ```
pub fn stabilized_reciprocal_mesh(
    mesh: [i32; 3],
    is_shift: [i32; 3],
    is_time_reversal: bool,
    rotations: &[Mat3I],
    qpoints: &[[f64; 3]],
) -> Result<StabilizedMesh, SymError> {
    let total = crate::kgrid::validate_mesh(&mesh)?;
    crate::kgrid::validate_shift(&is_shift)?;
    use crate::mathfunc::MatINT;
    let mut rot = MatINT::new(rotations.len());
    for (i, r) in rotations.iter().enumerate() {
        rot.mat[i] = *r;
    }
    let mut grid_address = vec![[0i32; 3]; total];
    let mut mapping_table = vec![0usize; total];
    let num_ir = crate::kpoint::kpt_get_stabilized_reciprocal_mesh(
        &mut grid_address,
        &mut mapping_table,
        &mesh,
        &is_shift,
        if is_time_reversal { 1 } else { 0 },
        &rot,
        qpoints,
    )?;
    Ok(StabilizedMesh { grid_addresses: grid_address, mapping_table, num_ir })
}

/// Apply rotations to a grid address, returning rotated grid point indices.
///
/// # Examples
///
/// ```
/// use cryspglib::{SymmetryOps, dense_grid_points_by_rotations};
///
/// let ops = SymmetryOps::from_database(517).unwrap(); // Pm-3m
/// let rots: Vec<_> = ops.operations.iter().map(|op| op.rotation).collect();
///
/// let points = dense_grid_points_by_rotations(
///     [0, 0, 0],     // Γ point
///     &rots,
///     [4, 4, 4],
///     [0, 0, 0],
/// ).unwrap();
/// // All rotations of Γ point map to index 0
/// for p in &points {
///     assert_eq!(*p, 0);
/// }
/// ```
pub fn dense_grid_points_by_rotations(
    address_orig: [i32; 3],
    rot_reciprocal: &[Mat3I],
    mesh: [i32; 3],
    is_shift: [i32; 3],
) -> Result<Vec<usize>, SymError> {
    crate::kgrid::validate_mesh(&mesh)?;
    crate::kgrid::validate_shift(&is_shift)?;
    use crate::mathfunc::MatINT;
    let mut rot = MatINT::new(rot_reciprocal.len());
    for (i, r) in rot_reciprocal.iter().enumerate() {
        rot.mat[i] = *r;
    }
    let mut rot_grid_points = vec![0usize; rot_reciprocal.len()];
    crate::kpoint::kpt_get_dense_grid_points_by_rotations(
        &mut rot_grid_points,
        &address_orig,
        &rot,
        &mesh,
        &is_shift,
    )?;
    Ok(rot_grid_points)
}

/// Apply rotations to a grid address, returning BZ-mapped grid point indices.
///
/// # Examples
///
/// ```no_run
/// use cryspglib::{SymmetryOps, dense_bz_grid_points_by_rotations};
///
/// let ops = SymmetryOps::from_database(517).unwrap();
/// let rots: Vec<_> = ops.operations.iter().map(|op| op.rotation).collect();
/// // BZ map for 2x2x2 mesh (8 points → 64 double-mesh points)
/// let bz_map: Vec<usize> = (0..64).map(|i| i % 8).collect();
///
/// let points = dense_bz_grid_points_by_rotations(
///     [0, 0, 0],
///     &rots,
///     [2, 2, 2],
///     [0, 0, 0],
///     &bz_map,
/// ).unwrap();
/// assert_eq!(points.len(), rots.len());
/// ```
pub fn dense_bz_grid_points_by_rotations(
    address_orig: [i32; 3],
    rot_reciprocal: &[Mat3I],
    mesh: [i32; 3],
    is_shift: [i32; 3],
    bz_map: &[usize],
) -> Result<Vec<usize>, SymError> {
    crate::kgrid::validate_mesh(&mesh)?;
    crate::kgrid::validate_shift(&is_shift)?;
    use crate::mathfunc::MatINT;
    let mut rot = MatINT::new(rot_reciprocal.len());
    for (i, r) in rot_reciprocal.iter().enumerate() {
        rot.mat[i] = *r;
    }
    let mut rot_grid_points = vec![0usize; rot_reciprocal.len()];
    crate::kpoint::kpt_get_dense_bz_grid_points_by_rotations(
        &mut rot_grid_points,
        &address_orig,
        &rot,
        &mesh,
        &is_shift,
        bz_map,
    )?;
    Ok(rot_grid_points)
}

/// Relocate grid addresses into the first Brillouin zone.
///
/// # Examples
///
/// ```
/// use cryspglib::{Crystal, relocate_bz_grid_address};
///
/// // Γ-centered 2x2x2 mesh in a cubic lattice
/// let cry = Crystal::new(
///     [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
///     vec![[0.0, 0.0, 0.0]],
///     vec![14],
/// );
/// let ds = cry.analyze().symprec(1e-5).dataset().unwrap();
/// let im = cry.analyze().irreducible_mesh([2, 2, 2], [0, 0, 0], true).unwrap();
///
/// // Relocate all mesh addresses into the first BZ
/// let reciprocal_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
/// let bz = relocate_bz_grid_address(
///     &im.grid_addresses,
///     [2, 2, 2],
///     &reciprocal_lattice,
///     [0, 0, 0],
/// ).unwrap();
/// assert!(bz.num_bz > 0);
/// ```
pub fn relocate_bz_grid_address(
    grid_address: &[[i32; 3]],
    mesh: [i32; 3],
    rec_lattice: &Mat3,
    is_shift: [i32; 3],
) -> Result<BzMesh, SymError> {
    let total = crate::kgrid::validate_mesh(&mesh)?;
    crate::kgrid::validate_shift(&is_shift)?;
    let num_bz_map = total.checked_mul(8).ok_or(SymError::ArraySizeShortage)?;
    let mut bz_grid_address = vec![[0i32; 3]; num_bz_map];
    let mut bz_map = vec![0usize; num_bz_map];
    let num_bz = crate::kpoint::kpt_relocate_bz_grid_address(
        &mut bz_grid_address,
        &mut bz_map,
        grid_address,
        &mesh,
        rec_lattice,
        &is_shift,
    )?;
    Ok(BzMesh { grid_addresses: bz_grid_address, bz_map, num_bz })
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
    let container = det_determine_all(cell, hall_number, symprec, angle_tolerance)?;

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

    for (i, &mapped) in mapping_table.iter().take(primitive.size).enumerate() {
        if mapped != i {
            debug::warning_print(format_args!(
                "spglib: spa_transform_to_primitive failed ({} != {})\n",
                mapped, i
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

        for (&mapped, &expected) in mapping_table
            .iter()
            .zip(&dataset.mapping_to_primitive)
        {
            if mapped != expected as usize {
                debug::warning_print(format_args!(
                    "spglib: spa_transform_to_primitive failed ({} != {})\n",
                    mapped, expected
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

        for (&mapped, &expected) in mapping_table
            .iter()
            .zip(&dataset.mapping_to_primitive)
        {
            if mapped != expected as usize {
                debug::warning_print(format_args!(
                    "spglib: spa_transform_to_primitive failed ({} != {})\n",
                    mapped, expected
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

#[cfg(test)]
mod ordinary_input_contract_tests {
    use super::*;

    fn cubic_lattice() -> Mat3 {
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }

    #[test]
    fn ordinary_analysis_rejects_empty_crystal_across_terminal_methods() {
        let empty = Crystal::new(cubic_lattice(), vec![], vec![]);
        let analysis = empty.analyze();

        assert!(matches!(analysis.dataset(), Err(SymError::InvalidInput)));
        assert!(matches!(analysis.symmetry(), Err(SymError::InvalidInput)));
        assert!(matches!(analysis.primitive_cell(), Err(SymError::InvalidInput)));
        assert!(matches!(
            analysis.standardize(true, false),
            Err(SymError::InvalidInput)
        ));
        assert!(matches!(analysis.hall_number(), Err(SymError::InvalidInput)));
        assert!(matches!(analysis.international(), Err(SymError::InvalidInput)));
        assert!(matches!(
            analysis.irreducible_mesh([2, 2, 2], [0, 0, 0], true),
            Err(SymError::InvalidInput)
        ));

        let empty_layer = Crystal::new(cubic_lattice(), vec![], vec![])
            .with_layer(AperiodicAxis::Z);
        assert!(matches!(
            empty_layer.analyze().dataset(),
            Err(SymError::InvalidInput)
        ));
    }

    #[test]
    fn ordinary_analysis_rejects_mutated_parallel_fields() {
        let mut empty_positions = Crystal::new(
            cubic_lattice(),
            vec![[0.0, 0.0, 0.0]],
            vec![14],
        );
        empty_positions.positions.clear();
        assert!(matches!(
            empty_positions.analyze().dataset(),
            Err(SymError::InvalidInput)
        ));

        let mut missing_type = Crystal::new(
            cubic_lattice(),
            vec![[0.0, 0.0, 0.0]],
            vec![14],
        );
        missing_type.types.clear();
        assert!(matches!(
            missing_type.analyze().primitive_cell(),
            Err(SymError::InvalidInput)
        ));

        let mut short_moments = Crystal::new(
            cubic_lattice(),
            vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
            vec![14, 14],
        );
        short_moments.moments = Some(vec![[0.0; 3]]);
        assert!(matches!(
            short_moments.analyze().standardize(false, false),
            Err(SymError::InvalidInput)
        ));
    }

    #[test]
    fn ordinary_analysis_valid_control_still_succeeds() {
        let crystal = Crystal::new(
            cubic_lattice(),
            vec![[0.0, 0.0, 0.0]],
            vec![14],
        );
        let dataset = crystal.analyze().dataset().unwrap();
        assert_eq!(dataset.spacegroup_number, 221);
    }

    #[test]
    fn mesh_apis_reject_invalid_meshes_and_shifts() {
        let identity = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]];
        for mesh in [[0, 2, 2], [-1, 2, 2], [i32::MIN, 1, 1]] {
            assert!(matches!(
                grid_point_from_address([0, 0, 0], mesh),
                Err(SymError::InvalidInput)
            ));
            assert!(matches!(
                stabilized_reciprocal_mesh(mesh, [0, 0, 0], true, &identity, &[]),
                Err(SymError::InvalidInput)
            ));
        }

        assert!(matches!(
            stabilized_reciprocal_mesh([2, 2, 2], [0, 0, 2], true, &identity, &[]),
            Err(SymError::InvalidInput)
        ));
        assert!(matches!(
            dense_grid_points_by_rotations(
                [0, 0, 0],
                &identity,
                [2, 2, 2],
                [0, -1, 0],
            ),
            Err(SymError::InvalidInput)
        ));
    }

    #[test]
    fn allocating_mesh_apis_reject_unsafe_sizes_and_short_buffers() {
        let identity = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]];
        assert!(matches!(
            stabilized_reciprocal_mesh(
                [1291, 1291, 1291],
                [0, 0, 0],
                true,
                &identity,
                &[],
            ),
            Err(SymError::ArraySizeShortage)
        ));
        assert!(matches!(
            dense_bz_grid_points_by_rotations(
                [0, 0, 0],
                &identity,
                [2, 2, 2],
                [0, 0, 0],
                &[],
            ),
            Err(SymError::ArraySizeShortage)
        ));
        assert!(matches!(
            relocate_bz_grid_address(
                &[],
                [2, 2, 2],
                &cubic_lattice(),
                [0, 0, 0],
            ),
            Err(SymError::ArraySizeShortage)
        ));
    }

    #[test]
    fn mesh_valid_controls_keep_existing_results() {
        let identity = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]];
        assert_eq!(
            grid_point_from_address([1, 2, 3], [4, 4, 4]).unwrap(),
            57
        );
        assert_eq!(
            grid_point_from_address([-1, 0, 0], [4, 4, 4]).unwrap(),
            3
        );
        assert_eq!(
            grid_point_from_address([i32::MIN, 0, 0], [4, 4, 4]).unwrap(),
            0
        );

        let points = dense_grid_points_by_rotations(
            [0, 0, 0],
            &identity,
            [4, 4, 4],
            [0, 0, 0],
        )
        .unwrap();
        assert_eq!(points, vec![0]);

        let stabilized = stabilized_reciprocal_mesh(
            [2, 2, 2],
            [0, 0, 0],
            true,
            &identity,
            &[[0.0, 0.0, 0.0]],
        )
        .unwrap();
        assert_eq!(stabilized.num_ir, 8);
        assert_eq!(stabilized.grid_addresses.len(), 8);
    }
}
