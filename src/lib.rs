//! cryspglib — Rust 晶体空间群识别库。
//!
//! 基于 [spglib](https://github.com/spglib/spglib) 的纯 Rust 移植，
//! 提供晶体对称性分析、空间群识别、标准晶胞构造和 k 点网格生成。
//!
//! # 快速开始
//!
//! ## 非磁空间群
//!
//! ```no_run
//! use cryspglib::Crystal;
//!
//! // FCC Al (space group Fm-3m, #225)
//! let al = Crystal::new(
//!     [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
//!     vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]],
//!     vec![13, 13, 13, 13],
//! );
//! let ds = al.analyze().symprec(1e-5).dataset()?;
//! println!("Space group #{}: {}", ds.spacegroup_number, ds.international_symbol);
//! # Ok::<(), cryspglib::SymError>(())
//! ```
//!
//! ## 磁性空间群
//!
//! ```no_run
//! use cryspglib::Crystal;
//!
//! // BCC AFM [111]: Fe at [0,0,0] and [0.5,0.5,0.5], opposite spins
//! let n = (3.0_f64).sqrt();
//! let fe = Crystal::new(
//!     [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
//!     vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
//!     vec![26, 26],
//! ).with_magnetic(vec![
//!     [1.0/n, 1.0/n, 1.0/n],
//!     [-1.0/n, -1.0/n, -1.0/n],
//! ]);
//!
//! let r = fe.analyze().symprec(1e-5).magnetic_dataset().unwrap();
//! println!("SG #{} → UNI={}, BNS={}", r.spacegroup_number, r.uni_number, r.bns_number);
//! ```
//!
//! # 主要类型
//!
//! | 类型 | 说明 |
//! |------|------|
//! | [`Crystal`] | 晶体结构（晶格 + 原子位置 + 可选磁矩） |
//! | [`SymmetryAnalysis`] | 对称性分析构建器（`.symprec()`, `.dataset()`, `.magnetic_dataset()` 等） |
//! | [`SpaceGroup`] | 空间群完整数据（编号、符号、Wyckoff 位置、对称操作等） |
//! | [`SymmetryOp`] | 单个对称操作 `{R\|t}` + `time_reversal` |
//! | [`SymmetryOps`] | 对称操作集合，支持 [`SymmetryOps::from_database`] |
//! | [`MagneticSymmetry`] | 磁空间群分析结果，实现 `Display` trait |
//! | [`MagneticSpaceGroupType`] | 磁空间群类型，支持返回 `Result` 的 `.from_uni()` 和 `.classify()` |
//! | [`SpaceGroupType`] | 空间群类型信息，支持 `.from_hall()` |
//! | [`IrMesh`] | 不可约 k 点网格 |
//! | [`StabilizedMesh`] | 稳定化倒易网格（含 q 点） |
//! | [`BzMesh`] | 第一布里渊区重定位结果 |
//!
//! # 晶格矩阵约定
//!
//! 所有 3x3 矩阵采用 `lattice[cart][vec]` 布局（行=笛卡尔分量，列=晶格矢量）。
//! 详见 [`mathfunc`] 模块文档。
//!
//! The public surface is Rust-native: owned outputs, typed errors, methods,
//! and domain types replace C-style output parameters and sentinel values.

pub mod api;
pub mod arithmetic;
pub mod cell;
pub mod debug;
pub mod delaunay;
pub mod determination;
pub mod hall_symbol;
pub mod irrep;
pub mod kgrid;
pub mod kpoint;
pub mod magnetic_spacegroup;
pub mod mathfunc;
pub mod msg_database;
pub mod niggli;
pub mod operation_group;
pub mod overlap;
pub mod parser;
pub mod pointgroup;
pub mod primitive;
pub mod refinement;
pub mod site_symmetry;
pub mod sitesym_database;
pub mod spacegroup;
pub mod spg_database;
pub mod spin;
pub mod symmetry;

use crate::mathfunc::{Mat3, Mat3I, Vec3, mat_inverse_matrix_d3, mat_multiply_matrix_d3};
use crate::pointgroup::ptg_get_pointgroup;
use crate::primitive::prm_get_primitive_symmetry;
use crate::spacegroup::spa_search_spacegroup_with_symmetry;
use crate::spg_database::spgdb_get_spacegroup_type;
use crate::symmetry::Symmetry;

// Re-export the new Rust-idiomatic API
pub use api::{
    BzMesh, Crystal, ExternalFields, IrMesh, StabilizedMesh, SymmetryAnalysis, SymmetryOp,
    SymmetryOps, dense_bz_grid_points_by_rotations, dense_grid_points_by_rotations,
    grid_point_from_address, relocate_bz_grid_address, stabilized_reciprocal_mesh,
};
pub use operation_group::{
    MagneticGroupIdentification, MagneticOperationSetError, SpinLiftError,
    ValidatedMagneticOperationSet, axial_spin_half_lift,
};
pub use pointgroup::pointgroup_from_rotations;

// ---------------------------------------------------------------------------
// Version constants
// ---------------------------------------------------------------------------
/// Library version.
pub const VERSION: &str = "0.2.0";

// ---------------------------------------------------------------------------
// Error codes
// ---------------------------------------------------------------------------
/// Symmetry analysis error codes.
#[derive(thiserror::Error, Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum SymError {
    /// 无错误
    #[error("no error")]
    Success = 0,
    /// 空间群搜索失败
    #[error("spacegroup search failed")]
    SpacegroupSearchFailed = 1,
    /// 晶胞标准化失败
    #[error("cell standardization failed")]
    CellStandardizationFailed = 2,
    /// 对称操作搜索失败
    #[error("symmetry operation search failed")]
    SymmetryOperationSearchFailed = 3,
    /// 原子间距过近
    #[error("too close distance between atoms")]
    AtomsTooClose = 4,
    /// 点群未找到
    #[error("pointgroup not found")]
    PointgroupNotFound = 5,
    /// Niggli 约化失败
    #[error("Niggli reduction failed")]
    NiggliFailed = 6,
    /// Delaunay 约化失败
    #[error("Delaunay reduction failed")]
    DelaunayFailed = 7,
    /// 数组大小不足
    #[error("array size shortage")]
    ArraySizeShortage = 8,
    /// 输入格式无效
    #[error("invalid input format")]
    InvalidInput = 9,
    /// 数学运算失败
    #[error("math operation failed")]
    MathFailed = 10,
    /// 磁操作生成失败
    #[error("magnetic operation generation failed")]
    MagneticOpGenerationFailed = 11,
    /// 磁参考空间群搜索失败
    #[error("magnetic reference space group search failed")]
    MagneticReferenceGroupFailed = 12,
    /// 磁群 fallback 参考设置失败
    #[error("magnetic fallback reference setting failed")]
    MagneticFallbackReferenceFailed = 13,
    /// Hall 编号无 UNI 候选
    #[error("no UNI candidates for Hall number")]
    MagneticUniCandidatesNotFound = 14,
    /// UNI 候选匹配全部失败
    #[error("all UNI candidates failed full magnetic Seitz-set matching")]
    MagneticUniMatchFailed = 15,
    /// 磁原胞晶格确定失败
    #[error("magnetic primitive lattice determination failed")]
    MagneticPrimitiveLatticeFailed = 16,
    /// 完整磁对称操作对应多个不同的 UNI；需要母群 Hall 信息消歧
    #[error("magnetic symmetry has multiple UNI candidates; parent Hall number is required")]
    MagneticUniAmbiguous = 17,
}

// ---------------------------------------------------------------------------
// Public data structures
// ---------------------------------------------------------------------------

/// 空间群数据集的完整结构。
///
/// 包含标准晶胞、对称操作、Wyckoff 位置标记和映射信息。
/// 所有动态数据由 Rust 的 [`Vec`] 所有权管理，无需手动释放。
#[derive(Debug, Clone)]
pub struct SpaceGroup {
    /// 空间群编号 (1–230)
    pub spacegroup_number: usize,
    /// Hall 编号 (1–530)
    pub hall_number: usize,
    /// 国际符号 (最多 11 字符)
    pub international_symbol: String,
    /// Hall 符号 (最多 17 字符)
    pub hall_symbol: String,
    /// 选择 (最多 6 字符)
    pub choice: String,
    /// 变换矩阵 (Bravais → 原始晶胞)
    pub transformation_matrix: Mat3,
    /// 原点平移
    pub origin_shift: Vec3,
    /// 对称操作数量
    pub n_operations: usize,
    /// 旋转矩阵 [n_operations][3][3]
    pub rotations: Vec<Mat3I>,
    /// 平移矢量 [n_operations][3]
    pub translations: Vec<Vec3>,
    /// 原子数
    pub n_atoms: usize,
    /// Wyckoff 字母编码 (0=a, 1=b, ..., 26=z)
    pub wyckoffs: Vec<i32>,
    /// 位点对称性符号
    pub site_symmetry_symbols: Vec<String>,
    /// 对等原子映射
    pub equivalent_atoms: Vec<i32>,
    /// 晶体学轨道
    pub crystallographic_orbits: Vec<i32>,
    /// 原子 → 原胞映射
    pub mapping_to_primitive: Vec<i32>,
    /// 标准晶胞原子数
    pub n_std_atoms: usize,
    /// 标准晶胞晶格
    pub std_lattice: Mat3,
    /// 标准晶胞原子位置
    pub std_positions: Vec<Vec3>,
    /// 标准晶胞原子类型
    pub std_types: Vec<i32>,
    /// 标准晶胞旋转矩阵
    pub std_rotation_matrix: Mat3,
    /// 标准晶胞 → 原胞映射
    pub std_mapping_to_primitive: Vec<i32>,
    /// 原胞晶格
    pub primitive_lattice: Mat3,
    /// 点群符号 (最多 6 字符)
    pub pointgroup_symbol: String,
}

/// 空间群类型信息（从数据库查询）。
#[derive(Debug, Clone)]
pub struct SpaceGroupType {
    /// 空间群编号 (1–230)
    pub number: usize,
    /// Hall 编号
    pub hall_number: usize,
    /// Schoenflies 符号
    pub schoenflies: String,
    /// Hall 符号
    pub hall_symbol: String,
    /// 选择
    pub choice: String,
    /// 国际符号（完整）
    pub international: String,
    /// 国际符号（完整，多行格式）
    pub international_full: String,
    /// 国际符号（短格式）
    pub international_short: String,
    /// 点群国际符号
    pub pointgroup_international: String,
    /// 点群 Schoenflies 符号
    pub pointgroup_schoenflies: String,
    /// 算术晶体类编号
    pub arithmetic_crystal_class_number: i32,
    /// 算术晶体类符号
    pub arithmetic_crystal_class_symbol: String,
}

impl SpaceGroupType {
    /// Look up space group type by Hall number (1–530).
    ///
    /// # Examples
    ///
    /// ```
    /// use cryspglib::SpaceGroupType;
    ///
    /// // Pm-3m (Hall number 517)
    /// let sg = SpaceGroupType::from_hall(517).unwrap();
    /// assert_eq!(sg.number, 221);
    /// assert_eq!(sg.international_short.trim(), "Pm-3m");
    /// assert_eq!(sg.schoenflies.trim(), "Oh^1");
    ///
    /// // Invalid Hall number
    /// assert!(SpaceGroupType::from_hall(999).is_err());
    /// ```
    pub fn from_hall(hall_number: usize) -> Result<Self, SymError> {
        if hall_number > 0 && hall_number < 531 {
            get_spacegroup_type(hall_number)
        } else {
            Err(SymError::SpacegroupSearchFailed)
        }
    }
}

/// 磁性空间群类型。
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MagneticType {
    /// 非磁 (UNI=0)
    NonMagnetic = 0,
    /// Type-1: 普通磁结构，无时间反演对称性
    Ordinary = 1,
    /// Type-2: 灰色磁结构（顺磁），含纯时间反演操作
    Grey = 2,
    /// Type-3: 黑白磁结构，反转动
    BlackWhite = 3,
    /// Type-4: 黑白磁结构，反平移
    AntiTranslation = 4,
}

/// 磁性空间群类型（从数据库查询）。
#[derive(Debug, Clone)]
pub struct MagneticSpaceGroupType {
    /// UNI 编号
    pub uni_number: usize,
    /// Litvin 编号
    pub litvin_number: usize,
    /// BNS 符号 (最多 8 字符)
    pub bns_number: String,
    /// OG 符号 (最多 12 字符)
    pub og_number: String,
    /// 晶体学编号 (1–230)
    pub number: usize,
    /// 磁性类型 (1-4)
    pub type_: MagneticType,
}

impl MagneticSpaceGroupType {
    /// Look up a magnetic space group type by UNI number (1–1651).
    ///
    /// # Errors
    ///
    /// Returns [`SymError::InvalidInput`] when `uni_number` is outside
    /// `1..=1651`.
    ///
    /// # Examples
    ///
    /// ```
    /// use cryspglib::{MagneticSpaceGroupType, MagneticType};
    ///
    /// // UNI 1331 = BNS 166.101 (Type-3 black-white, parent R-3m)
    /// let msg = MagneticSpaceGroupType::from_uni(1331).unwrap();
    /// assert_eq!(msg.uni_number, 1331);
    /// assert_eq!(msg.bns_number.trim(), "166.101");
    /// assert_eq!(msg.type_, MagneticType::BlackWhite);
    ///
    /// // UNI 1005 = BNS 123.345 (Type-3, ferromagnetic along [001])
    /// let msg = MagneticSpaceGroupType::from_uni(1005).unwrap();
    /// assert_eq!(msg.bns_number.trim(), "123.345");
    ///
    /// assert!(MagneticSpaceGroupType::from_uni(0).is_err());
    /// assert!(MagneticSpaceGroupType::from_uni(1652).is_err());
    /// ```
    pub fn from_uni(uni_number: usize) -> Result<Self, SymError> {
        if !(1..=1651).contains(&uni_number) {
            return Err(SymError::InvalidInput);
        }
        let msgtype = crate::msg_database::msgdb_get_magnetic_spacegroup_type(uni_number);
        Ok(MagneticSpaceGroupType {
            uni_number: msgtype.uni_number,
            litvin_number: msgtype.litvin_number,
            bns_number: msgtype.bns_number.to_string(),
            og_number: msgtype.og_number.to_string(),
            number: msgtype.number,
            type_: msgtype.type_,
        })
    }

    /// Classify magnetic space group type from a set of symmetry operations.
    ///
    /// `time_reversals` can be `None` (treated as all-false / ordinary operations).
    ///
    /// Returns the identification error, including
    /// [`SymError::MagneticUniAmbiguous`], instead of disguising failure as a
    /// default `UNI=0` non-magnetic result.
    ///
    /// Use [`crate::magnetic_spacegroup::msg_identify_with_parent_hall`] when
    /// the non-magnetic parent Hall number is known.
    ///
    /// # Examples
    ///
    /// ```
    /// use cryspglib::{MagneticSpaceGroupType, SymmetryOps, MagneticType};
    ///
    /// // Get Pm-3m symmetry operations from database
    /// let ops = SymmetryOps::from_database(517).unwrap();
    /// let rots: Vec<_> = ops.operations.iter().map(|op| op.rotation).collect();
    /// let trans: Vec<_> = ops.operations.iter().map(|op| op.translation).collect();
    /// let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    ///
    /// // Without time reversal → Type-1 (ordinary)
    /// let msg = MagneticSpaceGroupType::classify(
    ///     &rots, &trans, None, &lattice, 1e-5,
    /// ).unwrap();
    /// assert_eq!(msg.type_, MagneticType::Ordinary);
    /// assert!(msg.uni_number > 0);
    /// ```
    pub fn classify(
        rotations: &[Mat3I],
        translations: &[Vec3],
        time_reversals: Option<&[bool]>,
        lattice: &Mat3,
        symprec: f64,
    ) -> Result<Self, SymError> {
        let n_ops = rotations.len();
        if n_ops == 0
            || translations.len() != n_ops
            || time_reversals.is_some_and(|values| values.len() != n_ops)
        {
            return Err(SymError::InvalidInput);
        }

        let mut magnetic_symmetry = crate::symmetry::MagneticSymmetry::new(n_ops);
        for i in 0..n_ops {
            magnetic_symmetry.rot[i] = rotations[i];
            magnetic_symmetry.trans[i] = translations[i];
            magnetic_symmetry.timerev[i] = time_reversals.is_some_and(|values| values[i]);
        }

        let dataset = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
            lattice,
            &magnetic_symmetry,
            symprec,
        )?;
        Self::from_uni(dataset.uni_number)
    }
}

/// 从对称操作确定 Hall 编号。
///
/// 给定一组旋转和平移操作，搜索匹配的空间群 Hall 编号。
///
/// # Errors
///
/// `rotations` 和 `translations` 必须非空且长度相等；否则返回
/// [`SymError::InvalidInput`]。操作集合无法匹配到空间群时返回
/// [`SymError::SpacegroupSearchFailed`]。
pub(crate) fn hall_number_from_symmetry(
    rotations: &[Mat3I],
    translations: &[Vec3],
    symprec: f64,
) -> Result<usize, SymError> {
    let lattice: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let hall_number = identify_hall_number(rotations, translations, &lattice, false, symprec)?;
    if hall_number > 0 {
        Ok(hall_number)
    } else {
        Err(SymError::SpacegroupSearchFailed)
    }
}

/// 磁空间群 + 对称操作的完整分析结果。
pub struct MagneticSymmetry {
    /// 磁矩约束之前检测到的非磁结构母空间群编号 (1-230)。
    ///
    /// 这不一定等于磁群的 family-space-group 编号；后者应从
    /// [`ValidatedMagneticOperationSet::identify`] 的
    /// `MagneticGroupIdentification::spacegroup_number` 读取。
    pub spacegroup_number: usize,
    /// 国际符号（短）
    pub international_short: String,
    /// 非磁结构母空间群的 Hall setting (1-530)，用于记录输入坐标系。
    pub hall_number: usize,
    /// Hall 符号
    pub hall_symbol: String,
    /// 磁空间群 UNI 编号 (0 表示未找到)
    pub uni_number: usize,
    /// 磁性类型: 0=非磁, 1=ordinary, 2=grey, 3=black-white, 4=anti-translation
    pub magnetic_type: MagneticType,
    /// BNS 符号（如 "221.93"）
    pub bns_number: String,
    /// OG 符号（如 "221.2.1595"）
    pub og_number: String,
    /// 对称操作数
    pub num_operations: usize,
    /// 旋转矩阵 (整数 3x3)
    pub rotations: Vec<Mat3I>,
    /// 平移向量 (分数坐标)
    pub translations: Vec<Vec3>,
    /// 时间反演标记 (false=ordinary, true=anti)
    pub time_reversals: Vec<bool>,
}

/// Analyze magnetic symmetry from a validated crystal representation.
///
/// `magnetic_moments` 为 `None` 时不考虑磁性，仅返回非磁空间群。
/// 每个原子的磁矩为 3 分量 `[mx, my, mz]`。
///
/// 返回包含非磁空间群、磁空间群、对称操作的结构。
///
/// # Errors
///
/// Returns [`SymError::InvalidInput`] for an empty structure or inconsistent
/// positions/types/moments lengths. Other symmetry and magnetic-space-group
/// identification errors, including [`SymError::MagneticUniMatchFailed`], are
/// propagated unchanged.
pub(crate) fn magnetic_symmetry_from_crystal(
    lattice: &Mat3,
    positions: &[Vec3],
    types: &[i32],
    magnetic_moments: Option<&[[f64; 3]]>,
    symprec: f64,
) -> Result<MagneticSymmetry, SymError> {
    let n_atoms = positions.len();
    if n_atoms == 0 || types.len() != n_atoms {
        return Err(SymError::InvalidInput);
    }
    if magnetic_moments.is_some_and(|moments| moments.len() != n_atoms) {
        return Err(SymError::InvalidInput);
    }

    // --- 构建 Cell ---
    let has_mag = magnetic_moments.is_some();
    let tensor_rank = if has_mag {
        crate::cell::TensorRank::NonCollinear
    } else {
        crate::cell::TensorRank::NoSpin
    };

    let mut cell = crate::cell::Cell::new(n_atoms, tensor_rank);
    cell.set_cell(lattice, positions, types);

    if has_mag {
        let moments = magnetic_moments.unwrap();
        for (i, moment) in moments.iter().enumerate().take(n_atoms) {
            cell.tensors[i * 3] = moment[0];
            cell.tensors[i * 3 + 1] = moment[1];
            cell.tensors[i * 3 + 2] = moment[2];
        }
    }
    cell.aperiodic_axis = None;

    // --- 1. 非磁空间群 ---
    let primitive = crate::primitive::prm_get_primitive(&cell, symprec, -1.0)?;
    let spg = crate::spacegroup::spa_search_spacegroup(&primitive, 0, symprec, -1.0)?;
    let hall_number = spg.hall_number;

    // --- 2. 非磁对称操作 (用常规晶胞获取, 保证基矢正确) ---
    let nonspin_sym = crate::symmetry::sym_get_operation(&cell, symprec, -1.0)?;

    if !has_mag {
        // 无磁矩: 只返回非磁结果
        let rot = (0..nonspin_sym.len()).map(|i| nonspin_sym.rot[i]).collect();
        let trans = (0..nonspin_sym.len())
            .map(|i| nonspin_sym.trans[i])
            .collect();
        let timerev = vec![false; nonspin_sym.len()];
        let spg_type = crate::spg_database::spgdb_get_spacegroup_type(hall_number);
        return Ok(MagneticSymmetry {
            spacegroup_number: spg.number,
            international_short: spg.international_short.trim().to_string(),
            hall_number,
            hall_symbol: spg_type.hall_symbol.trim().to_string(),
            uni_number: 0,
            magnetic_type: MagneticType::NonMagnetic,
            bns_number: String::new(),
            og_number: String::new(),
            num_operations: nonspin_sym.len(),
            rotations: rot,
            translations: trans,
            time_reversals: timerev,
        });
    }

    // --- 3. 磁对称操作 (从磁矩计算 timerev 标记) ---
    let mag_sym =
        crate::spin::operations_with_site_tensors(crate::spin::MagneticOperationSearch {
            symmetry: &nonspin_sym,
            cell: &cell,
            with_time_reversal: true,
            is_axial: true,
            symprec,
            angle_tolerance: -1.0,
            magnetic_symprec: -1.0,
        })?;

    // 如果有磁矩但磁对称操作数为 0，尝试用简单方法
    // (operations_with_site_tensors 可能因原胞匹配失败)
    let (final_mag_sym, _used_fallback) = if mag_sym.is_empty() {
        // fallback: 手动计算 timerev
        let crystal_ops: Vec<(Mat3I, Vec3)> = (0..nonspin_sym.len())
            .map(|i| (nonspin_sym.rot[i], nonspin_sym.trans[i]))
            .collect();
        let moments = magnetic_moments.unwrap();
        let tr = manual_compute_timerev(positions, moments, &crystal_ops, symprec);
        let valid: Vec<usize> = tr
            .iter()
            .cloned()
            .enumerate()
            .filter(|(_, t)| *t != -1)
            .map(|(i, _)| i)
            .collect();
        let n = valid.len();
        if n == 0 {
            return Err(SymError::MagneticOpGenerationFailed);
        }
        let mut fallback = crate::symmetry::MagneticSymmetry::new(n);
        for (j, &idx) in valid.iter().enumerate() {
            fallback.rot[j] = nonspin_sym.rot[idx];
            fallback.trans[j] = nonspin_sym.trans[idx];
            fallback.timerev[j] = tr[idx] != 0;
        }
        (fallback, true)
    } else {
        (mag_sym, false)
    };

    // --- 4. 磁空间群识别 ---
    // 用已求得的非磁 Hall 编号作为 parent_hall_number fallback。
    // 当 FSG 空间群搜索失败时（如多原子磁细胞的原胞约化限制），
    // fallback 直接使用非磁母空间群来搜索 UNI 候选。
    let identification = crate::magnetic_spacegroup::msg_identify_with_parent_hall(
        lattice,
        &final_mag_sym,
        Some(hall_number),
        symprec,
    );
    let (uni_number, magnetic_type, bns_number, og_number) =
        magnetic_identification_metadata(identification)?;

    let spg_type = crate::spg_database::spgdb_get_spacegroup_type(hall_number);
    let rot_out = (0..final_mag_sym.len())
        .map(|i| final_mag_sym.rot[i])
        .collect();
    let trans_out = (0..final_mag_sym.len())
        .map(|i| final_mag_sym.trans[i])
        .collect();
    let tr_out = (0..final_mag_sym.len())
        .map(|i| final_mag_sym.timerev[i])
        .collect();

    Ok(MagneticSymmetry {
        spacegroup_number: spg.number,
        international_short: spg.international_short.trim().to_string(),
        hall_number,
        hall_symbol: spg_type.hall_symbol.trim().to_string(),
        uni_number,
        magnetic_type,
        bns_number,
        og_number,
        num_operations: final_mag_sym.len(),
        rotations: rot_out,
        translations: trans_out,
        time_reversals: tr_out,
    })
}

type MagneticIdentificationMetadata = (usize, MagneticType, String, String);

fn magnetic_identification_metadata(
    identification: Result<crate::magnetic_spacegroup::MagneticDataset, SymError>,
) -> Result<MagneticIdentificationMetadata, SymError> {
    let dataset = identification?;
    let msg_type = crate::msg_database::msgdb_get_magnetic_spacegroup_type(dataset.uni_number);
    Ok((
        dataset.uni_number,
        msg_type.type_,
        msg_type.bns_number.to_string(),
        msg_type.og_number.to_string(),
    ))
}

#[cfg(test)]
mod magnetic_dataset_contract_tests {
    use super::*;

    fn cubic_lattice() -> Mat3 {
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }

    #[test]
    fn magnetic_dataset_propagates_uni_match_failed() {
        let result = magnetic_identification_metadata(Err(SymError::MagneticUniMatchFailed));

        assert_eq!(result.unwrap_err(), SymError::MagneticUniMatchFailed);
    }

    #[test]
    fn magnetic_dataset_rejects_empty_structure() {
        let lattice = cubic_lattice();

        assert!(matches!(
            magnetic_symmetry_from_crystal(&lattice, &[], &[], None, 1e-5),
            Err(SymError::InvalidInput)
        ));
        assert!(matches!(
            magnetic_symmetry_from_crystal(&lattice, &[], &[], Some(&[]), 1e-5),
            Err(SymError::InvalidInput)
        ));
    }

    #[test]
    fn magnetic_dataset_rejects_parallel_array_length_mismatches() {
        let lattice = cubic_lattice();
        let positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];
        let types = [26, 26];
        let one_moment = [[1.0, 0.0, 0.0]];
        let three_moments = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];

        for bad_types in [&[26][..], &[26, 26, 26][..]] {
            assert!(matches!(
                magnetic_symmetry_from_crystal(&lattice, &positions, bad_types, None, 1e-5),
                Err(SymError::InvalidInput)
            ));
        }
        for bad_moments in [&[][..], &one_moment[..], &three_moments[..]] {
            assert!(matches!(
                magnetic_symmetry_from_crystal(
                    &lattice,
                    &positions,
                    &types,
                    Some(bad_moments),
                    1e-5,
                ),
                Err(SymError::InvalidInput)
            ));
        }
    }
}

#[cfg(test)]
mod hall_number_from_symmetry_contract_tests {
    use super::*;

    const IDENTITY: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

    #[test]
    fn hall_number_rejects_empty_and_mismatched_parallel_slices() {
        let one_rotation = [IDENTITY];
        let two_rotations = [IDENTITY, IDENTITY];
        let one_translation = [[0.0; 3]];
        let two_translations = [[0.0; 3], [0.5, 0.0, 0.0]];

        for (rotations, translations) in [
            (&[][..], &[][..]),
            (&[][..], &one_translation[..]),
            (&one_rotation[..], &[][..]),
            (&one_rotation[..], &two_translations[..]),
            (&two_rotations[..], &one_translation[..]),
        ] {
            assert!(matches!(
                hall_number_from_symmetry(rotations, translations, 1e-5),
                Err(SymError::InvalidInput)
            ));
        }
    }

    #[test]
    fn hall_number_identifies_valid_database_operations() {
        let operations = SymmetryOps::from_database(517).unwrap();
        let rotations: Vec<_> = operations
            .operations
            .iter()
            .map(|operation| operation.rotation)
            .collect();
        let translations: Vec<_> = operations
            .operations
            .iter()
            .map(|operation| operation.translation)
            .collect();

        let hall = hall_number_from_symmetry(&rotations, &translations, 1e-5).unwrap();
        assert_eq!(hall, 517);
        assert_eq!(SpaceGroupType::from_hall(hall).unwrap().number, 221);
    }
}

/// 手动计算磁矩变换的 timerev 标记 (fallback)。
fn manual_compute_timerev(
    positions: &[Vec3],
    moments: &[[f64; 3]],
    ops: &[(Mat3I, Vec3)],
    symprec: f64,
) -> Vec<i32> {
    use crate::mathfunc::mat_get_determinant_i3;
    let snap = |x: f64| (x * 2.0).round() / 2.0;
    let snapped_pos: Vec<_> = positions
        .iter()
        .map(|p| [snap(p[0]), snap(p[1]), snap(p[2])])
        .collect();

    ops.iter()
        .map(|(rot, trans)| {
            let det = mat_get_determinant_i3(rot);
            let mut global_tr: Option<i32> = None;

            for i in 0..positions.len() {
                let p_new = [
                    snap(
                        (rot[0][0] as f64 * positions[i][0]
                            + rot[0][1] as f64 * positions[i][1]
                            + rot[0][2] as f64 * positions[i][2]
                            + trans[0])
                            .rem_euclid(1.0),
                    ),
                    snap(
                        (rot[1][0] as f64 * positions[i][0]
                            + rot[1][1] as f64 * positions[i][1]
                            + rot[1][2] as f64 * positions[i][2]
                            + trans[1])
                            .rem_euclid(1.0),
                    ),
                    snap(
                        (rot[2][0] as f64 * positions[i][0]
                            + rot[2][1] as f64 * positions[i][1]
                            + rot[2][2] as f64 * positions[i][2]
                            + trans[2])
                            .rem_euclid(1.0),
                    ),
                ];

                let j = snapped_pos.iter().position(|sp| {
                    (sp[0] - p_new[0]).abs() < 0.01
                        && (sp[1] - p_new[1]).abs() < 0.01
                        && (sp[2] - p_new[2]).abs() < 0.01
                });
                let j = match j {
                    Some(j) => j,
                    None => return -1,
                };

                let m_new = [
                    (det as f64)
                        * (rot[0][0] as f64 * moments[i][0]
                            + rot[0][1] as f64 * moments[i][1]
                            + rot[0][2] as f64 * moments[i][2]),
                    (det as f64)
                        * (rot[1][0] as f64 * moments[i][0]
                            + rot[1][1] as f64 * moments[i][1]
                            + rot[1][2] as f64 * moments[i][2]),
                    (det as f64)
                        * (rot[2][0] as f64 * moments[i][0]
                            + rot[2][1] as f64 * moments[i][1]
                            + rot[2][2] as f64 * moments[i][2]),
                ];

                let preserved = (m_new[0] - moments[j][0]).abs() < symprec
                    && (m_new[1] - moments[j][1]).abs() < symprec
                    && (m_new[2] - moments[j][2]).abs() < symprec;
                let reversed = (m_new[0] + moments[j][0]).abs() < symprec
                    && (m_new[1] + moments[j][1]).abs() < symprec
                    && (m_new[2] + moments[j][2]).abs() < symprec;

                let this_tr = if preserved {
                    0
                } else if reversed {
                    1
                } else {
                    return -1;
                };

                match global_tr {
                    Some(tr) if tr != this_tr => return -1,
                    _ => global_tr = Some(this_tr),
                }
            }
            global_tr.unwrap_or(-1)
        })
        .collect()
}

// ========================================================================
// Internal functions
// ========================================================================

/// 从对称操作获取 Hall 编号。
pub(crate) fn identify_hall_number(
    rotations: &[Mat3I],
    translations: &[Vec3],
    lattice: &Mat3,
    transform_lattice_by_tmat: bool,
    symprec: f64,
) -> Result<usize, SymError> {
    let num_ops = rotations.len();
    if num_ops == 0 || translations.len() != num_ops {
        return Err(SymError::InvalidInput);
    }
    let mut symmetry = Symmetry::new(num_ops);
    symmetry.rot[..num_ops].copy_from_slice(&rotations[..num_ops]);
    symmetry.trans[..num_ops].copy_from_slice(&translations[..num_ops]);

    let (t_mat, prim_sym) =
        prm_get_primitive_symmetry(&symmetry, symprec).ok_or(SymError::SpacegroupSearchFailed)?;

    let prim_lat = if transform_lattice_by_tmat {
        let t_mat_inv = mat_inverse_matrix_d3(&t_mat, symprec)
            .ok()
            .ok_or(SymError::SpacegroupSearchFailed)?;
        mat_multiply_matrix_d3(lattice, &t_mat_inv)
    } else {
        *lattice
    };

    let spacegroup = spa_search_spacegroup_with_symmetry(&prim_sym, &prim_lat, symprec)?;
    Ok(spacegroup.hall_number)
}

/// 获取 SpaceGroupType。
fn get_spacegroup_type(hall_number: usize) -> Result<SpaceGroupType, SymError> {
    if hall_number == 0 || hall_number >= 531 {
        return Err(SymError::SpacegroupSearchFailed);
    }

    let spgtype = spgdb_get_spacegroup_type(hall_number);
    let pointgroup = ptg_get_pointgroup(spgtype.pointgroup_number);

    Ok(SpaceGroupType {
        number: spgtype.number,
        hall_number,
        schoenflies: spgtype.schoenflies,
        hall_symbol: spgtype.hall_symbol,
        choice: spgtype.choice,
        international: spgtype.international,
        international_full: spgtype.international_full,
        international_short: spgtype.international_short,
        pointgroup_international: pointgroup.symbol.to_string(),
        pointgroup_schoenflies: pointgroup.schoenflies.to_string(),
        arithmetic_crystal_class_number: 0, // TODO: arth_get_symbol
        arithmetic_crystal_class_symbol: String::new(),
    })
}
