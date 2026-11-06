//! Explicit physical line sources, separate from the original 73 W tables.
//!
//! The certified domain is the eight scalar beta lines in each of SG3 and SG4:
//! LD1LE1/LD2LE2,
//! U1UA1/U2UA2, V1VA1/V2VA2 and W1WA1/W2WA2. They have transverse offsets
//! (0,0), (1/2,1/2), (0,1/2), (1/2,0), respectively, and C2 point sign +1/-1.
//! Each has ONE actual scalar CIR (LD1/2, U1/2, V1/2, W1/2), not a stored
//! partner row. Its physical carrier pairs that CIR at +k with its derived
//! conjugate at -k, including the entire offset and point sign. Both component
//! arms enter one common folding, block solve and full character reconstruction.
//! Assigning dimension two to the +k table would instead describe two copies of
//! the same CIR and gives wrong translation/screw characters.
//!
//! At b in (1/2)Z the conjugate stars coincide; the two component identities are
//! retained. This is NOT spatial little-group enhancement: both little groups
//! are C2 for every b. In SG4 the half-beta screw characters are distinct +/-i,
//! and the full labelled parent period is TWO rather than one.
//! The result describes the continued two-dimensional carrier,
//! including its reducible specializations, not a new endpoint irrep label.
//!
//! This entry restricts to an explicitly supplied, validated fixed isotropy
//! context. It does not claim any line owns that context. There is no automatic
//! high-level routing, W-row insertion or general affine-source support. The
//! separate exact partition certifies geometry/character equivalence, not labels.
//! All pointwise evaluation uses exact rational b and the existing
//! Bloch convention; the result certifies characters, not matrix intertwiners.
//! Native archive hashes and the independent exact source oracle are documented
//! in [the contract](../../docs/physical-line-subduction.md).

use super::*;
use crate::irrep::w_little_characters_data::LittleOperation;

/// An immutable, independently audited physical line source.
///
/// Obtained only through [`physical_line_source`]. The frozen seed is a complex
/// component; its mathematical conjugate is not a second stored CIR row.
#[derive(Debug)]
pub struct PhysicalLineSource {
    label: &'static str,
    pir_irnumber: u32,
    seed_irnumber: u32,
    dimension: u32,
    // Point/direction storage only: this is NOT a complete offset CIR table.
    seed: &'static LittleCharacterTable,
    offset: [(i128, i128); 3],
}

impl PhysicalLineSource {
    /// Parent space group in the archive's conventional frame.
    pub const fn parent_sg(&self) -> u8 {
        self.seed.space_group
    }
    /// Exact physical PIR label, not a fabricated compiled IrrepRecord.
    pub const fn label(&self) -> &'static str {
        self.label
    }
    /// Authoritative PIR source number.
    pub const fn pir_irnumber(&self) -> u32 {
        self.pir_irnumber
    }
    /// The ONE authoritative CIR seed source number.
    pub const fn seed_irnumber(&self) -> u32 {
        self.seed_irnumber
    }
    /// Dimension of the continued physical carrier, including specializations.
    pub const fn dimension(&self) -> u32 {
        self.dimension
    }
    /// A complete ZERO-OFFSET complex seed table, not the physical carrier.
    ///
    /// Returns None for an affine source: LittleCharacterTable cannot express
    /// its offset. Callers of the original reference-returning accessor must
    /// now handle this Option; no incomplete U1 table escapes as a valid source.
    pub const fn seed_table(&self) -> Option<&'static LittleCharacterTable> {
        if self.offset[0].0 == 0 && self.offset[1].0 == 0 && self.offset[2].0 == 0 {
            Some(self.seed)
        } else {
            None
        }
    }
    /// Authoritative complex CIR seed label, including for an affine source.
    pub const fn seed_label(&self) -> &'static str {
        self.seed.label
    }
    /// Exact affine offset k0 in the archive's parent conventional frame.
    pub fn offset(&self) -> Vec3R {
        Vec3R::new(
            self.offset
                .map(|(n, d)| Rat::new(n, d).expect("certified source offset")),
        )
    }
    /// Exact coefficient of b, separate from the affine offset.
    pub fn direction(&self) -> Vec3R {
        line_direction(self.seed).expect("certified source direction")
    }
    /// Full labelled character period in the parent, not just the momentum period.
    /// P2 has period one; the P21 screw forces period two and swaps the two
    /// native source labels after a unit parameter shift.
    pub fn parent_character_period(&self) -> Rat {
        Rat::new(if self.parent_sg() == 4 { 2 } else { 1 }, 1)
            .expect("certified parent character period")
    }
    /// Existing provenance tags keep the same CIR number on both components.
    pub const fn components(&self) -> [SubductionComponent; 2] {
        [
            SubductionComponent::RealificationSeed {
                irnumber: self.seed_irnumber,
            },
            SubductionComponent::RealificationConjugate {
                irnumber: self.seed_irnumber,
            },
        ]
    }
}

const SG3_SCALAR_POINT_OPERATIONS: &[LittleOperation] = &[
    LittleOperation {
        element: "x,y,z",
        rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        translation: ["0", "0", "0"],
        character: [1, 0],
    },
    LittleOperation {
        element: "-x,y,-z",
        rotation: [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
        translation: ["0", "0", "0"],
        character: [1, 0],
    },
];

const SG3_NEGATIVE_POINT_OPERATIONS: &[LittleOperation] = &[
    SG3_SCALAR_POINT_OPERATIONS[0],
    LittleOperation {
        character: [-1, 0],
        ..SG3_SCALAR_POINT_OPERATIONS[1]
    },
];

// P21 has irtranslation equal to the stored screw translation. Thus the
// native t-tau+irtranslation is the full t already used by LineArmSource.
const SG4_SCALAR_POINT_OPERATIONS: &[LittleOperation] = &[
    SG3_SCALAR_POINT_OPERATIONS[0],
    LittleOperation {
        element: "-x,y+1/2,-z",
        translation: ["0", "1/2", "0"],
        ..SG3_SCALAR_POINT_OPERATIONS[1]
    },
];
const SG4_NEGATIVE_POINT_OPERATIONS: &[LittleOperation] = &[
    SG4_SCALAR_POINT_OPERATIONS[0],
    LittleOperation {
        character: [-1, 0],
        ..SG4_SCALAR_POINT_OPERATIONS[1]
    },
];

// Eight pinned SG3 scalar beta sources: irtype=3, one stored k arm, P centering,
// E/C2y point matrices I2/(epsilon I2), zero Seitz AND irtranslations.
// Hence the existing seed character D(R) exp(+2pi i k.t) is the native CIR.
// The original 73 W tables and their source/context lookup are not modified.
static SOURCES: &[PhysicalLineSource] = &[
    PhysicalLineSource {
        label: "LD1LE1",
        pir_irnumber: 29,
        seed_irnumber: 29,
        dimension: 2,
        offset: [(0, 1); 3],
        seed: &LittleCharacterTable {
            space_group: 3,
            label: "LD1",
            k_label: "LD",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG3_SCALAR_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "U1UA1",
        pir_irnumber: 45,
        seed_irnumber: 45,
        dimension: 2,
        offset: [(1, 2), (0, 1), (1, 2)],
        // Native PIR/CIR #45 has the same E/C2y point values and zero stored
        // translations/irtranslations. C2y*k-k=(-1,0,-1) is parent reciprocal;
        // thus its generic spatial star still has ONE arm. Only the complete
        // affine wave vector, never this private point table alone, evaluates U1.
        seed: &LittleCharacterTable {
            space_group: 3,
            label: "U1",
            k_label: "U",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG3_SCALAR_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "LD2LE2",
        pir_irnumber: 30,
        seed_irnumber: 30,
        dimension: 2,
        offset: [(0, 1), (0, 1), (0, 1)],
        seed: &LittleCharacterTable {
            space_group: 3,
            label: "LD2",
            k_label: "LD",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG3_NEGATIVE_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "U2UA2",
        pir_irnumber: 46,
        seed_irnumber: 46,
        dimension: 2,
        offset: [(1, 2), (0, 1), (1, 2)],
        seed: &LittleCharacterTable {
            space_group: 3,
            label: "U2",
            k_label: "U",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG3_NEGATIVE_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "V1VA1",
        pir_irnumber: 47,
        seed_irnumber: 47,
        dimension: 2,
        offset: [(0, 1), (0, 1), (1, 2)],
        seed: &LittleCharacterTable {
            space_group: 3,
            label: "V1",
            k_label: "V",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG3_SCALAR_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "V2VA2",
        pir_irnumber: 48,
        seed_irnumber: 48,
        dimension: 2,
        offset: [(0, 1), (0, 1), (1, 2)],
        seed: &LittleCharacterTable {
            space_group: 3,
            label: "V2",
            k_label: "V",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG3_NEGATIVE_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "W1WA1",
        pir_irnumber: 49,
        seed_irnumber: 49,
        dimension: 2,
        offset: [(1, 2), (0, 1), (0, 1)],
        seed: &LittleCharacterTable {
            space_group: 3,
            label: "W1",
            k_label: "W",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG3_SCALAR_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "W2WA2",
        pir_irnumber: 50,
        seed_irnumber: 50,
        dimension: 2,
        offset: [(1, 2), (0, 1), (0, 1)],
        seed: &LittleCharacterTable {
            space_group: 3,
            label: "W2",
            k_label: "W",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG3_NEGATIVE_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "LD1LE1",
        pir_irnumber: 56,
        seed_irnumber: 56,
        dimension: 2,
        offset: [(0, 1), (0, 1), (0, 1)],
        seed: &LittleCharacterTable {
            space_group: 4,
            label: "LD1",
            k_label: "LD",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG4_SCALAR_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "LD2LE2",
        pir_irnumber: 57,
        seed_irnumber: 57,
        dimension: 2,
        offset: [(0, 1), (0, 1), (0, 1)],
        seed: &LittleCharacterTable {
            space_group: 4,
            label: "LD2",
            k_label: "LD",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG4_NEGATIVE_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "U1UA1",
        pir_irnumber: 68,
        seed_irnumber: 72,
        dimension: 2,
        offset: [(1, 2), (0, 1), (1, 2)],
        seed: &LittleCharacterTable {
            space_group: 4,
            label: "U1",
            k_label: "U",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG4_SCALAR_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "U2UA2",
        pir_irnumber: 69,
        seed_irnumber: 73,
        dimension: 2,
        offset: [(1, 2), (0, 1), (1, 2)],
        seed: &LittleCharacterTable {
            space_group: 4,
            label: "U2",
            k_label: "U",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG4_NEGATIVE_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "V1VA1",
        pir_irnumber: 70,
        seed_irnumber: 74,
        dimension: 2,
        offset: [(0, 1), (0, 1), (1, 2)],
        seed: &LittleCharacterTable {
            space_group: 4,
            label: "V1",
            k_label: "V",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG4_SCALAR_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "V2VA2",
        pir_irnumber: 71,
        seed_irnumber: 75,
        dimension: 2,
        offset: [(0, 1), (0, 1), (1, 2)],
        seed: &LittleCharacterTable {
            space_group: 4,
            label: "V2",
            k_label: "V",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG4_NEGATIVE_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "W1WA1",
        pir_irnumber: 72,
        seed_irnumber: 76,
        dimension: 2,
        offset: [(1, 2), (0, 1), (0, 1)],
        seed: &LittleCharacterTable {
            space_group: 4,
            label: "W1",
            k_label: "W",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG4_SCALAR_POINT_OPERATIONS,
        },
    },
    PhysicalLineSource {
        label: "W2WA2",
        pir_irnumber: 73,
        seed_irnumber: 77,
        dimension: 2,
        offset: [(1, 2), (0, 1), (0, 1)],
        seed: &LittleCharacterTable {
            space_group: 4,
            label: "W2",
            k_label: "W",
            direction: ["0", "1", "0"],
            dimension: 1,
            operations: SG4_NEGATIVE_POINT_OPERATIONS,
        },
    },
];

/// Look up an explicitly certified physical line by its exact PIR label.
/// The SG3 and SG4 scalar beta cohorts retain distinct PIR/CIR namespaces,
/// their full affine offsets and native point signs.
pub fn physical_line_source(parent_sg: u8, label: &str) -> Option<&'static PhysicalLineSource> {
    SOURCES
        .iter()
        .find(|s| s.parent_sg() == parent_sg && s.label == label)
}

/// Flat arm identities are seed arms followed by their derived conjugates.
/// Conjugation applies to the WHOLE evaluated character, including Bloch phase.
pub(super) struct RealifiedLineSource<'a> {
    seed: LineArmSource<'a>,
}

impl RealifiedLineSource<'_> {
    fn seed_index(&self, index: usize) -> Result<(usize, bool), StarError> {
        let n = self.seed.arms.len();
        if index >= 2 * n {
            return Err(StarError::UnknownArmIndex { index, arms: 2 * n });
        }
        Ok((index % n, index >= n))
    }
    pub(super) fn q_block_dimension(&self, indices: &[usize]) -> Result<u32, StarError> {
        for &index in indices {
            self.seed_index(index)?;
        }
        self.seed.q_block_dimension(indices)
    }
    pub(super) fn q_block_character(
        &self,
        indices: &[usize],
        op: &ExactSeitz,
    ) -> Result<Complex64, StarError> {
        let mut value = Complex64::new(0.0, 0.0);
        for &index in indices {
            let (index, conjugate) = self.seed_index(index)?;
            let ch = self.seed.q_block_character(&[index], op)?;
            value += if conjugate { ch.conj() } else { ch };
        }
        Ok(value)
    }
    pub(super) fn full_character(&self, op: &ExactSeitz) -> Result<Complex64, StarError> {
        self.q_block_character(&(0..2 * self.seed.arms.len()).collect::<Vec<_>>(), op)
    }
}

/// Complete ordinary character decomposition of one continued physical line.
/// The borrowed record and embedding retain the exact restriction context.
#[derive(Debug)]
pub struct PhysicalLineSubduction<'a> {
    source: &'static PhysicalLineSource,
    subgroup: &'a IsotropySubgroup,
    embedding: &'a SubgroupEmbedding,
    parameter: Rat,
    wave_vectors: [Vec3R; 2],
    conjugate_stars_coincide: bool,
    blocks: Vec<FullStarBlock>,
    verification_operations: Vec<ExactSeitz>,
    parent_characters: Vec<Complex64>,
    reconstructed: Vec<Complex64>,
}
impl<'a> PhysicalLineSubduction<'a> {
    /// Immutable physical and single-CIR provenance.
    pub const fn source(&self) -> &'static PhysicalLineSource {
        self.source
    }
    /// Fixed record defining the restriction, NOT a record owned by this line.
    pub const fn subgroup(&self) -> &'a IsotropySubgroup {
        self.subgroup
    }
    /// Validated exact geometry, including setting and origin.
    pub const fn embedding(&self) -> &'a SubgroupEmbedding {
        self.embedding
    }
    /// Explicit rational b, unreduced.
    pub const fn parameter(&self) -> &Rat {
        &self.parameter
    }
    /// Seed and derived-conjugate momenta; neither is silently reduced.
    pub const fn wave_vectors(&self) -> &[Vec3R; 2] {
        &self.wave_vectors
    }
    /// Whether the two parent stars coincide modulo the parent reciprocal lattice.
    /// This is a geometric collision, not spatial little-group enhancement.
    pub const fn conjugate_stars_coincide(&self) -> bool {
        self.conjugate_stars_coincide
    }
    /// Joint child-star blocks. Coincident component arms remain counted twice.
    pub fn blocks(&self) -> &[FullStarBlock] {
        &self.blocks
    }
    /// Dimension accounted for by all blocks; always the physical dimension.
    pub fn covered_dimension(&self) -> u32 {
        self.blocks.iter().map(FullStarBlock::block_dimension).sum()
    }
    /// Parent-frame verification operations: every representative translated by
    /// coefficients 0,1,2 of all three subgroup lattice generators. This includes
    /// pure translations, their squares and mixed translations; point characters
    /// alone would miss an incorrectly doubled +k component.
    pub fn verification_operations(&self) -> &[ExactSeitz] {
        &self.verification_operations
    }
    /// Parent and independently reconstructed characters on the verification operations.
    pub fn reconstruction(&self) -> (&[Complex64], &[Complex64]) {
        (&self.parent_characters, &self.reconstructed)
    }
    /// Existing ordinary-character comparison tolerance.
    pub const fn tolerance(&self) -> f64 {
        SUBDUCTION_TOLERANCE
    }
}

// Keep the same subgroup representatives; add exact L_H translates without
// changing the child setting or origin. Degree-two mixed moments also distinguish
// opposite Bloch eigenvalues whose first trace vanishes (e.g. b=1/4).
fn verification_operations(
    embedding: &SubgroupEmbedding,
) -> Result<Vec<ExactSeitz>, SubductionError> {
    let basis = embedding.subgroup_lattice().rows().rows();
    let mut operations = Vec::new();
    for representative in embedding.representatives() {
        for a in 0..=2 {
            for b in 0..=2 {
                for c in 0..=2 {
                    let mut translation = *representative.translation().as_array();
                    for (row, coefficient) in basis.iter().zip([a, b, c]) {
                        for (entry, value) in translation.iter_mut().zip(row) {
                            *entry =
                                entry.checked_add(value.checked_mul(Rat::new(coefficient, 1)?)?)?;
                        }
                    }
                    operations.push(ExactSeitz::new(
                        representative.rotation(),
                        Vec3R::new(translation),
                    ));
                }
            }
        }
    }
    Ok(operations)
}

/// Restrict a certified physical line to an explicit fixed record/embedding at b.
///
/// Uses one shared fold/solve/reconstruct, never concatenates already completed
/// component reports. Special parameters retain the full physical carrier;
/// stored target labels are reported only when the existing catalogue proves them.
/// Context, missing data and exact-arithmetic errors retain their original causes.
///
/// `Z2` owns this fixed geometry; `LD1LE1` is the separate source being restricted.
/// `b=1/4` folds both momenta to one SG4 point, but the screw has opposite
/// complex eigenvalues, so two different stored CIR targets must survive.
///
/// ```
/// use cryspglib::irrep::{isotropy, LabelConvention};
/// use cryspglib::irrep::subduction::{Rat, SubgroupEmbedding};
/// use cryspglib::irrep::subduction::star::decompose::physical_line::{
///     physical_line_source, subduce_physical_line_with_embedding,
/// };
/// let context = isotropy::isotropy_subgroups(3, "Z2", LabelConvention::Cdml)?[0];
/// let embedding = SubgroupEmbedding::from_isotropy_subgroup(&context)?;
/// let source = physical_line_source(3, "LD1LE1").unwrap();
/// let report = subduce_physical_line_with_embedding(
///     &context, &embedding, source, Rat::new(1, 4)?,
/// )?;
/// assert_eq!(report.covered_dimension(), 2);
/// assert_eq!(report.blocks().len(), 1);
/// assert_eq!(report.blocks()[0].targets().len(), 2);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn subduce_physical_line_with_embedding<'a>(
    subgroup: &'a IsotropySubgroup,
    embedding: &'a SubgroupEmbedding,
    source: &'static PhysicalLineSource,
    parameter: Rat,
) -> Result<PhysicalLineSubduction<'a>, FullStarError> {
    validate_line_context(subgroup, embedding, source.seed)?;
    let direction = line_direction(source.seed)?;
    let arms = line_arms(subgroup, embedding, source.seed, &direction)?;
    let k = line_wave_vector(&direction, &parameter)?.checked_add(&source.offset())?;
    let negative_k = k.checked_neg()?;
    let paired = RealifiedLineSource {
        seed: LineArmSource {
            table: source.seed,
            direction,
            wave_vector: k,
            arms: &arms,
        },
    };
    let mut seed_arms = Vec::with_capacity(arms.len());
    for (_, rotation) in &arms {
        seed_arms.push(
            Mat3R::from_ints(*rotation)
                .inverse()?
                .transpose()
                .checked_mul_vector(&k)?,
        );
    }
    let reciprocal = embedding.parent_lattice().reciprocal()?;
    let mut coincident = false;
    for a in &seed_arms {
        for b in &seed_arms {
            coincident |= reciprocal.same_mod(a, &b.checked_neg()?)?;
        }
    }
    let mut folded_inputs = Vec::with_capacity(2 * arms.len());
    for conjugate in [false, true] {
        for arm in &seed_arms {
            folded_inputs.push(FoldArm {
                wave_vector: if conjugate { arm.checked_neg()? } else { *arm },
                dimension: usize::from(source.seed.dimension),
            });
        }
    }
    let folded = fold_arms(embedding, source.parent_sg(), &folded_inputs)?;
    let child_cell = Lattice::new(exact_primitive_basis(embedding.subgroup_sg())?)?;
    let child_reciprocal = child_cell.reciprocal()?;
    let evaluator = ArmCharacterSource::RealifiedLine(&paired);
    let mut blocks = Vec::with_capacity(folded.len());
    let mut covered = 0u32;
    for star in &folded {
        let block = build_block(embedding, &evaluator, star, &child_cell, &child_reciprocal)?;
        covered = covered.checked_add(block.block_dimension()).ok_or(
            SubductionError::RationalOverflow {
                operation: "physical line dimension",
            },
        )?;
        blocks.push(block);
    }
    if covered != source.dimension {
        return Err(FullStarError::TotalDimensionMismatch {
            expected: source.dimension,
            found: covered,
        });
    }
    let verification_operations = verification_operations(embedding)?;
    let (parent_characters, reconstructed) =
        reconstruct_on(embedding, &evaluator, &blocks, &verification_operations)?;
    Ok(PhysicalLineSubduction {
        source,
        subgroup,
        embedding,
        parameter,
        wave_vectors: [k, negative_k],
        conjugate_stars_coincide: coincident,
        blocks,
        verification_operations,
        parent_characters,
        reconstructed,
    })
}

/// The exact complex-character content of the restricted physical pair.
///
/// Both components are scalar for every registered P2/P21 source. A momentum
/// collision is not necessarily equality of their characters on the whole
/// subgroup: the SG4 screw distinguishes them at quarter parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhysicalLineKind {
    /// Two different singleton child stars, each with one scalar target.
    DistinctMomenta,
    /// One child momentum, with two inequivalent scalar targets of multiplicity one.
    DistinctCharactersAtOneMomentum,
    /// One scalar target with multiplicity two; both source copies survive.
    RepeatedCharacter,
}

/// One exact point cell within the source's full parent character period.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PhysicalLineBoundary {
    /// Parameter in [0,period). Only the full-period endpoint is represented by zero.
    pub parameter: Rat,
    /// Content at this point, not the content of an adjacent open interval.
    pub kind: PhysicalLineKind,
    /// Whether the two momenta also coincide modulo the *parent* reciprocal lattice.
    /// This is not a change of the parent's little co-group (always C2 here).
    pub parent_conjugate_stars_coincide: bool,
}

/// Exact geometric and scalar-character partition of the registered physical line.
///
/// Every open interval between the boundaries has two different singleton child
/// stars, one scalar target in each. The point cells retain both source copies.
/// This is a proof from exact congruences, not an inference from sampled calls.
/// It does not assign stored labels or continuous catalogue indices to an interval:
/// use the pointwise subduction to obtain target identities and numerical characters.
#[derive(Debug, Clone)]
pub struct PhysicalLinePartition<'a> {
    subgroup: &'a IsotropySubgroup,
    embedding: &'a SubgroupEmbedding,
    source: &'static PhysicalLineSource,
    parent_collision_step: Rat,
    momentum_collision_step: Rat,
    character_equivalence_step: Rat,
    boundaries: Vec<PhysicalLineBoundary>,
}

impl<'a> PhysicalLinePartition<'a> {
    /// The fixed record defining restriction geometry, not the physical source's owner.
    pub const fn subgroup(&self) -> &'a IsotropySubgroup {
        self.subgroup
    }
    /// The validated restriction embedding.
    pub const fn embedding(&self) -> &'a SubgroupEmbedding {
        self.embedding
    }
    /// Source identity, including the seed and derived-conjugate identities.
    pub const fn source(&self) -> &'static PhysicalLineSource {
        self.source
    }
    /// Positive generator of the parameters where parent momenta coincide (1/2).
    pub const fn parent_collision_step(&self) -> Rat {
        self.parent_collision_step
    }
    /// Positive generator of child-momentum collisions. Twice this is the period
    /// of each labelled component's translation character, not necessarily its
    /// full subgroup character (SG4 is the counterexample).
    pub const fn momentum_collision_step(&self) -> Rat {
        self.momentum_collision_step
    }
    /// Positive generator of parameters where the two full restricted characters
    /// agree. Twice this is their labelled full-character period. Doubled SG3
    /// has step 1/4 (full period 1/2); SG4 has step 1/2 (full period 1), despite
    /// the same momentum period 1/2.
    pub const fn character_equivalence_step(&self) -> Rat {
        self.character_equivalence_step
    }
    /// Point cells in ascending order within the source's full parent period.
    pub fn boundaries(&self) -> &[PhysicalLineBoundary] {
        &self.boundaries
    }
    /// Find a point cell modulo the full parent character period (one or two).
    /// This exact classification does not evaluate a potentially overflowing
    /// unreduced Bloch vector, nor change the pointwise evaluator's arithmetic contract.
    pub fn boundary(&self, parameter: &Rat) -> Option<&PhysicalLineBoundary> {
        let period = self.source.parent_character_period().numerator();
        // A canonical mod-two Rat need not fit i128: -1/MAX is a witness.
        // Keep the integer turn and fractional residue separate; never form 2*d.
        let turn = parameter
            .numerator()
            .div_euclid(parameter.denominator())
            .rem_euclid(period);
        let reduced = super::super::line_domain::reduce_modulo_one(*parameter);
        self.boundaries.iter().find(|b| {
            b.parameter
                .numerator()
                .div_euclid(b.parameter.denominator())
                == turn
                && super::super::line_domain::reduce_modulo_one(b.parameter) == reduced
        })
    }
    /// Exact content at any representable rational parameter, including negative
    /// and period-translated values. No tolerance or near-boundary snapping.
    pub fn kind_at(&self, parameter: &Rat) -> PhysicalLineKind {
        self.boundary(parameter)
            .map_or(PhysicalLineKind::DistinctMomenta, |b| b.kind)
    }
    /// Open interval endpoints in [0,period]. All interiors have DistinctMomenta content;
    /// endpoints must be evaluated as point cells, not by an interior sample.
    pub fn intervals(&self) -> Vec<(Rat, Rat)> {
        self.boundaries
            .iter()
            .map(|b| b.parameter)
            .zip(
                self.boundaries
                    .iter()
                    .skip(1)
                    .map(|b| b.parameter)
                    .chain(std::iter::once(self.source.parent_character_period())),
            )
            .collect()
    }
}

/// Partition any certified P2/P21 scalar physical pair over its full parent period.
///
/// Here D±({R|t})=epsilon(R) exp(±2πi(k0.t+b t_y)), epsilon(R) in {1,-1}.
/// The same real point sign cancels when comparing the two characters. Since
/// 2k0.t is integral, the offset cancels too (not from their values). Every subgroup
/// rotation fixes y. Thus each restricted component is already a scalar irrep.
/// Child momenta coincide iff 2b e_y is in L_H*, computed with the existing exact
/// lattice-coordinate step solver. The full characters coincide only if, in
/// addition, 2b t_y is integral for every coset representative. Intersecting
/// these exact cyclic subgroups of Q gives all point cells and multiplicities.
/// No little-group enhancement, fake parent rotation for conjugation, sampled
/// completeness argument, or single-CIR ParameterKind is used.
///
/// The original 73-source partition and all pointwise routes are unchanged.
/// This function returns content/geometry, not a general affine-line partition,
/// stored-label interpolation, numerical matrices or gauge/monodromy transport.
///
/// ```
/// use cryspglib::irrep::{isotropy, LabelConvention};
/// use cryspglib::irrep::subduction::{Rat, SubgroupEmbedding};
/// use cryspglib::irrep::subduction::star::decompose::physical_line::{
///     physical_line_source, partition_physical_line_with_embedding, PhysicalLineKind,
/// };
/// let record = isotropy::isotropy_subgroups(3, "Z2", LabelConvention::Cdml)?[0];
/// let embedding = SubgroupEmbedding::from_isotropy_subgroup(&record)?;
/// let partition = partition_physical_line_with_embedding(
///     &record, &embedding, physical_line_source(3, "LD1LE1").unwrap(),
/// )?;
/// assert_eq!(partition.boundaries().len(), 4);
/// assert_eq!(partition.kind_at(&Rat::new(1, 4)?),
///            PhysicalLineKind::DistinctCharactersAtOneMomentum);
/// assert_eq!(partition.kind_at(&Rat::new(1, 2)?), PhysicalLineKind::RepeatedCharacter);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn partition_physical_line_with_embedding<'a>(
    subgroup: &'a IsotropySubgroup,
    embedding: &'a SubgroupEmbedding,
    source: &'static PhysicalLineSource,
) -> Result<PhysicalLinePartition<'a>, FullStarError> {
    use super::super::line_domain as domain;
    validate_line_context(subgroup, embedding, source.seed)?;
    let period = source.parent_character_period();
    // Certified P2/P21 sources have real scalar point coefficients and half-
    // reciprocal transverse k0. Integer x/z translations cancel k0 from BOTH
    // comparison equations, even when the screw has half-integral y.
    // This is not a general affine congruence solver.
    let direction = line_direction(source.seed)?;
    let difference = direction.checked_add(&direction)?;
    let parent_reciprocal = embedding.parent_lattice().reciprocal()?;
    let child_reciprocal_in_parent = embedding.subgroup_lattice().reciprocal()?;
    if !parent_reciprocal.contains(&source.offset().checked_add(&source.offset())?)? {
        return Err(SubductionError::DomainCensusInconsistent {
            reason: "physical P2 offset must be half reciprocal",
        }
        .into());
    }
    domain::require_reciprocal_direction(
        &parent_reciprocal,
        &direction,
        source.parent_sg(),
        source.label(),
    )?;
    domain::require_reciprocal_direction(
        &child_reciprocal_in_parent,
        &direction,
        embedding.subgroup_sg(),
        source.label(),
    )?;
    let missing_step = || SubductionError::DomainCensusInconsistent {
        reason: "physical line direction must be nonzero",
    };
    let parent_collision_step =
        domain::minimal_parameter_step_via_coordinates(&parent_reciprocal, &difference)?
            .ok_or_else(missing_step)?;
    let momentum_collision_step =
        domain::minimal_parameter_step_via_coordinates(&child_reciprocal_in_parent, &difference)?
            .ok_or_else(missing_step)?;
    let mut character_equivalence_step = momentum_collision_step;
    for operation in embedding.representatives() {
        // y is an additive homomorphism on this P2/P21 subgroup. Translations and
        // one representative of each point class therefore test the entire group.
        let y = operation.translation.get(1);
        if operation.translation.get(0).denominator() != 1
            || operation.translation.get(2).denominator() != 1
            || period.numerator() % y.denominator() != 0
        {
            return Err(SubductionError::DomainCensusInconsistent {
                reason: if source.parent_sg() == 3 {
                    "physical P2 character must have parameter period one"
                } else {
                    "physical P21 requires integer transverse and half-integer axial translations"
                },
            }
            .into());
        }
        let constraints = Vec3R::new([
            Rat::ONE.checked_div(character_equivalence_step)?,
            y.checked_add(y)?,
            Rat::ZERO,
        ]);
        character_equivalence_step =
            domain::minimal_parameter_step_via_coordinates(&Lattice::integer(), &constraints)?
                .ok_or_else(missing_step)?;
    }
    // Every step divides the certified parent period. Both parent collisions and
    // character equivalence imply child collision, but in P21 a parent collision
    // need NOT imply character equivalence (the half-beta screw has +/-i).
    if [
        parent_collision_step,
        momentum_collision_step,
        character_equivalence_step,
    ]
    .iter()
    .any(|s| s.numerator() != 1)
        || parent_collision_step
            .checked_div(momentum_collision_step)?
            .denominator()
            != 1
        || character_equivalence_step
            .checked_div(momentum_collision_step)?
            .denominator()
            != 1
    {
        return Err(SubductionError::DomainCensusInconsistent {
            reason: "physical-pair collision steps must nest within the certified period",
        }
        .into());
    }
    let residues = domain::residues_of(momentum_collision_step)?;
    let mut boundaries = Vec::new();
    for turn in 0..period.numerator() {
        for residue in &residues {
            let parameter = residue.checked_add(Rat::new(turn, 1)?)?;
            boundaries.push(PhysicalLineBoundary {
                parameter,
                kind: if parameter
                    .checked_div(character_equivalence_step)?
                    .denominator()
                    == 1
                {
                    PhysicalLineKind::RepeatedCharacter
                } else {
                    PhysicalLineKind::DistinctCharactersAtOneMomentum
                },
                parent_conjugate_stars_coincide: parameter
                    .checked_div(parent_collision_step)?
                    .denominator()
                    == 1,
            });
        }
    }
    Ok(PhysicalLinePartition {
        subgroup,
        embedding,
        source,
        parent_collision_step,
        momentum_collision_step,
        character_equivalence_step,
        boundaries,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::irrep::{LabelConvention, isotropy};

    #[test]
    fn translated_checks_reject_wrong_pair_geometry_even_when_point_checks_pass() {
        let subgroup = isotropy::isotropy_subgroups(3, "GM1", LabelConvention::Cdml).unwrap()[0];
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
        let table = physical_line_source(3, "LD1LE1").unwrap().seed;
        let direction = line_direction(table).unwrap();
        let arms = line_arms(&subgroup, &embedding, table, &direction).unwrap();
        let cell = Lattice::new(exact_primitive_basis(3).unwrap()).unwrap();
        let reciprocal = cell.reciprocal().unwrap();
        for (b, wrong_parameters) in [
            (Rat::new(1, 7).unwrap(), [Rat::new(1, 7).unwrap(); 2]),
            (
                Rat::new(1, 4).unwrap(),
                [Rat::ZERO, Rat::new(1, 2).unwrap()],
            ),
        ] {
            let paired = RealifiedLineSource {
                seed: LineArmSource {
                    table,
                    direction,
                    wave_vector: line_wave_vector(&direction, &b).unwrap(),
                    arms: &arms,
                },
            };
            let source = ArmCharacterSource::RealifiedLine(&paired);
            let wrong: Vec<_> = wrong_parameters
                .iter()
                .map(|p| FoldArm {
                    wave_vector: line_wave_vector(&direction, p).unwrap(),
                    dimension: 1,
                })
                .collect();
            let folded = fold_arms(&embedding, 3, &wrong).unwrap();
            let blocks: Vec<_> = folded
                .iter()
                .map(|s| build_block(&embedding, &source, s, &cell, &reciprocal).unwrap())
                .collect();
            // Both malformed geometries pass dimension and pure-point characters.
            reconstruct(&embedding, &source, &blocks).unwrap();
            if b == Rat::new(1, 4).unwrap() {
                // Even the first translation moment aliases: i + (-i) = 1 + (-1).
                let ty = ExactSeitz::new(IDENTITY_ROTATION, direction);
                reconstruct_on(&embedding, &source, &blocks, &[ty]).unwrap();
            }
            let ops = verification_operations(&embedding).unwrap();
            let error = reconstruct_on(&embedding, &source, &blocks, &ops).unwrap_err();
            let FullStarError::ReconstructionMismatch { index, .. } = error else {
                panic!("{error:?}");
            };
            if b == Rat::new(1, 4).unwrap() {
                assert_eq!(ops[index].translation().get(1), Rat::new(2, 1).unwrap());
            }
        }
    }

    #[test]
    fn derived_component_conjugates_the_entire_coefficient_and_rejects_bad_indices() {
        // Synthetic valid C4 character, NOT another registered/archive source.
        static TABLE: LittleCharacterTable = LittleCharacterTable {
            space_group: 75,
            label: "synthetic-C4",
            k_label: "synthetic",
            direction: ["0", "0", "1"],
            dimension: 1,
            operations: &[
                LittleOperation {
                    element: "x,y,z",
                    rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    translation: ["0"; 3],
                    character: [1, 0],
                },
                LittleOperation {
                    element: "-y,x,z",
                    rotation: [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
                    translation: ["0"; 3],
                    character: [0, 1],
                },
                LittleOperation {
                    element: "-x,-y,z",
                    rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
                    translation: ["0"; 3],
                    character: [-1, 0],
                },
                LittleOperation {
                    element: "y,-x,z",
                    rotation: [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
                    translation: ["0"; 3],
                    character: [0, -1],
                },
            ],
        };
        let direction = line_direction(&TABLE).unwrap();
        let arms = [(direction, IDENTITY_ROTATION)];
        let source = RealifiedLineSource {
            seed: LineArmSource {
                table: &TABLE,
                direction,
                wave_vector: line_wave_vector(&direction, &Rat::new(1, 8).unwrap()).unwrap(),
                arms: &arms,
            },
        };
        let op = ExactSeitz::new([[0, -1, 0], [1, 0, 0], [0, 0, 1]], direction);
        let z = source.q_block_character(&[0], &op).unwrap();
        let conjugate = source.q_block_character(&[1], &op).unwrap();
        assert!((z - Complex64::new(-1.0, 1.0) / 2.0_f64.sqrt()).norm() < 1e-12);
        assert!((conjugate - z.conj()).norm() < 1e-12);
        assert!(
            (source.full_character(&op).unwrap() - Complex64::new(-2.0_f64.sqrt(), 0.0)).norm()
                < 1e-12
        );
        assert_eq!(source.q_block_dimension(&[0, 1]).unwrap(), 2);
        assert!(matches!(
            source.q_block_dimension(&[2]),
            Err(StarError::UnknownArmIndex { index: 2, arms: 2 })
        ));
        assert!(matches!(
            source.q_block_character(&[2], &op),
            Err(StarError::UnknownArmIndex { index: 2, arms: 2 })
        ));
    }
}
