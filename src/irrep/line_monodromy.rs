//! Monodromy of the frozen parametric-k line sources under a reciprocal shift.
//!
//! A parametric-k row lives on a line `k(t) = t * v` whose direction `v` is a
//! **parent reciprocal lattice vector**. So `t` and `t + 1` name the same point
//! of the parent's Brillouin zone.
//!
//! They do **not** name the same frozen table.  The frozen
//! [`LittleCharacterTable`]s store the little co-group matrices `D(R)` solved at
//! `k = Gamma` (no Bloch factor); the character of the band `alpha` at the
//! parameter `t` is
//!
//! ```text
//! chi_alpha^t(R, T_R) = D_alpha(R) * exp(2 pi i t (v . T_R)).
//! ```
//!
//! For a reciprocal shift `K`, multiplying by the **reciprocal-shift twist**
//! `Phi_K(R) = exp(2 pi i K.T_R)` gives the frozen character at `k + K`. When
//! `K = delta * v` for a source's direction `v`, this is the parameter shift
//! `t -> t + delta`. The twist must be a one-dimensional character of that
//! source's little group; when it is, the shifted character is another label:
//!
//! ```text
//! rho_{k+K, M_K(alpha)}  ~  rho_{k, alpha},        K in L*_G.
//! ```
//!
//! `M_K` is the **monodromy map** of a line little group. For operations
//! `g=(R_g,t_g)` and `h=(R_h,t_h)`, the exact condition is
//! `((I-R_g^T)K)·t_h ∈ Z` for every pair; requiring `R_g^{-T}K = K` is
//! sufficient but unnecessarily restrictive. A reciprocal shift that fails
//! the exact character test is returned as [`LabelImage::UnsupportedShift`]
//! rather than guessed from the phase formula. On a nonsymmorphic line, valid
//! monodromy need not be the identity. Along each source's own frozen
//! direction, SG 203/210 swap `DT1 <-> DT2` and `DT3 <-> DT4`, while SG
//! 227/228 swap `DT1 <-> DT3` and `DT2 <-> DT4`; SG 209 and the `SM` sources
//! are unchanged. The `DT3`/`DT4` complex-conjugation map in SG 209/210 is a
//! separate operation. Two consequences are used by the audit:
//!
//! * when the line step `v` is reciprocal and its twist is a character, the
//!   character data of label `alpha` at `t + 1` **is** the character data of
//!   label `M_v(alpha)` at `t`, so a decomposition transported by one parameter
//!   step must be compared against `M_v(alpha)` and not against `alpha` again;
//! * the trivial content is a gauge-invariant number, so it is constant on
//!   monodromy orbits of the *pinned* table -- that is the checkable half of the
//!   contract and it is what [`LineMonodromy::orbit`] feeds.
//!
//! # How the map is computed
//!
//! Never from label names. The shifted character vector `D_alpha(R) * Phi_K(R)`
//! is matched exactly against every frozen table of the same parent, and the
//! unique match is `M_K(alpha)`. The rational phase is reduced modulo one and
//! checked as an exact Gaussian-integer character; there is no floating-point
//! matching tolerance. That covers real sources, the `DT3`/`DT4` swap and any
//! cycle representable by the frozen Gaussian phases with one mechanism. A
//! label whose image is not unique (two
//! frozen tables sharing a fingerprint) or absent is reported as
//! [`LabelImage::Ambiguous`] / [`LabelImage::Missing`] rather than guessed.
//!
//! The fingerprints are keyed by the exact little-group **operation** (rotation
//! plus fractional translation), so two operations that share a rotation stay
//! distinct and the map is exact rather than tolerance-limited in its keys.
//!
//! See `docs/subduction-conventions.md` section 16 for the contract and its
//! verified parameter domain.

use crate::irrep::subduction::{
    Lattice, Mat3R, Rat, SubductionError, Vec3R, exact_primitive_basis,
};
use crate::irrep::w_little_characters_data::{
    LittleCharacterTable, LittleOperation, W_LITTLE_CHARACTERS,
};
use std::collections::{BTreeMap, BTreeSet};

/// Exact key of one frozen little-group operation: the rotation plus the
/// translation as exact `(numerator, denominator)` pairs.
///
/// A rotation alone is *not* a key (two coset representatives of one little
/// group can share it); the translation is what makes the key the operation.
type OperationKey = ([[i8; 3]; 3], [(i128, i128); 3]);

/// A frozen character value together with its exact rational phase in turns.
#[derive(Debug, Clone, Copy)]
struct PhasedCharacter {
    value: [i32; 2],
    phase: Rat,
}

/// A character vector of one frozen table, keyed by operation.
type Fingerprint = BTreeMap<OperationKey, PhasedCharacter>;

/// Where one label of a frozen table goes under a monodromy map.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LabelImage {
    /// Exactly one frozen table of the parent carries the shifted fingerprint.
    Unique(&'static str),
    /// Several do: the frozen data does not distinguish them, and the map is not
    /// guessed.
    Ambiguous(Vec<&'static str>),
    /// No frozen table of the parent carries the shifted fingerprint.
    Missing,
    /// The shift is reciprocal, but its exact Bloch phase does not define a
    /// one-dimensional character on this source's little group.
    UnsupportedShift,
}

impl LabelImage {
    /// The image when it is uniquely determined.
    pub fn unique(&self) -> Option<&'static str> {
        match self {
            Self::Unique(label) => Some(label),
            Self::Ambiguous(_) | Self::Missing | Self::UnsupportedShift => None,
        }
    }

    /// Whether the frozen data determines this image.
    pub const fn is_unique(&self) -> bool {
        matches!(self, Self::Unique(_))
    }
}

/// The monodromy map `M_K` of one parent, one shift.
///
/// Built by [`monodromy`] (a reciprocal shift in the parent's conventional frame)
/// or by [`complex_conjugation`] (the label map induced by `chi -> conj(chi)`).
#[derive(Debug, Clone)]
pub struct LineMonodromy {
    parent: u8,
    shift: Vec3R,
    conjugate: bool,
    images: BTreeMap<&'static str, LabelImage>,
}

impl LineMonodromy {
    /// Parent space group number.
    pub const fn parent(&self) -> u8 {
        self.parent
    }

    /// The shift `K` this map was built for, in the parent's conventional
    /// reciprocal basis (the frame the frozen directions live in).
    pub const fn shift(&self) -> &Vec3R {
        &self.shift
    }

    /// Whether this is the conjugation map `C` rather than a shift `M_K`.
    pub const fn is_conjugation(&self) -> bool {
        self.conjugate
    }

    /// Every frozen label of the parent and its image.
    pub fn images(&self) -> impl Iterator<Item = (&'static str, &LabelImage)> {
        self.images.iter().map(|(label, image)| (*label, image))
    }

    /// The image of one label, or `None` if the parent has no such frozen
    /// source.
    pub fn image(&self, label: &str) -> Option<&LabelImage> {
        self.images.get(label)
    }

    /// The image of one label when the frozen data determines it.
    pub fn unique_image(&self, label: &str) -> Option<&'static str> {
        self.images.get(label).and_then(LabelImage::unique)
    }

    /// The labels whose image the frozen data does not determine.
    pub fn undetermined(&self) -> Vec<(&'static str, &LabelImage)> {
        self.images
            .iter()
            .filter(|(_, image)| !image.is_unique())
            .map(|(label, image)| (*label, image))
            .collect()
    }

    /// The orbit `label, M(label), M^2(label), ...` of `steps` applications.
    ///
    /// `None` if any step is not uniquely determined, so a caller never marches
    /// through a guessed image.  The returned vector has `steps + 1` entries and
    /// starts with `label` itself.
    pub fn orbit(&self, label: &'static str, steps: usize) -> Option<Vec<&'static str>> {
        let mut orbit = Vec::with_capacity(steps + 1);
        let mut current = label;
        orbit.push(label);
        for _ in 0..steps {
            current = match self.images.get(current)? {
                LabelImage::Unique(image) => image,
                LabelImage::Ambiguous(_) | LabelImage::Missing | LabelImage::UnsupportedShift => {
                    return None;
                }
            };
            orbit.push(current);
        }
        Some(orbit)
    }
}

/// The parents that carry frozen parametric-k sources.
pub fn line_parents() -> Vec<u8> {
    let mut parents: Vec<u8> = W_LITTLE_CHARACTERS
        .iter()
        .map(|table| table.space_group)
        .collect();
    parents.sort_unstable();
    parents.dedup();
    parents
}

/// Every frozen source of one parent, in frozen order.
pub fn line_sources(parent: u8) -> Vec<&'static LittleCharacterTable> {
    W_LITTLE_CHARACTERS
        .iter()
        .filter(|table| table.space_group == parent)
        .collect()
}

/// One frozen source by parent and label.
pub fn line_table(parent: u8, label: &str) -> Option<&'static LittleCharacterTable> {
    W_LITTLE_CHARACTERS
        .iter()
        .find(|table| table.space_group == parent && table.label == label)
}

/// A frozen fraction such as `"-1/2"`, `"0"` or `"3"` as an exact rational.
pub fn parse_fraction(text: &str) -> Option<Rat> {
    let text = text.trim();
    let (numerator, denominator) = text.split_once('/').unwrap_or((text, "1"));
    let numerator = numerator.trim().parse::<i128>().ok()?;
    let denominator = denominator.trim().parse::<i128>().ok()?;
    Rat::new(numerator, denominator).ok()
}

/// The frozen direction of one source, in the parent's conventional reciprocal
/// basis.
pub fn line_direction(table: &LittleCharacterTable) -> Option<Vec3R> {
    let mut values = [Rat::ZERO; 3];
    for (axis, text) in table.direction.iter().enumerate() {
        values[axis] = parse_fraction(text)?;
    }
    Some(Vec3R::new(values))
}

/// The fractional translation of one frozen operation, in the parent's
/// conventional cell.
pub fn operation_translation(operation: &LittleOperation) -> Option<Vec3R> {
    let mut values = [Rat::ZERO; 3];
    for (axis, text) in operation.translation.iter().enumerate() {
        values[axis] = parse_fraction(text)?;
    }
    Some(Vec3R::new(values))
}

/// The exact operation key of one frozen operation.
fn operation_key(operation: &LittleOperation) -> Option<OperationKey> {
    let translation = operation_translation(operation)?;
    let mut key = [(0i128, 1i128); 3];
    for (axis, value) in key.iter_mut().enumerate() {
        *value = (
            translation.get(axis).numerator(),
            translation.get(axis).denominator(),
        );
    }
    Some((operation.rotation, key))
}

/// `shift . translation`, the exact phase in turns for one operation.
fn phase_turns(shift: &Vec3R, translation: &Vec3R) -> Result<Rat, SubductionError> {
    let mut turns = Rat::ZERO;
    for axis in 0..3 {
        turns = turns.checked_add(shift.get(axis).checked_mul(translation.get(axis))?)?;
    }
    Ok(turns)
}

/// One frozen table's character vector, optionally shifted by the twist `K` and
/// optionally conjugated.
fn fingerprint(
    table: &LittleCharacterTable,
    shift: Option<&Vec3R>,
    conjugate: bool,
) -> Result<Option<Fingerprint>, SubductionError> {
    let mut out = Fingerprint::new();
    for operation in table.operations {
        let Some(key) = operation_key(operation) else {
            return Ok(None);
        };
        let re = operation.character[0];
        let im = if conjugate {
            let Some(imaginary) = operation.character[1].checked_neg() else {
                return Ok(None);
            };
            imaginary
        } else {
            operation.character[1]
        };
        let value = [re, im];
        let phase = if value == [0, 0] {
            Rat::ZERO
        } else if let Some(shift) = shift {
            let Some(translation) = operation_translation(operation) else {
                return Ok(None);
            };
            phase_turns(shift, &translation)?
        } else {
            Rat::ZERO
        };
        out.insert(key, PhasedCharacter { value, phase });
    }
    Ok(Some(out))
}

/// The exact quarter-turn phase, reduced modulo one, if it is a Gaussian unit.
fn quarter_turn(phase: Rat) -> Option<u8> {
    let residue = Rat::new(
        phase.numerator().rem_euclid(phase.denominator()),
        phase.denominator(),
    )
    .ok()?;
    (0..4).find(|turn| Rat::new(i128::from(*turn), 4).ok() == Some(residue))
}

/// Apply an exact root-of-unity phase to a frozen Gaussian-integer character.
fn apply_phase(value: [i32; 2], phase: Rat) -> Option<[i32; 2]> {
    if value == [0, 0] {
        return Some(value);
    }
    match quarter_turn(phase)? {
        0 => Some(value),
        1 => Some([value[1].checked_neg()?, value[0]]),
        2 => Some([value[0].checked_neg()?, value[1].checked_neg()?]),
        3 => Some([value[1], value[0].checked_neg()?]),
        _ => None,
    }
}

fn fingerprint_matches(wanted: &Fingerprint, actual: &Fingerprint) -> bool {
    wanted.len() == actual.len()
        && wanted.iter().all(|(key, phased)| {
            actual.get(key).is_some_and(|candidate| {
                apply_phase(phased.value, phased.phase) == Some(candidate.value)
            })
        })
}

/// Whether `exp(2 pi i shift . t_g)` is a one-dimensional character on the
/// frozen little-group operations.
///
/// If `g h` is represented by `(R_g R_h, t_g + R_g t_h)` up to a parent
/// lattice vector, multiplicativity reduces to
/// `((I - R_g^T) shift) . t_h` being integral. The lattice-vector remainder is
/// integral because `shift` has already been checked against the parent's
/// reciprocal lattice. This includes cases where the shift is not fixed as an
/// exact vector but its residual pairs integrally with every representative
/// translation.
fn shift_twist_is_character(
    table: &LittleCharacterTable,
    shift: &Vec3R,
) -> Result<bool, SubductionError> {
    let transposed_images = table
        .operations
        .iter()
        .map(|operation| {
            let rotation = operation.rotation.map(|row| row.map(i32::from));
            Mat3R::from_ints(rotation)
                .transpose()
                .checked_mul_vector(shift)
        })
        .collect::<Result<Vec<_>, _>>()?;

    for image in transposed_images {
        let residual = shift.checked_sub(&image)?;
        for other in table.operations {
            let Some(translation) = operation_translation(other) else {
                return Ok(false);
            };
            let mut pairing = Rat::ZERO;
            for axis in 0..3 {
                pairing =
                    pairing.checked_add(residual.get(axis).checked_mul(translation.get(axis))?)?;
            }
            if !pairing.is_integer() {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn match_fingerprint(
    parent: u8,
    wanted: &Fingerprint,
    eligible: Option<&BTreeSet<&'static str>>,
) -> LabelImage {
    let mut found: Vec<&'static str> = Vec::new();
    for candidate in line_sources(parent) {
        if eligible.is_some_and(|eligible| !eligible.contains(candidate.label)) {
            continue;
        }
        let Ok(Some(actual)) = fingerprint(candidate, None, false) else {
            continue;
        };
        if fingerprint_matches(wanted, &actual) {
            found.push(candidate.label);
        }
    }
    match found.len() {
        0 => LabelImage::Missing,
        1 => LabelImage::Unique(found[0]),
        _ => LabelImage::Ambiguous(found),
    }
}

/// Build the monodromy map `M_K` of `parent` for a reciprocal shift `K`.
///
/// `K` is expressed in the parent's conventional reciprocal basis. The exact
/// reciprocal-lattice membership is checked here. Since one parent may carry
/// several line little groups, each label is independently marked
/// [`LabelImage::UnsupportedShift`] when the exact reciprocal-shift phase does
/// not define a one-dimensional character on that source's little group.
pub fn monodromy(parent: u8, shift: &Vec3R) -> Result<LineMonodromy, SubductionError> {
    let direct = Lattice::new(exact_primitive_basis(parent)?)?;
    if !direct.reciprocal()?.contains(shift)? {
        return Err(SubductionError::NonReciprocalShift {
            sg: parent,
            shift: *shift,
        });
    }
    let sources = line_sources(parent);
    let mut eligible = BTreeSet::new();
    for source in &sources {
        if shift_twist_is_character(source, shift)? {
            eligible.insert(source.label);
        }
    }
    let mut images = BTreeMap::new();
    for source in sources {
        let image = if !eligible.contains(source.label) {
            LabelImage::UnsupportedShift
        } else {
            match fingerprint(source, Some(shift), false)? {
                Some(shifted) => match_fingerprint(parent, &shifted, Some(&eligible)),
                // A frozen table whose own operations cannot be parsed has no
                // computable image; reporting `Missing` keeps the caller from
                // marching through a guess.
                None => LabelImage::Missing,
            }
        };
        images.insert(source.label, image);
    }
    Ok(LineMonodromy {
        parent,
        shift: *shift,
        conjugate: false,
        images,
    })
}

/// The map induced by complex conjugation of the character: `chi -> conj(chi)`.
///
/// For SG 209 and 210 this exchanges `DT3 <-> DT4`; it is the same as SG 210's
/// along-line shift, but SG 209's along-line shift is the identity. Keeping the
/// maps distinct lets the contract test check `C . M_K = M_-K . C` rather than
/// assuming they are the same operation.
pub fn complex_conjugation(parent: u8) -> LineMonodromy {
    let mut images = BTreeMap::new();
    for source in line_sources(parent) {
        let image = match fingerprint(source, None, true) {
            Ok(Some(conjugated)) => match_fingerprint(parent, &conjugated, None),
            Ok(None) => LabelImage::Missing,
            Err(_) => LabelImage::Missing,
        };
        images.insert(source.label, image);
    }
    LineMonodromy {
        parent,
        shift: Vec3R::zero(),
        conjugate: true,
        images,
    }
}

#[cfg(test)]
mod tests {
    use super::{apply_phase, parse_fraction, phase_turns};
    use crate::irrep::subduction::{Rat, Vec3R};

    #[test]
    fn gaussian_phase_turns_keep_the_positive_bloch_sign() {
        let translation = Vec3R::from_ints([1, 0, 0]);
        let positive_quarter = Vec3R::new([parse_fraction("1/4").unwrap(), Rat::ZERO, Rat::ZERO]);
        let negative_quarter = Vec3R::new([parse_fraction("-1/4").unwrap(), Rat::ZERO, Rat::ZERO]);

        assert_eq!(
            phase_turns(&positive_quarter, &translation).unwrap(),
            parse_fraction("1/4").unwrap()
        );
        assert_eq!(
            phase_turns(&negative_quarter, &translation).unwrap(),
            parse_fraction("-1/4").unwrap()
        );
        assert_eq!(
            apply_phase([1, 0], parse_fraction("1/4").unwrap()),
            Some([0, 1])
        );
        assert_eq!(
            apply_phase([1, 0], parse_fraction("-1/4").unwrap()),
            Some([0, -1])
        );
        assert_eq!(
            apply_phase([1, 0], parse_fraction("5/4").unwrap()),
            Some([0, 1])
        );

        // The zero character is unchanged by any phase; nonzero values outside
        // the frozen Gaussian units cannot be represented by this table format.
        assert_eq!(
            apply_phase([0, 0], parse_fraction("1/3").unwrap()),
            Some([0, 0])
        );
        assert_eq!(apply_phase([1, 0], parse_fraction("1/3").unwrap()), None);
    }
}
