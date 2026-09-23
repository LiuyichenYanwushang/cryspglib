//! Monodromy of the frozen parametric-k line sources under a reciprocal shift.
//!
//! A parametric-k row lives on a line `k(t) = t * v` whose direction `v` is a
//! **parent reciprocal lattice vector** (for the 73 frozen sources: `(0,2,0)`
//! and `(2,2,0)` in the parent's conventional reciprocal basis, both all-even
//! and therefore in the reciprocal lattice of every F-centred parent).  So `t`
//! and `t + 1` name the same point of the parent's Brillouin zone.
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
//! Shifting the parameter by one multiplies that character by the
//! **reciprocal-shift twist** `Phi_v(R) = exp(2 pi i v.T_R)`, which is a genuine
//! one-dimensional character of the line's little group.  The twisted character
//! is again an irreducible character of the same little group, so it is the
//! frozen table of *another* label:
//!
//! ```text
//! rho_{k+K, M_K(alpha)}  ~  rho_{k, alpha},        K in L*_G.
//! ```
//!
//! `M_K` is the **monodromy map** of the parent.  On a nonsymmorphic line it
//! need not be the identity: SG 209/210/227/228 swap their one-dimensional `DT`
//! and `SM` conjugate pairs.  Two consequences are used by the audit:
//!
//! * the character data of label `alpha` at `t + 1` **is** the character data of
//!   label `M_v(alpha)` at `t`, so a decomposition transported by one parameter
//!   step must be compared against `M_v(alpha)` and not against `alpha` again;
//! * the trivial content is a gauge-invariant number, so it is constant on
//!   monodromy orbits of the *pinned* table -- that is the checkable half of the
//!   contract and it is what [`LineMonodromy::orbit`] feeds.
//!
//! # How the map is computed
//!
//! Never from label names.  The shifted character vector `D_alpha(R) * Phi_K(R)`
//! is matched against every frozen table of the same parent, and the unique
//! match is `M_K(alpha)`.  That covers real sources, the `DT3`/`DT4` swap and any
//! longer cycle with one mechanism.  A label whose image is not unique (two
//! frozen tables sharing a fingerprint) or absent is reported as
//! [`LabelImage::Ambiguous`] / [`LabelImage::Missing`] rather than guessed.
//!
//! The fingerprints are keyed by the exact little-group **operation** (rotation
//! plus fractional translation), so two operations that share a rotation stay
//! distinct and the map is exact rather than tolerance-limited in its keys.
//!
//! See `docs/subduction-conventions.md` section 16 for the contract and the
//! parameter domain in which it has been measured.

use crate::irrep::subduction::{Rat, Vec3R};
use crate::irrep::w_little_characters_data::{
    LittleCharacterTable, LittleOperation, W_LITTLE_CHARACTERS,
};
use std::collections::BTreeMap;

/// Exact key of one frozen little-group operation: the rotation plus the
/// translation as exact `(numerator, denominator)` pairs.
///
/// A rotation alone is *not* a key (two coset representatives of one little
/// group can share it); the translation is what makes the key the operation.
type OperationKey = ([[i8; 3]; 3], [(i128, i128); 3]);

/// A character vector of one frozen table, keyed by operation.
type Fingerprint = BTreeMap<OperationKey, (f64, f64)>;

/// Tolerance of the fingerprint comparison.
///
/// The frozen characters are exact integers in `{-1, 0, 1}` (real and imaginary
/// parts apart), and the twist is a phase computed in `f64`, so an exact match
/// after the phase is reproduced to `1e-9`.  The tolerance is far below the
/// `sqrt(2)` that separates any two distinct characters of these tables.
const FINGERPRINT_TOLERANCE: f64 = 1e-9;

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
}

impl LabelImage {
    /// The image when it is uniquely determined.
    pub fn unique(&self) -> Option<&'static str> {
        match self {
            Self::Unique(label) => Some(label),
            Self::Ambiguous(_) | Self::Missing => None,
        }
    }

    /// Whether the frozen data determines this image.
    pub const fn is_unique(&self) -> bool {
        matches!(self, Self::Unique(_))
    }
}

/// The monodromy map `M_K` of one parent, one shift.
///
/// Built by [`monodromy`] (shift by a rational multiple of a frozen direction)
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
                LabelImage::Ambiguous(_) | LabelImage::Missing => return None,
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
        *value = (translation.get(axis).numerator(), translation.get(axis).denominator());
    }
    Some((operation.rotation, key))
}

/// `exp(2 pi i shift . translation)`, the reciprocal-shift twist on one
/// operation.
fn twist_phase(shift: &Vec3R, translation: &Vec3R) -> (f64, f64) {
    let mut angle = 0.0f64;
    for axis in 0..3 {
        angle += shift.get(axis).to_f64() * translation.get(axis).to_f64();
    }
    let angle = std::f64::consts::TAU * angle;
    (angle.cos(), angle.sin())
}

/// One frozen table's character vector, optionally shifted by the twist `K` and
/// optionally conjugated.
fn fingerprint(
    table: &LittleCharacterTable,
    shift: Option<&Vec3R>,
    conjugate: bool,
) -> Option<Fingerprint> {
    let mut out = Fingerprint::new();
    for operation in table.operations {
        let key = operation_key(operation)?;
        let (mut re, mut im) = (
            f64::from(operation.character[0]),
            f64::from(operation.character[1]),
        );
        if conjugate {
            im = -im;
        }
        if let Some(shift) = shift {
            let translation = operation_translation(operation)?;
            let (cos, sin) = twist_phase(shift, &translation);
            let (a, b) = (re, im);
            re = a * cos - b * sin;
            im = a * sin + b * cos;
        }
        out.insert(key, (re, im));
    }
    Some(out)
}

/// Match a fingerprint against every frozen table of one parent.
fn match_fingerprint(parent: u8, wanted: &Fingerprint) -> LabelImage {
    let mut found: Vec<&'static str> = Vec::new();
    for candidate in line_sources(parent) {
        let Some(actual) = fingerprint(candidate, None, false) else {
            continue;
        };
        if actual.len() != wanted.len() {
            continue;
        }
        let equal = actual.iter().all(|(key, value)| {
            wanted.get(key).is_some_and(|other| {
                (value.0 - other.0).abs() <= FINGERPRINT_TOLERANCE
                    && (value.1 - other.1).abs() <= FINGERPRINT_TOLERANCE
            })
        });
        if equal {
            found.push(candidate.label);
        }
    }
    match found.len() {
        0 => LabelImage::Missing,
        1 => LabelImage::Unique(found[0]),
        _ => LabelImage::Ambiguous(found),
    }
}

/// Build the monodromy map `M_K` of `parent` for the shift `K`.
///
/// The shift is expressed in the parent's **conventional reciprocal basis**, the
/// frame the frozen `direction` fields live in.  The map is a statement about
/// the frozen tables; whether `K` is a genuine parent reciprocal lattice vector
/// (the premise of the contract) is the caller's check -- see
/// `SubgroupEmbedding::parent_lattice` and `Lattice::contains`.
pub fn monodromy(parent: u8, shift: &Vec3R) -> LineMonodromy {
    let mut images = BTreeMap::new();
    for source in line_sources(parent) {
        let image = match fingerprint(source, Some(shift), false) {
            Some(shifted) => match_fingerprint(parent, &shifted),
            // A frozen table whose own operations cannot be parsed has no
            // computable image; reporting `Missing` keeps the caller from
            // marching through a guess.
            None => LabelImage::Missing,
        };
        images.insert(source.label, image);
    }
    LineMonodromy {
        parent,
        shift: *shift,
        conjugate: false,
        images,
    }
}

/// The map induced by complex conjugation of the character: `chi -> conj(chi)`.
///
/// For the frozen sources of SG 209/210 this is the `DT3 <-> DT4` pair, i.e. the
/// same two-cycle the shift produces; keeping the two maps apart is what lets
/// the contract test check `C . M_K = M_-K . C` instead of assuming it.
pub fn complex_conjugation(parent: u8) -> LineMonodromy {
    let mut images = BTreeMap::new();
    for source in line_sources(parent) {
        let image = match fingerprint(source, None, true) {
            Some(conjugated) => match_fingerprint(parent, &conjugated),
            None => LabelImage::Missing,
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
