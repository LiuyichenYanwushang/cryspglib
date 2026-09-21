//! Explicit CDML and Bradley–Cracknell label conventions.
//!
//! Generated records store two independent label spellings:
//!
//! * [`IrrepRecord::ml`] — the CDML / Miller–Love label (`"GM4+"`), plain
//!   ASCII, available for every record including spinors.
//! * [`IrrepRecord::bc`] — the **raw legacy** Bradley–Cracknell LaTeX string
//!   (`"\\Gamma_{4}^+"`).  This field is not rewritten; new code should treat
//!   it as source data and ask [`IrrepRecord::label`] for a comparable label.
//!
//! The BC field has two honest gaps:
//!
//! * a handful of scalar records store the literal placeholder `"***"`;
//! * all spinor `.bc` values are synthesized from `.ml` by
//!   `scripts/generate_irrep_data.py`, so they are not verified
//!   Bradley–Cracknell mappings.
//!
//! [`IrrepRecord::label`], [`IrrepRecord::labels`] and
//! [`IrrepRecord::k_label_with_convention`] report those gaps as `None`
//! instead of pretending a mapping exists.  CDML remains available for every
//! record.

use super::types::IrrepRecord;

/// Placeholder stored in `IrrepRecord::bc` when no BC label is known.
const MISSING_BC: &str = "***";

/// Which irrep label spelling a query should use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LabelConvention {
    /// CDML / Miller–Love label, taken literally from [`IrrepRecord::ml`]
    /// (the Γ point stays `"GM"`).
    Cdml,
    /// Bradley–Cracknell label, normalized from the stored raw LaTeX
    /// [`IrrepRecord::bc`] (`"\\Gamma_{3}^+"` → `"Γ3+"`).
    ///
    /// The stored field itself is never modified.  Spinor records and rows
    /// holding `"***"` have no BC label and are reported as unavailable.
    Bc,
}

/// Both label spellings of one irrep or one k-point.
///
/// The two spellings are independent source data, never derived from each
/// other: SG 213 CDML `X2` is BC `X1`, while SG 1 CDML `GM1` and `Z1` are both
/// BC `Γ1`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IrrepLabels {
    /// CDML / Miller–Love label, stored literally (`"GM4+"`, `"X2"`).
    pub cdml: &'static str,
    /// Normalized Bradley–Cracknell label (`"Γ4+"`, `"X1"`), or `None` when
    /// the stored `.bc` field is the `"***"` placeholder or the record is a
    /// double-valued (spinor) row whose `.bc` was synthesized from `.ml`.
    pub bc: Option<String>,
}

impl IrrepRecord {
    /// Label of this record under `convention`, or `None` when that convention
    /// has no genuine label for it.
    ///
    /// CDML is always available.  BC is unavailable for spinor records (whose
    /// stored `.bc` is synthesized from `.ml`) and for the `"***"` placeholder
    /// rows; it is never approximated by the CDML label.
    pub fn label(&self, convention: LabelConvention) -> Option<String> {
        match convention {
            LabelConvention::Cdml => Some(self.ml.to_string()),
            LabelConvention::Bc => self.bc_label(),
        }
    }

    /// Full label of this record in both conventions at once.
    ///
    /// CDML is always present; BC is present exactly when the stored `.bc`
    /// field is a genuine Bradley–Cracknell mapping.
    pub fn labels(&self) -> IrrepLabels {
        IrrepLabels {
            cdml: self.ml,
            bc: self.bc_label(),
        }
    }

    /// k-point prefix of this record's label in both conventions at once
    /// (`"GM4+"` → `"GM"` / `"Γ"`, `"X2"` → `"X"` / `"X"`).
    ///
    /// The BC half is `None` on records without a BC label.  The two prefixes
    /// are read independently: a BC prefix is never inferred from the CDML
    /// k-point symbol.
    pub fn k_labels(&self) -> IrrepLabels {
        IrrepLabels {
            cdml: self.k_label(),
            bc: self.k_label_with_convention(LabelConvention::Bc),
        }
    }

    /// k-point prefix of this record's label under `convention`
    /// (`"Γ3+"` → `"Γ"`, `"B1"` → `"B"`), or `None` when the convention has no
    /// label for the record.
    ///
    /// For CDML this is the existing [`IrrepRecord::k_label`] prefix, so the
    /// Γ point stays `"GM"`.
    pub fn k_label_with_convention(&self, convention: LabelConvention) -> Option<String> {
        match convention {
            LabelConvention::Cdml => Some(self.k_label().to_string()),
            LabelConvention::Bc => self
                .bc_label()
                .map(|label| label_prefix(&label).to_string()),
        }
    }

    /// Normalized BC label, or `None` when the stored BC field is not a
    /// genuine mapping (spinor synthesis or the `"***"` placeholder).
    fn bc_label(&self) -> Option<String> {
        if self.spinor || self.bc == MISSING_BC {
            return None;
        }
        Some(normalize_bc_label(self.bc))
    }
}

/// Prefix of a normalized label before its first index.
///
/// `"Γ3+"` → `"Γ"`, `"B1"` → `"B"`, `"H1H2"` → `"H"`.
fn label_prefix(label: &str) -> &str {
    let body = label.trim_end_matches(['+', '-']);
    match body.find(|c: char| c.is_ascii_digit()) {
        Some(end) => &body[..end],
        None => body,
    }
}

/// Normalize a stored Bradley–Cracknell label spelling to the plain form used
/// for comparisons.
///
/// Notation markup is removed: `_`, `^`, `{`, `}`, `$` and whitespace.  The
/// LaTeX k-point commands emitted by the generator (`\Gamma`, `\Delta`,
/// `\Lambda`, `\Sigma`) become their Unicode letters, and Unicode subscript
/// digits / superscript signs become ASCII.  Every other character — including
/// unknown LaTeX commands, signs and compound label components — is preserved
/// verbatim, so `"+"` never collapses into `"-"` and unidentified text never
/// turns into a match.
///
/// The generator writes only the **first** component of a physical compound
/// label as LaTeX; repeated Γ components keep the CDML `GM` spelling, e.g.
/// stored `"\\Gamma_{2}GM3"`.  A repeated `GM` token therefore normalizes to
/// `Γ`, so the stored raw LaTeX and its Unicode compound form are one query:
/// `"\\Gamma_{2}GM3"` and `"Γ2Γ3"` both become `"Γ2Γ3"`.  A **leading** `GM`
/// is the CDML spelling of Γ and is deliberately left alone, so a CDML label
/// never matches a BC query by accident.
pub(crate) fn normalize_bc_label(raw: &str) -> String {
    let mut normalized = String::with_capacity(raw.len());
    let mut characters = raw.chars().peekable();
    while let Some(character) = characters.next() {
        match character {
            '_' | '^' | '{' | '}' | '$' => {}
            '\\' => {
                let mut command = String::new();
                while let Some(&next) = characters.peek() {
                    if next.is_ascii_alphabetic() {
                        command.push(next);
                        characters.next();
                    } else {
                        break;
                    }
                }
                match latex_kpoint_command(&command) {
                    Some(letter) => normalized.push(letter),
                    // Not notation markup we understand: keep it verbatim
                    // rather than dropping unidentified text.
                    None => {
                        normalized.push('\\');
                        normalized.push_str(&command);
                    }
                }
            }
            // `GM` repeating a Γ k-point inside a compound BC label.
            'G' if !normalized.is_empty() && characters.peek() == Some(&'M') => {
                characters.next();
                normalized.push('Γ');
            }
            '₀'..='₉' => normalized.push(subscript_digit(character)),
            '⁺' => normalized.push('+'),
            '⁻' => normalized.push('-'),
            _ if character.is_whitespace() => {}
            _ => normalized.push(character),
        }
    }
    normalized
}

/// LaTeX k-point commands used by the generator's `KPOINT_LATEX` map.
fn latex_kpoint_command(command: &str) -> Option<char> {
    match command {
        "Gamma" => Some('Γ'),
        "Delta" => Some('Δ'),
        "Lambda" => Some('Λ'),
        "Sigma" => Some('Σ'),
        _ => None,
    }
}

/// Map `'₀'..='₉'` (U+2080..U+2089, contiguous) to `'0'..='9'`.
fn subscript_digit(character: char) -> char {
    let offset = character as u32 - '₀' as u32;
    char::from_u32('0' as u32 + offset).unwrap_or(character)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalization_removes_only_known_markup() {
        assert_eq!(normalize_bc_label("\\Gamma_{3}^+"), "Γ3+");
        assert_eq!(normalize_bc_label("$ \\Gamma_{3}^{+} $"), "Γ3+");
        assert_eq!(normalize_bc_label("Γ₃⁺"), "Γ3+");
        assert_eq!(normalize_bc_label("\\Delta_{10}^-"), "Δ10-");
        assert_eq!(normalize_bc_label("H_{1}H2"), "H1H2");
        assert_eq!(normalize_bc_label("B_{1}"), "B1");
    }

    #[test]
    fn repeated_gamma_tokens_normalize_but_a_leading_gm_stays_cdml() {
        // The generator leaves the second compound component in the CDML
        // spelling; BC must still spell every component as Γ.
        assert_eq!(normalize_bc_label("\\Gamma_{2}GM3"), "Γ2Γ3");
        assert_eq!(normalize_bc_label("\\Gamma_{2}^+GM3+"), "Γ2+Γ3+");
        assert_eq!(normalize_bc_label("Γ2Γ3"), "Γ2Γ3");
        assert_eq!(normalize_bc_label("\\Gamma_{5}^-GM6-"), "Γ5-Γ6-");
        // A leading GM is the CDML spelling of Γ, not a BC label.
        assert_eq!(normalize_bc_label("GM3"), "GM3");
        assert_eq!(normalize_bc_label("GM3GM4"), "GM3Γ4");
    }

    #[test]
    fn normalization_preserves_signs_compounds_and_unknown_text() {
        assert_ne!(
            normalize_bc_label("\\Gamma_{3}^+"),
            normalize_bc_label("\\Gamma_{3}^-")
        );
        assert_eq!(normalize_bc_label("M_{1}M2"), "M1M2");
        assert_eq!(normalize_bc_label("\\Unknown_{1}"), "\\Unknown1");
        assert_eq!(normalize_bc_label("\\UnknownGM_{1}"), "\\UnknownGM1");
        assert_eq!(normalize_bc_label("***"), "***");
    }

    #[test]
    fn prefix_stops_at_the_first_index() {
        assert_eq!(label_prefix("Γ3+"), "Γ");
        assert_eq!(label_prefix("Γ2Γ3"), "Γ");
        assert_eq!(label_prefix("B1"), "B");
        assert_eq!(label_prefix("H1H2"), "H");
        assert_eq!(label_prefix("R1R2"), "R");
    }
}
