//! Generated per-ordinal embedding settings for the full-subduction engine.
//!
//! `DO NOT EDIT`: regenerate with
//! `python3 scripts/generate_subduction_settings.py --write`; verify with
//! `python3 scripts/generate_subduction_settings.py --check`.
//!
//! Scope: 69 task-8 records plus six task-9 label/shear witnesses,
//! keyed by their **0-based raw isotropy ordinal** (the index
//! into `ISOTROPY_SUBGROUPS`).  Other records keep the engine's candidate
//! search; this table deliberately does not grow beyond that scope.
//!
//! Provenance, re-derived by every `--check` run:
//! * pinned archive `isotropy_subgroup/iso.zip`
//!   (SHA-256 `568667bfc8027095537d642297b319c872d00016b868143c666f90d5931d9f7b`); the machine isotropy tables are read
//!   from that ZIP, never from an extracted checkout. The binary and its
//!   data files are also extracted from this ZIP into a private directory;
//! * official `iso` 9.6.1 oracle, one isolated process per record and ITA
//!   setting, commands `SET I ALL OR 1|2`, `VALUE PARENT/IRREP/DIRECTION`,
//!   `SHOW SUBGROUP/DIRECTION/BASIS/ORIGIN/SIZE/ELEMENTS`, `DISPLAY ISOTROPY`;
//! * `setting` is `U = W . P_parent . (P_sub . B_oracle)^-1`, an exact integer
//!   matrix with `|det| = 1` (no signed-permutation restriction); the stored
//!   origin converted through `P_parent` must equal the printed origin exactly;
//! * `child_shift` is the **origin difference** between the recorded ITA
//!   setting (`SET I ALL OR 1`, the setting the tables are recorded in) and the
//!   printed setting that the shipped Hall row provably belongs to, mapped into
//!   the child cell.  The Hall row is `SG_DATA_HALL[subgroup]` in the
//!   digest-pinned `scripts/hall_operations.json`
//!   (SHA-256 `ebd1cf36668fb8c0efd633b2d7728c51ca1b404a3cc02ed871ece47b46a0d1c8`), the frame the engine loads;
//!   the setting is selected by exact equality of the full expanded operation
//!   set (never fitted to one matching operation) and a row matching no setting
//!   fails the generator.

use crate::mathfunc::Mat3I;

/// Identity setting transform.
const IDENTITY: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// Cyclic axis permutation used by the rhombohedral and C-centred records.
const CYCLIC: Mat3I = [[0, 0, 1], [1, 0, 0], [0, 1, 0]];

/// No child origin shift: the isotropy record already uses the Hall frame.
const NO_SHIFT: [i32; 4] = [0, 0, 0, 1];

/// `(ordinal, parent SG, subgroup SG, U numerator, U denominator, child shift)`.
///
/// The convention is the exact rational matrix `U numerator / U denominator`;
/// the denominator is one for every entry of this table.
///
/// `child_shift` is applied in the **subgroup's own conventional cell**: an
/// operation `(R | t)` of the shipped Hall row becomes
/// `(R | t + delta - R delta)` before the affine map into the parent.  It is
/// non-zero only where the isotropy table and the Hall row use different ITA
/// origin choices (#126).
pub type FrozenEmbeddingSetting = (usize, u8, u8, Mat3I, i32, [i32; 4]);

/// One entry per covered isotropy record ordinal.
pub static FROZEN_EMBEDDING_SETTINGS: &[FrozenEmbeddingSetting] = &[
    // 16 R1 P1 -> #22 (size 2)
    (314, 16, 22, IDENTITY, 1, NO_SHIFT),
    // 16 R2 P1 -> #22 (size 2)
    (315, 16, 22, IDENTITY, 1, NO_SHIFT),
    // 16 R3 P1 -> #22 (size 2)
    (316, 16, 22, IDENTITY, 1, NO_SHIFT),
    // 16 R4 P1 -> #22 (size 2)
    (317, 16, 22, IDENTITY, 1, NO_SHIFT),
    // 19 GM1 P1 -> #19 (size 1)
    (370, 19, 19, IDENTITY, 1, NO_SHIFT),
    // 21 R1 P1 -> #22 (size 4)
    (427, 21, 22, [[1, -1, 0], [1, 0, -1], [1, 0, 0]], 1, NO_SHIFT),
    // 23 GM1 P1 -> #23 (size 1)
    (474, 23, 23, IDENTITY, 1, NO_SHIFT),
    // 43 GM1 P1 -> #43 (size 1)
    (984, 43, 43, IDENTITY, 1, NO_SHIFT),
    // 45 GM1 P1 -> #45 (size 1)
    (1037, 45, 45, IDENTITY, 1, NO_SHIFT),
    // 64 GM2+ P1 -> #14 (size 1)
    (1942, 64, 14, [[2, 0, 1], [-1, 0, -1], [0, 1, 0]], 1, NO_SHIFT),
    // 67 GM2+ P1 -> #13 (size 1)
    (2102, 67, 13, [[2, 0, 1], [-1, 0, -1], [0, 1, 0]], 1, NO_SHIFT),
    // 83 GM1+ P1 -> #83 (size 1)
    (2793, 83, 83, IDENTITY, 1, NO_SHIFT),
    // 83 M1+ P1 -> #83 (size 2)
    (2799, 83, 83, IDENTITY, 1, NO_SHIFT),
    // 83 M2+ P1 -> #83 (size 2)
    (2800, 83, 83, IDENTITY, 1, NO_SHIFT),
    // 83 X1+ P1 -> #83 (size 4)
    (2806, 83, 83, IDENTITY, 1, NO_SHIFT),
    // 83 X2- P1 -> #83 (size 4)
    (2815, 83, 83, IDENTITY, 1, NO_SHIFT),
    // 83 Z1+ P1 -> #83 (size 2)
    (2835, 83, 83, IDENTITY, 1, NO_SHIFT),
    // 83 Z1- P1 -> #83 (size 2)
    (2838, 83, 83, IDENTITY, 1, NO_SHIFT),
    // 139 M1- P1 -> #126 (size 2), delta = ITA origin choice 1 vs 2
    (6213, 139, 126, IDENTITY, 1, [1, 1, 1, 4]),
    // 167 GM3+ P1 -> #15 (size 1)
    (7651, 167, 15, CYCLIC, 1, NO_SHIFT),
    // 167 F1+ C2 -> #15 (size 4)
    (7660, 167, 15, CYCLIC, 1, NO_SHIFT),
    // 167 F2+ P2 -> #15 (size 4)
    (7664, 167, 15, CYCLIC, 1, NO_SHIFT),
    // 167 F1- P2 -> #15 (size 4)
    (7668, 167, 15, CYCLIC, 1, NO_SHIFT),
    // 167 F2- P2 -> #15 (size 4)
    (7674, 167, 15, CYCLIC, 1, NO_SHIFT),
    // 167 L1 P4 -> #15 (size 4)
    (7681, 167, 15, CYCLIC, 1, NO_SHIFT),
    // 167 L1 P5 -> #15 (size 4)
    (7682, 167, 15, CYCLIC, 1, NO_SHIFT),
    // 221 GM3+ P1 -> #123 (size 1)
    (12397, 221, 123, IDENTITY, 1, NO_SHIFT),
    // 221 GM3+ C1 -> #47 (size 1)
    (12398, 221, 47, IDENTITY, 1, NO_SHIFT),
    // 221 GM4+ P3 -> #148 (size 1)
    (12399, 221, 148, IDENTITY, 1, NO_SHIFT),
    // 221 GM4+ P1 -> #83 (size 1)
    (12400, 221, 83, IDENTITY, 1, NO_SHIFT),
    // 221 GM4+ P2 -> #12 (size 1)
    (12401, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 GM5+ C2 -> #12 (size 1)
    (12405, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 R4+ C1 -> #12 (size 2)
    (12433, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 R5+ C2 -> #12 (size 2)
    (12438, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 R5+ C1 -> #12 (size 2)
    (12439, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 R4- C2 -> #12 (size 2)
    (12449, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 R4- C1 -> #12 (size 2)
    (12450, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 R5- C1 -> #12 (size 2)
    (12456, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 X1+ P1 -> #123 (size 2)
    (12458, 221, 123, IDENTITY, 1, NO_SHIFT),
    // 221 X1+ P2 -> #123 (size 4)
    (12459, 221, 123, IDENTITY, 1, NO_SHIFT),
    // 221 X1+ C1 -> #47 (size 4)
    (12460, 221, 47, IDENTITY, 1, NO_SHIFT),
    // 221 X1+ C2 -> #123 (size 8)
    (12462, 221, 123, IDENTITY, 1, NO_SHIFT),
    // 221 X1+ S1 -> #47 (size 8)
    (12463, 221, 47, IDENTITY, 1, NO_SHIFT),
    // 221 X2+ P2 -> #123 (size 4)
    (12465, 221, 123, IDENTITY, 1, NO_SHIFT),
    // 221 X2+ C1 -> #47 (size 4)
    (12466, 221, 47, IDENTITY, 1, NO_SHIFT),
    // 221 X2+ S1 -> #47 (size 8)
    (12469, 221, 47, IDENTITY, 1, NO_SHIFT),
    // 221 X5+ C15 -> #12 (size 4)
    (12490, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 X5+ C23 -> #148 (size 8)
    (12497, 221, 148, IDENTITY, 1, NO_SHIFT),
    // 221 X5+ S13 -> #12 (size 8)
    (12503, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 X3- P1 -> #123 (size 2)
    (12520, 221, 123, IDENTITY, 1, NO_SHIFT),
    // 221 X3- P2 -> #123 (size 4)
    (12521, 221, 123, IDENTITY, 1, NO_SHIFT),
    // 221 X3- C1 -> #47 (size 4)
    (12522, 221, 47, IDENTITY, 1, NO_SHIFT),
    // 221 X3- C2 -> #123 (size 8)
    (12524, 221, 123, IDENTITY, 1, NO_SHIFT),
    // 221 X3- S1 -> #47 (size 8)
    (12525, 221, 47, IDENTITY, 1, NO_SHIFT),
    // 221 X4- P2 -> #123 (size 4)
    (12527, 221, 123, IDENTITY, 1, NO_SHIFT),
    // 221 X4- C1 -> #47 (size 4)
    (12528, 221, 47, IDENTITY, 1, NO_SHIFT),
    // 221 X4- S1 -> #47 (size 8)
    (12531, 221, 47, IDENTITY, 1, NO_SHIFT),
    // 221 X5- C15 -> #12 (size 4)
    (12540, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 X5- C23 -> #148 (size 8)
    (12547, 221, 148, IDENTITY, 1, NO_SHIFT),
    // 221 X5- S13 -> #12 (size 8)
    (12553, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 M1+ P1 -> #123 (size 2)
    (12558, 221, 123, IDENTITY, 1, NO_SHIFT),
    // 221 M4+ P1 -> #123 (size 2)
    (12570, 221, 123, IDENTITY, 1, NO_SHIFT),
    // 221 M5+ C23 -> #148 (size 4)
    (12581, 221, 148, IDENTITY, 1, NO_SHIFT),
    // 221 M5+ S12 -> #12 (size 4)
    (12585, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 M5- C15 -> #12 (size 4)
    (12626, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 M5- S7 -> #12 (size 4)
    (12629, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 221 M5- S6 -> #12 (size 4)
    (12630, 221, 12, CYCLIC, 1, NO_SHIFT),
    // 225 GM4- C2 -> #8 (size 1)
    (13345, 225, 8, CYCLIC, 1, NO_SHIFT),
    // 225 GM4- C1 -> #8 (size 1)
    (13346, 225, 8, CYCLIC, 1, NO_SHIFT),
    // 225 GM5- C1 -> #8 (size 1)
    (13351, 225, 8, CYCLIC, 1, NO_SHIFT),
    // 225 X5- S13 -> #8 (size 4)
    (13533, 225, 8, CYCLIC, 1, NO_SHIFT),
    // 225 W5 4D22 -> #8 (size 16)
    (13729, 225, 8, CYCLIC, 1, NO_SHIFT),
    // 225 W5 6D23 -> #8 (size 32)
    (13813, 225, 8, CYCLIC, 1, NO_SHIFT),
    // 230 GM4- P2 -> #43 (size 1)
    (15125, 230, 43, IDENTITY, 1, NO_SHIFT),
    // 230 GM5- P2 -> #43 (size 1)
    (15131, 230, 43, IDENTITY, 1, NO_SHIFT),
];
