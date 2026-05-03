//! # Tetragonal space groups (#75–#142)
//!
//! Tetragonal crystals have a = b ≠ c, α = β = γ = 90°.  The Brillouin
//! zone is a tetragonal prism.
//!
//! ## Common k-point labels (simple tetragonal)
//!
//! | Label | Coords (fractional) | Little group |
//! |-------|---------------------|--------------|
//! | Γ | (0, 0, 0) | 4/mmm (D₄ₕ) |
//! | X | (½, 0, 0) | 4/mmm (D₄ₕ) |
//! | M | (½, ½, 0) | 4/mmm (D₄ₕ) |
//! | Z | (0, 0, ½) | 4/mmm (D₄ₕ) |
//! | R | (0, ½, ½) | 4/mmm (D₄ₕ) |
//! | A | (½, ½, ½) | 4/mmm (D₄ₕ) |
//!
//! ---

/// # 75 P4 (C₄¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3GM4 | $\Gamma_{2}GM4$ | $k_{17}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM4 ($\Gamma_{2}GM4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{2}$ | $k_{15}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{P1} | 75 | P4 | C4^1 | 4 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{P1} | 75 | P4 | C4^1 | 4 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3M4 | $M_{2}M4$ | $k_{18}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M3M4 ($M_{2}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 2 | B8a | no |
/// | R2 | $R_{2}$ | $k_{16}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 79 | I4 | C4^5 | 4 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 79 | I4 | C4^5 | 4 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{20}t1$ | 1 | A2a | no |
/// | A2 | $A_{3}$ | $k_{20}t3$ | 1 | A2a | no |
/// | A3A4 | $A_{2}A4$ | $k_{20}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 79 | I4 | C4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 79 | I4 | C4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for A3A4 ($A_{2}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z3Z4 | $Z_{2}Z4$ | $k_{19}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{2}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
///
pub struct Sg75;

/// # 76 P4₁ (C₄²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3GM4 | $\Gamma_{2}GM4$ | $k_{17}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM4 ($\Gamma_{2}GM4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{2}$ | $k_{15}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{P1} | 76 | P4_1 | C4^2 | 4 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{P1} | 76 | P4_1 | C4^2 | 4 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3M4 | $M_{2}M4$ | $k_{18}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M3M4 ($M_{2}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{16}t1t2$ | 4 | D16b | yes |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A3 | $A_{1}A4$ | $k_{20}t1t2$ | 2 | B8b | yes |
/// | A2A4 | $A_{2}A3$ | $k_{20}t3t4$ | 2 | B8b | yes |
///
/// ### Isotropy subgroups for A1A3 ($A_{1}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ### Isotropy subgroups for A2A4 ($A_{2}A3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z3 | $Z_{1}Z4$ | $k_{19}t1t2$ | 2 | B8b | yes |
/// | Z2Z4 | $Z_{2}Z3$ | $k_{19}t3t4$ | 2 | B8b | yes |
///
/// ### Isotropy subgroups for Z1Z3 ($Z_{1}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2Z4 ($Z_{2}Z3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
pub struct Sg76;

/// # 77 P4₂ (C₄³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3GM4 | $\Gamma_{2}GM4$ | $k_{17}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM4 ($\Gamma_{2}GM4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{2}$ | $k_{15}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{P1} | 77 | P4_2 | C4^3 | 4 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{P1} | 77 | P4_2 | C4^3 | 4 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3M4 | $M_{2}M4$ | $k_{18}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M3M4 ($M_{2}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 2 | B8a | no |
/// | R2 | $R_{2}$ | $k_{16}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 4 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 4 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{2}A4$ | $k_{20}t1t3$ | 2 | B4a | yes |
/// | A3 | $A_{1}$ | $k_{20}t2$ | 1 | A2a | no |
/// | A4 | $A_{3}$ | $k_{20}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for A1A2 ($A_{2}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z2 | $Z_{2}Z4$ | $k_{19}t1t3$ | 2 | B4a | yes |
/// | Z3 | $Z_{1}$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z4 | $Z_{3}$ | $k_{19}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1Z2 ($Z_{2}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 78 | P4_3 | C4^4 | 2 | 1 |
///
pub struct Sg77;

/// # 78 P4₃ (C₄⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3GM4 | $\Gamma_{2}GM4$ | $k_{17}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 78 | P4_3 | C4^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM4 ($\Gamma_{2}GM4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{2}$ | $k_{15}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{P1} | 78 | P4_3 | C4^4 | 4 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{P1} | 78 | P4_3 | C4^4 | 4 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3M4 | $M_{2}M4$ | $k_{18}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 78 | P4_3 | C4^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 78 | P4_3 | C4^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M3M4 ($M_{2}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{16}t1t2$ | 4 | D16b | yes |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A3 | $A_{1}A4$ | $k_{20}t2t3$ | 2 | B8b | yes |
/// | A2A4 | $A_{2}A3$ | $k_{20}t1t4$ | 2 | B8b | yes |
///
/// ### Isotropy subgroups for A1A3 ($A_{1}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ### Isotropy subgroups for A2A4 ($A_{2}A3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z3 | $Z_{1}Z4$ | $k_{19}t2t3$ | 2 | B8b | yes |
/// | Z2Z4 | $Z_{2}Z3$ | $k_{19}t1t4$ | 2 | B8b | yes |
///
/// ### Isotropy subgroups for Z1Z3 ($Z_{1}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2Z4 ($Z_{2}Z3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
pub struct Sg78;

/// # 79 I4 (C₄⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3GM4 | $\Gamma_{2}GM4$ | $k_{14}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 79 | I4 | C4^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM4 ($\Gamma_{2}GM4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | X2 | $X_{2}$ | $k_{13}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 75 | P4 | C4^1 | 4 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 77 | P4_2 | C4^3 | 4 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $Z_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | M2 | $Z_{3}$ | $k_{15}t3$ | 1 | A2a | no |
/// | M3M4 | $Z_{2}Z4$ | $k_{15}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for M1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M3M4 ($Z_{2}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{11}t1$ | 4 | D32b | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 8 | 4 |
/// | dir{P1} | 79 | I4 | C4^5 | 8 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1 | $P_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | P2P2 | $P_{2}P2$ | $k_{12}t2t2$ | 4 | D8a | no |
///
/// ### Isotropy subgroups for P1 ($P_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 4 | 2 |
/// | dir{P1} | 79 | I4 | C4^5 | 4 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for P2P2 ($P_{2}P2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 2 |
///
pub struct Sg79;

/// # 80 I4₁ (C₄⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3GM4 | $\Gamma_{2}GM4$ | $k_{14}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM4 ($\Gamma_{2}GM4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | X2 | $X_{2}$ | $k_{13}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 76 | P4_1 | C4^2 | 4 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 78 | P4_3 | C4^4 | 4 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M2 | $Z_{2}Z4$ | $k_{15}t1t3$ | 2 | B4a | yes |
/// | M3 | $Z_{1}$ | $k_{15}t2$ | 1 | A2a | no |
/// | M4 | $Z_{3}$ | $k_{15}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for M1M2 ($Z_{2}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ### Isotropy subgroups for M3 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 78 | P4_3 | C4^4 | 2 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{11}t1$ | 4 | D32a | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P2 | $P_{1}P2$ | $k_{12}t1t2$ | 4 | D16b | no |
///
/// ### Isotropy subgroups for P1P2 ($P_{1}P2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 2 |
///
pub struct Sg80;

/// # 81 P-4 (S₄¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3GM4 | $\Gamma_{2}GM4$ | $k_{17}t2t4$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM4 ($\Gamma_{2}GM4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{2}$ | $k_{15}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3M4 | $M_{2}M4$ | $k_{18}t2t4$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M3M4 ($M_{2}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 2 | B8a | no |
/// | R2 | $R_{2}$ | $k_{16}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 4 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 4 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{20}t1$ | 1 | A2a | no |
/// | A2 | $A_{3}$ | $k_{20}t3$ | 1 | A2a | no |
/// | A3A4 | $A_{2}A4$ | $k_{20}t2t4$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A3A4 ($A_{2}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z3Z4 | $Z_{2}Z4$ | $k_{19}t2t4$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{2}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
///
pub struct Sg81;

/// # 82 I-4 (S₄²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3GM4 | $\Gamma_{2}GM4$ | $k_{14}t2t4$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM4 ($\Gamma_{2}GM4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | X2 | $X_{2}$ | $k_{13}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $Z_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | M2 | $Z_{3}$ | $k_{15}t3$ | 1 | A2a | no |
/// | M3M4 | $Z_{2}Z4$ | $k_{15}t2t4$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for M1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M3M4 ($Z_{2}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{11}t1$ | 4 | D32b | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 4 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P1 | $P_{1}P1$ | $k_{12}t1t1*$ | 2 | B4a | no |
/// | P2P2 | $P_{3}P3$ | $k_{12}t3t3*$ | 2 | B4a | no |
/// | P3P3 | $P_{2}P2$ | $k_{12}t2t2*$ | 2 | B4a | no |
/// | P4P4 | $P_{4}P4$ | $k_{12}t4t4*$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for P1P1 ($P_{1}P1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 4 | 1 |
///
/// ### Isotropy subgroups for P2P2 ($P_{3}P3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 4 | 1 |
///
/// ### Isotropy subgroups for P3P3 ($P_{2}P2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 4 | 1 |
///
/// ### Isotropy subgroups for P4P4 ($P_{4}P4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 4 | 1 |
///
pub struct Sg82;

/// # 83 P4/m (C₄ₕ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+GM4+ | $\Gamma_{3}^+GM4+$ | $k_{17}t3t7$ | 2 | B4a | no |
/// | GM3-GM4- | $\Gamma_{3}^-GM4-$ | $k_{17}t4t8$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 83 | P4/m | C4h^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+GM4+ ($\Gamma_{3}^+GM4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 4 | 1 |
///
/// ### Isotropy subgroups for GM3-GM4- ($\Gamma_{3}^-GM4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{15}t1$ | 2 | B8a | no |
/// | X1- | $X_{1}^-$ | $k_{15}t2$ | 2 | B8a | no |
/// | X2+ | $X_{2}^+$ | $k_{15}t3$ | 2 | B8a | no |
/// | X2- | $X_{2}^-$ | $k_{15}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{P1} | 83 | P4/m | C4h^1 | 4 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 85 | P4/n | C4h^3 | 4 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X2+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 85 | P4/n | C4h^3 | 4 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X2- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{P1} | 83 | P4/m | C4h^1 | 4 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{18}t1$ | 1 | A2a | no |
/// | M1- | $M_{1}^-$ | $k_{18}t2$ | 1 | A2a | no |
/// | M2+ | $M_{2}^+$ | $k_{18}t5$ | 1 | A2a | no |
/// | M2- | $M_{2}^-$ | $k_{18}t6$ | 1 | A2a | no |
/// | M3+M4+ | $M_{3}^+M4+$ | $k_{18}t3t7$ | 2 | B4a | no |
/// | M3-M4- | $M_{3}^-M4-$ | $k_{18}t4t8$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 83 | P4/m | C4h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 85 | P4/n | C4h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M2+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 83 | P4/m | C4h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M2- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 85 | P4/n | C4h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M3+M4+ ($M_{3}^+M4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
///
/// ### Isotropy subgroups for M3-M4- ($M_{3}^-M4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{16}t1$ | 2 | B8a | no |
/// | R1- | $R_{1}^-$ | $k_{16}t2$ | 2 | B8a | no |
/// | R2+ | $R_{2}^+$ | $k_{16}t3$ | 2 | B8a | no |
/// | R2- | $R_{2}^-$ | $k_{16}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 4 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 2 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 4 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 2 |
///
/// ### Isotropy subgroups for R2+ ($R_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 4 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 2 |
///
/// ### Isotropy subgroups for R2- ($R_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 4 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1+ | $A_{1}^+$ | $k_{20}t1$ | 1 | A2a | no |
/// | A1- | $A_{1}^-$ | $k_{20}t2$ | 1 | A2a | no |
/// | A2+ | $A_{2}^+$ | $k_{20}t5$ | 1 | A2a | no |
/// | A2- | $A_{2}^-$ | $k_{20}t6$ | 1 | A2a | no |
/// | A3+A4+ | $A_{3}^+A4+$ | $k_{20}t3t7$ | 2 | B4a | no |
/// | A3-A4- | $A_{3}^-A4-$ | $k_{20}t4t8$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for A1+ ($A_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 87 | I4/m | C4h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for A1- ($A_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 87 | I4/m | C4h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for A2+ ($A_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 87 | I4/m | C4h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for A2- ($A_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 87 | I4/m | C4h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for A3+A4+ ($A_{3}^+A4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
///
/// ### Isotropy subgroups for A3-A4- ($A_{3}^-A4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{2}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | Z2- | $Z_{2}^-$ | $k_{19}t6$ | 1 | A2a | no |
/// | Z3+Z4+ | $Z_{3}^+Z4+$ | $k_{19}t3t7$ | 2 | B4a | no |
/// | Z3-Z4- | $Z_{3}^-Z4-$ | $k_{19}t4t8$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 83 | P4/m | C4h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 83 | P4/m | C4h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3+Z4+ ($Z_{3}^+Z4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 4 | 1 |
///
/// ### Isotropy subgroups for Z3-Z4- ($Z_{3}^-Z4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 4 | 1 |
///
pub struct Sg83;

/// # 84 P4₂/m (C₄ₕ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+GM4+ | $\Gamma_{3}^+GM4+$ | $k_{17}t3t7$ | 2 | B4a | no |
/// | GM3-GM4- | $\Gamma_{3}^-GM4-$ | $k_{17}t4t8$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+GM4+ ($\Gamma_{3}^+GM4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 4 | 1 |
///
/// ### Isotropy subgroups for GM3-GM4- ($\Gamma_{3}^-GM4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{15}t1$ | 2 | B8a | no |
/// | X1- | $X_{1}^-$ | $k_{15}t2$ | 2 | B8a | no |
/// | X2+ | $X_{2}^+$ | $k_{15}t3$ | 2 | B8a | no |
/// | X2- | $X_{2}^-$ | $k_{15}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 4 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 4 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X2+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 4 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X2- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 4 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{18}t1$ | 1 | A2a | no |
/// | M1- | $M_{1}^-$ | $k_{18}t2$ | 1 | A2a | no |
/// | M2+ | $M_{2}^+$ | $k_{18}t5$ | 1 | A2a | no |
/// | M2- | $M_{2}^-$ | $k_{18}t6$ | 1 | A2a | no |
/// | M3+M4+ | $M_{3}^+M4+$ | $k_{18}t3t7$ | 2 | B4a | no |
/// | M3-M4- | $M_{3}^-M4-$ | $k_{18}t4t8$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M2+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M2- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M3+M4+ ($M_{3}^+M4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
///
/// ### Isotropy subgroups for M3-M4- ($M_{3}^-M4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{16}t1$ | 2 | B8a | no |
/// | R1- | $R_{1}^-$ | $k_{16}t2$ | 2 | B8a | no |
/// | R2+ | $R_{2}^+$ | $k_{16}t3$ | 2 | B8a | no |
/// | R2- | $R_{2}^-$ | $k_{16}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 4 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 2 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 4 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 2 |
///
/// ### Isotropy subgroups for R2+ ($R_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 4 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 2 |
///
/// ### Isotropy subgroups for R2- ($R_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 4 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{2}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | A2 | $A_{1}$ | $k_{20}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 4 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 4 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{2}$ | $k_{19}t1$ | 2 | B8a | yes |
/// | Z2 | $Z_{1}$ | $k_{19}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 4 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg84;

/// # 85 P4/n (C₄ₕ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+GM4+ | $\Gamma_{3}^+GM4+$ | $k_{17}t3t7$ | 2 | B4a | no |
/// | GM3-GM4- | $\Gamma_{3}^-GM4-$ | $k_{17}t4t8$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 85 | P4/n | C4h^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+GM4+ ($\Gamma_{3}^+GM4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 4 | 1 |
///
/// ### Isotropy subgroups for GM3-GM4- ($\Gamma_{3}^-GM4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 8 | 2 |
/// | dir{P1} | 75 | P4 | C4^1 | 8 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{2}$ | $k_{18}t2$ | 2 | B8a | yes |
/// | M2 | $M_{1}$ | $k_{18}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32b | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
/// | dir{P1} | 79 | I4 | C4^5 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{2}$ | $k_{20}t2$ | 2 | B8a | yes |
/// | A2 | $A_{1}$ | $k_{20}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 79 | I4 | C4^5 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{2}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | Z2- | $Z_{2}^-$ | $k_{19}t6$ | 1 | A2a | no |
/// | Z3+Z4+ | $Z_{3}^+Z4+$ | $k_{19}t3t7$ | 2 | B4a | no |
/// | Z3-Z4- | $Z_{3}^-Z4-$ | $k_{19}t4t8$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 85 | P4/n | C4h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 85 | P4/n | C4h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3+Z4+ ($Z_{3}^+Z4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 4 | 1 |
///
/// ### Isotropy subgroups for Z3-Z4- ($Z_{3}^-Z4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 4 | 1 |
///
pub struct Sg85;

/// # 86 P4₂/n (C₄ₕ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+GM4+ | $\Gamma_{3}^+GM4+$ | $k_{17}t3t7$ | 2 | B4a | no |
/// | GM3-GM4- | $\Gamma_{3}^-GM4-$ | $k_{17}t4t8$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+GM4+ ($\Gamma_{3}^+GM4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 4 | 1 |
///
/// ### Isotropy subgroups for GM3-GM4- ($\Gamma_{3}^-GM4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 8 | 2 |
/// | dir{P1} | 77 | P4_2 | C4^3 | 8 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{2}$ | $k_{18}t2$ | 2 | B8a | yes |
/// | M2 | $M_{1}$ | $k_{18}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32b | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1+ | $A_{1}^+$ | $k_{20}t6$ | 1 | A2a | no |
/// | A1- | $A_{1}^-$ | $k_{20}t5$ | 1 | A2a | no |
/// | A2+ | $A_{2}^+$ | $k_{20}t2$ | 1 | A2a | no |
/// | A2- | $A_{2}^-$ | $k_{20}t1$ | 1 | A2a | no |
/// | A3+A4+ | $A_{3}^+A4+$ | $k_{20}t4t8$ | 2 | B4a | no |
/// | A3-A4- | $A_{3}^-A4-$ | $k_{20}t3t7$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for A1+ ($A_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for A1- ($A_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for A2+ ($A_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for A2- ($A_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for A3+A4+ ($A_{3}^+A4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
///
/// ### Isotropy subgroups for A3-A4- ($A_{3}^-A4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{2}$ | $k_{19}t1$ | 2 | B8a | yes |
/// | Z2 | $Z_{1}$ | $k_{19}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 4 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg86;

/// # 87 I4/m (C₄ₕ⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{14}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{14}t6$ | 1 | A2a | no |
/// | GM3+GM4+ | $\Gamma_{3}^+GM4+$ | $k_{14}t3t7$ | 2 | B4a | no |
/// | GM3-GM4- | $\Gamma_{3}^-GM4-$ | $k_{14}t4t8$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 87 | I4/m | C4h^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 79 | I4 | C4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+GM4+ ($\Gamma_{3}^+GM4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 4 | 1 |
///
/// ### Isotropy subgroups for GM3-GM4- ($\Gamma_{3}^-GM4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{13}t1$ | 2 | B8a | no |
/// | X1- | $X_{1}^-$ | $k_{13}t2$ | 2 | B8a | no |
/// | X2+ | $X_{2}^+$ | $k_{13}t3$ | 2 | B8a | no |
/// | X2- | $X_{2}^-$ | $k_{13}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 83 | P4/m | C4h^1 | 4 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 85 | P4/n | C4h^3 | 4 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X2+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 4 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X2- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 4 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $Z_{1}^+$ | $k_{15}t1$ | 1 | A2a | no |
/// | M1- | $Z_{1}^-$ | $k_{15}t2$ | 1 | A2a | no |
/// | M2+ | $Z_{2}^+$ | $k_{15}t5$ | 1 | A2a | no |
/// | M2- | $Z_{2}^-$ | $k_{15}t6$ | 1 | A2a | no |
/// | M3+M4+ | $Z_{3}^+Z4+$ | $k_{15}t3t7$ | 2 | B4a | no |
/// | M3-M4- | $Z_{3}^-Z4-$ | $k_{15}t4t8$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for M1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 83 | P4/m | C4h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 85 | P4/n | C4h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M2+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M2- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M3+M4+ ($Z_{3}^+Z4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 4 | 1 |
///
/// ### Isotropy subgroups for M3-M4- ($Z_{3}^-Z4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 4 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1+ | $N_{1}^+$ | $k_{11}t1$ | 4 | D32b | no |
/// | N1- | $N_{1}^-$ | $k_{11}t2$ | 4 | D32b | no |
///
/// ### Isotropy subgroups for N1+ ($N_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 8 | 4 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 8 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 4 |
///
/// ### Isotropy subgroups for N1- ($N_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 8 | 4 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 8 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1 | $P_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | P2 | $P_{3}$ | $k_{12}t3$ | 2 | B8a | no |
/// | P3P4 | $P_{2}P4$ | $k_{12}t2t4$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for P1 ($P_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 4 | 2 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 4 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
///
/// ### Isotropy subgroups for P2 ($P_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 4 | 2 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 4 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
///
/// ### Isotropy subgroups for P3P4 ($P_{2}P4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
///
pub struct Sg87;

/// # 88 I4₁/a (C₄ₕ⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{14}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{14}t6$ | 1 | A2a | no |
/// | GM3+GM4+ | $\Gamma_{3}^+GM4+$ | $k_{14}t3t7$ | 2 | B4a | no |
/// | GM3-GM4- | $\Gamma_{3}^-GM4-$ | $k_{14}t4t8$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+GM4+ ($\Gamma_{3}^+GM4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 4 | 1 |
///
/// ### Isotropy subgroups for GM3-GM4- ($\Gamma_{3}^-GM4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 4 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 4 | D32b | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 8 | 2 |
/// | dir{P1} | 78 | P4_3 | C4^4 | 8 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $Z_{2}$ | $k_{15}t1$ | 2 | B8a | yes |
/// | M2 | $Z_{1}$ | $k_{15}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M2 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 4 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1+ | $N_{1}^+$ | $k_{11}t1$ | 4 | D64a | no |
/// | N1- | $N_{1}^-$ | $k_{11}t2$ | 4 | D64a | no |
///
/// ### Isotropy subgroups for N1+ ($N_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 16 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{C1} | 2 | P-1 | Ci^1 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ### Isotropy subgroups for N1- ($N_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 16 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{C1} | 2 | P-1 | Ci^1 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P4 | $P_{1}P4$ | $k_{12}t1t4$ | 4 | D32e | no |
/// | P2P3 | $P_{2}P3$ | $k_{12}t2t3$ | 4 | D32e | no |
///
/// ### Isotropy subgroups for P1P4 ($P_{1}P4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ### Isotropy subgroups for P2P3 ($P_{2}P3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
pub struct Sg88;

/// # 89 P422 (D₄¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{15}t4$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{15}t2$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{15}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{P1} | 89 | P422 | D4^1 | 4 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 90 | P42_12 | D4^2 | 4 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 90 | P42_12 | D4^2 | 4 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{P1} | 89 | P422 | D4^1 | 4 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3 | $M_{2}$ | $k_{18}t2$ | 1 | A2a | no |
/// | M4 | $M_{4}$ | $k_{18}t4$ | 1 | A2a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 2 | B8a | no |
/// | R2 | $R_{3}$ | $k_{16}t4$ | 2 | B8a | no |
/// | R3 | $R_{4}$ | $k_{16}t2$ | 2 | B8a | no |
/// | R4 | $R_{2}$ | $k_{16}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 97 | I422 | D4^9 | 4 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 97 | I422 | D4^9 | 4 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 2 |
///
/// ### Isotropy subgroups for R3 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 97 | I422 | D4^9 | 4 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 2 |
///
/// ### Isotropy subgroups for R4 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 97 | I422 | D4^9 | 4 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{20}t1$ | 1 | A2a | no |
/// | A2 | $A_{3}$ | $k_{20}t3$ | 1 | A2a | no |
/// | A3 | $A_{2}$ | $k_{20}t2$ | 1 | A2a | no |
/// | A4 | $A_{4}$ | $k_{20}t4$ | 1 | A2a | no |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 2 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 4 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z3 | $Z_{2}$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z4 | $Z_{4}$ | $k_{19}t4$ | 1 | A2a | no |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg89;

/// # 90 P42₁2 (D₄²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 75 | P4 | C4^1 | 8 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M4 | $M_{1}M2$ | $k_{18}t2t3$ | 2 | B4a | no |
/// | M2M3 | $M_{3}M4$ | $k_{18}t1t4$ | 2 | B4a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1M4 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 4 | 1 |
///
/// ### Isotropy subgroups for M2M3 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 79 | I4 | C4^5 | 8 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A4 | $A_{1}A2$ | $k_{20}t2t3$ | 2 | B4a | no |
/// | A2A3 | $A_{3}A4$ | $k_{20}t1t4$ | 2 | B4a | no |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1A4 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 4 | 1 |
///
/// ### Isotropy subgroups for A2A3 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 4 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 79 | I4 | C4^5 | 4 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z3 | $Z_{2}$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z4 | $Z_{4}$ | $k_{19}t4$ | 1 | A2a | no |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg90;

/// # 91 P4₁22 (D₄³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{15}t4$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{15}t2$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{15}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 91 | P4_122 | D4^3 | 4 | 2 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 4 | 2 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 4 | 2 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 91 | P4_122 | D4^3 | 4 | 2 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3 | $M_{2}$ | $k_{18}t2$ | 1 | A2a | no |
/// | M4 | $M_{4}$ | $k_{18}t4$ | 1 | A2a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32d | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{20}t1$ | 2 | B16a | yes |
/// | A2 | $A_{2}$ | $k_{20}t2$ | 2 | B16a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t1$ | 2 | B16a | yes |
/// | Z2 | $Z_{2}$ | $k_{19}t2$ | 2 | B16a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
///
pub struct Sg91;

/// # 92 P4₁2₁2 (D₄⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 76 | P4_1 | C4^2 | 8 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M4 | $M_{1}M2$ | $k_{18}t2t3$ | 2 | B4a | no |
/// | M2M3 | $M_{3}M4$ | $k_{18}t1t4$ | 2 | B4a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1M4 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
///
/// ### Isotropy subgroups for M2M3 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R3 | $R_{3}R4$ | $k_{16}t1t2$ | 4 | D32e | no |
/// | R2R4 | $R_{1}R2$ | $k_{16}t3t4$ | 4 | D32e | no |
///
/// ### Isotropy subgroups for R1R3 ($R_{3}R4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ### Isotropy subgroups for R2R4 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{20}t1t2$ | 4 | D16a | no |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t1$ | 2 | B16a | yes |
/// | Z2 | $Z_{2}$ | $k_{19}t2$ | 2 | B16a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
///
pub struct Sg92;

/// # 93 P4₂22 (D₄⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{15}t4$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{15}t2$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{15}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 4 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 4 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 4 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 4 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3 | $M_{2}$ | $k_{18}t2$ | 1 | A2a | no |
/// | M4 | $M_{4}$ | $k_{18}t4$ | 1 | A2a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 2 | B8a | no |
/// | R2 | $R_{3}$ | $k_{16}t4$ | 2 | B8a | no |
/// | R3 | $R_{4}$ | $k_{16}t2$ | 2 | B8a | no |
/// | R4 | $R_{2}$ | $k_{16}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 4 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 4 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 2 |
///
/// ### Isotropy subgroups for R3 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 4 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 2 |
///
/// ### Isotropy subgroups for R4 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 4 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{20}t3$ | 1 | A2a | no |
/// | A2 | $A_{3}$ | $k_{20}t1$ | 1 | A2a | no |
/// | A3 | $A_{2}$ | $k_{20}t4$ | 1 | A2a | no |
/// | A4 | $A_{4}$ | $k_{20}t2$ | 1 | A2a | no |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 2 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 4 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z3 | $Z_{2}$ | $k_{19}t4$ | 1 | A2a | no |
/// | Z4 | $Z_{4}$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
pub struct Sg93;

/// # 94 P4₂2₁2 (D₄⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 77 | P4_2 | C4^3 | 8 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M4 | $M_{1}M2$ | $k_{18}t2t3$ | 2 | B4a | no |
/// | M2M3 | $M_{3}M4$ | $k_{18}t1t4$ | 2 | B4a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1M4 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 4 | 1 |
///
/// ### Isotropy subgroups for M2M3 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 8 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A4 | $A_{1}A2$ | $k_{20}t1t4$ | 2 | B4a | no |
/// | A2A3 | $A_{3}A4$ | $k_{20}t2t3$ | 2 | B4a | no |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1A4 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 4 | 1 |
///
/// ### Isotropy subgroups for A2A3 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 4 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 4 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{2}$ | $k_{19}t4$ | 1 | A2a | no |
/// | Z2 | $Z_{4}$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z3 | $Z_{1}$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z4 | $Z_{3}$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
pub struct Sg94;

/// # 95 P4₃22 (D₄⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 78 | P4_3 | C4^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{15}t4$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{15}t2$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{15}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 95 | P4_322 | D4^7 | 4 | 2 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 4 | 2 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 4 | 2 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 95 | P4_322 | D4^7 | 4 | 2 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3 | $M_{2}$ | $k_{18}t2$ | 1 | A2a | no |
/// | M4 | $M_{4}$ | $k_{18}t4$ | 1 | A2a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32d | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{20}t2$ | 2 | B16a | yes |
/// | A2 | $A_{2}$ | $k_{20}t1$ | 2 | B16a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t2$ | 2 | B16a | yes |
/// | Z2 | $Z_{2}$ | $k_{19}t1$ | 2 | B16a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
///
pub struct Sg95;

/// # 96 P4₃2₁2 (D₄⁸)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 78 | P4_3 | C4^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 78 | P4_3 | C4^4 | 8 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M4 | $M_{1}M2$ | $k_{18}t2t3$ | 2 | B4a | no |
/// | M2M3 | $M_{3}M4$ | $k_{18}t1t4$ | 2 | B4a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1M4 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
///
/// ### Isotropy subgroups for M2M3 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 78 | P4_3 | C4^4 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R3 | $R_{3}R4$ | $k_{16}t3t4$ | 4 | D32e | no |
/// | R2R4 | $R_{1}R2$ | $k_{16}t1t2$ | 4 | D32e | no |
///
/// ### Isotropy subgroups for R1R3 ($R_{3}R4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ### Isotropy subgroups for R2R4 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{20}t1t2$ | 4 | D16a | no |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t2$ | 2 | B16a | yes |
/// | Z2 | $Z_{2}$ | $k_{19}t1$ | 2 | B16a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
///
pub struct Sg96;

/// # 97 I422 (D₄⁹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{14}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 79 | I4 | C4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 22 | F222 | D2^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{13}t4$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{13}t3$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{13}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 89 | P422 | D4^1 | 4 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 90 | P42_12 | D4^2 | 4 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 4 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 4 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $Z_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | M2 | $Z_{3}$ | $k_{15}t3$ | 1 | A2a | no |
/// | M3 | $Z_{2}$ | $k_{15}t2$ | 1 | A2a | no |
/// | M4 | $Z_{4}$ | $k_{15}t4$ | 1 | A2a | no |
/// | M5 | $Z_{5}$ | $k_{15}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for M1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{11}t1$ | 4 | D64c | no |
/// | N2 | $N_{2}$ | $k_{11}t2$ | 4 | D64c | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 8 | 4 |
/// | dir{P1} | 97 | I422 | D4^9 | 8 | 4 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 4 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ### Isotropy subgroups for N2 ($N_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 8 | 4 |
/// | dir{P1} | 97 | I422 | D4^9 | 8 | 4 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 4 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1 | $P_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | P2 | $P_{3}$ | $k_{12}t4$ | 2 | B8a | no |
/// | P3P4 | $P_{2}P4$ | $k_{12}t2t3$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for P1 ($P_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 4 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 4 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 2 |
///
/// ### Isotropy subgroups for P2 ($P_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 4 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 4 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 2 |
///
/// ### Isotropy subgroups for P3P4 ($P_{2}P4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
///
pub struct Sg97;

/// # 98 I4₁22 (D₄¹⁰)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{14}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 22 | F222 | D2^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{13}t4$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{13}t3$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{13}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 91 | P4_122 | D4^3 | 4 | 2 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 4 | 2 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 95 | P4_322 | D4^7 | 4 | 2 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 4 | 2 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $Z_{2}$ | $k_{15}t1$ | 1 | A2a | no |
/// | M2 | $Z_{4}$ | $k_{15}t3$ | 1 | A2a | no |
/// | M3 | $Z_{1}$ | $k_{15}t2$ | 1 | A2a | no |
/// | M4 | $Z_{3}$ | $k_{15}t4$ | 1 | A2a | no |
/// | M5 | $Z_{5}$ | $k_{15}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for M1 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{11}t1$ | 4 | D64d | no |
/// | N2 | $N_{2}$ | $k_{11}t2$ | 4 | D64d | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 16 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ### Isotropy subgroups for N2 ($N_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 16 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1 | $P_{1}$ | $k_{12}t1$ | 4 | D32d | no |
///
/// ### Isotropy subgroups for P1 ($P_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
pub struct Sg98;

/// # 99 P4mm (C₄ᵥ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 99 | P4mm | C4v^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 4 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{15}t2$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{15}t3$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{15}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 4 | 1 |
/// | dir{P1} | 99 | P4mm | C4v^1 | 4 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{P1} | 100 | P4bm | C4v^2 | 4 | 2 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 4 | 1 |
/// | dir{P1} | 99 | P4mm | C4v^1 | 4 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{P1} | 100 | P4bm | C4v^2 | 4 | 2 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3 | $M_{4}$ | $k_{18}t4$ | 1 | A2a | no |
/// | M4 | $M_{2}$ | $k_{18}t2$ | 1 | A2a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 99 | P4mm | C4v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 100 | P4bm | C4v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 99 | P4mm | C4v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 100 | P4bm | C4v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 4 | 1 |
/// | dir{P1} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 2 | B8a | no |
/// | R2 | $R_{3}$ | $k_{16}t2$ | 2 | B8a | no |
/// | R3 | $R_{4}$ | $k_{16}t3$ | 2 | B8a | no |
/// | R4 | $R_{2}$ | $k_{16}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{P1} | 107 | I4mm | C4v^9 | 4 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{P1} | 108 | I4cm | C4v^10 | 4 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 8 | 2 |
///
/// ### Isotropy subgroups for R3 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{P1} | 107 | I4mm | C4v^9 | 4 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 8 | 2 |
///
/// ### Isotropy subgroups for R4 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{P1} | 108 | I4cm | C4v^10 | 4 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 8 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{20}t1$ | 1 | A2a | no |
/// | A2 | $A_{3}$ | $k_{20}t3$ | 1 | A2a | no |
/// | A3 | $A_{4}$ | $k_{20}t4$ | 1 | A2a | no |
/// | A4 | $A_{2}$ | $k_{20}t2$ | 1 | A2a | no |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 107 | I4mm | C4v^9 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 108 | I4cm | C4v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 107 | I4mm | C4v^9 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 108 | I4cm | C4v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 46 | Ima2 | C2v^22 | 4 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z3 | $Z_{4}$ | $k_{19}t4$ | 1 | A2a | no |
/// | Z4 | $Z_{2}$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 99 | P4mm | C4v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 105 | P4_2mc | C4v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 101 | P4_2cm | C4v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 103 | P4cc | C4v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg99;

/// # 100 P4bm (C₄ᵥ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 100 | P4bm | C4v^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 32 | Pba2 | C2v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 75 | P4 | C4^1 | 8 | 2 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 8 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M3 | $M_{3}M4$ | $k_{18}t1t4$ | 2 | B4a | yes |
/// | M2M4 | $M_{1}M2$ | $k_{18}t2t3$ | 2 | B4a | yes |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1M3 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 4 | 1 |
///
/// ### Isotropy subgroups for M2M4 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 32 | Pba2 | C2v^8 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 4 | 1 |
/// | dir{P1} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 79 | I4 | C4^5 | 8 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 8 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A3 | $A_{3}A4$ | $k_{20}t1t4$ | 2 | B4a | yes |
/// | A2A4 | $A_{1}A2$ | $k_{20}t2t3$ | 2 | B4a | yes |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1A3 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 44 | Imm2 | C2v^20 | 4 | 1 |
///
/// ### Isotropy subgroups for A2A4 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 45 | Iba2 | C2v^21 | 4 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 79 | I4 | C4^5 | 4 | 1 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z3 | $Z_{4}$ | $k_{19}t4$ | 1 | A2a | no |
/// | Z4 | $Z_{2}$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 100 | P4bm | C4v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 106 | P4_2bc | C4v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 102 | P4_2nm | C4v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 104 | P4nc | C4v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg100;

/// # 101 P4₂cm (C₄ᵥ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 101 | P4_2cm | C4v^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{15}t2$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{15}t3$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{15}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 4 | 1 |
/// | dir{P1} | 101 | P4_2cm | C4v^3 | 4 | 2 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{P1} | 102 | P4_2nm | C4v^4 | 4 | 2 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 4 | 1 |
/// | dir{P1} | 101 | P4_2cm | C4v^3 | 4 | 2 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{P1} | 102 | P4_2nm | C4v^4 | 4 | 2 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3 | $M_{4}$ | $k_{18}t4$ | 1 | A2a | no |
/// | M4 | $M_{2}$ | $k_{18}t2$ | 1 | A2a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 105 | P4_2mc | C4v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 106 | P4_2bc | C4v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 105 | P4_2mc | C4v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 106 | P4_2bc | C4v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 4 | 1 |
/// | dir{P1} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{16}t1t2$ | 4 | D16c | yes |
/// | R3R4 | $R_{3}R4$ | $k_{16}t3t4$ | 4 | D16c | yes |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 8 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
///
/// ### Isotropy subgroups for R3R4 ($R_{3}R4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 8 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A3 | $A_{1}A2$ | $k_{20}t1t4$ | 2 | B4a | yes |
/// | A2A4 | $A_{3}A4$ | $k_{20}t2t3$ | 2 | B4a | yes |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1A3 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 44 | Imm2 | C2v^20 | 4 | 1 |
///
/// ### Isotropy subgroups for A2A4 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 45 | Iba2 | C2v^21 | 4 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 4 | 1 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z3 | $Z_{1}Z2$ | $k_{19}t1t4$ | 2 | B4a | yes |
/// | Z2Z4 | $Z_{3}Z4$ | $k_{19}t2t3$ | 2 | B4a | yes |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1Z3 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 4 | 1 |
///
/// ### Isotropy subgroups for Z2Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 4 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 78 | P4_3 | C4^4 | 4 | 1 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg101;

/// # 102 P4₂nm (C₄ᵥ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 102 | P4_2nm | C4v^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 77 | P4_2 | C4^3 | 8 | 2 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 8 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M3 | $M_{3}M4$ | $k_{18}t1t4$ | 2 | B4a | yes |
/// | M2M4 | $M_{1}M2$ | $k_{18}t2t3$ | 2 | B4a | yes |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1M3 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 4 | 1 |
///
/// ### Isotropy subgroups for M2M4 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 32 | Pba2 | C2v^8 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 4 | 1 |
/// | dir{P1} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 8 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 8 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{4}$ | $k_{20}t4$ | 1 | A2a | no |
/// | A2 | $A_{1}$ | $k_{20}t2$ | 1 | A2a | no |
/// | A3 | $A_{2}$ | $k_{20}t1$ | 1 | A2a | no |
/// | A4 | $A_{3}$ | $k_{20}t3$ | 1 | A2a | no |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 109 | I4_1md | C4v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 110 | I4_1cd | C4v^12 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 109 | I4_1md | C4v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 110 | I4_1cd | C4v^12 | 2 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 46 | Ima2 | C2v^22 | 4 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z3 | $Z_{1}Z2$ | $k_{19}t1t4$ | 2 | B4a | yes |
/// | Z2Z4 | $Z_{3}Z4$ | $k_{19}t2t3$ | 2 | B4a | yes |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1Z3 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 4 | 1 |
///
/// ### Isotropy subgroups for Z2Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 4 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 78 | P4_3 | C4^4 | 4 | 1 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg102;

/// # 103 P4cc (C₄ᵥ⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 103 | P4cc | C4v^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{15}t2$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{15}t3$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{15}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 4 | 1 |
/// | dir{P1} | 103 | P4cc | C4v^5 | 4 | 2 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{P1} | 104 | P4nc | C4v^6 | 4 | 2 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 4 | 1 |
/// | dir{P1} | 103 | P4cc | C4v^5 | 4 | 2 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{P1} | 104 | P4nc | C4v^6 | 4 | 2 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3 | $M_{4}$ | $k_{18}t4$ | 1 | A2a | no |
/// | M4 | $M_{2}$ | $k_{18}t2$ | 1 | A2a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 103 | P4cc | C4v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 104 | P4nc | C4v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 103 | P4cc | C4v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 104 | P4nc | C4v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 4 | 1 |
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{16}t1t2$ | 4 | D16c | yes |
/// | R3R4 | $R_{3}R4$ | $k_{16}t3t4$ | 4 | D16c | yes |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 79 | I4 | C4^5 | 8 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
///
/// ### Isotropy subgroups for R3R4 ($R_{3}R4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 79 | I4 | C4^5 | 8 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A4 | $A_{3}A4$ | $k_{20}t1t2$ | 2 | B4a | yes |
/// | A2A3 | $A_{1}A2$ | $k_{20}t3t4$ | 2 | B4a | yes |
/// | A5A5 | $A_{5}A5$ | $k_{20}t5t5$ | 4 | D8a | no |
///
/// ### Isotropy subgroups for A1A4 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 79 | I4 | C4^5 | 4 | 1 |
///
/// ### Isotropy subgroups for A2A3 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 79 | I4 | C4^5 | 4 | 1 |
///
/// ### Isotropy subgroups for A5A5 ($A_{5}A5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z4 | $Z_{3}Z4$ | $k_{19}t1t2$ | 2 | B4a | yes |
/// | Z2Z3 | $Z_{1}Z2$ | $k_{19}t3t4$ | 2 | B4a | yes |
/// | Z5Z5 | $Z_{5}Z5$ | $k_{19}t5t5$ | 4 | D8a | no |
///
/// ### Isotropy subgroups for Z1Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 4 | 1 |
///
/// ### Isotropy subgroups for Z2Z3 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 4 | 1 |
///
/// ### Isotropy subgroups for Z5Z5 ($Z_{5}Z5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg103;

/// # 104 P4nc (C₄ᵥ⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 104 | P4nc | C4v^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 75 | P4 | C4^1 | 8 | 2 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 8 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M3 | $M_{3}M4$ | $k_{18}t1t4$ | 2 | B4a | yes |
/// | M2M4 | $M_{1}M2$ | $k_{18}t2t3$ | 2 | B4a | yes |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1M3 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 4 | 1 |
///
/// ### Isotropy subgroups for M2M4 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 4 | 1 |
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 79 | I4 | C4^5 | 8 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 8 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{20}t2t4$ | 2 | B4a | yes |
/// | A3A4 | $A_{3}A4$ | $k_{20}t1t3$ | 2 | B4a | yes |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 4 | 1 |
///
/// ### Isotropy subgroups for A3A4 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 4 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 79 | I4 | C4^5 | 4 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z4 | $Z_{3}Z4$ | $k_{19}t1t2$ | 2 | B4a | yes |
/// | Z2Z3 | $Z_{1}Z2$ | $k_{19}t3t4$ | 2 | B4a | yes |
/// | Z5Z5 | $Z_{5}Z5$ | $k_{19}t5t5$ | 4 | D8a | no |
///
/// ### Isotropy subgroups for Z1Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 75 | P4 | C4^1 | 4 | 1 |
///
/// ### Isotropy subgroups for Z2Z3 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 4 | 1 |
///
/// ### Isotropy subgroups for Z5Z5 ($Z_{5}Z5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg104;

/// # 105 P4₂mc (C₄ᵥ⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 105 | P4_2mc | C4v^7 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 4 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{15}t2$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{15}t3$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{15}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 4 | 1 |
/// | dir{P1} | 105 | P4_2mc | C4v^7 | 4 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{P1} | 106 | P4_2bc | C4v^8 | 4 | 2 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 4 | 1 |
/// | dir{P1} | 105 | P4_2mc | C4v^7 | 4 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{P1} | 106 | P4_2bc | C4v^8 | 4 | 2 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3 | $M_{4}$ | $k_{18}t4$ | 1 | A2a | no |
/// | M4 | $M_{2}$ | $k_{18}t2$ | 1 | A2a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 101 | P4_2cm | C4v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 102 | P4_2nm | C4v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 101 | P4_2cm | C4v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 102 | P4_2nm | C4v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 4 | 1 |
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 2 | B8a | no |
/// | R2 | $R_{3}$ | $k_{16}t2$ | 2 | B8a | no |
/// | R3 | $R_{4}$ | $k_{16}t3$ | 2 | B8a | no |
/// | R4 | $R_{2}$ | $k_{16}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{P1} | 109 | I4_1md | C4v^11 | 4 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{P1} | 110 | I4_1cd | C4v^12 | 4 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 8 | 2 |
///
/// ### Isotropy subgroups for R3 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{P1} | 109 | I4_1md | C4v^11 | 4 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 8 | 2 |
///
/// ### Isotropy subgroups for R4 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{P1} | 110 | I4_1cd | C4v^12 | 4 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 8 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{20}t1t3$ | 2 | B4a | yes |
/// | A3A4 | $A_{3}A4$ | $k_{20}t2t4$ | 2 | B4a | yes |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 4 | 1 |
///
/// ### Isotropy subgroups for A3A4 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 4 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 4 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z2 | $Z_{1}Z2$ | $k_{19}t1t3$ | 2 | B4a | yes |
/// | Z3Z4 | $Z_{3}Z4$ | $k_{19}t2t4$ | 2 | B4a | yes |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1Z2 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 4 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 4 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 4 | 1 |
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg105;

/// # 106 P4₂bc (C₄ᵥ⁸)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 106 | P4_2bc | C4v^8 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 32 | Pba2 | C2v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 77 | P4_2 | C4^3 | 8 | 2 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 8 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M3 | $M_{3}M4$ | $k_{18}t1t4$ | 2 | B4a | yes |
/// | M2M4 | $M_{1}M2$ | $k_{18}t2t3$ | 2 | B4a | yes |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1M3 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 4 | 1 |
///
/// ### Isotropy subgroups for M2M4 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 77 | P4_2 | C4^3 | 4 | 1 |
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 8 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 8 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A4 | $A_{1}A2$ | $k_{20}t1t2$ | 2 | B4a | yes |
/// | A2A3 | $A_{3}A4$ | $k_{20}t3t4$ | 2 | B4a | yes |
/// | A5A5 | $A_{5}A5$ | $k_{20}t5t5$ | 4 | D8a | no |
///
/// ### Isotropy subgroups for A1A4 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 4 | 1 |
///
/// ### Isotropy subgroups for A2A3 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 4 | 1 |
///
/// ### Isotropy subgroups for A5A5 ($A_{5}A5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z2 | $Z_{1}Z2$ | $k_{19}t1t3$ | 2 | B4a | yes |
/// | Z3Z4 | $Z_{3}Z4$ | $k_{19}t2t4$ | 2 | B4a | yes |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1Z2 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 32 | Pba2 | C2v^8 | 4 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 4 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 4 | 1 |
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg106;

/// # 107 I4mm (C₄ᵥ⁹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{14}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 107 | I4mm | C4v^9 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 44 | Imm2 | C2v^20 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 79 | I4 | C4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{13}t2$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{13}t4$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{13}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{P1} | 99 | P4mm | C4v^1 | 4 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 41 | Aba2 | C2v^17 | 4 | 1 |
/// | dir{P1} | 100 | P4bm | C4v^2 | 4 | 2 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 40 | Ama2 | C2v^16 | 4 | 1 |
/// | dir{P1} | 102 | P4_2nm | C4v^4 | 4 | 2 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{P1} | 101 | P4_2cm | C4v^3 | 4 | 2 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $Z_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | M2 | $Z_{3}$ | $k_{15}t3$ | 1 | A2a | no |
/// | M3 | $Z_{4}$ | $k_{15}t4$ | 1 | A2a | no |
/// | M4 | $Z_{2}$ | $k_{15}t2$ | 1 | A2a | no |
/// | M5 | $Z_{5}$ | $k_{15}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 99 | P4mm | C4v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 105 | P4_2mc | C4v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 102 | P4_2nm | C4v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 104 | P4nc | C4v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{11}t1$ | 4 | D64c | no |
/// | N2 | $N_{2}$ | $k_{11}t2$ | 4 | D64c | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 8 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 8 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 109 | I4_1md | C4v^11 | 8 | 4 |
/// | dir{P1} | 107 | I4mm | C4v^9 | 8 | 4 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 16 | 4 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 16 | 4 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 4 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 4 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ### Isotropy subgroups for N2 ($N_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 8 | 1 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 8 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 110 | I4_1cd | C4v^12 | 8 | 4 |
/// | dir{P1} | 108 | I4cm | C4v^10 | 8 | 4 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 16 | 4 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 16 | 4 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 4 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 4 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1 | $P_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | P2 | $P_{3}$ | $k_{12}t2$ | 2 | B8a | no |
/// | P3P4 | $P_{2}P4$ | $k_{12}t3t4$ | 4 | D16c | yes |
///
/// ### Isotropy subgroups for P1 ($P_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 109 | I4_1md | C4v^11 | 4 | 2 |
/// | dir{P1} | 107 | I4mm | C4v^9 | 4 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 8 | 2 |
///
/// ### Isotropy subgroups for P2 ($P_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 110 | I4_1cd | C4v^12 | 4 | 2 |
/// | dir{P1} | 108 | I4cm | C4v^10 | 4 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 8 | 2 |
///
/// ### Isotropy subgroups for P3P4 ($P_{2}P4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 46 | Ima2 | C2v^22 | 8 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 8 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
///
pub struct Sg107;

/// # 108 I4cm (C₄ᵥ¹⁰)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{14}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 108 | I4cm | C4v^10 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 45 | Iba2 | C2v^21 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 79 | I4 | C4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{13}t2$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{13}t4$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{13}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{P1} | 103 | P4cc | C4v^5 | 4 | 2 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 40 | Ama2 | C2v^16 | 4 | 1 |
/// | dir{P1} | 104 | P4nc | C4v^6 | 4 | 2 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 41 | Aba2 | C2v^17 | 4 | 1 |
/// | dir{P1} | 106 | P4_2bc | C4v^8 | 4 | 2 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{P1} | 105 | P4_2mc | C4v^7 | 4 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $Z_{1}$ | $k_{15}t2$ | 1 | A2a | no |
/// | M2 | $Z_{3}$ | $k_{15}t4$ | 1 | A2a | no |
/// | M3 | $Z_{4}$ | $k_{15}t3$ | 1 | A2a | no |
/// | M4 | $Z_{2}$ | $k_{15}t1$ | 1 | A2a | no |
/// | M5 | $Z_{5}$ | $k_{15}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 103 | P4cc | C4v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 101 | P4_2cm | C4v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 106 | P4_2bc | C4v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 100 | P4bm | C4v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1N2 | $N_{1}N2$ | $k_{11}t1t2$ | 8 | F64a | no |
///
/// ### Isotropy subgroups for N1N2 ($N_{1}N2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P3} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P3} | 80 | I4_1 | C4^6 | 16 | 4 |
/// | dir{C1} | 79 | I4 | C4^5 | 16 | 4 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 16 | 4 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 16 | 4 |
/// | dir{P3} | 9 | Cc | Cs^4 | 32 | 4 |
/// | dir{C1} | 8 | Cm | Cs^3 | 32 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P2 | $P_{1}P2$ | $k_{12}t1t2$ | 4 | D16c | yes |
/// | P3P3 | $P_{3}P3$ | $k_{12}t4t4$ | 4 | D8a | no |
/// | P4P4 | $P_{4}P4$ | $k_{12}t3t3$ | 4 | D8a | no |
///
/// ### Isotropy subgroups for P1P2 ($P_{1}P2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 8 | 2 |
/// | dir{P1} | 79 | I4 | C4^5 | 8 | 2 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
///
/// ### Isotropy subgroups for P3P3 ($P_{3}P3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 45 | Iba2 | C2v^21 | 8 | 2 |
///
/// ### Isotropy subgroups for P4P4 ($P_{4}P4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 44 | Imm2 | C2v^20 | 8 | 2 |
///
pub struct Sg108;

/// # 109 I4₁md (C₄ᵥ¹¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{14}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 109 | I4_1md | C4v^11 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 44 | Imm2 | C2v^20 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 76 | P4_1 | C4^2 | 8 | 2 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 8 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M2 | $Z_{1}Z2$ | $k_{15}t1t3$ | 2 | B4a | yes |
/// | M3M4 | $Z_{3}Z4$ | $k_{15}t2t4$ | 2 | B4a | yes |
/// | M5 | $Z_{5}$ | $k_{15}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1M2 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 4 | 1 |
///
/// ### Isotropy subgroups for M3M4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 4 | 1 |
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{11}t1$ | 4 | D64a | no |
/// | N2 | $N_{2}$ | $k_{11}t2$ | 4 | D64a | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 8 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 8 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 16 | 4 |
/// | dir{C1} | 8 | Cm | Cs^3 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ### Isotropy subgroups for N2 ($N_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 8 | 1 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 8 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 16 | 4 |
/// | dir{C1} | 9 | Cc | Cs^4 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1 | $P_{1}$ | $k_{12}t1$ | 4 | D32a | yes |
///
/// ### Isotropy subgroups for P1 ($P_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 8 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 8 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
pub struct Sg109;

/// # 110 I4₁cd (C₄ᵥ¹²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{14}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 110 | I4_1cd | C4v^12 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 45 | Iba2 | C2v^21 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 78 | P4_3 | C4^4 | 8 | 2 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 8 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M2 | $Z_{1}Z2$ | $k_{15}t2t4$ | 2 | B4a | yes |
/// | M3M4 | $Z_{3}Z4$ | $k_{15}t1t3$ | 2 | B4a | yes |
/// | M5 | $Z_{5}$ | $k_{15}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1M2 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 4 | 1 |
///
/// ### Isotropy subgroups for M3M4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 32 | Pba2 | C2v^8 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 78 | P4_3 | C4^4 | 4 | 1 |
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1N2 | $N_{1}N2$ | $k_{11}t1t2$ | 8 | F64b | no |
///
/// ### Isotropy subgroups for N1N2 ($N_{1}N2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P1 | $P_{1}P1$ | $k_{12}t1t1$ | 8 | F32a | no |
///
/// ### Isotropy subgroups for P1P1 ($P_{1}P1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
pub struct Sg110;

/// # 111 P-42m (D₂ₔ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{15}t4$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{15}t2$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{15}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 4 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 4 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 4 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 4 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t4$ | 1 | A2a | no |
/// | M3 | $M_{4}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M4 | $M_{2}$ | $k_{18}t2$ | 1 | A2a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 2 | B8a | no |
/// | R2 | $R_{3}$ | $k_{16}t4$ | 2 | B8a | no |
/// | R3 | $R_{4}$ | $k_{16}t2$ | 2 | B8a | no |
/// | R4 | $R_{2}$ | $k_{16}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 4 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 4 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 2 |
///
/// ### Isotropy subgroups for R3 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 4 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 2 |
///
/// ### Isotropy subgroups for R4 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 4 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{20}t1$ | 1 | A2a | no |
/// | A2 | $A_{3}$ | $k_{20}t4$ | 1 | A2a | no |
/// | A3 | $A_{4}$ | $k_{20}t3$ | 1 | A2a | no |
/// | A4 | $A_{2}$ | $k_{20}t2$ | 1 | A2a | no |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 2 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 46 | Ima2 | C2v^22 | 4 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{19}t4$ | 1 | A2a | no |
/// | Z3 | $Z_{4}$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z4 | $Z_{2}$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg111;

/// # 112 P-42c (D₂ₔ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{15}t4$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{15}t2$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{15}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 4 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 4 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 4 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 4 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{3}$ | $k_{18}t4$ | 1 | A2a | no |
/// | M3 | $M_{4}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M4 | $M_{2}$ | $k_{18}t2$ | 1 | A2a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 2 | B8a | no |
/// | R2 | $R_{3}$ | $k_{16}t2$ | 2 | B8a | no |
/// | R3 | $R_{4}$ | $k_{16}t3$ | 2 | B8a | no |
/// | R4 | $R_{2}$ | $k_{16}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 4 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 4 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 2 |
///
/// ### Isotropy subgroups for R3 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 4 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 2 |
///
/// ### Isotropy subgroups for R4 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 4 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{20}t1t4$ | 2 | B4a | no |
/// | A3A4 | $A_{3}A4$ | $k_{20}t2t3$ | 2 | B4a | no |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 22 | F222 | D2^7 | 4 | 1 |
///
/// ### Isotropy subgroups for A3A4 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 22 | F222 | D2^7 | 4 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 4 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z2 | $Z_{1}Z2$ | $k_{19}t1t4$ | 2 | B4a | no |
/// | Z3Z4 | $Z_{3}Z4$ | $k_{19}t2t3$ | 2 | B4a | no |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1Z2 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
pub struct Sg112;

/// # 113 P-42₁m (D₂ₔ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 8 | 2 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 8 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M3 | $M_{3}M4$ | $k_{18}t1t3$ | 2 | B4a | no |
/// | M2M4 | $M_{1}M2$ | $k_{18}t2t4$ | 2 | B4a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1M3 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 4 | 1 |
///
/// ### Isotropy subgroups for M2M4 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 32 | Pba2 | C2v^8 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 1 |
/// | dir{P1} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 8 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A3 | $A_{3}A4$ | $k_{20}t1t3$ | 2 | B4a | no |
/// | A2A4 | $A_{1}A2$ | $k_{20}t2t4$ | 2 | B4a | no |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1A3 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 44 | Imm2 | C2v^20 | 4 | 1 |
///
/// ### Isotropy subgroups for A2A4 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 45 | Iba2 | C2v^21 | 4 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 4 | 1 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{19}t4$ | 1 | A2a | no |
/// | Z3 | $Z_{4}$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z4 | $Z_{2}$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg113;

/// # 114 P-42₁c (D₂ₔ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 8 | 2 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 8 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M3 | $M_{3}M4$ | $k_{18}t1t3$ | 2 | B4a | no |
/// | M2M4 | $M_{1}M2$ | $k_{18}t2t4$ | 2 | B4a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1M3 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 4 | 1 |
///
/// ### Isotropy subgroups for M2M4 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 1 |
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 8 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A4 | $A_{3}A4$ | $k_{20}t1t2$ | 2 | B4a | no |
/// | A2A3 | $A_{1}A2$ | $k_{20}t3t4$ | 2 | B4a | no |
/// | A5A5 | $A_{5}A5$ | $k_{20}t5t5$ | 4 | D8a | yes |
///
/// ### Isotropy subgroups for A1A4 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 4 | 1 |
///
/// ### Isotropy subgroups for A2A3 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 4 | 1 |
///
/// ### Isotropy subgroups for A5A5 ($A_{5}A5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z2 | $Z_{1}Z2$ | $k_{19}t1t4$ | 2 | B4a | no |
/// | Z3Z4 | $Z_{3}Z4$ | $k_{19}t2t3$ | 2 | B4a | no |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1Z2 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 4 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 4 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
pub struct Sg114;

/// # 115 P-4m2 (D₂ₔ⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{15}t2$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{15}t3$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{15}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 4 | 1 |
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 4 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 4 | 2 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 4 | 1 |
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 4 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 4 | 2 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{4}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3 | $M_{2}$ | $k_{18}t2$ | 1 | A2a | no |
/// | M4 | $M_{3}$ | $k_{18}t4$ | 1 | A2a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 2 | B8a | no |
/// | R2 | $R_{3}$ | $k_{16}t2$ | 2 | B8a | no |
/// | R3 | $R_{4}$ | $k_{16}t3$ | 2 | B8a | no |
/// | R4 | $R_{2}$ | $k_{16}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 4 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 4 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 8 | 2 |
///
/// ### Isotropy subgroups for R3 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 4 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 8 | 2 |
///
/// ### Isotropy subgroups for R4 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 4 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 8 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{20}t1$ | 1 | A2a | no |
/// | A2 | $A_{4}$ | $k_{20}t3$ | 1 | A2a | no |
/// | A3 | $A_{2}$ | $k_{20}t2$ | 1 | A2a | no |
/// | A4 | $A_{3}$ | $k_{20}t4$ | 1 | A2a | no |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 2 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 4 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{4}$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z3 | $Z_{2}$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z4 | $Z_{3}$ | $k_{19}t4$ | 1 | A2a | no |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 4 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg115;

/// # 116 P-4c2 (D₂ₔ⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{15}t2$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{15}t3$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{15}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 4 | 1 |
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 4 | 2 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 4 | 2 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 4 | 1 |
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 4 | 2 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 4 | 2 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | M2 | $M_{4}$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3 | $M_{2}$ | $k_{18}t2$ | 1 | A2a | no |
/// | M4 | $M_{3}$ | $k_{18}t4$ | 1 | A2a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{16}t1t2$ | 4 | D16c | yes |
/// | R3R4 | $R_{3}R4$ | $k_{16}t3t4$ | 4 | D16c | yes |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
///
/// ### Isotropy subgroups for R3R4 ($R_{3}R4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A4 | $A_{1}A2$ | $k_{20}t2t3$ | 2 | B4a | no |
/// | A2A3 | $A_{3}A4$ | $k_{20}t1t4$ | 2 | B4a | no |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for A1A4 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 4 | 1 |
///
/// ### Isotropy subgroups for A2A3 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 4 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 4 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z4 | $Z_{1}Z2$ | $k_{19}t2t3$ | 2 | B4a | no |
/// | Z2Z3 | $Z_{3}Z4$ | $k_{19}t1t4$ | 2 | B4a | no |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1Z4 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
///
/// ### Isotropy subgroups for Z2Z3 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
pub struct Sg116;

/// # 117 P-4b2 (D₂ₔ⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 32 | Pba2 | C2v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 8 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M4 | $M_{1}M2$ | $k_{18}t2t3$ | 2 | B4a | no |
/// | M2M3 | $M_{3}M4$ | $k_{18}t1t4$ | 2 | B4a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for M1M4 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
///
/// ### Isotropy subgroups for M2M3 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A4 | $A_{1}A2$ | $k_{20}t2t3$ | 2 | B4a | no |
/// | A2A3 | $A_{3}A4$ | $k_{20}t1t4$ | 2 | B4a | no |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for A1A4 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 4 | 1 |
///
/// ### Isotropy subgroups for A2A3 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 4 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 4 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{4}$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z3 | $Z_{2}$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z4 | $Z_{3}$ | $k_{19}t4$ | 1 | A2a | no |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 4 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg117;

/// # 118 P-4n2 (D₂ₔ⁸)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{17}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 8 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M4 | $M_{1}M2$ | $k_{18}t2t3$ | 2 | B4a | no |
/// | M2M3 | $M_{3}M4$ | $k_{18}t1t4$ | 2 | B4a | no |
/// | M5 | $M_{5}$ | $k_{18}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for M1M4 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
///
/// ### Isotropy subgroups for M2M3 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{20}t2$ | 1 | A2a | no |
/// | A2 | $A_{4}$ | $k_{20}t4$ | 1 | A2a | no |
/// | A3 | $A_{2}$ | $k_{20}t1$ | 1 | A2a | no |
/// | A4 | $A_{3}$ | $k_{20}t3$ | 1 | A2a | no |
/// | A5 | $A_{5}$ | $k_{20}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 2 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 4 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z4 | $Z_{1}Z2$ | $k_{19}t2t3$ | 2 | B4a | no |
/// | Z2Z3 | $Z_{3}Z4$ | $k_{19}t1t4$ | 2 | B4a | no |
/// | Z5 | $Z_{5}$ | $k_{19}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1Z4 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
///
/// ### Isotropy subgroups for Z2Z3 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
///
/// ### Isotropy subgroups for Z5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
pub struct Sg118;

/// # 119 I-4m2 (D₂ₔ⁹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{4}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{14}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 44 | Imm2 | C2v^20 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 22 | F222 | D2^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{13}t4$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{13}t3$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{13}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 4 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 4 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 4 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 4 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $Z_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | M2 | $Z_{4}$ | $k_{15}t3$ | 1 | A2a | no |
/// | M3 | $Z_{2}$ | $k_{15}t2$ | 1 | A2a | no |
/// | M4 | $Z_{3}$ | $k_{15}t4$ | 1 | A2a | no |
/// | M5 | $Z_{5}$ | $k_{15}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 4 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{11}t1$ | 4 | D64c | no |
/// | N2 | $N_{2}$ | $k_{11}t2$ | 4 | D64c | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 8 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 8 | 4 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 8 | 4 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 16 | 4 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 16 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ### Isotropy subgroups for N2 ($N_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 8 | 1 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 8 | 4 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 8 | 4 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 16 | 4 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 16 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1 | $P_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | P2 | $P_{3}$ | $k_{12}t3$ | 2 | B8a | no |
/// | P3 | $P_{2}$ | $k_{12}t2$ | 2 | B8a | no |
/// | P4 | $P_{4}$ | $k_{12}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for P1 ($P_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 4 | 2 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 4 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
///
/// ### Isotropy subgroups for P2 ($P_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 4 | 2 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 4 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
///
/// ### Isotropy subgroups for P3 ($P_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 4 | 2 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 4 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
///
/// ### Isotropy subgroups for P4 ($P_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 4 | 2 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 4 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
///
pub struct Sg119;

/// # 120 I-4c2 (D₂ₔ¹⁰)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{4}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{14}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 45 | Iba2 | C2v^21 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 22 | F222 | D2^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{13}t4$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{13}t3$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{13}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 4 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 4 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 4 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 4 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $Z_{1}$ | $k_{15}t2$ | 1 | A2a | no |
/// | M2 | $Z_{4}$ | $k_{15}t4$ | 1 | A2a | no |
/// | M3 | $Z_{2}$ | $k_{15}t1$ | 1 | A2a | no |
/// | M4 | $Z_{3}$ | $k_{15}t3$ | 1 | A2a | no |
/// | M5 | $Z_{5}$ | $k_{15}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 4 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1N2 | $N_{1}N2$ | $k_{11}t1t2$ | 8 | F64a | no |
///
/// ### Isotropy subgroups for N1N2 ($N_{1}N2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{C1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P3} | 82 | I-4 | S4^2 | 16 | 4 |
/// | dir{P1} | 82 | I-4 | S4^2 | 16 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 16 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 16 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P2 | $P_{1}P3$ | $k_{12}t1t3$ | 4 | D16c | no |
/// | P3P4 | $P_{2}P4$ | $k_{12}t2t4$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for P1P2 ($P_{1}P3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
///
/// ### Isotropy subgroups for P3P4 ($P_{2}P4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
///
pub struct Sg120;

/// # 121 I-42m (D₂ₔ¹¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{14}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | X2 | $X_{3}$ | $k_{13}t2$ | 2 | B8a | no |
/// | X3 | $X_{4}$ | $k_{13}t4$ | 2 | B8a | no |
/// | X4 | $X_{2}$ | $k_{13}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 4 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 41 | Aba2 | C2v^17 | 4 | 1 |
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 4 | 2 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 2 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 40 | Ama2 | C2v^16 | 4 | 1 |
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 4 | 2 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 2 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 4 | 2 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $Z_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | M2 | $Z_{3}$ | $k_{15}t4$ | 1 | A2a | no |
/// | M3 | $Z_{4}$ | $k_{15}t3$ | 1 | A2a | no |
/// | M4 | $Z_{2}$ | $k_{15}t2$ | 1 | A2a | no |
/// | M5 | $Z_{5}$ | $k_{15}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{11}t1$ | 4 | D64c | no |
/// | N2 | $N_{2}$ | $k_{11}t2$ | 4 | D64c | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 8 | 4 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 8 | 4 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 4 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 4 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 4 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ### Isotropy subgroups for N2 ($N_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 8 | 4 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 8 | 4 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 4 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 4 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 4 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P1 | $P_{1}P1$ | $k_{12}t1t1*$ | 2 | B4a | no |
/// | P2P2 | $P_{3}P3$ | $k_{12}t4t4*$ | 2 | B4a | no |
/// | P3P3 | $P_{4}P4$ | $k_{12}t3t3*$ | 2 | B4a | no |
/// | P4P4 | $P_{2}P2$ | $k_{12}t2t2*$ | 2 | B4a | no |
/// | P5P5 | $P_{5}P5$ | $k_{12}t5t5*$ | 4 | D16c | yes |
///
/// ### Isotropy subgroups for P1P1 ($P_{1}P1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 4 | 1 |
///
/// ### Isotropy subgroups for P2P2 ($P_{3}P3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 4 | 1 |
///
/// ### Isotropy subgroups for P3P3 ($P_{4}P4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 4 | 1 |
///
/// ### Isotropy subgroups for P4P4 ($P_{2}P2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 4 | 1 |
///
/// ### Isotropy subgroups for P5P5 ($P_{5}P5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 1 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 8 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
///
pub struct Sg121;

/// # 122 I-42d (D₂ₔ¹²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{14}t5$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 8 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M2 | $Z_{3}Z4$ | $k_{15}t2t3$ | 2 | B4a | no |
/// | M3M4 | $Z_{1}Z2$ | $k_{15}t1t4$ | 2 | B4a | no |
/// | M5 | $Z_{5}$ | $k_{15}t5$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for M1M2 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
///
/// ### Isotropy subgroups for M3M4 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 4 | 1 |
///
/// ### Isotropy subgroups for M5 ($Z_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{2}$ | $k_{11}t2$ | 4 | D64b | no |
/// | N2 | $N_{1}$ | $k_{11}t1$ | 4 | D64b | no |
///
/// ### Isotropy subgroups for N1 ($N_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 16 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ### Isotropy subgroups for N2 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 16 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P1 | $P_{1}P1$ | $k_{12}t1t1*$ | 4 | D32e | yes |
/// | P2P2 | $P_{2}P2$ | $k_{12}t2t2*$ | 4 | D32e | yes |
///
/// ### Isotropy subgroups for P1P1 ($P_{1}P1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 1 |
///
/// ### Isotropy subgroups for P2P2 ($P_{2}P2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 1 |
///
pub struct Sg122;

/// # 123 P4/mmm (D₄ₕ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 47 | Pmmm | D2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 83 | P4/m | C4h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 99 | P4mm | C4v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 4 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{15}t1$ | 2 | B8a | no |
/// | X1- | $X_{1}^-$ | $k_{15}t2$ | 2 | B8a | no |
/// | X2+ | $X_{3}^+$ | $k_{15}t7$ | 2 | B8a | no |
/// | X2- | $X_{3}^-$ | $k_{15}t8$ | 2 | B8a | no |
/// | X3+ | $X_{4}^+$ | $k_{15}t3$ | 2 | B8a | no |
/// | X3- | $X_{4}^-$ | $k_{15}t4$ | 2 | B8a | no |
/// | X4+ | $X_{2}^+$ | $k_{15}t5$ | 2 | B8a | no |
/// | X4- | $X_{2}^-$ | $k_{15}t6$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 47 | Pmmm | D2h^1 | 4 | 1 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 4 | 2 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 4 | 2 |
/// | dir{P1} | 50 | Pban | D2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 127 | P4/mbm | D4h^5 | 4 | 2 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 8 | 2 |
///
/// ### Isotropy subgroups for X2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 4 | 2 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 8 | 2 |
///
/// ### Isotropy subgroups for X3+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 4 | 2 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 8 | 2 |
///
/// ### Isotropy subgroups for X3- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 127 | P4/mbm | D4h^5 | 4 | 2 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 8 | 2 |
///
/// ### Isotropy subgroups for X4+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 4 | 2 |
/// | dir{P1} | 50 | Pban | D2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X4- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 47 | Pmmm | D2h^1 | 4 | 1 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 4 | 2 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{18}t1$ | 1 | A2a | no |
/// | M1- | $M_{1}^-$ | $k_{18}t2$ | 1 | A2a | no |
/// | M2+ | $M_{3}^+$ | $k_{18}t5$ | 1 | A2a | no |
/// | M2- | $M_{3}^-$ | $k_{18}t6$ | 1 | A2a | no |
/// | M3+ | $M_{2}^+$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3- | $M_{2}^-$ | $k_{18}t4$ | 1 | A2a | no |
/// | M4+ | $M_{4}^+$ | $k_{18}t7$ | 1 | A2a | no |
/// | M4- | $M_{4}^-$ | $k_{18}t8$ | 1 | A2a | no |
/// | M5+ | $M_{5}^+$ | $k_{18}t9$ | 2 | B8a | no |
/// | M5- | $M_{5}^-$ | $k_{18}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M2+ ($M_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 127 | P4/mbm | D4h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for M2- ($M_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for M3+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 127 | P4/mbm | D4h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for M3- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for M4+ ($M_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M4- ($M_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M5+ ($M_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 4 | 1 |
/// | dir{P1} | 53 | Pmna | D2h^7 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
///
/// ### Isotropy subgroups for M5- ($M_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 4 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{16}t1$ | 2 | B8a | no |
/// | R1- | $R_{1}^-$ | $k_{16}t2$ | 2 | B8a | no |
/// | R2+ | $R_{3}^+$ | $k_{16}t7$ | 2 | B8a | no |
/// | R2- | $R_{3}^-$ | $k_{16}t8$ | 2 | B8a | no |
/// | R3+ | $R_{4}^+$ | $k_{16}t3$ | 2 | B8a | no |
/// | R3- | $R_{4}^-$ | $k_{16}t4$ | 2 | B8a | no |
/// | R4+ | $R_{2}^+$ | $k_{16}t5$ | 2 | B8a | no |
/// | R4- | $R_{2}^-$ | $k_{16}t6$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 4 | 1 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 4 | 2 |
/// | dir{P1} | 71 | Immm | D2h^25 | 8 | 2 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 4 | 1 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 4 | 2 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 8 | 2 |
///
/// ### Isotropy subgroups for R2+ ($R_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 4 | 1 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 4 | 2 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 8 | 2 |
///
/// ### Isotropy subgroups for R2- ($R_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 4 | 1 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 4 | 2 |
/// | dir{P1} | 71 | Immm | D2h^25 | 8 | 2 |
///
/// ### Isotropy subgroups for R3+ ($R_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 4 | 1 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 4 | 2 |
/// | dir{P1} | 71 | Immm | D2h^25 | 8 | 2 |
///
/// ### Isotropy subgroups for R3- ($R_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 4 | 1 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 4 | 2 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 8 | 2 |
///
/// ### Isotropy subgroups for R4+ ($R_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 4 | 1 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 4 | 2 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 8 | 2 |
///
/// ### Isotropy subgroups for R4- ($R_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 4 | 1 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 4 | 2 |
/// | dir{P1} | 71 | Immm | D2h^25 | 8 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1+ | $A_{1}^+$ | $k_{20}t1$ | 1 | A2a | no |
/// | A1- | $A_{1}^-$ | $k_{20}t2$ | 1 | A2a | no |
/// | A2+ | $A_{3}^+$ | $k_{20}t5$ | 1 | A2a | no |
/// | A2- | $A_{3}^-$ | $k_{20}t6$ | 1 | A2a | no |
/// | A3+ | $A_{2}^+$ | $k_{20}t3$ | 1 | A2a | no |
/// | A3- | $A_{2}^-$ | $k_{20}t4$ | 1 | A2a | no |
/// | A4+ | $A_{4}^+$ | $k_{20}t7$ | 1 | A2a | no |
/// | A4- | $A_{4}^-$ | $k_{20}t8$ | 1 | A2a | no |
/// | A5+ | $A_{5}^+$ | $k_{20}t9$ | 2 | B8a | no |
/// | A5- | $A_{5}^-$ | $k_{20}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1+ ($A_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 2 | 1 |
///
/// ### Isotropy subgroups for A1- ($A_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for A2+ ($A_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for A2- ($A_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 2 | 1 |
///
/// ### Isotropy subgroups for A3+ ($A_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for A3- ($A_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 2 | 1 |
///
/// ### Isotropy subgroups for A4+ ($A_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 2 | 1 |
///
/// ### Isotropy subgroups for A4- ($A_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for A5+ ($A_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 74 | Imma | D2h^28 | 4 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 4 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
///
/// ### Isotropy subgroups for A5- ($A_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 74 | Imma | D2h^28 | 4 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 4 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{3}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | Z2- | $Z_{3}^-$ | $k_{19}t6$ | 1 | A2a | no |
/// | Z3+ | $Z_{2}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z3- | $Z_{2}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | Z4+ | $Z_{4}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | Z4- | $Z_{4}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | Z5+ | $Z_{5}^+$ | $k_{19}t9$ | 2 | B8a | no |
/// | Z5- | $Z_{5}^-$ | $k_{19}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 124 | P4/mcc | D4h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 124 | P4/mcc | D4h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4+ ($Z_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4- ($Z_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 2 | 1 |
///
/// ### Isotropy subgroups for Z5+ ($Z_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 4 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 8 | 1 |
///
/// ### Isotropy subgroups for Z5- ($Z_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 4 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 8 | 1 |
///
pub struct Sg123;

/// # 124 P4/mcc (D₄ₕ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 124 | P4/mcc | D4h^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 83 | P4/m | C4h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 103 | P4cc | C4v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 66 | Cccm | D2h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 4 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{15}t1$ | 2 | B8a | no |
/// | X1- | $X_{1}^-$ | $k_{15}t2$ | 2 | B8a | no |
/// | X2+ | $X_{3}^+$ | $k_{15}t7$ | 2 | B8a | no |
/// | X2- | $X_{3}^-$ | $k_{15}t8$ | 2 | B8a | no |
/// | X3+ | $X_{4}^+$ | $k_{15}t3$ | 2 | B8a | no |
/// | X3- | $X_{4}^-$ | $k_{15}t4$ | 2 | B8a | no |
/// | X4+ | $X_{2}^+$ | $k_{15}t5$ | 2 | B8a | no |
/// | X4- | $X_{2}^-$ | $k_{15}t6$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
/// | dir{P1} | 124 | P4/mcc | D4h^2 | 4 | 2 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 50 | Pban | D2h^4 | 4 | 1 |
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 4 | 2 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 8 | 2 |
///
/// ### Isotropy subgroups for X2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 53 | Pmna | D2h^7 | 4 | 1 |
/// | dir{P1} | 128 | P4/mnc | D4h^6 | 4 | 2 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 8 | 2 |
///
/// ### Isotropy subgroups for X2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 130 | P4/ncc | D4h^8 | 4 | 2 |
/// | dir{P1} | 56 | Pccn | D2h^10 | 8 | 2 |
///
/// ### Isotropy subgroups for X3+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 130 | P4/ncc | D4h^8 | 4 | 2 |
/// | dir{P1} | 56 | Pccn | D2h^10 | 8 | 2 |
///
/// ### Isotropy subgroups for X3- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 53 | Pmna | D2h^7 | 4 | 1 |
/// | dir{P1} | 128 | P4/mnc | D4h^6 | 4 | 2 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 8 | 2 |
///
/// ### Isotropy subgroups for X4+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 50 | Pban | D2h^4 | 4 | 1 |
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 4 | 2 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 8 | 2 |
///
/// ### Isotropy subgroups for X4- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
/// | dir{P1} | 124 | P4/mcc | D4h^2 | 4 | 2 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{18}t1$ | 1 | A2a | no |
/// | M1- | $M_{1}^-$ | $k_{18}t2$ | 1 | A2a | no |
/// | M2+ | $M_{3}^+$ | $k_{18}t5$ | 1 | A2a | no |
/// | M2- | $M_{3}^-$ | $k_{18}t6$ | 1 | A2a | no |
/// | M3+ | $M_{2}^+$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3- | $M_{2}^-$ | $k_{18}t4$ | 1 | A2a | no |
/// | M4+ | $M_{4}^+$ | $k_{18}t7$ | 1 | A2a | no |
/// | M4- | $M_{4}^-$ | $k_{18}t8$ | 1 | A2a | no |
/// | M5+ | $M_{5}^+$ | $k_{18}t9$ | 2 | B8a | no |
/// | M5- | $M_{5}^-$ | $k_{18}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 124 | P4/mcc | D4h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M2+ ($M_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 128 | P4/mnc | D4h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M2- ($M_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 130 | P4/ncc | D4h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for M3+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 128 | P4/mnc | D4h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M3- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 130 | P4/ncc | D4h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for M4+ ($M_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 124 | P4/mcc | D4h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M4- ($M_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M5+ ($M_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 68 | Ccca | D2h^22 | 4 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
///
/// ### Isotropy subgroups for M5- ($M_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 66 | Cccm | D2h^20 | 4 | 1 |
/// | dir{P1} | 53 | Pmna | D2h^7 | 4 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32c | yes |
/// | R2 | $R_{2}$ | $k_{16}t2$ | 4 | D32c | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 8 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 8 | 2 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 8 | 2 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 8 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 16 | 2 |
/// | dir{P1} | 79 | I4 | C4^5 | 16 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 8 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 8 | 2 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 8 | 2 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 8 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 16 | 2 |
/// | dir{P1} | 79 | I4 | C4^5 | 16 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{2}$ | $k_{20}t2$ | 2 | B8a | yes |
/// | A2 | $A_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | A3A4 | $A_{3}A4$ | $k_{20}t3t4$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for A1 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 4 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 4 | 1 |
/// | dir{P1} | 79 | I4 | C4^5 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 4 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 4 | 1 |
/// | dir{P1} | 79 | I4 | C4^5 | 8 | 1 |
///
/// ### Isotropy subgroups for A3A4 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{2}$ | $k_{19}t2$ | 2 | B8a | yes |
/// | Z2 | $Z_{1}$ | $k_{19}t1$ | 2 | B8a | yes |
/// | Z3Z4 | $Z_{3}Z4$ | $k_{19}t3t4$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 4 | 1 |
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 4 | 1 |
/// | dir{P1} | 77 | P4_2 | C4^3 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 4 | 1 |
/// | dir{P1} | 83 | P4/m | C4h^1 | 4 | 1 |
/// | dir{P1} | 75 | P4 | C4^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 1 |
///
pub struct Sg124;

/// # 125 P4/nbm (D₄ₕ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 50 | Pban | D2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 85 | P4/n | C4h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 100 | P4bm | C4v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D64c | yes |
/// | X2 | $X_{2}$ | $k_{15}t2$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 1 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 8 | 2 |
/// | dir{P1} | 89 | P422 | D4^1 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 16 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 1 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 8 | 2 |
/// | dir{P1} | 90 | P42_12 | D4^2 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 16 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{3}$ | $k_{18}t1$ | 2 | B8a | no |
/// | M2 | $M_{1}$ | $k_{18}t2$ | 2 | B8a | no |
/// | M3 | $M_{4}$ | $k_{18}t4$ | 2 | B8a | yes |
/// | M4 | $M_{2}$ | $k_{18}t3$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for M1 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 4 | 1 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 4 | 1 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 4 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 4 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D64c | yes |
/// | R2 | $R_{2}$ | $k_{16}t2$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 8 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 8 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{3}$ | $k_{20}t1$ | 2 | B8a | no |
/// | A2 | $A_{1}$ | $k_{20}t2$ | 2 | B8a | no |
/// | A3 | $A_{4}$ | $k_{20}t4$ | 2 | B8a | yes |
/// | A4 | $A_{2}$ | $k_{20}t3$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 4 | 1 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 4 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 4 | 1 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 4 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 4 | 1 |
/// | dir{P1} | 74 | Imma | D2h^28 | 4 | 1 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 8 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 4 | 1 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 4 | 1 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{3}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | Z2- | $Z_{3}^-$ | $k_{19}t6$ | 1 | A2a | no |
/// | Z3+ | $Z_{2}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z3- | $Z_{2}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | Z4+ | $Z_{4}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | Z4- | $Z_{4}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | Z5+ | $Z_{5}^+$ | $k_{19}t9$ | 2 | B8a | no |
/// | Z5- | $Z_{5}^-$ | $k_{19}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 133 | P4_2/nbc | D4h^11 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4+ ($Z_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4- ($Z_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 133 | P4_2/nbc | D4h^11 | 2 | 1 |
///
/// ### Isotropy subgroups for Z5+ ($Z_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 4 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 4 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
///
/// ### Isotropy subgroups for Z5- ($Z_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 4 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 4 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
///
pub struct Sg125;

/// # 126 P4/nnc (D₄ₕ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 48 | Pnnn | D2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 85 | P4/n | C4h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 104 | P4nc | C4v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 68 | Ccca | D2h^22 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 4 | 1 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 4 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D64c | yes |
/// | X2 | $X_{2}$ | $k_{15}t2$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 1 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 8 | 2 |
/// | dir{P1} | 89 | P422 | D4^1 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 16 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 1 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 8 | 2 |
/// | dir{P1} | 90 | P42_12 | D4^2 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 16 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{3}$ | $k_{18}t1$ | 2 | B8a | no |
/// | M2 | $M_{1}$ | $k_{18}t2$ | 2 | B8a | no |
/// | M3 | $M_{4}$ | $k_{18}t4$ | 2 | B8a | yes |
/// | M4 | $M_{2}$ | $k_{18}t3$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for M1 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 4 | 1 |
/// | dir{P1} | 50 | Pban | D2h^4 | 4 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 4 | 1 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 4 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 4 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 4 | 1 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D64c | yes |
/// | R2 | $R_{2}$ | $k_{16}t2$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 1 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 8 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 1 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 8 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{3}$ | $k_{20}t1$ | 2 | B8a | no |
/// | A2 | $A_{1}$ | $k_{20}t2$ | 2 | B8a | no |
/// | A3 | $A_{2}$ | $k_{20}t4$ | 2 | B8a | yes |
/// | A4 | $A_{4}$ | $k_{20}t3$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 4 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 4 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 4 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 4 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 4 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 4 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 8 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 4 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 4 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{2}$ | $k_{19}t2$ | 2 | B8a | yes |
/// | Z2 | $Z_{1}$ | $k_{19}t1$ | 2 | B8a | yes |
/// | Z3Z4 | $Z_{3}Z4$ | $k_{19}t3t4$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 4 | 1 |
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 4 | 1 |
/// | dir{P1} | 77 | P4_2 | C4^3 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 4 | 1 |
/// | dir{P1} | 85 | P4/n | C4h^3 | 4 | 1 |
/// | dir{P1} | 75 | P4 | C4^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 1 |
///
pub struct Sg126;

/// # 127 P4/mbm (D₄ₕ⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 127 | P4/mbm | D4h^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 55 | Pbam | D2h^9 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 83 | P4/m | C4h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 100 | P4bm | C4v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 4 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{2}$ | $k_{15}t2$ | 4 | D32b | yes |
/// | X2 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 8 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 1 |
/// | dir{P1} | 85 | P4/n | C4h^3 | 8 | 2 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 8 | 2 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 16 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 16 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 16 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 32 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 8 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 16 | 1 |
/// | dir{P1} | 83 | P4/m | C4h^1 | 8 | 2 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 8 | 2 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 16 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 16 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 16 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+M4+ | $M_{1}^-M2-$ | $k_{18}t1t7$ | 2 | B4a | no |
/// | M1-M4- | $M_{3}^+M4+$ | $k_{18}t2t8$ | 2 | B4a | no |
/// | M2+M3+ | $M_{3}^-M4-$ | $k_{18}t3t5$ | 2 | B4a | no |
/// | M2-M3- | $M_{1}^+M2+$ | $k_{18}t4t6$ | 2 | B4a | no |
/// | M5+ | $M_{5}^-$ | $k_{18}t9$ | 2 | B8a | no |
/// | M5- | $M_{5}^+$ | $k_{18}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1+M4+ ($M_{1}^-M2-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 59 | Pmmn | D2h^13 | 4 | 1 |
///
/// ### Isotropy subgroups for M1-M4- ($M_{3}^+M4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 55 | Pbam | D2h^9 | 4 | 1 |
///
/// ### Isotropy subgroups for M2+M3+ ($M_{3}^-M4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 50 | Pban | D2h^4 | 4 | 1 |
///
/// ### Isotropy subgroups for M2-M3- ($M_{1}^+M2+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 47 | Pmmm | D2h^1 | 4 | 1 |
///
/// ### Isotropy subgroups for M5+ ($M_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 83 | P4/m | C4h^1 | 4 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M5- ($M_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 85 | P4/n | C4h^3 | 4 | 1 |
/// | dir{P1} | 53 | Pmna | D2h^7 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{2}$ | $k_{16}t2$ | 4 | D32b | yes |
/// | R2 | $R_{1}$ | $k_{16}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for R1 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 8 | 2 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 8 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 8 | 2 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 8 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1+A4+ | $A_{1}^-A2-$ | $k_{20}t1t7$ | 2 | B4a | no |
/// | A1-A4- | $A_{3}^+A4+$ | $k_{20}t2t8$ | 2 | B4a | no |
/// | A2+A3+ | $A_{3}^-A4-$ | $k_{20}t3t5$ | 2 | B4a | no |
/// | A2-A3- | $A_{1}^+A2+$ | $k_{20}t4t6$ | 2 | B4a | no |
/// | A5+ | $A_{5}^-$ | $k_{20}t9$ | 2 | B8a | no |
/// | A5- | $A_{5}^+$ | $k_{20}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1+A4+ ($A_{1}^-A2-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 71 | Immm | D2h^25 | 4 | 1 |
///
/// ### Isotropy subgroups for A1-A4- ($A_{3}^+A4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 72 | Ibam | D2h^26 | 4 | 1 |
///
/// ### Isotropy subgroups for A2+A3+ ($A_{3}^-A4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 72 | Ibam | D2h^26 | 4 | 1 |
///
/// ### Isotropy subgroups for A2-A3- ($A_{1}^+A2+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 71 | Immm | D2h^25 | 4 | 1 |
///
/// ### Isotropy subgroups for A5+ ($A_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 87 | I4/m | C4h^5 | 4 | 1 |
/// | dir{P1} | 74 | Imma | D2h^28 | 4 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
///
/// ### Isotropy subgroups for A5- ($A_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 87 | I4/m | C4h^5 | 4 | 1 |
/// | dir{P1} | 74 | Imma | D2h^28 | 4 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{3}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | Z2- | $Z_{3}^-$ | $k_{19}t6$ | 1 | A2a | no |
/// | Z3+ | $Z_{2}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z3- | $Z_{2}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | Z4+ | $Z_{4}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | Z4- | $Z_{4}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | Z5+ | $Z_{5}^+$ | $k_{19}t9$ | 2 | B8a | no |
/// | Z5- | $Z_{5}^-$ | $k_{19}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 127 | P4/mbm | D4h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 128 | P4/mnc | D4h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 135 | P4_2/mbc | D4h^13 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 128 | P4/mnc | D4h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 127 | P4/mbm | D4h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4+ ($Z_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4- ($Z_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 135 | P4_2/mbc | D4h^13 | 2 | 1 |
///
/// ### Isotropy subgroups for Z5+ ($Z_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 4 | 1 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 4 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 8 | 1 |
///
/// ### Isotropy subgroups for Z5- ($Z_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 4 | 1 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 4 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 8 | 1 |
///
pub struct Sg127;

/// # 128 P4/mnc (D₄ₕ⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 128 | P4/mnc | D4h^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 58 | Pnnm | D2h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 83 | P4/m | C4h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 104 | P4nc | C4v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 66 | Cccm | D2h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 4 | 1 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 4 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{2}$ | $k_{15}t2$ | 4 | D32b | yes |
/// | X2 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 8 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 1 |
/// | dir{P1} | 85 | P4/n | C4h^3 | 8 | 2 |
/// | dir{P1} | 68 | Ccca | D2h^22 | 8 | 2 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 16 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 16 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 16 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 32 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 8 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 16 | 1 |
/// | dir{P1} | 83 | P4/m | C4h^1 | 8 | 2 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 8 | 2 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 16 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 16 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 16 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+M4+ | $M_{1}^-M2-$ | $k_{18}t1t7$ | 2 | B4a | no |
/// | M1-M4- | $M_{3}^+M4+$ | $k_{18}t2t8$ | 2 | B4a | no |
/// | M2+M3+ | $M_{3}^-M4-$ | $k_{18}t3t5$ | 2 | B4a | no |
/// | M2-M3- | $M_{1}^+M2+$ | $k_{18}t4t6$ | 2 | B4a | no |
/// | M5+ | $M_{5}^-$ | $k_{18}t9$ | 2 | B8a | no |
/// | M5- | $M_{5}^+$ | $k_{18}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1+M4+ ($M_{1}^-M2-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 56 | Pccn | D2h^10 | 4 | 1 |
///
/// ### Isotropy subgroups for M1-M4- ($M_{3}^+M4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 58 | Pnnm | D2h^12 | 4 | 1 |
///
/// ### Isotropy subgroups for M2+M3+ ($M_{3}^-M4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 48 | Pnnn | D2h^2 | 4 | 1 |
///
/// ### Isotropy subgroups for M2-M3- ($M_{1}^+M2+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
///
/// ### Isotropy subgroups for M5+ ($M_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 83 | P4/m | C4h^1 | 4 | 1 |
/// | dir{P1} | 53 | Pmna | D2h^7 | 4 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M5- ($M_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 85 | P4/n | C4h^3 | 4 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{16}t1$ | 4 | D32b | no |
/// | R1- | $R_{1}^-$ | $k_{16}t2$ | 4 | D32b | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 8 | 2 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 8 | 2 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{2}$ | $k_{20}t3$ | 2 | B8a | no |
/// | A2 | $A_{1}$ | $k_{20}t4$ | 2 | B8a | no |
/// | A3A4 | $A_{3}A4$ | $k_{20}t1t2$ | 4 | D16c | yes |
///
/// ### Isotropy subgroups for A1 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 4 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 4 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 4 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 4 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 1 |
///
/// ### Isotropy subgroups for A3A4 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 8 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{2}$ | $k_{19}t2$ | 2 | B8a | yes |
/// | Z2 | $Z_{1}$ | $k_{19}t1$ | 2 | B8a | yes |
/// | Z3Z4 | $Z_{3}Z4$ | $k_{19}t3t4$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 4 | 1 |
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 4 | 1 |
/// | dir{P1} | 77 | P4_2 | C4^3 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 4 | 1 |
/// | dir{P1} | 83 | P4/m | C4h^1 | 4 | 1 |
/// | dir{P1} | 75 | P4 | C4^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 8 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 1 |
///
pub struct Sg128;

/// # 129 P4/nmm (D₄ₕ⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 59 | Pmmn | D2h^13 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 85 | P4/n | C4h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 99 | P4mm | C4v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 4 | 1 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D64c | yes |
/// | X2 | $X_{2}$ | $k_{15}t2$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 8 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 16 | 1 |
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 8 | 2 |
/// | dir{P1} | 99 | P4mm | C4v^1 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 16 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 1 |
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 8 | 2 |
/// | dir{P1} | 100 | P4bm | C4v^2 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 16 | 2 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{4}$ | $k_{18}t2$ | 2 | B8a | no |
/// | M2 | $M_{2}$ | $k_{18}t1$ | 2 | B8a | no |
/// | M3 | $M_{1}$ | $k_{18}t4$ | 2 | B8a | no |
/// | M4 | $M_{3}$ | $k_{18}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 4 | 1 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 4 | 1 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 99 | P4mm | C4v^1 | 4 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 100 | P4bm | C4v^2 | 4 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D64c | yes |
/// | R2 | $R_{2}$ | $k_{16}t2$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 8 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 8 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 16 | 1 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 8 | 2 |
/// | dir{P1} | 107 | I4mm | C4v^9 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 16 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 39 | Abm2 | C2v^15 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 1 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 8 | 2 |
/// | dir{P1} | 108 | I4cm | C4v^10 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 16 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{4}$ | $k_{20}t2$ | 2 | B8a | no |
/// | A2 | $A_{2}$ | $k_{20}t1$ | 2 | B8a | no |
/// | A3 | $A_{1}$ | $k_{20}t4$ | 2 | B8a | no |
/// | A4 | $A_{3}$ | $k_{20}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 4 | 1 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 4 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 4 | 1 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 4 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 107 | I4mm | C4v^9 | 4 | 1 |
/// | dir{P1} | 74 | Imma | D2h^28 | 4 | 1 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 8 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 108 | I4cm | C4v^10 | 4 | 1 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 4 | 1 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{19}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{3}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | Z2- | $Z_{3}^-$ | $k_{19}t6$ | 1 | A2a | no |
/// | Z3+ | $Z_{2}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | Z3- | $Z_{2}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | Z4+ | $Z_{4}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | Z4- | $Z_{4}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | Z5+ | $Z_{5}^+$ | $k_{19}t9$ | 2 | B8a | no |
/// | Z5- | $Z_{5}^-$ | $k_{19}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 130 | P4/ncc | D4h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 137 | P4_2/nmc | D4h^15 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 138 | P4_2/ncm | D4h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 130 | P4/ncc | D4h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4+ ($Z_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 138 | P4_2/ncm | D4h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4- ($Z_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 137 | P4_2/nmc | D4h^15 | 2 | 1 |
///
/// ### Isotropy subgroups for Z5+ ($Z_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 4 | 1 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 4 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
///
/// ### Isotropy subgroups for Z5- ($Z_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 4 | 1 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 4 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
///
pub struct Sg129;

/// # 130 P4/ncc (D₄ₕ⁸)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 130 | P4/ncc | D4h^8 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 56 | Pccn | D2h^10 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 85 | P4/n | C4h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 103 | P4cc | C4v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 68 | Ccca | D2h^22 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 4 | 1 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 4 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D64c | yes |
/// | X2 | $X_{2}$ | $k_{15}t2$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 1 |
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 8 | 2 |
/// | dir{P1} | 103 | P4cc | C4v^5 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 16 | 2 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 1 |
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 8 | 2 |
/// | dir{P1} | 104 | P4nc | C4v^6 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 16 | 2 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{4}$ | $k_{18}t2$ | 2 | B8a | no |
/// | M2 | $M_{2}$ | $k_{18}t1$ | 2 | B8a | no |
/// | M3 | $M_{1}$ | $k_{18}t4$ | 2 | B8a | no |
/// | M4 | $M_{3}$ | $k_{18}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 4 | 1 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 4 | 1 |
/// | dir{P1} | 50 | Pban | D2h^4 | 4 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 103 | P4cc | C4v^5 | 4 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 104 | P4nc | C4v^6 | 4 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 4 | 1 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{16}t1t2$ | 8 | F64a | no |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 1 |
/// | dir{P3} | 82 | I-4 | S4^2 | 16 | 2 |
/// | dir{C1} | 79 | I4 | C4^5 | 16 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 16 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 16 | 2 |
/// | dir{P3} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P3} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{3}A4$ | $k_{20}t1t2$ | 4 | D16c | no |
/// | A3A4 | $A_{1}A2$ | $k_{20}t3t4$ | 4 | D16c | yes |
///
/// ### Isotropy subgroups for A1A2 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
///
/// ### Isotropy subgroups for A3A4 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 79 | I4 | C4^5 | 8 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{2}$ | $k_{19}t2$ | 2 | B8a | yes |
/// | Z2 | $Z_{1}$ | $k_{19}t1$ | 2 | B8a | yes |
/// | Z3Z4 | $Z_{3}Z4$ | $k_{19}t3t4$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 4 | 1 |
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 4 | 1 |
/// | dir{P1} | 77 | P4_2 | C4^3 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 4 | 1 |
/// | dir{P1} | 85 | P4/n | C4h^3 | 4 | 1 |
/// | dir{P1} | 75 | P4 | C4^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 1 |
///
pub struct Sg130;

/// # 131 P4₂/mmc (D₄ₕ⁹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 47 | Pmmm | D2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 105 | P4_2mc | C4v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 66 | Cccm | D2h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 4 | 1 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 4 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{15}t1$ | 2 | B8a | no |
/// | X1- | $X_{1}^-$ | $k_{15}t2$ | 2 | B8a | no |
/// | X2+ | $X_{3}^+$ | $k_{15}t7$ | 2 | B8a | no |
/// | X2- | $X_{3}^-$ | $k_{15}t8$ | 2 | B8a | no |
/// | X3+ | $X_{4}^+$ | $k_{15}t3$ | 2 | B8a | no |
/// | X3- | $X_{4}^-$ | $k_{15}t4$ | 2 | B8a | no |
/// | X4+ | $X_{2}^+$ | $k_{15}t5$ | 2 | B8a | no |
/// | X4- | $X_{2}^-$ | $k_{15}t6$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 47 | Pmmm | D2h^1 | 4 | 1 |
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 4 | 2 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
/// | dir{P1} | 133 | P4_2/nbc | D4h^11 | 4 | 2 |
/// | dir{P1} | 50 | Pban | D2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 135 | P4_2/mbc | D4h^13 | 4 | 2 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 8 | 2 |
///
/// ### Isotropy subgroups for X2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 137 | P4_2/nmc | D4h^15 | 4 | 2 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 8 | 2 |
///
/// ### Isotropy subgroups for X3+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 137 | P4_2/nmc | D4h^15 | 4 | 2 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 8 | 2 |
///
/// ### Isotropy subgroups for X3- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 135 | P4_2/mbc | D4h^13 | 4 | 2 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 8 | 2 |
///
/// ### Isotropy subgroups for X4+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
/// | dir{P1} | 133 | P4_2/nbc | D4h^11 | 4 | 2 |
/// | dir{P1} | 50 | Pban | D2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X4- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 47 | Pmmm | D2h^1 | 4 | 1 |
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 4 | 2 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{18}t1$ | 1 | A2a | no |
/// | M1- | $M_{1}^-$ | $k_{18}t2$ | 1 | A2a | no |
/// | M2+ | $M_{3}^+$ | $k_{18}t5$ | 1 | A2a | no |
/// | M2- | $M_{3}^-$ | $k_{18}t6$ | 1 | A2a | no |
/// | M3+ | $M_{2}^+$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3- | $M_{2}^-$ | $k_{18}t4$ | 1 | A2a | no |
/// | M4+ | $M_{4}^+$ | $k_{18}t7$ | 1 | A2a | no |
/// | M4- | $M_{4}^-$ | $k_{18}t8$ | 1 | A2a | no |
/// | M5+ | $M_{5}^+$ | $k_{18}t9$ | 2 | B8a | no |
/// | M5- | $M_{5}^-$ | $k_{18}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 2 | 1 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for M2+ ($M_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for M2- ($M_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 138 | P4_2/ncm | D4h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for M3+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for M3- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 138 | P4_2/ncm | D4h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for M4+ ($M_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 2 | 1 |
///
/// ### Isotropy subgroups for M4- ($M_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for M5+ ($M_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 4 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
///
/// ### Isotropy subgroups for M5- ($M_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 4 | 1 |
/// | dir{P1} | 53 | Pmna | D2h^7 | 4 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{16}t1$ | 2 | B8a | no |
/// | R1- | $R_{1}^-$ | $k_{16}t2$ | 2 | B8a | no |
/// | R2+ | $R_{3}^+$ | $k_{16}t7$ | 2 | B8a | no |
/// | R2- | $R_{3}^-$ | $k_{16}t8$ | 2 | B8a | no |
/// | R3+ | $R_{4}^+$ | $k_{16}t3$ | 2 | B8a | no |
/// | R3- | $R_{4}^-$ | $k_{16}t4$ | 2 | B8a | no |
/// | R4+ | $R_{2}^+$ | $k_{16}t5$ | 2 | B8a | no |
/// | R4- | $R_{2}^-$ | $k_{16}t6$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 4 | 1 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 4 | 2 |
/// | dir{P1} | 74 | Imma | D2h^28 | 8 | 2 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 4 | 1 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 4 | 2 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 8 | 2 |
///
/// ### Isotropy subgroups for R2+ ($R_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 4 | 1 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 4 | 2 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 8 | 2 |
///
/// ### Isotropy subgroups for R2- ($R_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 4 | 1 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 4 | 2 |
/// | dir{P1} | 74 | Imma | D2h^28 | 8 | 2 |
///
/// ### Isotropy subgroups for R3+ ($R_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 4 | 1 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 4 | 2 |
/// | dir{P1} | 74 | Imma | D2h^28 | 8 | 2 |
///
/// ### Isotropy subgroups for R3- ($R_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 4 | 1 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 4 | 2 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 8 | 2 |
///
/// ### Isotropy subgroups for R4+ ($R_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 4 | 1 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 4 | 2 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 8 | 2 |
///
/// ### Isotropy subgroups for R4- ($R_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 4 | 1 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 4 | 2 |
/// | dir{P1} | 74 | Imma | D2h^28 | 8 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{4}$ | $k_{20}t3$ | 2 | B8a | yes |
/// | A2 | $A_{2}$ | $k_{20}t4$ | 2 | B8a | yes |
/// | A3 | $A_{3}$ | $k_{20}t1$ | 2 | B8a | no |
/// | A4 | $A_{1}$ | $k_{20}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 4 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 4 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 4 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 4 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 8 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 4 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 4 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 4 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 4 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{4}$ | $k_{19}t3$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{19}t4$ | 2 | B8a | yes |
/// | Z3 | $Z_{3}$ | $k_{19}t1$ | 2 | B8a | no |
/// | Z4 | $Z_{1}$ | $k_{19}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 4 | 1 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 4 | 1 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 4 | 1 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 4 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 4 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 1 |
///
pub struct Sg131;

/// # 132 P4₂/mcm (D₄ₕ¹⁰)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 101 | P4_2cm | C4v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{15}t1$ | 2 | B8a | no |
/// | X1- | $X_{1}^-$ | $k_{15}t2$ | 2 | B8a | no |
/// | X2+ | $X_{3}^+$ | $k_{15}t7$ | 2 | B8a | no |
/// | X2- | $X_{3}^-$ | $k_{15}t8$ | 2 | B8a | no |
/// | X3+ | $X_{4}^+$ | $k_{15}t3$ | 2 | B8a | no |
/// | X3- | $X_{4}^-$ | $k_{15}t4$ | 2 | B8a | no |
/// | X4+ | $X_{2}^+$ | $k_{15}t5$ | 2 | B8a | no |
/// | X4- | $X_{2}^-$ | $k_{15}t6$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 4 | 2 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 50 | Pban | D2h^4 | 4 | 1 |
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 4 | 2 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 8 | 2 |
///
/// ### Isotropy subgroups for X2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 53 | Pmna | D2h^7 | 4 | 1 |
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 4 | 2 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 8 | 2 |
///
/// ### Isotropy subgroups for X2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 138 | P4_2/ncm | D4h^16 | 4 | 2 |
/// | dir{P1} | 56 | Pccn | D2h^10 | 8 | 2 |
///
/// ### Isotropy subgroups for X3+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 138 | P4_2/ncm | D4h^16 | 4 | 2 |
/// | dir{P1} | 56 | Pccn | D2h^10 | 8 | 2 |
///
/// ### Isotropy subgroups for X3- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 53 | Pmna | D2h^7 | 4 | 1 |
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 4 | 2 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 8 | 2 |
///
/// ### Isotropy subgroups for X4+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 50 | Pban | D2h^4 | 4 | 1 |
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 4 | 2 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 8 | 2 |
///
/// ### Isotropy subgroups for X4- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 4 | 2 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{18}t1$ | 1 | A2a | no |
/// | M1- | $M_{1}^-$ | $k_{18}t2$ | 1 | A2a | no |
/// | M2+ | $M_{3}^+$ | $k_{18}t5$ | 1 | A2a | no |
/// | M2- | $M_{3}^-$ | $k_{18}t6$ | 1 | A2a | no |
/// | M3+ | $M_{2}^+$ | $k_{18}t3$ | 1 | A2a | no |
/// | M3- | $M_{2}^-$ | $k_{18}t4$ | 1 | A2a | no |
/// | M4+ | $M_{4}^+$ | $k_{18}t7$ | 1 | A2a | no |
/// | M4- | $M_{4}^-$ | $k_{18}t8$ | 1 | A2a | no |
/// | M5+ | $M_{5}^+$ | $k_{18}t9$ | 2 | B8a | no |
/// | M5- | $M_{5}^-$ | $k_{18}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 2 | 1 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 133 | P4_2/nbc | D4h^11 | 2 | 1 |
///
/// ### Isotropy subgroups for M2+ ($M_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 135 | P4_2/mbc | D4h^13 | 2 | 1 |
///
/// ### Isotropy subgroups for M2- ($M_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 137 | P4_2/nmc | D4h^15 | 2 | 1 |
///
/// ### Isotropy subgroups for M3+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 135 | P4_2/mbc | D4h^13 | 2 | 1 |
///
/// ### Isotropy subgroups for M3- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 137 | P4_2/nmc | D4h^15 | 2 | 1 |
///
/// ### Isotropy subgroups for M4+ ($M_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 2 | 1 |
///
/// ### Isotropy subgroups for M4- ($M_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 133 | P4_2/nbc | D4h^11 | 2 | 1 |
///
/// ### Isotropy subgroups for M5+ ($M_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 68 | Ccca | D2h^22 | 4 | 1 |
/// | dir{P1} | 53 | Pmna | D2h^7 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
///
/// ### Isotropy subgroups for M5- ($M_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 66 | Cccm | D2h^20 | 4 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D32c | yes |
/// | R2 | $R_{2}$ | $k_{16}t2$ | 4 | D32c | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 8 | 2 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 8 | 2 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 8 | 2 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 8 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 16 | 2 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 16 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 8 | 2 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 8 | 2 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 8 | 2 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 8 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 16 | 2 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 16 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{4}$ | $k_{20}t3$ | 2 | B8a | yes |
/// | A2 | $A_{2}$ | $k_{20}t4$ | 2 | B8a | yes |
/// | A3 | $A_{3}$ | $k_{20}t2$ | 2 | B8a | no |
/// | A4 | $A_{1}$ | $k_{20}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 4 | 1 |
/// | dir{P1} | 71 | Immm | D2h^25 | 4 | 1 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 4 | 1 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 4 | 1 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 8 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 4 | 1 |
/// | dir{P1} | 74 | Imma | D2h^28 | 4 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 4 | 1 |
/// | dir{P1} | 74 | Imma | D2h^28 | 4 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{4}$ | $k_{19}t3$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{19}t4$ | 2 | B8a | yes |
/// | Z3 | $Z_{3}$ | $k_{19}t2$ | 2 | B8a | no |
/// | Z4 | $Z_{1}$ | $k_{19}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 4 | 1 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 4 | 1 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 4 | 1 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 4 | 1 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 8 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 4 | 1 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 4 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 4 | 1 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 4 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 1 |
///
pub struct Sg132;

/// # 133 P4₂/nbc (D₄ₕ¹¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 133 | P4_2/nbc | D4h^11 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 50 | Pban | D2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 106 | P4_2bc | C4v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 68 | Ccca | D2h^22 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 4 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D64c | yes |
/// | X2 | $X_{2}$ | $k_{15}t2$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 1 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 8 | 2 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 16 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 1 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 8 | 2 |
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 16 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{4}$ | $k_{18}t4$ | 2 | B8a | yes |
/// | M2 | $M_{2}$ | $k_{18}t3$ | 2 | B8a | yes |
/// | M3 | $M_{1}$ | $k_{18}t2$ | 2 | B8a | no |
/// | M4 | $M_{3}$ | $k_{18}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 4 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 4 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 4 | 1 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 4 | 1 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 4 | 1 |
/// | dir{P1} | 50 | Pban | D2h^4 | 4 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{2}$ | $k_{16}t2$ | 4 | D64c | yes |
/// | R2 | $R_{1}$ | $k_{16}t1$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for R1 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 8 | 2 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 8 | 2 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | A2 | $A_{2}$ | $k_{20}t2$ | 2 | B8a | yes |
/// | A3A4 | $A_{3}A4$ | $k_{20}t3t4$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 4 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 4 | 1 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 4 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 4 | 1 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 8 | 1 |
///
/// ### Isotropy subgroups for A3A4 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 8 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{4}$ | $k_{19}t3$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{19}t4$ | 2 | B8a | yes |
/// | Z3 | $Z_{1}$ | $k_{19}t2$ | 2 | B8a | no |
/// | Z4 | $Z_{3}$ | $k_{19}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 4 | 1 |
/// | dir{P1} | 50 | Pban | D2h^4 | 4 | 1 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 4 | 1 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 4 | 1 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 4 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 4 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 1 |
///
pub struct Sg133;

/// # 134 P4₂/nnm (D₄ₕ¹²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 48 | Pnnn | D2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 102 | P4_2nm | C4v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 4 | 1 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D64c | yes |
/// | X2 | $X_{2}$ | $k_{15}t2$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 1 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 8 | 2 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 16 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 1 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 8 | 2 |
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 16 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{4}$ | $k_{18}t4$ | 2 | B8a | yes |
/// | M2 | $M_{2}$ | $k_{18}t3$ | 2 | B8a | yes |
/// | M3 | $M_{1}$ | $k_{18}t2$ | 2 | B8a | no |
/// | M4 | $M_{3}$ | $k_{18}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 4 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 4 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 4 | 1 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 4 | 1 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t2$ | 4 | D64c | yes |
/// | R2 | $R_{2}$ | $k_{16}t1$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 1 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 8 | 2 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 1 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 8 | 2 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1+ | $A_{1}^+$ | $k_{20}t3$ | 1 | A2a | no |
/// | A1- | $A_{1}^-$ | $k_{20}t4$ | 1 | A2a | no |
/// | A2+ | $A_{3}^+$ | $k_{20}t7$ | 1 | A2a | no |
/// | A2- | $A_{3}^-$ | $k_{20}t8$ | 1 | A2a | no |
/// | A3+ | $A_{2}^+$ | $k_{20}t1$ | 1 | A2a | no |
/// | A3- | $A_{2}^-$ | $k_{20}t2$ | 1 | A2a | no |
/// | A4+ | $A_{4}^+$ | $k_{20}t5$ | 1 | A2a | no |
/// | A4- | $A_{4}^-$ | $k_{20}t6$ | 1 | A2a | no |
/// | A5+ | $A_{5}^+$ | $k_{20}t9$ | 2 | B8a | no |
/// | A5- | $A_{5}^-$ | $k_{20}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1+ ($A_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for A1- ($A_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for A2+ ($A_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for A2- ($A_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for A3+ ($A_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for A3- ($A_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for A4+ ($A_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for A4- ($A_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for A5+ ($A_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 72 | Ibam | D2h^26 | 4 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
///
/// ### Isotropy subgroups for A5- ($A_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 72 | Ibam | D2h^26 | 4 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{4}$ | $k_{19}t3$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{19}t4$ | 2 | B8a | yes |
/// | Z3 | $Z_{1}$ | $k_{19}t2$ | 2 | B8a | no |
/// | Z4 | $Z_{3}$ | $k_{19}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 4 | 1 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 4 | 1 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 4 | 1 |
/// | dir{P1} | 68 | Ccca | D2h^22 | 4 | 1 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 8 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 4 | 1 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 4 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 4 | 1 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 4 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 1 |
///
pub struct Sg134;

/// # 135 P4₂/mbc (D₄ₕ¹³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 135 | P4_2/mbc | D4h^13 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 55 | Pbam | D2h^9 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 106 | P4_2bc | C4v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 66 | Cccm | D2h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 4 | 1 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 4 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{2}$ | $k_{15}t2$ | 4 | D32b | yes |
/// | X2 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 8 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 1 |
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 8 | 2 |
/// | dir{P1} | 68 | Ccca | D2h^22 | 8 | 2 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 16 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 16 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 16 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 32 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 8 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 16 | 1 |
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 8 | 2 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 8 | 2 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 16 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 16 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 16 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+M4+ | $M_{1}^-M2-$ | $k_{18}t1t7$ | 2 | B4a | no |
/// | M1-M4- | $M_{3}^+M4+$ | $k_{18}t2t8$ | 2 | B4a | no |
/// | M2+M3+ | $M_{3}^-M4-$ | $k_{18}t3t5$ | 2 | B4a | no |
/// | M2-M3- | $M_{1}^+M2+$ | $k_{18}t4t6$ | 2 | B4a | no |
/// | M5+ | $M_{5}^-$ | $k_{18}t9$ | 2 | B8a | no |
/// | M5- | $M_{5}^+$ | $k_{18}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1+M4+ ($M_{1}^-M2-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 56 | Pccn | D2h^10 | 4 | 1 |
///
/// ### Isotropy subgroups for M1-M4- ($M_{3}^+M4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 58 | Pnnm | D2h^12 | 4 | 1 |
///
/// ### Isotropy subgroups for M2+M3+ ($M_{3}^-M4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 48 | Pnnn | D2h^2 | 4 | 1 |
///
/// ### Isotropy subgroups for M2-M3- ($M_{1}^+M2+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
///
/// ### Isotropy subgroups for M5+ ($M_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 4 | 1 |
/// | dir{P1} | 53 | Pmna | D2h^7 | 4 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M5- ($M_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 4 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{2}$ | $k_{16}t2$ | 4 | D32b | yes |
/// | R2 | $R_{1}$ | $k_{16}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for R1 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 8 | 2 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 8 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 16 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 8 | 2 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 8 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{3}A4$ | $k_{20}t1t2$ | 4 | D16c | no |
/// | A3A4 | $A_{1}A2$ | $k_{20}t3t4$ | 4 | D16c | yes |
///
/// ### Isotropy subgroups for A1A2 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
///
/// ### Isotropy subgroups for A3A4 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 80 | I4_1 | C4^6 | 8 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{4}$ | $k_{19}t3$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{19}t4$ | 2 | B8a | yes |
/// | Z3 | $Z_{3}$ | $k_{19}t1$ | 2 | B8a | no |
/// | Z4 | $Z_{1}$ | $k_{19}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 4 | 1 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 4 | 1 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 4 | 1 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 4 | 1 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 4 | 1 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 4 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 8 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 4 | 1 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 4 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 8 | 1 |
///
pub struct Sg135;

/// # 136 P4₂/mnm (D₄ₕ¹⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 58 | Pnnm | D2h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 102 | P4_2nm | C4v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 4 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{2}$ | $k_{15}t2$ | 4 | D32b | yes |
/// | X2 | $X_{1}$ | $k_{15}t1$ | 4 | D32b | yes |
///
/// ### Isotropy subgroups for X1 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 8 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 1 |
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 8 | 2 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 8 | 2 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 16 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 16 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 16 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 32 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 8 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 16 | 1 |
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 8 | 2 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 8 | 2 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 16 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 16 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 16 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+M4+ | $M_{1}^+M2+$ | $k_{18}t4t6$ | 2 | B4a | no |
/// | M1-M4- | $M_{3}^-M4-$ | $k_{18}t3t5$ | 2 | B4a | no |
/// | M2+M3+ | $M_{3}^+M4+$ | $k_{18}t2t8$ | 2 | B4a | no |
/// | M2-M3- | $M_{1}^-M2-$ | $k_{18}t1t7$ | 2 | B4a | no |
/// | M5+ | $M_{5}^+$ | $k_{18}t10$ | 2 | B8a | no |
/// | M5- | $M_{5}^-$ | $k_{18}t9$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1+M4+ ($M_{1}^+M2+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 47 | Pmmm | D2h^1 | 4 | 1 |
///
/// ### Isotropy subgroups for M1-M4- ($M_{3}^-M4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 50 | Pban | D2h^4 | 4 | 1 |
///
/// ### Isotropy subgroups for M2+M3+ ($M_{3}^+M4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 55 | Pbam | D2h^9 | 4 | 1 |
///
/// ### Isotropy subgroups for M2-M3- ($M_{1}^-M2-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 59 | Pmmn | D2h^13 | 4 | 1 |
///
/// ### Isotropy subgroups for M5+ ($M_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 4 | 1 |
/// | dir{P1} | 53 | Pmna | D2h^7 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
///
/// ### Isotropy subgroups for M5- ($M_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 4 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{16}t2$ | 4 | D32b | no |
/// | R1- | $R_{1}^-$ | $k_{16}t1$ | 4 | D32b | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 8 | 2 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 8 | 2 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{20}t4$ | 2 | B8a | no |
/// | A2 | $A_{3}$ | $k_{20}t3$ | 2 | B8a | no |
/// | A3 | $A_{2}$ | $k_{20}t1$ | 2 | B8a | no |
/// | A4 | $A_{4}$ | $k_{20}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 109 | I4_1md | C4v^11 | 4 | 1 |
/// | dir{P1} | 71 | Immm | D2h^25 | 4 | 1 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 110 | I4_1cd | C4v^12 | 4 | 1 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 4 | 1 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 8 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 4 | 1 |
/// | dir{P1} | 74 | Imma | D2h^28 | 4 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 4 | 1 |
/// | dir{P1} | 74 | Imma | D2h^28 | 4 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{4}$ | $k_{19}t3$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{19}t4$ | 2 | B8a | yes |
/// | Z3 | $Z_{3}$ | $k_{19}t1$ | 2 | B8a | no |
/// | Z4 | $Z_{1}$ | $k_{19}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 4 | 1 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 4 | 1 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 4 | 1 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 4 | 1 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 8 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 4 | 1 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 4 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 4 | 1 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 4 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 1 |
///
pub struct Sg136;

/// # 137 P4₂/nmc (D₄ₕ¹⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 137 | P4_2/nmc | D4h^15 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 59 | Pmmn | D2h^13 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 105 | P4_2mc | C4v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 68 | Ccca | D2h^22 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 4 | 1 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 4 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D64c | yes |
/// | X2 | $X_{2}$ | $k_{15}t2$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 8 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 16 | 1 |
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 8 | 2 |
/// | dir{P1} | 105 | P4_2mc | C4v^7 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 16 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 1 |
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 8 | 2 |
/// | dir{P1} | 106 | P4_2bc | C4v^8 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 16 | 2 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{4}$ | $k_{18}t2$ | 2 | B8a | no |
/// | M2 | $M_{2}$ | $k_{18}t1$ | 2 | B8a | no |
/// | M3 | $M_{1}$ | $k_{18}t4$ | 2 | B8a | no |
/// | M4 | $M_{3}$ | $k_{18}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 4 | 1 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 4 | 1 |
/// | dir{P1} | 50 | Pban | D2h^4 | 4 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 101 | P4_2cm | C4v^3 | 4 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 102 | P4_2nm | C4v^4 | 4 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 4 | 1 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{16}t1$ | 4 | D64c | yes |
/// | R2 | $R_{2}$ | $k_{16}t2$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 8 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 8 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 16 | 1 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 8 | 2 |
/// | dir{P1} | 109 | I4_1md | C4v^11 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 16 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 39 | Abm2 | C2v^15 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 1 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 8 | 2 |
/// | dir{P1} | 110 | I4_1cd | C4v^12 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 16 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{2}$ | $k_{20}t4$ | 2 | B8a | no |
/// | A2 | $A_{1}$ | $k_{20}t3$ | 2 | B8a | no |
/// | A3A4 | $A_{3}A4$ | $k_{20}t1t2$ | 4 | D16c | yes |
///
/// ### Isotropy subgroups for A1 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 4 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 4 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 4 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 4 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 8 | 1 |
///
/// ### Isotropy subgroups for A3A4 ($A_{3}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 8 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 8 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{4}$ | $k_{19}t3$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{19}t4$ | 2 | B8a | yes |
/// | Z3 | $Z_{1}$ | $k_{19}t2$ | 2 | B8a | no |
/// | Z4 | $Z_{3}$ | $k_{19}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 4 | 1 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 4 | 1 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 4 | 1 |
/// | dir{P1} | 56 | Pccn | D2h^10 | 4 | 1 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 4 | 1 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 4 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 8 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 4 | 1 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 4 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 8 | 1 |
///
pub struct Sg137;

/// # 138 P4₂/ncm (D₄ₕ¹⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{17}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{17}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 138 | P4_2/ncm | D4h^16 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 56 | Pccn | D2h^10 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 101 | P4_2cm | C4v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 4 | 1 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{15}t1$ | 4 | D64c | yes |
/// | X2 | $X_{2}$ | $k_{15}t2$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 1 |
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 8 | 2 |
/// | dir{P1} | 101 | P4_2cm | C4v^3 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 16 | 2 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 16 | 1 |
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 8 | 2 |
/// | dir{P1} | 102 | P4_2nm | C4v^4 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 16 | 2 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{4}$ | $k_{18}t2$ | 2 | B8a | no |
/// | M2 | $M_{2}$ | $k_{18}t1$ | 2 | B8a | no |
/// | M3 | $M_{1}$ | $k_{18}t4$ | 2 | B8a | no |
/// | M4 | $M_{3}$ | $k_{18}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 4 | 1 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 4 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 4 | 1 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 4 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M3 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 105 | P4_2mc | C4v^7 | 4 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M4 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 106 | P4_2bc | C4v^8 | 4 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{16}t1t2$ | 8 | F64a | no |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 1 |
/// | dir{P3} | 82 | I-4 | S4^2 | 16 | 2 |
/// | dir{C1} | 80 | I4_1 | C4^6 | 16 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 16 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 16 | 2 |
/// | dir{P3} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P3} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1+A4+ | $A_{1}^-A2-$ | $k_{20}t3t5$ | 2 | B4a | no |
/// | A1-A4- | $A_{3}^+A4+$ | $k_{20}t4t6$ | 2 | B4a | no |
/// | A2+A3+ | $A_{3}^-A4-$ | $k_{20}t1t7$ | 2 | B4a | no |
/// | A2-A3- | $A_{1}^+A2+$ | $k_{20}t2t8$ | 2 | B4a | no |
/// | A5+ | $A_{5}^-$ | $k_{20}t9$ | 2 | B8a | no |
/// | A5- | $A_{5}^+$ | $k_{20}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1+A4+ ($A_{1}^-A2-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 74 | Imma | D2h^28 | 4 | 1 |
///
/// ### Isotropy subgroups for A1-A4- ($A_{3}^+A4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 73 | Ibca | D2h^27 | 4 | 1 |
///
/// ### Isotropy subgroups for A2+A3+ ($A_{3}^-A4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 73 | Ibca | D2h^27 | 4 | 1 |
///
/// ### Isotropy subgroups for A2-A3- ($A_{1}^+A2+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 74 | Imma | D2h^28 | 4 | 1 |
///
/// ### Isotropy subgroups for A5+ ($A_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 4 | 1 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
///
/// ### Isotropy subgroups for A5- ($A_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 4 | 1 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{4}$ | $k_{19}t3$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{19}t4$ | 2 | B8a | yes |
/// | Z3 | $Z_{1}$ | $k_{19}t2$ | 2 | B8a | no |
/// | Z4 | $Z_{3}$ | $k_{19}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 4 | 1 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 4 | 1 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 4 | 1 |
/// | dir{P1} | 68 | Ccca | D2h^22 | 4 | 1 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 8 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 4 | 1 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 4 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 4 | 1 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 4 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 1 |
///
pub struct Sg138;

/// # 139 I4/mmm (D₄ₕ¹⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{14}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{14}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{14}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{14}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{14}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{14}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 71 | Immm | D2h^25 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 87 | I4/m | C4h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 107 | I4mm | C4v^9 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 69 | Fmmm | D2h^23 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 4 | 1 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 4 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{13}t1$ | 2 | B8a | no |
/// | X1- | $X_{1}^-$ | $k_{13}t2$ | 2 | B8a | no |
/// | X2+ | $X_{3}^+$ | $k_{13}t7$ | 2 | B8a | no |
/// | X2- | $X_{3}^-$ | $k_{13}t8$ | 2 | B8a | no |
/// | X3+ | $X_{4}^+$ | $k_{13}t3$ | 2 | B8a | no |
/// | X3- | $X_{4}^-$ | $k_{13}t4$ | 2 | B8a | no |
/// | X4+ | $X_{2}^+$ | $k_{13}t5$ | 2 | B8a | no |
/// | X4- | $X_{2}^-$ | $k_{13}t6$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 4 | 1 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 4 | 2 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 8 | 2 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 68 | Ccca | D2h^22 | 4 | 1 |
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 4 | 2 |
/// | dir{P1} | 50 | Pban | D2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 4 | 1 |
/// | dir{P1} | 127 | P4/mbm | D4h^5 | 4 | 2 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 8 | 2 |
///
/// ### Isotropy subgroups for X2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 4 | 1 |
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 4 | 2 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 8 | 2 |
///
/// ### Isotropy subgroups for X3+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 4 | 1 |
/// | dir{P1} | 138 | P4_2/ncm | D4h^16 | 4 | 2 |
/// | dir{P1} | 56 | Pccn | D2h^10 | 8 | 2 |
///
/// ### Isotropy subgroups for X3- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 4 | 1 |
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 4 | 2 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 8 | 2 |
///
/// ### Isotropy subgroups for X4+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 66 | Cccm | D2h^20 | 4 | 1 |
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 4 | 2 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 8 | 2 |
///
/// ### Isotropy subgroups for X4- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 4 | 1 |
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 4 | 2 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $Z_{1}^+$ | $k_{15}t1$ | 1 | A2a | no |
/// | M1- | $Z_{1}^-$ | $k_{15}t2$ | 1 | A2a | no |
/// | M2+ | $Z_{3}^+$ | $k_{15}t5$ | 1 | A2a | no |
/// | M2- | $Z_{3}^-$ | $k_{15}t6$ | 1 | A2a | no |
/// | M3+ | $Z_{2}^+$ | $k_{15}t3$ | 1 | A2a | no |
/// | M3- | $Z_{2}^-$ | $k_{15}t4$ | 1 | A2a | no |
/// | M4+ | $Z_{4}^+$ | $k_{15}t7$ | 1 | A2a | no |
/// | M4- | $Z_{4}^-$ | $k_{15}t8$ | 1 | A2a | no |
/// | M5+ | $Z_{5}^+$ | $k_{15}t9$ | 2 | B8a | no |
/// | M5- | $Z_{5}^-$ | $k_{15}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for M1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for M2+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 2 | 1 |
///
/// ### Isotropy subgroups for M2- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for M3+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 128 | P4/mnc | D4h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M3- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for M4+ ($Z_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for M4- ($Z_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 137 | P4_2/nmc | D4h^15 | 2 | 1 |
///
/// ### Isotropy subgroups for M5+ ($Z_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 4 | 1 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 4 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
///
/// ### Isotropy subgroups for M5- ($Z_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 4 | 1 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 4 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1+ | $N_{1}^+$ | $k_{11}t1$ | 4 | D64c | no |
/// | N1- | $N_{1}^-$ | $k_{11}t2$ | 4 | D64c | no |
/// | N2+ | $N_{2}^+$ | $k_{11}t3$ | 4 | D64c | no |
/// | N2- | $N_{2}^-$ | $k_{11}t4$ | 4 | D64c | no |
///
/// ### Isotropy subgroups for N1+ ($N_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 8 | 4 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 8 | 4 |
/// | dir{P1} | 74 | Imma | D2h^28 | 16 | 4 |
/// | dir{P1} | 71 | Immm | D2h^25 | 16 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 32 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 32 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 32 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 64 | 4 |
///
/// ### Isotropy subgroups for N1- ($N_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 8 | 4 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 8 | 4 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 16 | 4 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 16 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 32 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 32 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 32 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 64 | 4 |
///
/// ### Isotropy subgroups for N2+ ($N_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 8 | 4 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 8 | 4 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 16 | 4 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 16 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 32 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 32 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 32 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 64 | 4 |
///
/// ### Isotropy subgroups for N2- ($N_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 8 | 4 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 8 | 4 |
/// | dir{P1} | 74 | Imma | D2h^28 | 16 | 4 |
/// | dir{P1} | 71 | Immm | D2h^25 | 16 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 32 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 32 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 32 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 64 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1 | $P_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | P2 | $P_{3}$ | $k_{12}t4$ | 2 | B8a | no |
/// | P3 | $P_{4}$ | $k_{12}t3$ | 2 | B8a | no |
/// | P4 | $P_{2}$ | $k_{12}t2$ | 2 | B8a | no |
/// | P5 | $P_{5}$ | $k_{12}t5$ | 4 | D32c | yes |
///
/// ### Isotropy subgroups for P1 ($P_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 4 | 2 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 4 | 2 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 8 | 2 |
///
/// ### Isotropy subgroups for P2 ($P_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 4 | 2 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 4 | 2 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 8 | 2 |
///
/// ### Isotropy subgroups for P3 ($P_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 4 | 2 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 4 | 2 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 8 | 2 |
///
/// ### Isotropy subgroups for P4 ($P_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 4 | 2 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 4 | 2 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 8 | 2 |
///
/// ### Isotropy subgroups for P5 ($P_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 8 | 2 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 8 | 2 |
/// | dir{P1} | 74 | Imma | D2h^28 | 8 | 2 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 8 | 2 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 8 | 2 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 8 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 16 | 2 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 16 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 16 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
///
pub struct Sg139;

/// # 140 I4/mcm (D₄ₕ¹⁸)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{14}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{14}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{14}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{14}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{14}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{14}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 72 | Ibam | D2h^26 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 87 | I4/m | C4h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 108 | I4cm | C4v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 69 | Fmmm | D2h^23 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 4 | 1 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 4 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{2}^-$ | $k_{13}t1$ | 2 | B8a | no |
/// | X1- | $X_{2}^+$ | $k_{13}t2$ | 2 | B8a | no |
/// | X2+ | $X_{4}^-$ | $k_{13}t7$ | 2 | B8a | no |
/// | X2- | $X_{4}^+$ | $k_{13}t8$ | 2 | B8a | no |
/// | X3+ | $X_{3}^-$ | $k_{13}t3$ | 2 | B8a | no |
/// | X3- | $X_{3}^+$ | $k_{13}t4$ | 2 | B8a | no |
/// | X4+ | $X_{1}^-$ | $k_{13}t5$ | 2 | B8a | no |
/// | X4- | $X_{1}^+$ | $k_{13}t6$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 4 | 1 |
/// | dir{P1} | 124 | P4/mcc | D4h^2 | 4 | 2 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 8 | 2 |
///
/// ### Isotropy subgroups for X1- ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 66 | Cccm | D2h^20 | 4 | 1 |
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 4 | 2 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 8 | 2 |
///
/// ### Isotropy subgroups for X2+ ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 4 | 1 |
/// | dir{P1} | 128 | P4/mnc | D4h^6 | 4 | 2 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 8 | 2 |
///
/// ### Isotropy subgroups for X2- ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 4 | 1 |
/// | dir{P1} | 130 | P4/ncc | D4h^8 | 4 | 2 |
/// | dir{P1} | 56 | Pccn | D2h^10 | 8 | 2 |
///
/// ### Isotropy subgroups for X3+ ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 4 | 1 |
/// | dir{P1} | 137 | P4_2/nmc | D4h^15 | 4 | 2 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 8 | 2 |
///
/// ### Isotropy subgroups for X3- ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 4 | 1 |
/// | dir{P1} | 135 | P4_2/mbc | D4h^13 | 4 | 2 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 8 | 2 |
///
/// ### Isotropy subgroups for X4+ ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 68 | Ccca | D2h^22 | 4 | 1 |
/// | dir{P1} | 133 | P4_2/nbc | D4h^11 | 4 | 2 |
/// | dir{P1} | 50 | Pban | D2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for X4- ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 4 | 1 |
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 4 | 2 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $Z_{3}^+$ | $k_{15}t2$ | 1 | A2a | no |
/// | M1- | $Z_{3}^-$ | $k_{15}t1$ | 1 | A2a | no |
/// | M2+ | $Z_{1}^+$ | $k_{15}t6$ | 1 | A2a | no |
/// | M2- | $Z_{1}^-$ | $k_{15}t5$ | 1 | A2a | no |
/// | M3+ | $Z_{4}^+$ | $k_{15}t4$ | 1 | A2a | no |
/// | M3- | $Z_{4}^-$ | $k_{15}t3$ | 1 | A2a | no |
/// | M4+ | $Z_{2}^+$ | $k_{15}t8$ | 1 | A2a | no |
/// | M4- | $Z_{2}^-$ | $k_{15}t7$ | 1 | A2a | no |
/// | M5+ | $Z_{5}^+$ | $k_{15}t10$ | 2 | B8a | no |
/// | M5- | $Z_{5}^-$ | $k_{15}t9$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 124 | P4/mcc | D4h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for M1- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M2+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 2 | 1 |
///
/// ### Isotropy subgroups for M2- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 133 | P4_2/nbc | D4h^11 | 2 | 1 |
///
/// ### Isotropy subgroups for M3+ ($Z_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 127 | P4/mbm | D4h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for M3- ($Z_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 130 | P4/ncc | D4h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for M4+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 135 | P4_2/mbc | D4h^13 | 2 | 1 |
///
/// ### Isotropy subgroups for M4- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 138 | P4_2/ncm | D4h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for M5+ ($Z_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 4 | 1 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 4 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
///
/// ### Isotropy subgroups for M5- ($Z_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 4 | 1 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 4 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{11}t1$ | 8 | F128a | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 1 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 16 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P3} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P3} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 32 | 2 |
/// | dir{P3} | 3 | P2 | C2^1 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{C1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{C1} | 1 | P1 | C1^1 | 64 | 2 |
/// | dir{C1} | 1 | P1 | C1^1 | 64 | 2 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 16 | 4 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 16 | 4 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 16 | 4 |
/// | dir{P1} | 97 | I422 | D4^9 | 16 | 4 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 16 | 4 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 16 | 4 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 16 | 4 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 16 | 4 |
/// | dir{P1} | 82 | I-4 | S4^2 | 32 | 4 |
/// | dir{P1} | 82 | I-4 | S4^2 | 32 | 4 |
/// | dir{C1} | 80 | I4_1 | C4^6 | 32 | 4 |
/// | dir{P1} | 79 | I4 | C4^5 | 32 | 4 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 32 | 4 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 32 | 4 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 32 | 4 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 32 | 4 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 32 | 4 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 32 | 4 |
/// | dir{P1} | 23 | I222 | D2^8 | 32 | 4 |
/// | dir{P1} | 23 | I222 | D2^8 | 32 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 32 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 32 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 32 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 32 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 32 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 32 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 32 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 32 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 32 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 32 | 4 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 32 | 4 |
/// | dir{P3} | 12 | C2/m | C2h^3 | 32 | 4 |
/// | dir{P1} | 9 | Cc | Cs^4 | 64 | 4 |
/// | dir{P1} | 9 | Cc | Cs^4 | 64 | 4 |
/// | dir{P1} | 8 | Cm | Cs^3 | 64 | 4 |
/// | dir{P3} | 8 | Cm | Cs^3 | 64 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 64 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 64 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 64 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 64 | 4 |
/// | dir{P3} | 5 | C2 | C2^3 | 64 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 64 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 64 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 64 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 128 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P3 | $P_{2}P3$ | $k_{12}t2t4$ | 4 | D16c | no |
/// | P2P4 | $P_{1}P4$ | $k_{12}t1t3$ | 4 | D16c | no |
/// | P5 | $P_{5}$ | $k_{12}t5$ | 4 | D32c | yes |
///
/// ### Isotropy subgroups for P1P3 ($P_{2}P3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 8 | 2 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 8 | 2 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 8 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 16 | 2 |
///
/// ### Isotropy subgroups for P2P4 ($P_{1}P4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 8 | 2 |
/// | dir{P1} | 74 | Imma | D2h^28 | 8 | 2 |
/// | dir{P1} | 71 | Immm | D2h^25 | 8 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 16 | 2 |
///
/// ### Isotropy subgroups for P5 ($P_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 8 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 8 | 2 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 8 | 2 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 8 | 2 |
/// | dir{P1} | 74 | Imma | D2h^28 | 8 | 2 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 8 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 16 | 2 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 16 | 2 |
/// | dir{P1} | 79 | I4 | C4^5 | 16 | 2 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
///
pub struct Sg140;

/// # 141 I4₁/amd (D₄ₕ¹⁹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{14}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{14}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{14}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{14}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{14}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{14}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 74 | Imma | D2h^28 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 109 | I4_1md | C4v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 70 | Fddd | D2h^24 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 4 | 1 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 4 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 4 | D64c | yes |
/// | X2 | $X_{2}$ | $k_{13}t2$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 1 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 8 | 2 |
/// | dir{P1} | 91 | P4_122 | D4^3 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 16 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 1 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 8 | 2 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 8 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 16 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 16 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $Z_{2}$ | $k_{15}t3$ | 2 | B8a | yes |
/// | M2 | $Z_{4}$ | $k_{15}t4$ | 2 | B8a | yes |
/// | M3 | $Z_{3}$ | $k_{15}t1$ | 2 | B8a | no |
/// | M4 | $Z_{1}$ | $k_{15}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 4 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 4 | 1 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 8 | 1 |
///
/// ### Isotropy subgroups for M2 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 4 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 4 | 1 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 8 | 1 |
///
/// ### Isotropy subgroups for M3 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 4 | 1 |
/// | dir{P1} | 53 | Pmna | D2h^7 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 1 |
///
/// ### Isotropy subgroups for M4 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 4 | 1 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 4 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1+ | $N_{1}^+$ | $k_{11}t1$ | 4 | D128a | no |
/// | N1- | $N_{1}^-$ | $k_{11}t2$ | 4 | D128a | no |
/// | N2+ | $N_{2}^+$ | $k_{11}t3$ | 4 | D128a | no |
/// | N2- | $N_{2}^-$ | $k_{11}t4$ | 4 | D128a | no |
///
/// ### Isotropy subgroups for N1+ ($N_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 16 | 4 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 32 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 64 | 4 |
/// | dir{C1} | 8 | Cm | Cs^3 | 64 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 64 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 128 | 4 |
///
/// ### Isotropy subgroups for N1- ($N_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 16 | 4 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 32 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 64 | 4 |
/// | dir{C1} | 9 | Cc | Cs^4 | 64 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 64 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 128 | 4 |
///
/// ### Isotropy subgroups for N2+ ($N_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 8 | 1 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 16 | 4 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 32 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 64 | 4 |
/// | dir{C1} | 9 | Cc | Cs^4 | 64 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 64 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 128 | 4 |
///
/// ### Isotropy subgroups for N2- ($N_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 8 | 1 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 16 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 16 | 4 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 32 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 64 | 4 |
/// | dir{C1} | 8 | Cm | Cs^3 | 64 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 64 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 128 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1 | $P_{1}$ | $k_{12}t1$ | 4 | D64d | yes |
/// | P2 | $P_{2}$ | $k_{12}t2$ | 4 | D64d | yes |
///
/// ### Isotropy subgroups for P1 ($P_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 8 | 2 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for P2 ($P_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 8 | 2 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 16 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
pub struct Sg141;

/// # 142 I4₁/acd (D₄ₕ²⁰)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{14}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{14}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{14}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{14}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{14}t9$ | 2 | B8a | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{14}t10$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 73 | Ibca | D2h^27 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 110 | I4_1cd | C4v^12 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 70 | Fddd | D2h^24 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 4 | 1 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 4 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 8 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{13}t1$ | 4 | D64c | yes |
/// | X2 | $X_{2}$ | $k_{13}t2$ | 4 | D64c | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 8 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 8 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 16 | 1 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 8 | 2 |
/// | dir{P1} | 91 | P4_122 | D4^3 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 16 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 8 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 8 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 16 | 1 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 8 | 2 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 8 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 16 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $Z_{4}$ | $k_{15}t4$ | 2 | B8a | yes |
/// | M2 | $Z_{2}$ | $k_{15}t3$ | 2 | B8a | yes |
/// | M3 | $Z_{1}$ | $k_{15}t2$ | 2 | B8a | no |
/// | M4 | $Z_{3}$ | $k_{15}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 4 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 8 | 1 |
///
/// ### Isotropy subgroups for M2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 4 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 8 | 1 |
///
/// ### Isotropy subgroups for M3 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 4 | 1 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 4 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 8 | 1 |
///
/// ### Isotropy subgroups for M4 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 4 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 4 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 8 | 1 |
///
/// ## Irreps at $N$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{11}t1$ | 8 | F128b | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 16 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 16 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 32 | 1 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 16 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 16 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 32 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
/// | dir{P3} | 82 | I-4 | S4^2 | 32 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 32 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 64 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 64 | 4 |
/// | dir{P3} | 5 | C2 | C2^3 | 64 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 64 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 128 | 4 |
///
/// ## Irreps at $P$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P2 | $P_{1}P2$ | $k_{12}t1t2$ | 8 | F64c | no |
///
/// ### Isotropy subgroups for P1P2 ($P_{1}P2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 16 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 16 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 32 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 32 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 64 | 2 |
///
pub struct Sg142;
