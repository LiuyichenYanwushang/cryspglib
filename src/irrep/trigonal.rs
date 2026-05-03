//! # Trigonal space groups (#143–#167)
//!
//! Trigonal space groups use hexagonal axes (a = b ≠ c, γ = 120°).
//! The Brillouin zone is a hexagonal prism.
//!
//! ## Common k-point labels (hexagonal setting)
//!
//! | Label | Coords (fractional) | Little group |
//! |-------|---------------------|--------------|
//! | Γ | (0, 0, 0) | -3m (D₃ₔ) |
//! | A | (0, 0, ½) | -3m (D₃ₔ) |
//! | K | (⅓, ⅓, 0) | 3m (C₃ᵥ) |
//! | H | (⅓, ⅓, ½) | 3m (C₃ᵥ) |
//! | M | (½, 0, 0) | 2/m (C₂ₕ) |
//! | L | (½, 0, ½) | 2/m (C₂ₕ) |
//!
//! ---

/// # 143 P3 (C₃¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2GM3 | $\Gamma_{2}GM3$ | $k_{16}t2t3$ | 2 | B3a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2GM3 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 3 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C12a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 4 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2A3 | $A_{2}A3$ | $k_{17}t2t3$ | 2 | B6b | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A2A3 ($A_{2}A3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H1 | $H_{1}H1$ | $k_{15}t1t1*$ | 2 | B6b | yes |
/// | H2H2 | $H_{2}H2$ | $k_{15}t3t3*$ | 2 | B6b | yes |
/// | H3H3 | $H_{3}H3$ | $k_{15}t2t2*$ | 2 | B6b | yes |
///
/// ### Isotropy subgroups for H1H1 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 1 |
///
/// ### Isotropy subgroups for H2H2 ($H_{2}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 1 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1K1 | $K_{1}K1$ | $k_{13}t1t1*$ | 2 | B3a | yes |
/// | K2K2 | $K_{2}K2$ | $k_{13}t3t3*$ | 2 | B3a | yes |
/// | K3K3 | $K_{3}K3$ | $k_{13}t2t2*$ | 2 | B3a | yes |
///
/// ### Isotropy subgroups for K1K1 ($K_{1}K1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 3 | 1 |
///
/// ### Isotropy subgroups for K2K2 ($K_{2}K2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 3 | 1 |
///
/// ### Isotropy subgroups for K3K3 ($K_{3}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 3 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
pub struct Sg143;

/// # 144 P3₁ (C₃²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2GM3 | $\Gamma_{2}GM3$ | $k_{16}t2t3$ | 2 | B3a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2GM3 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 3 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C12a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 4 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A3 | $A_{2}A3$ | $k_{17}t1t2$ | 2 | B6b | yes |
/// | A2 | $A_{1}$ | $k_{17}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for A1A3 ($A_{2}A3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H1 | $H_{2}H2$ | $k_{15}t1t1*$ | 2 | B6b | yes |
/// | H2H2 | $H_{1}H1$ | $k_{15}t3t3*$ | 2 | B6b | yes |
/// | H3H3 | $H_{3}H3$ | $k_{15}t2t2*$ | 2 | B6b | yes |
///
/// ### Isotropy subgroups for H1H1 ($H_{2}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 6 | 1 |
///
/// ### Isotropy subgroups for H2H2 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 6 | 1 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 6 | 1 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1K1 | $K_{1}K1$ | $k_{13}t1t1*$ | 2 | B3a | yes |
/// | K2K2 | $K_{2}K2$ | $k_{13}t3t3*$ | 2 | B3a | yes |
/// | K3K3 | $K_{3}K3$ | $k_{13}t2t2*$ | 2 | B3a | yes |
///
/// ### Isotropy subgroups for K1K1 ($K_{1}K1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 3 | 1 |
///
/// ### Isotropy subgroups for K2K2 ($K_{2}K2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 3 | 1 |
///
/// ### Isotropy subgroups for K3K3 ($K_{3}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 3 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 8 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
pub struct Sg144;

/// # 145 P3₂ (C₃³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2GM3 | $\Gamma_{2}GM3$ | $k_{16}t2t3$ | 2 | B3a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2GM3 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 3 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C12a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 4 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A2A3 | $A_{2}A3$ | $k_{17}t1t3$ | 2 | B6b | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A2A3 ($A_{2}A3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H1 | $H_{1}H1$ | $k_{15}t2t2*$ | 2 | B6b | yes |
/// | H2H2 | $H_{2}H2$ | $k_{15}t1t1*$ | 2 | B6b | yes |
/// | H3H3 | $H_{3}H3$ | $k_{15}t3t3*$ | 2 | B6b | yes |
///
/// ### Isotropy subgroups for H1H1 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 6 | 1 |
///
/// ### Isotropy subgroups for H2H2 ($H_{2}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 6 | 1 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 6 | 1 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1K1 | $K_{1}K1$ | $k_{13}t1t1*$ | 2 | B3a | yes |
/// | K2K2 | $K_{2}K2$ | $k_{13}t3t3*$ | 2 | B3a | yes |
/// | K3K3 | $K_{3}K3$ | $k_{13}t2t2*$ | 2 | B3a | yes |
///
/// ### Isotropy subgroups for K1K1 ($K_{1}K1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 3 | 1 |
///
/// ### Isotropy subgroups for K2K2 ($K_{2}K2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 3 | 1 |
///
/// ### Isotropy subgroups for K3K3 ($K_{3}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 3 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 8 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
pub struct Sg145;

/// # 146 R3 (C₃⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM2GM3 | $\Gamma_{2}GM3$ | $k_{7}t2t3$ | 2 | B3a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2GM3 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 3 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{4}t1$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $T$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $Z_{1}$ | $k_{8}t1$ | 1 | A2a | no |
/// | T2T3 | $Z_{2}Z3$ | $k_{8}t2t3$ | 2 | B6b | yes |
///
/// ### Isotropy subgroups for T1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 2 | 1 |
///
/// ### Isotropy subgroups for T2T3 ($Z_{2}Z3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $F$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | F1 | $F_{1}$ | $k_{5}t1$ | 3 | C12a | no |
///
/// ### Isotropy subgroups for F1 ($F_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 4 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 3 |
///
pub struct Sg146;

/// # 147 P-3 (C₃ᵢ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM2+GM3+ | $\Gamma_{2}^+GM3+$ | $k_{16}t3t5$ | 2 | B3a | no |
/// | GM2-GM3- | $\Gamma_{2}^-GM3-$ | $k_{16}t4t6$ | 2 | B6b | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 147 | P-3 | C3i^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+GM3+ ($\Gamma_{2}^+GM3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 3 | 1 |
///
/// ### Isotropy subgroups for GM2-GM3- ($\Gamma_{2}^-GM3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{12}t1$ | 3 | C12a | no |
/// | M1- | $M_{1}^-$ | $k_{12}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 4 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 3 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1+ | $A_{1}^+$ | $k_{17}t1$ | 1 | A2a | no |
/// | A1- | $A_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | A2+A3+ | $A_{2}^+A3+$ | $k_{17}t3t5$ | 2 | B6b | no |
/// | A2-A3- | $A_{2}^-A3-$ | $k_{17}t4t6$ | 2 | B6b | no |
///
/// ### Isotropy subgroups for A1+ ($A_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 147 | P-3 | C3i^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A1- ($A_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 147 | P-3 | C3i^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A2+A3+ ($A_{2}^+A3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
///
/// ### Isotropy subgroups for A2-A3- ($A_{2}^-A3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t1$ | 2 | B12a | yes |
/// | H2H3 | $H_{2}H3$ | $k_{15}t2t3$ | 4 | D36b | no |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 147 | P-3 | C3i^1 | 6 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H2H3 ($H_{2}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 18 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | yes |
/// | K2K3 | $K_{2}K3$ | $k_{13}t2t3$ | 4 | D18a | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 147 | P-3 | C3i^1 | 3 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2K3 ($K_{2}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 9 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 18 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^+$ | $k_{14}t1$ | 3 | C24a | no |
/// | L1- | $L_{1}^-$ | $k_{14}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 8 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 8 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
pub struct Sg147;

/// # 148 R-3 (C₃ᵢ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{7}t2$ | 1 | A2a | no |
/// | GM2+GM3+ | $\Gamma_{2}^+GM3+$ | $k_{7}t3t5$ | 2 | B3a | no |
/// | GM2-GM3- | $\Gamma_{2}^-GM3-$ | $k_{7}t4t6$ | 2 | B6b | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+GM3+ ($\Gamma_{2}^+GM3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 3 | 1 |
///
/// ### Isotropy subgroups for GM2-GM3- ($\Gamma_{2}^-GM3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^+$ | $k_{4}t1$ | 3 | C24a | no |
/// | L1- | $L_{1}^-$ | $k_{4}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 2 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 2 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ## Irreps at $T$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1+ | $Z_{1}^+$ | $k_{8}t1$ | 1 | A2a | no |
/// | T1- | $Z_{1}^-$ | $k_{8}t2$ | 1 | A2a | no |
/// | T2+T3+ | $Z_{2}^+Z3+$ | $k_{8}t3t5$ | 2 | B6b | no |
/// | T2-T3- | $Z_{2}^-Z3-$ | $k_{8}t4t6$ | 2 | B6b | no |
///
/// ### Isotropy subgroups for T1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 2 | 1 |
///
/// ### Isotropy subgroups for T1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 2 | 1 |
///
/// ### Isotropy subgroups for T2+T3+ ($Z_{2}^+Z3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
///
/// ### Isotropy subgroups for T2-T3- ($Z_{2}^-Z3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
///
/// ## Irreps at $F$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | F1+ | $F_{1}^+$ | $k_{5}t1$ | 3 | C12a | no |
/// | F1- | $F_{1}^-$ | $k_{5}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for F1+ ($F_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 4 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 3 |
///
/// ### Isotropy subgroups for F1- ($F_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
pub struct Sg148;

/// # 149 P312 (D₃¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{16}t3$ | 2 | B6a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 3 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
/// | dir{P1} | 149 | P312 | D3^1 | 4 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A3 | $A_{3}$ | $k_{17}t3$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t1$ | 2 | B12a | yes |
/// | H2 | $H_{2}$ | $k_{15}t3$ | 2 | B12a | yes |
/// | H3 | $H_{3}$ | $k_{15}t2$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | yes |
/// | K2 | $K_{2}$ | $k_{13}t3$ | 2 | B6a | yes |
/// | K3 | $K_{3}$ | $k_{13}t2$ | 2 | B6a | yes |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 3 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 3 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 3 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C48a | no |
/// | L2 | $L_{2}$ | $k_{14}t2$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg149;

/// # 150 P321 (D₃²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{16}t3$ | 2 | B6a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 3 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
/// | dir{P1} | 150 | P321 | D3^2 | 4 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A3 | $A_{3}$ | $k_{17}t3$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H1 | $H_{1}H1$ | $k_{15}t1t1*$ | 2 | B6b | no |
/// | H2H2 | $H_{2}H2$ | $k_{15}t2t2*$ | 2 | B6b | no |
/// | H3H3 | $H_{3}H3$ | $k_{15}t3t3*$ | 4 | D36b | no |
///
/// ### Isotropy subgroups for H1H1 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 6 | 1 |
///
/// ### Isotropy subgroups for H2H2 ($H_{2}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 6 | 1 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 1 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1K1 | $K_{1}K1$ | $k_{13}t1t1*$ | 2 | B3a | no |
/// | K2K2 | $K_{2}K2$ | $k_{13}t2t2*$ | 2 | B6b | no |
/// | K3K3 | $K_{3}K3$ | $k_{13}t3t3*$ | 4 | D18a | no |
///
/// ### Isotropy subgroups for K1K1 ($K_{1}K1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 3 | 1 |
///
/// ### Isotropy subgroups for K2K2 ($K_{2}K2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 1 |
///
/// ### Isotropy subgroups for K3K3 ($K_{3}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 9 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 18 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C48a | no |
/// | L2 | $L_{2}$ | $k_{14}t2$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg150;

/// # 151 P3₁12 (D₃³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{16}t3$ | 2 | B6a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 151 | P3_112 | D3^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 3 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
/// | dir{P1} | 151 | P3_112 | D3^3 | 4 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A3 | $A_{3}$ | $k_{17}t3$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 153 | P3_212 | D3^5 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 153 | P3_212 | D3^5 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{2}$ | $k_{15}t1$ | 2 | B12a | yes |
/// | H2 | $H_{1}$ | $k_{15}t3$ | 2 | B12a | yes |
/// | H3 | $H_{3}$ | $k_{15}t2$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for H1 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 6 | 2 |
/// | dir{P1} | 154 | P3_221 | D3^6 | 6 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 12 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 6 | 2 |
/// | dir{P1} | 154 | P3_221 | D3^6 | 6 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 12 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 6 | 2 |
/// | dir{P1} | 154 | P3_221 | D3^6 | 6 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 12 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | yes |
/// | K2 | $K_{2}$ | $k_{13}t3$ | 2 | B6a | yes |
/// | K3 | $K_{3}$ | $k_{13}t2$ | 2 | B6a | yes |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 3 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 3 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 6 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 3 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 6 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t2$ | 3 | C48a | no |
/// | L2 | $L_{2}$ | $k_{14}t1$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 153 | P3_212 | D3^5 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 153 | P3_212 | D3^5 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg151;

/// # 152 P3₁21 (D₃⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{16}t3$ | 2 | B6a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 3 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
/// | dir{P1} | 152 | P3_121 | D3^4 | 4 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A3 | $A_{3}$ | $k_{17}t3$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H1 | $H_{1}H1$ | $k_{15}t1t1*$ | 2 | B6b | no |
/// | H2H2 | $H_{2}H2$ | $k_{15}t2t2*$ | 2 | B6b | no |
/// | H3H3 | $H_{3}H3$ | $k_{15}t3t3*$ | 4 | D36b | no |
///
/// ### Isotropy subgroups for H1H1 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 153 | P3_212 | D3^5 | 6 | 1 |
///
/// ### Isotropy subgroups for H2H2 ($H_{2}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 153 | P3_212 | D3^5 | 6 | 1 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 1 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1K1 | $K_{1}K1$ | $k_{13}t1t1*$ | 2 | B3a | no |
/// | K2K2 | $K_{2}K2$ | $k_{13}t2t2*$ | 2 | B6b | no |
/// | K3K3 | $K_{3}K3$ | $k_{13}t3t3*$ | 4 | D18a | no |
///
/// ### Isotropy subgroups for K1K1 ($K_{1}K1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 151 | P3_112 | D3^3 | 3 | 1 |
///
/// ### Isotropy subgroups for K2K2 ($K_{2}K2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 6 | 1 |
///
/// ### Isotropy subgroups for K3K3 ($K_{3}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 9 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 18 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{2}$ | $k_{14}t2$ | 3 | C48a | no |
/// | L2 | $L_{1}$ | $k_{14}t1$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 154 | P3_221 | D3^6 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 154 | P3_221 | D3^6 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg152;

/// # 153 P3₂12 (D₃⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{16}t3$ | 2 | B6a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 153 | P3_212 | D3^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 3 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
/// | dir{P1} | 153 | P3_212 | D3^5 | 4 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A2 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A3 | $A_{3}$ | $k_{17}t3$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 151 | P3_112 | D3^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 151 | P3_112 | D3^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t2$ | 2 | B12a | yes |
/// | H2 | $H_{2}$ | $k_{15}t1$ | 2 | B12a | yes |
/// | H3 | $H_{3}$ | $k_{15}t3$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 6 | 2 |
/// | dir{P1} | 152 | P3_121 | D3^4 | 6 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 12 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 6 | 2 |
/// | dir{P1} | 152 | P3_121 | D3^4 | 6 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 12 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 6 | 2 |
/// | dir{P1} | 152 | P3_121 | D3^4 | 6 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 12 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | yes |
/// | K2 | $K_{2}$ | $k_{13}t3$ | 2 | B6a | yes |
/// | K3 | $K_{3}$ | $k_{13}t2$ | 2 | B6a | yes |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 3 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 3 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 6 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 3 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 6 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C48a | no |
/// | L2 | $L_{2}$ | $k_{14}t2$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 151 | P3_112 | D3^3 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 151 | P3_112 | D3^3 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg153;

/// # 154 P3₂21 (D₃⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{16}t3$ | 2 | B6a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 3 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
/// | dir{P1} | 154 | P3_221 | D3^6 | 4 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{2}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2 | $A_{1}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A3 | $A_{3}$ | $k_{17}t3$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H1 | $H_{2}H2$ | $k_{15}t1t1*$ | 2 | B6b | no |
/// | H2H2 | $H_{1}H1$ | $k_{15}t2t2*$ | 2 | B6b | no |
/// | H3H3 | $H_{3}H3$ | $k_{15}t3t3*$ | 4 | D36b | no |
///
/// ### Isotropy subgroups for H1H1 ($H_{2}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 151 | P3_112 | D3^3 | 6 | 1 |
///
/// ### Isotropy subgroups for H2H2 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 151 | P3_112 | D3^3 | 6 | 1 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 1 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1K1 | $K_{1}K1$ | $k_{13}t1t1*$ | 2 | B3a | no |
/// | K2K2 | $K_{2}K2$ | $k_{13}t2t2*$ | 2 | B6b | no |
/// | K3K3 | $K_{3}K3$ | $k_{13}t3t3*$ | 4 | D18a | no |
///
/// ### Isotropy subgroups for K1K1 ($K_{1}K1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 153 | P3_212 | D3^5 | 3 | 1 |
///
/// ### Isotropy subgroups for K2K2 ($K_{2}K2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 6 | 1 |
///
/// ### Isotropy subgroups for K3K3 ($K_{3}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 9 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 18 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t2$ | 3 | C48a | no |
/// | L2 | $L_{2}$ | $k_{14}t1$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 152 | P3_121 | D3^4 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 152 | P3_121 | D3^4 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg154;

/// # 155 R32 (D₃⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{7}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{7}t3$ | 2 | B6a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 3 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{4}t1$ | 3 | C48a | no |
/// | L2 | $L_{2}$ | $k_{4}t2$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $T$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $Z_{1}$ | $k_{8}t1$ | 1 | A2a | no |
/// | T2 | $Z_{2}$ | $k_{8}t2$ | 1 | A2a | no |
/// | T3 | $Z_{3}$ | $k_{8}t3$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for T1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 2 | 1 |
///
/// ### Isotropy subgroups for T2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 2 | 1 |
///
/// ### Isotropy subgroups for T3 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $F$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | F1 | $F_{1}$ | $k_{5}t1$ | 3 | C24b | no |
/// | F2 | $F_{2}$ | $k_{5}t2$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for F1 ($F_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
/// | dir{P1} | 155 | R32 | D3^7 | 4 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ### Isotropy subgroups for F2 ($F_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
pub struct Sg155;

/// # 156 P3m1 (C₃ᵥ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{16}t3$ | 2 | B6a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 156 | P3m1 | C3v^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 3 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 6 | 1 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 4 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 6 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A3 | $A_{3}$ | $k_{17}t3$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 156 | P3m1 | C3v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 158 | P3c1 | C3v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t1$ | 2 | B12a | no |
/// | H2 | $H_{2}$ | $k_{15}t3$ | 2 | B12a | no |
/// | H3 | $H_{3}$ | $k_{15}t2$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 6 | 2 |
/// | dir{P1} | 157 | P31m | C3v^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 6 | 2 |
/// | dir{P1} | 157 | P31m | C3v^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 6 | 2 |
/// | dir{P1} | 157 | P31m | C3v^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{2}$ | $k_{13}t3$ | 2 | B6a | no |
/// | K3 | $K_{3}$ | $k_{13}t2$ | 2 | B6a | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 3 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 3 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 3 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C48a | no |
/// | L2 | $L_{2}$ | $k_{14}t2$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 8 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 8 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg156;

/// # 157 P31m (C₃ᵥ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{16}t3$ | 2 | B6a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 3 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 6 | 1 |
/// | dir{P1} | 157 | P31m | C3v^2 | 4 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 6 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A3 | $A_{3}$ | $k_{17}t3$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H1 | $H_{1}H1$ | $k_{15}t1t1*$ | 2 | B6b | yes |
/// | H2H2 | $H_{2}H2$ | $k_{15}t2t2*$ | 2 | B6b | yes |
/// | H3H3 | $H_{3}H3$ | $k_{15}t3t3*$ | 4 | D36b | no |
///
/// ### Isotropy subgroups for H1H1 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 156 | P3m1 | C3v^1 | 6 | 1 |
///
/// ### Isotropy subgroups for H2H2 ($H_{2}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 158 | P3c1 | C3v^3 | 6 | 1 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 18 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 18 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 1 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1K1 | $K_{1}K1$ | $k_{13}t1t1*$ | 2 | B3a | yes |
/// | K2K2 | $K_{2}K2$ | $k_{13}t2t2*$ | 2 | B6b | yes |
/// | K3K3 | $K_{3}K3$ | $k_{13}t3t3*$ | 4 | D18a | no |
///
/// ### Isotropy subgroups for K1K1 ($K_{1}K1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 156 | P3m1 | C3v^1 | 3 | 1 |
///
/// ### Isotropy subgroups for K2K2 ($K_{2}K2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 1 |
///
/// ### Isotropy subgroups for K3K3 ($K_{3}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 9 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 18 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C48a | no |
/// | L2 | $L_{2}$ | $k_{14}t2$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 157 | P31m | C3v^2 | 8 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 159 | P31c | C3v^4 | 8 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg157;

/// # 158 P3c1 (C₃ᵥ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{16}t3$ | 2 | B6a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 158 | P3c1 | C3v^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 3 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 6 | 1 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 4 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 6 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{17}t1t2$ | 2 | B4a | yes |
/// | A3A3 | $A_{3}A3$ | $k_{17}t3t3$ | 4 | D12a | no |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 4 | 1 |
///
/// ### Isotropy subgroups for A3A3 ($A_{3}A3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H1 | $H_{1}H1$ | $k_{15}t1t1$ | 4 | D12a | no |
/// | H2H2 | $H_{2}H2$ | $k_{15}t3t3$ | 4 | D12a | no |
/// | H3H3 | $H_{3}H3$ | $k_{15}t2t2$ | 4 | D12a | no |
///
/// ### Isotropy subgroups for H1H1 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H2H2 ($H_{2}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{2}$ | $k_{13}t3$ | 2 | B6a | no |
/// | K3 | $K_{3}$ | $k_{13}t2$ | 2 | B6a | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 3 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 3 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 3 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1L2 | $L_{1}L2$ | $k_{14}t1t2$ | 6 | E48b | no |
///
/// ### Isotropy subgroups for L1L2 ($L_{1}L2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 16 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg158;

/// # 159 P31c (C₃ᵥ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{16}t3$ | 2 | B6a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 3 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 6 | 1 |
/// | dir{P1} | 159 | P31c | C3v^4 | 4 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 6 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{17}t1t2$ | 2 | B4a | yes |
/// | A3A3 | $A_{3}A3$ | $k_{17}t3t3$ | 4 | D12a | no |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 4 | 1 |
///
/// ### Isotropy subgroups for A3A3 ($A_{3}A3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H1 | $H_{1}H1$ | $k_{15}t1t1*$ | 2 | B12b | yes |
/// | H2H2 | $H_{2}H2$ | $k_{15}t2t2*$ | 2 | B12b | yes |
/// | H3H3 | $H_{3}H3$ | $k_{15}t3t3*$ | 4 | D36a | no |
///
/// ### Isotropy subgroups for H1H1 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 1 |
///
/// ### Isotropy subgroups for H2H2 ($H_{2}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 1 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 1 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1K1 | $K_{1}K1$ | $k_{13}t1t1*$ | 2 | B3a | yes |
/// | K2K2 | $K_{2}K2$ | $k_{13}t2t2*$ | 2 | B6b | yes |
/// | K3K3 | $K_{3}K3$ | $k_{13}t3t3*$ | 4 | D18a | no |
///
/// ### Isotropy subgroups for K1K1 ($K_{1}K1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 158 | P3c1 | C3v^3 | 3 | 1 |
///
/// ### Isotropy subgroups for K2K2 ($K_{2}K2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 1 |
///
/// ### Isotropy subgroups for K3K3 ($K_{3}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 9 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 18 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1L2 | $L_{1}L2$ | $k_{14}t1t2$ | 6 | E48b | no |
///
/// ### Isotropy subgroups for L1L2 ($L_{1}L2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 16 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg159;

/// # 160 R3m (C₃ᵥ⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{7}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{7}t3$ | 2 | B6a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 3 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{4}t1$ | 3 | C48a | no |
/// | L2 | $L_{2}$ | $k_{4}t2$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 160 | R3m | C3v^5 | 8 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 161 | R3c | C3v^6 | 8 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $T$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $Z_{1}$ | $k_{8}t1$ | 1 | A2a | no |
/// | T2 | $Z_{2}$ | $k_{8}t2$ | 1 | A2a | no |
/// | T3 | $Z_{3}$ | $k_{8}t3$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for T1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for T2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for T3 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $F$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | F1 | $F_{1}$ | $k_{5}t1$ | 3 | C24b | no |
/// | F2 | $F_{2}$ | $k_{5}t2$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for F1 ($F_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 6 | 1 |
/// | dir{P1} | 160 | R3m | C3v^5 | 4 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ### Isotropy subgroups for F2 ($F_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 6 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
pub struct Sg160;

/// # 161 R3c (C₃ᵥ⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{7}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{7}t3$ | 2 | B6a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 3 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1L2 | $L_{1}L2$ | $k_{4}t1t2$ | 6 | E48b | no |
///
/// ### Isotropy subgroups for L1L2 ($L_{1}L2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $T$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1T2 | $Z_{1}Z2$ | $k_{8}t1t2$ | 2 | B4a | yes |
/// | T3T3 | $Z_{3}Z3$ | $k_{8}t3t3$ | 4 | D12a | no |
///
/// ### Isotropy subgroups for T1T2 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 4 | 1 |
///
/// ### Isotropy subgroups for T3T3 ($Z_{3}Z3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $F$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | F1 | $F_{1}$ | $k_{5}t2$ | 3 | C24b | no |
/// | F2 | $F_{2}$ | $k_{5}t1$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for F1 ($F_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 6 | 1 |
/// | dir{P1} | 161 | R3c | C3v^6 | 4 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ### Isotropy subgroups for F2 ($F_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 6 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
pub struct Sg161;

/// # 162 P-31m (D₃ₔ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{16}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{16}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 162 | P-31m | D3d^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 147 | P-3 | C3i^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 3 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{12}t1$ | 3 | C24b | no |
/// | M1- | $M_{1}^-$ | $k_{12}t2$ | 3 | C48a | no |
/// | M2+ | $M_{2}^+$ | $k_{12}t3$ | 3 | C24c | no |
/// | M2- | $M_{2}^-$ | $k_{12}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 6 | 1 |
/// | dir{P1} | 162 | P-31m | D3d^1 | 4 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 6 | 1 |
/// | dir{P1} | 149 | P312 | D3^1 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M2+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 6 | 1 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 6 | 1 |
/// | dir{P1} | 157 | P31m | C3v^2 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1+ | $A_{1}^+$ | $k_{17}t1$ | 1 | A2a | no |
/// | A1- | $A_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | A2+ | $A_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | A2- | $A_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | A3+ | $A_{3}^+$ | $k_{17}t5$ | 2 | B12a | no |
/// | A3- | $A_{3}^-$ | $k_{17}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for A1+ ($A_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 162 | P-31m | D3d^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A1- ($A_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 163 | P-31c | D3d^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A2+ ($A_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 163 | P-31c | D3d^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A2- ($A_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 162 | P-31m | D3d^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A3+ ($A_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A3- ($A_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t1$ | 2 | B12a | yes |
/// | H2 | $H_{2}$ | $k_{15}t2$ | 2 | B12a | yes |
/// | H3 | $H_{3}$ | $k_{15}t3$ | 4 | D72d | no |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 6 | 2 |
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 6 | 2 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 6 | 2 |
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 6 | 2 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 12 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 18 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 18 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 18 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 18 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 36 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 36 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | yes |
/// | K2 | $K_{2}$ | $k_{13}t2$ | 2 | B12a | yes |
/// | K3 | $K_{3}$ | $k_{13}t3$ | 4 | D36c | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 3 | 2 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 9 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 18 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^+$ | $k_{14}t1$ | 3 | C48a | no |
/// | L1- | $L_{1}^-$ | $k_{14}t2$ | 3 | C48a | no |
/// | L2+ | $L_{2}^+$ | $k_{14}t3$ | 3 | C48a | no |
/// | L2- | $L_{2}^-$ | $k_{14}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 162 | P-31m | D3d^1 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 163 | P-31c | D3d^2 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2+ ($L_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 163 | P-31c | D3d^2 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2- ($L_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 162 | P-31m | D3d^1 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
pub struct Sg162;

/// # 163 P-31c (D₃ₔ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{16}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{16}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 163 | P-31c | D3d^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 147 | P-3 | C3i^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 3 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{12}t1$ | 3 | C24b | no |
/// | M1- | $M_{1}^-$ | $k_{12}t2$ | 3 | C48a | no |
/// | M2+ | $M_{2}^+$ | $k_{12}t3$ | 3 | C24c | no |
/// | M2- | $M_{2}^-$ | $k_{12}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 6 | 1 |
/// | dir{P1} | 163 | P-31c | D3d^2 | 4 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 6 | 1 |
/// | dir{P1} | 149 | P312 | D3^1 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M2+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 6 | 1 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 6 | 1 |
/// | dir{P1} | 159 | P31c | C3v^4 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{17}t2t3$ | 4 | D24e | no |
/// | A3 | $A_{3}$ | $k_{17}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 4 | 1 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 4 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t1$ | 2 | B24a | yes |
/// | H2 | $H_{2}$ | $k_{15}t2$ | 2 | B24a | yes |
/// | H3 | $H_{3}$ | $k_{15}t3$ | 4 | D72c | no |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 12 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 12 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | yes |
/// | K2 | $K_{2}$ | $k_{13}t2$ | 2 | B12a | yes |
/// | K3 | $K_{3}$ | $k_{13}t3$ | 4 | D36c | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 3 | 2 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 9 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 18 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 6 | E96d | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{C1} | 149 | P312 | D3^1 | 16 | 3 |
/// | dir{P3} | 147 | P-3 | C3i^1 | 16 | 3 |
/// | dir{P1} | 143 | P3 | C3^1 | 32 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
pub struct Sg163;

/// # 164 P-3m1 (D₃ₔ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{16}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{16}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 147 | P-3 | C3i^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 156 | P3m1 | C3v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 3 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{12}t1$ | 3 | C24b | no |
/// | M1- | $M_{1}^-$ | $k_{12}t2$ | 3 | C48a | no |
/// | M2+ | $M_{2}^+$ | $k_{12}t3$ | 3 | C24c | no |
/// | M2- | $M_{2}^-$ | $k_{12}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 6 | 1 |
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 4 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 6 | 1 |
/// | dir{P1} | 150 | P321 | D3^2 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M2+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 6 | 1 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 6 | 1 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1+ | $A_{1}^+$ | $k_{17}t1$ | 1 | A2a | no |
/// | A1- | $A_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | A2+ | $A_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | A2- | $A_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | A3+ | $A_{3}^+$ | $k_{17}t5$ | 2 | B12a | no |
/// | A3- | $A_{3}^-$ | $k_{17}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for A1+ ($A_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A1- ($A_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A2+ ($A_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A2- ($A_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A3+ ($A_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A3- ($A_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t1$ | 2 | B12a | no |
/// | H2 | $H_{2}$ | $k_{15}t2$ | 2 | B12a | no |
/// | H3 | $H_{3}$ | $k_{15}t3$ | 4 | D72d | yes |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 163 | P-31c | D3d^2 | 6 | 2 |
/// | dir{P1} | 162 | P-31m | D3d^1 | 6 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 163 | P-31c | D3d^2 | 6 | 2 |
/// | dir{P1} | 162 | P-31m | D3d^1 | 6 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 12 | 2 |
/// | dir{P1} | 157 | P31m | C3v^2 | 12 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 18 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 18 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 18 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 18 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 36 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 36 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{2}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3 | $K_{3}$ | $k_{13}t3$ | 4 | D36c | yes |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 162 | P-31m | D3d^1 | 3 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 6 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 6 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 9 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 18 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^+$ | $k_{14}t1$ | 3 | C48a | no |
/// | L1- | $L_{1}^-$ | $k_{14}t2$ | 3 | C48a | no |
/// | L2+ | $L_{2}^+$ | $k_{14}t3$ | 3 | C48a | no |
/// | L2- | $L_{2}^-$ | $k_{14}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2+ ($L_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2- ($L_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
pub struct Sg164;

/// # 165 P-3c1 (D₃ₔ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{16}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{16}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 147 | P-3 | C3i^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 158 | P3c1 | C3v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 3 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{12}t1$ | 3 | C24b | no |
/// | M1- | $M_{1}^-$ | $k_{12}t2$ | 3 | C48a | no |
/// | M2+ | $M_{2}^+$ | $k_{12}t3$ | 3 | C24c | no |
/// | M2- | $M_{2}^-$ | $k_{12}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 6 | 1 |
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 4 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 6 | 1 |
/// | dir{P1} | 150 | P321 | D3^2 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M2+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 6 | 1 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 6 | 1 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{17}t2t3$ | 4 | D24e | no |
/// | A3 | $A_{3}$ | $k_{17}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 4 | 1 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 4 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H2 | $H_{1}H2$ | $k_{15}t1t2$ | 4 | D24e | no |
/// | H3H3 | $H_{3}H3$ | $k_{15}t3t3$ | 8 | F72b | no |
///
/// ### Isotropy subgroups for H1H2 ($H_{1}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 12 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 12 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{2}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3 | $K_{3}$ | $k_{13}t3$ | 4 | D36c | yes |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 163 | P-31c | D3d^2 | 3 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 6 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 6 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 9 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 18 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 6 | E96d | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{C1} | 150 | P321 | D3^2 | 16 | 3 |
/// | dir{P3} | 147 | P-3 | C3i^1 | 16 | 3 |
/// | dir{P1} | 143 | P3 | C3^1 | 32 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
pub struct Sg165;

/// # 166 R-3m (D₃ₔ⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{7}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{7}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{7}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{7}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{7}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 3 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^+$ | $k_{4}t1$ | 3 | C48a | no |
/// | L1- | $L_{1}^-$ | $k_{4}t2$ | 3 | C48a | no |
/// | L2+ | $L_{2}^+$ | $k_{4}t3$ | 3 | C48a | no |
/// | L2- | $L_{2}^-$ | $k_{4}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2+ ($L_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2- ($L_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ## Irreps at $T$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1+ | $Z_{1}^+$ | $k_{8}t1$ | 1 | A2a | no |
/// | T1- | $Z_{1}^-$ | $k_{8}t2$ | 1 | A2a | no |
/// | T2+ | $Z_{2}^+$ | $k_{8}t3$ | 1 | A2a | no |
/// | T2- | $Z_{2}^-$ | $k_{8}t4$ | 1 | A2a | no |
/// | T3+ | $Z_{3}^+$ | $k_{8}t5$ | 2 | B12a | no |
/// | T3- | $Z_{3}^-$ | $k_{8}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for T1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 2 | 1 |
///
/// ### Isotropy subgroups for T1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 2 | 1 |
///
/// ### Isotropy subgroups for T2+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 2 | 1 |
///
/// ### Isotropy subgroups for T2- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 2 | 1 |
///
/// ### Isotropy subgroups for T3+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ### Isotropy subgroups for T3- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ## Irreps at $F$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | F1+ | $F_{1}^+$ | $k_{5}t1$ | 3 | C24b | no |
/// | F1- | $F_{1}^-$ | $k_{5}t2$ | 3 | C48a | no |
/// | F2+ | $F_{2}^+$ | $k_{5}t3$ | 3 | C24c | no |
/// | F2- | $F_{2}^-$ | $k_{5}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for F1+ ($F_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 6 | 1 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 4 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for F1- ($F_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 6 | 1 |
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for F2+ ($F_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 6 | 1 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for F2- ($F_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 6 | 1 |
/// | dir{P1} | 160 | R3m | C3v^5 | 8 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg166;

/// # 167 R-3c (D₃ₔ⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{7}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{7}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{7}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{7}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{7}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 3 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{4}t1$ | 6 | E96d | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{C1} | 155 | R32 | D3^7 | 16 | 3 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 16 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $T$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1T2 | $Z_{1}Z2$ | $k_{8}t2t3$ | 4 | D24e | no |
/// | T3 | $Z_{3}$ | $k_{8}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1T2 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for T3 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 4 | 1 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 4 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
///
/// ## Irreps at $F$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | F1+ | $F_{1}^+$ | $k_{5}t2$ | 3 | C24b | no |
/// | F1- | $F_{1}^-$ | $k_{5}t1$ | 3 | C48a | no |
/// | F2+ | $F_{2}^+$ | $k_{5}t4$ | 3 | C24c | no |
/// | F2- | $F_{2}^-$ | $k_{5}t3$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for F1+ ($F_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 6 | 1 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 4 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for F1- ($F_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 6 | 1 |
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for F2+ ($F_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 6 | 1 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for F2- ($F_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 6 | 1 |
/// | dir{P1} | 161 | R3c | C3v^6 | 8 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg167;
