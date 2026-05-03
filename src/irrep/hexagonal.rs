//! # Hexagonal space groups (#168–#194)
//!
//! Hexagonal crystals have a = b ≠ c, α = β = 90°, γ = 120°.
//! The Brillouin zone is a hexagonal prism.
//!
//! ## Common k-point labels
//!
//! | Label | Coords (fractional) | Little group |
//! |-------|---------------------|--------------|
//! | Γ | (0, 0, 0) | 6/mmm (D₆ₕ) |
//! | A | (0, 0, ½) | 6/mmm (D₆ₕ) |
//! | K | (⅓, ⅓, 0) | -6m2 (D₃ₕ) |
//! | H | (⅓, ⅓, ½) | -6m2 (D₃ₕ) |
//! | M | (½, 0, 0) | mmm (D₂ₕ) |
//! | L | (½, 0, ½) | mmm (D₂ₕ) |
//!
//! ---

/// # 168 P6 (C₆¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{4}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM3GM5 | $\Gamma_{5}GM6$ | $k_{16}t3t5$ | 2 | B3a | yes |
/// | GM4GM6 | $\Gamma_{2}GM3$ | $k_{16}t2t6$ | 2 | B6b | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 168 | P6 | C6^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM5 ($\Gamma_{5}GM6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 3 | 1 |
///
/// ### Isotropy subgroups for GM4GM6 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C12a | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
/// | dir{P1} | 168 | P6 | C6^1 | 4 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 12 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2 | $A_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | A3A5 | $A_{5}A6$ | $k_{17}t3t5$ | 2 | B6b | yes |
/// | A4A6 | $A_{2}A3$ | $k_{17}t2t6$ | 2 | B6b | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 168 | P6 | C6^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 173 | P6_3 | C6^6 | 2 | 1 |
///
/// ### Isotropy subgroups for A3A5 ($A_{5}A6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for A4A6 ($A_{2}A3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t1$ | 2 | B12a | no |
/// | H2H3 | $H_{2}H3$ | $k_{15}t2t3$ | 4 | D36b | no |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 173 | P6_3 | C6^6 | 6 | 2 |
/// | dir{P1} | 168 | P6 | C6^1 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H2H3 ($H_{2}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 18 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2K3 | $K_{2}K3$ | $k_{13}t2t3$ | 4 | D18a | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 168 | P6 | C6^1 | 3 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2K3 ($K_{2}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 9 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 18 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C24a | no |
/// | L2 | $L_{2}$ | $k_{14}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 168 | P6 | C6^1 | 8 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 173 | P6_3 | C6^6 | 8 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
///
pub struct Sg168;

/// # 169 P6₁ (C₆²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{4}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM3GM5 | $\Gamma_{5}GM6$ | $k_{16}t3t5$ | 2 | B3a | yes |
/// | GM4GM6 | $\Gamma_{2}GM3$ | $k_{16}t2t6$ | 2 | B6b | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 169 | P6_1 | C6^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM5 ($\Gamma_{5}GM6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 3 | 1 |
///
/// ### Isotropy subgroups for GM4GM6 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C12a | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
/// | dir{P1} | 169 | P6_1 | C6^2 | 4 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 8 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A6 | $A_{3}A6$ | $k_{17}t1t2$ | 2 | B12b | yes |
/// | A2A5 | $A_{4}A5$ | $k_{17}t4t5$ | 2 | B12b | yes |
/// | A3A4 | $A_{1}A2$ | $k_{17}t3t6$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for A1A6 ($A_{3}A6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A2A5 ($A_{4}A5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A3A4 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H3 | $H_{2}H3$ | $k_{15}t1t2$ | 4 | D36a | yes |
/// | H2H2 | $H_{1}H1$ | $k_{15}t3t3$ | 4 | D12a | no |
///
/// ### Isotropy subgroups for H1H3 ($H_{2}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ### Isotropy subgroups for H2H2 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 12 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2K3 | $K_{2}K3$ | $k_{13}t2t3$ | 4 | D18a | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 169 | P6_1 | C6^2 | 3 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 6 | 2 |
///
/// ### Isotropy subgroups for K2K3 ($K_{2}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 6 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 9 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 18 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1L2 | $L_{1}L2$ | $k_{14}t1t2$ | 6 | E48a | yes |
///
/// ### Isotropy subgroups for L1L2 ($L_{1}L2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
/// | dir{P3} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 16 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg169;

/// # 170 P6₅ (C₆³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{4}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM3GM5 | $\Gamma_{5}GM6$ | $k_{16}t3t5$ | 2 | B3a | yes |
/// | GM4GM6 | $\Gamma_{2}GM3$ | $k_{16}t2t6$ | 2 | B6b | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 170 | P6_5 | C6^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM5 ($\Gamma_{5}GM6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 3 | 1 |
///
/// ### Isotropy subgroups for GM4GM6 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C12a | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
/// | dir{P1} | 170 | P6_5 | C6^3 | 4 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 8 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{17}t2t5$ | 2 | B4a | yes |
/// | A3A6 | $A_{4}A5$ | $k_{17}t1t6$ | 2 | B12b | yes |
/// | A4A5 | $A_{3}A6$ | $k_{17}t3t4$ | 2 | B12b | yes |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 4 | 1 |
///
/// ### Isotropy subgroups for A3A6 ($A_{4}A5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A4A5 ($A_{3}A6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H1 | $H_{1}H1$ | $k_{15}t2t2$ | 4 | D12a | no |
/// | H2H3 | $H_{2}H3$ | $k_{15}t1t3$ | 4 | D36a | yes |
///
/// ### Isotropy subgroups for H1H1 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 12 | 2 |
///
/// ### Isotropy subgroups for H2H3 ($H_{2}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2K3 | $K_{2}K3$ | $k_{13}t2t3$ | 4 | D18a | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 170 | P6_5 | C6^3 | 3 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 6 | 2 |
///
/// ### Isotropy subgroups for K2K3 ($K_{2}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 6 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 9 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 18 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1L2 | $L_{1}L2$ | $k_{14}t1t2$ | 6 | E48a | yes |
///
/// ### Isotropy subgroups for L1L2 ($L_{1}L2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
/// | dir{P3} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 16 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg170;

/// # 171 P6₂ (C₆⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{4}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM3GM5 | $\Gamma_{5}GM6$ | $k_{16}t3t5$ | 2 | B3a | yes |
/// | GM4GM6 | $\Gamma_{2}GM3$ | $k_{16}t2t6$ | 2 | B6b | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 171 | P6_2 | C6^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM5 ($\Gamma_{5}GM6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 3 | 1 |
///
/// ### Isotropy subgroups for GM4GM6 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C12a | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
/// | dir{P1} | 171 | P6_2 | C6^4 | 4 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 12 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 8 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{4}$ | $k_{17}t5$ | 1 | A2a | no |
/// | A2 | $A_{1}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A3A5 | $A_{2}A3$ | $k_{17}t1t3$ | 2 | B6b | yes |
/// | A4A6 | $A_{5}A6$ | $k_{17}t4t6$ | 2 | B6b | yes |
///
/// ### Isotropy subgroups for A1 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 172 | P6_4 | C6^5 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 169 | P6_1 | C6^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A3A5 ($A_{2}A3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for A4A6 ($A_{5}A6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t2$ | 2 | B12a | no |
/// | H2H3 | $H_{2}H3$ | $k_{15}t1t3$ | 4 | D36b | no |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 172 | P6_4 | C6^5 | 6 | 2 |
/// | dir{P1} | 169 | P6_1 | C6^2 | 6 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 12 | 2 |
///
/// ### Isotropy subgroups for H2H3 ($H_{2}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 12 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 18 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2K3 | $K_{2}K3$ | $k_{13}t2t3$ | 4 | D18a | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 171 | P6_2 | C6^4 | 3 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 6 | 2 |
///
/// ### Isotropy subgroups for K2K3 ($K_{2}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 6 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 9 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 18 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C24a | no |
/// | L2 | $L_{2}$ | $k_{14}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 172 | P6_4 | C6^5 | 8 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 169 | P6_1 | C6^2 | 8 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
///
pub struct Sg171;

/// # 172 P6₄ (C₆⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{4}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM3GM5 | $\Gamma_{5}GM6$ | $k_{16}t3t5$ | 2 | B3a | yes |
/// | GM4GM6 | $\Gamma_{2}GM3$ | $k_{16}t2t6$ | 2 | B6b | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 172 | P6_4 | C6^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM5 ($\Gamma_{5}GM6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 3 | 1 |
///
/// ### Isotropy subgroups for GM4GM6 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C12a | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
/// | dir{P1} | 172 | P6_4 | C6^5 | 4 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 12 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 8 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A5 | $A_{5}A6$ | $k_{17}t1t5$ | 2 | B6b | yes |
/// | A2A6 | $A_{2}A3$ | $k_{17}t2t4$ | 2 | B6b | yes |
/// | A3 | $A_{1}$ | $k_{17}t3$ | 1 | A2a | no |
/// | A4 | $A_{4}$ | $k_{17}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for A1A5 ($A_{5}A6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for A2A6 ($A_{2}A3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 171 | P6_2 | C6^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 170 | P6_5 | C6^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H3 | $H_{2}H3$ | $k_{15}t1t2$ | 4 | D36b | no |
/// | H2 | $H_{1}$ | $k_{15}t3$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for H1H3 ($H_{2}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 145 | P3_2 | C3^3 | 12 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 18 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 171 | P6_2 | C6^4 | 6 | 2 |
/// | dir{P1} | 170 | P6_5 | C6^3 | 6 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 12 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2K3 | $K_{2}K3$ | $k_{13}t2t3$ | 4 | D18a | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 172 | P6_4 | C6^5 | 3 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 6 | 2 |
///
/// ### Isotropy subgroups for K2K3 ($K_{2}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 144 | P3_1 | C3^2 | 6 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 9 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 18 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C24a | no |
/// | L2 | $L_{2}$ | $k_{14}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 171 | P6_2 | C6^4 | 8 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 2 |
/// | dir{P1} | 170 | P6_5 | C6^3 | 8 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
///
pub struct Sg172;

/// # 173 P6₃ (C₆⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{4}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM3GM5 | $\Gamma_{5}GM6$ | $k_{16}t3t5$ | 2 | B3a | yes |
/// | GM4GM6 | $\Gamma_{2}GM3$ | $k_{16}t2t6$ | 2 | B6b | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 173 | P6_3 | C6^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM5 ($\Gamma_{5}GM6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 3 | 1 |
///
/// ### Isotropy subgroups for GM4GM6 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C12a | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
/// | dir{P1} | 173 | P6_3 | C6^6 | 4 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{17}t1t4$ | 2 | B4a | yes |
/// | A3A6 | $A_{3}A6$ | $k_{17}t2t3$ | 2 | B12b | yes |
/// | A4A5 | $A_{4}A5$ | $k_{17}t5t6$ | 2 | B12b | yes |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 4 | 1 |
///
/// ### Isotropy subgroups for A3A6 ($A_{3}A6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A4A5 ($A_{4}A5$)
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
/// | H2H3 | $H_{2}H3$ | $k_{15}t2t3$ | 4 | D36a | yes |
///
/// ### Isotropy subgroups for H1H1 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H2H3 ($H_{2}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2K3 | $K_{2}K3$ | $k_{13}t2t3$ | 4 | D18a | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 173 | P6_3 | C6^6 | 3 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2K3 ($K_{2}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 9 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 18 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1L2 | $L_{1}L2$ | $k_{14}t1t2$ | 6 | E48a | yes |
///
/// ### Isotropy subgroups for L1L2 ($L_{1}L2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
/// | dir{P3} | 1 | P1 | C1^1 | 24 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 16 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
pub struct Sg173;

/// # 174 P-6 (C₃ₕ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{4}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3GM5 | $\Gamma_{2}GM3$ | $k_{16}t3t5$ | 2 | B3a | no |
/// | GM4GM6 | $\Gamma_{5}GM6$ | $k_{16}t4t6$ | 2 | B6b | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3GM5 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 3 | 1 |
///
/// ### Isotropy subgroups for GM4GM6 ($\Gamma_{5}GM6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C12a | no |
/// | M2 | $M_{2}$ | $k_{12}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 6 | 1 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 4 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 12 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 6 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 12 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2 | $A_{4}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A3A5 | $A_{2}A3$ | $k_{17}t3t5$ | 2 | B6b | no |
/// | A4A6 | $A_{5}A6$ | $k_{17}t4t6$ | 2 | B6b | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A3A5 ($A_{2}A3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 6 | 1 |
///
/// ### Isotropy subgroups for A4A6 ($A_{5}A6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H1 | $H_{1}H1$ | $k_{15}t1t1*$ | 2 | B6b | no |
/// | H2H2 | $H_{4}H4$ | $k_{15}t2t2*$ | 2 | B6b | no |
/// | H3H3 | $H_{3}H3$ | $k_{15}t5t5*$ | 2 | B6b | no |
/// | H4H4 | $H_{6}H6$ | $k_{15}t6t6*$ | 2 | B6b | no |
/// | H5H5 | $H_{2}H2$ | $k_{15}t3t3*$ | 2 | B6b | no |
/// | H6H6 | $H_{5}H5$ | $k_{15}t4t4*$ | 2 | B6b | no |
///
/// ### Isotropy subgroups for H1H1 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for H2H2 ($H_{4}H4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for H4H4 ($H_{6}H6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for H5H5 ($H_{2}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for H6H6 ($H_{5}H5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1K1 | $K_{1}K1$ | $k_{13}t1t1*$ | 2 | B3a | no |
/// | K2K2 | $K_{4}K4$ | $k_{13}t2t2*$ | 2 | B6b | no |
/// | K3K3 | $K_{3}K3$ | $k_{13}t5t5*$ | 2 | B3a | no |
/// | K4K4 | $K_{6}K6$ | $k_{13}t6t6*$ | 2 | B6b | no |
/// | K5K5 | $K_{2}K2$ | $k_{13}t3t3*$ | 2 | B3a | no |
/// | K6K6 | $K_{5}K5$ | $k_{13}t4t4*$ | 2 | B6b | no |
///
/// ### Isotropy subgroups for K1K1 ($K_{1}K1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 3 | 1 |
///
/// ### Isotropy subgroups for K2K2 ($K_{4}K4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 1 |
///
/// ### Isotropy subgroups for K3K3 ($K_{3}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 3 | 1 |
///
/// ### Isotropy subgroups for K4K4 ($K_{6}K6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 1 |
///
/// ### Isotropy subgroups for K5K5 ($K_{2}K2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 3 | 1 |
///
/// ### Isotropy subgroups for K6K6 ($K_{5}K5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C24a | no |
/// | L2 | $L_{2}$ | $k_{14}t2$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 8 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 8 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 3 |
///
pub struct Sg174;

/// # 175 P6/m (C₆ₕ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{4}^+$ | $k_{16}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{4}^-$ | $k_{16}t8$ | 1 | A2a | no |
/// | GM3+GM5+ | $\Gamma_{5}^+GM6+$ | $k_{16}t5t9$ | 2 | B3a | no |
/// | GM3-GM5- | $\Gamma_{5}^-GM6-$ | $k_{16}t6t10$ | 2 | B6b | no |
/// | GM4+GM6+ | $\Gamma_{2}^+GM3+$ | $k_{16}t3t11$ | 2 | B6b | no |
/// | GM4-GM6- | $\Gamma_{2}^-GM3-$ | $k_{16}t4t12$ | 2 | B6b | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 175 | P6/m | C6h^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 168 | P6 | C6^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 147 | P-3 | C3i^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+GM5+ ($\Gamma_{5}^+GM6+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 3 | 1 |
///
/// ### Isotropy subgroups for GM3-GM5- ($\Gamma_{5}^-GM6-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4+GM6+ ($\Gamma_{2}^+GM3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4-GM6- ($\Gamma_{2}^-GM3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{12}t1$ | 3 | C12a | no |
/// | M1- | $M_{1}^-$ | $k_{12}t2$ | 3 | C24a | no |
/// | M2+ | $M_{2}^+$ | $k_{12}t3$ | 3 | C24a | no |
/// | M2- | $M_{2}^-$ | $k_{12}t4$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 6 | 1 |
/// | dir{P1} | 175 | P6/m | C6h^1 | 4 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 12 | 3 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 6 | 1 |
/// | dir{P1} | 168 | P6 | C6^1 | 8 | 3 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 12 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 13 | P2/c | C2h^4 | 6 | 1 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 8 | 3 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 12 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 6 | 1 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 8 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 12 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1+ | $A_{1}^+$ | $k_{17}t1$ | 1 | A2a | no |
/// | A1- | $A_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | A2+ | $A_{4}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | A2- | $A_{4}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | A3+A5+ | $A_{5}^+A6+$ | $k_{17}t5t9$ | 2 | B6b | no |
/// | A3-A5- | $A_{5}^-A6-$ | $k_{17}t6t10$ | 2 | B6b | no |
/// | A4+A6+ | $A_{2}^+A3+$ | $k_{17}t3t11$ | 2 | B6b | no |
/// | A4-A6- | $A_{2}^-A3-$ | $k_{17}t4t12$ | 2 | B6b | no |
///
/// ### Isotropy subgroups for A1+ ($A_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 175 | P6/m | C6h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A1- ($A_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 175 | P6/m | C6h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A2+ ($A_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A2- ($A_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A3+A5+ ($A_{5}^+A6+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for A3-A5- ($A_{5}^-A6-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 10 | P2/m | C2h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for A4+A6+ ($A_{2}^+A3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 6 | 1 |
///
/// ### Isotropy subgroups for A4-A6- ($A_{2}^-A3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 6 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t1$ | 2 | B12a | no |
/// | H2 | $H_{4}$ | $k_{15}t2$ | 2 | B12a | no |
/// | H3H5 | $H_{2}H3$ | $k_{15}t3t5$ | 4 | D36b | no |
/// | H4H6 | $H_{5}H6$ | $k_{15}t4t6$ | 4 | D36b | no |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 6 | 2 |
/// | dir{P1} | 175 | P6/m | C6h^1 | 6 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 6 | 2 |
/// | dir{P1} | 175 | P6/m | C6h^1 | 6 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H3H5 ($H_{2}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 18 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 18 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 36 | 2 |
///
/// ### Isotropy subgroups for H4H6 ($H_{5}H6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 18 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 18 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{4}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3K5 | $K_{2}K3$ | $k_{13}t3t5$ | 4 | D18a | no |
/// | K4K6 | $K_{5}K6$ | $k_{13}t4t6$ | 4 | D36b | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 175 | P6/m | C6h^1 | 3 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 168 | P6 | C6^1 | 6 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3K5 ($K_{2}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 9 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 18 | 2 |
///
/// ### Isotropy subgroups for K4K6 ($K_{5}K6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 18 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^+$ | $k_{14}t1$ | 3 | C24a | no |
/// | L1- | $L_{1}^-$ | $k_{14}t2$ | 3 | C24a | no |
/// | L2+ | $L_{2}^+$ | $k_{14}t3$ | 3 | C24a | no |
/// | L2- | $L_{2}^-$ | $k_{14}t4$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 175 | P6/m | C6h^1 | 8 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 24 | 3 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 175 | P6/m | C6h^1 | 8 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 24 | 3 |
///
/// ### Isotropy subgroups for L2+ ($L_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 8 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 3 |
///
/// ### Isotropy subgroups for L2- ($L_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 2 |
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 8 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 3 |
///
pub struct Sg175;

/// # 176 P6₃/m (C₆ₕ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{4}^+$ | $k_{16}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{4}^-$ | $k_{16}t8$ | 1 | A2a | no |
/// | GM3+GM5+ | $\Gamma_{5}^+GM6+$ | $k_{16}t5t9$ | 2 | B3a | no |
/// | GM3-GM5- | $\Gamma_{5}^-GM6-$ | $k_{16}t6t10$ | 2 | B6b | no |
/// | GM4+GM6+ | $\Gamma_{2}^+GM3+$ | $k_{16}t3t11$ | 2 | B6b | no |
/// | GM4-GM6- | $\Gamma_{2}^-GM3-$ | $k_{16}t4t12$ | 2 | B6b | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 173 | P6_3 | C6^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 147 | P-3 | C3i^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+GM5+ ($\Gamma_{5}^+GM6+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 3 | 1 |
///
/// ### Isotropy subgroups for GM3-GM5- ($\Gamma_{5}^-GM6-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4+GM6+ ($\Gamma_{2}^+GM3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4-GM6- ($\Gamma_{2}^-GM3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 6 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{12}t1$ | 3 | C12a | no |
/// | M1- | $M_{1}^-$ | $k_{12}t2$ | 3 | C24a | no |
/// | M2+ | $M_{2}^+$ | $k_{12}t3$ | 3 | C24a | no |
/// | M2- | $M_{2}^-$ | $k_{12}t4$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 6 | 1 |
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 4 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 12 | 3 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 6 | 1 |
/// | dir{P1} | 173 | P6_3 | C6^6 | 8 | 3 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 12 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M2+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 6 | 1 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 8 | 3 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 12 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 6 | 1 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 8 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 12 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{2}$ | $k_{17}t1$ | 2 | B8a | yes |
/// | A2A3 | $A_{1}A3$ | $k_{17}t2t3$ | 4 | D24d | yes |
///
/// ### Isotropy subgroups for A1 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 4 | 1 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 4 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 1 |
///
/// ### Isotropy subgroups for A2A3 ($A_{1}A3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H2 | $H_{1}H4$ | $k_{15}t1t2$ | 4 | D24e | no |
/// | H3H6 | $H_{2}H6$ | $k_{15}t4t5$ | 4 | D72b | no |
/// | H4H5 | $H_{3}H5$ | $k_{15}t3t6$ | 4 | D72b | no |
///
/// ### Isotropy subgroups for H1H2 ($H_{1}H4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 12 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
///
/// ### Isotropy subgroups for H3H6 ($H_{2}H6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ### Isotropy subgroups for H4H5 ($H_{3}H5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{4}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3K5 | $K_{2}K3$ | $k_{13}t3t5$ | 4 | D18a | no |
/// | K4K6 | $K_{5}K6$ | $k_{13}t4t6$ | 4 | D36b | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 3 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 173 | P6_3 | C6^6 | 6 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3K5 ($K_{2}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 9 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 18 | 2 |
///
/// ### Isotropy subgroups for K4K6 ($K_{5}K6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 18 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 6 | E96a | yes |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{C1} | 174 | P-6 | C3h^1 | 16 | 3 |
/// | dir{P3} | 147 | P-3 | C3i^1 | 16 | 3 |
/// | dir{P1} | 143 | P3 | C3^1 | 32 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 48 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
pub struct Sg176;

/// # 177 P622 (D₆¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{16}t6$ | 2 | B6a | yes |
/// | GM6 | $\Gamma_{6}$ | $k_{16}t5$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 177 | P622 | D6^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 168 | P6 | C6^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 3 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM6 ($\Gamma_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{12}t2$ | 3 | C24c | no |
/// | M3 | $M_{2}$ | $k_{12}t3$ | 3 | C48a | no |
/// | M4 | $M_{4}$ | $k_{12}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 6 | 1 |
/// | dir{P1} | 177 | P622 | D6^1 | 4 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 6 | 1 |
/// | dir{P1} | 168 | P6 | C6^1 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 6 | 1 |
/// | dir{P1} | 150 | P321 | D3^2 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 6 | 1 |
/// | dir{P1} | 149 | P312 | D3^1 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A3 | $A_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | A4 | $A_{3}$ | $k_{17}t4$ | 1 | A2a | no |
/// | A5 | $A_{5}$ | $k_{17}t6$ | 2 | B12a | yes |
/// | A6 | $A_{6}$ | $k_{17}t5$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 177 | P622 | D6^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 177 | P622 | D6^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 182 | P6_322 | D6^6 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 182 | P6_322 | D6^6 | 2 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A6 ($A_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 6 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 6 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t1$ | 2 | B12a | no |
/// | H2 | $H_{2}$ | $k_{15}t2$ | 2 | B12a | no |
/// | H3 | $H_{3}$ | $k_{15}t3$ | 4 | D72d | no |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 182 | P6_322 | D6^6 | 6 | 2 |
/// | dir{P1} | 177 | P622 | D6^1 | 6 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 182 | P6_322 | D6^6 | 6 | 2 |
/// | dir{P1} | 177 | P622 | D6^1 | 6 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 18 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 18 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 18 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 18 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 36 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{2}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3 | $K_{3}$ | $k_{13}t3$ | 4 | D36c | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 177 | P622 | D6^1 | 3 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 168 | P6 | C6^1 | 6 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 9 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C48a | no |
/// | L2 | $L_{3}$ | $k_{14}t2$ | 3 | C48a | no |
/// | L3 | $L_{2}$ | $k_{14}t3$ | 3 | C48a | no |
/// | L4 | $L_{4}$ | $k_{14}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 177 | P622 | D6^1 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 177 | P622 | D6^1 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L3 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 182 | P6_322 | D6^6 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 3 |
///
/// ### Isotropy subgroups for L4 ($L_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 182 | P6_322 | D6^6 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 3 |
///
pub struct Sg177;

/// # 178 P6₁22 (D₆²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{16}t6$ | 2 | B6a | yes |
/// | GM6 | $\Gamma_{6}$ | $k_{16}t5$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 178 | P6_122 | D6^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 169 | P6_1 | C6^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 151 | P3_112 | D3^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 3 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
///
/// ### Isotropy subgroups for GM6 ($\Gamma_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{12}t2$ | 3 | C24c | no |
/// | M3 | $M_{2}$ | $k_{12}t3$ | 3 | C48a | no |
/// | M4 | $M_{4}$ | $k_{12}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 6 | 1 |
/// | dir{P1} | 178 | P6_122 | D6^2 | 4 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 6 | 1 |
/// | dir{P1} | 169 | P6_1 | C6^2 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 6 | 1 |
/// | dir{P1} | 152 | P3_121 | D3^4 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 6 | 1 |
/// | dir{P1} | 151 | P3_112 | D3^3 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 2 | B24a | yes |
/// | A2 | $A_{2}$ | $k_{17}t3$ | 2 | B24a | yes |
/// | A3 | $A_{3}$ | $k_{17}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 4 | 1 |
/// | dir{P1} | 153 | P3_212 | D3^5 | 4 | 1 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H2 | $H_{1}H2$ | $k_{15}t1t2$ | 4 | D24e | no |
/// | H3 | $H_{3}$ | $k_{15}t3$ | 4 | D72c | yes |
///
/// ### Isotropy subgroups for H1H2 ($H_{1}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 12 | 2 |
/// | dir{P1} | 153 | P3_212 | D3^5 | 12 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 24 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 12 | 2 |
/// | dir{P1} | 154 | P3_221 | D3^6 | 12 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{2}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3 | $K_{3}$ | $k_{13}t3$ | 4 | D36c | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 178 | P6_122 | D6^2 | 3 | 2 |
/// | dir{P1} | 151 | P3_112 | D3^3 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 169 | P6_1 | C6^2 | 6 | 2 |
/// | dir{P1} | 152 | P3_121 | D3^4 | 6 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 6 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 9 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 6 | E96b | yes |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
/// | dir{C1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{C1} | 154 | P3_221 | D3^6 | 16 | 3 |
/// | dir{P3} | 153 | P3_212 | D3^5 | 16 | 3 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 32 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
pub struct Sg178;

/// # 179 P6₅22 (D₆³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{16}t6$ | 2 | B6a | yes |
/// | GM6 | $\Gamma_{6}$ | $k_{16}t5$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 179 | P6_522 | D6^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 170 | P6_5 | C6^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 153 | P3_212 | D3^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 3 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
///
/// ### Isotropy subgroups for GM6 ($\Gamma_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{12}t2$ | 3 | C24c | no |
/// | M3 | $M_{2}$ | $k_{12}t3$ | 3 | C48a | no |
/// | M4 | $M_{4}$ | $k_{12}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 6 | 1 |
/// | dir{P1} | 179 | P6_522 | D6^3 | 4 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 6 | 1 |
/// | dir{P1} | 170 | P6_5 | C6^3 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 6 | 1 |
/// | dir{P1} | 154 | P3_221 | D3^6 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 6 | 1 |
/// | dir{P1} | 153 | P3_212 | D3^5 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{3}$ | $k_{17}t2$ | 2 | B8a | yes |
/// | A2 | $A_{2}$ | $k_{17}t1$ | 2 | B24a | yes |
/// | A3 | $A_{1}$ | $k_{17}t3$ | 2 | B24a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 4 | 1 |
/// | dir{P1} | 151 | P3_112 | D3^3 | 4 | 1 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H2 | $H_{1}H2$ | $k_{15}t1t2$ | 4 | D24e | no |
/// | H3 | $H_{3}$ | $k_{15}t3$ | 4 | D72c | yes |
///
/// ### Isotropy subgroups for H1H2 ($H_{1}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 12 | 2 |
/// | dir{P1} | 151 | P3_112 | D3^3 | 12 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 24 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 12 | 2 |
/// | dir{P1} | 152 | P3_121 | D3^4 | 12 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{2}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3 | $K_{3}$ | $k_{13}t3$ | 4 | D36c | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 179 | P6_522 | D6^3 | 3 | 2 |
/// | dir{P1} | 153 | P3_212 | D3^5 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 170 | P6_5 | C6^3 | 6 | 2 |
/// | dir{P1} | 154 | P3_221 | D3^6 | 6 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 6 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 9 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 6 | E96b | yes |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
/// | dir{C1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{C1} | 152 | P3_121 | D3^4 | 16 | 3 |
/// | dir{P3} | 151 | P3_112 | D3^3 | 16 | 3 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 32 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
pub struct Sg179;

/// # 180 P6₂22 (D₆⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{16}t6$ | 2 | B6a | yes |
/// | GM6 | $\Gamma_{6}$ | $k_{16}t5$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 180 | P6_222 | D6^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 171 | P6_2 | C6^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 153 | P3_212 | D3^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 3 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM6 ($\Gamma_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{12}t2$ | 3 | C24c | no |
/// | M3 | $M_{2}$ | $k_{12}t3$ | 3 | C48a | no |
/// | M4 | $M_{4}$ | $k_{12}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 6 | 1 |
/// | dir{P1} | 180 | P6_222 | D6^4 | 4 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 6 | 1 |
/// | dir{P1} | 171 | P6_2 | C6^4 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 6 | 1 |
/// | dir{P1} | 154 | P3_221 | D3^6 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 6 | 1 |
/// | dir{P1} | 153 | P3_212 | D3^5 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | A2 | $A_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | A3 | $A_{2}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A4 | $A_{1}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A5 | $A_{6}$ | $k_{17}t5$ | 2 | B12a | yes |
/// | A6 | $A_{5}$ | $k_{17}t6$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 181 | P6_422 | D6^5 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 181 | P6_422 | D6^5 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 178 | P6_122 | D6^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 178 | P6_122 | D6^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A6 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 6 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 6 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t1$ | 2 | B12a | no |
/// | H2 | $H_{2}$ | $k_{15}t2$ | 2 | B12a | no |
/// | H3 | $H_{3}$ | $k_{15}t3$ | 4 | D72d | no |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 181 | P6_422 | D6^5 | 6 | 2 |
/// | dir{P1} | 178 | P6_122 | D6^2 | 6 | 2 |
/// | dir{P1} | 151 | P3_112 | D3^3 | 12 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 181 | P6_422 | D6^5 | 6 | 2 |
/// | dir{P1} | 178 | P6_122 | D6^2 | 6 | 2 |
/// | dir{P1} | 151 | P3_112 | D3^3 | 12 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 12 | 2 |
/// | dir{P1} | 152 | P3_121 | D3^4 | 12 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 18 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 18 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 18 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 18 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 36 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{2}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3 | $K_{3}$ | $k_{13}t3$ | 4 | D36c | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 180 | P6_222 | D6^4 | 3 | 2 |
/// | dir{P1} | 153 | P3_212 | D3^5 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 171 | P6_2 | C6^4 | 6 | 2 |
/// | dir{P1} | 154 | P3_221 | D3^6 | 6 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 6 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 9 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{3}$ | $k_{14}t2$ | 3 | C48a | no |
/// | L2 | $L_{1}$ | $k_{14}t1$ | 3 | C48a | no |
/// | L3 | $L_{4}$ | $k_{14}t4$ | 3 | C48a | no |
/// | L4 | $L_{2}$ | $k_{14}t3$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 181 | P6_422 | D6^5 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 181 | P6_422 | D6^5 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L3 ($L_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 178 | P6_122 | D6^2 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 3 |
///
/// ### Isotropy subgroups for L4 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 178 | P6_122 | D6^2 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 3 |
///
pub struct Sg180;

/// # 181 P6₄22 (D₆⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{16}t6$ | 2 | B6a | yes |
/// | GM6 | $\Gamma_{6}$ | $k_{16}t5$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 181 | P6_422 | D6^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 172 | P6_4 | C6^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 151 | P3_112 | D3^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 3 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM6 ($\Gamma_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{12}t2$ | 3 | C24c | no |
/// | M3 | $M_{2}$ | $k_{12}t3$ | 3 | C48a | no |
/// | M4 | $M_{4}$ | $k_{12}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 6 | 1 |
/// | dir{P1} | 181 | P6_422 | D6^5 | 4 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 6 | 1 |
/// | dir{P1} | 172 | P6_4 | C6^5 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 6 | 1 |
/// | dir{P1} | 152 | P3_121 | D3^4 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 6 | 1 |
/// | dir{P1} | 151 | P3_112 | D3^3 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A3 | $A_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | A4 | $A_{3}$ | $k_{17}t4$ | 1 | A2a | no |
/// | A5 | $A_{5}$ | $k_{17}t6$ | 2 | B12a | yes |
/// | A6 | $A_{6}$ | $k_{17}t5$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 180 | P6_222 | D6^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 180 | P6_222 | D6^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 179 | P6_522 | D6^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 179 | P6_522 | D6^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A6 ($A_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 6 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 6 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t1$ | 2 | B12a | no |
/// | H2 | $H_{2}$ | $k_{15}t2$ | 2 | B12a | no |
/// | H3 | $H_{3}$ | $k_{15}t3$ | 4 | D72d | no |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 180 | P6_222 | D6^4 | 6 | 2 |
/// | dir{P1} | 179 | P6_522 | D6^3 | 6 | 2 |
/// | dir{P1} | 153 | P3_212 | D3^5 | 12 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 180 | P6_222 | D6^4 | 6 | 2 |
/// | dir{P1} | 179 | P6_522 | D6^3 | 6 | 2 |
/// | dir{P1} | 153 | P3_212 | D3^5 | 12 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 154 | P3_221 | D3^6 | 12 | 2 |
/// | dir{P1} | 154 | P3_221 | D3^6 | 12 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 18 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 18 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 18 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 18 | 2 |
/// | dir{P1} | 145 | P3_2 | C3^3 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 36 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{2}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3 | $K_{3}$ | $k_{13}t3$ | 4 | D36c | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 181 | P6_422 | D6^5 | 3 | 2 |
/// | dir{P1} | 151 | P3_112 | D3^3 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 172 | P6_4 | C6^5 | 6 | 2 |
/// | dir{P1} | 152 | P3_121 | D3^4 | 6 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 152 | P3_121 | D3^4 | 6 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 9 | 2 |
/// | dir{P1} | 144 | P3_1 | C3^2 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C48a | no |
/// | L2 | $L_{3}$ | $k_{14}t2$ | 3 | C48a | no |
/// | L3 | $L_{2}$ | $k_{14}t3$ | 3 | C48a | no |
/// | L4 | $L_{4}$ | $k_{14}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 180 | P6_222 | D6^4 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 180 | P6_222 | D6^4 | 8 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L3 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 179 | P6_522 | D6^3 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 3 |
///
/// ### Isotropy subgroups for L4 ($L_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 179 | P6_522 | D6^3 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 3 |
///
pub struct Sg181;

/// # 182 P6₃22 (D₆⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{16}t6$ | 2 | B6a | yes |
/// | GM6 | $\Gamma_{6}$ | $k_{16}t5$ | 2 | B12a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 182 | P6_322 | D6^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 173 | P6_3 | C6^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 3 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
///
/// ### Isotropy subgroups for GM6 ($\Gamma_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{12}t2$ | 3 | C24c | no |
/// | M3 | $M_{2}$ | $k_{12}t3$ | 3 | C48a | no |
/// | M4 | $M_{4}$ | $k_{12}t4$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 6 | 1 |
/// | dir{P1} | 182 | P6_322 | D6^6 | 4 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 6 | 1 |
/// | dir{P1} | 173 | P6_3 | C6^6 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 6 | 1 |
/// | dir{P1} | 150 | P321 | D3^2 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 6 | 1 |
/// | dir{P1} | 149 | P312 | D3^1 | 8 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{3}$ | $k_{17}t1$ | 2 | B8a | yes |
/// | A2 | $A_{1}$ | $k_{17}t2$ | 2 | B24a | yes |
/// | A3 | $A_{2}$ | $k_{17}t3$ | 2 | B24a | yes |
///
/// ### Isotropy subgroups for A1 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 4 | 1 |
/// | dir{P1} | 149 | P312 | D3^1 | 4 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H2 | $H_{1}H2$ | $k_{15}t1t2$ | 4 | D24e | no |
/// | H3 | $H_{3}$ | $k_{15}t3$ | 4 | D72c | yes |
///
/// ### Isotropy subgroups for H1H2 ($H_{1}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 12 | 2 |
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
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{2}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3 | $K_{3}$ | $k_{13}t3$ | 4 | D36c | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 182 | P6_322 | D6^6 | 3 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 173 | P6_3 | C6^6 | 6 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 9 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 6 | E96b | yes |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
/// | dir{C1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{C1} | 150 | P321 | D3^2 | 16 | 3 |
/// | dir{P3} | 149 | P312 | D3^1 | 16 | 3 |
/// | dir{P1} | 143 | P3 | C3^1 | 32 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
pub struct Sg182;

/// # 183 P6mm (C₆ᵥ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{16}t6$ | 2 | B6a | no |
/// | GM6 | $\Gamma_{6}$ | $k_{16}t5$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 183 | P6mm | C6v^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 168 | P6 | C6^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 156 | P3m1 | C3v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 3 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM6 ($\Gamma_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{12}t2$ | 3 | C24c | no |
/// | M3 | $M_{4}$ | $k_{12}t4$ | 3 | C48a | no |
/// | M4 | $M_{2}$ | $k_{12}t3$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 6 | 1 |
/// | dir{P1} | 183 | P6mm | C6v^1 | 4 | 3 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 12 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 32 | Pba2 | C2v^8 | 6 | 1 |
/// | dir{P1} | 168 | P6 | C6^1 | 8 | 3 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 12 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 6 | 1 |
/// | dir{P1} | 157 | P31m | C3v^2 | 8 | 3 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 12 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 6 | 1 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 8 | 3 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 12 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A3 | $A_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | A4 | $A_{3}$ | $k_{17}t3$ | 1 | A2a | no |
/// | A5 | $A_{5}$ | $k_{17}t6$ | 2 | B12a | no |
/// | A6 | $A_{6}$ | $k_{17}t5$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 183 | P6mm | C6v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 184 | P6cc | C6v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 185 | P6_3cm | C6v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 186 | P6_3mc | C6v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 6 | 1 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 6 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A6 ($A_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 6 | 1 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 6 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 1 |
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
/// | dir{P1} | 186 | P6_3mc | C6v^4 | 6 | 2 |
/// | dir{P1} | 183 | P6mm | C6v^1 | 6 | 2 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 185 | P6_3cm | C6v^3 | 6 | 2 |
/// | dir{P1} | 184 | P6cc | C6v^2 | 6 | 2 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 12 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 12 | 2 |
/// | dir{P1} | 157 | P31m | C3v^2 | 12 | 2 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 18 | 2 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 18 | 2 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 18 | 2 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 18 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 36 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 36 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 36 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 36 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 36 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 36 | 2 |
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
/// | dir{P1} | 183 | P6mm | C6v^1 | 3 | 2 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 168 | P6 | C6^1 | 6 | 2 |
/// | dir{P1} | 157 | P31m | C3v^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 6 | 2 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 9 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 18 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 18 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C48a | no |
/// | L2 | $L_{3}$ | $k_{14}t2$ | 3 | C48a | no |
/// | L3 | $L_{4}$ | $k_{14}t4$ | 3 | C48a | no |
/// | L4 | $L_{2}$ | $k_{14}t3$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 44 | Imm2 | C2v^20 | 6 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 183 | P6mm | C6v^1 | 8 | 3 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 45 | Iba2 | C2v^21 | 6 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 184 | P6cc | C6v^2 | 8 | 3 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L3 ($L_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 46 | Ima2 | C2v^22 | 6 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 185 | P6_3cm | C6v^3 | 8 | 3 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 3 |
///
/// ### Isotropy subgroups for L4 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 46 | Ima2 | C2v^22 | 6 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 186 | P6_3mc | C6v^4 | 8 | 3 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 3 |
///
pub struct Sg183;

/// # 184 P6cc (C₆ᵥ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{16}t6$ | 2 | B6a | no |
/// | GM6 | $\Gamma_{6}$ | $k_{16}t5$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 184 | P6cc | C6v^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 168 | P6 | C6^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 158 | P3c1 | C3v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 3 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM6 ($\Gamma_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{12}t2$ | 3 | C24c | no |
/// | M3 | $M_{4}$ | $k_{12}t4$ | 3 | C48a | no |
/// | M4 | $M_{2}$ | $k_{12}t3$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 6 | 1 |
/// | dir{P1} | 184 | P6cc | C6v^2 | 4 | 3 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 12 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 6 | 1 |
/// | dir{P1} | 168 | P6 | C6^1 | 8 | 3 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 12 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 6 | 1 |
/// | dir{P1} | 159 | P31c | C3v^4 | 8 | 3 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 12 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 6 | 1 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 8 | 3 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 12 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{17}t1t2$ | 2 | B4a | yes |
/// | A3A4 | $A_{4}A5$ | $k_{17}t3t4$ | 2 | B4a | yes |
/// | A5A5 | $A_{3}A3$ | $k_{17}t6t6$ | 4 | D12a | no |
/// | A6A6 | $A_{6}A6$ | $k_{17}t5t5$ | 4 | D12a | no |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 168 | P6 | C6^1 | 4 | 1 |
///
/// ### Isotropy subgroups for A3A4 ($A_{4}A5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 173 | P6_3 | C6^6 | 4 | 1 |
///
/// ### Isotropy subgroups for A5A5 ($A_{3}A3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 3 | P2 | C2^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A6A6 ($A_{6}A6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H2 | $H_{1}H2$ | $k_{15}t1t2$ | 4 | D24c | yes |
/// | H3H3 | $H_{3}H3$ | $k_{15}t3t3$ | 8 | F72a | no |
///
/// ### Isotropy subgroups for H1H2 ($H_{1}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 173 | P6_3 | C6^6 | 12 | 2 |
/// | dir{P1} | 168 | P6 | C6^1 | 12 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 36 | 2 |
/// | dir{C1} | 3 | P2 | C2^1 | 36 | 2 |
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
/// | dir{P1} | 184 | P6cc | C6v^2 | 3 | 2 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 168 | P6 | C6^1 | 6 | 2 |
/// | dir{P1} | 159 | P31c | C3v^4 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 6 | 2 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 9 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 18 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 18 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1L2 | $L_{1}L2$ | $k_{14}t1t2$ | 6 | E48b | yes |
/// | L3L4 | $L_{3}L4$ | $k_{14}t3t4$ | 6 | E48b | yes |
///
/// ### Isotropy subgroups for L1L2 ($L_{1}L2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 12 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 168 | P6 | C6^1 | 16 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L3L4 ($L_{3}L4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 12 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 173 | P6_3 | C6^6 | 16 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 3 |
///
pub struct Sg184;

/// # 185 P6₃cm (C₆ᵥ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{16}t6$ | 2 | B6a | no |
/// | GM6 | $\Gamma_{6}$ | $k_{16}t5$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 185 | P6_3cm | C6v^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 173 | P6_3 | C6^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 158 | P3c1 | C3v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 3 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
///
/// ### Isotropy subgroups for GM6 ($\Gamma_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{12}t2$ | 3 | C24c | no |
/// | M3 | $M_{4}$ | $k_{12}t4$ | 3 | C48a | no |
/// | M4 | $M_{2}$ | $k_{12}t3$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 6 | 1 |
/// | dir{P1} | 185 | P6_3cm | C6v^3 | 4 | 3 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 12 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 6 | 1 |
/// | dir{P1} | 173 | P6_3 | C6^6 | 8 | 3 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 12 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 6 | 1 |
/// | dir{P1} | 157 | P31m | C3v^2 | 8 | 3 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 12 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 6 | 1 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 8 | 3 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 12 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A3 | $A_{1}A4$ | $k_{17}t1t4$ | 2 | B4a | yes |
/// | A2A4 | $A_{2}A5$ | $k_{17}t2t3$ | 2 | B4a | yes |
/// | A5A6 | $A_{3}A6$ | $k_{17}t5t6$ | 4 | D24c | yes |
///
/// ### Isotropy subgroups for A1A3 ($A_{1}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 4 | 1 |
///
/// ### Isotropy subgroups for A2A4 ($A_{2}A5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 4 | 1 |
///
/// ### Isotropy subgroups for A5A6 ($A_{3}A6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H1 | $H_{1}H1$ | $k_{15}t1t1$ | 4 | D12a | no |
/// | H2H2 | $H_{2}H2$ | $k_{15}t2t2$ | 4 | D12a | no |
/// | H3H3 | $H_{3}H3$ | $k_{15}t3t3$ | 8 | F72a | no |
///
/// ### Isotropy subgroups for H1H1 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 156 | P3m1 | C3v^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H2H2 ($H_{2}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 158 | P3c1 | C3v^3 | 12 | 2 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 36 | 2 |
/// | dir{C1} | 8 | Cm | Cs^3 | 36 | 2 |
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
/// | dir{P1} | 186 | P6_3mc | C6v^4 | 3 | 2 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 173 | P6_3 | C6^6 | 6 | 2 |
/// | dir{P1} | 159 | P31c | C3v^4 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 6 | 2 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 9 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 18 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 18 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1L3 | $L_{1}L2$ | $k_{14}t1t4$ | 6 | E96c | yes |
/// | L2L4 | $L_{3}L4$ | $k_{14}t2t3$ | 6 | E96c | yes |
///
/// ### Isotropy subgroups for L1L3 ($L_{1}L2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{P1} | 157 | P31m | C3v^2 | 16 | 3 |
/// | dir{C1} | 8 | Cm | Cs^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ### Isotropy subgroups for L2L4 ($L_{3}L4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{P1} | 159 | P31c | C3v^4 | 16 | 3 |
/// | dir{C1} | 9 | Cc | Cs^4 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
pub struct Sg185;

/// # 186 P6₃mc (C₆ᵥ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{16}t6$ | 2 | B6a | no |
/// | GM6 | $\Gamma_{6}$ | $k_{16}t5$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 186 | P6_3mc | C6v^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 173 | P6_3 | C6^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 156 | P3m1 | C3v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 3 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
///
/// ### Isotropy subgroups for GM6 ($\Gamma_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{12}t2$ | 3 | C24c | no |
/// | M3 | $M_{4}$ | $k_{12}t4$ | 3 | C48a | no |
/// | M4 | $M_{2}$ | $k_{12}t3$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 6 | 1 |
/// | dir{P1} | 186 | P6_3mc | C6v^4 | 4 | 3 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 12 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 6 | 1 |
/// | dir{P1} | 173 | P6_3 | C6^6 | 8 | 3 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 12 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 6 | 1 |
/// | dir{P1} | 159 | P31c | C3v^4 | 8 | 3 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 12 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 6 | 1 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 8 | 3 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 12 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A4 | $A_{1}A4$ | $k_{17}t1t3$ | 2 | B4a | yes |
/// | A2A3 | $A_{2}A5$ | $k_{17}t2t4$ | 2 | B4a | yes |
/// | A5A6 | $A_{3}A6$ | $k_{17}t5t6$ | 4 | D24c | yes |
///
/// ### Isotropy subgroups for A1A4 ($A_{1}A4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 156 | P3m1 | C3v^1 | 4 | 1 |
///
/// ### Isotropy subgroups for A2A3 ($A_{2}A5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 158 | P3c1 | C3v^3 | 4 | 1 |
///
/// ### Isotropy subgroups for A5A6 ($A_{3}A6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H2 | $H_{1}H2$ | $k_{15}t1t2$ | 4 | D24c | yes |
/// | H3 | $H_{3}$ | $k_{15}t3$ | 4 | D72a | no |
///
/// ### Isotropy subgroups for H1H2 ($H_{1}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 12 | 2 |
/// | dir{P1} | 157 | P31m | C3v^2 | 12 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 12 | 2 |
/// | dir{P1} | 157 | P31m | C3v^2 | 12 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 36 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 36 | 2 |
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
/// | dir{P1} | 185 | P6_3cm | C6v^3 | 3 | 2 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 173 | P6_3 | C6^6 | 6 | 2 |
/// | dir{P1} | 157 | P31m | C3v^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 6 | 2 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 9 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 18 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 18 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 18 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1L4 | $L_{1}L2$ | $k_{14}t1t3$ | 6 | E96c | yes |
/// | L2L3 | $L_{3}L4$ | $k_{14}t2t4$ | 6 | E96c | yes |
///
/// ### Isotropy subgroups for L1L4 ($L_{1}L2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 16 | 3 |
/// | dir{C1} | 8 | Cm | Cs^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ### Isotropy subgroups for L2L3 ($L_{3}L4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 16 | 3 |
/// | dir{C1} | 9 | Cc | Cs^4 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
pub struct Sg186;

/// # 187 P-6m2 (D₃ₕ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{16}t5$ | 2 | B6a | no |
/// | GM6 | $\Gamma_{6}$ | $k_{16}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 156 | P3m1 | C3v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 3 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM6 ($\Gamma_{6}$)
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
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{4}$ | $k_{12}t3$ | 3 | C24c | no |
/// | M3 | $M_{2}$ | $k_{12}t4$ | 3 | C48a | no |
/// | M4 | $M_{3}$ | $k_{12}t2$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 6 | 1 |
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 4 | 3 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 12 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 6 | 1 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 8 | 3 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 12 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 6 | 1 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 8 | 3 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 12 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 6 | 1 |
/// | dir{P1} | 149 | P312 | D3^1 | 8 | 3 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 12 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{17}t3$ | 1 | A2a | no |
/// | A3 | $A_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | A4 | $A_{3}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A5 | $A_{5}$ | $k_{17}t5$ | 2 | B12a | no |
/// | A6 | $A_{6}$ | $k_{17}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 40 | Ama2 | C2v^16 | 6 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 6 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A6 ($A_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 40 | Ama2 | C2v^16 | 6 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 6 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t1$ | 2 | B12a | no |
/// | H2 | $H_{4}$ | $k_{15}t2$ | 2 | B12a | no |
/// | H3 | $H_{3}$ | $k_{15}t5$ | 2 | B12a | no |
/// | H4 | $H_{6}$ | $k_{15}t6$ | 2 | B12a | no |
/// | H5 | $H_{2}$ | $k_{15}t3$ | 2 | B12a | no |
/// | H6 | $H_{5}$ | $k_{15}t4$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 6 | 2 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 6 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 6 | 2 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 6 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 6 | 2 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 6 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H4 ($H_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 6 | 2 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 6 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H5 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 6 | 2 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 6 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H6 ($H_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 6 | 2 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 6 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{4}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3 | $K_{3}$ | $k_{13}t5$ | 2 | B6a | no |
/// | K4 | $K_{6}$ | $k_{13}t6$ | 2 | B12a | no |
/// | K5 | $K_{2}$ | $k_{13}t3$ | 2 | B6a | no |
/// | K6 | $K_{5}$ | $k_{13}t4$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 189 | P-62m | D3h^3 | 3 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 6 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 189 | P-62m | D3h^3 | 3 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K4 ($K_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 6 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K5 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 189 | P-62m | D3h^3 | 3 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K6 ($K_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 6 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C48a | no |
/// | L2 | $L_{4}$ | $k_{14}t3$ | 3 | C48a | no |
/// | L3 | $L_{2}$ | $k_{14}t4$ | 3 | C48a | no |
/// | L4 | $L_{3}$ | $k_{14}t2$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 44 | Imm2 | C2v^20 | 6 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 12 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 8 | 3 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 46 | Ima2 | C2v^22 | 6 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 12 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 8 | 3 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L3 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 44 | Imm2 | C2v^20 | 6 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 12 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 8 | 3 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L4 ($L_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 46 | Ima2 | C2v^22 | 6 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 12 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 8 | 3 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
pub struct Sg187;

/// # 188 P-6c2 (D₃ₕ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{3}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{16}t5$ | 2 | B6a | no |
/// | GM6 | $\Gamma_{6}$ | $k_{16}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 158 | P3c1 | C3v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 40 | Ama2 | C2v^16 | 3 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM6 ($\Gamma_{6}$)
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
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{4}$ | $k_{12}t3$ | 3 | C24c | no |
/// | M3 | $M_{2}$ | $k_{12}t4$ | 3 | C48a | no |
/// | M4 | $M_{3}$ | $k_{12}t2$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 6 | 1 |
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 4 | 3 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 12 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 6 | 1 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 8 | 3 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 12 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 6 | 1 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 8 | 3 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 12 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 6 | 1 |
/// | dir{P1} | 149 | P312 | D3^1 | 8 | 3 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 12 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{3}$ | $k_{17}t1$ | 2 | B8a | yes |
/// | A2A3 | $A_{1}A2$ | $k_{17}t2t3$ | 4 | D24e | no |
///
/// ### Isotropy subgroups for A1 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 4 | 1 |
/// | dir{P1} | 149 | P312 | D3^1 | 4 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 1 |
///
/// ### Isotropy subgroups for A2A3 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 6 | Pm | Cs^1 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H2 | $H_{1}H4$ | $k_{15}t1t2$ | 4 | D24e | no |
/// | H3H4 | $H_{3}H6$ | $k_{15}t5t6$ | 4 | D24e | no |
/// | H5H6 | $H_{2}H5$ | $k_{15}t3t4$ | 4 | D24e | no |
///
/// ### Isotropy subgroups for H1H2 ($H_{1}H4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
///
/// ### Isotropy subgroups for H3H4 ($H_{3}H6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
///
/// ### Isotropy subgroups for H5H6 ($H_{2}H5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{4}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3 | $K_{3}$ | $k_{13}t5$ | 2 | B6a | no |
/// | K4 | $K_{6}$ | $k_{13}t6$ | 2 | B12a | no |
/// | K5 | $K_{2}$ | $k_{13}t3$ | 2 | B6a | no |
/// | K6 | $K_{5}$ | $k_{13}t4$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 3 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 6 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 3 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K4 ($K_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 6 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K5 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 3 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K6 ($K_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 6 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 6 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 6 | E96d | yes |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 12 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 12 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{P3} | 174 | P-6 | C3h^1 | 16 | 3 |
/// | dir{C1} | 149 | P312 | D3^1 | 16 | 3 |
/// | dir{P1} | 143 | P3 | C3^1 | 32 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 48 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
pub struct Sg188;

/// # 189 P-62m (D₃ₕ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{16}t5$ | 2 | B6a | no |
/// | GM6 | $\Gamma_{6}$ | $k_{16}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 189 | P-62m | D3h^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 3 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM6 ($\Gamma_{6}$)
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
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{12}t2$ | 3 | C48a | no |
/// | M3 | $M_{2}$ | $k_{12}t4$ | 3 | C48a | no |
/// | M4 | $M_{4}$ | $k_{12}t3$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 6 | 1 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 4 | 3 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 12 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 6 | 1 |
/// | dir{P1} | 150 | P321 | D3^2 | 8 | 3 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 12 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 6 | 1 |
/// | dir{P1} | 157 | P31m | C3v^2 | 8 | 3 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 12 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 6 | 1 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 8 | 3 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 12 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | A2 | $A_{3}$ | $k_{17}t2$ | 1 | A2a | no |
/// | A3 | $A_{4}$ | $k_{17}t4$ | 1 | A2a | no |
/// | A4 | $A_{2}$ | $k_{17}t3$ | 1 | A2a | no |
/// | A5 | $A_{5}$ | $k_{17}t5$ | 2 | B12a | no |
/// | A6 | $A_{6}$ | $k_{17}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 189 | P-62m | D3h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 189 | P-62m | D3h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A4 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 40 | Ama2 | C2v^16 | 6 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 6 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A6 ($A_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 40 | Ama2 | C2v^16 | 6 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 6 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H1 | $H_{1}H1$ | $k_{15}t1t1*$ | 2 | B6b | no |
/// | H2H2 | $H_{3}H3$ | $k_{15}t2t2*$ | 2 | B6b | no |
/// | H3H3 | $H_{4}H4$ | $k_{15}t4t4*$ | 2 | B6b | no |
/// | H4H4 | $H_{2}H2$ | $k_{15}t3t3*$ | 2 | B6b | no |
/// | H5H5 | $H_{5}H5$ | $k_{15}t5t5*$ | 4 | D36b | yes |
/// | H6H6 | $H_{6}H6$ | $k_{15}t6t6*$ | 4 | D36b | yes |
///
/// ### Isotropy subgroups for H1H1 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for H2H2 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 6 | 1 |
///
/// ### Isotropy subgroups for H3H3 ($H_{4}H4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for H4H4 ($H_{2}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 6 | 1 |
///
/// ### Isotropy subgroups for H5H5 ($H_{5}H5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 1 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 18 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 18 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 36 | 1 |
///
/// ### Isotropy subgroups for H6H6 ($H_{6}H6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 1 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 18 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 18 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 36 | 1 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1K1 | $K_{1}K1$ | $k_{13}t1t1*$ | 2 | B3a | no |
/// | K2K2 | $K_{3}K3$ | $k_{13}t2t2*$ | 2 | B6b | no |
/// | K3K3 | $K_{4}K4$ | $k_{13}t4t4*$ | 2 | B6b | no |
/// | K4K4 | $K_{2}K2$ | $k_{13}t3t3*$ | 2 | B6b | no |
/// | K5K5 | $K_{5}K5$ | $k_{13}t5t5*$ | 4 | D18a | yes |
/// | K6K6 | $K_{6}K6$ | $k_{13}t6t6*$ | 4 | D36b | yes |
///
/// ### Isotropy subgroups for K1K1 ($K_{1}K1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 3 | 1 |
///
/// ### Isotropy subgroups for K2K2 ($K_{3}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 6 | 1 |
///
/// ### Isotropy subgroups for K3K3 ($K_{4}K4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 156 | P3m1 | C3v^1 | 6 | 1 |
///
/// ### Isotropy subgroups for K4K4 ($K_{2}K2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for K5K5 ($K_{5}K5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 9 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 18 | 1 |
///
/// ### Isotropy subgroups for K6K6 ($K_{6}K6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 18 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 3 | C48a | no |
/// | L2 | $L_{3}$ | $k_{14}t2$ | 3 | C48a | no |
/// | L3 | $L_{2}$ | $k_{14}t4$ | 3 | C48a | no |
/// | L4 | $L_{4}$ | $k_{14}t3$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 44 | Imm2 | C2v^20 | 6 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 12 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 8 | 3 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 46 | Ima2 | C2v^22 | 6 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 12 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{P1} | 190 | P-62c | D3h^4 | 8 | 3 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L3 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 44 | Imm2 | C2v^20 | 6 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 12 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 8 | 3 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L4 ($L_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 46 | Ima2 | C2v^22 | 6 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 12 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{P1} | 190 | P-62c | D3h^4 | 8 | 3 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
pub struct Sg189;

/// # 190 P-62c (D₃ₕ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM5 | $\Gamma_{5}$ | $k_{16}t5$ | 2 | B6a | no |
/// | GM6 | $\Gamma_{6}$ | $k_{16}t6$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 150 | P321 | D3^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 40 | Ama2 | C2v^16 | 3 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM6 ($\Gamma_{6}$)
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
/// | M1 | $M_{1}$ | $k_{12}t1$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{12}t2$ | 3 | C48a | no |
/// | M3 | $M_{2}$ | $k_{12}t4$ | 3 | C48a | no |
/// | M4 | $M_{4}$ | $k_{12}t3$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 6 | 1 |
/// | dir{P1} | 190 | P-62c | D3h^4 | 4 | 3 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 12 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 6 | 1 |
/// | dir{P1} | 150 | P321 | D3^2 | 8 | 3 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 12 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 6 | 1 |
/// | dir{P1} | 159 | P31c | C3v^4 | 8 | 3 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 12 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 6 | 1 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 8 | 3 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 12 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 3 |
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
/// | dir{P1} | 6 | Pm | Cs^1 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 4 | 1 |
/// | dir{P1} | 150 | P321 | D3^2 | 4 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H1 | $H_{1}H1$ | $k_{15}t3t3*$ | 4 | D72b | no |
/// | H2H2 | $H_{2}H2$ | $k_{15}t2t2*$ | 4 | D72b | no |
/// | H3H3 | $H_{3}H3$ | $k_{15}t1t1*$ | 4 | D24d | yes |
///
/// ### Isotropy subgroups for H1H1 ($H_{1}H1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 1 |
///
/// ### Isotropy subgroups for H2H2 ($H_{2}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 1 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 1 |
/// | dir{P1} | 149 | P312 | D3^1 | 12 | 1 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1K1 | $K_{1}K1$ | $k_{13}t1t1*$ | 2 | B3a | no |
/// | K2K2 | $K_{3}K3$ | $k_{13}t2t2*$ | 2 | B6b | no |
/// | K3K3 | $K_{4}K4$ | $k_{13}t4t4*$ | 2 | B6b | no |
/// | K4K4 | $K_{2}K2$ | $k_{13}t3t3*$ | 2 | B6b | no |
/// | K5K5 | $K_{5}K5$ | $k_{13}t5t5*$ | 4 | D18a | yes |
/// | K6K6 | $K_{6}K6$ | $k_{13}t6t6*$ | 4 | D36b | yes |
///
/// ### Isotropy subgroups for K1K1 ($K_{1}K1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 3 | 1 |
///
/// ### Isotropy subgroups for K2K2 ($K_{3}K3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 149 | P312 | D3^1 | 6 | 1 |
///
/// ### Isotropy subgroups for K3K3 ($K_{4}K4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 158 | P3c1 | C3v^3 | 6 | 1 |
///
/// ### Isotropy subgroups for K4K4 ($K_{2}K2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for K5K5 ($K_{5}K5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 6 | 1 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 9 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 18 | 1 |
///
/// ### Isotropy subgroups for K6K6 ($K_{6}K6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 143 | P3 | C3^1 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 18 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 18 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 36 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 6 | E96d | yes |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 12 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 12 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{P3} | 174 | P-6 | C3h^1 | 16 | 3 |
/// | dir{C1} | 150 | P321 | D3^2 | 16 | 3 |
/// | dir{P1} | 143 | P3 | C3^1 | 32 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 48 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
pub struct Sg190;

/// # 191 P6/mmm (D₆ₕ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{16}t5$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{16}t6$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{3}^+$ | $k_{16}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{3}^-$ | $k_{16}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{6}^+$ | $k_{16}t11$ | 2 | B6a | no |
/// | GM5- | $\Gamma_{6}^-$ | $k_{16}t12$ | 2 | B12a | no |
/// | GM6+ | $\Gamma_{5}^+$ | $k_{16}t9$ | 2 | B12a | no |
/// | GM6- | $\Gamma_{5}^-$ | $k_{16}t10$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 191 | P6/mmm | D6h^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 177 | P622 | D6^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 175 | P6/m | C6h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 183 | P6mm | C6v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 189 | P-62m | D3h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 162 | P-31m | D3d^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{6}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 3 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{6}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 6 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM6+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM6- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 6 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 6 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{12}t1$ | 3 | C24b | no |
/// | M1- | $M_{1}^-$ | $k_{12}t2$ | 3 | C48a | no |
/// | M2+ | $M_{3}^+$ | $k_{12}t3$ | 3 | C24c | no |
/// | M2- | $M_{3}^-$ | $k_{12}t4$ | 3 | C48a | no |
/// | M3+ | $M_{2}^+$ | $k_{12}t5$ | 3 | C48a | no |
/// | M3- | $M_{2}^-$ | $k_{12}t6$ | 3 | C48a | no |
/// | M4+ | $M_{4}^+$ | $k_{12}t7$ | 3 | C48a | no |
/// | M4- | $M_{4}^-$ | $k_{12}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 47 | Pmmm | D2h^1 | 6 | 1 |
/// | dir{P1} | 191 | P6/mmm | D6h^1 | 4 | 3 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 12 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 50 | Pban | D2h^4 | 6 | 1 |
/// | dir{P1} | 177 | P622 | D6^1 | 8 | 3 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 12 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 24 | 3 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M2+ ($M_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 55 | Pbam | D2h^9 | 6 | 1 |
/// | dir{P1} | 175 | P6/m | C6h^1 | 8 | 3 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 12 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2- ($M_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 59 | Pmmn | D2h^13 | 6 | 1 |
/// | dir{P1} | 183 | P6mm | C6v^1 | 8 | 3 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 12 | 3 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 24 | 3 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M3+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 53 | Pmna | D2h^7 | 6 | 1 |
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 8 | 3 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 12 | 3 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 24 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M3- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 6 | 1 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 8 | 3 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 12 | 3 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 24 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4+ ($M_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 53 | Pmna | D2h^7 | 6 | 1 |
/// | dir{P1} | 162 | P-31m | D3d^1 | 8 | 3 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 12 | 3 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 24 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4- ($M_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 6 | 1 |
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 8 | 3 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 12 | 3 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 24 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1+ | $A_{1}^+$ | $k_{17}t1$ | 1 | A2a | no |
/// | A1- | $A_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | A2+ | $A_{2}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | A2- | $A_{2}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | A3+ | $A_{4}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | A3- | $A_{4}^-$ | $k_{17}t6$ | 1 | A2a | no |
/// | A4+ | $A_{3}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | A4- | $A_{3}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | A5+ | $A_{6}^+$ | $k_{17}t11$ | 2 | B12a | no |
/// | A5- | $A_{6}^-$ | $k_{17}t12$ | 2 | B12a | no |
/// | A6+ | $A_{5}^+$ | $k_{17}t9$ | 2 | B12a | no |
/// | A6- | $A_{5}^-$ | $k_{17}t10$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for A1+ ($A_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 191 | P6/mmm | D6h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A1- ($A_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 192 | P6/mcc | D6h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A2+ ($A_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 192 | P6/mcc | D6h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A2- ($A_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 191 | P6/mmm | D6h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A3+ ($A_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 194 | P6_3/mmc | D6h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A3- ($A_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 193 | P6_3/mcm | D6h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A4+ ($A_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 193 | P6_3/mcm | D6h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A4- ($A_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 194 | P6_3/mmc | D6h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A5+ ($A_{6}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 66 | Cccm | D2h^20 | 6 | 1 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 6 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A5- ($A_{6}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 66 | Cccm | D2h^20 | 6 | 1 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 6 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 12 | 1 |
///
/// ### Isotropy subgroups for A6+ ($A_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 6 | 1 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 6 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 12 | 1 |
///
/// ### Isotropy subgroups for A6- ($A_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 6 | 1 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 6 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{15}t1$ | 2 | B12a | no |
/// | H2 | $H_{3}$ | $k_{15}t2$ | 2 | B12a | no |
/// | H3 | $H_{4}$ | $k_{15}t4$ | 2 | B12a | no |
/// | H4 | $H_{2}$ | $k_{15}t3$ | 2 | B12a | no |
/// | H5 | $H_{5}$ | $k_{15}t5$ | 4 | D72d | yes |
/// | H6 | $H_{6}$ | $k_{15}t6$ | 4 | D72d | yes |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 194 | P6_3/mmc | D6h^4 | 6 | 2 |
/// | dir{P1} | 191 | P6/mmm | D6h^1 | 6 | 2 |
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 193 | P6_3/mcm | D6h^3 | 6 | 2 |
/// | dir{P1} | 192 | P6/mcc | D6h^2 | 6 | 2 |
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 12 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 194 | P6_3/mmc | D6h^4 | 6 | 2 |
/// | dir{P1} | 191 | P6/mmm | D6h^1 | 6 | 2 |
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 12 | 2 |
///
/// ### Isotropy subgroups for H4 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 193 | P6_3/mcm | D6h^3 | 6 | 2 |
/// | dir{P1} | 192 | P6/mcc | D6h^2 | 6 | 2 |
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 12 | 2 |
///
/// ### Isotropy subgroups for H5 ($H_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 12 | 2 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 12 | 2 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 18 | 2 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 18 | 2 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 18 | 2 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 18 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 24 | 2 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 36 | 2 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 36 | 2 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 36 | 2 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 36 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 36 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 36 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 72 | 2 |
///
/// ### Isotropy subgroups for H6 ($H_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 12 | 2 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 12 | 2 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 18 | 2 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 18 | 2 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 18 | 2 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 18 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 24 | 2 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 36 | 2 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 36 | 2 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 36 | 2 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 36 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 36 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 36 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{3}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3 | $K_{4}$ | $k_{13}t4$ | 2 | B12a | no |
/// | K4 | $K_{2}$ | $k_{13}t3$ | 2 | B12a | no |
/// | K5 | $K_{5}$ | $k_{13}t5$ | 4 | D36c | yes |
/// | K6 | $K_{6}$ | $k_{13}t6$ | 4 | D72d | yes |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 191 | P6/mmm | D6h^1 | 3 | 2 |
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 177 | P622 | D6^1 | 6 | 2 |
/// | dir{P1} | 162 | P-31m | D3d^1 | 6 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 183 | P6mm | C6v^1 | 6 | 2 |
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 6 | 2 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K4 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 189 | P-62m | D3h^3 | 6 | 2 |
/// | dir{P1} | 175 | P6/m | C6h^1 | 6 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K5 ($K_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 189 | P-62m | D3h^3 | 6 | 2 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 9 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 18 | 2 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 18 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 18 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 36 | 2 |
///
/// ### Isotropy subgroups for K6 ($K_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 12 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 18 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 18 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 18 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 18 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 36 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 36 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^+$ | $k_{14}t1$ | 3 | C48a | no |
/// | L1- | $L_{1}^-$ | $k_{14}t2$ | 3 | C48a | no |
/// | L2+ | $L_{3}^+$ | $k_{14}t3$ | 3 | C48a | no |
/// | L2- | $L_{3}^-$ | $k_{14}t4$ | 3 | C48a | no |
/// | L3+ | $L_{2}^+$ | $k_{14}t5$ | 3 | C48a | no |
/// | L3- | $L_{2}^-$ | $k_{14}t6$ | 3 | C48a | no |
/// | L4+ | $L_{4}^+$ | $k_{14}t7$ | 3 | C48a | no |
/// | L4- | $L_{4}^-$ | $k_{14}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 71 | Immm | D2h^25 | 6 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 12 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 191 | P6/mmm | D6h^1 | 8 | 3 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 24 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 72 | Ibam | D2h^26 | 6 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 12 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 192 | P6/mcc | D6h^2 | 8 | 3 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 24 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2+ ($L_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 72 | Ibam | D2h^26 | 6 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 12 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 192 | P6/mcc | D6h^2 | 8 | 3 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 24 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L2- ($L_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 71 | Immm | D2h^25 | 6 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 12 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 191 | P6/mmm | D6h^1 | 8 | 3 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 24 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 48 | 3 |
///
/// ### Isotropy subgroups for L3+ ($L_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 74 | Imma | D2h^28 | 6 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 12 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 194 | P6_3/mmc | D6h^4 | 8 | 3 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 24 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 48 | 3 |
///
/// ### Isotropy subgroups for L3- ($L_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 74 | Imma | D2h^28 | 6 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 12 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 193 | P6_3/mcm | D6h^3 | 8 | 3 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 24 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 48 | 3 |
///
/// ### Isotropy subgroups for L4+ ($L_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 74 | Imma | D2h^28 | 6 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 12 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 193 | P6_3/mcm | D6h^3 | 8 | 3 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 24 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 48 | 3 |
///
/// ### Isotropy subgroups for L4- ($L_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 74 | Imma | D2h^28 | 6 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 12 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 194 | P6_3/mmc | D6h^4 | 8 | 3 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 24 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 48 | 3 |
///
pub struct Sg191;

/// # 192 P6/mcc (D₆ₕ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{16}t5$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{16}t6$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{3}^+$ | $k_{16}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{3}^-$ | $k_{16}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{6}^+$ | $k_{16}t11$ | 2 | B6a | no |
/// | GM5- | $\Gamma_{6}^-$ | $k_{16}t12$ | 2 | B12a | no |
/// | GM6+ | $\Gamma_{5}^+$ | $k_{16}t9$ | 2 | B12a | no |
/// | GM6- | $\Gamma_{5}^-$ | $k_{16}t10$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 192 | P6/mcc | D6h^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 177 | P622 | D6^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 175 | P6/m | C6h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 184 | P6cc | C6v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 163 | P-31c | D3d^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{6}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 66 | Cccm | D2h^20 | 3 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{6}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 6 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM6+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM6- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 40 | Ama2 | C2v^16 | 6 | 1 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 6 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{12}t1$ | 3 | C24b | no |
/// | M1- | $M_{1}^-$ | $k_{12}t2$ | 3 | C48a | no |
/// | M2+ | $M_{3}^+$ | $k_{12}t3$ | 3 | C24c | no |
/// | M2- | $M_{3}^-$ | $k_{12}t4$ | 3 | C48a | no |
/// | M3+ | $M_{2}^+$ | $k_{12}t5$ | 3 | C48a | no |
/// | M3- | $M_{2}^-$ | $k_{12}t6$ | 3 | C48a | no |
/// | M4+ | $M_{4}^+$ | $k_{12}t7$ | 3 | C48a | no |
/// | M4- | $M_{4}^-$ | $k_{12}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 6 | 1 |
/// | dir{P1} | 192 | P6/mcc | D6h^2 | 4 | 3 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 12 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 48 | Pnnn | D2h^2 | 6 | 1 |
/// | dir{P1} | 177 | P622 | D6^1 | 8 | 3 |
/// | dir{P1} | 68 | Ccca | D2h^22 | 12 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 24 | 3 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M2+ ($M_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 58 | Pnnm | D2h^12 | 6 | 1 |
/// | dir{P1} | 175 | P6/m | C6h^1 | 8 | 3 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 12 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 24 | 3 |
///
/// ### Isotropy subgroups for M2- ($M_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 56 | Pccn | D2h^10 | 6 | 1 |
/// | dir{P1} | 184 | P6cc | C6v^2 | 8 | 3 |
/// | dir{P1} | 68 | Ccca | D2h^22 | 12 | 3 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 24 | 3 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 24 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M3+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 52 | Pnna | D2h^6 | 6 | 1 |
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 8 | 3 |
/// | dir{P1} | 68 | Ccca | D2h^22 | 12 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 3 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M3- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 53 | Pmna | D2h^7 | 6 | 1 |
/// | dir{P1} | 190 | P-62c | D3h^4 | 8 | 3 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 12 | 3 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 24 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4+ ($M_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 52 | Pnna | D2h^6 | 6 | 1 |
/// | dir{P1} | 163 | P-31c | D3d^2 | 8 | 3 |
/// | dir{P1} | 68 | Ccca | D2h^22 | 12 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 3 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4- ($M_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 53 | Pmna | D2h^7 | 6 | 1 |
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 8 | 3 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 12 | 3 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 24 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{17}t5t6$ | 4 | D24e | no |
/// | A3A4 | $A_{4}A5$ | $k_{17}t3t4$ | 4 | D24e | no |
/// | A5 | $A_{3}$ | $k_{17}t1$ | 2 | B8a | yes |
/// | A6 | $A_{6}$ | $k_{17}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 12 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 1 |
///
/// ### Isotropy subgroups for A3A4 ($A_{4}A5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 12 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
///
/// ### Isotropy subgroups for A5 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 177 | P622 | D6^1 | 4 | 1 |
/// | dir{P1} | 175 | P6/m | C6h^1 | 4 | 1 |
/// | dir{P1} | 168 | P6 | C6^1 | 8 | 1 |
///
/// ### Isotropy subgroups for A6 ($A_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 182 | P6_322 | D6^6 | 4 | 1 |
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 4 | 1 |
/// | dir{P1} | 173 | P6_3 | C6^6 | 8 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H2 | $H_{1}H2$ | $k_{15}t2t3$ | 8 | F144a | no |
/// | H3 | $H_{3}$ | $k_{15}t1$ | 4 | D48e | yes |
///
/// ### Isotropy subgroups for H1H2 ($H_{1}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 24 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 24 | 2 |
/// | dir{P3} | 21 | C222 | D2^6 | 36 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 36 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 36 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 36 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 48 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 72 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 72 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 72 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 72 | 2 |
/// | dir{C1} | 3 | P2 | C2^1 | 72 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 72 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 144 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 182 | P6_322 | D6^6 | 12 | 2 |
/// | dir{P1} | 177 | P622 | D6^1 | 12 | 2 |
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 12 | 2 |
/// | dir{P1} | 175 | P6/m | C6h^1 | 12 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 24 | 2 |
/// | dir{P1} | 173 | P6_3 | C6^6 | 24 | 2 |
/// | dir{P1} | 168 | P6 | C6^1 | 24 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 24 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 24 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 24 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 48 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{3}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3 | $K_{4}$ | $k_{13}t4$ | 2 | B12a | no |
/// | K4 | $K_{2}$ | $k_{13}t3$ | 2 | B12a | no |
/// | K5 | $K_{5}$ | $k_{13}t5$ | 4 | D36c | yes |
/// | K6 | $K_{6}$ | $k_{13}t6$ | 4 | D72d | yes |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 192 | P6/mcc | D6h^2 | 3 | 2 |
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 177 | P622 | D6^1 | 6 | 2 |
/// | dir{P1} | 163 | P-31c | D3d^2 | 6 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 184 | P6cc | C6v^2 | 6 | 2 |
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 6 | 2 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 12 | 2 |
///
/// ### Isotropy subgroups for K4 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 6 | 2 |
/// | dir{P1} | 175 | P6/m | C6h^1 | 6 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K5 ($K_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 6 | 2 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 9 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 18 | 2 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 18 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 18 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 36 | 2 |
///
/// ### Isotropy subgroups for K6 ($K_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 12 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 18 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 18 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 18 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 18 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 36 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 36 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 6 | E96d | yes |
/// | L2 | $L_{2}$ | $k_{14}t2$ | 6 | E96d | yes |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 12 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 12 | 2 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 12 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 24 | 2 |
/// | dir{C1} | 22 | F222 | D2^7 | 24 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{C1} | 177 | P622 | D6^1 | 16 | 3 |
/// | dir{P3} | 175 | P6/m | C6h^1 | 16 | 3 |
/// | dir{P1} | 168 | P6 | C6^1 | 32 | 3 |
/// | dir{P1} | 21 | C222 | D2^6 | 48 | 3 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 48 | 3 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 48 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 96 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 12 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 12 | 2 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 12 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 24 | 2 |
/// | dir{C1} | 22 | F222 | D2^7 | 24 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{C1} | 182 | P6_322 | D6^6 | 16 | 3 |
/// | dir{P3} | 176 | P6_3/m | C6h^2 | 16 | 3 |
/// | dir{P1} | 173 | P6_3 | C6^6 | 32 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 48 | 3 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 48 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 48 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 96 | 3 |
///
pub struct Sg192;

/// # 193 P6₃/mcm (D₆ₕ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{16}t5$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{16}t6$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{3}^+$ | $k_{16}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{3}^-$ | $k_{16}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{6}^+$ | $k_{16}t11$ | 2 | B6a | no |
/// | GM5- | $\Gamma_{6}^-$ | $k_{16}t12$ | 2 | B12a | no |
/// | GM6+ | $\Gamma_{5}^+$ | $k_{16}t9$ | 2 | B12a | no |
/// | GM6- | $\Gamma_{5}^-$ | $k_{16}t10$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 193 | P6_3/mcm | D6h^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 182 | P6_322 | D6^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 185 | P6_3cm | C6v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 189 | P-62m | D3h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 162 | P-31m | D3d^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{6}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 3 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 6 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{6}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 6 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 6 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 1 |
///
/// ### Isotropy subgroups for GM6+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM6- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 6 | 1 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 6 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{12}t1$ | 3 | C24b | no |
/// | M1- | $M_{1}^-$ | $k_{12}t2$ | 3 | C48a | no |
/// | M2+ | $M_{3}^+$ | $k_{12}t3$ | 3 | C24c | no |
/// | M2- | $M_{3}^-$ | $k_{12}t4$ | 3 | C48a | no |
/// | M3+ | $M_{2}^+$ | $k_{12}t5$ | 3 | C48a | no |
/// | M3- | $M_{2}^-$ | $k_{12}t6$ | 3 | C48a | no |
/// | M4+ | $M_{4}^+$ | $k_{12}t7$ | 3 | C48a | no |
/// | M4- | $M_{4}^-$ | $k_{12}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 6 | 1 |
/// | dir{P1} | 193 | P6_3/mcm | D6h^3 | 4 | 3 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 52 | Pnna | D2h^6 | 6 | 1 |
/// | dir{P1} | 182 | P6_322 | D6^6 | 8 | 3 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 3 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 3 |
///
/// ### Isotropy subgroups for M2+ ($M_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 62 | Pnma | D2h^16 | 6 | 1 |
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 8 | 3 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M2- ($M_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 62 | Pnma | D2h^16 | 6 | 1 |
/// | dir{P1} | 185 | P6_3cm | C6v^3 | 8 | 3 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 3 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 24 | 3 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 3 |
///
/// ### Isotropy subgroups for M3+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 60 | Pbcn | D2h^14 | 6 | 1 |
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 8 | 3 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 3 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M3- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 59 | Pmmn | D2h^13 | 6 | 1 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 8 | 3 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 3 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 24 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4+ ($M_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 58 | Pnnm | D2h^12 | 6 | 1 |
/// | dir{P1} | 162 | P-31m | D3d^1 | 8 | 3 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 3 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 24 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4- ($M_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 57 | Pbcm | D2h^11 | 6 | 1 |
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 8 | 3 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 3 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 24 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 2 | B8a | yes |
/// | A2 | $A_{2}$ | $k_{17}t2$ | 2 | B8a | yes |
/// | A3 | $A_{3}$ | $k_{17}t3$ | 4 | D48e | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 189 | P-62m | D3h^3 | 4 | 1 |
/// | dir{P1} | 162 | P-31m | D3d^1 | 4 | 1 |
/// | dir{P1} | 157 | P31m | C3v^2 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 4 | 1 |
/// | dir{P1} | 163 | P-31c | D3d^2 | 4 | 1 |
/// | dir{P1} | 159 | P31c | C3v^4 | 8 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 40 | Ama2 | C2v^16 | 12 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 12 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H3 | $H_{1}H4$ | $k_{15}t1t4$ | 4 | D24e | no |
/// | H2H4 | $H_{2}H3$ | $k_{15}t2t3$ | 4 | D24e | no |
/// | H5H6 | $H_{5}H6$ | $k_{15}t5t6$ | 8 | F144a | no |
///
/// ### Isotropy subgroups for H1H3 ($H_{1}H4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 12 | 2 |
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 12 | 2 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 24 | 2 |
///
/// ### Isotropy subgroups for H2H4 ($H_{2}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 12 | 2 |
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 12 | 2 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 24 | 2 |
///
/// ### Isotropy subgroups for H5H6 ($H_{5}H6$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 174 | P-6 | C3h^1 | 24 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 24 | 2 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 36 | 2 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 36 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 36 | 2 |
/// | dir{P3} | 12 | C2/m | C2h^3 | 36 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 72 | 2 |
/// | dir{C1} | 8 | Cm | Cs^3 | 72 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 72 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 72 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 72 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 72 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 144 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{3}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3 | $K_{4}$ | $k_{13}t4$ | 2 | B12a | no |
/// | K4 | $K_{2}$ | $k_{13}t3$ | 2 | B12a | no |
/// | K5 | $K_{5}$ | $k_{13}t5$ | 4 | D36c | yes |
/// | K6 | $K_{6}$ | $k_{13}t6$ | 4 | D72d | yes |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 194 | P6_3/mmc | D6h^4 | 3 | 2 |
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 182 | P6_322 | D6^6 | 6 | 2 |
/// | dir{P1} | 163 | P-31c | D3d^2 | 6 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 186 | P6_3mc | C6v^4 | 6 | 2 |
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 6 | 2 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K4 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 6 | 2 |
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 6 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K5 ($K_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 6 | 2 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 9 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 18 | 2 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 18 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 18 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 36 | 2 |
///
/// ### Isotropy subgroups for K6 ($K_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 159 | P31c | C3v^4 | 12 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 18 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 18 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 18 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 18 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 36 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 36 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 6 | E192a | yes |
/// | L2 | $L_{2}$ | $k_{14}t2$ | 6 | E192a | yes |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 44 | Imm2 | C2v^20 | 12 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P3} | 43 | Fdd2 | C2v^19 | 24 | 2 |
/// | dir{C1} | 42 | Fmm2 | C2v^18 | 24 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{C1} | 189 | P-62m | D3h^3 | 16 | 3 |
/// | dir{P3} | 162 | P-31m | D3d^1 | 16 | 3 |
/// | dir{P1} | 157 | P31m | C3v^2 | 32 | 3 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 48 | 3 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 48 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 96 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 96 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 46 | Ima2 | C2v^22 | 12 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P3} | 43 | Fdd2 | C2v^19 | 24 | 2 |
/// | dir{C1} | 42 | Fmm2 | C2v^18 | 24 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{C1} | 190 | P-62c | D3h^4 | 16 | 3 |
/// | dir{P3} | 163 | P-31c | D3d^2 | 16 | 3 |
/// | dir{P1} | 159 | P31c | C3v^4 | 32 | 3 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 48 | 3 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 48 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 96 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 96 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
pub struct Sg193;

/// # 194 P6₃/mmc (D₆ₕ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{16}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{16}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{16}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{16}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{16}t5$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{16}t6$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{3}^+$ | $k_{16}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{3}^-$ | $k_{16}t8$ | 1 | A2a | no |
/// | GM5+ | $\Gamma_{6}^+$ | $k_{16}t11$ | 2 | B6a | no |
/// | GM5- | $\Gamma_{6}^-$ | $k_{16}t12$ | 2 | B12a | no |
/// | GM6+ | $\Gamma_{5}^+$ | $k_{16}t9$ | 2 | B12a | no |
/// | GM6- | $\Gamma_{5}^-$ | $k_{16}t10$ | 2 | B12a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 194 | P6_3/mmc | D6h^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 182 | P6_322 | D6^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 186 | P6_3mc | C6v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 163 | P-31c | D3d^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{6}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 3 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 6 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{6}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 6 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 6 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 1 |
///
/// ### Isotropy subgroups for GM6+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM6- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 6 | 1 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 6 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{12}t1$ | 3 | C24b | no |
/// | M1- | $M_{1}^-$ | $k_{12}t2$ | 3 | C48a | no |
/// | M2+ | $M_{3}^+$ | $k_{12}t3$ | 3 | C24c | no |
/// | M2- | $M_{3}^-$ | $k_{12}t4$ | 3 | C48a | no |
/// | M3+ | $M_{2}^+$ | $k_{12}t5$ | 3 | C48a | no |
/// | M3- | $M_{2}^-$ | $k_{12}t6$ | 3 | C48a | no |
/// | M4+ | $M_{4}^+$ | $k_{12}t7$ | 3 | C48a | no |
/// | M4- | $M_{4}^-$ | $k_{12}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 6 | 1 |
/// | dir{P1} | 194 | P6_3/mmc | D6h^4 | 4 | 3 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 52 | Pnna | D2h^6 | 6 | 1 |
/// | dir{P1} | 182 | P6_322 | D6^6 | 8 | 3 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 3 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 3 |
///
/// ### Isotropy subgroups for M2+ ($M_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 62 | Pnma | D2h^16 | 6 | 1 |
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 8 | 3 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 3 |
///
/// ### Isotropy subgroups for M2- ($M_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 62 | Pnma | D2h^16 | 6 | 1 |
/// | dir{P1} | 186 | P6_3mc | C6v^4 | 8 | 3 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 3 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 24 | 3 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 24 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 3 |
///
/// ### Isotropy subgroups for M3+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 58 | Pnnm | D2h^12 | 6 | 1 |
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 8 | 3 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 3 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 24 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M3- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 57 | Pbcm | D2h^11 | 6 | 1 |
/// | dir{P1} | 190 | P-62c | D3h^4 | 8 | 3 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 3 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 24 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4+ ($M_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 60 | Pbcn | D2h^14 | 6 | 1 |
/// | dir{P1} | 163 | P-31c | D3d^2 | 8 | 3 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 3 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M4- ($M_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 59 | Pmmn | D2h^13 | 6 | 1 |
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 8 | 3 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 3 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 24 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{A}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{17}t1$ | 2 | B8a | yes |
/// | A2 | $A_{2}$ | $k_{17}t2$ | 2 | B8a | yes |
/// | A3 | $A_{3}$ | $k_{17}t3$ | 4 | D48e | yes |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 187 | P-6m2 | D3h^1 | 4 | 1 |
/// | dir{P1} | 164 | P-3m1 | D3d^3 | 4 | 1 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 8 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 4 | 1 |
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 4 | 1 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 8 | 1 |
///
/// ### Isotropy subgroups for A3 ($A_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 40 | Ama2 | C2v^16 | 12 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 12 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{H}\ (1/3,1/3,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{2}$ | $k_{15}t3$ | 4 | D144a | no |
/// | H2 | $H_{1}$ | $k_{15}t2$ | 4 | D144a | no |
/// | H3 | $H_{3}$ | $k_{15}t1$ | 4 | D48e | yes |
///
/// ### Isotropy subgroups for H1 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 12 | 2 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 12 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 36 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 36 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 72 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 72 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 72 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 72 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 72 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 144 | 2 |
///
/// ### Isotropy subgroups for H2 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 12 | 2 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 12 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 36 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 36 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 72 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 72 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 72 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 72 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 72 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 144 | 2 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 190 | P-62c | D3h^4 | 12 | 2 |
/// | dir{P1} | 189 | P-62m | D3h^3 | 12 | 2 |
/// | dir{P1} | 163 | P-31c | D3d^2 | 12 | 2 |
/// | dir{P1} | 162 | P-31m | D3d^1 | 12 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 24 | 2 |
/// | dir{P1} | 159 | P31c | C3v^4 | 24 | 2 |
/// | dir{P1} | 157 | P31m | C3v^2 | 24 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 24 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 24 | 2 |
/// | dir{P1} | 147 | P-3 | C3i^1 | 24 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 48 | 2 |
///
/// ## Irreps at $\mathrm{K}\ (1/3,1/3,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | K1 | $K_{1}$ | $k_{13}t1$ | 2 | B6a | no |
/// | K2 | $K_{3}$ | $k_{13}t2$ | 2 | B12a | no |
/// | K3 | $K_{4}$ | $k_{13}t4$ | 2 | B12a | no |
/// | K4 | $K_{2}$ | $k_{13}t3$ | 2 | B12a | no |
/// | K5 | $K_{5}$ | $k_{13}t5$ | 4 | D36c | yes |
/// | K6 | $K_{6}$ | $k_{13}t6$ | 4 | D72d | yes |
///
/// ### Isotropy subgroups for K1 ($K_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 193 | P6_3/mcm | D6h^3 | 3 | 2 |
/// | dir{P1} | 188 | P-6c2 | D3h^2 | 6 | 2 |
///
/// ### Isotropy subgroups for K2 ($K_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 182 | P6_322 | D6^6 | 6 | 2 |
/// | dir{P1} | 162 | P-31m | D3d^1 | 6 | 2 |
/// | dir{P1} | 149 | P312 | D3^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K3 ($K_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 185 | P6_3cm | C6v^3 | 6 | 2 |
/// | dir{P1} | 165 | P-3c1 | D3d^4 | 6 | 2 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 12 | 2 |
///
/// ### Isotropy subgroups for K4 ($K_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 189 | P-62m | D3h^3 | 6 | 2 |
/// | dir{P1} | 176 | P6_3/m | C6h^2 | 6 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
///
/// ### Isotropy subgroups for K5 ($K_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 189 | P-62m | D3h^3 | 6 | 2 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 9 | 2 |
/// | dir{P1} | 174 | P-6 | C3h^1 | 12 | 2 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 18 | 2 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 18 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 18 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 36 | 2 |
///
/// ### Isotropy subgroups for K6 ($K_{6}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 157 | P31m | C3v^2 | 12 | 2 |
/// | dir{P1} | 150 | P321 | D3^2 | 12 | 2 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 18 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 18 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 18 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 18 | 2 |
/// | dir{P1} | 143 | P3 | C3^1 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 36 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 36 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 36 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 36 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 72 | 2 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{14}t1$ | 6 | E192a | yes |
/// | L2 | $L_{2}$ | $k_{14}t2$ | 6 | E192a | yes |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 44 | Imm2 | C2v^20 | 12 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 24 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 24 | 2 |
/// | dir{P3} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P3} | 187 | P-6m2 | D3h^1 | 16 | 3 |
/// | dir{C1} | 164 | P-3m1 | D3d^3 | 16 | 3 |
/// | dir{P1} | 156 | P3m1 | C3v^1 | 32 | 3 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 48 | 3 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 48 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 96 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 96 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 46 | Ima2 | C2v^22 | 12 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 24 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 24 | 2 |
/// | dir{P3} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P3} | 188 | P-6c2 | D3h^2 | 16 | 3 |
/// | dir{C1} | 165 | P-3c1 | D3d^4 | 16 | 3 |
/// | dir{P1} | 158 | P3c1 | C3v^3 | 32 | 3 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 48 | 3 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 48 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 96 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 96 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
pub struct Sg194;
