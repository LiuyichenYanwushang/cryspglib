//! # Cubic space groups (#195–#230)
//!
//! Cubic crystals have a = b = c, α = β = γ = 90°.  The Brillouin zone
//! is a truncated octahedron for F- and I-centered lattices, and a cube
//! for primitive lattices.
//!
//! ## Common k-point labels (primitive cubic)
//!
//! | Label | Coords (fractional) | Little group |
//! |-------|---------------------|--------------|
//! | Γ | (0, 0, 0) | m-3m (Oₕ) |
//! | X | (0, ½, 0) | 4/mmm (D₄ₕ) |
//! | M | (½, ½, 0) | 4/mmm (D₄ₕ) |
//! | R | (½, ½, ½) | m-3m (Oₕ) |
//!
//! ---


/// # 195 P23 (T¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM2GM3 | $\Gamma_{2}GM3$ | $k_{12}t2t3$ | 2 | B3a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{12}t4$ | 3 | C12a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 195 | P23 | T^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2GM3 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 3 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 4 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{10}t1$ | 3 | C24a | no |
/// | X2 | $X_{3}$ | $k_{10}t2$ | 3 | C24a | no |
/// | X3 | $X_{4}$ | $k_{10}t3$ | 3 | C24a | no |
/// | X4 | $X_{2}$ | $k_{10}t4$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 6 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 12 | 2 |
/// | dir{P1} | 195 | P23 | T^1 | 8 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 6 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 12 | 2 |
/// | dir{P1} | 198 | P2_13 | T^4 | 8 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 3 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 6 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 12 | 2 |
/// | dir{P1} | 198 | P2_13 | T^4 | 8 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 3 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 6 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 12 | 2 |
/// | dir{P1} | 195 | P23 | T^1 | 8 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{11}t1$ | 3 | C12a | no |
/// | M2 | $M_{3}$ | $k_{11}t4$ | 3 | C12a | no |
/// | M3 | $M_{4}$ | $k_{11}t2$ | 3 | C12a | no |
/// | M4 | $M_{2}$ | $k_{11}t3$ | 3 | C12a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 197 | I23 | T^3 | 4 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 12 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 197 | I23 | T^3 | 4 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 12 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 199 | I2_13 | T^5 | 4 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 12 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 199 | I2_13 | T^5 | 4 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 12 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{13}t1$ | 1 | A2a | no |
/// | R2R3 | $R_{2}R3$ | $k_{13}t2t3$ | 2 | B6b | no |
/// | R4 | $R_{4}$ | $k_{13}t4$ | 3 | C24a | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 196 | F23 | T^2 | 2 | 1 |
///
/// ### Isotropy subgroups for R2R3 ($R_{2}R3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 22 | F222 | D2^7 | 6 | 1 |
///
/// ### Isotropy subgroups for R4 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
pub struct Sg195;

/// # 196 F23 (T²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM2GM3 | $\Gamma_{2}GM3$ | $k_{11}t2t3$ | 2 | B3a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{11}t4$ | 3 | C12a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 196 | F23 | T^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2GM3 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 22 | F222 | D2^7 | 3 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{10}t1$ | 3 | C12a | no |
/// | X2 | $X_{3}$ | $k_{10}t2$ | 3 | C12a | no |
/// | X3 | $X_{4}$ | $k_{10}t3$ | 3 | C12a | no |
/// | X4 | $X_{2}$ | $k_{10}t4$ | 3 | C12a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 195 | P23 | T^1 | 4 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 12 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 6 | 1 |
/// | dir{P1} | 198 | P2_13 | T^4 | 4 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 12 | 3 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 6 | 1 |
/// | dir{P1} | 198 | P2_13 | T^4 | 4 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 12 | 3 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 195 | P23 | T^1 | 4 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 12 | 3 |
///
/// ## Irreps at $\mathrm{L}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{9}t1$ | 4 | D96a | no |
/// | L2L3 | $L_{2}L3$ | $k_{9}t2t3$ | 8 | F96a | yes |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{P1} | 196 | F23 | T^2 | 8 | 4 |
/// | dir{P1} | 196 | F23 | T^2 | 8 | 4 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 4 |
///
/// ### Isotropy subgroups for L2L3 ($L_{2}L3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 24 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 24 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 4 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1 | $W_{1}$ | $k_{8}t1$ | 6 | E384d | yes |
/// | W2 | $W_{2}$ | $k_{8}t2$ | 6 | E384d | yes |
///
/// ### Isotropy subgroups for W1 ($W_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 12 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 96 | 4 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 96 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 96 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 96 | 4 |
/// | dir{P3} | 1 | P1 | C1^1 | 192 | 4 |
/// | dir{P1} | 198 | P2_13 | T^4 | 32 | 6 |
/// | dir{C1} | 198 | P2_13 | T^4 | 32 | 6 |
/// | dir{P3} | 195 | P23 | T^1 | 32 | 6 |
/// | dir{P1} | 195 | P23 | T^1 | 32 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 128 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 6 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 192 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 192 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 6 |
///
/// ### Isotropy subgroups for W2 ($W_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 12 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 12 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 96 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 96 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 96 | 4 |
/// | dir{P3} | 3 | P2 | C2^1 | 96 | 4 |
/// | dir{P3} | 1 | P1 | C1^1 | 192 | 4 |
/// | dir{P1} | 198 | P2_13 | T^4 | 32 | 6 |
/// | dir{P3} | 198 | P2_13 | T^4 | 32 | 6 |
/// | dir{C1} | 195 | P23 | T^1 | 32 | 6 |
/// | dir{P1} | 195 | P23 | T^1 | 32 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 128 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 6 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 6 |
/// | dir{C1} | 3 | P2 | C2^1 | 192 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 6 |
///
pub struct Sg196;

/// # 197 I23 (T³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM2GM3 | $\Gamma_{2}GM3$ | $k_{11}t2t3$ | 2 | B3a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{11}t4$ | 3 | C12a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 197 | I23 | T^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2GM3 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 3 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{12}t1$ | 1 | A2a | no |
/// | H2H3 | $H_{2}H3$ | $k_{12}t2t3$ | 2 | B6b | no |
/// | H4 | $H_{4}$ | $k_{12}t4$ | 3 | C24a | yes |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 195 | P23 | T^1 | 2 | 1 |
///
/// ### Isotropy subgroups for H2H3 ($H_{2}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for H4 ($H_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 6 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{N}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{9}t1$ | 6 | E96e | no |
/// | N2 | $N_{2}$ | $k_{9}t2$ | 6 | E96e | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
/// | dir{C1} | 197 | I23 | T^3 | 8 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 6 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 6 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 6 |
///
/// ### Isotropy subgroups for N2 ($N_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
/// | dir{C1} | 199 | I2_13 | T^5 | 8 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 6 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 6 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 6 |
///
/// ## Irreps at $\mathrm{P}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P1 | $P_{1}P1$ | $k_{10}t1t1*$ | 2 | B4a | no |
/// | P2P2 | $P_{2}P2$ | $k_{10}t2t2*$ | 2 | B12b | no |
/// | P3P3 | $P_{3}P3$ | $k_{10}t3t3*$ | 2 | B12b | no |
/// | P4P4 | $P_{4}P4$ | $k_{10}t4t4*$ | 6 | E48a | no |
///
/// ### Isotropy subgroups for P1P1 ($P_{1}P1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 196 | F23 | T^2 | 4 | 1 |
///
/// ### Isotropy subgroups for P2P2 ($P_{2}P2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 1 |
///
/// ### Isotropy subgroups for P3P3 ($P_{3}P3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 1 |
///
/// ### Isotropy subgroups for P4P4 ($P_{4}P4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 1 |
/// | dir{P3} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
pub struct Sg197;

/// # 198 P2₁3 (T⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM2GM3 | $\Gamma_{2}GM3$ | $k_{12}t2t3$ | 2 | B3a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{12}t4$ | 3 | C12a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 198 | P2_13 | T^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2GM3 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 3 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 4 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{10}t1$ | 6 | E96g | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M2 | $M_{1}M2$ | $k_{11}t1t4$ | 6 | E48c | no |
/// | M3M4 | $M_{3}M4$ | $k_{11}t2t3$ | 6 | E48c | no |
///
/// ### Isotropy subgroups for M1M2 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ### Isotropy subgroups for M3M4 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 4 | P2_1 | C2^2 | 12 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R3 | $R_{2}R3$ | $k_{13}t1t2$ | 4 | D24a | yes |
/// | R2R2 | $R_{1}R1$ | $k_{13}t3t3$ | 4 | D24b | yes |
///
/// ### Isotropy subgroups for R1R3 ($R_{2}R3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for R2R2 ($R_{1}R1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
pub struct Sg198;

/// # 199 I2₁3 (T⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM2GM3 | $\Gamma_{2}GM3$ | $k_{11}t2t3$ | 2 | B3a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{11}t4$ | 3 | C12a | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 199 | I2_13 | T^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2GM3 ($\Gamma_{2}GM3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 3 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 4 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 6 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 12 | 1 |
///
/// ## Irreps at $\mathrm{H}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{12}t1$ | 1 | A2a | no |
/// | H2H3 | $H_{2}H3$ | $k_{12}t2t3$ | 2 | B6b | no |
/// | H4 | $H_{4}$ | $k_{12}t4$ | 3 | C24a | yes |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 198 | P2_13 | T^4 | 2 | 1 |
///
/// ### Isotropy subgroups for H2H3 ($H_{2}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 6 | 1 |
///
/// ### Isotropy subgroups for H4 ($H_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 6 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{N}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{9}t1$ | 6 | E96j | no |
/// | N2 | $N_{2}$ | $k_{9}t2$ | 6 | E96j | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 6 |
///
/// ### Isotropy subgroups for N2 ($N_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 6 |
///
/// ## Irreps at $\mathrm{P}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P1 | $P_{2}P2$ | $k_{10}t2t2*$ | 4 | D48b | yes |
/// | P2P2 | $P_{1}P1$ | $k_{10}t3t3*$ | 4 | D48d | yes |
/// | P3P3 | $P_{3}P3$ | $k_{10}t1t1*$ | 4 | D48b | yes |
///
/// ### Isotropy subgroups for P1P1 ($P_{2}P2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 1 |
/// | dir{C1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for P2P2 ($P_{1}P1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for P3P3 ($P_{3}P3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 1 |
/// | dir{C1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
pub struct Sg199;

/// # 200 Pm-3 (Tₕ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{12}t2$ | 1 | A2a | no |
/// | GM2+GM3+ | $\Gamma_{2}^+GM3+$ | $k_{12}t3t5$ | 2 | B3a | no |
/// | GM2-GM3- | $\Gamma_{2}^-GM3-$ | $k_{12}t4t6$ | 2 | B6b | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{12}t7$ | 3 | C12a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{12}t8$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 200 | Pm-3 | Th^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 195 | P23 | T^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+GM3+ ($\Gamma_{2}^+GM3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 47 | Pmmm | D2h^1 | 3 | 1 |
///
/// ### Isotropy subgroups for GM2-GM3- ($\Gamma_{2}^-GM3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 4 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 6 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{10}t1$ | 3 | C24a | no |
/// | X1- | $X_{1}^-$ | $k_{10}t2$ | 3 | C24a | no |
/// | X2+ | $X_{3}^+$ | $k_{10}t3$ | 3 | C24a | no |
/// | X2- | $X_{3}^-$ | $k_{10}t4$ | 3 | C24a | no |
/// | X3+ | $X_{4}^+$ | $k_{10}t5$ | 3 | C24a | no |
/// | X3- | $X_{4}^-$ | $k_{10}t6$ | 3 | C24a | no |
/// | X4+ | $X_{2}^+$ | $k_{10}t7$ | 3 | C24a | no |
/// | X4- | $X_{2}^-$ | $k_{10}t8$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 47 | Pmmm | D2h^1 | 6 | 1 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 12 | 2 |
/// | dir{P1} | 200 | Pm-3 | Th^1 | 8 | 3 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 6 | 1 |
/// | dir{P1} | 50 | Pban | D2h^4 | 12 | 2 |
/// | dir{P1} | 201 | Pn-3 | Th^2 | 8 | 3 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 24 | 3 |
///
/// ### Isotropy subgroups for X2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 6 | 1 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 12 | 2 |
/// | dir{P1} | 205 | Pa-3 | Th^6 | 8 | 3 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 24 | 3 |
///
/// ### Isotropy subgroups for X2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 6 | 1 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 12 | 2 |
/// | dir{P1} | 205 | Pa-3 | Th^6 | 8 | 3 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 24 | 3 |
///
/// ### Isotropy subgroups for X3+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 6 | 1 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 12 | 2 |
/// | dir{P1} | 205 | Pa-3 | Th^6 | 8 | 3 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 24 | 3 |
///
/// ### Isotropy subgroups for X3- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 51 | Pmma | D2h^5 | 6 | 1 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 12 | 2 |
/// | dir{P1} | 205 | Pa-3 | Th^6 | 8 | 3 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 24 | 3 |
///
/// ### Isotropy subgroups for X4+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 49 | Pccm | D2h^3 | 6 | 1 |
/// | dir{P1} | 50 | Pban | D2h^4 | 12 | 2 |
/// | dir{P1} | 201 | Pn-3 | Th^2 | 8 | 3 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 24 | 3 |
///
/// ### Isotropy subgroups for X4- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 47 | Pmmm | D2h^1 | 6 | 1 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 12 | 2 |
/// | dir{P1} | 200 | Pm-3 | Th^1 | 8 | 3 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{11}t1$ | 3 | C12a | no |
/// | M1- | $M_{1}^-$ | $k_{11}t2$ | 3 | C24a | no |
/// | M2+ | $M_{3}^+$ | $k_{11}t7$ | 3 | C12a | no |
/// | M2- | $M_{3}^-$ | $k_{11}t8$ | 3 | C24a | no |
/// | M3+ | $M_{4}^+$ | $k_{11}t3$ | 3 | C12a | no |
/// | M3- | $M_{4}^-$ | $k_{11}t4$ | 3 | C24a | no |
/// | M4+ | $M_{2}^+$ | $k_{11}t5$ | 3 | C12a | no |
/// | M4- | $M_{2}^-$ | $k_{11}t6$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 6 | 1 |
/// | dir{P1} | 204 | Im-3 | Th^5 | 4 | 3 |
/// | dir{P1} | 71 | Immm | D2h^25 | 12 | 3 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 6 | 1 |
/// | dir{P1} | 197 | I23 | T^3 | 8 | 3 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 12 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
///
/// ### Isotropy subgroups for M2+ ($M_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 6 | 1 |
/// | dir{P1} | 204 | Im-3 | Th^5 | 4 | 3 |
/// | dir{P1} | 71 | Immm | D2h^25 | 12 | 3 |
///
/// ### Isotropy subgroups for M2- ($M_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 6 | 1 |
/// | dir{P1} | 197 | I23 | T^3 | 8 | 3 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 12 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
///
/// ### Isotropy subgroups for M3+ ($M_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 6 | 1 |
/// | dir{P1} | 206 | Ia-3 | Th^7 | 4 | 3 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 12 | 3 |
///
/// ### Isotropy subgroups for M3- ($M_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 6 | 1 |
/// | dir{P1} | 199 | I2_13 | T^5 | 8 | 3 |
/// | dir{P1} | 74 | Imma | D2h^28 | 12 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
///
/// ### Isotropy subgroups for M4+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 6 | 1 |
/// | dir{P1} | 206 | Ia-3 | Th^7 | 4 | 3 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 12 | 3 |
///
/// ### Isotropy subgroups for M4- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 6 | 1 |
/// | dir{P1} | 199 | I2_13 | T^5 | 8 | 3 |
/// | dir{P1} | 74 | Imma | D2h^28 | 12 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{13}t1$ | 1 | A2a | no |
/// | R1- | $R_{1}^-$ | $k_{13}t2$ | 1 | A2a | no |
/// | R2+R3+ | $R_{2}^+R3+$ | $k_{13}t3t5$ | 2 | B6b | no |
/// | R2-R3- | $R_{2}^-R3-$ | $k_{13}t4t6$ | 2 | B6b | no |
/// | R4+ | $R_{4}^+$ | $k_{13}t7$ | 3 | C24a | no |
/// | R4- | $R_{4}^-$ | $k_{13}t8$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 202 | Fm-3 | Th^3 | 2 | 1 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 202 | Fm-3 | Th^3 | 2 | 1 |
///
/// ### Isotropy subgroups for R2+R3+ ($R_{2}^+R3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 69 | Fmmm | D2h^23 | 6 | 1 |
///
/// ### Isotropy subgroups for R2-R3- ($R_{2}^-R3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 69 | Fmmm | D2h^23 | 6 | 1 |
///
/// ### Isotropy subgroups for R4+ ($R_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for R4- ($R_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
pub struct Sg200;

/// # 201 Pn-3 (Tₕ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{12}t2$ | 1 | A2a | no |
/// | GM2+GM3+ | $\Gamma_{2}^+GM3+$ | $k_{12}t3t5$ | 2 | B3a | no |
/// | GM2-GM3- | $\Gamma_{2}^-GM3-$ | $k_{12}t4t6$ | 2 | B6b | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{12}t7$ | 3 | C12a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{12}t8$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 201 | Pn-3 | Th^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 195 | P23 | T^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+GM3+ ($\Gamma_{2}^+GM3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 48 | Pnnn | D2h^2 | 3 | 1 |
///
/// ### Isotropy subgroups for GM2-GM3- ($\Gamma_{2}^-GM3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 4 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 6 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{10}t1$ | 6 | E192c | yes |
/// | X2 | $X_{2}$ | $k_{10}t2$ | 6 | E192c | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 16 | P222 | D2^1 | 12 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 12 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 2 |
/// | dir{P3} | 3 | P2 | C2^1 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{C1} | 195 | P23 | T^1 | 16 | 3 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P3} | 3 | P2 | C2^1 | 96 | 3 |
/// | dir{C1} | 3 | P2 | C2^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 12 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 12 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{C1} | 198 | P2_13 | T^4 | 16 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 96 | 3 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{11}t1$ | 6 | E96e | yes |
/// | M2 | $M_{2}$ | $k_{11}t2$ | 6 | E96e | yes |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 12 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 1 |
/// | dir{C1} | 197 | I23 | T^3 | 8 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 16 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 12 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 1 |
/// | dir{C1} | 199 | I2_13 | T^5 | 8 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 16 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{13}t2$ | 1 | A2a | no |
/// | R1- | $R_{1}^-$ | $k_{13}t1$ | 1 | A2a | no |
/// | R2+R3+ | $R_{2}^+R3+$ | $k_{13}t4t6$ | 2 | B6b | no |
/// | R2-R3- | $R_{2}^-R3-$ | $k_{13}t3t5$ | 2 | B6b | no |
/// | R4+ | $R_{4}^+$ | $k_{13}t8$ | 3 | C24a | no |
/// | R4- | $R_{4}^-$ | $k_{13}t7$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 203 | Fd-3 | Th^4 | 2 | 1 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 203 | Fd-3 | Th^4 | 2 | 1 |
///
/// ### Isotropy subgroups for R2+R3+ ($R_{2}^+R3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 70 | Fddd | D2h^24 | 6 | 1 |
///
/// ### Isotropy subgroups for R2-R3- ($R_{2}^-R3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 70 | Fddd | D2h^24 | 6 | 1 |
///
/// ### Isotropy subgroups for R4+ ($R_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for R4- ($R_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
pub struct Sg201;

/// # 202 Fm-3 (Tₕ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM2+GM3+ | $\Gamma_{2}^+GM3+$ | $k_{11}t3t5$ | 2 | B3a | no |
/// | GM2-GM3- | $\Gamma_{2}^-GM3-$ | $k_{11}t4t6$ | 2 | B6b | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{11}t7$ | 3 | C12a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{11}t8$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 202 | Fm-3 | Th^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 196 | F23 | T^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+GM3+ ($\Gamma_{2}^+GM3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 69 | Fmmm | D2h^23 | 3 | 1 |
///
/// ### Isotropy subgroups for GM2-GM3- ($\Gamma_{2}^-GM3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 22 | F222 | D2^7 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 4 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{10}t1$ | 3 | C12a | no |
/// | X1- | $X_{1}^-$ | $k_{10}t2$ | 3 | C24a | no |
/// | X2+ | $X_{3}^+$ | $k_{10}t3$ | 3 | C12a | no |
/// | X2- | $X_{3}^-$ | $k_{10}t4$ | 3 | C24a | no |
/// | X3+ | $X_{4}^+$ | $k_{10}t5$ | 3 | C12a | no |
/// | X3- | $X_{4}^-$ | $k_{10}t6$ | 3 | C24a | no |
/// | X4+ | $X_{2}^+$ | $k_{10}t7$ | 3 | C12a | no |
/// | X4- | $X_{2}^-$ | $k_{10}t8$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 6 | 1 |
/// | dir{P1} | 200 | Pm-3 | Th^1 | 4 | 3 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 12 | 3 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 68 | Ccca | D2h^22 | 6 | 1 |
/// | dir{P1} | 195 | P23 | T^1 | 8 | 3 |
/// | dir{P1} | 50 | Pban | D2h^4 | 12 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 6 | 1 |
/// | dir{P1} | 205 | Pa-3 | Th^6 | 4 | 3 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 12 | 3 |
///
/// ### Isotropy subgroups for X2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 6 | 1 |
/// | dir{P1} | 198 | P2_13 | T^4 | 8 | 3 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 12 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 3 |
///
/// ### Isotropy subgroups for X3+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 6 | 1 |
/// | dir{P1} | 205 | Pa-3 | Th^6 | 4 | 3 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 12 | 3 |
///
/// ### Isotropy subgroups for X3- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 6 | 1 |
/// | dir{P1} | 198 | P2_13 | T^4 | 8 | 3 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 12 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 3 |
///
/// ### Isotropy subgroups for X4+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 66 | Cccm | D2h^20 | 6 | 1 |
/// | dir{P1} | 201 | Pn-3 | Th^2 | 4 | 3 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 12 | 3 |
///
/// ### Isotropy subgroups for X4- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 6 | 1 |
/// | dir{P1} | 195 | P23 | T^1 | 8 | 3 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 12 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ## Irreps at $\mathrm{L}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^+$ | $k_{9}t1$ | 4 | D96a | no |
/// | L1- | $L_{1}^-$ | $k_{9}t2$ | 4 | D96a | no |
/// | L2+L3+ | $L_{2}^+L3+$ | $k_{9}t3t5$ | 8 | F96a | no |
/// | L2-L3- | $L_{2}^-L3-$ | $k_{9}t4t6$ | 8 | F96a | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 203 | Fd-3 | Th^4 | 8 | 4 |
/// | dir{P1} | 202 | Fm-3 | Th^3 | 8 | 4 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 4 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 203 | Fd-3 | Th^4 | 8 | 4 |
/// | dir{P1} | 202 | Fm-3 | Th^3 | 8 | 4 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 4 |
///
/// ### Isotropy subgroups for L2+L3+ ($L_{2}^+L3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 4 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 24 | 4 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 24 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 4 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 48 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 4 |
///
/// ### Isotropy subgroups for L2-L3- ($L_{2}^-L3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 4 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 24 | 4 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 24 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 4 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 48 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 4 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1 | $W_{1}$ | $k_{8}t1$ | 6 | E768a | yes |
/// | W2 | $W_{3}$ | $k_{8}t2$ | 6 | E768a | yes |
/// | W3 | $W_{4}$ | $k_{8}t3$ | 6 | E768a | yes |
/// | W4 | $W_{2}$ | $k_{8}t4$ | 6 | E768a | yes |
///
/// ### Isotropy subgroups for W1 ($W_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 74 | Imma | D2h^28 | 12 | 2 |
/// | dir{P1} | 71 | Immm | D2h^25 | 12 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 24 | 2 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 48 | 4 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 48 | 4 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 48 | 4 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 48 | 4 |
/// | dir{P3} | 31 | Pmn2_1 | C2v^7 | 96 | 4 |
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 96 | 4 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 96 | 4 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 96 | 4 |
/// | dir{P3} | 6 | Pm | Cs^1 | 192 | 4 |
/// | dir{P3} | 200 | Pm-3 | Th^1 | 32 | 6 |
/// | dir{P1} | 200 | Pm-3 | Th^1 | 32 | 6 |
/// | dir{C1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 256 | 6 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 96 | 6 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 96 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{C1} | 11 | P2_1/m | C2h^2 | 192 | 6 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 192 | 6 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 384 | 6 |
/// | dir{C1} | 6 | Pm | Cs^1 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
/// ### Isotropy subgroups for W2 ($W_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 73 | Ibca | D2h^27 | 12 | 2 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 12 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 24 | 2 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 48 | 4 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 48 | 4 |
/// | dir{P1} | 56 | Pccn | D2h^10 | 48 | 4 |
/// | dir{P1} | 50 | Pban | D2h^4 | 48 | 4 |
/// | dir{P3} | 33 | Pna2_1 | C2v^9 | 96 | 4 |
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 96 | 4 |
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 96 | 4 |
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 96 | 4 |
/// | dir{P3} | 7 | Pc | Cs^2 | 192 | 4 |
/// | dir{P3} | 201 | Pn-3 | Th^2 | 32 | 6 |
/// | dir{P1} | 201 | Pn-3 | Th^2 | 32 | 6 |
/// | dir{C1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 256 | 6 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 96 | 6 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 96 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{C1} | 14 | P2_1/c | C2h^5 | 192 | 6 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 192 | 6 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 384 | 6 |
/// | dir{C1} | 7 | Pc | Cs^2 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
/// ### Isotropy subgroups for W3 ($W_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 74 | Imma | D2h^28 | 12 | 2 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 12 | 2 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 24 | 2 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 48 | 4 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 48 | 4 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 48 | 4 |
/// | dir{P1} | 50 | Pban | D2h^4 | 48 | 4 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 96 | 4 |
/// | dir{P3} | 30 | Pnc2 | C2v^6 | 96 | 4 |
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 96 | 4 |
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 96 | 4 |
/// | dir{P3} | 7 | Pc | Cs^2 | 192 | 4 |
/// | dir{P1} | 205 | Pa-3 | Th^6 | 32 | 6 |
/// | dir{P3} | 205 | Pa-3 | Th^6 | 32 | 6 |
/// | dir{C1} | 195 | P23 | T^1 | 64 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 256 | 6 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 96 | 6 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 96 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{C1} | 13 | P2/c | C2h^4 | 192 | 6 |
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 192 | 6 |
/// | dir{C1} | 3 | P2 | C2^1 | 384 | 6 |
/// | dir{C1} | 7 | Pc | Cs^2 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
/// ### Isotropy subgroups for W4 ($W_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 74 | Imma | D2h^28 | 12 | 2 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 12 | 2 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 24 | 2 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 48 | 4 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 48 | 4 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 48 | 4 |
/// | dir{P1} | 50 | Pban | D2h^4 | 48 | 4 |
/// | dir{P3} | 32 | Pba2 | C2v^8 | 96 | 4 |
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 96 | 4 |
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 96 | 4 |
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 96 | 4 |
/// | dir{P3} | 7 | Pc | Cs^2 | 192 | 4 |
/// | dir{P1} | 205 | Pa-3 | Th^6 | 32 | 6 |
/// | dir{P3} | 205 | Pa-3 | Th^6 | 32 | 6 |
/// | dir{C1} | 195 | P23 | T^1 | 64 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 256 | 6 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 96 | 6 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 96 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{C1} | 13 | P2/c | C2h^4 | 192 | 6 |
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 192 | 6 |
/// | dir{C1} | 3 | P2 | C2^1 | 384 | 6 |
/// | dir{C1} | 7 | Pc | Cs^2 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
pub struct Sg202;

/// # 203 Fd-3 (Tₕ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM2+GM3+ | $\Gamma_{2}^+GM3+$ | $k_{11}t3t5$ | 2 | B3a | no |
/// | GM2-GM3- | $\Gamma_{2}^-GM3-$ | $k_{11}t4t6$ | 2 | B6b | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{11}t7$ | 3 | C12a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{11}t8$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 203 | Fd-3 | Th^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 196 | F23 | T^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+GM3+ ($\Gamma_{2}^+GM3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 70 | Fddd | D2h^24 | 3 | 1 |
///
/// ### Isotropy subgroups for GM2-GM3- ($\Gamma_{2}^-GM3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 22 | F222 | D2^7 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{10}t1$ | 6 | E96e | yes |
/// | X2 | $X_{2}$ | $k_{10}t2$ | 6 | E96e | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 12 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 1 |
/// | dir{C1} | 195 | P23 | T^1 | 8 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 16 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{C1} | 3 | P2 | C2^1 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 12 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{C1} | 198 | P2_13 | T^4 | 8 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 16 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $\mathrm{L}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^+$ | $k_{9}t1$ | 4 | D192a | no |
/// | L1- | $L_{1}^-$ | $k_{9}t2$ | 4 | D192a | no |
/// | L2+L3+ | $L_{2}^+L3+$ | $k_{9}t3t5$ | 8 | F192b | no |
/// | L2-L3- | $L_{2}^-L3-$ | $k_{9}t4t6$ | 8 | F192b | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 196 | F23 | T^2 | 16 | 4 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 4 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 196 | F23 | T^2 | 16 | 4 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 4 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L2+L3+ ($L_{2}^+L3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L2-L3- ($L_{2}^-L3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1 | $W_{1}$ | $k_{8}t1$ | 12 | G768c | no |
///
/// ### Isotropy subgroups for W1 ($W_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P3} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 96 | 4 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 96 | 4 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 96 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
/// | dir{P1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{C1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{P1} | 195 | P23 | T^1 | 64 | 6 |
/// | dir{P1} | 195 | P23 | T^1 | 64 | 6 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 128 | 6 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 128 | 6 |
/// | dir{C1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P3} | 18 | P2_12_12 | D2^3 | 192 | 6 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 192 | 6 |
/// | dir{C1} | 17 | P222_1 | D2^2 | 192 | 6 |
/// | dir{C1} | 17 | P222_1 | D2^2 | 192 | 6 |
/// | dir{C1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{C1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 256 | 6 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 384 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 384 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
pub struct Sg203;

/// # 204 Im-3 (Tₕ⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM2+GM3+ | $\Gamma_{2}^+GM3+$ | $k_{11}t3t5$ | 2 | B3a | no |
/// | GM2-GM3- | $\Gamma_{2}^-GM3-$ | $k_{11}t4t6$ | 2 | B6b | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{11}t7$ | 3 | C12a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{11}t8$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 204 | Im-3 | Th^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 197 | I23 | T^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+GM3+ ($\Gamma_{2}^+GM3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 71 | Immm | D2h^25 | 3 | 1 |
///
/// ### Isotropy subgroups for GM2-GM3- ($\Gamma_{2}^-GM3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 23 | I222 | D2^8 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 4 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{H}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1+ | $H_{1}^+$ | $k_{12}t1$ | 1 | A2a | no |
/// | H1- | $H_{1}^-$ | $k_{12}t2$ | 1 | A2a | no |
/// | H2+H3+ | $H_{2}^+H3+$ | $k_{12}t3t5$ | 2 | B6b | no |
/// | H2-H3- | $H_{2}^-H3-$ | $k_{12}t4t6$ | 2 | B6b | no |
/// | H4+ | $H_{4}^+$ | $k_{12}t7$ | 3 | C24a | no |
/// | H4- | $H_{4}^-$ | $k_{12}t8$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for H1+ ($H_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 200 | Pm-3 | Th^1 | 2 | 1 |
///
/// ### Isotropy subgroups for H1- ($H_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 201 | Pn-3 | Th^2 | 2 | 1 |
///
/// ### Isotropy subgroups for H2+H3+ ($H_{2}^+H3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 47 | Pmmm | D2h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for H2-H3- ($H_{2}^-H3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 48 | Pnnn | D2h^2 | 6 | 1 |
///
/// ### Isotropy subgroups for H4+ ($H_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 6 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for H4- ($H_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 6 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{N}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1+ | $N_{1}^+$ | $k_{9}t1$ | 6 | E96e | no |
/// | N1- | $N_{1}^-$ | $k_{9}t2$ | 6 | E192c | no |
/// | N2+ | $N_{2}^+$ | $k_{9}t3$ | 6 | E96e | no |
/// | N2- | $N_{2}^-$ | $k_{9}t4$ | 6 | E192c | no |
///
/// ### Isotropy subgroups for N1+ ($N_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 12 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 24 | 2 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 16 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{C1} | 204 | Im-3 | Th^5 | 8 | 6 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 6 |
/// | dir{P1} | 74 | Imma | D2h^28 | 24 | 6 |
/// | dir{P1} | 71 | Immm | D2h^25 | 24 | 6 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 48 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 6 |
///
/// ### Isotropy subgroups for N1- ($N_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 12 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{C1} | 197 | I23 | T^3 | 16 | 6 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 32 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 6 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 24 | 6 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 24 | 6 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 48 | 6 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 6 |
/// | dir{P3} | 15 | C2/c | C2h^6 | 48 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 6 |
/// | dir{P3} | 9 | Cc | Cs^4 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ### Isotropy subgroups for N2+ ($N_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 12 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 24 | 2 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 16 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{C1} | 206 | Ia-3 | Th^7 | 8 | 6 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 6 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 24 | 6 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 24 | 6 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 48 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 6 |
///
/// ### Isotropy subgroups for N2- ($N_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 12 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{C1} | 199 | I2_13 | T^5 | 16 | 6 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 32 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 6 |
/// | dir{P1} | 74 | Imma | D2h^28 | 24 | 6 |
/// | dir{P1} | 71 | Immm | D2h^25 | 24 | 6 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 48 | 6 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 6 |
/// | dir{P3} | 12 | C2/m | C2h^3 | 48 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 6 |
/// | dir{P3} | 8 | Cm | Cs^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ## Irreps at $\mathrm{P}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1 | $P_{1}$ | $k_{10}t1$ | 2 | B8a | no |
/// | P2P3 | $P_{2}P3$ | $k_{10}t2t3$ | 4 | D24d | no |
/// | P4 | $P_{4}$ | $k_{10}t4$ | 6 | E96a | yes |
///
/// ### Isotropy subgroups for P1 ($P_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 203 | Fd-3 | Th^4 | 4 | 2 |
/// | dir{P1} | 202 | Fm-3 | Th^3 | 4 | 2 |
/// | dir{P1} | 196 | F23 | T^2 | 8 | 2 |
///
/// ### Isotropy subgroups for P2P3 ($P_{2}P3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 70 | Fddd | D2h^24 | 12 | 2 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 12 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 24 | 2 |
///
/// ### Isotropy subgroups for P4 ($P_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1} | 148 | R-3 | C3i^2 | 16 | 2 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 16 | 2 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 12 | 2 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 12 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 24 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 24 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 24 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
///
pub struct Sg204;

/// # 205 Pa-3 (Tₕ⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{12}t2$ | 1 | A2a | no |
/// | GM2+GM3+ | $\Gamma_{2}^+GM3+$ | $k_{12}t3t5$ | 2 | B3a | no |
/// | GM2-GM3- | $\Gamma_{2}^-GM3-$ | $k_{12}t4t6$ | 2 | B6b | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{12}t7$ | 3 | C12a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{12}t8$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 205 | Pa-3 | Th^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 198 | P2_13 | T^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+GM3+ ($\Gamma_{2}^+GM3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 61 | Pbca | D2h^15 | 3 | 1 |
///
/// ### Isotropy subgroups for GM2-GM3- ($\Gamma_{2}^-GM3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 4 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 6 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{2}$ | $k_{10}t2$ | 6 | E192g | yes |
/// | X2 | $X_{1}$ | $k_{10}t1$ | 6 | E192g | yes |
///
/// ### Isotropy subgroups for X1 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 12 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 12 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P3} | 7 | Pc | Cs^2 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 12 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 12 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P3} | 7 | Pc | Cs^2 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M2 | $M_{1}M2$ | $k_{11}t1t2$ | 12 | G96a | no |
///
/// ### Isotropy subgroups for M1M2 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 16 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{4D1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+R3+ | $R_{2}^+R3+$ | $k_{13}t1t3$ | 4 | D24a | no |
/// | R1-R3- | $R_{2}^-R3-$ | $k_{13}t2t4$ | 4 | D24a | no |
/// | R2+R2+ | $R_{1}^+R1+$ | $k_{13}t5t5$ | 4 | D24b | no |
/// | R2-R2- | $R_{1}^-R1-$ | $k_{13}t6t6$ | 4 | D24b | no |
///
/// ### Isotropy subgroups for R1+R3+ ($R_{2}^+R3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for R1-R3- ($R_{2}^-R3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for R2+R2+ ($R_{1}^+R1+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for R2-R2- ($R_{1}^-R1-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
pub struct Sg205;

/// # 206 Ia-3 (Tₕ⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM2+GM3+ | $\Gamma_{2}^+GM3+$ | $k_{11}t3t5$ | 2 | B3a | no |
/// | GM2-GM3- | $\Gamma_{2}^-GM3-$ | $k_{11}t4t6$ | 2 | B6b | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{11}t7$ | 3 | C12a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{11}t8$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 206 | Ia-3 | Th^7 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 199 | I2_13 | T^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+GM3+ ($\Gamma_{2}^+GM3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 73 | Ibca | D2h^27 | 3 | 1 |
///
/// ### Isotropy subgroups for GM2-GM3- ($\Gamma_{2}^-GM3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 4 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 6 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{H}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1+ | $H_{1}^+$ | $k_{12}t1$ | 1 | A2a | no |
/// | H1- | $H_{1}^-$ | $k_{12}t2$ | 1 | A2a | no |
/// | H2+H3+ | $H_{2}^+H3+$ | $k_{12}t3t5$ | 2 | B6b | no |
/// | H2-H3- | $H_{2}^-H3-$ | $k_{12}t4t6$ | 2 | B6b | no |
/// | H4+ | $H_{4}^+$ | $k_{12}t7$ | 3 | C24a | no |
/// | H4- | $H_{4}^-$ | $k_{12}t8$ | 3 | C24a | no |
///
/// ### Isotropy subgroups for H1+ ($H_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 205 | Pa-3 | Th^6 | 2 | 1 |
///
/// ### Isotropy subgroups for H1- ($H_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 205 | Pa-3 | Th^6 | 2 | 1 |
///
/// ### Isotropy subgroups for H2+H3+ ($H_{2}^+H3+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 61 | Pbca | D2h^15 | 6 | 1 |
///
/// ### Isotropy subgroups for H2-H3- ($H_{2}^-H3-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 61 | Pbca | D2h^15 | 6 | 1 |
///
/// ### Isotropy subgroups for H4+ ($H_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 6 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for H4- ($H_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 6 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{N}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{9}t1$ | 12 | G192b | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 1 |
/// | dir{P3} | 37 | Ccc2 | C2v^13 | 24 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 2 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 24 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{C1} | 3 | P2 | C2^1 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 16 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 6 |
/// | dir{C1} | 148 | R-3 | C3i^2 | 32 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 6 |
/// | dir{P3} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{4D1} | 2 | P-1 | Ci^1 | 96 | 6 |
/// | dir{C1} | 2 | P-1 | Ci^1 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ## Irreps at $\mathrm{P}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P3 | $P_{2}P3$ | $k_{10}t1t2$ | 8 | F96c | no |
/// | P2P2 | $P_{1}P1$ | $k_{10}t3t3$ | 8 | F96d | no |
///
/// ### Isotropy subgroups for P1P3 ($P_{2}P3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P3} | 148 | R-3 | C3i^2 | 16 | 2 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 16 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
///
/// ### Isotropy subgroups for P2P2 ($P_{1}P1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
///
pub struct Sg206;

/// # 207 P432 (O¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{12}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{12}t3$ | 2 | B6a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{12}t5$ | 3 | C24c | yes |
/// | GM5 | $\Gamma_{5}$ | $k_{12}t4$ | 3 | C24b | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 207 | P432 | O^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 195 | P23 | T^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 3 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 75 | P4 | C4^1 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 4 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{10}t1$ | 3 | C48a | no |
/// | X2 | $X_{3}$ | $k_{10}t3$ | 3 | C48a | no |
/// | X3 | $X_{2}$ | $k_{10}t2$ | 3 | C48a | no |
/// | X4 | $X_{4}$ | $k_{10}t4$ | 3 | C48a | no |
/// | X5 | $X_{5}$ | $k_{10}t5$ | 6 | E192f | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 6 | 1 |
/// | dir{P1} | 89 | P422 | D4^1 | 12 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 207 | P432 | O^1 | 8 | 3 |
/// | dir{P1} | 89 | P422 | D4^1 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 6 | 1 |
/// | dir{P1} | 89 | P422 | D4^1 | 12 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 208 | P4_232 | O^2 | 8 | 3 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X3 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 6 | 1 |
/// | dir{P1} | 89 | P422 | D4^1 | 12 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 207 | P432 | O^1 | 8 | 3 |
/// | dir{P1} | 89 | P422 | D4^1 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X4 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 6 | 1 |
/// | dir{P1} | 89 | P422 | D4^1 | 12 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 208 | P4_232 | O^2 | 8 | 3 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X5 ($X_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 12 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{P1} | 90 | P42_12 | D4^2 | 12 | 2 |
/// | dir{P3} | 90 | P42_12 | D4^2 | 12 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{C1} | 198 | P2_13 | T^4 | 16 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{C1} | 20 | C222_1 | D2^5 | 48 | 3 |
/// | dir{C1} | 20 | C222_1 | D2^5 | 48 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 96 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 96 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{11}t1$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{11}t3$ | 3 | C24c | no |
/// | M3 | $M_{2}$ | $k_{11}t2$ | 3 | C24c | no |
/// | M4 | $M_{4}$ | $k_{11}t4$ | 3 | C24b | no |
/// | M5 | $M_{5}$ | $k_{11}t5$ | 6 | E96f | yes |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 6 | 1 |
/// | dir{P1} | 211 | I432 | O^5 | 4 | 3 |
/// | dir{P1} | 97 | I422 | D4^9 | 12 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 6 | 1 |
/// | dir{P1} | 197 | I23 | T^3 | 8 | 3 |
/// | dir{P1} | 97 | I422 | D4^9 | 12 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 6 | 1 |
/// | dir{P1} | 197 | I23 | T^3 | 8 | 3 |
/// | dir{P1} | 97 | I422 | D4^9 | 12 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 6 | 1 |
/// | dir{P1} | 211 | I432 | O^5 | 4 | 3 |
/// | dir{P1} | 97 | I422 | D4^9 | 12 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 12 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 1 |
/// | dir{C1} | 199 | I2_13 | T^5 | 8 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 16 | 3 |
/// | dir{P1} | 97 | I422 | D4^9 | 12 | 3 |
/// | dir{P3} | 97 | I422 | D4^9 | 12 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{13}t1$ | 1 | A2a | no |
/// | R2 | $R_{2}$ | $k_{13}t2$ | 1 | A2a | no |
/// | R3 | $R_{3}$ | $k_{13}t3$ | 2 | B12a | no |
/// | R4 | $R_{4}$ | $k_{13}t5$ | 3 | C48a | yes |
/// | R5 | $R_{5}$ | $k_{13}t4$ | 3 | C48a | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 209 | F432 | O^3 | 2 | 1 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 209 | F432 | O^3 | 2 | 1 |
///
/// ### Isotropy subgroups for R3 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 6 | 1 |
/// | dir{P1} | 97 | I422 | D4^9 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 1 |
///
/// ### Isotropy subgroups for R4 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 97 | I422 | D4^9 | 6 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for R5 ($R_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 97 | I422 | D4^9 | 6 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
pub struct Sg207;

/// # 208 P4₂32 (O²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{12}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{12}t3$ | 2 | B6a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{12}t5$ | 3 | C24c | yes |
/// | GM5 | $\Gamma_{5}$ | $k_{12}t4$ | 3 | C24b | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 208 | P4_232 | O^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 195 | P23 | T^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 3 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 77 | P4_2 | C4^3 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 4 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{10}t3$ | 3 | C48a | no |
/// | X2 | $X_{3}$ | $k_{10}t1$ | 3 | C48a | no |
/// | X3 | $X_{2}$ | $k_{10}t4$ | 3 | C48a | no |
/// | X4 | $X_{4}$ | $k_{10}t2$ | 3 | C48a | no |
/// | X5 | $X_{5}$ | $k_{10}t5$ | 6 | E192f | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 6 | 1 |
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 12 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{P1} | 213 | P4_132 | O^7 | 8 | 3 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 24 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 6 | 1 |
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 12 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{P1} | 212 | P4_332 | O^6 | 8 | 3 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 24 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
///
/// ### Isotropy subgroups for X3 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 6 | 1 |
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 12 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{P1} | 213 | P4_132 | O^7 | 8 | 3 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 24 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
///
/// ### Isotropy subgroups for X4 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 6 | 1 |
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 12 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{P1} | 212 | P4_332 | O^6 | 8 | 3 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 24 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
///
/// ### Isotropy subgroups for X5 ($X_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 12 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 1 |
/// | dir{P3} | 93 | P4_222 | D4^5 | 12 | 2 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 12 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P3} | 3 | P2 | C2^1 | 48 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{C1} | 195 | P23 | T^1 | 16 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{C1} | 21 | C222 | D2^6 | 48 | 3 |
/// | dir{C1} | 21 | C222 | D2^6 | 48 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 96 | 3 |
/// | dir{P3} | 3 | P2 | C2^1 | 96 | 3 |
/// | dir{C1} | 3 | P2 | C2^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{11}t1$ | 3 | C24c | no |
/// | M2 | $M_{3}$ | $k_{11}t3$ | 3 | C24b | no |
/// | M3 | $M_{2}$ | $k_{11}t2$ | 3 | C24b | no |
/// | M4 | $M_{4}$ | $k_{11}t4$ | 3 | C24c | no |
/// | M5 | $M_{5}$ | $k_{11}t5$ | 6 | E96f | yes |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 6 | 1 |
/// | dir{P1} | 199 | I2_13 | T^5 | 8 | 3 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 12 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 6 | 1 |
/// | dir{P1} | 214 | I4_132 | O^8 | 4 | 3 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 12 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 6 | 1 |
/// | dir{P1} | 214 | I4_132 | O^8 | 4 | 3 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 12 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 6 | 1 |
/// | dir{P1} | 199 | I2_13 | T^5 | 8 | 3 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 12 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 12 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 1 |
/// | dir{C1} | 197 | I23 | T^3 | 8 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 16 | 3 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 12 | 3 |
/// | dir{P3} | 98 | I4_122 | D4^10 | 12 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{13}t2$ | 1 | A2a | no |
/// | R2 | $R_{2}$ | $k_{13}t1$ | 1 | A2a | no |
/// | R3 | $R_{3}$ | $k_{13}t3$ | 2 | B12a | no |
/// | R4 | $R_{4}$ | $k_{13}t4$ | 3 | C48a | yes |
/// | R5 | $R_{5}$ | $k_{13}t5$ | 3 | C48a | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 210 | F4_132 | O^4 | 2 | 1 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 210 | F4_132 | O^4 | 2 | 1 |
///
/// ### Isotropy subgroups for R3 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 6 | 1 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 1 |
///
/// ### Isotropy subgroups for R4 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 6 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for R5 ($R_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 6 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
pub struct Sg208;

/// # 209 F432 (O³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{11}t3$ | 2 | B6a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{11}t5$ | 3 | C24c | yes |
/// | GM5 | $\Gamma_{5}$ | $k_{11}t4$ | 3 | C24b | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 209 | F432 | O^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 196 | F23 | T^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 3 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 79 | I4 | C4^5 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 4 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{10}t1$ | 3 | C24b | no |
/// | X2 | $X_{3}$ | $k_{10}t3$ | 3 | C24c | no |
/// | X3 | $X_{2}$ | $k_{10}t2$ | 3 | C24c | no |
/// | X4 | $X_{4}$ | $k_{10}t4$ | 3 | C24b | no |
/// | X5 | $X_{5}$ | $k_{10}t5$ | 6 | E96f | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 6 | 1 |
/// | dir{P1} | 207 | P432 | O^1 | 4 | 3 |
/// | dir{P1} | 89 | P422 | D4^1 | 12 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 6 | 1 |
/// | dir{P1} | 195 | P23 | T^1 | 8 | 3 |
/// | dir{P1} | 89 | P422 | D4^1 | 12 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X3 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 6 | 1 |
/// | dir{P1} | 195 | P23 | T^1 | 8 | 3 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 12 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X4 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 6 | 1 |
/// | dir{P1} | 208 | P4_232 | O^2 | 4 | 3 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 12 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X5 ($X_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 12 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{C1} | 198 | P2_13 | T^4 | 8 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 16 | 3 |
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 12 | 3 |
/// | dir{P3} | 90 | P42_12 | D4^2 | 12 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 3 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $\mathrm{L}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{9}t1$ | 4 | D192b | no |
/// | L2 | $L_{2}$ | $k_{9}t2$ | 4 | D192b | no |
/// | L3 | $L_{3}$ | $k_{9}t3$ | 8 | F192a | yes |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 24 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 2 |
/// | dir{P1} | 210 | F4_132 | O^4 | 8 | 4 |
/// | dir{P1} | 209 | F432 | O^3 | 8 | 4 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 4 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 4 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 24 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 2 |
/// | dir{P1} | 210 | F4_132 | O^4 | 8 | 4 |
/// | dir{P1} | 209 | F432 | O^3 | 8 | 4 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 4 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 4 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L3 ($L_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 24 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 4 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 4 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 24 | 4 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 24 | 4 |
/// | dir{P1} | 97 | I422 | D4^9 | 24 | 4 |
/// | dir{P1} | 97 | I422 | D4^9 | 24 | 4 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 4 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 48 | 4 |
/// | dir{P1} | 79 | I4 | C4^5 | 48 | 4 |
/// | dir{P3} | 24 | I2_12_12_1 | D2^9 | 48 | 4 |
/// | dir{C1} | 24 | I2_12_12_1 | D2^9 | 48 | 4 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 4 |
/// | dir{C1} | 23 | I222 | D2^8 | 48 | 4 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 4 |
/// | dir{P3} | 23 | I222 | D2^8 | 48 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P3} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1 | $W_{1}$ | $k_{8}t1$ | 6 | E768c | no |
/// | W2 | $W_{3}$ | $k_{8}t4$ | 6 | E768c | no |
/// | W3W4 | $W_{2}W4$ | $k_{8}t2t3$ | 12 | G768a | no |
///
/// ### Isotropy subgroups for W1 ($W_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 12 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 12 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 24 | 2 |
/// | dir{C1} | 94 | P4_22_12 | D4^6 | 48 | 4 |
/// | dir{P1} | 89 | P422 | D4^1 | 48 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P3} | 1 | P1 | C1^1 | 384 | 4 |
/// | dir{C1} | 213 | P4_132 | O^7 | 32 | 6 |
/// | dir{P1} | 212 | P4_332 | O^6 | 32 | 6 |
/// | dir{P3} | 208 | P4_232 | O^2 | 32 | 6 |
/// | dir{P1} | 207 | P432 | O^1 | 32 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 128 | 6 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 96 | 6 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 96 | 6 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 96 | 6 |
/// | dir{P1} | 89 | P422 | D4^1 | 96 | 6 |
/// | dir{P1} | 21 | C222 | D2^6 | 192 | 6 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 384 | 6 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 384 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
/// ### Isotropy subgroups for W2 ($W_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 12 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 12 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 24 | 2 |
/// | dir{C1} | 94 | P4_22_12 | D4^6 | 48 | 4 |
/// | dir{P1} | 89 | P422 | D4^1 | 48 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P3} | 1 | P1 | C1^1 | 384 | 4 |
/// | dir{P1} | 213 | P4_132 | O^7 | 32 | 6 |
/// | dir{C1} | 212 | P4_332 | O^6 | 32 | 6 |
/// | dir{P1} | 208 | P4_232 | O^2 | 32 | 6 |
/// | dir{P3} | 207 | P432 | O^1 | 32 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 128 | 6 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 96 | 6 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 96 | 6 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 96 | 6 |
/// | dir{P1} | 89 | P422 | D4^1 | 96 | 6 |
/// | dir{P1} | 21 | C222 | D2^6 | 192 | 6 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 384 | 6 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 384 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
/// ### Isotropy subgroups for W3W4 ($W_{2}W4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 2 |
/// | dir{C1} | 22 | F222 | D2^7 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 48 | 4 |
/// | dir{P3} | 90 | P42_12 | D4^2 | 48 | 4 |
/// | dir{C3} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{C1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 4 |
/// | dir{C2} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
/// | dir{C1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{C1} | 195 | P23 | T^1 | 64 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 128 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 128 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 256 | 6 |
/// | dir{C1} | 21 | C222 | D2^6 | 192 | 6 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 192 | 6 |
/// | dir{C1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 384 | 6 |
/// | dir{P1} | 5 | C2 | C2^3 | 384 | 6 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 384 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
pub struct Sg209;

/// # 210 F4₁32 (O⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{11}t3$ | 2 | B6a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{11}t5$ | 3 | C24c | yes |
/// | GM5 | $\Gamma_{5}$ | $k_{11}t4$ | 3 | C24b | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 210 | F4_132 | O^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 196 | F23 | T^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 3 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 4 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{3}$ | $k_{10}t3$ | 3 | C24c | no |
/// | X2 | $X_{1}$ | $k_{10}t1$ | 3 | C24b | no |
/// | X3 | $X_{4}$ | $k_{10}t4$ | 3 | C24b | no |
/// | X4 | $X_{2}$ | $k_{10}t2$ | 3 | C24c | no |
/// | X5 | $X_{5}$ | $k_{10}t5$ | 6 | E96f | yes |
///
/// ### Isotropy subgroups for X1 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 6 | 1 |
/// | dir{P1} | 198 | P2_13 | T^4 | 8 | 3 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 12 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 6 | 1 |
/// | dir{P1} | 212 | P4_332 | O^6 | 4 | 3 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 12 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 3 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 6 | 1 |
/// | dir{P1} | 213 | P4_132 | O^7 | 4 | 3 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 12 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 3 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 6 | 1 |
/// | dir{P1} | 198 | P2_13 | T^4 | 8 | 3 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 12 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 3 |
///
/// ### Isotropy subgroups for X5 ($X_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 12 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 1 |
/// | dir{C1} | 195 | P23 | T^1 | 8 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 16 | 3 |
/// | dir{P1} | 95 | P4_322 | D4^7 | 12 | 3 |
/// | dir{P3} | 91 | P4_122 | D4^3 | 12 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 24 | 3 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{C1} | 3 | P2 | C2^1 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $\mathrm{L}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{9}t1$ | 4 | D192c | no |
/// | L2 | $L_{2}$ | $k_{9}t2$ | 4 | D192c | no |
/// | L3 | $L_{3}$ | $k_{9}t3$ | 8 | F192c | yes |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P1} | 196 | F23 | T^2 | 16 | 4 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 4 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P1} | 196 | F23 | T^2 | 16 | 4 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 4 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L3 ($L_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 4 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 4 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 4 |
/// | dir{P3} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1 | $W_{1}$ | $k_{8}t1$ | 12 | G768d | no |
///
/// ### Isotropy subgroups for W1 ($W_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 48 | 4 |
/// | dir{C1} | 95 | P4_322 | D4^7 | 48 | 4 |
/// | dir{P3} | 92 | P4_12_12 | D4^4 | 48 | 4 |
/// | dir{C1} | 91 | P4_122 | D4^3 | 48 | 4 |
/// | dir{P1} | 78 | P4_3 | C4^4 | 96 | 4 |
/// | dir{P3} | 76 | P4_1 | C4^2 | 96 | 4 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 96 | 4 |
/// | dir{C1} | 20 | C222_1 | D2^5 | 96 | 4 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 4 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P3} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 96 | 4 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 96 | 4 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 96 | 4 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 96 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 4 |
/// | dir{P3} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
/// | dir{C1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{P1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{P1} | 195 | P23 | T^1 | 64 | 6 |
/// | dir{P1} | 195 | P23 | T^1 | 64 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 128 | 6 |
/// | dir{C1} | 155 | R32 | D3^7 | 128 | 6 |
/// | dir{C1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P3} | 18 | P2_12_12 | D2^3 | 192 | 6 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 192 | 6 |
/// | dir{C1} | 17 | P222_1 | D2^2 | 192 | 6 |
/// | dir{C1} | 17 | P222_1 | D2^2 | 192 | 6 |
/// | dir{C1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{C1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 256 | 6 |
/// | dir{P1} | 5 | C2 | C2^3 | 384 | 6 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 384 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
pub struct Sg210;

/// # 211 I432 (O⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{11}t3$ | 2 | B6a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{11}t5$ | 3 | C24c | yes |
/// | GM5 | $\Gamma_{5}$ | $k_{11}t4$ | 3 | C24b | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 211 | I432 | O^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 197 | I23 | T^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 3 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 79 | I4 | C4^5 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 4 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{H}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{12}t1$ | 1 | A2a | no |
/// | H2 | $H_{2}$ | $k_{12}t2$ | 1 | A2a | no |
/// | H3 | $H_{3}$ | $k_{12}t3$ | 2 | B12a | no |
/// | H4 | $H_{4}$ | $k_{12}t5$ | 3 | C48a | yes |
/// | H5 | $H_{5}$ | $k_{12}t4$ | 3 | C48a | yes |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 207 | P432 | O^1 | 2 | 1 |
///
/// ### Isotropy subgroups for H2 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 208 | P4_232 | O^2 | 2 | 1 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 6 | 1 |
/// | dir{P1} | 89 | P422 | D4^1 | 6 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 12 | 1 |
///
/// ### Isotropy subgroups for H4 ($H_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 90 | P42_12 | D4^2 | 6 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for H5 ($H_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 6 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{N}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{9}t1$ | 6 | E192d | no |
/// | N2 | $N_{3}$ | $k_{9}t4$ | 6 | E192e | no |
/// | N3 | $N_{4}$ | $k_{9}t3$ | 6 | E192d | no |
/// | N4 | $N_{2}$ | $k_{9}t2$ | 6 | E192e | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 89 | P422 | D4^1 | 12 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 16 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{C1} | 211 | I432 | O^5 | 8 | 6 |
/// | dir{P3} | 98 | I4_122 | D4^10 | 24 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 6 |
/// | dir{P1} | 97 | I422 | D4^9 | 24 | 6 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 6 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 6 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 6 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ### Isotropy subgroups for N2 ($N_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 1 |
/// | dir{P1} | 90 | P42_12 | D4^2 | 12 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{C1} | 197 | I23 | T^3 | 16 | 6 |
/// | dir{P3} | 155 | R32 | D3^7 | 32 | 6 |
/// | dir{P3} | 98 | I4_122 | D4^10 | 24 | 6 |
/// | dir{C1} | 97 | I422 | D4^9 | 24 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 6 |
/// | dir{P1} | 79 | I4 | C4^5 | 48 | 6 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 6 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ### Isotropy subgroups for N3 ($N_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 12 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 16 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{C1} | 214 | I4_132 | O^8 | 8 | 6 |
/// | dir{P3} | 97 | I422 | D4^9 | 24 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 6 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 24 | 6 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 6 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 6 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 6 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ### Isotropy subgroups for N4 ($N_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 1 |
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 12 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{C1} | 199 | I2_13 | T^5 | 16 | 6 |
/// | dir{P3} | 155 | R32 | D3^7 | 32 | 6 |
/// | dir{C1} | 98 | I4_122 | D4^10 | 24 | 6 |
/// | dir{P3} | 97 | I422 | D4^9 | 24 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 6 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 48 | 6 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 6 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ## Irreps at $\mathrm{P}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1 | $P_{1}$ | $k_{10}t1$ | 2 | B8a | no |
/// | P2 | $P_{2}$ | $k_{10}t2$ | 2 | B24a | no |
/// | P3 | $P_{3}$ | $k_{10}t3$ | 2 | B24a | no |
/// | P4 | $P_{4}$ | $k_{10}t4$ | 6 | E96b | no |
///
/// ### Isotropy subgroups for P1 ($P_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 210 | F4_132 | O^4 | 4 | 2 |
/// | dir{P1} | 209 | F432 | O^3 | 4 | 2 |
/// | dir{P1} | 196 | F23 | T^2 | 8 | 2 |
///
/// ### Isotropy subgroups for P2 ($P_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 12 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 12 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 24 | 2 |
///
/// ### Isotropy subgroups for P3 ($P_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 12 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 12 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 24 | 2 |
///
/// ### Isotropy subgroups for P4 ($P_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P3} | 155 | R32 | D3^7 | 16 | 2 |
/// | dir{C1} | 155 | R32 | D3^7 | 16 | 2 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 12 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 12 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 2 |
/// | dir{C1} | 24 | I2_12_12_1 | D2^9 | 24 | 2 |
/// | dir{C1} | 23 | I222 | D2^8 | 24 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 24 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
///
pub struct Sg211;

/// # 212 P4₃32 (O⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{12}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{12}t3$ | 2 | B6a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{12}t5$ | 3 | C24c | yes |
/// | GM5 | $\Gamma_{5}$ | $k_{12}t4$ | 3 | C24b | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 212 | P4_332 | O^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 198 | P2_13 | T^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 3 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 78 | P4_3 | C4^4 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 4 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{2}$ | $k_{10}t2$ | 6 | E192h | yes |
/// | X2 | $X_{1}$ | $k_{10}t1$ | 6 | E192h | yes |
///
/// ### Isotropy subgroups for X1 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
/// | dir{P1} | 78 | P4_3 | C4^4 | 24 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
/// | dir{P1} | 78 | P4_3 | C4^4 | 24 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M4 | $M_{3}M4$ | $k_{11}t1t4$ | 6 | E96h | no |
/// | M2M3 | $M_{1}M2$ | $k_{11}t2t3$ | 6 | E96k | no |
/// | M5 | $M_{5}$ | $k_{11}t5$ | 6 | E96i | no |
///
/// ### Isotropy subgroups for M1M4 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 12 | 1 |
/// | dir{P1} | 155 | R32 | D3^7 | 16 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ### Isotropy subgroups for M2M3 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 12 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 78 | P4_3 | C4^4 | 12 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 12 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{P3} | 155 | R32 | D3^7 | 16 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{13}t1t2$ | 4 | D48c | yes |
/// | R3 | $R_{3}$ | $k_{13}t3$ | 4 | D48a | yes |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for R3 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
pub struct Sg212;

/// # 213 P4₁32 (O⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{12}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{12}t3$ | 2 | B6a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{12}t5$ | 3 | C24c | yes |
/// | GM5 | $\Gamma_{5}$ | $k_{12}t4$ | 3 | C24b | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 213 | P4_132 | O^7 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 198 | P2_13 | T^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 3 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 76 | P4_1 | C4^2 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 4 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{2}$ | $k_{10}t1$ | 6 | E192h | yes |
/// | X2 | $X_{1}$ | $k_{10}t2$ | 6 | E192h | yes |
///
/// ### Isotropy subgroups for X1 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
/// | dir{P1} | 76 | P4_1 | C4^2 | 24 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
/// | dir{P1} | 76 | P4_1 | C4^2 | 24 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M4 | $M_{1}M2$ | $k_{11}t1t4$ | 6 | E96h | no |
/// | M2M3 | $M_{3}M4$ | $k_{11}t2t3$ | 6 | E96k | no |
/// | M5 | $M_{5}$ | $k_{11}t5$ | 6 | E96i | no |
///
/// ### Isotropy subgroups for M1M4 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 12 | 1 |
/// | dir{P1} | 155 | R32 | D3^7 | 16 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ### Isotropy subgroups for M2M3 ($M_{3}M4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 12 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 76 | P4_1 | C4^2 | 12 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 12 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{P3} | 155 | R32 | D3^7 | 16 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{13}t1t2$ | 4 | D48c | yes |
/// | R3 | $R_{3}$ | $k_{13}t3$ | 4 | D48a | yes |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for R3 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
pub struct Sg213;

/// # 214 I4₁32 (O⁸)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{11}t3$ | 2 | B6a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{11}t5$ | 3 | C24c | yes |
/// | GM5 | $\Gamma_{5}$ | $k_{11}t4$ | 3 | C24b | yes |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 214 | I4_132 | O^8 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 199 | I2_13 | T^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 3 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 4 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 6 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{H}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{2}$ | $k_{12}t2$ | 1 | A2a | no |
/// | H2 | $H_{1}$ | $k_{12}t1$ | 1 | A2a | no |
/// | H3 | $H_{3}$ | $k_{12}t3$ | 2 | B12a | no |
/// | H4 | $H_{4}$ | $k_{12}t4$ | 3 | C48a | yes |
/// | H5 | $H_{5}$ | $k_{12}t5$ | 3 | C48a | yes |
///
/// ### Isotropy subgroups for H1 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 213 | P4_132 | O^7 | 2 | 1 |
///
/// ### Isotropy subgroups for H2 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 212 | P4_332 | O^6 | 2 | 1 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 6 | 1 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 6 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 12 | 1 |
///
/// ### Isotropy subgroups for H4 ($H_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 91 | P4_122 | D4^3 | 6 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for H5 ($H_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 95 | P4_322 | D4^7 | 6 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{N}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{9}t1$ | 6 | E192i | no |
/// | N2 | $N_{4}$ | $k_{9}t4$ | 6 | E192j | no |
/// | N3 | $N_{3}$ | $k_{9}t3$ | 6 | E192i | no |
/// | N4 | $N_{2}$ | $k_{9}t2$ | 6 | E192j | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 91 | P4_122 | D4^3 | 12 | 2 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 24 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 16 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 6 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 6 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ### Isotropy subgroups for N2 ($N_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 1 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 12 | 2 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 32 | 6 |
/// | dir{C1} | 22 | F222 | D2^7 | 48 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 6 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ### Isotropy subgroups for N3 ($N_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 95 | P4_322 | D4^7 | 12 | 2 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 24 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 16 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 6 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 6 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ### Isotropy subgroups for N4 ($N_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 1 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 12 | 2 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 32 | 6 |
/// | dir{C1} | 22 | F222 | D2^7 | 48 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 6 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ## Irreps at $\mathrm{P}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1 | $P_{2}$ | $k_{10}t2$ | 4 | D96b | yes |
/// | P2 | $P_{1}$ | $k_{10}t3$ | 4 | D96c | yes |
/// | P3 | $P_{3}$ | $k_{10}t1$ | 4 | D96b | yes |
///
/// ### Isotropy subgroups for P1 ($P_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 16 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
///
/// ### Isotropy subgroups for P2 ($P_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
///
/// ### Isotropy subgroups for P3 ($P_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 16 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 16 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
///
pub struct Sg214;

/// # 215 P-43m (Tₔ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{12}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{12}t3$ | 2 | B6a | no |
/// | GM4 | $\Gamma_{5}$ | $k_{12}t4$ | 3 | C24b | no |
/// | GM5 | $\Gamma_{4}$ | $k_{12}t5$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 215 | P-43m | Td^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 195 | P23 | T^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 3 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 4 | 1 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{10}t1$ | 3 | C48a | no |
/// | X2 | $X_{3}$ | $k_{10}t4$ | 3 | C48a | no |
/// | X3 | $X_{4}$ | $k_{10}t3$ | 3 | C48a | no |
/// | X4 | $X_{2}$ | $k_{10}t2$ | 3 | C48a | no |
/// | X5 | $X_{5}$ | $k_{10}t5$ | 6 | E192f | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 6 | 1 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 12 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 215 | P-43m | Td^1 | 8 | 3 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 6 | 1 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 12 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 218 | P-43n | Td^4 | 8 | 3 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 6 | 1 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 12 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 215 | P-43m | Td^1 | 8 | 3 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 6 | 1 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 12 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 218 | P-43n | Td^4 | 8 | 3 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X5 ($X_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 12 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 12 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 12 | 2 |
/// | dir{P3} | 113 | P-42_1m | D2d^3 | 12 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{C1} | 8 | Cm | Cs^3 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{C1} | 198 | P2_13 | T^4 | 16 | 3 |
/// | dir{P3} | 161 | R3c | C3v^6 | 32 | 3 |
/// | dir{P1} | 160 | R3m | C3v^5 | 32 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{C1} | 36 | Cmc2_1 | C2v^12 | 48 | 3 |
/// | dir{C1} | 36 | Cmc2_1 | C2v^12 | 48 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
/// | dir{C1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 96 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 96 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{11}t1$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{11}t4$ | 3 | C24c | no |
/// | M3 | $M_{4}$ | $k_{11}t3$ | 3 | C24b | no |
/// | M4 | $M_{2}$ | $k_{11}t2$ | 3 | C24c | no |
/// | M5 | $M_{5}$ | $k_{11}t5$ | 6 | E96f | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 6 | 1 |
/// | dir{P1} | 217 | I-43m | Td^3 | 4 | 3 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 12 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 6 | 1 |
/// | dir{P1} | 197 | I23 | T^3 | 8 | 3 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 12 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 6 | 1 |
/// | dir{P1} | 217 | I-43m | Td^3 | 4 | 3 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 12 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 6 | 1 |
/// | dir{P1} | 197 | I23 | T^3 | 8 | 3 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 12 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 28 | Pma2 | C2v^4 | 12 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 1 |
/// | dir{C1} | 199 | I2_13 | T^5 | 8 | 3 |
/// | dir{P3} | 160 | R3m | C3v^5 | 16 | 3 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 12 | 3 |
/// | dir{P3} | 121 | I-42m | D2d^11 | 12 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
/// | dir{C1} | 8 | Cm | Cs^3 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{13}t1$ | 1 | A2a | no |
/// | R2 | $R_{2}$ | $k_{13}t2$ | 1 | A2a | no |
/// | R3 | $R_{3}$ | $k_{13}t3$ | 2 | B12a | no |
/// | R4 | $R_{5}$ | $k_{13}t4$ | 3 | C48a | no |
/// | R5 | $R_{4}$ | $k_{13}t5$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 216 | F-43m | Td^2 | 2 | 1 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 219 | F-43c | Td^5 | 2 | 1 |
///
/// ### Isotropy subgroups for R3 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 6 | 1 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 1 |
///
/// ### Isotropy subgroups for R4 ($R_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 8 | 1 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 6 | 1 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 12 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for R5 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 8 | 1 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 6 | 1 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
pub struct Sg215;

/// # 216 F-43m (Tₔ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{11}t3$ | 2 | B6a | no |
/// | GM4 | $\Gamma_{5}$ | $k_{11}t4$ | 3 | C24b | no |
/// | GM5 | $\Gamma_{4}$ | $k_{11}t5$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 216 | F-43m | Td^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 196 | F23 | T^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 3 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 4 | 1 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{10}t1$ | 3 | C24b | no |
/// | X2 | $X_{3}$ | $k_{10}t4$ | 3 | C24c | no |
/// | X3 | $X_{4}$ | $k_{10}t3$ | 3 | C24b | no |
/// | X4 | $X_{2}$ | $k_{10}t2$ | 3 | C24c | no |
/// | X5 | $X_{5}$ | $k_{10}t5$ | 6 | E96f | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 6 | 1 |
/// | dir{P1} | 215 | P-43m | Td^1 | 4 | 3 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 12 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 6 | 1 |
/// | dir{P1} | 195 | P23 | T^1 | 8 | 3 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 12 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 6 | 1 |
/// | dir{P1} | 215 | P-43m | Td^1 | 4 | 3 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 12 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 6 | 1 |
/// | dir{P1} | 195 | P23 | T^1 | 8 | 3 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 12 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X5 ($X_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 12 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{C1} | 198 | P2_13 | T^4 | 8 | 3 |
/// | dir{P3} | 160 | R3m | C3v^5 | 16 | 3 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 12 | 3 |
/// | dir{P3} | 113 | P-42_1m | D2d^3 | 12 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 3 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 3 |
/// | dir{C1} | 8 | Cm | Cs^3 | 48 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $\mathrm{L}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{9}t1$ | 4 | D192b | no |
/// | L2 | $L_{2}$ | $k_{9}t2$ | 4 | D192b | no |
/// | L3 | $L_{3}$ | $k_{9}t3$ | 8 | F192a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 8 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 24 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 2 |
/// | dir{P1} | 216 | F-43m | Td^2 | 8 | 4 |
/// | dir{P1} | 216 | F-43m | Td^2 | 8 | 4 |
/// | dir{P1} | 160 | R3m | C3v^5 | 32 | 4 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 48 | 4 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 48 | 4 |
/// | dir{C1} | 8 | Cm | Cs^3 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L2 ($L_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 8 | 1 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 24 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 48 | 2 |
/// | dir{P1} | 219 | F-43c | Td^5 | 8 | 4 |
/// | dir{P1} | 219 | F-43c | Td^5 | 8 | 4 |
/// | dir{P1} | 161 | R3c | C3v^6 | 32 | 4 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 48 | 4 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 48 | 4 |
/// | dir{C1} | 9 | Cc | Cs^4 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L3 ($L_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 24 | 2 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 24 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 48 | 2 |
/// | dir{P1} | 7 | Pc | Cs^2 | 48 | 2 |
/// | dir{P1} | 6 | Pm | Cs^1 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P1} | 161 | R3c | C3v^6 | 32 | 4 |
/// | dir{P1} | 160 | R3m | C3v^5 | 32 | 4 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 24 | 4 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 24 | 4 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 24 | 4 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 24 | 4 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 4 |
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 4 |
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 4 |
/// | dir{C1} | 46 | Ima2 | C2v^22 | 48 | 4 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 48 | 4 |
/// | dir{P3} | 45 | Iba2 | C2v^21 | 48 | 4 |
/// | dir{C1} | 45 | Iba2 | C2v^21 | 48 | 4 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 48 | 4 |
/// | dir{P3} | 44 | Imm2 | C2v^20 | 48 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 4 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 4 |
/// | dir{P3} | 8 | Cm | Cs^3 | 96 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1 | $W_{1}$ | $k_{8}t1$ | 6 | E768b | no |
/// | W2 | $W_{3}$ | $k_{8}t3$ | 6 | E768b | no |
/// | W3 | $W_{2}$ | $k_{8}t2$ | 6 | E768b | no |
/// | W4 | $W_{4}$ | $k_{8}t4$ | 6 | E768b | no |
///
/// ### Isotropy subgroups for W1 ($W_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 12 | 2 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 12 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 24 | 2 |
/// | dir{C1} | 113 | P-42_1m | D2d^3 | 48 | 4 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 48 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 4 |
/// | dir{P1} | 8 | Cm | Cs^3 | 192 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P3} | 1 | P1 | C1^1 | 384 | 4 |
/// | dir{P3} | 218 | P-43n | Td^4 | 32 | 6 |
/// | dir{P1} | 215 | P-43m | Td^1 | 32 | 6 |
/// | dir{C1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 256 | 6 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 96 | 6 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 96 | 6 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 192 | 6 |
/// | dir{P1} | 81 | P-4 | S4^1 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{P3} | 9 | Cc | Cs^4 | 384 | 6 |
/// | dir{P1} | 8 | Cm | Cs^3 | 384 | 6 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 384 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
/// ### Isotropy subgroups for W2 ($W_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 12 | 2 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 12 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 24 | 2 |
/// | dir{C1} | 113 | P-42_1m | D2d^3 | 48 | 4 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 48 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 4 |
/// | dir{P1} | 8 | Cm | Cs^3 | 192 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P3} | 1 | P1 | C1^1 | 384 | 4 |
/// | dir{P1} | 218 | P-43n | Td^4 | 32 | 6 |
/// | dir{P3} | 215 | P-43m | Td^1 | 32 | 6 |
/// | dir{C1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 256 | 6 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 96 | 6 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 96 | 6 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 192 | 6 |
/// | dir{P1} | 81 | P-4 | S4^1 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{P1} | 9 | Cc | Cs^4 | 384 | 6 |
/// | dir{P3} | 8 | Cm | Cs^3 | 384 | 6 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 384 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
/// ### Isotropy subgroups for W3 ($W_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 12 | 2 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 12 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 24 | 2 |
/// | dir{C1} | 113 | P-42_1m | D2d^3 | 48 | 4 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 48 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 4 |
/// | dir{P1} | 8 | Cm | Cs^3 | 192 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P3} | 1 | P1 | C1^1 | 384 | 4 |
/// | dir{P3} | 218 | P-43n | Td^4 | 32 | 6 |
/// | dir{P1} | 215 | P-43m | Td^1 | 32 | 6 |
/// | dir{C1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 256 | 6 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 96 | 6 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 96 | 6 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 192 | 6 |
/// | dir{P1} | 81 | P-4 | S4^1 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{P3} | 9 | Cc | Cs^4 | 384 | 6 |
/// | dir{P1} | 8 | Cm | Cs^3 | 384 | 6 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 384 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
/// ### Isotropy subgroups for W4 ($W_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 12 | 2 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 12 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 24 | 2 |
/// | dir{C1} | 113 | P-42_1m | D2d^3 | 48 | 4 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 48 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 4 |
/// | dir{P1} | 8 | Cm | Cs^3 | 192 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P3} | 1 | P1 | C1^1 | 384 | 4 |
/// | dir{P1} | 218 | P-43n | Td^4 | 32 | 6 |
/// | dir{P3} | 215 | P-43m | Td^1 | 32 | 6 |
/// | dir{C1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 256 | 6 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 96 | 6 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 96 | 6 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 192 | 6 |
/// | dir{P1} | 81 | P-4 | S4^1 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{P1} | 9 | Cc | Cs^4 | 384 | 6 |
/// | dir{P3} | 8 | Cm | Cs^3 | 384 | 6 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 384 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
pub struct Sg216;

/// # 217 I-43m (Tₔ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{11}t3$ | 2 | B6a | no |
/// | GM4 | $\Gamma_{5}$ | $k_{11}t4$ | 3 | C24b | no |
/// | GM5 | $\Gamma_{4}$ | $k_{11}t5$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 217 | I-43m | Td^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 197 | I23 | T^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 3 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 4 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 6 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{H}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{12}t1$ | 1 | A2a | no |
/// | H2 | $H_{2}$ | $k_{12}t2$ | 1 | A2a | no |
/// | H3 | $H_{3}$ | $k_{12}t3$ | 2 | B12a | no |
/// | H4 | $H_{5}$ | $k_{12}t4$ | 3 | C48a | no |
/// | H5 | $H_{4}$ | $k_{12}t5$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 215 | P-43m | Td^1 | 2 | 1 |
///
/// ### Isotropy subgroups for H2 ($H_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 218 | P-43n | Td^4 | 2 | 1 |
///
/// ### Isotropy subgroups for H3 ($H_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 6 | 1 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 6 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 12 | 1 |
///
/// ### Isotropy subgroups for H4 ($H_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 8 | 1 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 6 | 1 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 12 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for H5 ($H_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 8 | 1 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 6 | 1 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{N}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{9}t1$ | 6 | E192d | no |
/// | N2 | $N_{3}$ | $k_{9}t2$ | 6 | E192e | no |
/// | N3 | $N_{2}$ | $k_{9}t3$ | 6 | E192e | no |
/// | N4 | $N_{4}$ | $k_{9}t4$ | 6 | E192d | no |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 38 | Amm2 | C2v^14 | 12 | 1 |
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 12 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 24 | 2 |
/// | dir{P1} | 160 | R3m | C3v^5 | 16 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{C1} | 217 | I-43m | Td^3 | 8 | 6 |
/// | dir{P3} | 122 | I-42d | D2d^12 | 24 | 6 |
/// | dir{P1} | 160 | R3m | C3v^5 | 32 | 6 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 24 | 6 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 6 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 48 | 6 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 6 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 6 |
/// | dir{C1} | 8 | Cm | Cs^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ### Isotropy subgroups for N2 ($N_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 41 | Aba2 | C2v^17 | 12 | 1 |
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 12 | 2 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{C1} | 197 | I23 | T^3 | 16 | 6 |
/// | dir{P3} | 161 | R3c | C3v^6 | 32 | 6 |
/// | dir{P3} | 122 | I-42d | D2d^12 | 24 | 6 |
/// | dir{C1} | 121 | I-42m | D2d^11 | 24 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 6 |
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 6 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 6 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 6 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 6 |
/// | dir{C1} | 9 | Cc | Cs^4 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ### Isotropy subgroups for N3 ($N_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 39 | Abm2 | C2v^15 | 12 | 1 |
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 12 | 2 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 24 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{C1} | 199 | I2_13 | T^5 | 16 | 6 |
/// | dir{P3} | 160 | R3m | C3v^5 | 32 | 6 |
/// | dir{C1} | 122 | I-42d | D2d^12 | 24 | 6 |
/// | dir{P3} | 121 | I-42m | D2d^11 | 24 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 6 |
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 6 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 6 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 6 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 6 |
/// | dir{C1} | 8 | Cm | Cs^3 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ### Isotropy subgroups for N4 ($N_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 40 | Ama2 | C2v^16 | 12 | 1 |
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 12 | 2 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 24 | 2 |
/// | dir{P1} | 160 | R3m | C3v^5 | 16 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{C1} | 220 | I-43d | Td^6 | 8 | 6 |
/// | dir{P3} | 121 | I-42m | D2d^11 | 24 | 6 |
/// | dir{P1} | 161 | R3c | C3v^6 | 32 | 6 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 24 | 6 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 6 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 48 | 6 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 6 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 6 |
/// | dir{C1} | 9 | Cc | Cs^4 | 96 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ## Irreps at $\mathrm{P}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P1 | $P_{1}P1$ | $k_{10}t1t1*$ | 2 | B4a | no |
/// | P2P2 | $P_{2}P2$ | $k_{10}t2t2*$ | 2 | B4a | no |
/// | P3P3 | $P_{3}P3$ | $k_{10}t3t3*$ | 4 | D24c | no |
/// | P4P4 | $P_{5}P5$ | $k_{10}t4t4*$ | 6 | E96c | yes |
/// | P5P5 | $P_{4}P4$ | $k_{10}t5t5*$ | 6 | E96c | yes |
///
/// ### Isotropy subgroups for P1P1 ($P_{1}P1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 216 | F-43m | Td^2 | 4 | 1 |
///
/// ### Isotropy subgroups for P2P2 ($P_{2}P2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 219 | F-43c | Td^5 | 4 | 1 |
///
/// ### Isotropy subgroups for P3P3 ($P_{3}P3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 12 | 1 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 12 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 24 | 1 |
///
/// ### Isotropy subgroups for P4P4 ($P_{5}P5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 16 | 1 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 12 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 24 | 1 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 24 | 1 |
/// | dir{C1} | 8 | Cm | Cs^3 | 48 | 1 |
/// | dir{P3} | 5 | C2 | C2^3 | 48 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 1 |
///
/// ### Isotropy subgroups for P5P5 ($P_{4}P4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 16 | 1 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 12 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 24 | 1 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 24 | 1 |
/// | dir{C1} | 9 | Cc | Cs^4 | 48 | 1 |
/// | dir{P3} | 5 | C2 | C2^3 | 48 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 1 |
///
pub struct Sg217;

/// # 218 P-43n (Tₔ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{12}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{12}t3$ | 2 | B6a | no |
/// | GM4 | $\Gamma_{5}$ | $k_{12}t4$ | 3 | C24b | no |
/// | GM5 | $\Gamma_{4}$ | $k_{12}t5$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 218 | P-43n | Td^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 195 | P23 | T^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 3 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 4 | 1 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1X2 | $X_{1}X2$ | $k_{10}t1t4$ | 6 | E48b | no |
/// | X3X4 | $X_{3}X4$ | $k_{10}t2t3$ | 6 | E48b | no |
/// | X5 | $X_{5}$ | $k_{10}t5$ | 6 | E192b | yes |
///
/// ### Isotropy subgroups for X1X2 ($X_{1}X2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 12 | 1 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 12 | 2 |
/// | dir{P3} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{P1} | 198 | P2_13 | T^4 | 16 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
///
/// ### Isotropy subgroups for X3X4 ($X_{3}X4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 17 | P222_1 | D2^2 | 12 | 1 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 12 | 2 |
/// | dir{P3} | 18 | P2_12_12 | D2^3 | 24 | 2 |
/// | dir{P1} | 198 | P2_13 | T^4 | 16 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
///
/// ### Isotropy subgroups for X5 ($X_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 81 | P-4 | S4^1 | 12 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 12 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 1 |
/// | dir{P3} | 112 | P-42c | D2d^2 | 12 | 2 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 12 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 2 |
/// | dir{C1} | 9 | Cc | Cs^4 | 48 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 48 | 2 |
/// | dir{P3} | 3 | P2 | C2^1 | 48 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{C1} | 195 | P23 | T^1 | 16 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{C1} | 81 | P-4 | S4^1 | 48 | 3 |
/// | dir{C1} | 81 | P-4 | S4^1 | 48 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
/// | dir{C1} | 3 | P2 | C2^1 | 96 | 3 |
/// | dir{P1} | 3 | P2 | C2^1 | 96 | 3 |
/// | dir{P3} | 3 | P2 | C2^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{11}t4$ | 3 | C24b | no |
/// | M2 | $M_{3}$ | $k_{11}t1$ | 3 | C24c | no |
/// | M3 | $M_{4}$ | $k_{11}t2$ | 3 | C24b | no |
/// | M4 | $M_{2}$ | $k_{11}t3$ | 3 | C24c | no |
/// | M5 | $M_{5}$ | $k_{11}t5$ | 6 | E96f | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 6 | 1 |
/// | dir{P1} | 220 | I-43d | Td^6 | 4 | 3 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 12 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 6 | 1 |
/// | dir{P1} | 199 | I2_13 | T^5 | 8 | 3 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 12 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 6 | 1 |
/// | dir{P1} | 220 | I-43d | Td^6 | 4 | 3 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 12 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 6 | 1 |
/// | dir{P1} | 199 | I2_13 | T^5 | 8 | 3 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 12 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
///
/// ### Isotropy subgroups for M5 ($M_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 12 | 1 |
/// | dir{P1} | 21 | C222 | D2^6 | 12 | 1 |
/// | dir{P1} | 3 | P2 | C2^1 | 24 | 1 |
/// | dir{C1} | 197 | I23 | T^3 | 8 | 3 |
/// | dir{P3} | 161 | R3c | C3v^6 | 16 | 3 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 12 | 3 |
/// | dir{P3} | 122 | I-42d | D2d^12 | 12 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 3 |
/// | dir{C1} | 9 | Cc | Cs^4 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{13}t1t2$ | 2 | B4a | no |
/// | R3R3 | $R_{3}R3$ | $k_{13}t3t3$ | 4 | D12a | no |
/// | R4R5 | $R_{4}R5$ | $k_{13}t4t5$ | 6 | E48b | no |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 196 | F23 | T^2 | 4 | 1 |
///
/// ### Isotropy subgroups for R3R3 ($R_{3}R3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 1 |
///
/// ### Isotropy subgroups for R4R5 ($R_{4}R5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 12 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 1 |
/// | dir{P3} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
pub struct Sg218;

/// # 219 F-43c (Tₔ⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{11}t3$ | 2 | B6a | no |
/// | GM4 | $\Gamma_{5}$ | $k_{11}t4$ | 3 | C24b | no |
/// | GM5 | $\Gamma_{4}$ | $k_{11}t5$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 219 | F-43c | Td^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 196 | F23 | T^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 3 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 4 | 1 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{10}t4$ | 3 | C24b | no |
/// | X2 | $X_{3}$ | $k_{10}t1$ | 3 | C24c | no |
/// | X3 | $X_{4}$ | $k_{10}t2$ | 3 | C24b | no |
/// | X4 | $X_{2}$ | $k_{10}t3$ | 3 | C24c | no |
/// | X5 | $X_{5}$ | $k_{10}t5$ | 6 | E96f | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 6 | 1 |
/// | dir{P1} | 218 | P-43n | Td^4 | 4 | 3 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 12 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 6 | 1 |
/// | dir{P1} | 195 | P23 | T^1 | 8 | 3 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 12 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 6 | 1 |
/// | dir{P1} | 218 | P-43n | Td^4 | 4 | 3 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 12 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 6 | 1 |
/// | dir{P1} | 195 | P23 | T^1 | 8 | 3 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 12 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X5 ($X_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 12 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 12 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 24 | 1 |
/// | dir{C1} | 198 | P2_13 | T^4 | 8 | 3 |
/// | dir{P3} | 161 | R3c | C3v^6 | 16 | 3 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 12 | 3 |
/// | dir{P3} | 114 | P-42_1c | D2d^4 | 12 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 3 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 3 |
/// | dir{C1} | 9 | Cc | Cs^4 | 48 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
///
/// ## Irreps at $\mathrm{L}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1L2 | $L_{1}L2$ | $k_{9}t1t2$ | 8 | F192d | yes |
/// | L3L3 | $L_{3}L3$ | $k_{9}t3t3$ | 16 | H192a | no |
///
/// ### Isotropy subgroups for L1L2 ($L_{1}L2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P1} | 196 | F23 | T^2 | 16 | 4 |
/// | dir{P1} | 196 | F23 | T^2 | 16 | 4 |
/// | dir{P3} | 82 | I-4 | S4^2 | 48 | 4 |
/// | dir{C1} | 82 | I-4 | S4^2 | 48 | 4 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L3L3 ($L_{3}L3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 4 |
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 4 |
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 4 |
/// | dir{C1} | 22 | F222 | D2^7 | 48 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 4 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1W2 | $W_{1}W3$ | $k_{8}t1t3$ | 12 | G768b | no |
/// | W3W4 | $W_{2}W4$ | $k_{8}t2t4$ | 12 | G768b | no |
///
/// ### Isotropy subgroups for W1W2 ($W_{1}W3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1} | 82 | I-4 | S4^2 | 24 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 48 | 4 |
/// | dir{P3} | 112 | P-42c | D2d^2 | 48 | 4 |
/// | dir{C3} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{C1} | 16 | P222 | D2^1 | 96 | 4 |
/// | dir{C2} | 9 | Cc | Cs^4 | 192 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{C1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
/// | dir{P1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{C1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{C1} | 195 | P23 | T^1 | 64 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 256 | 6 |
/// | dir{P1} | 81 | P-4 | S4^1 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{C1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 384 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
/// ### Isotropy subgroups for W3W4 ($W_{2}W4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1} | 82 | I-4 | S4^2 | 24 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 24 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 48 | 4 |
/// | dir{P3} | 112 | P-42c | D2d^2 | 48 | 4 |
/// | dir{C3} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 96 | 4 |
/// | dir{C1} | 16 | P222 | D2^1 | 96 | 4 |
/// | dir{C2} | 9 | Cc | Cs^4 | 192 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 4 |
/// | dir{C1} | 3 | P2 | C2^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
/// | dir{C1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{P1} | 198 | P2_13 | T^4 | 64 | 6 |
/// | dir{C1} | 195 | P23 | T^1 | 64 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 256 | 6 |
/// | dir{P1} | 81 | P-4 | S4^1 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 6 |
/// | dir{C1} | 16 | P222 | D2^1 | 192 | 6 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 384 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 6 |
///
pub struct Sg219;

/// # 220 I-43d (Tₔ⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{3}$ | $k_{11}t3$ | 2 | B6a | no |
/// | GM4 | $\Gamma_{5}$ | $k_{11}t4$ | 3 | C24b | no |
/// | GM5 | $\Gamma_{4}$ | $k_{11}t5$ | 3 | C24c | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 220 | I-43d | Td^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 199 | I2_13 | T^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 3 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 6 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 4 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 8 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 6 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 12 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 24 | 1 |
///
/// ## Irreps at $\mathrm{H}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1H2 | $H_{1}H2$ | $k_{12}t1t2$ | 2 | B4a | no |
/// | H3H3 | $H_{3}H3$ | $k_{12}t3t3$ | 4 | D12a | no |
/// | H4H5 | $H_{4}H5$ | $k_{12}t4t5$ | 6 | E48b | no |
///
/// ### Isotropy subgroups for H1H2 ($H_{1}H2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 198 | P2_13 | T^4 | 4 | 1 |
///
/// ### Isotropy subgroups for H3H3 ($H_{3}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 12 | 1 |
///
/// ### Isotropy subgroups for H4H5 ($H_{4}H5$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 146 | R3 | C3^4 | 16 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 12 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 12 | 1 |
/// | dir{P3} | 3 | P2 | C2^1 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{N}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{1}$ | $k_{9}t1$ | 12 | G192a | yes |
///
/// ### Isotropy subgroups for N1 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P3} | 1 | P1 | C1^1 | 48 | 1 |
/// | dir{P3} | 81 | P-4 | S4^1 | 24 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 2 |
/// | dir{C1} | 3 | P2 | C2^1 | 48 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 2 |
/// | dir{P1} | 161 | R3c | C3v^6 | 16 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 3 |
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 6 |
/// | dir{P3} | 5 | C2 | C2^3 | 96 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 6 |
///
/// ## Irreps at $\mathrm{P}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P1 | $P_{1}P1$ | $k_{10}t2t2*$ | 4 | D96d | no |
/// | P2P2 | $P_{2}P2$ | $k_{10}t1t1*$ | 4 | D96d | no |
/// | P3P3 | $P_{3}P3$ | $k_{10}t3t3*$ | 8 | F96b | no |
///
/// ### Isotropy subgroups for P1P1 ($P_{1}P1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 1 |
///
/// ### Isotropy subgroups for P2P2 ($P_{2}P2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 1 |
///
/// ### Isotropy subgroups for P3P3 ($P_{3}P3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 24 | 1 |
/// | dir{P3} | 82 | I-4 | S4^2 | 24 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 1 |
///
pub struct Sg220;

/// # 221 Pm-3m (Oₕ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{12}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{12}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{12}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{12}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{12}t6$ | 2 | B12a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{12}t9$ | 3 | C24c | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{12}t10$ | 3 | C48a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{12}t7$ | 3 | C24b | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{12}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 221 | Pm-3m | Oh^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 207 | P432 | O^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 200 | Pm-3 | Th^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 215 | P-43m | Td^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 3 | 1 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 6 | 1 |
/// | dir{P1} | 89 | P422 | D4^1 | 6 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 83 | P4/m | C4h^1 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 8 | 1 |
/// | dir{P1} | 99 | P4mm | C4v^1 | 6 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 12 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 4 | 1 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 6 | 1 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 12 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{10}t1$ | 3 | C48a | no |
/// | X1- | $X_{1}^-$ | $k_{10}t2$ | 3 | C48a | no |
/// | X2+ | $X_{3}^+$ | $k_{10}t5$ | 3 | C48a | no |
/// | X2- | $X_{3}^-$ | $k_{10}t6$ | 3 | C48a | no |
/// | X3+ | $X_{2}^+$ | $k_{10}t3$ | 3 | C48a | no |
/// | X3- | $X_{2}^-$ | $k_{10}t4$ | 3 | C48a | no |
/// | X4+ | $X_{4}^+$ | $k_{10}t7$ | 3 | C48a | no |
/// | X4- | $X_{4}^-$ | $k_{10}t8$ | 3 | C48a | no |
/// | X5+ | $X_{5}^+$ | $k_{10}t9$ | 6 | E192f | no |
/// | X5- | $X_{5}^-$ | $k_{10}t10$ | 6 | E192f | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 6 | 1 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 12 | 2 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 24 | 2 |
/// | dir{P1} | 221 | Pm-3m | Oh^1 | 8 | 3 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 24 | 3 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 124 | P4/mcc | D4h^2 | 6 | 1 |
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 12 | 2 |
/// | dir{P1} | 50 | Pban | D2h^4 | 24 | 2 |
/// | dir{P1} | 222 | Pn-3n | Oh^2 | 8 | 3 |
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 24 | 3 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 48 | 3 |
///
/// ### Isotropy subgroups for X2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 6 | 1 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 12 | 2 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 24 | 2 |
/// | dir{P1} | 223 | Pm-3n | Oh^3 | 8 | 3 |
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 24 | 3 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 6 | 1 |
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 12 | 2 |
/// | dir{P1} | 50 | Pban | D2h^4 | 24 | 2 |
/// | dir{P1} | 224 | Pn-3m | Oh^4 | 8 | 3 |
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 24 | 3 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 48 | 3 |
///
/// ### Isotropy subgroups for X3+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 124 | P4/mcc | D4h^2 | 6 | 1 |
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 12 | 2 |
/// | dir{P1} | 50 | Pban | D2h^4 | 24 | 2 |
/// | dir{P1} | 222 | Pn-3n | Oh^2 | 8 | 3 |
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 24 | 3 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 48 | 3 |
///
/// ### Isotropy subgroups for X3- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 6 | 1 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 12 | 2 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 24 | 2 |
/// | dir{P1} | 221 | Pm-3m | Oh^1 | 8 | 3 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 24 | 3 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X4+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 6 | 1 |
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 12 | 2 |
/// | dir{P1} | 50 | Pban | D2h^4 | 24 | 2 |
/// | dir{P1} | 224 | Pn-3m | Oh^4 | 8 | 3 |
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 24 | 3 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 48 | 3 |
///
/// ### Isotropy subgroups for X4- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 6 | 1 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 12 | 2 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 24 | 2 |
/// | dir{P1} | 223 | Pm-3n | Oh^3 | 8 | 3 |
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 24 | 3 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X5+ ($X_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 12 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 1 |
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 12 | 2 |
/// | dir{P3} | 127 | P4/mbm | D4h^5 | 12 | 2 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 24 | 2 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 24 | 2 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 24 | 2 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 48 | 2 |
/// | dir{P3} | 11 | P2_1/m | C2h^2 | 48 | 2 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{C1} | 205 | Pa-3 | Th^6 | 16 | 3 |
/// | dir{P3} | 167 | R-3c | D3d^6 | 32 | 3 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 32 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 64 | 3 |
/// | dir{C1} | 64 | Cmca | D2h^18 | 48 | 3 |
/// | dir{C1} | 63 | Cmcm | D2h^17 | 48 | 3 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 48 | 3 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 48 | 3 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 96 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 3 |
/// | dir{P3} | 14 | P2_1/c | C2h^5 | 96 | 3 |
/// | dir{C1} | 14 | P2_1/c | C2h^5 | 96 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 3 |
///
/// ### Isotropy subgroups for X5- ($X_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 12 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 1 |
/// | dir{P3} | 129 | P4/nmm | D4h^7 | 12 | 2 |
/// | dir{P1} | 127 | P4/mbm | D4h^5 | 12 | 2 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 24 | 2 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 24 | 2 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 24 | 2 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P3} | 14 | P2_1/c | C2h^5 | 48 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 48 | 2 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{C1} | 205 | Pa-3 | Th^6 | 16 | 3 |
/// | dir{P3} | 167 | R-3c | D3d^6 | 32 | 3 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 32 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 64 | 3 |
/// | dir{C1} | 64 | Cmca | D2h^18 | 48 | 3 |
/// | dir{C1} | 63 | Cmcm | D2h^17 | 48 | 3 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 48 | 3 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 48 | 3 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 96 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 3 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 96 | 3 |
/// | dir{C1} | 14 | P2_1/c | C2h^5 | 96 | 3 |
/// | dir{P3} | 11 | P2_1/m | C2h^2 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{11}t1$ | 3 | C24b | no |
/// | M1- | $M_{1}^-$ | $k_{11}t2$ | 3 | C48a | no |
/// | M2+ | $M_{3}^+$ | $k_{11}t5$ | 3 | C24c | no |
/// | M2- | $M_{3}^-$ | $k_{11}t6$ | 3 | C48a | no |
/// | M3+ | $M_{2}^+$ | $k_{11}t3$ | 3 | C24c | no |
/// | M3- | $M_{2}^-$ | $k_{11}t4$ | 3 | C48a | no |
/// | M4+ | $M_{4}^+$ | $k_{11}t7$ | 3 | C24b | no |
/// | M4- | $M_{4}^-$ | $k_{11}t8$ | 3 | C48a | no |
/// | M5+ | $M_{5}^+$ | $k_{11}t9$ | 6 | E96f | no |
/// | M5- | $M_{5}^-$ | $k_{11}t10$ | 6 | E192f | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 6 | 1 |
/// | dir{P1} | 229 | Im-3m | Oh^9 | 4 | 3 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 12 | 3 |
/// | dir{P1} | 71 | Immm | D2h^25 | 24 | 3 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 6 | 1 |
/// | dir{P1} | 211 | I432 | O^5 | 8 | 3 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 12 | 3 |
/// | dir{P1} | 97 | I422 | D4^9 | 24 | 3 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 24 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 3 |
///
/// ### Isotropy subgroups for M2+ ($M_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 127 | P4/mbm | D4h^5 | 6 | 1 |
/// | dir{P1} | 204 | Im-3 | Th^5 | 8 | 3 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 12 | 3 |
/// | dir{P1} | 71 | Immm | D2h^25 | 24 | 3 |
///
/// ### Isotropy subgroups for M2- ($M_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 6 | 1 |
/// | dir{P1} | 217 | I-43m | Td^3 | 8 | 3 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 12 | 3 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 24 | 3 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 24 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 3 |
///
/// ### Isotropy subgroups for M3+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 127 | P4/mbm | D4h^5 | 6 | 1 |
/// | dir{P1} | 204 | Im-3 | Th^5 | 8 | 3 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 12 | 3 |
/// | dir{P1} | 71 | Immm | D2h^25 | 24 | 3 |
///
/// ### Isotropy subgroups for M3- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 6 | 1 |
/// | dir{P1} | 217 | I-43m | Td^3 | 8 | 3 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 12 | 3 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 24 | 3 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 24 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 3 |
///
/// ### Isotropy subgroups for M4+ ($M_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 6 | 1 |
/// | dir{P1} | 229 | Im-3m | Oh^9 | 4 | 3 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 12 | 3 |
/// | dir{P1} | 71 | Immm | D2h^25 | 24 | 3 |
///
/// ### Isotropy subgroups for M4- ($M_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 6 | 1 |
/// | dir{P1} | 211 | I432 | O^5 | 8 | 3 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 12 | 3 |
/// | dir{P1} | 97 | I422 | D4^9 | 24 | 3 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 24 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 3 |
///
/// ### Isotropy subgroups for M5+ ($M_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 12 | 1 |
/// | dir{P1} | 53 | Pmna | D2h^7 | 12 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 24 | 1 |
/// | dir{C1} | 206 | Ia-3 | Th^7 | 8 | 3 |
/// | dir{P3} | 166 | R-3m | D3d^5 | 16 | 3 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 12 | 3 |
/// | dir{P3} | 140 | I4/mcm | D4h^18 | 12 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 24 | 3 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 24 | 3 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 24 | 3 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
///
/// ### Isotropy subgroups for M5- ($M_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 12 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 12 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 24 | 1 |
/// | dir{C1} | 199 | I2_13 | T^5 | 16 | 3 |
/// | dir{P3} | 160 | R3m | C3v^5 | 32 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P3} | 139 | I4/mmm | D4h^17 | 12 | 3 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 12 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 74 | Imma | D2h^28 | 24 | 3 |
/// | dir{P1} | 71 | Immm | D2h^25 | 24 | 3 |
/// | dir{P1} | 71 | Immm | D2h^25 | 24 | 3 |
/// | dir{C1} | 42 | Fmm2 | C2v^18 | 48 | 3 |
/// | dir{C1} | 42 | Fmm2 | C2v^18 | 48 | 3 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 48 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 3 |
/// | dir{P3} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{C1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{P3} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{13}t1$ | 1 | A2a | no |
/// | R1- | $R_{1}^-$ | $k_{13}t2$ | 1 | A2a | no |
/// | R2+ | $R_{2}^+$ | $k_{13}t3$ | 1 | A2a | no |
/// | R2- | $R_{2}^-$ | $k_{13}t4$ | 1 | A2a | no |
/// | R3+ | $R_{3}^+$ | $k_{13}t5$ | 2 | B12a | no |
/// | R3- | $R_{3}^-$ | $k_{13}t6$ | 2 | B12a | no |
/// | R4+ | $R_{4}^+$ | $k_{13}t9$ | 3 | C48a | no |
/// | R4- | $R_{4}^-$ | $k_{13}t10$ | 3 | C48a | no |
/// | R5+ | $R_{5}^+$ | $k_{13}t7$ | 3 | C48a | no |
/// | R5- | $R_{5}^-$ | $k_{13}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 225 | Fm-3m | Oh^5 | 2 | 1 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 226 | Fm-3c | Oh^6 | 2 | 1 |
///
/// ### Isotropy subgroups for R2+ ($R_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 226 | Fm-3c | Oh^6 | 2 | 1 |
///
/// ### Isotropy subgroups for R2- ($R_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 225 | Fm-3m | Oh^5 | 2 | 1 |
///
/// ### Isotropy subgroups for R3+ ($R_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 6 | 1 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 6 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 12 | 1 |
///
/// ### Isotropy subgroups for R3- ($R_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 6 | 1 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 6 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 12 | 1 |
///
/// ### Isotropy subgroups for R4+ ($R_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 8 | 1 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 6 | 1 |
/// | dir{P1} | 74 | Imma | D2h^28 | 12 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
///
/// ### Isotropy subgroups for R4- ($R_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 8 | 1 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 6 | 1 |
/// | dir{P1} | 74 | Imma | D2h^28 | 12 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
///
/// ### Isotropy subgroups for R5+ ($R_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 8 | 1 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 6 | 1 |
/// | dir{P1} | 74 | Imma | D2h^28 | 12 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
///
/// ### Isotropy subgroups for R5- ($R_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 8 | 1 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 6 | 1 |
/// | dir{P1} | 74 | Imma | D2h^28 | 12 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
///
pub struct Sg221;

/// # 222 Pn-3n (Oₕ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{12}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{12}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{12}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{12}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{12}t6$ | 2 | B12a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{12}t9$ | 3 | C24c | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{12}t10$ | 3 | C48a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{12}t7$ | 3 | C24b | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{12}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 222 | Pn-3n | Oh^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 207 | P432 | O^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 201 | Pn-3 | Th^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 218 | P-43n | Td^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 3 | 1 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 6 | 1 |
/// | dir{P1} | 89 | P422 | D4^1 | 6 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 85 | P4/n | C4h^3 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 8 | 1 |
/// | dir{P1} | 104 | P4nc | C4v^6 | 6 | 1 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 4 | 1 |
/// | dir{P1} | 68 | Ccca | D2h^22 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 6 | 1 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 12 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{2}$ | $k_{10}t2$ | 6 | E384b | yes |
/// | X2 | $X_{1}$ | $k_{10}t1$ | 6 | E384b | yes |
/// | X3X4 | $X_{3}X4$ | $k_{10}t3t4$ | 12 | G384a | no |
///
/// ### Isotropy subgroups for X1 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 12 | 1 |
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 12 | 1 |
/// | dir{P1} | 77 | P4_2 | C4^3 | 24 | 1 |
/// | dir{P3} | 112 | P-42c | D2d^2 | 24 | 2 |
/// | dir{C1} | 89 | P422 | D4^1 | 24 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 2 |
/// | dir{C1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 96 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{C1} | 208 | P4_232 | O^2 | 16 | 3 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 64 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 128 | 3 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 48 | 3 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 48 | 3 |
/// | dir{C1} | 81 | P-4 | S4^1 | 96 | 3 |
/// | dir{P1} | 77 | P4_2 | C4^3 | 96 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 192 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 3 |
/// | dir{C1} | 3 | P2 | C2^1 | 192 | 3 |
/// | dir{P3} | 3 | P2 | C2^1 | 192 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 12 | 1 |
/// | dir{P1} | 85 | P4/n | C4h^3 | 12 | 1 |
/// | dir{P1} | 75 | P4 | C4^1 | 24 | 1 |
/// | dir{P3} | 112 | P-42c | D2d^2 | 24 | 2 |
/// | dir{C1} | 89 | P422 | D4^1 | 24 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 2 |
/// | dir{C1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 96 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{C1} | 207 | P432 | O^1 | 16 | 3 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 64 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 128 | 3 |
/// | dir{P1} | 89 | P422 | D4^1 | 48 | 3 |
/// | dir{P1} | 89 | P422 | D4^1 | 48 | 3 |
/// | dir{C1} | 81 | P-4 | S4^1 | 96 | 3 |
/// | dir{P1} | 75 | P4 | C4^1 | 96 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 192 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 3 |
/// | dir{C1} | 3 | P2 | C2^1 | 192 | 3 |
/// | dir{P3} | 3 | P2 | C2^1 | 192 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 3 |
///
/// ### Isotropy subgroups for X3X4 ($X_{3}X4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1} | 20 | C222_1 | D2^5 | 24 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 24 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 24 | 1 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 48 | 1 |
/// | dir{P3} | 114 | P-42_1c | D2d^4 | 24 | 2 |
/// | dir{P1} | 90 | P42_12 | D4^2 | 24 | 2 |
/// | dir{P3} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 2 |
/// | dir{P11} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{C2} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{C1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 96 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{C1} | 198 | P2_13 | T^4 | 32 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 64 | 3 |
/// | dir{C1} | 148 | R-3 | C3i^2 | 64 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 128 | 3 |
/// | dir{C1} | 20 | C222_1 | D2^5 | 96 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 3 |
/// | dir{C1} | 2 | P-1 | Ci^1 | 192 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 3 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 192 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{11}t1$ | 6 | E192d | no |
/// | M2 | $M_{3}$ | $k_{11}t2$ | 6 | E192e | no |
/// | M3 | $M_{4}$ | $k_{11}t3$ | 6 | E192d | yes |
/// | M4 | $M_{2}$ | $k_{11}t4$ | 6 | E192e | yes |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 89 | P422 | D4^1 | 12 | 1 |
/// | dir{P1} | 50 | Pban | D2h^4 | 12 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 1 |
/// | dir{C1} | 211 | I432 | O^5 | 8 | 3 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 16 | 3 |
/// | dir{P3} | 122 | I-42d | D2d^12 | 24 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 97 | I422 | D4^9 | 24 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 3 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 90 | P42_12 | D4^2 | 12 | 1 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 12 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 1 |
/// | dir{C1} | 197 | I23 | T^3 | 16 | 3 |
/// | dir{P3} | 161 | R3c | C3v^6 | 32 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P3} | 122 | I-42d | D2d^12 | 24 | 3 |
/// | dir{C1} | 97 | I422 | D4^9 | 24 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 79 | I4 | C4^5 | 48 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 3 |
/// | dir{C1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 12 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 12 | 1 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 24 | 1 |
/// | dir{C1} | 220 | I-43d | Td^6 | 8 | 3 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 16 | 3 |
/// | dir{P3} | 97 | I422 | D4^9 | 24 | 3 |
/// | dir{P1} | 161 | R3c | C3v^6 | 32 | 3 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 24 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 48 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 12 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 12 | 1 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 24 | 1 |
/// | dir{C1} | 199 | I2_13 | T^5 | 16 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{C1} | 122 | I-42d | D2d^12 | 24 | 3 |
/// | dir{P3} | 97 | I422 | D4^9 | 24 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | R2R3 | $R_{2}R3$ | $k_{13}t2t3$ | 4 | D24e | no |
/// | R4 | $R_{4}$ | $k_{13}t4$ | 6 | E96d | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 209 | F432 | O^3 | 4 | 1 |
/// | dir{P1} | 203 | Fd-3 | Th^4 | 4 | 1 |
/// | dir{P1} | 196 | F23 | T^2 | 8 | 1 |
///
/// ### Isotropy subgroups for R2R3 ($R_{2}R3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 97 | I422 | D4^9 | 12 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 12 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 24 | 1 |
///
/// ### Isotropy subgroups for R4 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1} | 155 | R32 | D3^7 | 16 | 1 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 16 | 1 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 12 | 1 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 12 | 1 |
/// | dir{P1} | 97 | I422 | D4^9 | 12 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 12 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 24 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 24 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 24 | 1 |
/// | dir{C1} | 24 | I2_12_12_1 | D2^9 | 24 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 24 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 48 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
/// | dir{P3} | 5 | C2 | C2^3 | 48 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 1 |
///
pub struct Sg222;

/// # 223 Pm-3n (Oₕ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{12}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{12}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{12}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{12}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{12}t6$ | 2 | B12a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{12}t9$ | 3 | C24c | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{12}t10$ | 3 | C48a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{12}t7$ | 3 | C24b | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{12}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 223 | Pm-3n | Oh^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 208 | P4_232 | O^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 200 | Pm-3 | Th^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 218 | P-43n | Td^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 3 | 1 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 6 | 1 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 6 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 8 | 1 |
/// | dir{P1} | 105 | P4_2mc | C4v^7 | 6 | 1 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 4 | 1 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 6 | 1 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 12 | 1 |
/// | dir{P1} | 6 | Pm | Cs^1 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{10}t3$ | 6 | E384a | yes |
/// | X2 | $X_{3}$ | $k_{10}t4$ | 6 | E384a | yes |
/// | X3 | $X_{2}$ | $k_{10}t1$ | 6 | E96d | no |
/// | X4 | $X_{4}$ | $k_{10}t2$ | 6 | E96d | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 12 | 1 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 12 | 1 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 24 | 1 |
/// | dir{P3} | 131 | P4_2/mmc | D4h^9 | 12 | 2 |
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 12 | 2 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 24 | 2 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 24 | 2 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 24 | 2 |
/// | dir{C1} | 40 | Ama2 | C2v^16 | 48 | 2 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 48 | 2 |
/// | dir{P3} | 25 | Pmm2 | C2v^1 | 48 | 2 |
/// | dir{P3} | 6 | Pm | Cs^1 | 96 | 2 |
/// | dir{C1} | 200 | Pm-3 | Th^1 | 16 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 64 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 128 | 3 |
/// | dir{C1} | 115 | P-4m2 | D2d^5 | 48 | 3 |
/// | dir{C1} | 115 | P-4m2 | D2d^5 | 48 | 3 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 48 | 3 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 3 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 96 | 3 |
/// | dir{P3} | 25 | Pmm2 | C2v^1 | 96 | 3 |
/// | dir{C1} | 25 | Pmm2 | C2v^1 | 96 | 3 |
/// | dir{P1} | 6 | Pm | Cs^1 | 192 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 12 | 1 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 12 | 1 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 24 | 1 |
/// | dir{P3} | 133 | P4_2/nbc | D4h^11 | 12 | 2 |
/// | dir{P1} | 133 | P4_2/nbc | D4h^11 | 12 | 2 |
/// | dir{P1} | 50 | Pban | D2h^4 | 24 | 2 |
/// | dir{P1} | 50 | Pban | D2h^4 | 24 | 2 |
/// | dir{P1} | 50 | Pban | D2h^4 | 24 | 2 |
/// | dir{C1} | 41 | Aba2 | C2v^17 | 48 | 2 |
/// | dir{P3} | 30 | Pnc2 | C2v^6 | 48 | 2 |
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 48 | 2 |
/// | dir{P3} | 7 | Pc | Cs^2 | 96 | 2 |
/// | dir{C1} | 201 | Pn-3 | Th^2 | 16 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 64 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 128 | 3 |
/// | dir{C1} | 118 | P-4n2 | D2d^8 | 48 | 3 |
/// | dir{C1} | 118 | P-4n2 | D2d^8 | 48 | 3 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 48 | 3 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 3 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 96 | 3 |
/// | dir{P3} | 34 | Pnn2 | C2v^10 | 96 | 3 |
/// | dir{C1} | 34 | Pnn2 | C2v^10 | 96 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 192 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 3 |
///
/// ### Isotropy subgroups for X3 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 12 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 12 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 24 | 1 |
/// | dir{P1} | 137 | P4_2/nmc | D4h^15 | 12 | 2 |
/// | dir{P1} | 135 | P4_2/mbc | D4h^13 | 12 | 2 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 24 | 2 |
/// | dir{C1} | 94 | P4_22_12 | D4^6 | 24 | 2 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 24 | 2 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 24 | 2 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 24 | 2 |
/// | dir{P3} | 18 | P2_12_12 | D2^3 | 48 | 2 |
/// | dir{C1} | 213 | P4_132 | O^7 | 16 | 3 |
/// | dir{P3} | 205 | Pa-3 | Th^6 | 16 | 3 |
/// | dir{P1} | 198 | P2_13 | T^4 | 32 | 3 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 48 | 3 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 48 | 3 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 48 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 3 |
///
/// ### Isotropy subgroups for X4 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 12 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 12 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 24 | 1 |
/// | dir{P1} | 137 | P4_2/nmc | D4h^15 | 12 | 2 |
/// | dir{P1} | 135 | P4_2/mbc | D4h^13 | 12 | 2 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 24 | 2 |
/// | dir{C1} | 94 | P4_22_12 | D4^6 | 24 | 2 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 24 | 2 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 24 | 2 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 24 | 2 |
/// | dir{P3} | 18 | P2_12_12 | D2^3 | 48 | 2 |
/// | dir{C1} | 213 | P4_132 | O^7 | 16 | 3 |
/// | dir{P3} | 205 | Pa-3 | Th^6 | 16 | 3 |
/// | dir{P1} | 198 | P2_13 | T^4 | 32 | 3 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 48 | 3 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 48 | 3 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 48 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{11}t5$ | 3 | C24b | no |
/// | M1- | $M_{1}^-$ | $k_{11}t6$ | 3 | C48a | no |
/// | M2+ | $M_{3}^+$ | $k_{11}t1$ | 3 | C24c | no |
/// | M2- | $M_{3}^-$ | $k_{11}t2$ | 3 | C48a | no |
/// | M3+ | $M_{2}^+$ | $k_{11}t7$ | 3 | C24c | no |
/// | M3- | $M_{2}^-$ | $k_{11}t8$ | 3 | C48a | no |
/// | M4+ | $M_{4}^+$ | $k_{11}t3$ | 3 | C24b | no |
/// | M4- | $M_{4}^-$ | $k_{11}t4$ | 3 | C48a | no |
/// | M5+ | $M_{5}^+$ | $k_{11}t9$ | 6 | E96f | no |
/// | M5- | $M_{5}^-$ | $k_{11}t10$ | 6 | E192f | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 6 | 1 |
/// | dir{P1} | 230 | Ia-3d | Oh^10 | 4 | 3 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 12 | 3 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 24 | 3 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 6 | 1 |
/// | dir{P1} | 214 | I4_132 | O^8 | 8 | 3 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 12 | 3 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 24 | 3 |
/// | dir{P1} | 74 | Imma | D2h^28 | 24 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 3 |
///
/// ### Isotropy subgroups for M2+ ($M_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 138 | P4_2/ncm | D4h^16 | 6 | 1 |
/// | dir{P1} | 206 | Ia-3 | Th^7 | 8 | 3 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 12 | 3 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 24 | 3 |
///
/// ### Isotropy subgroups for M2- ($M_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 6 | 1 |
/// | dir{P1} | 220 | I-43d | Td^6 | 8 | 3 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 12 | 3 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 24 | 3 |
/// | dir{P1} | 74 | Imma | D2h^28 | 24 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 3 |
///
/// ### Isotropy subgroups for M3+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 138 | P4_2/ncm | D4h^16 | 6 | 1 |
/// | dir{P1} | 206 | Ia-3 | Th^7 | 8 | 3 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 12 | 3 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 24 | 3 |
///
/// ### Isotropy subgroups for M3- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 6 | 1 |
/// | dir{P1} | 220 | I-43d | Td^6 | 8 | 3 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 12 | 3 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 24 | 3 |
/// | dir{P1} | 74 | Imma | D2h^28 | 24 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 3 |
///
/// ### Isotropy subgroups for M4+ ($M_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 6 | 1 |
/// | dir{P1} | 230 | Ia-3d | Oh^10 | 4 | 3 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 12 | 3 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 24 | 3 |
///
/// ### Isotropy subgroups for M4- ($M_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 6 | 1 |
/// | dir{P1} | 214 | I4_132 | O^8 | 8 | 3 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 12 | 3 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 24 | 3 |
/// | dir{P1} | 74 | Imma | D2h^28 | 24 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 3 |
///
/// ### Isotropy subgroups for M5+ ($M_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 12 | 1 |
/// | dir{P1} | 53 | Pmna | D2h^7 | 12 | 1 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 24 | 1 |
/// | dir{C1} | 204 | Im-3 | Th^5 | 8 | 3 |
/// | dir{P3} | 167 | R-3c | D3d^6 | 16 | 3 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 12 | 3 |
/// | dir{P3} | 141 | I4_1/amd | D4h^19 | 12 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P1} | 74 | Imma | D2h^28 | 24 | 3 |
/// | dir{P1} | 74 | Imma | D2h^28 | 24 | 3 |
/// | dir{P1} | 71 | Immm | D2h^25 | 24 | 3 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
///
/// ### Isotropy subgroups for M5- ($M_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 12 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 12 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 24 | 1 |
/// | dir{C1} | 197 | I23 | T^3 | 16 | 3 |
/// | dir{P3} | 161 | R3c | C3v^6 | 32 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P3} | 142 | I4_1/acd | D4h^20 | 12 | 3 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 12 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 24 | 3 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 24 | 3 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 24 | 3 |
/// | dir{C1} | 43 | Fdd2 | C2v^19 | 48 | 3 |
/// | dir{C1} | 43 | Fdd2 | C2v^19 | 48 | 3 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 48 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P3} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{C1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P3} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | R2R3 | $R_{2}R3$ | $k_{13}t2t3$ | 4 | D24e | no |
/// | R4 | $R_{4}$ | $k_{13}t4$ | 6 | E96d | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 210 | F4_132 | O^4 | 4 | 1 |
/// | dir{P1} | 202 | Fm-3 | Th^3 | 4 | 1 |
/// | dir{P1} | 196 | F23 | T^2 | 8 | 1 |
///
/// ### Isotropy subgroups for R2R3 ($R_{2}R3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 12 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 12 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 24 | 1 |
///
/// ### Isotropy subgroups for R4 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1} | 155 | R32 | D3^7 | 16 | 1 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 16 | 1 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 12 | 1 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 12 | 1 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 12 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 12 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 1 |
/// | dir{P1} | 82 | I-4 | S4^2 | 24 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 24 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 24 | 1 |
/// | dir{C1} | 23 | I222 | D2^8 | 24 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 24 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 48 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
/// | dir{P3} | 5 | C2 | C2^3 | 48 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 1 |
///
pub struct Sg223;

/// # 224 Pn-3m (Oₕ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{12}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{12}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{12}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{12}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{12}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{12}t6$ | 2 | B12a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{12}t9$ | 3 | C24c | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{12}t10$ | 3 | C48a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{12}t7$ | 3 | C24b | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{12}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 224 | Pn-3m | Oh^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 208 | P4_232 | O^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 201 | Pn-3 | Th^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 215 | P-43m | Td^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 3 | 1 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 6 | 1 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 6 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 86 | P4_2/n | C4h^4 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 8 | 1 |
/// | dir{P1} | 102 | P4_2nm | C4v^4 | 6 | 1 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 12 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 4 | 1 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 6 | 1 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 12 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{4}$ | $k_{10}t3$ | 6 | E384c | yes |
/// | X2 | $X_{2}$ | $k_{10}t4$ | 6 | E384c | yes |
/// | X3 | $X_{1}$ | $k_{10}t2$ | 6 | E384c | no |
/// | X4 | $X_{3}$ | $k_{10}t1$ | 6 | E384c | no |
///
/// ### Isotropy subgroups for X1 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 111 | P-42m | D2d^1 | 12 | 1 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 12 | 1 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 24 | 1 |
/// | dir{C1} | 111 | P-42m | D2d^1 | 24 | 2 |
/// | dir{P3} | 93 | P4_222 | D4^5 | 24 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 96 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{C1} | 215 | P-43m | Td^1 | 16 | 3 |
/// | dir{P3} | 166 | R-3m | D3d^5 | 32 | 3 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 32 | 3 |
/// | dir{P1} | 160 | R3m | C3v^5 | 64 | 3 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 48 | 3 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 48 | 3 |
/// | dir{C1} | 21 | C222 | D2^6 | 96 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 3 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 96 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 192 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 3 |
/// | dir{C1} | 8 | Cm | Cs^3 | 192 | 3 |
/// | dir{P3} | 3 | P2 | C2^1 | 192 | 3 |
/// | dir{C1} | 3 | P2 | C2^1 | 192 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 112 | P-42c | D2d^2 | 12 | 1 |
/// | dir{P1} | 68 | Ccca | D2h^22 | 12 | 1 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 24 | 1 |
/// | dir{C1} | 111 | P-42m | D2d^1 | 24 | 2 |
/// | dir{P3} | 93 | P4_222 | D4^5 | 24 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 2 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 96 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{C1} | 218 | P-43n | Td^4 | 16 | 3 |
/// | dir{P3} | 167 | R-3c | D3d^6 | 32 | 3 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 32 | 3 |
/// | dir{P1} | 161 | R3c | C3v^6 | 64 | 3 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 48 | 3 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 48 | 3 |
/// | dir{C1} | 21 | C222 | D2^6 | 96 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 3 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 96 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 192 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 3 |
/// | dir{C1} | 9 | Cc | Cs^4 | 192 | 3 |
/// | dir{P3} | 3 | P2 | C2^1 | 192 | 3 |
/// | dir{C1} | 3 | P2 | C2^1 | 192 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 3 |
///
/// ### Isotropy subgroups for X3 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 12 | 1 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 1 |
/// | dir{P3} | 113 | P-42_1m | D2d^3 | 24 | 2 |
/// | dir{C1} | 94 | P4_22_12 | D4^6 | 24 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 2 |
/// | dir{C1} | 8 | Cm | Cs^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 96 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{C1} | 212 | P4_332 | O^6 | 16 | 3 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 32 | 3 |
/// | dir{P3} | 166 | R-3m | D3d^5 | 32 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 64 | 3 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 48 | 3 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 48 | 3 |
/// | dir{C1} | 36 | Cmc2_1 | C2v^12 | 96 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 96 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 192 | 3 |
/// | dir{C1} | 8 | Cm | Cs^3 | 192 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 192 | 3 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 192 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 3 |
///
/// ### Isotropy subgroups for X4 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 91 | P4_122 | D4^3 | 12 | 1 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 1 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 1 |
/// | dir{P3} | 113 | P-42_1m | D2d^3 | 24 | 2 |
/// | dir{C1} | 94 | P4_22_12 | D4^6 | 24 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 2 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 2 |
/// | dir{C1} | 8 | Cm | Cs^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 96 | 2 |
/// | dir{P3} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{C1} | 213 | P4_132 | O^7 | 16 | 3 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 32 | 3 |
/// | dir{P3} | 166 | R-3m | D3d^5 | 32 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 64 | 3 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 48 | 3 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 48 | 3 |
/// | dir{C1} | 36 | Cmc2_1 | C2v^12 | 96 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 96 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 96 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 192 | 3 |
/// | dir{C1} | 8 | Cm | Cs^3 | 192 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 192 | 3 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 192 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 3 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{4}$ | $k_{11}t4$ | 6 | E192d | yes |
/// | M2 | $M_{2}$ | $k_{11}t3$ | 6 | E192e | yes |
/// | M3 | $M_{3}$ | $k_{11}t1$ | 6 | E192d | no |
/// | M4 | $M_{1}$ | $k_{11}t2$ | 6 | E192e | no |
///
/// ### Isotropy subgroups for M1 ($M_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 12 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 12 | 1 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 24 | 1 |
/// | dir{C1} | 217 | I-43m | Td^3 | 8 | 3 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 16 | 3 |
/// | dir{P3} | 98 | I4_122 | D4^10 | 24 | 3 |
/// | dir{P1} | 160 | R3m | C3v^5 | 32 | 3 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 24 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 48 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 12 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 12 | 1 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 24 | 1 |
/// | dir{C1} | 197 | I23 | T^3 | 16 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{C1} | 121 | I-42m | D2d^11 | 24 | 3 |
/// | dir{P3} | 98 | I4_122 | D4^10 | 24 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for M3 ($M_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 93 | P4_222 | D4^5 | 12 | 1 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 12 | 1 |
/// | dir{P1} | 16 | P222 | D2^1 | 24 | 1 |
/// | dir{C1} | 214 | I4_132 | O^8 | 8 | 3 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 16 | 3 |
/// | dir{P3} | 121 | I-42m | D2d^11 | 24 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 24 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 3 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for M4 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 12 | 1 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 12 | 1 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 24 | 1 |
/// | dir{C1} | 199 | I2_13 | T^5 | 16 | 3 |
/// | dir{P3} | 160 | R3m | C3v^5 | 32 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P3} | 121 | I-42m | D2d^11 | 24 | 3 |
/// | dir{C1} | 98 | I4_122 | D4^10 | 24 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 48 | 3 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 3 |
/// | dir{C1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{13}t4$ | 1 | A2a | no |
/// | R1- | $R_{1}^-$ | $k_{13}t3$ | 1 | A2a | no |
/// | R2+ | $R_{2}^+$ | $k_{13}t2$ | 1 | A2a | no |
/// | R2- | $R_{2}^-$ | $k_{13}t1$ | 1 | A2a | no |
/// | R3+ | $R_{3}^+$ | $k_{13}t6$ | 2 | B12a | no |
/// | R3- | $R_{3}^-$ | $k_{13}t5$ | 2 | B12a | no |
/// | R4+ | $R_{4}^+$ | $k_{13}t8$ | 3 | C48a | no |
/// | R4- | $R_{4}^-$ | $k_{13}t7$ | 3 | C48a | no |
/// | R5+ | $R_{5}^+$ | $k_{13}t10$ | 3 | C48a | no |
/// | R5- | $R_{5}^-$ | $k_{13}t9$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 227 | Fd-3m | Oh^7 | 2 | 1 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 228 | Fd-3c | Oh^8 | 2 | 1 |
///
/// ### Isotropy subgroups for R2+ ($R_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 228 | Fd-3c | Oh^8 | 2 | 1 |
///
/// ### Isotropy subgroups for R2- ($R_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 227 | Fd-3m | Oh^7 | 2 | 1 |
///
/// ### Isotropy subgroups for R3+ ($R_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 6 | 1 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 6 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 12 | 1 |
///
/// ### Isotropy subgroups for R3- ($R_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 6 | 1 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 6 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 12 | 1 |
///
/// ### Isotropy subgroups for R4+ ($R_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 8 | 1 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 6 | 1 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 12 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
///
/// ### Isotropy subgroups for R4- ($R_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 8 | 1 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 6 | 1 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 12 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
///
/// ### Isotropy subgroups for R5+ ($R_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 8 | 1 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 6 | 1 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 12 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
///
/// ### Isotropy subgroups for R5- ($R_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 8 | 1 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 6 | 1 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 12 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
///
pub struct Sg224;

/// # 225 Fm-3m (Oₕ⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{11}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{11}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{11}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{11}t6$ | 2 | B12a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{11}t9$ | 3 | C24c | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{11}t10$ | 3 | C48a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{11}t7$ | 3 | C24b | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{11}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 225 | Fm-3m | Oh^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 209 | F432 | O^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 202 | Fm-3 | Th^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 216 | F-43m | Td^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 3 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 6 | 1 |
/// | dir{P1} | 97 | I422 | D4^9 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 8 | 1 |
/// | dir{P1} | 107 | I4mm | C4v^9 | 6 | 1 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 12 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 4 | 1 |
/// | dir{P1} | 71 | Immm | D2h^25 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 6 | 1 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 12 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{10}t1$ | 3 | C24b | no |
/// | X1- | $X_{1}^-$ | $k_{10}t2$ | 3 | C48a | no |
/// | X2+ | $X_{3}^+$ | $k_{10}t5$ | 3 | C24c | no |
/// | X2- | $X_{3}^-$ | $k_{10}t6$ | 3 | C48a | no |
/// | X3+ | $X_{2}^+$ | $k_{10}t3$ | 3 | C24c | no |
/// | X3- | $X_{2}^-$ | $k_{10}t4$ | 3 | C48a | no |
/// | X4+ | $X_{4}^+$ | $k_{10}t7$ | 3 | C24b | no |
/// | X4- | $X_{4}^-$ | $k_{10}t8$ | 3 | C48a | no |
/// | X5+ | $X_{5}^+$ | $k_{10}t9$ | 6 | E96f | no |
/// | X5- | $X_{5}^-$ | $k_{10}t10$ | 6 | E192f | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 6 | 1 |
/// | dir{P1} | 221 | Pm-3m | Oh^1 | 4 | 3 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 12 | 3 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 6 | 1 |
/// | dir{P1} | 207 | P432 | O^1 | 8 | 3 |
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 12 | 3 |
/// | dir{P1} | 89 | P422 | D4^1 | 24 | 3 |
/// | dir{P1} | 50 | Pban | D2h^4 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 6 | 1 |
/// | dir{P1} | 200 | Pm-3 | Th^1 | 8 | 3 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 12 | 3 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 137 | P4_2/nmc | D4h^15 | 6 | 1 |
/// | dir{P1} | 215 | P-43m | Td^1 | 8 | 3 |
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 12 | 3 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 24 | 3 |
/// | dir{P1} | 50 | Pban | D2h^4 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X3+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 128 | P4/mnc | D4h^6 | 6 | 1 |
/// | dir{P1} | 201 | Pn-3 | Th^2 | 8 | 3 |
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 12 | 3 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 24 | 3 |
///
/// ### Isotropy subgroups for X3- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 6 | 1 |
/// | dir{P1} | 215 | P-43m | Td^1 | 8 | 3 |
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 12 | 3 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 24 | 3 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X4+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 6 | 1 |
/// | dir{P1} | 224 | Pn-3m | Oh^4 | 4 | 3 |
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 12 | 3 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 24 | 3 |
///
/// ### Isotropy subgroups for X4- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 6 | 1 |
/// | dir{P1} | 208 | P4_232 | O^2 | 8 | 3 |
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 12 | 3 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 24 | 3 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X5+ ($X_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 1 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 12 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 24 | 1 |
/// | dir{C1} | 205 | Pa-3 | Th^6 | 8 | 3 |
/// | dir{P3} | 166 | R-3m | D3d^5 | 16 | 3 |
/// | dir{P1} | 138 | P4_2/ncm | D4h^16 | 12 | 3 |
/// | dir{P3} | 127 | P4/mbm | D4h^5 | 12 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P1} | 56 | Pccn | D2h^10 | 24 | 3 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 24 | 3 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 24 | 3 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{C1} | 14 | P2_1/c | C2h^5 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
///
/// ### Isotropy subgroups for X5- ($X_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 1 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 12 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 1 |
/// | dir{C1} | 198 | P2_13 | T^4 | 16 | 3 |
/// | dir{P1} | 160 | R3m | C3v^5 | 32 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 12 | 3 |
/// | dir{P3} | 129 | P4/nmm | D4h^7 | 12 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 24 | 3 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 24 | 3 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 24 | 3 |
/// | dir{C1} | 39 | Abm2 | C2v^15 | 48 | 3 |
/// | dir{C1} | 38 | Amm2 | C2v^14 | 48 | 3 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 48 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
/// | dir{P3} | 14 | P2_1/c | C2h^5 | 48 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 48 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 96 | 3 |
/// | dir{P3} | 6 | Pm | Cs^1 | 96 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 96 | 3 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{L}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^+$ | $k_{9}t1$ | 4 | D192b | no |
/// | L1- | $L_{1}^-$ | $k_{9}t2$ | 4 | D192b | no |
/// | L2+ | $L_{2}^+$ | $k_{9}t3$ | 4 | D192b | no |
/// | L2- | $L_{2}^-$ | $k_{9}t4$ | 4 | D192b | no |
/// | L3+ | $L_{3}^+$ | $k_{9}t5$ | 8 | F192a | no |
/// | L3- | $L_{3}^-$ | $k_{9}t6$ | 8 | F192a | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 8 | 1 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 24 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 48 | 2 |
/// | dir{P1} | 227 | Fd-3m | Oh^7 | 8 | 4 |
/// | dir{P1} | 225 | Fm-3m | Oh^5 | 8 | 4 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 32 | 4 |
/// | dir{P1} | 74 | Imma | D2h^28 | 48 | 4 |
/// | dir{P1} | 71 | Immm | D2h^25 | 48 | 4 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 96 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 8 | 1 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 24 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 48 | 2 |
/// | dir{P1} | 228 | Fd-3c | Oh^8 | 8 | 4 |
/// | dir{P1} | 226 | Fm-3c | Oh^6 | 8 | 4 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 32 | 4 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 48 | 4 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 48 | 4 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 96 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L2+ ($L_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 8 | 1 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 24 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 48 | 2 |
/// | dir{P1} | 228 | Fd-3c | Oh^8 | 8 | 4 |
/// | dir{P1} | 226 | Fm-3c | Oh^6 | 8 | 4 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 32 | 4 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 48 | 4 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 48 | 4 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 96 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L2- ($L_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 8 | 1 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 24 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 48 | 2 |
/// | dir{P1} | 227 | Fd-3m | Oh^7 | 8 | 4 |
/// | dir{P1} | 225 | Fm-3m | Oh^5 | 8 | 4 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 32 | 4 |
/// | dir{P1} | 74 | Imma | D2h^28 | 48 | 4 |
/// | dir{P1} | 71 | Immm | D2h^25 | 48 | 4 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 96 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L3+ ($L_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 24 | 2 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 24 | 2 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 48 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 48 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 32 | 4 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 32 | 4 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 24 | 4 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 24 | 4 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 24 | 4 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 24 | 4 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 64 | 4 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 48 | 4 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 48 | 4 |
/// | dir{C1} | 74 | Imma | D2h^28 | 48 | 4 |
/// | dir{P1} | 74 | Imma | D2h^28 | 48 | 4 |
/// | dir{P3} | 73 | Ibca | D2h^27 | 48 | 4 |
/// | dir{C1} | 72 | Ibam | D2h^26 | 48 | 4 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 48 | 4 |
/// | dir{P3} | 71 | Immm | D2h^25 | 48 | 4 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 48 | 4 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 48 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 4 |
/// | dir{P3} | 12 | C2/m | C2h^3 | 96 | 4 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 96 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
///
/// ### Isotropy subgroups for L3- ($L_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
/// | dir{P1} | 67 | Cmma | D2h^21 | 24 | 2 |
/// | dir{P1} | 65 | Cmmm | D2h^19 | 24 | 2 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 48 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 48 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 32 | 4 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 32 | 4 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 24 | 4 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 24 | 4 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 24 | 4 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 24 | 4 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 64 | 4 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 48 | 4 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 48 | 4 |
/// | dir{C1} | 74 | Imma | D2h^28 | 48 | 4 |
/// | dir{P1} | 74 | Imma | D2h^28 | 48 | 4 |
/// | dir{P3} | 73 | Ibca | D2h^27 | 48 | 4 |
/// | dir{C1} | 72 | Ibam | D2h^26 | 48 | 4 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 48 | 4 |
/// | dir{P3} | 71 | Immm | D2h^25 | 48 | 4 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 48 | 4 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 48 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 4 |
/// | dir{P3} | 12 | C2/m | C2h^3 | 96 | 4 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 96 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1 | $W_{1}$ | $k_{8}t1$ | 6 | E1536a | no |
/// | W2 | $W_{4}$ | $k_{8}t3$ | 6 | E1536a | no |
/// | W3 | $W_{2}$ | $k_{8}t2$ | 6 | E1536a | no |
/// | W4 | $W_{3}$ | $k_{8}t4$ | 6 | E1536a | no |
/// | W5 | $W_{5}$ | $k_{8}t5$ | 12 | G1536a | yes |
///
/// ### Isotropy subgroups for W1 ($W_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 12 | 2 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 12 | 2 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 24 | 2 |
/// | dir{C1} | 136 | P4_2/mnm | D4h^14 | 48 | 4 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 48 | 4 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 96 | 4 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 96 | 4 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 96 | 4 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 192 | 4 |
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 192 | 4 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 192 | 4 |
/// | dir{P3} | 6 | Pm | Cs^1 | 384 | 4 |
/// | dir{P3} | 223 | Pm-3n | Oh^3 | 32 | 6 |
/// | dir{P1} | 221 | Pm-3m | Oh^1 | 32 | 6 |
/// | dir{C1} | 212 | P4_332 | O^6 | 64 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 256 | 6 |
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 96 | 6 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 96 | 6 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 192 | 6 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 192 | 6 |
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 192 | 6 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 192 | 6 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 192 | 6 |
/// | dir{P3} | 40 | Ama2 | C2v^16 | 384 | 6 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 384 | 6 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 384 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 384 | 6 |
/// | dir{C1} | 11 | P2_1/m | C2h^2 | 384 | 6 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 384 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 768 | 6 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 768 | 6 |
/// | dir{C1} | 6 | Pm | Cs^1 | 768 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 1536 | 6 |
///
/// ### Isotropy subgroups for W2 ($W_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 12 | 2 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 12 | 2 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 24 | 2 |
/// | dir{C1} | 136 | P4_2/mnm | D4h^14 | 48 | 4 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 48 | 4 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 96 | 4 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 96 | 4 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 96 | 4 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 192 | 4 |
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 192 | 4 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 192 | 4 |
/// | dir{P3} | 6 | Pm | Cs^1 | 384 | 4 |
/// | dir{P1} | 223 | Pm-3n | Oh^3 | 32 | 6 |
/// | dir{P3} | 221 | Pm-3m | Oh^1 | 32 | 6 |
/// | dir{C1} | 213 | P4_132 | O^7 | 64 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 256 | 6 |
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 96 | 6 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 96 | 6 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 192 | 6 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 192 | 6 |
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 192 | 6 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 192 | 6 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 192 | 6 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 384 | 6 |
/// | dir{P3} | 38 | Amm2 | C2v^14 | 384 | 6 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 384 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 384 | 6 |
/// | dir{C1} | 11 | P2_1/m | C2h^2 | 384 | 6 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 384 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 768 | 6 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 768 | 6 |
/// | dir{C1} | 6 | Pm | Cs^1 | 768 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 1536 | 6 |
///
/// ### Isotropy subgroups for W3 ($W_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 12 | 2 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 12 | 2 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 24 | 2 |
/// | dir{C1} | 138 | P4_2/ncm | D4h^16 | 48 | 4 |
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 48 | 4 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 96 | 4 |
/// | dir{P1} | 56 | Pccn | D2h^10 | 96 | 4 |
/// | dir{P1} | 50 | Pban | D2h^4 | 96 | 4 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 192 | 4 |
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 192 | 4 |
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 192 | 4 |
/// | dir{P3} | 7 | Pc | Cs^2 | 384 | 4 |
/// | dir{P1} | 224 | Pn-3m | Oh^4 | 32 | 6 |
/// | dir{P3} | 222 | Pn-3n | Oh^2 | 32 | 6 |
/// | dir{C1} | 213 | P4_132 | O^7 | 64 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 256 | 6 |
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 96 | 6 |
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 96 | 6 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 192 | 6 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 192 | 6 |
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 192 | 6 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 192 | 6 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 192 | 6 |
/// | dir{P3} | 41 | Aba2 | C2v^17 | 384 | 6 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 384 | 6 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 384 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 384 | 6 |
/// | dir{C1} | 14 | P2_1/c | C2h^5 | 384 | 6 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 384 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 768 | 6 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 768 | 6 |
/// | dir{C1} | 7 | Pc | Cs^2 | 768 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 1536 | 6 |
///
/// ### Isotropy subgroups for W4 ($W_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 12 | 2 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 12 | 2 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 24 | 2 |
/// | dir{C1} | 138 | P4_2/ncm | D4h^16 | 48 | 4 |
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 48 | 4 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 96 | 4 |
/// | dir{P1} | 56 | Pccn | D2h^10 | 96 | 4 |
/// | dir{P1} | 50 | Pban | D2h^4 | 96 | 4 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 192 | 4 |
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 192 | 4 |
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 192 | 4 |
/// | dir{P3} | 7 | Pc | Cs^2 | 384 | 4 |
/// | dir{P3} | 224 | Pn-3m | Oh^4 | 32 | 6 |
/// | dir{P1} | 222 | Pn-3n | Oh^2 | 32 | 6 |
/// | dir{C1} | 212 | P4_332 | O^6 | 64 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 256 | 6 |
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 96 | 6 |
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 96 | 6 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 192 | 6 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 192 | 6 |
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 192 | 6 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 192 | 6 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 192 | 6 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 384 | 6 |
/// | dir{P3} | 39 | Abm2 | C2v^15 | 384 | 6 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 384 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 384 | 6 |
/// | dir{C1} | 14 | P2_1/c | C2h^5 | 384 | 6 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 384 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 768 | 6 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 768 | 6 |
/// | dir{C1} | 7 | Pc | Cs^2 | 768 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 1536 | 6 |
///
/// ### Isotropy subgroups for W5 ($W_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 24 | 2 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 24 | 2 |
/// | dir{C1} | 74 | Imma | D2h^28 | 24 | 2 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 24 | 2 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 24 | 2 |
/// | dir{P3} | 69 | Fmmm | D2h^23 | 24 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 2 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 48 | 2 |
/// | dir{P3} | 43 | Fdd2 | C2v^19 | 48 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 48 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 2 |
/// | dir{C1} | 22 | F222 | D2^7 | 48 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 48 | 4 |
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 48 | 4 |
/// | dir{C1} | 129 | P4/nmm | D4h^7 | 48 | 4 |
/// | dir{P3} | 127 | P4/mbm | D4h^5 | 48 | 4 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 96 | 4 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 96 | 4 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 96 | 4 |
/// | dir{P3} | 90 | P42_12 | D4^2 | 96 | 4 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 96 | 4 |
/// | dir{C1} | 59 | Pmmn | D2h^13 | 96 | 4 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 96 | 4 |
/// | dir{P3} | 57 | Pbcm | D2h^11 | 96 | 4 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 96 | 4 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 96 | 4 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 96 | 4 |
/// | dir{P1} | 50 | Pban | D2h^4 | 96 | 4 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 96 | 4 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 96 | 4 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 192 | 4 |
/// | dir{P1} | 38 | Amm2 | C2v^14 | 192 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 192 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 192 | 4 |
/// | dir{C1} | 34 | Pnn2 | C2v^10 | 192 | 4 |
/// | dir{C11} | 32 | Pba2 | C2v^8 | 192 | 4 |
/// | dir{P3} | 31 | Pmn2_1 | C2v^7 | 192 | 4 |
/// | dir{C1} | 30 | Pnc2 | C2v^6 | 192 | 4 |
/// | dir{C1} | 29 | Pca2_1 | C2v^5 | 192 | 4 |
/// | dir{C12} | 28 | Pma2 | C2v^4 | 192 | 4 |
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 192 | 4 |
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 192 | 4 |
/// | dir{4D1} | 14 | P2_1/c | C2h^5 | 192 | 4 |
/// | dir{C1} | 14 | P2_1/c | C2h^5 | 192 | 4 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 192 | 4 |
/// | dir{P3} | 13 | P2/c | C2h^4 | 192 | 4 |
/// | dir{C1} | 13 | P2/c | C2h^4 | 192 | 4 |
/// | dir{C1} | 13 | P2/c | C2h^4 | 192 | 4 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 192 | 4 |
/// | dir{P1} | 10 | P2/m | C2h^1 | 192 | 4 |
/// | dir{C3} | 18 | P2_12_12 | D2^3 | 192 | 4 |
/// | dir{C1} | 18 | P2_12_12 | D2^3 | 192 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 4 |
/// | dir{P12} | 8 | Cm | Cs^3 | 384 | 4 |
/// | dir{P1} | 7 | Pc | Cs^2 | 384 | 4 |
/// | dir{4D1} | 7 | Pc | Cs^2 | 384 | 4 |
/// | dir{P1} | 7 | Pc | Cs^2 | 384 | 4 |
/// | dir{P1} | 7 | Pc | Cs^2 | 384 | 4 |
/// | dir{C11} | 7 | Pc | Cs^2 | 384 | 4 |
/// | dir{P1} | 6 | Pm | Cs^1 | 384 | 4 |
/// | dir{C8} | 6 | Pm | Cs^1 | 384 | 4 |
/// | dir{C2} | 5 | C2 | C2^3 | 384 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 384 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 384 | 4 |
/// | dir{C9} | 2 | P-1 | Ci^1 | 384 | 4 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 384 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 4 |
/// | dir{P1} | 218 | P-43n | Td^4 | 64 | 6 |
/// | dir{P3} | 215 | P-43m | Td^1 | 64 | 6 |
/// | dir{P1} | 205 | Pa-3 | Th^6 | 64 | 6 |
/// | dir{C1} | 205 | Pa-3 | Th^6 | 64 | 6 |
/// | dir{P3} | 167 | R-3c | D3d^6 | 128 | 6 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 128 | 6 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 128 | 6 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 128 | 6 |
/// | dir{C1} | 198 | P2_13 | T^4 | 128 | 6 |
/// | dir{C1} | 195 | P23 | T^1 | 128 | 6 |
/// | dir{C1} | 161 | R3c | C3v^6 | 256 | 6 |
/// | dir{P1} | 160 | R3m | C3v^5 | 256 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 256 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 256 | 6 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 256 | 6 |
/// | dir{C1} | 148 | R-3 | C3i^2 | 256 | 6 |
/// | dir{P1} | 68 | Ccca | D2h^22 | 192 | 6 |
/// | dir{P3} | 67 | Cmma | D2h^21 | 192 | 6 |
/// | dir{P1} | 66 | Cccm | D2h^20 | 192 | 6 |
/// | dir{C1} | 65 | Cmmm | D2h^19 | 192 | 6 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 192 | 6 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 192 | 6 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 192 | 6 |
/// | dir{C1} | 63 | Cmcm | D2h^17 | 192 | 6 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 192 | 6 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 192 | 6 |
/// | dir{P3} | 62 | Pnma | D2h^16 | 192 | 6 |
/// | dir{C1} | 62 | Pnma | D2h^16 | 192 | 6 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 192 | 6 |
/// | dir{C1} | 61 | Pbca | D2h^15 | 192 | 6 |
/// | dir{P3} | 41 | Aba2 | C2v^17 | 384 | 6 |
/// | dir{P3} | 41 | Aba2 | C2v^17 | 384 | 6 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 384 | 6 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 384 | 6 |
/// | dir{C1} | 39 | Abm2 | C2v^15 | 384 | 6 |
/// | dir{P1} | 39 | Abm2 | C2v^15 | 384 | 6 |
/// | dir{C1} | 38 | Amm2 | C2v^14 | 384 | 6 |
/// | dir{P3} | 38 | Amm2 | C2v^14 | 384 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 384 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 384 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 384 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 384 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 384 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 384 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 384 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 384 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 512 | 6 |
/// | dir{P1} | 81 | P-4 | S4^1 | 384 | 6 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 384 | 6 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 384 | 6 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 384 | 6 |
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 384 | 6 |
/// | dir{C1} | 31 | Pmn2_1 | C2v^7 | 384 | 6 |
/// | dir{C1} | 29 | Pca2_1 | C2v^5 | 384 | 6 |
/// | dir{P3} | 26 | Pmc2_1 | C2v^2 | 384 | 6 |
/// | dir{C1} | 21 | C222 | D2^6 | 384 | 6 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 384 | 6 |
/// | dir{C1} | 14 | P2_1/c | C2h^5 | 384 | 6 |
/// | dir{C1} | 14 | P2_1/c | C2h^5 | 384 | 6 |
/// | dir{C1} | 14 | P2_1/c | C2h^5 | 384 | 6 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 384 | 6 |
/// | dir{C10} | 13 | P2/c | C2h^4 | 384 | 6 |
/// | dir{C13} | 13 | P2/c | C2h^4 | 384 | 6 |
/// | dir{C1} | 11 | P2_1/m | C2h^2 | 384 | 6 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 384 | 6 |
/// | dir{C12} | 10 | P2/m | C2h^1 | 384 | 6 |
/// | dir{C1} | 19 | P2_12_12_1 | D2^4 | 384 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 384 | 6 |
/// | dir{P1} | 9 | Cc | Cs^4 | 768 | 6 |
/// | dir{P1} | 8 | Cm | Cs^3 | 768 | 6 |
/// | dir{P1} | 7 | Pc | Cs^2 | 768 | 6 |
/// | dir{P1} | 7 | Pc | Cs^2 | 768 | 6 |
/// | dir{P1} | 6 | Pm | Cs^1 | 768 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 768 | 6 |
/// | dir{P1} | 5 | C2 | C2^3 | 768 | 6 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 768 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 768 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 768 | 6 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 768 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 768 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 1536 | 6 |
///
pub struct Sg225;

/// # 226 Fm-3c (Oₕ⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{11}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{11}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{11}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{11}t6$ | 2 | B12a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{11}t9$ | 3 | C24c | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{11}t10$ | 3 | C48a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{11}t7$ | 3 | C24b | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{11}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 226 | Fm-3c | Oh^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 209 | F432 | O^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 202 | Fm-3 | Th^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 219 | F-43c | Td^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 3 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 6 | 1 |
/// | dir{P1} | 97 | I422 | D4^9 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 8 | 1 |
/// | dir{P1} | 108 | I4cm | C4v^10 | 6 | 1 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 4 | 1 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 6 | 1 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 12 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{4}^+$ | $k_{10}t8$ | 3 | C24b | no |
/// | X1- | $X_{4}^-$ | $k_{10}t7$ | 3 | C48a | no |
/// | X2+ | $X_{2}^+$ | $k_{10}t4$ | 3 | C24c | no |
/// | X2- | $X_{2}^-$ | $k_{10}t3$ | 3 | C48a | no |
/// | X3+ | $X_{3}^+$ | $k_{10}t6$ | 3 | C24c | no |
/// | X3- | $X_{3}^-$ | $k_{10}t5$ | 3 | C48a | no |
/// | X4+ | $X_{1}^+$ | $k_{10}t2$ | 3 | C24b | no |
/// | X4- | $X_{1}^-$ | $k_{10}t1$ | 3 | C48a | no |
/// | X5+ | $X_{5}^+$ | $k_{10}t10$ | 6 | E96f | no |
/// | X5- | $X_{5}^-$ | $k_{10}t9$ | 6 | E192f | no |
///
/// ### Isotropy subgroups for X1+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 6 | 1 |
/// | dir{P1} | 223 | Pm-3n | Oh^3 | 4 | 3 |
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 12 | 3 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X1- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 133 | P4_2/nbc | D4h^11 | 6 | 1 |
/// | dir{P1} | 208 | P4_232 | O^2 | 8 | 3 |
/// | dir{P1} | 133 | P4_2/nbc | D4h^11 | 12 | 3 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 24 | 3 |
/// | dir{P1} | 50 | Pban | D2h^4 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X2+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 127 | P4/mbm | D4h^5 | 6 | 1 |
/// | dir{P1} | 200 | Pm-3 | Th^1 | 8 | 3 |
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 12 | 3 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 24 | 3 |
///
/// ### Isotropy subgroups for X2- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 130 | P4/ncc | D4h^8 | 6 | 1 |
/// | dir{P1} | 218 | P-43n | Td^4 | 8 | 3 |
/// | dir{P1} | 133 | P4_2/nbc | D4h^11 | 12 | 3 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 24 | 3 |
/// | dir{P1} | 50 | Pban | D2h^4 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X3+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 135 | P4_2/mbc | D4h^13 | 6 | 1 |
/// | dir{P1} | 201 | Pn-3 | Th^2 | 8 | 3 |
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 12 | 3 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 24 | 3 |
///
/// ### Isotropy subgroups for X3- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 138 | P4_2/ncm | D4h^16 | 6 | 1 |
/// | dir{P1} | 218 | P-43n | Td^4 | 8 | 3 |
/// | dir{P1} | 124 | P4/mcc | D4h^2 | 12 | 3 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 24 | 3 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X4+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 124 | P4/mcc | D4h^2 | 6 | 1 |
/// | dir{P1} | 222 | Pn-3n | Oh^2 | 4 | 3 |
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 12 | 3 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 24 | 3 |
///
/// ### Isotropy subgroups for X4- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 6 | 1 |
/// | dir{P1} | 207 | P432 | O^1 | 8 | 3 |
/// | dir{P1} | 124 | P4/mcc | D4h^2 | 12 | 3 |
/// | dir{P1} | 89 | P422 | D4^1 | 24 | 3 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 24 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
///
/// ### Isotropy subgroups for X5+ ($X_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 1 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 12 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 24 | 1 |
/// | dir{C1} | 205 | Pa-3 | Th^6 | 8 | 3 |
/// | dir{P3} | 167 | R-3c | D3d^6 | 16 | 3 |
/// | dir{P3} | 135 | P4_2/mbc | D4h^13 | 12 | 3 |
/// | dir{P1} | 130 | P4/ncc | D4h^8 | 12 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P1} | 56 | Pccn | D2h^10 | 24 | 3 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 24 | 3 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 24 | 3 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{C1} | 14 | P2_1/c | C2h^5 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
///
/// ### Isotropy subgroups for X5- ($X_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 1 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 12 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 1 |
/// | dir{C1} | 198 | P2_13 | T^4 | 16 | 3 |
/// | dir{P1} | 161 | R3c | C3v^6 | 32 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P3} | 137 | P4_2/nmc | D4h^15 | 12 | 3 |
/// | dir{P1} | 128 | P4/mnc | D4h^6 | 12 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 24 | 3 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 24 | 3 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 24 | 3 |
/// | dir{C1} | 41 | Aba2 | C2v^17 | 48 | 3 |
/// | dir{C1} | 40 | Ama2 | C2v^16 | 48 | 3 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 48 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
/// | dir{P3} | 14 | P2_1/c | C2h^5 | 48 | 3 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 48 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 7 | Pc | Cs^2 | 96 | 3 |
/// | dir{P3} | 6 | Pm | Cs^1 | 96 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 96 | 3 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{L}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1L2 | $L_{1}L2$ | $k_{9}t2t3$ | 16 | H384a | no |
/// | L3 | $L_{3}$ | $k_{9}t1$ | 8 | F384b | yes |
///
/// ### Isotropy subgroups for L1L2 ($L_{1}L2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 1 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 48 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 48 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 96 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 64 | 4 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 64 | 4 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 48 | 4 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 48 | 4 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 48 | 4 |
/// | dir{P1} | 97 | I422 | D4^9 | 48 | 4 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 48 | 4 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 48 | 4 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 48 | 4 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 48 | 4 |
/// | dir{P1} | 146 | R3 | C3^4 | 128 | 4 |
/// | dir{P1} | 82 | I-4 | S4^2 | 96 | 4 |
/// | dir{P1} | 82 | I-4 | S4^2 | 96 | 4 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 96 | 4 |
/// | dir{P1} | 79 | I4 | C4^5 | 96 | 4 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 96 | 4 |
/// | dir{C1} | 42 | Fmm2 | C2v^18 | 96 | 4 |
/// | dir{C1} | 24 | I2_12_12_1 | D2^9 | 96 | 4 |
/// | dir{C1} | 24 | I2_12_12_1 | D2^9 | 96 | 4 |
/// | dir{C1} | 23 | I222 | D2^8 | 96 | 4 |
/// | dir{C1} | 23 | I222 | D2^8 | 96 | 4 |
/// | dir{C1} | 22 | F222 | D2^7 | 96 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 96 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 4 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 96 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 4 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 96 | 4 |
/// | dir{P1} | 9 | Cc | Cs^4 | 192 | 4 |
/// | dir{P1} | 8 | Cm | Cs^3 | 192 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
///
/// ### Isotropy subgroups for L3 ($L_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 16 | 1 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 16 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 1 |
/// | dir{P1} | 36 | Cmc2_1 | C2v^12 | 48 | 2 |
/// | dir{P1} | 21 | C222 | D2^6 | 48 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 96 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{P1} | 210 | F4_132 | O^4 | 16 | 4 |
/// | dir{P1} | 209 | F432 | O^3 | 16 | 4 |
/// | dir{P1} | 203 | Fd-3 | Th^4 | 16 | 4 |
/// | dir{P1} | 202 | Fm-3 | Th^3 | 16 | 4 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 48 | 4 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 48 | 4 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 48 | 4 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 48 | 4 |
/// | dir{P1} | 196 | F23 | T^2 | 32 | 4 |
/// | dir{P1} | 196 | F23 | T^2 | 32 | 4 |
/// | dir{P1} | 155 | R32 | D3^7 | 64 | 4 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 64 | 4 |
/// | dir{P3} | 82 | I-4 | S4^2 | 96 | 4 |
/// | dir{C1} | 82 | I-4 | S4^2 | 96 | 4 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 96 | 4 |
/// | dir{P1} | 79 | I4 | C4^5 | 96 | 4 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 96 | 4 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 96 | 4 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 96 | 4 |
/// | dir{P3} | 24 | I2_12_12_1 | D2^9 | 96 | 4 |
/// | dir{C1} | 23 | I222 | D2^8 | 96 | 4 |
/// | dir{P1} | 23 | I222 | D2^8 | 96 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 4 |
/// | dir{P1} | 146 | R3 | C3^4 | 128 | 4 |
/// | dir{P1} | 9 | Cc | Cs^4 | 192 | 4 |
/// | dir{P1} | 8 | Cm | Cs^3 | 192 | 4 |
/// | dir{P3} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1W2 | $W_{1}W4$ | $k_{8}t2t4$ | 12 | G1536b | no |
/// | W3W4 | $W_{2}W3$ | $k_{8}t1t3$ | 12 | G1536b | no |
/// | W5 | $W_{5}$ | $k_{8}t5$ | 12 | G1536c | yes |
///
/// ### Isotropy subgroups for W1W2 ($W_{1}W4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1} | 119 | I-4m2 | D2d^9 | 24 | 2 |
/// | dir{P1} | 74 | Imma | D2h^28 | 24 | 2 |
/// | dir{P1} | 71 | Immm | D2h^25 | 24 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 48 | 2 |
/// | dir{P3} | 131 | P4_2/mmc | D4h^9 | 48 | 4 |
/// | dir{P1} | 128 | P4/mnc | D4h^6 | 48 | 4 |
/// | dir{C3} | 59 | Pmmn | D2h^13 | 96 | 4 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 96 | 4 |
/// | dir{C1} | 47 | Pmmm | D2h^1 | 96 | 4 |
/// | dir{C2} | 40 | Ama2 | C2v^16 | 192 | 4 |
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 192 | 4 |
/// | dir{C1} | 25 | Pmm2 | C2v^1 | 192 | 4 |
/// | dir{P1} | 6 | Pm | Cs^1 | 384 | 4 |
/// | dir{C1} | 200 | Pm-3 | Th^1 | 64 | 6 |
/// | dir{C1} | 198 | P2_13 | T^4 | 128 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 256 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 512 | 6 |
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 192 | 6 |
/// | dir{C1} | 20 | C222_1 | D2^5 | 384 | 6 |
/// | dir{C1} | 47 | Pmmm | D2h^1 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 384 | 6 |
/// | dir{C11} | 11 | P2_1/m | C2h^2 | 384 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 768 | 6 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 384 | 6 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 768 | 6 |
/// | dir{P1} | 6 | Pm | Cs^1 | 768 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 1536 | 6 |
///
/// ### Isotropy subgroups for W3W4 ($W_{2}W3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1} | 120 | I-4c2 | D2d^10 | 24 | 2 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 24 | 2 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 24 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 48 | 2 |
/// | dir{P3} | 133 | P4_2/nbc | D4h^11 | 48 | 4 |
/// | dir{P1} | 130 | P4/ncc | D4h^8 | 48 | 4 |
/// | dir{C3} | 60 | Pbcn | D2h^14 | 96 | 4 |
/// | dir{P1} | 56 | Pccn | D2h^10 | 96 | 4 |
/// | dir{C1} | 50 | Pban | D2h^4 | 96 | 4 |
/// | dir{C2} | 41 | Aba2 | C2v^17 | 192 | 4 |
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 192 | 4 |
/// | dir{C1} | 30 | Pnc2 | C2v^6 | 192 | 4 |
/// | dir{P1} | 7 | Pc | Cs^2 | 384 | 4 |
/// | dir{C1} | 201 | Pn-3 | Th^2 | 64 | 6 |
/// | dir{C1} | 198 | P2_13 | T^4 | 128 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 256 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 512 | 6 |
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 192 | 6 |
/// | dir{C1} | 20 | C222_1 | D2^5 | 384 | 6 |
/// | dir{C1} | 48 | Pnnn | D2h^2 | 192 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 384 | 6 |
/// | dir{C11} | 14 | P2_1/c | C2h^5 | 384 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 768 | 6 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 384 | 6 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 768 | 6 |
/// | dir{P1} | 7 | Pc | Cs^2 | 768 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 1536 | 6 |
///
/// ### Isotropy subgroups for W5 ($W_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 98 | I4_122 | D4^10 | 24 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 24 | 2 |
/// | dir{P3} | 88 | I4_1/a | C4h^6 | 24 | 2 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 24 | 2 |
/// | dir{P1} | 74 | Imma | D2h^28 | 24 | 2 |
/// | dir{C1} | 72 | Ibam | D2h^26 | 24 | 2 |
/// | dir{C1} | 82 | I-4 | S4^2 | 48 | 2 |
/// | dir{P1} | 80 | I4_1 | C4^6 | 48 | 2 |
/// | dir{P3} | 79 | I4 | C4^5 | 48 | 2 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 48 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 137 | P4_2/nmc | D4h^15 | 48 | 4 |
/// | dir{P1} | 135 | P4_2/mbc | D4h^13 | 48 | 4 |
/// | dir{C1} | 126 | P4/nnc | D4h^4 | 48 | 4 |
/// | dir{P3} | 124 | P4/mcc | D4h^2 | 48 | 4 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 96 | 4 |
/// | dir{P3} | 112 | P-42c | D2d^2 | 96 | 4 |
/// | dir{P1} | 94 | P4_22_12 | D4^6 | 96 | 4 |
/// | dir{P1} | 89 | P422 | D4^1 | 96 | 4 |
/// | dir{P1} | 60 | Pbcn | D2h^14 | 96 | 4 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 96 | 4 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 96 | 4 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 96 | 4 |
/// | dir{P1} | 57 | Pbcm | D2h^11 | 96 | 4 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 96 | 4 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 96 | 4 |
/// | dir{P3} | 50 | Pban | D2h^4 | 96 | 4 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 96 | 4 |
/// | dir{C1} | 48 | Pnnn | D2h^2 | 96 | 4 |
/// | dir{P1} | 41 | Aba2 | C2v^17 | 192 | 4 |
/// | dir{P1} | 40 | Ama2 | C2v^16 | 192 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 192 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 192 | 4 |
/// | dir{P3} | 34 | Pnn2 | C2v^10 | 192 | 4 |
/// | dir{C1} | 32 | Pba2 | C2v^8 | 192 | 4 |
/// | dir{C1} | 31 | Pmn2_1 | C2v^7 | 192 | 4 |
/// | dir{P1} | 30 | Pnc2 | C2v^6 | 192 | 4 |
/// | dir{C11} | 29 | Pca2_1 | C2v^5 | 192 | 4 |
/// | dir{P1} | 28 | Pma2 | C2v^4 | 192 | 4 |
/// | dir{C1} | 26 | Pmc2_1 | C2v^2 | 192 | 4 |
/// | dir{C12} | 26 | Pmc2_1 | C2v^2 | 192 | 4 |
/// | dir{P3} | 14 | P2_1/c | C2h^5 | 192 | 4 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 192 | 4 |
/// | dir{C1} | 14 | P2_1/c | C2h^5 | 192 | 4 |
/// | dir{4D1} | 13 | P2/c | C2h^4 | 192 | 4 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 192 | 4 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 192 | 4 |
/// | dir{C1} | 11 | P2_1/m | C2h^2 | 192 | 4 |
/// | dir{C1} | 10 | P2/m | C2h^1 | 192 | 4 |
/// | dir{C3} | 18 | P2_12_12 | D2^3 | 192 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 192 | 4 |
/// | dir{C1} | 16 | P222 | D2^1 | 192 | 4 |
/// | dir{C2} | 9 | Cc | Cs^4 | 384 | 4 |
/// | dir{P1} | 7 | Pc | Cs^2 | 384 | 4 |
/// | dir{4D1} | 7 | Pc | Cs^2 | 384 | 4 |
/// | dir{P1} | 7 | Pc | Cs^2 | 384 | 4 |
/// | dir{C11} | 7 | Pc | Cs^2 | 384 | 4 |
/// | dir{C8} | 7 | Pc | Cs^2 | 384 | 4 |
/// | dir{P1} | 6 | Pm | Cs^1 | 384 | 4 |
/// | dir{P1} | 6 | Pm | Cs^1 | 384 | 4 |
/// | dir{P12} | 5 | C2 | C2^3 | 384 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 384 | 4 |
/// | dir{C9} | 2 | P-1 | Ci^1 | 384 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 384 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 384 | 4 |
/// | dir{C1} | 3 | P2 | C2^1 | 384 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 4 |
/// | dir{P3} | 213 | P4_132 | O^7 | 64 | 6 |
/// | dir{P1} | 212 | P4_332 | O^6 | 64 | 6 |
/// | dir{P1} | 208 | P4_232 | O^2 | 64 | 6 |
/// | dir{P1} | 207 | P432 | O^1 | 64 | 6 |
/// | dir{P1} | 205 | Pa-3 | Th^6 | 64 | 6 |
/// | dir{P1} | 205 | Pa-3 | Th^6 | 64 | 6 |
/// | dir{C1} | 198 | P2_13 | T^4 | 128 | 6 |
/// | dir{P1} | 198 | P2_13 | T^4 | 128 | 6 |
/// | dir{C1} | 195 | P23 | T^1 | 128 | 6 |
/// | dir{P3} | 155 | R32 | D3^7 | 256 | 6 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 256 | 6 |
/// | dir{P3} | 86 | P4_2/n | C4h^4 | 192 | 6 |
/// | dir{P1} | 85 | P4/n | C4h^3 | 192 | 6 |
/// | dir{P1} | 84 | P4_2/m | C4h^2 | 192 | 6 |
/// | dir{C1} | 83 | P4/m | C4h^1 | 192 | 6 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 192 | 6 |
/// | dir{P1} | 93 | P4_222 | D4^5 | 192 | 6 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 192 | 6 |
/// | dir{P1} | 89 | P422 | D4^1 | 192 | 6 |
/// | dir{P3} | 62 | Pnma | D2h^16 | 192 | 6 |
/// | dir{4D1} | 62 | Pnma | D2h^16 | 192 | 6 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 192 | 6 |
/// | dir{C1} | 61 | Pbca | D2h^15 | 192 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 512 | 6 |
/// | dir{P1} | 81 | P-4 | S4^1 | 384 | 6 |
/// | dir{P1} | 76 | P4_1 | C4^2 | 384 | 6 |
/// | dir{C1} | 77 | P4_2 | C4^3 | 384 | 6 |
/// | dir{P3} | 75 | P4 | C4^1 | 384 | 6 |
/// | dir{P1} | 33 | Pna2_1 | C2v^9 | 384 | 6 |
/// | dir{P1} | 31 | Pmn2_1 | C2v^7 | 384 | 6 |
/// | dir{P1} | 29 | Pca2_1 | C2v^5 | 384 | 6 |
/// | dir{P1} | 26 | Pmc2_1 | C2v^2 | 384 | 6 |
/// | dir{P3} | 21 | C222 | D2^6 | 384 | 6 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 384 | 6 |
/// | dir{P4} | 14 | P2_1/c | C2h^5 | 384 | 6 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 384 | 6 |
/// | dir{P5} | 14 | P2_1/c | C2h^5 | 384 | 6 |
/// | dir{C1} | 13 | P2/c | C2h^4 | 384 | 6 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 384 | 6 |
/// | dir{P3} | 11 | P2_1/m | C2h^2 | 384 | 6 |
/// | dir{P3} | 10 | P2/m | C2h^1 | 384 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 384 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 384 | 6 |
/// | dir{C1} | 16 | P222 | D2^1 | 384 | 6 |
/// | dir{C1} | 7 | Pc | Cs^2 | 768 | 6 |
/// | dir{P3} | 7 | Pc | Cs^2 | 768 | 6 |
/// | dir{P1} | 6 | Pm | Cs^1 | 768 | 6 |
/// | dir{C12} | 5 | C2 | C2^3 | 768 | 6 |
/// | dir{C1} | 2 | P-1 | Ci^1 | 768 | 6 |
/// | dir{C1} | 2 | P-1 | Ci^1 | 768 | 6 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 768 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 768 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 1536 | 6 |
///
pub struct Sg226;

/// # 227 Fd-3m (Oₕ⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{11}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{11}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{11}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{11}t6$ | 2 | B12a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{11}t9$ | 3 | C24c | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{11}t10$ | 3 | C48a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{11}t7$ | 3 | C24b | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{11}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 227 | Fd-3m | Oh^7 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 210 | F4_132 | O^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 203 | Fd-3 | Th^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 216 | F-43m | Td^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 3 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 6 | 1 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 8 | 1 |
/// | dir{P1} | 109 | I4_1md | C4v^11 | 6 | 1 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 4 | 1 |
/// | dir{P1} | 74 | Imma | D2h^28 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 6 | 1 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{10}t3$ | 6 | E192d | yes |
/// | X2 | $X_{3}$ | $k_{10}t4$ | 6 | E192e | yes |
/// | X3 | $X_{2}$ | $k_{10}t1$ | 6 | E192d | no |
/// | X4 | $X_{4}$ | $k_{10}t2$ | 6 | E192e | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 115 | P-4m2 | D2d^5 | 12 | 1 |
/// | dir{P1} | 51 | Pmma | D2h^5 | 12 | 1 |
/// | dir{P1} | 25 | Pmm2 | C2v^1 | 24 | 1 |
/// | dir{C1} | 215 | P-43m | Td^1 | 8 | 3 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 16 | 3 |
/// | dir{P3} | 91 | P4_122 | D4^3 | 24 | 3 |
/// | dir{P1} | 160 | R3m | C3v^5 | 32 | 3 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 24 | 3 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 48 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 48 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{C1} | 3 | P2 | C2^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 118 | P-4n2 | D2d^8 | 12 | 1 |
/// | dir{P1} | 52 | Pnna | D2h^6 | 12 | 1 |
/// | dir{P1} | 34 | Pnn2 | C2v^10 | 24 | 1 |
/// | dir{C1} | 195 | P23 | T^1 | 16 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{C1} | 111 | P-42m | D2d^1 | 24 | 3 |
/// | dir{P3} | 91 | P4_122 | D4^3 | 24 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 81 | P-4 | S4^1 | 48 | 3 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 48 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 3 | P2 | C2^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for X3 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 12 | 1 |
/// | dir{P1} | 53 | Pmna | D2h^7 | 12 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 24 | 1 |
/// | dir{C1} | 212 | P4_332 | O^6 | 8 | 3 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 16 | 3 |
/// | dir{P3} | 113 | P-42_1m | D2d^3 | 24 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 24 | 3 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 48 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for X4 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 12 | 1 |
/// | dir{P1} | 62 | Pnma | D2h^16 | 12 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 1 |
/// | dir{C1} | 198 | P2_13 | T^4 | 16 | 3 |
/// | dir{P3} | 160 | R3m | C3v^5 | 32 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P3} | 113 | P-42_1m | D2d^3 | 24 | 3 |
/// | dir{C1} | 92 | P4_12_12 | D4^4 | 24 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 76 | P4_1 | C4^2 | 48 | 3 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
/// | dir{C1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{L}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^+$ | $k_{9}t1$ | 4 | D384a | no |
/// | L1- | $L_{1}^-$ | $k_{9}t2$ | 4 | D384a | no |
/// | L2+ | $L_{2}^+$ | $k_{9}t3$ | 4 | D384a | no |
/// | L2- | $L_{2}^-$ | $k_{9}t4$ | 4 | D384a | no |
/// | L3+ | $L_{3}^+$ | $k_{9}t5$ | 8 | F384a | no |
/// | L3- | $L_{3}^-$ | $k_{9}t6$ | 8 | F384a | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 8 | 1 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 24 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 48 | 2 |
/// | dir{P1} | 216 | F-43m | Td^2 | 16 | 4 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 32 | 4 |
/// | dir{P1} | 160 | R3m | C3v^5 | 64 | 4 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 96 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 4 |
/// | dir{C1} | 8 | Cm | Cs^3 | 192 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 8 | 1 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 24 | 2 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 48 | 2 |
/// | dir{P1} | 219 | F-43c | Td^5 | 16 | 4 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 32 | 4 |
/// | dir{P1} | 161 | R3c | C3v^6 | 64 | 4 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 96 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 4 |
/// | dir{C1} | 9 | Cc | Cs^4 | 192 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
///
/// ### Isotropy subgroups for L2+ ($L_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 8 | 1 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 24 | 2 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 48 | 2 |
/// | dir{P1} | 219 | F-43c | Td^5 | 16 | 4 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 32 | 4 |
/// | dir{P1} | 161 | R3c | C3v^6 | 64 | 4 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 96 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 4 |
/// | dir{C1} | 9 | Cc | Cs^4 | 192 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
///
/// ### Isotropy subgroups for L2- ($L_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 8 | 1 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 24 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 48 | 2 |
/// | dir{P1} | 216 | F-43m | Td^2 | 16 | 4 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 32 | 4 |
/// | dir{P1} | 160 | R3m | C3v^5 | 64 | 4 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 96 | 4 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 4 |
/// | dir{C1} | 8 | Cm | Cs^3 | 192 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
///
/// ### Isotropy subgroups for L3+ ($L_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 24 | 2 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 24 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 48 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 48 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 32 | 4 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 32 | 4 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 48 | 4 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 48 | 4 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 64 | 4 |
/// | dir{P1} | 82 | I-4 | S4^2 | 96 | 4 |
/// | dir{C1} | 46 | Ima2 | C2v^22 | 96 | 4 |
/// | dir{C1} | 45 | Iba2 | C2v^21 | 96 | 4 |
/// | dir{P3} | 44 | Imm2 | C2v^20 | 96 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 96 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 4 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 96 | 4 |
/// | dir{P1} | 9 | Cc | Cs^4 | 192 | 4 |
/// | dir{P3} | 8 | Cm | Cs^3 | 192 | 4 |
/// | dir{P3} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
///
/// ### Isotropy subgroups for L3- ($L_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 24 | 2 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 24 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 48 | 2 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 48 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 48 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 32 | 4 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 32 | 4 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 48 | 4 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 48 | 4 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 64 | 4 |
/// | dir{P1} | 82 | I-4 | S4^2 | 96 | 4 |
/// | dir{C1} | 46 | Ima2 | C2v^22 | 96 | 4 |
/// | dir{C1} | 45 | Iba2 | C2v^21 | 96 | 4 |
/// | dir{P3} | 44 | Imm2 | C2v^20 | 96 | 4 |
/// | dir{P1} | 22 | F222 | D2^7 | 96 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 4 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 96 | 4 |
/// | dir{P1} | 9 | Cc | Cs^4 | 192 | 4 |
/// | dir{P3} | 8 | Cm | Cs^3 | 192 | 4 |
/// | dir{P3} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1 | $W_{1}$ | $k_{8}t1$ | 12 | G1536d | yes |
/// | W2 | $W_{2}$ | $k_{8}t2$ | 12 | G1536d | yes |
///
/// ### Isotropy subgroups for W1 ($W_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 24 | 2 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 24 | 2 |
/// | dir{P3} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 96 | 4 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 96 | 4 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 96 | 4 |
/// | dir{C1} | 95 | P4_322 | D4^7 | 96 | 4 |
/// | dir{P1} | 81 | P-4 | S4^1 | 192 | 4 |
/// | dir{P1} | 76 | P4_1 | C4^2 | 192 | 4 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 192 | 4 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 192 | 4 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 192 | 4 |
/// | dir{P3} | 18 | P2_12_12 | D2^3 | 192 | 4 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 192 | 4 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 192 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 4 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 192 | 4 |
/// | dir{C1} | 8 | Cm | Cs^3 | 384 | 4 |
/// | dir{P3} | 5 | C2 | C2^3 | 384 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 384 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 384 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 384 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 4 |
/// | dir{P1} | 218 | P-43n | Td^4 | 64 | 6 |
/// | dir{P1} | 215 | P-43m | Td^1 | 64 | 6 |
/// | dir{P1} | 198 | P2_13 | T^4 | 128 | 6 |
/// | dir{P3} | 167 | R-3c | D3d^6 | 128 | 6 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 128 | 6 |
/// | dir{C1} | 166 | R-3m | D3d^5 | 128 | 6 |
/// | dir{C1} | 166 | R-3m | D3d^5 | 128 | 6 |
/// | dir{P1} | 161 | R3c | C3v^6 | 256 | 6 |
/// | dir{P1} | 160 | R3m | C3v^5 | 256 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 256 | 6 |
/// | dir{C1} | 155 | R32 | D3^7 | 256 | 6 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 256 | 6 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 256 | 6 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 192 | 6 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 192 | 6 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 192 | 6 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 192 | 6 |
/// | dir{C1} | 37 | Ccc2 | C2v^13 | 384 | 6 |
/// | dir{P3} | 36 | Cmc2_1 | C2v^12 | 384 | 6 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 384 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 384 | 6 |
/// | dir{P3} | 18 | P2_12_12 | D2^3 | 384 | 6 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 384 | 6 |
/// | dir{C1} | 17 | P222_1 | D2^2 | 384 | 6 |
/// | dir{C1} | 16 | P222 | D2^1 | 384 | 6 |
/// | dir{C1} | 16 | P222 | D2^1 | 384 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 384 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 384 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 384 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 384 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 512 | 6 |
/// | dir{P3} | 81 | P-4 | S4^1 | 384 | 6 |
/// | dir{P1} | 9 | Cc | Cs^4 | 768 | 6 |
/// | dir{P3} | 8 | Cm | Cs^3 | 768 | 6 |
/// | dir{P1} | 5 | C2 | C2^3 | 768 | 6 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 768 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 768 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 768 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 768 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 1536 | 6 |
///
/// ### Isotropy subgroups for W2 ($W_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 24 | 2 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 24 | 2 |
/// | dir{P3} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 96 | 4 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 96 | 4 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 96 | 4 |
/// | dir{C1} | 95 | P4_322 | D4^7 | 96 | 4 |
/// | dir{P1} | 81 | P-4 | S4^1 | 192 | 4 |
/// | dir{P1} | 76 | P4_1 | C4^2 | 192 | 4 |
/// | dir{P1} | 35 | Cmm2 | C2v^11 | 192 | 4 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 192 | 4 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 192 | 4 |
/// | dir{P3} | 18 | P2_12_12 | D2^3 | 192 | 4 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 192 | 4 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 192 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 4 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 192 | 4 |
/// | dir{C1} | 8 | Cm | Cs^3 | 384 | 4 |
/// | dir{P3} | 5 | C2 | C2^3 | 384 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 384 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 384 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 384 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 4 |
/// | dir{P1} | 218 | P-43n | Td^4 | 64 | 6 |
/// | dir{P1} | 215 | P-43m | Td^1 | 64 | 6 |
/// | dir{P1} | 198 | P2_13 | T^4 | 128 | 6 |
/// | dir{C1} | 167 | R-3c | D3d^6 | 128 | 6 |
/// | dir{C1} | 167 | R-3c | D3d^6 | 128 | 6 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 128 | 6 |
/// | dir{P3} | 166 | R-3m | D3d^5 | 128 | 6 |
/// | dir{P1} | 161 | R3c | C3v^6 | 256 | 6 |
/// | dir{P1} | 160 | R3m | C3v^5 | 256 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 256 | 6 |
/// | dir{C1} | 155 | R32 | D3^7 | 256 | 6 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 256 | 6 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 256 | 6 |
/// | dir{P1} | 114 | P-42_1c | D2d^4 | 192 | 6 |
/// | dir{P1} | 113 | P-42_1m | D2d^3 | 192 | 6 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 192 | 6 |
/// | dir{P1} | 111 | P-42m | D2d^1 | 192 | 6 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 384 | 6 |
/// | dir{P3} | 36 | Cmc2_1 | C2v^12 | 384 | 6 |
/// | dir{C1} | 35 | Cmm2 | C2v^11 | 384 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 384 | 6 |
/// | dir{P3} | 18 | P2_12_12 | D2^3 | 384 | 6 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 384 | 6 |
/// | dir{C1} | 17 | P222_1 | D2^2 | 384 | 6 |
/// | dir{C1} | 16 | P222 | D2^1 | 384 | 6 |
/// | dir{C1} | 16 | P222 | D2^1 | 384 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 384 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 384 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 384 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 384 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 512 | 6 |
/// | dir{P3} | 81 | P-4 | S4^1 | 384 | 6 |
/// | dir{P3} | 9 | Cc | Cs^4 | 768 | 6 |
/// | dir{P1} | 8 | Cm | Cs^3 | 768 | 6 |
/// | dir{P1} | 5 | C2 | C2^3 | 768 | 6 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 768 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 768 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 768 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 768 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 1536 | 6 |
///
pub struct Sg227;

/// # 228 Fd-3c (Oₕ⁸)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{11}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{11}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{11}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{11}t6$ | 2 | B12a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{11}t9$ | 3 | C24c | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{11}t10$ | 3 | C48a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{11}t7$ | 3 | C24b | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{11}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 228 | Fd-3c | Oh^8 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 210 | F4_132 | O^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 203 | Fd-3 | Th^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 219 | F-43c | Td^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 3 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 6 | 1 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 6 | 1 |
/// | dir{P1} | 22 | F222 | D2^7 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 8 | 1 |
/// | dir{P1} | 110 | I4_1cd | C4v^12 | 6 | 1 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 4 | 1 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 6 | 1 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{4}$ | $k_{10}t4$ | 6 | E192d | yes |
/// | X2 | $X_{2}$ | $k_{10}t3$ | 6 | E192e | yes |
/// | X3 | $X_{3}$ | $k_{10}t1$ | 6 | E192d | no |
/// | X4 | $X_{1}$ | $k_{10}t2$ | 6 | E192e | no |
///
/// ### Isotropy subgroups for X1 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 12 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 12 | 1 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 24 | 1 |
/// | dir{C1} | 218 | P-43n | Td^4 | 8 | 3 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 16 | 3 |
/// | dir{P3} | 91 | P4_122 | D4^3 | 24 | 3 |
/// | dir{P1} | 161 | R3c | C3v^6 | 32 | 3 |
/// | dir{P1} | 112 | P-42c | D2d^2 | 24 | 3 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 48 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 48 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{C1} | 3 | P2 | C2^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 12 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 12 | 1 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 24 | 1 |
/// | dir{C1} | 195 | P23 | T^1 | 16 | 3 |
/// | dir{P3} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{C1} | 112 | P-42c | D2d^2 | 24 | 3 |
/// | dir{P3} | 91 | P4_122 | D4^3 | 24 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 81 | P-4 | S4^1 | 48 | 3 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 48 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{P1} | 16 | P222 | D2^1 | 48 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 3 | P2 | C2^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for X3 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 95 | P4_322 | D4^7 | 12 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 12 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 24 | 1 |
/// | dir{C1} | 212 | P4_332 | O^6 | 8 | 3 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 16 | 3 |
/// | dir{P3} | 114 | P-42_1c | D2d^4 | 24 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 24 | 3 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 48 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ### Isotropy subgroups for X4 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 96 | P4_32_12 | D4^8 | 12 | 1 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 12 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 1 |
/// | dir{C1} | 198 | P2_13 | T^4 | 16 | 3 |
/// | dir{P3} | 161 | R3c | C3v^6 | 32 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P3} | 114 | P-42_1c | D2d^4 | 24 | 3 |
/// | dir{C1} | 92 | P4_12_12 | D4^4 | 24 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 76 | P4_1 | C4^2 | 48 | 3 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 48 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 3 |
/// | dir{C1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 4 | P2_1 | C2^2 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
///
/// ## Irreps at $\mathrm{L}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1L2 | $L_{1}L2$ | $k_{9}t2t3$ | 16 | H384b | no |
/// | L3 | $L_{3}$ | $k_{9}t1$ | 8 | F384c | yes |
///
/// ### Isotropy subgroups for L1L2 ($L_{1}L2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 1 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 48 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 48 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 96 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{P1} | 155 | R32 | D3^7 | 64 | 4 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 64 | 4 |
/// | dir{P1} | 146 | R3 | C3^4 | 128 | 4 |
/// | dir{P1} | 82 | I-4 | S4^2 | 96 | 4 |
/// | dir{C1} | 22 | F222 | D2^7 | 96 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
///
/// ### Isotropy subgroups for L3 ($L_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 16 | 1 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 16 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 1 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 48 | 2 |
/// | dir{P1} | 20 | C222_1 | D2^5 | 48 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 96 | 2 |
/// | dir{P1} | 3 | P2 | C2^1 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{P1} | 196 | F23 | T^2 | 32 | 4 |
/// | dir{P1} | 155 | R32 | D3^7 | 64 | 4 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 64 | 4 |
/// | dir{P3} | 82 | I-4 | S4^2 | 96 | 4 |
/// | dir{P1} | 146 | R3 | C3^4 | 128 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 4 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 4 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1W2 | $W_{1}W2$ | $k_{8}t1t2$ | 1 | K1536a | no |
///
/// ### Isotropy subgroups for W1W2 ($W_{1}W2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1} | 82 | I-4 | S4^2 | 48 | 2 |
/// | dir{C1} | 24 | I2_12_12_1 | D2^9 | 48 | 2 |
/// | dir{C1} | 23 | I222 | D2^8 | 48 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{C1} | 114 | P-42_1c | D2d^4 | 96 | 4 |
/// | dir{C1} | 112 | P-42c | D2d^2 | 96 | 4 |
/// | dir{P1} | 95 | P4_322 | D4^7 | 96 | 4 |
/// | dir{C1} | 92 | P4_12_12 | D4^4 | 96 | 4 |
/// | dir{P1} | 81 | P-4 | S4^1 | 192 | 4 |
/// | dir{P3} | 78 | P4_3 | C4^4 | 192 | 4 |
/// | dir{P1} | 37 | Ccc2 | C2v^13 | 192 | 4 |
/// | dir{C1} | 20 | C222_1 | D2^5 | 192 | 4 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 192 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 192 | 4 |
/// | dir{P1} | 18 | P2_12_12 | D2^3 | 192 | 4 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 192 | 4 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 192 | 4 |
/// | dir{P1} | 16 | P222 | D2^1 | 192 | 4 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 192 | 4 |
/// | dir{C1} | 9 | Cc | Cs^4 | 384 | 4 |
/// | dir{P1} | 5 | C2 | C2^3 | 384 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 384 | 4 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 384 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 4 |
/// | dir{P1} | 3 | P2 | C2^1 | 384 | 4 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 384 | 4 |
/// | dir{P1} | 1 | P1 | C1^1 | 768 | 4 |
/// | dir{P1} | 198 | P2_13 | T^4 | 128 | 6 |
/// | dir{P1} | 198 | P2_13 | T^4 | 128 | 6 |
/// | dir{P1} | 195 | P23 | T^1 | 128 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 256 | 6 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 256 | 6 |
/// | dir{C1} | 81 | P-4 | S4^1 | 384 | 6 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 384 | 6 |
/// | dir{P3} | 19 | P2_12_12_1 | D2^4 | 384 | 6 |
/// | dir{C1} | 18 | P2_12_12 | D2^3 | 384 | 6 |
/// | dir{P3} | 17 | P222_1 | D2^2 | 384 | 6 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 384 | 6 |
/// | dir{P1} | 16 | P222 | D2^1 | 384 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 512 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 768 | 6 |
/// | dir{P1} | 4 | P2_1 | C2^2 | 768 | 6 |
/// | dir{P1} | 3 | P2 | C2^1 | 768 | 6 |
/// | dir{C1} | 2 | P-1 | Ci^1 | 768 | 6 |
/// | dir{P3} | 1 | P1 | C1^1 | 1536 | 6 |
///
pub struct Sg228;

/// # 229 Im-3m (Oₕ⁹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{11}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{11}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{11}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{11}t6$ | 2 | B12a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{11}t9$ | 3 | C24c | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{11}t10$ | 3 | C48a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{11}t7$ | 3 | C24b | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{11}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 229 | Im-3m | Oh^9 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 211 | I432 | O^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 204 | Im-3 | Th^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 217 | I-43m | Td^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 3 | 1 |
/// | dir{P1} | 71 | Immm | D2h^25 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 121 | I-42m | D2d^11 | 6 | 1 |
/// | dir{P1} | 97 | I422 | D4^9 | 6 | 1 |
/// | dir{P1} | 23 | I222 | D2^8 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 160 | R3m | C3v^5 | 8 | 1 |
/// | dir{P1} | 107 | I4mm | C4v^9 | 6 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 12 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 4 | 1 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 6 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 6 | 1 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 12 | 1 |
/// | dir{P1} | 8 | Cm | Cs^3 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{H}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1+ | $H_{1}^+$ | $k_{12}t1$ | 1 | A2a | no |
/// | H1- | $H_{1}^-$ | $k_{12}t2$ | 1 | A2a | no |
/// | H2+ | $H_{2}^+$ | $k_{12}t3$ | 1 | A2a | no |
/// | H2- | $H_{2}^-$ | $k_{12}t4$ | 1 | A2a | no |
/// | H3+ | $H_{3}^+$ | $k_{12}t5$ | 2 | B12a | no |
/// | H3- | $H_{3}^-$ | $k_{12}t6$ | 2 | B12a | no |
/// | H4+ | $H_{4}^+$ | $k_{12}t9$ | 3 | C48a | no |
/// | H4- | $H_{4}^-$ | $k_{12}t10$ | 3 | C48a | no |
/// | H5+ | $H_{5}^+$ | $k_{12}t7$ | 3 | C48a | no |
/// | H5- | $H_{5}^-$ | $k_{12}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for H1+ ($H_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 221 | Pm-3m | Oh^1 | 2 | 1 |
///
/// ### Isotropy subgroups for H1- ($H_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 222 | Pn-3n | Oh^2 | 2 | 1 |
///
/// ### Isotropy subgroups for H2+ ($H_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 223 | Pm-3n | Oh^3 | 2 | 1 |
///
/// ### Isotropy subgroups for H2- ($H_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 224 | Pn-3m | Oh^4 | 2 | 1 |
///
/// ### Isotropy subgroups for H3+ ($H_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 131 | P4_2/mmc | D4h^9 | 6 | 1 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 6 | 1 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 12 | 1 |
///
/// ### Isotropy subgroups for H3- ($H_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 6 | 1 |
/// | dir{P1} | 126 | P4/nnc | D4h^4 | 6 | 1 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 12 | 1 |
///
/// ### Isotropy subgroups for H4+ ($H_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 8 | 1 |
/// | dir{P1} | 128 | P4/mnc | D4h^6 | 6 | 1 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
///
/// ### Isotropy subgroups for H4- ($H_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 8 | 1 |
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 6 | 1 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
///
/// ### Isotropy subgroups for H5+ ($H_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 166 | R-3m | D3d^5 | 8 | 1 |
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 6 | 1 |
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 24 | 1 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
///
/// ### Isotropy subgroups for H5- ($H_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 8 | 1 |
/// | dir{P1} | 137 | P4_2/nmc | D4h^15 | 6 | 1 |
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 24 | 1 |
/// | dir{P1} | 11 | P2_1/m | C2h^2 | 24 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{N}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1+ | $N_{1}^+$ | $k_{9}t1$ | 6 | E192d | no |
/// | N1- | $N_{1}^-$ | $k_{9}t2$ | 6 | E384c | no |
/// | N2+ | $N_{3}^+$ | $k_{9}t7$ | 6 | E192e | no |
/// | N2- | $N_{3}^-$ | $k_{9}t8$ | 6 | E384c | no |
/// | N3+ | $N_{4}^+$ | $k_{9}t5$ | 6 | E192d | no |
/// | N3- | $N_{4}^-$ | $k_{9}t6$ | 6 | E384c | no |
/// | N4+ | $N_{2}^+$ | $k_{9}t3$ | 6 | E192e | no |
/// | N4- | $N_{2}^-$ | $k_{9}t4$ | 6 | E384c | no |
///
/// ### Isotropy subgroups for N1+ ($N_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 65 | Cmmm | D2h^19 | 12 | 1 |
/// | dir{P1} | 123 | P4/mmm | D4h^1 | 12 | 2 |
/// | dir{P1} | 47 | Pmmm | D2h^1 | 24 | 2 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 16 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 229 | Im-3m | Oh^9 | 8 | 6 |
/// | dir{P3} | 141 | I4_1/amd | D4h^19 | 24 | 6 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 32 | 6 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 24 | 6 |
/// | dir{P1} | 74 | Imma | D2h^28 | 48 | 6 |
/// | dir{P1} | 71 | Immm | D2h^25 | 48 | 6 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 48 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 6 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 96 | 6 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 96 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 6 |
///
/// ### Isotropy subgroups for N1- ($N_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 68 | Ccca | D2h^22 | 12 | 1 |
/// | dir{P1} | 125 | P4/nbm | D4h^3 | 12 | 2 |
/// | dir{P1} | 50 | Pban | D2h^4 | 24 | 2 |
/// | dir{P3} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
/// | dir{C1} | 211 | I432 | O^5 | 16 | 6 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 32 | 6 |
/// | dir{P3} | 142 | I4_1/acd | D4h^20 | 24 | 6 |
/// | dir{C1} | 140 | I4/mcm | D4h^18 | 24 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 64 | 6 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 48 | 6 |
/// | dir{P1} | 97 | I422 | D4^9 | 48 | 6 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 48 | 6 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 48 | 6 |
/// | dir{C1} | 43 | Fdd2 | C2v^19 | 96 | 6 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 96 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 6 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 96 | 6 |
/// | dir{P1} | 23 | I222 | D2^8 | 96 | 6 |
/// | dir{P1} | 22 | F222 | D2^7 | 96 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 6 |
/// | dir{P1} | 9 | Cc | Cs^4 | 192 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 6 |
/// | dir{P3} | 9 | Cc | Cs^4 | 192 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 6 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 192 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 6 |
///
/// ### Isotropy subgroups for N2+ ($N_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 1 |
/// | dir{P1} | 127 | P4/mbm | D4h^5 | 12 | 2 |
/// | dir{P1} | 55 | Pbam | D2h^9 | 24 | 2 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 204 | Im-3 | Th^5 | 16 | 6 |
/// | dir{P3} | 167 | R-3c | D3d^6 | 32 | 6 |
/// | dir{P3} | 141 | I4_1/amd | D4h^19 | 24 | 6 |
/// | dir{C1} | 139 | I4/mmm | D4h^17 | 24 | 6 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 64 | 6 |
/// | dir{P1} | 87 | I4/m | C4h^5 | 48 | 6 |
/// | dir{P1} | 74 | Imma | D2h^28 | 48 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 6 |
/// | dir{P1} | 71 | Immm | D2h^25 | 48 | 6 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 96 | 6 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 96 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 6 |
///
/// ### Isotropy subgroups for N2- ($N_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 1 |
/// | dir{P1} | 129 | P4/nmm | D4h^7 | 12 | 2 |
/// | dir{P1} | 59 | Pmmn | D2h^13 | 24 | 2 |
/// | dir{P1} | 160 | R3m | C3v^5 | 32 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
/// | dir{C1} | 217 | I-43m | Td^3 | 16 | 6 |
/// | dir{P3} | 166 | R-3m | D3d^5 | 32 | 6 |
/// | dir{P3} | 142 | I4_1/acd | D4h^20 | 24 | 6 |
/// | dir{C1} | 140 | I4/mcm | D4h^18 | 24 | 6 |
/// | dir{P1} | 160 | R3m | C3v^5 | 64 | 6 |
/// | dir{P1} | 121 | I-42m | D2d^11 | 48 | 6 |
/// | dir{P1} | 108 | I4cm | C4v^10 | 48 | 6 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 48 | 6 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 48 | 6 |
/// | dir{C1} | 43 | Fdd2 | C2v^19 | 96 | 6 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 96 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 6 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 96 | 6 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 96 | 6 |
/// | dir{P1} | 23 | I222 | D2^8 | 96 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 6 |
/// | dir{P1} | 9 | Cc | Cs^4 | 192 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 6 |
/// | dir{P3} | 9 | Cc | Cs^4 | 192 | 6 |
/// | dir{C1} | 8 | Cm | Cs^3 | 192 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 6 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 192 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 6 |
///
/// ### Isotropy subgroups for N3+ ($N_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 66 | Cccm | D2h^20 | 12 | 1 |
/// | dir{P1} | 134 | P4_2/nnm | D4h^12 | 12 | 2 |
/// | dir{P1} | 48 | Pnnn | D2h^2 | 24 | 2 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 16 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 230 | Ia-3d | Oh^10 | 8 | 6 |
/// | dir{P3} | 140 | I4/mcm | D4h^18 | 24 | 6 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 32 | 6 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 24 | 6 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 48 | 6 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 48 | 6 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 48 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 6 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 96 | 6 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 96 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 6 |
///
/// ### Isotropy subgroups for N3- ($N_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 67 | Cmma | D2h^21 | 12 | 1 |
/// | dir{P1} | 132 | P4_2/mcm | D4h^10 | 12 | 2 |
/// | dir{P1} | 49 | Pccm | D2h^3 | 24 | 2 |
/// | dir{P3} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
/// | dir{C1} | 214 | I4_132 | O^8 | 16 | 6 |
/// | dir{P1} | 166 | R-3m | D3d^5 | 32 | 6 |
/// | dir{C1} | 141 | I4_1/amd | D4h^19 | 24 | 6 |
/// | dir{P3} | 139 | I4/mmm | D4h^17 | 24 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 64 | 6 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 48 | 6 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 48 | 6 |
/// | dir{P1} | 74 | Imma | D2h^28 | 48 | 6 |
/// | dir{P1} | 71 | Immm | D2h^25 | 48 | 6 |
/// | dir{C1} | 42 | Fmm2 | C2v^18 | 96 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 6 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 96 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 6 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 96 | 6 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 96 | 6 |
/// | dir{P1} | 22 | F222 | D2^7 | 96 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 6 |
/// | dir{P1} | 8 | Cm | Cs^3 | 192 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 6 |
/// | dir{P3} | 8 | Cm | Cs^3 | 192 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 6 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 192 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 6 |
///
/// ### Isotropy subgroups for N4+ ($N_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 64 | Cmca | D2h^18 | 12 | 1 |
/// | dir{P1} | 138 | P4_2/ncm | D4h^16 | 12 | 2 |
/// | dir{P1} | 56 | Pccn | D2h^10 | 24 | 2 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{C1} | 206 | Ia-3 | Th^7 | 16 | 6 |
/// | dir{P3} | 166 | R-3m | D3d^5 | 32 | 6 |
/// | dir{C1} | 142 | I4_1/acd | D4h^20 | 24 | 6 |
/// | dir{P3} | 140 | I4/mcm | D4h^18 | 24 | 6 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 64 | 6 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 48 | 6 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 48 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 6 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 48 | 6 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 96 | 6 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 96 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 6 |
///
/// ### Isotropy subgroups for N4- ($N_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 63 | Cmcm | D2h^17 | 12 | 1 |
/// | dir{P1} | 136 | P4_2/mnm | D4h^14 | 12 | 2 |
/// | dir{P1} | 58 | Pnnm | D2h^12 | 24 | 2 |
/// | dir{P1} | 160 | R3m | C3v^5 | 32 | 3 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 3 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
/// | dir{C1} | 220 | I-43d | Td^6 | 16 | 6 |
/// | dir{P3} | 167 | R-3c | D3d^6 | 32 | 6 |
/// | dir{C1} | 141 | I4_1/amd | D4h^19 | 24 | 6 |
/// | dir{P3} | 139 | I4/mmm | D4h^17 | 24 | 6 |
/// | dir{P1} | 161 | R3c | C3v^6 | 64 | 6 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 48 | 6 |
/// | dir{P1} | 109 | I4_1md | C4v^11 | 48 | 6 |
/// | dir{P1} | 74 | Imma | D2h^28 | 48 | 6 |
/// | dir{P1} | 71 | Immm | D2h^25 | 48 | 6 |
/// | dir{C1} | 42 | Fmm2 | C2v^18 | 96 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 6 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 96 | 6 |
/// | dir{C1} | 12 | C2/m | C2h^3 | 96 | 6 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 96 | 6 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 96 | 6 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 96 | 6 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 96 | 6 |
/// | dir{P1} | 8 | Cm | Cs^3 | 192 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 6 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 192 | 6 |
/// | dir{C1} | 9 | Cc | Cs^4 | 192 | 6 |
/// | dir{P3} | 8 | Cm | Cs^3 | 192 | 6 |
/// | dir{C1} | 5 | C2 | C2^3 | 192 | 6 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 192 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 6 |
///
/// ## Irreps at $\mathrm{P}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1 | $P_{1}$ | $k_{10}t1$ | 2 | B8a | no |
/// | P2 | $P_{2}$ | $k_{10}t2$ | 2 | B8a | no |
/// | P3 | $P_{3}$ | $k_{10}t3$ | 4 | D48e | no |
/// | P4 | $P_{5}$ | $k_{10}t4$ | 6 | E192a | yes |
/// | P5 | $P_{4}$ | $k_{10}t5$ | 6 | E192a | yes |
///
/// ### Isotropy subgroups for P1 ($P_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 227 | Fd-3m | Oh^7 | 4 | 2 |
/// | dir{P1} | 225 | Fm-3m | Oh^5 | 4 | 2 |
/// | dir{P1} | 216 | F-43m | Td^2 | 8 | 2 |
///
/// ### Isotropy subgroups for P2 ($P_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 228 | Fd-3c | Oh^8 | 4 | 2 |
/// | dir{P1} | 226 | Fm-3c | Oh^6 | 4 | 2 |
/// | dir{P1} | 219 | F-43c | Td^5 | 8 | 2 |
///
/// ### Isotropy subgroups for P3 ($P_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 12 | 2 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 12 | 2 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 12 | 2 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 12 | 2 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 24 | 2 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 24 | 2 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 24 | 2 |
/// | dir{P1} | 97 | I422 | D4^9 | 24 | 2 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 24 | 2 |
/// | dir{P1} | 69 | Fmmm | D2h^23 | 24 | 2 |
/// | dir{P1} | 22 | F222 | D2^7 | 48 | 2 |
///
/// ### Isotropy subgroups for P4 ($P_{5}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P3} | 166 | R-3m | D3d^5 | 16 | 2 |
/// | dir{C1} | 166 | R-3m | D3d^5 | 16 | 2 |
/// | dir{P1} | 141 | I4_1/amd | D4h^19 | 12 | 2 |
/// | dir{P1} | 139 | I4/mmm | D4h^17 | 12 | 2 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 24 | 2 |
/// | dir{P3} | 121 | I-42m | D2d^11 | 24 | 2 |
/// | dir{C1} | 74 | Imma | D2h^28 | 24 | 2 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 24 | 2 |
/// | dir{P1} | 160 | R3m | C3v^5 | 32 | 2 |
/// | dir{P1} | 119 | I-4m2 | D2d^9 | 24 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 2 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 48 | 2 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 48 | 2 |
/// | dir{P1} | 44 | Imm2 | C2v^20 | 48 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 48 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 48 | 2 |
/// | dir{C1} | 24 | I2_12_12_1 | D2^9 | 48 | 2 |
/// | dir{C1} | 23 | I222 | D2^8 | 48 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{C1} | 8 | Cm | Cs^3 | 96 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 2 |
///
/// ### Isotropy subgroups for P5 ($P_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P3} | 167 | R-3c | D3d^6 | 16 | 2 |
/// | dir{C1} | 167 | R-3c | D3d^6 | 16 | 2 |
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 12 | 2 |
/// | dir{P1} | 140 | I4/mcm | D4h^18 | 12 | 2 |
/// | dir{P1} | 122 | I-42d | D2d^12 | 24 | 2 |
/// | dir{P3} | 121 | I-42m | D2d^11 | 24 | 2 |
/// | dir{C1} | 74 | Imma | D2h^28 | 24 | 2 |
/// | dir{P1} | 72 | Ibam | D2h^26 | 24 | 2 |
/// | dir{P1} | 161 | R3c | C3v^6 | 32 | 2 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 24 | 2 |
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 2 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 48 | 2 |
/// | dir{P1} | 46 | Ima2 | C2v^22 | 48 | 2 |
/// | dir{P1} | 45 | Iba2 | C2v^21 | 48 | 2 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 48 | 2 |
/// | dir{P1} | 42 | Fmm2 | C2v^18 | 48 | 2 |
/// | dir{C1} | 24 | I2_12_12_1 | D2^9 | 48 | 2 |
/// | dir{C1} | 23 | I222 | D2^8 | 48 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 12 | C2/m | C2h^3 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P1} | 8 | Cm | Cs^3 | 96 | 2 |
/// | dir{C1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{C1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P3} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 2 |
///
pub struct Sg229;

/// # 230 Ia-3d (Oₕ¹⁰)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{11}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{11}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{11}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{11}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{11}t5$ | 2 | B6a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{11}t6$ | 2 | B12a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{11}t9$ | 3 | C24c | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{11}t10$ | 3 | C48a | no |
/// | GM5+ | $\Gamma_{5}^+$ | $k_{11}t7$ | 3 | C24b | no |
/// | GM5- | $\Gamma_{5}^-$ | $k_{11}t8$ | 3 | C48a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 230 | Ia-3d | Oh^10 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 214 | I4_132 | O^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 206 | Ia-3 | Th^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 220 | I-43d | Td^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 142 | I4_1/acd | D4h^20 | 3 | 1 |
/// | dir{P1} | 73 | Ibca | D2h^27 | 6 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 122 | I-42d | D2d^12 | 6 | 1 |
/// | dir{P1} | 98 | I4_122 | D4^10 | 6 | 1 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 12 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 148 | R-3 | C3i^2 | 8 | 1 |
/// | dir{P1} | 88 | I4_1/a | C4h^6 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 161 | R3c | C3v^6 | 8 | 1 |
/// | dir{P1} | 110 | I4_1cd | C4v^12 | 6 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ### Isotropy subgroups for GM5+ ($\Gamma_{5}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 167 | R-3c | D3d^6 | 4 | 1 |
/// | dir{P1} | 70 | Fddd | D2h^24 | 6 | 1 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 12 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 24 | 1 |
///
/// ### Isotropy subgroups for GM5- ($\Gamma_{5}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 155 | R32 | D3^7 | 8 | 1 |
/// | dir{P1} | 120 | I-4c2 | D2d^10 | 6 | 1 |
/// | dir{P1} | 43 | Fdd2 | C2v^19 | 12 | 1 |
/// | dir{P1} | 9 | Cc | Cs^4 | 24 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 24 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 48 | 1 |
///
/// ## Irreps at $\mathrm{H}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | H1 | $H_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | H2H3 | $H_{2}H3$ | $k_{12}t2t3$ | 4 | D24e | no |
/// | H4 | $H_{4}$ | $k_{12}t4$ | 6 | E96d | no |
///
/// ### Isotropy subgroups for H1 ($H_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 213 | P4_132 | O^7 | 4 | 1 |
/// | dir{P1} | 205 | Pa-3 | Th^6 | 4 | 1 |
/// | dir{P1} | 198 | P2_13 | T^4 | 8 | 1 |
///
/// ### Isotropy subgroups for H2H3 ($H_{2}H3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 12 | 1 |
/// | dir{P1} | 61 | Pbca | D2h^15 | 12 | 1 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 24 | 1 |
///
/// ### Isotropy subgroups for H4 ($H_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1} | 155 | R32 | D3^7 | 16 | 1 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 16 | 1 |
/// | dir{P1} | 117 | P-4b2 | D2d^7 | 12 | 1 |
/// | dir{P1} | 116 | P-4c2 | D2d^6 | 12 | 1 |
/// | dir{P1} | 91 | P4_122 | D4^3 | 12 | 1 |
/// | dir{P1} | 54 | Pcca | D2h^8 | 12 | 1 |
/// | dir{P1} | 146 | R3 | C3^4 | 32 | 1 |
/// | dir{P1} | 81 | P-4 | S4^1 | 24 | 1 |
/// | dir{P1} | 32 | Pba2 | C2v^8 | 24 | 1 |
/// | dir{P1} | 27 | Pcc2 | C2v^3 | 24 | 1 |
/// | dir{C1} | 21 | C222 | D2^6 | 24 | 1 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 24 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 24 | 1 |
/// | dir{P1} | 7 | Pc | Cs^2 | 48 | 1 |
/// | dir{P1} | 5 | C2 | C2^3 | 48 | 1 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 48 | 1 |
/// | dir{P3} | 3 | P2 | C2^1 | 48 | 1 |
/// | dir{P1} | 1 | P1 | C1^1 | 96 | 1 |
///
/// ## Irreps at $\mathrm{N}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | N1 | $N_{2}$ | $k_{9}t2$ | 12 | G384b | yes |
/// | N2 | $N_{1}$ | $k_{9}t1$ | 12 | G384c | yes |
///
/// ### Isotropy subgroups for N1 ($N_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 20 | C222_1 | D2^5 | 24 | 1 |
/// | dir{P1} | 14 | P2_1/c | C2h^5 | 24 | 1 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 48 | 1 |
/// | dir{P3} | 114 | P-42_1c | D2d^4 | 24 | 2 |
/// | dir{P1} | 92 | P4_12_12 | D4^4 | 24 | 2 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 19 | P2_12_12_1 | D2^4 | 48 | 2 |
/// | dir{C1} | 18 | P2_12_12 | D2^3 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P3} | 4 | P2_1 | C2^2 | 96 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{P1} | 161 | R3c | C3v^6 | 32 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 32 | 3 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 64 | 6 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 64 | 6 |
/// | dir{C1} | 22 | F222 | D2^7 | 96 | 6 |
/// | dir{P1} | 82 | I-4 | S4^2 | 96 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 128 | 6 |
/// | dir{P3} | 5 | C2 | C2^3 | 192 | 6 |
/// | dir{P3} | 5 | C2 | C2^3 | 192 | 6 |
/// | dir{C1} | 2 | P-1 | Ci^1 | 192 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 6 |
///
/// ### Isotropy subgroups for N2 ($N_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 21 | C222 | D2^6 | 24 | 1 |
/// | dir{P1} | 13 | P2/c | C2h^4 | 24 | 1 |
/// | dir{P3} | 3 | P2 | C2^1 | 48 | 1 |
/// | dir{P3} | 112 | P-42c | D2d^2 | 24 | 2 |
/// | dir{P1} | 91 | P4_122 | D4^3 | 24 | 2 |
/// | dir{C1} | 15 | C2/c | C2h^6 | 48 | 2 |
/// | dir{P1} | 17 | P222_1 | D2^2 | 48 | 2 |
/// | dir{C1} | 16 | P222 | D2^1 | 48 | 2 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P3} | 3 | P2 | C2^1 | 96 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 2 |
/// | dir{P1} | 167 | R-3c | D3d^6 | 16 | 3 |
/// | dir{P1} | 155 | R32 | D3^7 | 32 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 15 | C2/c | C2h^6 | 48 | 3 |
/// | dir{P1} | 9 | Cc | Cs^4 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 2 | P-1 | Ci^1 | 96 | 3 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 3 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 3 |
/// | dir{P1} | 148 | R-3 | C3i^2 | 64 | 6 |
/// | dir{P1} | 155 | R32 | D3^7 | 64 | 6 |
/// | dir{P1} | 82 | I-4 | S4^2 | 96 | 6 |
/// | dir{P1} | 146 | R3 | C3^4 | 128 | 6 |
/// | dir{P1} | 22 | F222 | D2^7 | 96 | 6 |
/// | dir{P3} | 5 | C2 | C2^3 | 192 | 6 |
/// | dir{C1} | 2 | P-1 | Ci^1 | 192 | 6 |
/// | dir{P1} | 5 | C2 | C2^3 | 192 | 6 |
/// | dir{P1} | 1 | P1 | C1^1 | 384 | 6 |
///
/// ## Irreps at $\mathrm{P}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | P1P2 | $P_{1}P2$ | $k_{10}t1t2$ | 8 | F192f | no |
/// | P3 | $P_{3}$ | $k_{10}t3$ | 8 | F192e | no |
///
/// ### Isotropy subgroups for P1P2 ($P_{1}P2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1} | 82 | I-4 | S4^2 | 48 | 2 |
/// | dir{P1} | 24 | I2_12_12_1 | D2^9 | 48 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 2 |
///
/// ### Isotropy subgroups for P3 ($P_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1} | 155 | R32 | D3^7 | 32 | 2 |
/// | dir{P3} | 148 | R-3 | C3i^2 | 32 | 2 |
/// | dir{P3} | 82 | I-4 | S4^2 | 48 | 2 |
/// | dir{C1} | 24 | I2_12_12_1 | D2^9 | 48 | 2 |
/// | dir{P1} | 23 | I222 | D2^8 | 48 | 2 |
/// | dir{P1} | 146 | R3 | C3^4 | 64 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P1} | 5 | C2 | C2^3 | 96 | 2 |
/// | dir{P3} | 2 | P-1 | Ci^1 | 96 | 2 |
/// | dir{P1} | 1 | P1 | C1^1 | 192 | 2 |
///
pub struct Sg230;
