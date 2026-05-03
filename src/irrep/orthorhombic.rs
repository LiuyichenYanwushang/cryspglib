//! # Orthorhombic space groups (#16–#74)
//!
//! Orthorhombic crystals have three mutually perpendicular axes of
//! different lengths.  The Brillouin zone is an orthorhombic prism.
//!
//! ## Common k-point labels
//!
//! | Label | Coords (fractional) | Little group |
//! |-------|---------------------|--------------|
//! | Γ | (0, 0, 0) | mmm (D₂ₕ) |
//! | X | (½, 0, 0) | mmm (D₂ₕ) |
//! | Y | (0, ½, 0) | mmm (D₂ₕ) |
//! | Z | (0, 0, ½) | mmm (D₂ₕ) |
//! | U | (½, ½, 0) | mmm (D₂ₕ) |
//! | T | (0, ½, ½) | mmm (D₂ₕ) |
//! | S | (½, 0, ½) | mmm (D₂ₕ) |
//! | R | (½, ½, ½) | mmm (D₂ₕ) |
//!
//! ---


/// # 16 P222 (D₂¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{19}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 16 | P222 | D2^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{20}t1$ | 1 | A2a | no |
/// | X2 | $X_{3}$ | $k_{20}t4$ | 1 | A2a | no |
/// | X3 | $X_{4}$ | $k_{20}t2$ | 1 | A2a | no |
/// | X4 | $X_{2}$ | $k_{20}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{26}t1$ | 1 | A2a | no |
/// | R2 | $R_{3}$ | $k_{26}t4$ | 1 | A2a | no |
/// | R3 | $R_{4}$ | $k_{26}t2$ | 1 | A2a | no |
/// | R4 | $R_{2}$ | $k_{26}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 22 | F222 | D2^7 | 2 | 1 |
///
/// ### Isotropy subgroups for R2 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 22 | F222 | D2^7 | 2 | 1 |
///
/// ### Isotropy subgroups for R3 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 22 | F222 | D2^7 | 2 | 1 |
///
/// ### Isotropy subgroups for R4 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 22 | F222 | D2^7 | 2 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{23}t1$ | 1 | A2a | no |
/// | T2 | $T_{3}$ | $k_{23}t4$ | 1 | A2a | no |
/// | T3 | $T_{4}$ | $k_{23}t2$ | 1 | A2a | no |
/// | T4 | $T_{2}$ | $k_{23}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for T3 ($T_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for T4 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{22}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{22}t4$ | 1 | A2a | no |
/// | Z3 | $Z_{4}$ | $k_{22}t2$ | 1 | A2a | no |
/// | Z4 | $Z_{2}$ | $k_{22}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{25}t1$ | 1 | A2a | no |
/// | S2 | $S_{3}$ | $k_{25}t4$ | 1 | A2a | no |
/// | S3 | $S_{4}$ | $k_{25}t2$ | 1 | A2a | no |
/// | S4 | $S_{2}$ | $k_{25}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for S2 ($S_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for S3 ($S_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for S4 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $U_{1}$ | $k_{24}t1$ | 1 | A2a | no |
/// | U2 | $U_{3}$ | $k_{24}t4$ | 1 | A2a | no |
/// | U3 | $U_{4}$ | $k_{24}t2$ | 1 | A2a | no |
/// | U4 | $U_{2}$ | $k_{24}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for U1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for U2 ($U_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for U3 ($U_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for U4 ($U_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{21}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{3}$ | $k_{21}t4$ | 1 | A2a | no |
/// | Y3 | $Y_{4}$ | $k_{21}t2$ | 1 | A2a | no |
/// | Y4 | $Y_{2}$ | $k_{21}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 16 | P222 | D2^1 | 2 | 1 |
///
pub struct Sg16;

/// # 17 P222_1 (D₂²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{19}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{20}t1$ | 1 | A2a | no |
/// | X2 | $X_{3}$ | $k_{20}t4$ | 1 | A2a | no |
/// | X3 | $X_{4}$ | $k_{20}t2$ | 1 | A2a | no |
/// | X4 | $X_{2}$ | $k_{20}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{26}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{23}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{22}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,0)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{25}t1$ | 1 | A2a | no |
/// | S2 | $S_{3}$ | $k_{25}t4$ | 1 | A2a | no |
/// | S3 | $S_{4}$ | $k_{25}t2$ | 1 | A2a | no |
/// | S4 | $S_{2}$ | $k_{25}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ### Isotropy subgroups for S2 ($S_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ### Isotropy subgroups for S3 ($S_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ### Isotropy subgroups for S4 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $U_{1}$ | $k_{24}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for U1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{21}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{3}$ | $k_{21}t4$ | 1 | A2a | no |
/// | Y3 | $Y_{4}$ | $k_{21}t2$ | 1 | A2a | no |
/// | Y4 | $Y_{2}$ | $k_{21}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
pub struct Sg17;

/// # 18 P2_12_12 (D₂³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{19}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{26}t1t2$ | 2 | B4a | no |
/// | R3R4 | $R_{3}R4$ | $k_{26}t3t4$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ### Isotropy subgroups for R3R4 ($R_{3}R4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{23}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{22}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{22}t4$ | 1 | A2a | no |
/// | Z3 | $Z_{4}$ | $k_{22}t2$ | 1 | A2a | no |
/// | Z4 | $Z_{2}$ | $k_{22}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 19 | P2_12_12_1 | D2^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 19 | P2_12_12_1 | D2^4 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1S2 | $S_{1}S2$ | $k_{25}t1t2$ | 2 | B4a | no |
/// | S3S4 | $S_{3}S4$ | $k_{25}t3t4$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for S1S2 ($S_{1}S2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ### Isotropy subgroups for S3S4 ($S_{3}S4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $U_{1}$ | $k_{24}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for U1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
pub struct Sg18;

/// # 19 P2_12_12_1 (D₂⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{19}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 19 | P2_12_12_1 | D2^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $Y_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,0)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R1 | $R_{1}R1$ | $k_{26}t1t1$ | 4 | D8a | no |
///
/// ### Isotropy subgroups for R1R1 ($R_{1}R1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1T3 | $U_{3}U4$ | $k_{23}t1t2$ | 2 | B4a | no |
/// | T2T4 | $U_{1}U2$ | $k_{23}t3t4$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for T1T3 ($U_{3}U4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ### Isotropy subgroups for T2T4 ($U_{1}U2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{22}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,0)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1S2 | $S_{3}S4$ | $k_{25}t1t4$ | 2 | B4a | no |
/// | S3S4 | $S_{1}S2$ | $k_{25}t2t3$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for S1S2 ($S_{3}S4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ### Isotropy subgroups for S3S4 ($S_{1}S2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1U4 | $T_{3}T4$ | $k_{24}t3t4$ | 2 | B4a | no |
/// | U2U3 | $T_{1}T2$ | $k_{24}t1t2$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for U1U4 ($T_{3}T4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ### Isotropy subgroups for U2U3 ($T_{1}T2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $X_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,0)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
pub struct Sg19;

/// # 20 C222_1 (D₂⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{14}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{13}t1t2$ | 4 | D16c | yes |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{C12(2)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{17}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{16}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | S2 | $S_{2}$ | $k_{12}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,0)} | 20 | C222_1 | D2^5 | 4 | 2 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ### Isotropy subgroups for S2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,0)} | 20 | C222_1 | D2^5 | 4 | 2 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{3}$ | $k_{15}t4$ | 1 | A2a | no |
/// | Y3 | $Y_{4}$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y4 | $Y_{2}$ | $k_{15}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 19 | P2_12_12_1 | D2^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
pub struct Sg20;

/// # 21 C222 (D₂⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{14}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | R2 | $R_{2}$ | $k_{13}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 22 | F222 | D2^7 | 4 | 2 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 22 | F222 | D2^7 | 4 | 2 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | T2 | $T_{3}$ | $k_{17}t4$ | 1 | A2a | no |
/// | T3 | $T_{4}$ | $k_{17}t2$ | 1 | A2a | no |
/// | T4 | $T_{2}$ | $k_{17}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 23 | I222 | D2^8 | 2 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 23 | I222 | D2^8 | 2 | 1 |
///
/// ### Isotropy subgroups for T3 ($T_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 24 | I2_12_12_1 | D2^9 | 2 | 1 |
///
/// ### Isotropy subgroups for T4 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 24 | I2_12_12_1 | D2^9 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{16}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{16}t4$ | 1 | A2a | no |
/// | Z3 | $Z_{4}$ | $k_{16}t2$ | 1 | A2a | no |
/// | Z4 | $Z_{2}$ | $k_{16}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | S2 | $S_{2}$ | $k_{12}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,0)} | 21 | C222 | D2^6 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for S2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,0)} | 21 | C222 | D2^6 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{3}$ | $k_{15}t4$ | 1 | A2a | no |
/// | Y3 | $Y_{4}$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y4 | $Y_{2}$ | $k_{15}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
pub struct Sg21;

/// # 22 F222 (D₂⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{14}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 22 | F222 | D2^7 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{10}t1$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{P3(1)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{P5(1)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{P4(1)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{C3(2)/4D} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{C2(2)/4D} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{P12(1)/4D} | 22 | F222 | D2^7 | 8 | 4 |
/// | dir{P11(1)/4D} | 22 | F222 | D2^7 | 8 | 4 |
/// | dir{C9(2)/4D} | 5 | C2 | C2^3 | 16 | 4 |
/// | dir{C8(2)/4D} | 5 | C2 | C2^3 | 16 | 4 |
/// | dir{C13(2)/4D} | 5 | C2 | C2^3 | 16 | 4 |
/// | dir{C12(2)/4D} | 5 | C2 | C2^3 | 16 | 4 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 16 | 4 |
/// | dir{C10(2)/4D} | 5 | C2 | C2^3 | 16 | 4 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 4 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $X_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | T2 | $X_{3}$ | $k_{15}t4$ | 1 | A2a | no |
/// | T3 | $X_{4}$ | $k_{15}t2$ | 1 | A2a | no |
/// | T4 | $X_{2}$ | $k_{15}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for T2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ### Isotropy subgroups for T3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for T4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{17}t4$ | 1 | A2a | no |
/// | Z3 | $Z_{4}$ | $k_{17}t2$ | 1 | A2a | no |
/// | Z4 | $Z_{2}$ | $k_{17}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{16}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{3}$ | $k_{16}t4$ | 1 | A2a | no |
/// | Y3 | $Y_{4}$ | $k_{16}t2$ | 1 | A2a | no |
/// | Y4 | $Y_{2}$ | $k_{16}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
pub struct Sg22;

/// # 23 I222 (D₂⁸)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 23 | I222 | D2^8 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | X2 | $X_{3}$ | $k_{18}t4$ | 1 | A2a | no |
/// | X3 | $X_{4}$ | $k_{18}t2$ | 1 | A2a | no |
/// | X4 | $X_{2}$ | $k_{18}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{11}t1$ | 2 | B8a | no |
/// | R2 | $R_{2}$ | $k_{11}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 21 | C222 | D2^6 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 21 | C222 | D2^6 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | T2 | $T_{2}$ | $k_{12}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 21 | C222 | D2^6 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for T2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 21 | C222 | D2^6 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1W1 | $W_{1}W1$ | $k_{16}t1t1*$ | 2 | B4a | no |
/// | W2W2 | $W_{3}W3$ | $k_{16}t4t4*$ | 2 | B4a | no |
/// | W3W3 | $W_{4}W4$ | $k_{16}t2t2*$ | 2 | B4a | no |
/// | W4W4 | $W_{2}W2$ | $k_{16}t3t3*$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for W1W1 ($W_{1}W1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 22 | F222 | D2^7 | 4 | 1 |
///
/// ### Isotropy subgroups for W2W2 ($W_{3}W3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 22 | F222 | D2^7 | 4 | 1 |
///
/// ### Isotropy subgroups for W3W3 ($W_{4}W4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 22 | F222 | D2^7 | 4 | 1 |
///
/// ### Isotropy subgroups for W4W4 ($W_{2}W2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 22 | F222 | D2^7 | 4 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{10}t1$ | 2 | B8a | no |
/// | S2 | $S_{2}$ | $k_{10}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 21 | C222 | D2^6 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for S2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 21 | C222 | D2^6 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
pub struct Sg23;

/// # 24 I2_12_12_1 (D₂⁹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{17}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 24 | I2_12_12_1 | D2^9 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{18}t3$ | 1 | A2a | no |
/// | X2 | $X_{3}$ | $k_{18}t2$ | 1 | A2a | no |
/// | X3 | $X_{2}$ | $k_{18}t4$ | 1 | A2a | no |
/// | X4 | $X_{4}$ | $k_{18}t1$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for X3 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for X4 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 19 | P2_12_12_1 | D2^4 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $S_{2}$ | $k_{11}t2$ | 2 | B8a | no |
/// | R2 | $S_{1}$ | $k_{11}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 20 | C222_1 | D2^5 | 4 | 2 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 20 | C222_1 | D2^5 | 4 | 2 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{12}t2$ | 2 | B8a | no |
/// | T2 | $T_{2}$ | $k_{12}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 20 | C222_1 | D2^5 | 4 | 2 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ### Isotropy subgroups for T2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 20 | C222_1 | D2^5 | 4 | 2 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1W1 | $W_{1}W1$ | $k_{16}t1t1*$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for W1W1 ($W_{1}W1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{C12(2)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $R_{2}$ | $k_{10}t2$ | 2 | B8a | no |
/// | S2 | $R_{1}$ | $k_{10}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 20 | C222_1 | D2^5 | 4 | 2 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ### Isotropy subgroups for S2 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 20 | C222_1 | D2^5 | 4 | 2 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 2 |
///
pub struct Sg24;

/// # 25 Pmm2 (C₂ᵥ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{19}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 25 | Pmm2 | C2v^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{20}t1$ | 1 | A2a | no |
/// | X2 | $X_{3}$ | $k_{20}t2$ | 1 | A2a | no |
/// | X3 | $X_{4}$ | $k_{20}t3$ | 1 | A2a | no |
/// | X4 | $X_{2}$ | $k_{20}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{26}t1$ | 1 | A2a | no |
/// | R2 | $R_{3}$ | $k_{26}t2$ | 1 | A2a | no |
/// | R3 | $R_{4}$ | $k_{26}t3$ | 1 | A2a | no |
/// | R4 | $R_{2}$ | $k_{26}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 42 | Fmm2 | C2v^18 | 2 | 1 |
///
/// ### Isotropy subgroups for R2 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 42 | Fmm2 | C2v^18 | 2 | 1 |
///
/// ### Isotropy subgroups for R3 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 42 | Fmm2 | C2v^18 | 2 | 1 |
///
/// ### Isotropy subgroups for R4 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 42 | Fmm2 | C2v^18 | 2 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{23}t1$ | 1 | A2a | no |
/// | T2 | $T_{3}$ | $k_{23}t2$ | 1 | A2a | no |
/// | T3 | $T_{4}$ | $k_{23}t3$ | 1 | A2a | no |
/// | T4 | $T_{2}$ | $k_{23}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 38 | Amm2 | C2v^14 | 2 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 39 | Abm2 | C2v^15 | 2 | 1 |
///
/// ### Isotropy subgroups for T3 ($T_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 38 | Amm2 | C2v^14 | 2 | 1 |
///
/// ### Isotropy subgroups for T4 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 39 | Abm2 | C2v^15 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{22}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{22}t2$ | 1 | A2a | no |
/// | Z3 | $Z_{4}$ | $k_{22}t3$ | 1 | A2a | no |
/// | Z4 | $Z_{2}$ | $k_{22}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 27 | Pcc2 | C2v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{25}t1$ | 1 | A2a | no |
/// | S2 | $S_{3}$ | $k_{25}t2$ | 1 | A2a | no |
/// | S3 | $S_{4}$ | $k_{25}t3$ | 1 | A2a | no |
/// | S4 | $S_{2}$ | $k_{25}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 35 | Cmm2 | C2v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for S2 ($S_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 35 | Cmm2 | C2v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for S3 ($S_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 35 | Cmm2 | C2v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for S4 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 35 | Cmm2 | C2v^11 | 2 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $U_{1}$ | $k_{24}t1$ | 1 | A2a | no |
/// | U2 | $U_{3}$ | $k_{24}t2$ | 1 | A2a | no |
/// | U3 | $U_{4}$ | $k_{24}t3$ | 1 | A2a | no |
/// | U4 | $U_{2}$ | $k_{24}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for U1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 38 | Amm2 | C2v^14 | 2 | 1 |
///
/// ### Isotropy subgroups for U2 ($U_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 39 | Abm2 | C2v^15 | 2 | 1 |
///
/// ### Isotropy subgroups for U3 ($U_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 39 | Abm2 | C2v^15 | 2 | 1 |
///
/// ### Isotropy subgroups for U4 ($U_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 38 | Amm2 | C2v^14 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{21}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{3}$ | $k_{21}t2$ | 1 | A2a | no |
/// | Y3 | $Y_{4}$ | $k_{21}t3$ | 1 | A2a | no |
/// | Y4 | $Y_{2}$ | $k_{21}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
pub struct Sg25;

/// # 26 Pmc2_1 (C₂ᵥ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{19}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{21}t1$ | 1 | A2a | no |
/// | X2 | $X_{3}$ | $k_{21}t2$ | 1 | A2a | no |
/// | X3 | $X_{4}$ | $k_{21}t4$ | 1 | A2a | no |
/// | X4 | $X_{2}$ | $k_{21}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R3 | $R_{1}R2$ | $k_{26}t1t4$ | 2 | B4a | yes |
/// | R2R4 | $R_{3}R4$ | $k_{26}t2t3$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for R1R3 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 4 | 1 |
///
/// ### Isotropy subgroups for R2R4 ($R_{3}R4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1T3 | $T_{1}T2$ | $k_{24}t1t4$ | 2 | B4a | yes |
/// | T2T4 | $T_{3}T4$ | $k_{24}t2t3$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for T1T3 ($T_{1}T2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 4 | 1 |
///
/// ### Isotropy subgroups for T2T4 ($T_{3}T4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z3 | $Z_{1}Z2$ | $k_{22}t1t4$ | 2 | B4a | yes |
/// | Z2Z4 | $Z_{3}Z4$ | $k_{22}t2t3$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for Z1Z3 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 4 | 1 |
///
/// ### Isotropy subgroups for Z2Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{25}t1$ | 1 | A2a | no |
/// | S2 | $S_{3}$ | $k_{25}t2$ | 1 | A2a | no |
/// | S3 | $S_{4}$ | $k_{25}t4$ | 1 | A2a | no |
/// | S4 | $S_{2}$ | $k_{25}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 36 | Cmc2_1 | C2v^12 | 2 | 1 |
///
/// ### Isotropy subgroups for S2 ($S_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 36 | Cmc2_1 | C2v^12 | 2 | 1 |
///
/// ### Isotropy subgroups for S3 ($S_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 36 | Cmc2_1 | C2v^12 | 2 | 1 |
///
/// ### Isotropy subgroups for S4 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 36 | Cmc2_1 | C2v^12 | 2 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1U3 | $U_{1}U2$ | $k_{23}t1t4$ | 2 | B4a | yes |
/// | U2U4 | $U_{3}U4$ | $k_{23}t2t3$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for U1U3 ($U_{1}U2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 4 | 1 |
///
/// ### Isotropy subgroups for U2U4 ($U_{3}U4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{20}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{3}$ | $k_{20}t2$ | 1 | A2a | no |
/// | Y3 | $Y_{4}$ | $k_{20}t4$ | 1 | A2a | no |
/// | Y4 | $Y_{2}$ | $k_{20}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
pub struct Sg26;

/// # 27 Pcc2 (C₂ᵥ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{19}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 27 | Pcc2 | C2v^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{20}t1$ | 1 | A2a | no |
/// | X2 | $X_{3}$ | $k_{20}t2$ | 1 | A2a | no |
/// | X3 | $X_{4}$ | $k_{20}t3$ | 1 | A2a | no |
/// | X4 | $X_{2}$ | $k_{20}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 27 | Pcc2 | C2v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 27 | Pcc2 | C2v^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{26}t1t2$ | 2 | B4a | yes |
/// | R3R4 | $R_{3}R4$ | $k_{26}t3t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ### Isotropy subgroups for R3R4 ($R_{3}R4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1T2 | $T_{1}T2$ | $k_{23}t1t2$ | 2 | B4a | yes |
/// | T3T4 | $T_{3}T4$ | $k_{23}t3t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for T1T2 ($T_{1}T2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ### Isotropy subgroups for T3T4 ($T_{3}T4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z2 | $Z_{1}Z2$ | $k_{22}t1t2$ | 2 | B4a | yes |
/// | Z3Z4 | $Z_{3}Z4$ | $k_{22}t3t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for Z1Z2 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{25}t1$ | 1 | A2a | no |
/// | S2 | $S_{3}$ | $k_{25}t2$ | 1 | A2a | no |
/// | S3 | $S_{4}$ | $k_{25}t3$ | 1 | A2a | no |
/// | S4 | $S_{2}$ | $k_{25}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 37 | Ccc2 | C2v^13 | 2 | 1 |
///
/// ### Isotropy subgroups for S2 ($S_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 37 | Ccc2 | C2v^13 | 2 | 1 |
///
/// ### Isotropy subgroups for S3 ($S_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 37 | Ccc2 | C2v^13 | 2 | 1 |
///
/// ### Isotropy subgroups for S4 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 37 | Ccc2 | C2v^13 | 2 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1U2 | $U_{1}U2$ | $k_{24}t1t2$ | 2 | B4a | yes |
/// | U3U4 | $U_{3}U4$ | $k_{24}t3t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for U1U2 ($U_{1}U2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ### Isotropy subgroups for U3U4 ($U_{3}U4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{21}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{3}$ | $k_{21}t2$ | 1 | A2a | no |
/// | Y3 | $Y_{4}$ | $k_{21}t3$ | 1 | A2a | no |
/// | Y4 | $Y_{2}$ | $k_{21}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 27 | Pcc2 | C2v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 27 | Pcc2 | C2v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
pub struct Sg27;

/// # 28 Pma2 (C₂ᵥ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{19}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $Y_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 6 | Pm | Cs^1 | 4 | 1 |
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{26}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $U_{1}$ | $k_{23}t1$ | 1 | A2a | no |
/// | T2 | $U_{3}$ | $k_{23}t2$ | 1 | A2a | no |
/// | T3 | $U_{2}$ | $k_{23}t3$ | 1 | A2a | no |
/// | T4 | $U_{4}$ | $k_{23}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 40 | Ama2 | C2v^16 | 2 | 1 |
///
/// ### Isotropy subgroups for T2 ($U_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 41 | Aba2 | C2v^17 | 2 | 1 |
///
/// ### Isotropy subgroups for T3 ($U_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 40 | Ama2 | C2v^16 | 2 | 1 |
///
/// ### Isotropy subgroups for T4 ($U_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 41 | Aba2 | C2v^17 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{22}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{22}t2$ | 1 | A2a | no |
/// | Z3 | $Z_{2}$ | $k_{22}t3$ | 1 | A2a | no |
/// | Z4 | $Z_{4}$ | $k_{22}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{25}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $T_{1}$ | $k_{24}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for U1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $X_{1}$ | $k_{21}t1$ | 1 | A2a | no |
/// | Y2 | $X_{3}$ | $k_{21}t2$ | 1 | A2a | no |
/// | Y3 | $X_{2}$ | $k_{21}t3$ | 1 | A2a | no |
/// | Y4 | $X_{4}$ | $k_{21}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 32 | Pba2 | C2v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 32 | Pba2 | C2v^8 | 2 | 1 |
///
pub struct Sg28;

/// # 29 Pca2_1 (C₂ᵥ⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{19}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $Y_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R1 | $R_{1}R1$ | $k_{26}t1t1$ | 4 | D8a | no |
///
/// ### Isotropy subgroups for R1R1 ($R_{1}R1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1T4 | $U_{1}U2$ | $k_{23}t1t4$ | 2 | B4a | yes |
/// | T2T3 | $U_{3}U4$ | $k_{23}t2t3$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for T1T4 ($U_{1}U2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 4 | 1 |
///
/// ### Isotropy subgroups for T2T3 ($U_{3}U4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z4 | $Z_{1}Z2$ | $k_{22}t1t4$ | 2 | B4a | yes |
/// | Z2Z3 | $Z_{3}Z4$ | $k_{22}t2t3$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for Z1Z4 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 4 | 1 |
///
/// ### Isotropy subgroups for Z2Z3 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{25}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1U1 | $T_{1}T1$ | $k_{24}t1t1$ | 4 | D8a | no |
///
/// ### Isotropy subgroups for U1U1 ($T_{1}T1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $X_{1}$ | $k_{21}t1$ | 1 | A2a | no |
/// | Y2 | $X_{3}$ | $k_{21}t2$ | 1 | A2a | no |
/// | Y3 | $X_{2}$ | $k_{21}t3$ | 1 | A2a | no |
/// | Y4 | $X_{4}$ | $k_{21}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
pub struct Sg29;

/// # 30 Pnc2 (C₂ᵥ⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{19}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{21}t1$ | 1 | A2a | no |
/// | X2 | $X_{3}$ | $k_{21}t2$ | 1 | A2a | no |
/// | X3 | $X_{4}$ | $k_{21}t4$ | 1 | A2a | no |
/// | X4 | $X_{2}$ | $k_{21}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 34 | Pnn2 | C2v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 34 | Pnn2 | C2v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{26}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{24}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z2 | $Z_{1}Z2$ | $k_{22}t1t2$ | 2 | B4a | yes |
/// | Z3Z4 | $Z_{3}Z4$ | $k_{22}t3t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for Z1Z2 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{25}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1U2 | $U_{1}U2$ | $k_{23}t1t2$ | 2 | B4a | yes |
/// | U3U4 | $U_{3}U4$ | $k_{23}t3t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for U1U2 ($U_{1}U2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ### Isotropy subgroups for U3U4 ($U_{3}U4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
pub struct Sg30;

/// # 31 Pmn2_1 (C₂ᵥ⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{19}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $Y_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 6 | Pm | Cs^1 | 4 | 1 |
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{26}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{(a,a)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1T3 | $U_{1}U2$ | $k_{23}t1t3$ | 2 | B4a | yes |
/// | T2T4 | $U_{3}U4$ | $k_{23}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for T1T3 ($U_{1}U2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 4 | 1 |
///
/// ### Isotropy subgroups for T2T4 ($U_{3}U4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z3 | $Z_{1}Z2$ | $k_{22}t1t3$ | 2 | B4a | yes |
/// | Z2Z4 | $Z_{3}Z4$ | $k_{22}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for Z1Z3 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 4 | 1 |
///
/// ### Isotropy subgroups for Z2Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{25}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $T_{1}$ | $k_{24}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for U1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $X_{1}$ | $k_{21}t1$ | 1 | A2a | no |
/// | Y2 | $X_{3}$ | $k_{21}t2$ | 1 | A2a | no |
/// | Y3 | $X_{2}$ | $k_{21}t3$ | 1 | A2a | no |
/// | Y4 | $X_{4}$ | $k_{21}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
pub struct Sg31;

/// # 32 Pba2 (C₂ᵥ⁸)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{19}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 32 | Pba2 | C2v^8 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{3}R4$ | $k_{26}t3t4$ | 2 | B4a | yes |
/// | R3R4 | $R_{1}R2$ | $k_{26}t1t2$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for R1R2 ($R_{3}R4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ### Isotropy subgroups for R3R4 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{23}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{22}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{22}t2$ | 1 | A2a | no |
/// | Z3 | $Z_{4}$ | $k_{22}t3$ | 1 | A2a | no |
/// | Z4 | $Z_{2}$ | $k_{22}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 32 | Pba2 | C2v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 34 | Pnn2 | C2v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1S2 | $S_{3}S4$ | $k_{25}t3t4$ | 2 | B4a | yes |
/// | S3S4 | $S_{1}S2$ | $k_{25}t1t2$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for S1S2 ($S_{3}S4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ### Isotropy subgroups for S3S4 ($S_{1}S2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $U_{1}$ | $k_{24}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for U1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
pub struct Sg32;

/// # 33 Pna2_1 (C₂ᵥ⁹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{19}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $Y_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R3 | $R_{1}R2$ | $k_{26}t1t3$ | 2 | B4a | yes |
/// | R2R4 | $R_{3}R4$ | $k_{26}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for R1R3 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 4 | 1 |
///
/// ### Isotropy subgroups for R2R4 ($R_{3}R4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 4 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $U_{1}$ | $k_{23}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for T1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{(a,a)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z4 | $Z_{1}Z2$ | $k_{22}t1t4$ | 2 | B4a | yes |
/// | Z2Z3 | $Z_{3}Z4$ | $k_{22}t2t3$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for Z1Z4 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 4 | 1 |
///
/// ### Isotropy subgroups for Z2Z3 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1S2 | $S_{3}S4$ | $k_{25}t3t4$ | 2 | B4a | yes |
/// | S3S4 | $S_{1}S2$ | $k_{25}t1t2$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for S1S2 ($S_{3}S4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ### Isotropy subgroups for S3S4 ($S_{1}S2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1U1 | $T_{1}T1$ | $k_{24}t1t1$ | 4 | D8a | no |
///
/// ### Isotropy subgroups for U1U1 ($T_{1}T1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $X_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
pub struct Sg33;

/// # 34 Pnn2 (C₂ᵥ¹⁰)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{19}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 34 | Pnn2 | C2v^10 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{4}$ | $k_{26}t3$ | 1 | A2a | no |
/// | R2 | $R_{2}$ | $k_{26}t4$ | 1 | A2a | no |
/// | R3 | $R_{1}$ | $k_{26}t1$ | 1 | A2a | no |
/// | R4 | $R_{3}$ | $k_{26}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for R1 ($R_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 43 | Fdd2 | C2v^19 | 2 | 1 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 43 | Fdd2 | C2v^19 | 2 | 1 |
///
/// ### Isotropy subgroups for R3 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 43 | Fdd2 | C2v^19 | 2 | 1 |
///
/// ### Isotropy subgroups for R4 ($R_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 43 | Fdd2 | C2v^19 | 2 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{23}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z2 | $Z_{1}Z2$ | $k_{22}t1t2$ | 2 | B4a | yes |
/// | Z3Z4 | $Z_{3}Z4$ | $k_{22}t3t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for Z1Z2 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1S2 | $S_{3}S4$ | $k_{25}t3t4$ | 2 | B4a | yes |
/// | S3S4 | $S_{1}S2$ | $k_{25}t1t2$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for S1S2 ($S_{3}S4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ### Isotropy subgroups for S3S4 ($S_{1}S2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $U_{1}$ | $k_{24}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for U1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
pub struct Sg34;

/// # 35 Cmm2 (C₂ᵥ¹¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{14}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 35 | Cmm2 | C2v^11 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | R2 | $R_{2}$ | $k_{13}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 42 | Fmm2 | C2v^18 | 4 | 2 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 42 | Fmm2 | C2v^18 | 4 | 2 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | T2 | $T_{3}$ | $k_{17}t2$ | 1 | A2a | no |
/// | T3 | $T_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | T4 | $T_{2}$ | $k_{17}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 44 | Imm2 | C2v^20 | 2 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 45 | Iba2 | C2v^21 | 2 | 1 |
///
/// ### Isotropy subgroups for T3 ($T_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 46 | Ima2 | C2v^22 | 2 | 1 |
///
/// ### Isotropy subgroups for T4 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 46 | Ima2 | C2v^22 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{16}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{16}t2$ | 1 | A2a | no |
/// | Z3 | $Z_{4}$ | $k_{16}t3$ | 1 | A2a | no |
/// | Z4 | $Z_{2}$ | $k_{16}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 35 | Cmm2 | C2v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 37 | Ccc2 | C2v^13 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 36 | Cmc2_1 | C2v^12 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 36 | Cmc2_1 | C2v^12 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | S2 | $S_{2}$ | $k_{12}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,0)} | 35 | Cmm2 | C2v^11 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for S2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,0)} | 35 | Cmm2 | C2v^11 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{3}$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y3 | $Y_{4}$ | $k_{15}t3$ | 1 | A2a | no |
/// | Y4 | $Y_{2}$ | $k_{15}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 32 | Pba2 | C2v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
pub struct Sg35;

/// # 36 Cmc2_1 (C₂ᵥ¹²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{14}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 36 | Cmc2_1 | C2v^12 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{13}t1t2$ | 4 | D16c | yes |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{C12(2)/4D} | 9 | Cc | Cs^4 | 8 | 2 |
/// | dir{C11(2)/4D} | 8 | Cm | Cs^3 | 8 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1T3 | $T_{1}T2$ | $k_{17}t1t3$ | 2 | B4a | yes |
/// | T2T4 | $T_{3}T4$ | $k_{17}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for T1T3 ($T_{1}T2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 4 | 1 |
///
/// ### Isotropy subgroups for T2T4 ($T_{3}T4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z3 | $Z_{1}Z2$ | $k_{16}t1t3$ | 2 | B4a | yes |
/// | Z2Z4 | $Z_{3}Z4$ | $k_{16}t2t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for Z1Z3 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 4 | 1 |
///
/// ### Isotropy subgroups for Z2Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 4 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | S2 | $S_{2}$ | $k_{12}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,0)} | 36 | Cmc2_1 | C2v^12 | 4 | 2 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ### Isotropy subgroups for S2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,0)} | 36 | Cmc2_1 | C2v^12 | 4 | 2 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{3}$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y3 | $Y_{2}$ | $k_{15}t3$ | 1 | A2a | no |
/// | Y4 | $Y_{4}$ | $k_{15}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
pub struct Sg36;

/// # 37 Ccc2 (C₂ᵥ¹³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{14}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 37 | Ccc2 | C2v^13 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | R2 | $R_{2}$ | $k_{13}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 43 | Fdd2 | C2v^19 | 4 | 2 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 43 | Fdd2 | C2v^19 | 4 | 2 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1T2 | $T_{1}T2$ | $k_{17}t1t2$ | 2 | B4a | yes |
/// | T3T4 | $T_{3}T4$ | $k_{17}t3t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for T1T2 ($T_{1}T2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ### Isotropy subgroups for T3T4 ($T_{3}T4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 5 | C2 | C2^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z2 | $Z_{1}Z2$ | $k_{16}t1t2$ | 2 | B4a | yes |
/// | Z3Z4 | $Z_{3}Z4$ | $k_{16}t3t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for Z1Z2 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | S2 | $S_{2}$ | $k_{12}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,0)} | 37 | Ccc2 | C2v^13 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for S2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,0)} | 37 | Ccc2 | C2v^13 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{3}$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y3 | $Y_{4}$ | $k_{15}t3$ | 1 | A2a | no |
/// | Y4 | $Y_{2}$ | $k_{15}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 27 | Pcc2 | C2v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 34 | Pnn2 | C2v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
pub struct Sg37;

/// # 38 Amm2 (C₂ᵥ¹⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{14}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 38 | Amm2 | C2v^14 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{13}t1$ | 2 | B8a | no |
/// | R2 | $R_{2}$ | $k_{13}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,0)} | 42 | Fmm2 | C2v^18 | 4 | 2 |
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,0)} | 42 | Fmm2 | C2v^18 | 4 | 2 |
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | T2 | $T_{3}$ | $k_{17}t2$ | 1 | A2a | no |
/// | T3 | $T_{2}$ | $k_{17}t3$ | 1 | A2a | no |
/// | T4 | $T_{4}$ | $k_{17}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 44 | Imm2 | C2v^20 | 2 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 46 | Ima2 | C2v^22 | 2 | 1 |
///
/// ### Isotropy subgroups for T3 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 44 | Imm2 | C2v^20 | 2 | 1 |
///
/// ### Isotropy subgroups for T4 ($T_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 46 | Ima2 | C2v^22 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{16}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{16}t2$ | 1 | A2a | no |
/// | Z3 | $Z_{2}$ | $k_{16}t3$ | 1 | A2a | no |
/// | Z4 | $Z_{4}$ | $k_{16}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 38 | Amm2 | C2v^14 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 40 | Ama2 | C2v^16 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 38 | Amm2 | C2v^14 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 40 | Ama2 | C2v^16 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | S2 | $S_{2}$ | $k_{12}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 6 | Pm | Cs^1 | 4 | 1 |
/// | dir{(a,0)} | 38 | Amm2 | C2v^14 | 4 | 2 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 2 |
///
/// ### Isotropy subgroups for S2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,0)} | 39 | Abm2 | C2v^15 | 4 | 2 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{3}$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y3 | $Y_{2}$ | $k_{15}t3$ | 1 | A2a | no |
/// | Y4 | $Y_{4}$ | $k_{15}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
pub struct Sg38;

/// # 39 Abm2 (C₂ᵥ¹⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{14}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 39 | Abm2 | C2v^15 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{13}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{C12(2)/4D} | 8 | Cm | Cs^3 | 8 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{3}$ | $k_{17}t1$ | 1 | A2a | no |
/// | T2 | $T_{1}$ | $k_{17}t2$ | 1 | A2a | no |
/// | T3 | $T_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | T4 | $T_{2}$ | $k_{17}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1 ($T_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 45 | Iba2 | C2v^21 | 2 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 46 | Ima2 | C2v^22 | 2 | 1 |
///
/// ### Isotropy subgroups for T3 ($T_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 45 | Iba2 | C2v^21 | 2 | 1 |
///
/// ### Isotropy subgroups for T4 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 46 | Ima2 | C2v^22 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{16}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{16}t2$ | 1 | A2a | no |
/// | Z3 | $Z_{2}$ | $k_{16}t3$ | 1 | A2a | no |
/// | Z4 | $Z_{4}$ | $k_{16}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 39 | Abm2 | C2v^15 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 41 | Aba2 | C2v^17 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 39 | Abm2 | C2v^15 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 41 | Aba2 | C2v^17 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1S2 | $S_{1}S2$ | $k_{12}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for S1S2 ($S_{1}S2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{C12(2)/4D} | 8 | Cm | Cs^3 | 8 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{3}$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{1}$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y3 | $Y_{4}$ | $k_{15}t3$ | 1 | A2a | no |
/// | Y4 | $Y_{2}$ | $k_{15}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 27 | Pcc2 | C2v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
pub struct Sg39;

/// # 40 Ama2 (C₂ᵥ¹⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{14}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 40 | Ama2 | C2v^16 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{2}$ | $k_{13}t1$ | 2 | B8a | no |
/// | R2 | $R_{1}$ | $k_{13}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,0)} | 43 | Fdd2 | C2v^19 | 4 | 2 |
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,0)} | 43 | Fdd2 | C2v^19 | 4 | 2 |
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 8 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{17}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{16}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 6 | Pm | Cs^1 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | S2 | $S_{2}$ | $k_{12}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 6 | Pm | Cs^1 | 4 | 1 |
/// | dir{(a,0)} | 40 | Ama2 | C2v^16 | 4 | 2 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 2 |
///
/// ### Isotropy subgroups for S2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,0)} | 41 | Aba2 | C2v^17 | 4 | 2 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{3}$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y3 | $Y_{2}$ | $k_{15}t3$ | 1 | A2a | no |
/// | Y4 | $Y_{4}$ | $k_{15}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 34 | Pnn2 | C2v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
pub struct Sg40;

/// # 41 Aba2 (C₂ᵥ¹⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{14}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 41 | Aba2 | C2v^17 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{13}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{C12(2)/4D} | 9 | Cc | Cs^4 | 8 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{17}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{16}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1S2 | $S_{1}S2$ | $k_{12}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for S1S2 ($S_{1}S2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{C12(2)/4D} | 9 | Cc | Cs^4 | 8 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{3}$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{1}$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y3 | $Y_{4}$ | $k_{15}t3$ | 1 | A2a | no |
/// | Y4 | $Y_{2}$ | $k_{15}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 32 | Pba2 | C2v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
pub struct Sg41;

/// # 42 Fmm2 (C₂ᵥ¹⁸)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{14}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 42 | Fmm2 | C2v^18 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{10}t1$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{P5(1)/4D} | 8 | Cm | Cs^3 | 8 | 2 |
/// | dir{P4(1)/4D} | 8 | Cm | Cs^3 | 8 | 2 |
/// | dir{P3(1)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{C3(2)/4D} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{C2(2)/4D} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{P12(1)/4D} | 43 | Fdd2 | C2v^19 | 8 | 4 |
/// | dir{P11(1)/4D} | 42 | Fmm2 | C2v^18 | 8 | 4 |
/// | dir{C13(2)/4D} | 9 | Cc | Cs^4 | 16 | 4 |
/// | dir{C11(2)/4D} | 9 | Cc | Cs^4 | 16 | 4 |
/// | dir{C12(2)/4D} | 8 | Cm | Cs^3 | 16 | 4 |
/// | dir{C10(2)/4D} | 8 | Cm | Cs^3 | 16 | 4 |
/// | dir{C9(2)/4D} | 5 | C2 | C2^3 | 16 | 4 |
/// | dir{C8(2)/4D} | 5 | C2 | C2^3 | 16 | 4 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 4 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $X_{1}$ | $k_{15}t1$ | 1 | A2a | no |
/// | T2 | $X_{3}$ | $k_{15}t2$ | 1 | A2a | no |
/// | T3 | $X_{4}$ | $k_{15}t3$ | 1 | A2a | no |
/// | T4 | $X_{2}$ | $k_{15}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 38 | Amm2 | C2v^14 | 2 | 1 |
///
/// ### Isotropy subgroups for T2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 41 | Aba2 | C2v^17 | 2 | 1 |
///
/// ### Isotropy subgroups for T3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 40 | Ama2 | C2v^16 | 2 | 1 |
///
/// ### Isotropy subgroups for T4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 39 | Abm2 | C2v^15 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{17}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{3}$ | $k_{17}t2$ | 1 | A2a | no |
/// | Z3 | $Z_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | Z4 | $Z_{2}$ | $k_{17}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 35 | Cmm2 | C2v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 37 | Ccc2 | C2v^13 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3 ($Z_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 36 | Cmc2_1 | C2v^12 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 36 | Cmc2_1 | C2v^12 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{16}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{3}$ | $k_{16}t2$ | 1 | A2a | no |
/// | Y3 | $Y_{4}$ | $k_{16}t3$ | 1 | A2a | no |
/// | Y4 | $Y_{2}$ | $k_{16}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 38 | Amm2 | C2v^14 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 41 | Aba2 | C2v^17 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3 ($Y_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 39 | Abm2 | C2v^15 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 40 | Ama2 | C2v^16 | 2 | 1 |
///
pub struct Sg42;

/// # 43 Fdd2 (C₂ᵥ¹⁹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{14}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 43 | Fdd2 | C2v^19 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{10}t1$ | 4 | D32b | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{P5(1)/4D} | 9 | Cc | Cs^4 | 8 | 2 |
/// | dir{P4(1)/4D} | 9 | Cc | Cs^4 | 8 | 2 |
/// | dir{P3(1)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{C3(2)/4D} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{C2(2)/4D} | 1 | P1 | C1^1 | 16 | 2 |
/// | dir{C8(2)/4D} | 5 | C2 | C2^3 | 16 | 4 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 4 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z2 | $Z_{1}Z2$ | $k_{17}t1t2$ | 2 | B4a | yes |
/// | Z3Z4 | $Z_{3}Z4$ | $k_{17}t3t4$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for Z1Z2 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 3 | P2 | C2^1 | 4 | 1 |
///
/// ### Isotropy subgroups for Z3Z4 ($Z_{3}Z4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{16}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
pub struct Sg43;

/// # 44 Imm2 (C₂ᵥ²⁰)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 44 | Imm2 | C2v^20 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{18}t1$ | 1 | A2a | no |
/// | X2 | $X_{3}$ | $k_{18}t2$ | 1 | A2a | no |
/// | X3 | $X_{4}$ | $k_{18}t3$ | 1 | A2a | no |
/// | X4 | $X_{2}$ | $k_{18}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 34 | Pnn2 | C2v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{11}t1$ | 2 | B8a | no |
/// | R2 | $R_{2}$ | $k_{11}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,0)} | 38 | Amm2 | C2v^14 | 4 | 2 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 2 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{(a,0)} | 39 | Abm2 | C2v^15 | 4 | 2 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | T2 | $T_{2}$ | $k_{12}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 35 | Cmm2 | C2v^11 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for T2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 35 | Cmm2 | C2v^11 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1 | $W_{1}$ | $k_{16}t1$ | 2 | B8a | no |
/// | W2 | $W_{2}$ | $k_{16}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for W1 ($W_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 43 | Fdd2 | C2v^19 | 4 | 2 |
/// | dir{(a,a)} | 42 | Fmm2 | C2v^18 | 4 | 2 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for W2 ($W_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 43 | Fdd2 | C2v^19 | 4 | 2 |
/// | dir{(a,a)} | 42 | Fmm2 | C2v^18 | 4 | 2 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{10}t1$ | 2 | B8a | no |
/// | S2 | $S_{2}$ | $k_{10}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,0)} | 38 | Amm2 | C2v^14 | 4 | 2 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 2 |
///
/// ### Isotropy subgroups for S2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{(a,0)} | 39 | Abm2 | C2v^15 | 4 | 2 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 2 |
///
pub struct Sg44;

/// # 45 Iba2 (C₂ᵥ²¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{4}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{2}$ | $k_{17}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 45 | Iba2 | C2v^21 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{18}t2$ | 1 | A2a | no |
/// | X2 | $X_{3}$ | $k_{18}t1$ | 1 | A2a | no |
/// | X3 | $X_{4}$ | $k_{18}t4$ | 1 | A2a | no |
/// | X4 | $X_{2}$ | $k_{18}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 27 | Pcc2 | C2v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 32 | Pba2 | C2v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for X3 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for X4 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{11}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{C12(2)/4D} | 9 | Cc | Cs^4 | 8 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{12}t1$ | 2 | B8a | no |
/// | T2 | $T_{2}$ | $k_{12}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 37 | Ccc2 | C2v^13 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ### Isotropy subgroups for T2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 37 | Ccc2 | C2v^13 | 4 | 2 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1W1 | $W_{1}W1$ | $k_{16}t1t1$ | 4 | D8a | no |
/// | W2W2 | $W_{2}W2$ | $k_{16}t2t2$ | 4 | D8a | no |
///
/// ### Isotropy subgroups for W1W1 ($W_{1}W1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{4D1(4)/4D} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ### Isotropy subgroups for W2W2 ($W_{2}W2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{4D1(4)/4D} | 5 | C2 | C2^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1S2 | $S_{1}S2$ | $k_{10}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for S1S2 ($S_{1}S2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{C12(2)/4D} | 9 | Cc | Cs^4 | 8 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 2 |
///
pub struct Sg45;

/// # 46 Ima2 (C₂ᵥ²²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{3}$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM3 | $\Gamma_{2}$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM4 | $\Gamma_{4}$ | $k_{17}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 46 | Ima2 | C2v^22 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4 ($\Gamma_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{18}t2$ | 1 | A2a | no |
/// | X2 | $X_{3}$ | $k_{18}t1$ | 1 | A2a | no |
/// | X3 | $X_{2}$ | $k_{18}t4$ | 1 | A2a | no |
/// | X4 | $X_{4}$ | $k_{18}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for X3 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for X4 ($X_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $S_{1}S2$ | $k_{11}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for R1R2 ($S_{1}S2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 8 | 1 |
/// | dir{C12(2)/4D} | 8 | Cm | Cs^3 | 8 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{2}$ | $k_{12}t2$ | 2 | B8a | no |
/// | T2 | $T_{1}$ | $k_{12}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for T1 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 36 | Cmc2_1 | C2v^12 | 4 | 2 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ### Isotropy subgroups for T2 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 36 | Cmc2_1 | C2v^12 | 4 | 2 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1W2 | $W_{1}W2$ | $k_{16}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for W1W2 ($W_{1}W2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C12(2)/4D} | 9 | Cc | Cs^4 | 8 | 2 |
/// | dir{C11(2)/4D} | 8 | Cm | Cs^3 | 8 | 2 |
/// | dir{C1(2)/4D} | 5 | C2 | C2^3 | 8 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 2 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $R_{1}$ | $k_{10}t2$ | 2 | B8a | no |
/// | S2 | $R_{2}$ | $k_{10}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,0)} | 40 | Ama2 | C2v^16 | 4 | 2 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 2 |
///
/// ### Isotropy subgroups for S2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{(a,0)} | 41 | Aba2 | C2v^17 | 4 | 2 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 2 |
///
pub struct Sg46;

/// # 47 Pmmm (D₂ₕ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{19}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 47 | Pmmm | D2h^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{20}t1$ | 1 | A2a | no |
/// | X1- | $X_{1}^-$ | $k_{20}t2$ | 1 | A2a | no |
/// | X2+ | $X_{3}^+$ | $k_{20}t7$ | 1 | A2a | no |
/// | X2- | $X_{3}^-$ | $k_{20}t8$ | 1 | A2a | no |
/// | X3+ | $X_{4}^+$ | $k_{20}t3$ | 1 | A2a | no |
/// | X3- | $X_{4}^-$ | $k_{20}t4$ | 1 | A2a | no |
/// | X4+ | $X_{2}^+$ | $k_{20}t5$ | 1 | A2a | no |
/// | X4- | $X_{2}^-$ | $k_{20}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 47 | Pmmm | D2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for X2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for X2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for X3+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for X3- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 47 | Pmmm | D2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for X4+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for X4- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{26}t1$ | 1 | A2a | no |
/// | R1- | $R_{1}^-$ | $k_{26}t2$ | 1 | A2a | no |
/// | R2+ | $R_{3}^+$ | $k_{26}t7$ | 1 | A2a | no |
/// | R2- | $R_{3}^-$ | $k_{26}t8$ | 1 | A2a | no |
/// | R3+ | $R_{4}^+$ | $k_{26}t3$ | 1 | A2a | no |
/// | R3- | $R_{4}^-$ | $k_{26}t4$ | 1 | A2a | no |
/// | R4+ | $R_{2}^+$ | $k_{26}t5$ | 1 | A2a | no |
/// | R4- | $R_{2}^-$ | $k_{26}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 69 | Fmmm | D2h^23 | 2 | 1 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 69 | Fmmm | D2h^23 | 2 | 1 |
///
/// ### Isotropy subgroups for R2+ ($R_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 69 | Fmmm | D2h^23 | 2 | 1 |
///
/// ### Isotropy subgroups for R2- ($R_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 69 | Fmmm | D2h^23 | 2 | 1 |
///
/// ### Isotropy subgroups for R3+ ($R_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 69 | Fmmm | D2h^23 | 2 | 1 |
///
/// ### Isotropy subgroups for R3- ($R_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 69 | Fmmm | D2h^23 | 2 | 1 |
///
/// ### Isotropy subgroups for R4+ ($R_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 69 | Fmmm | D2h^23 | 2 | 1 |
///
/// ### Isotropy subgroups for R4- ($R_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 69 | Fmmm | D2h^23 | 2 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1+ | $T_{1}^+$ | $k_{23}t1$ | 1 | A2a | no |
/// | T1- | $T_{1}^-$ | $k_{23}t2$ | 1 | A2a | no |
/// | T2+ | $T_{3}^+$ | $k_{23}t7$ | 1 | A2a | no |
/// | T2- | $T_{3}^-$ | $k_{23}t8$ | 1 | A2a | no |
/// | T3+ | $T_{4}^+$ | $k_{23}t3$ | 1 | A2a | no |
/// | T3- | $T_{4}^-$ | $k_{23}t4$ | 1 | A2a | no |
/// | T4+ | $T_{2}^+$ | $k_{23}t5$ | 1 | A2a | no |
/// | T4- | $T_{2}^-$ | $k_{23}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1+ ($T_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for T1- ($T_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for T2+ ($T_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for T2- ($T_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for T3+ ($T_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for T3- ($T_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for T4+ ($T_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for T4- ($T_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{22}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{22}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{3}^+$ | $k_{22}t7$ | 1 | A2a | no |
/// | Z2- | $Z_{3}^-$ | $k_{22}t8$ | 1 | A2a | no |
/// | Z3+ | $Z_{4}^+$ | $k_{22}t3$ | 1 | A2a | no |
/// | Z3- | $Z_{4}^-$ | $k_{22}t4$ | 1 | A2a | no |
/// | Z4+ | $Z_{2}^+$ | $k_{22}t5$ | 1 | A2a | no |
/// | Z4- | $Z_{2}^-$ | $k_{22}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 47 | Pmmm | D2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 47 | Pmmm | D2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3+ ($Z_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3- ($Z_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1+ | $S_{1}^+$ | $k_{25}t1$ | 1 | A2a | no |
/// | S1- | $S_{1}^-$ | $k_{25}t2$ | 1 | A2a | no |
/// | S2+ | $S_{3}^+$ | $k_{25}t7$ | 1 | A2a | no |
/// | S2- | $S_{3}^-$ | $k_{25}t8$ | 1 | A2a | no |
/// | S3+ | $S_{4}^+$ | $k_{25}t3$ | 1 | A2a | no |
/// | S3- | $S_{4}^-$ | $k_{25}t4$ | 1 | A2a | no |
/// | S4+ | $S_{2}^+$ | $k_{25}t5$ | 1 | A2a | no |
/// | S4- | $S_{2}^-$ | $k_{25}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for S1+ ($S_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for S1- ($S_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for S2+ ($S_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for S2- ($S_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for S3+ ($S_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for S3- ($S_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for S4+ ($S_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for S4- ($S_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1+ | $U_{1}^+$ | $k_{24}t1$ | 1 | A2a | no |
/// | U1- | $U_{1}^-$ | $k_{24}t2$ | 1 | A2a | no |
/// | U2+ | $U_{3}^+$ | $k_{24}t7$ | 1 | A2a | no |
/// | U2- | $U_{3}^-$ | $k_{24}t8$ | 1 | A2a | no |
/// | U3+ | $U_{4}^+$ | $k_{24}t3$ | 1 | A2a | no |
/// | U3- | $U_{4}^-$ | $k_{24}t4$ | 1 | A2a | no |
/// | U4+ | $U_{2}^+$ | $k_{24}t5$ | 1 | A2a | no |
/// | U4- | $U_{2}^-$ | $k_{24}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for U1+ ($U_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for U1- ($U_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for U2+ ($U_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for U2- ($U_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for U3+ ($U_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for U3- ($U_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for U4+ ($U_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for U4- ($U_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Y_{1}^+$ | $k_{21}t1$ | 1 | A2a | no |
/// | Y1- | $Y_{1}^-$ | $k_{21}t2$ | 1 | A2a | no |
/// | Y2+ | $Y_{3}^+$ | $k_{21}t7$ | 1 | A2a | no |
/// | Y2- | $Y_{3}^-$ | $k_{21}t8$ | 1 | A2a | no |
/// | Y3+ | $Y_{4}^+$ | $k_{21}t3$ | 1 | A2a | no |
/// | Y3- | $Y_{4}^-$ | $k_{21}t4$ | 1 | A2a | no |
/// | Y4+ | $Y_{2}^+$ | $k_{21}t5$ | 1 | A2a | no |
/// | Y4- | $Y_{2}^-$ | $k_{21}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 47 | Pmmm | D2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Y_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Y_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3+ ($Y_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3- ($Y_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 47 | Pmmm | D2h^1 | 2 | 1 |
///
pub struct Sg47;

/// # 48 Pnnn (D₂ₕ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{19}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 48 | Pnnn | D2h^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 34 | Pnn2 | C2v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 34 | Pnn2 | C2v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 34 | Pnn2 | C2v^10 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | X2 | $X_{2}$ | $k_{20}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^-$ | $k_{26}t2$ | 1 | A2a | no |
/// | R1- | $R_{1}^+$ | $k_{26}t1$ | 1 | A2a | no |
/// | R2+ | $R_{3}^-$ | $k_{26}t8$ | 1 | A2a | no |
/// | R2- | $R_{3}^+$ | $k_{26}t7$ | 1 | A2a | no |
/// | R3+ | $R_{4}^-$ | $k_{26}t4$ | 1 | A2a | no |
/// | R3- | $R_{4}^+$ | $k_{26}t3$ | 1 | A2a | no |
/// | R4+ | $R_{2}^-$ | $k_{26}t6$ | 1 | A2a | no |
/// | R4- | $R_{2}^+$ | $k_{26}t5$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 70 | Fddd | D2h^24 | 2 | 1 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 70 | Fddd | D2h^24 | 2 | 1 |
///
/// ### Isotropy subgroups for R2+ ($R_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 70 | Fddd | D2h^24 | 2 | 1 |
///
/// ### Isotropy subgroups for R2- ($R_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 70 | Fddd | D2h^24 | 2 | 1 |
///
/// ### Isotropy subgroups for R3+ ($R_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 70 | Fddd | D2h^24 | 2 | 1 |
///
/// ### Isotropy subgroups for R3- ($R_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 70 | Fddd | D2h^24 | 2 | 1 |
///
/// ### Isotropy subgroups for R4+ ($R_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 70 | Fddd | D2h^24 | 2 | 1 |
///
/// ### Isotropy subgroups for R4- ($R_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 70 | Fddd | D2h^24 | 2 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{23}t1$ | 2 | B8a | yes |
/// | T2 | $T_{2}$ | $k_{23}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{22}t1$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{22}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{25}t1$ | 2 | B8a | yes |
/// | S2 | $S_{2}$ | $k_{25}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for S2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $U_{1}$ | $k_{24}t1$ | 2 | B8a | yes |
/// | U2 | $U_{2}$ | $k_{24}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for U1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for U2 ($U_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
/// | Y2 | $Y_{2}$ | $k_{21}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg48;

/// # 49 Pccm (D₂ₕ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{19}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 27 | Pcc2 | C2v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{20}t1$ | 1 | A2a | no |
/// | X1- | $X_{1}^-$ | $k_{20}t2$ | 1 | A2a | no |
/// | X2+ | $X_{3}^+$ | $k_{20}t7$ | 1 | A2a | no |
/// | X2- | $X_{3}^-$ | $k_{20}t8$ | 1 | A2a | no |
/// | X3+ | $X_{4}^+$ | $k_{20}t3$ | 1 | A2a | no |
/// | X3- | $X_{4}^-$ | $k_{20}t4$ | 1 | A2a | no |
/// | X4+ | $X_{2}^+$ | $k_{20}t5$ | 1 | A2a | no |
/// | X4- | $X_{2}^-$ | $k_{20}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 50 | Pban | D2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for X2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for X2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for X3+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 50 | Pban | D2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for X3- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for X4+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for X4- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{26}t1$ | 2 | B8a | yes |
/// | R2 | $R_{2}$ | $k_{26}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 22 | F222 | D2^7 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 22 | F222 | D2^7 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{23}t1$ | 2 | B8a | yes |
/// | T2 | $T_{2}$ | $k_{23}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{22}t1$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{22}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{(a,0)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{(a,0)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1+ | $S_{1}^+$ | $k_{25}t1$ | 1 | A2a | no |
/// | S1- | $S_{1}^-$ | $k_{25}t2$ | 1 | A2a | no |
/// | S2+ | $S_{3}^+$ | $k_{25}t7$ | 1 | A2a | no |
/// | S2- | $S_{3}^-$ | $k_{25}t8$ | 1 | A2a | no |
/// | S3+ | $S_{4}^+$ | $k_{25}t3$ | 1 | A2a | no |
/// | S3- | $S_{4}^-$ | $k_{25}t4$ | 1 | A2a | no |
/// | S4+ | $S_{2}^+$ | $k_{25}t5$ | 1 | A2a | no |
/// | S4- | $S_{2}^-$ | $k_{25}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for S1+ ($S_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 66 | Cccm | D2h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for S1- ($S_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 68 | Ccca | D2h^22 | 2 | 1 |
///
/// ### Isotropy subgroups for S2+ ($S_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 66 | Cccm | D2h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for S2- ($S_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 68 | Ccca | D2h^22 | 2 | 1 |
///
/// ### Isotropy subgroups for S3+ ($S_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 68 | Ccca | D2h^22 | 2 | 1 |
///
/// ### Isotropy subgroups for S3- ($S_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 66 | Cccm | D2h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for S4+ ($S_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 68 | Ccca | D2h^22 | 2 | 1 |
///
/// ### Isotropy subgroups for S4- ($S_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 66 | Cccm | D2h^20 | 2 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $U_{1}$ | $k_{24}t1$ | 2 | B8a | yes |
/// | U2 | $U_{2}$ | $k_{24}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for U1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for U2 ($U_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Y_{1}^+$ | $k_{21}t1$ | 1 | A2a | no |
/// | Y1- | $Y_{1}^-$ | $k_{21}t2$ | 1 | A2a | no |
/// | Y2+ | $Y_{3}^+$ | $k_{21}t7$ | 1 | A2a | no |
/// | Y2- | $Y_{3}^-$ | $k_{21}t8$ | 1 | A2a | no |
/// | Y3+ | $Y_{4}^+$ | $k_{21}t3$ | 1 | A2a | no |
/// | Y3- | $Y_{4}^-$ | $k_{21}t4$ | 1 | A2a | no |
/// | Y4+ | $Y_{2}^+$ | $k_{21}t5$ | 1 | A2a | no |
/// | Y4- | $Y_{2}^-$ | $k_{21}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 50 | Pban | D2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Y_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Y_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3+ ($Y_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3- ($Y_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 50 | Pban | D2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 2 | 1 |
///
pub struct Sg49;

/// # 50 Pban (D₂ₕ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{19}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 50 | Pban | D2h^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 16 | P222 | D2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 32 | Pba2 | C2v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | X2 | $X_{2}$ | $k_{20}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{26}t1$ | 2 | B8a | yes |
/// | R2 | $R_{2}$ | $k_{26}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 22 | F222 | D2^7 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 22 | F222 | D2^7 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{23}t1$ | 2 | B8a | yes |
/// | T2 | $T_{2}$ | $k_{23}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{22}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{22}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{3}^+$ | $k_{22}t7$ | 1 | A2a | no |
/// | Z2- | $Z_{3}^-$ | $k_{22}t8$ | 1 | A2a | no |
/// | Z3+ | $Z_{4}^+$ | $k_{22}t3$ | 1 | A2a | no |
/// | Z3- | $Z_{4}^-$ | $k_{22}t4$ | 1 | A2a | no |
/// | Z4+ | $Z_{2}^+$ | $k_{22}t5$ | 1 | A2a | no |
/// | Z4- | $Z_{2}^-$ | $k_{22}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 50 | Pban | D2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 48 | Pnnn | D2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 48 | Pnnn | D2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 50 | Pban | D2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3+ ($Z_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3- ($Z_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{25}t1$ | 2 | B8a | yes |
/// | S2 | $S_{2}$ | $k_{25}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for S2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $U_{1}$ | $k_{24}t1$ | 2 | B8a | yes |
/// | U2 | $U_{2}$ | $k_{24}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for U1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for U2 ($U_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
/// | Y2 | $Y_{2}$ | $k_{21}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 16 | P222 | D2^1 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg50;

/// # 51 Pmma (D₂ₕ⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{4}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{4}^-$ | $k_{19}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{19}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $Z_{1}$ | $k_{22}t1$ | 2 | B8a | yes |
/// | X2 | $Z_{2}$ | $k_{22}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 25 | Pmm2 | C2v^1 | 4 | 1 |
/// | dir{(a,0)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ### Isotropy subgroups for X2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{26}t1$ | 2 | B8a | yes |
/// | R2 | $R_{2}$ | $k_{26}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 42 | Fmm2 | C2v^18 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 8 | 1 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 42 | Fmm2 | C2v^18 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1+ | $S_{1}^+$ | $k_{25}t1$ | 1 | A2a | no |
/// | T1- | $S_{1}^-$ | $k_{25}t2$ | 1 | A2a | no |
/// | T2+ | $S_{4}^+$ | $k_{25}t5$ | 1 | A2a | no |
/// | T2- | $S_{4}^-$ | $k_{25}t6$ | 1 | A2a | no |
/// | T3+ | $S_{3}^+$ | $k_{25}t7$ | 1 | A2a | no |
/// | T3- | $S_{3}^-$ | $k_{25}t8$ | 1 | A2a | no |
/// | T4+ | $S_{2}^+$ | $k_{25}t3$ | 1 | A2a | no |
/// | T4- | $S_{2}^-$ | $k_{25}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1+ ($S_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 2 | 1 |
///
/// ### Isotropy subgroups for T1- ($S_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for T2+ ($S_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for T2- ($S_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 2 | 1 |
///
/// ### Isotropy subgroups for T3+ ($S_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 2 | 1 |
///
/// ### Isotropy subgroups for T3- ($S_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for T4+ ($S_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for T4- ($S_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $X_{1}^+$ | $k_{21}t1$ | 1 | A2a | no |
/// | Z1- | $X_{1}^-$ | $k_{21}t2$ | 1 | A2a | no |
/// | Z2+ | $X_{4}^+$ | $k_{21}t5$ | 1 | A2a | no |
/// | Z2- | $X_{4}^-$ | $k_{21}t6$ | 1 | A2a | no |
/// | Z3+ | $X_{3}^+$ | $k_{21}t7$ | 1 | A2a | no |
/// | Z3- | $X_{3}^-$ | $k_{21}t8$ | 1 | A2a | no |
/// | Z4+ | $X_{2}^+$ | $k_{21}t3$ | 1 | A2a | no |
/// | Z4- | $X_{2}^-$ | $k_{21}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 57 | Pbcm | D2h^11 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 55 | Pbam | D2h^9 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 55 | Pbam | D2h^9 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 57 | Pbcm | D2h^11 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $T_{1}$ | $k_{24}t1$ | 2 | B8a | yes |
/// | S2 | $T_{2}$ | $k_{24}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for S1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 35 | Cmm2 | C2v^11 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 8 | 1 |
///
/// ### Isotropy subgroups for S2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 35 | Cmm2 | C2v^11 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $U_{1}$ | $k_{23}t1$ | 2 | B8a | yes |
/// | U2 | $U_{2}$ | $k_{23}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for U1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{(a,0)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ### Isotropy subgroups for U2 ($U_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Y_{1}^+$ | $k_{20}t1$ | 1 | A2a | no |
/// | Y1- | $Y_{1}^-$ | $k_{20}t2$ | 1 | A2a | no |
/// | Y2+ | $Y_{4}^+$ | $k_{20}t5$ | 1 | A2a | no |
/// | Y2- | $Y_{4}^-$ | $k_{20}t6$ | 1 | A2a | no |
/// | Y3+ | $Y_{3}^+$ | $k_{20}t7$ | 1 | A2a | no |
/// | Y3- | $Y_{3}^-$ | $k_{20}t8$ | 1 | A2a | no |
/// | Y4+ | $Y_{2}^+$ | $k_{20}t3$ | 1 | A2a | no |
/// | Y4- | $Y_{2}^-$ | $k_{20}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Y_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 57 | Pbcm | D2h^11 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Y_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 59 | Pmmn | D2h^13 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3+ ($Y_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 59 | Pmmn | D2h^13 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3- ($Y_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 57 | Pbcm | D2h^11 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
pub struct Sg51;

/// # 52 Pnna (D₂ₕ⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{19}t6$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 34 | Pnn2 | C2v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
/// | X2 | $X_{2}$ | $k_{21}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{26}t1$ | 2 | B8a | yes |
/// | R2 | $R_{2}$ | $k_{26}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 43 | Fdd2 | C2v^19 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 8 | 1 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 43 | Fdd2 | C2v^19 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 8 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1+ | $T_{1}^-$ | $k_{24}t1$ | 2 | B8a | no |
/// | T1- | $T_{1}^+$ | $k_{24}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for T1+ ($T_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for T1- ($T_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Y_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | Z2 | $Y_{2}$ | $k_{20}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1S2 | $U_{1}U2$ | $k_{23}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for S1S2 ($U_{1}U2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{C1(2)/4D} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{C12(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $S_{2}$ | $k_{25}t1$ | 2 | B8a | yes |
/// | U2 | $S_{1}$ | $k_{25}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for U1 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ### Isotropy subgroups for U2 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Z_{2}$ | $k_{22}t2$ | 2 | B8a | yes |
/// | Y2 | $Z_{1}$ | $k_{22}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 34 | Pnn2 | C2v^10 | 4 | 1 |
/// | dir{(a,a)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{(a,a)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
pub struct Sg52;

/// # 53 Pmna (D₂ₕ⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{19}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $Y_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | X2 | $Y_{2}$ | $k_{20}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{(a,0)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for X2 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{(a,0)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^-$ | $k_{26}t2$ | 2 | B8a | no |
/// | R1- | $R_{1}^+$ | $k_{26}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $U_{1}$ | $k_{23}t1$ | 2 | B8a | yes |
/// | T2 | $U_{2}$ | $k_{23}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 40 | Ama2 | C2v^16 | 4 | 1 |
/// | dir{(a,a)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ### Isotropy subgroups for T2 ($U_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 41 | Aba2 | C2v^17 | 4 | 1 |
/// | dir{(a,a)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{22}t1$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{22}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{(a,a)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 32 | Pba2 | C2v^8 | 4 | 1 |
/// | dir{(a,a)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{25}t1$ | 2 | B8a | yes |
/// | S2 | $S_{2}$ | $k_{25}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for S2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1+ | $T_{1}^-$ | $k_{24}t2$ | 2 | B8a | no |
/// | U1- | $T_{1}^+$ | $k_{24}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for U1+ ($T_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for U1- ($T_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $X_{1}^+$ | $k_{21}t1$ | 1 | A2a | no |
/// | Y1- | $X_{1}^-$ | $k_{21}t2$ | 1 | A2a | no |
/// | Y2+ | $X_{3}^+$ | $k_{21}t7$ | 1 | A2a | no |
/// | Y2- | $X_{3}^-$ | $k_{21}t8$ | 1 | A2a | no |
/// | Y3+ | $X_{2}^+$ | $k_{21}t3$ | 1 | A2a | no |
/// | Y3- | $X_{2}^-$ | $k_{21}t4$ | 1 | A2a | no |
/// | Y4+ | $X_{4}^+$ | $k_{21}t5$ | 1 | A2a | no |
/// | Y4- | $X_{4}^-$ | $k_{21}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 60 | Pbcn | D2h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 58 | Pnnm | D2h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 58 | Pnnm | D2h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 60 | Pbcn | D2h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
pub struct Sg53;

/// # 54 Pcca (D₂ₕ⁸)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{4}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{4}^-$ | $k_{19}t6$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{19}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 17 | P222_1 | D2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 27 | Pcc2 | C2v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 32 | Pba2 | C2v^8 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $Z_{1}$ | $k_{22}t1$ | 2 | B8a | yes |
/// | X2 | $Z_{2}$ | $k_{22}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 27 | Pcc2 | C2v^3 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for X2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{26}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{C1(2)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{C12(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $S_{1}$ | $k_{25}t1$ | 2 | B8a | yes |
/// | T2 | $S_{2}$ | $k_{25}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for T2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $X_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
/// | Z2 | $X_{2}$ | $k_{21}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 17 | P222_1 | D2^2 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $T_{1}$ | $k_{24}t1$ | 2 | B8a | yes |
/// | S2 | $T_{2}$ | $k_{24}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for S1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 37 | Ccc2 | C2v^13 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 8 | 1 |
///
/// ### Isotropy subgroups for S2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 37 | Ccc2 | C2v^13 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 8 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1U2 | $U_{1}U2$ | $k_{23}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for U1U2 ($U_{1}U2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{C11(2)/4D} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{C12(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Y_{1}^+$ | $k_{20}t1$ | 1 | A2a | no |
/// | Y1- | $Y_{1}^-$ | $k_{20}t2$ | 1 | A2a | no |
/// | Y2+ | $Y_{4}^+$ | $k_{20}t5$ | 1 | A2a | no |
/// | Y2- | $Y_{4}^-$ | $k_{20}t6$ | 1 | A2a | no |
/// | Y3+ | $Y_{3}^+$ | $k_{20}t7$ | 1 | A2a | no |
/// | Y3- | $Y_{3}^-$ | $k_{20}t8$ | 1 | A2a | no |
/// | Y4+ | $Y_{2}^+$ | $k_{20}t3$ | 1 | A2a | no |
/// | Y4- | $Y_{2}^-$ | $k_{20}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Y_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 60 | Pbcn | D2h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Y_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 56 | Pccn | D2h^10 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3+ ($Y_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 56 | Pccn | D2h^10 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3- ($Y_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 60 | Pbcn | D2h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
pub struct Sg54;

/// # 55 Pbam (D₂ₕ⁹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{19}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 55 | Pbam | D2h^9 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 32 | Pba2 | C2v^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | X2 | $X_{2}$ | $k_{20}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 26 | Pmc2_1 | C2v^2 | 4 | 1 |
/// | dir{(a,a)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 29 | Pca2_1 | C2v^5 | 4 | 1 |
/// | dir{(a,a)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+R2+ | $R_{1}^+R2+$ | $k_{26}t1t7$ | 2 | B4a | no |
/// | R1-R2- | $R_{1}^-R2-$ | $k_{26}t2t8$ | 2 | B4a | no |
/// | R3+R4+ | $R_{3}^+R4+$ | $k_{26}t3t5$ | 2 | B4a | no |
/// | R3-R4- | $R_{3}^-R4-$ | $k_{26}t4t6$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for R1+R2+ ($R_{1}^+R2+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 12 | C2/m | C2h^3 | 4 | 1 |
///
/// ### Isotropy subgroups for R1-R2- ($R_{1}^-R2-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 12 | C2/m | C2h^3 | 4 | 1 |
///
/// ### Isotropy subgroups for R3+R4+ ($R_{3}^+R4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 12 | C2/m | C2h^3 | 4 | 1 |
///
/// ### Isotropy subgroups for R3-R4- ($R_{3}^-R4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 12 | C2/m | C2h^3 | 4 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{2}$ | $k_{23}t2$ | 2 | B8a | yes |
/// | T2 | $T_{1}$ | $k_{23}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 8 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{22}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{22}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{3}^+$ | $k_{22}t7$ | 1 | A2a | no |
/// | Z2- | $Z_{3}^-$ | $k_{22}t8$ | 1 | A2a | no |
/// | Z3+ | $Z_{4}^+$ | $k_{22}t3$ | 1 | A2a | no |
/// | Z3- | $Z_{4}^-$ | $k_{22}t4$ | 1 | A2a | no |
/// | Z4+ | $Z_{2}^+$ | $k_{22}t5$ | 1 | A2a | no |
/// | Z4- | $Z_{2}^-$ | $k_{22}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 55 | Pbam | D2h^9 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 58 | Pnnm | D2h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 58 | Pnnm | D2h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 55 | Pbam | D2h^9 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3+ ($Z_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3- ($Z_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1+S2+ | $S_{1}^+S2+$ | $k_{25}t1t7$ | 2 | B4a | no |
/// | S1-S2- | $S_{1}^-S2-$ | $k_{25}t2t8$ | 2 | B4a | no |
/// | S3+S4+ | $S_{3}^+S4+$ | $k_{25}t3t5$ | 2 | B4a | no |
/// | S3-S4- | $S_{3}^-S4-$ | $k_{25}t4t6$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for S1+S2+ ($S_{1}^+S2+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 4 | 1 |
///
/// ### Isotropy subgroups for S1-S2- ($S_{1}^-S2-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 4 | 1 |
///
/// ### Isotropy subgroups for S3+S4+ ($S_{3}^+S4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 4 | 1 |
///
/// ### Isotropy subgroups for S3-S4- ($S_{3}^-S4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $U_{1}$ | $k_{24}t1$ | 2 | B8a | yes |
/// | U2 | $U_{2}$ | $k_{24}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for U1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 8 | 1 |
///
/// ### Isotropy subgroups for U2 ($U_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{2}$ | $k_{21}t2$ | 2 | B8a | yes |
/// | Y2 | $Y_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 29 | Pca2_1 | C2v^5 | 4 | 1 |
/// | dir{(a,a)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 26 | Pmc2_1 | C2v^2 | 4 | 1 |
/// | dir{(a,a)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
pub struct Sg55;

/// # 56 Pccn (D₂ₕ¹⁰)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{19}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 56 | Pccn | D2h^10 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 27 | Pcc2 | C2v^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | X2 | $X_{2}$ | $k_{20}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 27 | Pcc2 | C2v^3 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+R2+ | $R_{1}^-R2-$ | $k_{26}t2t8$ | 2 | B4a | no |
/// | R1-R2- | $R_{1}^+R2+$ | $k_{26}t1t7$ | 2 | B4a | no |
/// | R3+R4+ | $R_{3}^-R4-$ | $k_{26}t4t6$ | 2 | B4a | no |
/// | R3-R4- | $R_{3}^+R4+$ | $k_{26}t3t5$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for R1+R2+ ($R_{1}^-R2-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 15 | C2/c | C2h^6 | 4 | 1 |
///
/// ### Isotropy subgroups for R1-R2- ($R_{1}^+R2+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 15 | C2/c | C2h^6 | 4 | 1 |
///
/// ### Isotropy subgroups for R3+R4+ ($R_{3}^-R4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 15 | C2/c | C2h^6 | 4 | 1 |
///
/// ### Isotropy subgroups for R3-R4- ($R_{3}^+R4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 15 | C2/c | C2h^6 | 4 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1T2 | $T_{1}T2$ | $k_{23}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for T1T2 ($T_{1}T2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{C11(2)/4D} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{C12(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{22}t1$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{22}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 19 | P2_12_12_1 | D2^4 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{25}t1$ | 2 | B8a | no |
/// | S2 | $S_{2}$ | $k_{25}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 37 | Ccc2 | C2v^13 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for S2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 37 | Ccc2 | C2v^13 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1U2 | $U_{1}U2$ | $k_{24}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for U1U2 ($U_{1}U2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{C11(2)/4D} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{C12(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
/// | Y2 | $Y_{2}$ | $k_{21}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 27 | Pcc2 | C2v^3 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
pub struct Sg56;

/// # 57 Pbcm (D₂ₕ¹¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{19}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 57 | Pbcm | D2h^11 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 28 | Pma2 | C2v^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $Z_{1}^+$ | $k_{22}t1$ | 1 | A2a | no |
/// | X1- | $Z_{1}^-$ | $k_{22}t2$ | 1 | A2a | no |
/// | X2+ | $Z_{2}^+$ | $k_{22}t3$ | 1 | A2a | no |
/// | X2- | $Z_{2}^-$ | $k_{22}t4$ | 1 | A2a | no |
/// | X3+ | $Z_{3}^+$ | $k_{22}t7$ | 1 | A2a | no |
/// | X3- | $Z_{3}^-$ | $k_{22}t8$ | 1 | A2a | no |
/// | X4+ | $Z_{4}^+$ | $k_{22}t5$ | 1 | A2a | no |
/// | X4- | $Z_{4}^-$ | $k_{22}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 57 | Pbcm | D2h^11 | 2 | 1 |
///
/// ### Isotropy subgroups for X1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 60 | Pbcn | D2h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for X2+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for X2- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 61 | Pbca | D2h^15 | 2 | 1 |
///
/// ### Isotropy subgroups for X3+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 60 | Pbcn | D2h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for X3- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 57 | Pbcm | D2h^11 | 2 | 1 |
///
/// ### Isotropy subgroups for X4+ ($Z_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 61 | Pbca | D2h^15 | 2 | 1 |
///
/// ### Isotropy subgroups for X4- ($Z_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{26}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C12(2)/4D} | 8 | Cm | Cs^3 | 8 | 1 |
/// | dir{C1(2)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{C11(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1T2 | $S_{1}S2$ | $k_{25}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for T1T2 ($S_{1}S2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C12(2)/4D} | 8 | Cm | Cs^3 | 8 | 1 |
/// | dir{C1(2)/4D} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{C11(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Y_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | Z2 | $Y_{2}$ | $k_{20}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 26 | Pmc2_1 | C2v^2 | 4 | 1 |
/// | dir{(a,a)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 31 | Pmn2_1 | C2v^7 | 4 | 1 |
/// | dir{(a,a)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $U_{2}$ | $k_{23}t2$ | 2 | B8a | yes |
/// | S2 | $U_{1}$ | $k_{23}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for S1 ($U_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 41 | Aba2 | C2v^17 | 4 | 1 |
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for S2 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 40 | Ama2 | C2v^16 | 4 | 1 |
/// | dir{(a,a)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $T_{1}$ | $k_{24}t1$ | 2 | B8a | yes |
/// | U2 | $T_{2}$ | $k_{24}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for U1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 8 | 1 |
///
/// ### Isotropy subgroups for U2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $X_{2}$ | $k_{21}t2$ | 2 | B8a | yes |
/// | Y2 | $X_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 32 | Pba2 | C2v^8 | 4 | 1 |
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for Y2 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{(a,a)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
pub struct Sg57;

/// # 58 Pnnm (D₂ₕ¹²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{19}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 58 | Pnnm | D2h^12 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 34 | Pnn2 | C2v^10 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | X2 | $X_{2}$ | $k_{20}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 31 | Pmn2_1 | C2v^7 | 4 | 1 |
/// | dir{(a,a)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 33 | Pna2_1 | C2v^9 | 4 | 1 |
/// | dir{(a,a)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{26}t1$ | 2 | B8a | no |
/// | R2 | $R_{2}$ | $k_{26}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 43 | Fdd2 | C2v^19 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 43 | Fdd2 | C2v^19 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1+ | $T_{1}^+$ | $k_{23}t2$ | 2 | B8a | no |
/// | T1- | $T_{1}^-$ | $k_{23}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for T1+ ($T_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for T1- ($T_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{22}t1$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{22}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{(a,0)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 19 | P2_12_12_1 | D2^4 | 4 | 1 |
/// | dir{(a,0)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1+S2+ | $S_{1}^+S2+$ | $k_{25}t1t7$ | 2 | B4a | no |
/// | S1-S2- | $S_{1}^-S2-$ | $k_{25}t2t8$ | 2 | B4a | no |
/// | S3+S4+ | $S_{3}^+S4+$ | $k_{25}t3t5$ | 2 | B4a | no |
/// | S3-S4- | $S_{3}^-S4-$ | $k_{25}t4t6$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for S1+S2+ ($S_{1}^+S2+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 4 | 1 |
///
/// ### Isotropy subgroups for S1-S2- ($S_{1}^-S2-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 4 | 1 |
///
/// ### Isotropy subgroups for S3+S4+ ($S_{3}^+S4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 4 | 1 |
///
/// ### Isotropy subgroups for S3-S4- ($S_{3}^-S4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1+ | $U_{1}^+$ | $k_{24}t2$ | 2 | B8a | no |
/// | U1- | $U_{1}^-$ | $k_{24}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for U1+ ($U_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for U1- ($U_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{2}$ | $k_{21}t2$ | 2 | B8a | yes |
/// | Y2 | $Y_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 33 | Pna2_1 | C2v^9 | 4 | 1 |
/// | dir{(a,a)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 31 | Pmn2_1 | C2v^7 | 4 | 1 |
/// | dir{(a,a)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
pub struct Sg58;

/// # 59 Pmmn (D₂ₕ¹³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{19}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 59 | Pmmn | D2h^13 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 25 | Pmm2 | C2v^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | X2 | $X_{2}$ | $k_{20}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 25 | Pmm2 | C2v^1 | 4 | 1 |
/// | dir{(a,0)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{26}t1$ | 2 | B8a | no |
/// | R2 | $R_{2}$ | $k_{26}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 42 | Fmm2 | C2v^18 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for R2 ($R_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 42 | Fmm2 | C2v^18 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{23}t1$ | 2 | B8a | yes |
/// | T2 | $T_{2}$ | $k_{23}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{(a,0)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{22}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{22}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{3}^+$ | $k_{22}t7$ | 1 | A2a | no |
/// | Z2- | $Z_{3}^-$ | $k_{22}t8$ | 1 | A2a | no |
/// | Z3+ | $Z_{4}^+$ | $k_{22}t3$ | 1 | A2a | no |
/// | Z3- | $Z_{4}^-$ | $k_{22}t4$ | 1 | A2a | no |
/// | Z4+ | $Z_{2}^+$ | $k_{22}t5$ | 1 | A2a | no |
/// | Z4- | $Z_{2}^-$ | $k_{22}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 59 | Pmmn | D2h^13 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 56 | Pccn | D2h^10 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 56 | Pccn | D2h^10 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 59 | Pmmn | D2h^13 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3+ ($Z_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3- ($Z_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{25}t1$ | 2 | B8a | no |
/// | S2 | $S_{2}$ | $k_{25}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 35 | Cmm2 | C2v^11 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for S2 ($S_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 35 | Cmm2 | C2v^11 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $U_{1}$ | $k_{24}t1$ | 2 | B8a | yes |
/// | U2 | $U_{2}$ | $k_{24}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for U1 ($U_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{(a,0)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ### Isotropy subgroups for U2 ($U_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
/// | Y2 | $Y_{2}$ | $k_{21}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 25 | Pmm2 | C2v^1 | 4 | 1 |
/// | dir{(a,0)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 28 | Pma2 | C2v^4 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
pub struct Sg59;

/// # 60 Pbcn (D₂ₕ¹⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{19}t6$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 60 | Pbcn | D2h^14 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 18 | P2_12_12 | D2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 30 | Pnc2 | C2v^6 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
/// | X2 | $X_{2}$ | $k_{21}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 29 | Pca2_1 | C2v^5 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 33 | Pna2_1 | C2v^9 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{26}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C11(2)/4D} | 9 | Cc | Cs^4 | 8 | 1 |
/// | dir{C1(2)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{C12(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1T2 | $T_{1}T2$ | $k_{24}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for T1T2 ($T_{1}T2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C12(2)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{C11(2)/4D} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{C1(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Y_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | Z2 | $Y_{2}$ | $k_{20}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 30 | Pnc2 | C2v^6 | 4 | 1 |
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 34 | Pnn2 | C2v^10 | 4 | 1 |
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1+ | $U_{1}^+$ | $k_{23}t1$ | 2 | B8a | no |
/// | S1- | $U_{1}^-$ | $k_{23}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1+ ($U_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for S1- ($U_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1U2 | $S_{1}S2$ | $k_{25}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for U1U2 ($S_{1}S2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C11(2)/4D} | 9 | Cc | Cs^4 | 8 | 1 |
/// | dir{C1(2)/4D} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{C12(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Z_{1}$ | $k_{22}t1$ | 2 | B8a | yes |
/// | Y2 | $Z_{2}$ | $k_{22}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 18 | P2_12_12 | D2^3 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 19 | P2_12_12_1 | D2^4 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg60;

/// # 61 Pbca (D₂ₕ¹⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{2}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{2}^-$ | $k_{19}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{4}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{4}^-$ | $k_{19}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 61 | Pbca | D2h^15 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 19 | P2_12_12_1 | D2^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 29 | Pca2_1 | C2v^5 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $Y_{1}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | X2 | $Y_{2}$ | $k_{20}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 29 | Pca2_1 | C2v^5 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for X2 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 33 | Pna2_1 | C2v^9 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+R1+ | $R_{1}^+R1+$ | $k_{26}t1t1$ | 4 | D8a | no |
/// | R1-R1- | $R_{1}^-R1-$ | $k_{26}t2t2$ | 4 | D8a | no |
///
/// ### Isotropy subgroups for R1+R1+ ($R_{1}^+R1+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{4D1(4)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ### Isotropy subgroups for R1-R1- ($R_{1}^-R1-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{4D1(4)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1T2 | $U_{1}U2$ | $k_{23}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for T1T2 ($U_{1}U2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C12(2)/4D} | 9 | Cc | Cs^4 | 8 | 1 |
/// | dir{C1(2)/4D} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{C11(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{22}t1$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{22}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 29 | Pca2_1 | C2v^5 | 4 | 1 |
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 33 | Pna2_1 | C2v^9 | 4 | 1 |
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1S2 | $S_{1}S2$ | $k_{25}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for S1S2 ($S_{1}S2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C11(2)/4D} | 9 | Cc | Cs^4 | 8 | 1 |
/// | dir{C1(2)/4D} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{C12(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1U2 | $T_{1}T2$ | $k_{24}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for U1U2 ($T_{1}T2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C11(2)/4D} | 9 | Cc | Cs^4 | 8 | 1 |
/// | dir{C1(2)/4D} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{C12(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $X_{2}$ | $k_{21}t2$ | 2 | B8a | yes |
/// | Y2 | $X_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 33 | Pna2_1 | C2v^9 | 4 | 1 |
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for Y2 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 29 | Pca2_1 | C2v^5 | 4 | 1 |
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
pub struct Sg61;

/// # 62 Pnma (D₂ₕ¹⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{19}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{19}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{19}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{19}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{19}t5$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{19}t6$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{19}t3$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{19}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 19 | P2_12_12_1 | D2^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 31 | Pmn2_1 | C2v^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 26 | Pmc2_1 | C2v^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 33 | Pna2_1 | C2v^9 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $X_{1}$ | $k_{21}t1$ | 2 | B8a | yes |
/// | X2 | $X_{2}$ | $k_{21}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for X1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 31 | Pmn2_1 | C2v^7 | 4 | 1 |
/// | dir{(a,0)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ### Isotropy subgroups for X2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 33 | Pna2_1 | C2v^9 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1R2 | $R_{1}R2$ | $k_{26}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for R1R2 ($R_{1}R2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C11(2)/4D} | 9 | Cc | Cs^4 | 8 | 1 |
/// | dir{C12(2)/4D} | 8 | Cm | Cs^3 | 8 | 1 |
/// | dir{C1(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{2}$ | $k_{24}t1$ | 2 | B8a | no |
/// | T2 | $T_{1}$ | $k_{24}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for T1 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 36 | Cmc2_1 | C2v^12 | 4 | 1 |
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{22}t1$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{22}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 26 | Pmc2_1 | C2v^2 | 4 | 1 |
/// | dir{(a,0)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{(a,b)} | 6 | Pm | Cs^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 29 | Pca2_1 | C2v^5 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1S2 | $S_{1}S2$ | $k_{25}t1t2$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for S1S2 ($S_{1}S2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C11(2)/4D} | 8 | Cm | Cs^3 | 8 | 1 |
/// | dir{C1(2)/4D} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{C12(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 16 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1+U4+ | $U_{1}^+U2+$ | $k_{23}t1t3$ | 2 | B4a | no |
/// | U1-U4- | $U_{1}^-U2-$ | $k_{23}t2t4$ | 2 | B4a | no |
/// | U2+U3+ | $U_{3}^+U4+$ | $k_{23}t5t7$ | 2 | B4a | no |
/// | U2-U3- | $U_{3}^-U4-$ | $k_{23}t6t8$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for U1+U4+ ($U_{1}^+U2+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
///
/// ### Isotropy subgroups for U1-U4- ($U_{1}^-U2-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
///
/// ### Isotropy subgroups for U2+U3+ ($U_{3}^+U4+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
///
/// ### Isotropy subgroups for U2-U3- ($U_{3}^-U4-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{2}$ | $k_{20}t1$ | 2 | B8a | yes |
/// | Y2 | $Y_{1}$ | $k_{20}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 31 | Pmn2_1 | C2v^7 | 4 | 1 |
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 26 | Pmc2_1 | C2v^2 | 4 | 1 |
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 7 | Pc | Cs^2 | 8 | 1 |
///
pub struct Sg62;

/// # 63 Cmcm (D₂ₕ¹⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{14}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{14}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{14}t5$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{14}t6$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{14}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 36 | Cmc2_1 | C2v^12 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 40 | Ama2 | C2v^16 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 38 | Amm2 | C2v^14 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{13}t1$ | 4 | D32c | yes |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P12(1)/4D} | 8 | Cm | Cs^3 | 8 | 1 |
/// | dir{P5(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{C13(2)/4D} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P3(1)/4D} | 43 | Fdd2 | C2v^19 | 8 | 2 |
/// | dir{P4(1)/4D} | 42 | Fmm2 | C2v^18 | 8 | 2 |
/// | dir{P1(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P11(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{C9(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C1(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C11(2)/4D} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{C10(2)/4D} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{C2(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C8(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C3(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C12(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{17}t1$ | 2 | B8a | yes |
/// | T2 | $T_{2}$ | $k_{17}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 44 | Imm2 | C2v^20 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 8 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 46 | Ima2 | C2v^22 | 4 | 1 |
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{16}t1$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{16}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 38 | Amm2 | C2v^14 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 40 | Ama2 | C2v^16 | 4 | 1 |
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1+ | $S_{1}^+$ | $k_{12}t1$ | 2 | B8a | no |
/// | S1- | $S_{1}^-$ | $k_{12}t2$ | 2 | B8a | no |
/// | S2+ | $S_{2}^+$ | $k_{12}t3$ | 2 | B8a | no |
/// | S2- | $S_{2}^-$ | $k_{12}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1+ ($S_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{(a,0)} | 63 | Cmcm | D2h^17 | 4 | 2 |
/// | dir{(a,b)} | 11 | P2_1/m | C2h^2 | 8 | 2 |
///
/// ### Isotropy subgroups for S1- ($S_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,0)} | 64 | Cmca | D2h^18 | 4 | 2 |
/// | dir{(a,b)} | 14 | P2_1/c | C2h^5 | 8 | 2 |
///
/// ### Isotropy subgroups for S2+ ($S_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,0)} | 64 | Cmca | D2h^18 | 4 | 2 |
/// | dir{(a,b)} | 14 | P2_1/c | C2h^5 | 8 | 2 |
///
/// ### Isotropy subgroups for S2- ($S_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{(a,0)} | 63 | Cmcm | D2h^17 | 4 | 2 |
/// | dir{(a,b)} | 11 | P2_1/m | C2h^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Y_{1}^+$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y1- | $Y_{1}^-$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y2+ | $Y_{3}^+$ | $k_{15}t7$ | 1 | A2a | no |
/// | Y2- | $Y_{3}^-$ | $k_{15}t8$ | 1 | A2a | no |
/// | Y3+ | $Y_{4}^+$ | $k_{15}t5$ | 1 | A2a | no |
/// | Y3- | $Y_{4}^-$ | $k_{15}t6$ | 1 | A2a | no |
/// | Y4+ | $Y_{2}^+$ | $k_{15}t3$ | 1 | A2a | no |
/// | Y4- | $Y_{2}^-$ | $k_{15}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Y_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Y_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3+ ($Y_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 58 | Pnnm | D2h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3- ($Y_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 57 | Pbcm | D2h^11 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 60 | Pbcn | D2h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 59 | Pmmn | D2h^13 | 2 | 1 |
///
pub struct Sg63;

/// # 64 Cmca (D₂ₕ¹⁸)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{14}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{14}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{14}t5$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{14}t6$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{14}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 20 | C222_1 | D2^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 36 | Cmc2_1 | C2v^12 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 41 | Aba2 | C2v^17 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 39 | Abm2 | C2v^15 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+R2+ | $R_{1}^+R2+$ | $k_{13}t1t3$ | 4 | D16c | no |
/// | R1-R2- | $R_{1}^-R2-$ | $k_{13}t2t4$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for R1+R2+ ($R_{1}^+R2+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{C12(2)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{C11(2)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{4D1(4)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
///
/// ### Isotropy subgroups for R1-R2- ($R_{1}^-R2-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{C12(2)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{C11(2)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{4D1(4)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{17}t2$ | 2 | B8a | yes |
/// | T2 | $T_{2}$ | $k_{17}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 46 | Ima2 | C2v^22 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 8 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 45 | Iba2 | C2v^21 | 4 | 1 |
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{16}t1$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{16}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 39 | Abm2 | C2v^15 | 4 | 1 |
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 8 | Cm | Cs^3 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 41 | Aba2 | C2v^17 | 4 | 1 |
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 9 | Cc | Cs^4 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{12}t1$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 4 | P2_1 | C2^2 | 8 | 1 |
/// | dir{P3(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P5(1)/4D} | 36 | Cmc2_1 | C2v^12 | 8 | 2 |
/// | dir{P4(1)/4D} | 20 | C222_1 | D2^5 | 8 | 2 |
/// | dir{P12(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P11(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{C13(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C12(2)/4D} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C10(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C3(2)/4D} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{C2(2)/4D} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{C9(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C8(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Y_{1}^+$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y1- | $Y_{1}^-$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y2+ | $Y_{3}^+$ | $k_{15}t8$ | 1 | A2a | no |
/// | Y2- | $Y_{3}^-$ | $k_{15}t7$ | 1 | A2a | no |
/// | Y3+ | $Y_{4}^+$ | $k_{15}t6$ | 1 | A2a | no |
/// | Y3- | $Y_{4}^-$ | $k_{15}t5$ | 1 | A2a | no |
/// | Y4+ | $Y_{2}^+$ | $k_{15}t4$ | 1 | A2a | no |
/// | Y4- | $Y_{2}^-$ | $k_{15}t3$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Y_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 61 | Pbca | D2h^15 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Y_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3+ ($Y_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 55 | Pbam | D2h^9 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3- ($Y_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 60 | Pbcn | D2h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 56 | Pccn | D2h^10 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 57 | Pbcm | D2h^11 | 2 | 1 |
///
pub struct Sg64;

/// # 65 Cmmm (D₂ₕ¹⁹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{14}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{14}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{14}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{14}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 35 | Cmm2 | C2v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 38 | Amm2 | C2v^14 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 38 | Amm2 | C2v^14 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{13}t1$ | 2 | B8a | no |
/// | R1- | $R_{1}^-$ | $k_{13}t2$ | 2 | B8a | no |
/// | R2+ | $R_{2}^+$ | $k_{13}t3$ | 2 | B8a | no |
/// | R2- | $R_{2}^-$ | $k_{13}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 69 | Fmmm | D2h^23 | 4 | 2 |
/// | dir{(a,b)} | 12 | C2/m | C2h^3 | 8 | 2 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 69 | Fmmm | D2h^23 | 4 | 2 |
/// | dir{(a,b)} | 12 | C2/m | C2h^3 | 8 | 2 |
///
/// ### Isotropy subgroups for R2+ ($R_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 69 | Fmmm | D2h^23 | 4 | 2 |
/// | dir{(a,b)} | 12 | C2/m | C2h^3 | 8 | 2 |
///
/// ### Isotropy subgroups for R2- ($R_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 69 | Fmmm | D2h^23 | 4 | 2 |
/// | dir{(a,b)} | 12 | C2/m | C2h^3 | 8 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1+ | $T_{1}^+$ | $k_{17}t1$ | 1 | A2a | no |
/// | T1- | $T_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | T2+ | $T_{3}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | T2- | $T_{3}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | T3+ | $T_{4}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | T3- | $T_{4}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | T4+ | $T_{2}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | T4- | $T_{2}^-$ | $k_{17}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1+ ($T_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 71 | Immm | D2h^25 | 2 | 1 |
///
/// ### Isotropy subgroups for T1- ($T_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 72 | Ibam | D2h^26 | 2 | 1 |
///
/// ### Isotropy subgroups for T2+ ($T_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 72 | Ibam | D2h^26 | 2 | 1 |
///
/// ### Isotropy subgroups for T2- ($T_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 71 | Immm | D2h^25 | 2 | 1 |
///
/// ### Isotropy subgroups for T3+ ($T_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 74 | Imma | D2h^28 | 2 | 1 |
///
/// ### Isotropy subgroups for T3- ($T_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 74 | Imma | D2h^28 | 2 | 1 |
///
/// ### Isotropy subgroups for T4+ ($T_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 74 | Imma | D2h^28 | 2 | 1 |
///
/// ### Isotropy subgroups for T4- ($T_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 74 | Imma | D2h^28 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{16}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{16}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{3}^+$ | $k_{16}t7$ | 1 | A2a | no |
/// | Z2- | $Z_{3}^-$ | $k_{16}t8$ | 1 | A2a | no |
/// | Z3+ | $Z_{4}^+$ | $k_{16}t3$ | 1 | A2a | no |
/// | Z3- | $Z_{4}^-$ | $k_{16}t4$ | 1 | A2a | no |
/// | Z4+ | $Z_{2}^+$ | $k_{16}t5$ | 1 | A2a | no |
/// | Z4- | $Z_{2}^-$ | $k_{16}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 66 | Cccm | D2h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 66 | Cccm | D2h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3+ ($Z_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3- ($Z_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1+ | $S_{1}^+$ | $k_{12}t1$ | 2 | B8a | no |
/// | S1- | $S_{1}^-$ | $k_{12}t2$ | 2 | B8a | no |
/// | S2+ | $S_{2}^+$ | $k_{12}t3$ | 2 | B8a | no |
/// | S2- | $S_{2}^-$ | $k_{12}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1+ ($S_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,0)} | 65 | Cmmm | D2h^19 | 4 | 2 |
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ### Isotropy subgroups for S1- ($S_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,0)} | 67 | Cmma | D2h^21 | 4 | 2 |
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for S2+ ($S_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,0)} | 67 | Cmma | D2h^21 | 4 | 2 |
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for S2- ($S_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,0)} | 65 | Cmmm | D2h^19 | 4 | 2 |
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Y_{1}^+$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y1- | $Y_{1}^-$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y2+ | $Y_{3}^+$ | $k_{15}t7$ | 1 | A2a | no |
/// | Y2- | $Y_{3}^-$ | $k_{15}t8$ | 1 | A2a | no |
/// | Y3+ | $Y_{4}^+$ | $k_{15}t3$ | 1 | A2a | no |
/// | Y3- | $Y_{4}^-$ | $k_{15}t4$ | 1 | A2a | no |
/// | Y4+ | $Y_{2}^+$ | $k_{15}t5$ | 1 | A2a | no |
/// | Y4- | $Y_{2}^-$ | $k_{15}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 47 | Pmmm | D2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 50 | Pban | D2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Y_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 55 | Pbam | D2h^9 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Y_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 59 | Pmmn | D2h^13 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3+ ($Y_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3- ($Y_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
pub struct Sg65;

/// # 66 Cccm (D₂ₕ²⁰)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{14}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{14}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{14}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{14}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 66 | Cccm | D2h^20 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 37 | Ccc2 | C2v^13 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 40 | Ama2 | C2v^16 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 40 | Ama2 | C2v^16 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{13}t1$ | 2 | B8a | no |
/// | R1- | $R_{1}^-$ | $k_{13}t2$ | 2 | B8a | no |
/// | R2+ | $R_{2}^+$ | $k_{13}t3$ | 2 | B8a | no |
/// | R2- | $R_{2}^-$ | $k_{13}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 70 | Fddd | D2h^24 | 4 | 2 |
/// | dir{(a,b)} | 15 | C2/c | C2h^6 | 8 | 2 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 70 | Fddd | D2h^24 | 4 | 2 |
/// | dir{(a,b)} | 15 | C2/c | C2h^6 | 8 | 2 |
///
/// ### Isotropy subgroups for R2+ ($R_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 70 | Fddd | D2h^24 | 4 | 2 |
/// | dir{(a,b)} | 15 | C2/c | C2h^6 | 8 | 2 |
///
/// ### Isotropy subgroups for R2- ($R_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 70 | Fddd | D2h^24 | 4 | 2 |
/// | dir{(a,b)} | 15 | C2/c | C2h^6 | 8 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{17}t1$ | 2 | B8a | yes |
/// | T2 | $T_{2}$ | $k_{17}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 23 | I222 | D2^8 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 24 | I2_12_12_1 | D2^9 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{16}t1$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{16}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{(a,0)} | 11 | P2_1/m | C2h^2 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1+ | $S_{1}^+$ | $k_{12}t1$ | 2 | B8a | no |
/// | S1- | $S_{1}^-$ | $k_{12}t2$ | 2 | B8a | no |
/// | S2+ | $S_{2}^+$ | $k_{12}t3$ | 2 | B8a | no |
/// | S2- | $S_{2}^-$ | $k_{12}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1+ ($S_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,0)} | 66 | Cccm | D2h^20 | 4 | 2 |
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ### Isotropy subgroups for S1- ($S_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,0)} | 68 | Ccca | D2h^22 | 4 | 2 |
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for S2+ ($S_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,0)} | 68 | Ccca | D2h^22 | 4 | 2 |
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for S2- ($S_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 10 | P2/m | C2h^1 | 4 | 1 |
/// | dir{(a,0)} | 66 | Cccm | D2h^20 | 4 | 2 |
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Y_{1}^+$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y1- | $Y_{1}^-$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y2+ | $Y_{3}^+$ | $k_{15}t7$ | 1 | A2a | no |
/// | Y2- | $Y_{3}^-$ | $k_{15}t8$ | 1 | A2a | no |
/// | Y3+ | $Y_{4}^+$ | $k_{15}t3$ | 1 | A2a | no |
/// | Y3- | $Y_{4}^-$ | $k_{15}t4$ | 1 | A2a | no |
/// | Y4+ | $Y_{2}^+$ | $k_{15}t5$ | 1 | A2a | no |
/// | Y4- | $Y_{2}^-$ | $k_{15}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 48 | Pnnn | D2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Y_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 58 | Pnnm | D2h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Y_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 56 | Pccn | D2h^10 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3+ ($Y_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3- ($Y_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
pub struct Sg66;

/// # 67 Cmma (D₂ₕ²¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{14}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{14}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{14}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{14}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 35 | Cmm2 | C2v^11 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 39 | Abm2 | C2v^15 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 39 | Abm2 | C2v^15 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{13}t1$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P3(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P5(1)/4D} | 42 | Fmm2 | C2v^18 | 8 | 2 |
/// | dir{P4(1)/4D} | 22 | F222 | D2^7 | 8 | 2 |
/// | dir{P12(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{P11(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{C13(2)/4D} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{C12(2)/4D} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{C3(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C2(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C10(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C9(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C8(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1+ | $T_{1}^+$ | $k_{17}t1$ | 1 | A2a | no |
/// | T1- | $T_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | T2+ | $T_{3}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | T2- | $T_{3}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | T3+ | $T_{4}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | T3- | $T_{4}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | T4+ | $T_{2}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | T4- | $T_{2}^-$ | $k_{17}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1+ ($T_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 72 | Ibam | D2h^26 | 2 | 1 |
///
/// ### Isotropy subgroups for T1- ($T_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 72 | Ibam | D2h^26 | 2 | 1 |
///
/// ### Isotropy subgroups for T2+ ($T_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 72 | Ibam | D2h^26 | 2 | 1 |
///
/// ### Isotropy subgroups for T2- ($T_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 72 | Ibam | D2h^26 | 2 | 1 |
///
/// ### Isotropy subgroups for T3+ ($T_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 74 | Imma | D2h^28 | 2 | 1 |
///
/// ### Isotropy subgroups for T3- ($T_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 73 | Ibca | D2h^27 | 2 | 1 |
///
/// ### Isotropy subgroups for T4+ ($T_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 73 | Ibca | D2h^27 | 2 | 1 |
///
/// ### Isotropy subgroups for T4- ($T_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 74 | Imma | D2h^28 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{16}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{16}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{3}^+$ | $k_{16}t7$ | 1 | A2a | no |
/// | Z2- | $Z_{3}^-$ | $k_{16}t8$ | 1 | A2a | no |
/// | Z3+ | $Z_{4}^+$ | $k_{16}t3$ | 1 | A2a | no |
/// | Z3- | $Z_{4}^-$ | $k_{16}t4$ | 1 | A2a | no |
/// | Z4+ | $Z_{2}^+$ | $k_{16}t5$ | 1 | A2a | no |
/// | Z4- | $Z_{2}^-$ | $k_{16}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 68 | Ccca | D2h^22 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 68 | Ccca | D2h^22 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3+ ($Z_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3- ($Z_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 2 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{12}t1$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P3(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P5(1)/4D} | 35 | Cmm2 | C2v^11 | 8 | 2 |
/// | dir{P4(1)/4D} | 21 | C222 | D2^6 | 8 | 2 |
/// | dir{P12(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{P11(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{C13(2)/4D} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{C12(2)/4D} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C10(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C3(2)/4D} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{C2(2)/4D} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{C9(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C8(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Y_{1}^+$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y1- | $Y_{1}^-$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y2+ | $Y_{3}^+$ | $k_{15}t7$ | 1 | A2a | no |
/// | Y2- | $Y_{3}^-$ | $k_{15}t8$ | 1 | A2a | no |
/// | Y3+ | $Y_{4}^+$ | $k_{15}t3$ | 1 | A2a | no |
/// | Y3- | $Y_{4}^-$ | $k_{15}t4$ | 1 | A2a | no |
/// | Y4+ | $Y_{2}^+$ | $k_{15}t5$ | 1 | A2a | no |
/// | Y4- | $Y_{2}^-$ | $k_{15}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Y_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 57 | Pbcm | D2h^11 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Y_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 57 | Pbcm | D2h^11 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3+ ($Y_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3- ($Y_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
pub struct Sg67;

/// # 68 Ccca (D₂ₕ²²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{14}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{14}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{14}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{14}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 68 | Ccca | D2h^22 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 21 | C222 | D2^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 37 | Ccc2 | C2v^13 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 41 | Aba2 | C2v^17 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 41 | Aba2 | C2v^17 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{13}t1$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P3(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P5(1)/4D} | 43 | Fdd2 | C2v^19 | 8 | 2 |
/// | dir{P4(1)/4D} | 22 | F222 | D2^7 | 8 | 2 |
/// | dir{P11(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P12(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{C12(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C13(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C3(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C2(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C10(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C9(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C8(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{17}t1$ | 2 | B8a | yes |
/// | T2 | $T_{2}$ | $k_{17}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 23 | I222 | D2^8 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ### Isotropy subgroups for T2 ($T_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 24 | I2_12_12_1 | D2^9 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,b)} | 5 | C2 | C2^3 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{16}t1$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{16}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{12}t1$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 3 | P2 | C2^1 | 8 | 1 |
/// | dir{P3(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P5(1)/4D} | 37 | Ccc2 | C2v^13 | 8 | 2 |
/// | dir{P4(1)/4D} | 21 | C222 | D2^6 | 8 | 2 |
/// | dir{P11(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P12(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{C12(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C13(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C10(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C3(2)/4D} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{C2(2)/4D} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{C8(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C9(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Y_{1}^+$ | $k_{15}t1$ | 1 | A2a | no |
/// | Y1- | $Y_{1}^-$ | $k_{15}t2$ | 1 | A2a | no |
/// | Y2+ | $Y_{3}^+$ | $k_{15}t7$ | 1 | A2a | no |
/// | Y2- | $Y_{3}^-$ | $k_{15}t8$ | 1 | A2a | no |
/// | Y3+ | $Y_{4}^+$ | $k_{15}t3$ | 1 | A2a | no |
/// | Y3- | $Y_{4}^-$ | $k_{15}t4$ | 1 | A2a | no |
/// | Y4+ | $Y_{2}^+$ | $k_{15}t5$ | 1 | A2a | no |
/// | Y4- | $Y_{2}^-$ | $k_{15}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 50 | Pban | D2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 50 | Pban | D2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Y_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 60 | Pbcn | D2h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Y_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 60 | Pbcn | D2h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3+ ($Y_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3- ($Y_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
pub struct Sg68;

/// # 69 Fmmm (D₂ₕ²³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{14}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{14}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{14}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{14}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 69 | Fmmm | D2h^23 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 22 | F222 | D2^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 42 | Fmm2 | C2v^18 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 42 | Fmm2 | C2v^18 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 42 | Fmm2 | C2v^18 | 2 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^+$ | $k_{10}t1$ | 4 | D32c | no |
/// | L1- | $L_{1}^-$ | $k_{10}t2$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{P3(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{P5(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{P4(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{C1(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C3(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C2(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P12(1)/4D} | 70 | Fddd | D2h^24 | 8 | 4 |
/// | dir{P11(1)/4D} | 69 | Fmmm | D2h^23 | 8 | 4 |
/// | dir{C9(2)/4D} | 15 | C2/c | C2h^6 | 16 | 4 |
/// | dir{C13(2)/4D} | 15 | C2/c | C2h^6 | 16 | 4 |
/// | dir{C11(2)/4D} | 15 | C2/c | C2h^6 | 16 | 4 |
/// | dir{C8(2)/4D} | 12 | C2/m | C2h^3 | 16 | 4 |
/// | dir{C12(2)/4D} | 12 | C2/m | C2h^3 | 16 | 4 |
/// | dir{C10(2)/4D} | 12 | C2/m | C2h^3 | 16 | 4 |
/// | dir{4D1(4)/4D} | 2 | P-1 | Ci^1 | 32 | 4 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{P3(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{P5(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{P4(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{C1(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C3(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C2(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P12(1)/4D} | 70 | Fddd | D2h^24 | 8 | 4 |
/// | dir{P11(1)/4D} | 69 | Fmmm | D2h^23 | 8 | 4 |
/// | dir{C9(2)/4D} | 15 | C2/c | C2h^6 | 16 | 4 |
/// | dir{C13(2)/4D} | 15 | C2/c | C2h^6 | 16 | 4 |
/// | dir{C11(2)/4D} | 15 | C2/c | C2h^6 | 16 | 4 |
/// | dir{C8(2)/4D} | 12 | C2/m | C2h^3 | 16 | 4 |
/// | dir{C12(2)/4D} | 12 | C2/m | C2h^3 | 16 | 4 |
/// | dir{C10(2)/4D} | 12 | C2/m | C2h^3 | 16 | 4 |
/// | dir{4D1(4)/4D} | 2 | P-1 | Ci^1 | 32 | 4 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1+ | $X_{1}^+$ | $k_{15}t1$ | 1 | A2a | no |
/// | T1- | $X_{1}^-$ | $k_{15}t2$ | 1 | A2a | no |
/// | T2+ | $X_{3}^+$ | $k_{15}t7$ | 1 | A2a | no |
/// | T2- | $X_{3}^-$ | $k_{15}t8$ | 1 | A2a | no |
/// | T3+ | $X_{4}^+$ | $k_{15}t3$ | 1 | A2a | no |
/// | T3- | $X_{4}^-$ | $k_{15}t4$ | 1 | A2a | no |
/// | T4+ | $X_{2}^+$ | $k_{15}t5$ | 1 | A2a | no |
/// | T4- | $X_{2}^-$ | $k_{15}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for T1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 68 | Ccca | D2h^22 | 2 | 1 |
///
/// ### Isotropy subgroups for T2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for T2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 2 | 1 |
///
/// ### Isotropy subgroups for T3+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 66 | Cccm | D2h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for T3- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for T4+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for T4- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{17}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{3}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | Z2- | $Z_{3}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | Z3+ | $Z_{4}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | Z3- | $Z_{4}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | Z4+ | $Z_{2}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | Z4- | $Z_{2}^-$ | $k_{17}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 68 | Ccca | D2h^22 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 66 | Cccm | D2h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3+ ($Z_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for Z3- ($Z_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for Z4- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Y_{1}^+$ | $k_{16}t1$ | 1 | A2a | no |
/// | Y1- | $Y_{1}^-$ | $k_{16}t2$ | 1 | A2a | no |
/// | Y2+ | $Y_{3}^+$ | $k_{16}t7$ | 1 | A2a | no |
/// | Y2- | $Y_{3}^-$ | $k_{16}t8$ | 1 | A2a | no |
/// | Y3+ | $Y_{4}^+$ | $k_{16}t3$ | 1 | A2a | no |
/// | Y3- | $Y_{4}^-$ | $k_{16}t4$ | 1 | A2a | no |
/// | Y4+ | $Y_{2}^+$ | $k_{16}t5$ | 1 | A2a | no |
/// | Y4- | $Y_{2}^-$ | $k_{16}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 65 | Cmmm | D2h^19 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 68 | Ccca | D2h^22 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Y_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Y_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3+ ($Y_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 64 | Cmca | D2h^18 | 2 | 1 |
///
/// ### Isotropy subgroups for Y3- ($Y_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 63 | Cmcm | D2h^17 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 66 | Cccm | D2h^20 | 2 | 1 |
///
/// ### Isotropy subgroups for Y4- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 67 | Cmma | D2h^21 | 2 | 1 |
///
pub struct Sg69;

/// # 70 Fddd (D₂ₕ²⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{14}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{14}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{14}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{14}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{14}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{14}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{14}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 70 | Fddd | D2h^24 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 22 | F222 | D2^7 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 43 | Fdd2 | C2v^19 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 43 | Fdd2 | C2v^19 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 43 | Fdd2 | C2v^19 | 2 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^+$ | $k_{10}t1$ | 4 | D64c | no |
/// | L1- | $L_{1}^-$ | $k_{10}t2$ | 4 | D64c | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{P3(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P4(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P5(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{C1(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C2(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C3(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P11(1)/4D} | 22 | F222 | D2^7 | 16 | 4 |
/// | dir{C8(2)/4D} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{C10(2)/4D} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{C12(2)/4D} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{S1(3)/4D} | 2 | P-1 | Ci^1 | 32 | 4 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{P3(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P4(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P5(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{C1(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C2(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C3(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{P11(1)/4D} | 22 | F222 | D2^7 | 16 | 4 |
/// | dir{C8(2)/4D} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{C10(2)/4D} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{C12(2)/4D} | 5 | C2 | C2^3 | 32 | 4 |
/// | dir{S1(3)/4D} | 2 | P-1 | Ci^1 | 32 | 4 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 64 | 4 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $X_{1}$ | $k_{15}t1$ | 2 | B8a | yes |
/// | T2 | $X_{2}$ | $k_{15}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for T1 ($X_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for T2 ($X_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{17}t1$ | 2 | B8a | yes |
/// | Z2 | $Z_{2}$ | $k_{17}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{16}t1$ | 2 | B8a | yes |
/// | Y2 | $Y_{2}$ | $k_{16}t2$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 21 | C222 | D2^6 | 4 | 1 |
/// | dir{(a,0)} | 13 | P2/c | C2h^4 | 4 | 1 |
/// | dir{(a,b)} | 3 | P2 | C2^1 | 8 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 20 | C222_1 | D2^5 | 4 | 1 |
/// | dir{(a,0)} | 14 | P2_1/c | C2h^5 | 4 | 1 |
/// | dir{(a,b)} | 4 | P2_1 | C2^2 | 8 | 1 |
///
pub struct Sg70;

/// # 71 Immm (D₂ₕ²⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{17}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 71 | Immm | D2h^25 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 23 | I222 | D2^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 44 | Imm2 | C2v^20 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 44 | Imm2 | C2v^20 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 44 | Imm2 | C2v^20 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{18}t1$ | 1 | A2a | no |
/// | X1- | $X_{1}^-$ | $k_{18}t2$ | 1 | A2a | no |
/// | X2+ | $X_{3}^+$ | $k_{18}t7$ | 1 | A2a | no |
/// | X2- | $X_{3}^-$ | $k_{18}t8$ | 1 | A2a | no |
/// | X3+ | $X_{4}^+$ | $k_{18}t3$ | 1 | A2a | no |
/// | X3- | $X_{4}^-$ | $k_{18}t4$ | 1 | A2a | no |
/// | X4+ | $X_{2}^+$ | $k_{18}t5$ | 1 | A2a | no |
/// | X4- | $X_{2}^-$ | $k_{18}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 47 | Pmmm | D2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 48 | Pnnn | D2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for X2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 58 | Pnnm | D2h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for X2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 59 | Pmmn | D2h^13 | 2 | 1 |
///
/// ### Isotropy subgroups for X3+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 58 | Pnnm | D2h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for X3- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 59 | Pmmn | D2h^13 | 2 | 1 |
///
/// ### Isotropy subgroups for X4+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 58 | Pnnm | D2h^12 | 2 | 1 |
///
/// ### Isotropy subgroups for X4- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 59 | Pmmn | D2h^13 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{11}t1$ | 2 | B8a | no |
/// | R1- | $R_{1}^-$ | $k_{11}t2$ | 2 | B8a | no |
/// | R2+ | $R_{2}^+$ | $k_{11}t3$ | 2 | B8a | no |
/// | R2- | $R_{2}^-$ | $k_{11}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 65 | Cmmm | D2h^19 | 4 | 2 |
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,0)} | 67 | Cmma | D2h^21 | 4 | 2 |
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for R2+ ($R_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,0)} | 67 | Cmma | D2h^21 | 4 | 2 |
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for R2- ($R_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 65 | Cmmm | D2h^19 | 4 | 2 |
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1+ | $T_{1}^+$ | $k_{12}t1$ | 2 | B8a | no |
/// | T1- | $T_{1}^-$ | $k_{12}t2$ | 2 | B8a | no |
/// | T2+ | $T_{2}^+$ | $k_{12}t3$ | 2 | B8a | no |
/// | T2- | $T_{2}^-$ | $k_{12}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for T1+ ($T_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 65 | Cmmm | D2h^19 | 4 | 2 |
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ### Isotropy subgroups for T1- ($T_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,0)} | 67 | Cmma | D2h^21 | 4 | 2 |
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for T2+ ($T_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,0)} | 67 | Cmma | D2h^21 | 4 | 2 |
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for T2- ($T_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 65 | Cmmm | D2h^19 | 4 | 2 |
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1 | $W_{1}$ | $k_{16}t1$ | 2 | B8a | no |
/// | W2 | $W_{3}$ | $k_{16}t4$ | 2 | B8a | no |
/// | W3 | $W_{4}$ | $k_{16}t2$ | 2 | B8a | no |
/// | W4 | $W_{2}$ | $k_{16}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for W1 ($W_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 70 | Fddd | D2h^24 | 4 | 2 |
/// | dir{(a,a)} | 69 | Fmmm | D2h^23 | 4 | 2 |
/// | dir{(a,b)} | 22 | F222 | D2^7 | 8 | 2 |
///
/// ### Isotropy subgroups for W2 ($W_{3}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 70 | Fddd | D2h^24 | 4 | 2 |
/// | dir{(a,a)} | 69 | Fmmm | D2h^23 | 4 | 2 |
/// | dir{(a,b)} | 22 | F222 | D2^7 | 8 | 2 |
///
/// ### Isotropy subgroups for W3 ($W_{4}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 70 | Fddd | D2h^24 | 4 | 2 |
/// | dir{(a,a)} | 69 | Fmmm | D2h^23 | 4 | 2 |
/// | dir{(a,b)} | 22 | F222 | D2^7 | 8 | 2 |
///
/// ### Isotropy subgroups for W4 ($W_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,0)} | 70 | Fddd | D2h^24 | 4 | 2 |
/// | dir{(a,a)} | 69 | Fmmm | D2h^23 | 4 | 2 |
/// | dir{(a,b)} | 22 | F222 | D2^7 | 8 | 2 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1+ | $S_{1}^+$ | $k_{10}t1$ | 2 | B8a | no |
/// | S1- | $S_{1}^-$ | $k_{10}t2$ | 2 | B8a | no |
/// | S2+ | $S_{2}^+$ | $k_{10}t3$ | 2 | B8a | no |
/// | S2- | $S_{2}^-$ | $k_{10}t4$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1+ ($S_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 65 | Cmmm | D2h^19 | 4 | 2 |
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ### Isotropy subgroups for S1- ($S_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,0)} | 67 | Cmma | D2h^21 | 4 | 2 |
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for S2+ ($S_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,0)} | 67 | Cmma | D2h^21 | 4 | 2 |
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for S2- ($S_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 65 | Cmmm | D2h^19 | 4 | 2 |
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 8 | 2 |
///
pub struct Sg71;

/// # 72 Ibam (D₂ₕ²⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{17}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 72 | Ibam | D2h^26 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 23 | I222 | D2^8 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 45 | Iba2 | C2v^21 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 46 | Ima2 | C2v^22 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 46 | Ima2 | C2v^22 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{1}^+$ | $k_{18}t2$ | 1 | A2a | no |
/// | X1- | $X_{1}^-$ | $k_{18}t1$ | 1 | A2a | no |
/// | X2+ | $X_{3}^+$ | $k_{18}t8$ | 1 | A2a | no |
/// | X2- | $X_{3}^-$ | $k_{18}t7$ | 1 | A2a | no |
/// | X3+ | $X_{4}^+$ | $k_{18}t4$ | 1 | A2a | no |
/// | X3- | $X_{4}^-$ | $k_{18}t3$ | 1 | A2a | no |
/// | X4+ | $X_{2}^+$ | $k_{18}t6$ | 1 | A2a | no |
/// | X4- | $X_{2}^-$ | $k_{18}t5$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 49 | Pccm | D2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for X1- ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 50 | Pban | D2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for X2+ ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 55 | Pbam | D2h^9 | 2 | 1 |
///
/// ### Isotropy subgroups for X2- ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 56 | Pccn | D2h^10 | 2 | 1 |
///
/// ### Isotropy subgroups for X3+ ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 60 | Pbcn | D2h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for X3- ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 57 | Pbcm | D2h^11 | 2 | 1 |
///
/// ### Isotropy subgroups for X4+ ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 60 | Pbcn | D2h^14 | 2 | 1 |
///
/// ### Isotropy subgroups for X4- ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 57 | Pbcm | D2h^11 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{11}t1$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P3(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P5(1)/4D} | 36 | Cmc2_1 | C2v^12 | 8 | 2 |
/// | dir{P4(1)/4D} | 21 | C222 | D2^6 | 8 | 2 |
/// | dir{P12(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P11(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{C13(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C12(2)/4D} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{C10(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C3(2)/4D} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{C2(2)/4D} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{C9(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C8(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1+ | $T_{1}^+$ | $k_{12}t2$ | 2 | B8a | no |
/// | T1- | $T_{1}^-$ | $k_{12}t1$ | 2 | B8a | no |
/// | T2+ | $T_{2}^+$ | $k_{12}t4$ | 2 | B8a | no |
/// | T2- | $T_{2}^-$ | $k_{12}t3$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for T1+ ($T_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 66 | Cccm | D2h^20 | 4 | 2 |
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ### Isotropy subgroups for T1- ($T_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,0)} | 68 | Ccca | D2h^22 | 4 | 2 |
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for T2+ ($T_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,0)} | 68 | Ccca | D2h^22 | 4 | 2 |
/// | dir{(a,b)} | 13 | P2/c | C2h^4 | 8 | 2 |
///
/// ### Isotropy subgroups for T2- ($T_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 66 | Cccm | D2h^20 | 4 | 2 |
/// | dir{(a,b)} | 10 | P2/m | C2h^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1W2 | $W_{1}W3$ | $k_{16}t1t4$ | 4 | D16c | no |
/// | W3W4 | $W_{2}W4$ | $k_{16}t2t3$ | 4 | D16c | no |
///
/// ### Isotropy subgroups for W1W2 ($W_{1}W3$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 22 | F222 | D2^7 | 8 | 2 |
/// | dir{C12(2)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{C11(2)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{4D1(4)/4D} | 5 | C2 | C2^3 | 16 | 2 |
///
/// ### Isotropy subgroups for W3W4 ($W_{2}W4$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{C1(2)/4D} | 22 | F222 | D2^7 | 8 | 2 |
/// | dir{C12(2)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{C11(2)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{4D1(4)/4D} | 5 | C2 | C2^3 | 16 | 2 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{10}t1$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P3(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P5(1)/4D} | 36 | Cmc2_1 | C2v^12 | 8 | 2 |
/// | dir{P4(1)/4D} | 21 | C222 | D2^6 | 8 | 2 |
/// | dir{P12(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P11(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{C13(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C12(2)/4D} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{C10(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C3(2)/4D} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{C2(2)/4D} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{C9(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C8(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 2 |
///
pub struct Sg72;

/// # 73 Ibca (D₂ₕ²⁷)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{17}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 73 | Ibca | D2h^27 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 24 | I2_12_12_1 | D2^9 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 45 | Iba2 | C2v^21 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 45 | Iba2 | C2v^21 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 45 | Iba2 | C2v^21 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{3}^-$ | $k_{18}t5$ | 1 | A2a | no |
/// | X1- | $X_{3}^+$ | $k_{18}t6$ | 1 | A2a | no |
/// | X2+ | $X_{1}^-$ | $k_{18}t3$ | 1 | A2a | no |
/// | X2- | $X_{1}^+$ | $k_{18}t4$ | 1 | A2a | no |
/// | X3+ | $X_{2}^-$ | $k_{18}t7$ | 1 | A2a | no |
/// | X3- | $X_{2}^+$ | $k_{18}t8$ | 1 | A2a | no |
/// | X4+ | $X_{4}^-$ | $k_{18}t1$ | 1 | A2a | no |
/// | X4- | $X_{4}^+$ | $k_{18}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for X1- ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for X2+ ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for X2- ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for X3+ ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for X3- ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 54 | Pcca | D2h^8 | 2 | 1 |
///
/// ### Isotropy subgroups for X4+ ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 61 | Pbca | D2h^15 | 2 | 1 |
///
/// ### Isotropy subgroups for X4- ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 61 | Pbca | D2h^15 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $R_{1}$ | $k_{11}t1$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for R1 ($R_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P3(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P5(1)/4D} | 37 | Ccc2 | C2v^13 | 8 | 2 |
/// | dir{P4(1)/4D} | 20 | C222_1 | D2^5 | 8 | 2 |
/// | dir{P11(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P12(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{C12(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C13(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C10(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C2(2)/4D} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{C3(2)/4D} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{C8(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C9(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{12}t1$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P3(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P5(1)/4D} | 37 | Ccc2 | C2v^13 | 8 | 2 |
/// | dir{P4(1)/4D} | 20 | C222_1 | D2^5 | 8 | 2 |
/// | dir{P12(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P11(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{C13(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C12(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C10(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C2(2)/4D} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{C3(2)/4D} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{C8(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C9(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1W1 | $W_{1}W1$ | $k_{16}t1t1$ | 8 | F32b | no |
///
/// ### Isotropy subgroups for W1W1 ($W_{1}W1$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{4D3(4)/8D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{4D24(4)/8D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{4D10(4)/8D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{4D18(4)/8D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{4D25(4)/8D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{8D1(8)/8D} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1 | $S_{1}$ | $k_{10}t1$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for S1 ($S_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P3(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P5(1)/4D} | 37 | Ccc2 | C2v^13 | 8 | 2 |
/// | dir{P4(1)/4D} | 20 | C222_1 | D2^5 | 8 | 2 |
/// | dir{P11(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P12(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{C12(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C13(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C10(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C2(2)/4D} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{C3(2)/4D} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{C9(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C8(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 2 |
///
pub struct Sg73;

/// # 74 Imma (D₂ₕ²⁸)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{17}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{17}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{3}^+$ | $k_{17}t7$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{3}^-$ | $k_{17}t8$ | 1 | A2a | no |
/// | GM3+ | $\Gamma_{4}^+$ | $k_{17}t3$ | 1 | A2a | no |
/// | GM3- | $\Gamma_{4}^-$ | $k_{17}t4$ | 1 | A2a | no |
/// | GM4+ | $\Gamma_{2}^+$ | $k_{17}t5$ | 1 | A2a | no |
/// | GM4- | $\Gamma_{2}^-$ | $k_{17}t6$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 74 | Imma | D2h^28 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 24 | I2_12_12_1 | D2^9 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 44 | Imm2 | C2v^20 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3+ ($\Gamma_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM3- ($\Gamma_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 46 | Ima2 | C2v^22 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM4- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 46 | Ima2 | C2v^22 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $X_{3}^-$ | $k_{18}t6$ | 1 | A2a | no |
/// | X1- | $X_{3}^+$ | $k_{18}t5$ | 1 | A2a | no |
/// | X2+ | $X_{1}^-$ | $k_{18}t4$ | 1 | A2a | no |
/// | X2- | $X_{1}^+$ | $k_{18}t3$ | 1 | A2a | no |
/// | X3+ | $X_{2}^-$ | $k_{18}t8$ | 1 | A2a | no |
/// | X3- | $X_{2}^+$ | $k_{18}t7$ | 1 | A2a | no |
/// | X4+ | $X_{4}^-$ | $k_{18}t2$ | 1 | A2a | no |
/// | X4- | $X_{4}^+$ | $k_{18}t1$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1+ ($X_{3}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for X1- ($X_{3}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for X2+ ($X_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 52 | Pnna | D2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for X2- ($X_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 51 | Pmma | D2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for X3+ ($X_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for X3- ($X_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 53 | Pmna | D2h^7 | 2 | 1 |
///
/// ### Isotropy subgroups for X4+ ($X_{4}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ### Isotropy subgroups for X4- ($X_{4}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 62 | Pnma | D2h^16 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $R_{1}^+$ | $k_{11}t3$ | 2 | B8a | no |
/// | R1- | $R_{1}^-$ | $k_{11}t4$ | 2 | B8a | no |
/// | R2+ | $R_{2}^+$ | $k_{11}t1$ | 2 | B8a | no |
/// | R2- | $R_{2}^-$ | $k_{11}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for R1+ ($R_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 63 | Cmcm | D2h^17 | 4 | 2 |
/// | dir{(a,b)} | 11 | P2_1/m | C2h^2 | 8 | 2 |
///
/// ### Isotropy subgroups for R1- ($R_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,0)} | 64 | Cmca | D2h^18 | 4 | 2 |
/// | dir{(a,b)} | 14 | P2_1/c | C2h^5 | 8 | 2 |
///
/// ### Isotropy subgroups for R2+ ($R_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,0)} | 64 | Cmca | D2h^18 | 4 | 2 |
/// | dir{(a,b)} | 14 | P2_1/c | C2h^5 | 8 | 2 |
///
/// ### Isotropy subgroups for R2- ($R_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 63 | Cmcm | D2h^17 | 4 | 2 |
/// | dir{(a,b)} | 11 | P2_1/m | C2h^2 | 8 | 2 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $T_{1}$ | $k_{12}t1$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for T1 ($T_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P1(1)/4D} | 5 | C2 | C2^3 | 8 | 1 |
/// | dir{P3(1)/4D} | 2 | P-1 | Ci^1 | 8 | 1 |
/// | dir{C1(2)/4D} | 1 | P1 | C1^1 | 16 | 1 |
/// | dir{P5(1)/4D} | 35 | Cmm2 | C2v^11 | 8 | 2 |
/// | dir{P4(1)/4D} | 20 | C222_1 | D2^5 | 8 | 2 |
/// | dir{P12(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{P11(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{C13(2)/4D} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{C12(2)/4D} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{C11(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C10(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C2(2)/4D} | 4 | P2_1 | C2^2 | 16 | 2 |
/// | dir{C3(2)/4D} | 3 | P2 | C2^1 | 16 | 2 |
/// | dir{C8(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C9(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{W}$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | W1 | $W_{1}$ | $k_{16}t1$ | 4 | D32c | no |
///
/// ### Isotropy subgroups for W1 ($W_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{P3(1)/4D} | 43 | Fdd2 | C2v^19 | 8 | 2 |
/// | dir{P11(1)/4D} | 42 | Fmm2 | C2v^18 | 8 | 2 |
/// | dir{P12(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P1(1)/4D} | 15 | C2/c | C2h^6 | 8 | 2 |
/// | dir{P5(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{P4(1)/4D} | 12 | C2/m | C2h^3 | 8 | 2 |
/// | dir{C9(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C1(2)/4D} | 9 | Cc | Cs^4 | 16 | 2 |
/// | dir{C12(2)/4D} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{C10(2)/4D} | 8 | Cm | Cs^3 | 16 | 2 |
/// | dir{C8(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C13(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C2(2)/4D} | 5 | C2 | C2^3 | 16 | 2 |
/// | dir{C11(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{C3(2)/4D} | 2 | P-1 | Ci^1 | 16 | 2 |
/// | dir{4D1(4)/4D} | 1 | P1 | C1^1 | 32 | 2 |
///
/// ## Irreps at $\mathrm{S}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | S1+ | $S_{2}^+$ | $k_{10}t4$ | 2 | B8a | no |
/// | S1- | $S_{2}^-$ | $k_{10}t3$ | 2 | B8a | no |
/// | S2+ | $S_{1}^+$ | $k_{10}t2$ | 2 | B8a | no |
/// | S2- | $S_{1}^-$ | $k_{10}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for S1+ ($S_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 63 | Cmcm | D2h^17 | 4 | 2 |
/// | dir{(a,b)} | 11 | P2_1/m | C2h^2 | 8 | 2 |
///
/// ### Isotropy subgroups for S1- ($S_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,0)} | 64 | Cmca | D2h^18 | 4 | 2 |
/// | dir{(a,b)} | 14 | P2_1/c | C2h^5 | 8 | 2 |
///
/// ### Isotropy subgroups for S2+ ($S_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 15 | C2/c | C2h^6 | 4 | 1 |
/// | dir{(a,0)} | 64 | Cmca | D2h^18 | 4 | 2 |
/// | dir{(a,b)} | 14 | P2_1/c | C2h^5 | 8 | 2 |
///
/// ### Isotropy subgroups for S2- ($S_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 12 | C2/m | C2h^3 | 4 | 1 |
/// | dir{(a,0)} | 63 | Cmcm | D2h^17 | 4 | 2 |
/// | dir{(a,b)} | 11 | P2_1/m | C2h^2 | 8 | 2 |
///
pub struct Sg74;
