//! # Monoclinic space groups (#3–#15)
//!
//! Unique axis **b**, cell choice 1 (International Tables).  The Brillouin
//! zone is a monoclinic prism.  High-symmetry k-points include **Γ**, plus
//! points on the unique-axis zone boundary.
//!
//! ## Common k-point labels (monoclinic, b-unique)
//!
//! | Label | Coords (fractional) | Little group |
//! |-------|---------------------|--------------|
//! | Γ | (0, 0, 0) | 2/m (C₂ₕ) |
//! | Y | (0, ½, 0) | 2/m (C₂ₕ) |
//! | A | (½, ½, 0) | 2/m (C₂ₕ) |
//! | M | (½, ½, ½) | 2/m (C₂ₕ) |
//!
//! ---


/// # 3 P2 (C₂¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{7}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{14}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{14}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{11}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{2}$ | $k_{11}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ## Irreps at $D$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | D1 | $C_{1}$ | $k_{8}t1$ | 1 | A2a | no |
/// | D2 | $C_{2}$ | $k_{8}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for D1 ($C_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for D2 ($C_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $B_{1}$ | $k_{13}t1$ | 1 | A2a | no |
/// | Y2 | $B_{2}$ | $k_{13}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($B_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($B_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ## Irreps at $B$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | B1 | $Y_{1}$ | $k_{12}t1$ | 1 | A2a | no |
/// | B2 | $Y_{2}$ | $k_{12}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for B1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for B2 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ## Irreps at $C$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | C1 | $D_{1}$ | $k_{9}t1$ | 1 | A2a | no |
/// | C2 | $D_{2}$ | $k_{9}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for C1 ($D_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for C2 ($D_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ## Irreps at $E$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | E1 | $E_{1}$ | $k_{10}t1$ | 1 | A2a | no |
/// | E2 | $E_{2}$ | $k_{10}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for E1 ($E_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for E2 ($E_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
pub struct Sg3;

/// # 4 P2_1 (C₂²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{7}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{14}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{14}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1Z2 | $Z_{1}Z2$ | $k_{11}t1t2$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for Z1Z2 ($Z_{1}Z2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $D$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | D1D2 | $C_{1}C2$ | $k_{8}t1t2$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for D1D2 ($C_{1}C2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $B_{1}$ | $k_{13}t1$ | 1 | A2a | no |
/// | Y2 | $B_{2}$ | $k_{13}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($B_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($B_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ## Irreps at $B$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | B1 | $Y_{1}$ | $k_{12}t1$ | 1 | A2a | no |
/// | B2 | $Y_{2}$ | $k_{12}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for B1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for B2 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ## Irreps at $C$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | C1C2 | $D_{1}D2$ | $k_{9}t1t2$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for C1C2 ($D_{1}D2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $E$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | E1E2 | $E_{1}E2$ | $k_{10}t1t2$ | 2 | B4a | yes |
///
/// ### Isotropy subgroups for E1E2 ($E_{1}E2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 1 | P1 | C1^1 | 4 | 1 |
///
pub struct Sg4;

/// # 5 C2 (C₂³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{6}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{6}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{9}t1$ | 1 | A2a | no |
/// | M2 | $M_{2}$ | $k_{9}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{7}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{7}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{5}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 1 | P1 | C1^1 | 4 | 1 |
/// | dir{(a,0)} | 5 | C2 | C2^3 | 4 | 2 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{V}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | V1 | $V_{1}$ | $k_{4}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for V1 ($V_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 1 | P1 | C1^1 | 4 | 1 |
/// | dir{(a,0)} | 5 | C2 | C2^3 | 4 | 2 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Z_{1}$ | $k_{8}t1$ | 1 | A2a | no |
/// | Y2 | $Z_{2}$ | $k_{8}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
pub struct Sg5;

/// # 6 Pm (Cₛ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{7}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{14}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{14}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{11}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{2}$ | $k_{11}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ## Irreps at $D$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | D1 | $C_{1}$ | $k_{8}t1$ | 1 | A2a | no |
/// | D2 | $C_{2}$ | $k_{8}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for D1 ($C_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ### Isotropy subgroups for D2 ($C_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $B_{1}$ | $k_{13}t1$ | 1 | A2a | no |
/// | Y2 | $B_{2}$ | $k_{13}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($B_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($B_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $B$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | B1 | $Y_{1}$ | $k_{12}t1$ | 1 | A2a | no |
/// | B2 | $Y_{2}$ | $k_{12}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for B1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ### Isotropy subgroups for B2 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $C$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | C1 | $D_{1}$ | $k_{9}t1$ | 1 | A2a | no |
/// | C2 | $D_{2}$ | $k_{9}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for C1 ($D_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ### Isotropy subgroups for C2 ($D_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ## Irreps at $E$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | E1 | $E_{1}$ | $k_{10}t1$ | 1 | A2a | no |
/// | E2 | $E_{2}$ | $k_{10}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for E1 ($E_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ### Isotropy subgroups for E2 ($E_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
pub struct Sg6;

/// # 7 Pc (Cₛ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{7}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{14}t1t2$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{11}t1$ | 1 | A2a | no |
/// | Z2 | $Z_{2}$ | $k_{11}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $D$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | D1D2 | $D_{1}D2$ | $k_{8}t1t2$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for D1D2 ($D_{1}D2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Y_{1}$ | $k_{13}t1$ | 1 | A2a | no |
/// | Y2 | $Y_{2}$ | $k_{13}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Y_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Y_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $B$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | B1B2 | $B_{1}B2$ | $k_{12}t1t2$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for B1B2 ($B_{1}B2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $C$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | C1 | $C_{1}$ | $k_{9}t1$ | 1 | A2a | no |
/// | C2 | $C_{2}$ | $k_{9}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for C1 ($C_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ### Isotropy subgroups for C2 ($C_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ## Irreps at $E$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | E1E2 | $E_{1}E2$ | $k_{10}t1t2$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for E1E2 ($E_{1}E2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 1 | P1 | C1^1 | 4 | 1 |
///
pub struct Sg7;

/// # 8 Cm (Cₛ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{6}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{6}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{9}t1$ | 1 | A2a | no |
/// | M2 | $M_{2}$ | $k_{9}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M2 ($M_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{7}t1$ | 1 | A2a | no |
/// | A2 | $A_{2}$ | $k_{7}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A2 ($A_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{5}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 1 | P1 | C1^1 | 4 | 1 |
/// | dir{(a,0)} | 8 | Cm | Cs^3 | 4 | 2 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{V}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | V1 | $V_{1}$ | $k_{4}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for V1 ($V_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 1 | P1 | C1^1 | 4 | 1 |
/// | dir{(a,0)} | 8 | Cm | Cs^3 | 4 | 2 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Z_{1}$ | $k_{8}t1$ | 1 | A2a | no |
/// | Y2 | $Z_{2}$ | $k_{8}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
pub struct Sg8;

/// # 9 Cc (Cₛ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{6}t1$ | 1 | A1a | no |
/// | GM2 | $\Gamma_{2}$ | $k_{6}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM2 ($\Gamma_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1M2 | $M_{1}M2$ | $k_{9}t1t2$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for M1M2 ($M_{1}M2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1A2 | $A_{1}A2$ | $k_{7}t1t2$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for A1A2 ($A_{1}A2$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 1 | P1 | C1^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1 | $L_{1}$ | $k_{5}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for L1 ($L_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 1 | P1 | C1^1 | 4 | 1 |
/// | dir{(a,0)} | 9 | Cc | Cs^4 | 4 | 2 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{V}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | V1 | $V_{1}$ | $k_{4}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for V1 ($V_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 1 | P1 | C1^1 | 4 | 1 |
/// | dir{(a,0)} | 9 | Cc | Cs^4 | 4 | 2 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $Z_{1}$ | $k_{8}t1$ | 1 | A2a | no |
/// | Y2 | $Z_{2}$ | $k_{8}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2 ($Z_{2}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
pub struct Sg9;

/// # 10 P2/m (C₂ₕ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{7}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{7}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{7}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1+ | $A_{1}^+$ | $k_{14}t1$ | 1 | A2a | no |
/// | A1- | $A_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | A2+ | $A_{2}^+$ | $k_{14}t3$ | 1 | A2a | no |
/// | A2- | $A_{2}^-$ | $k_{14}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for A1+ ($A_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for A1- ($A_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A2+ ($A_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for A2- ($A_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{11}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{11}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{2}^+$ | $k_{11}t3$ | 1 | A2a | no |
/// | Z2- | $Z_{2}^-$ | $k_{11}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
/// ## Irreps at $D$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | D1+ | $C_{1}^+$ | $k_{8}t1$ | 1 | A2a | no |
/// | D1- | $C_{1}^-$ | $k_{8}t2$ | 1 | A2a | no |
/// | D2+ | $C_{2}^+$ | $k_{8}t3$ | 1 | A2a | no |
/// | D2- | $C_{2}^-$ | $k_{8}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for D1+ ($C_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for D1- ($C_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for D2+ ($C_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for D2- ($C_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $B_{1}^+$ | $k_{13}t1$ | 1 | A2a | no |
/// | Y1- | $B_{1}^-$ | $k_{13}t2$ | 1 | A2a | no |
/// | Y2+ | $B_{2}^+$ | $k_{13}t3$ | 1 | A2a | no |
/// | Y2- | $B_{2}^-$ | $k_{13}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($B_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($B_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($B_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($B_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ## Irreps at $B$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | B1+ | $Y_{1}^+$ | $k_{12}t1$ | 1 | A2a | no |
/// | B1- | $Y_{1}^-$ | $k_{12}t2$ | 1 | A2a | no |
/// | B2+ | $Y_{2}^+$ | $k_{12}t3$ | 1 | A2a | no |
/// | B2- | $Y_{2}^-$ | $k_{12}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for B1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for B1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for B2+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for B2- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ## Irreps at $C$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | C1+ | $D_{1}^+$ | $k_{9}t1$ | 1 | A2a | no |
/// | C1- | $D_{1}^-$ | $k_{9}t2$ | 1 | A2a | no |
/// | C2+ | $D_{2}^+$ | $k_{9}t3$ | 1 | A2a | no |
/// | C2- | $D_{2}^-$ | $k_{9}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for C1+ ($D_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for C1- ($D_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for C2+ ($D_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for C2- ($D_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ## Irreps at $E$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | E1+ | $E_{1}^+$ | $k_{10}t1$ | 1 | A2a | no |
/// | E1- | $E_{1}^-$ | $k_{10}t2$ | 1 | A2a | no |
/// | E2+ | $E_{2}^+$ | $k_{10}t3$ | 1 | A2a | no |
/// | E2- | $E_{2}^-$ | $k_{10}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for E1+ ($E_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for E1- ($E_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for E2+ ($E_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for E2- ($E_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
pub struct Sg10;

/// # 11 P2_1/m (C₂ₕ²)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{7}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{7}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{7}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 6 | Pm | Cs^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1+ | $A_{1}^+$ | $k_{14}t1$ | 1 | A2a | no |
/// | A1- | $A_{1}^-$ | $k_{14}t2$ | 1 | A2a | no |
/// | A2+ | $A_{2}^+$ | $k_{14}t3$ | 1 | A2a | no |
/// | A2- | $A_{2}^-$ | $k_{14}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for A1+ ($A_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for A1- ($A_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for A2+ ($A_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for A2- ($A_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{11}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 6 | Pm | Cs^1 | 4 | 1 |
/// | dir{(a,0)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $D$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | D1 | $C_{1}$ | $k_{8}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for D1 ($C_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,0)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $B_{1}^+$ | $k_{13}t1$ | 1 | A2a | no |
/// | Y1- | $B_{1}^-$ | $k_{13}t2$ | 1 | A2a | no |
/// | Y2+ | $B_{2}^+$ | $k_{13}t3$ | 1 | A2a | no |
/// | Y2- | $B_{2}^-$ | $k_{13}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($B_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($B_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($B_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($B_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
/// ## Irreps at $B$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | B1+ | $Y_{1}^+$ | $k_{12}t1$ | 1 | A2a | no |
/// | B1- | $Y_{1}^-$ | $k_{12}t2$ | 1 | A2a | no |
/// | B2+ | $Y_{2}^+$ | $k_{12}t3$ | 1 | A2a | no |
/// | B2- | $Y_{2}^-$ | $k_{12}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for B1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
/// ### Isotropy subgroups for B1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for B2+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for B2- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
/// ## Irreps at $C$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | C1 | $D_{1}$ | $k_{9}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for C1 ($D_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,0)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $E$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | E1 | $E_{1}$ | $k_{10}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for E1 ($E_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 8 | Cm | Cs^3 | 4 | 1 |
/// | dir{(a,0)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
pub struct Sg11;

/// # 12 C2/m (C₂ₕ³)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{6}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{6}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{6}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{6}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 8 | Cm | Cs^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1+ | $M_{1}^+$ | $k_{9}t1$ | 1 | A2a | no |
/// | M1- | $M_{1}^-$ | $k_{9}t2$ | 1 | A2a | no |
/// | M2+ | $M_{2}^+$ | $k_{9}t3$ | 1 | A2a | no |
/// | M2- | $M_{2}^-$ | $k_{9}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for M1+ ($M_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for M1- ($M_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M2+ ($M_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for M2- ($M_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1+ | $A_{1}^+$ | $k_{7}t1$ | 1 | A2a | no |
/// | A1- | $A_{1}^-$ | $k_{7}t2$ | 1 | A2a | no |
/// | A2+ | $A_{2}^+$ | $k_{7}t3$ | 1 | A2a | no |
/// | A2- | $A_{2}^-$ | $k_{7}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for A1+ ($A_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ### Isotropy subgroups for A1- ($A_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for A2+ ($A_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for A2- ($A_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 12 | C2/m | C2h^3 | 2 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^+$ | $k_{5}t1$ | 2 | B8a | no |
/// | L1- | $L_{1}^-$ | $k_{5}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 2 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 2 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 2 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{V}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | V1+ | $V_{1}^+$ | $k_{4}t1$ | 2 | B8a | no |
/// | V1- | $V_{1}^-$ | $k_{4}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for V1+ ($V_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 2 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 2 |
///
/// ### Isotropy subgroups for V1- ($V_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,0)} | 12 | C2/m | C2h^3 | 4 | 2 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Z_{1}^+$ | $k_{8}t1$ | 1 | A2a | no |
/// | Y1- | $Z_{1}^-$ | $k_{8}t2$ | 1 | A2a | no |
/// | Y2+ | $Z_{2}^+$ | $k_{8}t3$ | 1 | A2a | no |
/// | Y2- | $Z_{2}^-$ | $k_{8}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 10 | P2/m | C2h^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 11 | P2_1/m | C2h^2 | 2 | 1 |
///
pub struct Sg12;

/// # 13 P2/c (C₂ₕ⁴)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{7}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{7}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{7}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 3 | P2 | C2^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{14}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,0)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $Z_{1}^+$ | $k_{11}t1$ | 1 | A2a | no |
/// | Z1- | $Z_{1}^-$ | $k_{11}t2$ | 1 | A2a | no |
/// | Z2+ | $Z_{2}^+$ | $k_{11}t3$ | 1 | A2a | no |
/// | Z2- | $Z_{2}^-$ | $k_{11}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Z2- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ## Irreps at $D$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | D1 | $D_{1}$ | $k_{8}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for D1 ($D_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Y_{1}^+$ | $k_{13}t1$ | 1 | A2a | no |
/// | Y1- | $Y_{1}^-$ | $k_{13}t2$ | 1 | A2a | no |
/// | Y2+ | $Y_{2}^+$ | $k_{13}t3$ | 1 | A2a | no |
/// | Y2- | $Y_{2}^-$ | $k_{13}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ## Irreps at $B$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | B1 | $B_{1}$ | $k_{12}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for B1 ($B_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 3 | P2 | C2^1 | 4 | 1 |
/// | dir{(a,0)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $C$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | C1+ | $C_{1}^+$ | $k_{9}t1$ | 1 | A2a | no |
/// | C1- | $C_{1}^-$ | $k_{9}t2$ | 1 | A2a | no |
/// | C2+ | $C_{2}^+$ | $k_{9}t3$ | 1 | A2a | no |
/// | C2- | $C_{2}^-$ | $k_{9}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for C1+ ($C_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for C1- ($C_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for C2+ ($C_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ### Isotropy subgroups for C2- ($C_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 2 | 1 |
///
/// ## Irreps at $E$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | E1 | $E_{1}$ | $k_{10}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for E1 ($E_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
pub struct Sg13;

/// # 14 P2_1/c (C₂ₕ⁵)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{7}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{7}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{7}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{7}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 4 | P2_1 | C2^2 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 7 | Pc | Cs^2 | 2 | 1 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{14}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,0)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $Z_{1}$ | $k_{11}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for Z1 ($Z_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 7 | Pc | Cs^2 | 4 | 1 |
/// | dir{(a,0)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $D$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | D1+D2+ | $D_{1}^-D2-$ | $k_{8}t2t4$ | 2 | B4a | no |
/// | D1-D2- | $D_{1}^+D2+$ | $k_{8}t1t3$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for D1+D2+ ($D_{1}^-D2-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 4 | 1 |
///
/// ### Isotropy subgroups for D1-D2- ($D_{1}^+D2+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 4 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Y_{1}^+$ | $k_{13}t1$ | 1 | A2a | no |
/// | Y1- | $Y_{1}^-$ | $k_{13}t2$ | 1 | A2a | no |
/// | Y2+ | $Y_{2}^+$ | $k_{13}t3$ | 1 | A2a | no |
/// | Y2- | $Y_{2}^-$ | $k_{13}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Y_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Y_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Y_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Y_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ## Irreps at $B$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | B1 | $B_{1}$ | $k_{12}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for B1 ($B_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 4 | P2_1 | C2^2 | 4 | 1 |
/// | dir{(a,0)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $C$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | C1 | $C_{1}$ | $k_{9}t1$ | 2 | B8a | yes |
///
/// ### Isotropy subgroups for C1 ($C_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 9 | Cc | Cs^4 | 4 | 1 |
/// | dir{(a,0)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $E$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | E1+E2+ | $E_{1}^-E2-$ | $k_{10}t2t4$ | 2 | B4a | no |
/// | E1-E2- | $E_{1}^+E2+$ | $k_{10}t1t3$ | 2 | B4a | no |
///
/// ### Isotropy subgroups for E1+E2+ ($E_{1}^-E2-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 4 | 1 |
///
/// ### Isotropy subgroups for E1-E2- ($E_{1}^+E2+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 4 | 1 |
///
pub struct Sg14;

/// # 15 C2/c (C₂ₕ⁶)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{6}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{6}t2$ | 1 | A2a | no |
/// | GM2+ | $\Gamma_{2}^+$ | $k_{6}t3$ | 1 | A2a | no |
/// | GM2- | $\Gamma_{2}^-$ | $k_{6}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 15 | C2/c | C2h^6 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 5 | C2 | C2^3 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2+ ($\Gamma_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ### Isotropy subgroups for GM2- ($\Gamma_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 9 | Cc | Cs^4 | 2 | 1 |
///
/// ## Irreps at $\mathrm{M}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | M1 | $M_{1}$ | $k_{9}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for M1 ($M_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{A}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | A1 | $A_{1}$ | $k_{7}t1$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for A1 ($A_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 5 | C2 | C2^3 | 4 | 1 |
/// | dir{(a,0)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,b)} | 1 | P1 | C1^1 | 8 | 1 |
///
/// ## Irreps at $\mathrm{L}\ (1/2,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | L1+ | $L_{1}^-$ | $k_{5}t1$ | 2 | B8a | no |
/// | L1- | $L_{1}^+$ | $k_{5}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for L1+ ($L_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 2 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 2 |
///
/// ### Isotropy subgroups for L1- ($L_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 2 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{V}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | V1+ | $V_{1}^+$ | $k_{4}t1$ | 2 | B8a | no |
/// | V1- | $V_{1}^-$ | $k_{4}t2$ | 2 | B8a | no |
///
/// ### Isotropy subgroups for V1+ ($V_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 2 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 2 |
///
/// ### Isotropy subgroups for V1- ($V_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a,a)} | 2 | P-1 | Ci^1 | 4 | 1 |
/// | dir{(a,0)} | 15 | C2/c | C2h^6 | 4 | 2 |
/// | dir{(a,b)} | 2 | P-1 | Ci^1 | 8 | 2 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $Z_{1}^+$ | $k_{8}t1$ | 1 | A2a | no |
/// | Y1- | $Z_{1}^-$ | $k_{8}t2$ | 1 | A2a | no |
/// | Y2+ | $Z_{2}^+$ | $k_{8}t3$ | 1 | A2a | no |
/// | Y2- | $Z_{2}^-$ | $k_{8}t4$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($Z_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($Z_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 13 | P2/c | C2h^4 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2+ ($Z_{2}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
/// ### Isotropy subgroups for Y2- ($Z_{2}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 14 | P2_1/c | C2h^5 | 2 | 1 |
///
pub struct Sg15;
