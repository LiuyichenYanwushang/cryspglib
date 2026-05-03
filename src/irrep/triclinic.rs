//! # Triclinic space groups (#1–#2)
//!
//! Triclinic crystals have no symmetry-imposed constraints on the lattice
//! parameters (a ≠ b ≠ c, α ≠ β ≠ γ ≠ 90°).  The Brillouin zone is a
//! triclinic parallelepiped with only the **Γ** point (0,0,0) possessing
//! non-trivial little-group symmetry.
//!
//! ## Brillouin zone
//!
//! | Label | Coords (fractional) | Little group |
//! |-------|---------------------|--------------|
//! | Γ | (0, 0, 0) | -1 (Ci) |
//!
//! ---


/// # 1 P1 (C₁¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1 | $\Gamma_{1}$ | $k_{8}t1$ | 1 | A1a | no |
///
/// ### Isotropy subgroups for GM1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 1 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1 | $B_{1}$ | $k_{5}t1$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1 ($B_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1 | $***$ | $k_{1}t1$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for R1 ($***$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1 | $***$ | $k_{2}t1$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1 ($***$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1 | $\Gamma_{1}$ | $k_{7}t1$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1 ($\Gamma_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1 | $***$ | $k_{3}t1$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for U1 ($***$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{V}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | V1 | $***$ | $k_{4}t1$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for V1 ($***$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1 | $F_{1}$ | $k_{6}t1$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1 ($F_{1}$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
pub struct Sg1;

/// # 2 P-1 (Cᵢ¹)
///
/// ## Irreps at $\Gamma\ (0,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | GM1+ | $\Gamma_{1}^+$ | $k_{8}t1$ | 1 | A1a | no |
/// | GM1- | $\Gamma_{1}^-$ | $k_{8}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for GM1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 1 | 1 |
///
/// ### Isotropy subgroups for GM1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 1 | P1 | C1^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{X}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | X1+ | $B_{1}^+$ | $k_{5}t1$ | 1 | A2a | no |
/// | X1- | $B_{1}^-$ | $k_{5}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for X1+ ($B_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ### Isotropy subgroups for X1- ($B_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{R}\ (1/2,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | R1+ | $***$ | $k_{1}t1$ | 1 | A2a | no |
/// | R1- | $***$ | $k_{1}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for R1+ ($***$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ### Isotropy subgroups for R1- ($***$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{T}\ (0,1/2,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | T1+ | $***$ | $k_{2}t1$ | 1 | A2a | no |
/// | T1- | $***$ | $k_{2}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for T1+ ($***$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ### Isotropy subgroups for T1- ($***$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Z}\ (0,0,1/2)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Z1+ | $\Gamma_{1}^+$ | $k_{7}t1$ | 1 | A2a | no |
/// | Z1- | $\Gamma_{1}^-$ | $k_{7}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Z1+ ($\Gamma_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Z1- ($\Gamma_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{U}\ (1/2,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | U1+ | $***$ | $k_{3}t1$ | 1 | A2a | no |
/// | U1- | $***$ | $k_{3}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for U1+ ($***$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ### Isotropy subgroups for U1- ($***$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{V}\ (1/2,0,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | V1+ | $***$ | $k_{4}t1$ | 1 | A2a | no |
/// | V1- | $***$ | $k_{4}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for V1+ ($***$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ### Isotropy subgroups for V1- ($***$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ## Irreps at $\mathrm{Y}\ (0,1/2,0)$
///
/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |
/// |-----|-----|---------|-----|-------|----------|
/// | Y1+ | $F_{1}^+$ | $k_{6}t1$ | 1 | A2a | no |
/// | Y1- | $F_{1}^-$ | $k_{6}t2$ | 1 | A2a | no |
///
/// ### Isotropy subgroups for Y1+ ($F_{1}^+$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
/// ### Isotropy subgroups for Y1- ($F_{1}^-$)
///
/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |
/// |-----------|-----|--------|-------------|---------|------|
/// | dir{(a)} | 2 | P-1 | Ci^1 | 2 | 1 |
///
pub struct Sg2;
