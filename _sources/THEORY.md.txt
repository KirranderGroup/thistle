# THISTLE: what is computed, and how

*Authors: Andrés Moreno Carrascosa, Mats Simmermacher, Patrick Yuheng Wang, Eirik Moen Liane, Joseph Cooper, Nikola Zotev, Thomas Northey and Adam Kirrander.*

THISTLE (with its kernel library `pyxscat`) computes the scattering of X-rays and fast electrons
by molecules directly from their electronic wavefunction, given as a Gaussian basis set and one-
and two-particle density matrices (1RDM, 2RDM). It does not depend on any electronic structure
program: PySCF objects, Molden files, OpenMolcas RASSCF and Molpro MULTI outputs are read by thin
adapters, and any other source only needs to provide a basis and density matrices.

Atomic units throughout unless stated (q in bohr⁻¹, lengths in bohr). The user interface
converts to Å⁻¹ and Å² by default.

---

## 1. X-ray scattering

### 1.1 Total, elastic and inelastic scattering

For a molecule in electronic state Ψ, the X-ray scattering in electron units (Thomson cross
section = 1), summed over all final electronic states (the closure or Waller–Hartree
approximation, valid when the photon energy greatly exceeds the electronic excitation energies),
is

    I_tot(q) = ⟨Ψ| Σ_ij exp(iq·(r_i − r_j)) |Ψ⟩ = N + ∫∫ Γ(r₁, r₂) exp(iq·(r₁ − r₂)) dr₁ dr₂,

with N the number of electrons and Γ the spin-summed pair density (normalised to N(N−1)).
The elastic part is the square of the form factor F(q) = ∫ ρ(r) exp(iq·r) dr,

    I_el(q) = |F(q)|²,         I_inel(q) = I_tot(q) − I_el(q).

In the Gaussian basis {χ_μ}, with ρ = Σ D_μν χ_μ χ_ν and Γ = Σ Γ_μνλσ χ_μ χ_ν (1) χ_λ χ_σ (2)
(chemists' order, PySCF's `make_rdm2` convention),

    F(q)  = Σ_μν D_μν F_μν(q),       F_μν(q) = ∫ χ_μ χ_ν exp(iq·r) dr,
    I₂(q) = Σ_μνλσ Γ_μνλσ ⟨μν|exp(iq·r₁₂)|λσ⟩.

### 1.2 Decomposition

Writing Γ = D⊗D − ½ D⊗D (exchange) + Λ, with Λ the cumulant (the genuinely correlated part),

    I_tot = N + I_Hartree + I_exchange + I_cumulant,     I_el = I_Hartree,

where the first three terms are built from the 1RDM alone. For a closed-shell Hartree–Fock
determinant Λ = 0, so total scattering follows from the orbitals (`hf_total`, `--hf`). For
natural orbitals of correlated wavefunctions it does not, which is why a Molden file alone only
gives elastic scattering unless the user declares it Hartree–Fock.

**Coupled cluster without the four-index density.** The CCSD 2RDM of PySCF is an n_MO⁴ then
n_AO⁴ array, several copies of it: benzene/cc-pVDZ (120 Cartesian functions) needs more than
6 GB, 500 functions would need 360 GB for the array alone. Its cumulant has structure: every
block carries occupied indices except the all-virtual one, and that block is exactly a sum over
occupied pairs of products,

    Λ_abcd (vvvv) = Σ_ij λ₂^{ij}_{ac} (τ^{ij} − ½ τ^{ij,T})_{bd} (+ images),    τ = t₂ + t₁t₁,

and the blocks with one occupied index (ovvv) are sums of products of the λ- and t-amplitudes
with one summed occupied index. The whole density is therefore written, exactly, as

    Γ_μνλσ = Σ_k A^C_k(μν) B^C_k(λσ)  +  Σ_k A^X_k(μλ) B^X_k(νσ),

with K_C = 2 n_occ² + 2 Coulomb-paired and K_X = 5 n_occ² + 2 exchange-paired terms (the
one-particle part is 4 P⊗(D − P) − 2 P⊗(D − P)_exchange with P the occupied projector), each
an n_AO × n_AO matrix: memory 2(K_C + K_X) n_AO² instead of n_AO⁴ (benzene 0.7 GB; 500 functions
and 40 occupied orbitals 45 GB). The intermediates are those of PySCF's `ccsd_rdm`
(passes 1 and 2, at most n_occ² n_vir²); pass 3, which builds the n_vir⁴ and n_occ n_vir³
blocks, is replaced by the factorisation. The kernels contract this form directly (mode 3): the
Coulomb-paired terms in the canonical shell-quartet loop (their block is 8-fold symmetric by
itself), the exchange-paired ones in a pass ordered by the (μ,λ) shell pair, where the block of
a batch of shell quartets is built from eight matrix products over k (BLAS `dgemm`, taken from
SciPy at run time) of factor blocks shared across the batch. It reproduces the dense route to
10⁻¹³ (every block of the MO density, and total, j₂, fixed-orientation and electron scattering;
frozen core included) and is the default beyond 60 basis functions. Timings on 2 cores:
water/cc-pVTZ 4.3 s against 2.2 s dense; benzene/cc-pVDZ 120 s in 2.4 GB, where the dense route
does not fit in 6 GB. The same representation serves any method whose 2RDM is given as products
of one-particle matrices.

### 1.3 Rotational average (gas phase)

For randomly oriented molecules, exp(iq·r) averages to j₀(qr) = sin(qr)/(qr):

    ⟨I(q)⟩ = N + Σ Γ_μνλσ ⟨μν|j₀(q r₁₂)|λσ⟩.

**Evaluation.** Each primitive product χ_a χ_b is expanded in Hermite Gaussians at its centre
P with exponent p (McMurchie–Davidson; p' and Q for the second pair). The Hermite coefficients come from the usual recurrence, valid for any
angular momentum. A primitive quartet then reduces to derivatives of j₀(q|P − Q|) with respect to
the vector between the two pair centres:

    ⟨ab|j₀(q r₁₂)|cd⟩ = (π²/pp')^{3/2} exp(−q²/4p − q²/4p') Σ E^{ab}_{tuv} E^{cd}_{τνφ} (−1)^{τ+ν+φ}
                        ∂^{t+τ}_x ∂^{u+ν}_y ∂^{v+φ}_z j₀(q|h|) at h = P − Q.

The Cartesian derivatives of j₀ are combinations Σ_b c_b j_b(qh)(q/h)^b = Σ_b c_b q^{2b} g_b(qh) with
g_b(x) = j_b(x)/x^b (b up to 4l_max + 2). Each order must be accurate for every x, and no single
recurrence direction achieves that: the upward recurrence is stable only for orders b < x, the
downward one for b > x. The radial functions are therefore evaluated as follows:

| x = qh | method |
|---|---|
| x < 1 | power series g_b = Σ_k (−x²/2)^k / [k!(2b+2k+1)!!] for the two highest orders, then the downward recurrence g_{b−1} = (2b+1) g_b − x² g_{b+1} |
| 1 ≤ x ≤ b_max + 1 | upward recurrence from j₀, j₁ for b ≤ x; Miller's backward recurrence for b > x, started where the neglected minimal solution is below 10⁻¹⁷ and matched to the upward values at the junction |
| x > b_max + 1 | upward recurrence for all orders |

Against 40-digit reference values every order up to 18 is accurate to 2·10⁻¹⁵ for
10⁻³ ≤ x ≤ 10³ (`tests/test_radial.py`). The earlier routine used the closed forms of j₁, j₂, j₃
for quartets of total order ≤ 3, which lose 2, 4 and 6 powers of x to cancellation (7·10⁻³
relative error at x = 10⁻³), and Miller's recurrence above, with its start capped at order 500
(wrong for x ≳ 480). Coincident centres need no special treatment; the closed form

    ∂^{L}_x ∂^{M}_y ∂^{N}_z j₀(q|h|) at h = 0 = (−1)^{S/2} q^S (L−1)!!(M−1)!!(N−1)!!/(S+1)!!,  S = L+M+N even,

is still used for h < `cutoffcentre` (10⁻¹⁰ bohr), where it is exact to machine precision.

**Density providers.** The quartet loop is the same for every quantity; only the weights
change: an explicit 2RDM, or products built on the fly from one-particle matrices
(Z = w_c/2 (A⊗B + B⊗A) + w_x A⊗A exchange). This gives elastic, transition and coherent cross
terms (w_c = 1, w_x = 0), exchange (w_c = 0, w_x = −1/4) and the Hartree–Fock mean field
(w_c = 1, w_x = −1/4) without storing any four-index array.

**Efficiency.** This is the "one-loop" algorithm of the earlier code (the latest version,
`Total_scattering_one_loop_opt.f90`): a canonical loop over contracted-shell quartets with the
primitive quartets inside, so each block of the pair density is gathered once per shell
quartet; quartets are screened by bounds on Gaussian prefactors and density weights
(`screen_tol`, 10⁻¹² by default); the outer loop is parallelised with OpenMP. The cost scales as
the number of significant primitive quartets × the number of q points.

**Elastic scattering by density-weighted pairs.** For the elastic term (and coherent/transition
terms) the pair density is a product, D⊗D, and a much cheaper route exists. The density is folded
into every primitive-shell pair first, turning it into one Hermite charge distribution at its
centre P,

    W_P(tuv) = Σ_{a∈P, b∈P} c_a c_b D_ab E^{ab}_{tuv},

and the elastic scattering becomes a generalised Debye sum over pairs of such distributions,

    ⟨|F(q)|²⟩ = Σ_{P,Q} pref_P pref_Q e^{−q²/4p − q²/4p'} Σ_m S_m(P − Q) (−q²)^m j_m(qh)/(qh)^m,

    S_m = Σ_{TUV} U_{TUV} Σ_{ijk: T+U+V−i−j−k=m} a_{Ti} a_{Uj} a_{Vk} h_x^{T−2i} h_y^{U−2j} h_z^{V−2k},
    U_{TUV} = Σ W_P(tuv) W_Q(t'u'v') (−1)^{t'+u'+v'}  (t+t' = T, …),   a_{Ti} = T!/(i!(T−2i)! 2^i),

where the second line is the closed form of the Cartesian derivatives of a radial function
(f_n = (1/h d/dh)^n f). Per pair of pairs the work is one convolution of two Hermite tables, a
q-independent radial sum and, for each q, one set of j_m(x)/x^m (evaluated as above). Because |∂^{tuv} j₀(q|h|)| ≤ q^{t+u+v}, the contribution at each q
is rigorously bounded by pref_P pref_Q e^{−q²/4p−q²/4p'} Σ_m |U_m| q^m, which drops the q tail of
diffuse pairs. Nothing is special at coincident or nearly coincident pair centres, so there is
no `cutoffcentre` approximation. This route is the default for elastic scattering:

| Elastic, 2 threads, 200 q | basis functions | pairs route | quartet route | speed-up |
|---|---|---|---|---|
| benzene / 6-31G | 66 | 9.8 s | 12.8 s | 1.3 |
| benzene / cc-pVDZ | 120 | 15.1 s | 25.6 s | 1.7 |
| benzene / aug-cc-pVDZ | 204 | 66 s | 132 s | 2.0 |
| CO / cc-pVQZ (7 q) | 110 | 0.57 s | 36 s | 63 |

The gain grows with angular momentum, because the quartet route loops over every Cartesian
component of four functions. Accuracy also improves: at the same tolerance (10⁻¹²) the
screening error is 10⁻⁹ instead of 10⁻⁶, I(0) = N² holds to machine precision, and the j₂ moment
is exact where the coincident-centre approximation of the earlier code (`cutoffcentre` = 0.01)
costs 2.5·10⁻⁴.

### 1.4 Partially aligned molecules: the j₂ moment

For molecules aligned along a laboratory axis (laser alignment, with axial symmetry), the signal
at a scattering vector making the angle θ_q with that axis is

    I(q, θ_q) = I₀(q) + 5 ⟨P₂(cos β)⟩ P₂(cos θ_q) S₂(q) + (l = 4, …),

with ⟨P₂⟩ the degree of alignment and

    S₂(q) = i² Σ Γ ⟨μν| j₂(q r₁₂) P₂(cos θ₁₂) |λσ⟩ = −Σ Γ ⟨μν| j₂(q r₁₂) P₂(cos θ₁₂) |λσ⟩,

θ₁₂ being the angle of r₁₂ with the molecular z axis. The factor i² = −1 is the phase of the
l = 2 term of the plane-wave expansion, e^{iq·r} = Σ_l (2l+1) iˡ j_l(qr) P_l(q̂·r̂); with it
included every moment enters with a plus sign, I = Σ_l (2l+1) ⟨P_l⟩ P_l(cos θ_q) S_l, and S₂ is
directly the projection S₂ = ⟨I(q) P₂(cos θ_q)⟩ over directions of q for perfect alignment
(= a₂/5, a₂ the l = 2 Legendre coefficient). (Before 8 October 2026 pyxscat returned −S₂, with
I = I₀ − 5⟨P₂⟩P₂S₂; the numbers are the same up to that sign.) The j₂P₂ kernel is reduced to j₀ by
Hobson's theorem,

    j₂(qh) P₂(cos θ_h) = [3 ∂²/∂h_z² + q²] j₀(qh) / (2q²),

so the same derivative machinery serves (two more orders), with closed-form coincident limits.
The isotropic self term N does not contribute to S₂. The convention of the earlier code (with
Y₂₀ and without i²) differs by a factor −√(5/4π).

### 1.5 Fixed orientation (oriented molecules, single-shot detector images)

Since exp(iq·(r₁ − r₂)) = exp(iq·r₁) exp(−iq·r₂), all fixed-orientation quantities follow from
the analytic pair transforms

    F_ab(q) = (π/p)^{3/2} exp(−q²/4p) exp(iq·P) Σ E^{ab}_{tuv} (iq_x)^t (iq_y)^u (iq_z)^v,

contracted with density matrices: elastic |Σ D F|², two-electron Σ Γ F F* (one matrix product
per q), exchange and mean field from the 1RDM at O(n³) per q. There is no Bessel expansion and no
coincident-centre special case. Detector geometry: q vectors on the Ewald sphere for a given
wavelength (beam along z), and molecular orientations as z-y-z Euler angles.

### 1.6 Detector grids and aligned ensembles in the laboratory

The beam runs along z with wavenumber k₀ (k₀ = E·0.26817 bohr⁻¹/keV for photons; relativistic
2π/λ for electrons), and the scattered wave vectors lie on the Ewald sphere. Besides a square
(q_x, q_y) image, results can be computed on the **polar grid of the earlier `grid.f90`
(`qcoord3`)**, reproduced point for point: n₁ rings of transverse momentum
q⊥,i = i·k₀ sin θ_max / n₁, n_φ equally spaced azimuths per ring (φ from the y axis:
q_x = q⊥ sin φ, q_y = q⊥ cos φ), q_z = k₀[1 − (1 − q⊥²/k₀²)^{1/2}], and q = 0 first
(N = n₁ n_φ + 1 points; `grid.f90`'s n₂ is n_φ + 1). Its sign is q = k_in − k_out; since
I(−q) = I(q) for real wavefunctions and the alignment weight below is even, the patterns are the
same in either convention. Any other list of laboratory vectors can be read from a file
(columns q_x q_y q_z in bohr⁻¹).

*Fixed orientation.* For Euler angles (α, β, γ), q_mol = Rᵀ q_lab and §1.5 applies point by point.

*Partially aligned ensemble.* For an alignment axis â (e.g. the laser polarisation) and a degree
of alignment ⟨P₂(cos β)⟩ of the molecular z axis, §1.4 gives at every grid point

    I(q) = I₀(|q|) + 5 ⟨P₂⟩ P₂(q̂·â) S₂(|q|),

with I₀ and S₂ computed once per distinct |q| (once per ring on the polar grid). Terms l ≥ 4
(⟨P₄⟩ S₄, …) are neglected; the relation is exact for any distribution with no components beyond
l = 2, e.g. f(β) ∝ 1 + 5⟨P₂⟩P₂(cos β). The same holds for electron cross sections.

*IAM on the same grid.* At fixed orientation the IAM is |Σ_A f_A(q) e^{iq·R_A}|² plus the atomic
inelastic terms (electrons: amplitudes Z_A − f_A); for the aligned ensemble its I₀ and S₂ are the
Debye forms of §3 at l = 0 and 2. Molecule and IAM are thus compared point by point, as
100 (I − I_IAM)/I_IAM.


---

## 2. Electron scattering

### 2.1 Cross sections

A fast electron is scattered by the whole charge distribution, nuclei and electrons. In the
first Born approximation, with Z(q) = Σ_A Z_A exp(iq·R_A) the nuclear form factor,

    dσ/dΩ = (4 / a₀² q⁴) [ |Z(q) − F(q)|² + I_inel(q) ],

where I_inel is the same inelastic term as for X-rays (nuclei do not scatter inelastically at
these energies). Expanding the square,

    dσ/dΩ = (4 / a₀² q⁴) [ I_nn − 2 I_ne + I_ee ],

    I_nn = Σ_AB Z_A Z_B exp(iq·R_AB)                (nucleus–nucleus)
    I_ne = Σ_A Z_A Re[exp(−iq·R_A) F(q)]             (electron–nucleus)
    I_ee = I_tot (X-ray, total) or |F|² (elastic)    (electron–electron)

The code reports the bracket (as the earlier code did) at the kernel level and dσ/dΩ in the
interface (Å²/sr for q in Å⁻¹; bohr²/sr for q in bohr⁻¹). For a neutral molecule the bracket
vanishes at q = 0; for a polar molecule the elastic bracket starts as q² μ²/3 (μ the dipole
moment), so dσ/dΩ diverges as 1/q² and the q = 0 point is undefined. For fast beams, the
relativistic correction multiplies the cross section by γ² = (1 + E/mc²)²; the wavelength
λ = h/p is relativistic as well (`electron_wavelength`).

### 2.2 Rotational average and the electron–nucleus integrals

On rotational averaging, I_nn = Σ Z_A Z_B j₀(qR_AB) and

    I_ne(q) = Σ_A Z_A Σ_μν D_μν ∫ χ_μ χ_ν j₀(q|r − R_A|) dr.

This one-electron integral is evaluated analytically. With the Hermite expansion of χ_a χ_b at P,

    ∫ Λ_tuv(r − P) j₀(q|r − A|) dr = (π/p)^{3/2} exp(−q²/4p) ∂^t_x ∂^u_y ∂^v_z j₀(q|h|),  h = P − A,

and the derivatives of the radial function follow from the McMurchie–Davidson recursion

    R^n_{000} = f_n(h),   R^n_{t+1,u,v} = t R^{n+1}_{t−1,u,v} + h_x R^{n+1}_{t,u,v}   (and y, z),

with the auxiliary functions

    f_n(h) = (1/h d/dh)^n j₀(qh) = (−q²)^n j_n(qh)/(qh)^n.

The functions j_n(x)/x^n are evaluated by their power series for x < 1, by Miller's backward
recurrence for 1 ≤ x ≤ n_max + 1 and by upward recurrence beyond, so they are accurate for every
x including x = 0. Nothing is singular when a nucleus sits on a pair centre. The j₂P₂ moment for
aligned molecules uses Hobson's theorem as in §1.4, and I_nn at l = 2 is
−Σ Z_A Z_B j₂(qR_AB) P₂(cos θ_AB). The relation I(q, θ_q) = I₀ + 5⟨P₂⟩P₂(cos θ_q) S₂ holds term by
term, so it holds for electrons as for X-rays.

The cost is one pass over primitive pairs × atoms × q points, negligible next to the
two-electron terms (0.7 s for benzene/cc-pVDZ at 200 q on two cores, against 30 s for the X-ray
elastic term).

### 2.3 Fixed orientation

At fixed orientation, Z(q) is a sum of plane waves and F(q) comes from the pair transforms of
§1.5, so |Z − F|² and the inelastic term are exact for any list of q vectors (detector images
for UED, given the beam energy).

### 2.4 Modified scattering intensity sM(s)

Electron diffraction experiments report the modified intensity

    sM(s) = s (I − I_at) / I_at,      s = q,

where I_at is the incoherent sum of atomic scattering. THISTLE gives it for the ab initio
molecule (with I_at from the independent atom model below) and for the IAM itself
(s I_mol/I_at, the usual simulation). The difference is the effect of chemical bonding and
electron correlation on the diffraction pattern.

---

## 3. Independent atom model (IAM)

The reference against which molecular effects are measured: free, spherically averaged atoms at
the nuclear positions.

    X-rays     I_el = Σ_AB f_A f_B j₀(qR_AB),     I_tot = Σ_A S_A + Σ_{A≠B} f_A f_B j₀(qR_AB)
    electrons  I_el = Σ_AB (Z_A − f_A)(Z_B − f_B) j₀(qR_AB),   + Σ_A (S_A − f_A²) inelastic

(the electron amplitude Z_A − f_A is the Mott–Bethe formula). At fixed orientation j₀(qR) becomes
exp(iq·R_AB); for l = 2, j₂(qR_AB) P₂(cos θ_AB).

The atomic form factors f_A and total scattering S_A = N_A + I₂,A come either from ab initio
free-atom calculations through the same kernels (UHF or closed-shell CCSD, in a chosen basis per
element; the spherical average is exact on a Lebedev grid of degree 17 because the angular
structure of an atomic density is a polynomial of low degree), or from the Waasmaier–Kirfel
tables (elastic only). Using the molecule's own method and basis isolates the effect of bonding.

---

## 4. Wavefunctions and conventions

| Source | Read through | Gives |
|---|---|---|
| PySCF, in memory | `interfaces.pyscf` | RHF, CCSD, CISD, FCI, CASSCF/CASCI (state-averaged: any root) |
| Molden file (Molpro, OpenMolcas, BAGEL, ORCA, PySCF, …) | `interfaces.molden` | basis, orbitals, occupations |
| OpenMolcas RASSCF output + `.rasscf.molden` | `interfaces.molcas` | CI vector of any root (GUGA CSFs → determinants) |
| OpenMolcas MRCI output (uncontracted SDCI/ACPF, `THRPrint = 0.0`) + `.rasscf.molden` | `interfaces.molcas` | the complete MRCI vector of any state: CSFs over inactive, active and virtual orbitals (≤ 64 correlated orbitals) |
| Molpro MULTI output (`config,det`) or punch file + Molden | `interfaces.molpro` | CI vectors of all states, transition 1RDMs |
| Your own CI vector | `pyxscat.ci` | spin-summed 1RDM, 2RDM, transition 1RDM |

* Spherical basis functions (Molden `[5D]`, `[7F]`, `[9G]`) are converted to Cartesian ones by
  real solid harmonics; the Cartesian normalisation convention of the Molden file is detected
  from MO orthonormality. Molpro orbital order is restored from the symmetry labels.
* Molpro's orbital-ordered determinants carry the phase that converts them to α-string-first
  order; truncated CI printouts are detected (and the punch file used instead when present).
* Active-space RDMs are embedded with the doubly occupied core. Beyond 60 basis functions (or
  with `factored=True`) the embedded 2RDM is never formed as an n⁴ array but kept as products of
  one-particle matrices, Γ = Σ_k w_k V_k ⊗ V_k + D_c⊗D_c + D_c⊗D_a + D_a⊗D_c − ½(exchange
  pairings of the same three), where V_k, w_k diagonalise the active Γ as a symmetric
  (pq)×(rs) matrix (at most n_act(n_act+1)/2 terms: 21 for CAS(6,6), 136 for CAS(16,16)) and
  D_c, D_a are the core and active 1RDMs. It is exact (same scattering to 10⁻¹⁴) and is
  contracted by the factored kernel (mode 3, as for CCSD). Memory 2K n_AO² instead of n_AO⁴:
  benzene CAS(6,6)/cc-pVTZ needs 39 MB of factors instead of a 65 GB array (peak 0.33 GB for the
  whole run instead of > 200 GB). Molpro, Molcas, PySCF CASSCF/CASCI and FCI all use it.
* **RDMs from determinant expansions** (`pyxscat.ci.make_rdm12`), three exact builders, chosen
  automatically by an estimate of their cost (all agree to 10⁻¹⁵ and with PySCF):

  | Builder | How | Cost | Suited to |
  |---|---|---|---|
  | dense | Γ_pqrs = ⟨E_qp Ψ\|E_rs Ψ⟩ − δ_qr γ_ps, intermediates E_rs\|Ψ⟩ over the complete active space (α × β strings, in blocks), one matrix product | norb⁴ N_FCI (BLAS) | CASSCF/CASCI vectors: CAS(12,12) 1.6 s, CAS(13,13) 6 s |
  | bitwise | the determinant-pair algorithm of the earlier code (Fortran, OpenMP): pairs screened by popcnt(ieor), Slater–Condon rules for ≤ 2 differences; multiword bit strings, so any number of orbitals (the RDMs accumulate in thread-private copies while they are small, otherwise in the shared array with atomic updates) | ≈ 1 ns per pair and thread, N_det²/2 pairs; independent of N_FCI | truncated and selected expansions (MCCI, print thresholds, RAS) in large orbital spaces: 10⁵ determinants in 27 orbitals 2.9 s on 2 cores; a CAS vector on scattered orbitals of a 128-orbital space reproduces PySCF to 10⁻¹⁴ |
  | external (`make_rdm12_external`) | for expansions with at most two electrons in an "external" set of orbitals (uncontracted MRCI, CISD-like): determinants grouped by their external label; the pair loop runs only over label pairs that can couple (differing by ≤ 3 spin-orbitals), and the external-external-external-external block, the only one fed by labels differing by 4, is a matrix product over the internal parts, Γ_vvvv[c,a,d,b] = Σ_K Σ_στ T_K^στ[c,d] T_K^στ[a,b] | pairs ∝ N_det × (neighbouring labels), plus GEMMs | MRCI: water/cc-pVDZ (50 034 CSFs) and a 120-orbital test agree with the full pair loop to 10⁻¹³; 4× faster at 70 orbitals with an MRCI-like label distribution, and the gain grows with the virtual space |
  | sparse | the same intermediates as dense, but only on the determinants reached (bit-string excitations, sparse product, bounded memory) | linear in N_det | very long expansions on few cores (4·10⁵ determinants: 55 s, against 264 s for the pair loop on 2 cores) |

  Measured on 2 cores (`benchmarks/ci_rdm_scaling`): a full CAS(14,14) vector (1.2·10⁷
  determinants) takes 39 s with the dense builder (peak 3.2 GB); FCI of water/6-31G (1.7·10⁶
  determinants) 4.9 s. The pair loop needs N_det²/2 pair tests whatever the language: in Fortran
  ~1 ns per pair (10⁵ determinants in 3 s), in pure Python 170 ns per pair for the *screening
  alone* (0.1 h at 6·10⁴ determinants, 64 h at 1.7·10⁶); that is why the earlier
  pure-Python version stopped at ~10³ determinants. None of the builders is pure Python: the
  dense and sparse ones do their work in BLAS/SciPy, the bitwise one in Fortran.
* Nuclear charges for electron scattering are the atomic numbers (pass `charges=` to override,
  e.g. with effective core potentials).

---

## 5. Validation

| What | Against | Agreement |
|---|---|---|
| Rotationally averaged kernel (elastic, exchange, CCSD two-electron) | independent route: analytic pair transforms (equal to PySCF `ft_aopair`) averaged on a 131-point Lebedev grid | 10⁻⁸ |
| Radial functions j_b(x)/x^b, b ≤ 18, 10⁻⁴ ≤ x ≤ 10³ | 40-digit values (mpmath); SciPy | 2·10⁻¹⁵ |
| Elastic by density-weighted pairs (l = 0 and 2, q from 10⁻³; transition densities) | the same Lebedev route; the quartet kernel up to f functions | 10⁻¹¹; 10⁻⁹ |
| j₂ at small q through the quartet kernel (exchange, CCSD two-electron, elastic) | Lebedev projection of the fixed-orientation terms; the pair route | 10⁻⁹ for q ≥ 0.01 bohr⁻¹, 10⁻⁶ at 10⁻³ |
| The earlier code (BigMemoryScattering kernel) | stored reference curves (HF and CCSD, four molecules) | 10⁻⁹ (after its single-precision π) |
| j₂ moment | Lebedev projection of the fixed-orientation signal; the earlier code (Y₂₀ convention); IAM for separated atoms (Ne₂ at 14 bohr) | 5·10⁻⁸; 10⁻⁶; 10⁻⁸ |
| Electron–nucleus integrals (l = 0 and 2) | Lebedev average of Re[Z* F] from pair transforms, up to g functions, Cartesian and spherical, q up to 12 bohr⁻¹ | 10⁻¹⁰ |
| Electron–nucleus integrals | the earlier code (`nuclei_electron_integration`) | 4·10⁻⁸ (its single-precision π) |
| Electron scattering of one atom | Mott–Bethe: elastic = (Z − f)² with f² the X-ray elastic | 10⁻⁸ |
| H atom, elastic electron scattering | exact [1 − 16/(4 + q²)²]² | 2·10⁻⁴ (basis set) |
| Polar molecule at small q | elastic bracket → q² μ²/3 with the SCF dipole | 2·10⁻³ at q = 2·10⁻³ |
| Fixed-orientation electron and X-ray scattering | rotational kernels after Lebedev averaging (l = 0 and l = 2) | 10⁻⁹–10⁻⁸ |
| Partially aligned ensemble on a detector grid (X-ray elastic and HF total, electron elastic) | explicit orientational average of fixed-orientation patterns with f(β) ∝ 1 + 5⟨P₂⟩P₂(cos β) (Lebedev order 59 × 24 γ) | 10⁻¹²–10⁻¹⁰ |
| Polar detector grid | `qgrid.dat` written by the original `grid.f90` | identical file |
| IAM on detector grids: fixed orientation (electrons), aligned ensemble (X-rays, electrons) | Lebedev average to the Debye form; explicit orientational average as above | 10⁻¹⁰; 10⁻⁹ |
| CCSD density as a sum of products (kernel mode 3), all routes | PySCF `make_rdm2` block by block (10⁻¹⁵); the dense route for total, cumulant, j₂, fixed orientation, electrons (10⁻¹³) | 10⁻¹³ |
| RDMs from determinants: dense, bitwise and sparse builders | PySCF FCI (`direct_spin1`) and selected CI (`selected_ci`, 27 orbitals); each other up to 64 orbitals | 10⁻¹⁵ |
| Molden reader | PySCF round trip, up to g functions | density 10⁻¹¹, scattering 10⁻⁹ |
| h functions (l = 5): elastic (pairs, quartets), exchange, l = 0 and 2; electron–nucleus | Lebedev route | 10⁻¹¹ (pairs, electron–nucleus); 10⁻⁸ (quartets) |
| Molcas and Molpro adapters | state energies rebuilt from our RDMs | 10⁻⁸ (Molpro), CASCI eigenvector (Molcas) |
| IAM | a closed-shell atom through the molecular kernels; Waasmaier–Kirfel tables (H to S, aug-cc-pVQZ) | 10⁻⁸; within 1 % of Z |

Exact limits are tested throughout: I_tot(0) = I_el(0) = N², I_ne(0) = Σ Z_A · N, the
decompositions sum to the totals, and electron brackets vanish at q = 0 for neutral molecules.
The suite (`python -m pytest`, about 5 minutes) runs 217 tests; `PYXSCAT_SLOW=1` adds 5 on large basis sets.

---

## 6. Limits and plans

* Angular momentum: tested up to h functions (l = 5) in every kernel. The quartet kernel's
  derivative tables are dimensioned for 4 l_max + 3 ≤ 32 (l ≤ 7); the pair, electron–nucleus and
  fixed-orientation kernels have no fixed limit. The Molden format defines shells only up to g, so
  h functions come in through PySCF (or any adapter that supplies the basis directly).
* Closure (Waller–Hartree) and first Born approximations; no Compton energy resolution, no
  multiple scattering.
* The j₂ moment is computed through Hobson's theorem, [3∂²_z + q²] j₀/(2q²), which cancels as
  q → 0: with the stable radial functions the relative error is 10⁻¹¹ at q = 0.01 bohr⁻¹ and
  4·10⁻⁸ at 10⁻³ (the elastic j₂ through the pair route: 10⁻¹³).
* j₂ / aligned ensembles: S₂ is the moment about the **molecular z axis** of the input frame
  (the laboratory alignment axis is free). A molecule whose aligned axis is not z must be
  reoriented before the quantum-chemistry run; nothing warns about it. To do: a
  `molecular_axis` option (rotate basis and densities, exact for Cartesian Gaussians) or the full
  l = 2 tensor T (S₂(n) = n·T·n, five components, any axis afterwards). Only one-dimensional
  alignment (⟨P₂⟩ of one axis); 3D alignment of asymmetric tops would need higher terms.
* Planned: j₁ terms; BAGEL and NEXMD adapters. (Factored CCSD and CAS/CI densities: done.)
