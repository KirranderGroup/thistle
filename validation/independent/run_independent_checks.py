# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Independent checks of every quantity THISTLE computes.

Each check compares a THISTLE result with a reference that does not go through THISTLE's
kernels: PySCF's own analytic Fourier transforms of basis-function pairs (ft_aopair) and its own
density matrices, numerical integration on PySCF's DFT grids, Lebedev averages over directions,
mpmath, the xraydb tables, closed-form limits, the energies printed by Molpro and OpenMolcas, the
earlier production kernels and thistle-dev's Fortran engine.  Every check prints one line:

  id | quantity | system | reference | max abs dev | max rel dev | tolerance | PASS/FAIL

Relative deviations are max|ours - ref| / max|ref| over the q points listed, unless stated.

python validation/independent/run_independent_checks.py   (about 10 minutes on 4 cores)
-> validation/independent/results.md and results.txt
"""
import os, sys, time, traceback
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))
DATA = os.path.join(ROOT, "tests", "data")

import pyxscat as xs
import thistle as th
from pyxscat import aligned as al, electron as el, iam, grids, ci as xci
from pyxscat.interfaces import pyscf as xp
from pyscf import gto, scf, cc, ci, mcscf, fci, dft
from pyscf.gto import ft_ao
from pyscf.dft import gen_grid

ROWS = []
T0 = time.time()


def record(cid, area, quantity, system, reference, ours, ref, tol, note="", absolute=False):
    ours = np.atleast_1d(np.asarray(ours, dtype=float)); ref = np.atleast_1d(np.asarray(ref, dtype=float))
    dabs = float(np.max(np.abs(ours - ref)))
    scale = float(np.max(np.abs(ref)))
    drel = dabs if absolute else dabs / scale if scale > 0 else dabs
    ok = drel <= tol
    ROWS.append(dict(id=cid, area=area, quantity=quantity, system=system, reference=reference,
                     dabs=dabs, drel=drel, tol=tol, ok=ok, note=note, absolute=absolute))
    print(f"{cid:6s} {'PASS' if ok else 'FAIL'}  {quantity:46.46s} {system:30.30s} abs {dabs:8.1e}  "
          f"{'abs' if absolute else 'rel'} {drel:8.1e}  tol {tol:7.0e}  | {reference}", flush=True)


def error(cid, area, quantity, system, exc):
    ROWS.append(dict(id=cid, area=area, quantity=quantity, system=system, reference="ERROR",
                     dabs=np.nan, drel=np.nan, tol=np.nan, ok=False, note=str(exc)[:200], absolute=False))
    print(f"{cid:6s} ERROR {quantity}: {exc}", flush=True)
    traceback.print_exc()


def lebedev(n):
    g = gen_grid.MakeAngularGrid(n)
    return g[:, :3], g[:, 3] / g[:, 3].sum()


def pair_ft(mol, qvec):
    """PySCF's analytic F_mn(q) = int chi_m chi_n exp(+i q.r)  (its convention is exp(-i G.r))."""
    return ft_ao.ft_aopair(mol, -np.asarray(qvec, dtype=float), aosym="s1")      # (nq, nao, nao)


def ref_elastic_fixed(mol, dm, qvec):
    F = pair_ft(mol, qvec)
    f = np.einsum("qmn,mn->q", F, dm)
    return np.abs(f) ** 2


def ref_two_electron_fixed(mol, dm2, qvec, chunk=600):
    n = mol.nao
    G = dm2.reshape(n * n, n * n)
    out = np.empty(len(qvec))
    for s in range(0, len(qvec), chunk):
        F = pair_ft(mol, qvec[s:s + chunk]).reshape(-1, n * n)
        out[s:s + chunk] = np.real(np.einsum("qa,qa->q", F @ G, F.conj()))
    return out


def sphere_avg(fun, q, ndir):
    """<fun(q n)> over directions n on a Lebedev grid, for each |q|."""
    n, w = lebedev(ndir)
    out = []
    for qq in q:
        out.append(w @ fun(qq * n))
    return np.array(out)


# =========================================================================== systems
print("building wavefunctions ...", flush=True)
WATER = "O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692"
wmol = gto.M(atom=WATER, basis="cc-pvdz", cart=True, verbose=0)
wmf = scf.RHF(wmol).run(conv_tol=1e-12)
wcc = cc.CCSD(wmf).run(conv_tol=1e-10); wcc.solve_lambda()
wci = ci.CISD(wmf).run(conv_tol=1e-10)
wcas = mcscf.CASSCF(wmf, 6, 6); wcas.conv_tol = 1e-11; wcas.kernel()
wsph = gto.M(atom=WATER, basis="cc-pvdz", cart=False, verbose=0)
wsph_mf = scf.RHF(wsph).run(conv_tol=1e-12)
comol = gto.M(atom="C 0 0 0; O 0 0 1.128", basis="cc-pvtz", cart=True, verbose=0)
comf = scf.RHF(comol).run(conv_tol=1e-12)
N_W = wmol.nelectron
Q = np.array([0.5, 1.5, 3.0, 5.0])            # bohr^-1, for the Lebedev checks
NDIR = 2030                                   # Lebedev degree 77
NDIR2 = 3074                                  # degree 95, to show the reference is converged
print(f"  done ({time.time()-T0:.0f} s)", flush=True)

AO = lambda mf: mf.make_rdm1()
cc_dm2_ao = wcc.make_rdm2(ao_repr=True)
ci_dm2_ao = wci.make_rdm2(ao_repr=True)
cas_dm1_ao, cas_dm2_ao = mcscf.addons.make_rdm12(wcas)

# =========================================================================== A. X-ray, rotational average
area = "A. X-ray, rotationally averaged"
try:
    m = th.load(wmf)
    ours = m.scattering(q=Q, unit="1/bohr", parts="all")
    ref = sphere_avg(lambda qv: ref_elastic_fixed(wmol, AO(wmf), qv), Q, NDIR)
    ref2 = sphere_avg(lambda qv: ref_elastic_fixed(wmol, AO(wmf), qv), Q, NDIR2)
    record("A1", area, "elastic <|F|^2>, RHF", "water/cc-pVDZ", "PySCF ft_aopair + Lebedev 2030", ours["elastic"], ref, 1e-9,
           note=f"reference converged: Lebedev 2030 vs 3074 differ by {np.max(np.abs(ref-ref2))/np.max(ref):.1e}")
    # Waller-Hartree for a closed-shell determinant: I_inel = N - 2 sum_ij <|F_ij|^2>, occupied MOs
    C = wmf.mo_coeff[:, wmf.mo_occ > 0]
    def wh(qv):
        F = pair_ft(wmol, qv); Fij = np.einsum("mi,qmn,nj->qij", C, F, C)
        return N_W - 2.0 * np.sum(np.abs(Fij) ** 2, axis=(1, 2))
    ref = sphere_avg(wh, Q, NDIR)
    record("A2", area, "inelastic, RHF (Waller-Hartree)", "water/cc-pVDZ", "N - 2 sum_ij <|F_ij|^2>, PySCF ft", ours["inelastic"], ref, 1e-9)
    record("A3", area, "total, RHF", "water/cc-pVDZ", "elastic + Waller-Hartree inelastic, PySCF", ours["total"],
           sphere_avg(lambda qv: ref_elastic_fixed(wmol, AO(wmf), qv), Q, NDIR) + ref, 1e-9)
    # f functions
    mco = th.load(comf); oc = mco.scattering(q=Q, unit="1/bohr", parts="all")
    C2 = comf.mo_coeff[:, comf.mo_occ > 0]
    def wh2(qv):
        F = pair_ft(comol, qv); Fij = np.einsum("mi,qmn,nj->qij", C2, F, C2)
        return comol.nelectron - 2.0 * np.sum(np.abs(Fij) ** 2, axis=(1, 2))
    record("A4", area, "elastic, RHF (f functions)", "CO/cc-pVTZ", "PySCF ft_aopair + Lebedev", oc["elastic"],
           sphere_avg(lambda qv: ref_elastic_fixed(comol, AO(comf), qv), Q, NDIR), 1e-9)
    record("A5", area, "inelastic, RHF (f functions)", "CO/cc-pVTZ", "Waller-Hartree, PySCF ft", oc["inelastic"], sphere_avg(wh2, Q, NDIR), 1e-9)
    # spherical basis (conversion to Cartesian functions inside THISTLE)
    ms = th.load(wsph_mf); os_ = ms.scattering(q=Q, unit="1/bohr", parts="all")
    record("A6", area, "elastic, RHF, spherical basis", "water/cc-pVDZ (5d)", "PySCF ft_aopair on spherical AOs", os_["elastic"],
           sphere_avg(lambda qv: ref_elastic_fixed(wsph, AO(wsph_mf), qv), Q, NDIR), 1e-9)
except Exception as e:
    error("A1-6", area, "RHF block", "", e)

try:
    for cid, label, obj, dm1, dm2, extra in (("A7", "CCSD", wcc, wcc.make_rdm1(ao_repr=True), cc_dm2_ao, {}),
                                              ("A8", "CISD", wci, wci.make_rdm1(ao_repr=True), ci_dm2_ao, {}),
                                              ("A9", "CASSCF(6,6)", wcas, cas_dm1_ao, cas_dm2_ao, {})):
        mm = th.load(obj)                  # factored or dense as THISTLE decides (here: dense, 25 functions)
        mf_ = th.load(obj, factored=True)  # the factored route explicitly
        o = mm.scattering(q=Q, unit="1/bohr", parts="all")
        of = mf_.scattering(q=Q, unit="1/bohr", parts="all")
        ref_tot = N_W + sphere_avg(lambda qv: ref_two_electron_fixed(wmol, dm2, qv), Q, NDIR)
        ref_el = sphere_avg(lambda qv: ref_elastic_fixed(wmol, dm1, qv), Q, NDIR)
        record(cid + "a", area, f"total, {label}", "water/cc-pVDZ", "N + <sum Gamma F F*>, PySCF Gamma + ft", o["total"], ref_tot, 1e-9)
        record(cid + "b", area, f"total, {label}, factored density", "water/cc-pVDZ", "same reference", of["total"], ref_tot, 1e-9)
        record(cid + "c", area, f"elastic, {label}", "water/cc-pVDZ", "<|sum D F|^2>, PySCF D + ft", o["elastic"], ref_el, 1e-9)
        record(cid + "d", area, f"inelastic, {label}", "water/cc-pVDZ", "difference of the two references", o["inelastic"], ref_tot - ref_el, 1e-9)
        # the cumulant against the reference with the PySCF density split by hand
        hart = ref_el
        if "cumulant" in o.curves:
            # exchange reference: -1/2 sum D_ps D_rq F_pq F*_rs
            def exch(qv):
                F = pair_ft(wmol, qv)
                return -0.5 * np.real(np.einsum("ps,rq,xpq,xrs->x", dm1, dm1, F, F.conj(), optimize=True))
            ref_x = sphere_avg(exch, Q, NDIR)
            record(cid + "e", area, f"exchange part, {label}", "water/cc-pVDZ", "-1/2 <sum D D F F*>, PySCF", o["exchange"], ref_x, 1e-9)
            record(cid + "f", area, f"cumulant part, {label}", "water/cc-pVDZ", "total - N - Hartree - exchange (refs)",
                   o["cumulant"], ref_tot - N_W - hart - ref_x, 1e-7, note="small quantity: relative to its own size")
except Exception as e:
    error("A7-9", area, "correlated block", "", e)

# =========================================================================== B. limits and sum rules
area = "B. Exact limits and sum rules"
try:
    q0 = np.array([0.0, 1e-6])
    mcc = th.load(wcc)
    r0 = mcc.scattering(q=q0, unit="1/bohr", parts="all")
    record("B1", area, "I_tot(0) = N^2", "water CCSD", "N^2 = 100", r0["total"], [100.0, 100.0], 1e-10)
    record("B2", area, "I_el(0) = N^2", "water CCSD", "N^2 = 100", r0["elastic"], [100.0, 100.0], 1e-10)
    record("B3", area, "I_inel(0) = 0", "water CCSD", "0", r0["inelastic"], [0.0, 0.0], 1e-9, absolute=True)
    qbig = np.array([8.0, 12.0, 16.0])
    rbig = mcc.scattering(q=qbig, unit="1/bohr", parts="all")
    ref_big = N_W + sphere_avg(lambda qv: ref_two_electron_fixed(wmol, cc_dm2_ao, qv), qbig, 5810)
    record("B4", area, "total at large q (8, 12, 16 bohr^-1)", "water CCSD", "N + <sum Gamma F F*>, PySCF, Lebedev 5810",
           rbig["total"], ref_big, 1e-9, note=f"two-electron term there: {', '.join(f'{x:.1e}' for x in ref_big - N_W)} "
           "(tight core pairs keep it finite in a Gaussian basis; it does not vanish)")
    d2 = th.load(wcc, factored=False).wf.dm2
    S = wmol.intor("int1e_ovlp")
    tr2 = np.einsum("pqrs,pq,rs->", wcc.make_rdm2(ao_repr=True), S, S)
    record("B5", area, "sum Gamma S S = N(N-1)", "water CCSD (PySCF Gamma)", "90", [tr2], [90.0], 1e-10)
except Exception as e:
    error("B1-5", area, "limits", "", e)

# =========================================================================== C. j2 moment / partial alignment
area = "C. j2 moment (aligned molecules)"
try:
    def p2(z): return 0.5 * (3 * z * z - 1)
    n2, w2 = lebedev(NDIR)
    for cid, label, obj, dm1, dm2 in (("C1", "RHF", wmf, AO(wmf), None), ("C2", "CCSD", wcc, wcc.make_rdm1(ao_repr=True), cc_dm2_ao)):
        a = th.load(obj).anisotropy(q=Q, unit="1/bohr")
        ref_el, ref_tot = [], []
        for qq in Q:
            qv = qq * n2
            iel = ref_elastic_fixed(wmol, dm1, qv)
            ref_el.append(w2 @ (iel * p2(n2[:, 2])))
            if dm2 is not None:
                it = N_W + ref_two_electron_fixed(wmol, dm2, qv)
            else:
                C = wmf.mo_coeff[:, wmf.mo_occ > 0]
                F = pair_ft(wmol, qv); Fij = np.einsum("mi,qmn,nj->qij", C, F, C)
                it = iel + N_W - 2.0 * np.sum(np.abs(Fij) ** 2, axis=(1, 2))
            ref_tot.append(w2 @ (it * p2(n2[:, 2])))
        record(cid + "a", area, f"S2 elastic, {label}", "water/cc-pVDZ", "<I_el(q n) P2(n_z)>, PySCF ft + Lebedev", a["S2 elastic"], ref_el, 1e-8)
        record(cid + "b", area, f"S2 total, {label}", "water/cc-pVDZ", "<I_tot(q n) P2(n_z)>, PySCF", a["S2 total"], ref_tot, 1e-8)
except Exception as e:
    error("C1-2", area, "j2", "", e)

# =========================================================================== D. fixed orientation and detector grids
area = "D. Fixed orientation, detector grids"
try:
    rng = np.random.default_rng(1)
    qv = rng.normal(size=(40, 3)); qv *= (rng.uniform(0.2, 5.0, size=40) / np.linalg.norm(qv, axis=1))[:, None]
    wcc_f = th.load(wcc).wf
    record("D1", area, "elastic |F(q)|^2, 40 random q vectors", "water CCSD", "PySCF ft_aopair", al.scattering(wcc_f, qv, "elastic"),
           ref_elastic_fixed(wmol, wcc.make_rdm1(ao_repr=True), qv), 1e-10)
    record("D2", area, "total, 40 random q vectors", "water CCSD", "N + sum Gamma F F*, PySCF", al.scattering(wcc_f, qv, "total"),
           N_W + ref_two_electron_fixed(wmol, cc_dm2_ao, qv), 1e-10)
    # Ewald sphere and photon wavelength
    from scipy import constants as sc
    E = 10.0
    k0_ref = E * 1e3 * sc.e / (sc.hbar * sc.c) * sc.physical_constants["Bohr radius"][0]
    record("D3", area, "photon k0 at 10 keV", "-", "E/(hbar c), CODATA (scipy)", [grids.photon_k0(E)], [k0_ref], 1e-6)
    tth = np.linspace(0.05, 1.0, 9); phi = np.linspace(0, 2 * np.pi, 9)
    qe = al.ewald_q(k0_ref, tth, phi)
    record("D4", area, "|q| on the Ewald sphere", "9 directions", "2 k0 sin(2theta/2)", np.linalg.norm(qe, axis=-1),
           2 * k0_ref * np.sin(tth / 2), 1e-12)
    # polar grid of the earlier grid.f90
    ref_grid = grids.read_qgrid(os.path.join(DATA, "grid", "qgrid_10keV_60deg_49_73.dat"))
    ours_grid = grids.qcoord3(10.0, 60.0, 49, 73)
    record("D5", area, "polar detector grid (qcoord3)", "10 keV, 60 deg, 49 x 73", "qgrid.dat written by the original grid.f90",
           np.asarray(ours_grid).ravel(), np.asarray(ref_grid).ravel(), 1e-11, absolute=True, note=f"{np.asarray(ref_grid).shape[0]} points; file has 12 decimals")
except Exception as e:
    error("D1-5", area, "fixed orientation", "", e)

try:
    # aligned ensemble: explicit orientational average with f(beta) = 1 + 5 <P2> P2(cos beta)
    P2v = 0.4
    m = th.load(wmf)
    a_ax = np.array([0.0, 0.0, 1.0])
    qlab = np.array([[1.2, 0.0, 0.3], [0.0, 2.0, 1.1], [0.7, -0.4, -2.5], [3.0, 1.0, 0.0]])
    qn = np.linalg.norm(qlab, axis=1)
    sc_ = m.scattering(q=qn, unit="1/bohr", parts="standard"); an = m.anisotropy(q=qn, unit="1/bohr")
    ours = grids.aligned_ensemble(qlab, a_ax, sc_["elastic"], an["S2 elastic"], p2=P2v)
    nd, wd = lebedev(590)
    gam = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    ref = np.zeros(len(qlab))
    D = AO(wmf)
    for ni, wi in zip(nd, wd):
        f = 1.0 + 5.0 * P2v * 0.5 * (3 * (ni @ a_ax) ** 2 - 1)
        t = np.array([1.0, 0, 0]) if abs(ni[0]) < 0.9 else np.array([0, 1.0, 0])
        e1 = np.cross(ni, t); e1 /= np.linalg.norm(e1); e2 = np.cross(ni, e1)
        qm = []
        for g_ in gam:
            u = np.cos(g_) * e1 + np.sin(g_) * e2; v = np.cross(ni, u)
            R = np.column_stack([u, v, ni])          # molecule frame axes in the lab: mol z along n
            qm.append(qlab @ R)                      # q in the molecular frame
        qm = np.concatenate(qm)
        I = ref_elastic_fixed(wmol, D, qm).reshape(len(gam), len(qlab)).mean(axis=0)
        ref += wi * f * I
    record("D6", area, "aligned ensemble, elastic, <P2> = 0.4", "water RHF, 4 lab q vectors",
           "explicit average over 590 x 24 orientations, PySCF ft", ours, ref, 1e-7)
except Exception as e:
    error("D6", area, "aligned ensemble", "", e)

# =========================================================================== E. electron scattering
area = "E. Electron scattering"
try:
    # electron-nucleus integral on PySCF's numerical DFT grid: no Gaussian integrals at all
    grid = dft.gen_grid.Grids(wmol); grid.level = 6; grid.build()
    ao = dft.numint.eval_ao(wmol, grid.coords)
    rho = dft.numint.eval_rho(wmol, ao, AO(wmf))
    Z, R = wmol.atom_charges().astype(float), wmol.atom_coords()
    Qe = np.array([0.0, 0.5, 1.5, 3.0, 5.0])
    def j0(x): return np.sinc(x / np.pi)
    ref_ne = np.array([sum(Z[a] * np.sum(grid.weights * rho * j0(qq * np.linalg.norm(grid.coords - R[a], axis=1)))
                           for a in range(len(Z))) for qq in Qe])
    wf = th.load(wmf).wf
    ours_ne = el.electron_nucleus(wf.basis, Qe, wf.dm1, Z, R)
    record("E1", area, "electron-nucleus term I_ne(q)", "water RHF", "sum_A Z_A int rho j0(q|r-R_A|), DFT grid (level 6)", ours_ne, ref_ne, 1e-7,
           note=f"grid integrates N to {abs(np.sum(grid.weights*rho)-N_W):.1e}")
    record("E2", area, "I_ne(0) = N sum Z", "water RHF", "10 x 10", [ours_ne[0]], [N_W * Z.sum()], 1e-10)
    # nucleus-nucleus: closed form
    Rab = np.linalg.norm(R[:, None] - R[None], axis=-1)
    ref_nn = np.array([np.sum(np.outer(Z, Z) * j0(qq * Rab)) for qq in Qe])
    record("E3", area, "nucleus-nucleus term", "water", "sum Z_A Z_B j0(q R_AB)", el.nucleus_nucleus(Z, R, Qe), ref_nn, 1e-14)
    # the whole bracket, against the X-ray references of section A
    d = el.decompose(th.load(wcc).wf, Q)
    ref_ee_tot = N_W + sphere_avg(lambda qv: ref_two_electron_fixed(wmol, cc_dm2_ao, qv), Q, NDIR)
    Dcc = wcc.make_rdm1(ao_repr=True)
    ref_ne_cc = sphere_avg(lambda qv: np.real(np.einsum("qmn,mn->q", pair_ft(wmol, qv), Dcc) *
                                              sum(Z[a] * np.exp(-1j * qv @ R[a]) for a in range(len(Z)))), Q, NDIR)
    ref_nn_q = np.array([np.sum(np.outer(Z, Z) * j0(qq * Rab)) for qq in Q])
    record("E4", area, "electron bracket, total, CCSD", "water/cc-pVDZ", "I_nn - 2 I_ne + I_ee, PySCF ft + Lebedev", d["total"],
           ref_nn_q - 2 * ref_ne_cc + ref_ee_tot, 1e-9)
    # small-q limit of the elastic bracket of a polar molecule: q^2 mu^2 / 3
    mu = np.asarray(wmf.dip_moment(unit="au", verbose=0))
    qs = np.array([2e-3])
    de = el.decompose(th.load(wmf).wf, qs, hf=True)
    record("E5", area, "elastic bracket -> q^2 mu^2/3 (q = 2e-3)", "water RHF", "PySCF SCF dipole", de["elastic"], qs ** 2 * (mu @ mu) / 3, 3e-3,
           note="leading term only; next order is O(q^4)")
    # cross section and electron wavelength
    from scipy import constants as sc
    Ek = 100.0                                   # keV
    mc2 = sc.m_e * sc.c ** 2
    p = np.sqrt((Ek * 1e3 * sc.e) ** 2 + 2 * Ek * 1e3 * sc.e * mc2) / sc.c
    k_ref = p / sc.hbar * sc.physical_constants["Bohr radius"][0]
    record("E6", area, "electron k0 at 100 keV (relativistic)", "-", "p/hbar, CODATA (scipy)", [grids.electron_k0(Ek)], [k_ref], 1e-6)
    br = np.array([0.3, 1.2]); qq = np.array([0.7, 2.0])
    record("E7", area, "cross section from the bracket", "-", "4 bracket / (a0^2 q^4) (bohr^2, no gamma)", el.cross_section(br, qq),
           4 * br / qq ** 4, 1e-14, note="units of the interface: bohr^2/sr for q in bohr^-1")
except Exception as e:
    error("E1-7", area, "electron", "", e)

# =========================================================================== F. hydrogen atom, exact
area = "F. Exact one-electron results"
try:
    hmol = gto.M(atom="H 0 0 0", basis="aug-cc-pv5z", spin=1, verbose=0)
    hmf = scf.UHF(hmol).run(conv_tol=1e-12)
    qh = np.array([0.5, 1.0, 2.0, 4.0])
    Fex = 16.0 / (4.0 + qh ** 2) ** 2
    dmh = np.asarray(hmf.make_rdm1()); dmh = dmh[0] + dmh[1]            # spin-summed AO density
    wfh = xp.from_density(hmol, dmh)                                       # open shell: no th.load yet
    el_h = xs.elastic(wfh.basis, qh, wfh.dm1)
    record("F1", area, "H atom elastic F(q)^2", "H/aug-cc-pV5Z UHF", "exact [16/(4+q^2)^2]^2", el_h, Fex ** 2, 3e-4,
           note=f"UHF energy error {hmf.e_tot + 0.5:.1e} Eh: basis-set limit")
    record("F2", area, "H atom inelastic 1 - F^2 (one electron: I_tot = 1)", "H/aug-cc-pV5Z UHF", "exact 1 - F(q)^2", 1.0 - el_h, 1 - Fex ** 2, 3e-4)
    Zh, Rh = np.array([1.0]), np.zeros((1, 3))
    br = el.nucleus_nucleus(Zh, Rh, qh) - 2 * el.electron_nucleus(wfh.basis, qh, wfh.dm1, Zh, Rh) + el_h
    record("F3", area, "H atom electron elastic bracket", "H/aug-cc-pV5Z UHF", "exact (1 - F)^2 (Mott-Bethe)", br, (1 - Fex) ** 2, 3e-4)
except Exception as e:
    error("F1-3", area, "H atom", "", e)

# =========================================================================== G. IAM
area = "G. Independent atom model"
try:
    import xraydb
    qi = np.array([0.5, 1.0, 2.0, 4.0, 6.0])
    s = qi / 0.52917721092 / (4 * np.pi)          # iam.py uses the CODATA 2010 Bohr radius
    devs = []
    for sym in ("H", "C", "N", "O", "F", "Ne"):
        ours = iam.atomic_form_factor(sym, qi, source="ab initio", basis="aug-cc-pvqz")
        ref = np.asarray(xraydb.f0(sym, s))
        z = gto.charge(sym)
        devs.append((sym, np.max(np.abs(ours - ref)) / z))
        record(f"G1{sym}", area, f"f(q) of the free atom {sym}", "HF/aug-cc-pVQZ", "xraydb f0 (Waasmaier-Kirfel)", ours, ref, 1e-2,
               note="method difference (HF in a finite basis vs relativistic HF tables): agreement to ~1% of Z expected")
    # Debye sum of the IAM against a sum written here with xraydb
    Zs = wmol.atom_charges(); R = wmol.atom_coords(); syms = [wmol.atom_symbol(i) for i in range(3)]
    f = np.array([xraydb.f0(sy, s) for sy in syms])
    Rab = np.linalg.norm(R[:, None] - R[None], axis=-1)
    ref = np.array([np.sum(np.outer(f[:, k], f[:, k]) * np.sinc(qq * Rab / np.pi)) for k, qq in enumerate(qi)])
    ours = iam.scattering(Zs, R, qi, source="tables", total=False)["elastic"]
    record("G2", area, "IAM elastic (Debye sum)", "water, xraydb f", "sum f_A f_B j0(q R_AB) written here", ours, ref, 1e-12)
except Exception as e:
    error("G1-2", area, "IAM", "", e)

# =========================================================================== H. radial functions
area = "H. Radial functions j_n(x)/x^n"
try:
    import mpmath as mp, ctypes
    from pyxscat._lib import lib
    mp.mp.dps = 40
    xs_ = np.logspace(-3, 3, 25); worst = 0.0
    for nmax in (4, 10, 18):
        for b in range(nmax + 1):
            coef = np.zeros(nmax + 1); coef[b] = 1.0
            out = np.zeros(xs_.size)
            lib.pyxscat_radial_series(xs_.size, np.ascontiguousarray(xs_), 1.0, nmax, coef, out)
            ref = np.array([float(mp.mpf(x) ** (2 * b) * mp.sqrt(mp.pi / (2 * mp.mpf(x))) * mp.besselj(b + mp.mpf(1) / 2, x) / mp.mpf(x) ** b) for x in xs_])
            worst = max(worst, float(np.max(np.abs(out - ref) / np.maximum(np.abs(ref), 1e-300))))
    record("H1", area, "q^2b j_b(x)/x^b, b <= 18, x = 1e-3 .. 1e3", "-", "mpmath, 40 digits", [worst], [0.0], 1e-13, absolute=True,
           note="pointwise relative error, worst over all orders and x")
except Exception as e:
    error("H1", area, "radial", "", e)

# =========================================================================== I. density matrices
area = "I. Density matrices"
try:
    from pyxscat.core import symmetrize_dm2
    # factored CCSD / CISD against PySCF's dense 2RDM (8-fold averaged, which is all the kernels see)
    for cid, label, obj, ref2 in (("I1", "CCSD", wcc, wcc.make_rdm2()), ("I2", "CISD", wci, wci.make_rdm2())):
        f = xp.from_cc(obj, factored=True) if label == "CCSD" else xp.from_ci(obj, factored=True)
        C = obj.mo_coeff; Ci = C.T @ wmol.intor("int1e_ovlp")
        mo = np.einsum("pqrs,ap,bq,cr,ds->abcd", f.dm2.dense(), Ci, Ci, Ci, Ci, optimize=True)
        record(cid, area, f"factored 2RDM, {label} (MO basis)", "water/cc-pVDZ", "PySCF make_rdm2", symmetrize_dm2(mo), symmetrize_dm2(ref2), 1e-12, absolute=True,
               note="max element difference")
    fC = th.load(wcas, factored=True).wf.dm2.dense()
    record("I3", area, "factored 2RDM, CASSCF(6,6) (AO basis)", "water/cc-pVDZ", "PySCF mcscf.addons.make_rdm12",
           symmetrize_dm2(fC), symmetrize_dm2(cas_dm2_ao), 1e-12, absolute=True)
    # RDM builders against PySCF FCI on a CAS(8,8) vector
    from pyscf.fci import cistring
    hmol8 = gto.M(atom="N 0 0 0; N 0 0 1.1", basis="6-31g", verbose=0)
    m8 = scf.RHF(hmol8).run(); c8 = mcscf.CASCI(m8, 8, 8).run()
    na, nb = c8.nelecas
    sa, sb = cistring.make_strings(range(8), na), cistring.make_strings(range(8), nb)
    dets = [(int(a), int(b)) for a in sa for b in sb]
    coeffs = np.asarray(c8.ci).ravel()
    d1r, d2r = fci.direct_spin1.make_rdm12(c8.ci, 8, (na, nb))
    for meth in ("dense", "bitwise", "sparse"):
        d1, d2 = xci.make_rdm12(8, na, nb, dets, coeffs, method=meth)
        record(f"I4{meth[0]}", area, f"RDM builder '{meth}', CAS(8,8) FCI vector", "N2/6-31G", "PySCF fci.direct_spin1.make_rdm12",
               np.concatenate([d1.ravel(), d2.ravel()]), np.concatenate([d1r.ravel(), d2r.ravel()]), 1e-12, absolute=True)
    # standalone FCI solver through th.load
    smol = gto.M(atom=WATER, basis="sto-3g", verbose=0); smf = scf.RHF(smol).run(conv_tol=1e-12)
    fs = fci.FCI(smf); fs.nroots = 2; fs.kernel()
    from pyscf.fci import direct_spin1
    mfci = th.load(fs, state=2, mo_coeff=smf)
    d1f, d2f = direct_spin1.make_rdm12(fs.ci[1], fs.norb, fs.nelec)
    Cm = smf.mo_coeff; hmol8 = smol
    w_ = mfci.wf
    ref_wf = xp.from_density(hmol8, d1f, d2f, mo_coeff=Cm)
    q4 = np.array([0.5, 2.0])
    record("I5", area, "two-electron term, FCI root 2 via th.load(fci.FCI)", "water/STO-3G", "same kernel on PySCF's root-2 RDMs",
           xs.two_electron(w_.basis, q4, w_.dm2), xs.two_electron(ref_wf.basis, q4, ref_wf.dm2), 1e-11,
           note="checks the loader (orbitals, root selection) with the kernel held fixed")
except Exception as e:
    error("I1-5", area, "density matrices", "", e)


try:
    # DMRG (block2) through th.load, against PySCF's exact CASCI density in the same orbitals
    sys.path.insert(0, os.path.join(ROOT, "benchmarks", "interface_comparison"))
    import importlib.util
    spec = importlib.util.spec_from_file_location("ci_", os.path.join(ROOT, "benchmarks", "interface_comparison", "compare_interfaces.py"))
    src = open(spec.origin).read()
    ns = {}
    exec(src[src.index("class DMRGCI"):src.index("# ------------------------------------------------------------------ wavefunctions")], {"np": np, "tempfile": __import__("tempfile")}, ns)
    dm = mcscf.CASCI(wmf, 6, 6); dm.fcisolver = ns["DMRGCI"](wmol, bond_dim=500); dm.kernel()
    ex = mcscf.CASCI(wmf, 6, 6); ex.kernel()
    _, ex_dm2_ao = mcscf.addons.make_rdm12(ex)
    od = th.load(dm, factored=True).scattering(q=Q, unit="1/bohr", parts="all")
    ref = N_W + sphere_avg(lambda qv: ref_two_electron_fixed(wmol, ex_dm2_ao, qv), Q, NDIR)
    a1d, a2d = dm.fcisolver.make_rdm12(0, 6, (3, 3)); a1e, a2e = ex.fcisolver.make_rdm12(ex.ci, 6, (3, 3))
    record("I6", area, "total, DMRG-CASCI(6,6) (block2, M = 500), factored", "water/cc-pVDZ",
           "PySCF exact CASCI Gamma + ft + Lebedev", od["total"], ref, 1e-8,
           note=f"DMRG - CASCI: energy {dm.e_tot - ex.e_tot:.1e} Eh, active 2RDM {np.abs(a2d - a2e).max():.1e} "
                "(DMRG convergence of the density, not THISTLE)")
except Exception as e:
    error("I6", area, "DMRG", "", e)

# =========================================================================== J. readers
area = "J. Readers (energies printed by the programs)"
try:
    from pyscf import ao2mo
    from pyscf.tools import molden as pm
    from pyxscat.interfaces import molpro, molcas
    info = molpro.read_multi_output(os.path.join(DATA, "molpro", "lif_sa4_casscf.out"))
    mol, _, moc, _, _, _ = pm.load(os.path.join(DATA, "molpro", "lif_sa4_casscf.mld")); mol.verbose = 0
    c = moc[:, :7]; h = c.T @ (mol.intor("int1e_kin") + mol.intor("int1e_nuc")) @ c
    eri = ao2mo.restore(1, ao2mo.kernel(mol, c), 7)
    es, er = [], []
    for k, s_ in sorted(info["states"].items()):
        d1, d2 = molpro.multi_rdms(info, k)
        es.append(mol.energy_nuc() + np.einsum("pq,pq", h, d1) + 0.5 * np.einsum("pqrs,pqrs", eri, d2)); er.append(s_["energy"])
    record("J1", area, "Molpro MULTI: 4 state energies from our RDMs", "LiF SA4-CASSCF(6,4)", "energies printed by Molpro", es, er, 2e-7, absolute=True,
           note="Hartree; limited by the 8 printed decimals")
    log = os.path.join(DATA, "molcas_mrci", "water_mrci.log"); mld = os.path.join(DATA, "molcas_mrci", "water_mrci.rasscf.molden")
    info = molcas.read_mrci_output(log); dm1, dm2 = molcas.mrci_rdms(info, 1)
    mol, _, moc, _, _, _ = pm.load(mld); mol.verbose = 0; c = moc[:, :24]
    h = c.T @ (mol.intor("int1e_kin") + mol.intor("int1e_nuc")) @ c
    eri = ao2mo.restore(1, ao2mo.kernel(mol, c), 24)
    e = np.einsum("pq,pq", h, dm1) + 0.5 * np.einsum("pqrs,pqrs", eri, dm2) + mol.energy_nuc()
    record("J2", area, "OpenMolcas MRCI: energy from our RDMs (12626 CSFs)", "water/cc-pVDZ", "CI energy printed by OpenMolcas",
           [e], [-76.23719469], 3e-7, absolute=True, note="limited by the 6 printed decimals of the CSF coefficients")
    # OpenMolcas RASSCF (C7H8 CASSCF(8,8), 1764 CSFs): Molcas used RICD integrals, so its energy cannot be
    # rebuilt with exact integrals (an exact PySCF CASCI in the same orbitals is 0.024 Eh lower); the
    # integral-free check is the natural occupations it prints
    logc = os.path.join(DATA, "molcas", "casscf_c7h8.log")
    infc = molcas.read_rasscf_output(logc); d1c, d2c = molcas.rasscf_rdms(infc)
    act = d1c[21:, 21:]
    record("J4", area, "OpenMolcas RASSCF: active natural occupations from our 1RDM", "C7H8/ANO-L-MB CASSCF(8,8)",
           "occupations printed by OpenMolcas", np.diag(act), infc["roots"][1]["occupations"], 1e-5, absolute=True,
           note="CI vector of 1764 CSFs -> determinants; printed to 4 decimals. Energy not usable: Molcas ran with RICD integrals")
    # Molden reader: elastic scattering against PySCF's own reader
    mldf = os.path.join(DATA, "molden", "co_uhf_pyscf.molden")
    molm, _, mo_c, mo_o, _, _ = pm.load(mldf); molm.verbose = 0
    if isinstance(mo_c, tuple) or (np.asarray(mo_c).ndim == 3):
        dmr = sum(np.asarray(c_) * np.asarray(o_) @ np.asarray(c_).T for c_, o_ in zip(mo_c, mo_o))
    else:
        dmr = (np.asarray(mo_c) * np.asarray(mo_o)) @ np.asarray(mo_c).T
    mm = th.load(mldf)
    record("J3", area, "Molden reader: elastic scattering", "CO UHF (PySCF Molden file)", "PySCF molden.load + ft_aopair + Lebedev",
           mm.scattering(q=Q, unit="1/bohr")["elastic"], sphere_avg(lambda qv: ref_elastic_fixed(molm, dmr, qv), Q, NDIR), 1e-9)
except Exception as e:
    error("J1-3", area, "readers", "", e)

# =========================================================================== K. other implementations
area = "K. Other implementations"
try:
    PIFAC = (float(np.arccos(np.float32(-1.0))) / np.pi) ** 3
    import glob
    from pyscf import cc as pcc
    for path in sorted(glob.glob(os.path.join(DATA, "legacy_*.npz"))):
        refd = np.load(path); name = os.path.basename(path)[7:-4]
        mol = gto.M(atom=str(refd["atom"]), basis=str(refd["basis"]), cart=True, verbose=0)
        mf = scf.RHF(mol).run(conv_tol=1e-12); myc = pcc.CCSD(mf).run(conv_tol=1e-10)
        basis = xp.basis_from_mole(mol)
        OLD = dict(screen_tol=0.0, cutoffcentre=0.01, cutoffz=1e-20, cutoffmd=1e-20)
        record(f"K1{name[:6]}", area, "earlier production kernel, HF, its own cut-off", name, "stored output of the earlier code",
               xs.mean_field(basis, refd["q"], mf.make_rdm1(), **OLD) * PIFAC, refd["hf"], 1e-9)
        wfc = xp.from_cc(myc)
        record(f"K2{name[:6]}", area, "earlier production kernel, CCSD 2-e, its own cut-off", name, "stored output of the earlier code",
               xs.two_electron(basis, refd["q"], wfc.dm2, **OLD) * PIFAC, refd["ccsd"], 1e-9)
    from pyxscat.interfaces import deck
    for nm in ("ch4_total", "ch4_total_zcotr"):
        wfd, info = deck.read(os.path.join(DATA, "thistle_dev", nm + ".h5"))
        gold = np.loadtxt(os.path.join(DATA, "thistle_dev", "ch4_total_golden.dat"))
        keep = np.abs(gold[:, 1]) > 1e-8 * np.abs(gold[:, 1]).max()
        ours = xs.two_electron(wfd.basis, gold[:, 0], wfd.dm2, screen_tol=1e-16, cutoffcentre=1e-2) * PIFAC
        record(f"K3{nm[-5:]}", area, f"thistle-dev Fortran engine ({'AO' if 'zcotr' in nm else 'MO'} deck)", "CH4/STO-3G",
               "thistle-dev golden output (ifort engine)", ours[keep], gold[keep, 1], 1e-11)
    data = np.loadtxt(os.path.join(DATA, "eirik_elastic_j2.txt"))
    molj = gto.M(atom="O 0.05 0.02 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", basis="6-31g", cart=True, verbose=0)
    wfj = xp.from_scf(scf.RHF(molj).run())
    j2 = -np.sqrt(5 / (4 * np.pi)) * xs.elastic(wfj.basis, data[:, 0], wfj.dm1, l=2, screen_tol=0.0, cutoffcentre=0.01, method="quartets")
    record("K4", area, "earlier j2 code (Y20 convention), its own cut-off", "water/6-31G", "stored output (E. Moen Liane)", j2, data[:, 1], 1e-6)
    tp = os.path.join(ROOT, "paper", "tab_three_programs_water_cas.txt")
    if os.path.exists(tp):
        print("K5    (three programs, from the stored run) " + open(tp).read().replace("\n", " | "))
except Exception as e:
    error("K1-5", area, "other implementations", "", e)

# =========================================================================== write
print(f"\n{sum(r['ok'] for r in ROWS)} of {len(ROWS)} checks pass ({time.time()-T0:.0f} s)")
with open(os.path.join(HERE, "results.txt"), "w") as fh:
    for r in ROWS:
        fh.write(f"{r['id']}\t{'PASS' if r['ok'] else 'FAIL'}\t{r['area']}\t{r['quantity']}\t{r['system']}\t{r['reference']}\t"
                 f"{r['dabs']:.2e}\t{r['drel']:.2e}\t{r['tol']:.0e}\t{'abs' if r['absolute'] else 'rel'}\t{r['note']}\n")
import json
json.dump(ROWS, open(os.path.join(HERE, "results.json"), "w"), indent=1, default=float)
