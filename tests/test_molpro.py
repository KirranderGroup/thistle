"""Molpro MULTI interface on a stored state-averaged CASSCF(6,4) of LiF, four singlet states
(aug-cc-pVTZ on F, cc-pVTZ on Li, Molpro 2022, config,det)."""
import os

import numpy as np
import pytest

import pyxscat as xs
from pyxscat import aligned
from pyxscat.interfaces import molpro

D = os.path.join(os.path.dirname(__file__), "data", "molpro")
OUT = os.path.join(D, "lif_sa4_casscf.out")
MLD = os.path.join(D, "lif_sa4_casscf.mld")


@pytest.fixture(scope="module")
def info():
    return molpro.read_multi_output(OUT)


@pytest.fixture(scope="module")
def ints():
    from pyscf import ao2mo
    from pyscf.tools import molden as pm
    mol, _, moc, _, _, _ = pm.load(MLD)
    c = moc[:, :7]
    h = c.T @ (mol.intor("int1e_kin") + mol.intor("int1e_nuc")) @ c
    eri = ao2mo.restore(1, ao2mo.kernel(mol, c), 7)
    r = np.array([c.T @ x @ c for x in mol.intor("int1e_r")])
    return mol, h, eri, r


def test_parsing(info):
    assert (info["nclosed"], info["nactive"], info["nactive_electrons"], info["ms2"]) == (3, 4, 6, 0)
    assert sorted(info["states"]) == [1, 2, 3, 4]
    assert info["states"][1]["energy"] == pytest.approx(-106.90987819)


def test_state_energies_reproduced(info, ints):
    """E = h.D + 1/2 (pq|rs) Gamma + E_nuc with PySCF integrals in the Molden orbitals matches
    Molpro for every state: checks the determinant phases, the RDMs and the reader together."""
    mol, h, eri, _ = ints
    for k, s in info["states"].items():
        d1, d2 = molpro.multi_rdms(info, k)
        e = mol.energy_nuc() + np.einsum("pq,pq", h, d1) + 0.5 * np.einsum("pqrs,pqrs", eri, d2)
        assert e == pytest.approx(s["energy"], abs=2e-7)


def test_state_averaged_occupations(info):
    avg = sum(np.diag(molpro.multi_rdms(info, k)[0])[3:] for k in info["states"]) / 4
    assert np.allclose(avg, [1.75, 1.75, 1.75, 0.75], atol=1e-6)      # Molpro natural occupations


def test_dipoles_and_transition_dipoles(info, ints):
    mol, _, _, r = ints
    d1, _ = molpro.multi_rdms(info, 1)
    mu = -np.einsum("pq,pq", r[2], d1) + mol.atom_charges() @ mol.atom_coords()[:, 2]
    assert mu == pytest.approx(-2.69917604, abs=1e-7)
    for (a, b, comp), ref in {(1, 4, 2): -0.595447349279, (1, 3, 0): -0.223917907042,
                              (1, 2, 1): -0.223917907042}.items():
        t = molpro.transition_dm1(info, a, b)
        assert -np.einsum("pq,pq", r[comp], t) == pytest.approx(ref, abs=1e-7)


def test_wavefunction_and_fixed_orientation_scattering(info):
    wf = molpro.from_multi(OUT, MLD, state=1)
    assert wf.nelec == 12
    q = np.array([[0.0, 0.0, 0.0], [0.3, 0.1, 1.2]])
    d = aligned.decompose(wf, q)
    assert d["total"][0] == pytest.approx(144, rel=1e-8)
    assert d["elastic"][0] == pytest.approx(144, rel=1e-8)
    assert d["inelastic"][1] > 0


def test_transition_scattering_vanishes_at_q0(info):
    """Orthogonal states: the transition form factor is zero at q = 0 and finite elsewhere."""
    from pyxscat.interfaces.molden import read_molden
    data = read_molden(MLD)
    c = data.mo_coeff[:, :7]
    t = c @ molpro.transition_dm1(info, 1, 4) @ c.T
    f = aligned.form_factor(data.basis, np.array([[0, 0, 0.0], [0, 0, 1.0]]), t)
    assert abs(f[0]) < 1e-7 and abs(f[1]) > 1e-3


def test_csf_printout_is_rejected(tmp_path):
    p = tmp_path / "x.out"
    p.write_text(" PROGRAM * MULTI\n Number of closed-shell orbitals:   1 (1)\n Number of active  orbitals: 2 (2)\n"
                 " Number of active electrons:   2    Spin symmetry=Singlet\n CI vector\n 2/\\0 0.9\n")
    with pytest.raises(ValueError):
        molpro.read_multi_output(str(p))


SF6P = os.path.join(D, "sf6p_cas56_truncated.out")
SF6P_MLD = os.path.join(D, "sf6p_cas56_truncated.molden")


def test_older_printout_doublet_and_truncation_warning():
    """SF6+ doublet CASSCF(5,6) in the older 'CI vector ... TOTAL ENERGIES' format, printed with the
    default threshold: 9 of 300 determinants.  The reader warns, and the printed coefficients are
    exactly those of the CASCI eigenvector in the Molden orbitals (overlap = printed weight),
    which checks the determinant phases for an open shell."""
    from pyscf import mcscf, scf
    from pyscf.fci import cistring
    from pyscf.tools import molden as pm
    from pyxscat import ci
    with pytest.warns(UserWarning, match="truncated"):
        info = molpro.read_multi_output(SF6P)
    assert (info["nclosed"], info["nactive"], info["nactive_electrons"], info["ms2"]) == (32, 6, 5, 1)
    s = info["states"][1]
    assert len(s["dets"]) == 9 and s["printed_weight"] == pytest.approx(0.9988, abs=1e-4)
    mol, _, moc, _, _, _ = pm.load(SF6P_MLD)
    mol.charge, mol.spin, mol.verbose = 1, 1, 0
    mol.build()
    mc = mcscf.CASCI(scf.ROHF(mol), 6, (3, 2))
    mc.verbose = 0
    mc.kernel(moc)
    assert mc.e_tot == pytest.approx(s["energy"], abs=1e-6)
    c, (sa, ia), (sb, ib) = ci.ci_matrix(6, 3, 2, s["dets"], s["coeffs"])
    cp = c[np.ix_([ia[int(x)] for x in cistring.make_strings(range(6), 3)],
                  [ib[int(x)] for x in cistring.make_strings(range(6), 2)])]
    assert abs(np.einsum("ij,ij", cp, mc.ci)) == pytest.approx(s["printed_weight"], abs=2e-7)


def test_consistency_check(info):
    from pyxscat.interfaces.molden import read_molden
    assert molpro.check_molden_consistency(info, read_molden(MLD)) < 1e-6


def _energy_in_our_orbitals(mld, info, state):
    """Energy from our RDMs with exact integrals, in the orbitals as OUR Molden reader returns them
    (PySCF provides only the basis function space, matched by function values)."""
    from pyscf import ao2mo
    from pyscf.tools import molden as pm
    from pyxscat.interfaces.molden import read_molden
    data = read_molden(mld)
    mol = pm.load(mld)[0]
    mol.cart, mol.verbose = True, 0
    mol.build()
    pts = np.random.default_rng(0).normal(scale=1.5, size=(3000, 3)) + data.coords[0]
    n = info["nclosed"] + info["nactive"]
    c = np.linalg.lstsq(mol.eval_gto("GTOval_cart", pts), data.basis.eval(pts) @ data.mo_coeff[:, :n], rcond=None)[0]
    h = c.T @ (mol.intor("int1e_kin") + mol.intor("int1e_nuc")) @ c
    eri = ao2mo.restore(1, ao2mo.kernel(mol, c), n)
    d1, d2 = molpro.multi_rdms(info, state)
    return mol.energy_nuc() + np.einsum("pq,pq", h, d1) + 0.5 * np.einsum("pqrs,pqrs", eri, d2)


def test_punch_file_restores_full_vector():
    """The SF6+ punch file holds 299 of 300 determinants (weight 0.9999999) where the output printed
    9; the state energy is then reproduced to 1e-7 Hartree."""
    info = molpro.read_punch(os.path.join(D, "sf6p_cas56.pun"), molden=SF6P_MLD)
    assert (info["nclosed"], info["nactive"], info["ms2"]) == (32, 6, 1)
    s = info["states"][1]
    assert len(s["dets"]) == 299 and s["printed_weight"] > 0.9999
    assert _energy_in_our_orbitals(SF6P_MLD, info, 1) == pytest.approx(s["energy"], abs=1e-7)


def test_molpro_orbital_order_restored():
    """He SA-CASSCF(2,10): Molpro lists the natural orbitals out of its internal order in the
    Molden file (Sym= 1.1 2.1 6.1 3.1 ...).  Read in label order, both state energies match Molpro;
    read in file order they are off by 30-40 mEh."""
    out, pun, mld = (os.path.join(D, f"he_sa2_casscf.{e}") for e in ("out", "pun", "mld"))
    info = molpro.read_punch(pun, output=out)
    for k in (1, 2):
        assert _energy_in_our_orbitals(mld, info, k) == pytest.approx(info["states"][k]["energy"], abs=1e-7)
    from pyxscat.interfaces.molden import read_molden
    assert molpro.check_molden_consistency(info, read_molden(mld)) < 1e-3
    assert molpro.check_molden_consistency(info, read_molden(mld, reorder=False)) > 1e-3


# ----------------------------------------------------------------------------- punch file without output
def test_punch_and_molden_without_output(tmp_path):
    """Only the punch and Molden files: the same wavefunction as with the output; a missing
    output falls back to the punch file; a Molden file with punch= is the MULTI calculation."""
    import shutil
    import warnings
    import thistle as th
    for ext in ("pun", "mld"):
        shutil.copy(os.path.join(D, f"he_sa2_casscf.{ext}"), tmp_path)
    ref = th.load(os.path.join(D, "he_sa2_casscf.out"), state=2).wf
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = th.load(str(tmp_path / "he_sa2_casscf.pun"), state=2).wf
        c = th.load(str(tmp_path / "he_sa2_casscf.mld"), punch=str(tmp_path / "he_sa2_casscf.pun"), state=2).wf
    with pytest.warns(UserWarning, match="not found"):
        b = th.load(str(tmp_path / "he_sa2_casscf.out"), state=2).wf
    for w in (a, b, c):
        assert np.abs(w.dm1 - ref.dm1).max() < 1e-12
        assert np.abs(w.dm2 - ref.dm2).max() < 1e-12
    with pytest.warns(UserWarning, match="elastic scattering follows"):
        m = th.load(str(tmp_path / "he_sa2_casscf.mld"))
    assert m.wf.dm2 is None
    with pytest.raises(FileNotFoundError):
        th.load(str(tmp_path / "nothing.out"))
