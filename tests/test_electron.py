"""Electron scattering: the electron-nucleus integrals, the assembled cross sections and the IAM."""
import os

import numpy as np
import pytest

import pyxscat as xs
from pyxscat import aligned as al
from pyxscat import electron as el
from pyxscat import iam

DATA = os.path.join(os.path.dirname(__file__), "data")
WATER = "O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692"


def _sphere(order=131):
    from scipy.integrate import lebedev_rule
    xyz, w = lebedev_rule(order)
    return xyz.T, w / w.sum()


def _scf(atom, basis, cart=False, unit="angstrom", conv_tol=1e-11):
    from pyscf import gto, scf
    mol = gto.M(atom=atom, basis=basis, cart=cart, unit=unit, verbose=0)
    return scf.RHF(mol).run(conv_tol=conv_tol)


@pytest.fixture(scope="module")
def water():
    return xs.interfaces.pyscf.from_scf(_scf(WATER, "cc-pvdz", cart=True))


@pytest.fixture(scope="module")
def water_ccsd():
    from pyscf import cc
    mf = _scf(WATER, "6-31g")
    return xs.interfaces.pyscf.from_cc(cc.CCSD(mf).run(conv_tol=1e-10))


def _lebedev_ne(wf, q, order=131):
    """<Re[Z(q)^* F(q)]> and its P2 projection by quadrature over directions of q."""
    n, w = _sphere(order)
    p2 = 1.5 * n[:, 2] ** 2 - 0.5
    r0, r2 = [], []
    for qq in q:
        qv = qq * n
        zq = np.exp(1j * qv @ wf.coords.T) @ wf.atom_numbers
        v = (np.conj(zq) * al.form_factor(wf.basis, qv, wf.dm1)).real
        r0.append(w @ v)
        r2.append(w @ (v * p2))            # S2 = <I P2(n_z)>, i.e. -j2(qh) P2(cos theta_h) (i^2 = -1)
    return np.array(r0), np.array(r2)


@pytest.mark.parametrize("basis,cart", [("cc-pvdz", True), ("cc-pvqz", False)])
def test_electron_nucleus_against_quadrature(basis, cart):
    # CO off the axes, up to g functions; q spans the series, Miller and upward Bessel regimes
    wf = xs.interfaces.pyscf.from_scf(_scf("C 0 0.1 -0.2; O 0.3 0 1.9", basis, cart=cart, unit="bohr"))
    q = np.array([0.0, 0.05, 1.9, 7.5, 12.0])
    r0, r2 = _lebedev_ne(wf, q)
    a0 = el.electron_nucleus(wf.basis, q, wf.dm1, wf.atom_numbers, wf.coords)
    a2 = el.electron_nucleus(wf.basis, q, wf.dm1, wf.atom_numbers, wf.coords, l=2)
    np.testing.assert_allclose(a0, r0, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(a2, r2, rtol=1e-9, atol=1e-10)


def test_electron_nucleus_at_q0_is_ZN(water):
    ne = el.electron_nucleus(water.basis, [0.0], water.dm1, water.atom_numbers, water.coords)
    assert abs(ne[0] - 10 * 10) < 1e-8


def test_against_old_code(water):
    ref = np.loadtxt(os.path.join(DATA, "legacy_electron_nucleus_h2o_ccpvdz.txt"))
    new = el.electron_nucleus(water.basis, ref[:, 0], water.dm1, water.atom_numbers, water.coords)
    np.testing.assert_allclose(new, ref[:, 1], rtol=1e-7)     # single-precision pi in the old code


def test_neutral_molecule_and_dipole(water):
    """Z(q) - F(q) = -i q.mu + O(q^2): the elastic bracket starts as q^2 mu^2 / 3."""
    from pyscf import gto, scf
    mol = gto.M(atom=WATER, basis="cc-pvdz", cart=True, verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-11)
    mu = np.linalg.norm(mf.dip_moment(unit="au", verbose=0))
    q = np.array([0.0, 2e-3, 4e-3])
    d = el.decompose(water, q, hf=True)
    assert abs(d["elastic"][0]) < 1e-8 and abs(d["total"][0]) < 1e-8
    np.testing.assert_allclose(d["elastic"][1:] / q[1:] ** 2, mu ** 2 / 3, rtol=2e-3)


def test_single_atom_mott_bethe():
    """One spherical atom: elastic = (Z - f)^2 with f^2 the X-ray elastic scattering."""
    wf = xs.interfaces.pyscf.from_scf(_scf("Ne 0.3 -0.2 0.5", "aug-cc-pvtz"))
    q = np.linspace(0.0, 8.0, 9)
    f = np.sqrt(xs.elastic(wf.basis, q, wf.dm1))
    d = el.decompose(wf, q, hf=True)
    np.testing.assert_allclose(d["elastic"], (10 - f) ** 2, atol=1e-8)
    np.testing.assert_allclose(d["ne"], 10 * f, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(d["inelastic"], xs.scattering(wf, q, "hf_total") - f ** 2, atol=1e-10)


def test_hydrogen_atom_exact():
    """Elastic electron scattering of H(1s): [1 - 16/(4 + q^2)^2]^2 (Mott-Bethe, exact)."""
    from pyscf import gto, scf
    mol = gto.M(atom="H 0 0 0", basis="aug-cc-pv5z", spin=1, verbose=0)
    mf = scf.UHF(mol).run(conv_tol=1e-12)
    da, db = mf.make_rdm1()
    wf = xs.interfaces.pyscf.from_density(mol, da + db)
    q = np.array([0.0, 0.5, 1.0, 2.0, 4.0])
    exact = (1 - 16 / (4 + q ** 2) ** 2) ** 2
    got = el.scattering(wf, q, "elastic")
    np.testing.assert_allclose(got, exact, atol=2e-4)


def test_aligned_average_is_rotational(water_ccsd):
    q = np.array([0.4, 1.3, 3.0])
    n, w = _sphere(59)
    ref = el.decompose(water_ccsd, q)
    ref2 = el.decompose(water_ccsd, q, l=2)
    for i, qq in enumerate(q):
        a = el.aligned(water_ccsd, qq * n)
        for k in ("elastic", "inelastic", "total"):
            assert abs(w @ a[k] - ref[k][i]) < 1e-9 * max(1, abs(ref[k][i]))
        p2 = 1.5 * n[:, 2] ** 2 - 0.5
        assert abs((w @ (a["total"] * p2)) - ref2["total"][i]) < 1e-8      # quadrature error


def test_parts_and_kinds(water_ccsd):
    q = np.linspace(0.0, 6.0, 7)
    d = el.decompose(water_ccsd, q)
    np.testing.assert_allclose(d["total"], d["nn"] - 2 * d["ne"] + d["ee_total"], rtol=1e-12)
    np.testing.assert_allclose(d["inelastic"], xs.scattering(water_ccsd, q, "inelastic"), rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(el.scattering(water_ccsd, q, "elastic"), d["elastic"], rtol=1e-10, atol=1e-9)
    orbitals_only = xs.Wavefunction(water_ccsd.basis, water_ccsd.dm1, nelec=10,
                                    atom_numbers=water_ccsd.atom_numbers, coords=water_ccsd.coords)
    with pytest.raises(ValueError, match="2RDM"):
        el.scattering(orbitals_only, q, "total")
    cs = el.cross_section(d["total"], q)
    assert np.isnan(cs[0]) and np.allclose(cs[1:], 4 * d["total"][1:] / q[1:] ** 4)


def test_iam_electron(water):
    q = np.linspace(0.0, 8.0, 9)
    r = iam.electron_for_wavefunction(water, q, basis="cc-pvdz")
    x = iam.for_wavefunction(water, q, basis="cc-pvdz")
    z = np.array([8.0, 1.0, 1.0])
    # Mott-Bethe from the X-ray atomic form factors
    f = {zz: iam.atomic_form_factor(int(zz), q, basis="cc-pvdz") for zz in set(z)}
    amp = np.array([zz - f[zz] for zz in z])
    rr = np.linalg.norm(water.coords[:, None] - water.coords[None], axis=2)
    qr = q[None, None] * rr[..., None]
    j0 = np.where(qr > 0, np.sin(qr) / np.where(qr > 0, qr, 1), 1)
    np.testing.assert_allclose(r["elastic"], np.einsum("aq,bq,abq->q", amp, amp, j0), atol=1e-10)
    np.testing.assert_allclose(r["inelastic"], x["inelastic"], atol=1e-10)
    np.testing.assert_allclose(r["atomic"] + r["molecular"], r["total"], atol=1e-10)
    assert abs(r["total"][0]) < 1e-6
    # bonding moves the molecule away from the free atoms, but not by a lot at large q
    mol = el.decompose(water, q, hf=True)
    assert abs(mol["total"][-1] / r["total"][-1] - 1) < 0.02


def test_iam_electron_l2_tables():
    """l = 2 IAM: -sum_{A != B} (Z_A - f_A)(Z_B - f_B) j2(q R) P2(cos theta), with tabulated f."""
    import xraydb
    z, c = [8, 6], np.array([[0.3, 0.0, -1.1], [0.0, 0.2, 1.0]])
    q = np.array([0.5, 2.0, 5.0])
    r = iam.electron_scattering(z, c, q, source="tables", total=False, l=2)
    s = q / iam.BOHR_IN_ANGSTROM / (4 * np.pi)
    amp = [8 - np.asarray(xraydb.f0("O", s)), 6 - np.asarray(xraydb.f0("C", s))]
    d = c[0] - c[1]
    h = np.linalg.norm(d)
    p2 = 1.5 * (d[2] / h) ** 2 - 0.5
    np.testing.assert_allclose(r["elastic"], -2 * amp[0] * amp[1] * p2 * iam._j2(q * h), rtol=1e-12)
    # point nuclei: the same with f = 0
    np.testing.assert_allclose(el.nucleus_nucleus(z, c, q, l=2), -2 * 48 * p2 * iam._j2(q * h), rtol=1e-12)
