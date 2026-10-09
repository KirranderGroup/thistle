"""Fixed-orientation (aligned) scattering."""
import numpy as np
import pytest

import pyxscat as xs
from pyxscat import aligned as al
from pyxscat.interfaces import pyscf as xp
from conftest import get_case, H2O, CO


def _rand_q(n, seed=0, scale=2.0):
    q = np.random.default_rng(seed).normal(scale=scale, size=(n, 3))
    q[0] = 0.0
    return q


@pytest.mark.parametrize("atom,basis", [(CO, "6-31g*"), (H2O, "cc-pvtz"), ("Ne 0.1 0.2 0.3", "cc-pvqz")])
def test_pair_transform_matches_pyscf(atom, basis):
    """PySCF's ft_aopair uses exp(-i q.r); ours is its complex conjugate (exp(+i q.r)). Up to g."""
    from pyscf import gto
    from pyscf.gto import ft_ao
    mol = gto.M(atom=atom, basis=basis, cart=True, verbose=0)
    q = _rand_q(9)
    f = al.pair_transform(xp.basis_from_mole(mol), q)
    ref = ft_ao.ft_aopair(mol, q).reshape(len(q), mol.nao, mol.nao)
    assert np.abs(f - ref.conj()).max() < 1e-12 * max(1.0, np.abs(ref).max())


def test_spherical_basis_form_factor():
    from pyscf.gto import ft_ao
    case = get_case(H2O, "cc-pvdz", False)
    wf = xp.from_scf(case.mf)
    q = _rand_q(6, seed=3)
    ref = np.einsum("qij,ij->q", ft_ao.ft_aopair(case.mol, q).reshape(len(q), case.mol.nao, case.mol.nao),
                    case.mf.make_rdm1()).conj()
    assert np.allclose(al.form_factor(wf.basis, q, wf.dm1), ref, atol=1e-11)


def test_zero_q_limits(h2o_631g):
    wf = xp.from_cc(h2o_631g.cc)
    n = wf.nelec
    d = al.decompose(wf, np.zeros((1, 3)))
    assert d["elastic"][0] == pytest.approx(n * n, rel=1e-10)
    assert d["two_electron"][0] == pytest.approx(n * (n - 1), rel=1e-10)
    assert d["total"][0] == pytest.approx(n * n, rel=1e-10)


@pytest.mark.parametrize("kind", ["elastic", "exchange", "two_electron"])
def test_orientational_average_equals_j0_kernel(h2o_631g, kind):
    """Averaging the fixed-orientation result over a Lebedev grid reproduces the rotationally
    averaged kernel (two independent code paths)."""
    from scipy.integrate import lebedev_rule
    wf = xp.from_cc(h2o_631g.cc)
    xyz, w = lebedev_rule(131)
    w = w / w.sum()
    qs = np.array([0.5, 1.5, 3.0, 6.0])
    avg = []
    for q in qs:
        qv = q * xyz.T
        if kind == "elastic":
            v = al.elastic(wf.basis, qv, wf.dm1)
        elif kind == "exchange":
            v = al.exchange(wf.basis, qv, wf.dm1)
        else:
            v = al.two_electron(wf.basis, qv, wf.dm2)
        avg.append(w @ v)
    fn = {"elastic": xs.elastic, "exchange": xs.exchange, "two_electron": xs.two_electron}[kind]
    arg = wf.dm2 if kind == "two_electron" else wf.dm1
    ref = fn(wf.basis, qs, arg, screen_tol=0.0, cutoffcentre=1e-6)
    assert np.allclose(avg, ref, rtol=1e-8, atol=1e-8)


def test_mean_field_equals_explicit_hf_dm2(h2o_631g):
    wf = xp.from_scf(h2o_631g.mf, with_dm2=True)
    q = _rand_q(8, seed=5)
    assert np.allclose(al.two_electron(wf.basis, q, wf.dm2), al.mean_field(wf.basis, q, wf.dm1), atol=1e-10)


def test_decomposition_and_scattering_agree(h2o_631g):
    wf = xp.from_cc(h2o_631g.cc)
    q = _rand_q(8, seed=7)
    d = al.decompose(wf, q)
    assert np.allclose(d["self"] + d["hartree"] + d["exchange"] + d["cumulant"], d["total"], atol=1e-10)
    assert np.allclose(d["total"], al.scattering(wf, q, "total"), atol=1e-10)
    assert np.allclose(d["inelastic"], al.scattering(wf, q, "inelastic"), atol=1e-10)
    assert np.allclose(d["exchange"], al.exchange(wf.basis, q, wf.dm1), atol=1e-12)


def test_inversion_symmetry_of_intensities(h2o_631g):
    """Real wavefunctions: I(q) = I(-q) (Friedel's law), even for a non-centrosymmetric molecule."""
    wf = xp.from_cc(h2o_631g.cc)
    q = _rand_q(6, seed=9)
    assert np.allclose(al.scattering(wf, q, "total"), al.scattering(wf, -q, "total"), atol=1e-11)
    assert np.allclose(al.elastic(wf.basis, q, wf.dm1), al.elastic(wf.basis, -q, wf.dm1), atol=1e-11)


def test_rotating_q_equals_rotating_the_molecule():
    from pyscf import gto, scf
    r = al.rotation_matrix(0.3, 1.1, -0.7)
    m1 = gto.M(atom=CO, basis="6-31g", verbose=0)
    m1.build()
    coords = m1.atom_coords() @ r.T
    m2 = gto.M(atom=[(m1.atom_symbol(i), coords[i]) for i in range(m1.natm)], unit="Bohr", basis="6-31g",
               verbose=0)
    f1 = scf.RHF(m1).run()
    f2 = scf.RHF(m2).run()
    q_lab = _rand_q(7, seed=11)
    lab = al.elastic(xp.basis_from_mole(m2), q_lab, xp.from_scf(f2).dm1)
    mol_frame = al.elastic(xp.basis_from_mole(m1), al.rotate_q(q_lab, r), xp.from_scf(f1).dm1)
    assert np.allclose(lab, mol_frame, rtol=1e-7, atol=1e-7)


def test_ewald_geometry():
    k0 = 2 * np.pi / 1.5
    tt = np.linspace(0, np.pi, 7)
    q = al.ewald_q(k0, tt, 0.3)
    assert np.allclose(np.linalg.norm(q, axis=1), 2 * k0 * np.sin(tt / 2))
    qd, ok = al.detector_q(k0, 0.8 * k0, 5)
    assert qd.shape == (25, 3)
    assert ok.sum() == 21                       # the four corners lie beyond |q_perp| = k0
    assert np.allclose(np.linalg.norm(qd[ok] + np.array([0, 0, k0]), axis=1), k0)   # on the Ewald sphere
