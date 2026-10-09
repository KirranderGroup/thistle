"""Second Legendre moment (j2) of the scattering about the molecular z axis.

S2(q) = -sum Gamma <j2(q r12) P2(cos theta_12)> (the i^2 of the plane-wave expansion included);
for molecules aligned along z, I(q, theta_q) = I0(q) + 5 <P2> P2(cos theta_q) S2(q) + ...
"""
import os

import numpy as np
import pytest

import pyxscat as xs
from pyxscat import aligned as al, iam
from pyxscat.interfaces import pyscf as xp
from conftest import get_case

WATER = "O 0.05 0.02 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692"      # no symmetry plane along z
Q = np.array([0.5, 1.0, 2.0, 4.0])


def _projection(fn, qs, order=131):
    """S2 from fixed-orientation results: S2(q) = <P2(cos theta_q) I(q n)>_n (Lebedev)."""
    from scipy.integrate import lebedev_rule
    xyz, w = lebedev_rule(order)
    w = w / w.sum()
    p2 = 1.5 * xyz[2] ** 2 - 0.5
    return np.array([(w * p2) @ fn(q * xyz.T) for q in qs])


@pytest.fixture(scope="module")
def water():
    case = get_case(WATER, "6-31g*", False)
    return xp.from_cc(case.cc)


@pytest.mark.parametrize("kind", ["elastic", "exchange", "two_electron"])
def test_j2_kernel_equals_projection_of_aligned_results(water, kind):
    """Two independent routes: the rotational kernel (j0 derivatives, Hobson's theorem) and the
    P2 projection of the fixed-orientation plane-wave results."""
    b = water.basis
    opts = dict(screen_tol=0.0, cutoffcentre=1e-6)
    if kind == "elastic":
        got = xs.elastic(b, Q, water.dm1, l=2, **opts)
        ref = _projection(lambda qv: al.elastic(b, qv, water.dm1), Q)
    elif kind == "exchange":
        got = xs.exchange(b, Q, water.dm1, l=2, **opts)
        ref = _projection(lambda qv: al.exchange(b, qv, water.dm1), Q)
    else:
        got = xs.two_electron(b, Q, water.dm2, l=2, **opts)
        ref = _projection(lambda qv: al.two_electron(b, qv, water.dm2), Q)
    assert np.allclose(got, ref, atol=5e-8)
    assert np.abs(ref).max() > 1e-3


def test_matches_eirik_elastic_j2():
    """The new kernel reproduces Eirik's elastic j2 code, which uses Y20 instead of P2."""
    from pyscf import gto, scf
    data = np.loadtxt(os.path.join(os.path.dirname(__file__), "data", "eirik_elastic_j2.txt"))
    mol = gto.M(atom=WATER, basis="6-31g", cart=True, verbose=0)
    wf = xp.from_scf(scf.RHF(mol).run())          # PySCF default convergence, as for the reference
    # the quartet kernel with Eirik's settings (pair centres closer than 0.01 bohr treated as
    # coincident) reproduces his numbers
    # Eirik's j2 has no i^2 = -1: ours = -his / sqrt(5/(4 pi))
    ours = -xs.elastic(wf.basis, data[:, 0], wf.dm1, l=2, screen_tol=0.0, cutoffcentre=0.01, method="quartets")
    assert np.allclose(np.sqrt(5 / (4 * np.pi)) * ours, data[:, 1], rtol=1e-6, atol=1e-9)
    # the default (pair) route is exact; the 0.01 bohr approximation costs 2.5e-4 here
    exact = -xs.elastic(wf.basis, data[:, 0], wf.dm1, l=2)
    assert np.allclose(np.sqrt(5 / (4 * np.pi)) * exact, data[:, 1], rtol=5e-4)
    assert not np.allclose(np.sqrt(5 / (4 * np.pi)) * exact, data[:, 1], rtol=1e-5)


def test_separated_atoms_equal_iam():
    """Two Ne atoms 14 bohr apart on an oblique axis: the molecule IS the IAM, and its j2 moment is
    the interference term f f j2(qR) P2(cos theta_R) alone."""
    from pyscf import gto, scf
    v = 14.0 * np.array([0.48, 0.36, 0.8])
    mol = gto.M(atom=f"Ne 0 0 0; Ne {v[0]} {v[1]} {v[2]}", unit="Bohr", basis="aug-cc-pvtz", verbose=0)
    wf = xp.from_scf(scf.RHF(mol).run(conv_tol=1e-11))
    q = np.array([0.2, 0.5, 1.0, 2.0, 4.0, 8.0])
    ref = iam.for_wavefunction(wf, q, basis="aug-cc-pvtz", l=2)
    opts = dict(screen_tol=0.0, cutoffcentre=1e-6)
    assert np.allclose(xs.elastic(wf.basis, q, wf.dm1, l=2, **opts), ref["elastic"], rtol=1e-6, atol=1e-8)
    assert np.allclose(xs.scattering(wf, q, "hf_total", l=2, **opts), ref["total"], rtol=1e-6, atol=1e-8)


def test_spherical_atom_has_no_j2():
    from pyscf import gto, scf
    mol = gto.M(atom="Ne 0.3 -0.2 0.1", basis="cc-pvtz", verbose=0)
    wf = xp.from_scf(scf.RHF(mol).run())
    assert np.abs(xs.scattering(wf, Q, "hf_total", l=2)).max() < 1e-9


def test_iam_j2_equals_projection_of_aligned_iam(h2o_631g):
    wf = xp.from_scf(h2o_631g.mf)
    q = np.array([0.5, 1.5, 3.0, 6.0])
    r = iam.for_wavefunction(wf, q, basis="6-31g", l=2)
    ref = _projection(lambda qv: iam.for_wavefunction(wf, qv, aligned=True, basis="6-31g")["total"], q)
    assert np.allclose(r["total"], ref, atol=1e-10)
    assert np.allclose(r["inelastic"], 0.0)


def test_small_q_and_decomposition(water):
    d = xs.decompose(water, [0.01, 1.0, 3.0], l=2)
    assert abs(d["total"][0]) < 1e-4                                   # S2 ~ q^2 at small q
    assert np.allclose(d["self"], 0.0)
    assert np.allclose(d["hartree"] + d["exchange"] + d["cumulant"], d["total"], atol=1e-12)
    assert np.allclose(d["total"], xs.scattering(water, [0.01, 1.0, 3.0], "total", l=2), atol=1e-12)


def test_bad_order_is_refused(h2o_631g):
    wf = xp.from_scf(h2o_631g.mf)
    with pytest.raises(RuntimeError):
        xs.elastic(wf.basis, Q, wf.dm1, l=4)


def test_small_q_quartet_route():
    """The j2 moment from the quartet kernel at very small q, where [3 d2/dz2 + q^2] j0 / (2 q^2)
    cancels: exchange and CCSD two-electron terms against the Lebedev projection of the
    fixed-orientation terms, elastic against the exact pair route.  (Before the radial functions
    were evaluated stably, the quartet route was off by 36 % at q = 0.01 bohr^-1.)"""
    from pyscf import gto, scf, cc
    from scipy.integrate import lebedev_rule
    mol = gto.M(atom=WATER, basis="6-31g", verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-11)
    wf = xp.from_cc(cc.CCSD(mf).run(conv_tol=1e-10))
    xyz, w = lebedev_rule(59)
    n, w = xyz.T, w / w.sum()
    p2 = 1.5 * n[:, 2] ** 2 - 0.5
    q = np.array([1e-3, 1e-2, 0.1, 1.0])
    ex = xs.exchange(wf.basis, q, wf.dm1, l=2, screen_tol=0.0)
    i2 = xs.two_electron(wf.basis, q, wf.dm2, l=2, screen_tol=0.0)
    el_q = xs.elastic(wf.basis, q, wf.dm1, l=2, screen_tol=0.0, method="quartets")
    el_p = xs.elastic(wf.basis, q, wf.dm1, l=2, screen_tol=0.0)
    ref_ex = np.array([(w @ (al.exchange(wf.basis, qq * n, wf.dm1) * p2)) for qq in q])
    ref_i2 = np.array([(w @ (al.two_electron(wf.basis, qq * n, wf.dm2) * p2)) for qq in q])
    for got, ref in ((ex, ref_ex), (i2, ref_i2), (el_q, el_p)):
        np.testing.assert_allclose(got[1:], ref[1:], rtol=1e-9)          # q >= 0.01
        np.testing.assert_allclose(got[0], ref[0], rtol=1e-6)            # q = 0.001
