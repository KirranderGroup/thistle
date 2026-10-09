"""Physical limits and comparison with an independent evaluation."""
import numpy as np
import pytest

import pyxscat as xs
from pyxscat.interfaces import pyscf as xp
from conftest import get_case, rotational_average_reference, H2O, CO

Q = np.array([0.0, 0.3, 1.0, 2.5, 5.0, 8.0])


def test_zero_q_limits(h2o_631g):
    wf = xp.from_cc(h2o_631g.cc)
    n = wf.nelec
    d = xs.decompose(wf, [0.0])
    assert d["two_electron"][0] == pytest.approx(n * (n - 1), rel=1e-9)
    assert d["elastic"][0] == pytest.approx(n * n, rel=1e-9)
    assert d["total"][0] == pytest.approx(n * n, rel=1e-9)
    assert d["inelastic"][0] == pytest.approx(0.0, abs=1e-7)
    # exchange at q = 0 is -1/2 Tr(D S D S)
    assert d["exchange"][0] == pytest.approx(-0.5 * np.trace(_dsds(h2o_631g)), rel=1e-9)


def _dsds(case):
    """D S D S in the AO basis (trace = 2N for an idempotent closed-shell density)."""
    s = case.mol.intor("int1e_ovlp")
    d = case.cc.mo_coeff @ case.cc.make_rdm1() @ case.cc.mo_coeff.T
    return d @ s @ d @ s


def test_large_q_limit(h2o_631g):
    """Total scattering tends to N: the two-electron term dies off."""
    wf = xp.from_cc(h2o_631g.cc)
    i2 = xs.two_electron(wf.basis, [40.0, 80.0], wf.dm2)
    assert abs(i2[0]) < 5e-3 and abs(i2[1]) < 5e-4


@pytest.mark.parametrize("atom,basis,cart", [
    (H2O, "6-31g", True),
    (CO, "6-31g*", True),
    (H2O, "cc-pvdz", False),      # spherical d functions, general contraction
])
def test_elastic_against_independent_reference(atom, basis, cart):
    case = get_case(atom, basis, cart)
    wf = xp.from_scf(case.mf)
    ref = rotational_average_reference(case.mol, Q, case.mf.make_rdm1())
    got = xs.elastic(wf.basis, Q, wf.dm1, screen_tol=0.0, cutoffcentre=1e-6)
    assert np.allclose(got, ref, rtol=1e-8, atol=1e-8)


@pytest.mark.parametrize("atom,basis,cart", [
    (H2O, "6-31g", True),
    (CO, "6-31g*", True),
    (H2O, "cc-pvdz", False),
])
def test_ccsd_two_electron_against_independent_reference(atom, basis, cart):
    case = get_case(atom, basis, cart)
    wf = xp.from_cc(case.cc)
    c = case.cc.mo_coeff
    dm2_ao = np.einsum("pqrs,ap,bq,cr,ds->abcd", case.cc.make_rdm2(), c, c, c, c, optimize=True)
    ref = rotational_average_reference(case.mol, Q, None, dm2=dm2_ao)
    got = xs.two_electron(wf.basis, Q, wf.dm2, screen_tol=0.0, cutoffcentre=1e-6)
    assert np.allclose(got, ref, rtol=1e-8, atol=1e-8)


def test_f_functions():
    """f functions (cc-pVTZ on a small system) against the independent reference."""
    case = get_case("Li 0 0 0; H 0.2 0.1 1.6", "cc-pvtz", False)
    wf = xp.from_scf(case.mf)
    q = np.array([0.0, 1.0, 4.0])
    ref = rotational_average_reference(case.mol, q, case.mf.make_rdm1())
    got = xs.elastic(wf.basis, q, wf.dm1, screen_tol=0.0, cutoffcentre=1e-6)
    assert np.allclose(got, ref, rtol=1e-8, atol=1e-8)


def test_transition_like_cross_term(h2o_631g):
    """Two different one-particle matrices: rotational average of F_a F_b^*."""
    case = h2o_631g
    da = case.mf.make_rdm1()
    c = case.mf.mo_coeff
    db = np.outer(c[:, 3], c[:, 5]) + np.outer(c[:, 5], c[:, 3])       # a transition-density-like matrix
    wf = xp.from_density(case.mol, da)
    ref = rotational_average_reference(case.mol, Q, da, db)
    got = xs.elastic(wf.basis, Q, da, db, screen_tol=0.0, cutoffcentre=1e-6)
    assert np.allclose(got, ref, rtol=1e-8, atol=1e-8)


def test_near_coincident_pair_centres():
    """Known limitation (inherited): pair centres closer than ~1e-3 bohr are handled either by the
    exact h = 0 limit or by the general expansion, which loses digits as h -> 0.  With the default
    switch-over (cutoffcentre = 1e-4) the relative error stays below 1e-4; see docs/STATUS.md."""
    from pyscf import gto
    rng = np.random.default_rng(0)
    q = np.array([0.0, 0.5, 2.0, 6.0, 12.0])
    for d in (3e-2, 3e-3, 3e-4, 3e-5):
        mol = gto.M(atom=f"O 0 0 0; O {d} {d/2} {d/3}; H 0.3 0.9 1.2", basis="cc-pvdz", cart=True,
                    spin=1, verbose=0)
        dm = rng.normal(size=(mol.nao, mol.nao)); dm = dm + dm.T
        ref = rotational_average_reference(mol, q, dm)
        got = xs.elastic(xp.basis_from_mole(mol), q, dm, screen_tol=0.0)
        assert (np.abs(got - ref) / np.abs(ref).max()).max() < 1e-4
