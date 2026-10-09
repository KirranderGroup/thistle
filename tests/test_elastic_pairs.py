"""The pair-factorised elastic kernel (xs_elastic.f90), the default route for elastic scattering."""
import numpy as np
import pytest

import pyxscat as xs
from pyxscat import aligned as al
from pyxscat.interfaces import pyscf as xp

WATER = "O 0.05 0.02 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692"


def _wf(atom, basis, cart=False, unit="angstrom"):
    from pyscf import gto, scf
    mol = gto.M(atom=atom, basis=basis, cart=cart, unit=unit, verbose=0)
    return xp.from_scf(scf.RHF(mol).run(conv_tol=1e-10))


def _lebedev(wf, q, l, da=None, db=None, order=131):
    from scipy.integrate import lebedev_rule
    xyz, w = lebedev_rule(order)
    n, w = xyz.T, w / w.sum()
    da = wf.dm1 if da is None else da
    db = da if db is None else db
    weight = np.ones(len(w)) if l == 0 else (1.5 * n[:, 2] ** 2 - 0.5)    # S2 = <I P2>: includes i^2 = -1
    out = []
    for qq in q:
        fa = al.form_factor(wf.basis, qq * n, da)
        fb = al.form_factor(wf.basis, qq * n, db)
        out.append(w @ ((fa * fb.conj()).real * weight))
    return np.array(out)


@pytest.fixture(scope="module")
def water():
    return _wf(WATER, "6-31g", cart=True)


@pytest.mark.parametrize("l", [0, 2])
def test_exact_against_quadrature(water, l):
    """Includes very small q (where the j2 combination cancels) and pair centres a few hundredths of
    a bohr apart (which the quartet kernel treats as coincident below cutoffcentre)."""
    q = np.array([1e-3, 0.05, 0.3, 1.0, 3.0, 5.0])
    got = xs.elastic(water.basis, q, water.dm1, l=l, screen_tol=0.0)
    ref = _lebedev(water, q, l)
    np.testing.assert_allclose(got, ref, rtol=1e-11, atol=1e-14)


@pytest.mark.parametrize("basis,cart", [("cc-pvdz", True), ("cc-pvtz", False)])
def test_same_as_quartet_kernel(basis, cart):
    wf = _wf("C 0 0.1 -0.2; O 0.3 0 1.9", basis, cart=cart, unit="bohr")
    q = np.array([0.0, 0.4, 1.7, 4.5, 9.0])
    a = xs.elastic(wf.basis, q, wf.dm1, screen_tol=0.0)
    b = xs.elastic(wf.basis, q, wf.dm1, screen_tol=0.0, cutoffcentre=1e-8, method="quartets")
    np.testing.assert_allclose(a, b, rtol=1e-9)
    assert a[0] == pytest.approx(14.0 ** 2, rel=1e-12)          # N^2 exactly at q = 0


def test_transition_densities(water):
    """Mixed terms <F_a F_b^*> with two different (non-symmetric) one-particle matrices."""
    rng = np.random.default_rng(1)
    n = water.basis.nao
    da = water.dm1 + 0.05 * rng.standard_normal((n, n))
    db = 0.3 * rng.standard_normal((n, n))
    q = np.array([0.0, 0.8, 2.5])
    got = xs.elastic(water.basis, q, da, db, screen_tol=0.0)
    quart = xs.elastic(water.basis, q, da, db, screen_tol=0.0, cutoffcentre=1e-8, method="quartets")
    np.testing.assert_allclose(got, quart, rtol=1e-8, atol=1e-10)
    sym = lambda m: 0.5 * (m + m.T)
    np.testing.assert_allclose(got, _lebedev(water, q, 0, sym(da), sym(db)), rtol=1e-10, atol=1e-12)


def test_screening_threads_and_order(water):
    q = np.linspace(0.0, 12.0, 41)
    ref = xs.elastic(water.basis, q, water.dm1, screen_tol=0.0)
    got, st = xs.elastic(water.basis, q, water.dm1, return_stats=True)      # default screening
    assert np.abs(got - ref).max() < 1e-9
    assert 0 < st["pair_pairs_computed"] <= st["pair_pairs"]
    np.testing.assert_allclose(xs.elastic(water.basis, q, water.dm1, nthreads=1), got, rtol=1e-12, atol=1e-13)
    perm = np.random.default_rng(0).permutation(q.size)                   # q in any order
    np.testing.assert_allclose(xs.elastic(water.basis, q[perm], water.dm1), got[perm], rtol=1e-12, atol=1e-13)


def test_tighter_than_quartet_screening():
    """At the same tolerance the pair route's (rigorous, q-dependent) screening is far tighter."""
    wf = _wf(WATER, "aug-cc-pvdz")
    q = np.array([0.0, 1.0, 4.0])
    exact = xs.elastic(wf.basis, q, wf.dm1, screen_tol=0.0)
    pairs = xs.elastic(wf.basis, q, wf.dm1, screen_tol=1e-10)
    assert np.abs(pairs - exact).max() < 1e-8
