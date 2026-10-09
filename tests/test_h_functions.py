"""h functions (l = 5): every kernel against the fixed-orientation route averaged on a Lebedev grid."""
import numpy as np
import pytest

import pyxscat as xs
from pyxscat import aligned as al
from pyxscat import electron as el
from pyxscat.interfaces import pyscf as xp


@pytest.fixture(scope="module")
def co_h():
    from pyscf import gto, scf
    c = gto.basis.load("6-31g", "C") + [[5, [0.9, 1.0]]]           # one h shell on carbon
    mol = gto.M(atom="C 0 0.1 -0.2; O 0.3 0 1.9", unit="bohr", basis={"C": c, "O": "sto-3g"}, verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-10)
    return xp.from_scf(mf)


def _sphere():
    from scipy.integrate import lebedev_rule
    xyz, w = lebedev_rule(131)
    return xyz.T, w / w.sum()


def test_basis_has_h(co_h):
    assert co_h.basis.lmn.sum(axis=1).max() == 5


@pytest.mark.parametrize("l", [0, 2])
def test_elastic_and_exchange(co_h, l):
    n, w = _sphere()
    wt = np.ones(len(w)) if l == 0 else (1.5 * n[:, 2] ** 2 - 0.5)
    q = np.array([0.05, 0.6, 1.5])
    ref_el = np.array([w @ (np.abs(al.form_factor(co_h.basis, qq * n, co_h.dm1)) ** 2 * wt) for qq in q])
    ref_ex = np.array([w @ (al.exchange(co_h.basis, qq * n, co_h.dm1) * wt) for qq in q])
    np.testing.assert_allclose(xs.elastic(co_h.basis, q, co_h.dm1, l=l, screen_tol=0), ref_el, rtol=1e-11)
    np.testing.assert_allclose(xs.elastic(co_h.basis, q, co_h.dm1, l=l, screen_tol=0, method="quartets"),
                               ref_el, rtol=1e-8)
    np.testing.assert_allclose(xs.exchange(co_h.basis, q, co_h.dm1, l=l, screen_tol=0), ref_ex, rtol=1e-8)


def test_electron_nucleus(co_h):
    n, w = _sphere()
    q = np.array([0.0, 0.6, 2.0])
    ref = []
    for qq in q:
        zq = np.exp(1j * qq * n @ co_h.coords.T) @ co_h.atom_numbers
        ref.append(w @ (np.conj(zq) * al.form_factor(co_h.basis, qq * n, co_h.dm1)).real)
    got = el.electron_nucleus(co_h.basis, q, co_h.dm1, co_h.atom_numbers, co_h.coords)
    np.testing.assert_allclose(got, ref, rtol=1e-11)
