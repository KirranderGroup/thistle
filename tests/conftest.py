import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

pyscf = pytest.importorskip("pyscf")

# Tests on large primitive basis sets (ANO) take minutes; run them with PYXSCAT_SLOW=1.
slow = pytest.mark.skipif(not os.environ.get("PYXSCAT_SLOW"), reason="slow; set PYXSCAT_SLOW=1")
from pyscf import gto, scf, cc  # noqa: E402

H2O = "O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692"
CO = "C 0.1 -0.2 0; O 0.3 0.4 1.1"          # deliberately off-axis


class Case:
    def __init__(self, atom, basis, cart):
        self.mol = gto.M(atom=atom, basis=basis, cart=cart, verbose=0)
        self.mf = scf.RHF(self.mol).run(conv_tol=1e-12)
        self._cc = None

    @property
    def cc(self):
        if self._cc is None:
            self._cc = cc.CCSD(self.mf).run(conv_tol=1e-10)
        return self._cc


_cache = {}


def get_case(atom, basis, cart):
    key = (atom, basis, cart)
    if key not in _cache:
        _cache[key] = Case(atom, basis, cart)
    return _cache[key]


@pytest.fixture(scope="session")
def h2o_631g():
    return get_case(H2O, "6-31g", True)


@pytest.fixture(scope="session")
def h2o_dz_sph():
    return get_case(H2O, "cc-pvdz", False)


def rotational_average_reference(mol, q, dm1_a, dm1_b=None, dm2=None, order=131):
    """Independent reference: Fourier transforms of AO pairs from PySCF (analytic), averaged over
    orientations with a Lebedev quadrature.  dm1/dm2 in the AO basis of `mol`.

    Returns <F_a F_b^*> if dm2 is None, else <sum dm2 f f^*>.
    """
    from pyscf.gto import ft_ao
    from scipy.integrate import lebedev_rule

    xyz, w = lebedev_rule(order)
    w = w / w.sum()
    out = []
    n = mol.nao
    for qq in np.atleast_1d(q):
        f = ft_ao.ft_aopair(mol, qq * xyz.T).reshape(-1, n * n)      # (npts, n^2)
        if dm2 is None:
            fa = f @ dm1_a.reshape(-1)
            fb = fa if dm1_b is None else f @ dm1_b.reshape(-1)
            out.append(np.sum(w * (fa * fb.conj()).real))
        else:
            g = dm2.reshape(n * n, n * n)
            out.append(np.sum(w * np.einsum("pi,pi->p", f @ g, f.conj()).real))
    return np.array(out)
