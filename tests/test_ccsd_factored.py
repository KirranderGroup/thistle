"""The CCSD two-particle density as a sum of products (pyxscat.ccsd, kernel mode 3): the same
numbers as PySCF's four-index array through every route, without ever forming it."""
import numpy as np
import pytest

import pyxscat as xs
from pyxscat import aligned as al, ccsd, electron as xe
from pyxscat.core import FactoredDM2

pyscf = pytest.importorskip("pyscf")
from pyscf import cc, gto, scf  # noqa: E402

WATER = "O 0.05 0.02 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692"
TIGHT = dict(screen_tol=1e-16)         # the screening error of either route (3e-10 by default) is not the point


@pytest.fixture(scope="module", params=[("cc-pvdz", None), ("6-31g", 1)], ids=["cc-pVDZ", "6-31G/frozen-core"])
def case(request):
    basis, frozen = request.param
    mol = gto.M(atom=WATER, basis=basis, verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-12)
    mycc = cc.CCSD(mf, frozen=frozen).run(conv_tol=1e-10)
    mycc.solve_lambda()
    return mycc


def test_factors_reproduce_make_rdm2_blockwise(case):
    """Every block of the symmetrised MO density agrees to machine precision."""
    g = ccsd.factored_rdm2(case)
    dense = xs.core.symmetrize_dm2(case.make_rdm2())
    gd = g.dense()
    nmo = case.mo_occ.size
    nocc = int(np.count_nonzero(case.mo_occ > 0))
    o, v = slice(0, nocc), slice(nocc, nmo)
    for name, sl in {"oooo": (o, o, o, o), "ooov": (o, o, o, v), "oovv": (o, o, v, v), "ovov": (o, v, o, v),
                     "ovvo": (o, v, v, o), "ovvv": (o, v, v, v), "vvvv": (v, v, v, v)}.items():
        assert np.abs(gd[sl] - dense[sl]).max() < 1e-13, name
    assert np.abs(gd - dense).max() < 1e-13
    assert g.trace() == pytest.approx(case.mol.nelectron * (case.mol.nelectron - 1), abs=1e-10)
    assert g.nterms == (2 * case.t1.shape[0] ** 2 + 2, 5 * case.t1.shape[0] ** 2 + 2)


def test_random_kernel_contraction(case):
    """Contracted with a random 8-fold symmetric tensor (what every kernel is)."""
    g = ccsd.factored_rdm2(case)
    rng = np.random.default_rng(1)
    k = xs.core.symmetrize_dm2(rng.standard_normal(g.shape))
    ref = np.einsum("abcd,abcd->", case.make_rdm2(), k)
    assert abs(ccsd.contract_dense(g, k) - ref) < 1e-12 * abs(ref)


def test_all_routes_agree_with_the_dense_density(case):
    wd = xs.interfaces.pyscf.from_cc(case, factored=False)
    wf = xs.interfaces.pyscf.from_cc(case, factored=True)
    assert isinstance(wf.dm2, FactoredDM2) and not isinstance(wd.dm2, FactoredDM2)
    assert np.abs(wd.dm1 - wf.dm1).max() < 1e-14
    q = np.linspace(0.0, 6.0, 7)
    a, b = xs.decompose(wd, q, **TIGHT), xs.decompose(wf, q, **TIGHT)
    for k in ("total", "cumulant", "inelastic"):
        assert np.abs(a[k] - b[k]).max() < 1e-11 * np.abs(a[k]).max(), k
    assert b["total"][0] == pytest.approx(case.mol.nelectron ** 2, rel=1e-10)
    # j2 moment
    a2 = xs.scattering(wd, q[1:], "total", l=2, **TIGHT)
    b2 = xs.scattering(wf, q[1:], "total", l=2, **TIGHT)
    assert np.abs(a2 - b2).max() < 1e-11 * np.abs(a2).max()
    # fixed orientation
    qv = np.array([[0.3, 0.1, 0.5], [1.5, -0.7, 0.2], [0.0, 0.0, 3.0], [2.0, 2.0, 1.0]])
    da, db = al.decompose(wd, qv), al.decompose(wf, qv)
    assert np.abs(da["total"] - db["total"]).max() < 1e-12 * np.abs(da["total"]).max()
    assert np.abs(da["cumulant"] - db["cumulant"]).max() < 1e-12
    assert np.abs(al.two_electron(wd.basis, qv, wd.dm2) - al.two_electron(wf.basis, qv, wf.dm2)).max() < 1e-11
    # electrons
    ea, eb = xe.decompose(wd, q[1:], **TIGHT), xe.decompose(wf, q[1:], **TIGHT)
    assert np.abs(ea["total"] - eb["total"]).max() < 1e-11 * np.abs(ea["total"]).max()


def test_kernel_mode3_with_random_factors():
    """Kernel mode 3 against mode 1 on arbitrary (non-symmetric) factor matrices, l = 0 and 2,
    with and without the BLAS dgemm."""
    from pyxscat import _lib
    mol = gto.M(atom=WATER, basis="6-31g", verbose=0)
    wf = xs.interfaces.pyscf.from_scf(scf.RHF(mol).run())
    n = wf.basis.nao
    rng = np.random.default_rng(0)
    g = FactoredDM2(rng.standard_normal((5, n, n)), rng.standard_normal((5, n, n)),
                    rng.standard_normal((3, n, n)), rng.standard_normal((3, n, n)))
    q = np.linspace(0.0, 4.0, 5)
    ref = xs.two_electron(wf.basis, q, g.dense(), **TIGHT)
    ref2 = xs.two_electron(wf.basis, q[1:], g.dense(), l=2, **TIGHT)
    saved = _lib.DGEMM
    try:
        for ptr in (saved, None):
            _lib.DGEMM = ptr
            assert np.abs(xs.two_electron(wf.basis, q, g, **TIGHT) - ref).max() < 1e-12 * np.abs(ref).max()
            assert np.abs(xs.two_electron(wf.basis, q[1:], g, l=2, **TIGHT) - ref2).max() < 1e-12 * np.abs(ref).max()
    finally:
        _lib.DGEMM = saved
    # one thread and all threads
    assert np.abs(xs.two_electron(wf.basis, q, g, nthreads=1, **TIGHT) - xs.two_electron(wf.basis, q, g, **TIGHT)).max() \
        < 1e-12 * np.abs(ref).max()
    # the Hartree-Fock pair density as factors
    d = wf.dm1
    hf = FactoredDM2(ac=d[None], bc=d[None], ax=(-0.5 * d)[None], bx=d[None])
    assert np.abs(xs.two_electron(wf.basis, q, hf, **TIGHT) - xs.mean_field(wf.basis, q, d, **TIGHT)).max() < 1e-11


def test_thistle_uses_the_factored_density(case):
    import thistle as th
    from pyxscat import core as xp            # the threshold lives in pyxscat.core
    saved = xp.FACTORED_ABOVE
    try:
        xp.FACTORED_ABOVE = 0
        m = th.load(case)
        assert isinstance(m.wf.dm2, FactoredDM2) and m.correlated
        r = m.scattering(qmax=3, points=4, parts="all")
        xp.FACTORED_ABOVE = 10 ** 6
        r2 = th.load(case).scattering(qmax=3, points=4, parts="all")
    finally:
        xp.FACTORED_ABOVE = saved
    for k in r.curves:
        assert np.abs(r[k] - r2[k]).max() < 1e-8 * max(1.0, np.abs(r2[k]).max()), k
