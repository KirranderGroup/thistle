"""Factored CAS/CI two-particle density (pyxscat.ci.embed_core(factored=True)): the same
scattering as the four-index array, from PySCF, Molpro and Molcas inputs."""
import os

import numpy as np
import pytest

import pyxscat as xs
from pyxscat import ci, aligned
from pyxscat.core import FactoredDM2, symmetrize_dm2

Q = np.linspace(0.0, 6.0, 13)
KW = dict(screen_tol=1e-16)


def _rel(a, b):
    return float(np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-3 * np.abs(a).max())))


@pytest.fixture(scope="module")
def co():
    from pyscf import gto, scf, mcscf
    mol = gto.M(atom="C 0 0 0; O 0 0 1.13", basis="6-31g*", verbose=0)
    mf = scf.RHF(mol).run()
    mc = mcscf.CASSCF(mf, 6, 6).state_average_([0.5, 0.5]).run()
    return mol, mf, mc


@pytest.mark.parametrize("ncore", [0, 3])
def test_embedding_matches_dense(co, ncore):
    _, _, mc = co
    from pyscf import fci
    d1, d2 = fci.direct_spin1.make_rdm12(mc.ci[0], 6, (3, 3))
    D1, D2 = ci.embed_core(d1, d2, ncore)
    F1, F2 = ci.embed_core(d1, d2, ncore, factored=True)
    assert isinstance(F2, FactoredDM2)
    assert F2.nterms[0] <= 6 * 7 // 2 + 3 and F2.nterms[1] == (3 if ncore else 0)
    assert np.abs(D1 - F1).max() == 0.0
    assert np.abs(symmetrize_dm2(D2) - F2.dense()).max() < 1e-13
    n = 6 + 2 * ncore
    assert F2.trace() == pytest.approx(n * (n - 1), rel=1e-13)


def test_active_factorisation_rank(co):
    _, _, mc = co
    from pyscf import fci
    _, d2 = fci.direct_spin1.make_rdm12(mc.ci[1], 6, (3, 3))
    w, v = ci.factor_active(d2)
    assert len(w) <= 21
    g = np.einsum("k,kpq,krs->pqrs", w, v, v)
    assert np.abs(g - symmetrize_dm2(d2)).max() < 1e-13


@pytest.mark.parametrize("state", [0, 1])
def test_pyscf_casscf_scattering(co, state):
    from pyxscat.interfaces import pyscf as xp
    _, _, mc = co
    wd = xp.from_casscf(mc, state=state, factored=False)
    wf = xp.from_casscf(mc, state=state, factored=True)
    assert isinstance(wf.dm2, FactoredDM2) and isinstance(wd.dm2, np.ndarray)
    for kind, tol in (("total", 1e-12), ("inelastic", 1e-9)):        # inelastic = total - elastic
        assert _rel(xs.scattering(wd, Q, kind=kind, **KW), xs.scattering(wf, Q, kind=kind, **KW)) < tol
    # j2 (partial alignment) through the factored density as well
    a = xs.scattering(wd, Q[1:], l=2, **KW)
    b = xs.scattering(wf, Q[1:], l=2, **KW)
    assert np.max(np.abs(a - b)) < 1e-10 * np.max(np.abs(a))


def test_pyscf_reference_array(co):
    """Against PySCF's own AO array (mcscf.addons.make_rdm12, all MOs)."""
    from pyscf import mcscf
    from pyscf.mcscf import addons
    from pyxscat.interfaces import pyscf as xp
    mol, mf, _ = co
    mc = mcscf.CASSCF(mf, 6, 6).run()
    dm1, dm2 = addons.make_rdm12(mc)
    ref = xs.scattering(xp.from_density(mol, dm1, dm2), Q, **KW)
    assert _rel(ref, xs.scattering(xp.from_casscf(mc, factored=True), Q, **KW)) < 1e-12


def test_pyscf_fci(co):
    from pyscf import gto, scf, fci
    from pyxscat.interfaces import pyscf as xp
    mol = gto.M(atom="O 0 0 0; H 0 .757 .587; H 0 -.757 .587", basis="sto-3g", verbose=0)
    mf = scf.RHF(mol).run()
    _, c = fci.FCI(mf).kernel()
    a = xs.scattering(xp.from_fci(mol, mf.mo_coeff, c, 7, (5, 5), factored=False), Q, **KW)
    b = xs.scattering(xp.from_fci(mol, mf.mo_coeff, c, 7, (5, 5), factored=True), Q, **KW)
    assert _rel(a, b) < 1e-12


def test_molpro_factored():
    from pyxscat.interfaces import molpro
    d = os.path.join(os.path.dirname(__file__), "data", "molpro")
    out, mld = os.path.join(d, "lif_sa4_casscf.out"), os.path.join(d, "lif_sa4_casscf.mld")
    wd = molpro.from_multi(out, mld, state=2, factored=False)
    wf = molpro.from_multi(out, mld, state=2, factored=True)
    assert isinstance(wf.dm2, FactoredDM2)
    assert _rel(xs.scattering(wd, Q, **KW), xs.scattering(wf, Q, **KW)) < 1e-12


@pytest.mark.parametrize("frozen", [None, 1])
def test_cisd_factored_density(frozen):
    """CISD: the factored density against PySCF's make_rdm2, blockwise and through the kernels."""
    from pyscf import gto, scf, ci as pci
    from pyxscat import ccsd
    from pyxscat.interfaces import pyscf as xp
    mol = gto.M(atom="O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", basis="6-31g*", verbose=0)
    myci = pci.CISD(scf.RHF(mol).run(conv_tol=1e-12))
    myci.frozen = frozen
    myci.run(conv_tol=1e-10)
    f = ccsd.factored_rdm2_cisd(myci)
    assert np.abs(symmetrize_dm2(myci.make_rdm2()) - f.dense()).max() < 1e-13
    assert sum(f.nterms) == 4 * myci.nocc ** 2 + 2              # nocc: active occupied
    a = xs.scattering(xp.from_ci(myci, factored=False), Q, **KW)
    b = xs.scattering(xp.from_ci(myci, factored=True), Q, **KW)
    assert _rel(a, b) < 1e-12


def test_factored_switch():
    """thistle.load(factored=...) and the input-file / command-line values."""
    import thistle as th
    from thistle.inputfile import _factored
    out = os.path.join(os.path.dirname(__file__), "data", "molpro", "lif_sa4_casscf.out")
    assert isinstance(th.load(out).wf.dm2, FactoredDM2)                  # auto: 90 functions
    assert isinstance(th.load(out, factored=False).wf.dm2, np.ndarray)
    assert "factored" in repr(th.load(out))
    assert [_factored(v) for v in ("auto", "TRUE", "false", True, False)] == ["auto", True, False, True, False]
    with pytest.raises(ValueError):
        _factored("maybe")
