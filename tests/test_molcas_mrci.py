"""OpenMolcas uncontracted MRCI (water/cc-pVDZ, CASSCF(6,5) reference, O 1s frozen, THRPrint = 0;
OpenMolcas 24.10).  The log keeps only the CSFs with non-zero coefficients."""
import os

import numpy as np
import pytest

import pyxscat as xs
from pyxscat.core import FactoredDM2
from pyxscat.interfaces import molcas

D = os.path.join(os.path.dirname(__file__), "data", "molcas_mrci")
LOG, MLD = os.path.join(D, "water_mrci.log"), os.path.join(D, "water_mrci.rasscf.molden")
E_MRCI = -76.23719469                      # "CI ENERGY" in the log
OCC = [2.0, 1.98511, 1.97382, 1.96615, 1.96324, 0.02804, 0.02595, 0.01623, 0.01082, 0.00563]


@pytest.fixture(scope="module")
def info():
    return molcas.read_mrci_output(LOG)


@pytest.fixture(scope="module")
def rdms(info):
    return molcas.mrci_rdms(info, 1)


def test_parsing(info):
    assert (info["nfrozen"], info["norb_ci"], info["ninternal"], info["nelec"], info["spin"]) == (1, 23, 6, 8, 1)
    s = info["states"][1]
    assert s["energy"] == pytest.approx(E_MRCI)
    assert len(s["csfs"]) == 12626
    assert s["coeffs"] @ s["coeffs"] == pytest.approx(1.0, abs=1e-5)       # six printed decimals
    assert molcas.has_mrci(LOG)


def test_natural_occupations(rdms):
    dm1, dm2 = rdms
    w = np.sort(np.linalg.eigvalsh(dm1))[::-1]
    assert np.allclose(w[:10], OCC, atol=6e-6)
    assert np.trace(dm1) == pytest.approx(10.0)
    assert np.einsum("iijj->", dm2) == pytest.approx(90.0)


def test_energy_reproduced(rdms):
    """The MRCI energy from our RDMs with PySCF's integrals in the Molden orbitals."""
    from pyscf import ao2mo
    from pyscf.tools import molden as pm
    dm1, dm2 = rdms
    mol, _, moc, _, _, _ = pm.load(MLD)
    mol.verbose = 0
    c = moc[:, :24]
    h = c.T @ (mol.intor("int1e_kin") + mol.intor("int1e_nuc")) @ c
    eri = ao2mo.restore(1, ao2mo.kernel(mol, c), 24)
    e = np.einsum("pq,pq", h, dm1) + 0.5 * np.einsum("pqrs,pqrs", eri, dm2) + mol.energy_nuc()
    assert e == pytest.approx(E_MRCI, abs=3e-7)


def test_scattering_and_load():
    import thistle as th
    wd = molcas.from_mrci(LOG, MLD, factored=False)
    wf = molcas.from_mrci(LOG, MLD, factored=True)
    assert isinstance(wf.dm2, FactoredDM2) and wd.nelec == 10
    q = np.linspace(0.0, 6.0, 7)
    a, b = xs.scattering(wd, q, screen_tol=1e-16), xs.scattering(wf, q, screen_tol=1e-16)
    assert a[0] == pytest.approx(100.0, rel=1e-8)
    assert np.max(np.abs(a - b) / a) < 1e-10
    # the MRCI inelastic signal differs from the CASSCF one (dynamic correlation of the 2s pair)
    wc = molcas.from_rasscf(LOG, MLD)
    assert abs(xs.scattering(wd, [2.0], "inelastic")[0] - xs.scattering(wc, [2.0], "inelastic")[0]) > 0.1
    m = th.load(LOG)
    assert m.origin["kind"] == "molcas_mrci" and m.correlated


def test_structured_builder_equals_pair_loop(info):
    """The restricted pair loop + products (any number of orbitals) against the full pair loop,
    also with the virtuals spread over 120 orbitals (two 64-bit words, atomic accumulation)."""
    from pyxscat import ci
    a1, a2 = molcas.mrci_rdms(info, 1, method="bitwise")
    b1, b2 = molcas.mrci_rdms(info, 1, method="external")
    assert np.abs(a1 - b1).max() < 1e-12 and np.abs(a2 - b2).max() < 1e-12
    s = info["states"][1]
    acc = {}
    for (orbs, step), cv in zip(s["csfs"], s["coeffs"]):
        pos = [o - 2 for o in orbs]
        for v, ab, bb in ci.csf_to_determinants(step, 0):
            a, b, sign = molcas._place(ab, bb, pos)
            acc[a, b] = acc.get((a, b), 0.0) + sign * v * cv
    dets, co = list(acc), np.array(list(acc.values()))
    r1, r2 = ci.make_rdm12(23, 4, 4, dets, co, method="bitwise")
    vmap = list(range(6)) + sorted(np.random.default_rng(1).choice(range(6, 120), 17, replace=False).tolist())
    remap = lambda x: sum(1 << vmap[k] for k in range(23) if x >> k & 1)
    e1, e2 = ci.make_rdm12_external(120, 6, [(remap(a), remap(b)) for a, b in dets], co)
    assert np.abs(e1[np.ix_(vmap, vmap)] - r1).max() < 1e-12
    assert np.abs(e2[np.ix_(vmap, vmap, vmap, vmap)] - r2).max() < 1e-12
    assert abs(np.abs(e2).sum() - np.abs(r2).sum()) < 1e-10
