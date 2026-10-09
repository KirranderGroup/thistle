"""Independent atom model."""
import numpy as np
import pytest

import pyxscat as xs
from pyxscat import aligned, iam
from pyxscat.interfaces import pyscf as xp
from conftest import get_case, H2O

Q = np.array([0.0, 0.5, 1.0, 2.0, 4.0, 8.0])


@pytest.mark.parametrize("el", ["H", "C", "O", "F", "Ne", "S"])
def test_ab_initio_form_factors_match_tables(el):
    """Free-atom Hartree-Fock form factors (aug-cc-pVQZ) against Waasmaier-Kirfel: within 1 %."""
    pytest.importorskip("xraydb")
    fa = iam.atomic_form_factor(el, Q, basis="aug-cc-pvqz")
    ft = iam.atomic_form_factor(el, Q, source="tables")
    z = iam.ELEMENTS.index(el)
    assert fa[0] == pytest.approx(z, rel=1e-9)
    assert np.all(np.abs(fa - ft) < 0.01 * z + 1e-3)


def test_single_atom_iam_equals_molecular_code():
    """A closed-shell atom is its own IAM: elastic and total through the molecular kernels agree."""
    from pyscf import gto, scf
    mol = gto.M(atom="Ne 0 0 0", basis="aug-cc-pvtz", verbose=0)
    wf = xp.from_scf(scf.RHF(mol).run(conv_tol=1e-11))
    ref = iam.for_wavefunction(wf, Q, basis="aug-cc-pvtz")
    assert np.allclose(xs.elastic(wf.basis, Q, wf.dm1), ref["elastic"], rtol=1e-8, atol=1e-8)
    assert np.allclose(xs.scattering(wf, Q, "hf_total"), ref["total"], rtol=1e-8, atol=1e-8)


def test_limits_and_partition(h2o_631g):
    wf = xp.from_scf(h2o_631g.mf)
    r = iam.for_wavefunction(wf, Q, basis="6-31g")
    n = wf.nelec
    assert r["elastic"][0] == pytest.approx(n * n, rel=1e-9)
    assert r["total"][0] == pytest.approx(n * n, rel=1e-9)
    assert np.allclose(r["total"], r["elastic"] + r["inelastic"])
    assert np.all(r["inelastic"][1:] > 0)


def test_aligned_iam_averages_to_debye(h2o_631g):
    from scipy.integrate import lebedev_rule
    wf = xp.from_scf(h2o_631g.mf)
    xyz, w = lebedev_rule(41)
    w = w / w.sum()
    qs = np.array([0.7, 2.0, 5.0])
    deb = iam.for_wavefunction(wf, qs, basis="6-31g")
    for i, q in enumerate(qs):
        al = iam.for_wavefunction(wf, q * xyz.T, aligned=True, basis="6-31g")
        assert w @ al["elastic"] == pytest.approx(deb["elastic"][i], rel=1e-9)
        assert w @ al["total"] == pytest.approx(deb["total"][i], rel=1e-9)


def test_molecule_approaches_iam_at_high_q():
    """Bonding and correlation are valence effects: the molecular total approaches the IAM at large
    q, while at intermediate q the difference is of the order of per cent."""
    case = get_case(H2O, "aug-cc-pvtz", True)
    wf = xp.from_cc(case.cc)
    q = np.array([1.0, 2.0, 10.0, 14.0])
    mol = xs.scattering(wf, q, "total")
    ref = iam.for_wavefunction(wf, q, basis="aug-cc-pvtz")["total"]
    pd = iam.percent_difference(mol, ref)
    assert abs(pd[-1]) < 0.5 and abs(pd[-2]) < 1.0
    assert np.abs(pd[:2]).max() > 0.5


def test_tables_refuse_total():
    pytest.importorskip("xraydb")
    with pytest.raises(ValueError):
        iam.scattering([8, 1, 1], np.zeros((3, 3)), Q, source="tables")
    r = iam.scattering([8, 1, 1], [[0, 0, 0], [0, 1.4, -1.1], [0, -1.4, -1.1]], Q, source="tables", total=False)
    assert r["elastic"][0] == pytest.approx(100, rel=1e-3)
