"""Different routes through the code must give the same numbers."""
import numpy as np
import pytest

import pyxscat as xs
from pyxscat.interfaces import pyscf as xp

Q = np.linspace(0.0, 8.0, 17)


def test_mean_field_on_the_fly_equals_explicit_dm2(h2o_631g):
    wf = xp.from_scf(h2o_631g.mf, with_dm2=True)
    explicit = xs.two_electron(wf.basis, Q, wf.dm2, screen_tol=0.0)
    on_the_fly = xs.mean_field(wf.basis, Q, wf.dm1, screen_tol=0.0)
    assert np.allclose(on_the_fly, explicit, rtol=1e-11, atol=1e-11)


def test_hartree_plus_exchange_is_mean_field(h2o_631g):
    wf = xp.from_scf(h2o_631g.mf)
    h = xs.elastic(wf.basis, Q, wf.dm1, screen_tol=0.0)
    x = xs.exchange(wf.basis, Q, wf.dm1, screen_tol=0.0)
    mf = xs.mean_field(wf.basis, Q, wf.dm1, screen_tol=0.0)
    assert np.allclose(h + x, mf, rtol=1e-11, atol=1e-11)


def test_decomposition_sums_to_total(h2o_631g):
    wf = xp.from_cc(h2o_631g.cc)
    d = xs.decompose(wf, Q)
    assert np.allclose(d["self"] + d["hartree"] + d["exchange"] + d["cumulant"], d["total"], atol=1e-10)
    assert np.allclose(d["total"], xs.scattering(wf, Q, "total"), atol=1e-12)
    assert np.allclose(d["inelastic"], xs.scattering(wf, Q, "inelastic"), atol=1e-12)
    # electron correlation lowers the pair density at short range: the cumulant is not zero
    assert np.abs(d["cumulant"]).max() > 1e-3


def test_hf_cumulant_vanishes(h2o_631g):
    wf = xp.from_scf(h2o_631g.mf, with_dm2=True)
    d = xs.decompose(wf, Q)
    assert np.abs(d["cumulant"]).max() < 1e-9


def test_cartesian_and_spherical_route_agree():
    """The same wavefunction through a spherical Mole and through its Cartesian representation."""
    from conftest import get_case, H2O
    case = get_case(H2O, "cc-pvdz", False)
    wf = xp.from_scf(case.mf)
    c2s = case.mol.cart2sph_coeff(normalized="sp")
    assert np.allclose(wf.dm1, c2s @ case.mf.make_rdm1() @ c2s.T)
    n = case.mol.nelectron
    assert xs.elastic(wf.basis, [0.0], wf.dm1)[0] == pytest.approx(n * n, rel=1e-9)


def test_symmetrisation_does_not_change_result(h2o_631g):
    wf = xp.from_cc(h2o_631g.cc)
    a = xs.two_electron(wf.basis, Q, wf.dm2)
    b = xs.two_electron(wf.basis, Q, xs.symmetrize_dm2(wf.dm2), symmetrize=False)
    assert np.allclose(a, b, rtol=1e-12, atol=1e-12)


def test_screening_converges(h2o_631g):
    wf = xp.from_cc(h2o_631g.cc)
    exact, s0 = xs.two_electron(wf.basis, Q, wf.dm2, screen_tol=0.0, return_stats=True)
    tight = xs.two_electron(wf.basis, Q, wf.dm2, screen_tol=1e-12)
    loose, s1 = xs.two_electron(wf.basis, Q, wf.dm2, screen_tol=1e-6, return_stats=True)
    assert np.abs(tight - exact).max() < 1e-8
    assert np.abs(loose - exact).max() < 1e-2
    assert s1["primitive_quartets_computed"] < s0["primitive_quartets_computed"]


def test_result_independent_of_thread_count(h2o_631g):
    wf = xp.from_cc(h2o_631g.cc)
    one = xs.two_electron(wf.basis, Q, wf.dm2, nthreads=1)
    two = xs.two_electron(wf.basis, Q, wf.dm2, nthreads=2)
    assert np.allclose(one, two, rtol=1e-11, atol=1e-11)


def test_bad_input_is_reported(h2o_631g):
    wf = xp.from_scf(h2o_631g.mf)
    with pytest.raises(ValueError):
        xs.scattering(wf, Q, "total")            # no dm2
    with pytest.raises(ValueError):
        xs.elastic(wf.basis, Q, wf.dm1[:-1, :-1])
    with pytest.raises(TypeError):
        xs.elastic(wf.basis, Q, wf.dm1, cutoff=1.0)
