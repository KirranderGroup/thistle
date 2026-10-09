import numpy as np
import pytest

import pyxscat as xs
from pyxscat.interfaces import pyscf as xp
from conftest import get_case, H2O


@pytest.mark.parametrize("basis", ["sto-3g", "6-31g*", "cc-pvdz", "cc-pvtz"])
def test_cartesian_functions_match_pyscf(basis):
    """The primitive expansion reproduces PySCF's Cartesian basis functions (s to f, general contraction)."""
    from pyscf import gto
    mol = gto.M(atom=H2O, basis=basis, cart=True, verbose=0)
    b = xp.basis_from_mole(mol)
    pts = np.random.default_rng(1).normal(scale=1.5, size=(40, 3))
    assert b.nao == mol.nao
    assert np.allclose(b.eval(pts), mol.eval_gto("GTOval_cart", pts), atol=1e-13)


def test_spherical_density_is_mapped_to_cartesian_functions():
    """For a pure (spherical) basis the density built in Cartesian functions is the same function of r."""
    case = get_case(H2O, "cc-pvdz", False)
    wf = xp.from_scf(case.mf)
    pts = np.random.default_rng(2).normal(scale=1.5, size=(40, 3))
    ao_sph = case.mol.eval_gto("GTOval_sph", pts)
    rho_ref = np.einsum("pi,ij,pj->p", ao_sph, case.mf.make_rdm1(), ao_sph)
    ao = wf.basis.eval(pts)
    rho = np.einsum("pi,ij,pj->p", ao, wf.dm1, ao)
    assert np.allclose(rho, rho_ref, atol=1e-12)


def test_from_shells_numbering():
    shells = [((0.0, 0.0, 0.0), 0, [3.0, 1.0], [[0.5, 0.1], [0.6, 0.9]]),     # two s contractions
              ((0.0, 0.0, 1.0), 1, [0.8], [1.0])]                             # one p shell
    b = xs.Basis.from_shells(shells)
    assert b.nao == 2 + 3
    assert b.nprim == 2 * 2 + 3
    assert list(b.group) == [1, 1, 2, 2, 3, 3, 3]
    assert list(b.contr) == [1, 2, 1, 2, 3, 4, 5]


def test_general_contraction_written_as_segments_is_merged():
    seg = [((0.0, 0.0, 0.0), 0, [3.0, 1.0], [0.5, 0.6]),
           ((0.0, 0.0, 0.0), 0, [3.0, 1.0], [0.1, 0.9]),
           ((0.0, 0.0, 1.0), 1, [0.8], [1.0])]
    gen = [((0.0, 0.0, 0.0), 0, [3.0, 1.0], [[0.5, 0.1], [0.6, 0.9]]),
           ((0.0, 0.0, 1.0), 1, [0.8], [1.0])]
    a = xs.Basis.from_shells(seg)
    b = xs.Basis.from_shells(gen)
    assert a.group.max() == b.group.max() == 3
    pts = np.random.default_rng(0).normal(size=(10, 3))
    assert np.allclose(a.eval(pts), b.eval(pts))
    c = xs.Basis.from_shells(seg, merge=False)
    assert c.group.max() == 5 and np.allclose(c.eval(pts), a.eval(pts))
