"""Molden reader: conventions, round trips through PySCF and the example files."""
import os

import numpy as np
import pytest

import pyxscat as xs
from pyxscat.interfaces import molden as M
from pyxscat.interfaces import pyscf as xp
from conftest import slow

DATA = os.path.join(os.path.dirname(__file__), "data")
H2O = "O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692"


def test_solid_harmonics():
    assert M.solid_harmonic(2, 0) == {(0, 0, 2): 1.0, (2, 0, 0): -0.5, (0, 2, 0): -0.5}
    assert M.solid_harmonic(2, 2) == {(2, 0, 0): 1.0, (0, 2, 0): -1.0}
    assert M.solid_harmonic(2, -2) == {(1, 1, 0): 2.0}
    assert M.solid_harmonic(2, 1) == {(1, 0, 1): 2.0}
    assert M.solid_harmonic(2, -1) == {(0, 1, 1): 2.0}


@pytest.mark.parametrize("basis,cart", [("cc-pvtz", False), ("cc-pvtz", True), ("cc-pvqz", False),
                                        ("6-31g*", True), ("def2-svp", False)])
def test_pyscf_round_trip(tmp_path, basis, cart):
    """PySCF -> Molden -> pyxscat gives the same density and the same scattering as the in-memory adapter."""
    from pyscf import gto, scf
    from pyscf.tools import molden as pm
    mol = gto.M(atom=H2O, basis=basis, cart=cart, verbose=0)
    mf = scf.RHF(mol).run()
    f = tmp_path / "x.molden"
    pm.from_scf(mf, str(f))
    data = M.read_molden(str(f))
    assert M.check_orbitals(data) < 1e-10
    pts = np.random.default_rng(0).normal(scale=1.3, size=(50, 3))
    ao = mol.eval_gto("GTOval_cart" if cart else "GTOval_sph", pts)
    rho_ref = np.einsum("pi,ij,pj->p", ao, mf.make_rdm1(), ao)
    a2 = data.basis.eval(pts)
    assert np.allclose(np.einsum("pi,ij,pj->p", a2, data.density(), a2), rho_ref, atol=1e-11)
    q = np.array([0.0, 1.0, 4.0])
    wf_m = data.wavefunction()
    wf_p = xp.from_scf(mf)
    assert np.allclose(xs.elastic(wf_m.basis, q, wf_m.dm1), xs.elastic(wf_p.basis, q, wf_p.dm1), rtol=1e-9)


EXAMPLES = [("molden/c5h6_orbitals.molden", 32), ("molden/co_uhf_pyscf.molden", 14),
            ("molcas/casscf_c7h8.rasscf.molden", 50)]


@pytest.mark.parametrize("name,nelec", EXAMPLES)
def test_example_files(name, nelec):
    data = M.read_molden(os.path.join(DATA, name))
    # Molcas and Molpro print 8-10 significant digits
    assert M.check_orbitals(data) < 1e-7
    assert data.nelec == pytest.approx(nelec)
    s = M.overlap_matrix(data.basis)
    assert np.einsum("ij,ji", data.density(), s) == pytest.approx(nelec, rel=1e-7)


@slow
@pytest.mark.parametrize("name,nelec", EXAMPLES)
def test_example_files_scattering(name, nelec):
    wf = M.read_molden(os.path.join(DATA, name)).wavefunction()
    assert xs.elastic(wf.basis, [0.0], wf.dm1)[0] == pytest.approx(nelec ** 2, rel=1e-6)


def test_unrestricted_molden_density():
    data = M.read_molden(os.path.join(DATA, "molden/co_uhf_pyscf.molden"))
    assert set(data.mo_spin) == {"Alpha", "Beta"}
    d = data.density()
    assert np.allclose(d, data.density("Alpha") + data.density("Beta"))


@pytest.mark.parametrize("basis,cart", [("cc-pvtz", False), ("cc-pvqz", False), ("6-31g*", True)])
def test_scf_energy_check(tmp_path, basis, cart):
    """The Molden convention check: the SCF energy rebuilt from the file equals the program's."""
    from pyscf import gto, scf
    from pyscf.tools import molden as pm
    mol = gto.M(atom=H2O, basis=basis, cart=cart, verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-11)
    f = tmp_path / "x.molden"
    pm.from_scf(mf, str(f))
    e, resid = M.scf_energy(str(f))
    assert resid < 1e-10
    assert e == pytest.approx(mf.e_tot, abs=1e-8)
