"""Molcas RASSCF interface on a stored CASSCF(8,8) calculation (C7H8, ANO-L-MB, OpenMolcas 22.10)."""
import os

import numpy as np
import pytest

import pyxscat as xs
from pyxscat import ci
from pyxscat.interfaces import molcas
from conftest import slow

DATA = os.path.join(os.path.dirname(__file__), "data", "molcas")
LOG = os.path.join(DATA, "casscf_c7h8.log")
MLD = os.path.join(DATA, "casscf_c7h8.rasscf.molden")


@pytest.fixture(scope="module")
def info():
    return molcas.read_rasscf_output(LOG)


def test_output_parsing(info):
    assert (info["ninactive"], info["nactive"], info["nactive_electrons"], info["spin"]) == (21, 8, 8, 0.0)
    r = info["roots"][1]
    assert len(r["csfs"]) == 1764
    assert r["energy"] == pytest.approx(-269.24529814)
    assert r["coeffs"] @ r["coeffs"] == pytest.approx(1.0, abs=1e-4)      # five printed decimals


def test_natural_occupations_reproduced(info):
    """The CI vector is printed in the natural-orbital basis: our 1RDM must be diagonal there,
    with the occupations Molcas prints."""
    dm1, dm2 = molcas.rasscf_rdms(info)
    act = dm1[21:, 21:]
    assert np.abs(act - np.diag(np.diag(act))).max() < 2e-5
    assert np.allclose(np.diag(act), info["roots"][1]["occupations"], atol=5e-6)
    n = 50
    assert np.trace(dm1) == pytest.approx(n)
    assert np.einsum("iijj->", dm2) == pytest.approx(n * (n - 1))


def test_ci_vector_is_the_casci_eigenvector(info):
    """In the same orbitals and basis, the converted determinant vector is PySCF's CASCI ground state.
    (Absolute energies differ by about 24 mEh, as do the SCF energies of the two programs: the
    Molcas run used RI two-electron integrals (RICD).)"""
    from pyscf import fci, mcscf, scf
    from pyscf.fci import cistring
    from pyscf.tools import molden as pm
    mol, moe, moc, occ, _, _ = pm.load(MLD)
    mol.verbose = 0
    mf = scf.RHF(mol)
    mc = mcscf.CASCI(mf, 8, 8)
    mc.verbose = 0
    mc.kernel(moc)
    r = info["roots"][1]
    dets, co = ci.csfs_to_determinants(r["csfs"], r["coeffs"], 0)
    c, (sa, ia), _ = ci.ci_matrix(8, 4, 4, dets, co)
    p = [ia[int(s)] for s in cistring.make_strings(range(8), 4)]
    cp = c[np.ix_(p, p)] / np.linalg.norm(c)
    assert abs(np.einsum("ij,ij", cp, mc.ci)) > 0.99999


def test_wavefunction_in_ao_basis():
    from pyxscat.interfaces.molden import overlap_matrix
    wf = molcas.from_rasscf(LOG, MLD)
    s = overlap_matrix(wf.basis)
    assert np.einsum("ij,ji", wf.dm1, s) == pytest.approx(50, rel=1e-7)
    assert np.einsum("ijkl,ji,lk", wf.dm2, s, s) == pytest.approx(50 * 49, rel=1e-7)
    assert wf.nelec == 50


@slow
def test_scattering_limits():
    wf = molcas.from_rasscf(LOG, MLD)
    d = xs.decompose(wf, [0.0, 2.0], screen_tol=1e-10)
    assert d["total"][0] == pytest.approx(2500, rel=1e-6)
    assert d["elastic"][0] == pytest.approx(2500, rel=1e-6)
    assert d["inelastic"][1] > 0


@slow
def test_command_line_molcas(tmp_path):
    from pyxscat import cli
    out = tmp_path / "res.dat"
    cli.main(["--molden", MLD, "--molcas-log", LOG, "--kind", "decompose", "--qmax", "4", "--npoints", "3",
              "--screen-tol", "1e-10", "-o", str(out)])
    data = np.loadtxt(out)
    assert data.shape == (3, 7)
    assert data[0, 1] == pytest.approx(2500, rel=1e-6)
