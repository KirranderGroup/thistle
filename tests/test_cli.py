import numpy as np
import pytest

from pyxscat import cli
from conftest import H2O


def test_command_line_elastic_from_molden(tmp_path):
    from pyscf import gto, scf
    from pyscf.tools import molden as pm
    mol = gto.M(atom=H2O, basis="6-31g", verbose=0)
    mf = scf.RHF(mol).run()
    f = tmp_path / "w.molden"
    pm.from_scf(mf, str(f))
    out = tmp_path / "el.dat"
    cli.main(["--molden", str(f), "--kind", "elastic", "--qmax", "5", "--npoints", "6", "-o", str(out)])
    d = np.loadtxt(out)
    assert d.shape == (6, 2)
    assert d[0, 1] == pytest.approx(100.0, rel=1e-9)
    out2 = tmp_path / "hf.dat"
    cli.main(["--molden", str(f), "--kind", "hf_total", "--qmax", "5", "--npoints", "6", "-o", str(out2)])
    assert np.loadtxt(out2)[0, 1] == pytest.approx(100.0, rel=1e-9)


def test_command_line_refuses_total_without_dm2(tmp_path):
    with pytest.raises(SystemExit):
        cli.main(["--molden", "x.molden", "--kind", "total"])


def test_command_line_aligned(tmp_path):
    from pyscf import gto, scf
    from pyscf.tools import molden as pm
    from pyxscat import aligned
    from pyxscat.interfaces import pyscf as xp
    mol = gto.M(atom=H2O, basis="6-31g", verbose=0)
    mf = scf.RHF(mol).run()
    f = tmp_path / "w.molden"
    pm.from_scf(mf, str(f))
    q = np.array([[0.0, 0.0, 0.0], [0.5, -0.2, 1.0], [2.0, 1.0, -0.3]])
    np.savetxt(tmp_path / "q.txt", q)
    out = tmp_path / "al.dat"
    cli.main(["--molden", str(f), "--qvec", str(tmp_path / "q.txt"), "--orientation", "0.2", "0.4", "0.6",
              "--kind", "elastic", "-o", str(out)])
    d = np.loadtxt(out)
    ref = aligned.elastic(xp.basis_from_mole(mol), aligned.rotate_q(q, aligned.rotation_matrix(0.2, 0.4, 0.6)),
                          mf.make_rdm1() if mol.cart else xp.from_scf(mf).dm1)
    assert np.allclose(d[:, 3], ref, rtol=1e-8)
    assert d[0, 3] == pytest.approx(100.0, rel=1e-9)
