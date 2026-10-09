"""The simple interface (``import thistle``) and the ``thistle`` command."""
import os
import shutil

import numpy as np
import pytest

import pyxscat as xs
import thistle as th
from thistle import inputfile

DATA = os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture(scope="module")
def water_scf():
    from pyscf import gto, scf
    mol = gto.M(atom="O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", basis="6-31g", verbose=0)
    return scf.RHF(mol).run(conv_tol=1e-11)


@pytest.fixture(scope="module")
def water_cc(water_scf):
    from pyscf import cc
    return cc.CCSD(water_scf).run(conv_tol=1e-9)


def test_sniff():
    assert th.core._sniff(os.path.join(DATA, "molden", "co_uhf_pyscf.molden")) == "molden"
    assert th.core._sniff(os.path.join(DATA, "molpro", "lif_sa4_casscf.out")) == "molpro_output"
    assert th.core._sniff(os.path.join(DATA, "molpro", "sf6p_cas56.pun")) == "molpro_punch"
    assert th.core._sniff(os.path.join(DATA, "molcas", "casscf_c7h8.log")) == "molcas_output"


def test_units_and_grid(water_scf):
    m = th.load(water_scf)
    a = m.scattering(q=[0.0, 1.0, 2.0], unit="1/A")
    b = m.scattering(q=np.array([0.0, 1.0, 2.0]) * th.BOHR, unit="1/bohr")
    for k in a.curves:
        np.testing.assert_allclose(a[k], b[k], rtol=1e-12, atol=1e-12)
    assert a.unit == "1/A" and list(a.q) == [0.0, 1.0, 2.0]


def test_scf_gives_total_and_matches_pyxscat(water_scf):
    m = th.load(water_scf)
    assert not m.correlated
    r = m.scattering(qmax=3, points=7)
    assert set(r.curves) == {"elastic", "total", "inelastic"}
    qb = r.q * th.BOHR
    np.testing.assert_allclose(r.total, xs.scattering(m.wf, qb, "hf_total"), rtol=1e-10)
    assert abs(r.elastic[0] - 100.0) < 1e-6 and abs(r.total[0] - 100.0) < 1e-6


def test_natural_orbitals_give_elastic_only():
    m = th.load(os.path.join(DATA, "molpro", "lif_sa4_casscf.mld"))
    assert set(m.scattering(q=[0.5]).curves) == {"elastic"}
    with pytest.raises(ValueError, match="closed-shell Hartree-Fock"):
        th.load(os.path.join(DATA, "molpro", "lif_sa4_casscf.mld"), hf=True)
    with pytest.raises(ValueError, match="closed-shell Hartree-Fock"):      # UHF orbitals
        th.load(os.path.join(DATA, "molden", "co_uhf_pyscf.molden"), hf=True)


def test_hf_molden_needs_the_option(tmp_path, water_scf):
    from pyscf.tools import molden
    path = str(tmp_path / "water_rhf.molden")
    molden.from_scf(water_scf, path)
    plain = th.load(path).scattering(q=[0.0, 1.0, 2.0])
    assert set(plain.curves) == {"elastic"}
    hf = th.load(path, hf=True).scattering(q=[0.0, 1.0, 2.0], parts="all")
    ref = th.load(water_scf).scattering(q=[0.0, 1.0, 2.0])
    for k in ("total", "elastic", "inelastic"):
        np.testing.assert_allclose(hf[k], ref[k], rtol=1e-7, atol=1e-7)
    np.testing.assert_allclose(hf.total, 10 + hf.hartree + hf.exchange, rtol=1e-12)
    prefix = str(tmp_path / "cli")
    inputfile.main(["run", path, "--hf", "--qmax", "2", "--points", "3", "--out", prefix])
    head = open(prefix + "_scattering.csv").read().splitlines()[1]
    assert head.split(",")[1:] == ["total", "elastic", "inelastic"]


def test_fci_solver():
    """A standalone fci.FCI(mf) solver: any root, orbitals given (PySCF keeps only integrals)."""
    from pyscf import gto, scf, fci
    from pyxscat.interfaces import pyscf as xp
    mol = gto.M(atom="O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", basis="sto-3g", verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-12)
    s = fci.FCI(mf); s.nroots = 2; s.kernel()
    with pytest.raises(ValueError, match="mo_coeff"):
        th.load(s, state=2)
    m = th.load(s, state=2, mo_coeff=mf)
    d1, d2 = fci.direct_spin1.make_rdm12(s.ci[1], s.norb, s.nelec)
    ref = xp.from_density(mol, d1, d2, mo_coeff=mf.mo_coeff)
    q = np.array([1e-8, 0.5, 2.0])
    np.testing.assert_allclose(xs.two_electron(m.wf.basis, q, m.wf.dm2),
                               xs.two_electron(ref.basis, q, ref.dm2), rtol=1e-12)
    r = th.load(s, state=1, mo_coeff=mf.mo_coeff).scattering(qmax=1, points=3)
    assert abs(r.total[0] - 100) < 1e-6


def test_ccsd_parts(water_cc):
    m = th.load(water_cc)
    assert m.correlated
    r = m.scattering(qmax=3, points=5, parts="all")
    np.testing.assert_allclose(r.total, r.elastic + r.inelastic, rtol=1e-12)
    np.testing.assert_allclose(r.total, 10 + r.hartree + r.exchange + r.cumulant, rtol=1e-10)
    assert abs(r.total[0] - 100) < 1e-5


def test_molpro_output_finds_molden(recwarn):
    m = th.load(os.path.join(DATA, "molpro", "lif_sa4_casscf.out"), state=2)
    assert m.correlated and m.origin["molden"].endswith("lif_sa4_casscf.mld")
    assert m.origin["state"] == 2


def test_molpro_output_prefers_punch():
    m = th.load(os.path.join(DATA, "molpro", "he_sa2_casscf.out"))
    ref = xs.interfaces.molpro.from_multi(os.path.join(DATA, "molpro", "he_sa2_casscf.out"),
                                          os.path.join(DATA, "molpro", "he_sa2_casscf.mld"),
                                          punch=os.path.join(DATA, "molpro", "he_sa2_casscf.pun"))
    np.testing.assert_allclose(m.wf.dm2, ref.dm2, atol=1e-14)


def test_ambiguous_molden_is_an_error(tmp_path):
    for f in ("he_sa2_casscf.out", "he_sa2_casscf.mld", "lif_sa4_casscf.mld"):
        shutil.copy(os.path.join(DATA, "molpro", f), tmp_path / f.replace("he_sa2_casscf.mld", "orbs.mld"))
    with pytest.raises(FileNotFoundError, match="several Molden"):
        th.load(str(tmp_path / "he_sa2_casscf.out"))
    # the punch file is deliberately not copied, so the reader falls back on the truncated
    # printout of the output and must say so
    with pytest.warns(UserWarning, match="truncated"):
        m = th.load(str(tmp_path / "he_sa2_casscf.out"), molden=str(tmp_path / "orbs.mld"))
    assert m.correlated


def test_iam_and_compare(water_scf):
    m = th.load(water_scf)
    r, ref = m.scattering(qmax=4, points=5), m.iam(qmax=4, points=5, basis="6-31g")
    d = th.compare(r, ref)
    assert d.curves.keys() == {"%diff total"}
    assert abs(d["%diff total"][0]) < 1e-6           # both N^2 at q = 0
    assert np.all(np.abs(d["%diff total"]) < 20)


def test_anisotropy_matches_l2(water_cc):
    m = th.load(water_cc)
    a = m.anisotropy(q=[0.5, 1.5])
    np.testing.assert_allclose(a["S2 total"], xs.scattering(m.wf, np.array([0.5, 1.5]) * th.BOHR, "total", l=2))


def test_detector(water_cc):
    m = th.load(water_cc)
    img = m.detector(wavelength=1.0, pixels=9, orientation=(0, 90, 0))
    assert img["total"].shape == (9, 9)
    centre = img["total"][4, 4]
    assert abs(centre - 100) < 1e-5                   # forward direction, q = 0


def test_result_save_and_plot(tmp_path, water_scf):
    r = th.load(water_scf).scattering(qmax=2, points=5)
    path = r.save(str(tmp_path / "w.csv"))
    arr = np.loadtxt(path, delimiter=",", skiprows=2)
    cols = open(path).read().splitlines()[1].split(",")
    np.testing.assert_allclose(arr[:, cols.index("elastic")], r.elastic, rtol=1e-10)
    pytest.importorskip("matplotlib")
    r.plot(str(tmp_path / "w.png"))
    assert os.path.getsize(tmp_path / "w.png") > 1000


def test_input_file_workflow(tmp_path, capsys):
    path = str(tmp_path / "lif.toml")
    th.template(path)
    with pytest.raises(FileExistsError):
        th.template(path)
    for f in ("lif_sa4_casscf.out", "lif_sa4_casscf.mld"):
        shutil.copy(os.path.join(DATA, "molpro", f), tmp_path / f)
    text = open(path).read()
    text = text.replace('file  = "my_calculation.out"', 'file  = "lif_sa4_casscf.out"')
    text = text.replace("max    = 4.0", "max    = 2.0").replace("points = 201", "points = 5")
    text = text.replace("iam        = true", "iam        = false")
    text = text.replace("plots  = true", "plots  = false")
    open(path, "w").write(text)
    written = inputfile.run(path, verbose=False)
    assert os.path.exists(tmp_path / "results" / "lif_scattering.csv")
    assert set(written) == {"scattering", "plot_script", "run_log"}          # plots = false
    inputfile.main(["submit", path])
    script = open(tmp_path / "lif.slurm").read()
    assert "#SBATCH --cpus-per-task=16" in script and "thistle run" in script


def test_check_command(capsys):
    inputfile.main(["check", os.path.join(DATA, "molpro", "he_sa2_casscf.out")])
    out = capsys.readouterr().out
    assert out.count("OK") == 2 and "PROBLEM" not in out


def test_run_calculation_file_directly(tmp_path):
    for f in ("lif_sa4_casscf.out", "lif_sa4_casscf.mld"):
        shutil.copy(os.path.join(DATA, "molpro", f), tmp_path / f)
    prefix = str(tmp_path / "s2")
    inputfile.main(["run", str(tmp_path / "lif_sa4_casscf.out"), "--state", "2", "--qmax", "1",
                    "--points", "3", "--out", prefix])
    arr = np.loadtxt(prefix + "_scattering.csv", delimiter=",", skiprows=2)
    ref = th.load(str(tmp_path / "lif_sa4_casscf.out"), state=2).scattering(qmax=1, points=3)
    np.testing.assert_allclose(arr[:, 1], ref.total, rtol=1e-10)


def test_slurm_setup_lines(tmp_path):
    path = tmp_path / "x.toml"
    path.write_text('[input]\nfile = "a.out"\n[cluster]\ncpus = 4\nsetup = ["module load python"]\n')
    s = inputfile.slurm_script(str(path))
    assert "module load python" in s and f"cd {tmp_path}" in s and "OMP_NUM_THREADS=4" in s


def test_electron_interface(water_cc, water_scf):
    m = th.load(water_cc)
    r = m.scattering(q=[0.0, 1.0, 3.0], probe="electron", parts="all")
    assert np.isnan(r.total[0])
    qb = np.array([1.0, 3.0]) * th.BOHR
    d = xs.electron.decompose(m.wf, qb)
    np.testing.assert_allclose(r.total[1:], 4 * d["total"] / qb ** 4 * th.BOHR ** 2, rtol=1e-10)   # Angstrom^2
    rb = m.scattering(q=qb, unit="1/bohr", probe="electron")
    np.testing.assert_allclose(rb.total, 4 * d["total"] / qb ** 4, rtol=1e-10)                       # bohr^2
    np.testing.assert_allclose(r["nucleus-nucleus"] + r["electron-nucleus"] + r["electron-electron"], r.total,
                               rtol=1e-8)
    fast = m.scattering(q=[1.0], probe="electron", energy_kev=3000)
    np.testing.assert_allclose(fast.total, r.total[1:2] * (1 + 3000 / 510.99895) ** 2, rtol=1e-10)
    # HF Molden-like input: elastic only unless hf
    plain = th.Molecule(xs.Wavefunction(m.wf.basis, m.wf.dm1, nelec=10, atom_numbers=m.wf.atom_numbers,
                                        coords=m.wf.coords))
    assert set(plain.scattering(q=[1.0], probe="electron").curves) == {"elastic"}
    assert set(th.load(water_scf).scattering(q=[1.0], probe="e").curves) == {"total", "elastic", "inelastic"}


def test_electron_iam_and_sm(water_scf):
    m = th.load(water_scf)
    r = m.scattering(q=[0.5, 2.0, 4.0], probe="electron")
    ref = m.iam(q=[0.5, 2.0, 4.0], probe="electron", basis="6-31g")
    np.testing.assert_allclose(ref.atomic + ref.molecular, ref.total, rtol=1e-12)
    s = th.sm(r, ref)
    np.testing.assert_allclose(s["sM (total)"], r.q * (r.total - ref.atomic) / ref.atomic)
    np.testing.assert_allclose(s["sM (IAM)"], r.q * ref.molecular / ref.atomic)
    with pytest.raises(ValueError, match="electron IAM"):
        th.sm(r, m.iam(q=[0.5, 2.0, 4.0], basis="6-31g"))


def test_electron_detector_and_wavelength(water_cc):
    assert abs(th.electron_wavelength(100) - 0.037014) < 1e-5          # 100 keV TEM
    assert abs(th.electron_wavelength(3000) - 0.0035693) < 1e-6        # 3 MeV UED
    img = th.load(water_cc).detector(probe="electron", energy_kev=3000, pixels=9, qmax=5.0)
    assert img["total"].shape == (9, 9) and np.isnan(img["total"][4, 4])   # q = 0 in the centre
    assert np.all(img["total"][np.isfinite(img["total"])] > 0)


def test_cli_electrons(tmp_path):
    for f in ("lif_sa4_casscf.out", "lif_sa4_casscf.mld"):
        shutil.copy(os.path.join(DATA, "molpro", f), tmp_path / f)
    prefix = str(tmp_path / "e")
    w = inputfile.run(str(tmp_path / "lif_sa4_casscf.out"), verbose=False, probe="both", qmax=2.0, points=5,
                      prefix=prefix, plots=False)
    assert {"scattering", "electron_scattering"} <= set(w)
    assert os.path.exists(prefix + "_electron_scattering.csv")
    inputfile.main(["run", str(tmp_path / "lif_sa4_casscf.out"), "--electrons", "--energy", "3000",
                    "--qmax", "2", "--points", "3", "--out", prefix + "2"])
    head = open(prefix + "2_electron_scattering.csv").readline()
    assert '"energy_kev": 3000.0' in head and '"probe": "electron"' in head


def test_lif_sum_rules_through_the_interface():
    """The case a tester reported: total(0) = elastic(0) = N^2 for the Molpro LiF wavefunction."""
    m = th.load(os.path.join(DATA, "molpro", "lif_sa4_casscf.out"))
    r = m.scattering(q=[0.0, 4.0])
    assert r.total[0] == pytest.approx(144.0, rel=1e-6)
    assert r.elastic[0] == pytest.approx(144.0, rel=1e-9)
    assert abs(r.inelastic[0]) < 1e-5
    assert r.total[1] == pytest.approx(34.12804, rel=1e-5) and r.elastic[1] == pytest.approx(28.01261, rel=1e-5)
