"""The polar detector grid of grid.f90 (qcoord3), fixed-orientation patterns on it, and partially
aligned ensembles from the j2 moment."""
import os
import warnings

import numpy as np
import pytest

from pyxscat import grids, aligned as al

HERE = os.path.dirname(__file__)
QGRID = os.path.join(HERE, "data", "grid", "qgrid_10keV_60deg_49_73.dat")


def test_qcoord3_matches_grid_f90():
    ref = grids.read_qgrid(QGRID)
    q = grids.qcoord3(10.0, 60.0, 49, 73)
    assert q.shape == ref.shape == (49 * 72 + 1, 3)
    assert np.max(np.abs(q - ref)) < 1e-11                    # the file has 12 decimals


def test_write_qgrid_is_byte_identical(tmp_path):
    out = grids.write_qgrid(tmp_path / "qgrid.dat", grids.qcoord3(10.0, 60.0, 49, 73))
    with open(out, "rb") as a, open(QGRID, "rb") as b:
        assert a.read() == b.read()


def test_polar_grid_geometry():
    k0 = grids.photon_k0(10.0)
    q = grids.polar_grid(k0, 60.0, 5, 8)
    assert np.allclose(q[0], 0.0)
    kout = np.array([0, 0, k0]) - q                        # q = k_in - k_out: all on the Ewald sphere
    assert np.allclose(np.linalg.norm(kout, axis=1), k0, rtol=0, atol=1e-12)
    two_theta = np.degrees(np.arccos(kout[-1, 2] / k0))
    assert two_theta == pytest.approx(60.0)                # outer ring at the largest angle
    qn, uq, inv = grids.unique_norms(q)
    assert len(uq) == 6 and np.allclose(uq[inv], qn)        # one |q| per ring (+ the origin)


def test_electron_k0():
    import thistle as th
    lam = th.electron_wavelength(3000.0)                   # Angstrom
    assert grids.electron_k0(3000.0) == pytest.approx(2 * np.pi / (lam / th.BOHR), rel=1e-12)


@pytest.fixture(scope="module")
def water():
    from pyscf import gto, scf
    import thistle as th
    mol = gto.M(atom="O 0.05 0.02 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", basis="6-31g", verbose=0)
    return th.load(scf.RHF(mol).run(conv_tol=1e-11))


def test_polar_detector_equals_fixed_orientation(water, tmp_path):
    euler = (10.0, 40.0, 70.0)
    pat = water.detector(grid="polar", energy_kev=10, rings=6, azimuths=8, orientation=euler)
    qlab = grids.polar_grid(grids.photon_k0(10), 60.0, 6, 8)
    assert np.allclose(pat.qvec, qlab)
    qm = al.rotate_q(qlab, al.rotation_matrix(*np.radians(euler)))
    wf = water.wf
    assert np.allclose(pat["elastic"], np.abs(al.form_factor(wf.basis, qm, wf.dm1)) ** 2, rtol=1e-12)
    assert pat["total"][0] == pytest.approx(wf.nelec ** 2, rel=1e-10)
    # a grid file gives the same pattern
    grids.write_qgrid(tmp_path / "g.dat", qlab)
    pat2 = water.detector(grid=str(tmp_path / "g.dat"), orientation=euler)
    assert np.allclose(pat2["total"], pat["total"], rtol=1e-10)
    # text and npz output
    pat.save(tmp_path / "p.dat")
    data = np.loadtxt(tmp_path / "p.dat")
    assert data.shape == (len(qlab), 4 + len(pat.values))
    assert np.allclose(data[:, :3], qlab, atol=1e-11)
    pat.save(tmp_path / "p.npz")


def _orientation_average(wf, qlab, axis, p2, order=59, ngamma=24, electron=False):
    """Explicit average over orientations with the distribution 1 + 5 p2 P2(cos beta) of the
    molecular z axis about ``axis`` (exactly <P2> = p2, nothing beyond l = 2)."""
    from scipy.integrate import lebedev_rule
    from pyxscat import electron as xe
    xyz, w = lebedev_rule(order)
    w = w / w.sum()
    acc = 0.0
    for n, wn in zip(xyz.T, w):
        beta, alpha = np.arccos(np.clip(n[2], -1, 1)), np.arctan2(n[1], n[0])
        f = 1 + 5 * p2 * (1.5 * (n @ axis) ** 2 - 0.5)
        for g in np.arange(ngamma) * 2 * np.pi / ngamma:
            qm = al.rotate_q(qlab, al.rotation_matrix(alpha, beta, g))
            v = xe.aligned(wf, qm, hf=True)["elastic"] if electron else np.abs(al.form_factor(wf.basis, qm, wf.dm1)) ** 2
            acc = acc + wn / ngamma * f * v
    return acc


def test_aligned_ensemble_equals_orientation_average(water):
    """I0 - 5 <P2> P2 S2 is exact for a distribution with no l > 2 components."""
    pytest.importorskip("scipy.integrate", minversion=None)
    axis = np.array([0.3, -0.5, 0.8])
    axis /= np.linalg.norm(axis)
    p2 = 0.6
    qlab = grids.polar_grid(grids.photon_k0(10), 60.0, 3, 4)
    ens = water.aligned_ensemble(axis=axis, p2=p2, grid=qlab, unit="1/bohr")
    ref = _orientation_average(water.wf, qlab, axis, p2)
    iso = water.scattering(q=np.linalg.norm(qlab, axis=1), unit="1/bohr")["elastic"]
    assert np.max(np.abs(ens["elastic"] - iso) / iso) > 1e-3          # alignment is visible
    assert np.max(np.abs(ens["elastic"] - ref) / ref) < 1e-9


def test_aligned_ensemble_electrons(water):
    from pyxscat import electron as xe
    import thistle as th
    axis = np.array([0.0, 1.0, 0.0])
    p2 = -0.3
    qlab = grids.polar_grid(grids.electron_k0(100.0), 1.2, 2, 3)[1:]
    ens = water.aligned_ensemble(axis=axis, p2=p2, grid=qlab, unit="1/bohr", probe="electron", energy_kev=100)
    ref = _orientation_average(water.wf, qlab, axis, p2, order=59, ngamma=24, electron=True)
    ref = xe.cross_section(ref, np.linalg.norm(qlab, axis=1)) * th.core._gamma2(100)
    assert np.max(np.abs(ens["elastic"] - ref) / ref) < 1e-8


def test_aligned_ensemble_limits(water):
    pat0 = water.aligned_ensemble(p2=0.0, rings=4, azimuths=6, energy_kev=10)
    pat1 = water.aligned_ensemble(p2=1.0, rings=4, azimuths=6, energy_kev=10)
    iso = water.scattering(q=pat0.q, unit="1/A")
    for k in ("total", "elastic"):
        assert np.allclose(pat0[k], iso[k], rtol=1e-12)                 # no alignment: isotropic
    # full alignment changes the pattern but not the forward value N^2
    assert pat1["elastic"][0] == pytest.approx(water.wf.nelec ** 2, rel=1e-10)
    assert not np.allclose(pat1["total"], pat0["total"])
    sq = water.aligned_ensemble(grid="square", pixels=11, p2=0.5)
    assert sq.images["total"].shape == (11, 11)


def test_input_file_detector_and_alignment(tmp_path):
    import shutil
    import thistle as th
    src = os.path.join(HERE, "data", "molpro", "lif_sa4_casscf.out")
    for f in os.listdir(os.path.dirname(src)):
        if f.startswith("lif_sa4_casscf"):
            shutil.copy(os.path.join(os.path.dirname(src), f), tmp_path)
    shutil.copy(QGRID, tmp_path / "qgrid.dat")
    inp = tmp_path / "lif.toml"
    inp.write_text('''[input]
file = "lif_sa4_casscf.out"
[q]
max = 2.0
points = 5
[compute]
scattering = false
iam = false
[detector]
grid = "polar"
photon_kev = 10
rings = 3
azimuths = 4
orientation = [0, 90, 0]
[alignment]
axis = [0, 1, 0]
p2 = 0.4
grid = "qgrid.dat"
[output]
prefix = "out/lif"
plots = false
''')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w = th.run(str(inp), verbose=False)
    det = np.loadtxt(w["detector"])
    assert det.shape[0] == 3 * 4 + 1
    al_ = np.loadtxt(w["aligned"])
    assert al_.shape[0] == 49 * 72 + 1
    assert np.allclose(al_[:, :3], grids.read_qgrid(QGRID), atol=1e-11)


# ----------------------------------------------------------------------------- IAM on the same grids
def test_iam_detector_matches_iam_module(water):
    from pyxscat import iam
    import thistle as th
    euler = (20.0, 50.0, 10.0)
    kw = dict(grid="polar", energy_kev=10, rings=4, azimuths=6, orientation=euler, iam_basis="6-31g")
    pat = water.detector(model="iam", **kw)
    qm = al.rotate_q(pat.qvec, al.rotation_matrix(*np.radians(euler)))
    ref = iam.for_wavefunction(water.wf, qm, aligned=True, basis="6-31g")
    for k in ("total", "elastic", "inelastic"):
        assert np.allclose(pat[k], ref[k], rtol=1e-12, atol=1e-12)
    assert pat.name == "IAM"
    diff = th.compare(water.detector(**{k: v for k, v in kw.items() if k != "iam_basis"}), pat, "total")
    assert abs(diff["%diff total"][0]) < 1e-8                     # both N^2 at q = 0
    assert np.all(np.isfinite(diff["%diff total"]))


def test_electron_iam_fixed_orientation_averages_to_rotational(water):
    """The new fixed-orientation electron IAM, averaged over orientations, is the Debye form."""
    from scipy.integrate import lebedev_rule
    from pyxscat import iam
    wf = water.wf
    q = np.array([0.7, 2.5])
    xyz, w = lebedev_rule(59)
    w = w / w.sum()
    avg = {k: np.array([w @ iam.electron_for_wavefunction(wf, qq * xyz.T, aligned=True, basis="6-31g")[k]
                        for qq in q]) for k in ("total", "elastic")}
    rot = iam.electron_for_wavefunction(wf, q, basis="6-31g")
    for k in ("total", "elastic"):
        assert np.allclose(avg[k], rot[k], rtol=1e-10)


@pytest.mark.parametrize("probe", ["xray", "electron"])
def test_iam_aligned_ensemble_equals_orientation_average(water, probe, monkeypatch):
    from scipy.integrate import lebedev_rule
    from pyxscat import iam, electron as xe
    import thistle as th
    # the free-atom terms depend on |q| only, the same for every orientation: compute them once
    memo, terms = {}, iam._atom_terms

    def cached(symbols, q, *args):
        key = (tuple(symbols), np.asarray(q).tobytes(), repr(args))
        if key not in memo:
            memo[key] = terms(symbols, q, *args)
        return memo[key]
    monkeypatch.setattr(iam, "_atom_terms", cached)
    axis = np.array([0.6, 0.0, 0.8])
    p2 = 0.45
    k0 = grids.photon_k0(10) if probe == "xray" else grids.electron_k0(100.0)
    qlab = grids.polar_grid(k0, 60.0 if probe == "xray" else 1.5, 3, 5)[1:]
    ens = water.aligned_ensemble(axis=axis, p2=p2, grid=qlab, unit="1/bohr", model="iam", iam_basis="6-31g",
                                 probe=probe, energy_kev=100 if probe == "electron" else None)
    fn = iam.electron_for_wavefunction if probe == "electron" else iam.for_wavefunction
    xyz, w = lebedev_rule(59)
    w = w / w.sum()
    ng = 24
    acc = 0.0
    for n, wn in zip(xyz.T, w):
        beta, alpha = np.arccos(np.clip(n[2], -1, 1)), np.arctan2(n[1], n[0])
        f = 1 + 5 * p2 * (1.5 * (n @ axis) ** 2 - 0.5)
        for g in np.arange(ng) * 2 * np.pi / ng:
            qm = al.rotate_q(qlab, al.rotation_matrix(alpha, beta, g))
            acc = acc + wn / ng * f * fn(water.wf, qm, aligned=True, basis="6-31g")["total"]
    if probe == "electron":
        acc = xe.cross_section(acc, np.linalg.norm(qlab, axis=1)) * th.core._gamma2(100)
    assert np.max(np.abs(ens["total"] - acc) / np.abs(acc)) < 1e-9


def test_input_file_writes_iam_patterns(tmp_path):
    import shutil
    import thistle as th
    src = os.path.join(HERE, "data", "molpro", "lif_sa4_casscf.out")
    for f in os.listdir(os.path.dirname(src)):
        if f.startswith("lif_sa4_casscf"):
            shutil.copy(os.path.join(os.path.dirname(src), f), tmp_path)
    inp = tmp_path / "lif.toml"
    inp.write_text('''[input]
file = "lif_sa4_casscf.out"
[q]
max = 2.0
points = 3
[compute]
scattering = false
iam = true
iam_basis = "cc-pvdz"
[detector]
grid = "polar"
photon_kev = 10
rings = 2
azimuths = 4
orientation = [90, 90, 0]
[alignment]
axis = [0, 1, 0]
p2 = 0.4
rings = 3
[output]
prefix = "out/lif"
plots = true
''')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w = th.run(str(inp), verbose=False)
    for k in ("detector", "detector_iam", "detector_vs_iam", "aligned", "aligned_iam", "aligned_vs_iam",
              "detector_vs_iam_plot", "aligned_vs_iam_plot"):
        assert os.path.exists(w[k]), k
    d = np.loadtxt(w["aligned_vs_iam"])
    assert d.shape == (3 * 4 + 1, 5) and abs(d[0, 4]) < 1e-6
