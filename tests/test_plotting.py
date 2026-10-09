"""Plotters for every calculation, reading back the files THISTLE writes, the run report and
the generated plot script."""
import json
import os
import subprocess
import sys
import warnings

import numpy as np
import pytest

import thistle as th
import thistle.plotting as tp

HERE = os.path.dirname(__file__)


@pytest.fixture(scope="module")
def water():
    from pyscf import gto, scf
    mol = gto.M(atom="O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", basis="6-31g", verbose=0)
    return th.load(scf.RHF(mol).run())


def _png(path):
    assert os.path.getsize(path) > 5000
    with open(path, "rb") as f:
        assert f.read(8) == b"\x89PNG\r\n\x1a\n"


def test_every_plotter(water, tmp_path):
    q = dict(qmax=4, points=21)
    res = water.scattering(parts="all", **q)
    iam = water.iam(basis="6-31g", **q)
    e = water.scattering(probe="electron", **q)
    ie = water.iam(probe="electron", basis="6-31g", **q)
    pat = water.detector(grid="polar", energy_kev=10, rings=6, azimuths=12)
    ipat = water.detector(grid="polar", energy_kev=10, rings=6, azimuths=12, model="iam", iam_basis="6-31g")
    img = water.detector(pixels=15)
    figs = {
        "scattering": tp.scattering(res, path=tmp_path / "a.png"),
        "versus_iam": tp.versus_iam(res, iam, path=tmp_path / "b.png"),
        "electron_vs_iam": tp.versus_iam(e, ie, path=tmp_path / "c.png"),
        "sm": tp.sm(th.sm(e, ie), path=tmp_path / "d.png"),
        "anisotropy": tp.anisotropy(water.anisotropy(**q), path=tmp_path / "e.png"),
        "detector_polar": tp.detector(pat, path=tmp_path / "f.png"),
        "detector_square": tp.detector(img, path=tmp_path / "g.png"),
        "detector_vs_iam": tp.detector_versus_iam(pat, ipat, path=tmp_path / "h.png"),
        "overlay": tp.overlay([res, iam], "elastic", labels=["HF", "IAM"], path=tmp_path / "i.png"),
        "difference": tp.difference(th.compare(res, iam), path=tmp_path / "j.png"),
        "aligned": tp.detector(water.aligned_ensemble(p2=0.5, rings=5, azimuths=8, energy_kev=10),
                               path=tmp_path / "k.png"),
    }
    for f in "abcdefghijk":
        _png(tmp_path / f"{f}.png")
    assert len(figs["versus_iam"].axes) == 2 and len(figs["detector_vs_iam"].axes) >= 3
    # the same quantity has the same colour everywhere; identical curves are drawn once
    labels = [ln.get_label() for ln in figs["scattering"].axes[0].get_lines()]
    assert "elastic = hartree" in labels
    total = [ln for ln in figs["scattering"].axes[0].get_lines() if ln.get_label() == "total"][0]
    assert total.get_color() == tp.colour("total") == tp.colour("S2 total") == tp.colour("%diff total")
    # the methods of the result objects use the same plotters
    assert water.scattering(**q).plot(str(tmp_path / "m.png")) is not None
    _png(tmp_path / "m.png")


def test_files_read_back(water, tmp_path):
    res = water.scattering(qmax=3, points=11, parts="all")
    back = tp.read(res.save(str(tmp_path / "r.csv")))
    assert back.unit == res.unit and back.name == res.name and set(back.curves) == set(res.curves)
    for k in res.curves:
        assert np.allclose(back[k], res[k], rtol=1e-11)
    pat = water.detector(grid="polar", energy_kev=10, rings=4, azimuths=6)
    for ext in (".dat", ".npz"):
        b = tp.read(pat.save(str(tmp_path / f"p{ext}")))
        assert b.unit == pat.unit and b.meta["rings"] == 4
        assert np.allclose(b.qvec, pat.qvec, atol=1e-11) and np.allclose(b["total"], pat["total"], rtol=1e-12)
    img = water.detector(pixels=9, unit="1/bohr")
    b = tp.read(img.save(str(tmp_path / "i.npz")))
    assert b.unit == "1/bohr" and np.allclose(b.qx, img.qx) and np.allclose(b["total"], img["total"], equal_nan=True)
    for f in ("r.csv", "p.dat", "i.npz"):
        tp.plot_file(str(tmp_path / f))
        _png(tmp_path / (os.path.splitext(f)[0] + ".png"))


def test_run_report_files_and_plot_script(tmp_path, capsys):
    import shutil
    src = os.path.join(HERE, "data", "molpro")
    for f in os.listdir(src):
        if f.startswith("lif_sa4_casscf"):
            shutil.copy(os.path.join(src, f), tmp_path)
    inp = tmp_path / "lif.toml"
    inp.write_text('''[input]
file = "lif_sa4_casscf.out"
[q]
max = 3.0
points = 7
[compute]
iam = true
iam_basis = "cc-pvdz"
anisotropy = true
[detector]
grid = "polar"
photon_kev = 10
rings = 3
azimuths = 6
[output]
prefix = "out/lif"
''')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        written = th.run(str(inp))
    text = capsys.readouterr().out
    assert "Finished in" in text and "Files written in" in text
    for p in written.values():
        assert os.path.exists(p), p
        assert os.path.basename(p) in text                      # every file is named in the report
    info = json.load(open(written["run_log"]))
    assert info["wall_time_s"] > 0 and any(s["step"].startswith("X-ray scattering") for s in info["steps"])
    assert set(info["files"]) >= {"scattering", "iam", "vs_iam", "anisotropy", "detector", "plot_script"}
    # the generated script redraws every plot from the files alone
    pngs = [p for p in written.values() if p.endswith(".png")]
    for p in pngs:
        os.remove(p)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(sys.path))
    r = subprocess.run([sys.executable, written["plot_script"]], capture_output=True, text=True, env=env,
                       cwd=str(tmp_path.parent))
    assert r.returncode == 0, r.stderr
    for p in pngs:
        _png(p)


def test_plot_command(water, tmp_path):
    from thistle.inputfile import main
    a = water.scattering(qmax=3, points=9).save(str(tmp_path / "a.csv"))
    b = water.iam(qmax=3, points=9, basis="6-31g").save(str(tmp_path / "b.csv"))
    main(["plot", a, "--iam", b, "-o", str(tmp_path / "v.png")])
    main(["plot", a, b, "--curve", "elastic", "-o", str(tmp_path / "o.png")])
    main(["plot", a, "--log"])
    for f in ("v.png", "o.png", "a.png"):
        _png(tmp_path / f)
