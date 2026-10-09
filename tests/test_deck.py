# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""thistle-dev input decks (HDF5 schema v1): read in both representations and reproduce
thistle-dev's reference curves.  The decks and curves live in thistle-dev; this test uses
the copies in tests/data/thistle_dev (CH4 STO-3G total, MO and dense AO)."""
import os

import numpy as np
import pytest

import pyxscat as xs
import thistle as th
from pyxscat.interfaces import deck

h5py = pytest.importorskip("h5py")
DATA = os.path.join(os.path.dirname(__file__), "data", "thistle_dev")
PIFAC = 1.0 + 8.3483e-8                 # thistle-dev's references carry the legacy single-precision pi


@pytest.mark.parametrize("name", ["ch4_total", "ch4_total_zcotr"])
def test_deck_reproduces_thistle_dev(name):
    path = os.path.join(DATA, f"{name}.h5")
    wf, info = deck.read(path)
    assert info["rdm_mode"] == ("AO" if "zcotr" in name else "MO")
    assert wf.nelec == pytest.approx(10.0, abs=1e-9)
    assert xs.elastic(wf.basis, [0.0], wf.dm1)[0] == pytest.approx(100.0, rel=1e-10)
    gold = np.loadtxt(os.path.join(DATA, "ch4_total_golden.dat"))
    q = gold[:, 0]
    keep = gold[:, 1] > 1e-8 * gold[:, 1].max()
    # the reference is the two-electron term, computed with cutoffcentre = 0.01
    i2 = xs.two_electron(wf.basis, q, wf.dm2, screen_tol=1e-16, cutoffcentre=info["cutoffcentre"])
    assert np.max(np.abs(i2 * PIFAC - gold[:, 1])[keep] / gold[keep, 1]) < 1e-11
    # and our exact default stays within the size of that approximation
    exact = xs.two_electron(wf.basis, q, wf.dm2)
    assert np.max(np.abs(exact * PIFAC - gold[:, 1])[keep] / gold[keep, 1]) < 5e-5


def test_thistle_load_deck():
    m = th.load(os.path.join(DATA, "ch4_total.h5"))
    assert m.correlated
    r = m.scattering(qmax=2.0, points=3, unit="1/bohr")
    assert r["total"][0] == pytest.approx(100.0, rel=1e-9)
    ao = th.load(os.path.join(DATA, "ch4_total_zcotr.h5"))
    assert np.allclose(ao.scattering(qmax=2.0, points=3, unit="1/bohr")["total"], r["total"], rtol=1e-9)


def test_writer_round_trip(tmp_path):
    """Our deck writer reproduces thistle-dev's own file (every dataset), and a deck written
    from a factored CASSCF density reads back to the same scattering."""
    import h5py
    from pyxscat.interfaces import deck, molpro
    src = os.path.join(DATA, "ch4_total_zcotr.h5")
    wf, info = deck.read(src)
    out = str(tmp_path / "rt.h5")
    deck.write(out, wf, q=info["q"])
    with h5py.File(src) as a, h5py.File(out) as b:
        for k in ("basis/contr", "basis/group", "basis/l", "basis/m", "basis/n", "basis/x", "basis/y",
                  "basis/z", "basis/exponent", "contraction/coeffs", "contraction/gs", "contraction/gf",
                  "contraction/gc", "molecule/coords", "molecule/atomic_numbers", "rdm/gamma_ao"):
            assert a[k].shape == b[k].shape and np.abs(a[k][()] - b[k][()]).max() < 1e-14, k
        assert dict(a["params"].attrs) == dict(b["params"].attrs)
        assert b.attrs["rdm_mode"] == "AO" and list(b.attrs["states"]) == [1, 1]
        # orbitals over primitives = contraction coefficient x contracted orbitals (as theirs)
        mm, c = b["mos/mmod"][()], b["mos2/C"][()].T
        cf, ct = b["contraction/coeffs"][()], b["basis/contr"][()]
        assert np.abs(mm - cf[:, None] * c[ct - 1, :]).max() < 1e-14
    d = os.path.join(os.path.dirname(__file__), "data", "molpro")
    w = molpro.from_multi(os.path.join(d, "lif_sa4_casscf.out"), os.path.join(d, "lif_sa4_casscf.mld"),
                          state=2, factored=True)
    deck.write(str(tmp_path / "lif.h5"), w)
    w2, _ = deck.read(str(tmp_path / "lif.h5"))
    q = np.linspace(0.0, 5.0, 6)
    a, b = xs.scattering(w, q, screen_tol=1e-16), xs.scattering(w2, q, screen_tol=1e-16)
    assert np.max(np.abs(a - b) / a) < 1e-10
