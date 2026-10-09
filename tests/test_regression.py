"""Regression against the pre-existing production code (see tests/data/README.md)."""
import glob
import os

import numpy as np
import pytest

import pyxscat as xs
from pyxscat.interfaces import pyscf as xp
from conftest import get_case

DATA = sorted(glob.glob(os.path.join(os.path.dirname(__file__), "data", "legacy_*.npz")))
# the legacy kernel uses a single-precision pi in the pair prefactors
PI_FACTOR = (float(np.arccos(np.float32(-1.0))) / np.pi) ** 3


@pytest.mark.parametrize("path", DATA, ids=[os.path.basename(p)[7:-4] for p in DATA])
def test_matches_legacy_kernel(path):
    ref = np.load(path)
    case = get_case(str(ref["atom"]), str(ref["basis"]), True)
    q = ref["q"]
    opts = dict(screen_tol=0.0, cutoffcentre=0.01, cutoffz=1e-20, cutoffmd=1e-20)

    hf = xs.mean_field(xp.basis_from_mole(case.mol), q, case.mf.make_rdm1(), **opts)
    assert np.allclose(hf * PI_FACTOR, ref["hf"], rtol=1e-9, atol=1e-9)

    wf = xp.from_cc(case.cc)
    ccsd = xs.two_electron(wf.basis, q, wf.dm2, **opts)
    assert np.allclose(ccsd * PI_FACTOR, ref["ccsd"], rtol=1e-8, atol=1e-8)
