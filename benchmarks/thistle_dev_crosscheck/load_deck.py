# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Read a thistle-dev schema-v1 deck.h5 into a pyxscat Wavefunction (each primitive its own
function; the MO coefficients over primitives carry contraction and normalisation)."""
import h5py, numpy as np
from pyxscat.basis import Basis
from pyxscat.core import Wavefunction, transform_dm2

def load(path):
    f = h5py.File(path, "r")
    b = f["basis"]
    x, y, z = b["x"][()], b["y"][()], b["z"][()]
    ga = b["exponent"][()]
    lmn = np.column_stack([b["l"][()], b["m"][()], b["n"][()]])
    grp = b["group"][()].astype(int)
    n = ga.size
    basis = Basis(np.column_stack([x, y, z]), ga, lmn, grp, np.arange(1, n + 1), np.ones(n), n)
    M = f["mos/mmod"][()]                      # stored (ngto, norbs)
    if M.shape[0] != n:
        M = M.T
    attrs = dict(f.attrs)
    if attrs["rdm_mode"] == "MO":
        d1 = f["rdm/onerdm_mo"][()]
        g = f["rdm/mo2rdm"][()]
        dm1 = M @ d1 @ M.T
        dm2 = transform_dm2(g, M)
    else:
        dm1 = dm2 = None
    zat = f["molecule/atomic_numbers"][()]
    xyz = f["molecule/coords"][()]
    xyz = xyz.T if xyz.shape[0] == 3 else xyz
    wf = Wavefunction(basis, dm1, dm2, nelec=float(np.trace(f["rdm/onerdm_mo"][()])) if "rdm/onerdm_mo" in f else None,
                      atom_numbers=zat, coords=xyz)
    return wf, attrs, f
