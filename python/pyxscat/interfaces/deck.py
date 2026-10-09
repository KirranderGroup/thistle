# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Read the HDF5 "deck" of thistle-dev (KirranderGroup/thistle-dev, schema v1).

A deck holds a primitive Cartesian basis, molecular orbitals over the primitives (``mos/mmod``,
contraction and normalisation included), the geometry, the q grid and cutoffs (``params``), and
the density in one of two representations (attribute ``rdm_mode``):

``MO``  ``rdm/onerdm_mo`` and ``rdm/mo2rdm`` in the orbital basis (chemists' notation, trace
        N(N-1), PySCF's convention).  Here every primitive becomes its own function and the
        densities are transformed with ``mmod``.
``AO``  ``rdm/gamma_ao`` (total) or ``rdm/onerdm_ao`` (elastic) over contracted functions,
        with ``basis/contr`` and ``contraction/coeffs`` defining them: exactly pyxscat's Basis.
        Without a one-particle density (total decks) it follows from the two-particle one,
        D = sum_rs Gamma[:, :, r, s] S_rs / (N - 1), with N(N - 1) = sum Gamma S S.

Arrays were written in Fortran order; read back from Python they come with reversed axes,
which is harmless here: the 1RDM is symmetric and the 2RDM is invariant under reversing all
four indices.  Checked against the five test decks of thistle-dev to 4e-13 (with the deck's
cutoffcentre); see benchmarks/thistle_dev_crosscheck.
"""
import numpy as np

from ..basis import Basis
from ..core import Wavefunction, transform_dm2


def is_deck(path):
    """True for an HDF5 file with the attributes of a thistle-dev deck."""
    try:
        import h5py
        with h5py.File(path, "r") as f:
            return "rdm_mode" in f.attrs and "basis" in f and "mos" in f
    except Exception:
        return False


def read(path):
    """``(Wavefunction, info)`` from a deck; ``info`` holds the q grid, cutoffs, observable,
    states and the representation."""
    import h5py
    with h5py.File(path, "r") as f:
        a = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in f.attrs.items()}
        b = f["basis"]
        xyz = np.column_stack([b["x"][()], b["y"][()], b["z"][()]])
        ga = b["exponent"][()]
        lmn = np.column_stack([b["l"][()], b["m"][()], b["n"][()]])
        group = b["group"][()].astype(int)
        n = ga.size
        mode = str(a["rdm_mode"]).upper()
        zat = np.asarray(f["molecule/atomic_numbers"][()], dtype=float)
        coords = np.asarray(f["molecule/coords"][()], dtype=float).T   # stored (natoms, 3) in Fortran order
        params = {k: (v.item() if hasattr(v, "item") else v) for k, v in f["params"].attrs.items()} \
            if "params" in f else {}

        if mode == "MO":
            basis = Basis(xyz, ga, lmn, group, np.arange(1, n + 1), np.ones(n), n)
            mmod = f["mos/mmod"][()]
            mmod = mmod if mmod.shape[0] == n else mmod.T             # (ngto, norbs)
            d1 = f["rdm/onerdm_mo"][()] if "rdm/onerdm_mo" in f else None
            g = f["rdm/mo2rdm"][()] if "rdm/mo2rdm" in f else None
            dm1 = None if d1 is None else mmod @ d1 @ mmod.T
            dm2 = None if g is None else transform_dm2(g, mmod)
            nelec = None if d1 is None else float(np.trace(d1))
            if dm1 is None and g is not None:                         # 1RDM from the 2RDM
                nn = float(np.einsum("ppqq->", g))
                nelec = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * nn))
                dm1 = mmod @ (np.einsum("pqrr->pq", g) / (nelec - 1.0)) @ mmod.T
        elif mode == "AO":
            contr = b["contr"][()].astype(int)
            coeff = f["contraction/coeffs"][()]
            nao = int(contr.max())
            basis = Basis(xyz, ga, lmn, group, contr, coeff, nao)
            dm1 = f["rdm/onerdm_ao"][()] if "rdm/onerdm_ao" in f else None
            dm2 = f["rdm/gamma_ao"][()] if "rdm/gamma_ao" in f else None
            nelec = None
            if dm2 is not None:
                from .molden import overlap_matrix
                s = overlap_matrix(basis)
                nn = float(np.einsum("pqrs,pq,rs->", dm2, s, s))
                nelec = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * nn))
                if dm1 is None:
                    dm1 = np.einsum("pqrs,rs->pq", dm2, s) / (nelec - 1.0)
            elif dm1 is not None:
                from .molden import overlap_matrix
                nelec = float(np.einsum("pq,qp->", dm1, overlap_matrix(basis)))
        else:
            raise ValueError(f"{path}: unknown rdm_mode {mode!r}")

        q = None
        if {"qmin", "qmax", "npoints"} <= set(params):
            q = np.linspace(float(params["qmin"]), float(params["qmax"]), int(params["npoints"]))
        info = {"rdm_mode": mode, "observable": a.get("observable"), "typec": a.get("typec"),
                "states": np.asarray(a.get("states", [])).tolist(), "q": q,
                "cutoffcentre": params.get("cutoffcentre"), "cutoffz": params.get("cutoffz"),
                "cutoffmd": params.get("cutoffmd"), "file_out": a.get("file_out")}
    wf = Wavefunction(basis, dm1, dm2, nelec=nelec, label=f"thistle-dev deck ({mode})",
                      atom_numbers=zat, coords=coords)
    return wf, info


# ----------------------------------------------------------------------------- writer
def _native(grp, name, a, dtype):
    """Store ``a`` (the engines' native, column-major logical shape) as thistle-dev does: the
    Fortran-order transpose, so that a Fortran/Julia read recovers it index for index."""
    a = np.asarray(a, dtype=dtype)
    return grp.create_dataset(name, data=np.ascontiguousarray(np.asfortranarray(a).T))


def write(path, wf, observable="total", q=None, cutoffcentre=0.01, cutoffz=1e-20, cutoffmd=1e-20,
          states=None, typec=None, file_out="scattering_out.dat", max_gb=8.0):
    """Write a wavefunction as a thistle-dev deck (HDF5 schema v1, AO representation), the
    input of its Fortran and Julia engines.

    observable
        ``"total"`` writes the two-particle density (``rdm/gamma_ao``, chemists' notation, trace
        N(N-1)); ``"elastic"`` the one-particle density (``rdm/onerdm_ao``).
    q
        the q grid for the deck's ``params`` (bohr^-1, evenly spaced; default 1e-10 to 10, 200
        points).  ``cutoffcentre`` defaults to thistle-dev's 0.01 (THISTLE itself uses the exact
        limit), so that both codes compute the same numbers.

    The deck format holds the four-index array: a factored density is expanded here, up to
    ``max_gb`` gigabytes.  The orbitals in ``mos`` are natural orbitals of the 1RDM.
    """
    import h5py
    if observable not in ("total", "elastic"):
        raise ValueError("observable must be 'total' or 'elastic'")
    b = wf.basis
    nao, ngto = b.nao, b.ga.size
    from .molden import overlap_matrix
    s = overlap_matrix(b)
    rdm = {}
    if observable == "total":
        g = wf.dm2
        if g is None:
            raise ValueError("a total deck needs a two-particle density (a correlated wavefunction)")
        from ..core import FactoredDM2
        if isinstance(g, FactoredDM2):
            need = 8.0 * nao ** 4 / 1e9
            if need > max_gb:
                raise MemoryError(f"the deck needs the {nao}^4 array ({need:.0f} GB > max_gb={max_gb})")
            g = g.dense()
        rdm["gamma_ao"] = np.asarray(g, dtype=float)
    else:
        rdm["onerdm_ao"] = np.asarray(wf.dm1, dtype=float)
    # natural orbitals: S^1/2 D S^1/2 = U n U^T, C = S^-1/2 U (contracted functions x orbitals)
    w, v = np.linalg.eigh(s)
    sh, smh = (v * np.sqrt(w)) @ v.T, (v / np.sqrt(w)) @ v.T
    occ, u = np.linalg.eigh(sh @ wf.dm1 @ sh)
    order = np.argsort(-occ)
    c = smh @ u[:, order]                                     # (ncontr, norbs)
    mmod = b.coeff[:, None] * c[b.contr - 1, :]               # (ngto, norbs): primitives x orbitals
    counts = np.bincount(b.contr - 1, minlength=nao).astype(np.int32)
    gf = np.cumsum(counts).astype(np.int32)
    gs = (gf - counts + 1).astype(np.int32)
    if q is None:
        q = np.linspace(1e-10, 10.0, 200)
    q = np.asarray(q, dtype=float)
    if q.size > 2 and not np.allclose(np.diff(q), q[1] - q[0], rtol=1e-8, atol=1e-12):
        raise ValueError("a deck stores qmin, qmax and npoints: q must be evenly spaced")
    if typec is None:
        typec = 1 if observable == "total" else 2
    states = [1, 1] if states is None else list(states)
    with h5py.File(path, "w", libver=("earliest", "v110")) as f:
        f.attrs["deck_format"] = 1
        f.attrs["rdm_mode"] = "AO"
        f.attrs["observable"] = observable
        f.attrs["typec"] = int(typec)
        f.attrs["states"] = np.asarray(states, dtype=np.int32)
        f.attrs["file_out"] = np.bytes_(file_out)
        f.attrs["source_var"] = "readtwordm"
        f.attrs["source_rdm_files"] = [b"THISTLE"]
        mol = f.create_group("molecule")
        z = wf.atom_numbers if wf.atom_numbers is not None else np.zeros(0)
        _native(mol, "atomic_numbers", z, np.float64)
        _native(mol, "coords", wf.coords if wf.coords is not None else np.zeros((0, 3)), np.float64)
        prm = f.create_group("params")
        for k, v in (("cutoffcentre", cutoffcentre), ("cutoffz", cutoffz), ("cutoffmd", cutoffmd),
                     ("qmin", q[0]), ("qmax", q[-1])):
            prm.attrs[k] = float(v)
        prm.attrs["npoints"] = int(q.size)
        bas = f.create_group("basis")
        for k, a in (("x", b.xyz[:, 0]), ("y", b.xyz[:, 1]), ("z", b.xyz[:, 2]), ("exponent", b.ga)):
            _native(bas, k, a, np.float64)
        for k, a in (("l", b.lmn[:, 0]), ("m", b.lmn[:, 1]), ("n", b.lmn[:, 2]), ("group", b.group),
                     ("contr", b.contr)):
            _native(bas, k, a, np.int32)
        mos = f.create_group("mos")
        _native(mos, "mmod", mmod.T, np.float64)              # native (norbs, ngto)
        mos.attrs["norbs"] = int(nao)
        mos.attrs["ngto"] = int(ngto)
        con = f.create_group("contraction")
        con.attrs["ncontr"] = int(nao)
        _native(con, "coeffs", b.coeff, np.float64)
        for k, a in (("gs", gs), ("gf", gf), ("gc", counts)):
            _native(con, k, a, np.int32)
        m2 = f.create_group("mos2")
        _native(m2, "C", c, np.float64)                       # native (ncontr, norbs)
        r = f.create_group("rdm")
        for k, a in rdm.items():
            _native(r, k, a, np.float64)
    return path
