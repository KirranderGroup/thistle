# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Scattering from molecules with a fixed orientation (aligned or oriented molecules, single
detector images), at arbitrary scattering vectors q = (qx, qy, qz) in the molecular frame.

Because exp(i q.(r1 - r2)) = exp(i q.r1) exp(-i q.r2), every quantity follows from the Fourier
transforms of basis-function products::

    F_mu,nu(q) = int chi_mu(r) chi_nu(r) exp(i q.r) dr,

computed analytically (``pair_transform``) and contracted with the density matrices::

    form factor     F(q)   = sum D_mu,nu F_mu,nu(q)                  elastic = |F(q)|^2
    two-electron    I2(q)  = sum Gamma_{mu nu, la si} F_mu,nu F_la,si^*
    total           I(q)   = N + I2(q)

The cost is one pass over primitive pairs plus dense linear algebra; there is no Bessel
expansion, so close pair centres need no special treatment.

Detector geometry: ``ewald_q`` gives the scattering vectors for given scattering and azimuthal
angles, and ``rotate_q`` turns laboratory-frame vectors into the molecular frame for a given
molecular orientation.
"""
import ctypes

import numpy as np

from ._lib import lib, ERRORS
from .core import _need

DEFAULT_SCREEN = 1e-14


def pair_transform(basis, qvec, screen=DEFAULT_SCREEN, nthreads=None):
    """F[q, mu, nu] = int chi_mu chi_nu exp(i q.r) dr for qvec of shape (nq, 3) (bohr^-1)."""
    qvec = np.atleast_2d(np.asarray(qvec, dtype=np.float64))
    if qvec.shape[1] != 3:
        raise ValueError("qvec must have shape (nq, 3)")
    nq, n = qvec.shape[0], basis.nao
    f = np.zeros((nq, n, n), dtype=np.complex128, order="F")
    ierr = ctypes.c_int(0)
    lmn = basis.lmn
    lib.pyxscat_pair_transform(
        basis.nprim, n, nq,
        np.asfortranarray(basis.xyz[:, 0]), np.asfortranarray(basis.xyz[:, 1]), np.asfortranarray(basis.xyz[:, 2]),
        basis.ga, np.asfortranarray(lmn[:, 0]), np.asfortranarray(lmn[:, 1]), np.asfortranarray(lmn[:, 2]),
        basis.group, basis.contr, basis.coeff, np.asfortranarray(qvec.T), float(screen), int(nthreads or 0),
        f, ctypes.byref(ierr))
    if ierr.value:
        raise RuntimeError("PyXSCAT pair transform: " + ERRORS.get(ierr.value, f"error {ierr.value}"))
    return f


def _chunks(nq, n, max_bytes=5e8):
    step = max(1, int(max_bytes // (16 * n * n)))
    for i in range(0, nq, step):
        yield slice(i, min(nq, i + step))


def form_factor(basis, qvec, dm1, **kw):
    """Complex form factor F(q) = sum D F_mu,nu(q) (exp(+i q.r) convention)."""
    qvec = np.atleast_2d(np.asarray(qvec, dtype=float))
    d = np.asarray(dm1, dtype=float)
    out = np.empty(qvec.shape[0], dtype=np.complex128)
    for sl in _chunks(len(qvec), basis.nao):
        f = pair_transform(basis, qvec[sl], **kw)
        out[sl] = np.einsum("qij,ij->q", f, d)
    return out


def elastic(basis, qvec, dm1, dm1_b=None, **kw):
    """``|F(q)|^2``, or ``Re[F_a(q) F_b(q)^*]`` for two one-particle matrices (transition/mixed terms)."""
    fa = form_factor(basis, qvec, dm1, **kw)
    fb = fa if dm1_b is None else form_factor(basis, qvec, dm1_b, **kw)
    return (fa * fb.conj()).real


def _factored_two_electron(f, g):
    """sum Gamma F F^* for a FactoredDM2 and pair transforms f (nq, n, n), F symmetric:
    Coulomb  sum_k (A_k . F)(B_k . F^*);  exchange  sum_k sum_ac A_k[a,c] (F B_k F^*)[a,c],
    averaged with the transposed term (the images under a<->b, c<->d coincide for symmetric F)."""
    out = np.zeros(len(f))
    if g.ac.size:
        pa = np.einsum("qab,kab->qk", f, g.ac)
        pb = np.einsum("qab,kab->qk", f.conj(), g.bc)
        out += 0.5 * (np.einsum("qk,qk->q", pa, pb) + np.einsum("qk,qk->q", pa.conj(), pb.conj())).real
    for a, b in zip(g.ax, g.bx):
        fb = f @ b                                       # (nq, n, n) = F B
        t = np.einsum("qab,qcb,ac->q", fb, f.conj(), a)                 # sum_ac A[a,c] (F B F*^T)[a,c]
        fbt = f @ b.T
        t2 = np.einsum("qab,qcb,ca->q", fbt, f.conj(), a)               # the transposed image
        out += 0.5 * (t + t2).real
    return out


def two_electron(basis, qvec, dm2, **kw):
    """Two-electron term sum Gamma F F^* at fixed q (tends to N(N-1) at q = 0)."""
    from .core import FactoredDM2
    qvec = np.atleast_2d(np.asarray(qvec, dtype=float))
    n = basis.nao
    if isinstance(dm2, FactoredDM2):
        out = np.empty(qvec.shape[0])
        for sl in _chunks(len(qvec), n):
            out[sl] = _factored_two_electron(pair_transform(basis, qvec[sl], **kw), dm2)
        return out
    g = np.asarray(dm2, dtype=float).reshape(n * n, n * n)
    out = np.empty(qvec.shape[0])
    for sl in _chunks(len(qvec), n):
        f = pair_transform(basis, qvec[sl], **kw).reshape(-1, n * n, order="C")
        out[sl] = np.einsum("qa,qa->q", f @ g, f.conj()).real
    return out


def exchange(basis, qvec, dm1, **kw):
    """Exchange term of the mean-field pair density, -1/2 sum D_ps D_rq F_pq F_rs^*."""
    qvec = np.atleast_2d(np.asarray(qvec, dtype=float))
    d = np.asarray(dm1, dtype=float)
    out = np.empty(qvec.shape[0])
    for sl in _chunks(len(qvec), basis.nao):
        f = pair_transform(basis, qvec[sl], **kw)
        # sum_pq F_pq (D F^* D)_pq
        t = np.einsum("ps,qrs->qpr", d, f.conj())          # (D F*)  [q, p, r] = sum_s D_ps F*_rs  (F symmetric)
        t = np.einsum("qpr,rk->qpk", t, d)                  # (D F* D)
        out[sl] = -0.5 * np.einsum("qpk,qpk->q", f, t).real
    return out


def mean_field(basis, qvec, dm1, **kw):
    """Two-electron term of D D - 1/2 D D (closed-shell Hartree-Fock) without a 2RDM."""
    return elastic(basis, qvec, dm1, **kw) + exchange(basis, qvec, dm1, **kw)


def scattering(wf, qvec, kind="total", **kw):
    """Fixed-orientation counterpart of :func:`pyxscat.scattering` (same ``kind`` values)."""
    b = wf.basis
    if kind == "elastic":
        return elastic(b, qvec, wf.dm1, **kw)
    if kind == "hf_total":
        _need(wf.nelec is not None, "nelec")
        return wf.nelec + mean_field(b, qvec, wf.dm1, **kw)
    if kind not in ("total", "inelastic", "two_electron"):
        raise ValueError(f"unknown kind {kind!r}")
    _need(wf.dm2 is not None, "dm2")
    i2 = two_electron(b, qvec, wf.dm2, **kw)
    if kind == "two_electron":
        return i2
    _need(wf.nelec is not None, "nelec")
    if kind == "total":
        return wf.nelec + i2
    return wf.nelec + i2 - elastic(b, qvec, wf.dm1, **kw)


def decompose(wf, qvec, **kw):
    """Fixed-orientation counterpart of :func:`pyxscat.decompose`."""
    _need(wf.dm2 is not None, "dm2")
    _need(wf.nelec is not None, "nelec")
    b = wf.basis
    qvec = np.atleast_2d(np.asarray(qvec, dtype=float))
    n = b.nao
    from .core import FactoredDM2
    d = np.asarray(wf.dm1, dtype=float)
    factored = isinstance(wf.dm2, FactoredDM2)
    g = wf.dm2 if factored else np.asarray(wf.dm2, dtype=float).reshape(n * n, n * n)
    res = {k: np.empty(len(qvec)) for k in ("hartree", "exchange", "two_electron")}
    for sl in _chunks(len(qvec), n):
        f = pair_transform(b, qvec[sl], **kw)
        ff = np.einsum("qij,ij->q", f, d)
        res["hartree"][sl] = (ff * ff.conj()).real
        t = np.einsum("qpr,rk->qpk", np.einsum("ps,qrs->qpr", d, f.conj()), d)
        res["exchange"][sl] = -0.5 * np.einsum("qpk,qpk->q", f, t).real
        if factored:
            res["two_electron"][sl] = _factored_two_electron(f, g)
        else:
            fl = f.reshape(-1, n * n)
            res["two_electron"][sl] = np.einsum("qa,qa->q", fl @ g, fl.conj()).real
    nel = float(wf.nelec)
    res["self"] = np.full(len(qvec), nel)
    res["cumulant"] = res["two_electron"] - res["hartree"] - res["exchange"]
    res["elastic"] = res["hartree"]
    res["total"] = nel + res["two_electron"]
    res["inelastic"] = res["total"] - res["hartree"]
    res["qvec"] = qvec
    return res


# ----------------------------------------------------------------------------- geometry helpers
def ewald_q(k0, two_theta, phi):
    """Scattering vectors q = k_out - k_in on the Ewald sphere, incident beam along +z.

    k0 = 2 pi / wavelength (bohr^-1); two_theta = scattering angle, phi = azimuth (radians,
    broadcast against each other).  Returns an array (..., 3); ``|q| = 2 k0 sin(two_theta / 2)``.
    """
    tt, ph = np.broadcast_arrays(np.asarray(two_theta, float), np.asarray(phi, float))
    return k0 * np.stack([np.sin(tt) * np.cos(ph), np.sin(tt) * np.sin(ph), np.cos(tt) - 1.0], axis=-1)


def detector_q(k0, qmax, npix):
    """A square detector image in q space: (npix x npix) pixels covering ``|q_perp| <= qmax`` on the
    Ewald sphere, beam along +z.  Returns (qvec (npix*npix, 3), mask of pixels with q <= 2 k0)."""
    s = np.linspace(-qmax, qmax, npix)
    qx, qy = np.meshgrid(s, s, indexing="xy")
    qperp = np.hypot(qx, qy)
    ok = qperp < k0
    tt = np.arcsin(np.clip(qperp / k0, 0, 1))       # q_perp = k0 sin(2 theta)
    ph = np.arctan2(qy, qx)
    q = ewald_q(k0, tt, ph)
    return q.reshape(-1, 3), ok.reshape(-1)


def rotation_matrix(alpha, beta, gamma):
    """Rotation R(alpha, beta, gamma) = Rz(alpha) Ry(beta) Rz(gamma) (active, radians)."""
    def rz(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    def ry(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    return rz(alpha) @ ry(beta) @ rz(gamma)


def rotate_q(qvec_lab, rotation):
    """Laboratory-frame q vectors -> molecular frame, for a molecule rotated by ``rotation``
    (3x3, active) from its input orientation: q_mol = R^T q_lab."""
    return np.asarray(qvec_lab) @ np.asarray(rotation)
