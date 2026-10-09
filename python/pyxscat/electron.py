# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Electron scattering (ultrafast electron diffraction, gas electron diffraction).

In the first Born approximation a fast electron is scattered by the whole charge distribution
of the molecule.  With Z(q) = sum_A Z_A exp(i q.R_A) the nuclear and F(q) the electronic form
factor, the differential cross section is::

    dsigma/dOmega = (4 / a0^2 q^4) [ <|Z(q) - F(q)|^2> + S_inel(q) ]
                  = (4 / a0^2 q^4) [ I_nn - 2 I_ne + I_ee ],

where I_ee is the total X-ray scattering (N + two-electron term) and::

    I_nn(q) = sum_AB Z_A Z_B j0(q R_AB)                               nucleus-nucleus
    I_ne(q) = sum_A Z_A int rho(r) j0(q |r - R_A|) dr                 electron-nucleus
    elastic   I_nn - 2 I_ne + <|F|^2>
    inelastic I_ee - <|F|^2>  (the same as for X-rays: the nuclei do not scatter inelastically)

The functions here return the bracket [...] ("scattering in electron units", the quantity of the
earlier code, 0 at q = 0 for a neutral molecule); :func:`cross_section` multiplies by 4/q^4.
Relativistic beam electrons: multiply by gamma^2 = (1 + E/mc^2)^2.

The electron-nucleus integrals are analytic (McMurchie-Davidson with the auxiliary functions
(1/h d/dh)^n j0(qh)), in Fortran.  ``l=2`` gives the second Legendre moment about z, for
aligned molecules: I(q, theta_q) = I0 + 5 <P2> P2(cos theta_q) I2, as for X-rays.

q in bohr^-1, coordinates in bohr.
"""
import ctypes

import numpy as np

from ._lib import lib, ERRORS
from . import core

DEFAULT_SCREEN = 1e-14


def _charges(wf, charges):
    if charges is not None:
        return np.asarray(charges, dtype=float)
    core._need(wf.atom_numbers is not None and wf.coords is not None, "atom_numbers and coords")
    return np.asarray(wf.atom_numbers, dtype=float)


def nucleus_nucleus(charges, coords, q, l=0):
    """sum_AB Z_A Z_B j0(q R_AB), or with l=2 -sum_AB Z_A Z_B j2(q R_AB) P2(cos theta_AB) (i^2 = -1)."""
    from .iam import _j2
    z = np.asarray(charges, dtype=float)
    c = np.asarray(coords, dtype=float)
    q = np.atleast_1d(np.asarray(q, dtype=float))
    d = c[:, None, :] - c[None, :, :]
    r = np.linalg.norm(d, axis=2)
    qr = q[None, None, :] * r[:, :, None]
    zz = z[:, None] * z[None, :]
    if l == 0:
        j0 = np.where(qr > 1e-12, np.sin(qr) / np.where(qr > 1e-12, qr, 1.0), 1.0)
        return np.einsum("ab,abq->q", zz, j0)
    if l == 2:
        with np.errstate(invalid="ignore", divide="ignore"):
            p2 = np.where(r > 1e-12, 1.5 * (d[:, :, 2] / np.where(r > 1e-12, r, 1.0)) ** 2 - 0.5, 0.0)
        return -np.einsum("ab,ab,abq->q", zz, p2, _j2(qr))          # i^2 = -1 (see core.lsign)
    raise ValueError("only l = 0 and l = 2 are implemented")


def electron_nucleus(basis, q, dm1, charges, coords, l=0, nthreads=None, screen=DEFAULT_SCREEN):
    """sum_A Z_A sum D_mu,nu int chi_mu chi_nu K_l(q, r - R_A) dr, K_0 = j0, K_2 = j2 P2."""
    q = np.atleast_1d(np.asarray(q, dtype=np.float64))
    if q.ndim != 1 or np.any(q < 0):
        raise ValueError("q must be a one-dimensional array of non-negative values")
    z = np.asfortranarray(np.atleast_1d(np.asarray(charges, dtype=np.float64)))
    r = np.asfortranarray(np.asarray(coords, dtype=np.float64).reshape(-1, 3).T)
    n = basis.nao
    d = np.asarray(dm1, dtype=np.float64)
    if d.shape != (n, n):
        raise ValueError(f"dm1 must have shape {(n, n)}")
    res = np.zeros(q.size, order="F")
    ierr = ctypes.c_int(0)
    lmn = basis.lmn
    lib.pyxscat_electron_nucleus(
        basis.nprim, n, q.size,
        np.asfortranarray(basis.xyz[:, 0]), np.asfortranarray(basis.xyz[:, 1]), np.asfortranarray(basis.xyz[:, 2]),
        basis.ga, np.asfortranarray(lmn[:, 0]), np.asfortranarray(lmn[:, 1]), np.asfortranarray(lmn[:, 2]),
        basis.group, basis.contr, basis.coeff, np.asfortranarray(q),
        z.size, z, r, np.asfortranarray(d), int(l), float(screen), int(nthreads or 0), res, ctypes.byref(ierr))
    if ierr.value:
        raise RuntimeError("PyXSCAT electron-nucleus kernel: " + ERRORS.get(ierr.value, f"error {ierr.value}"))
    return core.lsign(l) * res


def decompose(wf, q, l=0, hf=False, charges=None, nthreads=None, **cutoffs):
    """All parts of the electron scattering of a wavefunction (bracket units, see module doc).

    Returns a dict over q with ``nn``, ``ne``, ``ee_elastic`` (= X-ray elastic) and
    ``elastic``; with a 2RDM (or ``hf=True`` for a closed-shell determinant) also ``ee_total``
    (= X-ray total), ``inelastic`` and ``total``.
    """
    q = np.atleast_1d(np.asarray(q, dtype=float))
    z = _charges(wf, charges)
    b = wf.basis
    nn = nucleus_nucleus(z, wf.coords, q, l=l)
    ne = electron_nucleus(b, q, wf.dm1, z, wf.coords, l=l, nthreads=nthreads)
    ee_el = core.elastic(b, q, wf.dm1, nthreads=nthreads, l=l, **cutoffs)
    out = {"q": q, "nn": nn, "ne": ne, "ee_elastic": ee_el, "elastic": nn - 2.0 * ne + ee_el}
    if wf.dm2 is not None:
        ee_tot = core.scattering(wf, q, "total", nthreads=nthreads, l=l, **cutoffs)
    elif hf:
        ee_tot = core.scattering(wf, q, "hf_total", nthreads=nthreads, l=l, **cutoffs)
    else:
        return out
    out["ee_total"] = ee_tot
    out["inelastic"] = ee_tot - ee_el
    out["total"] = out["elastic"] + out["inelastic"]
    return out


def scattering(wf, q, kind="total", l=0, hf=False, charges=None, nthreads=None, **cutoffs):
    """Electron scattering of kind ``total``, ``elastic`` or ``inelastic`` (bracket units)."""
    d = decompose(wf, q, l=l, hf=hf, charges=charges, nthreads=nthreads, **cutoffs)
    if kind not in d:
        raise ValueError(f"{kind!r} needs a 2RDM (or hf=True for a closed-shell determinant)"
                         if kind in ("total", "inelastic") else f"unknown kind {kind!r}")
    return d[kind]


def cross_section(bracket, q):
    """dsigma/dOmega = 4 [...] / q^4 (bohr^2 per steradian with q in bohr^-1); NaN at q = 0."""
    q = np.asarray(q, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(q > 0, 4.0 * np.asarray(bracket) / np.where(q > 0, q, 1.0) ** 4, np.nan)


def aligned(wf, qvec, hf=False, charges=None, **kw):
    """Fixed orientation: ``|Z(q) - F(q)|^2`` (+ inelastic) for qvec (nq, 3) in the molecular frame."""
    from . import aligned as al
    z = _charges(wf, charges)
    qvec = np.atleast_2d(np.asarray(qvec, dtype=float))
    zq = np.exp(1j * qvec @ np.asarray(wf.coords, dtype=float).T) @ z
    fq = al.form_factor(wf.basis, qvec, wf.dm1, **kw)
    amp = zq - fq
    out = {"qvec": qvec, "elastic": (amp * amp.conj()).real, "ee_elastic": (fq * fq.conj()).real}
    if wf.dm2 is not None:
        ee_tot = al.scattering(wf, qvec, "total", **kw)
    elif hf:
        ee_tot = al.scattering(wf, qvec, "hf_total", **kw)
    else:
        return out
    out["inelastic"] = ee_tot - out["ee_elastic"]
    out["total"] = out["elastic"] + out["inelastic"]
    return out
