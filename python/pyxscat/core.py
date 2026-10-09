# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Scattering from a basis set and density matrices.

Nothing in this module depends on an electronic structure package.  All density matrices are
expressed in the contracted functions of the :class:`~pyxscat.basis.Basis` (atomic units; q in
inverse bohr).  Conventions:

* ``dm1[mu, nu]`` is the spin-summed one-particle density matrix,
  rho(r) = sum dm1[mu, nu] chi_mu(r) chi_nu(r).
* ``dm2[mu, nu, la, si]`` is the spin-summed two-particle density matrix in chemists' notation,
  Gamma(r1, r2) = sum dm2[mu, nu, la, si] chi_mu(r1) chi_nu(r1) chi_la(r2) chi_si(r2),
  normalised to N(N-1).  This is the convention of PySCF's ``make_rdm2``.
"""
import ctypes
from dataclasses import dataclass
from typing import Optional

import numpy as np

from ._lib import lib, ERRORS
from . import _lib
from .basis import Basis

DEFAULTS = dict(cutoffcentre=1e-10, cutoffz=1e-20, cutoffmd=1e-20, screen_tol=1e-12)
Q_FLOOR = 1e-10      # the Bessel recurrences are evaluated at q >= Q_FLOOR


FACTORED_ABOVE = 60   # basis functions: beyond this correlated 2RDMs are kept factored by default


def use_factored(factored, nao):
    """Resolve a ``factored`` option ("auto", True, False) for a basis of ``nao`` functions."""
    if factored == "auto":
        return nao > FACTORED_ABOVE
    return bool(factored)


class FactoredDM2:
    """A two-particle density given as a sum of products of one-particle matrices,

        Gamma[mu,nu,la,si] = sum_k Ac_k[mu,nu] Bc_k[la,si]  +  sum_k Ax_k[mu,la] Bx_k[nu,si]

    (Coulomb and exchange pairings).  The kernels symmetrise over the 8 index permutations, so
    only contractions with symmetric kernels are defined.  Memory 2 K n^2 instead of n^4: the
    CCSD density of :func:`pyxscat.ccsd.factored_rdm2` has K = 7 nocc^2 + 4 terms.  Weights
    are folded into the A matrices.  ``transform(c)`` changes the basis of every matrix
    (M -> c M c^T); ``dense()`` builds the four-index array (small systems only)."""

    def __init__(self, ac=None, bc=None, ax=None, bx=None):
        self.ac, self.bc = self._arr(ac), self._arr(bc)
        self.ax, self.bx = self._arr(ax), self._arr(bx)
        n = {a.shape[1] for a in (self.ac, self.bc, self.ax, self.bx) if a.size}
        if len(n) != 1 or self.ac.shape != self.bc.shape or self.ax.shape != self.bx.shape:
            raise ValueError("Ac/Bc and Ax/Bx must be (K, n, n) arrays over the same n")
        self.nao = n.pop()

    @staticmethod
    def _arr(a):
        if a is None:
            return np.zeros((0, 0, 0))
        a = np.asarray(a, dtype=float)
        return a.reshape(-1, a.shape[-2], a.shape[-1]) if a.size else np.zeros((0, 0, 0))

    @property
    def shape(self):
        return (self.nao,) * 4

    @property
    def nterms(self):
        return (len(self.ac), len(self.ax))

    @property
    def nbytes(self):
        return sum(a.nbytes for a in (self.ac, self.bc, self.ax, self.bx))

    def transform(self, c):
        """Every matrix M -> c M c^T (e.g. molecular orbitals -> atomic orbitals), in chunks so
        that the temporaries stay small."""
        c = np.asarray(c, dtype=float)

        def f(a):
            if not a.size:
                return None
            out = np.empty((a.shape[0], c.shape[0], c.shape[0]))
            step = max(1, int(2e8 // (8 * c.shape[0] * max(a.shape[1], c.shape[0]))))
            for k0 in range(0, a.shape[0], step):
                t = np.matmul(c, a[k0:k0 + step])                 # (chunk, p, j)
                out[k0:k0 + step] = np.matmul(t, c.T)
            return out
        return FactoredDM2(f(self.ac), f(self.bc), f(self.ax), f(self.bx))

    def dense(self):
        """The four-index array (chemists' notation, 8-fold symmetrised)."""
        n = self.nao
        g = np.einsum("kab,kcd->abcd", self.ac, self.bc) if self.ac.size else np.zeros((n,) * 4)
        if self.ax.size:
            g = g + np.einsum("kac,kbd->abcd", self.ax, self.bx)
        return symmetrize_dm2(g)

    def trace(self):
        """sum_{pq} Gamma[p,p,q,q] (= N(N-1) for a normalised density)."""
        t = float(np.einsum("kaa,kcc->", self.ac, self.bc)) if self.ac.size else 0.0
        if self.ax.size:
            t += float(np.einsum("kac,kac->", self.ax, self.bx))
        return t

    def _kernel_arrays(self):
        """(nkc, Ac, Bc, nkx, Ax, Bx) in the layout the Fortran kernel reads (k fastest); the
        Coulomb matrices symmetrised (allowed: the kernel is symmetric within each pair)."""
        kf = lambda a: np.asfortranarray(a) if a.size else None
        ac = 0.5 * (self.ac + self.ac.transpose(0, 2, 1))
        bc = 0.5 * (self.bc + self.bc.transpose(0, 2, 1))
        return len(self.ac), kf(ac), kf(bc), len(self.ax), kf(self.ax), kf(self.bx)


@dataclass
class Wavefunction:
    """Everything the scattering kernels need: a basis, densities and the electron count.
    ``dm2`` is the four-index array or a :class:`FactoredDM2`."""

    basis: Basis
    dm1: np.ndarray
    dm2: Optional[object] = None
    nelec: Optional[float] = None
    label: str = ""
    atom_numbers: Optional[np.ndarray] = None      # nuclear charges, for the IAM and electron scattering
    coords: Optional[np.ndarray] = None            # nuclear positions (bohr)


def symmetrize_dm2(dm2, inplace=False):
    """Average a two-particle density matrix over its 8 index permutations.

    The integrals are invariant under these permutations, so the scattering is unchanged; the
    kernel loops over symmetry-unique quartets and needs the symmetric form.
    """
    g = dm2 if inplace else np.array(dm2, dtype=np.float64)
    g += g.transpose(1, 0, 2, 3).copy()
    g += g.transpose(0, 1, 3, 2).copy()
    g += g.transpose(2, 3, 0, 1).copy()
    g *= 0.125
    return g


def transform_dm2(g, c):
    """g[pqrs] -> sum_pqrs c[a,p] c[b,q] c[c,r] c[d,s] g[pqrs] (e.g. MO -> AO), one index at a time
    (a :class:`FactoredDM2` transforms its matrices instead)."""
    if isinstance(g, FactoredDM2):
        return g.transform(c)
    g = np.tensordot(c, g, axes=(1, 0))
    g = np.tensordot(c, g, axes=(1, 1))
    g = np.tensordot(c, g, axes=(1, 2))
    g = np.tensordot(c, g, axes=(1, 3))
    return np.ascontiguousarray(g.transpose(3, 2, 1, 0))


def lsign(l):
    """i^l for even l: the phase of the l-th term of the plane-wave expansion
    exp(i q.r) = sum_l (2l+1) i^l j_l(qr) P_l(q^.r^).  Every l-moment returned by pyxscat includes
    it, so S2 = -sum Gamma <j2 P2>, and the moments add up with plus signs:
    I(q, theta_q) = sum_l (2l+1) <P_l(cos beta)> P_l(cos theta_q) S_l(q)."""
    if l % 2:
        raise ValueError("odd l is not implemented")
    return -1.0 if (l // 2) % 2 else 1.0


def _run(basis, q, mode, dm2, da, db, wc, wx, nthreads, return_stats, l=0, factors=None, **cut):
    opts = dict(DEFAULTS)
    unknown = set(cut) - set(opts)
    if unknown:
        raise TypeError(f"unknown options: {sorted(unknown)}")
    opts.update(cut)
    q = np.atleast_1d(np.asarray(q, dtype=np.float64))
    if q.ndim != 1 or np.any(q < 0):
        raise ValueError("q must be a one-dimensional array of non-negative values")
    qk = np.asfortranarray(np.maximum(q, Q_FLOOR))
    n = basis.nao
    res = np.zeros(q.size, order="F")
    stats = np.zeros(4, order="F")
    ierr = ctypes.c_int(0)
    if dm2 is not None:
        if dm2.shape != (n, n, n, n):
            raise ValueError(f"dm2 must have shape {(n,) * 4}")
        dm2 = np.ascontiguousarray(dm2, dtype=np.float64)   # 8-fold symmetric: memory order is irrelevant
        ptr = dm2.ctypes.data_as(ctypes.c_void_p)
    else:
        ptr = None
    nkc = nkx = 0
    pc = [None] * 4
    if factors is not None:
        if factors.nao != n:
            raise ValueError(f"the factored dm2 is over {factors.nao} functions, the basis has {n}")
        nkc, ac, bc, nkx, ax, bx = factors._kernel_arrays()
        keep = (ac, bc, ax, bx)                                  # alive during the call
        pc = [a.ctypes.data_as(ctypes.c_void_p) if a is not None else None for a in keep]
    if da is None:
        da_f = db_f = np.zeros((n, n), order="F")       # not referenced by the explicit-dm2 provider
    else:
        if da.shape != (n, n) or db.shape != (n, n):
            raise ValueError(f"one-particle matrices must have shape {(n, n)}")
        da_f = np.asfortranarray(0.5 * (da + da.T), dtype=np.float64)
        db_f = np.asfortranarray(0.5 * (db + db.T), dtype=np.float64)
    lmn = basis.lmn
    lib.pyxscat_two_electron_j0(
        basis.nprim, n, q.size,
        np.asfortranarray(basis.xyz[:, 0]), np.asfortranarray(basis.xyz[:, 1]), np.asfortranarray(basis.xyz[:, 2]),
        basis.ga,
        np.asfortranarray(lmn[:, 0]), np.asfortranarray(lmn[:, 1]), np.asfortranarray(lmn[:, 2]),
        basis.group, basis.contr, basis.coeff, qk,
        mode, ptr, da_f, db_f, wc, wx, nkc, pc[0], pc[1], nkx, pc[2], pc[3], _lib.DGEMM if factors is not None else None,
        opts["cutoffcentre"], opts["cutoffz"], opts["cutoffmd"], opts["screen_tol"],
        int(l), int(nthreads or 0), res, stats, ctypes.byref(ierr))
    if ierr.value != 0:
        raise RuntimeError("PyXSCAT kernel: " + ERRORS.get(ierr.value, f"error {ierr.value}"))
    res = lsign(l) * res
    if return_stats:
        keys = ("shell_quartets", "shell_quartets_computed", "primitive_quartets", "primitive_quartets_computed")
        return res, dict(zip(keys, (int(s) for s in stats)))
    return res


def two_electron(basis, q, dm2, symmetrize=True, nthreads=None, return_stats=False, l=0, **cutoffs):
    """Two-electron part of the total scattering, sum Gamma <ij|j0(q r12)|kr>.

    Tends to N(N-1) at q = 0 and to zero at large q.  Set ``symmetrize=False`` only if ``dm2``
    already has the full 8-fold permutational symmetry (saves one copy of the array).

    ``l=2`` gives the second Legendre moment about the molecular z axis,
    sum Gamma <ij| j2(q r12) P2(cos theta_12) |kr> (see :func:`scattering`).
    """
    if isinstance(dm2, FactoredDM2):              # sum of products: no four-index array at all
        return _run(basis, q, 3, None, None, None, 0.0, 0.0, nthreads, return_stats, l=l, factors=dm2, **cutoffs)
    g = symmetrize_dm2(dm2) if symmetrize else dm2
    return _run(basis, q, 1, g, None, None, 0.0, 0.0, nthreads, return_stats, l=l, **cutoffs)


def elastic(basis, q, dm1, dm1_b=None, nthreads=None, return_stats=False, l=0, method="pairs", **cutoffs):
    """Rotationally averaged elastic scattering, <|F(q)|^2>, from a one-particle density matrix.

    With two different matrices this is the rotational average of F_a(q) F_b(q)^*, the building
    block of coherent mixed and inelastic terms from transition densities.

    method
        ``"pairs"`` (default): the elastic pair density is a product, so the density is folded
        into every primitive pair first and pairs of Hermite charge distributions are combined
        (``xs_elastic.f90``).  Exact at coincident pair centres, rigorous q-dependent screening.
        ``"quartets"``: the general shell-quartet kernel with the pair density built on the fly
        (the one-loop algorithm; kept for comparison).
    """
    db = dm1 if dm1_b is None else dm1_b
    if method == "quartets":
        return _run(basis, q, 2, None, np.asarray(dm1), np.asarray(db), 1.0, 0.0, nthreads, return_stats, l=l,
                    **cutoffs)
    if method != "pairs":
        raise ValueError("method must be 'pairs' or 'quartets'")
    unknown = set(cutoffs) - set(DEFAULTS)
    if unknown:
        raise TypeError(f"unknown options: {sorted(unknown)}")
    screen = cutoffs.get("screen_tol", DEFAULTS["screen_tol"])
    q = np.atleast_1d(np.asarray(q, dtype=np.float64))
    if q.ndim != 1 or np.any(q < 0):
        raise ValueError("q must be a one-dimensional array of non-negative values")
    n = basis.nao
    da, dbb = np.asarray(dm1, dtype=np.float64), np.asarray(db, dtype=np.float64)
    if da.shape != (n, n) or dbb.shape != (n, n):
        raise ValueError(f"one-particle matrices must have shape {(n, n)}")
    order = np.argsort(q, kind="stable")             # ascending q makes the q screening effective
    qs = np.asfortranarray(q[order])
    res = np.zeros(q.size, order="F")
    stats = np.zeros(2, order="F")
    ierr = ctypes.c_int(0)
    lmn = basis.lmn
    lib.pyxscat_elastic_pairs(
        basis.nprim, n, q.size,
        np.asfortranarray(basis.xyz[:, 0]), np.asfortranarray(basis.xyz[:, 1]), np.asfortranarray(basis.xyz[:, 2]),
        basis.ga, np.asfortranarray(lmn[:, 0]), np.asfortranarray(lmn[:, 1]), np.asfortranarray(lmn[:, 2]),
        basis.group, basis.contr, basis.coeff, qs,
        np.asfortranarray(0.5 * (da + da.T)), np.asfortranarray(0.5 * (dbb + dbb.T)),
        int(l), float(screen), int(nthreads or 0), res, stats, ctypes.byref(ierr))
    if ierr.value != 0:
        raise RuntimeError("PyXSCAT elastic kernel: " + ERRORS.get(ierr.value, f"error {ierr.value}"))
    out = np.empty_like(res)
    out[order] = lsign(l) * res
    if return_stats:
        return out, {"pair_pairs": int(stats[0]), "pair_pairs_computed": int(stats[1])}
    return out


def exchange(basis, q, dm1, nthreads=None, return_stats=False, l=0, **cutoffs):
    """Exchange part of the two-electron scattering of a closed-shell single determinant,
    built from -1/2 D(1,2) D(2,1) without forming a four-index array."""
    d = np.asarray(dm1)
    return _run(basis, q, 2, None, d, d, 0.0, -0.25, nthreads, return_stats, l=l, **cutoffs)


def mean_field(basis, q, dm1, nthreads=None, return_stats=False, l=0, **cutoffs):
    """Two-electron scattering of the mean-field pair density D D - 1/2 D D (exact for
    closed-shell Hartree-Fock), without forming a four-index array."""
    d = np.asarray(dm1)
    return _run(basis, q, 2, None, d, d, 1.0, -0.25, nthreads, return_stats, l=l, **cutoffs)


def scattering(wf, q, kind="total", nthreads=None, l=0, **cutoffs):
    """Rotationally averaged X-ray scattering of a :class:`Wavefunction`.

    kind
        ``"total"``         N + two-electron term (needs ``wf.dm2`` and ``wf.nelec``)
        ``"elastic"``       <|F(q)|^2> from ``wf.dm1``
        ``"inelastic"``     total - elastic
        ``"two_electron"``  the two-electron term alone
        ``"hf_total"``      N + mean-field two-electron term from ``wf.dm1`` only

    l
        ``0`` (default): the isotropic signal, the average over all orientations.
        ``2``: the second Legendre moment about the molecular z axis,
        S2(q) = -sum Gamma <j2(q r12) P2(cos theta_12)> (the i^2 = -1 of the plane-wave expansion
        included, see :func:`lsign`; S2 = <I(q) P2(cos theta_q)> for perfect alignment).  For molecules whose z axis is aligned
        along a laboratory axis (free rotation about it), the signal at a scattering vector
        making the angle theta_q with that axis is

            I(q, theta_q) = I0(q) + 5 <P2(cos beta)> P2(cos theta_q) S2(q) + (l = 4, ...)

        with <P2(cos beta)> the degree of alignment (1 for perfect alignment).  The isotropic
        one-electron term N does not contribute at l = 2.
    """
    b = wf.basis
    self_term = 0.0 if l else None
    if kind == "elastic":
        return elastic(b, q, wf.dm1, nthreads=nthreads, l=l, **cutoffs)
    if kind == "hf_total":
        _need(wf.nelec is not None, "nelec")
        n = 0.0 if l else wf.nelec
        return n + mean_field(b, q, wf.dm1, nthreads=nthreads, l=l, **cutoffs)
    if kind not in ("total", "inelastic", "two_electron"):
        raise ValueError(f"unknown kind {kind!r}")
    _need(wf.dm2 is not None, "dm2")
    i2 = two_electron(b, q, wf.dm2, nthreads=nthreads, l=l, **cutoffs)
    if kind == "two_electron":
        return i2
    _need(wf.nelec is not None, "nelec")
    n = wf.nelec if self_term is None else self_term
    if kind == "total":
        return n + i2
    return n + i2 - elastic(b, q, wf.dm1, nthreads=nthreads, l=l, **cutoffs)


def decompose(wf, q, nthreads=None, l=0, **cutoffs):
    """Split the total scattering into mean-field and correlation contributions.

    Gamma = D D  -  1/2 D D (exchange)  +  Lambda (cumulant), so that

        total = N + hartree + exchange + cumulant,     elastic = hartree.

    Returns a dict of arrays over q.  The mean-field parts are built from the 1RDM on the fly;
    only the full two-electron term uses the four-index array.  ``l=2``: the same split of the
    second Legendre moment (see :func:`scattering`); its self term is zero.
    """
    _need(wf.dm2 is not None, "dm2")
    _need(wf.nelec is not None, "nelec")
    b = wf.basis
    hartree = elastic(b, q, wf.dm1, nthreads=nthreads, l=l, **cutoffs)
    exch = exchange(b, q, wf.dm1, nthreads=nthreads, l=l, **cutoffs)
    i2 = two_electron(b, q, wf.dm2, nthreads=nthreads, l=l, **cutoffs)
    q = np.atleast_1d(np.asarray(q, dtype=float))
    n = 0.0 if l else float(wf.nelec)
    return {
        "q": q,
        "self": np.full(q.shape, n),
        "hartree": hartree,
        "exchange": exch,
        "cumulant": i2 - hartree - exch,
        "two_electron": i2,
        "elastic": hartree,
        "total": n + i2,
        "inelastic": n + i2 - hartree,
    }


def _need(ok, what):
    if not ok:
        raise ValueError(f"this calculation needs `{what}` in the Wavefunction")
