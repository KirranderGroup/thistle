# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Independent atom model (IAM), for X-rays and electrons: the reference against which molecular (bonding, correlation)
effects are measured.

The molecule is modelled as free, spherically averaged atoms at the nuclear positions:

    elastic     I_el(q)  = sum_AB f_A(q) f_B(q) j0(q R_AB)                 (rotational average)
                         = |sum_A f_A(q) exp(i q.R_A)|^2                   (fixed orientation)
    total       I_tot(q) = sum_A S_A(q) + sum_{A != B} f_A f_B j0(q R_AB)   (or cos(q.R_AB))
    inelastic   I_tot - I_el = sum_A [S_A(q) - f_A(q)^2]

f_A is the form factor of the spherically averaged atomic density and S_A = N_A + I2_A the total
(elastic + inelastic) scattering of the free atom.  Two sources:

* ``source="ab initio"`` (default): free-atom calculations (``basis`` is a PySCF basis name, or
  a dict element -> basis name) with PySCF at a chosen level
  (UHF by default, or CCSD for closed-shell atoms) and basis, through the same kernels as the
  molecule.  Elastic and total.  With the molecule's own method and basis, molecule minus IAM
  isolates the effect of bonding.
* ``source="tables"``: Waasmaier-Kirfel form factors (through the ``xraydb`` package).  Elastic
  only; the tables give no incoherent scattering.

q is in inverse bohr throughout.
"""
from functools import lru_cache

import numpy as np

from . import aligned as _aligned
from . import core

ELEMENTS = ["X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
            "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
            "Br", "Kr"]
BOHR_IN_ANGSTROM = 0.52917721092


def _symbol(z_or_symbol):
    if isinstance(z_or_symbol, str):
        return z_or_symbol.capitalize()
    return ELEMENTS[int(z_or_symbol)]


# ----------------------------------------------------------------------------- free atoms
_GROUND_STATE_SPIN = {"H": 1, "He": 0, "Li": 1, "Be": 0, "B": 1, "C": 2, "N": 3, "O": 2, "F": 1, "Ne": 0,
                      "Na": 1, "Mg": 0, "Al": 1, "Si": 2, "P": 3, "S": 2, "Cl": 1, "Ar": 0}


@lru_cache(maxsize=None)
def _atom_wavefunction(symbol, basis, method):
    from pyscf import gto, scf
    from .interfaces import pyscf as xp
    spin = _GROUND_STATE_SPIN.get(symbol)
    if spin is None:
        raise ValueError(f"no ground-state multiplicity tabulated for {symbol}; pass source='tables'")
    mol = gto.M(atom=f"{symbol} 0 0 0", basis=basis, spin=spin, verbose=0, symmetry=False)
    if method == "hf":
        mf = scf.UHF(mol).run(conv_tol=1e-11) if spin else scf.RHF(mol).run(conv_tol=1e-11)
        if spin:
            da, db = mf.make_rdm1()
            wf = xp.from_density(mol, da + db)
            sa, sb = xp.from_density(mol, da).dm1, xp.from_density(mol, db).dm1
            return wf, (sa, sb)
        return xp.from_scf(mf), None
    if method == "ccsd":
        if spin:
            raise NotImplementedError("CCSD atoms are implemented for closed shells only; use method='hf'")
        from pyscf import cc
        mf = scf.RHF(mol).run(conv_tol=1e-11)
        return xp.from_cc(cc.CCSD(mf).run(conv_tol=1e-10)), None
    raise ValueError("method must be 'hf' or 'ccsd'")


def _sphere(order=17):
    from scipy.integrate import lebedev_rule
    xyz, w = lebedev_rule(order)
    return xyz.T, w / w.sum()


def atomic_form_factor(symbol, q, source="ab initio", basis="aug-cc-pvtz", method="hf"):
    """Form factor f(q) of the spherically averaged free atom (q in bohr^-1)."""
    symbol = _symbol(symbol)
    q = np.atleast_1d(np.asarray(q, dtype=float))
    if source == "tables":
        import xraydb
        s = q / BOHR_IN_ANGSTROM / (4 * np.pi)            # sin(theta)/lambda in Angstrom^-1
        return np.asarray(xraydb.f0(symbol, s), dtype=float)
    wf, _ = _atom_wavefunction(symbol, basis, method)
    # <F(q n)> over directions n = spherical average of the density, exactly (Lebedev grid of
    # degree 17 integrates the angular structure of densities up to g functions)
    n, w = _sphere()
    out = np.empty(q.size)
    for i, qq in enumerate(q):
        out[i] = (w @ _aligned.form_factor(wf.basis, qq * n, wf.dm1)).real
    return out


def atomic_total(symbol, q, basis="aug-cc-pvtz", method="hf", **cutoffs):
    """Total (elastic + inelastic) rotationally averaged scattering S(q) = N + I2(q) of a free atom."""
    symbol = _symbol(symbol)
    q = np.atleast_1d(np.asarray(q, dtype=float))
    wf, spins = _atom_wavefunction(symbol, basis, method)
    if method == "ccsd":
        return core.scattering(wf, q, "total", **cutoffs)
    if spins is None:
        return wf.nelec + core.mean_field(wf.basis, q, wf.dm1, **cutoffs)
    da, db = spins
    # UHF pair density: D(1)D(2) - Da Da - Db Db  (exchange within each spin)
    i2 = (core.elastic(wf.basis, q, wf.dm1, **cutoffs)
          + 0.5 * core.exchange(wf.basis, q, 2 * da, **cutoffs) + 0.5 * core.exchange(wf.basis, q, 2 * db, **cutoffs))
    return wf.nelec + i2


def _atom_terms(symbols, q, source, basis, method, need_total, cutoffs):
    uniq = sorted(set(symbols))
    bas = {s: (basis[s] if isinstance(basis, dict) else basis) for s in uniq}
    f = {s: atomic_form_factor(s, q, source, bas[s], method) for s in uniq}
    tot = None
    if need_total:
        if source == "tables":
            raise ValueError("the tables give no inelastic scattering: use source='ab initio' for total/inelastic")
        tot = {s: atomic_total(s, q, bas[s], method, **cutoffs) for s in uniq}
    return f, tot


# ----------------------------------------------------------------------------- molecules
def scattering(atom_numbers, coords, q, source="ab initio", basis="aug-cc-pvtz", method="hf", total=True,
               l=0, **cutoffs):
    """Rotationally averaged IAM (Debye formula).  ``coords`` in bohr.

    Returns a dict over q: elastic, and with ``total=True`` also total, inelastic and the atomic
    and interference parts.

    ``l=2``: the second Legendre moment about the z axis (the IAM counterpart of
    ``pyxscat.scattering(..., l=2)``),

        S2(q) = -sum_{A != B} f_A f_B j2(q R_AB) P2(cos theta_AB)   (i^2 = -1, see core.lsign),

    where theta_AB is the angle of R_A - R_B with z.  Free atoms are spherical, so their own
    elastic and inelastic scattering is isotropic: atomic and inelastic parts vanish at l = 2.
    """
    if l == 2:
        return _scattering_l2(atom_numbers, coords, q, source, basis, method, total)
    if l != 0:
        raise ValueError("only l = 0 and l = 2 are implemented")
    symbols = [_symbol(z) for z in atom_numbers]
    coords = np.asarray(coords, dtype=float)
    q = np.atleast_1d(np.asarray(q, dtype=float))
    f, tot = _atom_terms(symbols, q, source, basis, method, total, cutoffs)
    fa = np.array([f[s] for s in symbols])                         # (natom, nq)
    r = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
    qr = q[None, None, :] * r[:, :, None]
    j0 = np.where(qr > 1e-12, np.sin(qr) / np.where(qr > 1e-12, qr, 1.0), 1.0)
    elastic = np.einsum("aq,bq,abq->q", fa, fa, j0)
    out = {"q": q, "elastic": elastic}
    if total:
        atomic = sum(tot[s] for s in symbols)
        self_el = (fa ** 2).sum(axis=0)
        out["atomic"] = atomic
        out["interference"] = elastic - self_el
        out["total"] = atomic + elastic - self_el
        out["inelastic"] = out["total"] - elastic
    return out


def _j2(x):
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-3
    xs = x[small]
    out[small] = xs ** 2 / 15.0 - xs ** 4 / 210.0
    xl = x[~small]
    out[~small] = (3.0 / xl ** 3 - 1.0 / xl) * np.sin(xl) - 3.0 * np.cos(xl) / xl ** 2
    return out


def _scattering_l2(atom_numbers, coords, q, source, basis, method, total):
    symbols = [_symbol(z) for z in atom_numbers]
    coords = np.asarray(coords, dtype=float)
    q = np.atleast_1d(np.asarray(q, dtype=float))
    f, _ = _atom_terms(symbols, q, source, basis, method, False, {})
    fa = np.array([f[s] for s in symbols])
    d = coords[:, None, :] - coords[None, :, :]
    r = np.linalg.norm(d, axis=2)
    with np.errstate(invalid="ignore", divide="ignore"):
        p2 = np.where(r > 1e-12, 1.5 * (d[:, :, 2] / np.where(r > 1e-12, r, 1.0)) ** 2 - 0.5, 0.0)
    j2 = _j2(q[None, None, :] * r[:, :, None])
    elastic = -np.einsum("aq,bq,ab,abq->q", fa, fa, p2, j2)        # i^2 = -1 (core.lsign)
    out = {"q": q, "elastic": elastic}
    if total:
        zero = np.zeros_like(q)
        out.update(atomic=zero, interference=elastic, total=elastic.copy(), inelastic=zero)
    return out


def electron_scattering(atom_numbers, coords, q, source="ab initio", basis="aug-cc-pvtz", method="hf",
                        total=True, l=0, **cutoffs):
    """Electron scattering of independent atoms (rotational average), in the bracket units of
    :mod:`pyxscat.electron` (multiply by 4/q^4 for the cross section).

    Each atom scatters with the Mott-Bethe amplitude Z_A - f_A(q):

        elastic    sum_AB (Z_A - f_A)(Z_B - f_B) j0(q R_AB)
        atomic     sum_A  (Z_A - f_A)^2  (+ inelastic sum_A [S_A - f_A^2] with ``total``)
        molecular  elastic - sum_A (Z_A - f_A)^2      (the interference, as in sM(s))

    ``l=2``: -sum_{A != B} (Z_A - f_A)(Z_B - f_B) j2(q R_AB) P2(cos theta_AB)  (i^2 = -1).
    """
    symbols = [_symbol(z) for z in atom_numbers]
    zat = np.array([ELEMENTS.index(s) for s in symbols], dtype=float)
    coords = np.asarray(coords, dtype=float)
    q = np.atleast_1d(np.asarray(q, dtype=float))
    f, tot = _atom_terms(symbols, q, source, basis, method, total and l == 0, cutoffs)
    amp = zat[:, None] - np.array([f[s] for s in symbols])        # (natom, nq)
    d = coords[:, None, :] - coords[None, :, :]
    r = np.linalg.norm(d, axis=2)
    qr = q[None, None, :] * r[:, :, None]
    if l == 2:
        with np.errstate(invalid="ignore", divide="ignore"):
            p2 = np.where(r > 1e-12, 1.5 * (d[:, :, 2] / np.where(r > 1e-12, r, 1.0)) ** 2 - 0.5, 0.0)
        el = -np.einsum("aq,bq,ab,abq->q", amp, amp, p2, _j2(qr))      # i^2 = -1 (core.lsign)
        zero = np.zeros_like(q)
        out = {"q": q, "elastic": el, "molecular": el, "atomic": zero}
        if total:
            out.update(total=el.copy(), inelastic=zero)
        return out
    if l != 0:
        raise ValueError("only l = 0 and l = 2 are implemented")
    j0 = np.where(qr > 1e-12, np.sin(qr) / np.where(qr > 1e-12, qr, 1.0), 1.0)
    el = np.einsum("aq,bq,abq->q", amp, amp, j0)
    self_el = (amp ** 2).sum(axis=0)
    out = {"q": q, "elastic": el, "molecular": el - self_el, "atomic": self_el}
    if total:
        inel = sum(tot[s] - f[s] ** 2 for s in symbols)
        out["inelastic"] = inel
        out["total"] = el + inel
        out["atomic"] = self_el + inel
    return out


def electron_for_wavefunction(wf, q, aligned=False, **kw):
    """Electron IAM for the geometry of a :class:`~pyxscat.Wavefunction` (``aligned=True``: q
    vectors in the molecular frame)."""
    if wf.atom_numbers is None or wf.coords is None:
        raise ValueError("this Wavefunction carries no nuclear positions; call iam.electron_scattering directly")
    fn = electron_scattering_aligned if aligned else electron_scattering
    return fn(wf.atom_numbers, wf.coords, q, **kw)


def scattering_aligned(atom_numbers, coords, qvec, source="ab initio", basis="aug-cc-pvtz", method="hf",
                       total=True, **cutoffs):
    """IAM at fixed orientation; ``qvec`` (nq, 3) in the molecular frame (bohr^-1)."""
    symbols = [_symbol(z) for z in atom_numbers]
    coords = np.asarray(coords, dtype=float)
    qvec = np.atleast_2d(np.asarray(qvec, dtype=float))
    qn = np.linalg.norm(qvec, axis=1)
    uq, inv = np.unique(np.round(qn, 12), return_inverse=True)
    f, tot = _atom_terms(symbols, uq, source, basis, method, total, cutoffs)
    fa = np.array([f[s][inv] for s in symbols])                    # (natom, nq)
    amp = np.einsum("aq,aq->q", fa, np.exp(1j * coords @ qvec.T))
    elastic = (amp * amp.conj()).real
    out = {"qvec": qvec, "elastic": elastic}
    if total:
        atomic = sum(tot[s][inv] for s in symbols)
        self_el = (fa ** 2).sum(axis=0)
        out["atomic"] = atomic
        out["interference"] = elastic - self_el
        out["total"] = atomic + elastic - self_el
        out["inelastic"] = out["total"] - elastic
    return out


def electron_scattering_aligned(atom_numbers, coords, qvec, source="ab initio", basis="aug-cc-pvtz",
                                method="hf", total=True, **cutoffs):
    """Electron IAM at fixed orientation (bracket units; multiply by 4/q^4 for the cross section):
    elastic |sum_A (Z_A - f_A) exp(i q.R_A)|^2, plus the atomic inelastic sum_A [S_A - f_A^2]."""
    symbols = [_symbol(z) for z in atom_numbers]
    zat = np.array([ELEMENTS.index(s) for s in symbols], dtype=float)
    coords = np.asarray(coords, dtype=float)
    qvec = np.atleast_2d(np.asarray(qvec, dtype=float))
    qn = np.linalg.norm(qvec, axis=1)
    uq, inv = np.unique(np.round(qn, 12), return_inverse=True)
    f, tot = _atom_terms(symbols, uq, source, basis, method, total, cutoffs)
    amp = zat[:, None] - np.array([f[s][inv] for s in symbols])     # (natom, nq)
    phase = np.einsum("aq,aq->q", amp, np.exp(1j * coords @ qvec.T))
    el = (phase * phase.conj()).real
    self_el = (amp ** 2).sum(axis=0)
    out = {"qvec": qvec, "elastic": el, "molecular": el - self_el, "atomic": self_el}
    if total:
        inel = sum(tot[s][inv] - f[s][inv] ** 2 for s in symbols)
        out.update(inelastic=inel, total=el + inel, atomic=self_el + inel)
    return out


def for_wavefunction(wf, q, aligned=False, **kw):
    """IAM for the geometry of a :class:`~pyxscat.Wavefunction` (atoms taken from the adapter).
    ``l=2`` gives the second Legendre moment (rotational form only)."""
    if wf.atom_numbers is None or wf.coords is None:
        raise ValueError("this Wavefunction carries no nuclear positions; call iam.scattering directly")
    fn = scattering_aligned if aligned else scattering
    return fn(wf.atom_numbers, wf.coords, q, **kw)


def percent_difference(molecular, reference):
    """100 (I_molecule - I_IAM) / I_IAM, the usual way of showing molecular effects."""
    return 100.0 * (np.asarray(molecular) - np.asarray(reference)) / np.asarray(reference)
