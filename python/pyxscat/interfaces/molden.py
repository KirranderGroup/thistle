# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Molden files: basis set, geometry and molecular orbitals from any program that writes them
(Molpro, Molcas/OpenMolcas, BAGEL, ORCA, PySCF, Psi4, ...).

The reader does not depend on any electronic structure package.  Spherical (pure) shells are
expressed exactly in Cartesian functions, so the MO coefficients it returns refer to the
Cartesian :class:`~pyxscat.basis.Basis` it builds.

Conventions handled
-------------------
* ``[Atoms] AU`` or ``Angs``; Fortran ``D`` exponents.
* ``[5D]``, ``[5D7F]``, ``[5D10F]``, ``[7F]``, ``[9G]`` flags (default: Cartesian 6D/10F/15G).
* Molden component order: d xx yy zz xy xz yz;  f xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz;
  g xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy;
  spherical m = 0, +1, -1, +2, -2, ...
* Contraction coefficients refer to normalised primitives; the contracted function is
  renormalised.
* Cartesian components of d, f, g shells are normalised either each to one (``"unit"``) or all
  with the factor of the x^l component (``"axis"``).  Programs differ; with the default
  ``cart_norm="auto"`` the convention that makes the occupied MOs orthonormal is chosen.
"""
import re
from dataclasses import dataclass, field
from math import comb, factorial, gamma, pi, sqrt

import numpy as np

from ..basis import Basis, cartesian_powers
from ..core import Wavefunction, transform_dm2

ANGSTROM = 1.0 / 0.52917721092
_LCHAR = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}

_MOLDEN_CART = {
    0: [(0, 0, 0)],
    1: [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    2: [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)],
    3: [(3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0), (2, 0, 1), (1, 0, 2), (0, 1, 2), (0, 2, 1),
        (1, 1, 1)],
    4: [(4, 0, 0), (0, 4, 0), (0, 0, 4), (3, 1, 0), (3, 0, 1), (1, 3, 0), (0, 3, 1), (1, 0, 3), (0, 1, 3),
        (2, 2, 0), (2, 0, 2), (0, 2, 2), (2, 1, 1), (1, 2, 1), (1, 1, 2)],
}


def _dfact(n):
    """Double factorial, (-1)!! = 1."""
    out = 1
    while n > 1:
        out *= n
        n -= 2
    return out


def _gauss_moment(n, p):
    """int x^n exp(-p x^2) dx over the real line."""
    if n % 2:
        return 0.0
    return gamma((n + 1) / 2) / p ** ((n + 1) / 2)


def _same_centre_overlap(pa, pb, a, b):
    """Overlap of x^pa e^{-a r^2} and x^pb e^{-b r^2} (powers are 3-tuples) on the same centre."""
    p = a + b
    return (_gauss_moment(pa[0] + pb[0], p) * _gauss_moment(pa[1] + pb[1], p)
            * _gauss_moment(pa[2] + pb[2], p))


def _prim_norm(a, l):
    """Normalisation of x^l exp(-a r^2)."""
    return (2 * a / pi) ** 0.75 * (4 * a) ** (l / 2) / sqrt(_dfact(2 * l - 1))


def solid_harmonic(l, m):
    """Real solid harmonic of degree l and order m as {(lx, ly, lz): coefficient} (unnormalised,
    Molden/Helgaker sign convention: d+1 = xz, d-1 = yz, d+2 = x^2-y^2, d-2 = xy)."""
    am = abs(m)
    vm = 0.0 if m >= 0 else 0.5
    out = {}
    for t in range((l - am) // 2 + 1):
        for u in range(t + 1):
            v = vm
            while v <= (am / 2 - vm) // 1 + vm + 1e-9:
                c = ((-1) ** int(round(t + v - vm)) * 0.25 ** t * comb(l, t) * comb(l - t, am + t)
                     * comb(t, u) * comb(am, int(round(2 * v))))
                key = (int(round(2 * t + am - 2 * (u + v))), int(round(2 * (u + v))), l - 2 * t - am)
                out[key] = out.get(key, 0.0) + c
                v += 1.0
    return {k: c for k, c in out.items() if c != 0.0}


@dataclass
class MoldenData:
    """Contents of a Molden file, with orbitals expressed in the Cartesian ``basis``."""

    basis: Basis
    atom_numbers: np.ndarray
    coords: np.ndarray                  # bohr
    mo_coeff: np.ndarray                # (nao_cartesian, nmo)
    mo_occ: np.ndarray
    mo_energy: np.ndarray
    mo_spin: np.ndarray                 # "Alpha" / "Beta"
    mo_sym: list
    cart_norm: str
    file_to_cart: np.ndarray = field(repr=False)   # chi_file = chi_cart @ file_to_cart

    @property
    def nelec(self):
        return float(self.mo_occ.sum())

    def density(self, spin=None):
        """One-particle density matrix sum_i occ_i C_i C_i^T (spin-summed unless spin is given)."""
        sel = np.ones(self.mo_occ.size, bool) if spin is None else (self.mo_spin == spin)
        c = self.mo_coeff[:, sel]
        return (c * self.mo_occ[sel]) @ c.T

    def wavefunction(self, dm1_mo=None, dm2_mo=None, orbitals=None, label="Molden"):
        """Build a :class:`~pyxscat.Wavefunction`.

        Without arguments the 1RDM is built from the occupations in the file (enough for
        elastic scattering and for the mean-field total).  Correlated densities in the MO basis
        (e.g. from a CI vector) are passed as ``dm1_mo`` / ``dm2_mo``; ``orbitals`` selects the
        MO columns they refer to (default: the first ``dm1_mo.shape[0]``).
        """
        if dm1_mo is None:
            return Wavefunction(self.basis, self.density(), None, nelec=self.nelec, label=label,
                                atom_numbers=self.atom_numbers, coords=self.coords)
        n = dm1_mo.shape[0]
        c = self.mo_coeff[:, :n] if orbitals is None else self.mo_coeff[:, orbitals]
        dm1 = c @ dm1_mo @ c.T
        dm2 = None
        if dm2_mo is not None:
            dm2 = transform_dm2(dm2_mo, c)
        nel = float(np.trace(dm1_mo))
        if abs(nel - round(nel)) < 1e-6:
            nel = float(round(nel))
        return Wavefunction(self.basis, dm1, dm2, nelec=nel, label=label,
                            atom_numbers=self.atom_numbers, coords=self.coords)


def _float(s):
    return float(s.replace("D", "E").replace("d", "e"))


def _sections(lines):
    out = {}
    current = None
    for line in lines:
        s = line.strip()
        m = re.match(r"^\[([^\]]+)\](.*)$", s)
        if m:
            current = m.group(1).strip().upper()
            out[current] = {"opts": m.group(2).strip(), "lines": []}
        elif current is not None:
            out[current]["lines"].append(line.rstrip("\n"))
    return out


def _parse_atoms(sec):
    unit = ANGSTROM if "ANG" in sec["opts"].upper() else 1.0
    z, xyz = [], []
    for line in sec["lines"]:
        w = line.split()
        if len(w) < 6:
            continue
        z.append(int(w[2]))
        xyz.append([_float(v) * unit for v in w[3:6]])
    return np.array(z), np.array(xyz)


def _parse_gto(sec):
    """Shells as (atom_index (0-based), l, exponents, coefficients)."""
    shells = []
    lines = [l for l in sec["lines"]]
    i = 0
    atom = None
    while i < len(lines):
        w = lines[i].split()
        i += 1
        if not w:
            atom = None
            continue
        if atom is None:
            atom = int(w[0]) - 1
            continue
        kind = w[0].lower()
        nprim = int(w[1])
        scale = _float(w[2]) if len(w) > 2 else 1.0
        rows = [lines[i + k].split() for k in range(nprim)]
        i += nprim
        exps = np.array([_float(r[0]) for r in rows]) * scale ** 2
        if kind == "sp":
            shells.append((atom, 0, exps, np.array([_float(r[1]) for r in rows])))
            shells.append((atom, 1, exps, np.array([_float(r[2]) for r in rows])))
        else:
            if kind not in _LCHAR:
                raise ValueError(f"unsupported shell type {kind!r} in Molden file")
            shells.append((atom, _LCHAR[kind], exps, np.array([_float(r[1]) for r in rows])))
    return shells


def _parse_mo(sec):
    mos = []
    cur = None
    for line in sec["lines"]:
        s = line.strip()
        if not s:
            continue
        if "=" in s:
            key, val = [x.strip() for x in s.split("=", 1)]
            key = key.lower()
            if cur is None or cur["coef"]:
                cur = {"sym": "", "ene": 0.0, "spin": "Alpha", "occ": 0.0, "coef": {}}
                mos.append(cur)
            if key == "sym":
                cur["sym"] = val
            elif key == "ene":
                cur["ene"] = _float(val)
            elif key == "spin":
                cur["spin"] = "Beta" if val.lower().startswith("b") else "Alpha"
            elif key == "occup":
                cur["occ"] = _float(val)
        else:
            w = s.split()
            cur["coef"][int(w[0]) - 1] = _float(w[1])
    return mos


def _molpro_order(mos):
    """Molpro writes ``Sym= n.irrep`` labels and may list the orbitals out of its internal order
    (natural orbitals sorted for output).  CI vectors refer to the internal order, so with C1
    labels the orbitals are put back in label order.  Returns the permutation, or None."""
    labels = [m["sym"].strip() for m in mos]
    if not labels or not all(re.fullmatch(r"\d+\.1", lab) for lab in labels):
        return None
    for spin in ("Alpha", "Beta"):
        nums = [int(m["sym"].split(".")[0]) for m in mos if m["spin"] == spin]
        if nums and sorted(nums) != list(range(1, len(nums) + 1)):
            return None
    order = sorted(range(len(mos)), key=lambda i: (mos[i]["spin"] == "Beta", int(mos[i]["sym"].split(".")[0])))
    return None if order == list(range(len(mos))) else order


def read_molden(path, cart_norm="auto", reorder="auto"):
    """Read a Molden file.  Returns :class:`MoldenData`.

    ``reorder="auto"`` puts Molpro orbitals (``Sym= n.1`` labels) back in Molpro's internal order,
    the order its CI vectors use; ``reorder=False`` keeps the file order.
    """
    with open(path) as f:
        lines = f.readlines()
    secs = _sections(lines)
    if "ATOMS" not in secs or "GTO" not in secs or "MO" not in secs:
        raise ValueError(f"{path}: [Atoms], [GTO] and [MO] sections are required")
    znum, coords = _parse_atoms(secs["ATOMS"])
    shells = _parse_gto(secs["GTO"])
    keys = set(secs)
    pure = {2: False, 3: False, 4: False}
    if "5D" in keys or "5D7F" in keys:
        pure[2] = pure[3] = True
    if "5D10F" in keys:
        pure[2] = True
    if "7F" in keys:
        pure[3] = True
    if "9G" in keys:
        pure[4] = True
    mos = _parse_mo(secs["MO"])
    if reorder == "auto":
        perm = _molpro_order(mos)
        if perm is not None:
            mos = [mos[i] for i in perm]

    if cart_norm == "auto":
        cands = ["unit", "axis"] if any(l >= 2 and not pure.get(l, False) for _, l, _, _ in shells) else ["unit"]
        best = None
        for cn in cands:
            data = _build(znum, coords, shells, pure, mos, cn)
            err = _orthonormality_error(data)
            if best is None or err < best[0]:
                best = (err, data)
        return best[1]
    return _build(znum, coords, shells, pure, mos, cart_norm)


def _build(znum, coords, shells, pure, mos, cart_norm):
    bshells = []          # for Basis.from_shells (Cartesian, one contraction per shell)
    blocks = []           # per shell: (n_file_functions, transformation block ncart x nfile)
    for atom, l, exps, cs in shells:
        powers = cartesian_powers(l)
        # bare-primitive coefficients for the x^l component, normalised contraction
        bare = cs * np.array([_prim_norm(a, l) for a in exps])
        axis = (l, 0, 0)
        s = sum(bi * bj * _same_centre_overlap(axis, axis, ai, aj)
                for ai, bi in zip(exps, bare) for aj, bj in zip(exps, bare))
        bare = bare / sqrt(s)
        if pure.get(l, False):
            comp = np.ones(len(powers))                       # all components share the x^l factor
            bshells.append((coords[atom], l, exps, bare, comp))
            # overlap of the Cartesian monomials of this contracted shell (same centre)
            smono = np.array([[sum(bi * bj * _same_centre_overlap(pa, pb, ai, aj)
                                   for ai, bi in zip(exps, bare) for aj, bj in zip(exps, bare))
                               for pb in powers] for pa in powers])
            ms = [0] + [s * k for k in range(1, l + 1) for s in (1, -1)]
            t = np.zeros((len(powers), 2 * l + 1))
            for col, m in enumerate(ms):
                for p, c in solid_harmonic(l, m).items():
                    t[powers.index(p), col] = c
                t[:, col] /= sqrt(t[:, col] @ smono @ t[:, col])
            blocks.append(t)
        else:
            if cart_norm == "unit":
                comp = np.array([sqrt(_dfact(2 * l - 1) / (_dfact(2 * a - 1) * _dfact(2 * b - 1) * _dfact(2 * c - 1)))
                                 for a, b, c in powers])
            elif cart_norm == "axis":
                comp = np.ones(len(powers))
            else:
                raise ValueError("cart_norm must be 'auto', 'unit' or 'axis'")
            bshells.append((coords[atom], l, exps, bare, comp))
            t = np.zeros((len(powers), len(powers)))
            for col, p in enumerate(_MOLDEN_CART[l]):
                t[powers.index(p), col] = 1.0
            blocks.append(t)

    basis = Basis.from_shells(bshells)
    nfile = sum(b.shape[1] for b in blocks)
    tmat = np.zeros((basis.nao, nfile))
    r = c = 0
    for b in blocks:
        tmat[r:r + b.shape[0], c:c + b.shape[1]] = b
        r += b.shape[0]
        c += b.shape[1]

    cfile = np.zeros((nfile, len(mos)))
    for k, mo in enumerate(mos):
        for i, v in mo["coef"].items():
            if i >= nfile:
                raise ValueError("MO coefficient index beyond the number of basis functions")
            cfile[i, k] = v
    return MoldenData(basis=basis, atom_numbers=znum, coords=coords, mo_coeff=tmat @ cfile,
                      mo_occ=np.array([m["occ"] for m in mos]), mo_energy=np.array([m["ene"] for m in mos]),
                      mo_spin=np.array([m["spin"] for m in mos]), mo_sym=[m["sym"] for m in mos],
                      cart_norm=cart_norm, file_to_cart=tmat)


def overlap_matrix(basis):
    """Overlap matrix of the contracted functions of ``basis`` (pure numpy, closed form)."""
    from numpy.polynomial.hermite import hermgauss
    t, w = hermgauss(12)                      # exact up to polynomial degree 23
    xyz, ga, lmn, cf = basis.xyz, basis.ga, basis.lmn, basis.coeff
    n = basis.nprim
    sp = np.zeros((n, n))
    for i in range(n):
        p = ga[i] + ga
        P = (ga[i] * xyz[i] + ga[:, None] * xyz) / p[:, None]
        K = np.exp(-ga[i] * ga / p * ((xyz[i] - xyz) ** 2).sum(1))
        val = K / p ** 1.5
        for d in range(3):
            x = P[:, d][:, None] + t[None, :] / np.sqrt(p)[:, None]
            val = val * ((x - xyz[i, d]) ** lmn[i, d] * (x - xyz[:, d][:, None]) ** lmn[:, d][:, None] * w).sum(1)
        sp[i] = val * cf[i] * cf
    a = np.zeros((n, basis.nao))
    a[np.arange(n), basis.contr - 1] = 1.0
    return a.T @ sp @ a


def _orthonormality_error(data):
    occ = data.mo_occ > 1e-8
    if not occ.any():
        occ[:] = True
    s = overlap_matrix(data.basis)
    err = 0.0
    for spin in ("Alpha", "Beta"):
        sel = occ & (data.mo_spin == spin)
        if sel.any():
            c = data.mo_coeff[:, sel]
            err = max(err, np.abs(c.T @ s @ c - np.eye(sel.sum())).max())
    return err


def check_orbitals(data):
    """Largest deviation of the MOs from orthonormality (a diagnostic for the reader)."""
    s = overlap_matrix(data.basis)
    err = 0.0
    for spin in ("Alpha", "Beta"):
        sel = data.mo_spin == spin
        if sel.any():
            c = data.mo_coeff[:, sel]
            err = max(err, np.abs(c.T @ s @ c - np.eye(sel.sum())).max())
    return err


def from_molden(path, cart_norm="auto"):
    """Single-determinant / natural-orbital wavefunction from the occupations in a Molden file."""
    return read_molden(path, cart_norm).wavefunction()


def scf_energy(path, cart_norm="auto", reorder="auto", npoints=4000, seed=0):
    """Energy of the single determinant (or natural-orbital 1RDM, Hartree-Fock-like) defined by the
    occupations in a Molden file, with exact integrals from PySCF.

    Used to validate a program's Molden conventions: for an SCF Molden file the result must equal
    the program's SCF energy.  It matches only if the basis, the contraction normalisation, the
    Cartesian/spherical conversion, the component order and every sign are read correctly.
    Closed-shell (or spin-summed restricted) files only.  Returns (energy, fit_residual).
    """
    from pyscf import scf
    from pyscf.tools import molden as pm
    data = read_molden(path, cart_norm, reorder)
    mol = pm.load(path)[0]                      # used only for the basis function space
    mol.cart = True
    mol.verbose = 0
    mol.build()
    rng = np.random.default_rng(seed)
    pts = (data.coords[rng.integers(0, len(data.coords), npoints)]
           + rng.normal(scale=1.0, size=(npoints, 3)) * rng.choice([0.2, 1.0, 3.0], size=(npoints, 1)))
    a = mol.eval_gto("GTOval_cart", pts)
    occ = data.mo_occ > 1e-10
    mo_vals = data.basis.eval(pts) @ data.mo_coeff[:, occ]
    c = np.linalg.lstsq(a, mo_vals, rcond=None)[0]
    resid = float(np.abs(a @ c - mo_vals).max())
    dm = (c * data.mo_occ[occ]) @ c.T
    return float(scf.RHF(mol).energy_tot(dm)), resid
