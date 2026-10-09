# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""User-facing objects: load anything, get a Molecule, ask it for scattering."""
import csv
import json
import os
import time
import warnings
from dataclasses import dataclass, field

import numpy as np

import pyxscat as xs
from pyxscat import aligned as _aligned
from pyxscat import iam as _iam
from pyxscat import electron as _electron
from pyxscat import grids as _grids

BOHR = 0.529177210903                       # Angstrom
Q_UNITS = {"1/A": 1.0 / BOHR, "1/bohr": 1.0}  # factor from bohr^-1 to the unit


# ----------------------------------------------------------------------------- loading
def _sniff(path):
    """Guess what a file is from its name and first lines."""
    low = path.lower()
    if low.endswith(".pun"):
        return "molpro_punch"
    if low.endswith((".h5", ".hdf5")):
        from pyxscat.interfaces import deck
        if deck.is_deck(path):
            return "deck"
        raise ValueError(f"{path} is an HDF5 file but not a thistle-dev deck")
    with open(path, errors="replace") as f:
        head = f.read(200_000)
    if "[Molden Format]" in head or "[MOLDEN FORMAT]" in head.upper()[:200]:
        return "molden"
    up = head.upper()
    if "&RASSCF" in up or "&MRCI" in up or ("MOLCAS" in up and "RASSCF" in up):
        return "molcas_output"
    if "PROGRAM * MULTI" in head or ("Molpro" in head and "MULTI" in head):
        return "molpro_output"
    if "PROGRAM SYSTEM MOLPRO" in head or "***  " in head[:10]:
        return "molpro_output"
    raise ValueError(f"cannot tell what {path} is; supported: Molden, Molcas/Molpro output, Molpro punch")


def _find_molden(path):
    """The Molden file written by the same run, if it sits next to the output."""
    base = os.path.splitext(path)[0]
    folder = os.path.dirname(path) or "."
    for cand in (base + ".rasscf.molden", base + ".molden", base + ".mld",
                 os.path.join(folder, "molden.mld"), os.path.join(folder, "molpro.mld")):
        if os.path.exists(cand):
            return cand
    others = [n for n in sorted(os.listdir(folder)) if n.lower().endswith((".molden", ".mld"))]
    if len(others) == 1:                      # only one candidate: it must be this one
        return os.path.join(folder, others[0])
    if others:
        raise FileNotFoundError(f"several Molden files next to {path} ({', '.join(others)}); "
                                f"say which one with molden='...'")
    return None


def load(source, molden=None, state=1, punch=None, name=None, hf=False, factored="auto", mo_coeff=None):
    """Load a wavefunction from (almost) anything.

    source
        * a PySCF object: SCF (RHF/RKS), CCSD, CISD, CASSCF/CASCI (any active-space solver,
          DMRG and SHCI included) or a standalone FCI / selected-CI solver (``fci.FCI(mf)``)
        * a Molden file (orbitals and occupations: elastic and mean-field scattering)
        * an OpenMolcas RASSCF output   (its ``.rasscf.molden`` is found automatically); if the
          run went on to the uncontracted MRCI with ``THRPrint = 0.0``, the MRCI state is read
        * a thistle-dev input deck (``deck.h5``, MO or dense-AO representation)
        * a Molpro MULTI output or punch file (the Molden file next to it is found).  The output
          is not needed: a punch file and the Molden file are enough, and a missing output is
          replaced by the punch file of the same name.  A Molden file together with
          ``punch=`` is read as that MULTI calculation.
    state
        root / state number for multi-state calculations (1 = lowest)
    hf
        the orbitals are those of a closed-shell Hartree-Fock calculation.  A Molden file holds
        only orbitals and occupations, which give the elastic scattering; for a closed-shell
        determinant they also fix the pair density, hence the total and inelastic scattering.
        This is not true for natural orbitals of a correlated calculation, so it has to be said.
        (A PySCF RHF object is recognised by itself.)
    factored
        how a correlated two-particle density is held: ``"auto"`` (default) as products of
        one-particle matrices beyond 60 basis functions (no n^4 array), ``True`` always,
        ``False`` never.  Same numbers either way (to ~1e-12); see THEORY §1.2 and §4.
    mo_coeff
        orbitals of a standalone PySCF FCI solver, which keeps only its integrals: the SCF
        object or the orbital array (``th.load(fci.FCI(mf).run(), mo_coeff=mf)``)
    """
    from pyxscat.interfaces import molcas, molden as mld, molpro
    if not isinstance(source, (str, os.PathLike)):
        cls = type(source).__name__
        auto_hf = "HF" in cls.upper() and "KS" not in cls.upper() and not hasattr(source, "ncas")
        return Molecule(_from_pyscf(source, state, factored, mo_coeff), name=name or cls,
                        origin={"source": "PySCF " + cls, "state": state}, hf=hf or auto_hf)
    path = os.fspath(source)
    if not os.path.exists(path):
        pun = punch or os.path.splitext(path)[0] + ".pun"
        if not os.path.exists(pun):
            raise FileNotFoundError(f"{path} not found (nor a punch file {os.path.basename(pun)} next to it)")
        warnings.warn(f"{path} not found: reading the CI vectors from {pun} and the Molden file instead",
                      stacklevel=2)
        path, punch = pun, None
    kind = _sniff(path)
    label = name or os.path.splitext(os.path.basename(path))[0]
    if kind == "deck":
        from pyxscat.interfaces import deck
        wf, info = deck.read(path)
        return Molecule(wf, name=label, origin={"source": path, "kind": "thistle-dev deck",
                                                  "rdm_mode": info["rdm_mode"], "states": info["states"]})
    if kind == "molden":
        if punch is not None:                      # Molden + punch: the MULTI calculation itself
            wf = molpro.from_multi(None, path, state=state, punch=punch, factored=factored)
            return Molecule(wf, name=label, origin={"source": punch, "molden": path, "kind": "molpro_punch",
                                                    "state": state})
        sib = os.path.splitext(path)[0] + ".pun"
        if os.path.exists(sib) and not hf:
            warnings.warn(f"a Molden file holds orbitals and occupations only (for a state-averaged MULTI: the "
                          f"state-averaged density), so only the elastic scattering follows from it. "
                          f"{os.path.basename(sib)} is next to it: run on that file (or pass punch=...) for the "
                          f"total, inelastic and state-specific scattering", stacklevel=2)
        return Molecule(mld.from_molden(path), name=label, origin={"source": path, "kind": "Molden"}, hf=hf)
    molden = molden or _find_molden(path)
    if molden is None:
        raise FileNotFoundError(f"no Molden file found next to {path}; pass molden='...'")
    if kind == "molcas_output":
        if molcas.has_mrci(path):                       # an MRCI after the RASSCF: the correlated state
            wf = molcas.from_mrci(path, molden, state=state, factored=factored)
            kind = "molcas_mrci"
        else:
            wf = molcas.from_rasscf(path, molden, root=state, factored=factored)
    elif kind == "molpro_punch":
        out = os.path.splitext(path)[0] + ".out"
        wf = molpro.from_multi(out if os.path.exists(out) else None, molden, state=state, punch=path,
                               factored=factored)
    else:
        pun = punch or (os.path.splitext(path)[0] + ".pun")
        wf = molpro.from_multi(path, molden, state=state, punch=pun if os.path.exists(pun) else None,
                               factored=factored)
    return Molecule(wf, name=label, origin={"source": path, "molden": molden, "kind": kind, "state": state})


def _from_pyscf(obj, state, factored="auto", mo_coeff=None):
    from pyxscat.interfaces import pyscf as xp
    cls = type(obj).__name__.upper()
    if hasattr(obj, "fcisolver") and hasattr(obj, "ncas"):
        st = None
        if isinstance(obj.ci, (list, tuple)):
            st = state - 1
        return xp.from_casscf(obj, state=st, factored=factored)
    if (type(obj).__module__ or "").startswith("pyscf.fci") and hasattr(obj, "make_rdm12"):
        st = state - 1 if isinstance(obj.ci, (list, tuple)) else None
        return xp.from_fci_solver(obj, mo_coeff, state=st, factored=factored)
    if "CCSD" in cls:
        return xp.from_cc(obj, factored=factored)
    if "CISD" in cls:
        return xp.from_ci(obj, factored=factored)
    if hasattr(obj, "make_rdm1") and hasattr(obj, "mo_coeff"):
        return xp.from_scf(obj)
    raise TypeError(f"don't know how to read a {type(obj).__name__}")


# ----------------------------------------------------------------------------- q grids
def _q_grid(q, qmax, points, unit):
    if unit not in Q_UNITS:
        raise ValueError(f"unit must be one of {list(Q_UNITS)}")
    if q is None:
        q = np.linspace(0.0, qmax, points)
    q = np.atleast_1d(np.asarray(q, dtype=float))
    return q, q / Q_UNITS[unit]          # (in the user's unit, in bohr^-1)


# ----------------------------------------------------------------------------- results
@dataclass
class Result:
    """Scattering curves on one q grid.  ``res["total"]``, ``res.q``, ``res.save()``, ``res.plot()``."""

    q: np.ndarray
    unit: str
    curves: dict
    name: str = ""
    meta: dict = field(default_factory=dict)

    def __getitem__(self, key):
        return self.curves[key]

    def __getattr__(self, key):
        curves = self.__dict__.get("curves", {})
        if key in curves:
            return curves[key]
        raise AttributeError(key)

    def __repr__(self):
        rows = [f"Scattering of {self.name or 'molecule'}: {len(self.q)} q points, "
                f"{self.q.min():g} to {self.q.max():g} {self.unit}"]
        ok = np.flatnonzero(np.all(np.isfinite(np.array(list(self.curves.values()))), axis=0))
        i0 = ok[0] if ok.size else 0
        for k, v in self.curves.items():
            rows.append(f"  {k:<17} I(q={self.q[i0]:.3g}) = {v[i0]:12.6g}   I(qmax) = {v[-1]:12.6g}")
        return "\n".join(rows)

    def save(self, path):
        """Write a CSV (one column per curve) and a JSON file with how it was computed."""
        with open(path, "w", newline="") as f:
            f.write(f"# {self.name}  q in {self.unit}  " + json.dumps(self.meta) + "\n")
            w = csv.writer(f)
            w.writerow([f"q ({self.unit})"] + list(self.curves))
            for i, qq in enumerate(self.q):
                w.writerow([f"{qq:.8g}"] + [f"{v[i]:.12g}" for v in self.curves.values()])
        return path

    def plot(self, path=None, curves=None, ax=None):
        """Plot the curves (all by default) in the THISTLE style (:mod:`thistle.plotting` has a
        plotter for every kind of result).  Saves to ``path`` if given, returns the axes."""
        from . import plotting
        names = " ".join(self.curves)
        if "sM" in names and curves is None:
            fig = plotting.sm(self, path=path, ax=ax)
        elif "%diff" in names and curves is None:
            fig = plotting.difference(self, path=path, ax=ax)
        else:
            fig = plotting.scattering(self, curves=curves, path=path, ax=ax)
        return ax if ax is not None else fig.axes[0]


def compare(a, b, curve="total"):
    """Percent difference 100 (a - b) / b of one curve of two results on the same grid: curves
    (:class:`Result`), detector images or detector patterns (e.g. ab initio vs IAM)."""
    if isinstance(a, (DetectorPattern, DetectorImage)):
        va, vb = (a.values, b.values) if isinstance(a, DetectorPattern) else (a.images, b.images)
        if curve not in va or curve not in vb:
            curve = "elastic"
        with np.errstate(divide="ignore", invalid="ignore"):
            pd = 100.0 * (va[curve] - vb[curve]) / vb[curve]
        meta = dict(a.meta, compared_with=b.name, ylabel=f"% difference ({curve})")
        name = f"{a.name} vs {b.name}"
        if isinstance(a, DetectorPattern):
            if a.qvec.shape != b.qvec.shape or not np.allclose(a.qvec, b.qvec):
                raise ValueError("the two patterns are on different q grids")
            return DetectorPattern(a.qvec, {f"%diff {curve}": pd}, a.unit, name=name, meta=meta)
        if not (np.allclose(a.qx, b.qx) and np.allclose(a.qy, b.qy)):
            raise ValueError("the two images are on different q grids")
        return DetectorImage(a.qx, a.qy, {f"%diff {curve}": pd}, a.unit, name=name, meta=meta)
    if not np.allclose(a.q, b.q):
        raise ValueError("the two results have different q grids")
    with np.errstate(divide="ignore", invalid="ignore"):
        pd = 100.0 * (a[curve] - b[curve]) / b[curve]
    return Result(a.q, a.unit, {f"%diff {curve}": pd}, name=f"{a.name} vs {b.name}",
                  meta={**a.meta, "ylabel": f"% difference ({curve})", "yscale": None, "plot_curves": None,
                        "compared_with": b.name})


def sm(molecule, iam, curve="total"):
    """Modified scattering intensity sM(s) = s (I - I_at) / I_at of electron diffraction, with
    I_at the atomic part of the independent atom model (s = q, in the unit of the results).

    Two curves: the ab initio molecule, and the IAM itself (s I_mol / I_at, the usual
    simulation); their difference is the effect of bonding and electron correlation."""
    if not np.allclose(molecule.q, iam.q):
        raise ValueError("the two results have different q grids")
    if "atomic" not in iam.curves:
        raise ValueError("sM needs an electron IAM: mol.iam(probe='electron')")
    with np.errstate(divide="ignore", invalid="ignore"):
        val = molecule.q * (molecule[curve] - iam["atomic"]) / iam["atomic"]
        ref = molecule.q * iam["molecular"] / iam["atomic"]
    return Result(molecule.q, molecule.unit, {f"sM ({curve})": val, "sM (IAM)": ref},
                  name=f"{molecule.name}: sM(s)",
                  meta={**molecule.meta, "ylabel": f"sM(s) ({molecule.unit})", "probe": "electron",
                        "kind": "modified scattering sM(s)", "plot_curves": None, "yscale": None})


def _colour(values, log=True):
    """Colour map and normalisation: a symmetric diverging scale for signed data (differences),
    a log scale for positive intensities."""
    from matplotlib.colors import LogNorm, Normalize
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if v.size and v.min() < 0:
        m = float(np.max(np.abs(v))) or 1.0
        return "RdBu_r", Normalize(-m, m)
    if log and v.size and v.min() > 0:
        return "viridis", LogNorm()
    return "viridis", None


@dataclass
class DetectorImage:
    """A detector image at fixed orientation: ``img["total"]`` is a 2D array (NaN off the detector)."""

    qx: np.ndarray
    qy: np.ndarray
    images: dict
    unit: str
    name: str = ""
    meta: dict = field(default_factory=dict)

    def __getitem__(self, key):
        return self.images[key]

    def plot(self, path=None, curve="total", log=True):
        """Show one image (log colour scale by default: the structure sits far below the centre)."""
        from . import plotting
        return plotting.detector(self, curve=curve if curve in self.images else None, path=path, log=log).axes[0]

    def save(self, path):
        meta = dict(self.meta, q_unit=self.unit, name=self.name)
        np.savez(path, qx=self.qx, qy=self.qy, **self.images, meta=json.dumps(meta))
        return path


@dataclass
class DetectorPattern:
    """Values at an arbitrary list of laboratory scattering vectors (polar detector grid, a q-grid
    file, ...): ``pat["total"]`` is an array over the points, ``pat.qvec`` the vectors (bohr^-1)."""

    qvec: np.ndarray
    values: dict
    unit: str
    name: str = ""
    meta: dict = field(default_factory=dict)

    def __getitem__(self, key):
        return self.values[key]

    @property
    def q(self):
        """|q| of every point, in the result's unit."""
        return np.linalg.norm(self.qvec, axis=1) * Q_UNITS[self.unit]

    def plot(self, path=None, curve="total", log=True):
        """Colour map over (q_x, q_y) of the detector (log scale by default)."""
        from . import plotting
        return plotting.detector(self, curve=curve if curve in self.values else None, path=path, log=log).axes[0]

    def save(self, path):
        """``.npz``: arrays; any other extension: text, columns qx qy qz |q| (bohr^-1, as in
        qgrid.dat) followed by one column per curve."""
        meta = dict(self.meta, q_unit=self.unit, name=self.name)
        if str(path).endswith(".npz"):
            np.savez(path, qvec=self.qvec, **self.values, meta=json.dumps(meta))
            return path
        qn = np.linalg.norm(self.qvec, axis=1)
        cols = np.column_stack([self.qvec, qn] + [np.asarray(v) for v in self.values.values()])
        head = (f"{self.name}  " + json.dumps(meta) + "\nqx qy qz |q| (bohr^-1)  "
                + "  ".join(self.values))
        fmt = ["%19.12f"] * 4 + ["%22.14e"] * len(self.values)
        np.savetxt(path, cols, fmt=fmt, header=head)
        return path


def _lab_k0(probe, wavelength, energy_kev):
    """Beam wavenumber (bohr^-1) from a wavelength (Angstrom) or an energy (keV; photons or
    electrons according to the probe).  Defaults: 1 Angstrom X-rays, 100 keV electrons."""
    if wavelength is not None:
        return 2 * np.pi / (wavelength / BOHR)
    if _probe(probe) == "electron":
        return _grids.electron_k0(energy_kev or 100.0)
    if energy_kev is not None:
        return _grids.photon_k0(energy_kev)
    return 2 * np.pi / (1.0 / BOHR)


def _theta(theta_max, k0, qmax_b, electron):
    """Largest scattering angle of the polar grid (degrees): as given; X-rays 60 degrees (as in
    grid.f90); electrons the angle that reaches q_perp = qmax (default 12 1/A), since fast
    electrons scatter into a few degrees."""
    if theta_max is not None:
        return float(theta_max)
    if not electron:
        return 60.0
    q = qmax_b if qmax_b is not None else 12.0 * BOHR
    return float(np.degrees(np.arcsin(min(q / k0, 0.999))))


def _lab_grid(grid, k0, theta_max, rings, azimuths):
    """Laboratory scattering vectors (bohr^-1) of a 'polar' grid, a q-grid file or an array."""
    if isinstance(grid, str) and grid.lower() == "polar":
        return _grids.polar_grid(k0, theta_max, rings, azimuths)
    if isinstance(grid, (str, os.PathLike)):
        return _grids.read_qgrid(grid)
    q = np.atleast_2d(np.asarray(grid, float))
    if q.shape[1] != 3:
        raise ValueError("grid must be 'square', 'polar', a q-grid file or an (n, 3) array")
    return q


# ----------------------------------------------------------------------------- the molecule
PROBES = ("xray", "electron")
ELECTRON_REST_KEV = 510.99895


def _model(model):
    m = str(model).lower().replace(" ", "").replace("-", "").replace("_", "")
    if m in ("abinitio", "molecule", "wavefunction"):
        return "ab initio"
    if m in ("iam", "independentatoms", "independentatommodel"):
        return "iam"
    raise ValueError(f"model must be 'ab initio' or 'iam', not {model!r}")


def _probe(probe):
    p = str(probe).lower().replace("-", "").replace(" ", "")
    if p in ("xray", "xrays", "x"):
        return "xray"
    if p in ("electron", "electrons", "e", "ued", "ged"):
        return "electron"
    raise ValueError(f"probe must be 'xray' or 'electron', not {probe!r}")


def electron_wavelength(energy_kev):
    """Relativistic de Broglie wavelength (Angstrom) of an electron of kinetic energy E (keV)."""
    e = float(energy_kev) * 1e3                                  # eV
    return 12.264259 / np.sqrt(e * (1.0 + 0.978475e-6 * e))


def _gamma2(energy_kev):
    return 1.0 if energy_kev is None else (1.0 + float(energy_kev) / ELECTRON_REST_KEV) ** 2


def _cross_section(bracket, qb, unit, energy_kev=None):
    """4 [...] / q^4 in bohr^2, or Angstrom^2 when q is in 1/A; times gamma^2 if E is given."""
    cs = _electron.cross_section(bracket, qb) * _gamma2(energy_kev)
    return cs * (BOHR ** 2 if unit == "1/A" else 1.0)


def _area(unit):
    return "\u00c5\u00b2" if unit == "1/A" else "bohr\u00b2"


class Molecule:
    """A loaded wavefunction.  Ask it for ``scattering()``, ``iam()``, ``anisotropy()`` or
    ``detector()``, for X-rays (default) or electrons (``probe="electron"``).
    ``mol.wf`` is the underlying :class:`pyxscat.Wavefunction`."""

    def __init__(self, wf, name="", origin=None, hf=False):
        self.wf = wf
        self.name = name
        self.origin = origin or {}
        self.hf = bool(hf) and wf.dm2 is None
        if self.hf and not _single_determinant(wf):
            raise ValueError(f"{name}: hf=True (--hf) needs the orbitals of a closed-shell Hartree-Fock run "
                             "(occupations 0 or 2); these are not, so only the elastic scattering can be "
                             "computed from them")
        if self.hf:
            self.origin["hf"] = True

    @property
    def correlated(self):
        return self.wf.dm2 is not None

    def write_deck(self, path, observable="total", qmax=None, points=200, unit="1/A", **kw):
        """Write this wavefunction as a thistle-dev input deck (HDF5 schema v1, AO
        representation) for its Fortran and Julia engines; ``observable`` "total" or
        "elastic".  The q grid runs from ~0 to ``qmax`` (default 4 1/A).  Other options:
        :func:`pyxscat.interfaces.deck.write`."""
        from pyxscat.interfaces import deck
        qmax = (4.0 if unit == "1/A" else 4.0 * BOHR) if qmax is None else qmax
        q = np.linspace(1e-10, qmax / Q_UNITS[unit], points)
        return deck.write(path, self.wf, observable=observable, q=q, **kw)

    @property
    def has_total(self):
        """Total and inelastic scattering are available (2RDM, or closed-shell HF)."""
        return self.correlated or self.hf

    def __repr__(self):
        b = self.wf.basis
        from pyxscat.core import FactoredDM2
        g = self.wf.dm2
        two = (f"2RDM factored, {sum(g.nterms)} terms, {g.nbytes / 1e6:.1f} MB" if isinstance(g, FactoredDM2)
               else f"2RDM array, {g.nbytes / 1e6:.0f} MB" if g is not None else "")
        kind = (f"correlated (1RDM + {two})" if self.correlated else
                "closed-shell Hartree-Fock" if self.hf else "orbitals only (1RDM): elastic only, --hf for HF")
        atoms = ""
        if self.wf.atom_numbers is not None:
            from collections import Counter
            from pyxscat.iam import ELEMENTS
            c = Counter(ELEMENTS[int(z)] for z in self.wf.atom_numbers)
            atoms = "".join(f"{k}{v if v > 1 else ''}" for k, v in c.items()) + ", "
        return (f"Molecule {self.name!r}: {atoms}{self.wf.nelec:g} electrons, {b.nao} basis functions, {kind}\n"
                f"  from {self.origin.get('source', '?')}")

    def _meta(self, **kw):
        return {"thistle": "0.2.0.dev0", "from": self.origin, "created": time.strftime("%Y-%m-%d %H:%M"), **kw}

    # ------------------------------------------------------------------ rotational average
    def scattering(self, q=None, qmax=4.0, points=201, unit="1/A", parts="standard", threads=None,
                   probe="xray", energy_kev=None):
        """Rotationally averaged scattering (gas phase, isotropic sample).

        probe
            ``"xray"``: intensity in electron units (Thomson cross section = 1).
            ``"electron"``: differential cross section dsigma/dOmega (first Born approximation)
            in Angstrom^2/sr (bohr^2/sr with ``unit="1/bohr"``); NaN at q = 0.  ``energy_kev``
            adds the relativistic factor gamma^2 of the beam electrons.
        parts
            ``"standard"``: total, elastic, inelastic (elastic only if the input has no 2RDM)
            ``"all"``: X-rays: also Hartree, exchange and correlation (cumulant) parts;
            electrons: also the nucleus-nucleus, electron-nucleus and electron-electron terms.
        """
        qu, qb = _q_grid(q, qmax, points, unit)
        if _probe(probe) == "electron":
            return self._electron(qu, qb, unit, parts, threads, energy_kev)
        wf = self.wf
        if self.correlated:
            d = xs.decompose(wf, qb, nthreads=threads)
            keys = ["total", "elastic", "inelastic"]
            if parts == "all":
                keys += ["hartree", "exchange", "cumulant"]
            curves = {k: d[k] for k in keys}
        elif self.hf:
            # closed-shell determinant: the pair density follows from the 1RDM exactly
            el = xs.elastic(wf.basis, qb, wf.dm1, nthreads=threads)
            ex = xs.exchange(wf.basis, qb, wf.dm1, nthreads=threads)
            curves = {"total": wf.nelec + el + ex, "elastic": el, "inelastic": wf.nelec + ex}
            if parts == "all":
                curves.update(hartree=el, exchange=ex)
        else:
            curves = {"elastic": xs.elastic(wf.basis, qb, wf.dm1, nthreads=threads)}
        return Result(qu, unit, curves, name=self.name,
                      meta=self._meta(q_unit=unit, probe="xray", kind="rotational average"))

    def _electron(self, qu, qb, unit, parts, threads, energy_kev, l=0):
        d = _electron.decompose(self.wf, qb, l=l, hf=self.hf, nthreads=threads)
        keys = ["total", "elastic", "inelastic"] if "total" in d else ["elastic"]
        brackets = {k: d[k] for k in keys}
        if parts == "all":
            brackets["nucleus-nucleus"] = d["nn"]
            brackets["electron-nucleus"] = -2.0 * d["ne"]
            brackets["electron-electron"] = d.get("ee_total", d["ee_elastic"])
        curves = {k: _cross_section(v, qb, unit, energy_kev) for k, v in brackets.items()}
        if l == 2:
            curves = {f"S2 {k}": v for k, v in curves.items()}
        return Result(qu, unit, curves, name=self.name,
                      meta=self._meta(q_unit=unit, probe="electron", energy_kev=energy_kev,
                                      kind="rotational average" if l == 0 else "j2 moment about z",
                                      ylabel=("S2: " if l else "") + f"d\u03c3/d\u03a9 ({_area(unit)}/sr)",
                                      yscale="symlog" if l else "log",
                                      # the nuclear and electronic parts cancel to 1 part in 10^4 at
                                      # small q: they are written to the file but not plotted
                                      plot_curves=[k for k in curves if "nucle" not in k and "electron-" not in k]))

    def anisotropy(self, q=None, qmax=4.0, points=201, unit="1/A", threads=None, probe="xray",
                   energy_kev=None):
        """Second Legendre moment S2(q) about the molecular z axis (partially aligned molecules):
        I(q, theta) = I0(q) + 5 <P2> P2(cos theta) S2(q) + ...  (the same relation for electrons,
        in cross-section units)."""
        qu, qb = _q_grid(q, qmax, points, unit)
        if _probe(probe) == "electron":
            return self._electron(qu, qb, unit, "standard", threads, energy_kev, l=2)
        wf = self.wf
        if self.correlated:
            curves = {"S2 total": xs.scattering(wf, qb, "total", nthreads=threads, l=2),
                      "S2 elastic": xs.elastic(wf.basis, qb, wf.dm1, nthreads=threads, l=2)}
        elif self.hf:
            curves = {"S2 total": xs.scattering(wf, qb, "hf_total", nthreads=threads, l=2),
                      "S2 elastic": xs.elastic(wf.basis, qb, wf.dm1, nthreads=threads, l=2)}
        else:
            curves = {"S2 elastic": xs.elastic(wf.basis, qb, wf.dm1, nthreads=threads, l=2)}
        return Result(qu, unit, curves, name=self.name,
                      meta=self._meta(q_unit=unit, probe="xray", kind="j2 moment about z", ylabel="S2(q)"))

    def iam(self, q=None, qmax=4.0, points=201, unit="1/A", basis="aug-cc-pvtz", source="ab initio",
            probe="xray", energy_kev=None):
        """Independent atom model on the same grid, for :func:`compare` (and :func:`sm` for electrons).

        For electrons the curves also contain ``atomic`` (sum over atoms) and ``molecular``
        (interference) parts."""
        qu, qb = _q_grid(q, qmax, points, unit)
        total = source != "tables"
        meta = self._meta(q_unit=unit, kind="independent atom model", basis=basis, source=source)
        if _probe(probe) == "electron":
            r = _iam.electron_for_wavefunction(self.wf, qb, basis=basis, source=source, total=total)
            keys = (["total", "elastic", "inelastic"] if total else ["elastic"]) + ["atomic", "molecular"]
            curves = {k: _cross_section(r[k], qb, unit, energy_kev) for k in keys}
            meta.update(probe="electron", energy_kev=energy_kev, yscale="log",
                        ylabel=f"d\u03c3/d\u03a9 ({_area(unit)}/sr)")
            return Result(qu, unit, curves, name="IAM", meta=meta)
        r = _iam.for_wavefunction(self.wf, qb, basis=basis, source=source, total=total)
        keys = ["total", "elastic", "inelastic"] if total else ["elastic"]
        meta.update(probe="xray")
        return Result(qu, unit, {k: r[k] for k in keys}, name="IAM", meta=meta)

    # ------------------------------------------------------------------ fixed orientation
    def detector(self, wavelength=None, qmax=None, pixels=101, orientation=(0.0, 0.0, 0.0), unit="1/A",
                 parts="standard", threads=None, probe="xray", energy_kev=None, grid="square",
                 theta_max=None, rings=49, azimuths=72, model="ab initio", iam_basis="aug-cc-pvtz",
                 iam_source="ab initio"):
        """Scattering for one molecular orientation (oriented molecules, single detector images).

        Beam: ``wavelength`` (Angstrom) or ``energy_kev`` (photon energy for X-rays, kinetic
        energy for electrons).  Defaults: 1 Angstrom X-rays, 100 keV electrons.  The beam runs
        along z; ``orientation`` = Euler angles (degrees, z-y-z) of the molecule.

        grid
            ``"square"``: a pixels x pixels image in (q_x, q_y) up to ``qmax`` (in ``unit``;
            default: as far as the detector reaches, at most 12 1/A) -> :class:`DetectorImage`.
            ``"polar"``: the qgrid3d/qcoord3 grid of the earlier code: ``rings`` rings up to the
            scattering angle ``theta_max`` (degrees; default 60 for X-rays, for electrons the angle
            reaching ``qmax``, default 12 1/A), ``azimuths`` points per ring, plus q = 0
            (49 rings, 72 azimuths, 60 degrees by default) -> :class:`DetectorPattern`.
            a file name (columns qx qy qz in bohr^-1, e.g. qgrid.dat) or an (n, 3) array of
            laboratory q vectors (bohr^-1) -> :class:`DetectorPattern`.
        model
            ``"ab initio"`` (the wavefunction) or ``"iam"`` (independent atoms on the same grid,
            free atoms in ``iam_basis``; compare the two with :func:`compare`).
        """
        electron = _probe(probe) == "electron"
        k0 = _lab_k0(probe, wavelength, energy_kev)
        model = _model(model)
        meta = self._meta(probe="electron" if electron else "xray", k0_bohr=k0, energy_kev=energy_kev,
                          wavelength_A=2 * np.pi / k0 * BOHR, orientation_deg=list(orientation), model=model)
        iam_kw = dict(basis=iam_basis, source=iam_source) if model == "iam" else None
        if iam_kw:
            meta.update(iam_basis=iam_basis, iam_source=iam_source)
        rot = _aligned.rotation_matrix(*np.radians(orientation))
        if isinstance(grid, str) and grid.lower() == "square":
            qmax_b = min(k0 * 0.999, 12.0 * BOHR) if qmax is None else qmax / Q_UNITS[unit]
            qv, on = _aligned.detector_q(k0, qmax_b, pixels)
            vals = self._fixed(_aligned.rotate_q(qv[on], rot), electron, unit, parts, threads, energy_kev, iam_kw)
            images = {}
            for k, v in vals.items():
                img = np.full(on.size, np.nan)
                img[on] = v
                images[k] = img.reshape(pixels, pixels)
            axis = np.linspace(-qmax_b, qmax_b, pixels) * Q_UNITS[unit]
            return DetectorImage(axis, axis, images, unit, name=self._label(model), meta=meta)
        theta_max = _theta(theta_max, k0, None if qmax is None else qmax / Q_UNITS[unit], electron)
        qlab = _lab_grid(grid, k0, theta_max, rings, azimuths)
        meta.update(grid=grid if isinstance(grid, str) else "array", theta_max_deg=theta_max, rings=rings,
                    azimuths=azimuths)
        vals = self._fixed(_aligned.rotate_q(qlab, rot), electron, unit, parts, threads, energy_kev, iam_kw)
        return DetectorPattern(qlab, vals, unit, name=self._label(model), meta=meta)

    def _label(self, model):
        return "IAM" if model == "iam" else self.name

    def _fixed(self, qm, electron, unit, parts, threads, energy_kev, iam_kw=None):
        """Fixed-orientation curves at molecular-frame vectors qm (bohr^-1)."""
        wf = self.wf
        if iam_kw:
            total = iam_kw["source"] != "tables"
            keys = ["total", "elastic", "inelastic"] if total else ["elastic"]
            if electron:
                r = _iam.electron_for_wavefunction(wf, qm, aligned=True, total=total, **iam_kw)
                qn = np.linalg.norm(qm, axis=1)
                return {k: _cross_section(r[k], qn, unit, energy_kev) for k in keys}
            r = _iam.for_wavefunction(wf, qm, aligned=True, total=total, **iam_kw)
            return {k: r[k] for k in keys}
        if electron:
            d = _electron.aligned(wf, qm, hf=self.hf, nthreads=threads)
            qn = np.linalg.norm(qm, axis=1)
            return {k: _cross_section(d[k], qn, unit, energy_kev) for k in ("total", "elastic", "inelastic") if k in d}
        if self.correlated:
            d = _aligned.decompose(wf, qm, nthreads=threads)
            keys = ["total", "elastic", "inelastic"] + (["cumulant"] if parts == "all" else [])
            return {k: d[k] for k in keys}
        if self.hf:
            el = _aligned.elastic(wf.basis, qm, wf.dm1, nthreads=threads)
            tot = _aligned.scattering(wf, qm, "hf_total", nthreads=threads)
            return {"total": tot, "elastic": el, "inelastic": tot - el}
        return {"elastic": _aligned.elastic(wf.basis, qm, wf.dm1, nthreads=threads)}

    # ------------------------------------------------------------------ partially aligned ensembles
    def aligned_ensemble(self, axis=(0.0, 0.0, 1.0), p2=1.0, grid="polar", wavelength=None, energy_kev=None,
                         theta_max=None, rings=49, azimuths=72, qmax=None, pixels=101, unit="1/A",
                         threads=None, probe="xray", model="ab initio", iam_basis="aug-cc-pvtz",
                         iam_source="ab initio"):
        """Detector pattern of molecules whose z axis is aligned along ``axis`` (laboratory frame,
        e.g. the laser polarisation; the beam runs along z) with degree of alignment
        ``p2`` = <P2(cos beta)>, from the isotropic signal and the j2 moment:

            I(q) = I0(|q|) + 5 p2 P2(q^ . axis^) S2(|q|)     (terms l >= 4 neglected)

        I0 and S2 are rotational averages, computed once per distinct |q| (once per ring on the
        polar grid).  Grids and ``model`` (``"ab initio"`` or ``"iam"``) as in :meth:`detector`.
        """
        electron = _probe(probe) == "electron"
        model = _model(model)
        k0 = _lab_k0(probe, wavelength, energy_kev)
        square = isinstance(grid, str) and grid.lower() == "square"
        if square:
            qmax_b = min(k0 * 0.999, 12.0 * BOHR) if qmax is None else qmax / Q_UNITS[unit]
            qv, on = _aligned.detector_q(k0, qmax_b, pixels)
            qlab = qv[on]
        else:
            theta_max = _theta(theta_max, k0, None if qmax is None else qmax / Q_UNITS[unit], electron)
            qlab = _lab_grid(grid, k0, theta_max, rings, azimuths)
        qn, uq, inv = _grids.unique_norms(qlab)
        pk = dict(probe=probe, energy_kev=energy_kev) if electron else {}
        qu = uq * Q_UNITS[unit]                    # same points; cross sections then come in unit^-2
        if model == "iam":
            i0c = self.iam(q=qu, unit=unit, basis=iam_basis, source=iam_source, **pk).curves
            s2c = self._iam_s2(qu, unit, electron, energy_kev, iam_basis, iam_source)
        else:
            i0c = self.scattering(q=qu, unit=unit, threads=threads, **pk).curves
            s2c = self.anisotropy(q=qu, unit=unit, threads=threads, **pk).curves
        vals = {}
        for k, v in i0c.items():
            if f"S2 {k}" in s2c:
                vals[k] = _grids.aligned_ensemble(qlab, axis, v[inv], s2c[f"S2 {k}"][inv], p2)
        meta = self._meta(probe="electron" if electron else "xray", k0_bohr=k0, energy_kev=energy_kev,
                          alignment_axis=list(map(float, axis)), p2=p2, kind="aligned ensemble (j2)", model=model,
                          grid=grid if isinstance(grid, str) else "array", theta_max_deg=theta_max,
                          rings=rings, azimuths=azimuths)
        if square:
            images = {}
            for k, v in vals.items():
                img = np.full(on.size, np.nan)
                img[on] = v
                images[k] = img.reshape(pixels, pixels)
            ax_ = np.linspace(-qmax_b, qmax_b, pixels) * Q_UNITS[unit]
            return DetectorImage(ax_, ax_, images, unit, name=self._label(model), meta=meta)
        return DetectorPattern(qlab, vals, unit, name=self._label(model), meta=meta)

    def _iam_s2(self, qu, unit, electron, energy_kev, basis, source):
        """j2 moment of the IAM about the molecular z axis, in the units of :meth:`anisotropy`."""
        qb = qu / Q_UNITS[unit]
        total = source != "tables"
        if electron:
            r = _iam.electron_for_wavefunction(self.wf, qb, basis=basis, source=source, total=total, l=2)
            return {f"S2 {k}": _cross_section(r[k], qb, unit, energy_kev)
                    for k in ("total", "elastic") if k in r}
        r = _iam.for_wavefunction(self.wf, qb, basis=basis, source=source, total=total, l=2)
        return {f"S2 {k}": r[k] for k in ("total", "elastic") if k in r}


def _single_determinant(wf):
    """True if the 1RDM is that of a closed-shell determinant (natural occupations 0 or 2)."""
    from pyxscat.interfaces.molden import overlap_matrix
    if wf.nelec is None:
        return False
    s = overlap_matrix(wf.basis)
    w, v = np.linalg.eigh(s)
    sh = (v * np.sqrt(np.clip(w, 0, None))) @ v.T
    occ = np.linalg.eigvalsh(sh @ np.asarray(wf.dm1) @ sh)
    return bool(np.all(np.minimum(np.abs(occ), np.abs(occ - 2.0)) < 1e-6))
