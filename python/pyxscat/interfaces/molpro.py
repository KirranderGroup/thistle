# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Molpro MULTI (CASSCF/MCSCF) wavefunctions: CI vectors from the output, orbitals from a
Molden file written after the MULTI step (``put,molden,file.mld``).

Requirements on the Molpro input:
  * ``gprint,civector`` and ``config,det`` inside ``{multi ...}`` (determinant printout:
    strings of 2/a/b/0 over the active orbitals);
  * no point-group symmetry (``symmetry,nosym``);
  * the Molden file must hold the orbitals the CI vector refers to (the MULTI orbitals that
    ``put,molden`` writes right after the MULTI step).

Molpro writes each determinant with the creation operators in orbital order
(1alpha 1beta 2alpha ...); the phase to the alpha-first convention of :mod:`pyxscat.ci` is
applied on reading (verified: the energies of all states are reproduced to 1e-8 Hartree).

Molpro prints CI coefficients with eight decimals and only those above its print threshold.
Older versions default to 0.05, which can drop several per cent of the wavefunction: add
``gthresh,printci=0`` to the input.  A truncated vector triggers a warning (its printed weight
is stored as ``printed_weight``); each vector is renormalised after reading.

Example::

    from pyxscat.interfaces import molpro
    wf = molpro.from_multi("lif.out", "lif.mld", state=1)
"""
import re
import warnings

import numpy as np

from ..core import use_factored

from .. import ci
from .molden import read_molden

PRINT_WEIGHT_TOL = 1e-4


def read_multi_output(path):
    """Parse the last MULTI block of a Molpro output.

    Returns a dict with ``nclosed``, ``nactive``, ``nactive_electrons``, ``ms2`` (2 S from the
    spin label) and ``states``: {state number (1-based): {"dets": [(a, b)], "coeffs", "energy"}}.
    """
    with open(path) as f:
        text = f.read()
    starts = [m.start() for m in re.finditer(r"PROGRAM \* MULTI", text)]
    if not starts:
        raise ValueError(f"{path}: no MULTI program found")
    block = text[starts[-1]:]
    nxt = re.search(r"\n\s*PROGRAM \*", block[20:])
    if nxt:
        block = block[:nxt.start() + 20]

    def grab(pattern, cast=int, default=None):
        m = re.search(pattern, block)
        if not m:
            if default is not None:
                return default
            raise ValueError(f"{path}: could not find {pattern!r} in the MULTI block")
        return cast(m.group(1))

    if re.search(r"Point group\s+(?!C1\b)\S+", text) and not re.search(r"Point group\s+C1\b", text):
        raise NotImplementedError("Molpro calculations with point-group symmetry are not supported; use symmetry,nosym")
    info = {
        "nclosed": grab(r"Number of closed-shell orbitals:\s+(\d+)", default=0),
        "nfrozen": grab(r"Number of frozen orbitals:\s+(\d+)", default=0),
        "nactive": grab(r"Number of active\s+orbitals:\s+(\d+)"),
        "nactive_electrons": grab(r"Number of (?:active )?electrons:\s+(\d+)"),
    }
    spin = grab(r"Spin symmetry=(\w+)", str)
    names = {"Singlet": 0, "Doublet": 1, "Triplet": 2, "Quartet": 3, "Quintet": 4, "Sextet": 5}
    if spin not in names:
        raise ValueError(f"unknown spin label {spin!r}")
    info["ms2"] = names[spin]

    # recent Molpro: "CI Coefficients of symmetry 1 ... Energy:"; older: "CI vector ... TOTAL ENERGIES"
    blocks = list(re.finditer(r"CI Coefficients of symmetry\s+1\s*\n\s*=+\s*\n(.*?)\n\s*Energy:([^\n]*)", block, re.S))
    if not blocks:
        blocks = list(re.finditer(r"CI vector\s*\n\s*=+\s*\n(.*?)\n\s*TOTAL ENERGIES([^\n]*)", block, re.S))
    if not blocks:
        raise ValueError(f"{path}: no CI coefficients printed (use gprint,civector)")
    body, energies = blocks[-1].group(1), blocks[-1].group(2)
    energies = [float(e) for e in energies.split()]
    dets, rows = [], []
    for line in body.splitlines():
        w = line.split()
        if not w:
            continue
        string = w[0]
        if not set(string) <= set("2ab0"):
            raise ValueError(f"unexpected determinant string {string!r}: CSF printout (rerun MULTI with "
                             "config,det) or point-group symmetry")
        if len(string) != info["nactive"]:
            raise ValueError(f"determinant {string!r} does not have {info['nactive']} active orbitals")
        abits = sum(1 << i for i, c in enumerate(string) if c in "2a")
        bbits = sum(1 << i for i, c in enumerate(string) if c in "2b")
        dets.append((abits, bbits))
        # Molpro orders the creation operators orbital by orbital (alpha, beta); convert
        rows.append([ci.interleaved_phase(abits, bbits) * float(v) for v in w[1:]])
    rows = np.array(rows)
    info["states"] = {k + 1: {"dets": dets, "coeffs": rows[:, k], "energy": energies[k] if k < len(energies) else None,
                              "printed_weight": float(rows[:, k] @ rows[:, k])}
                      for k in range(rows.shape[1])}
    low = [k for k, s in info["states"].items() if s["printed_weight"] < 1 - PRINT_WEIGHT_TOL]
    if low:
        warnings.warn(f"{path}: the printed CI vector of state(s) {low} is truncated (weight "
                      f"{min(info['states'][k]['printed_weight'] for k in low):.4f}); Molpro prints only "
                      "coefficients above its threshold. Use gthresh,printci=0 (or 1d-8) for exact densities.",
                      stacklevel=2)
    return info


_DET_LINE = re.compile(r"^\s*([2ab0]+)((?:\s+-?\d+\.\d+)+)\s*$")


def read_punch(path, nclosed=None, molden=None, output=None):
    """CI vectors from a Molpro punch file (``punch,file.pun``).

    Older Molpro versions write the MULTI CI vector to the punch file with seven decimals and a
    much lower threshold than the output printout, so it is the better source when the output
    vector is truncated.  The punch file does not record the closed-shell count; it is taken from
    ``nclosed``, or from a Molpro ``output``, or from the electron count of the ``molden`` file
    (sum of its occupations).  Returns the same dict as :func:`read_multi_output`.
    """
    with open(path) as f:
        lines = f.read().splitlines()
    blocks, cur = [], []
    energies = {}
    for line in lines:
        m = _DET_LINE.match(line)
        if m:
            cur.append((m.group(1), [float(v) for v in m.group(2).split()]))
            continue
        if cur:
            blocks.append(cur)
            cur = []
        e = re.match(r"\s*MCSCF STATE\s+(\d+)\.\d+\s+Energy\s+(-?[\d.]+)", line)
        if e:
            energies[int(e.group(1))] = float(e.group(2))
    if cur:
        blocks.append(cur)
    if not blocks:
        raise ValueError(f"{path}: no CI vector in the punch file")
    rows = blocks[-1]
    norb = len(rows[0][0])
    if any(len(s) != norb for s, _ in rows):
        raise ValueError(f"{path}: determinant strings of different lengths (point-group symmetry?)")
    first = rows[0][0]
    na = sum(c in "2a" for c in first)
    nb = sum(c in "2b" for c in first)
    info = {"nactive": norb, "nactive_electrons": na + nb, "ms2": na - nb, "nfrozen": 0}
    if nclosed is None and output is not None:
        with warnings.catch_warnings():          # only nclosed is wanted from the output
            warnings.simplefilter("ignore")
            nclosed = read_multi_output(output)["nclosed"]
    if nclosed is None and molden is not None:
        nel = read_molden(molden).mo_occ.sum()
        nclosed = int(round((nel - (na + nb)) / 2))
    if nclosed is None:
        raise ValueError("give nclosed, a Molpro output or the Molden file to fix the number of closed shells")
    info["nclosed"] = nclosed
    dets, coeffs = [], []
    for string, vals in rows:
        abits = sum(1 << i for i, c in enumerate(string) if c in "2a")
        bbits = sum(1 << i for i, c in enumerate(string) if c in "2b")
        dets.append((abits, bbits))
        coeffs.append([ci.interleaved_phase(abits, bbits) * v for v in vals])
    coeffs = np.array(coeffs)
    info["states"] = {k + 1: {"dets": dets, "coeffs": coeffs[:, k], "energy": energies.get(k + 1),
                              "printed_weight": float(coeffs[:, k] @ coeffs[:, k])}
                      for k in range(coeffs.shape[1])}
    low = [k for k, st in info["states"].items() if st["printed_weight"] < 1 - PRINT_WEIGHT_TOL]
    if low:
        warnings.warn(f"{path}: the punched CI vector of state(s) {low} is truncated", stacklevel=2)
    return info


def multi_rdms(info, state=1, factored=False):
    """Spin-summed 1- and 2-RDMs of one MULTI state (MO basis: closed shells, then active).
    ``factored=True``: the 2RDM as a FactoredDM2 (see :func:`pyxscat.ci.embed_core`)."""
    s = info["states"][state]
    nel = info["nactive_electrons"]
    na = sum(bin(a).count("1") for a, _ in s["dets"][:1])
    nb = nel - na
    dm1a, dm2a = ci.make_rdm12(info["nactive"], na, nb, s["dets"], s["coeffs"], normalize=True)
    return ci.embed_core(dm1a, dm2a, info["nclosed"], factored=factored)


def transition_dm1(info, state_a, state_b):
    """Spin-summed transition 1-RDM <a|E_pq|b> between two MULTI states (MO basis)."""
    sa, sb = info["states"][state_a], info["states"][state_b]
    nel = info["nactive_electrons"]
    na = bin(sa["dets"][0][0]).count("1")
    t = ci.make_transition_rdm1(info["nactive"], na, nel - na, sa["dets"], sa["coeffs"], sb["dets"], sb["coeffs"])
    n = info["nclosed"] + info["nactive"]
    out = np.zeros((n, n))
    out[info["nclosed"]:, info["nclosed"]:] = t
    return out


def check_molden_consistency(info, data):
    """Largest difference between the state-averaged active 1RDM built from the printed CI
    vectors and the occupations in the Molden file (natural orbitals).  Large values mean the CI
    vector and the Molden orbitals do not refer to the same orbitals (or the printout is
    truncated)."""
    nc, na = info["nclosed"], info["nactive"]
    avg = sum(multi_rdms(info, k)[0] for k in info["states"]) / len(info["states"])
    act = avg[nc:nc + na, nc:nc + na]
    return float(np.abs(act - np.diag(data.mo_occ[nc:nc + na])).max())


def from_multi(output, molden, state=1, cart_norm="auto", check=True, punch=None, nclosed=None,
               factored="auto"):
    """Wavefunction (basis + 1RDM + 2RDM) of one MULTI state.

    ``punch``: take the CI vector from this punch file instead of the output printout (more
    determinants and digits in older Molpro versions).  ``output`` may then be None.
    ``factored``: keep the 2RDM as products of one-particle matrices (no n_AO^4 array); "auto"
    beyond ``pyxscat.core.FACTORED_ABOVE`` basis functions.
    """
    if punch is not None:
        info = read_punch(punch, nclosed=nclosed, molden=molden, output=output)
    else:
        info = read_multi_output(output)
    if state not in info["states"]:
        raise ValueError(f"state {state} not printed; available: {sorted(info['states'])}")
    data = read_molden(molden, cart_norm)
    dm1, dm2 = multi_rdms(info, state, factored=use_factored(factored, data.basis.nao))
    if check:
        dev = check_molden_consistency(info, data)
        if dev > 1e-3:
            warnings.warn(f"state-averaged 1RDM from the printed CI vectors differs from the Molden occupations "
                          f"by {dev:.3g}: truncated printout, unequal state weights, or orbitals that do not "
                          "match the CI vector", stacklevel=2)
    return data.wavefunction(dm1, dm2, label=f"MULTI state {state}")
