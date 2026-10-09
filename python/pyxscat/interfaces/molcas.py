# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""(Open)Molcas RASSCF/CASSCF wavefunctions: CI vector from the output file, orbitals from the
Molden file written by the same run (``*.rasscf.molden``).

The RASSCF input must print the full CI vector, ``PRWF = -1`` (or a threshold small enough
for the precision you need).  Molcas prints the coefficients with five decimals; the vector is
renormalised after reading.  Only calculations without point-group symmetry (C1) are supported
so far.

Example::

    from pyxscat.interfaces import molcas
    wf = molcas.from_rasscf("molcas.log", "molcas.rasscf.molden", root=1)
"""
import re

import numpy as np

from ..core import use_factored

from .. import ci
from .molden import read_molden


def read_rasscf_output(path):
    """Parse a RASSCF output.  Returns a dict with the orbital spaces, spin, and for every root
    printed: the CSF strings, coefficients, energy and natural occupations."""
    with open(path) as f:
        text = f.read()
    lines = text.splitlines()

    def last_int(pattern):
        vals = [int(m.group(1)) for m in re.finditer(pattern, text)]
        if not vals:
            raise ValueError(f"{path}: could not find {pattern!r}")
        return vals[-1]

    sym_lines = [l for l in lines if re.match(r"^\s*Symmetry species\s+\d+(\s+\d+)+\s*$", l)]
    if sym_lines:
        raise NotImplementedError("RASSCF with point-group symmetry is not supported yet; rerun in C1")

    info = {
        "nclosed_electrons": last_int(r"Number of closed shell electrons\s+(\d+)"),
        "nactive_electrons": last_int(r"Number of electrons in active shells\s+(\d+)"),
        "ninactive": last_int(r"Number of inactive orbitals\s+(\d+)"),
        "nfrozen": _optional_int(text, r"Number of frozen shell orbitals\s+(\d+)"),
    }
    m = re.findall(r"Spin quantum number\s+([\d.]+)", text)
    info["spin"] = float(m[-1]) if m else 0.0
    nras = [_optional_int(text, rf"Number of {k} orbitals\s+(\d+)") for k in ("RAS1", "RAS2", "RAS3")]
    info["nactive"] = sum(nras)
    if info["nactive"] == 0:
        info["nactive"] = _optional_int(text, r"Number of active orbitals\s+(\d+)")

    roots = {}
    i = 0
    while i < len(lines):
        mm = re.search(r"printout of CI-coefficients larger than\s+\S+\s+for root\s+(\d+)", lines[i])
        if not mm:
            i += 1
            continue
        root = int(mm.group(1))
        energy = None
        csfs, coeffs = [], []
        i += 1
        while i < len(lines) and not lines[i].strip().startswith("conf/sym"):
            e = re.search(r"energy=\s*(-?[\d.]+)", lines[i])
            if e:
                energy = float(e.group(1))
            i += 1
        i += 1
        while i < len(lines) and lines[i].strip():
            w = lines[i].split()
            if not w[0].isdigit():
                break
            csfs.append("".join(w[1:-2]))
            coeffs.append(float(w[-2]))
            i += 1
        roots[root] = {"csfs": csfs, "coeffs": np.array(coeffs), "energy": energy}   # last printout wins
    if not roots:
        raise ValueError(f"{path}: no CI printout found (use PRWF = -1 in the RASSCF input)")
    for mm in re.finditer(r"Natural orbitals and occupation numbers for root\s+(\d+)\s*\n\s*sym 1:([^\n]*(?:\n\s{10,}[\d.\s]+)*)", text):
        r = int(mm.group(1))
        if r in roots:
            roots[r]["occupations"] = np.array([float(v) for v in mm.group(2).split()])
    for mm in re.finditer(r"RASSCF root number\s+(\d+)\s+Total energy:\s+(-?[\d.]+)", text):
        r = int(mm.group(1))
        if r in roots:
            roots[r]["energy"] = float(mm.group(2))
    info["roots"] = roots
    return info


def _optional_int(text, pattern):
    vals = re.findall(pattern, text)
    return int(vals[-1]) if vals else 0


def rasscf_rdms(info, root=1, factored=False):
    """Spin-summed 1- and 2-RDMs (MO basis: closed shells, then active orbitals)."""
    r = info["roots"][root]
    nact = info["nactive"]
    nel = info["nactive_electrons"]
    ms2 = int(round(2 * info["spin"]))
    for csf in r["csfs"]:
        if len(csf) != nact:
            raise ValueError(f"CSF {csf!r} does not have {nact} active orbitals")
    dets, coeffs = ci.csfs_to_determinants(r["csfs"], r["coeffs"], ms2)
    na = (nel + ms2) // 2
    nb = nel - na
    dm1a, dm2a = ci.make_rdm12(nact, na, nb, dets, coeffs, normalize=True)
    ncore = info["nfrozen"] + info["ninactive"]
    if 2 * ncore != info["nclosed_electrons"]:
        raise ValueError("closed-shell electrons and inactive orbitals do not match")
    return ci.embed_core(dm1a, dm2a, ncore, factored=factored)


def from_rasscf(output, molden, root=1, cart_norm="auto", factored="auto"):
    """Wavefunction (basis + 1RDM + 2RDM) of a RASSCF/CASSCF root."""
    info = read_rasscf_output(output)
    if root not in info["roots"]:
        raise ValueError(f"root {root} not printed in {output}; available: {sorted(info['roots'])}")
    data = read_molden(molden, cart_norm)
    dm1, dm2 = rasscf_rdms(info, root, factored=use_factored(factored, data.basis.nao))
    wf = data.wavefunction(dm1, dm2, label=f"RASSCF root {root}")
    return wf


# ----------------------------------------------------------------------------- MRCI (uncontracted)
_MRCI_CASE = {"0": "0", "1": "u", "2": "d", "3": "2"}     # GUGA case numbers -> step letters


def has_mrci(path):
    """True if the output holds an MRCI section with printed CI coefficients."""
    with open(path, errors="replace") as f:
        t = f.read()
    return "&MRCI" in t and "FINAL RESULTS FOR STATE NR" in t and "CI-COEFFICIENTS LARGER THAN" in t


def read_mrci_output(path):
    """Parse the MRCI (uncontracted multireference SDCI / ACPF) section of an OpenMolcas output
    run with ``THRPrint = 0.0`` (the complete CI vector printed).

    Returns a dict: ``nfrozen`` (orbitals frozen at the integral transformation, doubly occupied
    and not shown), ``norb_ci`` (correlated orbitals = inactive + active + virtual), ``nelec``
    (correlated electrons), ``spin`` (2S+1), and ``states`` {k: dict(energy, csfs, coeffs)} where
    each CSF is ``(orbitals, step)``: the 1-based orbital indices (within the symmetry; with
    ``Group = NoSym`` these are global) in the order of the spin coupling (virtual, inactive,
    active), and the step string over them (``0``, ``u``, ``d``, ``2``).  Orbitals not listed are
    empty.  Only ``NoSym`` runs are supported (one irrep).
    """
    with open(path) as f:
        text = f.read()
    i0 = text.find("&MRCI")
    if i0 < 0:
        raise ValueError(f"{path}: no MRCI section")
    t = text[i0:]
    m = re.search(r"POINT GROUP ORDER\s+(\d+)", t)
    if m and int(m.group(1)) != 1:
        raise NotImplementedError("MRCI reader: run with Group = NoSym (one irrep)")
    nel = int(re.search(r"CORRELATED ELECTRONS\s+(\d+)", t).group(1))
    spin2 = float(re.search(r"SPIN QUANTUM NUMBER\s+([\d.]+)", t).group(1))
    norb_ci = int(re.search(r"SUM:ORBITALS IN CI\s+(\d+)", t).group(1))
    nfrozen = int(re.search(r"PRE-FROZEN ORBITALS\s+(\d+)", t).group(1))
    ninact = int(re.search(r"\n\s*INACTIVE ORBITALS([ \d]+)\n", t).group(1).split()[-1])   # last column: the sum
    nact = int(re.search(r"\n\s*ACTIVE ORBITALS([ \d]+)\n", t).group(1).split()[-1])
    states = {}
    blocks = re.split(r"FINAL RESULTS FOR STATE NR\s+(\d+)", t)
    for k in range(1, len(blocks), 2):
        st = int(blocks[k]); body = blocks[k + 1]
        energy = float(re.search(r"\n\s+CI ENERGY:\s+(-?[\d.]+)", body).group(1))
        end = body.find("NATURAL ORBITALS")
        body = body[:end] if end > 0 else body
        csfs, coeffs = [], []
        for m in re.finditer(r"CONFIGURATION\s+\d+\s+COEFFICIENT\s+(-?[\d.]+)\s+\w+\s*\n"
                             r"\s*SYMMETRY([ \d]+)\n\s*ORBITALS([ \d]+)\n\s*OCCUPATION([ \d]+)\n"
                             r"\s*SPIN-COUPLING([ \d]+)", body):
            c = float(m.group(1))
            if c == 0.0:
                continue
            orbs = [int(x) for x in m.group(3).split()]
            step = "".join(_MRCI_CASE[x] for x in m.group(5).split())
            csfs.append((orbs, step)); coeffs.append(c)
        states[st] = {"energy": energy, "csfs": csfs, "coeffs": np.array(coeffs)}
    if not states:
        raise ValueError(f"{path}: no MRCI states with CI coefficients (use THRPrint = 0.0)")
    return {"nfrozen": nfrozen, "norb_ci": norb_ci, "ninternal": ninact + nact, "nelec": nel,
            "spin": int(round(2 * spin2 + 1)), "states": states}


def _place(ab, bb, pos):
    """Determinant bits over the positions of a step string -> bits over orbital indices ``pos``,
    with the sign of reordering the creation operators into increasing orbital order (the
    convention of the determinant builders; the coupling order of the string is not monotonic:
    virtual, inactive, active)."""
    a = b = 0
    sign = 1
    for bits in (ab, bb):
        occ = [pos[j] for j in range(len(pos)) if bits >> j & 1]
        inv = sum(1 for i in range(len(occ)) for k in range(i + 1, len(occ)) if occ[i] > occ[k])
        if inv % 2:
            sign = -sign
    for j, p in enumerate(pos):
        if ab >> j & 1:
            a |= 1 << p
        if bb >> j & 1:
            b |= 1 << p
    return a, b, sign


def mrci_rdms(info, state=1, factored=False, method="external"):
    """Spin-summed 1- and 2-RDMs of one MRCI state in the MO basis (frozen orbitals, then the
    correlated ones, as numbered by Molcas).  The CSFs are expanded in determinants, the
    determinant bits placed on the orbital indices, and the RDMs built over the correlated
    orbitals (inactive + active + virtual, any number): ``method="external"`` (default) uses
    the at-most-two-electrons-in-the-virtuals structure (:func:`pyxscat.ci.make_rdm12_external`),
    any other value the general builders of :func:`pyxscat.ci.make_rdm12`.  The frozen
    orbitals are added as a doubly occupied core.  ``factored=True`` returns a FactoredDM2
    (feasible up to ~80 orbitals in all: beyond that the eigendecomposition of the active
    block is too large, keep the four-index array)."""
    s = info["states"][state]
    n = info["norb_ci"]
    nf = info["nfrozen"]
    ms2 = (info["spin"] - 1) % 2
    acc = {}
    for (orbs, step), cv in zip(s["csfs"], s["coeffs"]):
        pos = [o - nf - 1 for o in orbs]                   # 0-based index among the correlated orbitals
        for v, ab, bb in ci.csf_to_determinants(step, ms2):
            a, b, sign = _place(ab, bb, pos)
            acc[a, b] = acc.get((a, b), 0.0) + sign * v * cv
    dets = list(acc)
    coeffs = np.array([acc[d] for d in dets])
    na = bin(dets[0][0]).count("1"); nb = bin(dets[0][1]).count("1")
    if na + nb != info["nelec"]:
        raise ValueError("electron count of the determinants does not match the output")
    if method == "external":
        d1, d2 = ci.make_rdm12_external(n, info["ninternal"], dets, coeffs, normalize=True)
    else:
        d1, d2 = ci.make_rdm12(n, na, nb, dets, coeffs, normalize=True, method=method)
    return ci.embed_core(d1, d2, nf, factored=factored)


def from_mrci(output, molden, state=1, cart_norm="auto", factored="auto"):
    """Wavefunction (basis + 1RDM + 2RDM) of an OpenMolcas MRCI state.  ``molden`` is the
    RASSCF Molden file of the same run (the orbitals MOTRA transformed)."""
    info = read_mrci_output(output)
    data = read_molden(molden, cart_norm)
    f = use_factored(factored, data.basis.nao) and info["nfrozen"] + info["norb_ci"] <= 80
    dm1, dm2 = mrci_rdms(info, state, factored=f)
    return data.wavefunction(dm1, dm2, label=f"MRCI state {state}")
