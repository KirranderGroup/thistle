# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Factored and four-index two-particle densities on the same wavefunction: LiF, state-averaged
CASSCF(6,4), four singlet states (Molpro, tests/data/molpro), 90 Cartesian basis functions.

Computes the rotationally averaged total and inelastic scattering and the j2 moment S2 with both
routes, with the default cutoffs and with screening switched off (screen_tol = 1e-16), and plots
the percent difference 100 (factored - array) / array.  Writes factored_vs_dense_lif.png/.npz.
"""
import os
import time

import numpy as np

import pyxscat as xs
from pyxscat.interfaces import molpro

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "tests", "data", "molpro")
OUT, MLD = os.path.join(D, "lif_sa4_casscf.out"), os.path.join(D, "lif_sa4_casscf.mld")
BOHR = 0.529177210903
QA = np.linspace(0.0, 8.0, 81)                      # 1/Angstrom
Q = QA * BOHR                                       # 1/bohr for the kernels
STATES = (1, 2, 3, 4)


def compute():
    res = {"q": QA}
    for st in STATES:
        wd = molpro.from_multi(OUT, MLD, state=st, factored=False)
        wf = molpro.from_multi(OUT, MLD, state=st, factored=True)
        el = xs.scattering(wd, Q, kind="elastic")
        res[f"{st}_elastic"] = el
        for tag, kw in (("default", {}), ("tight", {"screen_tol": 1e-16})):
            for route, w in (("array", wd), ("factored", wf)):
                t = time.perf_counter()
                res[f"{st}_{tag}_{route}_total"] = xs.scattering(w, Q, **kw)
                res[f"{st}_{tag}_{route}_S2"] = xs.scattering(w, Q[1:], l=2, **kw)
                res[f"{st}_{tag}_{route}_s"] = time.perf_counter() - t
        print(f"state {st}: {wf.dm2.nterms} factor terms "
              f"({wf.dm2.nbytes / 1e6:.2f} MB vs {wd.dm2.nbytes / 1e6:.1f} MB)", flush=True)
    res["nterms"] = np.array(wf.dm2.nterms)
    res["MB"] = np.array([wf.dm2.nbytes / 1e6, wd.dm2.nbytes / 1e6])
    return res


def pct(a, b):
    return 100.0 * (a - b) / b


def plot(r, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from thistle.plotting import COLOURS, INK2, GRID, SLOTS
    q = r["q"]
    sc = {1: SLOTS[0], 2: SLOTS[1], 3: SLOTS[2], 4: SLOTS[6]}
    fig, ax = plt.subplots(2, 2, figsize=(10.5, 7.2), constrained_layout=True)

    a = ax[0, 0]                                     # the curves, state 1
    tot_a, tot_f = r["1_default_array_total"], r["1_default_factored_total"]
    el = r["1_elastic"]
    for name, ya, yf in (("total", tot_a, tot_f), ("inelastic", tot_a - el, tot_f - el)):
        a.plot(q, ya, color=COLOURS[name], lw=2.2, label=f"{name}, four-index array")
        a.plot(q, yf, color="k", lw=1.0, ls=(0, (4, 3)), label=f"{name}, factored")
    a.plot(q, el, color=COLOURS["elastic"], lw=1.6, label="elastic (1RDM only)")
    a.set(xlabel="q (Å⁻¹)", ylabel="I(q) (e.u.)", title="(a) LiF state 1: scattering, both routes")
    a.legend(fontsize=8, frameon=False)

    a = ax[0, 1]                                     # S2 curves
    for st in STATES:
        a.plot(q[1:], r[f"{st}_default_array_S2"], color=sc[st], lw=2.0, label=f"state {st}")
        a.plot(q[1:], r[f"{st}_default_factored_S2"], color="k", lw=0.9, ls=(0, (4, 3)))
    a.axhline(0, color=INK2, lw=0.6)
    a.set(xlabel="q (Å⁻¹)", ylabel="S₂(q) (e.u.)", title="(b) j₂ moment, all states (dashed: factored)")
    a.legend(fontsize=8, frameon=False)

    floor = 1e-16
    for col, (key, lab, sl) in enumerate((("total", "total I(q)", slice(None)), ("S2", "S₂(q)", slice(None)))):
        a = ax[1, col]
        qq = q if key == "total" else q[1:]
        for st in STATES:
            for tag, ls, lw in (("default", "-", 1.8), ("tight", ":", 1.4)):
                d = np.abs(pct(r[f"{st}_{tag}_factored_{key}"], r[f"{st}_{tag}_array_{key}"]))
                if key == "S2":                      # S2 crosses zero: relative to its largest value
                    s = r[f"{st}_{tag}_array_{key}"]
                    d = 100.0 * np.abs(r[f"{st}_{tag}_factored_{key}"] - s) / np.abs(s).max()
                a.semilogy(qq, np.maximum(d, floor), color=sc[st], ls=ls, lw=lw,
                           label=f"state {st}, {'default cutoffs' if tag == 'default' else 'no screening'}")
        a.axhline(100 * 2.2e-16, color=INK2, lw=0.8, ls="--")
        a.text(q[-1], 100 * 2.2e-16 / 2.2, "machine ε", ha="right", va="top", fontsize=8, color=INK2,
               bbox=dict(fc="white", ec="none", pad=1))
        a.set_ylim(floor, 1e-5)
        a.set(xlabel="q (Å⁻¹)",
              ylabel="|difference| (%)" + (" of max |S₂|" if key == "S2" else ""),
              title=f"({'c' if col == 0 else 'd'}) factored vs array: {lab}")
        a.grid(True, which="major", color=GRID, lw=0.6)
        if col == 0:
            a.legend(fontsize=7, ncol=2, frameon=False, loc="upper left")
        a.text(0.99, 0.03, "states 2 and 3 (degenerate Π pair) coincide", transform=a.transAxes,
               ha="right", fontsize=7.5, color=INK2)
    nk = r["nterms"]
    fig.suptitle(f"LiF SA-CASSCF(6,4), 90 basis functions: factored 2RDM ({nk[0]}+{nk[1]} terms, "
                 f"{r['MB'][0]:.2f} MB) vs four-index array ({r['MB'][1]:.0f} MB)", fontsize=11)
    fig.savefig(path, dpi=160)


if __name__ == "__main__":
    import sys
    npz = os.path.join(HERE, "factored_vs_dense_lif.npz")
    if "--plot-only" in sys.argv:
        plot(dict(np.load(npz)), os.path.join(HERE, "factored_vs_dense_lif.png"))
        sys.exit()
    r = compute()
    np.savez(npz, **r)
    for st in STATES:
        for tag in ("default", "tight"):
            t = np.abs(pct(r[f"{st}_{tag}_factored_total"], r[f"{st}_{tag}_array_total"])).max()
            s = r[f"{st}_{tag}_array_S2"]
            s2 = 100 * np.abs(r[f"{st}_{tag}_factored_S2"] - s).max() / np.abs(s).max()
            print(f"state {st} {tag:8s} max |%diff| total {t:.2e}  S2 {s2:.2e}   "
                  f"time array {r[f'{st}_{tag}_array_s']:.1f} s factored {r[f'{st}_{tag}_factored_s']:.1f} s")
    plot(r, os.path.join(HERE, "factored_vs_dense_lif.png"))
