# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Turn results.json of run_independent_checks.py into the line-by-line report results.md."""
import json, os
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
rows = json.load(open(os.path.join(HERE, "results.json")))
n_ok = sum(r["ok"] for r in rows)


def e(x):
    return "—" if x != x else f"{x:.1e}"


def c(t):
    return str(t).replace("|", "\\|")


for r in rows:
    if r["area"].startswith("J. Readers"):
        r["area"] = "J. Readers (numbers printed by Molpro and OpenMolcas)"


out = [f"# THISTLE: independent checks, line by line",
       "",
       f"Run of `validation/independent/run_independent_checks.py`, 9 October 2026: **{n_ok} of {len(rows)} checks pass.**",
       "",
       "Every line compares a THISTLE result with a reference that does not go through THISTLE's kernels. "
       "The references are:",
       "",
       "- **PySCF analytic Fourier transforms** of basis-function pairs (`ft_aopair`) with PySCF's own density matrices "
       "(`make_rdm1`, `make_rdm2`, `mcscf.addons.make_rdm12`), averaged over directions on Lebedev grids "
       "(2030 points, degree 77; 5810 points above q = 8 bohr⁻¹). Shares no code with THISTLE.",
       "- **Numerical integration** of the electron density on PySCF's DFT grids (no Gaussian integrals at all).",
       "- **mpmath** at 40 digits; **xraydb** atomic form factors (Waasmaier–Kirfel); **CODATA** constants (scipy).",
       "- **Closed forms and exact limits**: N², N(N−1), N ΣZ, the dipole limit q²μ²/3, the hydrogen atom.",
       "- **Numbers printed by other programs** (Molpro energies, OpenMolcas energies and occupations), the stored output of the "
       "earlier production kernels, and thistle-dev's Fortran engine.",
       "",
       "Columns: *abs* = max |THISTLE − reference| over the points listed; *rel* = abs / max |reference| "
       "(for lines marked *abs*, the tolerance applies to the absolute value). q in bohr⁻¹; "
       "Lebedev checks at q = 0.5, 1.5, 3.0 and 5.0 bohr⁻¹ unless stated.",
       ""]
areas = OrderedDict()
for r in rows:
    areas.setdefault(r["area"], []).append(r)
out += ["## Summary", "", "| Area | Checks | Pass | Largest deviation (rel or abs) |", "|---|---|---|---|"]
for a, rs in areas.items():
    out.append(f"| {a} | {len(rs)} | {sum(x['ok'] for x in rs)} | {e(max(x['drel'] for x in rs))} |")
out.append("")
for a, rs in areas.items():
    out += [f"## {a}", "", "| ID | Quantity | System | Independent reference | abs | rel / abs | Tolerance | Result |",
            "|---|---|---|---|---|---|---|---|"]
    for r in rs:
        kind = "abs" if r["absolute"] else "rel"
        out.append(f"| {r['id']} | {c(r['quantity'])} | {c(r['system'])} | {c(r['reference'])} | {e(r['dabs'])} | "
                   f"{e(r['drel'])} ({kind}) | {e(r['tol'])} | {'PASS' if r['ok'] else 'FAIL'} |")
    notes = [r for r in rs if r["note"]]
    if notes:
        out.append("")
        for r in notes:
            out.append(f"- **{r['id']}**: {r['note']}")
    if a.startswith("K."):
        tp = os.path.join(HERE, "..", "..", "paper", "tab_three_programs_water_cas.txt")
        if os.path.exists(tp):
            out.append("- **Three programs** (stored run of `paper/fig_three_programs.py`, not repeated here): OpenMolcas "
                       "RASSCF read by THISTLE against PySCF CASCI in the same orbitals, water CASSCF(6,5)/cc-pVDZ: "
                       "total 1.6e-08, elastic 3.7e-08, inelastic 2.9e-07, S2 4.0e-06 (limited by the 8 decimals of the Molden file).")
    out.append("")
out += ["## What the checks found besides the numbers",
        "",
        "- **Open-shell SCF objects** are not accepted by `th.load` (`from_scf` handles closed shells only). The hydrogen-atom "
        "checks (F1–F3) therefore pass the spin-summed UHF density through `from_density`.",
        "- **Two values of the Bohr radius** are in use: CODATA 2018 (0.529177210903 Å) in `thistle/core.py` and `grids.py`, "
        "CODATA 2010 (0.52917721092 Å) in `iam.py` and the Molden reader. They differ by 1.7×10⁻¹¹ relative; harmless, "
        "but worth unifying.",
        "- **The OpenMolcas C₇H₈ test file was computed with RICD integrals.** Its printed energy cannot be rebuilt with exact "
        "integrals: an exact PySCF CASCI in the same orbitals is 0.024 Eh lower, and our densities reproduce that PySCF energy "
        "to 5×10⁻⁶. The reader is checked instead on the natural occupations Molcas prints (J4).",
        "- **Large-q behaviour.** In a Gaussian basis the two-electron term does not vanish at q = 60–80 bohr⁻¹ "
        "(3.5×10⁻⁵ for water): tight core pairs keep e^{−q²/4p} near 1. So I_tot → N is not a usable test there; "
        "B4 checks the actual values at q = 8–16 bohr⁻¹ against PySCF instead.",
        "- **DMRG.** The block2 density differs from exact CASCI by 2.6×10⁻⁷ (convergence of the density, although the "
        "energies agree to 10⁻¹³); THISTLE carries that into the scattering at the 10⁻¹⁰ level.",
        "",
        "## Rerun",
        "",
        "```",
        "python validation/independent/run_independent_checks.py   # about 5 minutes on 4 cores",
        "python validation/independent/make_report.py              # -> results.md",
        "```",
        ""]
open(os.path.join(HERE, "results.md"), "w").write("\n".join(out))
print(f"wrote results.md ({n_ok}/{len(rows)} pass)")
