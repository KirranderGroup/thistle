# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Command-line driver (for batch and cluster jobs).

Examples::

    pyxscat --molden scf.molden --kind elastic --qmax 10 --npoints 201 -o elastic.dat
    pyxscat --molcas-log casscf.log --molden casscf.rasscf.molden --root 1 --kind decompose -o water.dat

Output: plain text, first column q (inverse bohr), then the requested quantities.
"""
import argparse
import sys
import time

import numpy as np

from . import core
from .interfaces import molden as molden_io
from .interfaces import molcas


def build_parser():
    p = argparse.ArgumentParser(prog="pyxscat", description="Ab initio X-ray and electron scattering from wavefunctions.")
    src = p.add_argument_group("wavefunction")
    src.add_argument("--molden", required=True, help="Molden file with basis set and orbitals")
    src.add_argument("--molcas-log", help="(Open)Molcas RASSCF output with the CI vector (PRWF = -1)")
    src.add_argument("--molpro-out", help="Molpro output with a MULTI CI vector (config,det; gprint,civector)")
    src.add_argument("--molpro-punch", help="Molpro punch file with the MULTI CI vector (more complete than the output)")
    src.add_argument("--root", type=int, default=1, help="RASSCF root or MULTI state (default 1)")
    src.add_argument("--cart-norm", default="auto", choices=["auto", "unit", "axis"],
                     help="normalisation of Cartesian d/f/g components in the Molden file")
    calc = p.add_argument_group("calculation")
    calc.add_argument("--kind", default="elastic",
                      choices=["elastic", "total", "inelastic", "two_electron", "hf_total", "decompose"])
    calc.add_argument("--probe", default="xray", choices=["xray", "electron"],
                      help="electron: Born-approximation electron scattering, written as the bracket "
                           "I_nn - 2 I_ne + I_ee (multiply by 4/q^4 for dsigma/dOmega in bohr^2)")
    calc.add_argument("--qvec", help="fixed orientation: text file with columns qx qy qz (bohr^-1, molecular frame)")
    calc.add_argument("--polar-grid", type=float, nargs=4, metavar=("E0_KEV", "THSMAX", "N1", "N2"),
                      help="fixed orientation on the polar detector grid of grid.f90 (qcoord3): photon energy "
                           "(keV), largest scattering angle (degrees), N1 rings, N2-1 azimuths, plus q = 0")
    calc.add_argument("--orientation", type=float, nargs=3, metavar=("ALPHA", "BETA", "GAMMA"),
                      help="with --qvec/--polar-grid: Euler angles (radians, zyz) of the molecule; q is then "
                           "in the lab frame")
    calc.add_argument("--legendre", type=int, default=0, choices=[0, 2],
                      help="0: isotropic signal; 2: second Legendre moment about the molecular z axis")
    calc.add_argument("--qmin", type=float, default=0.0)
    calc.add_argument("--qmax", type=float, default=10.0, help="inverse bohr")
    calc.add_argument("--npoints", type=int, default=201)
    calc.add_argument("--screen-tol", type=float, default=core.DEFAULTS["screen_tol"])
    calc.add_argument("--cutoffcentre", type=float, default=core.DEFAULTS["cutoffcentre"])
    calc.add_argument("--threads", type=int, default=0, help="OpenMP threads (0: OMP_NUM_THREADS)")
    p.add_argument("-o", "--output", default="scattering.dat")
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    t0 = time.time()
    if a.molcas_log and (a.molpro_out or a.molpro_punch):
        sys.exit("give either --molcas-log or --molpro-out, not both")
    if a.molcas_log:
        wf = molcas.from_rasscf(a.molcas_log, a.molden, root=a.root, cart_norm=a.cart_norm)
    elif a.molpro_out or a.molpro_punch:
        from .interfaces import molpro
        wf = molpro.from_multi(a.molpro_out, a.molden, state=a.root, cart_norm=a.cart_norm, punch=a.molpro_punch)
    else:
        if a.kind in ("total", "inelastic", "two_electron", "decompose"):
            sys.exit(f"--kind {a.kind} needs a correlated two-particle density (e.g. --molcas-log); "
                     "(--molcas-log or --molpro-out); with a Molden file alone use elastic or hf_total")
        wf = molden_io.from_molden(a.molden, a.cart_norm)
    if a.qvec and a.polar_grid:
        sys.exit("give either --qvec or --polar-grid, not both")
    if a.qvec or a.polar_grid:
        return _aligned(a, wf, t0)
    q = np.linspace(a.qmin, a.qmax, a.npoints)
    opts = dict(screen_tol=a.screen_tol, cutoffcentre=a.cutoffcentre, l=a.legendre)
    if a.probe == "electron":
        from . import electron
        if a.kind == "two_electron":
            sys.exit("--kind two_electron is an X-ray quantity; use total, elastic, inelastic or decompose")
        d = electron.decompose(wf, q, hf=(a.kind == "hf_total"), nthreads=a.threads, **opts)
        kind = "total" if a.kind == "hf_total" else a.kind
        if a.kind == "decompose":
            cols = [c for c in ("total", "elastic", "inelastic", "nn", "ne", "ee_elastic", "ee_total") if c in d]
        else:
            cols = [kind]
        data = np.column_stack([q] + [d[c] for c in cols])
        cols = ["electron_" + c for c in cols]
    elif a.kind == "decompose":
        d = core.decompose(wf, q, nthreads=a.threads, **opts)
        cols = ["total", "elastic", "inelastic", "hartree", "exchange", "cumulant"]
        data = np.column_stack([q] + [d[c] for c in cols])
    else:
        cols = [a.kind]
        data = np.column_stack([q, core.scattering(wf, q, a.kind, nthreads=a.threads, **opts)])
    header = (f"pyxscat  molden={a.molden}" + (f"  molcas={a.molcas_log} root={a.root}" if a.molcas_log else "")
              + f"  N={wf.nelec:g}  legendre l={a.legendre}  probe={a.probe}\nq " + " ".join(cols))
    np.savetxt(a.output, data, header=header)
    print(f"{a.output} written ({a.npoints} q points, {time.time() - t0:.1f} s)")


def _aligned(a, wf, t0):
    from . import aligned
    from . import grids
    if a.polar_grid:
        e0, ths, n1, n2 = a.polar_grid
        qlab = grids.qcoord3(e0, ths, int(n1), int(n2))
    else:
        qlab = grids.read_qgrid(a.qvec)
    qmol = qlab if a.orientation is None else aligned.rotate_q(qlab, aligned.rotation_matrix(*a.orientation))
    if a.probe == "electron":
        from . import electron
        d = electron.aligned(wf, qmol, hf=(a.kind == "hf_total"), nthreads=a.threads)
        cols = [c for c in ("total", "elastic", "inelastic") if c in d] if a.kind == "decompose" else \
            ["total" if a.kind == "hf_total" else a.kind]
        data = np.column_stack([qlab] + [d[c] for c in cols])
        cols = ["electron_" + c for c in cols]
    elif a.kind == "decompose":
        d = aligned.decompose(wf, qmol, nthreads=a.threads)
        cols = ["total", "elastic", "inelastic", "hartree", "exchange", "cumulant"]
        data = np.column_stack([qlab] + [d[c] for c in cols])
    else:
        cols = [a.kind]
        data = np.column_stack([qlab, aligned.scattering(wf, qmol, a.kind, nthreads=a.threads)])
    header = (f"pyxscat (fixed orientation)  molden={a.molden}"
              + (f"  euler={a.orientation}" if a.orientation is not None else "") + f"  N={wf.nelec:g}\n"
              + "qx qy qz " + " ".join(cols))
    np.savetxt(a.output, data, header=header)
    print(f"{a.output} written ({len(qlab)} q vectors, {time.time() - t0:.1f} s)")


if __name__ == "__main__":
    main()
