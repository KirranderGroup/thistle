"""Checks for a new machine or compiler, beyond `thistle selftest`.

    python benchmarks/cluster_check.py reference ifx.npz      # compute a set of curves, save them
    python benchmarks/cluster_check.py compare ifx.npz gnu.npz  # two builds (or machines) agree?
    python benchmarks/cluster_check.py scaling                # timings against the thread count

`reference` covers every kernel and both probes: the Molpro LiF wavefunction of the test data
(X-ray total/elastic/inelastic and all parts, j2 moment, electrons, a detector image, the
polar grid of grid.f90 and an aligned ensemble) and a PySCF
CCSD water calculation (2RDM through the quartet kernel).  Run it from the unpacked source folder.
"""
import os
import sys
import time
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "tests", "data", "molpro", "lif_sa4_casscf.out")


def reference(path):
    import pyxscat as xs
    import thistle as th
    from pyxscat.selftest import build_info
    warnings.simplefilter("ignore")
    out = {}
    t0 = time.time()
    lif = th.load(DATA)
    q = np.linspace(0.0, 8.0, 41)
    for k, v in lif.scattering(q=q, parts="all").curves.items():
        out[f"lif_xray_{k}"] = v
    for k, v in lif.anisotropy(q=q[1:]).curves.items():
        out[f"lif_j2_{k}"] = v
    for k, v in lif.scattering(q=q[1:], probe="electron", parts="all").curves.items():
        out[f"lif_electron_{k}"] = v
    out["lif_detector_total"] = lif.detector(pixels=21, orientation=(10, 70, 30))["total"]
    out["lif_polar_total"] = lif.detector(grid="polar", energy_kev=10, rings=8, azimuths=12,
                                          orientation=(10, 70, 30))["total"]
    out["lif_aligned_total"] = lif.aligned_ensemble(axis=(0, 1, 0), p2=0.5, rings=8, azimuths=12,
                                                    energy_kev=10)["total"]
    from pyscf import gto, scf, cc
    mol = gto.M(atom="O 0.05 0.02 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", basis="cc-pvdz", verbose=0)
    mycc = cc.CCSD(scf.RHF(mol).run(conv_tol=1e-12)).run(conv_tol=1e-10)
    wf = xs.interfaces.pyscf.from_cc(mycc)
    d = xs.decompose(wf, q)
    for k in ("total", "elastic", "hartree", "exchange", "cumulant"):
        out[f"h2o_ccsd_{k}"] = d[k]
    out["h2o_ccsd_total_j2"] = xs.scattering(wf, q[1:], "total", l=2)
    wff = xs.interfaces.pyscf.from_cc(mycc, factored=True)             # the same without the 4-index array
    out["h2o_ccsd_total_factored"] = xs.scattering(wff, q, "total")
    info = build_info()
    np.savez(path, **out, info=str(info))
    print(f"{path}: {len(out)} curves in {time.time() - t0:.0f} s, built with {info['compiler'].split('|')[0]}")


def compare(a, b):
    A, B = np.load(a), np.load(b)
    print(f"A: {str(A['info'])[:160]}\nB: {str(B['info'])[:160]}\n")
    worst = 0.0
    for k in sorted(set(A.files) - {"info"}):
        x, y = A[k], B[k]
        scale = np.nanmax(np.abs(y))
        err = np.nanmax(np.abs(x - y)) / scale if scale > 0 else np.nanmax(np.abs(x - y))
        worst = max(worst, err)
        print(f"  {k:34s} max |A - B| / max |B| = {err:8.1e}  {'OK' if err < 1e-8 else 'DIFFERENT'}")
    print(f"\nworst {worst:.1e}: {'the two builds agree' if worst < 1e-8 else 'the builds DISAGREE - please report'}")


def scaling():
    import pyxscat as xs
    from pyxscat._lib import lib
    from pyscf import gto, scf
    warnings.simplefilter("ignore")
    bz = """C 0.000 1.396 0; C 1.209 0.698 0; C 1.209 -0.698 0; C 0 -1.396 0; C -1.209 -0.698 0;
            C -1.209 0.698 0; H 0 2.479 0; H 2.147 1.240 0; H 2.147 -1.240 0; H 0 -2.479 0;
            H -2.147 -1.240 0; H -2.147 1.240 0"""
    mol = gto.M(atom=bz, basis="cc-pvdz", verbose=0)
    mf = scf.RHF(mol).run()
    wf = xs.interfaces.pyscf.from_scf(mf)
    q = np.linspace(0.0, 10.0, 100)
    nmax = int(lib.pyxscat_max_threads())
    print(f"benzene/cc-pVDZ ({wf.basis.nao} functions), 100 q points; up to {nmax} threads")
    print(f"{'threads':>8} {'elastic (pairs)':>16} {'HF total (quartets)':>20}")
    base = None
    for nt in sorted({1, 2, 4, 8, 16, 32, 64, nmax}):
        if nt > nmax:
            continue
        t = time.time(); xs.elastic(wf.basis, q, wf.dm1, nthreads=nt); te = time.time() - t
        t = time.time(); xs.mean_field(wf.basis, q, wf.dm1, nthreads=nt); tm = time.time() - t
        base = base or (te, tm)
        print(f"{nt:8d} {te:9.2f} s ({base[0] / te:4.1f}x) {tm:11.2f} s ({base[1] / tm:4.1f}x)", flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "reference" and len(sys.argv) == 3:
        reference(sys.argv[2])
    elif cmd == "compare" and len(sys.argv) == 4:
        compare(sys.argv[2], sys.argv[3])
    elif cmd == "scaling":
        scaling()
    else:
        print(__doc__)
