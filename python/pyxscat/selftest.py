# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""A one-minute check of an installation: every kernel against exact sum rules and against an
independent route, on water/6-31G (PySCF), plus the build information needed in a bug report.

    thistle selftest            (or: python -m pyxscat.selftest)
"""
import ctypes
import os
import platform
import sys
import time

import numpy as np


def build_info():
    from ._lib import lib, _find_library
    buf = ctypes.create_string_buffer(2048)
    try:
        lib.pyxscat_build_info(buf, 2048)
        comp = buf.value.decode(errors="replace")
    except AttributeError:
        comp = "unknown (library predates pyxscat_build_info)"
    return {"library": _find_library(), "compiler": comp, "threads": int(lib.pyxscat_max_threads()),
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", "unset"),
            "python": sys.version.split()[0], "numpy": np.__version__, "platform": platform.platform()}


def _checks():
    from pyscf import gto, scf
    from . import aligned as al, core, electron
    from .interfaces import pyscf as xp

    mol = gto.M(atom="O 0.05 0.02 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", basis="6-31g", verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-11)
    wf = xp.from_scf(mf, with_dm2=True)                       # explicit Hartree-Fock 2RDM as well
    b, d, n = wf.basis, wf.dm1, float(wf.nelec)
    q0 = np.array([0.0])
    q = np.array([0.5, 1.5, 4.0])
    rel = lambda a, r: float(np.max(np.abs(np.asarray(a) - r) / np.maximum(np.abs(r), 1e-300)))

    def lebedev(fn, l=0):
        from scipy.integrate import lebedev_rule
        xyz, w = lebedev_rule(131)
        nn, w = xyz.T, w / w.sum()
        wt = np.ones(len(w)) if l == 0 else (1.5 * nn[:, 2] ** 2 - 0.5)
        return np.array([w @ (fn(qq * nn) * wt) for qq in q])

    yield "elastic I(0) = N^2 (pair route)", rel(core.elastic(b, q0, d), n * n), 1e-10
    yield "elastic I(0) = N^2 (quartet route)", rel(core.elastic(b, q0, d, method="quartets"), n * n), 1e-8
    yield "explicit 2RDM: I2(0) = N(N-1)", rel(core.two_electron(b, q0, wf.dm2), n * (n - 1)), 1e-8
    yield "mean field: I2(0) = N(N-1)", rel(core.mean_field(b, q0, d), n * (n - 1)), 1e-8
    yield "exchange(0) = -N", rel(core.exchange(b, q0, d), -n), 1e-8
    ref_el = lebedev(lambda qv: np.abs(al.form_factor(b, qv, d)) ** 2)
    yield "elastic, pair route vs fixed-orientation average", rel(core.elastic(b, q, d), ref_el), 1e-9
    yield "elastic, quartet route vs fixed-orientation average", rel(core.elastic(b, q, d, method="quartets"),
                                                                       ref_el), 1e-8
    ref_ex = lebedev(lambda qv: al.exchange(b, qv, d))
    yield "exchange vs fixed-orientation average", rel(core.exchange(b, q, d), ref_ex), 1e-8
    mfq = core.mean_field(b, q, d)
    yield "explicit 2RDM (mode 1) vs on-the-fly mean field (mode 2)", rel(core.two_electron(b, q, wf.dm2), mfq), 1e-9
    yield "one thread vs all threads (2RDM)", rel(core.two_electron(b, q, wf.dm2, nthreads=1),
                                                  core.two_electron(b, q, wf.dm2)), 1e-11
    ref2 = lebedev(lambda qv: np.abs(al.form_factor(b, qv, d)) ** 2, l=2)
    yield "j2 elastic, pair route vs projection", rel(core.elastic(b, q, d, l=2), ref2), 1e-9
    yield "j2 elastic, quartet route vs projection", rel(core.elastic(b, q, d, l=2, method="quartets"), ref2), 1e-8
    from . import ci
    rng = np.random.default_rng(7)
    strs = [sum(1 << k for k in c) for c in __import__("itertools").combinations(range(8), 4)]
    dets = [(a, b) for a in strs for b in strs]
    cv = rng.standard_normal(len(dets)) * (rng.random(len(dets)) < 0.3)
    dd = ci.make_rdm12(8, 4, 4, dets, cv, method="dense")[1]
    yield "2RDM from determinants: bitwise (Fortran) vs dense", float(np.abs(ci.make_rdm12(8, 4, 4, dets, cv, method="bitwise")[1] - dd).max()), 1e-12
    yield "2RDM from determinants: sparse vs dense", float(np.abs(ci.make_rdm12(8, 4, 4, dets, cv, method="sparse")[1] - dd).max()), 1e-12
    z = np.asarray(wf.atom_numbers, float)
    yield "electron-nucleus(0) = Z N", rel(electron.electron_nucleus(b, q0, d, z, wf.coords), z.sum() * n), 1e-10


def _file_checks(path, state=1):
    """Sum rules for a user's own calculation, separating the Python side (reading, RDMs) from
    the Fortran kernels."""
    import thistle as th
    from . import core
    from .interfaces.molden import overlap_matrix
    mol = th.load(path, state=state)
    wf = mol.wf
    n = float(wf.nelec)
    s = overlap_matrix(wf.basis)
    rel = lambda a, r: float(abs(a - r) / abs(r))
    yield f"[{os.path.basename(path)}] Python: tr(D S) = N", rel(np.einsum("ij,ji->", wf.dm1, s), n), 1e-6
    yield f"[{os.path.basename(path)}] Fortran: elastic I(0) = N^2", rel(core.elastic(wf.basis, [0.0], wf.dm1)[0], n * n), 1e-6
    if wf.dm2 is not None:
        tr2 = np.einsum("ijkl,ij,kl->", wf.dm2, s, s)
        yield f"[{os.path.basename(path)}] Python: sum Gamma S S = N(N-1)", rel(tr2, n * (n - 1)), 1e-6
        i2 = core.two_electron(wf.basis, [0.0], wf.dm2)[0]
        yield f"[{os.path.basename(path)}] Fortran: I2(0) = sum Gamma S S", rel(i2, tr2), 1e-8
        i2b = core.two_electron(wf.basis, [0.0], wf.dm2, nthreads=1)[0]
        yield f"[{os.path.basename(path)}] Fortran: I2(0), one thread vs all", rel(i2b, i2), 1e-11


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    print("THISTLE / pyxscat self-test")
    for k, v in build_info().items():
        print(f"  {k:16s} {v}")
    print()
    t0 = time.time()
    failed = 0
    try:
        checks = list(_checks()) if not argv else []
        for path in argv:
            checks += list(_file_checks(path))
        for name, err, tol in checks:
            ok = err <= tol
            failed += not ok
            print(f"  {'PASS' if ok else 'FAIL'}  {name:58s} {err:9.1e}  (tolerance {tol:.0e})", flush=True)
    except ImportError as exc:
        print(f"  cannot run the checks: {exc} (install with pip install '.[test]')")
        return 2
    print(f"\n{'all checks passed' if not failed else f'{failed} check(s) FAILED'} in {time.time() - t0:.0f} s")
    if failed:
        print("Please send this whole output with the bug report.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
