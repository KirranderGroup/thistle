# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Time and peak memory of the CI -> active-space RDM step, one case per process.

    python builders.py <norb> <nalpha> <nbeta> <method> [ndet]

Full-CAS vectors are random (the cost does not depend on the values); with ``ndet`` a random
subset of that many determinants is kept (a truncated / selected expansion).
Prints one JSON line: time, peak RSS, sizes.
"""
import json, resource, sys, time
from math import comb
import numpy as np
from pyscf.fci import cistring
from pyxscat import ci

norb, na, nb, method = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
ndet = int(sys.argv[5]) if len(sys.argv) > 5 else None
rng = np.random.default_rng(1)
sa = cistring.make_strings(range(norb), na)
sb = cistring.make_strings(range(norb), nb)
if ndet is None:
    ia, ib = np.divmod(np.arange(sa.size * sb.size), sb.size)
else:
    pick = rng.choice(sa.size * sb.size, ndet, replace=False)
    ia, ib = np.divmod(pick, sb.size)
dets = np.column_stack([sa[ia], sb[ib]]).astype(np.int64)
c = rng.standard_normal(len(dets))
rss0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
t = time.perf_counter()
dm1, dm2 = ci.make_rdm12(norb, na, nb, [tuple(map(int, d)) for d in dets] if method == "dense" else dets,
                         c, method=method)
t = time.perf_counter() - t
rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print(json.dumps({"norb": norb, "nel": na + nb, "ndet": len(dets), "nfci": comb(norb, na) * comb(norb, nb),
                  "method": method, "seconds": round(t, 2), "peak_GB": round(rss / 1e6, 2),
                  "before_GB": round(rss0 / 1e6, 2), "trace": float(np.trace(dm1)),
                  "dm2_MB": dm2.nbytes / 1e6}))
