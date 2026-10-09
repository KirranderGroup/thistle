# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
import numpy as np, pyxscat as xs
from load_deck import load
import sys
D = (sys.argv[1] if len(sys.argv) > 1 else "../../..") + "/tests/decks/"
PIFAC = 1.0 + 8.3483e-8
for d, kind in (("ch4_total", "two"), ("ch4_elastic", "el")):
    wf, a, f = load(D + d + "/deck.h5")
    gold = np.loadtxt(D + d + "/golden.dat"); q = gold[:, 0]
    keep = np.abs(gold[:, 1]) > 1e-8 * np.abs(gold[:, 1]).max()
    for cc in (1e-10, 1e-2):
        if kind == "two":
            ours = xs.two_electron(wf.basis, q, wf.dm2, screen_tol=1e-16, cutoffcentre=cc)
        else:
            ours = xs.elastic(wf.basis, q, wf.dm1, screen_tol=1e-16, cutoffcentre=cc, method="quartets")
        rel = np.abs(ours * PIFAC - gold[:, 1])[keep] / np.abs(gold[:, 1][keep])
        i = np.argmax(rel)
        print(f"{d:12s} cutoffcentre={cc:7.0e}: max rel diff {rel.max():.2e} (at q = {q[keep][i]:.2f})")
