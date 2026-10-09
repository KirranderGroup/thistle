# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
import numpy as np, pyxscat as xs
from pyxscat import aligned as al
from scipy.integrate import lebedev_rule
from load_deck import load
wf, a, f = load((__import__("sys").argv[1] if len(__import__("sys").argv) > 1 else "../../..") + "/tests/decks/ch4_total/deck.h5")
q = np.array([0.5, 1.0, 2.0, 3.0, 4.0, 4.8, 6.0])
s2 = xs.two_electron(wf.basis, q, wf.dm2, l=2)
xyz, w = lebedev_rule(131); w = w / w.sum(); n = xyz.T
p2 = 1.5 * n[:, 2] ** 2 - 0.5
ref = np.array([-(w @ (al.two_electron(wf.basis, qq * n, wf.dm2) * p2)) for qq in q])
i0 = xs.two_electron(wf.basis, q, wf.dm2)
print(" q      S2 two-electron (j2 kernel)   independent check (Lebedev)   rel diff    I2 j0 (for scale)")
for r in zip(q, s2, ref, np.abs(s2 - ref) / np.abs(ref), i0):
    print(" %.1f   %16.8e          %16.8e        %.1e    %10.4f" % r)
print("Eirik's Y20 convention = ours x sqrt(5/(4 pi)) =", np.sqrt(5 / (4 * np.pi)))
