import ctypes

import numpy as np
from numpy.polynomial.hermite import hermgauss

from pyxscat._lib import lib


def test_hermite_coefficients_reproduce_overlaps():
    """sqrt(pi/p) K_ab E(0; i, j) is the 1D overlap of two Cartesian Gaussians, up to h functions."""
    rng = np.random.default_rng(3)
    ng, maxl = 5, 5
    lg = np.full(ng, maxl, dtype=np.int32, order="F")
    xg = np.asfortranarray(rng.normal(size=ng))
    ag = np.asfortranarray(rng.uniform(0.3, 3.0, size=ng))
    dd = np.zeros((2 * maxl + 1, maxl + 1, maxl + 1, ng, ng), order="F")
    lib.pyxscat_hermite_table(ng, maxl, lg, xg, ag, dd)

    t, w = hermgauss(40)
    for a in range(ng):
        for b in range(ng):
            p = ag[a] + ag[b]
            P = (ag[a] * xg[a] + ag[b] * xg[b]) / p
            K = np.exp(-ag[a] * ag[b] / p * (xg[a] - xg[b]) ** 2)
            x = P + t / np.sqrt(p)
            for i in range(maxl + 1):
                for j in range(maxl + 1):
                    ref = K * np.sum(w * (x - xg[a]) ** i * (x - xg[b]) ** j) / np.sqrt(p)
                    got = np.sqrt(np.pi / p) * K * dd[0, i, j, b, a]
                    assert abs(got - ref) < 1e-11 * max(1.0, abs(ref))
