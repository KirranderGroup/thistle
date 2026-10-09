# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Real-space transforms of scattering curves."""
import numpy as np


def sine_transform(q, sq, r):
    """(2/pi) int q S(q) sin(q r) dq, by the trapezoidal rule on the given q grid."""
    q = np.asarray(q, dtype=float)
    sq = np.asarray(sq, dtype=float)
    r = np.atleast_1d(np.asarray(r, dtype=float))
    integrand = q[None, :] * sq[None, :] * np.sin(q[None, :] * r[:, None])
    trap = getattr(np, "trapezoid", None) or np.trapz          # numpy < 2 has only trapz
    return (2.0 / np.pi) * trap(integrand, q, axis=1)


def to_real_space(q, sq, r):
    """r times the sine transform: the radial pair-distribution-like function."""
    return np.asarray(r, dtype=float) * sine_transform(q, sq, r)
