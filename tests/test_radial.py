"""The radial functions behind every j0 kernel: sum_b c_b j_b(qh) (q/h)^b = sum_b c_b q^(2b) g_b(qh),
g_b(x) = j_b(x)/x^b.  Each order must be accurate for every x: small x (series and downward
recurrence), orders above x (Miller's downward recurrence; the upward recurrence is unstable
there) and large x (upward recurrence)."""
import math

import numpy as np
import pytest

from pyxscat._lib import lib


def _radial(x, order, b):
    """g_b(x) through the library: unit coefficient on order b, h = 1 so that q = x."""
    x = np.asfortranarray(np.atleast_1d(np.asarray(x, dtype=float)))
    coef = np.zeros(order + 1, order="F")
    coef[b] = 1.0
    out = np.zeros(x.size, order="F")
    lib.pyxscat_radial_series(x.size, x, 1.0, order, coef, out)
    return out / x ** (2 * b)


def _reference(b, x):
    from scipy.special import spherical_jn
    if x >= 0.1:
        return spherical_jn(b, x) / x ** b
    s, term, k = 0.0, 1.0 / math.prod(range(1, 2 * b + 2, 2)), 0
    while abs(term) > 1e-20 * max(abs(s), 1e-300):
        s += term
        k += 1
        term *= -x * x / (2.0 * k * (2 * b + 2 * k + 1))
    return s


def _envelope(b, x):
    """|g_b| without its zeros: the small-x value or the large-x amplitude, whichever is smaller."""
    return min(1.0 / math.prod(range(1, 2 * b + 2, 2)), x ** -(b + 1))


@pytest.mark.parametrize("order", [2, 6, 10, 14, 18])
def test_every_order_every_x(order):
    xs = np.logspace(-4, 3, 281)
    worst = 0.0
    for b in range(order + 1):
        got = _radial(xs, order, b)
        for x, g in zip(xs, got):
            ref = _reference(b, x)
            worst = max(worst, abs(g - ref) / max(abs(ref), _envelope(b, x)))
    assert worst < 1e-12, worst


def test_orders_above_x():
    """The region that breaks the upward recurrence: 1 < x < order."""
    order = 16
    xs = np.linspace(1.05, 15.9, 60)
    for b in (8, 12, 16):
        got = _radial(xs, order, b)
        ref = np.array([_reference(b, x) for x in xs])
        env = np.array([_envelope(b, x) for x in xs])
        assert np.max(np.abs(got - ref) / np.maximum(np.abs(ref), env)) < 1e-12
