# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Plots for every THISTLE calculation, from results in memory or from the files that
``thistle run`` writes.

    import thistle as th
    import thistle.plotting as tp

    res = mol.scattering()                       # or tp.read("results/lif_scattering.csv")
    tp.scattering(res, path="lif.png")           # total, elastic, inelastic (and parts)
    tp.versus_iam(res, mol.iam())                # molecule and IAM, and the % difference below
    tp.sm(th.sm(e, mol.iam(probe="electron")))   # electrons: sM(s)
    tp.anisotropy(mol.anisotropy())              # j2 moment
    tp.detector(mol.detector(grid="polar"))      # a detector image or a pattern on any grid
    tp.detector_versus_iam(pat, iam_pat)         # molecule, IAM, % difference
    tp.overlay([res1, res2], "total", labels=["state 1", "state 2"])
    tp.plot_file("results/lif_vs_iam.csv")       # the right plot for any file THISTLE wrote

Every function returns the matplotlib Figure (``fig.axes`` for further changes) and saves it
when ``path`` is given.  Curves keep the same colour in every plot (total blue, elastic orange,
inelastic green, ...), from a palette checked for colour-blind readers.
"""
import json
import os

import numpy as np

# Fixed colour per quantity (colour follows the entity, never its position in a plot).
# Categorical slots of the validated default palette, in order.
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
COLOURS = {
    "total": SLOTS[0], "elastic": SLOTS[1], "inelastic": SLOTS[2],
    "hartree": SLOTS[3], "exchange": SLOTS[4], "cumulant": SLOTS[5],
    "nucleus-nucleus": SLOTS[6], "electron-nucleus": SLOTS[7], "electron-electron": SLOTS[3],
    "atomic": SLOTS[6], "molecular": SLOTS[5], "interference": SLOTS[5],
}
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
REFERENCE = dict(color=INK2, ls=(0, (5, 3)), lw=1.6)          # the IAM / reference curve



# ----------------------------------------------------------------------------- helpers
def _base(name):
    """The quantity behind a curve label: 'S2 total', '%diff total', 'sM (total)' -> 'total'."""
    n = name.lower().replace("%diff", "").replace("s2 ", "").replace("sm (", "").replace(")", "")
    return n.replace("electron_", "").strip()


def colour(name, i=0):
    """The fixed colour of a curve (unknown names take the free slots in order)."""
    b = _base(name)
    if b in COLOURS:
        return COLOURS[b]
    free = [c for c in SLOTS if c not in COLOURS.values()] or SLOTS
    return free[i % len(free)]


def _figure(ncols=1, nrows=1, size=(6.4, 4.0), path=None, **kw):
    """A figure with the THISTLE style; without pyplot when it only goes to a file."""
    if path is not None:
        from matplotlib.figure import Figure
        fig = Figure(figsize=size, layout="constrained")
        axes = fig.subplots(nrows, ncols, squeeze=False, **kw)
    else:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(nrows, ncols, figsize=size, squeeze=False, layout="constrained", **kw)
    for ax in axes.flat:
        _style(ax)
    return fig, axes


def _style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelcolor=INK2, labelsize=8.5)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(INK2)
    ax.yaxis.label.set_color(INK2)
    ax.title.set_color(INK)


def _axes(ax, path):
    """The given axes (and its figure), or a new single-panel figure."""
    if ax is not None:
        _style(ax)
        return ax.figure, ax
    fig, axes = _figure(path=path)
    return fig, axes[0, 0]


def _finish(fig, path):
    if path:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fig.savefig(path, dpi=160, facecolor="white")
    return fig


def _legend(ax, n):
    if n >= 2:
        ax.legend(frameon=False, loc="best", handlelength=2.4)


def _yscale(ax, ys, scale):
    if scale == "log":
        ax.set_yscale("log")
    elif scale == "symlog":
        vals = np.concatenate([np.abs(y)[np.isfinite(y)] for y in ys]) if ys else np.zeros(0)
        if vals.size and vals.max() > 0:
            ax.set_yscale("symlog", linthresh=1e-3 * vals.max())


def _zero(ax):
    ax.axhline(0.0, color=AXIS, lw=1.0, zorder=1)


def _ylabel(res, default="intensity (electron units)"):
    return res.meta.get("ylabel", default)


def _title(ax, text, sub=None):
    """Bold title over the axes; the provenance line goes in the figure's footer."""
    ax.set_title(text, loc="left", fontweight="bold", fontsize=10.5)
    if sub:
        _footer(ax.figure, sub)


def _footer(fig, text):
    import textwrap
    if getattr(fig, "_thistle_footer", None) is not None:
        return
    lines = textwrap.wrap(text, int(fig.get_figwidth() * 17))
    fig._thistle_footer = fig.text(0.01, 0.008, "\n".join(lines), ha="left", va="bottom", color=MUTED,
                                   fontsize=7.5)
    fig._thistle_footer.set_in_layout(False)
    eng = fig.get_layout_engine()
    if eng is not None and hasattr(eng, "set"):              # keep the footer clear of the axes
        bottom = (len(lines) * 7.5 * 1.35 + 8.0) / 72.0 / fig.get_figheight()
        eng.set(rect=(0.0, bottom, 1.0, 1.0 - bottom))


def _u(unit):
    return {"1/A": "\u00c5$^{-1}$", "1/bohr": "bohr$^{-1}$"}.get(unit, unit)


def _subtitle(meta):
    """One line on where the numbers came from (file, state, probe, beam)."""
    src = meta.get("from", {}) or {}
    bits = []
    if src.get("source"):
        bits.append(os.path.basename(str(src["source"])))
    if src.get("state") not in (None, 1) or src.get("kind") in ("molpro_output", "molpro_punch", "molcas_output"):
        if src.get("state") is not None:
            bits.append(f"state {src['state']}")
    if meta.get("probe") == "electron":
        bits.append("electrons" + (f", {meta['energy_kev']:g} keV" if meta.get("energy_kev") else ""))
    elif meta.get("probe") == "xray":
        bits.append("X-rays")
    return "  ·  ".join(bits) or None


# ----------------------------------------------------------------------------- curves
def scattering(res, curves=None, path=None, yscale=None, ax=None, title=None):
    """Rotationally averaged curves (any :class:`thistle.Result`): one line per quantity."""
    fig, ax = _axes(ax, path)
    names = list(curves or res.meta.get("plot_curves") or res.curves)
    drawn = []
    for i, k in enumerate(names):
        same = [d for d in drawn if np.allclose(res.curves[k], res.curves[d], rtol=1e-9, atol=0, equal_nan=True)]
        if same:                                  # identical curves (Hartree = elastic, ...): one line
            for line in ax.get_lines():
                if line.get_label() == same[0] or line.get_label().startswith(same[0] + " ="):
                    line.set_label(line.get_label() + f" = {k}")
            continue
        ax.plot(res.q, res.curves[k], color=colour(k, i), label=k, zorder=3)
        drawn.append(k)
    ax.set_xlabel(f"q ({_u(res.unit)})")
    ax.set_ylabel(_ylabel(res))
    _yscale(ax, [res.curves[k] for k in names], yscale or res.meta.get("yscale"))
    if any(np.nanmin(res.curves[k]) < 0 for k in names) and ax.get_yscale() != "log":
        _zero(ax)
    _legend(ax, len(names))
    _title(ax, title or res.name, _subtitle(res.meta))
    return _finish(fig, path)


def anisotropy(res, curves=None, path=None, ax=None, title=None):
    """The j2 moment S2(q) (sign matters: the zero line is drawn)."""
    fig = scattering(res, curves=curves, path=None, ax=ax, title=title or f"{res.name}: j₂ moment S₂(q)")
    _zero(fig.axes[0])
    return _finish(fig, path)


def versus_iam(res, iam, curve=None, path=None, title=None):
    """Molecule and independent atom model on the same grid (top) and the percent difference
    100 (I - I_IAM) / I_IAM (bottom).  ``curve``: total if both have it, else elastic."""
    if curve is None:
        curve = "total" if ("total" in res.curves and "total" in iam.curves) else "elastic"
    if not np.allclose(res.q, iam.q):
        raise ValueError("the two results have different q grids")
    fig, axes = _figure(nrows=2, size=(6.4, 5.6), path=path, sharex=True,
                        gridspec_kw={"height_ratios": [1.6, 1.0]})
    top, bot = axes[0, 0], axes[1, 0]
    c = colour(curve)
    top.plot(res.q, res[curve], color=c, label=f"{res.name or 'molecule'}", zorder=3)
    top.plot(iam.q, iam[curve], label="IAM", zorder=2, **REFERENCE)
    top.set_ylabel(_ylabel(res).replace("intensity", f"{curve} intensity"))
    if res.meta.get("yscale") == "log":
        top.set_yscale("log")
    _legend(top, 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        pd = 100.0 * (res[curve] - iam[curve]) / iam[curve]
    bot.plot(res.q, pd, color=c, zorder=3)
    _zero(bot)
    bot.set_ylabel("100 (I − I$_\\mathrm{IAM}$) / I$_\\mathrm{IAM}$ (%)")
    bot.set_xlabel(f"q ({_u(res.unit)})")
    _title(top, title or f"{res.name}: {curve} scattering vs independent atoms", _subtitle(res.meta))
    return _finish(fig, path)


def difference(res, path=None, ax=None, title=None):
    """A percent-difference curve as written by ``compare`` (``*_vs_iam.csv``)."""
    fig = scattering(res, path=None, ax=ax, title=title or res.name)
    ax = fig.axes[0]
    _zero(ax)
    ax.set_ylabel(res.meta.get("ylabel", "% difference"))
    return _finish(fig, path)


def sm(res, path=None, ax=None, title=None):
    """Modified scattering intensity sM(s) of electron diffraction: molecule (solid) and the IAM
    simulation (dashed)."""
    fig, ax = _axes(ax, path)
    for i, k in enumerate(res.curves):
        if "iam" in k.lower():
            ax.plot(res.q, res[k], label=k, zorder=2, **REFERENCE)
        else:
            ax.plot(res.q, res[k], color=colour(k, i), label=k, zorder=3)
    _zero(ax)
    ax.set_xlabel(f"s ({_u(res.unit)})")
    ax.set_ylabel(f"sM(s) ({_u(res.unit)})")
    _legend(ax, len(res.curves))
    _title(ax, title or res.name, _subtitle(res.meta))
    return _finish(fig, path)


def overlay(results, curve="total", labels=None, path=None, title=None, yscale=None):
    """One curve from several results (states, methods, basis sets) on the same axes; each
    result takes the next palette slot."""
    fig, axes = _figure(path=path)
    ax = axes[0, 0]
    labels = labels or [r.name for r in results]
    for i, (r, lab) in enumerate(zip(results, labels)):
        ax.plot(r.q, r[curve], color=SLOTS[i % len(SLOTS)], label=lab, zorder=3)
    ax.set_xlabel(f"q ({_u(results[0].unit)})")
    ax.set_ylabel(_ylabel(results[0]))
    _yscale(ax, [r[curve] for r in results], yscale or results[0].meta.get("yscale"))
    _legend(ax, len(results))
    _title(ax, title or curve)
    return _finish(fig, path)


# ----------------------------------------------------------------------------- detector
def _xy(obj):
    from .core import Q_UNITS
    if hasattr(obj, "images"):
        return obj.qx, obj.qy, obj.images
    f = Q_UNITS[obj.unit]
    return obj.qvec[:, 0] * f, obj.qvec[:, 1] * f, obj.values


def _cbar_label(obj, curve):
    if "%" in curve:
        return "100 (I − I$_\\mathrm{IAM}$) / I$_\\mathrm{IAM}$ (%)"
    if obj.meta.get("probe") == "electron":
        return f"{curve}: dσ/dΩ ({'Å²' if obj.unit == '1/A' else 'bohr²'}/sr)"
    return f"{curve} (electron units)"


def _draw(ax, obj, curve, norm=None, cmap=None, log=True):
    from .core import _colour
    x, y, vals = _xy(obj)
    v = np.asarray(vals[curve], float)
    if cmap is None:
        cmap, norm0 = _colour(v, log)
        norm = norm or norm0
        fin = v[np.isfinite(v)]
        if norm0 is not None and type(norm0).__name__ == "LogNorm" and fin.size and fin.max() > 1e3 * fin.min():
            from matplotlib.colors import LogNorm   # a few points next to q = 0 (electrons: 1/q^4)
            norm = LogNorm(np.percentile(fin, 0.5), np.percentile(fin, 99.0))   # would take the whole scale
    if hasattr(obj, "images"):
        ext = [x.min(), x.max(), y.min(), y.max()]
        im = ax.imshow(v, origin="lower", extent=ext, cmap=cmap, norm=norm, interpolation="nearest")
    else:
        tri = _triangulation(obj, x, y, v)
        im = ax.tripcolor(tri, np.where(np.isfinite(v), v, np.nanmin(v)), shading="gouraud", cmap=cmap, norm=norm)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel(f"q$_x$ ({_u(obj.unit)})")
    ax.set_ylabel(f"q$_y$ ({_u(obj.unit)})")
    return im


def _triangulation(obj, x, y, v):
    """Triangles over the detector points.  The polar grid is triangulated along its own rings
    (points on a ring are co-circular, so a Delaunay triangulation would be arbitrary and change
    with the last digit); other grids by Delaunay on rounded coordinates.  Triangles touching a
    point without a value (q = 0 for electrons) are masked."""
    from matplotlib.tri import Triangulation
    n = len(x)
    rings, naz = obj.meta.get("rings"), obj.meta.get("azimuths")
    if obj.meta.get("grid") == "polar" and rings and naz and n == rings * naz + 1:
        idx = lambda i, j: 1 + (i * naz) + (j % naz)              # ring i (0-based), azimuth j
        t = [(0, idx(0, j), idx(0, j + 1)) for j in range(naz)]
        for i in range(rings - 1):
            for j in range(naz):
                a, b, c, d = idx(i, j), idx(i, j + 1), idx(i + 1, j), idx(i + 1, j + 1)
                t += [(a, c, d), (a, d, b)]
        tri = Triangulation(x, y, np.array(t))
    else:
        tri = Triangulation(np.round(x, 9), np.round(y, 9))
    bad = ~np.isfinite(v)
    if bad.any():
        tri.set_mask(bad[tri.triangles].any(axis=1))
    return tri


def detector(obj, curve=None, path=None, log=True, title=None):
    """A detector image (square grid) or pattern (polar grid, q-grid file): log colour scale for
    intensities, a diverging scale for signed data (differences)."""
    vals = obj.images if hasattr(obj, "images") else obj.values
    curve = curve or ("total" if "total" in vals else next(iter(vals)))
    fig, axes = _figure(size=(5.6, 4.9), path=path)
    ax = axes[0, 0]
    im = _draw(ax, obj, curve, log=log)
    fig.colorbar(im, ax=ax, label=_cbar_label(obj, curve), shrink=0.9)
    _title(ax, title or f"{obj.name}: {curve}", _detector_sub(obj.meta))
    return _finish(fig, path)


def _detector_sub(meta):
    bits = [_subtitle(meta) or ""]
    if meta.get("kind", "").startswith("aligned"):
        axis = ", ".join(f"{float(v):g}" for v in meta.get("alignment_axis", []))
        bits.append(f"aligned along ({axis}), <P\u2082> = {float(meta.get('p2', 1)):g}")
    elif meta.get("orientation_deg") is not None:
        eu = ", ".join(f"{float(v):g}" for v in meta["orientation_deg"])
        bits.append(f"Euler angles ({eu})\u00b0")
    if meta.get("grid") == "polar":
        bits.append(f"polar grid {meta.get('rings')} \u00d7 {meta.get('azimuths')}, "
                    f"2\u03b8 \u2264 {float(meta.get('theta_max_deg') or 0):.3g}\u00b0")
    return "  \u00b7  ".join(b for b in bits if b)


def detector_versus_iam(obj, iam, curve="total", path=None, title=None):
    """Three panels on the same grid: molecule, IAM (one shared log scale) and the percent
    difference (diverging scale centred on zero)."""
    from matplotlib.colors import LogNorm, Normalize
    va = obj.images if hasattr(obj, "images") else obj.values
    vb = iam.images if hasattr(iam, "images") else iam.values
    if curve not in va or curve not in vb:
        curve = "elastic"
    a, b = np.asarray(va[curve], float), np.asarray(vb[curve], float)
    with np.errstate(divide="ignore", invalid="ignore"):
        d = 100.0 * (a - b) / b
    fig, axes = _figure(ncols=3, size=(13.0, 4.3), path=path)
    both = np.concatenate([a[np.isfinite(a)].ravel(), b[np.isfinite(b)].ravel()])
    norm = LogNorm(both.min(), both.max()) if both.size and both.min() > 0 else None
    for ax, o, lab in ((axes[0, 0], obj, obj.name or "molecule"), (axes[0, 1], iam, "IAM")):
        im = _draw(ax, o, curve, norm=norm, cmap="viridis")
        _title(ax, f"{lab}: {curve}")
    fig.colorbar(im, ax=axes[0, :2].tolist(), label=_cbar_label(obj, curve), shrink=0.9)
    diff = type(obj)(**{**_fields(obj), ("images" if hasattr(obj, "images") else "values"): {"%diff": d}})
    m = float(np.nanmax(np.abs(d))) or 1.0
    im = _draw(axes[0, 2], diff, "%diff", norm=Normalize(-m, m), cmap="RdBu_r")
    fig.colorbar(im, ax=axes[0, 2], label=_cbar_label(obj, "%"), shrink=0.9)
    _title(axes[0, 2], "difference from IAM")
    fig.suptitle(title or _detector_sub(obj.meta), color=MUTED, fontsize=9)
    return _finish(fig, path)


def _fields(obj):
    from dataclasses import fields
    return {f.name: getattr(obj, f.name) for f in fields(obj)}


# ----------------------------------------------------------------------------- files
def read(path):
    """Any file ``thistle run`` writes, back as the object that wrote it: a ``.csv`` gives a
    :class:`~thistle.Result`, a ``.npz`` a :class:`~thistle.DetectorImage` (or pattern), a
    ``.dat`` a :class:`~thistle.DetectorPattern`."""
    from .core import Result, DetectorImage, DetectorPattern
    path = os.fspath(path)
    if path.endswith(".csv"):
        with open(path) as f:
            first = f.readline()
            header = f.readline().strip().split(",")
        name, unit, meta = _parse_first(first)
        data = np.atleast_2d(np.genfromtxt(path, delimiter=",", skip_header=2))
        unit = unit or header[0][header[0].find("(") + 1:header[0].rfind(")")]
        return Result(data[:, 0], unit, {k: data[:, i + 1] for i, k in enumerate(header[1:])},
                      name=name, meta=meta)
    if path.endswith(".npz"):
        z = np.load(path, allow_pickle=False)
        meta = json.loads(str(z["meta"])) if "meta" in z.files else {}
        unit = meta.get("q_unit", "1/A")
        if "qx" in z.files:
            imgs = {k: z[k] for k in z.files if k not in ("qx", "qy", "meta")}
            return DetectorImage(z["qx"], z["qy"], imgs, unit, name=meta.get("name", ""), meta=meta)
        vals = {k: z[k] for k in z.files if k not in ("qvec", "meta")}
        return DetectorPattern(z["qvec"], vals, unit, name=meta.get("name", ""), meta=meta)
    # text pattern: '# name  {json}' / '# qx qy qz |q| (bohr^-1)  curve curve ...'
    with open(path) as f:
        first, second = f.readline(), f.readline()
    name, _, meta = _parse_first(first)
    cols = second.lstrip("#").split("(bohr^-1)")[-1].split()
    data = np.atleast_2d(np.loadtxt(path))
    return DetectorPattern(data[:, :3], {k: data[:, 4 + i] for i, k in enumerate(cols)},
                           meta.get("q_unit", "1/A"), name=name, meta=meta)


def _parse_first(line):
    line = line.lstrip("#").strip()
    j = line.find("{")
    meta = json.loads(line[j:]) if j >= 0 else {}
    head = line[:j].strip() if j >= 0 else line
    unit = None
    if "  q in " in head:
        head, unit = head.split("  q in ", 1)
        unit = unit.strip()
    return head.strip(), unit, meta


def plot_file(path, out=None, **kw):
    """The natural plot of any file THISTLE wrote (saved next to it as .png unless ``out``)."""
    obj = read(path)
    out = out or os.path.splitext(path)[0] + ".png"
    if hasattr(obj, "curves"):
        lab = obj.meta.get("ylabel", "")
        names = " ".join(obj.curves)
        if "sM" in names:
            return sm(obj, path=out, **kw)
        if "%" in lab or "%diff" in names:
            return difference(obj, path=out, **kw)
        if "j2" in obj.meta.get("kind", "") or names.startswith("S2"):
            return anisotropy(obj, path=out, **kw)
        return scattering(obj, path=out, **kw)
    return detector(obj, path=out, **kw)
