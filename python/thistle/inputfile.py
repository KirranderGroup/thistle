# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Plain-text input files and the ``thistle`` command.

    thistle new water.toml        # write a commented template to fill in
    thistle run water.toml        # compute everything it asks for, write CSV + PNG
    thistle check orbitals.molden # is this file read correctly?
    thistle submit water.toml     # write (and optionally send) a Slurm job script
    thistle plot results/water_scattering.csv [--iam results/water_iam.csv]   # (re)draw a plot
"""
import argparse
import os
import sys
import time

try:
    import tomllib
except ModuleNotFoundError:                       # Python < 3.11
    import tomli as tomllib

TEMPLATE = '''# THISTLE input file.  Lines starting with # are comments.
# Run with:   thistle run {name}

[input]
# What to read.  THISTLE recognises the file type by itself:
#   a Molden file, an OpenMolcas RASSCF output, a Molpro MULTI output or punch file.
# For Molcas/Molpro outputs the Molden file of the same run is found automatically;
# give it here only if it has a different name.
file  = "my_calculation.out"
# molden = "my_calculation.molden"
# punch  = "my_calculation.pun"   # Molpro: the CI vectors from the punch file (no output needed)
state = 1                 # which root / state (1 = lowest)
# hf = true               # a Molden file from closed-shell Hartree-Fock: also total and inelastic
                          # (without it a Molden file gives the elastic scattering only)
# factored = "auto"       # correlated 2RDM as products of 1-particle matrices (no n^4 array):
                          # "auto" beyond 60 basis functions, true always, false never

[q]
unit   = "1/A"            # "1/A" (inverse Angstrom) or "1/bohr"
max    = 4.0
points = 201

[compute]
probe      = "xray"       # "xray", "electron" (UED / gas electron diffraction) or "both"
# energy_kev = 3000       # electron beam energy: relativistic factor and detector wavelength
scattering = true         # rotationally averaged total, elastic, inelastic
parts      = "standard"   # "all" adds Hartree, exchange and correlation parts
iam        = true         # independent atom model, the % difference from it (and sM(s) for electrons);
                          # also on the [detector] and [alignment] grids
anisotropy = false        # j2 moment about the molecular z axis (aligned molecules)
# iam_basis = "aug-cc-pvtz"                        # basis of the free atoms, or per element:
# iam_basis = { Li = "cc-pvtz", F = "aug-cc-pvtz" }  # (use the molecule's basis)

# [detector]              # uncomment for one fixed orientation of the molecule
# orientation = [0, 90, 0] # Euler angles in degrees (z-y-z); the beam runs along z
# grid        = "square"   # "square": an image in (qx, qy);  "polar": the rings of grid.f90 (qcoord3);
#                          # or a q-grid file: columns qx qy qz in bohr^-1, e.g. qgrid.dat of grid.f90
# pixels      = 101        # square grid: pixels per side
# theta_max   = 60         # polar grid: largest scattering angle (degrees; X-rays 60 by default,
#                          #   electrons: the angle that reaches the [q] max)
# rings       = 49         #   number of rings (n1 of grid.f90)
# azimuths    = 72         #   points per ring (n2 - 1 of grid.f90)
# photon_kev  = 10         # X-ray photon energy (or wavelength = 1.0, Angstrom); electrons use energy_kev

# [alignment]             # partially aligned molecules (laser alignment): from the j2 moment,
#                         #   I(q) = I0(|q|) + 5 <P2> P2(q^.axis) S2(|q|)
# axis = [0, 1, 0]         # laboratory direction of the alignment (e.g. the laser polarisation)
# p2   = 0.5               # degree of alignment <P2(cos beta)> of the molecular z axis (1 = perfect)
# grid = "polar"           # the detector grid; grid, theta_max, rings, azimuths, photon_kev, pixels
#                          # as in [detector] (taken from there when not given here)

[output]
prefix = "results/{stem}"   # files: <prefix>_scattering.csv, .png, ...
plots  = true

[run]
threads = 0               # 0 = use all cores of this machine

# Only for `thistle submit` (clusters with Slurm):
[cluster]
time      = "02:00:00"
cpus      = 16
memory    = "16G"
# partition = "standard"
# account   = "my_project"
# setup     = ["module load python", "source ~/thistle-env/bin/activate"]  # run before thistle
'''


def template(path):
    """Write a commented input file to fill in."""
    if os.path.exists(path):
        raise FileExistsError(f"{path} already exists")
    text = TEMPLATE.replace("{name}", os.path.basename(path))
    text = text.replace("{stem}", os.path.splitext(os.path.basename(path))[0])
    with open(path, "w") as f:
        f.write(text)
    return path


def _read(path):
    """Settings from an input file, or defaults for a calculation file given directly."""
    if path.lower().endswith(".toml"):
        with open(path, "rb") as f:
            cfg = tomllib.load(f)
    else:                                         # `thistle run lif.out`: no input file needed
        stem = os.path.splitext(os.path.basename(path))[0]
        cfg = {"input": {"file": path}, "compute": {"scattering": True, "iam": False},
               "output": {"prefix": stem}}
        return cfg
    base = os.path.dirname(os.path.abspath(path))
    inp = cfg.get("input", {})
    if "file" not in inp:
        raise ValueError(f"{path}: [input] needs file = '...'")
    for k in ("file", "molden", "punch"):
        if k in inp and not os.path.isabs(inp[k]):
            inp[k] = os.path.join(base, inp[k])
    for sec in ("detector", "alignment"):
        g = cfg.get(sec, {}).get("grid")
        if isinstance(g, str) and g.lower() not in ("square", "polar") and not os.path.isabs(g):
            cfg[sec]["grid"] = os.path.join(base, g)
    out = cfg.setdefault("output", {})
    out.setdefault("prefix", os.path.splitext(path)[0])
    if not os.path.isabs(out["prefix"]):
        out["prefix"] = os.path.join(base, out["prefix"])
    return cfg


class _Log:
    """Steps with their wall time, the files written with a line on each, and the plot script."""

    def __init__(self, say):
        self.say = say
        self.steps = []
        self.files = {}                 # key -> path
        self.notes = {}                 # path -> description
        self.script = []                # lines of <prefix>_plots.py
        self.plot_time, self.nplots = 0.0, 0
        self.t0 = time.perf_counter()
        self.started = time.strftime("%Y-%m-%d %H:%M:%S")

    def step(self, label):
        log = self

        class _Timer:
            def __enter__(self):
                self.t = time.perf_counter()
                log.say(f"[{time.strftime('%H:%M:%S')}] {label} ...", flush=True)

            def __exit__(self, *exc):
                dt = time.perf_counter() - self.t
                log.steps.append((label, dt))
                if exc[0] is None:
                    log.say(f"           done in {_fmt_time(dt)}", flush=True)
        return _Timer()

    def file(self, key, path, note):
        self.files[key] = path
        self.notes[path] = note
        return path

    def plot(self, comment, code, png):
        self.script += ["", f"# {comment}", *code]


def _fmt_time(t):
    if t < 60:
        return f"{t:.1f} s"
    m, s_ = divmod(t, 60)
    if m < 60:
        return f"{int(m)} min {s_:.0f} s"
    h, m = divmod(m, 60)
    return f"{int(h)} h {int(m)} min"


def run(path, verbose=True, **overrides):
    """Do everything an input file asks for.  Returns {name: written file}.

    ``path`` may also be a calculation file (Molden, Molcas/Molpro output): then the default
    settings are used.  ``overrides`` (state, qmax, points, unit, parts, iam, anisotropy, prefix)
    replace the corresponding settings; ``detector={...}`` and ``alignment={...}`` add/update those sections.

    Besides the data (CSV, .dat, .npz) and plots (PNG), every run writes ``<prefix>_plots.py``,
    a short script that redraws every plot from the data files (edit it to change the plots),
    and ``<prefix>_run.json`` with the settings, the time of each step and the file list.
    """
    from . import core, plotting
    say = print if verbose else (lambda *a, **k: None)
    log = _Log(say)
    cfg = _read(path)
    _apply(cfg, overrides)
    inp, qc, comp, out = cfg["input"], cfg.get("q", {}), cfg.get("compute", {}), cfg["output"]
    threads = cfg.get("run", {}).get("threads", 0) or None
    say(f"THISTLE {_version()}: {path}  ({time.strftime('%Y-%m-%d %H:%M:%S')}, "
        f"{threads or _max_threads()} threads)")
    with log.step("read the wavefunction"):
        mol = core.load(inp["file"], molden=inp.get("molden"), state=inp.get("state", 1), punch=inp.get("punch"),
                        hf=inp.get("hf", False), factored=_factored(inp.get("factored", "auto")))
    say(mol)
    folder = os.path.dirname(os.path.abspath(out["prefix"]))
    os.makedirs(folder, exist_ok=True)
    grid = dict(qmax=qc.get("max", 4.0), points=qc.get("points", 201), unit=qc.get("unit", "1/A"))
    plots = out.get("plots", True)
    probes = _probes(comp.get("probe", "xray"))
    energy = comp.get("energy_kev")
    iam_kw = dict(iam_basis=comp.get("iam_basis", "aug-cc-pvtz"), iam_source=comp.get("iam_source", "ab initio"))
    ib = iam_kw["iam_basis"]
    iam_atoms = ", ".join(f"{k} {v}" for k, v in ib.items()) if isinstance(ib, dict) else ib
    rel = lambda p: os.path.basename(p)                     # names inside the plot script

    def draw(fn, png, name, *args, note="plot", **kw):
        if plots:
            t = time.perf_counter()
            fn(*args, path=png, **kw)
            log.plot_time += time.perf_counter() - t
            log.nplots += 1
            log.file(name + "_plot", png, note)

    for probe in probes:
        # X-ray files keep the plain names; electron files carry "_electron"
        tag = "" if probe == "xray" else "_electron"
        key = "" if probe == "xray" else "electron_"
        what = "X-ray" if probe == "xray" else "electron"
        pre = out["prefix"] + tag
        pk = dict(probe=probe) if probe == "xray" else dict(probe=probe, energy_kev=energy)

        def csv_(name, result, note):
            p = log.file(key + name, result.save(f"{pre}_{name}.csv"), note)
            return p

        res = ref = None
        if comp.get("scattering", True):
            with log.step(f"{what} scattering, rotational average"):
                res = mol.scattering(parts=comp.get("parts", "standard"), threads=threads, **grid, **pk)
            say(res)
            f = csv_("scattering", res, f"{what} scattering: {', '.join(res.curves)} ({len(res.q)} q points)")
            png = f"{pre}_scattering.png"
            draw(plotting.scattering, png, key + "scattering", res)
            log.plot(f"{what} scattering", [f'res = tp.read(here("{rel(f)}"))',
                                            f'tp.scattering(res, path=here("{rel(png)}"))'
                                            '   # curves=["total", "elastic"], yscale="log"'], png)
        if comp.get("iam", False):
            with log.step(f"{what} independent atom model"):
                ref = mol.iam(**grid, basis=iam_kw["iam_basis"], source=iam_kw["iam_source"], **pk)
            fi = csv_("iam", ref, f"{what} IAM on the same q grid (free atoms: {iam_atoms})")
            if res is not None:
                curve = "total" if ("total" in res.curves and "total" in ref.curves) else "elastic"
                csv_("vs_iam", core.compare(res, ref, curve), f"100 (I - I_IAM) / I_IAM, {curve}")
                png = f"{pre}_vs_iam.png"
                draw(plotting.versus_iam, png, key + "vs_iam", res, ref, curve=curve, note="plot: molecule and IAM, and the % difference")
                log.plot(f"{what}: molecule against the IAM", [
                    f'res, iam = tp.read(here("{rel(f)}")), tp.read(here("{rel(fi)}"))',
                    f'tp.versus_iam(res, iam, curve="{curve}", path=here("{rel(png)}"))'], png)
                if probe == "electron":
                    smr = core.sm(res, ref, curve)
                    fs = csv_("sM", smr, "modified scattering sM(s): molecule and IAM")
                    png = f"{pre}_sM.png"
                    draw(plotting.sm, png, key + "sM", smr)
                    log.plot("electrons: modified scattering intensity sM(s)", [
                        f'tp.sm(tp.read(here("{rel(fs)}")), path=here("{rel(png)}"))'], png)
        if comp.get("anisotropy", False):
            with log.step(f"{what} j2 moment (aligned molecules)"):
                s2 = mol.anisotropy(threads=threads, **grid, **pk)
            fa = csv_("anisotropy", s2, f"{what} j2 moment S2(q) about the molecular z axis")
            png = f"{pre}_anisotropy.png"
            draw(plotting.anisotropy, png, key + "anisotropy", s2)
            log.plot(f"{what} j2 moment", [f'tp.anisotropy(tp.read(here("{rel(fa)}")), path=here("{rel(png)}"))'],
                     png)
        for sec, label in (("detector", "fixed orientation"), ("alignment", "aligned ensemble")):
            if sec not in cfg:
                continue
            sc = cfg[sec]
            name = "detector" if sec == "detector" else "aligned"
            stem = f"{pre}_{name}"
            qmax_sec = sc.get("qmax", grid["qmax"] if probe == "electron" else None)
            if sec == "detector":
                kw = dict(orientation=tuple(sc.get("orientation", (0, 0, 0))), unit=grid["unit"],
                          qmax=qmax_sec, threads=threads, **_lab(sc, {}, probe, energy))
                fn = mol.detector
            else:
                kw = dict(axis=tuple(sc.get("axis", (0, 0, 1))), p2=sc.get("p2", 1.0), unit=grid["unit"],
                          qmax=qmax_sec, threads=threads,
                          **_lab(sc, cfg.get("detector", {}), probe, energy, default_grid="polar"))
                fn = mol.aligned_ensemble
            with log.step(f"{what} {label} ({kw['grid'] if isinstance(kw['grid'], str) else 'q grid'})"):
                pat = fn(**kw)
            fd = _save_pattern(pat, stem, key + name, log, f"{what} {label}: {', '.join(_curves(pat))} "
                               f"on {_npoints(pat)} points")
            png = stem + ".png"
            draw(plotting.detector, png, key + name, pat)
            log.plot(f"{what} {label}", [f'tp.detector(tp.read(here("{rel(fd)}")), path=here("{rel(png)}"))'
                                         '   # curve="elastic", log=False'], png)
            if comp.get("iam", False):
                with log.step(f"{what} {label}, IAM"):
                    ipat = fn(model="iam", **kw, **iam_kw)
                fdi = _save_pattern(ipat, stem + "_iam", key + name + "_iam", log, f"{what} {label}, IAM")
                _save_pattern(core.compare(pat, ipat, "total"), stem + "_vs_iam", key + name + "_vs_iam", log,
                              "100 (I - I_IAM) / I_IAM")
                png = stem + "_vs_iam.png"
                draw(plotting.detector_versus_iam, png, key + name + "_vs_iam", pat, ipat, note="plot: molecule, IAM, % difference")
                log.plot(f"{what} {label}: molecule against the IAM", [
                    f'pat, iam = tp.read(here("{rel(fd)}")), tp.read(here("{rel(fdi)}"))',
                    f'tp.detector_versus_iam(pat, iam, path=here("{rel(png)}"))'], png)
    if log.nplots:
        log.steps.append((f"plots ({log.nplots} figures)", log.plot_time))
    _write_extras(cfg, path, mol, out["prefix"], log)
    _report(log, folder)
    return log.files


def _version():
    from . import __version__
    return __version__


def _max_threads():
    from pyxscat._lib import lib
    return int(lib.pyxscat_max_threads())


def _curves(pat):
    return list(pat.images if hasattr(pat, "images") else pat.values)


def _npoints(pat):
    if hasattr(pat, "images"):
        v = next(iter(pat.images.values()))
        return f"{v.shape[0]} x {v.shape[1]} pixels"
    return f"{len(pat.qvec)}"


def _write_extras(cfg, path, mol, prefix, log):
    """<prefix>_plots.py (redraws every plot) and <prefix>_run.json (settings, timings, files)."""
    import json
    script = f"{prefix}_plots.py"
    head = [f'"""Plots of the THISTLE run {os.path.basename(str(path))} ({time.strftime("%Y-%m-%d %H:%M")}).',
            "",
            f"    python {os.path.basename(script)}",
            "",
            "redraws every plot from the data files in this folder.  Change the calls below (curves,",
            "scales, titles) or work with the arrays: res.q, res[\"total\"], pat.qvec, pat[\"total\"] ...",
            "Each tp.* function returns the matplotlib Figure (fig.axes[0].set_xlim(0, 2), ...).",
            '"""',
            "import os",
            "",
            "import thistle.plotting as tp",
            "",
            "HERE = os.path.dirname(os.path.abspath(__file__))",
            "",
            "",
            "def here(name):",
            "    return os.path.join(HERE, name)",
            ""]
    with open(script, "w") as f:
        f.write("\n".join(head + log.script + ["", 'print("plots redrawn in", HERE)', ""]))
    log.file("plot_script", script, f"redraws every plot from these files: python {os.path.basename(script)}")
    info = {"thistle": _version(), "input": os.path.abspath(str(path)), "settings": cfg,
            "wavefunction": str(mol).splitlines()[0], "started": log.started,
            "wall_time_s": round(time.perf_counter() - log.t0, 3),
            "steps": [{"step": k, "seconds": round(t, 3)} for k, t in log.steps],
            "files": {k: os.path.abspath(v) for k, v in log.files.items()}}
    runj = f"{prefix}_run.json"
    log.file("run_log", runj, "settings, time of each step, file list")
    info["files"]["run_log"] = os.path.abspath(runj)
    with open(runj, "w") as f:
        json.dump(info, f, indent=1, default=str)


def _report(log, folder):
    say = log.say
    total = time.perf_counter() - log.t0
    say("")
    say(f"Finished in {_fmt_time(total)}.")
    w = max(len(k) for k, _ in log.steps) if log.steps else 10
    for k, t in log.steps:
        say(f"  {k:<{w}}  {_fmt_time(t):>12}")
    say("")
    paths = list(dict.fromkeys(log.files.values()))
    folders = {os.path.dirname(os.path.abspath(p)) for p in paths}
    say(f"Files written in {folder}{os.sep}" if folders == {folder} else "Files written:")
    w = max(len(os.path.basename(p)) for p in paths)
    for p in paths:
        shown = os.path.basename(p) if folders == {folder} else os.path.abspath(p)
        say(f"  {shown:<{w}}  {log.notes.get(p, '')}")


def _lab(sec, fallback, probe, energy, default_grid="square"):
    """Beam and grid keywords of Molecule.detector / aligned_ensemble from an input-file section
    (missing keys from ``fallback``, e.g. [alignment] from [detector])."""
    get = lambda k, d=None: sec.get(k, fallback.get(k, d))
    kw = dict(grid=get("grid", default_grid), pixels=get("pixels", 101), theta_max=get("theta_max"),
              rings=get("rings", 49), azimuths=get("azimuths", 72), probe=probe)
    if probe == "xray":
        kw.update(wavelength=get("wavelength"), energy_kev=get("photon_kev"))
    else:                                 # electrons: the beam energy of [compute] unless given here
        kw.update(energy_kev=get("energy_kev", energy))
    return kw


def _save_pattern(pat, stem, key, log, note):
    """A square image as .npz, a pattern on a polar/file grid as text (qx qy qz |q| curves...)."""
    images = getattr(pat, "images", None)
    return log.file(key, pat.save(stem + (".npz" if images is not None else ".dat")), note)


def _probes(value):
    from .core import _probe
    if isinstance(value, str) and value.lower() == "both":
        return ["xray", "electron"]
    vals = [value] if isinstance(value, str) else list(value)
    return [_probe(v) for v in vals]


_WHERE = {"state": ("input", "state"), "molden": ("input", "molden"), "punch": ("input", "punch"), "hf": ("input", "hf"), "factored": ("input", "factored"), "qmax": ("q", "max"),
          "points": ("q", "points"), "unit": ("q", "unit"), "parts": ("compute", "parts"),
          "iam": ("compute", "iam"), "probe": ("compute", "probe"), "energy_kev": ("compute", "energy_kev"), "anisotropy": ("compute", "anisotropy"), "prefix": ("output", "prefix"),
          "plots": ("output", "plots"), "threads": ("run", "threads")}


def _factored(v):
    """"auto" | true | false (also as strings from the command line)."""
    if isinstance(v, str):
        k = v.strip().lower()
        if k in ("auto", ""):
            return "auto"
        if k in ("true", "yes", "on", "1"):
            return True
        if k in ("false", "no", "off", "0"):
            return False
        raise ValueError(f"factored must be auto, true or false (got {v!r})")
    return bool(v)


def _apply(cfg, overrides):
    for k, v in overrides.items():
        if v is None:
            continue
        if k in ("detector", "alignment"):         # whole sections (from the command line)
            cfg.setdefault(k, {}).update(v)
            continue
        if k not in _WHERE:
            raise TypeError(f"unknown setting {k!r}; known: {', '.join(_WHERE)}")
        sec, key = _WHERE[k]
        cfg.setdefault(sec, {})[key] = v


SLURM = """#!/bin/bash
#SBATCH --job-name={job}
#SBATCH --time={time}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={memory}
#SBATCH --output={job}.%j.log
{extra}
cd {folder}
{setup}
export OMP_NUM_THREADS={cpus}
thistle run {input}
"""


def slurm_script(path):
    cfg = _read(path)
    cl = cfg.get("cluster", {})
    extra = "\n".join(f"#SBATCH --{k}={cl[k]}" for k in ("partition", "account", "qos") if k in cl)
    job = os.path.splitext(os.path.basename(path))[0]
    setup = cl.get("setup", [])
    setup = "\n".join([setup] if isinstance(setup, str) else setup)
    return SLURM.format(job=job, time=cl.get("time", "02:00:00"), cpus=cl.get("cpus", 16),
                        memory=cl.get("memory", "16G"), extra=extra, setup=setup,
                        folder=os.path.dirname(os.path.abspath(path)), input=os.path.basename(path))


def main(argv=None):
    p = argparse.ArgumentParser(prog="thistle", description="X-ray scattering from wavefunctions.")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("new", help="write a commented input file to fill in")
    s.add_argument("file")
    s = sub.add_parser("run", help="run an input file (.toml), or a calculation file with default settings")
    s.add_argument("file")
    s.add_argument("--state", type=int, help="root / state number (1 = lowest)")
    s.add_argument("--molden", help="Molden file, if not found automatically")
    s.add_argument("--punch", help="Molpro punch file with the CI vectors (e.g. with a Molden file as FILE)")
    s.add_argument("--hf", "--scf", dest="hf", action="store_const", const=True,
                   help="the Molden file is from closed-shell Hartree-Fock: also total and inelastic")
    s.add_argument("--factored", choices=["auto", "true", "false"],
                   help="correlated 2RDM as products of one-particle matrices (default auto: beyond 60 "
                        "basis functions)")
    s.add_argument("--qmax", type=float)
    s.add_argument("--points", type=int)
    s.add_argument("--unit", choices=["1/A", "1/bohr"])
    s.add_argument("--all-parts", dest="parts", action="store_const", const="all",
                   help="also Hartree, exchange and correlation parts")
    s.add_argument("--iam", action="store_const", const=True, help="also the independent atom model")
    s.add_argument("--electrons", dest="probe", action="store_const", const="electron",
                   help="electron scattering (UED/GED) instead of X-rays")
    s.add_argument("--probe", dest="probe", choices=["xray", "electron", "both"], help="xray (default), electron or both")
    s.add_argument("--energy", dest="energy_kev", type=float,
                   help="electron beam energy in keV (adds the relativistic factor; detector wavelength)")
    s.add_argument("--anisotropy", action="store_const", const=True, help="also the j2 moment")
    g = s.add_argument_group("fixed orientation and aligned molecules (detector grids)")
    g.add_argument("--orientation", type=float, nargs=3, metavar=("A", "B", "G"),
                   help="one fixed orientation: Euler angles in degrees (z-y-z); writes <prefix>_detector")
    g.add_argument("--align", type=float, nargs=3, metavar=("X", "Y", "Z"),
                   help="partially aligned ensemble along this laboratory axis; writes <prefix>_aligned")
    g.add_argument("--p2", type=float, help="degree of alignment <P2(cos beta)> (default 1)")
    g.add_argument("--grid", help="square, polar (the grid of grid.f90) or a q-grid file (qx qy qz, bohr^-1); "
                                  "default: square for --orientation, polar for --align")
    g.add_argument("--photon-kev", type=float, help="X-ray photon energy in keV (default: 1 Angstrom)")
    g.add_argument("--theta-max", type=float, help="polar grid: largest scattering angle, degrees (60)")
    g.add_argument("--rings", type=int, help="polar grid: number of rings (49)")
    g.add_argument("--azimuths", type=int, help="polar grid: points per ring (72)")
    s.add_argument("--out", dest="prefix", help="prefix of the output files")
    s = sub.add_parser("plot", help="plot result files written by thistle run (CSV, .dat, .npz)")
    s.add_argument("files", nargs="+")
    s.add_argument("-o", "--out", help="output image (default: next to the first file, .png); .pdf/.svg work too")
    s.add_argument("--iam", help="IAM file on the same grid: molecule vs IAM with the %% difference")
    s.add_argument("--curve", help="the curve to show (overlay of several files, detector, vs IAM)")
    s.add_argument("--curves", nargs="+", help="curves of a scattering file to show")
    s.add_argument("--log", dest="yscale", action="store_const", const="log", help="logarithmic intensity axis")
    s.add_argument("--linear", dest="yscale", action="store_const", const="linear", help="linear axis")
    s = sub.add_parser("check", help="check that an input file (Molden, output) is read correctly")
    s.add_argument("file")
    s.add_argument("--scf-energy", type=float, help="the program's SCF energy, for a Molden convention check")
    s = sub.add_parser("selftest", help="check this installation (seconds; send the output with bug reports)")
    s.add_argument("files", nargs="*", help="also (instead) check sum rules for these calculation files")
    s = sub.add_parser("submit", help="write a Slurm script for an input file (and send it with --now)")
    s.add_argument("file")
    s.add_argument("--now", action="store_true", help="submit with sbatch straight away")
    a = p.parse_args(argv)

    if a.cmd == "new":
        template(a.file)
        print(f"wrote {a.file}: fill in [input] file = ..., then `thistle run {a.file}`")
    elif a.cmd == "run":
        over = {k: getattr(a, k) for k in _WHERE if hasattr(a, k)}
        lab = {k: getattr(a, k) for k in ("grid", "photon_kev", "theta_max", "rings", "azimuths")
               if getattr(a, k) is not None}
        if a.orientation is not None or (a.grid and a.align is None):
            over["detector"] = dict(lab, **({"orientation": a.orientation} if a.orientation is not None else {}))
        if a.align is not None:
            over["alignment"] = dict(lab, axis=a.align, **({"p2": a.p2} if a.p2 is not None else {}))
        run(a.file, **over)
    elif a.cmd == "plot":
        _plot_command(a)
    elif a.cmd == "selftest":
        from pyxscat import selftest
        sys.exit(selftest.main(a.files))
    elif a.cmd == "check":
        _check(a.file, a.scf_energy)
    elif a.cmd == "submit":
        script = os.path.splitext(a.file)[0] + ".slurm"
        with open(script, "w") as f:
            f.write(slurm_script(a.file))
        if a.now:
            import subprocess
            r = subprocess.run(["sbatch", script], capture_output=True, text=True)
            print(r.stdout or r.stderr)
        else:
            print(f"wrote {script}; send it with `sbatch {script}` (or `thistle submit {a.file} --now`)")


def _plot_command(a):
    from . import plotting as tp
    objs = [tp.read(f) for f in a.files]
    out = a.out or os.path.splitext(a.files[0])[0] + ("_vs_iam" if a.iam else "_overlay" if len(objs) > 1 else "") + ".png"
    first = objs[0]
    if a.iam:
        ref = tp.read(a.iam)
        if hasattr(first, "curves"):
            tp.versus_iam(first, ref, curve=a.curve, path=out)
        else:
            tp.detector_versus_iam(first, ref, curve=a.curve or "total", path=out)
    elif len(objs) > 1:
        curve = a.curve or ("total" if "total" in first.curves else next(iter(first.curves)))
        tp.overlay(objs, curve, labels=[os.path.splitext(os.path.basename(f))[0] for f in a.files], path=out,
                   yscale=a.yscale)
    elif hasattr(first, "curves") and (a.curves or a.yscale):
        tp.scattering(first, curves=a.curves, yscale=a.yscale, path=out)
    elif not hasattr(first, "curves"):
        tp.detector(first, curve=a.curve, path=out, log=a.yscale != "linear")
    else:
        tp.plot_file(a.files[0], out=out)
    print(f"wrote {os.path.abspath(out)}")


def _check(path, scf_energy=None):
    from . import core
    from pyxscat.interfaces import molden as mld
    mol = core.load(path)
    print(mol)
    if core._sniff(path) == "molden":
        data = mld.read_molden(path)
        print(f"  MO orthonormality error {mld.check_orbitals(data):.1e} "
              f"(Cartesian normalisation: {data.cart_norm})")
        if scf_energy is not None:
            e, _ = mld.scf_energy(path)
            print(f"  SCF energy from the file {e:.9f}, program {scf_energy:.9f}, difference {e - scf_energy:.1e}")
    import pyxscat as xs
    n = mol.wf.nelec
    el0 = float(xs.elastic(mol.wf.basis, [0.0], mol.wf.dm1)[0])
    ok = abs(el0 - n * n) < 1e-4 * n * n
    print(f"  elastic scattering at q = 0: {el0:.6f} (must be N^2 = {n * n:g})  {'OK' if ok else 'PROBLEM'}")
    if mol.correlated:
        tot0 = float(xs.scattering(mol.wf, [0.0], "total")[0])
        ok2 = abs(tot0 - n * n) < 1e-4 * n * n
        print(f"  total scattering at q = 0:   {tot0:.6f} (must be N^2)  {'OK' if ok2 else 'PROBLEM'}")


if __name__ == "__main__":
    main()
