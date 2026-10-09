# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""THISTLE: X-ray and electron scattering from wavefunctions, the simple way.

    import thistle as th

    wf  = th.load("water.molden")                         # any supported input, auto-detected
    res = wf.scattering(qmax=4)                           # q in 1/Angstrom by default
    res.plot("water.png"); res.save("water.csv")

    th.compare(wf.scattering(), wf.iam()).plot()          # molecule vs independent atom model
    e = wf.scattering(probe="electron")                   # electron diffraction cross section
    th.sm(e, wf.iam(probe="electron")).plot()             # modified scattering sM(s)
    wf.detector(wavelength=1.0).plot("image.png")         # fixed orientation, detector image
    wf.detector(grid="polar", energy_kev=10)              # ... on the polar grid of grid.f90
    wf.aligned_ensemble(axis=(0, 1, 0), p2=0.5)           # partially aligned molecules (j2)

Everything here is a thin layer over :mod:`pyxscat`; the full control (cutoffs, density
providers, raw arrays) stays available there.

Authors: Andrés Moreno Carrascosa, Mats Simmermacher, Patrick Yuheng Wang, Eirik Moen Liane, Joseph Cooper, Nikola Zotev, Thomas Northey and Adam Kirrander.
"""
from .core import (load, Molecule, Result, DetectorImage, DetectorPattern, compare, sm, electron_wavelength, BOHR, Q_UNITS)
from .inputfile import run, template
from pyxscat.grids import qcoord3, polar_grid, read_qgrid, write_qgrid

__all__ = ["load", "Molecule", "Result", "DetectorImage", "DetectorPattern", "compare", "sm", "electron_wavelength", "run", "template",
           "BOHR", "Q_UNITS", "qcoord3", "polar_grid", "read_qgrid", "write_qgrid"]
__version__ = "0.2.0.dev0"
__author__ = "A. Moreno Carrascosa, M. Simmermacher, P. Y. Wang, E. Moen Liane, J. Cooper, N. Zotev, T. Northey, A. Kirrander"
