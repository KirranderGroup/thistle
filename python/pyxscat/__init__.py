# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""PyXSCAT: ab initio X-ray and electron scattering from wavefunctions.

    import pyxscat as xs
    wf = xs.interfaces.pyscf.from_cc(mycc)
    total = xs.scattering(wf, q, "total")

Authors: Andrés Moreno Carrascosa, Mats Simmermacher, Patrick Yuheng Wang, Eirik Moen Liane, Joseph Cooper, Nikola Zotev, Thomas Northey and Adam Kirrander.
"""
from .basis import Basis
from .core import (Wavefunction, scattering, decompose, two_electron, elastic, exchange,
                   mean_field, symmetrize_dm2)
from . import interfaces, transforms, ci, aligned, iam, electron

__version__ = "0.2.0.dev0"
__author__ = "A. Moreno Carrascosa, M. Simmermacher, P. Y. Wang, E. Moen Liane, J. Cooper, N. Zotev, T. Northey, A. Kirrander"
__all__ = ["Basis", "Wavefunction", "scattering", "decompose", "two_electron", "elastic",
           "exchange", "mean_field", "symmetrize_dm2", "interfaces", "transforms", "ci", "aligned", "iam", "electron"]
