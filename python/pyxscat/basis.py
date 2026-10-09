# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Package-independent description of a Gaussian basis set.

The scattering kernels only ever see this object.  Every contracted basis function is a
linear combination of primitive Cartesian Gaussians::

    chi_mu(r) = sum_p  coeff_p (x-X)^l (y-Y)^m (z-Z)^n exp(-ga |r-R|^2),

where ``coeff`` already contains the contraction coefficient and all normalisation factors of
the electronic structure package the basis came from.  Density matrices handed to the kernels
must be expressed in the functions chi_mu, whatever their normalisation.
"""
from dataclasses import dataclass

import numpy as np


def cartesian_powers(l):
    """Cartesian powers (lx, ly, lz) of a shell, in the usual order xx, xy, xz, yy, yz, zz, ..."""
    return [(lx, ly, l - lx - ly) for lx in range(l, -1, -1) for ly in range(l - lx, -1, -1)]


@dataclass
class Basis:
    """Primitive Cartesian Gaussians grouped into primitive shells (every array has one entry
    per primitive Cartesian function).

    Attributes
    ----------
    xyz : ndarray, shape (nprim, 3)
        Centres in bohr.
    ga : ndarray
        Exponents.
    lmn : ndarray, shape (nprim, 3)
        Cartesian powers.
    group : ndarray
        Primitive shell index, 1-based, contiguous.  A primitive shell collects all functions
        that share centre, exponent and l+m+n.
    contr : ndarray
        Index (1-based) of the contracted function the primitive belongs to.
    coeff : ndarray
        Contraction coefficient times normalisation.
    nao : int
        Number of contracted functions.
    """

    xyz: np.ndarray
    ga: np.ndarray
    lmn: np.ndarray
    group: np.ndarray
    contr: np.ndarray
    coeff: np.ndarray
    nao: int

    def __post_init__(self):
        self.xyz = np.asfortranarray(self.xyz, dtype=np.float64)
        self.ga = np.asfortranarray(self.ga, dtype=np.float64)
        self.lmn = np.asfortranarray(self.lmn, dtype=np.int32)
        self.group = np.asfortranarray(self.group, dtype=np.int32)
        self.contr = np.asfortranarray(self.contr, dtype=np.int32)
        self.coeff = np.asfortranarray(self.coeff, dtype=np.float64)
        self.nao = int(self.nao)
        n = self.ga.size
        if not (self.xyz.shape == (n, 3) and self.lmn.shape == (n, 3)
                and self.group.shape == (n,) and self.contr.shape == (n,) and self.coeff.shape == (n,)):
            raise ValueError("inconsistent array sizes in Basis")

    @property
    def nprim(self):
        return self.ga.size

    @property
    def lmax(self):
        return int(self.lmn.sum(axis=1).max())

    @classmethod
    def from_shells(cls, shells, merge=True):
        """Build a basis from a list of shells.

        Each shell is a tuple ``(centre, l, exponents, coefficients)`` with ``centre`` in bohr and
        ``coefficients`` of shape (nprim, nctr) (or (nprim,) for a single contraction).  The
        coefficients must multiply the *bare* primitives x^l y^m z^n exp(-a r^2), i.e. they
        include every normalisation factor.  An optional fifth element gives per-component
        factors of length (l+1)(l+2)/2 for packages that normalise Cartesian components
        individually.  Contracted functions are numbered shell by shell, contraction by
        contraction, Cartesian component by Cartesian component.
        """
        # members[g] collects (lmn, contracted index, coefficient) of primitive shell g.  Primitive
        # shells that appear in several contractions (general contraction written as separate
        # segmented shells, as in most Molden files) are merged, so that their Hermite and
        # Bessel factors are computed once.
        keys = {}
        members = []
        centres, exps_g, ls = [], [], []
        iao = 0
        for shell in shells:
            centre, l, exps, cs = shell[:4]
            centre = np.asarray(centre, dtype=float)
            cs = np.asarray(cs, dtype=float)
            if cs.ndim == 1:
                cs = cs[:, None]
            exps = np.asarray(exps, dtype=float)
            powers = cartesian_powers(l)
            comp = np.ones(len(powers)) if len(shell) < 5 else np.asarray(shell[4], dtype=float)
            nctr = cs.shape[1]
            for ip, a in enumerate(exps):
                key = (tuple(np.round(centre, 10)), float(a), int(l)) if merge else len(members)
                if key not in keys:
                    keys[key] = len(members)
                    members.append([])
                    centres.append(centre)
                    exps_g.append(a)
                    ls.append(l)
                g = keys[key]
                for ic in range(nctr):
                    for icart, p in enumerate(powers):
                        if cs[ip, ic] != 0.0:
                            members[g].append((p, iao + ic * len(powers) + icart + 1, cs[ip, ic] * comp[icart]))
            iao += nctr * len(powers)
        xyz, ga, lmn, group, contr, coeff = [], [], [], [], [], []
        ig = 0
        for g, mem in enumerate(members):
            if not mem:
                continue
            ig += 1
            for p, c, v in mem:
                xyz.append(centres[g])
                ga.append(exps_g[g])
                lmn.append(p)
                group.append(ig)
                contr.append(c)
                coeff.append(v)
        return cls(np.array(xyz), np.array(ga), np.array(lmn), np.array(group), np.array(contr),
                   np.array(coeff), iao)

    def eval(self, points):
        """Values of the contracted functions at ``points`` (npts, 3); returns (npts, nao)."""
        points = np.atleast_2d(points)
        d = points[:, None, :] - self.xyz[None, :, :]
        prim = (d[:, :, 0] ** self.lmn[:, 0] * d[:, :, 1] ** self.lmn[:, 1] * d[:, :, 2] ** self.lmn[:, 2]
                * np.exp(-self.ga * (d ** 2).sum(axis=2))) * self.coeff
        out = np.zeros((points.shape[0], self.nao))
        np.add.at(out.T, self.contr - 1, prim.T)
        return out
