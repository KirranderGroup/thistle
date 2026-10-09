# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Scattering-vector grids for fixed-orientation and partially aligned calculations.

The polar detector grid is the one of ``qgrid3d.f90`` (``qcoord3``) of the earlier code, kept
point for point so that results line up with existing files:

* ``n_rings`` rings of equally spaced transverse momentum q_perp = i * qxymax / n_rings,
  i = 1 .. n_rings, with qxymax = k0 sin(theta_max) (theta_max = largest scattering angle);
* ``n_azimuth`` equally spaced azimuths per ring, phi_j = 2 pi j / n_azimuth, measured from the
  y axis: q_x = q_perp sin(phi), q_y = q_perp cos(phi);
* q_z = k0 [1 - sqrt(1 - q_perp^2 / k0^2)] on the Ewald sphere (beam along z);
* the point q = 0 first, then ring by ring, azimuth fastest:
  npts = n_rings * n_azimuth + 1  (in qcoord3's arguments: n1 = n_rings, n2 = n_azimuth + 1).

Its sign convention is q = k_in - k_out (q_z >= 0); THISTLE's square detector uses k_out - k_in.
Scattering intensities are the same at q and -q (I(-q) = I(q) for real wavefunctions), and the
P2 weight of an aligned ensemble is even, so the two conventions give identical patterns.

All q in bohr^-1 (atomic units).
"""
import numpy as np

KEV_TO_AU = 0.2681727606190511          # photon wavenumber per keV (bohr^-1), as in qgrid3d.f90
BOHR_ANGSTROM = 0.529177210903


def photon_k0(energy_kev):
    """Wavenumber of a photon of energy E (keV), bohr^-1."""
    return float(energy_kev) * KEV_TO_AU


def electron_k0(energy_kev):
    """Wavenumber of an electron of kinetic energy E (keV), relativistic, bohr^-1."""
    e = float(energy_kev) * 1e3
    lam = 12.264259 / np.sqrt(e * (1.0 + 0.978475e-6 * e))       # Angstrom
    return 2.0 * np.pi / (lam / BOHR_ANGSTROM)


def polar_grid(k0, theta_max_deg, n_rings, n_azimuth):
    """The qcoord3 polar detector grid (see module doc).  Returns qvec (npts, 3)."""
    k0 = float(k0)
    qxymax = k0 * np.sin(np.radians(theta_max_deg))
    qperp = np.arange(1, n_rings + 1) * (qxymax / n_rings)
    phi = np.arange(n_azimuth) * (2.0 * np.pi / n_azimuth)
    qp, ph = np.meshgrid(qperp, phi, indexing="ij")                    # ring-major, azimuth fastest
    qx = (np.sin(ph) * qp).ravel()
    qy = (np.cos(ph) * qp).ravel()
    qz = k0 * (1.0 - np.sqrt(1.0 - (qx ** 2 + qy ** 2) / k0 ** 2))
    return np.vstack([np.zeros((1, 3)), np.column_stack([qx, qy, qz])])


def qcoord3(e0_kev, thsmax_deg, n1, n2):
    """Exactly ``call qcoord3(e0, thsmax, n1, n2, qxyz)`` of qgrid3d.f90 (photons of e0 keV)."""
    return polar_grid(photon_k0(e0_kev), thsmax_deg, n1, n2 - 1)


def read_qgrid(path):
    """A q grid from a text file with columns qx qy qz (bohr^-1), e.g. qgrid.dat; extra columns
    (such as ``|q|``) are ignored."""
    q = np.atleast_2d(np.loadtxt(path, comments="#"))
    if q.shape[1] < 3:
        raise ValueError(f"{path}: need at least three columns qx qy qz")
    return q[:, :3]


def write_qgrid(path, qvec):
    """Write qx qy qz in the format of grid.f90 (qgrid.dat, 3f19.12)."""
    np.savetxt(path, np.asarray(qvec)[:, :3], fmt="%19.12f", delimiter="")
    return path


def unique_norms(qvec, decimals=12):
    """``|q|`` of every point, the distinct values, and the index of each point into them (so that
    rotationally averaged quantities are computed once per ring)."""
    qn = np.linalg.norm(np.asarray(qvec, float), axis=1)
    uq, inv = np.unique(np.round(qn, decimals), return_inverse=True)
    return qn, uq, inv


def aligned_ensemble(qvec, axis, i0, s2, p2=1.0):
    """Partially aligned ensemble from the isotropic signal and the j2 moment::

        I(q) = I0(|q|) + 5 <P2(cos beta)> P2(q^ . a^) S2(|q|)        (+ l = 4, ... neglected)

    ``axis``: alignment axis in the laboratory frame (e.g. the laser polarisation); ``p2``: the
    degree of alignment <P2(cos beta)> (1 = perfect); ``i0``, ``s2``: values at every point.
    The molecular axis that aligns is the z axis of the molecular frame in which S2 was computed.
    """
    qvec = np.asarray(qvec, float)
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    qn = np.linalg.norm(qvec, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.where(qn > 0, qvec @ a / np.where(qn > 0, qn, 1.0), 0.0)
    p2q = 1.5 * c * c - 0.5
    return np.asarray(i0) + 5.0 * p2 * p2q * np.asarray(s2)
