# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""PySCF adapter: turns PySCF objects into a :class:`pyxscat.Wavefunction`.

The basis is taken from the ``Mole`` object in memory (no Molden file).  Spherical (pure) basis
sets are handled by expressing the density matrices in the underlying Cartesian functions.
"""
import numpy as np

from ..basis import Basis, cartesian_powers
from ..core import FACTORED_ABOVE, Wavefunction, transform_dm2 as _transform4, use_factored

# libcint folds the angular normalisation of s and p functions into the Cartesian functions
_ANG = {0: 0.282094791773878143, 1: 0.488602511902919921}


def basis_from_mole(mol):
    """Cartesian :class:`Basis` of a PySCF ``Mole`` (functions ordered as in ``mol.cart = True``)."""
    shells = []
    for ib in range(mol.nbas):
        l = int(mol.bas_angular(ib))
        cs = np.asarray(mol._libcint_ctr_coeff(ib)) * _ANG.get(l, 1.0)
        shells.append((mol.bas_coord(ib), l, mol.bas_exp(ib), cs))
    return Basis.from_shells(shells)


def _to_cartesian(mol, dm1, dm2):
    """Density matrices in mol's AO basis -> Cartesian functions of basis_from_mole."""
    if mol.cart:
        return dm1, dm2
    c = mol.cart2sph_coeff(normalized="sp")          # chi_sph = chi_cart @ c
    dm1 = c @ dm1 @ c.T
    if dm2 is not None:
        dm2 = _transform4(dm2, c)                      # four-index array or FactoredDM2
    return dm1, dm2


def from_density(mol, dm1, dm2=None, mo_coeff=None, label=""):
    """Wavefunction from spin-summed density matrices.

    If ``mo_coeff`` is given, ``dm1`` and ``dm2`` are in that molecular-orbital basis (as returned
    by PySCF's ``make_rdm1`` / ``make_rdm2`` of post-HF methods); otherwise they are in the AO
    basis of ``mol``.  ``dm2`` follows PySCF's convention (chemists' notation, trace N(N-1)).
    """
    from ..core import FactoredDM2
    dm1 = np.asarray(dm1, dtype=float)
    if dm2 is not None and not isinstance(dm2, FactoredDM2):
        dm2 = np.asarray(dm2, dtype=float)
    if mo_coeff is not None:
        c = np.asarray(mo_coeff)
        dm1 = c @ dm1 @ c.T
        if dm2 is not None:
            dm2 = _transform4(dm2, c)
    dm1, dm2 = _to_cartesian(mol, dm1, dm2)
    return Wavefunction(basis_from_mole(mol), dm1, dm2, nelec=float(mol.nelectron), label=label,
                       atom_numbers=mol.atom_charges().astype(int), coords=mol.atom_coords())


def from_scf(mf, with_dm2=False):
    """Closed-shell (RHF/RKS) wavefunction.  The mean-field two-electron term does not need the
    four-index array (see ``pyxscat.mean_field``); set ``with_dm2=True`` to build it anyway."""
    mol = mf.mol
    dm1 = mf.make_rdm1()
    if dm1.ndim != 2:
        raise NotImplementedError("only closed-shell mean-field objects are supported so far")
    dm2 = None
    if with_dm2:
        dm2 = np.einsum("pq,rs->pqrs", dm1, dm1) - 0.5 * np.einsum("ps,rq->pqrs", dm1, dm1)
    return from_density(mol, dm1, dm2, label="SCF")




def from_cc(mycc, factored="auto"):
    """CCSD wavefunction (unrelaxed lambda densities from PySCF).

    factored
        ``True``: the two-particle density as a sum of 7 nocc^2 + 4 products of one-particle
        matrices (:mod:`pyxscat.ccsd`; memory 2 K nao^2 instead of nao^4, the kernels contract
        it directly, same numbers to 1e-13).  ``False``: PySCF's four-index array.  ``"auto"``
        (default): factored beyond ``FACTORED_ABOVE`` basis functions.
    """
    factored = use_factored(factored, _nao(mycc.mol))
    dm1 = mycc.make_rdm1()
    if factored:
        from .. import ccsd
        dm2 = ccsd.factored_rdm2(mycc, dm1=dm1)
    else:
        dm2 = mycc.make_rdm2()
    return from_density(mycc.mol, dm1, dm2, mo_coeff=mycc.mo_coeff, label="CCSD")


def from_ci(myci, factored="auto"):
    """CISD wavefunction.  ``factored``: as in :func:`from_cc` (K = 4 nocc^2 + 2 products of
    one-particle matrices instead of the n^4 array; :func:`pyxscat.ccsd.factored_rdm2_cisd`)."""
    dm1 = myci.make_rdm1()
    if use_factored(factored, _nao(myci.mol)):
        from .. import ccsd
        dm2 = ccsd.factored_rdm2_cisd(myci, dm1=dm1)
    else:
        dm2 = myci.make_rdm2()
    return from_density(myci.mol, dm1, dm2, mo_coeff=myci.mo_coeff, label="CISD")


def from_fci(mol, mo_coeff, civec, norb, nelec, factored="auto"):
    """Full CI wavefunction in the orbitals ``mo_coeff[:, :norb]`` (no frozen core).
    ``factored``: as in :func:`from_casscf`."""
    from pyscf import fci
    dm1, dm2 = fci.direct_spin1.make_rdm12(civec, norb, nelec)
    if use_factored(factored, _nao(mol)):
        from .. import ci as _ci
        dm1, dm2 = _ci.embed_core(dm1, dm2, 0, factored=True)
    return from_density(mol, dm1, dm2, mo_coeff=np.asarray(mo_coeff)[:, :norb], label="FCI")


def from_fci_solver(solver, mo_coeff, state=None, factored="auto"):
    """A standalone PySCF FCI or selected-CI solver after ``kernel()`` (``fci.FCI(mf)``,
    ``fci.selected_ci.SCI``).  PySCF does not keep the orbitals in the solver (only the
    integrals), so they are given: ``mo_coeff`` is the array or the SCF object.  ``state``
    selects a root (0-based) when ``nroots > 1``; the densities come from the solver's own
    ``make_rdm12``."""
    mo = getattr(mo_coeff, "mo_coeff", mo_coeff)
    if mo is None:
        raise ValueError("a PySCF FCI solver does not keep its orbitals: "
                         "load it with mo_coeff=mf (or the orbital array)")
    ci = solver.ci
    if isinstance(ci, (list, tuple)):
        if state is None:
            raise ValueError("several roots: choose one with `state`")
        ci = ci[state]
    elif state not in (None, 0):
        raise ValueError(f"state {state} requested but the solver holds a single root")
    norb, nelec = solver.norb, solver.nelec
    dm1, dm2 = solver.make_rdm12(ci, norb, nelec)
    if use_factored(factored, _nao(solver.mol)):
        from .. import ci as _ci
        dm1, dm2 = _ci.embed_core(dm1, dm2, 0, factored=True)
    return from_density(solver.mol, dm1, dm2, mo_coeff=np.asarray(mo)[:, :norb], label="FCI")


def from_casscf(mc, state=None, factored="auto"):
    """CASSCF/CASCI wavefunction; for state-averaged calculations pick the root with ``state``.

    factored
        ``True``: the 2RDM as products of one-particle matrices (core-core and core-active
        blocks from the core and active 1RDMs, the active block diagonalised as a (pq) x (rs)
        matrix, at most n_act (n_act + 1) / 2 + 3 terms; :func:`pyxscat.ci.embed_core`), never
        the n_AO^4 array.  ``False``: the four-index array over the core and active orbitals,
        transformed to AOs.
        ``"auto"`` (default): factored beyond ``FACTORED_ABOVE`` basis functions.
    """
    from .. import ci as _ci
    d1a, d2a = _active_rdm12(mc, state)
    f = use_factored(factored, _nao(mc.mol))
    dm1, dm2 = _ci.embed_core(d1a, d2a, mc.ncore, factored=f)
    c = np.asarray(mc.mo_coeff)[:, :mc.ncore + mc.ncas]       # core and active orbitals only
    return from_density(mc.mol, dm1, dm2, mo_coeff=c, label="CASSCF")


def _active_rdm12(mc, state=None):
    """Active-space 1- and 2-RDMs (PySCF convention) of one root; state-averaged solvers return
    every root through ``states_make_rdm12``."""
    ci = mc.ci
    if isinstance(ci, (list, tuple)):
        if state is None:
            raise ValueError("state-averaged calculation: choose the root with `state`")
        solver = mc.fcisolver
        if hasattr(solver, "states_make_rdm12"):
            d1, d2 = solver.states_make_rdm12(ci, mc.ncas, mc.nelecas)
            return np.asarray(d1[state]), np.asarray(d2[state])
        ci = ci[state]
    elif state not in (None, 0):
        raise ValueError(f"state {state} requested but the calculation has a single root")
    return mc.fcisolver.make_rdm12(ci, mc.ncas, mc.nelecas)


def _nao(mol):
    return mol.nao if mol.cart else mol.nao_cart()
