# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
# Parts derived from PySCF (pyscf/cc/ccsd_rdm.py), Copyright 2014-2021 The PySCF Developers,
# Apache License 2.0; see NOTICE.
"""The CCSD two-particle density as a sum of products, without the four-index array.

PySCF's ``make_rdm2`` builds the full nmo^4 (then nao^4) tensor, several copies of it: a
molecule of 114 basis functions already needs 6 GB.  The density has structure, though: every
block of the cumulant carries occupied indices except the virtual-virtual-virtual-virtual
one, and that one is exactly a sum over occupied pairs of products,

    Gamma_abcd (vvvv) = sum_ij lambda2_ij,ac tau_ij,bd  (+ permutations),  tau = t2 + t1 t1,

so the whole density can be written as :class:`pyxscat.core.FactoredDM2`,

    Gamma[mu,nu,la,si] = sum_k Ac_k[mu,nu] Bc_k[la,si] + sum_k Ax_k[mu,la] Bx_k[nu,si],

with K = 7 nocc^2 + 4 terms: memory 2 K nao^2 instead of nao^4 (benzene/cc-pVDZ: 0.6 GB
against 1.35 GB for the array alone; 500 functions, 40 occupied: 45 GB against 360 GB).  The
scattering kernels contract such a density directly (kernel mode 3).

The terms reproduce PySCF's ``ccsd_rdm._make_rdm2`` (with the one-particle terms, in
chemists' notation, trace N(N-1)) exactly when contracted with any kernel that is symmetric
under the 8 index permutations, which all scattering kernels are; the intermediates of
passes 1 and 2 of PySCF's ``_gamma2_outcore`` (BSD licence, Q. Sun et al.) are built here in
memory (they are at most nocc^2 nvir^2), its pass 3 (the nvir^4 and nocc nvir^3 blocks) is
replaced by the factorisation.  Checked against ``make_rdm2`` to 1e-13 (tests/test_ccsd_factored.py).
"""
import numpy as np

from .core import FactoredDM2


def _intermediates(t1, t2, l1, l2):
    """doooo, dooov, dovov, doovv, dovvo (as PySCF defines them) and the pass-1 intermediates
    mvOOv, mvoOV, mvv.  A transcription of passes 1-2 of pyscf.cc.ccsd_rdm._gamma2_outcore."""
    einsum = np.einsum
    nocc, nvir = t1.shape
    pvOOv = einsum("ikca,jkcb->aijb", l2, t2, optimize=True)
    moo = einsum("dljd->jl", pvOOv) * 2
    mvv = einsum("blld->db", pvOOv) * 2
    gooov = einsum("kc,cija->jkia", t1, pvOOv, optimize=True)
    mvOOv = pvOOv
    pvoOV = -einsum("ikca,jkbc->aijb", l2, t2, optimize=True)
    theta = t2 * 2 - t2.transpose(0, 1, 3, 2)
    pvoOV += einsum("ikac,jkbc->aijb", l2, theta, optimize=True)
    moo += einsum("dljd->jl", pvoOV)
    mvv += einsum("blld->db", pvoOV)
    gooov -= einsum("jc,cika->jkia", t1, pvoOV, optimize=True)
    mvoOV = pvoOV
    mia = einsum("kc,ikac->ia", l1, t2) * 2 - einsum("kc,ikca->ia", l1, t2)
    mab = einsum("kc,kb->cb", l1, t1)
    mij = einsum("kc,jc->jk", l1, t1) + moo * .5
    tau = einsum("ia,jb->ijab", t1, t1) + t2
    goooo = einsum("ijab,klab->ijkl", tau, l2, optimize=True) * .5
    doooo = goooo.transpose(0, 2, 1, 3) * 2 - goooo.transpose(0, 3, 1, 2)
    gooov += einsum("ji,ka->jkia", -.5 * moo, t1)
    gooov += einsum("la,jkil->jkia", 2 * t1, goooo, optimize=True)
    gooov -= einsum("ib,jkba->jkia", l1, tau, optimize=True)
    gooov -= einsum("jkba,ib->jkia", l2, t1, optimize=True)
    dooov = gooov.transpose(0, 2, 1, 3) * 2 - gooov.transpose(1, 2, 0, 3)
    goovv = einsum("ia,jb->ijab", mia, t1)
    tmpoovv = einsum("ijkl,klab->ijab", goooo, tau, optimize=True)
    tmpoovv -= einsum("jk,ikab->ijab", mij, tau, optimize=True)
    tmpoovv -= einsum("cb,ijac->ijab", mab, t2, optimize=True)
    tmpoovv -= einsum("bd,ijad->ijab", mvv * .5, tau, optimize=True)
    tmpoovv += .5 * tau
    tmpoovv += .5 * l2
    goovv += tmpoovv
    gOvvO = einsum("kiac,jc,kb->iabj", l2, t1, t1, optimize=True)
    gOvvO += einsum("aijb->iabj", pvOOv)
    govVO = einsum("ia,jb->iabj", l1, t1)
    govVO -= einsum("ikac,jc,kb->iabj", l2, t1, t1, optimize=True)
    govVO += einsum("aijb->iabj", pvoOV)
    dovvo = 2 * govVO + gOvvO
    doovv = (-2 * gOvvO - govVO).transpose(3, 0, 1, 2)
    tau2 = tau - t2 * .5
    goovv += einsum("dlib,jlda->ijab", pvOOv, tau2, optimize=True)
    goovv -= einsum("dlia,jldb->ijab", pvoOV, tau2, optimize=True)
    goovv += einsum("dlia,jlbd->ijab", pvoOV + pvOOv * .5, t2, optimize=True)
    dovov = goovv.transpose(0, 2, 1, 3) * 2 - goovv.transpose(1, 2, 0, 3)
    return dict(doooo=doooo, dooov=dooov, dovov=dovov, doovv=doovv, dovvo=dovvo,
                mvOOv=mvOOv, mvoOV=mvoOV, mvv=mvv, tau=tau)


def factored_rdm2(mycc, t1=None, t2=None, l1=None, l2=None, dm1=None):
    """CCSD spin-summed two-particle density as a :class:`~pyxscat.core.FactoredDM2` in the
    molecular-orbital basis of ``mycc`` (all orbitals, frozen ones included; chemists'
    notation, trace N(N-1), the same contraction values as ``mycc.make_rdm2()``).

    ``dm1``: the MO one-particle density (``mycc.make_rdm1()`` by default).  The amplitudes
    default to those stored in ``mycc`` (``solve_lambda`` is run if needed).
    """
    t1 = mycc.t1 if t1 is None else t1
    t2 = mycc.t2 if t2 is None else t2
    if l1 is None or l2 is None:
        if mycc.l1 is None or mycc.l2 is None:
            mycc.solve_lambda()
        l1, l2 = mycc.l1, mycc.l2
    if dm1 is None:
        dm1 = mycc.make_rdm1()
    nmo = mycc.mo_occ.size
    mask = mycc.get_frozen_mask()
    act = np.flatnonzero(mask)
    nocc_act = t1.shape[0]
    occ_act, vir_act = act[:nocc_act], act[nocc_act:]
    occ_all = np.flatnonzero(mycc.mo_occ > 0)
    return _build(nmo, occ_all, occ_act, vir_act, np.asarray(dm1, float), t1, t2, l1, l2)


def _build(nmo, occ_all, occ, vir, dm1, t1, t2, l1, l2):
    nocc, nvir = t1.shape
    X = _intermediates(t1, t2, l1, l2)
    tau = X["tau"]
    zero = lambda: np.zeros((nmo, nmo))

    def ov(m):                      # an (nocc_act, nvir_act) block into an nmo x nmo matrix
        out = zero(); out[np.ix_(occ, vir)] = m; return out

    def vv(m):
        out = zero(); out[np.ix_(vir, vir)] = m; return out

    def oo(m):
        out = zero(); out[np.ix_(occ, occ)] = m; return out

    def e_pair(i, j):               # e_i e_j^T
        out = zero(); out[occ[i], occ[j]] = 1.0; return out

    def e_vec(i, v):                # e_i v^T with v over the active virtuals
        out = zero(); out[occ[i], vir] = v; return out

    # ---- one-particle part: 4 P (D - P) Coulomb, -2 P (D - P) exchange, P = occupied projector
    P = zero(); P[occ_all, occ_all] = 1.0
    DP = dm1 - P
    ac, bc, ax, bx = [4.0 * P], [DP], [-2.0 * P], [DP]

    # ---- cumulant blocks with two occupied indices (weight 4 each: the four images of each
    #      block in the full tensor are related by symmetries of the kernel)
    M, N, L, T = X["mvoOV"], X["mvOOv"], l2, t2          # M[a,k,i,c], N[a,k,i,b], L[j,i,b,c], T[j,i,b,c]
    for i in range(nocc):
        for j in range(nocc):
            # Coulomb with A = e_i e_j^T: doovv (vv), doooo (oo), dooov (ov)
            ac.append(4.0 * e_pair(i, j))
            bc.append(vv(X["doovv"][i, j]) + oo(X["doooo"][i, j]) + ov(X["dooov"][i, j]))
            # exchange with A = e_i e_j^T: dovov[i,a,j,b] on (a,b); dovvo[i,a,b,j] on (a,b)
            ax.append(4.0 * e_pair(i, j))
            bx.append(vv(X["dovov"][i, :, j, :] + X["dovvo"][i, :, :, j]))
            # exchange, vvvv block: A = lambda2_ij on (a, y), B = tau_ij - 1/2 tau_ij^T on (x, z)
            ax.append(4.0 * vv(l2[i, j]))
            bx.append(vv(tau[i, j] - 0.5 * tau[i, j].T))
            # ---- ovvv block: dovvv[i,a,b,c] = 2 g[b,i,c,a] - g[b,i,a,c] with g = gvovv ----
            # (1) -sum_d gvvvv[a',d,b',c'] t1[i',d]: A[i',a'] = (t1 lambda2_jk^T)[i',a'], B = tau_jk
            ax.append(ov(t1 @ l2[i, j].T))
            bx.append(vv(4.0 * tau[i, j] - 8.0 * tau[i, j].T))
            # (2),(3),(7): A = e_i (t1_k)^T  [here k -> j]
            ac.append(e_vec(i, t1[j]))
            bc.append(vv(-4.0 * M[:, j, i, :] - 8.0 * N[:, j, i, :]))
            ax.append(e_vec(i, t1[j]))
            bx.append(vv(8.0 * M[:, j, i, :].T + 4.0 * N[:, j, i, :].T + 8.0 * L[j, i].T - 4.0 * L[j, i]))
            # (4): A = e_i (lambda1_j)^T, B = 8 t2_ji^T - 4 t2_ji
            ax.append(e_vec(i, l1[j]))
            bx.append(vv(8.0 * T[j, i].T - 4.0 * T[j, i]))
    # (5), (6): single terms with A = t1 (occupied x virtual)
    ac.append(ov(t1))
    bc.append(vv(8.0 * l1.T @ t1 + 4.0 * X["mvv"].T))
    ax.append(ov(t1))
    bx.append(vv(-4.0 * t1.T @ l1 - 2.0 * X["mvv"]))
    return FactoredDM2(np.array(ac), np.array(bc), np.array(ax), np.array(bx))


def contract_dense(g, k4):
    """sum Gamma K for a four-index or factored density (tests)."""
    if isinstance(g, FactoredDM2):
        return float(np.einsum("abcd,abcd->", g.dense(), k4))
    return float(np.einsum("abcd,abcd->", g, k4))


# ----------------------------------------------------------------------------------- CISD
def _cisd_intermediates(c0, c1, c2):
    """doooo, dooov, dovov, doovv, dovvo of a CISD vector (as PySCF's cisd._gamma2_outcore)."""
    einsum = np.einsum
    dovov = 2 * c0 * c2.transpose(0, 2, 1, 3) - c0 * c2.transpose(1, 2, 0, 3)
    g = einsum("ijab,klab->ijkl", c2, c2, optimize=True)
    doooo = g.transpose(0, 2, 1, 3) - g.transpose(1, 2, 0, 3) * .5
    g = -einsum("ia,klac->klic", c1 * 2, c2, optimize=True)
    dooov = g.transpose(0, 2, 1, 3) * 2 - g.transpose(1, 2, 0, 3)
    theta = c2 * 2 - c2.transpose(1, 0, 2, 3)
    doovv = einsum("ia,kc->ikca", c1, -c1)
    doovv -= einsum("kjcb,kica->jiab", c2, theta, optimize=True)
    doovv -= einsum("ikcb,jkca->ijab", c2, theta, optimize=True)
    dovvo = einsum("ikac,jkbc->iabj", theta, theta, optimize=True)
    dovvo += einsum("ia,kc->iack", c1, c1) * 2
    return dict(doooo=doooo, dooov=dooov, dovov=dovov, doovv=doovv, dovvo=dovvo)


def factored_rdm2_cisd(myci, civec=None, dm1=None):
    """CISD spin-summed two-particle density as a :class:`~pyxscat.core.FactoredDM2` in the MO
    basis of ``myci`` (all orbitals, frozen ones included; same contraction values as
    ``myci.make_rdm2()``).  The structure is that of the CCSD density with lambda2 = tau = c2:

        Gamma_abcd (vvvv) = sum_ij c2_ij,ac (c2_ij,bd - 1/2 c2_ij,db),
        Gamma_iabc (ovvv) from 2 sum_j c1_ja c2_ji,bc,

    K = 4 nocc^2 + 2 terms (with nocc active occupied orbitals)."""
    civec = myci.ci if civec is None else civec
    nmo = myci.mo_occ.size
    mask = myci.get_frozen_mask()
    act = np.flatnonzero(mask)
    c0, c1, c2 = myci.cisdvec_to_amplitudes(civec, mask.sum(), myci.nocc, copy=False)
    if dm1 is None:
        dm1 = myci.make_rdm1(civec)
    nocc = c1.shape[0]
    occ, vir = act[:nocc], act[nocc:]
    occ_all = np.flatnonzero(myci.mo_occ > 0)
    X = _cisd_intermediates(c0, c1, c2)
    dm1 = np.asarray(dm1, float)
    zero = lambda: np.zeros((nmo, nmo))

    def ov(m):
        out = zero(); out[np.ix_(occ, vir)] = m; return out

    def vv(m):
        out = zero(); out[np.ix_(vir, vir)] = m; return out

    def oo(m):
        out = zero(); out[np.ix_(occ, occ)] = m; return out

    def e_pair(i, j):
        out = zero(); out[occ[i], occ[j]] = 1.0; return out

    def e_vec(i, v):
        out = zero(); out[occ[i], vir] = v; return out

    P = zero(); P[occ_all, occ_all] = 1.0
    DP = dm1 - P
    ac, bc, ax, bx = [4.0 * P], [DP], [-2.0 * P], [DP]
    t1 = 2.0 * c1                        # the ovvv block has the CCSD form with t1 -> 2 c1, l2 -> c2
    for i in range(nocc):
        for j in range(nocc):
            ac.append(4.0 * e_pair(i, j))
            bc.append(vv(X["doovv"][i, j]) + oo(X["doooo"][i, j]) + ov(X["dooov"][i, j]))
            ax.append(4.0 * e_pair(i, j))
            bx.append(vv(X["dovov"][i, :, j, :] + X["dovvo"][i, :, :, j]))
            ax.append(4.0 * vv(c2[i, j]))
            bx.append(vv(c2[i, j] - 0.5 * c2[i, j].T))
            ax.append(e_vec(i, t1[j]))
            bx.append(vv(8.0 * c2[j, i].T - 4.0 * c2[j, i]))
    return FactoredDM2(np.array(ac), np.array(bc), np.array(ax), np.array(bx))
