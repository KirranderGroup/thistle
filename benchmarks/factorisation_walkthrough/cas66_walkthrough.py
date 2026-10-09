# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""The factorisation of a CASSCF two-particle density, step by step, on water CASSCF(6,6)/cc-pVDZ.

Every step prints what it produces; the numbers go into the note cas66_walkthrough.tex.
python cas66_walkthrough.py   (repository root or this folder)
"""
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python"))
import pyxscat as xs
import thistle as th
from pyxscat.core import symmetrize_dm2, FactoredDM2
from pyxscat.ci import factor_active, embed_core
from pyxscat.interfaces import pyscf as xp
from pyscf import gto, scf, mcscf

out = []
def log(s=""):
    print(s, flush=True); out.append(s)

# ---------------------------------------------------------------- 0. the wavefunction
mol = gto.M(atom="O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", basis="cc-pvdz", cart=True, verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-12)
mc = mcscf.CASSCF(mf, 6, 6); mc.conv_tol = 1e-11; mc.kernel()
n, nc, na = mol.nao, mc.ncore, mc.ncas
N = mol.nelectron
log(f"step 0  water CASSCF(6,6)/cc-pVDZ (Cartesian): n = {n} basis functions, {nc} core, {na} active, "
    f"{n - nc - na} virtual orbitals; N = {N} electrons; E = {mc.e_tot:.8f} Eh")

# ---------------------------------------------------------------- 1. active densities
d1a, d2a = mc.fcisolver.make_rdm12(mc.ci, na, mc.nelecas)
log(f"step 1  active 1RDM {d1a.shape}, trace {np.trace(d1a):.10f} (= 6 active electrons); "
    f"natural occupations {np.round(np.sort(np.linalg.eigvalsh(d1a))[::-1], 4).tolist()}")
log(f"        active 2RDM {d2a.shape} = {d2a.size} numbers, sum_pq Gamma_ppqq = {np.einsum('ppqq->', d2a):.10f} (= 6 x 5)")

# ---------------------------------------------------------------- 2. 8-fold average
g = symmetrize_dm2(d2a)
log(f"step 2  8-fold average: largest change {np.abs(g - d2a).max():.3f} (of max |Gamma| {np.abs(d2a).max():.3f}); "
    f"trace unchanged: {np.einsum('ppqq->', g):.10f}")

# ---------------------------------------------------------------- 3. eigen-decomposition of the (pq) x (rs) matrix
M = g.reshape(na * na, na * na)
w_all, v_all = np.linalg.eigh(M)
order = np.argsort(-np.abs(w_all))
w_all, v_all = w_all[order], v_all[:, order]
nz = np.abs(w_all) > 1e-12 * np.abs(w_all).max()
log(f"step 3  M = Gamma reshaped to {M.shape}; {nz.sum()} non-zero eigenvalues of {M.shape[0]} "
    f"(bound n_act(n_act+1)/2 = {na*(na+1)//2}); the other {(~nz).sum()} are {np.abs(w_all[~nz]).max():.0e} or less")
log("        eigenvalues w_k, largest |w| first: " + ", ".join(f"{x:.4f}" for x in w_all[nz]))
w, V = factor_active(d2a)
symm = max(np.abs(Vk - Vk.T).max() for Vk in V)
log(f"        factor_active: {len(w)} terms, each V_k a {V.shape[1]}x{V.shape[2]} matrix, symmetric to {symm:.0e}")
rec = np.einsum("k,kpq,krs->pqrs", w, V, V)
log(f"        sum_k w_k V_k(pq) V_k(rs) vs averaged Gamma: max error {np.abs(rec - g).max():.1e}")

# ---------------------------------------------------------------- 4. core embedding in the MO basis
dm1_mo, f_mo = embed_core(d1a, d2a, nc, factored=True)
nmo = nc + na
log(f"step 4  core embedding (core + active = {nmo} orbitals): {f_mo.ac.shape[0]} Coulomb terms "
    f"({len(w)} active + 3 core), {f_mo.ax.shape[0]} exchange terms, each {nmo}x{nmo}")
dc = np.zeros((nmo, nmo)); dc[:nc, :nc] = 2 * np.eye(nc)
da = np.zeros((nmo, nmo)); da[nc:, nc:] = d1a
tr = lambda a: float(np.trace(a))
parts = {
    "core Hartree  Dc x Dc": tr(dc) * tr(dc),
    "core exchange -1/2 Dc.Dc": -0.5 * float(np.sum(dc * dc)),
    "core-active Hartree  Dc x Da + Da x Dc": 2 * tr(dc) * tr(da),
    "core-active exchange -1/2 (Dc.Da + Da.Dc)": -1.0 * float(np.sum(dc * da)),
    "active block sum_k w_k tr(V_k)^2": float(np.sum(w * np.array([np.trace(v) ** 2 for v in V]))),
}
for k, v in parts.items():
    log(f"        {k:44s} contributes {v:8.4f} to N(N-1)")
log(f"        total {sum(parts.values()):.10f} = N(N-1) = {N*(N-1)}")
_, d_mo_dense = embed_core(d1a, d2a, nc, factored=False)
log(f"        rebuilt MO tensor vs dense embedding (both averaged): "
    f"{np.abs(symmetrize_dm2(f_mo.dense()) - symmetrize_dm2(d_mo_dense)).max():.1e}")

# ---------------------------------------------------------------- 5. to the AO basis
C = np.asarray(mc.mo_coeff)[:, :nmo]                       # core + active columns only
f_ao = f_mo.transform(C) if hasattr(f_mo, "transform") else None
wf_f = th.load(mc, factored=True).wf
wf_d = th.load(mc, factored=False).wf
f_ao = wf_f.dm2
log(f"step 5  A_AO = C A C^T with C the {C.shape[0]}x{C.shape[1]} core+active coefficients: "
    f"{f_ao.ac.shape[0]} + {f_ao.ax.shape[0]} terms, each {f_ao.nao}x{f_ao.nao}")
dense = wf_d.dm2
nb_f = f_ao.ac.nbytes + f_ao.bc.nbytes + f_ao.ax.nbytes + f_ao.bx.nbytes
log(f"        storage: factors {2*(f_ao.ac.shape[0]+f_ao.ax.shape[0])} matrices = "
    f"{2*(f_ao.ac.shape[0]+f_ao.ax.shape[0])*n*n:,} numbers ({nb_f/1e6:.3f} MB) vs dense {dense.size:,} ({dense.nbytes/1e6:.2f} MB)")
err_ao = np.abs(symmetrize_dm2(f_ao.dense()) - symmetrize_dm2(dense)).max()
log(f"        rebuilt AO tensor vs dense AO tensor (averaged): {err_ao:.1e}")

# ---------------------------------------------------------------- 6. one block as the kernel builds it
# shells: the kernel works shell quartet by shell quartet; take the O 2p shell (3 functions) with itself
# and the first H s shell with itself, and build Gamma(I,J,K,L) from the factors by one matrix product
lbl = mol.ao_labels()
I = [i for i, s in enumerate(lbl) if s.split()[:2] == ["0", "O"] and s.split()[2].startswith("2p")]
K = [i for i, s in enumerate(lbl) if s.split()[:3] == ["1", "H", "1s"]]
Ac = f_ao.ac[:, I][:, :, I].reshape(f_ao.ac.shape[0], -1)          # (K_C, |I||J|)
Bc = f_ao.bc[:, K][:, :, K].reshape(f_ao.bc.shape[0], -1)          # (K_C, |K||L|)
blk_c = Ac.T @ Bc                                                  # the DGEMM over k
Ax = f_ao.ax[:, I][:, :, K]; Bx = f_ao.bx[:, I][:, :, K]          # exchange pairing (mu,la)(nu,si)
blk_x = np.einsum("kac,kbd->abcd", Ax, Bx).reshape(len(I) * len(I), len(K) * len(K))
blk = blk_c + blk_x
ref = symmetrize_dm2(dense)[np.ix_(I, I, K, K)].reshape(len(I) ** 2, len(K) ** 2)
blk_s = symmetrize_dm2(f_ao.dense())[np.ix_(I, I, K, K)].reshape(len(I) ** 2, len(K) ** 2)
log(f"step 6  block Gamma(O 2p, O 2p | H 1s, H 1s), {len(I)}x{len(I)}x{len(K)}x{len(K)}: Coulomb part = A^T B "
    f"({Ac.shape[1]}x{Ac.shape[0]} times {Bc.shape[0]}x{Bc.shape[1]}), exchange part the same over {Ax.shape[0]} terms")
raw = f_ao.dense()[np.ix_(I, I, K, K)].reshape(len(I) ** 2, len(K) ** 2)
log(f"        block by the matrix products vs the same block of the rebuilt tensor: {np.abs(blk - raw).max():.1e}; "
    f"8-fold averaged, vs the dense tensor: {np.abs(blk_s - ref).max():.1e}")

# ---------------------------------------------------------------- 7. scattering
q = np.linspace(1e-8, 10.0, 41)
t0 = time.time(); a = xs.two_electron(wf_d.basis, q, dense, screen_tol=0.0); ta = time.time() - t0
t0 = time.time(); b = xs.two_electron(wf_f.basis, q, f_ao, screen_tol=0.0); tb = time.time() - t0
log(f"step 7  two-electron scattering, 41 q to 10 bohr^-1, unscreened: dense vs factored max rel {np.max(np.abs(a-b))/np.max(np.abs(a)):.1e}; "
    f"I2(0) = {b[0]:.8f} (N(N-1) = {N*(N-1)}); {ta:.2f} s dense, {tb:.2f} s factored")
m = th.load(mc, factored=True)
r = m.scattering(q=q, unit="1/bohr", parts="all")
log(f"        parts at q = 2.5 bohr^-1: total {r['total'][10]:.6f} = N {N} + Hartree {r['hartree'][10]:.6f} "
    f"+ exchange {r['exchange'][10]:.6f} + cumulant {r['cumulant'][10]:.6f}")

open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "cas66_walkthrough.txt"), "w").write("\n".join(out) + "\n")
