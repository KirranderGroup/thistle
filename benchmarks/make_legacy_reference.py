# Regenerates tests/data/legacy_*.npz with the pre-existing code (needs the old Fortran executable and cumulant_scattering).
"""Produce reference curves with the CURRENT code path (cumulant_scattering + Fortran one-loop kernel)."""
import sys, os, subprocess, numpy as np
sys.path.insert(0, '/home/claude/legacy/py')
from pyscf import gto, scf, cc
import pyscf.tools.molden
import cumulant_scattering.rdm.mo2ao as mo2ao
from cumulant_scattering.interfaces.File_just_pyscf import prepare_inputs_from_pyscf_molden, PySCFScatteringConfig

EXE = '/home/claude/legacy/Main_ref.exe'

def write_molden(mf, fname='X.mld'):
    pyscf.tools.molden.from_scf(mf, fname, ignore_h=True)
    with open(fname) as f:
        lines = [l for l in f if l.strip().upper() not in ('[6D]', '[10F]', '[15G]', '[5D]', '[7F]', '[9G]')]
    out = []
    for l in lines:
        if l.strip().startswith('[MO]'):
            while out and not out[-1].strip():
                out.pop()
        if l.lstrip().startswith('Ene='):
            l = f" Ene= {float(l.split('=')[1]):.10f}\n"
        if l.lstrip().startswith('Occup='):
            l = f" Occup= {float(l.split('=')[1]):.10f}\n"
        out.append(l)
    with open(fname, 'w') as f:
        f.writelines(out)

def legacy_total(mf, mol, Gamma_mo, tag, qmax=8.0, npoints=41):
    d = f'work_{tag}'; os.makedirs(d, exist_ok=True); cwd = os.getcwd(); os.chdir(d)
    try:
        write_molden(mf, 'X.mld')
        mo2ao.create_Zcotr(mf, mol, Gamma_mo)
        cfg = PySCFScatteringConfig(outfile='out.dat', qmin=1e-10, Type=1, qmax=qmax, npoints=npoints,
                                    cutoffcentre=0.01, cutoffz=1e-20, cutoffmd=1e-20)
        prepare_inputs_from_pyscf_molden(cfg)
        env = dict(os.environ, SCREEN_TOL='0', OMP_NUM_THREADS='2')
        r = subprocess.run([EXE], env=env, capture_output=True, text=True)
        open('log.txt', 'w').write(r.stdout + r.stderr)
        return np.loadtxt('out.dat')
    finally:
        os.chdir(cwd)

cases = {
 'h2o_sto3g': ('O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692', 'sto-3g'),
 'h2o_631g':  ('O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692', '6-31g'),
 'h2o_ccpvdz':('O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692', 'cc-pvdz'),
 'co_631gs':  ('C 0.1 -0.2 0; O 0.3 0.4 1.1', '6-31g*'),
}
only = sys.argv[1:] or list(cases)
for tag in only:
    atom, basis = cases[tag]
    mol = gto.M(atom=atom, basis=basis, cart=True, verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-12)
    mycc = cc.CCSD(mf).run(conv_tol=1e-10)
    g1 = mycc.make_rdm1(); g2 = mycc.make_rdm2()
    N = mol.nelectron
    # HF 2RDM in MO basis
    occ = np.diag(mf.mo_occ)
    ghf = np.einsum('pq,rs->pqrs', occ, occ) - 0.5*np.einsum('ps,qr->pqrs', occ, occ)
    res = {}
    for name, G in (('hf', ghf), ('ccsd', g2)):
        out = legacy_total(mf, mol, G, f'{tag}_{name}')
        res[name] = out
        print(tag, name, 'I(q0)=', out[0,1], 'N^2=', N*N, ' I(qmax)=', out[-1,1], flush=True)
    np.savez(f'ref_{tag}.npz', q=res['hf'][:,0], hf=res['hf'][:,1], ccsd=res['ccsd'][:,1],
             atom=atom, basis=basis)
