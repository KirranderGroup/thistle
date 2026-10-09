"""python check.py file.molden REFERENCE_SCF_ENERGY"""
import sys

from pyxscat.interfaces.molden import read_molden, scf_energy, check_orbitals

path, ref = sys.argv[1], float(sys.argv[2])
data = read_molden(path)
e, resid = scf_energy(path)
print(f"basis functions (Cartesian): {data.basis.nao}   MOs: {data.mo_coeff.shape[1]}   "
      f"Cartesian normalisation: {data.cart_norm}")
print(f"MO orthonormality error: {check_orbitals(data):.1e}   basis-fit residual: {resid:.1e}")
print(f"SCF energy from the Molden file: {e:.9f}   program: {ref:.9f}   difference: {e - ref:.2e} Eh")
