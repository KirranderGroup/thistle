#!/bin/bash
# Build thistle-dev's Fortran engine (Main.exe) with gfortran + OpenBLAS, no Intel compiler or MKL,
# linking the HDF5 C libraries bundled with h5py.  Apply thistle_dev_gfortran.patch first
# (casts and one `use ifport` only; the five goldens then pass to 4e-14).
#   cd THISTLE_DEV/fortran/src && patch -p2 < .../thistle_dev_gfortran.patch && bash .../build_thistle_dev_gfortran.sh THISTLE_DEV
set -e
DEV="$(cd "$1" && pwd)"; cd "$DEV/fortran/src"
L="$(python3 -c 'import h5py,os;print(os.path.join(os.path.dirname(os.path.dirname(h5py.__file__)),"h5py.libs"))')"
FLAGS="-O3 -march=native -fPIC -fopenmp -ffree-line-length-none -fno-range-check"
rm -f ./*.o ./*.mod
for f in Logging calculate_form_factors unique_module MD Bessels_j0 p0_cases twordms one_rdm Integrals Variables \
         Reader ReaderHDF5 Total_scattering_j0_groups_fast main_calculation linspace; do
    gfortran $FLAGS -c "$f.f90"
done
gfortran $FLAGS Main.f90 ./*.o -lopenblas "$L"/libhdf5_hl-*.so* "$L"/libhdf5-*.so* -Wl,-rpath,"$L" -o "$DEV/fortran/Main.exe"
echo "built $DEV/fortran/Main.exe"
