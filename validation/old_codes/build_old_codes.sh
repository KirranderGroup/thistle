#!/bin/bash
# Build the earlier group codes with gfortran (OpenBLAS) in scratch copies; the originals are not touched.
#   build_old_codes.sh BIGMEMORY_SRC EIRIK_J2_SRC OUTDIR
# BIGMEMORY_SRC = .../BigMemoryScattering, EIRIK_J2_SRC = .../Scattering_eirik/Eirik_j2
# Produces OUTDIR/bigmem/Main_opt.exe, OUTDIR/eirik/Main.exe (as-is kernels) and
# OUTDIR/eirik_repaired/Main.exe (eirik_repairs.patch on top).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BM=$(cd "$1" && pwd); EI=$(cd "$2" && pwd); OUT=$3
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
F="gfortran -O3 -march=native -fPIC -fopenmp -ffree-line-length-none -fallow-argument-mismatch -fno-range-check -c"

# BigMemoryScattering: its own recipe (script_compilation_fortran.sh) with gfortran.
# Patch: the kernel gate `if (2.eqv..False.)` is .true. under ifort (2 has low bit 0, i.e. .false.)
# and is written out as .True.; integer kinds in one_rdm.f90.
rm -rf "$OUT/bigmem"; mkdir -p "$OUT/bigmem"; cp "$BM"/*.f90 "$OUT/bigmem/"; cd "$OUT/bigmem"
patch -s -p1 < "$HERE/bigmemory_gfortran.patch"
for f in Reader calculate_form_factors unique_module twordms one_rdm p0_cases MD Bessels_opt Total_scattering_one_loop_opt linspace; do $F $f.f90; done
gfortran -O3 -march=native -fopenmp -ffree-line-length-none -fallow-argument-mismatch Main.f90 *.o -o Main_opt.exe -lopenblas

# Eirik_j2: its recipe with gfortran. Patch: integer kinds and logical comparisons; Main reads the
# 2RDM file with the number of lines counted (was 4549 fixed); Type 7/8 read the 2RDM file.
build_eirik () {
  rm -rf "$1"; mkdir -p "$1"; cp "$EI"/*.f90 "$1/"; cd "$1"
  patch -s -p1 < "$HERE/eirik_gfortran.patch"
  [ -n "$2" ] && patch -s -p1 < "$2"
  for f in calculate_form_factors unique_module twordms one_rdm Integrals Variables main_calculation linspace; do $F $f.f90; done
  gfortran -O3 -march=native -fopenmp -ffree-line-length-none Main.f90 *.o -o Main.exe -lopenblas
}
build_eirik "$OUT/eirik" ""
build_eirik "$OUT/eirik_repaired" "$HERE/eirik_repairs.patch"
echo "built: $OUT/bigmem/Main_opt.exe $OUT/eirik/Main.exe $OUT/eirik_repaired/Main.exe"
