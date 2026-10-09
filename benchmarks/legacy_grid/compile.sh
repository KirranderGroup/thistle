#/bin/bash 

rm -rf *.mod *.o

ifx -O3 -g -traceback -qopenmp -xCORE-AVX512 -heap-arrays 0 -qmkl -c qgrid3d.f90

ifx -O3 -g -traceback -qopenmp -xCORE-AVX512 -heap-arrays 0 -qmkl grid.f90 *.o -o grid.exe
