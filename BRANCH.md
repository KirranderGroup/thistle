# This branch: `amc/thistle-core`

Andrés Moreno Carrascosa's THISTLE (the refactored PyXSCAT kernels, `pyxscat` + `thistle`), added
next to thistle-dev for comparison and a later merge. Nothing outside `thistle_core/` is changed.

* Install: `pip install ./thistle_core` (gfortran or ifx; no MKL or HDF5 library needed;
  `h5py` to read decks). Check: `thistle selftest`.
* thistle-dev decks run directly: `import thistle as th; th.load("tests/decks/ch4_total/deck.h5")`.
* Cross-check against this repository's goldens (agreement 4e-13 with the deck's
  cutoffcentre = 0.01): `cd thistle_core/benchmarks/thistle_dev_crosscheck && python compare_goldens.py`.
  The comparison of the two codes is in that folder's README.md.
* Licence: Apache 2.0 (`thistle_core/LICENSE`); thistle-dev's own files keep their licence.
