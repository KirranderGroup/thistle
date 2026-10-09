# grid.f90 / qgrid3d.f90 (provenance)

The polar detector grid of the earlier code, kept here unchanged.  THISTLE reproduces it in
`pyxscat.grids.qcoord3` (point for point; `tests/test_grids.py` compares with the `qgrid.dat` this
program writes for e0 = 10 keV, thsmax = 60 degrees, n1 = 49, n2 = 73, stored in
`tests/data/grid/qgrid_10keV_60deg_49_73.dat`).

Note on the arguments: the grid has `n1` rings and **n2 - 1** points per ring (the azimuth loop runs
over n2 - 1 values), plus q = 0; npts = n1 (n2 - 1) + 1.  In THISTLE: `rings = n1`, `azimuths = n2 - 1`.

    thistle run lif.out --orientation 0 90 0 --grid polar --photon-kev 10 --rings 49 --azimuths 72
    pyxscat --molden h2o.molden --kind hf_total --polar-grid 10 60 49 73 -o pattern.dat
