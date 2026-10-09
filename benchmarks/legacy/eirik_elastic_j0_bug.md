# Elastic j0 in the old code calls the j1 Bessel sum

`elastic_j0.f90`, line 744 (identical in `Scattering_eirik/`, `Scattering_eirik/Eirik_j2/` and
`BigMemoryScattering/`):

    call bessels1rr(F,LLmax,mu,H,h_saved)

`bessels1rr` is the j1 routine. It sums mu^(beta-1)/h^beta instead of (mu/h)^beta and loops
beta = 0…16 over an `h_saved` of size 13 (out of bounds; gfortran `-fcheck=bounds` stops there).
The elastic j0 result is right only at q = 1 bohr⁻¹, where mu^(beta-1) = mu^beta.

Run on H₂O RHF/6-31G (Cartesian), with the routines called directly from Fortran:

| q (bohr⁻¹) | old elastic_integration | pyxscat (validated against PySCF) | ratio |
|---|---|---|---|
| 0.3 | 161.057 | 94.492 | 1.704 |
| 0.7 | 82.557 | 74.227 | 1.112 |
| 1.0 | 55.894 | 55.894 | 1.000 |
| 2.0 | 16.803 | 16.877 | 0.996 |
| 3.0 | 8.143 | 6.575 | 1.238 |
| 5.0 | 53651.7 | 2.647 | 2.0e4 |

The total j0 kernels (`total_j0.f90`, `Total_scattering_one_loop_opt.f90`) call the j0 routine
and are not affected (the optimised one is reproduced to 1e-9 by the regression tests). The
elastic j2 routine (`elastic_j2.f90`) is correct: it equals pyxscat's j2 times sqrt(5/4π) to 1e-6
(`tests/test_j2.py::test_matches_eirik_elastic_j2`).
