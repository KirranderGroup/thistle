# Reference data

`legacy_*.npz` hold the two-electron j0 scattering (Hartree-Fock and CCSD pair densities) computed with
the pre-existing production path: `cumulant_scattering` (Molden reader, `create_Zcotr`) feeding the
Fortran one-loop kernel `Total_scattering_one_loop_opt.f90`, with screening switched off
(`SCREEN_TOL=0`), `cutoffcentre=0.01`, `cutoffz=cutoffmd=1e-20`, Cartesian basis functions, q in
41 points from 1e-10 to 8 a.u.

That kernel evaluates pi in single precision (`acos(-1.0000)`), which scales every result by
(pi_single/pi_double)^3 = 1 + 8.35e-8.  The regression test applies this factor before comparing.
