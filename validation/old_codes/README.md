# The earlier group codes against THISTLE

9 October 2026. Two earlier Fortran codes were compiled with gfortran from the group's sources,
given the same wavefunctions through the `cumulant_scattering` package, and compared with THISTLE
(`run_bigmemory.py`, `run_eirik_j2.py`; results in `bigmemory_results.txt`, `eirik_j2_results.txt`).
The original sources are not changed. The build script copies them and applies the patches here.

| code | what it computes | route |
|---|---|---|
| `BigMemoryScattering` (`Main_opt.exe`) | total two-electron term, j₀ | dense AO 2RDM (`2rdmAO`), one-loop kernel |
| `Scattering_eirik/Eirik_j2` (`Main.exe`) | total j₀ (Type 1), total j₂ (Type 7), elastic j₂ (Type 8) | sparse MO 2RDM (`2rdm.txt`), MO engine |

`cumulant_scattering` (`run_scattering`) writes the Molden file, `basis.dat`, `MOs.dat` and
`options.dat`. It also splits each 2RDM into Hartree γ⊗γ, exchange γ⊗ₓγ and cumulant Λ, so the
old codes are run on every part separately. THISTLE is run with the quartet kernel at the old codes'
`cutoffcentre = 0.01` bohr, and also at the exact coincident-centre limit.

## Results

Max |old − THISTLE| / max |THISTLE| over q = 0–10 bohr⁻¹.

**BigMemoryScattering** (total j₀, 41 q):

| case | part | cut-off 0.01 | exact |
|---|---|---|---|
| water/6-31G CCSD | Total / H / X / Λ | 2.0e-14 / 2.4e-14 / 2.9e-14 / 2.1e-12 | 1.7e-7 / 8.8e-9 / 1.5e-6 / 3.6e-7 |
| water/cc-pVDZ CASSCF(6,6) | Total / H / X / Λ | 1.3e-14 / 1.4e-14 / 2.6e-14 / 2.1e-12 | 2.3e-7 / 1.5e-7 / 2.1e-6 / 2.1e-6 |
| CO/cc-pVDZ CCSD | Total / H / X / Λ | 7.0e-14 / 2.2e-14 / 8.5e-13 / 2.3e-12 | 4.1e-7 / 3.1e-7 / 1.6e-6 / 4.0e-6 |

**Eirik_j2** (100 q). The j₂ values are compared in Eirik's convention, −√(5/4π) × THISTLE: his
code has no i² = −1 and normalises with Y₂₀ instead of P₂.

| case | observable | as compiled | repaired |
|---|---|---|---|
| water/6-31G CCSD | total j₀ / total j₂ / elastic j₂ | 2.1e-14 / 0.55 / 1.6 | — / 6.5e-14 / 7.4e-14 |
| water/cc-pVDZ CASSCF(6,6) | total j₀ / total j₂ / elastic j₂ | 1.1e-14 / 1.1 / 2.2 | — / 4.4e-14 / 1.3e-13 |
| H₄ chain/6-31G CCSD (s functions only) | total j₀ / total j₂ / elastic j₂ | 9.3e-15 / 9.5e-15 / 1.0 | — / 9.3e-15 / 1.3e-14 |

Parts, repaired build, total j₀ / total j₂: water CCSD H 2.2e-14 / 5.7e-14, X 2.2e-14 / 3.5e-13,
Λ 1.7e-12 / 1.4e-12; water CASSCF H 1.3e-14 / 6.0e-14, X 1.9e-14 / 7.3e-13, Λ 2.5e-12 / 2.0e-13.

Λ is small (|Λ| ≤ 0.4 at q = 0), so its relative error sits at ~10⁻¹². Every other number is at
round-off. The "exact" column is the 0.01-bohr cut-off error of the old codes: 10⁻⁸–10⁻⁶ at every q.

## What the comparison found

* **Single-precision π.** Both codes set `pi = acos(-1.0000)`, which is 3.14159274. Every result
  carries a factor (π_sp/π)³ = 1 + 8.35×10⁻⁸. The drivers divide it out. Without that, the
  agreement stops at 8×10⁻⁸.
* **BigMemory's kernel gate.** `Main.f90` calls the kernel under `if (2.eqv..False.)`. Under ifort
  the integer 2 is `.false.` (only the low bit is tested), so the test is true. A literal port to
  `.False.` skips the kernel and writes an unset array. The gfortran patch writes `.True.`.
* **Eirik's total j₂ (Type 7) drops every coincident-centre quartet.** The table call is commented
  out (`!CALL set_P0_j2(P0matrix, …)`, after `P0matrix = 0`), so quartets with |P − Q| < 0.01 bohr add
  nothing. That is 55–110 % of the signal for water. With s functions only, those quartets carry no
  j₂ moment, and the code as compiled is exact (H₄ row). Restoring the call gives 10⁻¹³ everywhere.
* **Eirik's elastic j₂ (Type 8) in `Integrals.f90` is a debugging state.** It sets `f = 0` after
  every two-centre integral, so only the coincident-centre quartets are summed. It writes
  `signal.dat` (0.7 M lines) and stops after the first loop, writing `j2.dat` instead of the output
  file. The correct routine is the module `elastic_j2.f90`: the earlier THISTLE test
  (`tests/data/eirik_elastic_j2.txt`) called it directly, but it is not in `Main.exe`'s build.
  Removing the `f = 0` line gives 10⁻¹³.
* **Input handling in Eirik's `Main`.** It reads three state indices (`cumulant_scattering` writes
  two). The number of 2RDM lines is fixed at 4549. Types 7 and 8 never read the 2RDM file:
  `read2rdm2 = .False.` is hard-coded, and Type 8 builds the 1RDM from CI bitstrings. Type 2
  (elastic j₀) divides by a fixed N − 1 = 23, and its j₀ routine calls the j₁ Bessel sum
  (`benchmarks/legacy/eirik_elastic_j0_bug.md`). The debug writes in the j₂ kernel need a
  `prueba_am.dat` and at least 90 q points.
* **MO-route memory.** Eirik's engine stores four arrays of (number of non-zero Γ_MO) × n_prim².
  The drivers stay at water size; CO/cc-pVDZ CCSD (3.4×10⁵ non-zero elements, 70 primitives) would need ~53 GB.
* **gfortran.** Both codes were written for ifort. The patches add integer-kind casts (`dphrec`,
  `combo`, `integer2binary_orbs`, mixed-kind `iand`) and use `.eqv.` for logical comparisons.

## Rerun

```
./build_old_codes.sh ~/BigMemoryScattering ~/Scattering_eirik/Eirik_j2 /tmp/oldbuild
python run_bigmemory.py /tmp/oldbuild/bigmem/Main_opt.exe ~/cummulant_scattering/src
python run_eirik_j2.py  /tmp/oldbuild/eirik/Main.exe /tmp/oldbuild/eirik_repaired/Main.exe ~/cummulant_scattering/src
```

About 5 minutes on 2 cores. Needs gfortran, OpenBLAS and PySCF. BigMemory starts 16 threads
regardless of `OMP_NUM_THREADS`.
