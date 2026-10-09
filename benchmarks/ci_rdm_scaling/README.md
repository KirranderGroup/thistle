# CI wavefunctions: where the memory goes

Question (7 October): does a large CI wavefunction (FCI, say) explode the Python side of the
2RDM, which is why the earliest version did it in Fortran?

## 1. CI vector → active-space RDMs (`builders.py`)

Random full-CAS vectors (cost does not depend on the values); 2 cores, 7 GB.

| Active space | Determinants | dense | sparse | bitwise (Fortran pair loop) |
|---|---|---|---|---|
| CAS(10,10) | 6.4·10⁴ | 0.17 s / 0.16 GB | 1.0 s / 0.4 GB | 1.7 s / 0.11 GB |
| CAS(12,12) | 8.5·10⁵ | 2.7 s / 1.2 GB | 37 s / 3.1 GB | 216 s / 0.22 GB |
| FCI H₂O/6-31G (10e,13o) | 1.7·10⁶ | 4.9 s / 1.4 GB | 94 s / 3.3 GB | ≈ 13 min (N_det² estimate) |
| CAS(14,14) | 1.2·10⁷ | 39 s / 3.2 GB | — | ≈ 11 h |
| CAS(16,16) | 1.7·10⁸ | killed at 7 GB: the vector (1.3 GB) plus the determinant list of this script | | |

A pure-Python pair loop costs 170 ns per pair for the screening alone (`popcount(xor)`, before
any Slater–Condon work): 0.1 h at 6·10⁴ determinants, 17 h at 8.5·10⁵, 64 h at 1.7·10⁶. That is
why the earliest Python version was limited to ~10³ determinants. None of the present builders
loops in Python: dense and sparse are BLAS/SciPy, bitwise is Fortran. The result is tiny
(n_act⁴: 0.3 MB for 14 orbitals).

## 2. Active RDMs → what the kernels read (`casci_ao.py`)

Benzene CASCI(6,6), total scattering at 30 q:

| Basis (Cartesian AO) | Route | 2RDM held | Peak memory | Build | Kernel |
|---|---|---|---|---|---|
| cc-pVDZ (120) | PySCF `make_rdm12` (all MOs, then AO) | 1.66 GB | 5.9 GB | 10 s | 19.9 s |
| cc-pVDZ (120) | four-index array over core + active, then AO | 1.66 GB | 5.9 GB | 8 s | 19.5 s |
| cc-pVDZ (120) | **factored** (24 + 3 terms) | 6 MB | 0.14 GB | 0.03 s | 28 s |
| cc-pVTZ (300) | four-index array | 65 GB | > 200 GB (est.) | — | — |
| cc-pVTZ (300) | **factored** | 39 MB | 0.33 GB | 0.1 s | 650 s |
| cc-pVQZ (630) | four-index array | 1.26 TB | — | — | — |
| cc-pVQZ (630) | factored (est.) | 0.17 GB | < 1 GB | | |

CAS(12,12) on cc-pVDZ: 81 + 3 terms, kernel 28.9 s, the same as with 24 terms: the cost is in
the integrals, not in the number of terms.

So the explosion was never the CI → RDM step; it was the n_AO⁴ array built afterwards. The
factored embedding (`pyxscat.ci.embed_core(..., factored=True)`) removes it and is the default
beyond 60 basis functions.

Run: `python builders.py 12 6 6 dense`, `python casci_ao.py prepare cc-pvdz`,
`python casci_ao.py run cc-pvdz factored` (with pyxscat importable).
