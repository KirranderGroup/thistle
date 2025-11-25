# Configuration interaction helpers

Thistle's CI utilities translate the Fortran reduced density matrix builders
into a compact Julia interface. Configurations are represented as matrices with
one row per configuration and one column per spin orbital. A value of `0`
indicates an empty spin orbital while `1` and `2` encode spin-up and spin-down
occupations respectively. For a two-orbital singlet, for example, a row such as
`[1 0 2 0]` means `α` occupancy on orbital 1 and `β` occupancy on orbital 2.

## Phase and difference bookkeeping

The Fortran drivers allocate two integer matrices to track relative phases and
occupation differences between configurations:

- `ep2` records the sign changes produced when reordering occupied spin
  orbitals between configurations.
- `ndiff` counts the number of spin-orbital positions that differ.

Both are created in Julia by [`maxcoincidence`](@ref) and returned as regular
`Matrix{Int}` arrays so they can be passed directly to the density matrix
routines.

```jldoctest
julia> using CIUtils

julia> confs = [
           1 0 2 0;
           0 1 2 0;
       ];

julia> ep2, ndiff = maxcoincidence(confs);

julia> ep2
2×2 Matrix{Int64}:
 1  1
 1  1

julia> ndiff
2×2 Matrix{Int64}:
 0  2
 2  0
```

## Reduced density matrices

[`create_onerdm`](@ref) and [`create_twordm`](@ref) mirror the sign and phase
logic of the legacy Fortran subroutines while relying on Julia's explicit array
allocations. They accept CI vectors column-wise via the `civs` matrix and the
phase/occupancy arrays computed above.

```jldoctest
julia> using CIUtils

julia> confs = [
           1 0 1 0;  # |α₁ α₂⟩
           0 1 1 0;  # |α₂ α₂⟩
       ];

julia> civs = [0.8 0.1; 0.6 0.2];

julia> ep2, ndiff = maxcoincidence(confs);

julia> onerdm = create_onerdm(confs, civs, ndiff, ep2);

julia> onerdm
2×2 Matrix{Float64}:
 1.96  0.0
 0.0   1.0

julia> twordm = create_twordm(confs, civs, ndiff, ep2);

julia> twordm[1, 1, 1, 1], twordm[2, 1, 1, 2]
(0.0, 1.0)
```

The examples above exercise single and double excitations in a minimal CI space
and serve as regression tests for the ported sign conventions.
