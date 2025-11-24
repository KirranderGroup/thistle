# Numeric conventions

Thistle's Julia tools standardize numeric types through the `Types` module located in `src/Types.jl`:

- `DP` is the default floating-point precision (an alias for `Float64`).
- `IK` is the default integer kind (an alias for `Int64`).

All Julia entry points import these aliases to keep type choices consistent between scripts. The aliases keep the codebase concise while making it easy to experiment with different numeric widths.

## Adjusting precision

If you want to run calculations at a different precision, update the aliases in `src/Types.jl` and rebuild any compiled artifacts:

```julia
module Types
    const DP = BigFloat   # or Float32, etc.
    const IK = Int64      # keep or switch to a narrower integer type
end
```

Using a lower precision (e.g., `Float32`) can improve performance at the cost of accuracy, while higher precision (`BigFloat`) increases numerical stability but may slow execution. Because the aliases are imported in each Julia script, changing them in one place updates the numeric behavior throughout the Julia tooling.

## Integral kernels

`Integrals.tot_integral_k_ijkr` mirrors the nested control flow of the Fortran
`tot_integral_k_ijkr` routine: it builds the Hermite-like derivative grid on the
simplex `ll + mm + nn ≤ LLmax` once and reuses it for every `(i, j, k, r)`
quartet. Unlike the higher-level Fortran drivers, the Julia entry point expects
the caller to provide derivative tables (`dx*`, `dy*`, `dz*`) and group metadata
explicitly rather than consulting COMMON blocks. This keeps the symmetry
between `zcontr` and `zcontr2` intact while remaining easy to exercise in unit
tests.
