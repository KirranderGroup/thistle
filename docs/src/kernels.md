# Kernel guides

This page groups the Julia counterparts of the numerical kernels into **total**
and **elastic** families. Examples emphasize preserving the q-grid dimensionality
and the angular-momentum bookkeeping inherited from the Fortran drivers.

## Total routines

`Integrals.tot_integral_k_ijkr` mirrors the nested Fortran loops while keeping
all inputs explicit. Its internal `h_pre2` cache preserves the symmetry between
`zcontr` and `zcontr2` by reusing the same derivative grid for both branches.

```jldoctest
julia> include("../src/Integrals.jl"); using .Integrals

julia> l = m = n = fill(0, 2);

julia> groups = (group_start = [1, 2], group_count = [1, 1]);

julia> dx = dy = dz = ones(1, 1, 1);

julia> zc = fill(0.5, 1, 1, 1, 1);

julia> mu = [0.05, 0.15, 0.25];

julia> res = tot_integral_k_ijkr(mu, l, m, n, groups.group_start, groups.group_count,
                                0.1, 0.0, 0.0, 0.1,
                                dx, dy, dz, dx, dy, dz,
                                1, 1, 2, 2, zc, zc);

julia> size(res)
(3,)
```

When a cached set of derivative prefactors is already available, it can be
projected on a q-grid without touching the geometry by calling
`Integrals.bessel_sum`, which uses `SpecialFunctions.sphericalbesselj` under the
hood:

```jldoctest
julia> h_saved = [1.0, 0.1, 0.01];

julia> bvals = bessel_sum(mu, 0.4, h_saved);

julia> all(size(bvals) .== size(mu))
true
```

## Elastic routines

The elastic kernels `bessel_deriv1j` and `bessel_deriv2j` translate the
Fortran `BesselDeriv1j`/`2j` routines while maintaining the angular-momentum
ranges `(ℓ, m, n)` and their symmetry-optimized iteration order.

```jldoctest
julia> a = b = c = rrdj0(2, 0.05);

julia> ap = bp = cp = rrdj0(2, 0.1);

julia> bd1 = bessel_deriv1j(0, 0, 0, 0.0, 0.1, a, b, c);

julia> bd2 = bessel_deriv2j(0, 0, 0, 0.0, 0.1, a, b, c, ap, bp, cp);

julia> length(bd1) >= 3 && length(bd2) >= 4
true
```

Because the derivative arrays reuse the `rrdj0` output, the symmetry between
coordinates is retained even when the geometry is perturbed (e.g., different
`hx`, `hy`, and `hz` values for `a`, `b`, and `c`). The resulting Bessel
prefactors can be combined with `bessel_sum` to evaluate elastic scattering
factors on arbitrary q-grids without reshaping intermediate arrays.
