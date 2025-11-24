module Integrals

"""
    IntegralGeometry(px, py, pz)

Lightweight container holding precomputed Cartesian offsets between centers.
Use it instead of relying on implicit module state; everything is passed
explicitly to the integration kernels for clarity and easier testing.
"""
struct IntegralGeometry
    px::Matrix{Float64}
    py::Matrix{Float64}
    pz::Matrix{Float64}
end

# -- Basic special functions -------------------------------------------------

"""
    doublefactorial(n::Integer) -> Int64

Compute the double factorial `(2n+1)!!` used in the small-argument expansion
of the spherical Bessel functions. The implementation is iterative to avoid
recursion overhead and overflows for modest orders.
"""
function doublefactorial(n::Integer)
    n < 0 && throw(ArgumentError("n must be non-negative"))
    prod = 1
    for k in n:-2:1
        prod *= k
    end
    return prod
end

"""
    spherical_besselj_series(order::Integer, x::Float64) -> Vector{Float64}

Return the spherical Bessel functions `j_0(x) … j_order(x)` using a numerically
stable recurrence. For small `x` the series


is used to avoid catastrophic cancellation. Larger arguments fall back to the
stable upward recurrence `j_{n+1}(x) = ((2n+1)/x) j_n(x) - j_{n-1}(x)` seeded
by the analytic `j₀` and `j₁` definitions.
"""
function spherical_besselj_series(order::Integer, x::Float64)
    order < 0 && throw(ArgumentError("order must be non-negative"))
    vals = zeros(Float64, order + 1)
    ax = abs(x)
    if ax < 1e-6
        x2 = x * x
        vals[1] = 1 - x2 / 6 + x2 * x2 / 120
        if order >= 1
            vals[2] = x / 3 - x2 * x / 30 + x2 * x2 * x / 840
        end
        for n in 2:order
            vals[n + 1] = (ax == 0) ? (n == 0 ? 1.0 : 0.0) : x^n / doublefactorial(2n + 1)
        end
        return vals
    end

    vals[1] = sin(x) / x
    if order == 0
        return vals
    end

    vals[2] = (sin(x) / (x * x)) - (cos(x) / x)
    for n in 1:(order - 1)
        vals[n + 2] = ((2n + 1) / x) * vals[n + 1] - vals[n]
    end
    return vals
end

"""
    hermite_like_coeffs(lmax::Integer, h::Float64) -> Matrix{Float64}

Compute the Hermite-like coefficient table `a` used in the derivative expansion
of the translation operator. The recurrence mirrors the Fortran implementation
but uses explicit iteration and allocation-free broadcasting where possible.
"""
function hermite_like_coeffs(lmax::Integer, h::Float64)
    lmax < 0 && throw(ArgumentError("lmax must be non-negative"))
    a = zeros(Float64, lmax + 1, lmax + 1)
    a[1, 1] = 1.0
    if lmax > 0
        a[2, 2] = -h
    end
    for L in 2:lmax
        a[L + 1, 1] = -(L - 1) * a[L - 1, 1]
        @inbounds for k in 1:lmax
            a[L + 1, k + 1] = -h * a[L, k] - (L - 1) * a[L - 1, k + 1]
        end
    end
    return a
end

"""
    rrdj0(lmax::Integer, x::Float64) -> Matrix{Float64}

Wrapper around [`hermite_like_coeffs`](@ref) for parity with the Fortran
`rrdj0` routine. The output has size `(lmax+1, lmax+1)` and can be reused
across many integral evaluations.
"""
rrdj0(lmax::Integer, x::Float64) = hermite_like_coeffs(lmax, x)

"""
    bessel_deriv(ll::Integer, mm::Integer, nn::Integer, a, b, c) -> Vector{Float64}

Evaluate the derivative coefficients `BD` described in equation

````
BDₕ = Σ_{i=0}^{ℓ} Σ_{j=0}^{m} Σ_{k=0}^{n} a_{ℓ,i} b_{m,j} c_{n,k}
       with h-order ceil((ℓ+m+n-i-j-k)/2) + i + j + k
````

`a`, `b`, and `c` are the Hermite-like tables returned by [`rrdj0`](@ref) for
`h_x`, `h_y`, and `h_z` respectively. The computation short-circuits whenever
an intermediate coefficient is smaller than `1e-30` to keep the accumulation
numerically stable.
"""
function bessel_deriv(ll::Integer, mm::Integer, nn::Integer, a, b, c)
    ll < 0 && throw(ArgumentError("orders must be non-negative"))
    mm < 0 && throw(ArgumentError("orders must be non-negative"))
    nn < 0 && throw(ArgumentError("orders must be non-negative"))
    llmax = size(a, 1) - 1
    BD = zeros(Float64, llmax + 1)
    for ii in 0:ll
        c1 = a[ll + 1, ii + 1]
        abs(c1) < 1e-30 && continue
        for jj in 0:mm
            ct2 = b[mm + 1, jj + 1]
            abs(ct2) < 1e-30 && continue
            c2 = c1 * ct2
            for kk in 0:nn
                ct3 = c[nn + 1, kk + 1]
                abs(ct3) < 1e-30 && continue
                temp = ll + mm + nn - ii - jj - kk
                horder = temp ÷ 2 + (temp % 2) + ii + jj + kk
                BD[horder + 1] += c2 * ct3
            end
        end
    end
    return BD
end

"""
    bessel_series_sum(mu::AbstractVector{<:Real}, h::Real, h_saved::AbstractVector{<:Real})

Vectorized version of the `bessels0rr` summation. For each `μᵢ` the routine
computes

````
F(μᵢ) = Σ_{β=0}^{L} (μᵢ/h)^β j_β(μᵢ h) h_saved[β+1]
````

where `j_β` denotes the spherical Bessel function of the first kind. The
computation reuses the recurrence from [`spherical_besselj_series`](@ref) and
avoids division near the origin through the tailored series expansion.
"""
function bessel_series_sum(mu::AbstractVector{<:Real}, h::Real, h_saved::AbstractVector{<:Real})
    order = length(h_saved) - 1
    h <= 0 && throw(ArgumentError("h must be positive"))
    f = zeros(Float64, length(mu))
    for (idx, μ) in pairs(mu)
        pmu = float(μ) * h
        bvals = spherical_besselj_series(order, pmu)
        scale = (float(μ) / h) .^ collect(0:order)
        f[idx] = sum(scale .* bvals .* h_saved)
    end
    return f
end

# -- Geometry helpers --------------------------------------------------------

"""
    pairwise_offsets(geom::IntegralGeometry, i, j, k, r)

Return `hx`, `hy`, `hz`, and their Euclidean norm `h` for a quartet of centers
indexed by `i`, `j`, `k`, and `r` within the geometry arrays. This helper keeps
call sites free from COMMON-like implicit access patterns.
"""
function pairwise_offsets(geom::IntegralGeometry, i, j, k, r)
    hx = geom.px[k, r] - geom.px[i, j]
    hy = geom.py[k, r] - geom.py[i, j]
    hz = geom.pz[k, r] - geom.pz[i, j]
    h = hypot(hx, hypot(hy, hz))
    return hx, hy, hz, h
end

# -- Main integral kernel ----------------------------------------------------

"""
    tot_integral_k_ijkr(mu, l, m, n, group_start, group_count,
                       hx, hy, hz, h, dx1, dy1, dz1, dx2, dy2, dz2,
                       gi, gj, gk, gr, zcontr, zcontr2;
                       cutoff1=1e-12, cutoff2=1e-12)

Julia counterpart of the `tot_integral_k_ijkr` Fortran subroutine. The routine
accumulates mixed-derivative prefactors (stored in `dx*`, `dy*`, `dz*`) against
contracted density matrices `zcontr` and `zcontr2` before projecting them onto
a radial Bessel series. All inputs are explicit; no COMMON blocks or global
state are used.

The returned vector has the same length as `mu` and is computed via the
spherical-Bessel expansion

````
F(μ) = Σ_{β} (μ/h)^β j_β(μh) h_saved[β]
````

where `h_saved` collects every valid derivative contribution exceeding
`cutoff2`. Small `h` values remain safe because the Bessel evaluation falls
back to a series expansion.

# Examples

```julia
julia> using .Integrals
julia> l = m = n = fill(0, 2);
julia> groups = (group_start = [1, 2], group_count = [1, 1]);
julia> dx = dy = dz = ones(1, 1, 1);
julia> zc = fill(0.5, 1, 1, 1, 1);
julia> res = tot_integral_k_ijkr([0.1, 0.2], l, m, n, groups.group_start, groups.group_count,
                                0.05, 0.02, 0.01, 0.06,
                                dx, dy, dz, dx, dy, dz,
                                1, 1, 2, 2, zc, zc);
julia> all(res .> 0)
true
```
"""
function tot_integral_k_ijkr(mu, l, m, n, group_start, group_count,
                             hx, hy, hz, h, dx1, dy1, dz1, dx2, dy2, dz2,
                             gi::Integer, gj::Integer, gk::Integer, gr::Integer,
                             zcontr, zcontr2; cutoff1=1e-12, cutoff2=1e-12)
    # --- Input validation ---
    # Check group indices
    for idx in (gi, gj, gk, gr)
        if idx < 1 || idx > length(group_start)
            throw(ArgumentError("Group index $idx out of bounds for group_start of length $(length(group_start))"))
        end
        if idx > length(group_count)
            throw(ArgumentError("Group index $idx out of bounds for group_count of length $(length(group_count))"))
        end
    end
    # Check l, m, n arrays are long enough for all group_start + group_count - 1
    maxidx = maximum([group_start[gi] + group_count[gi] - 1,
                      group_start[gj] + group_count[gj] - 1,
                      group_start[gk] + group_count[gk] - 1,
                      group_start[gr] + group_count[gr] - 1])
    for arr in (l, m, n)
        if maxidx > length(arr)
            throw(ArgumentError("Basis index $maxidx out of bounds for array of length $(length(arr))"))
        end
    end
    # Check derivative arrays have sufficient size
    # (Assume they are at least 3D arrays; check size along each dimension >= max needed)
    for (name, arr) in zip(("dx1", "dy1", "dz1", "dx2", "dy2", "dz2"), (dx1, dy1, dz1, dx2, dy2, dz2))
        if ndims(arr) < 3
            throw(ArgumentError("Array $name must have at least 3 dimensions, got $(ndims(arr))"))
        end
    end
    # Check zcontr and zcontr2 have compatible shapes
    if size(zcontr) != size(zcontr2)
        throw(ArgumentError("zcontr and zcontr2 must have the same shape, got $(size(zcontr)) and $(size(zcontr2))"))
    end
    LLmax = l[group_start[gi]] + m[group_start[gi]] + n[group_start[gi]] +
            l[group_start[gj]] + m[group_start[gj]] + n[group_start[gj]] +
            l[group_start[gk]] + m[group_start[gk]] + n[group_start[gk]] +
            l[group_start[gr]] + m[group_start[gr]] + n[group_start[gr]]
    a = rrdj0(LLmax, hx)
    b = rrdj0(LLmax, hy)
    c = rrdj0(LLmax, hz)

    h_saved = zeros(Float64, LLmax + 1)

    posI = 1
    for i in group_start[gi]:(group_start[gi] + group_count[gi] - 1)
        posJ = 1
        for j in group_start[gj]:(group_start[gj] + group_count[gj] - 1)
            posK = 1
            for k in group_start[gk]:(group_start[gk] + group_count[gk] - 1)
                posR = 1
                for r in group_start[gr]:(group_start[gr] + group_count[gr] - 1)
                    ztot = zcontr[posJ, posK, posR, posI] + zcontr2[posR, posI, posJ, posK]
                    if abs(ztot) >= cutoff1
                        for ll in 0:(l[i] + l[j])
                            mdl = dx1[ll + 1, l[i] + 1, l[j] + 1]
                            prod1 = mdl * ztot
                            for mm in 0:(m[i] + m[j])
                                mdm = dy1[mm + 1, m[i] + 1, m[j] + 1]
                                prod2 = mdm * prod1
                                for nn in 0:(n[i] + n[j])
                                    mdn = dz1[nn + 1, n[i] + 1, n[j] + 1]
                                    phase = isodd(ll + mm + nn) ? -1.0 : 1.0
                                    prod3 = mdn * prod2 * phase
                                    for llp in 0:(l[k] + l[r])
                                        mdlp = dx2[llp + 1, l[k] + 1, l[r] + 1]
                                        prod4 = mdlp * prod3
                                        for mmp in 0:(m[k] + m[r])
                                            mdmp = dy2[mmp + 1, m[k] + 1, m[r] + 1]
                                            prod5 = mdmp * prod4
                                            for nnp in 0:(n[k] + n[r])
                                                mdnp = dz2[nnp + 1, n[k] + 1, n[r] + 1]
                                                prod6 = mdnp * prod5
                                                if abs(prod6) > cutoff2
                                                    bd = bessel_deriv(ll + llp, mm + mmp, nn + nnp, a, b, c)
                                                    h_saved .+= prod6 .* bd
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    posR += 1
                end
                posK += 1
            end
            posJ += 1
        end
        posI += 1
    end

    return bessel_series_sum(mu, h, h_saved)
end

"""
    _selftest()

Lightweight runtime checks used by CI scripts and developers to verify that the
Julia counterparts behave sensibly. The assertions deliberately reuse the
examples in the docstrings.
"""
function _selftest()
    l = m = n = fill(0, 2)
    groups = (group_start = [1, 2], group_count = [1, 1])
    dx = dy = dz = ones(1, 1, 1)
    zc = fill(0.5, 1, 1, 1, 1)
    res = tot_integral_k_ijkr([0.1, 0.2], l, m, n, groups.group_start, groups.group_count,
                              0.05, 0.02, 0.01, 0.06,
                              dx, dy, dz, dx, dy, dz,
                              1, 1, 2, 2, zc, zc)
    @assert all(res .> 0)
    bvals = spherical_besselj_series(3, 1e-7)
    @assert isapprox(bvals[1], 1; atol=1e-12)
    @assert isapprox(bvals[2], (1e-7 / 3); atol=1e-12)
    return true
end

end # module
