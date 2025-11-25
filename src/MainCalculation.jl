module MainCalculation

if !isdefined(Main, :Types)
    Base.include(Main, joinpath(@__DIR__, "Types.jl"))
end
if !isdefined(Main, :CIUtils)
    Base.include(Main, joinpath(@__DIR__, "CIUtils.jl"))
end
if !isdefined(Main, :Integrals)
    Base.include(Main, joinpath(@__DIR__, "Integrals.jl"))
end

const TypesMod = Main.Types
const CIUtilsMod = Main.CIUtils
const IntegralsMod = Main.Integrals

using .TypesMod: DP, IK
using .CIUtilsMod
using .IntegralsMod: IntegralGeometry, form_factor, pairwise_offsets, tot_integral_k_ijkr
using Base.Threads

export total_scattering_calculation, group_metadata

"""
    group_metadata(group::AbstractVector{<:Integer}) -> (group_start, group_count)

Recreate the `group_start`/`group_count` bookkeeping arrays produced by the
Fortran drivers. `group_start[i]` stores the first index whose group label
matches `i` and `group_count[i]` records how many consecutive entries share that
label. Missing groups are represented with zeros, mirroring the implicit COMMON
layout in the original routines.
"""
function group_metadata(group::AbstractVector{<:Integer})
    max_group = maximum(group)
    starts = zeros(IK, max_group)
    counts = zeros(IK, max_group)

    for g in 1:max_group
        idxs = findall(==(g), group)
        if !isempty(idxs)
            starts[g] = first(idxs)
            counts[g] = length(idxs)
        end
    end

    return (group_start = starts, group_count = counts)
end

# -- Density helpers ----------------------------------------------------------

"""
    density_from_ci(type, confs, civs; state1=1, state2=1)

Build the contracted density tensors used by the scattering kernels. For
`type == 1` (total scattering) the spin-free 2-RDM is returned twice to feed
both branches of the total integral. For `type == 2` (elastic scattering) the
1-RDM is expanded to a four-index tensor by placing the `(p, q)` element at
`z[p, 1, 1, q]`, matching the `zcontr[posJ, posK, posR, posI]` indexing in the
integration kernel.
"""
function density_from_ci(type::Integer, confs, civs; state1::Integer=1, state2::Integer=1)
    ep2, ndiff = maxcoincidence(confs)

    if type == 1
        z = create_twordm(confs, civs, ndiff, ep2; state1=state1, state2=state2)
        return z, z
    elseif type == 2
        onerdm = create_onerdm(confs, civs, ndiff, ep2; state1=state1, state2=state2)
        norb = size(onerdm, 1)
        z = zeros(DP, norb, 1, 1, norb)
        @inbounds for i in 1:norb, j in 1:norb
            z[j, 1, 1, i] = onerdm[i, j]
        end
        return z, z
    else
        throw(ArgumentError("type must be 1 (total) or 2 (elastic), got $type"))
    end
end

# -- Core accumulation --------------------------------------------------------

function accumulate_quartets!(accum, mu, geom, l, m, n, group_start, group_count,
                              ddx, ddy, ddz, zcontr, zcontr2;
                              unit_scale=1.0, cutoff1=1e-12, cutoff2=1e-12)
    ng = length(group_start)
    partials = [zeros(eltype(accum), length(accum)) for _ in 1:nthreads()]

    @threads for gi in 1:ng
        tid = threadid()
        local_accum = partials[tid]
        for gj in 1:ng, gk in 1:ng, gr in 1:ng
            hx, hy, hz, h = pairwise_offsets(geom, gi, gj, gk, gr)
            h_scaled = unit_scale * h
            if h_scaled <= 0
                continue
            end

            local_accum .+= tot_integral_k_ijkr(mu, l, m, n, group_start, group_count,
                                                unit_scale * hx, unit_scale * hy, unit_scale * hz, h_scaled,
                                                ddx, ddy, ddz, ddx, ddy, ddz,
                                                gi, gj, gk, gr, zcontr, zcontr2;
                                                cutoff1=cutoff1, cutoff2=cutoff2)
        end
    end

    for part in partials
        accum .+= part
    end

    return accum
end

# -- Public API ---------------------------------------------------------------

"""
    total_scattering_calculation(type, q, geom, l, m, n, group, ddx, ddy, ddz;
                                 zcontr=nothing, zcontr2=nothing,
                                 confs=nothing, civs=nothing,
                                 state1::Integer=1, state2::Integer=1,
                                 cutoff1=1e-12, cutoff2=1e-12, unit_scale=1.0)

Julia translation of the Fortran `total_scattering_calculation` driver. All
implicit COMMON state is replaced by explicit arguments:

* `type` selects the scattering channel (`1` for total, `2` for elastic).
* `q` is the momentum-transfer grid (one-dimensional).
* `geom` holds the precomputed pairwise offsets in an [`IntegralGeometry`](@ref).
* `l`, `m`, `n` store angular momentum exponents for each basis function.
* `group` labels each basis function with its symmetry group.
* `ddx`, `ddy`, `ddz` are the derivative prefactor tables shared by both
  branches of the total integral.

Density tensors can be provided directly through `zcontr`/`zcontr2`. When they
are omitted, the routine builds them from CI configurations (`confs`) and
coefficients (`civs`) using [`CIUtils.create_twordm`](@ref) for total scattering
or [`CIUtils.create_onerdm`](@ref) for elastic scattering. The return value is a
vector with the same length as `q`.
"""
function total_scattering_calculation(type::Integer, q, geom::IntegralGeometry,
                                      l, m, n, group, ddx, ddy, ddz;
                                      zcontr=nothing, zcontr2=nothing,
                                      confs=nothing, civs=nothing,
                                      state1::Integer=1, state2::Integer=1,
                                      cutoff1=1e-12, cutoff2=1e-12, unit_scale=1.0)
    groups = group_metadata(group)

    dens1 = zcontr
    dens2 = zcontr2
    if dens1 === nothing || dens2 === nothing
        (dens1, dens2) = density_from_ci(type, confs, civs; state1=state1, state2=state2)
    end

    result = zeros(DP, length(q))
    return accumulate_quartets!(result, q, geom, l, m, n, groups.group_start, groups.group_count,
                                ddx, ddy, ddz, dens1, dens2;
                                unit_scale=unit_scale, cutoff1=cutoff1, cutoff2=cutoff2)
end

end
