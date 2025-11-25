module CIUtils

include("Types.jl")
using .Types: DP, IK

export excitations, maxcoincidence, create_onerdm, create_twordm

"""
    excitations(det1, det2)

Count the number of orbital excitations required to transform determinant
bitstrings `det1` into `det2`. Both determinants are represented as
`Nint × 2` integer arrays where each column stores the alpha and beta spin
bit strings respectively.

The implementation mirrors the Fortran `excitations` helper by summing the
population count of the pairwise XOR for each spin string and dividing by
`2` to obtain the number of distinct orbital substitutions.

# Examples
```jldoctest
julia> using CIUtils

julia> det1 = reshape(Int[0x3, 0x0], 1, 2); # |αα⟩

julia> det2 = reshape(Int[0x1, 0x2], 1, 2); # |αβ⟩

julia> excitations(det1, det2)
1
```
"""
function excitations(det1::AbstractMatrix{<:Integer}, det2::AbstractMatrix{<:Integer})
    size(det1) == size(det2) || throw(ArgumentError("determinants must have the same shape"))

    n_ex = 0
    for l in axes(det1, 1)
        n_ex += count_ones(det1[l, 1] ⊻ det2[l, 1])
        n_ex += count_ones(det1[l, 2] ⊻ det2[l, 2])
    end

    return n_ex >>> 1
end

# helper to convert a spin-orbital column index into the spatial orbital number
orbital_index(col::Integer) = (col + 1) >>> 1

"""
    maxcoincidence(confs) -> (ep2, ndiff)

Recreate the phase-tracking matrices produced by the Fortran CI drivers. The
`confs` matrix lists configurations as rows and spin-orbitals as columns using
`0` for unoccupied, `1` for spin-up and `2` for spin-down entries. The return
values are:

* `ep2` – a symmetric integer matrix containing the relative phase between
  configurations induced by reordering spin-orbitals into a shared ordering.
* `ndiff` – the number of differing spin-orbital occupations for every pair of
  configurations.

Both arrays are allocated explicitly so callers can write directly into the
buffers without worrying about Julia's lazy array semantics.

# Examples
```jldoctest
julia> using CIUtils

julia> confs = [
           1 0 2 0;  # |α₁ β₂⟩
           0 1 2 0;  # |α₂ β₂⟩
       ];

julia> ep2, ndiff = maxcoincidence(confs);

julia> ep2
2×2 Matrix{Int64}:
 1   1
 1   1

julia> ndiff
2×2 Matrix{Int64}:
 0  2
 2  0
```
"""
function maxcoincidence(confs::AbstractMatrix{<:Integer})
    nconf, _ = size(confs)
    ep2 = ones(Int, nconf, nconf)
    ndiff = zeros(Int, nconf, nconf)

    occupied = [findall(!iszero, confs[i, :]) for i in 1:nconf]

    for c1 in 1:nconf
        for c2 in (c1 + 1):nconf
            perm = collect(occupied[c2])
            ep = 1

            for (idx, pos1) in enumerate(occupied[c1])
                if idx > length(perm)
                    break
                end

                if perm[idx] != pos1
                    swap_at = findfirst(==(pos1), perm)
                    if swap_at !== nothing
                        ep = -ep
                        perm[idx], perm[swap_at] = perm[swap_at], perm[idx]
                    end
                end
            end

            ep2[c1, c2] = ep2[c2, c1] = ep
            diff = sum(@view(confs[c1, :]) .!= @view(confs[c2, :]))
            ndiff[c1, c2] = ndiff[c2, c1] = diff
        end
    end

    return ep2, ndiff
end

# Internal helper carrying the per-pair bookkeeping used by both density matrix
# builders. The signature mirrors the Fortran routine `maxc_individual`.
function _maxc_individual(c1::AbstractVector{<:Integer}, c2::AbstractVector{<:Integer})
    occ1 = findall(!iszero, c1)
    occ2 = findall(!iszero, c2)

    perm = collect(occ2)
    ep = 1
    for (idx, pos1) in enumerate(occ1)
        if idx > length(perm)
            break
        end

        if perm[idx] != pos1
            swap_at = findfirst(==(pos1), perm)
            if swap_at !== nothing
                ep = -ep
                perm[idx], perm[swap_at] = perm[swap_at], perm[idx]
            end
        end
    end

    ndiff = sum(c1 .!= c2)
    diffs1 = Int[]
    diffs2 = Int[]
    spin1 = Int[]
    spin2 = Int[]

    for idx in eachindex(c1)
        if c1[idx] != c2[idx]
            if c1[idx] != 0
                push!(diffs1, orbital_index(idx))
                push!(spin1, c1[idx])
            else
                push!(diffs2, orbital_index(idx))
                push!(spin2, c2[idx])
            end
        end
    end

    return ep, ndiff, diffs1, diffs2, spin1, spin2
end

"""
    create_onerdm(confs, civs, ndiff, ep2; state1=1, state2=1)

Build the spin-free one-particle reduced density matrix (1-RDM) from CI
configurations and coefficients. `confs` stores configurations as rows while
`civs` holds CI vectors column-wise, indexed by the `state1`/`state2`
arguments. The `ndiff` and `ep2` arrays are typically produced by
[`maxcoincidence`](@ref).

The routine keeps the original Fortran sign bookkeeping: double differences
(`ndiff == 2`) accumulate off-diagonal elements with the relative phase from
`ep2`, while identical configurations contribute to the diagonal.

# Examples
```jldoctest
julia> using CIUtils

julia> confs = [
           1 0 1 0;  # |α₁ α₂⟩
           0 1 1 0;  # |α₂ α₂⟩ (single excitation)
       ];

julia> civs = [0.8 0.1; 0.6 0.2];

julia> ep2, ndiff = maxcoincidence(confs);

julia> onerdm = create_onerdm(confs, civs, ndiff, ep2);

julia> onerdm
2×2 Matrix{Float64}:
 1.96  0.0
 0.0   1.0
```
"""
function create_onerdm(
    confs::AbstractMatrix{<:Integer},
    civs::AbstractMatrix{<:Real},
    ndiff::AbstractMatrix{<:Integer},
    ep2::AbstractMatrix{<:Integer};
    state1::Integer = 1,
    state2::Integer = 1,
)
    nconf, norb2 = size(confs)
    maxnmo = cld(norb2, 2)
    onerdm = zeros(DP, maxnmo, maxnmo)

    for c1 in 1:nconf
        for c2 in 1:nconf
            nd = ndiff[c1, c2]
            coeff = civs[c1, state1] * civs[c2, state2]
            ep = ep2[c1, c2]

            if nd == 2
                diffs1 = Int[]
                diffs2 = Int[]
                spin1 = Int[]
                spin2 = Int[]
                for n in 1:norb2
                    if confs[c1, n] != confs[c2, n]
                        if confs[c1, n] != 0
                            push!(diffs1, orbital_index(n))
                            push!(spin1, confs[c1, n])
                        elseif confs[c2, n] != 0
                            push!(diffs2, orbital_index(n))
                            push!(spin2, confs[c2, n])
                        end
                    end
                end

                if length(diffs1) == 1 && length(diffs2) == 1 && spin1 == spin2
                    onerdm[diffs1[1], diffs2[1]] += coeff * ep
                end
            elseif nd == 0
                for i1 in 1:norb2
                    if confs[c1, i1] != 0
                        orb = orbital_index(i1)
                        onerdm[orb, orb] += coeff
                    end
                end
            end
        end
    end

    return onerdm
end

"""
    create_twordm(confs, civs, ndiff, ep2; state1=1, state2=1)

Construct the spin-free two-particle reduced density matrix (2-RDM) following
the legacy Fortran control flow. The function handles three relevant cases:

* `ndiff == 4` – double excitations contribute four signed permutations based on
  the spin pattern of the departing and arriving electrons.
* `ndiff == 2` – single excitations pair the differing orbital with every
  matching-spin occupied orbital shared by both configurations.
* `ndiff == 0` – diagonal contributions accumulate both Coulomb and exchange
  terms for each occupied orbital pair.

The resulting tensor has shape `(norbs, norbs, norbs, norbs)` where `norbs` is
half the number of columns in `confs`.
"""
function create_twordm(
    confs::AbstractMatrix{<:Integer},
    civs::AbstractMatrix{<:Real},
    ndiff::AbstractMatrix{<:Integer},
    ep2::AbstractMatrix{<:Integer};
    state1::Integer = 1,
    state2::Integer = 1,
)
    nconf, norb2 = size(confs)
    norbs = cld(norb2, 2)
    twordm = zeros(DP, norbs, norbs, norbs, norbs)

    for c1 in 1:nconf
        for c2 in 1:nconf
            nd = ndiff[c1, c2]
            ep_pair = ep2[c1, c2]
            coeff = civs[c1, state1] * civs[c2, state2]

            if nd == 4
                ep, _, diffs1, diffs2, spin1, spin2 = _maxc_individual(
                    @view(confs[c1, :]), @view(confs[c2, :]),
                )

                if length(diffs1) == 2 && length(diffs2) == 2
                    # direct term
                    if spin2[1] == spin1[1] && spin2[2] == spin1[2]
                        twordm[diffs1[2], diffs1[1], diffs2[1], diffs2[2]] +=
                            coeff * ep * ep_pair
                    end
                    # exchange permutations
                    if spin2[1] == spin1[2] && spin2[2] == spin1[1]
                        twordm[diffs1[1], diffs1[2], diffs2[1], diffs2[2]] -=
                            coeff * ep * ep_pair
                        twordm[diffs1[2], diffs1[1], diffs2[2], diffs2[1]] -=
                            coeff * ep * ep_pair
                    end
                    if spin2[1] == spin1[1] && spin2[2] == spin1[2]
                        twordm[diffs1[1], diffs1[2], diffs2[2], diffs2[1]] +=
                            coeff * ep * ep_pair
                    end
                end
            elseif nd == 2
                ep, _, diffs1, diffs2, spin1, spin2 = _maxc_individual(
                    @view(confs[c1, :]), @view(confs[c2, :]),
                )

                if length(diffs1) == 1 && length(diffs2) == 1
                    removed = diffs1[1]
                    added = diffs2[1]
                    spin_removed = spin1[1]
                    spin_added = spin2[1]

                    # pair the differing orbital with every like-spin shared occupation
                    for idx in 1:norb2
                        if confs[c1, idx] != 0 && confs[c2, idx] != 0
                            common_orb = orbital_index(idx)
                            spin_common = confs[c1, idx]

                            if spin_common == spin_removed && spin_added == spin_removed
                                twordm[removed, common_orb, common_orb, added] +=
                                    coeff * ep * ep_pair
                                twordm[common_orb, removed, common_orb, added] -=
                                    coeff * ep * ep_pair
                                twordm[common_orb, common_orb, removed, added] +=
                                    coeff * ep * ep_pair
                            end
                        end
                    end
                end
            elseif nd == 0
                ep = ep_pair
                for i1 in 1:norb2
                    if confs[c1, i1] == 0
                        continue
                    end
                    for i2 in 1:norb2
                        if i1 == i2 || confs[c2, i2] == 0
                            continue
                        end
                        sorb = orbital_index(i1)
                        qorb = orbital_index(i2)

                        porb = qorb
                        rorb = sorb
                        twordm[porb, rorb, sorb, qorb] += coeff * ep

                        if confs[c1, i1] == confs[c1, i2]
                            porb = sorb
                            rorb = qorb
                            twordm[porb, rorb, sorb, qorb] -= coeff * ep
                        end
                    end
                end
            end
        end
    end

    return twordm
end

end
