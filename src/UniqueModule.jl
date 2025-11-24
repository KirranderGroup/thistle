module UniqueModule

export unique_integer, unique_real, unique_total

# Internal helper to build stable unique row information.
function _unique_rows(matrix::AbstractMatrix, key_fn)
    nrows, _ = size(matrix)
    irec = zeros(Int, nrows)
    iuni = Int[]
    row_to_index = Dict{Any, Int}()

    for i in 1:nrows
        row_view = @view matrix[i, :]
        key = key_fn(row_view)
        group = get(row_to_index, key, 0)

        if group == 0
            group = length(iuni) + 1
            row_to_index[key] = group
            push!(iuni, i)
        end

        irec[i] = group
    end

    rmat = matrix[iuni, :]
    return rmat, iuni, irec
end

"""
    unique_integer(matrix::AbstractMatrix{<:Integer}) -> (rmat, iuni, irec)

Return unique rows from an integer matrix while preserving the first occurrence of
each row.

The algorithm scans rows from top to bottom, keeping the earliest occurrence of each
row as the representative. `iuni` lists the original row indices for the retained
rows (in the order they first appear), `rmat` contains those rows in the same order,
and `irec[i]` marks which representative row each original row belongs to.

Ordering is stable because rows are recorded the first time they are encountered and
never reordered afterwards.

# Examples
```jldoctest
julia> using UniqueModule

julia> mat = [1 2; 3 4; 1 2; 3 4];

julia> rmat, iuni, irec = unique_integer(mat);

julia> rmat
2×2 Matrix{Int64}:
 1  2
 3  4

julia> iuni
2-element Vector{Int64}:
 1
 2

julia> irec
4-element Vector{Int64}:
 1
 2
 1
 2
```
"""
function unique_integer(matrix::AbstractMatrix{<:Integer})
    _unique_rows(matrix, row -> Tuple(row))
end

"""
    unique_real(matrix::AbstractMatrix{<:Real}) -> (rmat, iuni, irec)

Return unique rows from a real-valued matrix with stable ordering.

Rows are considered identical only when every element matches exactly (no tolerance).
`iuni` records the first index of each distinct row, `rmat` contains the
corresponding rows in their discovery order, and `irec` maps each original row to its
representative.

# Examples
```jldoctest
julia> using UniqueModule

julia> mat = [1.0 2.0; 1.0 2.0; 1.0 3.0];

julia> rmat, iuni, irec = unique_real(mat);

julia> rmat
2×2 Matrix{Float64}:
 1.0  2.0
 1.0  3.0

julia> iuni
2-element Vector{Int64}:
 1
 3

julia> irec
3-element Vector{Int64}:
 1
 1
 2
```
"""
function unique_real(matrix::AbstractMatrix{<:Real})
    _unique_rows(matrix, row -> Tuple(row))
end

"""
    unique_total(matrix::AbstractMatrix{<:Integer}, total::AbstractVector{<:Real})
        -> (rmat, rtot, iuni, irec)

Collapse duplicate integer rows while summing associated totals.

Rows are de-duplicated in the order they first appear, mirroring `unique_integer`.
`rtot` stores the sum of `total` values for each unique row, aligned with the rows in
`rmat`. The vectors `iuni` and `irec` record the first-occurrence indices and group
membership, respectively.

# Examples
```jldoctest
julia> using UniqueModule

julia> mat = [1 2; 3 4; 1 2];

julia> totals = [10.0, 20.0, 1.5];

julia> rmat, rtot, iuni, irec = unique_total(mat, totals);

julia> rmat
2×2 Matrix{Int64}:
 1  2
 3  4

julia> rtot
2-element Vector{Float64}:
 11.5
 20.0

julia> iuni
2-element Vector{Int64}:
 1
 2

julia> irec
3-element Vector{Int64}:
 1
 2
 1
```
"""
function unique_total(matrix::AbstractMatrix{<:Integer}, total::AbstractVector{<:Real})
    nrows, _ = size(matrix)
    length(total) == nrows || throw(ArgumentError("total must have one entry per row"))

    rmat, iuni, irec = _unique_rows(matrix, row -> Tuple(row))

    totals = zeros(promote_type(eltype(total), Float64), length(iuni))
    for (row_idx, group_idx) in enumerate(irec)
        totals[group_idx] += total[row_idx]
    end

    return rmat, totals, iuni, irec
end

end
