"""
    linspace(from::Real, to::Real, n::Integer) -> Vector{Float64}

Return `n` linearly spaced values between `from` and `to`, inclusive.

# Arguments
- `from::Real`: Starting value of the range.
- `to::Real`: Ending value of the range.
- `n::Integer`: Number of points to generate (must be at least 2).

# Returns
A `Vector{Float64}` containing `n` evenly spaced samples from `from` to `to`.

# Examples
```julia
julia> linspace(0.0, 1.0, 5)
5-element Vector{Float64}:
 0.0
 0.25
 0.5
 0.75
 1.0
```
"""
function linspace(from::Real, to::Real, n::Integer)
    n < 2 && throw(ArgumentError("n must be at least 2 to include both endpoints"))
    collect(range(start = float(from), stop = float(to), length = n))
end

