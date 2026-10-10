"""
    DiscreteIntegration

Numerical integration of discrete data, given as sampled values `y` on an uneven grid `x` or an even grid with spacing `Δx`.

# Exports:
- `integrate`: definite integral over all points
- `cumul_integrate`: cumulative integral, returns a new vector
- `cumul_integrate!`: cumulative integral, inplace
- `Trapezoidal`, `Simpson`: integration methods
"""
module DiscreteIntegration

export integrate, cumul_integrate, cumul_integrate!, Trapezoidal, Simpson


# define integration methods
############################

abstract type IntegrationMethod end

"""
    Trapezoidal()

Integration method: trapezoidal rule. Exact for polynomials up to degree 1.
"""
struct Trapezoidal  <: IntegrationMethod end


"""
    Simpson()

Integration method: Simpson rule. Exact for polynomials up to degree 2 on uneven grids (`x`, `y`) and up to degree 3 on even grids (`Δx`, `y`).
If the number of points is even, the last interval is treated separately:
- uneven grid: Cartwright correction
- even grid: Simpson 3/8 rule

Requires at least 3 points.
"""
struct Simpson      <: IntegrationMethod end


# implementation
################

# integration type trapezoidal
include("trapezoidal.jl")

# integration type simpson
include("simpson.jl")


# function wrapper
###############

"""
    integrate(x, y, method = Trapezoidal())

Computes the numerical integral of y(x) over all points. The grid may be uneven.

# Arguments:
- `x::AbstractVector{<:Real}`: vector or range of x values
- `y::AbstractVector{<:Real}`: vector of corresponding y values
- `method::IntegrationMethod`: integration method, `Trapezoidal()` (default) or `Simpson()`

# Returns:
- scalar of at least Float64

# Throws:
- `DimensionMismatch`: if `x` and `y` differ in length
- `ArgumentError`: if fewer than 2 points (`Trapezoidal`) or 3 points (`Simpson`) are given

# Notes:
- `x` may be descending, the integral then changes sign.
- For `Simpson`, `x` must be strictly monotonic (no duplicate x values).
- Elements of `x` and `y` are loaded only once, so lazy arrays are supported.
- Vectors must have contiguous 1D indices.
"""
function integrate(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    method::IntegrationMethod = Trapezoidal()
) 
    _integrate(x, y, method)
end

"""
    integrate(Δx, y, method = Trapezoidal())

Computes the numerical integral of `y` over all points for an even grid spacing.

# Arguments:
- `Δx::Real`: constant spacing of the y values (may be negative)
- `y::AbstractVector{<:Real}`: vector of y values
- `method::IntegrationMethod`: integration method, `Trapezoidal()` (default) or `Simpson()`

# Returns:
- scalar of at least Float64

# Throws:
- `ArgumentError`: if fewer than 2 points (`Trapezoidal`) or 3 points (`Simpson`) are given

# Notes:
- `Δx = 0` returns 0.
- Elements of `y` are loaded only once, so lazy arrays are supported.
- Vectors must have contiguous 1D indices.
"""
function integrate(
    Δx::Real,
    y::AbstractVector{<:Real},
    method::IntegrationMethod = Trapezoidal()
)
    _integrate(Δx, y, method)
end

"""
    cumul_integrate(x, y, method = Trapezoidal())

Computes the cumulative numerical integral of y(x). The grid may be uneven.

# Arguments:
- `x::AbstractVector{<:Real}`: vector or range of x values
- `y::AbstractVector{<:Real}`: vector of corresponding y values
- `method::IntegrationMethod`: integration method, `Trapezoidal()` (default) or `Simpson()`

# Returns:
- `Vector` of at least Float64 with the same length as `y`. The first element is zero, element `k` is the integral from the first point to point `k`.

# Throws:
- `DimensionMismatch`: if `x` and `y` differ in length
- `ArgumentError`: if fewer than 2 points (`Trapezoidal`) or 3 points (`Simpson`) are given

# Notes:
- For `Simpson`, `x` must be strictly monotonic.
"""
function cumul_integrate(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    method::IntegrationMethod = Trapezoidal()
)   
    # get return type
    T = promote_type(eltype(x), eltype(y), Float64)

    # create dest vector
    dest = Vector{T}(undef, length(y))

    # calculate integral inplace
    _cumul_integrate!(dest, x, y, method)
end

"""
    cumul_integrate!(dest, x, y, method = Trapezoidal())

Computes the cumulative numerical integral of y(x) inplace. The grid may be uneven.

# Arguments:
- `dest::AbstractVector{<:AbstractFloat}`: output vector, same length as `y`
- `x::AbstractVector{<:Real}`: vector or range of x values
- `y::AbstractVector{<:Real}`: vector of corresponding y values
- `method::IntegrationMethod`: integration method, `Trapezoidal()` (default) or `Simpson()`

# Returns:
- `dest`, with `dest[1] = 0` and element `k` as the integral from the first point to point `k`

# Throws:
- `DimensionMismatch`: if `dest`, `x` and `y` differ in length
- `ArgumentError`: if fewer than 2 points (`Trapezoidal`) or 3 points (`Simpson`) are given

# Notes:
- Computation uses at least Float64 and is rounded to `eltype(dest)` on write. The caller is responsible for choosing a suitable `dest` type.
- `dest` may be identical to `x` or `y`. Partially overlapping views are not allowed.
- For `Simpson`, `x` must be strictly monotonic.
- Vectors must have contiguous 1D indices.
"""
function cumul_integrate!(
    dest::AbstractVector{<:AbstractFloat},
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    method::IntegrationMethod = Trapezoidal()
) 
    _cumul_integrate!(dest, x, y, method)
end

"""
    cumul_integrate(Δx, y, method = Trapezoidal())

Computes the cumulative numerical integral of `y` for an even grid spacing.

# Arguments:
- `Δx::Real`: constant spacing of the y values (may be negative)
- `y::AbstractVector{<:Real}`: vector of y values
- `method::IntegrationMethod`: integration method, `Trapezoidal()` (default) or `Simpson()`

# Returns:
- `Vector` of at least Float64 with the same length as `y`. The first element is zero, element `k` is the integral from the first point to point `k`.

# Throws:
- `ArgumentError`: if fewer than 2 points (`Trapezoidal`) or 3 points (`Simpson`) are given

# Notes:
- For `Simpson`, the second element is computed with the trapezoidal rule.
"""
function cumul_integrate(
    Δx::Real,
    y::AbstractVector{<:Real},
    method::IntegrationMethod = Trapezoidal()
)   
    # get return type
    T = promote_type(typeof(Δx), eltype(y), Float64)

    # create dest vector
    dest = Vector{T}(undef, length(y))

    # calculate integral inplace
    _cumul_integrate!(dest, Δx, y, method)
end

"""
    cumul_integrate!(dest, Δx, y, method = Trapezoidal())

Computes the cumulative numerical integral of `y` inplace for an even grid spacing.

# Arguments:
- `dest::AbstractVector{<:AbstractFloat}`: output vector, same length as `y`
- `Δx::Real`: constant spacing of the y values (may be negative)
- `y::AbstractVector{<:Real}`: vector of y values
- `method::IntegrationMethod`: integration method, `Trapezoidal()` (default) or `Simpson()`

# Returns:
- `dest`, with `dest[1] = 0` and element `k` as the integral from the first point to point `k`

# Throws:
- `DimensionMismatch`: if `dest` and `y` differ in length
- `ArgumentError`: if fewer than 2 points (`Trapezoidal`) or 3 points (`Simpson`) are given

# Notes:
- Computation uses at least Float64 and is rounded to `eltype(dest)` on write. The caller is responsible for choosing a suitable `dest` type.
- `dest` may be identical to `y`. Partially overlapping views are not allowed.
- For `Simpson`, the second element is computed with the trapezoidal rule.
- Vectors must have contiguous 1D indices.
"""
function cumul_integrate!(
    dest::AbstractVector{<:AbstractFloat},
    Δx::Real,
    y::AbstractVector{<:Real},
    method::IntegrationMethod = Trapezoidal()
) 
    _cumul_integrate!(dest, Δx, y, method)
end

end