
# src/trapezoidal.jl

"""
    _integrate(
        x::AbstractVector{<:Real}, 
        y::AbstractVector{<:Real},
        ::Trapezoidal
    )

Uses Trapezoidal rule, to compute numerical integral of y(x) over x.

# Arguments:
- `x::AbstractVector{<:Real}`: vector or range of x values
- `y::AbstractVector{<:Real}`: vector of corresponding y values
- `::Trapezoidal`: integration method

# Returns:
- scalar of at least Float64
"""
function _integrate(
    x::AbstractVector{<:Real}, 
    y::AbstractVector{<:Real},
    ::Trapezoidal
)

    # check size
    len_y = length(y)
    length(x) == len_y || throw(DimensionMismatch("x and y must have same length, got $(length(x)) and $len_y"))
    len_y >= 2 || throw(ArgumentError("at least 2 points required, got $len_y"))

    # get type
    T = promote_type(eltype(x), eltype(y), Float64)
    retval = zero(T)

    # get index offsets
    ix = firstindex(x) - 1
    iy = firstindex(y) - 1

    # cache first x and y value
    x1, y1 = T(x[ix + 1]), T(y[iy + 1]) 

    @inbounds @fastmath for i in 2:len_y
        
        # cache next x and y value
        x2, y2 = T(x[ix + i]), T(y[iy + i])

        # calc integral
        retval += (x2 - x1) * (y2 + y1)
        
        # transfer cached x and y value
        x1, y1 = x2, y2
    end

    return retval/2
end

"""
    _integrate(
        Δx::Real, 
        y::AbstractVector{<:Real},
        ::Trapezoidal
    )

Uses Trapezoidal rule, to compute numerical integral of y(x) with even spacing Δx.

# Arguments:
- `Δx::Real`: even spacing of y values
- `y::AbstractVector{<:Real}`: vector of corresponding y values
- `::Trapezoidal`: integration method

# Returns:
- scalar of at least Float64
"""
function _integrate(
    Δx::Real, 
    y::AbstractVector{<:Real},
    ::Trapezoidal
)   
    
    # check size
    len_y = length(y)
    len_y >= 2 || throw(ArgumentError("at least 2 points required, got $len_y"))

    # get type
    T = promote_type(eltype(y), typeof(Δx), Float64)

    # get index offset
    iy = firstindex(y) - 1

    @inbounds begin
    
        # cache first half y value
        retval = T(y[iy + 1]) / 2

        # inner values are full
        @fastmath for i in 2:(len_y-1)

            # calc integral
            retval += T(y[iy + i])
        end

        # add last half value
        retval += T(y[iy + len_y]) / 2
    end

    return retval * Δx
end

"""
    _cumul_integrate!(
        dest::AbstractVector{<:AbstractFloat},
        x::AbstractVector{<:Real}, 
        y::AbstractVector{<:Real},
        ::Trapezoidal   
    )

Uses Trapezoidal rule, to compute numerical cumulative integral of y(x) over x inplace in dest.

# Arguments:
- `dest::AbstractVector{<:AbstractFloat}`: destination vector for cumulative integral
- `x::AbstractVector{<:Real}`: vector or range of x values
- `y::AbstractVector{<:Real}`: vector of corresponding y values
- `::Trapezoidal`: integration method
"""
function _cumul_integrate!(
    dest::AbstractVector{<:AbstractFloat},
    x::AbstractVector{<:Real}, 
    y::AbstractVector{<:Real},
    ::Trapezoidal   
)

    # check size
    len_y = length(y)
    ((length(dest) == len_y) && (length(x) == len_y)) || throw(DimensionMismatch("dest, x and y must have same length, got $(length(dest)), $(length(x)) and $len_y"))
    len_y >= 2 || throw(ArgumentError("at least 2 points required, got $len_y"))

    # get type
    T = promote_type(eltype(x), eltype(y), Float64)
    accum = zero(T)
    
    # get index offsets
    id = firstindex(dest) - 1
    ix = firstindex(x) - 1
    iy = firstindex(y) - 1
    
    @inbounds begin
        # cache first x and y value
        x1, y1 = T(x[ix + 1]), T(y[iy + 1])

        # first value in dest is zero
        dest[id + 1] = zero(eltype(dest))

        @fastmath for i in 2:len_y
            
            # cache next value
            x2, y2 = T(x[ix + i]), T(y[iy + i])

            # calc next cumul integral value
            accum += (x2 - x1) * (y2 + y1) / 2

            # write accum to dest vector
            dest[id + i] = accum

            # transfer cached x and y value
            x1, y1 = x2, y2
        end
    end
    
    return dest
end

"""
    _cumul_integrate!(
        dest::AbstractVector{<:AbstractFloat},
        Δx::Real, 
        y::AbstractVector{<:Real},
        ::Trapezoidal   
    )

Uses Trapezoidal rule, to compute numerical cumulative integral of y(x) with even spacing Δx inplace in dest.

# Arguments:
- `dest::AbstractVector{<:AbstractFloat}`: destination vector for cumulative integral
- `Δx::Real`: vector or range of x values
- `y::AbstractVector{<:Real}`: vector of corresponding y values
- `::Trapezoidal`: integration method
"""
function _cumul_integrate!(
    dest::AbstractVector{<:AbstractFloat},
    Δx::Real, 
    y::AbstractVector{<:Real},
    ::Trapezoidal   
)

    # check size
    len_y = length(y)
    len_y >= 2 || throw(ArgumentError("at least 2 points required, got $len_y"))
    length(dest) == len_y || throw(DimensionMismatch("dest and y must have same length, got $(length(dest)) and $len_y"))

    # get type
    T = promote_type(typeof(Δx), eltype(y), Float64)
    accum = zero(T)

    # get index offsets
    id = firstindex(dest) - 1
    iy = firstindex(y) - 1
    
    # convert type
    Δx_T_half = T(Δx) / 2
    
    @inbounds begin
    
        # cache first y value
        y1 = T(y[iy + 1])

        # first value in dest is zero
        dest[id + 1] = zero(eltype(dest))
        
        @fastmath for i in 2:len_y
            
            # cache next value
            y2 = T(y[iy + i])

            # calc next cumul integral value
            accum += Δx_T_half * (y1 + y2)

            # write accum to dest vector
            dest[id + i] = accum

            # transfer cached y value
            y1 = y2
        end
    end

    return dest
end
