# src/simpson.jl


# compound simpson funcitons
############################

"""
    _simpson13(Δx, y, iy_start, N)

Uses Simpson rule, to compute numerical integral over N elements of y from index iy_start.
N must be an odd integer.

# Arguments:
- `Δx`: even spacing of y values
- `y`: vector of corresponding y values
- `iy_start`: start index
- `N`: number of elements for integration (odd)

# Returns:
- `retval`: integral value
- `y1`: cached value for _simpson38 function
"""
function _simpson13(Δx, y, iy_start, N)
    
    # get types
    T = promote_type(eltype(y), typeof(Δx), Float64)
    
    # cach sums
    s2 = zero(T)
    s4 = zero(T)
    
    @inbounds begin
        
        # simpson sums
        @fastmath for i in (iy_start + 1):2:(iy_start + N - 4)
            s4 += y[i]
            s2 += y[i+1]
        end
        
        # add missing element
        s4 += T(y[iy_start + N - 2])

        # cache and return last element
        y1 = T(y[iy_start + N - 1])

        retval =  (T(y[iy_start]) + 4 * s4 + 2 * s2 + y1) * (T(Δx) / 3)
    end

    return retval, y1
end


"""
    _simpson(x, y, ix_start, iy_start, N)

Uses extended Simpson rule, to compute numerical integral over N elements of y(x) from index ix_start and iy_start.
N must be an odd integer.

# Arguments
- `x`: vector of x values
- `y`: vector of corresponding y values
- `ix_start, iy_start`: start indices
- `N`: number of elements for integration (odd)

# Returns
- `retval`: integral value
- `x1, x2, y1, y2`: cached values for _simpson_correction function
"""
function _simpson(x, y, ix_start, iy_start, N)

    # get types
    T = promote_type(eltype(y), eltype(x), Float64)
    retval = zero(T)

    @inbounds begin
        
        # cache first values
        x_first = T(x[ix_start])
        y_first = T(y[iy_start])

        # init x_mid and y_mid to be available after the for loop
        x_mid = y_mid = zero(T)

        for k in 1:div(N - 1, 2)

            # cache other interval values
            x_mid = T(x[2k + ix_start - 1])
            y_mid = T(y[2k + iy_start - 1])

            x_last = T(x[2k + ix_start])
            y_last = T(y[2k + iy_start])
            
            retval += _simpson_interval(x_first, x_mid, x_last, y_first, y_mid, y_last)

            # transfer cached variables
            x_first, y_first = x_last, y_last
        end
    end

    # return last cached elements for correction calculation
    x1, x2, y1, y2 = x_mid, x_first, y_mid, y_first

    return retval, x1, x2, y1, y2
end



# simpson single
################

"""
    _simpson13_interval(Δx, y1, y2, y3)

Uses Simpson 1/3 rule, to compute numerical integral over y1, y2, y3.

# Arguments:
- `Δx`: even spacing of y values
- `y1, y2, y3`: y values for integration

# Returns:
- `retval`: integral value
"""
@inline function _simpson13_interval(Δx, y1, y2, y3)
    
    retval = (Δx / 3) * (y1 + 4 * y2 + y3)
    
    return retval
end

"""
    _simpson38_interval(Δx, y1, y2, y3, y4)

Uses Simpson 3/8 rule, to compute numerical integral over y1, y2, y3, y4.

# Arguments:
- `Δx`: even spacing of y values
- `y1, y2, y3, y4`: y values for integration

# Returns:
- `retval`: integral value
"""
@inline function _simpson38_interval(Δx, y1, y2, y3, y4)
    
    retval = (3 * Δx / 8) * (y1 + 3 * y2 + 3 * y3 + y4) 
   
    return retval
end

"""
    _simpson_interval(x1, x2, x3, y1, y2, y3)

Uses extended simpson rule to compute numercal integral over y1, y2, y3 with the corresponding x values x1, x2, x3

# Arguments:
- `x1, x2, x3`: x values
- `y1, y2, y3`: corresponding y values

# Returns:
- `retval`: integral value
"""
@inline function _simpson_interval(x1, x2, x3, y1, y2, y3)
    
    @fastmath begin

        h1 = x2 - x1
        h2 = x3 - x2
        H = h1 + h2

        α = 2 - h2 / h1
        β = H^2 / (h1 * h2)
        η = 2 - h1 / h2
        
        retval = H / 6 * (α * y1 + β * y2 + η * y3)
    end

    return retval
end

"""
    _simpson_correction(x1, x2, x3, y1, y2, y3)

Correction for `_simpson`, if length of y is even.

- `x1, x2, x3`: x_values, where x3 is the missing x value with even index
- `y1, y2, y3`: y values, where y3 is the missing y value with even index
"""
@inline function _simpson_correction(x1, x2, x3, y1, y2, y3)
    
    @fastmath begin
        
        h1 = x2 - x1
        h2 = x3 - x2
        H = h1 + h2

        α = (2 * h2^2 + 3 * h2 * h1) / (6 * H)
        β = (h2^2 + 3 * h2 * h1) / (6 * h1)
        η = h2^3 / (6 * h1 * H)

        retval = α * y3 + β * y2 - η * y1
    end

    return retval
end



# simpson integration
#####################

"""
    _integrate(
        x::AbstractVector{<:Real}, 
        y::AbstractVector{<:Real},
        ::Simpson
    )

Uses Simpson rule, to compute numerical integral of y(x) over x.

# Arguments:
- `x::AbstractVector{<:Real}`: vector or range of x values
- `y::AbstractVector{<:Real}`: vector of corresponding y values
- `::Simpson`: integration method

# Returns:
- scalar of same type as input
"""
function _integrate(
    x::AbstractVector{<:Real}, 
    y::AbstractVector{<:Real},
    ::Simpson
)

    # check size
    len_y = length(y)
    length(x) == len_y || throw(DimensionMismatch("x and y must have same length, got $(length(x)) and $len_y"))
    len_y >= 3 || throw(ArgumentError("at least 3 points required, got $len_y"))

    # get types
    T = promote_type(eltype(x), eltype(y), Float64)

    # get index offsets
    ix = firstindex(x) - 1
    iy = firstindex(y) - 1

    if isodd(len_y)
        retval, _, _, _, _ = _simpson(x, y, ix + 1, iy + 1, len_y)
    else
        retval, x1, x2, y1, y2 = _simpson(x, y, ix + 1, iy + 1, len_y - 1)

        # correction
        @inbounds x3, y3 = T(x[ix + len_y]), T(y[iy + len_y])

        retval += _simpson_correction(x1, x2, x3, y1, y2, y3)
    end

    return retval
end


"""
    _integrate(
        Δx::Real, 
        y::AbstractVector{<:Real},
        ::Simpson
    )

Uses Simpson rule, to compute numerical integral of y(x) with even spacing Δx.

# Arguments:
- `Δx::Real`: even spacing of y values
- `y::AbstractVector{<:Real}`: vector of corresponding y values
- `::Simpson`: integration method

# Returns:
- scalar of at least Float64
"""
function _integrate(
    Δx::Real, 
    y::AbstractVector{<:Real},
    ::Simpson
)   
    
    # check size
    len_y = length(y)
    len_y >= 3 || throw(ArgumentError("at least 3 points required, got $len_y"))
    
    # get types
    T = promote_type(typeof(Δx), eltype(y), Float64)

    # get index offset
    iy = firstindex(y) - 1

    if isodd(len_y)
        retval, _ = _simpson13(Δx, y, iy + 1, len_y)
    elseif (len_y >= 6) # and even
        retval, y1 = _simpson13(Δx, y, iy + 1, len_y - 3)

        # get elements
        y2 = T(y[iy + len_y - 2])
        y3 = T(y[iy + len_y - 1])
        y4 = T(y[iy + len_y])

        # add last part with simpson 3/8
        retval += _simpson38_interval(Δx, y1, y2, y3, y4)
    else # (len_y == 4)
        
        # get elements
        y1 = T(y[iy + 1])
        y2 = T(y[iy + 2])
        y3 = T(y[iy + 3])
        y4 = T(y[iy + 4])
        
        # calc integral
        retval = _simpson38_interval(Δx, y1, y2, y3, y4)
    end

    return retval
end

"""
    _cumul_integrate!(
        dest::AbstractVector{<:AbstractFloat},
        x::AbstractVector{<:Real}, 
        y::AbstractVector{<:Real},
        ::Simpson   
    )

Uses Simpson rule, to compute numerical cumulative integral of y(x) over x inplace in dest.

# Arguments:
- `dest::AbstractVector{<:AbstractFloat`: destination vector for cumulative integral
- `x::AbstractVector{<:Real}`: vector or range of x values
- `y::AbstractVector{<:Real}`: vector of corresponding y values
- `::Simpson`: integration method
"""
function _cumul_integrate!(
    dest::AbstractVector{<:AbstractFloat},
    x::AbstractVector{<:Real}, 
    y::AbstractVector{<:Real},
    ::Simpson   
)

    # check size
    len_y = length(y)
    (length(x) == len_y) && (length(dest) == len_y) || throw(DimensionMismatch("x, y and dest must have same lenth, got $(length(x)), $len_y and $(length(dest))"))
    len_y >= 3 || throw(ArgumentError("at least 3 points required, got $len_y"))

    # get type
    T = promote_type(eltype(dest), eltype(x), eltype(y), Float64)

    # get index offsets
    id = firstindex(dest) - 1
    iy = firstindex(y) - 1
    ix = firstindex(x) - 1

    @inbounds begin

        # cache first elements
        x1, y1 = T(x[ix + 1]), T(y[iy + 1])
        x2, y2 = T(x[ix + 2]), T(y[iy + 2])
        x3, y3 = T(x[ix + 3]), T(y[iy + 3])

        # calc first elements
        dest[id + 1] = zero(eltype(dest))
        accum = -_simpson_correction(x3, x2, x1, y3, y2, y1)
        dest[id + 2] = accum
        dest[id + 3] = _simpson_interval(x1, x2, x3, y1, y2, y3)

        # transfer cached values
        x1, y1 = x2, y2
        x2, y2 = x3, y3

        # calculate further elements in loop
        @fastmath for i in 4:2:(len_y - 1)

            # cache next elements
            x3, y3 = T(x[ix + i]), T(y[iy + i])
            x4, y4 = T(x[ix + i + 1]), T(y[iy + i + 1])

            # accum next integral value
            accum += _simpson_interval(x1, x2, x3, y1, y2, y3)

            correction = _simpson_correction(x2, x3, x4, y2, y3, y4)

            # write results to dest
            dest[id + i] = accum
            dest[id + i + 1] = accum + correction

            # transfer cached variables
            x1, y1 = x3, y3
            x2, y2 = x4, y4
        end

        # add last element, if length is even
        if iseven(len_y)
            
            # cache elements
            x3, y3 = T(x[ix + len_y ]), T(y[iy + len_y])

            accum += _simpson_interval(x1, x2, x3, y1, y2, y3)
            dest[id + len_y] = accum
        end
    end

    return dest
end

"""
    _cumul_integrate!(
        dest::AbstractVector{<:AbstractFloat},
        Δx::Real, 
        y::AbstractVector{<:Real},
        ::Simpson   
    )

Uses Simpson rule, to compute numerical cumulative integral of y(x) with even spacing Δx inplace in dest.

# Arguments:
- `dest::AbstractVector{<:AbstractFloat}`: destination vector for cumulative integral
- `Δx::Real`: even spacing of y values
- `y::AbstractVector{<:Real}`: vector of corresponding y values
- `::Simpson`: integration method
"""
function _cumul_integrate!(
    dest::AbstractVector{<:AbstractFloat},
    Δx::Real, 
    y::AbstractVector{<:Real},
    ::Simpson   
)

    # check size
    len_y = length(y)
    length(dest) == len_y || throw(DimensionMismatch("dest and y must have same length, got $(length(dest)) and $len_y"))
    len_y >= 3 || throw(ArgumentError("at least 3 points required, got $len_y"))

    # get type
    T = promote_type(typeof(Δx), eltype(y), Float64)
    Δx_T = T(Δx)

    # get index offsets
    id = firstindex(dest) - 1
    iy = firstindex(y) - 1
    
    @inbounds begin

        # calc first value
        dest[id + 1] = zero(eltype(dest))

        # cache first three y values
        y1 = T(y[iy + 1])
        y2 = T(y[iy + 2])
        y3 = T(y[iy + 3])

        # if length is 3, calculate with quadratic fit 
        if len_y == 3
            dest[id + 2] = @fastmath Δx_T * ((T(5)/12) * y1 + (T(2)/3) * y2 - (T(1)/12) * y3)
            dest[id + 3] = _simpson13_interval(Δx_T, y1, y2, y3)
        else
            # cache one more y value
            y4 = T(y[iy + 4])
            
            s38 = _simpson38_interval(Δx_T, y1, y2, y3, y4)
            dest[id + 2] = s38 - _simpson13_interval(Δx_T, y2, y3, y4)
            accum = _simpson13_interval(Δx_T, y1, y2, y3)
            dest[id + 3] = accum
            dest[id + 4] = s38

            # transfer cached values
            y1, y2 = y3, y4

            # sum up following values
            for i in 5:2:(len_y - 1)

                # cache values
                y3 = T(y[iy + i])
                y4 = T(y[iy + i + 1])
                
                # calc every second value with simpson 3/8
                dest[id + i + 1] = _simpson38_interval(Δx_T, y1, y2, y3, y4) + accum
                
                # and accumulate with simpson 1/3
                accum += _simpson13_interval(Δx_T, y1, y2, y3)
                dest[id + i] = accum

                # transfer cached values
                y1, y2 = y3, y4
            end

            # if len_y is odd, last value is missing
            if isodd(len_y)
                dest[id + len_y] = accum + _simpson13_interval(Δx_T, y1, y2, T(y[iy + len_y]))
            end
        end
    end

    return dest
end
