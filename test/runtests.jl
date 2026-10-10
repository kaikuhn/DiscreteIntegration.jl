using DiscreteIntegration
using Test
using LazyArrays
using OffsetArrays

function alloc(f, args...)
    f(args...)
    return @allocated f(args...)
end

@testset "DiscreteIntegration.jl" begin

    # define test functions
    linear          = x -> x
    linear_int      = x -> 0.5 * x^2

    quadratic       = x -> x^2
    quadratic_int   = x -> 1/3 * x^3

    # test range
    x_3 = range(1, 2; length=3)
    x_4 = range(1, 2; length=4)
    x_odd = range(1, 2; length=101)
    x_even = range(1, 2; length=100)

    for x in (x_3, x_4, x_even, x_odd)
    
        # Trapezoidal
        #############

        for (f, F) in [(linear, linear_int)]
            
            # calc y_array
            y = f.(x)

            @test isapprox(integrate(x, y) - (F(x[end]) - F(x[begin])), 0, atol=1e-12)
            @test isapprox(integrate(step(x), y) - (F(x[end]) - F(x[begin])), 0, atol=1e-12)

            @test isapprox(cumul_integrate(x, y) .- (F.(x) .- F(x[begin])), zeros(length(x)), atol= 1e-12)
            @test isapprox(cumul_integrate(step(x), y) .- (F.(x) .- F(x[begin])), zeros(length(x)), atol= 1e-12)

        end

        # Simpson
        #########
        
        for (f, F) in [(linear, linear_int), (quadratic, quadratic_int)]

            
            # calc y_array
            y = f.(x)

            @test isapprox(integrate(x, y, Simpson()) - (F(x[end]) - F(x[begin])), 0, atol=1e-12)
            @test isapprox(integrate(step(x), y, Simpson()) - (F(x[end]) - F(x[begin])), 0,  atol=1e-12)

            @test isapprox(cumul_integrate(x, y, Simpson()) .- (F.(x) .- F(x[begin])), zeros(length(x)), atol= 1e-8)
            @test isapprox(cumul_integrate(step(x), y, Simpson()) .- (F.(x) .- F(x[begin])), zeros(length(x)), atol= 1e-12)
        end
    end
end


@testset "LazyArrays" begin
    
    # define test functions
    f(x) = x
    F(x) = 0.5 * x * x

    # test range
    x = range(1, 2; length=101)
    y = BroadcastArray(@~ f.(x))

    # Trapezoidal
    #############
    @test isapprox(integrate(x, y) - (F(x[end]) - F(x[begin])), 0, atol=1e-12)
    @test isapprox(integrate(step(x), y) - (F(x[end]) - F(x[begin])), 0, atol=1e-12)

    @test alloc(integrate, x, y) == 0
    @test alloc(integrate, step(x), y) == 0

    @test isapprox(cumul_integrate(x, y) .- (F.(x) .- F(x[begin])), zeros(length(x)), atol= 1e-12)
    @test isapprox(cumul_integrate(step(x), y) - (F.(x) .- F(x[begin])), zeros(length(x)), atol= 1e-12)

    # Simpson
    #########
    @test isapprox(integrate(x, y, Simpson()) - (F(x[end]) - F(x[begin])), 0, atol=1e-12)
    @test isapprox(integrate(step(x), y, Simpson()) - (F(x[end]) - F(x[begin])), 0, atol=1e-12)

    @test alloc(integrate, x, y, Simpson()) == 0
    @test alloc(integrate, step(x), y, Simpson()) == 0

    @test isapprox(cumul_integrate(x, y, Simpson()) .- (F.(x) .- F(x[begin])), zeros(length(x)), atol= 1e-12)
    @test isapprox(cumul_integrate(step(x), y, Simpson()) .- (F.(x) .- F(x[begin])), zeros(length(x)), atol= 1e-12)
end


@testset "OffsetArrays" begin
    
    # test functions
    f(x) = x
    F(x) = 0.5 * x * x

    # define test functions
    x = range(1, 2; length=101)

    # x and y as OffsetArrays
    y_off = OffsetArray(f.(x), -50:50)
    x_off = OffsetArray(collect(x), 10:110)

    for m in (Trapezoidal(), Simpson())
        # nur y verschoben
        @test integrate(x, y_off, m) ≈ 1.5 atol=1e-12
        @test integrate(step(x), y_off, m) ≈ 1.5 atol=1e-12

        # x und y verschoben
        @test integrate(x_off, y_off, m) ≈ 1.5 atol=1e-12

        # kumulativ (Rückgabe ist ein normaler Vector)
        @test isapprox(cumul_integrate(x, y_off, m), F.(x) .- F(x[1]); atol=1e-12)
        @test isapprox(cumul_integrate(x_off, y_off, m), F.(x) .- F(x[1]); atol=1e-12)
        @test isapprox(cumul_integrate(step(x), y_off, m), F.(x) .- F(x[1]); atol=1e-12)

        # inplace mit verschobenem dest
        dest = OffsetArray(zeros(101), 3:103)
        cumul_integrate!(dest, x_off, y_off, m)
        @test isapprox(parent(dest), F.(x) .- F(x[1]); atol=1e-12)

        # Allokationen
        @test alloc(integrate, x_off, y_off, m) == 0
        @test alloc(integrate, step(x), y_off, m) == 0
    end
end
