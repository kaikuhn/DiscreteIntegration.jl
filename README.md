# DiscreteIntegration

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kaikuhn.github.io/DiscreteIntegration.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kaikuhn.github.io/DiscreteIntegration.jl/dev/)
[![Build Status](https://github.com/kaikuhn/DiscreteIntegration.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kaikuhn/DiscreteIntegration.jl/actions/workflows/CI.yml?query=branch%3Amain)

Numerical integration of discrete data on uniform (`Δx`) or non-uniform (`x`) grids, using the trapezoidal rule or Simpson's rule. Definite and cumulative integrals are supported.

## Features

- Definite integral (`integrate`) and cumulative integral (`cumul_integrate`, `cumul_integrate!`)
- Uniform grids (spacing `Δx`) and non-uniform grids (vector or range `x`)
- Methods: `Trapezoidal()` (default) and `Simpson()`
- Every element of `x` and `y` is loaded only once, so lazy arrays (e.g. from [LazyArrays.jl](https://github.com/JuliaArrays/LazyArrays.jl)) are supported
- Computation in at least `Float64`, also for `Float32` and integer input
- Works with non-1-based indices (e.g. OffsetArrays)

## Relation to NumericalIntegration.jl

DiscreteIntegration.jl and [NumericalIntegration.jl](https://github.com/JuliaMath/NumericalIntegration.jl) both integrate sampled data, so their basic use case overlaps. Main differences:

- **Single pass:** `x` and `y` are each accessed exactly once per element, so lazy arrays (e.g. from [LazyArrays.jl](https://github.com/JuliaArrays/LazyArrays.jl)) are evaluated only once and no temporary arrays are needed.
- **Simpson's rule on uniform grids:** the composite Simpson 1/3 rule is used, with the 3/8 rule for the last interval if the number of points is even.
- **Simpson's rule on non-uniform grids:** a composite Simpson's rule for irregularly spaced data, with a correction for the last interval if the number of points is even.
- **In-place cumulative integrals:** `cumul_integrate!` writes directly into a provided `dest` vector, which may be identical to `x` or `y`. Computation uses at least `Float64`, the element type of `dest` is chosen by the caller.
- **No low-order start value:** all elements of the cumulative integral are computed with Simpson's rule.
- **Uniform grids without `x`:** pass only the spacing `Δx`.
- **Non-1-based indices** (e.g. OffsetArrays) are supported.

Use DiscreteIntegration.jl if you work with LazyArrays.jl, need the classical 1/3 Simpson's rule or a Simpson's rule on non-uniform grids. Use NumericalIntegration.jl for maximum speed (using `TrapezoidalFast()` or `SimpsonEvenFast()`), or if you need the methods it offers, that are not available here (for example Romberg).

Currently available are the Trapezoidal rule and Simpson's rule.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/kaikuhn/DiscreteIntegration.jl")
```

(After registration: `Pkg.add("DiscreteIntegration")`.)

## Usage

```julia
using DiscreteIntegration

x = range(0, 1; length = 101)
y = x .^ 2

# definite integral, trapezoidal rule (default)
integrate(x, y)

# Simpson's rule
integrate(x, y, Simpson())

# uniform grid: pass the spacing instead of x
integrate(step(x), y, Simpson())

# cumulative integral (first element is 0)
cumul_integrate(x, y, Simpson())

# cumulative integral inplace
dest = similar(y)
cumul_integrate!(dest, x, y, Simpson())
```

## Methods

| Method | Exact for polynomials up to degree |
|---|---|
| `Trapezoidal()` | 1 |
| `Simpson()`, non-uniform grid `x` | 2 |
| `Simpson()`, uniform grid `Δx` | 3 |

`Simpson()` requires at least 3 points. For an even number of points, the last interval is handled separately (non-uniform grid: Cartwright correction, uniform grid: Simpson 3/8 rule).

## Notes

- `x` must be strictly monotonic for `Simpson()` (duplicate x values lead to a division by zero).
- Vectors must have contiguous 1D indices.
- Internally `@fastmath` is used, results may differ from strict IEEE arithmetic in the last digits, and `NaN`/`Inf` in the input give undefined results.
- For `cumul_integrate!`, the caller chooses the element type of `dest`. Computation uses at least `Float64` and is rounded on write.

## Documentation

See the [documentation](https://kaikuhn.github.io/DiscreteIntegration.jl/dev/) for the full API.

## License

See [LICENSE](LICENSE).