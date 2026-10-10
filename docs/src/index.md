```@meta
CurrentModule = DiscreteIntegration
```

# DiscreteIntegration.jl

Numerical integration of discrete data: sampled values `y` on a uniform grid (spacing `Δx`) or a non-uniform grid `x`.
Supports the trapezoidal rule and Simpson's rule, as definite and as cumulative integrals.
Elements are loaded only once, so lazy arrays are supported.

## Installation

```julia
using Pkg
Pkg.add("DiscreteIntegration")
```

## Quick start

```@example
using DiscreteIntegration

x = range(0, 1; length = 101)
y = x .^ 2

integrate(x, y)                  # trapezoidal rule (default)
```

```@example
using DiscreteIntegration

x = range(0, 1; length = 101)
y = x .^ 2

integrate(step(x), y, Simpson())   # uniform grid, Simpson
```

```@example
using DiscreteIntegration

x = range(0, 1; length = 5)
y = x .^ 2

cumul_integrate(x, y, Simpson())   # cumulative integral
```

## Methods

| Method | Exact for polynomials up to degree |
|---|---|
| `Trapezoidal()` | 1 |
| `Simpson()`, non-uniform grid `x` | 2 |
| `Simpson()`, uniform grid `Δx` | 3 |

## API

```@index
```

```@autodocs
Modules = [DiscreteIntegration]
Private = false
```
