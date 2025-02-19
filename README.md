# DiffusionX.jl
English | [简体中文](README-zh.md)

> [!NOTE]
> Development is in progress. DiffusionX.jl is a Julia-based library for stochastic process simulation. Rust/Python versions are also under development, see [DiffusionX](https://github.com/tangxiangong/diffusionx).

## Usage
```julia
using DiffusionX

# Simulate standard Brownian motion path
B = Bm()              
Bₜ = B(10)            
t, x = simulate(Bₜ; τ=0.01)

# Calculate moments
T = collect(10:10:100)
@. 𝔼(B(T))
@. 𝔼((B(T))^2)

# Simulate first passage time
fpt = FPT((-1, 1), B)
simulate(fpt) 
𝔼(fpt)

# Calculate time-averaged mean square displacement
𝔼(δ̄²(B, T=100, Δ=0.1))
```

## Extensibility

DiffusionX.jl leverages Julia's type system and multiple dispatch to make extending new stochastic processes intuitive and straightforward:

1. **Custom Processes**: Define new stochastic processes by inheriting from the `StochasticProcess` abstract type:
```julia
struct MyProcess <: StochasticProcess
    # Define process parameters
    x₀::Float64  # e.g., initial position
end
```

2. **Simulation Implementation**: Implement the `simulate` method for your new process:
```julia
function simulate(p::MyProcess, T::Union{Int,Float64}; τ::Float64=0.01)
    # Implement simulation logic
    # Return time series t and corresponding trajectory x
    t, x
end
```

3. **Functional Calculations**: Utilize built-in functional calculation tools:
```julia
# First passage time
fpt = FPT((-1, 1), process)  # Time to reach boundaries of interval [-1,1]
simulate(fpt)
# Occupation time
ot = OccupationTime(T, (a, b), process)  # Time spent in interval [a,b]
simulate(ot)
```

4. **Statistical Analysis**: All processes automatically support moment calculations and expectation analysis:
```julia
# Calculate expectations
𝔼(process(T))  # Expectation at time T
𝔼(process(T)^2)  # Second moment at time T
𝔼(fpt)  # Expected first passage time
𝔼(ot)  # Expected occupation time
𝔼(δ̄²(process, 100, 0.1))  # Time-averaged mean square displacement
```

All extended functionalities seamlessly integrate into the existing framework, automatically supporting trajectory simulation, functional calculations, and statistical analysis capabilities.

## License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).
