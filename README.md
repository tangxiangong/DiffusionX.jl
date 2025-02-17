# DiffusionX.jl
English | [简体中文](README-zh.md)

> [!NOTE]
> Development is in progress. DiffusionX.jl is a Julia-based library for stochastic process simulation. Rust/Python versions are also under development, see [DiffusionX](https://github.com/tangxiangong/diffusionx).

## Usage Examples
```julia
using DiffusionX.Simulation

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
𝔼(δ̄²(x; T=100, Δ=0.1))
```

## License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).
