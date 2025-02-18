@doc raw"""
    `StochasticProcess` 
    
    The abstract type of stochastic processes.
"""
abstract type StochasticProcess end

@doc raw"""
    `Trajectory(sp::StochasticProcess, T::Float64)` 
    
    The type of stochastic process trajectories.

    # Arguments
    - `sp`: The stochastic process.
    - `T`: The length of the trajectory.

    # Examples
    ```julia
    using DiffusionX
    sp = Bm()
    traj = sp(10)
    ```
"""
struct Trajectory{SP<:StochasticProcess}
    sp::SP
    T::Float64
    function Trajectory(sp::SP, T) where {SP<:StochasticProcess}
        @assert T > 0
        if T isa Integer
            T = Float64(T)
        end
        new{SP}(sp, T)
    end
end

function Base.show(io::IO, traj::Trajectory)
    print(io, "The trajectory of $(traj.sp) with length $(traj.T)")
end

@doc raw"""
    `(sp::StochasticProcess)(T)`

    The callable constructor of `Trajectory`.

    # Arguments
    - `T`: The length of the trajectory.
"""
function (sp::StochasticProcess)(T)
    Trajectory(sp, T)
end

simulate_method(sp::SP) where {SP<:StochasticProcess} = nothing

@doc raw"""
    `simulate(traj::Trajectory, τ=1e-2)`

    The function to simulate the trajectory.
    
    **Need the method `simulate_method` for the stochastic process.**

    # Arguments
    - `traj`: The trajectory.
    - `τ`: The time step.
"""
function simulate(traj::Trajectory, τ=1e-2)
    method = simulate_method(traj.sp)
    isnothing(method) && throw(ArgumentError("No simulation method for $(traj.sp)"))
    method(traj.T, args(traj.sp)...; τ=τ)
end

@doc raw"""
   `PowerTrajectory(traj::Trajectory, order::Int)`

   The type of power trajectories.

   # Arguments
   - `traj`: The trajectory.
   - `order`: The order of the power trajectory.
"""
struct PowerTrajectory
    traj::Trajectory
    order::Int
end

@doc raw"""
    `^(traj::Trajectory, order::Int)`

    The callable constructor of `PowerTrajectory`.

    # Arguments
    - `traj`: The trajectory.
    - `order`: The order of the power trajectory.
"""
function Base.:^(traj::Trajectory, order::Int)
    PowerTrajectory(traj, order)
end

@doc raw"""
    `moment(traj::Trajectory, N::Int; τ=0.01, order::Int=1)`

    The function to calculate the moment of the trajectory.

    # Arguments
    - `traj`: The trajectory.   
    - `N`: The number of samples.
    - `τ`: The time step.
    - `order`: The order of the moment.

    # Examples
    ```julia
    using DiffusionX
    sp = Bm()
    traj = sp(10)
    moment(traj, 1000; τ=0.01, order=2)
    ```
"""
function moment(traj::Trajectory; N::Int=1000, τ=0.01, order::Int=1)
    moment = zeros(nthreads())
    method = simulate_method(traj.sp)
    isnothing(method) && throw(ArgumentError("No simulation method for $(traj.sp)"))
    T = traj.T
    vargs = args(traj.sp)
    @threads for _ in 1:N
        __, x = method(T, vargs...; τ=τ)
        @inbounds moment[threadid()] += x[end]^order
    end
    sum(moment) / N
end

moment(ptraj::PowerTrajectory; N::Int=1_000, τ=0.01) = moment(ptraj.traj; N=N, τ=τ, order=ptraj.order)

𝔼(traj::Trajectory; N::Int=1_000, τ=0.01) = moment(traj; N=N, τ=τ)
𝔼(ptraj::PowerTrajectory; N::Int=1_000, τ=0.01) = moment(ptraj; N=N, τ=τ)

mean(traj::Trajectory; N::Int=1_000, τ=0.01) = moment(traj; N=N, τ=τ)
msd(traj::Trajectory; N::Int=1_000, τ=0.01) = moment(traj; N=N, τ=τ, order=2) - moment(traj; N=N, τ=τ)^2

@doc raw"""
    `Functional{SP<:StochasticProcess}`

    The abstract type of functionals.
"""
abstract type Functional{SP<:StochasticProcess} end

@doc raw"""
    `^(functional::Functional, order::Int)`

    The callable constructor of `FunctionalPower`.
"""
struct FunctionalPower{F<:Functional}
    functional::F
    order::Int
end

moment(fp::FunctionalPower; N::Int=1_000, τ=0.01) = moment(fp.functional; N=N, τ=τ, order=fp.order)

𝔼(fp::FunctionalPower; N::Int=1_000, τ=0.01) = moment(fp; N=N, τ=τ)


@doc raw"""
    `FPT(domain::NTuple{2,Float64}, sp::StochasticProcess)`

    The type of first passage time functionals.
"""
struct FPT{SP<:StochasticProcess} <: Functional{SP}
    domain::NTuple{2,Float64}
    sp::SP
    function FPT(domain, sp::SP) where {SP<:StochasticProcess}
        domain[1] >= domain[2] && throw(ArgumentError("domain[1] must be less than domain[2]"))
        if domain[1] isa Integer || domain[2] isa Integer
            domain = (Float64(domain[1]), Float64(domain[2]))
        end
        new{SP}(domain, sp)
    end
end

function Base.show(io::IO, f::FPT)
    print(io, "The first passage time of $(f.sp) on the interval $(f.domain)")
end

@doc raw"""
    `^(functional::Functional, order::Int)`

    The callable constructor of `FunctionalPower`.
"""
function Base.:^(functional::F, order::Int) where {F<:Functional}
    FunctionalPower(functional, order)
end

function firstpassagetime(domain, method, vargs...; τ=1e-2)
    a, b = domain
    counter = 1
    while true
        t, path = method(2^counter, vargs...; τ=τ)
        index = findfirst(path) do x
            x <= a || x >= b
        end
        !isnothing(index) && return t[index]
        counter += 1
        counter == 60 && error("NOT PASS!")
    end
end

@doc raw"""
    `simulate(f::FPT, τ=1e-2)`

    The function to simulate the first passage time.
"""
function simulate(f::FPT; τ=1e-2)
    method = simulate_method(f.sp)
    isnothing(method) && throw(ArgumentError("No simulation method for $(f.sp)"))
    firstpassagetime(f.domain, method, args(f.sp)...; τ=τ)
end

function simulate_uncheck(f::FPT; τ=1e-2)
    method = simulate_method(f.sp)
    firstpassagetime(f.domain, method, args(f.sp)...; τ=τ)
end

@doc raw"""
    `moment(functional::F, N::Int; τ=1e-2, order::Int=1) where {F<:Functional}` 

    The function to calculate the moment of the functional.

    # Arguments
    - `functional`: The functional.
    - `N`: The number of samples.
    - `τ`: The time step.
    - `order`: The order of the moment.
"""
function moment(functional::F; N::Int=1000, τ=1e-2, order::Int=1) where {F<:Functional}
    moment = zeros(nthreads())
    method = simulate_method(functional.sp)
    isnothing(method) && throw(ArgumentError("No simulation method for $(functional.sp)"))
    @threads for _ in 1:N
        x = simulate_uncheck(functional; τ=τ)
        @inbounds moment[threadid()] += x^order
    end
    sum(moment) / N
end

@doc raw"""
    `OccupationTime(T::Float64, domain::NTuple{2, Float64}, sp::StochasticProcess)`

    The type of occupation time functionals.
"""
struct OccupationTime{SP<:StochasticProcess} <: Functional{SP}
    T::Float64
    domain::NTuple{2,Float64}
    sp::SP
    function OccupationTime(T, domain, sp::SP) where {SP<:StochasticProcess}
        domain[1] >= domain[2] && throw(ArgumentError("$(domain[1]) must be less than $(domain[2])"))
        T <= 0 && throw(ArgumentError("T must be positive"))
        if T isa Integer
            T = Float64(T)
        end
        if domain[1] isa Integer || domain[2] isa Integer
            domain = (Float64(domain[1]), Float64(domain[2]))
        end
        new{SP}(T, domain, sp)
    end
end

function Base.show(io::IO, ot::OccupationTime)
    print(io, "The occupation time of $(ot.sp) on the interval $(ot.domain) in the time interval [0, $(ot.T)]")
end

function occupationtime(domain, method, T, vargs...; τ=1e-2)
    a, b = domain
    t, x = method(T, vargs...; τ=τ)
    indices = findall(x -> a <= x <= b, x)
    isnothing(indices) && return 0
    length(indices) == length(x) && return T
    temp = diff(indices)
    isempty(temp) && return 0
    jump = findall(x -> x != 1, temp)
    isnothing(jump) && return t[indices[end]] - t[indices[begin]]
    isempty(jump) && return t[indices[end]] - t[indices[begin]]
    start_index = 1
    end_index = jump[1]
    dural = t[indices[end_index]] - t[indices[start_index]]
    @inbounds for k in firstindex(jump)+1:lastindex(jump)
        start_index = end_index + 2
        end_index = jump[k]
        start_index < end_index && (dural += t[indices[end_index]] - t[indices[start_index]])
    end
    dural
end

@doc raw"""
    `simulate(oc::OccupationTime; τ=1e-2)`

    The function to simulate the occupation time.
"""
function simulate(oc::OccupationTime; τ=1e-2)
    method = simulate_method(oc.sp)
    isnothing(method) && throw(ArgumentError("No simulation method for $(oc.sp)"))
    occupationtime(oc.domain, method, oc.T, args(oc.sp)...; τ=τ)
end

function simulate_uncheck(oc::OccupationTime; τ=1e-2)
    method = simulate_method(oc.sp)
    occupationtime(oc.domain, method, oc.T, args(oc.sp)...; τ=τ)
end

𝔼(functional::F; τ=1e-2, N::Int=1000) where {F<:Functional} = moment(functional; N=N, τ=τ)
𝔼(fp::FunctionalPower; τ=1e-2, N::Int=1000) = moment(fp.functional; N=N, τ=τ, order=fp.order)

@doc raw"""
    `TimeAverage(sp::StochasticProcess, T::Real, Δ::Real)`

    The type of time average functionals.
"""
struct TimeAverage{SP<:StochasticProcess}
    sp::SP
    T::Float64
    Δ::Float64
    function TimeAverage(sp::SP, T, Δ) where {SP<:StochasticProcess}
        T <= 0 && throw(ArgumentError("T must be positive"))
        Δ <= 0 && throw(ArgumentError("Δ must be positive"))
        if T isa Integer
            T = Float64(T)
        end
        if Δ isa Integer
            Δ = Float64(Δ)
        end
        new{SP}(sp, T, Δ)
    end
end
δ̄² = TimeAverage

# x(t+Δ)*x(t)
struct TrajMultiplication{SP<:StochasticProcess}
    sp::SP
    T::Float64
    Δ::Float64
end

Base.:*(xt::Trajectory{SP}, xs::Trajectory{SP}) where {SP} = TrajMultiplication(xt.sp, xs.T, xt.T - xs.T)

function trajmulmean(tm::TrajMultiplication; N::Int=1000, τ=0.01)
    T = tm.T
    Δ = tm.Δ
    slag = round(Int, Δ / τ)
    means = zeros(nthreads())
    method = simulate_method(tm.sp)
    isnothing(method) && throw(ArgumentError("No simulation method for $(tm.sp)"))
    @threads for _ in 1:N
        _, x = method(T + Δ, args(tm.sp)...; τ=τ)
        @inbounds means[threadid()] += x[end] * x[end-slag]
    end
    sum(means) / N
end

𝔼(tm::TrajMultiplication; N::Int=1000, τ=0.01) = trajmulmean(tm; N=N, τ=τ)

import FastGaussQuadrature: gausslegendre

function get_weights_nodes(a, b, order)
    nodes_unit, weights_unit = gausslegendre(order)
    weights = @. (b - a) * weights_unit / 2
    nodes = @. (b - a) * nodes_unit / 2 + (b + a) / 2
    weights, nodes
end

"""
    TAMSD(TA::TimeAverage, order::Int, τ::Float64, N::Int)
    𝔼(TA::TimeAverage; τ::Float64=1e-2, N::Int=100_000, order::Int=10)

Calculate the time average mean square displacement (TAMSD).

# Arguments
- `TA` : The time average, constructed by the constructor `δ̄²(::StochasticProcess, T::Real, Δ::Real)`
- `τ` : The time step of the Euler method
- `N` : The number of particles in the Monte Carlo simulation
- `order` : The number of terms in the Gaussian-Legendre quadrature

# Usage
```julia
T = 100; Δ = 1
x::StochasticProcess = ....  # Define a stochastic process instance
𝔼(δ̄²(x, T, Δ))  # Calculate TAMSD
```
"""
function TAMSD(TA::TimeAverage; N::Int=1000, τ=0.01, order::Int=10)
    T, Δ = TA.T, TA.Δ
    x = TA.sp
    kwargs = (τ=τ, N=N)
    ω, t = get_weights_nodes(0, T - Δ, order)
    m = 0.0
    @inbounds @simd for k in eachindex(t)
        m += (𝔼((x(t[k] + Δ))^2; kwargs...) + 𝔼((x(t[k])^2); kwargs...) - 2𝔼(x(t[k] + Δ) * x(t[k]); kwargs...)) * ω[k]
    end
    m / (T - Δ)
end

𝔼(TA::TimeAverage; τ::Float64=1e-2, N::Int=1000, order::Int=10) = TAMSD(TA; τ=τ, N=N, order=order)