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

show(io::IO, traj::Trajectory) = print(io, "The trajectory of $(traj.sp) with length $(traj.T)")

@doc raw"""
    `(sp::StochasticProcess)(T)`

    The callable constructor of `Trajectory`.

    # Arguments
    - `T`: The length of the trajectory.
"""
(sp::StochasticProcess)(T) = Trajectory(sp, T)

simulate(traj::Trajectory, τ=1e-2) = traj.sp.method(traj.T, τ, args(traj.sp)...)

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
^(traj::Trajectory, order::Int) = PowerTrajectory(traj, order)

@doc raw"""
    `moments(traj::Trajectory, N::Int; τ=0.01, order::Int=1)`

    The function to calculate the moments of the trajectory.

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
    moments(traj, 1000; τ=0.01, order=2)
    ```
"""
function moments(traj::Trajectory; N::Int=1000, τ=0.01, order::Int=1)
    moment = zeros(nthreads())
    @threads for _ in 1:N
        __, x = simulate(traj, τ)
        @inbounds moment[threadid()] += x[end]^order
    end
    sum(moment) / N
end

@doc raw"""
    `Functional`

    The abstract type of functionals.
"""
abstract type Functional end

@doc raw"""
    `^(functional::Functional, order::Int)`

    The callable constructor of `FunctionalPower`.
"""
struct FunctionalPower{F<:Functional}
    functional::F
    order::Int
end

@doc raw"""
    `FPT(domain::NTuple{2,Float64}, sp::StochasticProcess)`

    The type of first passage time functionals.
"""
struct FPT <: Functional
    domain::NTuple{2,Float64}
    sp::StochasticProcess
    function FPT(domain, sp)
        domain[1] >= domain[2] && throw(ArgumentError("domain[1] must be less than domain[2]"))
        if domain[1] isa Integer || domain[2] isa Integer
            domain = (Float64(domain[1]), Float64(domain[2]))
        end
        new(domain, sp)
    end
end

show(io::IO, f::FPT) = print(io, "The first passage time of $(f.sp) on the interval $(f.domain)")

@doc raw"""
    `^(functional::Functional, order::Int)`

    The callable constructor of `FunctionalPower`.
"""
^(functional::Functional, order::Int) = FunctionalPower(functional, order)

@doc raw"""
    `simulate(f::FPT, τ=1e-2)`

    The function to simulate the first passage time.
"""
simulate(f::FPT, τ=1e-2) = firstpassagetime(f.domain, f.sp.method, τ, args(f.sp)...)


@doc raw"""
    `moments(functional::F, N::Int; τ=1e-2, order::Int=1) where {F<:Functional}` 

    The function to calculate the moments of the functional.

    # Arguments
    - `functional`: The functional.
    - `N`: The number of samples.
    - `τ`: The time step.
    - `order`: The order of the moment.
"""
function moments(functional::F, N::Int; τ=1e-2, order::Int=1) where {F<:Functional}
    moment = zeros(nthreads())
    @threads for _ in 1:N
        x = simulate(functional, τ)
        @inbounds moment[threadid()] += x^order
    end
    sum(moment) / N
end

@doc raw"""
    `OccupationTime(T::Float64, domain::NTuple{2, Float64}, sp::StochasticProcess)`

    The type of occupation time functionals.
"""
struct OccupationTime <: Functional
    T::Float64
    domain::NTuple{2,Float64}
    sp::StochasticProcess
    function OccupationTime(T, domain, sp)
        domain[1] >= domain[2] && throw(ArgumentError("$(domain[1]) must be less than $(domain[2])"))
        T <= 0 && throw(ArgumentError("T must be positive"))
        if T isa Integer
            T = Float64(T)
        end
        if domain[1] isa Integer || domain[2] isa Integer
            domain = (Float64(domain[1]), Float64(domain[2]))
        end
        new(T, domain, sp)
    end
end

show(io::IO, ot::OccupationTime) = print(io, "$(ot.sp) 在 [0, $(ot.T)] 内逗留在 $(ot.domain) 的时间")

@doc raw"""
    `simulate(oc::OccupationTime, τ=1e-2)`

    The function to simulate the occupation time.
"""
simulate(oc::OccupationTime, τ=1e-2) = occupationtime(oc.domain, oc.sp.method, oc.T, τ, args(oc.sp)...)



𝔼(traj::Trajectory; N::Int=100_000, τ=0.01) = moments(traj, N; τ=τ)
𝔼(ptraj::PowerTrajectory; N::Int=100_000, τ=0.01) = moments(ptraj.traj, N; τ=τ, order=ptraj.order)
𝔼(functional::F; τ=1e-2, N::Int=100_000) where {F<:Functional} = moments(functional, N; τ=τ)
𝔼(fp::FunctionalPower; τ=1e-2, N::Int=100_000) = moments(fp.functional, N; τ=τ, order=fp.order)

@doc raw"""
    `TimeAverage(sp::StochasticProcess, T::Real, Δ::Real)`

    The type of time average functionals.
"""
struct TimeAverage
    sp::StochasticProcess
    T::Float64
    Δ::Float64
end
δ̄² = TimeAverage

# x(t+Δ)*x(t)
struct TrajMultiplication
    sp::StochasticProcess
    T::Float64
    Δ::Float64
end

*(xt::Trajectory{SP,T1}, xs::Trajectory{SP,T2}) where {SP,T1,T2} = TrajMultiplication(xt.sp, xs.T, xt.T - xs.T)

function trajmulmean(tm::TrajMultiplication, τ, N)
    T = tm.T
    Δ = tm.Δ
    slag = round(Int, Δ / τ)
    means = zeros(nthreads())
    @threads for _ in 1:N
        _, x = tm.sp.method(T + Δ, τ, tm.sp.args...)
        @inbounds means[threadid()] += x[end] * x[end-slag]
    end
    sum(means) / N
end

𝔼(tm::TrajMultiplication; N::Int=100_000, τ=0.01) = trajmulmean(tm, τ, N)

"""
    TAMSD(TA::TimeAverage, order::Int, τ::Float64, N::Int)
    𝔼(TA::TimeAverage; τ::Float64=1e-2, N::Int=100_000, order::Int=10)

计算 Langevin 方程的时间平均均方位移 (TAMSD)

# Arguments
- `TA` : 时间平均， 由构造函数 `δ̄²(::StochasticProcess, T::Real, Δ::Real)` 构造
- `τ` : 模拟 Langevin 方程时所取的欧拉格式的步长 
- `N` : 蒙特卡罗模拟的粒子数量
- `order` : 高斯-勒让德数值积分所取正交多项式的次数

# 使用方法
```julia
T = 100; Δ = 1
x::StochasticProcess = ....  # 定义一个随机过程实例
𝔼(δ̄²(x, T, Δ))  # 计算 TAMSD
```
"""
function TAMSD(TA::TimeAverage, order::Int, τ::Float64, N::Int)
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

𝔼(TA::TimeAverage; τ::Float64=1e-2, N::Int=100_000, order::Int=10) = TAMSD(TA, order, τ, N)