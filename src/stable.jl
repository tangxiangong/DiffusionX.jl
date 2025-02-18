function sample_standard_alpha(α, β)
    half_pi = π / 2.0
    tmp = β * tan(α * half_pi)
    v = (rand() - 1 / 2) * 2 * half_pi
    w = randexp()
    b = atan(tmp) / α
    s = (1.0 + tmp * tmp)^(1.0 / (2.0 * α))
    c1 = α * sin(v + b) / cos(v)^(1.0 / α)
    c2 = cos(v - α * (v + b)) / w^(1.0 / α)
    s * c1 * c2
end

function sample_standard_alpha_one(β)
    half_pi = π / 2.0
    v = (rand() - 1 / 2) * 2 * half_pi
    w = randexp()
    c1 = (half_pi + β * v) * tan(v)
    c2 = ((half_pi * w * cos(v)) / log(half_pi + β * v)) * β
    2.0 * (c1 - c2) / π
end

function sample_stable_alpha(α, β, σ, μ)
    r = sample_standard_alpha(α, β)
    σ * r + μ
end

function sample_stable_alpha_one(α, β, σ, μ)
    r = sample_standard_alpha_one(β)
    σ * r + μ + 2.0 * β * σ * σ * log(σ) / π
end


function stable_rands(α, β, σ, μ, n)
    res = zeros(n)
    if α == 1.0
        gen = sample_stable_alpha_one
    else
        gen = sample_stable_alpha
    end

    # 使用多线程分块计算
    chunk_size = div(n, Threads.nthreads())
    Threads.@threads for tid in 1:Threads.nthreads()
        start_idx = (tid - 1) * chunk_size + 1
        end_idx = tid == Threads.nthreads() ? n : tid * chunk_size

        # 每个线程使用自己的随机数生成器
        rng = Random.default_rng()
        Random.seed!(rng, tid * time_ns())

        # 在每个块内使用 SIMD
        @simd for i in start_idx:end_idx
            @fastmath @inbounds res[i] = gen(α, β, σ, μ)
        end
    end
    res
end

function skewed_stable_rands(α, n)
    res = zeros(n)

    # 使用多线程分块计算
    chunk_size = div(n, Threads.nthreads())
    Threads.@threads for tid in 1:Threads.nthreads()
        start_idx = (tid - 1) * chunk_size + 1
        end_idx = tid == Threads.nthreads() ? n : tid * chunk_size

        # 每个线程使用自己的随机数生成器
        rng = Random.default_rng()
        Random.seed!(rng, tid * time_ns())

        # 在每个块内使用 SIMD
        @simd for i in start_idx:end_idx
            @fastmath @inbounds res[i] = sample_standard_stable(α, 1.0)
        end
    end
    res
end

@doc raw"""
    ``\alpha``-Stable distribution.

    # Parameters
    - `α`: Stability parameter.
    - `β`: Skewness parameter.
    - `σ`: Scale parameter.
    - `μ`: Location parameter.

    # Examples
    ```julia
    using DiffusionX
    S = Stable(1.0, 0.0, 1.0, 0.0)
    rand(S)
    ```
"""
struct Stable
    α::Float64
    β::Float64
    σ::Float64
    μ::Float64
    function Stable(α, β=0.0, σ=1.0, μ=0.0)
        !(0 < α < 2) && throw(ArgumentError("α must be in (0, 2)"))
        !(-1 <= β <= 1) && throw(ArgumentError("β must be in [-1, 1]"))
        !(σ > 0) && throw(ArgumentError("σ must be positive"))
        if α isa Integer
            α = Float64(α)
        end
        if β isa Integer
            β = Float64(β)
        end
        if σ isa Integer
            σ = Float64(σ)
        end
        if μ isa Integer
            μ = Float64(μ)
        end
        new(α, β, σ, μ)
    end
end

Base.show(io::IO, S::Stable) = print(io, "Stable distribution with stability index α=$(S.α), skewness β=$(S.β), scale σ=$(S.σ), location μ=$(S.μ)")

function Base.rand(S::Stable)
    if S.α ≈ 1.0
        sample_stable_alpha_one(S.α, S.β, S.σ, S.μ)
    else
        sample_stable_alpha(S.α, S.β, S.σ, S.μ)
    end
end

function Base.rand(S::Stable, n::Int)
    if n <= 0
        throw(ArgumentError("n must be positive"))
    end
    stable_rands(S.α, S.β, S.σ, S.μ, n)
end

function Base.rand(S::Stable, sz::NTuple{N,Int}) where {N}
    any(x -> x <= 0, sz) && throw(ArgumentError("sz must all be positive"))
    n = prod(sz)
    res = stable_rands(S.α, S.β, S.σ, S.μ, n)
    reshape(res, sz)
end

@doc raw"""
    ``\alpha``-SkewedStable distribution.

    # Parameters
    - `α`: Stability index. 

    # Examples
    ```julia
    using DiffusionX
    S = SkewedStable(0.5)
    rand(S)
    ```
"""
struct SkewedStable
    α::Float64
    function SkewedStable(α)
        !(0 < α < 1) && throw(ArgumentError("α must be in (0, 1)"))
        new(α)
    end
end

Base.show(io::IO, S::SkewedStable) = print(io, "totally skewed stable distribution with stability index α=$(S.α)")

function Base.rand(S::SkewedStable)
    sample_standard_stable(S.α, 1.0)
end

function Base.rand(S::SkewedStable, n::Int)
    if n <= 0
        throw(ArgumentError("n must be positive"))
    end
    skewed_stable_rands(S.α, n)
end

function Base.rand(S::SkewedStable, sz::NTuple{N,Int}) where {N}
    any(x -> x <= 0, sz) && throw(ArgumentError("sz must all be positive"))
    n = prod(sz)
    res = skewed_stable_rands(S.α, n)
    reshape(res, sz)
end

# using PoissonRandom
# import Random: randexp


# normal_poisson_rand(λ, μ, σ) = σ * randn() + μ + pois_rand(λ)

# function _rand_skewed_stable(α)
#     0 < α < 1 || throw(ArgumentError("参数在 (0,1) 之间"))
#     v = (rand() - 1 / 2) * π
#     w = randexp()
#     c₁ = (cospi(α / 2))^(-1 / α)
#     c₂ = π / 2
#     temp₁ = sin(α * (v + c₂))
#     temp₂ = (cos(v - α * (v + c₂)) / w)^(1 / α - 1)
#     temp₃ = (cos(v))^(1 / α)
#     c₁ * temp₁ * temp₂ / temp₃
# end


# function skewed_randtempered(α, γ)
#     U = rand()
#     V = _rand_skewed_stable(α)
#     U <= exp(-γ * V) && return V
#     return skewed_randtempered(α, γ)
# end

# function _unit_rand_stable(α)
#     @assert 0 < α < 1 || 1 < α < 2
#     U = π * rand() - π / 2
#     W = randexp()
#     ζ = -tanpi(α / 2)
#     ξ = -atan(-ζ) / α
#     temp₁ = (1 + ζ^2)^(1 / (2α))
#     temp₂ = sin(α * (U + ξ)) / (cos(U))^(1 / α)
#     temp₃ = cos(U - α * (U + ξ)) / W
#     temp₁ * temp₂ * temp₃^((1 - α) / α)
# end

# function _unit_randtempered(α, γ)
#     U = rand()
#     V = _unit_rand_stable(α)
#     constant = α * γ^(α - 1) / cospi(α / 2)
#     U <= exp(-γ * V) && return V - constant
#     return _unit_randtempered(α, γ)
# end

# function randtempered(α, γ)
#     Y⁺ = _unit_randtempered(α, γ)
#     Y⁻ = _unit_randtempered(α, γ)
#     2^(-1 / α) * (Y⁺ - Y⁻)
# end


# function _rand_power(α)
#     α > 0 || throw(ArgumentError("参数需为正"))
#     r = rand()
#     return (1 - r)^(-1 / α) - 1
# end

# randpow(α::Float64) = 0 < α < 1 ? _rand_skewed_stable(α) : _rand_power(α)

# tempered_normal_poisson_rand(β, γ, λ, μ, σ) = randtempered(β, γ) + normal_poisson_rand(λ, μ, σ)


# function vectorize(RNG::F, N::Integer, args...) where {F}
#     N > 0 || throw(ArgumentError("第二个参数需大于 1"))
#     N == 1 && return RNG(args...)
#     x = zeros(N)
#     Threads.@threads for k in eachindex(x)
#         @inbounds x[k] = RNG(args...)
#     end
#     x
# end

# normal_poisson_rand(λ, μ, σ, N::Int) = vectorize(normal_poisson_rand, N, λ, μ, σ)

# randpow(α, N::Int) = vectorize(randpow, N, α)

# skewed_randtempered(α, γ, N::Int) = vectorize(skewed_randtempered, N, α, γ)

# randtempered(α, γ, N::Int) = vectorize(randtempered, N, α, γ)

# tempered_normal_poisson_rand(β, γ, λ, μ, σ, N::Int) = vectorize(tempered_normal_poisson_rand, N, β, γ, λ, μ, σ)