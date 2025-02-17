using PoissonRandom
import Random: randexp


normal_poisson_rand(λ, μ, σ) = σ * randn() + μ + pois_rand(λ)

function _rand_skewed_stable(α)
    0 < α < 1 || throw(ArgumentError("参数在 (0,1) 之间"))
    v = (rand() - 1 / 2) * π
    w = randexp()
    c₁ = (cospi(α / 2))^(-1 / α)
    c₂ = π / 2
    temp₁ = sin(α * (v + c₂))
    temp₂ = (cos(v - α * (v + c₂)) / w)^(1 / α - 1)
    temp₃ = (cos(v))^(1 / α)
    c₁ * temp₁ * temp₂ / temp₃
end


function skewed_randtempered(α, γ)
    U = rand()
    V = _rand_skewed_stable(α)
    U <= exp(-γ * V) && return V
    return skewed_randtempered(α, γ)
end

function _unit_rand_stable(α)
    @assert 0 < α < 1 || 1 < α < 2
    U = π * rand() - π / 2
    W = randexp()
    ζ = -tanpi(α / 2)
    ξ = -atan(-ζ) / α
    temp₁ = (1 + ζ^2)^(1 / (2α))
    temp₂ = sin(α * (U + ξ)) / (cos(U))^(1 / α)
    temp₃ = cos(U - α * (U + ξ)) / W
    temp₁ * temp₂ * temp₃^((1 - α) / α)
end

function _unit_randtempered(α, γ)
    U = rand()
    V = _unit_rand_stable(α)
    constant = α * γ^(α - 1) / cospi(α / 2)
    U <= exp(-γ * V) && return V - constant
    return _unit_randtempered(α, γ)
end

function randtempered(α, γ)
    Y⁺ = _unit_randtempered(α, γ)
    Y⁻ = _unit_randtempered(α, γ)
    2^(-1 / α) * (Y⁺ - Y⁻)
end


function _rand_power(α)
    α > 0 || throw(ArgumentError("参数需为正"))
    r = rand()
    return (1 - r)^(-1 / α) - 1
end

randpow(α::Float64) = 0 < α < 1 ? _rand_skewed_stable(α) : _rand_power(α)

tempered_normal_poisson_rand(β, γ, λ, μ, σ) = randtempered(β, γ) + normal_poisson_rand(λ, μ, σ)


function vectorize(RNG::F, N::Integer, args...) where {F}
    N > 0 || throw(ArgumentError("第二个参数需大于 1"))
    N == 1 && return RNG(args...)
    x = zeros(N)
    Threads.@threads for k in eachindex(x)
        @inbounds x[k] = RNG(args...)
    end
    x
end

normal_poisson_rand(λ, μ, σ, N::Int) = vectorize(normal_poisson_rand, N, λ, μ, σ)

randpow(α, N::Int) = vectorize(randpow, N, α)

skewed_randtempered(α, γ, N::Int) = vectorize(skewed_randtempered, N, α, γ)

randtempered(α, γ, N::Int) = vectorize(randtempered, N, α, γ)

tempered_normal_poisson_rand(β, γ, λ, μ, σ, N::Int) = vectorize(tempered_normal_poisson_rand, N, β, γ, λ, μ, σ)