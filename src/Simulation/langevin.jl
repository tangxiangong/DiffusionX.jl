#  dx(s(t)) = (λ+μ)ds(t) + σdB(s(t))
function _langevin!(τs, x, λ, μ, σ, x₀)
    x[1] = x₀
    @inbounds for n in eachindex(τs)
        x[n+1] = x[n] + (λ + μ)τs[n] + σ * sqrt(τs[n]) * randn()
    end
end

function telomere_langevin(T, τ, x₀, α, λ, μ, σ)
    t = collect(0.0:τ:T)
    if α == 1
        τs = τ * ones(length(t) - 1)
    else
        _, s = inv_subordinator(T, α, τ)
        τs = diff(s)
    end
    x = similar(t)
    _langevin!(τs, x, λ, μ, σ, x₀)
    t, x
end