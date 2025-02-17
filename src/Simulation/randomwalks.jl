function _randomwalk(T, x₀, wait, jump, wait_args, jump_args)
    t = Vector{Float64}()
    x = Vector{Float64}()
    sizehint!(t, round(Int, T))
    sizehint!(x, round(Int, T))
    current_time = 0.0
    current_position = x₀
    while current_time <= T
        append!(t, current_time)
        append!(x, current_position)
        τ = wait(wait_args...)
        ξ = jump(jump_args...)
        current_time += τ
        current_position += ξ
    end
    if !(t[end] ≈ T) 
        append!(t, T)
        append!(x, x[end])
    end
    t, x
end

telomere_randomwalk(T, τ, x₀, α, λ, μ, σ) = _randomwalk(T, x₀, randpow, normal_poisson_rand, (α,), (λ, μ, σ))

telomere_temporal_tempered_randomwalk(T, τ, x₀, α, γ, λ, μ, σ) = _randomwalk(T, x₀, skewed_randtempered, normal_poisson_rand, (α, γ), (λ, μ, σ))

telomere_spatial_tempered_randomwalk(T, τ, x₀, α, β, γ, λ, μ, σ) = _randomwalk(T, x₀, randpow, tempered_normal_poisson_rand, (α,), (β, γ, λ, μ, σ))