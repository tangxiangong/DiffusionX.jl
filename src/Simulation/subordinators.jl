function subordinator(T, α, τ=0.01)
    t = collect(0:τ:T)
    x = similar(t)
    x[1] = 0.0
    @inbounds for n in firstindex(t)+1:lastindex(t)
        x[n] = x[n-1] + τ^(1/α)*randpow(α)
    end
    t, x
end

# E(t) = inf{u: S(u)>t}.
function inv_subordinator(T, α, τ)
    t_sub, s_sub = subordinator(T, α, τ)
    temp_counter = 1
    while s_sub[end] < T
        t_sub, s_sub = subordinator(T*(temp_counter+1), α, τ)
        temp_counter += 1
    end
    t = collect(0.0:τ:T)
    x = similar(t)
    @inbounds @simd for n in eachindex(t)
        index = findfirst(s->s>=t[n], s_sub)
        x[n] = t_sub[index]
    end   
    t, x
end