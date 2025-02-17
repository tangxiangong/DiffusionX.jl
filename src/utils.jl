import FastGaussQuadrature: gausslegendre


function get_weights_nodes(a, b, order)
    nodes_unit, weights_unit = gausslegendre(order)
    weights = @. (b-a) * weights_unit / 2
    nodes = @. (b-a)*nodes_unit/2 + (b+a)/ 2
    weights, nodes
end

function firstpassagetime(domain, method, args...)
    a, b = domain
    counter = 1
    while true
        t, path = method(2^counter, args...)
        index = findfirst(path) do x
            x<=a || x>= b
        end
        !isnothing(index) && return t[index]
        counter += 1
        counter == 60 && error("NOT PASS!")
    end
end

function occupationtime(domain, method, T, args...)
    a, b = domain
    t, x = method(T, args...)
    indices = findall(x->a<=x<=b, x)
    isnothing(indices) && return 0
    length(indices) == length(x) && return T
    temp = diff(indices)
    isempty(temp) && return 0
    jump = findall(x->x!=1, temp)
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

