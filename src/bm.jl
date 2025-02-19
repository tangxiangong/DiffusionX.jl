@doc raw"""
    `Bm(xₒ=0, D=1)`

    The Brownian motion process. 
    ``B(t) \sim \mathcal{N}(xₒ, 2Dt)``

    # Arguments
    - `xₒ`: The initial value.
    - `D`: The diffusion coefficient.
"""
struct Bm <: StochasticProcess
    xₒ::Float64
    D::Float64
    function Bm(xₒ=0, D=1)
        if D isa Integer
            D = Float64(D)
        end
        if xₒ isa Integer
            xₒ = Float64(xₒ)
        end
        D <= 0 && throw(ArgumentError("D must be positive"))
        new(xₒ, D)
    end
end

Base.show(io::IO, bm::Bm) = print(io, "Brownian motion with initial position $(bm.xₒ) and diffusion coefficient $(bm.D)")

function simulate(bm::Bm, T::Union{Int,Float64}; τ::Float64=0.01)
    x₀, D = bm.xₒ, bm.D
    n = ceil(Int, T / τ)
    x = zeros(Float64, n + 1)
    t = zeros(Float64, n + 1)
    x[1] = xₒ
    t[1] = 0.0
    @inbounds for i in 1:n
        x[i+1] = x[i] + sqrt(2.0 * D * τ) * randn()
        t[i+1] = t[i] + τ
    end
    t, x
end
