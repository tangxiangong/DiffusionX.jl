# include("_randomwalk1d.jl")
# include("_randomwalk2d.jl")
# include("_randomwalk3d.jl")

struct RandomWalk{N}
    step_num::Integer
    start_point::NTuple{N, Real}
    step_size::Real
    prob::NTuple
    dims::Integer
    function RandomWalk{N}(step_num::Integer,
        start_point::NTuple{N, Real},
        step_size::Real,
        prob::NTuple,
        dims=N::Integer) where N
        @assert step_num >= 1 && step_size > zero(typeof(step_size))
        @assert N ∈ [1, 2, 3]
        @assert length(prob) == 2dims
        @assert sum(prob) ≈ 1 && all(p->0<=p<=1, prob)
        new{N}(step_num, start_point, step_size, prob, dims)
    end
end

RandomWalk{N}(step_num, start_point=ntuple(i->0, Val(N)), step_size=1, prob=ntuple(i-> 1/(2N), Val(2N))) where N = RandomWalk{length(start_point)}(step_num, start_point, step_size, prob)

RandomWalk1D = RandomWalk{1}
RandomWalk2D = RandomWalk{2}
RandomWalk3D = RandomWalk{3}

include("_randomwalk1d.jl")
include("_randomwalk2d.jl")
include("_randomwalk3d.jl")
