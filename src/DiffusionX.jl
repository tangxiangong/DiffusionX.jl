module DiffusionX
using Random
using Base.Threads

export
    # interfaces
    Stable,
    StochasticProcess,
    Trajectory,
    Functional,
    OccupationTime,
    TimeAverage,
    FPT,
    TAMSD,
    simulate,
    moment,
    𝔼,
    δ̄²,
    # processes
    Bm


# files
include("stable.jl")
include("interfaces.jl")
include("bm.jl")

end # module DiffusionX
