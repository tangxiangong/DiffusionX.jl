module DiffusionX
using Random

export
    Stable,
    StochasticProcess,
    Trajectory,
    Functional,
    OccupationTime,
    TimeAverage,
    FPT

include("stable.jl")
include("interfaces.jl")
# include("utils.jl")
# include("simulation/simulation.jl")
# greet() = print("Hello World!")

end # module DiffusionX
