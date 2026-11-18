"""
    TethysDrift.jl

Lagrangian particle tracking and larval dispersal engine integrated with Tethys.jl.
"""
module TethysDrift

using
    Random,
    Statistics,
    LinearAlgebra,
    TOML,
    JSON3,
    Interpolations,
    Dates,
    Downloads,
    Tethys

using Tethys: available_data, DataSource

include("data/manifest.jl")
using .Manifest

include("config/configuration.jl")
include("data/store_reader.jl")
include("driver.jl")
include("biology/larval_behavior.jl")
include("analysis/larval_storage.jl")
include("analysis/empirical_analysis.jl")
include("output/visualization.jl")
include("output/biology_diagnostics.jl")

export TethysDrift

end # module
