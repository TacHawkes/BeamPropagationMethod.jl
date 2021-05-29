module BeamPropagationMethod

using   Arpack,
        Interpolations, 
        LinearAlgebra,
        Plots,
        Printf, 
        SparseArrays

include("types.jl")
include("util.jl")
include("propagator.jl")
include("modes.jl")
include("fd_bpm.jl")
include("fft_bpm.jl")

end
