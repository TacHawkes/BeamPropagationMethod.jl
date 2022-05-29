module BeamPropagationMethod

using Printf: Threads
using   Arpack,
        FFTW,
        CairoMakie,
        Dates,
        Interpolations,
        LinearAlgebra,
        ProgressMeter,
        Printf,
        SparseArrays

include("types.jl")
include("field.jl")
include("model.jl")
include("util.jl")
include("propagator.jl")
include("modes.jl")
include("fd_bpm.jl")
include("fft_bpm.jl")

export fdbpm!, fftbpm!, find_modes!, mode_superposition, PropagationParameters, propagation_length, Simple2x2FiberCoupler

end
