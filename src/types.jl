mutable struct PropagationParameters
    Nx::Int64
    Ny::Int64
    dx::Float32
    dy::Float32
    iz_start::Int64
    iz_end::Int64
    taperPerStep::Float32
    twistPerStep::Float32
    d::Float32
    n_cladding::Float32
    claddingAbsorption::Float32
    n_0::Float32
    Nshapes::Int64
    shapexs::Vector{Float32}
    shapeys::Vector{Float32}
    shapeRs::Vector{Float32}
    shapeTypes::Vector{Float32}
    shapeRIs::Vector{Float32}
    shapegs::Vector{Float32}
    shapeAbsorptions::Vector{Float32}
    shapexs_transformed::Vector{Float32}
    shapeys_transformed::Vector{Float32}
    shapeRs_transformed::Vector{Float32}
    Efinal::Matrix{ComplexF32}
    E1::Matrix{ComplexF32}
    E2::Matrix{ComplexF32}
    Eyx::Matrix{ComplexF32}
    n_out::Matrix{Float32}
    b::Vector{ComplexF32}
    multiplier::Matrix{ComplexF32}
    ax::ComplexF32
    ay::ComplexF32
    rho_e::Float32
    RoC::Float32
    sinBendDirection::Float32
    cosBendDirection::Float32
    precisePower::Float64
    EfieldPower::Float64
end

mutable struct EField
    field::Matrix{ComplexF32}
    Lx::Float64
    Ly::Float64
end

mutable struct Parameters
    name::String
    useAllCPUs::Bool
    useGPU::Bool

    figNum::Int64
    figTitle::String
    updates::Int64
    downsampleImages::Bool
    displayScaling::Int64
    disableStepsizeWarning::Bool

    Intensity_colormap::Int64
    Phase_colormap::Int64
    n_colormap::Int64

    Lx_main::Float64
    Ly_main::Float64
    Nx_main::Int64
    Ny_main::Int64
    padfactor::Int64
    dz_target::Float64
    alpha::Float64
    calcModeOverlaps::Bool

    lambda::Float64
    n_cladding::Float64
    n_0::Float64
    Lz::Float64

    z::Vector{Float64}

    shapes::Matrix{Float64}
    originalShapesInput::Matrix{Float64}
    
    E::EField
    Efunction::Union{Nothing, Function}
    originalEinput::Union{Nothing, EField}
    Einitial::Matrix{Float64}
    Eparameters::Dict{Any,Any}
    taperScaling::Int64
    twistRate::Int64
    rho_e::Float64
    bendingRoC::Float64
    bendDirection::Int64

    xzSlice::Vector{Matrix{Float64}}
    yzSlice::Vector{Matrix{Float64}}
end