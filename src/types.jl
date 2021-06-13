abstract type WaveguideComponent end

mutable struct PropagationParameters
    Nx::Int64
    Ny::Int64
    dx::Float32
    dy::Float32
    dz::Float64
    iz_start::Int64
    iz_end::Int64
    taperPerStep::Float32
    twistPerStep::Float32
    d::Float32    
    n_0::Float32
    n_in::Matrix{ComplexF32}
    n_func::Union{Nothing, WaveguideComponent}
    E1::Matrix{ComplexF32}
    E2::Matrix{ComplexF32}
    Eyx::Matrix{ComplexF32}
    n_out::Matrix{ComplexF32}
    b1::Vector{ComplexF32}
    b2::Vector{ComplexF32}
    w1::Vector{ComplexF32}
    w2::Vector{ComplexF32}
    multiplier::Matrix{Float32}
    amultiplier::Matrix{Float32}
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

mutable struct RefractiveIndex
    n::Matrix{ComplexF32}
    Lx::Float64
    Ly::Float64
end