@enum Symmetry begin
    NoSymmetry
    Symmetric
    AntiSymmetric
end

Base.@kwdef mutable struct ElectricFieldProfile
    lx::Float64 = 1.0
    ly::Float64 = 1.0
    field::Matrix{ComplexF32} = Matrix{ComplexF32}(undef, 0, 0)
    xsymmetry::Symmetry = NoSymmetry
    ysymmetry::Symmetry = NoSymmetry
    label::String = ""
    n_effective::Float64 = NaN
end

Base.@kwdef mutable struct RefractiveIndexProfile
    lx::Float64 = 1.0
    ly::Float64 = 1.0
    n::Array{ComplexF32, 3} = Array{ComplexF32}(undef, 0, 0, 0)
    xsymmetry::Symmetry = NoSymmetry
    ysymmetry::Symmetry = NoSymmetry
end

Base.@kwdef mutable struct Model
    # visualisation
    name::String = "BPM-Julia model " * (Dates.format(now(), dateformat"YYYY/m/d HH:MM:SS"))
    fig_title::String = ""
    fig::Union{Figure, Nothing} = nothing
    updates::UInt = 50
    plot_field_max::Float64 = 0.0
    plot_zoom::Float64 = 1.0
    store_field_3D::Bool = false
    save_video::Bool = false
    intensity_colormap::Symbol = :viridis
    phase_colormap::Symbol = :hsv
    refractive_index_colormap::Symbol = :viridis
    refractive_index_colorrange::Tuple{Int, Int} = (1, 0)
    calc_mode_overlaps::Bool = false
    disable_stepsize_warning::Bool = false
    disable_plot_time_warning::Bool = false
    disable_downsampling_warning::Bool = false

    # solver
    use_all_cpus::Bool = false
    use_gpu::Bool = false
    nx_main::UInt = 2
    ny_main::UInt = 2
    xsymmetry::Symmetry = NoSymmetry
    ysymmetry::Symmetry = NoSymmetry
    dz_target::Float64 = 1e-6
    padfactor::Float64 = 1.5
    alpha::Float64 = 3e14

    # geometry
    lx_main::Float64 = 1.0
    ly_main::Float64 = 1.0
    lz::Float64 = 1.0
    taper_scaling::Float64 = 1.0
    twist_rate::Float64 = 0.0
    bending_radius_of_curvature::Float64 = Inf
    bend_direction::Float64 = 0.0

    # optical and material parameters
    λ::Float64 = 1.0
    nbackground::Float64 = 1.0
    n0::Float64 = 1.0
    ρe::Float64 = 0.22

    # refractice index profile
    n::RefractiveIndexProfile = RefractiveIndexProfile()

    # electric field profile
    E::ElectricFieldProfile = ElectricFieldProfile()

    # internals
    modes::ElectricFieldProfile = ElectricFieldProfile()
    prior_data::Bool = false
    z::Vector{Float64} = Vector{Float64}(undef, 0)
    powers::Vector{Float64} = Vector{Float64}(undef, 0)
    xzslice::Vector{Matrix{ComplexF32}} = Matrix{ComplexF32}[]
    yzslice::Vector{Matrix{ComplexF32}} = Matrix{ComplexF32}[]
    E3D::Vector{Array{ComplexF32,3}} = Array{ComplexF32,3}[]
end

Base.@kwdef mutable struct Parameters
    m::Model
    nx::Int64
    ny::Int64
    dx::Float32
    dy::Float32
    dz::Float64
    iz_start::Int64
    iz_end::Int64
    xsymmetry::Symmetry
    ysymmetry::Symmetry
    taper_per_step::Float32
    twist_per_step::Float32
    d::Float32
    n0::Float32
    n_in::Array{ComplexF32,3}
    nx_n::Int64
    ny_n::Int64
    nz_n::Int64
    dz_n::Float32
    Efinal::Matrix{ComplexF32}
    E1::Matrix{ComplexF32}
    E2::Matrix{ComplexF32}
    Eyx::Matrix{ComplexF32}
    n_out::Array{ComplexF32,2}
    b1::Vector{ComplexF32}
    b2::Vector{ComplexF32}
    w1::Vector{ComplexF32}
    w2::Vector{ComplexF32}
    multiplier::Matrix{Float32}
    ax::ComplexF32
    ay::ComplexF32
    ρe::Float32
    radius_of_curvature::Float32
    sin_bend_direction::Float32
    cos_bend_direction::Float32
    precise_power::Float64
    precise_power_diff::Float32 = zero(Float32)
    E_field_power::Float64 = zero(Float64)
end
