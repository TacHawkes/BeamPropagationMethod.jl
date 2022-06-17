function fd_bpm!(m::Model)
    # Validations
    isempty(m.n.n) && @error("You have to initialize the refractive index object.")
    isempty(m.E.field) && @error("You have to initialize the electric field object.")

    if size(m.n.n, 3) > 1
        m.twist_rate != 0.0 && @error("You cannot specify both twisting and a 3D refractive index profile.")
        m.taper_scaling != 1.0 && @error("You cannot specify both tapering and a 3D refractive index profile.")
    end

    if (m.xsymmetry != NoSymmetry && !isinf(m.bending_radius_of_curvature) && sind(m.bend_direction) != 0.0
        || m.ysymmetry != NoSymmetry && !isinf(m.bending_radius_of_curvature) && cosd(m.bend_direction != 0.0))
        @error("The specified bending direction is inconsistent with the symmetry assumption")
    end

    k0 = 2π/m.λ

    _dx = dx(m)
    _dy = dy(m)
    _nx = nx(m)
    _ny = ny(m)
    _lx = lx(m)
    _ly = ly(m)

    # TODO: Modes overlap

    _x = x(m)
    _y = y(m)

    ix_plot = 1:_nx
    iy_plot = 1:_ny
    x_plot = _x[ix_plot]
    y_plot = _y[iy_plot]

    if !m.prior_data
        power_fraction = 1/(1 + Int(m.E.xsymmetry != NoSymmetry)) / (1 + Int(m.E.ysymmetry != NoSymmetry))
        m.E.field = m.E.field ./ √(sum(abs2.(m.E.field))/power_fraction)
    end
    power_fraction = 1/(1 + Int(m.xsymmetry != NoSymmetry)) / (1 + Int(m.ysymmetry != NoSymmetry))

    # interpolation
    nx_source, ny_source = size(m.E.field)
    dx_source = m.E.lx / nx_source
    dy_source = m.E.ly / ny_source
    x_source = get_grid_array(nx_source, dx_source, m.E.ysymmetry)
    y_source = get_grid_array(ny_source, dy_source, m.E.xsymmetry)
    x_source, y_source, E_source = calc_full_field(x_source, y_source, m.E.field)
    itp = Interpolations.LinearInterpolation((x_source, y_source), E_source; extrapolation_bc=zero(eltype(E_source)))
    E = Matrix{ComplexF32}(undef, length(_x), length(_y))
    for j in eachindex(_y)
        for i in eachindex(_x)
            E[i,j] = itp(_x[i], _y[j])
        end
    end
    E_source_norm = sum(abs2.(E_source))
    E_norm = sum(abs2.(E)/power_fraction)
    E .*= √(E_source_norm / E_norm)

    m.ysymmetry == AntiSymmetric && (E[_x .== 0, :] .= zero(eltype(E)))
    m.xsymmetry == AntiSymmetric && (E[:, _y .== 0] .= zero(eltype(E)))

    # refractive index init
    dz_n = 0
    nx_source, ny_source, nz_source = size(m.n.n)
    dx_source = m.n.lx / nx_source
    dy_source = m.n.ly / ny_source
    x_source = get_grid_array(nx_source, dx_source, m.n.ysymmetry)
    y_source = get_grid_array(ny_source, dy_source, m.n.xsymmetry)
    x_source, y_source, n_source = calc_full_refractive_index(x_source, y_source, m.n.n)
    if nz_source == 1
        # 2d
        itp = Interpolations.LinearInterpolation((x_source, y_source), @view(n_source[:,:,1]); extrapolation_bc=m.nbackground)
        n = Array{eltype(n_source), 3}(undef, length(_x), length(_y), 1)
        for j in eachindex(_y)
            for i in eachindex(_x)
                n[i,j,:] .= itp(_x[i], _y[j])
            end
        end
        n_slice = n
    else
        dz_n = m.lz/(nz_source-1)
        z_source = dz_n * (0:nz_source - 1)
        itp = Interpolations.LinearInterpolation((x_source, y_source, z_source), n_source; extrapolation_bc=m.nbackground)
        n = Array{eltype(n_source), 3}(undef, length(_x), length(_y), length(z_source))
        for k in eachindex(z_source)
            for j in eachindex(_y)
                for i in eachindex(_x)
                    n[i,j,k] = itp(_x[i], _y[j], z_source[k])
                end
            end
        end
        n_slice = @view(n[:,:,1])
        n = trim_refractive_index(n, m.nbackground)
    end

    nx_n, ny_n, nz_n = size(n)
    # TODO: Volumetric plot
    _nz = nz(m)
    _dz = dz(m)

    if !m.disable_stepsize_warning
        max_a = 5
        max_d = 15
        dz_max1 = max_a*4*_dx^2*k0*m.n0
        dz_max2 = max_a*4*_dy^2*k0*m.n0
        dz_max3 = max_d*2*m.n0/k0
        min_dz = minimum([dz_max1 dz_max2 dz_max3])
        if _dz > min_dz
            @warn("z step size is high (> $(min_dz) m), which may introduce numerical
            artifacts. Alleviate this by decreasing dz_target, decreasing nx_main and/or
            ny_main, or increasing lx_main and/or ly_main. You can disable this warning
            by setting m.disable_stepsize_warning = true.")
        end
    end

    z_update_idxs = round.(Int64, collect((1:m.updates)/m.updates*_nz))
    if m.prior_data
        m.z = [m.z; _dz*z_update_idxs .+ m.z[end]]
    else
        m.z = [0; _dz*z_update_idxs]
    end

    # porportionality factors
    ax = _dz/(4im*_dx^2*k0*m.n0)
    ay = _dz/(4im*_dy^2*k0*m.n0)
    d = -_dz*k0

    # edge absorption
    xedge = m.lx_main * (1 + Int(m.ysymmetry != NoSymmetry))/2
    yedge = m.ly_main * (1 + Int(m.xsymmetry != NoSymmetry))/2
    multiplier = Matrix{Float32}(undef, length(_x), length(_y))
    for j in eachindex(_y)
        for i in eachindex(_x)
            multiplier[i,j] = exp(-_dz*max(0, max(abs(_y[j]) - yedge, abs(_x[i]) - xedge))^2*m.alpha)
        end
    end

    # visualisation
    if isnothing(m.fig)
        m.fig = Figure()
    else
        empty!(m.fig)
    end

    if m.save_video && isnothing(m.video_stream)
        m.video_stream = VideoStream(m.fig; framerate=10)
    end

    xlimits = ([-1 1] .+ Int(m.ysymmetry != NoSymmetry)) .* _lx/(2*m.plot_zoom)
    ylimits = ([-1 1] .+ Int(m.xsymmetry != NoSymmetry)) .* _ly/(2*m.plot_zoom)

    redline_x, redline_y = redlines(m::Model)

    # TODO: Axis options
    axis1 = m.fig[1,1][1,1] = Axis(m.fig,
                            xlabel="x [μm]",
                            ylabel="y [μm]",
                            title="Real part of refractive index",
                            aspect=DataAspect(),
                            limits=Tuple(1e6.*[xlimits ylimits])
    )
    ri_data = Observable(real(n_slice[ix_plot, iy_plot]))

    hm1a = heatmap!(
        axis1,
        x_plot*1e6,
        y_plot*1e6,
        ri_data,
        colormap=m.refractive_index_colormap
    )
    Colorbar(m.fig[1,1][1,2], hm1a, tickformat="{:.3f}", ticks=LinearTicks(5))

    update_stored_fields!(m, _nx, _ny, E)

    axis2 = m.fig[1,2] = Axis(m.fig,
        xlabel="Propagation distance [mm]",
        ylabel="Relative power remaining",
        yminorticksvisible = true,
        yminorgridvisible = true
    )
    power_plot_data = Observable([Point2f(x_i,y_i) for (x_i, y_i) in zip(m.z*1e3, m.powers)])
    plot2 = lines!(axis2, power_plot_data, linewidth=2)
    xlims!(axis2, 0, m.z[end]*1e3)
    ylims!(axis2, 0, 1.1)

    axis3 = m.fig[2,1][1,1] = Axis(m.fig,
                                    xlabel="x [μm]",
                                    ylabel="y [μm]",
                                    title="Intensity [W/cm^2]",
                                    aspect=DataAspect(),
                                    limits=Tuple(1e6.*[xlimits ylimits])
    )
    field_plot_data = Observable(abs2.(E[ix_plot,iy_plot]).*1e4)
    hm3a = heatmap!(axis3,
                    x_plot*1e6,
                    y_plot*1e6,
                    field_plot_data,
                    colormap=m.intensity_colormap,
                    colorrange=(m.plot_field_max > 0.0 ? (0, 1e4*m.plot_field_max*maximum(abs2.(E))) : CairoMakie.Makie.Automatic())
    )
    lines!(axis3, vec(redline_x), vec(redline_y), color=:red, linestyle=:dash, linewidth=2)
    Colorbar(m.fig[2,1][1,2], hm3a, tickformat="{:.3f}", ticks=LinearTicks(5))

    axis4 = m.fig[2,2][1,1] = Axis(m.fig,
                                    xlabel="x [μm]",
                                    ylabel="x [μm]",
                                    title="Phase [rad]",
                                    limits=Tuple(1e6.*[xlimits ylimits]),
                                    aspect=DataAspect()

    )
    maxE0 = maximum(abs.(E))
    alpha_data = max.(0, @. (1 + log10(abs(E[ix_plot,iy_plot]/maxE0)^2)/3))
    Φ = angle.(E[ix_plot, iy_plot])
    Φ[alpha_data .== 0] .= NaN
    phase_plot_data = Observable(Φ)

    hm3b = heatmap!(axis4,
                    x_plot*1e6,
                    y_plot*1e6,
                    phase_plot_data,
                    colormap=m.phase_colormap,
                    colorrange=(-π,π)
    )
    Colorbar(m.fig[2,2][1,2], hm3b,ticks=(-π:π:π, ["-π", "0", "π"]))
    p = Parameters(;
            m,
            dx=_dx,
            dy=_dy,
            dz=_dz,
            taper_per_step=((1-m.taper_scaling)/_nz),
            twist_per_step=(m.twist_rate*m.lz/_nz),
            multiplier,
            n_in=n,
            n_out=zeros(ComplexF32, size(n, 1), size(n,2)),
            dz_n,
            d,
            n0=m.n0,
            ax,
            ay,
            radius_of_curvature=m.bending_radius_of_curvature,
            ρe=m.ρe,
            sin_bend_direction=sind(m.bend_direction),
            cos_bend_direction=cosd(m.bend_direction),
            precise_power=m.powers[end-length(z_update_idxs)]*power_fraction,
            xsymmetry=m.xsymmetry,
            ysymmetry=m.ysymmetry,
            E1=E,
            E2=similar(E),
            Efinal=similar(E),
            Eyx=similar(E),
            b1=zeros(ComplexF32, max(size(E)[1],size(E)[2])),
            b2=zeros(ComplexF32, max(size(E)[1],size(E)[2])),
            w1=zeros(ComplexF32, max(size(E)[1],size(E)[2])),
            w2=zeros(ComplexF32, max(size(E)[1],size(E)[2])),
            nx=size(E,1),
            ny=size(E,2),
            nx_n=size(n,1),
            ny_n=size(n,2),
            nz_n=size(n,3),
            iz_start=1,
            iz_end=z_update_idxs[1]
    )
    !isempty(m.fig_title) && Label(m.fig[0, :], text=m.fig_title)
    trim!(m.fig.layout)
    display(m.fig)

    m.save_video && recordframe!(m.video_stream)

    @showprogress "Calculating beam propagation..." for updidx in eachindex(z_update_idxs)
        if updidx > 1
            p.iz_start = z_update_idxs[updidx-1]
            p.iz_end = z_update_idxs[updidx]
        end

        fdbpm_propagator!(p)

        precise_power = p.precise_power
        m.powers[end-length(z_update_idxs) + updidx] = precise_power/power_fraction

        # update observables
        ri_data[] = real(p.n_out[ix_plot, iy_plot])
        field_plot_data[] = abs2.(p.E2[ix_plot,iy_plot]).*1e4
        alpha_data = max.(0, @. (1 + log10(abs(p.E2[ix_plot,iy_plot]/maxE0)^2)/3))
        maxE0 = maximum(abs.(p.E2))
        Φ = angle.(p.E2[ix_plot, iy_plot])
        Φ[alpha_data .== 0] .= NaN
        phase_plot_data[] = Φ
        power_plot_data[][end-length(z_update_idxs) + updidx] = Point2f(m.z[end-length(z_update_idxs) + updidx]*1e3, m.powers[end-length(z_update_idxs) + updidx])

        if m.store_field_3D
            m.E3D[end][:,:,end-length(z_update_idxs)+updidx] = p.E2
        end

        recordframe!(m.video_stream)

        swap_pointers!(p.E1, p.E2)
    end

    m.E.lx = lx(m)
    m.E.ly = ly(m)
    m.E.field = p.E2
    m.E.xsymmetry = m.xsymmetry
    m.E.ysymmetry = m.ysymmetry

    m.n.lx = lx(m)
    m.n.ly = ly(m)
    m.n.n = reshape(p.n_out, size(p.n_out)..., :)
    m.n.xsymmetry = m.xsymmetry
    m.n.ysymmetry = m.ysymmetry

    m.prior_data = true
    display(m.fig)
end
