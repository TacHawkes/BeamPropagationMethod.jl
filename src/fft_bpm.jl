function fft_bpm!(m::Model)
    m.dz_target = m.lz/m.updates
    nz = nz(m)

    _dx = dx(m)
    _dy = dy(m)
    _dz = dz(m)
    _nx = nx(m)
    _ny = ny(m)
    _lx = lx(m)
    _ly = ly(m)

    kx = 2π/_lx*[0:floor(Int64,(_nx-1)/2);floor(Int64,-(_nx-1)/2):-1]
    ky = 2π/_ly*[0:floor(Int64,(_ny-1)/2);floor(Int64,-(_ny-1)/2):-1]

    prop_kernel = [exp(im*_dz*(kxi^2 + kyi^2)*m.λ/(4π*m.n0)) for (kxi, kyj) in Iterators.product(kx,ky)]

    x = _dx*(-(_nx-1)/2:(_nx-1)/2)
    y = _dy*(-(_ny-1)/2:(_ny-1)/2)

    absorber = [exp(-_dz*max(0, max(abs(yj) - m.ly_main/2, abs(xi) - m.lx_main/2))^2*m.alpha) for (xi, yj) in Iterators.product(x,y)]

    if m.prior_data
        m.z = [m.z _dz*(1:_nz) + m.z[end]]
    else
        m.z = [0 _dz*(1:_nz)]
    end

    # beam initialization
    nx_source, ny_source = size(m.E.field)
    dx_source = m.E.lx / nx_source
    dy_source = m.E.ly / ny_source
    x_source = get_grid_array(nx_source, dx_source, m.E.ysymmetry)
    y_source = get_grid_array(ny_source, dy_source, m.E.xsymmetry)
    x_source, y_source, E_source = calc_full_field(x_source, y_source, m.E.field)
    itp = Interpolations.LinearInterpolation((x_source, y_source), E_source; extrapolation_bc=zero(eltype(E_source)))
    E = [itp(xi, yj) for (xi, yj) in Iterators.product(_x, _y)]
    E *= √(sum(abs2.(E_source))/sum(abs2.(E)))

    # visualisation
    if isnothing(m.fig)
        m.fig = Figure()
    else
        empty!(m.fig)
    end

    if m.save_video && isnothing(m.video_stream)
        m.video_stream = VideoStream(m.fig; framerate=10)
    end

    update_stored_fields!(m, _nx, _ny, E)
    redline_x, redline_y = redlines(m::Model)

    xlims = [-1 1]*_lx/(2*m.plot_zoom)
    ylims = [-1 1]*_ly/(2*m.plot_zoom)

    axis2 = m.fig[1,2] = Axis(m.fig,
                            xlabel = "Propagation distance [m]",
                            ylabel = "Relative power remaining",
                            yminorticksvisible = true,
                            yminorgridvisible = true
    )
    xlims!(axis2, 0, m.z[end]*1e3)
    ylims!(axis2, 0, 1.1)
    power_plot_data = Observable([Point2f(x_i,y_i) for (x_i, y_i) in zip(m.z*1e3, m.powers)])
    plot2 = lines!(axis2, power_plot_data, linewidth=2)

    axis3a = m.fig[2,1][1,1] = Axis(m.fig,
                            xlabel = "x [m]",
                            ylabel = "y [m]",
                            title = "Intensity [W/m^2]",
                            limits = Tuple([xlimits ylimits])
    )
    field_plot_data = Observable(abs2.(E[ix_plot,iy_plot]).*1e4)
    hm3a = heatmap!(axis3a,
                    x_plot*1e6,
                    y_plot*1e6,
                    field_plot_data,
                    colormap=m.intensity_colormap,
                    colorrange=(m.plot_field_max > 0.0 ? (0, 1e4*m.plot_field_max*maximum(abs2.(E))) : CairoMakie.Makie.Automatic())
    )
    lines!(axis3, vec(redline_x), vec(redline_y), color=:red, linestyle=:dash, linewidth=2)
    Colorbar(m.fig[2,1][1,2], hm3a, tickformat="{:.3f}", ticks=LinearTicks(5))

    axis3b = m.fig[2,2][1,1] = Axis(m.fig,
                             xlabel = "x [m]",
                             ylabel = "x [m]",
                             title = "Phase [rad]",
                             limits = Tuple([xlimits ylimits]),
                             aspect = DataAspect()
    )
    maxE0 = maximum(abs.(E))
    alpha_data = max.(0, @. (1 + log10(abs(E[ix_plot,iy_plot]/maxE0)^2)/3))
    Φ = angle.(E[ix_plot, iy_plot])
    Φ[alpha_data .== 0] .= NaN
    phase_plot_data = Observable(Φ)

    hm3b = heatmap!(axis3b,
                    x_plot*1e6,
                    y_plot*1e6,
                    phase_plot_data,
                    colormap=m.phase_colormap,
                    colorrange=(-π,π)
    )
    Colorbar(m.fig[2,2][1,2], hm3b,ticks=(-π:π:π, ["-π", "0", "π"]))

    display(m.fig)
    m.save_video && recordframe!(m.video_stream)

    Pfft = plan_fft!(similar(E); flags=FFTW.PATIENT, timelimit=5)
    Pifft = inv(Pfft)

    for zidx in 1:nz
        # propagate
        # forward-FFT
        Pfft * E
        Threads.@threads for i in eachindex(E)
            E[i] *= prop_kernel[i]
        end
        # inverse-FFT
        Pifft * E
        Threads.@threads for i in eachindex(E)
            E[i] *= absorber[i]
        end

        if m.store_field_3D
            m.E3D[end][:,:,end-nz+zidx] = E
        end
        m.powers[end-nz+zidx] = sub(abs2.(E))
        power_plot_data[][end-nz+zidx] = Point2f(m.z[end-nz+zidx]*1e3, m.powers[end-nz+zidx])
        field_plot_data[] = abs2.(E).*1e4
        alpha_data = max.(0, @. (1 + log10(abs(E/maxE0)^2)/3))
        maxE0 = maximum(abs.(E))
        Φ = angle.(E)
        Φ[alpha_data .== 0] .= NaN
        phase_plot_data[] = Φ

        recordframe!(m.video_stream)
    end

    @info "FFT-BPM iteration completed"
end
