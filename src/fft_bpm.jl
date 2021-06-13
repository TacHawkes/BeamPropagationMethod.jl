function fftbpm!(p)
    # check fields
    get!(p, :name, @__FILE__)    
    get!(p, :figNum, 1)
    get!(p, :figTitle, "")
    get!(p, :finalizeVideo, false)    
    get!(p, :saveVideo, false)
    get!(p, :finalizeVideo, false)    
    if p[:saveVideo] && !haskey(p, :videoName)
        p[:videoName] = p[:name] * ".mp4"
    end
    get!(p, :Intensity_colormap, 1)
    get!(p, :Phase_colormap, 3)        
    
    Nz = p[:updates]

    targetLx = p[:padfactor] * p[:Lx_main]
    targetLy = p[:padfactor] * p[:Ly_main]

    dx = p[:Lx_main] / p[:Nx_main]
    dy = p[:Ly_main] / p[:Ny_main]
    dz = p[:Lz]/Nz

    Nx = round(Int64, targetLx/dx)
    if Nx % 2 != p[:Nx_main] % 2
        Nx += 1
    end
    Ny = round(Int64, targetLy/dy)
    if Ny % 2 != p[:Ny_main] % 2
        Ny += 1
    end

    Lx = Nx*dx
    Ly = Ny*dy

    kx = 2π/Lx*[0:floor(Int64,(Nx-1)/2);floor(Int64,-(Nx-1)/2):-1]
    ky = 2π/Ly*[0:floor(Int64,(Ny-1)/2);floor(Int64,-(Ny-1)/2):-1]
    kX, kY = ndgrid(convert.(Float32, kx), convert.(Float32, ky))

    prop_kernel = @. exp(im*dz*(kX^2 + kY^2)*p[:lambda]/(4π*p[:n_0]))

    x = collect(dx*(-(Nx-1)/2:(Nx-1)/2))
    y = collect(dy*(-(Ny-1)/2:(Ny-1)/2))
    X, Y = ndgrid(convert.(Float32, x), convert.(Float32, y))

    absorber = @. exp(-dz*max(0, max(abs(Y) - p[:Ly_main]/2, abs(X) - p[:Lx_main]/2))^2*p[:alpha])

    if isa(p[:E], Function)
        E = p[:E](X, Y, p[:Eparameters])
        p[:E_0] = E        
    else
        Nx_source, Ny_source = size(p[:E].field)
        dx_source = p[:E].Lx/Nx_source
        dy_source = p[:E].Ly/Ny_source
        x_source = dx_source*(-(Nx_source - 1)/2:(Nx_source-1)/2)
        y_source = dy_source*(-(Ny_source - 1)/2:(Ny_source-1)/2)        
        E = Interpolations.LinearInterpolation((x_source, y_source), p[:E].field; extrapolation_bc=Line())            
        E = E.(X, Y)        
    end

    if p[:saveVideo]
        if haskey(p, :videoHandle)
            video = p[:videoHandle]
        else
            video = Animation()
        end
    end

    # init plot
    plots = plot(
        layout=(1,2),
        size=(1000, 800),
        framestyle=:box
    )
    plot_idx = 1

    heatmap!(x*1e6, y*1e6, abs.(E).^2, c=get_colormap(p[:Intensity_colormap]), sp=1)
    Efield_plot = plots.series_list[plot_idx]
    Efield_plot.plotattributes[:z] = abs.(E).^2

    plot_idx += 1
    plot!(
        xlims=([-1, 1]*1e6*Lx/(2*p[:displayScaling])),
        ylims=([-1, 1]*1e6*Ly/(2*p[:displayScaling])),
        xlabel="x [μm]",
        ylabel="y [μm]",
        title="Intensity [W/m^2]",
        legend=false,
        aspect_ratio=:equal,
        sp=1
    )    
    if haskey(p, :plotEmax)
        plot!(
            clims=(0, p[:plotEmax]*maximum(abs.(E).^2)),
            sp=1
        )    
    end
    plot!(1e6.*[-p[:Lx_main], p[:Lx_main], p[:Lx_main], -p[:Lx_main], -p[:Lx_main]]/2,1e6.*[p[:Ly_main], p[:Ly_main], -p[:Ly_main], -p[:Ly_main], p[:Ly_main]]/2,c=:red,linestyle=:dash,sp=1)
    plot_idx += 1
    heatmap!(x*1e6, y*1e6, angle.(E), c=get_colormap(p[:Phase_colormap]), sp=2)        
    phase_plot = plots.series_list[plot_idx]
    phase_plot.plotattributes[:z] = angle.(E)
    phase_plot.plotattributes[:background_color_inside] = 0.7*[1,1,1]    
    plot!(
        sp=2,
        xlims=([-1, 1]*1e6*Lx/(2*p[:displayScaling])),
        ylims=([-1, 1]*1e6*Ly/(2*p[:displayScaling])),
        clims=(-π,π),
        aspect_ratio=:equal,
        xlabel="x [μm]",
        ylabel="y [μm]",
        title="Phase [rad]",
        legend=false        
    )
    
    plot_idx += 1
    plot!(1e6.*[-p[:Lx_main], p[:Lx_main], p[:Lx_main], -p[:Lx_main], -p[:Lx_main]]/2,1e6.*[p[:Ly_main], p[:Ly_main], -p[:Ly_main], -p[:Ly_main], p[:Ly_main]]/2,c=:red,linestyle=:dash,sp=2)
    plot_idx += 1
    display(plot!())    

    if p[:saveVideo]
        frame(video)
    end

    updatesliceindices = unique(round.(collect(LinRange(1, Nz, min(Nz, p[:updates])))))
    nextupdatesliceindicesindex = 1
    #status information
    status = Progress(Nz; dt=0.5, showspeed=true, barglyphs=BarGlyphs("[=> ]"), color=:cyan)       
    @info "Pre-planning optimal FFT..."    
    Pfft = plan_fft!(similar(E); flags=FFTW.PATIENT, timelimit=5)
    Pifft = inv(Pfft)
    @info "Pre-planning FFT done"
    @info "Starting FFT-BPM iteration..."
    for zidx=1:Nz     
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

        if zidx == updatesliceindices[nextupdatesliceindicesindex]
            Efield_plot.plotattributes[:z] = abs.(E).^2
            phase_plot.plotattributes[:z] = angle.(E)

            nextupdatesliceindicesindex += 1
            display(plot!())
            sleep(0.01)

            if p[:saveVideo]
                frame(video)
            end

            ProgressMeter.next!(status)
        end
    end

    if p[:saveVideo]
        if p[:finalizeVideo]
            mp4(video, p[:videoName]; fps=20, loop=0)
        else
            p[:videoHandle] = video
        end
    end    
    p[:E] = EField(
        E,
        Lx,
        Ly
    )
    
    p[:x] = x;
    p[:y] = y;

    @info "FFT-BPM iteration completed"
end