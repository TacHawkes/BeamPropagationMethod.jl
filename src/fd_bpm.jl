function fdbpm!(p)
    k_0 = 2π/p[:lambda]    

    # check fields
    get!(p, :name, @__FILE__)
    get!(p, :shapes, [0 0 0 1 0])
    get!(p, :figNum, 1)
    get!(p, :figTitle, "")
    get!(p, :finalizeVideo, false)    
    get!(p, :saveVideo, false)
    get!(p, :finalizeVideo, false)
    get!(p, :saveData, false)
    get!(p, :useGPU, false)    
    get!(p, :downsampleImages, false)
    get!(p, :displayScaling, 1)
    get!(p, :disableStepsizeWarning, false)
    get!(p, :calcModeOverlaps, false)
    get!(p, :noShapeOutline, false)
    get!(p, :modes, nothing)
    get!(p, :Eparameters, Dict())
    get!(p, :taperScaling, 1)
    get!(p, :twistRate, 0)
    get!(p, :rho_e, 0.22)
    get!(p, :bendingRoC, Inf)
    get!(p, :bendDirection, 0)
    if size(p[:shapes], 2) == 5
        if any(p[:shapes][:,4] .== 4) || any(p[:shapes][:,4] .== 5)
            @error("Since you have a GRIN lens, you must define the gradient constant g in the shapes array")
        else
            new_shapes = fill(NaN, size(p[:shapes], 1), 6)
            new_shapes[:,1:5] = p[:shapes]
            p[:shapes] = new_shapes
        end 
    end
    if p[:saveVideo] && !haskey(p, :videoName)
        p[:videoName] = p[:name] * ".gif"
    end
    get!(p, :Intensity_colormap, 1)
    get!(p, :Phase_colormap, 3)
    get!(p, :n_colormap, 1)
    get!(p, :calcModeOverlaps, false)

    if p[:saveVideo]
        if haskey(p, :videoHandle)
            video = p[:videoHandle]
        else
            video = Animation()
        end
    end

    dx = p[:Lx_main] / p[:Nx_main]
    dy = p[:Ly_main] / p[:Ny_main]

    targetLx = p[:padfactor] * p[:Lx_main]
    targetLy = p[:padfactor] * p[:Ly_main]

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

    if p[:calcModeOverlaps] && (p[:modes][1][:Lx] != Lx || p[:modes][1][:Ly] != Ly || size(p[:modes][1][:field],1) != Nx || size(p[:modes][1][:field],2) != Ny)
        @error("The pre-calculated mode fields do not match the simulation Lx, Ly, Nx or Ny")
    end
    x = dx*(-(Nx-1)/2:(Nx-1)/2)
    y = dy*(-(Ny-1)/2:(Ny-1)/2)
    X, Y = ndgrid(x, y)   

    if p[:downsampleImages]
        if Nx > 500
            ix_plot = round.(Int64, LinRange(1, Nx, 500))
        else
            ix_plot = 1:Nx
        end
        x_plot = x[ix_plot]

        if Ny > 500
            iy_plot = round.(Int64, LinRange(1, Ny, 500))
        else
            iy_plot = 1:Ny
        end
        y_plot = y[iy_plot]
    end

    priorData = haskey(p, :originalEinput)
    if !priorData
        if !isa(p[:E], Function)
            p[:E].field = p[:E].field/sqrt(sum(abs.(p[:E].field).^2))
        end
        p[:originalEinput] = p[:E]
        p[:originalShapesInput] = p[:shapes]
    end

    if isa(p[:E], Function)
        E = p[:E](X, Y, p[:Eparameters])
        E ./= sqrt(sum(abs.(E).^2))
    else
        Nx_Esource, Ny_Esource = size(p[:E].field)
        if (Nx_Esource != Nx || Ny_Esource != Ny)        
            dx_Esource = p[:E].Lx/Nx_Esource
            dy_Esource = p[:E].Ly/Ny_Esource
            x_Esource = dx_Esource*(-(Nx_Esource - 1)/2:(Nx_Esource-1)/2)
            y_Esource = dy_Esource*(-(Ny_Esource - 1)/2:(Ny_Esource-1)/2)        
            E = Interpolations.interpolate((x_Esource, y_Esource), p[:E].field, Gridded(Linear()))
            E = E.(X, Y)
            E .*= √(sum(abs.(p[:E].field).^2)/sum(abs.(E).^2))        
        else
            E = p[:E].field
        end
    end

    E = convert.(ComplexF32, E)

    if !priorData
        p[:Einitial] = E
    end

    Nz = max(p[:updates], round(Int64, p[:Lz] / p[:dz_target]))
    dz = p[:Lz] / Nz

    if !p[:disableStepsizeWarning]
        max_a = 5
        max_d = 2.5
        dz_max1 = max_a*4*dx^2*k_0*p[:n_0]
        dz_max2 = max_a*4*dy^2*k_0*p[:n_0]
        dz_max3 = max_d*2*p[:n_0]/k_0
        if any(p[:shapes][:,4] .== 1) || any(p[:shapes][:,4] .== 2)
            if dz > min(dz_max1, dz_max2, dz_max3)
                @warn(@sprintf("z step size is high (> %.1e m), which may introduce numerical artifacts. You can disable this warning by setting p[:disableStepsizeWarning] = true.", min(dz_max1, dz_max2, dz_max3)))
            end
        else
            if dz > min(dz_max1, dz_max2)
                @warn(@sprintf("z step size is high (> %.1e m), which may introduce numerical artifacts. You can disable this warning by setting p[:disableStepsizeWarning] = true.", min(dz_max1, dz_max2)))
            end
        end
    end

    zUpdateIdxs = round.(Int64, collect((1:p[:updates])/p[:updates]*Nz))
    if priorData
        p[:z] = [p[:z];dz*zUpdateIdxs .+ p[:z][end]]
    else
        p[:z] = [0;dz*zUpdateIdxs]
    end

    ax = dz/(4im*dx^2*k_0*p[:n_0])
    ay = dz/(4im*dy^2*k_0*p[:n_0])
    d = -dz*k_0/(2*p[:n_0])

    absorber = @. exp(-dz*max(0, max(abs(Y) - p[:Ly_main]/2, abs(X) - p[:Lx_main]/2))^2*p[:alpha])
    multiplier = absorber

    shapeAbsorptions = @. exp(-dz*2*π*imag(p[:shapes][:,5])/p[:lambda])
    claddingAbsorption = @. exp(-dz*2*π*imag(p[:n_cladding])/p[:lambda])

    # init plot
    plots = plot(
        layout=(2,2),
        size=(1000, 800),
        framestyle=:box
    )
    plot_idx = 1
    if p[:downsampleImages]
        p1 = heatmap!(x_plot, y_plot, zeros(min(500, Ny), min(500, Nx)), c=get_colormap(p[:n_colormap]))
        n_plot = plots.series_list[plot_idx]
        n_plot.plotattributes[:z] = zeros(min(500, Ny), min(500, Nx))
    else
        p1 = heatmap!(x*1e6, y*1e6, zeros(Ny, Nx), c=get_colormap(p[:n_colormap]))
        n_plot = plots.series_list[plot_idx]
        n_plot.plotattributes[:z] = n_plot.plotattributes[:z] = zeros(min(500, Ny), min(500, Nx))
    end    
    plot_idx += 1
    plot!(
        sp=1,
        xlims=([-1 1]*1e6*Lx/(2*p[:displayScaling])),
        ylims=([-1 1]*1e6*Ly/(2*p[:displayScaling])),
        xlabel="x [μm]",
        ylabel="y [μm]",
        title="Real part of refractive index",
        aspect_ratio=:equal
    )    
    
    if priorData
        p[:powers] = [p[:powers];vec(fill(NaN, 1, p[:updates]))]
        push!(p[:xzSlice], fill(NaN, Nx, p[:updates]))
        push!(p[:yzSlice], fill(NaN, Ny, p[:updates]))
    else
        p[:powers] = vec(fill(NaN, 1, p[:updates] + 1))
        p[:powers][1] = 1
        p[:xzSlice] = Matrix{ComplexF32}[fill(NaN, Nx, p[:updates])]
        p[:yzSlice] = Matrix{ComplexF32}[fill(NaN, Ny, p[:updates])]
        p[:xzSlice][1][:,1] = E[:, round(Int64, (Nx-1)/2 + 1)]
        p[:yzSlice][1][:,1] = E[round(Int64, (Ny-1)/2 + 1), :]
    end

    plot!(p[:z]*1e3, p[:powers], lw=2, sp=2)
    powers_plot = plots.series_list[plot_idx]
    plot_idx += 1    
    plot!(
        xlims=(0, p[:z][end]*1e3),
        ylims=(0, 1.1),
        xlabel="Propagation distance [mm]",
        ylabel="Relative power remaining",    
        legend=false,    
        sp=2
    )
    
    if p[:downsampleImages]
        heatmap!(x_plot*1e6, y_plot*1e6, abs.(E[iy_plot, ix_plot]).^2, c=get_colormap(p[:Intensity_colormap]), sp=3)
        Efield_plot = plots.series_list[plot_idx]
        Efield_plot.plotattributes[:z] = abs.(E[iy_plot, ix_plot]).^2
    else
        heatmap!(x*1e6, y*1e6, abs.(E').^2, c=get_colormap(p[:Intensity_colormap]), sp=3)
        Efield_plot = plots.series_list[plot_idx]
        Efield_plot.plotattributes[:z] = abs.(E').^2
    end            
    plot_idx += 1        
    plot!(
        xlims=([-1 1]*1e6*Lx/(2*p[:displayScaling])),
        ylims=([-1 1]*1e6*Ly/(2*p[:displayScaling])),
        xlabel="x [μm]",
        ylabel="y [μm]",
        title="Intensity [W/m^2]",
        legend=false,
        aspect_ratio=:equal,
        sp=3
    )    
    plot!(1e6.*[-p[:Lx_main], p[:Lx_main], p[:Lx_main], -p[:Lx_main], -p[:Lx_main]]/2,1e6.*[p[:Ly_main], p[:Ly_main], -p[:Ly_main], -p[:Ly_main], p[:Ly_main]]/2,c=:red,linestyle=:dash,sp=3)
    plot_idx += 1
    if p[:noShapeOutline] == false && p[:taperScaling] == 1 && p[:twistRate] == 0
        for iShape=1:size(p[:shapes])[1]
            if p[:shapes][iShape,4] <= 3
                plot!(
                    1e6*(p[:shapes][iShape,1] .+ p[:shapes][iShape,3].*cos.(LinRange(0,2π,100))),
                    1e6*(p[:shapes][iShape,2] .+ p[:shapes][iShape,3].*sin.(LinRange(0,2π,100))),
                    color=:white,
                    linestyle=:dash,
                    sp=3
                )
                plot_idx += 1
            end
        end
    end

    if haskey(p, :plotEmax)
        plot!(
            clims=(0, p[:plotEmax]*maximum(abs.(E).^2)),
            sp=3
        )    
    end
    
    if p[:calcModeOverlaps]
        nModes = length(p[:modes])
        if priorData
            p[:modeOverlaps] = [p[:modeOverlaps] fill(NaN, nModes, p[:updates])]
        else
            p[:modeOverlaps] = fill(NaN, nModes, p[:updates]+1)
            for iMode=1:nModes
                p[:modeOverlaps][iMode,1] = abs(sum(E.*conj(p[:modes][iMode][:field])))^2
            end
        end

        mode_plots = []
        for iMode=1:nModes
            plot!(
                p[:z]*1e3, p[:modeOverlaps][iMode,:],
                lw=2,                
                sp=4
            )
            push!(mode_plots, plots.series_list[plot_idx])
            plot_idx += 1
        end        

        plot!(  xlims=(0, p[:z][end]*1e3),
                ylim=(1e-4, 2),
                xlabel="Propagation distance [mm]",
                ylabel="Mode overlaps",
                yaxis=:log,
                sp=4
            )
    else
        if p[:downsampleImages]
            heatmap!(x_plot*1e6, y_plot*1e6, angle(E[ix_plot, iy_plot]), c=get_colormap(p[:Phase_colormap]), sp=4)
        else
            heatmap!(x*1e6, y*1e6, angle.(E'), c=get_colormap(p[:Phase_colormap]), sp=4)
        end
        phase_plot = plots.series_list[plot_idx]
        phase_plot.plotattributes[:z] = angle.(E')
        phase_plot.plotattributes[:background_color_inside] = 0.7*[1,1,1]    
        plot!(
            sp=4,
            clims=(-π,π),
            aspect_ratio=:equal,
            xlabel="x [μm]",
            ylabel="y [μm]",
            title="Phase [rad]",
            legend=false        
        )
        plot_idx += 1
        plot!(1e6.*[-p[:Lx_main], p[:Lx_main], p[:Lx_main], -p[:Lx_main], -p[:Lx_main]]/2,1e6.*[p[:Ly_main], p[:Ly_main], -p[:Ly_main], -p[:Ly_main], p[:Ly_main]]/2,c=:red,linestyle=:dash,sp=4)
        plot_idx += 1
        if p[:noShapeOutline] == false && p[:taperScaling] == 1 && p[:twistRate] == 0
            for iShape=1:size(p[:shapes])[1]
                if p[:shapes][iShape,4] <= 3
                    plot!(
                        1e6*(p[:shapes][iShape,1] .+ p[:shapes][iShape,3].*cos.(LinRange(0,2π,100))),
                        1e6*(p[:shapes][iShape,2] .+ p[:shapes][iShape,3].*sin.(LinRange(0,2π,100))),
                        color=:white,
                        linestyle=:dash,
                        sp=4
                    )
                    plot_idx += 1
                end
            end
        end
    end
    display(plot!())
    sleep(1.0)

    if p[:saveVideo]
        frame(video)
    end
    
    prop_par = Dict(
        :dx => convert(Float32, dx),
        :dy => convert(Float32, dy),
        :taperPerStep => Float32((1-p[:taperScaling])/Nz),
        :twistPerStep => Float32(p[:twistRate]*p[:Lz]/Nz),
        :shapes => convert.(Float32, real(p[:shapes])),
        :shapeAbsorptions => convert.(Float32, real(shapeAbsorptions)),
        :n_cladding => Float32(p[:n_cladding]),
        :claddingAbsorption => Float32(real(claddingAbsorption)),
        :multiplier => convert.(ComplexF32, multiplier),
        :d => Float32(d),
        :n_0 => Float32(p[:n_0]),
        :ax => ComplexF32(ax),
        :ay => ComplexF32(ay),        
        :RoC => p[:bendingRoC],
        :rho_e => Float32(p[:rho_e]),
        :bendDirection => Float32(p[:bendDirection]),
        :inputPrecisePower => p[:powers][end-length(zUpdateIdxs)]
    )    
    prop_par[:iz_start] = Int32(0)
    prop_par[:iz_end] = Int32(zUpdateIdxs[1])
    for updidx=1:length(zUpdateIdxs)
        if updidx > 1
            prop_par[:iz_start] = Int32(zUpdateIdxs[updidx-1])
            prop_par[:iz_end] = Int32(zUpdateIdxs[updidx])
        end
        E, n, precisePower = fdbpm_propagator(E, prop_par)
        prop_par[:inputPrecisePower] = precisePower
        p[:powers][end-length(zUpdateIdxs) + updidx] = precisePower
        
        if p[:downsampleImages]
            heatmap!(z=n[ix_plot, iy_plot],sp=1)
        else
            n_plot.plotattributes[:z] = n'
            Efield_plot.plotattributes[:z] = abs.(E').^2
            if (!p[:calcModeOverlaps])
                phase_plot.plotattributes[:z] = angle.(E')
            end
        end
        
        powers_plot.plotattributes[:y] = p[:powers]
        
        p[:xzSlice][end][:,end-length(zUpdateIdxs)+updidx] = E[:, round(Int64,(Nx-1)/2+1)]
        p[:yzSlice][end][:,end-length(zUpdateIdxs)+updidx] = E[round(Int64,(Ny-1)/2+1),:]

        if p[:calcModeOverlaps]
            for iMode=1:nModes
                p[:modeOverlaps][iMode,end-length(zUpdateIdxs)+updidx] = abs(sum(E.*conj(p[:modes][iMode][:field])))^2                
                mode_plots[iMode].plotattributes[:y] = p[:modeOverlaps][iMode,:]
            end
        end

        display(plot!())
        sleep(0.01)     
        
        if p[:saveVideo]
            frame(video)
        end
    end

    if p[:saveVideo]
        if p[:finalizeVideo]
            mp4(video, p[:videoName]; fps=20, loop=0)
        else
            p[:videoHandle] = video
        end
    end

    shapesFinal = p[:shapes]
    shapesFinal[:,1] = p[:taperScaling]*(cos(p[:twistRate]*p[:Lz])*p[:shapes][:,1] - sin(p[:twistRate]*p[:Lz])*p[:shapes][:,2])
    shapesFinal[:,2] = p[:taperScaling]*(sin(p[:twistRate]*p[:Lz])*p[:shapes][:,1] + cos(p[:twistRate]*p[:Lz])*p[:shapes][:,2])
    #shapesFinal[:,3] = p[:taperScaling]*p[:shapes][:,3]
    shapesFinal[:,3] = p[:shapes][:,3]
    p[:shapes] = shapesFinal
    
    p[:E] = EField(
        E,
        Lx,
        Ly
    )

    p[:x] = x
    p[:y] = y
end