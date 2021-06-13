"""
    fdbpm!(p)

Performs an interative finite-differnece beam propagation simulation of the optical environment specified
in `p` using the Douglas-Gunn alternating direction implicit algorithm (DG-ADI). The `p` parameter is a Dict
containing the fields specified in TODO.
"""
function fdbpm!(p)
    k_0 = 2π/p[:lambda]    

    # check fields
    get!(p, :name, @__FILE__)
    if haskey(p, :n_cladding)
        error("Error: n_cladding has been renamed n_background")
    end
    if (haskey(p, :nFunc) && haskey(p, :shapes))
        error("You must specify exactly one of the fields 'shapes' and 'nFunc'")
    end
    if haskey(p, :shapes) && isempty(p[:shapes])
        p[:shapes] = [0 0 0 1 0]
    end
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
    get!(p, :nParameters, Dict())
    get!(p, :taperScaling, 1)
    get!(p, :twistRate, 0)
    get!(p, :rho_e, 0.22)
    get!(p, :bendingRoC, Inf)
    get!(p, :bendDirection, 0)
    if haskey(p, :shapes) && size(p[:shapes], 2) == 5
        if any(p[:shapes][:,4] .== 4) || any(p[:shapes][:,4] .== 5)
            @error("Since you have a GRIN lens, you must define the gradient constant g in the shapes array")
        else
            new_shapes = fill(NaN, size(p[:shapes], 1), 6)
            new_shapes[:,1:5] = p[:shapes]
            p[:shapes] = new_shapes
        end 
    end
    if p[:saveVideo] && !haskey(p, :videoName)
        p[:videoName] = p[:name] * ".mp4"
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
        if (haskey(p, :shapes))
            p[:originalShapesInput] = p[:shapes]
        end
    end

    if isa(p[:E], Function)
        E = p[:E](X, Y, p[:Eparameters])
        E ./= sqrt(sum(abs.(E).^2))
    else
        Ny_source, Nx_source = size(p[:E].field)
        if (Nx_source != Nx || Ny_source != Ny)        
            dx_source = p[:E].Lx/Nx_source
            dy_source = p[:E].Ly/Ny_source
            x_source = dx_source*(-(Nx_source - 1)/2:(Nx_source-1)/2)
            y_source = dy_source*(-(Ny_source - 1)/2:(Ny_source-1)/2)        
            E = Interpolations.LinearInterpolation((y_source, x_source), p[:E].field; extrapolation_bc=Line())
            E = E.(X, Y)
            E .*= √(sum(abs.(p[:E].field).^2)/sum(abs.(E).^2))
        else
            E = p[:E].field
        end
    end

    E = convert.(ComplexF32, E)

    if !priorData
        p[:Einitial] = E
        p[:I0] = maximum(abs.(E).^2)
    end

    if haskey(p, :shapes)
        n = p[:n_background] * ones(Float32, Nx, Ny)
        for iShape = 1:size(p[:shapes], 1)
            if p[:shapes][iShape,4] == 1
                n[@. ((X - p[:shapes][iShape,1])^2 + (Y - p[:shapes][iShape,2])^2 < p[:shapes][iShape,3]^2)] .= p[:shapes][iShape,5]
            elseif p[:shapes][iShape,4] == 2
                delta = max(dx, dy)                              
                r_diff = @. √((X - p[:shapes][iShape,1])^2 + (Y - p[:shapes][iShape,2])^2) - p[:shapes][iShape,3] + delta/2                
                n[r_diff .< 0] .= p[:shapes][iShape,5]
                n[(r_diff .>= 0) .& (r_diff .< delta)] = @. r_diff[(r_diff >= 0) & (r_diff < delta)]/delta*(p[:n_background] - p[:shapes][iShape,5]) + p[:shapes][iShape,5]
            elseif p[:shapes][iShape,4] == 3
                r_ratio_sqr = @. ((X - p[:shapes][iShape,1])^2 + (Y - p[:shapes][iShape,2])^2)/p[:shapes][iShape,3]^2
                n[r_ratio_sqr .< 1] = r_ratio_sqr[r_ratio_sqr .< 1] *(p[:n_background] - p[:shapes][iShape,5]) .+ p[:shapes][iShape,5]
            elseif p[:shapes][iShape,4] == 4
                r = @. √((X - p[:shapes][iShape,1])^2 + (Y - p[:shapes][iShape,2])^2)
                n[r .< p[:shapes][iShape,3]] = p[:shapes][iShape,5]*sech.(p[:shapes][iShape,6]*r[r .< p[:shapes][iShape,3]])
            elseif p[:shapes][iShape,4] == 5
                r = @. abs(Y - p[:shapes][iShape,2])
                n[r .< p[:shapes][iShape,3]] = p[:shapes][iShape,5]*sech.(p[:shapes][iShape,6]*r[r .< p[:shapes][iShape,3]])
            end
        end
    elseif isa(p[:n], Function)
        n = convert.(Float32, p[:n](X,Y,p[:n_background],p[:nParameters]))
    else
        Ny_source, Nx_source = size(p[:n].n)        
        dx_source = p[:n].Lx/Nx_source
        dy_source = p[:n].Ly/Ny_source
        x_source = dx_source*(-(Nx_source - 1)/2:(Nx_source-1)/2)
        y_source = dy_source*(-(Ny_source - 1)/2:(Ny_source-1)/2)       
        itp = Interpolations.LinearInterpolation((y_source, x_source), p[:n].n; extrapolation_bc=Flat())        
        n = convert.(Float32, itp.(Y, X))
    end        
    
    Nz = max(p[:updates], round(Int64, p[:Lz] / p[:dz_target]))
    dz = p[:Lz] / Nz

    if !p[:disableStepsizeWarning]
        max_a = 5
        max_d = 2.5
        dz_max1 = max_a*4*dx^2*k_0*p[:n_0]
        dz_max2 = max_a*4*dy^2*k_0*p[:n_0]
        dz_max3 = max_d*2*p[:n_0]/k_0
        if haskey(p, :shapes) && (any(p[:shapes][:,4] .== 1) || any(p[:shapes][:,4] .== 2))
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
    d = -dz*k_0

    multiplier = @. exp(-dz*max(0, max(abs(Y) - p[:Ly_main]/2, abs(X) - p[:Lx_main]/2))^2*p[:alpha])    
    
    # init plot
    plots = plot(
        layout=(2,2),
        size=(1000, 800),
        framestyle=:box
    )
    plot_idx = 1
    if p[:downsampleImages]
        p1 = heatmap!(x_plot, y_plot, real.(n[ix_plot, iy_plot]), c=get_colormap(p[:n_colormap]))
        n_plot = plots.series_list[plot_idx]        
    else
        p1 = heatmap!(x*1e6, y*1e6, real.(n), c=get_colormap(p[:n_colormap]))
        n_plot = plots.series_list[plot_idx]        
    end        
    plot_idx += 1
    plot!(
        sp=1,
        xlims=([-1, 1]*1e6*Lx/(2*p[:displayScaling])),
        ylims=([-1, 1]*1e6*Ly/(2*p[:displayScaling])),
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
        p[:xzSlice][1][:,1] = E[round(Int64, (Ny-1)/2 + 1),:]
        p[:yzSlice][1][:,1] = E[:,round(Int64, (Nx-1)/2 + 1)]
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
    else
        heatmap!(x*1e6, y*1e6, abs.(E).^2, c=get_colormap(p[:Intensity_colormap]), sp=3)
        Efield_plot = plots.series_list[plot_idx]        
    end                
    plot_idx += 1        
    plot!(
        xlims=([-1, 1]*1e6*Lx/(2*p[:displayScaling])),
        ylims=([-1, 1]*1e6*Ly/(2*p[:displayScaling])),
        xlabel="x [μm]",
        ylabel="y [μm]",
        title="Intensity [W/m^2]",        
        legend=false,
        aspect_ratio=:equal,
        sp=3
    )    
    plot!(1e6.*[-p[:Lx_main], p[:Lx_main], p[:Lx_main], -p[:Lx_main], -p[:Lx_main]]/2,1e6.*[p[:Ly_main], p[:Ly_main], -p[:Ly_main], -p[:Ly_main], p[:Ly_main]]/2,c=:red,linestyle=:dash,sp=3)
    plot_idx += 1
    if haskey(p, :shapes) && (p[:noShapeOutline] == false && p[:taperScaling] == 1 && p[:twistRate] == 0)
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

        plot!(  xlim=(0, p[:z][end]*1e3),
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
            heatmap!(x*1e6, y*1e6, angle.(E), c=get_colormap(p[:Phase_colormap]), sp=4)
        end
        phase_plot = plots.series_list[plot_idx]        
        phase_plot.plotattributes[:background_color_inside] = 0.7*[1,1,1]    
        plot!(
            sp=4,
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
        plot!(1e6.*[-p[:Lx_main], p[:Lx_main], p[:Lx_main], -p[:Lx_main], -p[:Lx_main]]/2,1e6.*[p[:Ly_main], p[:Ly_main], -p[:Ly_main], -p[:Ly_main], p[:Ly_main]]/2,c=:red,linestyle=:dash,sp=4)
        plot_idx += 1
        if haskey(p, :shapes) && (p[:noShapeOutline] == false && p[:taperScaling] == 1 && p[:twistRate] == 0)
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
    if p[:saveVideo]
        frame(video)
    end

    if haskey(p, :nFunc)
        nFunc = p[:nFunc]
    else
        nFunc = nothing
    end
        
    pp = PropagationParameters(
        size(E)[2],
        size(E)[1],
        dx,
        dy,
        dz,
        1,
        zUpdateIdxs[1],
        Float32((1-p[:taperScaling])/Nz),
        Float32(p[:twistRate]*p[:Lz]/Nz),
        d,        
        p[:n_0],
        convert.(ComplexF32, n),  
        nFunc,       
        E,
        similar(E),
        similar(E),
        zeros(ComplexF32, size(n)),
        zeros(ComplexF32, Threads.nthreads()*max(size(E)[1],size(E)[2])),
        zeros(ComplexF32, Threads.nthreads()*max(size(E)[1],size(E)[2])),
        zeros(ComplexF32, Threads.nthreads()*max(size(E)[1],size(E)[2])),
        zeros(ComplexF32, Threads.nthreads()*max(size(E)[1],size(E)[2])),
        multiplier,
        similar(multiplier),
        ax,
        ay,
        Float32(p[:rho_e]),
        p[:bendingRoC],
        sin(p[:bendDirection]/180*π),
        cos(p[:bendDirection]/180*π),
        p[:powers][end-length(zUpdateIdxs)],
        0.0
    ) 
    
    display(plot!())    
    
    #status information
    status = Progress(length(zUpdateIdxs); dt=0.5, showspeed=true, barglyphs=BarGlyphs("[=> ]"), color=:cyan)   
    if p[:updates] > 2
        @info "Starting FD-BPM iteration..."
    end

    for updidx=1:length(zUpdateIdxs)
        if updidx > 1
            pp.iz_start = Int32(zUpdateIdxs[updidx-1])
            pp.iz_end = Int32(zUpdateIdxs[updidx])
        end

        # iterate
        fdbpm_propagator!(pp)  
        
        p[:powers][end-length(zUpdateIdxs) + updidx] = pp.precisePower

        if p[:downsampleImages]
            heatmap!(z=pp.n_out[ix_plot, iy_plot],sp=1)
        else
            n_plot.plotattributes[:z] = Surface(real.(pp.n_out))
            Efield_plot.plotattributes[:z] = Surface(abs.(pp.E2).^2)
            if (!p[:calcModeOverlaps])
                phase_plot.plotattributes[:z] = Surface(angle.(pp.E2))
            end
        end        
        
        powers_plot.plotattributes[:y] = p[:powers]
        
        p[:xzSlice][end][:,end-length(zUpdateIdxs)+updidx] = pp.E2[round(Int64,(Ny-1)/2+1),:]
        p[:yzSlice][end][:,end-length(zUpdateIdxs)+updidx] = pp.E2[:,round(Int64,(Nx-1)/2+1)]

        if p[:calcModeOverlaps]
            for iMode=1:nModes
                p[:modeOverlaps][iMode,end-length(zUpdateIdxs)+updidx] = abs(sum(pp.E2.*conj(p[:modes][iMode][:field])))^2                
                mode_plots[iMode].plotattributes[:y] = p[:modeOverlaps][iMode,:]
            end
        end        
        display(plot!())        
        
        if p[:saveVideo]
            frame(video)
        end
        
        # prepare next iteration
        swap_pointers!(pp.E1, pp.E2)              

        ProgressMeter.next!(status)
    end

    if p[:saveVideo]
        if p[:finalizeVideo]
            mp4(video, p[:videoName]; fps=20, loop=0)
        else
            p[:videoHandle] = video
        end
    end

    if haskey(p, :shapes)
        shapesFinal = p[:shapes]
        shapesFinal[:,1] = p[:taperScaling]*(cos(p[:twistRate]*p[:Lz])*p[:shapes][:,1] - sin(p[:twistRate]*p[:Lz])*p[:shapes][:,2])
        shapesFinal[:,2] = p[:taperScaling]*(sin(p[:twistRate]*p[:Lz])*p[:shapes][:,1] + cos(p[:twistRate]*p[:Lz])*p[:shapes][:,2])
        shapesFinal[:,3] = p[:taperScaling]*p[:shapes][:,3]        
        p[:shapes] = shapesFinal
    end
    
    p[:E] = EField(
        pp.E2,
        Lx,
        Ly
    )
    p[:n] = RefractiveIndex(
        pp.n_out,
        Lx,
        Ly
    )

    p[:x] = x;
    p[:y] = y;
    if p[:updates] > 2
        @info "FD-BPM iteration completed"
    end
end