"""
    find_modes!(p, nModes, singleCoreModes, sortByLoss, plotModes)

A mode-solver which finds the eigenmodes of the given shape/refractive index geometry as 
given in `p`.
"""
function find_modes!(p, nModes, singleCoreModes, sortByLoss, plotModes)
    if haskey(p, :n_cladding)
        error("Error: n_cladding has been renamed n_background")
    end
    if (haskey(p, :nFunc) && haskey(p, :shapes))
        error("You must specify exactly one of the fields 'shapes' and 'nFunc'")
    end
    get!(p, :nParameters, Dict())
    get!(p, :rho_e, 0.22)
    get!(p, :bendingRoC, Inf)
    get!(p, :bendDirection, 0)
    get!(p, :Intensity_colormap, 1)
    get!(p, :Phase_colormap, 3)

    k0 = 2π/p[:lambda]
    dx = p[:Lx_main] / p[:Nx_main]
    dy = p[:Ly_main] / p[:Ny_main]

    targetLx = p[:padfactor]*p[:Lx_main]
    targetLy = p[:padfactor]*p[:Ly_main]

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

    x = dx*(-(Nx-1)/2:(Nx-1)/2)
    y = dy*(-(Ny-1)/2:(Ny-1)/2)    
    X, Y = ndgrid(x, y)        

    # init plot    
    if plotModes
        plots = plot(
            layout=nModes,
            size=(1000, 800),
            framestyle=:box
        )
    end    

    V = []
    D = []

    if !haskey(p, :shapes)
        shapesToInclude_2Darray = 1
    elseif singleCoreModes
        shapesToInclude_2Darray = (1:size(p[:shapes], 1))'
    else
        shapesToInclude_2Darray = 1:size(p[:shapes], 1)
    end

    dz = 1e-10

    @info "Finding modes..."
    for iModeFinderRun = 1:size(shapesToInclude_2Darray, 1)
        n = p[:n_background]*ones(Nx, Ny)
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
            Nx_source, Ny_source = size(p[:n].n)        
            dx_source = p[:n].Lx/Nx_source
            dy_source = p[:n].Ly/Ny_source
            x_source = dx_source*(-(Nx_source - 1)/2:(Nx_source-1)/2)
            y_source = dy_source*(-(Ny_source - 1)/2:(Ny_source-1)/2)        
            itp = Interpolations.interpolate((x_source, y_source), p[:n].n, Gridded(Linear()))
            n = convert.(Float32, itp.(X, Y))
        end 
        
        n_eff = @. real(n)*(1-(real(n)^2*(X*cosd(p[:bendDirection]) + Y*sind(p[:bendDirection]))*p[:rho_e]/(2*p[:bendingRoC]))).*exp((X*cosd(p[:bendDirection]) + Y * sind(p[:bendDirection]))/p[:bendingRoC])        
        
        delta_n_2 = @. real(n_eff)^2 - p[:n_0]^2
        
        absorber = @. exp(-dz*(max(0, max(abs(Y) - p[:Ly_main]/2, abs(X) - p[:Lx_main]/2))^2 *p[:alpha] + 2π*imag(n)/p[:lambda]))        
        # deliberate symmetry breaking to find the e/o modes aligned to the x/y grid
        ax = @. 1.00001*dz/(dx^2*2im*k0*p[:n_0])
        ay = @. dz/(dy^2*2im*k0*p[:n_0])
        
        N = Nx*Ny
        # matrix C
        M_rhs = sparse(1:N, 1:N, absorber[1:N] + delta_n_2[1:N]*dz*k0/(2im*p[:n_0]),N,N) +
                sparse(1:N-1,2:N,vec([repeat([fill(ax,1,Nx-1) 0],1,Ny-1) fill(ax,1,Nx-1)]),N,N) +
                sparse(2:N,1:N-1,vec([repeat([fill(ax,1,Nx-1) 0],1,Ny-1) fill(ax,1,Nx-1)]),N,N) +
                sparse(1:N-Nx,1+Nx:N,ay,N,N) +
                sparse(1+Nx:N,1:N-Nx,ay,N,N)        
        M_rhs[1:N+1:N*N] = M_rhs[1:N+1:N*N] - vec(repeat([ax fill(2*ax,1,Nx-2) ax],1,Ny))
        M_rhs[1:N+1:N*N] = M_rhs[1:N+1:N*N] - vec([fill(ay,1,Nx) fill(2*ay,1,Nx*(Ny-2)) fill(ay,1,Nx)])
        absorberdiag = sparse(1:N,1:N,absorber[1:N],N,N)
        M_rhs = M_rhs*absorberdiag
        
        Drun, Vrun = eigs(M_rhs; nev=ceil(Int, nModes/size(shapesToInclude_2Darray,1)), sigma=1.0, ncv=nModes*10)
        if iModeFinderRun == 1
            V = Vrun
            D = Drun
        else
            V = [V Vrun]
            D = [D Drun]
        end        
    end    
    @info "Done"
    if sortByLoss
        sortedidx = sortperm(vec(real.(D)); rev=true)
    else
        sortedidx = sortperm(vec(imag.(D)))
    end        
    p[:modes] = Vector{Dict}(undef, nModes)
    
    for iMode=1:nModes
        p[:modes][iMode] = Dict{Symbol,Any}()  
        p[:modes][iMode][:Lx] = Lx
        p[:modes][iMode][:Ly] = Ly        
        E = reshape(V[:,sortedidx[iMode]], (Nx, Ny))
        p[:modes][iMode][:field] = E.*exp(-im*angle(maximum(abs.(E))))
        p[:modes][iMode][:eigenval] = D[sortedidx[iMode]]
        if haskey(p, :shapes) && (isinf(p[:bendingRoC]) && (size(p[:shapes],1) == 1 || singleCoreModes))
            iMax = argmax(abs.(E))
            xMax = X[iMax]
            yMax = Y[iMax]
            theta = atan(yMax -p[:shapes][1,2], xMax - p[:shapes][1,1])
            itp = Interpolations.LinearInterpolation((x, y), E'; extrapolation_bc=Line())
            radialE = itp.(p[:shapes][1,1] .+ LinRange(0, max(Lx, Ly), 1000)*cos(theta), p[:shapes][1,2] .+ LinRange(0, max(Lx,Ly), 1000)*sin(theta))
            radialEpruned = radialE[abs.(radialE) .> 0.01*maximum(abs.(radialE))]
            m = sum(abs.(diff(angle.(radialEpruned) .> 0))) + 1

            R = √((xMax - p[:shapes][1,1])^2 + (yMax - p[:shapes][1,2])^2)            
            azimuthalE = itp.(p[:shapes][1,1] .+ R*cos.(theta .+ LinRange(0,2π,1000)),p[:shapes][1,2] .+ R*sin.(theta .+ LinRange(0, 2π, 1000)))
            azimuthalEpruned = azimuthalE[abs.(azimuthalE) .> 0.01*maximum(abs.(azimuthalE))]
            l = sum(abs.(diff(angle.(azimuthalEpruned) .> 0))) / 2

            if l > 0
                Emaxmirrored = Interpolations.interpolate((x, y), E', Gridded(Linear()))
                Emaxmirrored = Emaxmirrored(xMax, 2*p[:shapes][1,2] - yMax)
                if (real(E[iMax]/Emaxmirrored)) < 0
                    parity = "o"
                else
                    parity = "e"
                end
            else
                parity = ""
            end
            p[:modes][iMode][:label] = @sprintf("Mode %d, %s%d%d%s", iMode, "LP", l, m, parity)
        else
            p[:modes][iMode][:label] = @sprintf("Mode %d", iMode)
        end
        if plotModes
            mode_loss = -10*log10(real(D[sortedidx[iMode]]))/dz
            heatmap!(
                x*1e6, y*1e6, abs.(E').^2,
                aspect_ration=:equal,
                c=get_colormap(p[:Intensity_colormap]),
                colorbar=false,
                aspect_ratio=:equal,
                sp=iMode,
                title=p[:modes][iMode][:label],
                xlims=([-1, 1]*1e6*Lx/(2)),
                ylims=([-1, 1]*1e6*Ly/(2)),
            )
            annotate!((minimum(x)*1e6, maximum(y)*1e6, Plots.text(@sprintf("%.2f dB/m loss", mode_loss), 10, :white, :top, :left)), sp=iMode)
            annotate!((minimum(x)*1e6, minimum(y)*1e6, Plots.text(@sprintf("eigenvalue-1=%.3e", (real(D[sortedidx[iMode]]) - 1)), 10, :white, :bottom, :left)), sp=iMode)
            for iShape=1:size(p[:shapes], 1)
                if (p[:shapes][iShape,4]) <= 3
                    plot!(
                        1e6*(p[:shapes][iShape,1] .+ p[:shapes][iShape,3].*cos.(LinRange(0,2π,100))),
                        1e6*(p[:shapes][iShape,2] .+ p[:shapes][iShape,3].*sin.(LinRange(0,2π,100))),
                        color=:white,
                        linestyle=:dash,
                        legend=false,
                        sp=iMode
                    )
                end
            end            
            # TODO: Complete and adapt to Plots.jl with a single figure
        end
    end

    if plotModes
        display(plot!())
        sleep(0.01)
    end
end

"""
    mode_superposition(p, modeidxs, coeffs=ones(prod(size(modeidxs)), 1))

Calculates and returns a mode superposition of the pre-calulcated eigenmodes with the indices
given by `modeidxs`. The superposition coefficients (may be complex) can be passed using the `coeffs` parameter.
"""
function mode_superposition(p, modeidxs, coeffs=ones(prod(size(modeidxs)), 1))
    out = EField(
        zeros(ComplexF32, size(p[:modes][1][:field])),
        p[:modes][modeidxs[1]][:Lx],
        p[:modes][modeidxs[1]][:Ly],
    )    
    for modeIdx=1:prod(size(modeidxs))
        out.field += coeffs[modeIdx]*p[:modes][modeidxs[modeIdx]][:field]
    end

    return out
end