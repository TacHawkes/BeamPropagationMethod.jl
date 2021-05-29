function findModes!(p, nModes, singleCoreModes, sortByLoss, plotModes)
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

    if singleCoreModes
        shapesToInclude_2Darray = (1:size(p[:shapes], 1))'
    else
        shapesToInclude_2Darray = 1:size(p[:shapes], 1)
    end

    # init plot
    if plotModes
        plots = plot(
            layout=(2,1),
            size=(1000, 800),
            framestyle=:box
        )
    end

    V = []
    D = []
    println("Finding modes...")
    for iModeFinderRun = 1:size(shapesToInclude_2Darray, 1)
        n_mat = p[:n_cladding]*ones(Nx, Ny)
        for iShape in shapesToInclude_2Darray[iModeFinderRun,:]
            if p[:shapes][iShape, 4] == 1
                n_mat[@.((X-p[:shapes][iShape,1])^2 + (Y-p[:shapes][iShape,2])^2 < p[:shapes][iShape,3]^2)] .= p[:shapes][iShape,5]
            end
            if p[:shapes][iShape, 4] == 2
                delta = max(dx, dy)
                r_diff = @. √((X-p[:shapes][iShape,1])^2 + (Y-p[:shapes][iShape,2])^2) - p[:shapes][iShape,3] + delta/2
                n_mat[r_diff .< delta] = @. min(r_diff[r_diff < delta]/delta*(p[:n_cladding] - p[:shapes][iShape,5]) + p[:shapes][iShape,5], p[:shapes][iShape,5])                
            end
            if p[:shapes][iShape, 4] == 3
                r_ratio_sqr = @. ((X-p[:shapes][iShape,1]).^2 + (Y-p[:shapes][iShape,2]).^2)/p[:shapes][iShape,3]^2
                n_mat[r_ratio_sqr .< 1] = @. r_ratio_sqr[r_ratio_sqr .< 1]*(p[:n_cladding] - p[:shapes][iShape,5]) + p[:shapes][iShape,5]
            end
            if p[:shapes][iShape, 4] == 4
                r_ratio_sqr = @. ((X-p[:shapes][iShape,1]).^2 + (Y-p[:shapes][iShape,2]).^2)/p[:shapes][iShape,3]^2
                r_abs = .√((X.-p[:shapes][iShape,1]).^2 .+ (Y.-p[:shapes][iShape,2]).^2)
                n_mat[r_ratio_sqr .< 1] .= @. 2*p[:shapes][iShape,5]*exp(p[:shapes][iShape,6]*r_abs[r_ratio_sqr .< 1])/(exp(2*p[:shapes][iShape,6]*r_abs[r_ratio_sqr .< 1]) + 1)
            end
            if p[:shapes][iShape, 4] == 5
                r_ratio_sqr = @. (Y-p[:shapes][iShape,2]).^2/p[:shapes][iShape,3]^2
                r_abs = Y .- p[:shapes][iShape,2]
                n_mat[r_ratio_sqr .< 1] .= @. 2*p[:shapes][iShape,5]*exp(p[:shapes][iShape,6]*r_abs[r_ratio_sqr .< 1])/(exp(2*p[:shapes][iShape,6]*r_abs[r_ratio_sqr .< 1]) + 1)
            end
        end            
        n_mat = @. n_mat.*(1 - (n_mat.^2 .*(X*cosd(p[:bendDirection]) + Y*sind(p[:bendDirection]))*p[:rho_e]/(2*p[:bendingRoC]))).*exp((X*cosd(p[:bendDirection]) + Y * sind(p[:bendDirection]))/p[:bendingRoC])
        
        delta_n_2 = real(n_mat).^2 .- p[:n_0]^2

        dz = 1e-10
        absorber = @. exp(-dz*(max(0, max(abs(Y) - p[:Ly_main]/2, abs(X) - p[:Lx_main]/2))^2 *p[:alpha] + 2π*imag(n_mat)/p[:lambda]))        
        ax = dz/(dx^2*2im*k0*p[:n_0])
        #ax = dz / 2 / dx^2
        ay = dz/(dy^2*2im*k0*p[:n_0])
        #ay = dz / 2 / dy^2

        #c = @. 2im*k0*p[:n_0] - ax - ay + k0^2*dz/2 * delta_n_2

        #display(heatmap(x,y, n_mat', aspect_ratio=1,xlims=(-5e-6,5e-6),ylims=(-5e-6,5e-6)))
        #return
        
        N = Nx*Ny
        # matrix C
        M_rhs = sparse(1:N, 1:N, absorber[1:N] + delta_n_2[1:N]*dz*k0/(2im*p[:n_0]),N,N)
                + sparse(1:N-1,2:N,vec([repeat([fill(ax,1,Nx-1) 0],1,Ny-1) fill(ax,1,Nx-1)]),N,N)
                + sparse(2:N,1:N-1,vec([repeat([fill(ax,1,Nx-1) 0],1,Ny-1) fill(ax,1,Nx-1)]),N,N)
                + sparse(1:N-Nx,1+Nx:N,ay,N,N)
                + sparse(1+Nx:N,1:N-Nx,ay,N,N)        
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
    println("Done")        
    if sortByLoss
        sortedidx = [sortperm(vec(real.(D)); rev=true)]
    else
        sortedidx = sortperm(vec(imag.(D)))
    end        
    p[:modes] = Vector{Dict}(undef, nModes)
    cropN = 10    
    for iMode=nModes:-1:1      
        p[:modes][iMode] = Dict{Symbol,Any}()  
        p[:modes][iMode][:Lx] = Lx
        p[:modes][iMode][:Ly] = Ly
        E = reshape(V[:,sortedidx[iMode]], (Nx, Ny))
        # set border of field to zeros
        #E[1:cropN,:] .= 0.0 + 0.0*im
        #E[end-cropN:end,:] .= 0.0 + 0.0*im
        #E[:,1:cropN] .= 0.0 + 0.0*im
        #E[:,end-cropN:end] .= 0.0 + 0.0*im        

        p[:modes][iMode][:field] = E.*exp(-im*angle(maximum(abs.(E))))
        p[:modes][iMode][:eigenval] = D[sortedidx[iMode]]
        if isinf(p[:bendingRoC]) && (size(p[:shapes],1) == 1 || singleCoreModes)            
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
            p[:modes][iMode][:label] = @sprintf("%s%d%d%s", "LP", l, m, parity)
        end
        if plotModes
            heatmap!(
                x*1e6, y*1e6, abs.(E').^2,
                aspect_ration=:equal,
                c=get_colormap(p[:Intensity_colormap]),
                sp=1
            )
            for iShape=1:size(p[:shapes], 1)
                if (p[:shapes][iShape,4]) <= 3
                    plot!(
                        1e6*(p[:shapes][iShape,1] .+ p[:shapes][iShape,3].*cos.(LinRange(0,2π,100))),
                        1e6*(p[:shapes][iShape,2] .+ p[:shapes][iShape,3].*sin.(LinRange(0,2π,100))),
                        color=:white,
                        linestyle=:dash,
                        sp=3
                    )
                end
            end

            # TODO: Complete and adapt to Plots.jl with a single figure
        end
    end
end