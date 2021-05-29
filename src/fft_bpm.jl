function fftbpm!(p)
       # check fields
    get!(p, :name, @__FILE__)    
    get!(p, :figNum, 1)
    get!(p, :figTitle, "")
    get!(p, :finalizeVideo, false)    
    get!(p, :saveVideo, false)
    get!(p, :finalizeVideo, false)    
    if p[:saveVideo] && !haskey(p, :videoName)
        p[:videoName] = p[:name] * ".gif"
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

    kx = 2π/Lx*[0:floor(Int64,(Nx-1)/2), floor(Int64,-(Nx-1)/2:-1)]
    ky = 2π/Ly*[0:floor(Int64,(Ny-1)/2), floor(Int64,-(Ny-1)/2:-1)]
    kX, kY = ndgrid(convert.(Float32, kx), convert.(Float32, ky))

    prop_kernel = exp(im*dz*(kX.^2 + kY.^2)*p[:lambda]/(4π*p[:n_0]))

    x = dx*(-(Nx-1)/2:(Nx-1)/2)
    y = dy*(-(Ny-1)/2:(Ny-1)/2)
    X, Y = ndgrid(convert.(Float32, x), convert.(Float32, y))

    absorber = @. exp(-dz*max(0, max(abs(Y) - p[:Ly_main]/2, abs(X) - p[:Lx_main]/2))^2*p[:alpha])

    if isa(p[:E], Function)
        E = p[:E](X, Y, p[:Eparameters])
        p[:E_0] = E        
    else
        Nx_Esource, Ny_Esource = size(p[:E].field)
        dx_Esource = p[:E].Lx/Nx_Esource
        dy_Esource = p[:E].Ly/Ny_Esource
        x_Esource = dx_Esource*(-(Nx_Esource - 1)/2:(Nx_Esource-1)/2)
        y_Esource = dy_Esource*(-(Ny_Esource - 1)/2:(Ny_Esource-1)/2)        
        E = Interpolations.interpolate((x_Esource, y_Esource), p[:E].field, Gridded(Linear()))
        E = E.(X, Y)        
    end

    if p[:saveVideo]
        if haskey(p, :videoHandle)
            video = p[:videoHandle]
        else
            video = Animation()
        end
    end
end