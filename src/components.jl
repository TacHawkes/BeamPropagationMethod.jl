
struct Simple2x2FiberCoupler{T} <: WaveguideComponent
    r_cores::T
    n_core::Complex{T}
    n_background::Complex{T}
    d_cores_initial::T
    d_cores_final::T    
    z_mid::T
    d_mid::T
end

function (c::Simple2x2FiberCoupler{T})(p::PropagationParameters, iz) where T
    Ny, Nx = size(p.n_in)     
    z = p.dz*(iz-1)            
    if z < 2*c.z_mid
        # beginning section
        z_2 = z
        x_core_0 = 0.5*(c.d_cores_initial - (c.d_cores_initial-c.d_cores_final)*sigmoid(6/c.z_mid*(z-c.z_mid)))
    elseif z < 2*c.z_mid + c.d_mid
        # middle section
        x_core_0 = 0.5*c.d_cores_final        
    else
        # end splitting section
        z_2 = z - 2*c.z_mid - c.d_mid        
        x_core_0 = 0.5*(c.d_cores_initial - (c.d_cores_initial-c.d_cores_final)*(1-sigmoid(6/c.z_mid*(z_2-c.z_mid))))
    end
    x_core_1 = -x_core_0    
    calc_2x2_coupler!(p.n_in, p.dx, p.dy, Nx, Ny, c.n_core, c.n_background, [x_core_0, x_core_1], c.r_cores)
end

function calc_2x2_coupler!(M, dx, dy, Nx, Ny, n_core, n_background, x_cores, r_core)
    Threads.@threads for ix=1:Nx
        @inbounds for iy=1:Ny
            x = dx*((ix-1) - (Nx - 1)/2.0f0)        
            y = dy*((iy-1) - (Ny - 1)/2.0f0)
            if (((x-x_cores[1])^2 + y^2) < r_core^2 || ((x-x_cores[2])^2 + y^2) < r_core^2)
                M[iy,ix] = n_core
            else
                M[iy,ix] = n_background
            end
        end
    end
end

function (c::Simple2x2FiberCoupler{T})(Lx, Ly, Nx, Ny) where T    
    x_core_0 = 0.5*c.d_cores_initial    
    x_core_1 = -x_core_0    
    n_out = RefractiveIndex(
        zeros(ComplexF32, Ny, Nx),
        Lx,
        Ly
    )
    
    calc_2x2_coupler!(n_out.n, Lx/Nx, Ly/Ny, Nx, Ny, c.n_core, c.n_background, [x_core_0, x_core_1], c.r_cores)

    n_out
end  

propagation_length(c::Simple2x2FiberCoupler{T}) where T = 4*c.z_mid + c.d_mid

struct SimpleMachZehnderModulator{T} <: WaveguideComponent
    d_0::T
    r_core::T
    n_core::Complex{T}
    n_background::Complex{T}
    d_cores::T
    z_mid::T
    phase_mod::T    
    n_offset::T
    λ::T
    d_mid::T

    SimpleMachZehnderModulator{T}(
        d_0,
        r_core,
        n_core,
        n_background,
        d_cores,
        z_mid,
        phase_mod,
        n_offset,
        λ) where T =  new(d_0, r_core, n_core, n_background, d_cores, z_mid, phase_mod, n_offset, λ, phase_mod * λ / 2 / π / n_offset)    
end

function (c::SimpleMachZehnderModulator{T})(p::PropagationParameters, iz) where T
    Ny, Nx = size(p.n_in)     
    z = p.dz*(iz-1)
    n_cores = [c.n_core, c.n_core]      
    if z < c.d_0
        x_core_0 = 0.0
    elseif z < c.d_0 + 2*c.z_mid
        # splitting section
        z_2 = z - c.d_0
        x_core_0 = 0.5*(c.d_cores*sigmoid(6/c.z_mid*(z_2-c.z_mid)))
    elseif z < c.d_0 + 2*c.z_mid + c.d_mid
        # middle section
        x_core_0 = 0.5*c.d_cores
        n_cores[1] = c.n_core + c.n_offset
    elseif z < c.d_0 + 4*c.z_mid + c.d_mid
        # combinging section
        z_2 = z - 2*c.z_mid - c.d_mid - c.d_0      
        x_core_0 = 0.5*(c.d_cores*(1-sigmoid(6/c.z_mid*(z_2-c.z_mid))))
    else
        x_core_0 = 0.0
    end
    x_core_1 = -x_core_0    
    calc_simple_mzi!(p.n_in, p.dx, p.dy, Nx, Ny, n_cores, c.n_background, [x_core_0, x_core_1], c.r_core)
end

function calc_simple_mzi!(M, dx, dy, Nx, Ny, n_cores, n_background, x_cores, r_core)    
    Threads.@threads for ix=1:Nx
        @inbounds for iy=1:Ny
            x = dx*((ix-1) - (Nx - 1)/2.0f0)        
            y = dy*((iy-1) - (Ny - 1)/2.0f0)
            if ((x-x_cores[1])^2 + y^2) < r_core^2
                M[iy,ix] = n_cores[1]
            elseif (x-x_cores[2])^2 + y^2 < r_core^2
                M[iy,ix] = n_cores[2]
            elseif abs(x-x_cores[1]) < 2r_core
                M[iy,ix] = n_background + (n_cores[1] - n_cores[2])
            else
                M[iy,ix] = n_background
            end
        end
    end
end

function (c::SimpleMachZehnderModulator{T})(Lx, Ly, Nx, Ny) where T    
    x_core_0 = 0.0
    x_core_1 = 0.0
    n_out = RefractiveIndex(
        zeros(ComplexF32, Ny, Nx),
        Lx,
        Ly
    )
    
    calc_simple_mzi!(n_out.n, Lx/Nx, Ly/Ny, Nx, Ny, [c.n_core, c.n_core], c.n_background, [x_core_0, x_core_1], c.r_core)

    n_out
end 

propagation_length(c::SimpleMachZehnderModulator{T}) where T = 2*c.d_0 + 4*c.z_mid + c.d_mid

struct SimpleMultiCore6PortPhotonicLantern{T} <: WaveguideComponent
    d_0::T
    d_1::T   
    r_sm_core::T  
    r_mm_core::T
    n_core::Complex{T}
    n_cladd_i::Complex{T}
    n_cladd_o::Complex{T}    
    r_pentagon::T
    taper_factor::T
    taper_length::T
    reverse::Bool   
end 

function (c::SimpleMultiCore6PortPhotonicLantern{T})(p::PropagationParameters, iz) where T
    Ny, Nx = size(p.n_in)
    z = p.dz*(iz-1)    
    
    if c.reverse
        if z < c.d_1
            shapes = [0 0 c.r_mm_core*c.taper_factor c.n_cladd_i]            
        elseif z < c.d_1 + c.taper_length
            z_2 = z - c.d_1
            shapes = pentagon_lantern_base_shape(c)            
            taper_factor = c.taper_factor + ((1-c.taper_factor)*sigmoid(12/c.taper_length*(z_2-c.taper_length/2)))
            for i=1:3                
                shapes[:,i] .*= taper_factor
            end
        else
            shapes = pentagon_lantern_base_shape(c)
        end
    else
        if z < c.d_0
            # base multi-core shape without tapering
            shapes = pentagon_lantern_base_shape(c)
        elseif z < c.d_0 + c.taper_length
            # tapered base shape
            z_2 = z - c.d_0
            shapes = pentagon_lantern_base_shape(c)            
            taper_factor = c.taper_factor + ((1-c.taper_factor)*(1-sigmoid(12/c.taper_length*(z_2-c.taper_length/2))))
            for i=1:3                
                shapes[:,i] .*= taper_factor
            end
        else
            # shapes is now only the multimode core, no more sm cores
            shapes = [0 0 c.r_mm_core*c.taper_factor c.n_cladd_i]
        end
    end
    
    calc_lantern_shape!(p.n_in, p.dx, p.dy, Nx, Ny, c.n_cladd_o, shapes)
end

function pentagon_lantern_base_shape(c::SimpleMultiCore6PortPhotonicLantern{T}) where T
    shapes = zeros(7, 4)

    # multimode-core
    shapes[1,3] = c.r_mm_core
    shapes[1,4] = c.n_cladd_i

    # center core
    shapes[2,3] = c.r_sm_core
    shapes[2,4] = c.n_core

    # pentagon cores
    for i=3:7
        shapes[i,1] = sind((i-3)*72)*c.r_pentagon                
        shapes[i,2] = cosd((i-3)*72)c.r_pentagon
        shapes[i,3] = c.r_sm_core
        shapes[i,4] = c.n_core
    end    

    return shapes
end

function calc_lantern_shape!(M, dx, dy, Nx, Ny, n_outer_cladding, shapes)
    Threads.@threads for ix=1:Nx
        @inbounds for iy=1:Ny
            x = dx*((ix-1) - (Nx - 1)/2.0f0)        
            y = dy*((iy-1) - (Ny - 1)/2.0f0)
            M[iy,ix] = n_outer_cladding
            for row in eachrow(shapes)
                if real((x - row[1])^2 + (y - row[2])^2) < real(row[3]^2)
                    M[iy,ix] = row[4]             
                end
            end            
        end
    end
end

function (c::SimpleMultiCore6PortPhotonicLantern{T})(Lx, Ly, Nx, Ny) where T        
    n_out = RefractiveIndex(
        zeros(ComplexF32, Ny, Nx),
        Lx,
        Ly
    )
    if c.reverse
        # only mm-core
        shapes = [0 0 c.r_mm_core*c.taper_factor c.n_cladd_i]
    else
        shapes = pentagon_lantern_base_shape(c)
    end
    calc_lantern_shape!(n_out.n, Lx/Nx, Ly/Ny, Nx, Ny, c.n_cladd_o, shapes)

    n_out
end

propagation_length(c::SimpleMultiCore6PortPhotonicLantern{T}) where T = c.d_0 + c.d_1 + c.taper_length