@inline function substep1a!(p)    
    copyto!(p.E2, p.E1)

    # branch-free implementation for increased performance
    Threads.@threads for ix=2:(p.Nx-1)
        @inbounds @simd for iy=2:(p.Ny-1)            
            p.E2[iy,ix] +=  (p.E1[iy,ix-1] - p.E1[iy,ix])*p.ax +
                            (p.E1[iy,ix+1] - p.E1[iy,ix])*p.ax +
                            (p.E1[iy-1,ix] - p.E1[iy,ix])*p.ay*2.0f0 +
                            (p.E1[iy+1,ix] - p.E1[iy,ix])*p.ay*2.0f0
        end  
        # iy edge cases (iy=1 or iy=p.Ny)           
        p.E2[1,ix] +=   (p.E1[1,ix-1] - p.E1[1,ix])*p.ax +
                        (p.E1[1,ix+1] - p.E1[1,ix])*p.ax +
                        (p.E1[2,ix] - p.E1[1,ix])*p.ay*2.0f0
        p.E2[p.Ny,ix] +=    (p.E1[p.Ny,ix-1] - p.E1[p.Ny,ix])*p.ax +
                            (p.E1[p.Ny,ix+1] - p.E1[p.Ny,ix])*p.ax +
                            (p.E1[p.Ny-1,ix] - p.E1[p.Ny,ix])*p.ay*2.0f0                    
    end

    @inbounds @simd for iy=2:(p.Ny-1)
        # ix edge cases (ix=1 or ix=p.Nx)
        p.E2[iy,1] +=   (p.E1[iy,2] - p.E1[iy,1])*p.ax +
                        (p.E1[iy-1,1] - p.E1[iy,1])*p.ay*2.0f0 +
                        (p.E1[iy+1,1] - p.E1[iy,1])*p.ay*2.0f0
        p.E2[iy,p.Nx] +=    (p.E1[iy,p.Nx-1] - p.E1[iy,p.Nx])*p.ax +
                            (p.E1[iy-1,1] - p.E1[iy,1])*p.ay*2.0f0 +
                            (p.E1[iy+1,1] - p.E1[iy,1])*p.ay*2.0f0
    end

    # combined edge cases
    p.E2[1,1] += (p.E1[1,2] - p.E1[1,1])*p.ax +                
                 (p.E1[2,1] - p.E1[1,1])*p.ay*2.0f0    
    p.E2[p.Ny,p.Nx] +=  (p.E1[p.Ny,p.Nx-1] - p.E1[p.Ny,p.Nx])*p.ax +                 
                 (p.E1[p.Ny-1,p.Nx] - p.E1[p.Ny,p.Nx])*p.ay*2.0f0                 
    p.E2[p.Ny,1] +=  (p.E1[p.Ny,2] - p.E1[p.Ny,1])*p.ax +
                 (p.E1[p.Ny-1,1] - p.E1[p.Ny,1])*p.ay*2.0f0
    p.E2[1,p.Nx] +=  (p.E1[1,p.Nx-1] - p.E1[1,p.Nx])*p.ax +                                  
                 (p.E1[2,p.Nx] - p.E1[1,p.Nx])*p.ay*2.0f0              
end

@inline function substep1b!(p)    
    Threads.@threads for iy=1:p.Ny
        tid = Threads.threadid()
        @inbounds @simd for ix=2:p.Nx
            ib = ix + (tid-1)*p.Nx              
            p.E2[iy, ix] -= p.w1[ib]*p.E2[iy,ix-1]            
        end

        # edge case ix=p.Nx
        p.E2[iy,p.Nx] /= p.b1[p.Nx*tid]
        @inbounds @simd for ix=(p.Nx-1):-1:1
            p.E2[iy,ix] = (p.E2[iy,ix] + p.ax*p.E2[iy,ix+1])/p.b1[ix + (tid-1)*p.Nx]
        end        
    end
end

@inline function substep2a!(p)
    Threads.@threads for ix=1:p.Nx
        @inbounds @simd for iy=2:(p.Ny-1)
            p.E2[iy,ix] -=  (p.E1[iy-1,ix] - p.E1[iy,ix])*p.ay +
                            (p.E1[iy+1,ix] - p.E1[iy,ix])*p.ay            
        end

        # edge cases iy=1 and iy=p.Ny    
        p.E2[p.Ny,ix] -= (p.E1[p.Ny-1,ix] - p.E1[p.Ny,ix])*p.ay
        p.E2[1,ix] -= (p.E1[2,ix] - p.E1[1,ix])*p.ay        
    end    
end

@inline function substep2b!(p)      
    Threads.@threads for ix=1:p.Nx
        tid = Threads.threadid()
        @inbounds @simd for iy=2:p.Ny
            ib = iy + (tid-1)*p.Ny            
            p.E2[iy,ix] -= p.w2[ib]*p.E2[iy-1,ix]            
        end

        # edge case iy=p.Ny
        p.E2[p.Ny,ix] /= p.b2[p.Ny*tid]
        @inbounds @simd for iy=(p.Ny-1):-1:1            
            p.E2[iy,ix] = (p.E2[iy,ix] + p.ay*p.E2[iy+1,ix])/p.b2[iy + (tid-1)*p.Ny] 
        end
    end    
    p.EfieldPower += field_power(p.E2)
end

function apply_multiplier!(p, iz)    
    fieldCorrection = √(p.precisePower/p.EfieldPower)    
    cosvalue = cos(-p.twistPerStep*iz)
    sinvalue = sin(-p.twistPerStep*iz)
    scaling = 1/(1-p.taperPerStep*iz)    
    Threads.@threads for ix=1:p.Nx
        @inbounds @simd for iy=1:p.Ny
            x = p.dx*((ix-1) - (p.Nx - 1)/2.0f0)
            y = p.dy*((iy-1) - (p.Ny - 1)/2.0f0)

            n = zero(ComplexF32)
            if p.taperPerStep != 0 || p.twistPerStep != 0
                x_src = scaling*(cosvalue*x - sinvalue*y)
                y_src = scaling*(sinvalue*x + cosvalue*y)
                ix_src = min(max(0.0f0, x_src/p.dx + (p.Nx - 1)/2.0f0), (p.Nx - 1)*(1-eps(Float32)))
                iy_src = min(max(0.0f0, y_src/p.dy + (p.Ny - 1)/2.0f0), (p.Ny - 1)*(1-eps(Float32)))
                ix_low = floor(Int, ix_src) + 1
                iy_low = floor(Int, iy_src) + 1
                ix_frac = ix_src - floor(ix_src)
                iy_frac = iy_src - floor(iy_src)                
                n = p.n_in[iy_low, ix_low]*(1-ix_frac)*(1-iy_frac) +   
                    p.n_in[iy_low, ix_low+1]*(ix_frac)*(1-iy_frac) +
                    p.n_in[iy_low+1,ix_low]*(1-ix_frac)*(iy_frac) +
                    p.n_in[iy_low+1,ix_low+1]*(ix_frac)*(iy_frac)
            else
                n = p.n_in[iy, ix]
            end
            if iz == p.iz_end
                p.n_out[iy,ix] = n
            end                    
            n_eff = real(n)*(1-(real(n)^2*(x*p.cosBendDirection+y*p.sinBendDirection)/2/p.RoC*p.rho_e))*exp((x*p.cosBendDirection+y*p.sinBendDirection)/p.RoC)                        
            a = p.multiplier[iy,ix]*exp(p.d*(imag(n) + (n_eff^2 - p.n_0^2)*im/(2*p.n_0)))
            p.E2[iy,ix] *= fieldCorrection*a            
            anormsqr = abs2(a)
            if anormsqr > (1 - 10*eps(Float32))
                anormsqr = 1
            end
            p.amultiplier[iy,ix] = anormsqr            
        end
    end    

    Esum = 0.0
    @inbounds @simd for i in eachindex(p.E2)
        Esum += abs2(p.E2[i])*(1 - 1/p.amultiplier[i])
    end
    p.precisePower += Esum
end

function swap_field_pointers!(p, iz)
    if iz > p.iz_start
        temp = p.E1
        p.E1 = p.E2
        p.E2 = temp        
    elseif (p.iz_end - p.iz_start) % 2 != 0
        p.E1 = p.E2
        p.E2 = Matrix{ComplexF32}(undef, p.Ny, p.Nx)            
    end
end

transform_refractive_index!(n_func::WaveguideComponent, p, iz) = n_func(p, iz)
transform_refractive_index!(::Nothing, vargs...) = nothing

function initialize_bw!(p)
    for tid=1:Threads.nthreads()
        @inbounds for ix=1:p.Nx            
            ib = ix + (tid-1)*p.Nx
            p.b1[ib] = 1
            if (ix < p.Nx) 
                p.b1[ib] += p.ax
            end
            if (ix > 1)
                p.b1[ib] += p.ax
                p.w1[ib] = -p.ax/p.b1[ib-1]
                p.b1[ib] += p.w1[ib]*p.ax            
            end       
        end
    
        @inbounds for iy=1:p.Ny
            ib = iy + (tid-1)*p.Ny
            p.b2[ib] = 1
            if (iy < p.Ny)
                p.b2[ib] += p.ay
            end
            if (iy > 1)
                p.b2[ib] += p.ay
                p.w2[ib] = -p.ay / p.b2[ib-1]
                p.b2[ib] += p.w2[ib]*p.ay                                        
            end
        end
    end
end

"""
    fdbpm_propagator!(p)

Propagator function which translates the system a given number of z-steps.
The `p` parameter is of type `PropagationParameters`.
"""
function fdbpm_propagator!(p)
    # initial work
    initialize_bw!(p)

    for iz=p.iz_start:p.iz_end    
        # transform refractive index
        transform_refractive_index!(p.n_func, p, iz)

        # perform DG-ADI steps    
        substep1a!(p)
        substep1b!(p)
        substep2a!(p)
        substep2b!(p)
        # apply phase multiplier 
        apply_multiplier!(p, iz)        
        # reset the step power tracking
        p.EfieldPower = 0

        # swap the E1/E2 pointers for the next iteration step
        if (iz < p.iz_end)
            swap_field_pointers!(p, iz)
        end        
    end    
end