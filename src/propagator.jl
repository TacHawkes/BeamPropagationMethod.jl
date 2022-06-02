function substep1a!(p::Parameters)
    copyto!(p.E2, p.E1)
    x_is_anti_symm = p.xsymmetry == AntiSymmetric
    y_is_anti_symm = p.ysymmetry == AntiSymmetric

    # branch-free implementation for increased performance
    Threads.@threads for iy=2:(p.ny-1)
        @inbounds @simd for ix=2:(p.nx-1)
            p.E2[ix,iy] +=  (p.E1[ix-1,iy] - p.E1[ix,iy])*p.ax +
                            (p.E1[ix+1,iy] - p.E1[ix,iy])*p.ax +
                            (p.E1[ix,iy-1] - p.E1[ix,iy])*p.ay*2.0f0 +
                            (p.E1[ix,iy+1] - p.E1[ix,iy])*p.ay*2.0f0
        end
        # ix edge cases (ix=1 or ix=p.Nx)
        p.E2[1,iy] +=   ifelse(!y_is_anti_symm, (p.E1[2,iy] - p.E1[1,iy])*p.ax, zero(eltype(p.E1))) +
                        (p.E1[1,iy-1] - p.E1[1,iy])*p.ay*2.0f0 +
                        (p.E1[1,iy+1] - p.E1[1,iy])*p.ay*2.0f0

        p.E2[p.nx,iy] +=    (p.E1[p.nx-1,iy] - p.E1[p.nx,iy])*p.ax +
                            (p.E1[1,iy-1] - p.E1[1,iy])*p.ay*2.0f0 +
                            (p.E1[1,iy+1] - p.E1[1,iy])*p.ay*2.0f0
    end

    @inbounds @simd for ix=2:(p.nx-1)
        # iy edge cases (iy=1 or iy=p.Ny)
        p.E2[ix,1] +=   (p.E1[ix-1,1] - p.E1[ix,1])*p.ax +
                        (p.E1[ix+1,1] - p.E1[ix,1])*p.ax +
                        ifelse(!x_is_anti_symm, (p.E1[ix,2] - p.E1[ix,1])*p.ay*2.0f0, zero(eltype(p.E1)))
        p.E2[ix,p.ny] +=    (p.E1[ix-1,p.ny] - p.E1[ix,p.ny])*p.ax +
                            (p.E1[ix+1,p.ny] - p.E1[ix,p.ny])*p.ax +
                            (p.E1[ix,p.ny-1] - p.E1[ix,p.ny])*p.ay*2.0f0
    end

    # combined edge cases
    p.E2[1,1] += ifelse(!y_is_anti_symm, (p.E1[2,1] - p.E1[1,1])*p.ax, zero(eltype(p.E1))) +
                 ifelse(!x_is_anti_symm, (p.E1[1,2] - p.E1[1,1])*p.ay*2.0f0, zero(eltype(p.E1)))
    p.E2[p.nx,p.ny] +=  (p.E1[p.nx-1,p.ny] - p.E1[p.nx,p.ny])*p.ax +
                 (p.E1[p.nx,p.ny-1] - p.E1[p.ny,p.ny])*p.ay*2.0f0
    p.E2[1,p.ny] += ifelse(!y_is_anti_symm, (p.E1[2,p.ny] - p.E1[1,p.ny])*p.ax, zero(eltype(p.E1))) +
                 (p.E1[1,p.ny-1] - p.E1[1,p.ny])*p.ay*2.0f0
    p.E2[p.nx,1] +=  (p.E1[p.nx-1,1] - p.E1[p.nx,1])*p.ax +
                 ifelse(!x_is_anti_symm, (p.E1[p.nx,2] - p.E1[p.nx,1])*p.ay*2.0f0, zero(eltype(p.E1)))
end

function substep1b!(p::Parameters)
    y_is_anti_symm = p.ysymmetry == AntiSymmetric
    Threads.@threads for iy=1:p.ny
        @inbounds @simd for ix=2:p.nx
            p.E2[ix, iy] -= p.w1[ix]*p.E2[ix-1,iy]
        end

        # edge case ix=p.Nx
        p.E2[p.nx,iy] /= p.b1[p.nx]
        @inbounds @simd for ix=(p.nx-1):-1:(1 + Int(y_is_anti_symm))
            p.E2[ix,iy] = (p.E2[ix,iy] + p.ax*p.E2[ix+1,iy])/p.b1[ix]
        end
    end
end

function substep2a!(p::Parameters)
    x_is_anti_symm = p.xsymmetry == AntiSymmetric
    Threads.@threads for iy in 2:p.ny-1
        for ix in 1:p.nx
            @inbounds p.E2[ix,iy] -= (p.E1[ix,iy-1] - p.E1[ix,iy])*p.ay +
                                     (p.E1[ix,iy+1] - p.E1[ix,iy])*p.ay
        end
    end

    # edge cases iy=1 and iy=p.ny
    for ix in 1:p.nx
        @inbounds p.E2[ix,1] -= ifelse(!x_is_anti_symm, (p.E1[ix,2] - p.E1[ix,1])*p.ay, zero(eltype(p.E2)))
        @inbounds p.E2[ix,p.ny] -= (p.E1[ix,p.ny-1] - p.E1[ix,p.ny])*p.ay
    end
end

function substep2b!(p::Parameters)
    x_is_anti_symm = p.xsymmetry == AntiSymmetric
    Threads.@threads for ix=1:p.nx
        for iy=2:p.ny
            @inbounds p.E2[ix,iy] -= p.w2[iy]*p.E2[ix,iy-1]
        end

        # edge case iy=p.Ny
        @inbounds p.E2[ix,p.ny] /= p.b2[p.ny]
        for iy=(p.ny-1):-1:(1 + Int(x_is_anti_symm))
            @inbounds p.E2[ix,iy] = (p.E2[ix,iy] + p.ay*p.E2[ix,iy+1])/p.b2[iy]
        end
    end
    p.E_field_power += field_power(p.E2)
end

function apply_multiplier!(p::Parameters, iz)
    field_correction = √(p.precise_power/p.E_field_power)
    cosvalue = cos(-p.twist_per_step*iz)
    sinvalue = sin(-p.twist_per_step*iz)
    scaling = 1/(1-p.taper_per_step*iz)
    xnosymmetry = Int(p.xsymmetry == NoSymmetry)
    ynosymmetry = Int(p.ysymmetry == NoSymmetry)
    amultiplier = Array{Float64}(undef, size(p.E2))
    Threads.@threads for ix=1:p.nx
        for iy=1:p.ny
            x = p.dx*((ix-1) - (p.nx - 1)/2.0f0*ynosymmetry)
            y = p.dy*((iy-1) - (p.ny - 1)/2.0f0*xnosymmetry)

            n = zero(ComplexF32)
            if p.taper_per_step != 0 || p.twist_per_step != 0
                x_src = scaling*(cosvalue*x - sinvalue*y)
                y_src = scaling*(sinvalue*x + cosvalue*y)
                ix_src = min(max(0.0f0, x_src/p.dx + (p.nx - 1)/2.0f0*ynosymmetry), (p.nx - 1)*(1-eps(Float32)))
                iy_src = min(max(0.0f0, y_src/p.dy + (p.ny - 1)/2.0f0*xnosymmetry), (p.ny - 1)*(1-eps(Float32)))
                ix_low = floor(Int, ix_src) + 1
                iy_low = floor(Int, iy_src) + 1
                ix_frac = ix_src - floor(ix_src)
                iy_frac = iy_src - floor(iy_src)
                n = p.n_in[ix_low, iy_low]*(1-ix_frac)*(1-iy_frac) +
                    p.n_in[ix_low+1,iy_low]*(ix_frac)*(1-iy_frac) +
                    p.n_in[ix_low,iy_low+1]*(1-ix_frac)*(iy_frac) +
                    p.n_in[ix_low+1,iy_low+1]*(ix_frac)*(iy_frac)
            elseif p.nz_n == 1
                # 2D RI
                n = p.n_in[ix, iy, 1]
            else
                # 3D RI
                z = p.dz*(iz-1)
                ix_n = min(max(1, ix-1 - (p.nx - p.nx_n)/2), p.nx_n)
                iy_n = min(max(1, iy-1 - (p.ny - p.ny_n)/2), p.ny_n)
                iz_n = min(max(0.0f0, z/p.dz_n), (p.nz_n - 1)*(1-eps(Float32)))
                iz_n_low = floor(Int, iz_n) + 1
                iz_n_frac = iz_n - iz_n_low
                n = p.n_in[ix_n, iy_n, iz_n_low] * (1 - iz_n_frac) +
                    p.n_in[ix_n, iy_n, iz_n_low + 1] * (iz_n_frac)
            end
            if iz == p.iz_end
                p.n_out[ix,iy] = n
            end
            n_bend = real(n)*(1-(real(n)^2*(x*p.cos_bend_direction+y*p.sin_bend_direction)/2/p.radius_of_curvature*p.ρe))*exp((x*p.cos_bend_direction+y*p.sin_bend_direction)/p.radius_of_curvature)
            a = p.multiplier[ix,iy]*exp(p.d*(imag(n) + (n_bend^2 - p.n0^2)*im/(2*p.n0)))
            p.E2[ix,iy] *= field_correction*a
            anormsqr = abs2(a)
            if anormsqr > (1 - 10*eps(Float32)) && anormsqr < 1 + 10*eps(Float32)
                anormsqr = 1
            end
            amultiplier[ix,iy] = anormsqr
        end
    end

    Esum = 0.0
    @inbounds @simd for i in eachindex(p.E2)
        Esum += abs2(p.E2[i])*(1 - 1/amultiplier[i])
    end
    p.precise_power_diff += Esum
end

function swap_field_pointers!(p::Parameters, iz)
    if iz > p.iz_start
        temp = p.E1
        p.E1 = p.E2
        p.E2 = temp
    elseif (p.iz_end - p.iz_start) % 2 != 0
        p.E1 = p.E2
        p.E2 = Matrix{ComplexF32}(undef, p.nx, p.ny)
    else
        p.E1 = p.E2
        p.E2 = p.Efinal
    end
end

function update_precise_power!(p::Parameters)
    p.precise_power += p.precise_power_diff
    p.precise_power_diff = zero(typeof(p.precise_power_diff))
end

function initialize_bw!(p)
    x_is_anti_symm = p.xsymmetry == AntiSymmetric
    y_is_anti_symm = p.ysymmetry == AntiSymmetric
    for ix=1:p.nx
        if (ix == 1 && y_is_anti_symm)
            p.b1[ix] = 1.0f0
        elseif (ix == 1)
            p.b1[ix] = 1.0f0 + p.ax
        elseif (ix < p.nx)
            p.b1[ix] = 1.0f0 + 2.0f0*p.ax
        else
            p.b1[ix] = 1.0f0 + p.ax
        end
        if (ix > 1)
            p.w1[ix] = -p.ax/p.b1[ix-1]
            p.b1[ix] += p.w1[ix]*(ix == 2 && y_is_anti_symm ? zero(eltype(p.b1)) : p.ax)
        end
    end

    for iy=1:p.ny
        if (iy == 1 && x_is_anti_symm)
            p.b2[iy] = 1.0f0
        elseif (iy == 1)
            p.b2[iy] = 1.0f0 + p.ay
        elseif (iy < p.ny)
            p.b2[iy] = 1.0f0 + 2.0f0*p.ay
        else
            p.b2[iy] = 1.0f0 + p.ay
        end
        if (iy > 1)
            p.w2[iy] = -p.ay/p.b2[iy-1]
            p.b2[iy] += p.w2[iy]*(iy == 2 && x_is_anti_symm ? zero(eltype(p.b2)) : p.ay)
        end
    end
end

function fdbpm_propagator!(p)
    # initial work
    initialize_bw!(p)

    for iz=p.iz_start:p.iz_end
        # perform DG-ADI steps
        substep1a!(p)
        substep1b!(p)
        substep2a!(p)
        substep2b!(p)
        # apply phase multiplier
        apply_multiplier!(p, iz)
        # reset the step power tracking
        p.E_field_power = zero(typeof(p.E_field_power))

        # swap the E1/E2 pointers for the next iteration step
        if (iz < p.iz_end)
            swap_field_pointers!(p, iz)
        end
        update_precise_power!(p)
    end
end
