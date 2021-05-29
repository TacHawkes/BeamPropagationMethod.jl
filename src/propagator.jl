function substep1a(p)
    Threads.@threads for iy=1:p.Ny
        for ix=1:p.Nx
            p.E2[iy, ix] = p.E1[iy, ix]
            if ix != 1
                p.E2[iy, ix] += (p.E1[iy, ix-1] - p.E1[iy, ix])*p.ax
            end
            if ix != p.Nx
                p.E2[iy,ix] += (p.E1[iy,ix+1] - p.E1[iy,ix])*p.ax
            end
            if iy != 1
                p.E2[iy,ix] += (p.E1[iy-1,ix] - p.E1[iy,ix])*p.ay*2.0f0
            end
            if iy != p.Ny
                p.E2[iy,ix] += (p.E1[iy+1,ix] - p.E1[iy,ix])*p.ay*2.0f0
            end
        end
    end
end

function substep1b(p)        
    Threads.@threads for iy=1:p.Ny
        for ix=1:p.Nx
            ib = ix + (Threads.threadid()-1)*p.Nx

            p.b[ib] = 1
            if (ix < p.Nx) 
                p.b[ib] += p.ax
            end
            if (ix > 1)
                p.b[ib] += p.ax
                w = -p.ax/p.b[ib-1]
                p.b[ib] += w*p.ax                
                p.E2[iy, ix] -= w*p.E2[iy,ix - 1]
            end
        end

        for ix=p.Nx:-1:1
            ib = ix + (Threads.threadid()-1)*p.Nx            
            p.E2[iy,ix] = (p.E2[iy,ix] + (ix == p.Nx ? 0 : p.ax*p.E2[iy,ix+1]))/p.b[ib]
        end
    end
end


function substep2a(p)
    Threads.@threads for iy=1:p.Ny
        for ix=1:p.Nx            
            if (iy != 1)
                p.E2[iy, ix] -= (p.E1[iy-1,ix] - p.E1[iy,ix])*p.ay
            end
            if (iy != p.Ny)
                p.E2[iy,ix] -= (p.E1[iy+1,ix] - p.E1[iy,ix])*p.ay
            end
        end
    end
end

function substep2b(p)
    EfieldPowerThread = Threads.Atomic{Float64}(0.0)    
    Threads.@threads for ix=1:p.Nx
        for iy=1:p.Ny
            ib = iy + (Threads.threadid()-1)*p.Ny
            p.b[ib] = 1
            if (iy < p.Ny)
                p.b[ib] += p.ay
            end
            if (iy > 1)
                p.b[ib] += p.ay
                w = -p.ay / p.b[ib-1]
                p.b[ib] += w*p.ay                
                p.E2[iy,ix] -= w*p.E2[iy-1,ix]
            end
        end

        for iy=p.Ny:-1:1
            ib = iy + (Threads.threadid()-1)*p.Ny            
            p.E2[iy,ix] = (p.E2[iy,ix] + (iy == p.Ny ? 0 : p.ay*p.E2[iy+1,ix]))/p.b[ib]
            Threads.atomic_add!(EfieldPowerThread, convert(Float64, real(p.E2[iy,ix])^2 + imag(p.E2[iy,ix])^2))
        end
    end
    p.EfieldPower += EfieldPowerThread[]
end

function calcShapexyr(p, iz)
    cosvalue = cos(p.twistPerStep*iz)
    sinvalue = sin(p.twistPerStep*iz)
    scaling = 1 - p.taperPerStep*iz
    for iShape=1:p.Nshapes
        p.shapexs_transformed[iShape] = scaling*(cosvalue*p.shapexs[iShape] - sinvalue*p.shapeys[iShape])
        p.shapeys_transformed[iShape] = scaling*(sinvalue*p.shapexs[iShape] + cosvalue*p.shapeys[iShape])
        p.shapeRs_transformed[iShape] =  p.shapeRs[iShape]
    end
end

function applyMultiplier(p, iz)
    precisePowerDiffThread = Threads.Atomic{Float64}(0.0)
    fieldCorrection = √(p.precisePower/p.EfieldPower)
    if (p.taperPerStep != 0 || p.twistPerStep != 0)
        Threads.@threads for i=1:p.Nx*p.Ny
            ix = (i - 1) % p.Nx
            x = p.dx*(ix - (p.Nx - 1)/2.0f0)
            iy = (i - 1) ÷ p.Nx
            y = p.dy*(iy - (p.Ny - 1)/2.0f0)

            n = p.n_cladding
            absorption = p.claddingAbsorption
            for iShape=1:p.Nshapes
                if p.shapeTypes[iShape] == 1
                    if ((x - p.shapexs_transformed[iShape])^2 + (y - p.shapeys_transformed[iShape])^2 < (p.shapeRs_transformed[iShape])^2)
                        n = p.shapeRIs[iShape]
                        absorption = p.shapeAbsorptions[iShape]
                    end
                elseif p.shapeTypes[iShape] == 2
                    delta = max(p.dx, p.dy)
                    r_diff = √((x - p.shapexs_transformed[iShape])^2 + (y - p.shapeys_transformed[iShape])^2) - p.shapeRs_transformed[iShape] + delta/2.0f0
                    if r_diff < 0
                        n = p.shapeRIs[iShape]
                        absorption = p.shapeAbsorptions[iShape]
                    elseif r_diff < delta
                        n = r_diff/delta * (p.n_cladding - p.shapeRIs[iShape]) + p.shapeRIs[iShape]
                        absorption = r_diff/delta*(p.claddingAbsorption - p.shapeAbsorptions[iShape]) + p.shapeAbsorptions[iShape]
                    end
                elseif p.shapeTypes[iShape] == 3
                    r_ratio_sqr = ((x - p.shapexs_transformed[iShape])^2 + (y - p.shapeys_transformed[iShape])^2)/(p.shapeRs_transformed[iShape])^2
                    if r_ratio_sqr < 1
                        n = r_ratio_sqr * (p.n_cladding - p.shapeRIs[iShape]) + p.shapeRIs[iShape]
                        absorption = p.shapeAbsorptions[iShape]
                    end
                elseif p.shapeTypes[iShape] == 4
                    r_ratio_sqr = ((x - p.shapexs_transformed[iShape])^2 + (y - p.shapeys_transformed[iShape])^2)/(p.shapeRs_transformed[iShape])^2
                    r_abs = √((x - p.shapexs_transformed[iShape])^2 + (y - p.shapeys_transformed[iShape])^2)
                    if r_ratio_sqr < 1
                        n = 2*p.shapeRIs[iShape] * exp(p.shapegs[iShape]*r_abs) / (exp(2*p.shapegs[iShape]*r_abs)+1)
                        absorption = p.shapeAbsorptions[iShape]
                    end
                elseif p.shapeTypes[iShape] == 5
                    r_ratio_sqr = (y - p.shapeys_transformed[iShape])^2/(p.shapeRs_transformed[iShape])^2
                    r_abs = y - p.shapeys_transformed[iShape]
                    if r_ratio_sqr < 1
                        n = 2*p.shapeRIs[iShape] * exp(p.shapegs[iShape]*r_abs) / (exp(2*p.shapegs[iShape]*r_abs)+1)
                        absorption = p.shapeAbsorptions[iShape]
                    end
                end
            end
            n_eff = n*(1-(n^2*(x*p.cosBendDirection+y*p.sinBendDirection)/2/p.RoC*p.rho_e))*exp((x*p.cosBendDirection+y*p.sinBendDirection)/p.RoC)
            if iz == p.iz_end
                p.n_out[i] = n_eff
            end
            p.E2[i] *= fieldCorrection*p.multiplier[i]*exp(im*p.d*(n_eff^2 - p.n_0^2))*absorption
            multipliernormsqr = (real(p.multiplier[i])*absorption)^2 + (imag(p.multiplier[i])*absorption)^2
            if multipliernormsqr > 1 - 10*eps(Float32)
                multipliernormsqr = 1
            end
            Threads.atomic_add!(precisePowerDiffThread, convert(Float64, ((real(p.E2[i]))^2 + (imag(p.E2[i]))^2)*(1 - 1/multipliernormsqr)))
        end
    else
        Threads.@threads for i=1:(p.Nx*p.Ny)
            p.E2[i] *= fieldCorrection * p.multiplier[i]
            multipliernormsqr = real(p.multiplier[i])^2 + imag(p.multiplier[i])^2
            if multipliernormsqr > 1 - 10*eps(Float32)
                multipliernormsqr = 1
            end
            Threads.atomic_add!(precisePowerDiffThread, convert(Float64, ((real(p.E2[i]))^2 + (imag(p.E2[i]))^2)*(1 - 1/multipliernormsqr)))
        end
    end

    p.precisePower += precisePowerDiffThread[]
end

function swapEPointers(p, iz)
    if iz > p.iz_start
        temp = p.E1
        p.E1 = p.E2
        p.E2 = temp
    elseif (p.iz_end - p.iz_start) % 2 != 0
        p.E1 = p.E2
        p.E2 = Matrix{ComplexF32}(undef, p.Ny, p.Nx)
    else
        p.E1 = p.E2
        p.E2 = p.Efinal
    end
end

function fdbpm_propagator(E, pd)   
    multiplier = pd[:multiplier]  
    p = PropagationParameters(
        size(E)[2],
        size(E)[1],
        pd[:dx],
        pd[:dy],
        pd[:iz_start],
        pd[:iz_end],
        pd[:taperPerStep],
        pd[:twistPerStep],
        pd[:d],
        pd[:n_cladding],
        pd[:claddingAbsorption],
        pd[:n_0],
        size(pd[:shapes])[1],
        pd[:shapes][:,1],
        pd[:shapes][:,2],
        pd[:shapes][:,3],
        pd[:shapes][:,4],
        pd[:shapes][:,5],
        pd[:shapes][:,6],
        pd[:shapeAbsorptions],
        zeros(Float32, size(pd[:shapes])[1]),
        zeros(Float32, size(pd[:shapes])[1]),
        zeros(Float32, size(pd[:shapes])[1]),
        similar(E),
        E,
        similar(E),
        similar(E),
        zeros(Float32, size(E)),
        zeros(ComplexF32, Threads.nthreads()*max(size(E)[1],size(E)[2])),
        similar(E),
        pd[:ax],
        pd[:ay],
        pd[:rho_e],
        pd[:RoC],
        sin(pd[:bendDirection]/180*π),
        cos(pd[:bendDirection]/180*π),
        pd[:inputPrecisePower],
        0.0
    )
    
    if (p.taperPerStep != 0 || p.twistPerStep != 0)
        for i=1:(p.Nx*p.Ny)
            p.multiplier[i] = multiplier[i]
        end
    else
        for ix=1:p.Nx
            x = p.dx * ((ix-1) - (p.Nx-1)/2.0f0)
            for iy=1:p.Ny
                y = p.dy * ((iy - 1) - (p.Ny-1)/2.0f0)
                i = (ix-1)*p.Ny + iy
                n = p.n_cladding
                absorption = p.claddingAbsorption
                for iShape=1:p.Nshapes
                    if p.shapeTypes[iShape] == 1
                        # step-index disk
                        if ((x - p.shapexs[iShape])^2 + (y - p.shapeys[iShape])^2 < (p.shapeRs[iShape])^2)
                            n = p.shapeRIs[iShape]
                            absorption = p.shapeAbsorptions[iShape]
                        end
                    elseif p.shapeTypes[iShape] == 2
                        # anti-aliased step-index disk
                        delta = max(p.dx, p.dy)
                        r_diff = √((x - p.shapexs[iShape])^2 + (y - p.shapeys[iShape])^2) - p.shapeRs[iShape] + delta / 2.0f0
                        if r_diff < 0
                            n = p.shapeRIs[iShape]
                            absorption = p.shapeAbsorptions[iShape]
                        elseif (r_diff < delta)
                            n = r_diff / delta * (p.n_cladding - p.shapeRIs[iShape]) + p.shapeRIs[iShape]
                            absorption = r_diff / delta * (p.claddingAbsorption - p.shapeAbsorptions[iShape]) + p.shapeAbsorptions[iShape]
                        end
                    elseif p.shapeTypes[iShape] == 3
                        # parabolic graded index disk
                        r_ratio_sqr = ((x - p.shapexs[iShape])^2 + (y - p.shapeys[iShape])^2)/(p.shapeRs[iShape])^2
                        if r_ratio_sqr < 1
                            n = r_ratio_sqr*(p.n_cladding - p.shapeRIs[iShape]) + p.shapeRIs[iShape]
                            absorption = p.shapeAbsorptions[iShape]
                        end
                    elseif p.shapeTypes[iShape] == 4
                        # 2D hyberbolic GRIN lens
                        r_ratio_sqr = ((x - p.shapexs[iShape])^2 + (y - p.shapeys[iShape])^2)/(p.shapeRs[iShape])^2
                        r_abs = √((x - p.shapexs[iShape])^2 + (y - p.shapeys[iShape])^2)
                        if r_ratio_sqr < 1
                            n = 2 * p.shapeRIs[iShape] * exp(p.shapegs[iShape]*r_abs) / (exp(2*p.shapegs[iShape]*r_abs) + 1)
                            absorption = p.shapeAbsorptions[iShape]
                        end
                    elseif p.shapeTypes[iShape] == 5
                        # 1D (y) hyperbolic GRIN lens
                        r_ratio_sqr = (y - p.shapeys[iShape])^2 / (p.shapeRs[iShape])^2
                        r_abs = y - p.shapeys[iShape]
                        if r_ratio_sqr < 1
                            n = 2 * p.shapeRIs[iShape] * exp(p.shapegs[iShape]*r_abs) / (exp(2*p.shapegs[iShape]*r_abs) + 1)
                            absorption = p.shapeAbsorptions[iShape]
                        end
                    end
                end
                n_eff = n*(1-(n^2*(x*p.cosBendDirection + y*p.sinBendDirection)/2/p.RoC*p.rho_e))*exp((x*p.cosBendDirection+y*p.sinBendDirection)/p.RoC)
                p.n_out[i] = n_eff                
                p.multiplier[i] = multiplier[i]*exp(im*p.d*(n_eff^2 - p.n_0^2))*absorption                
            end
        end
    end

    p.EfieldPower = 0    
    for iz=p.iz_start:p.iz_end        
        substep1a(p)
        substep1b(p)
        substep2a(p)
        substep2b(p)
        if (p.taperPerStep != 0 || p.twistPerStep != 0)
            calcShapexyr(p, iz)
        end
        applyMultiplier(p, iz)        
        p.EfieldPower = 0

        if (iz + 1 < p.iz_end)
            swapEPointers(p, iz)
        end
    end    
    return p.Efinal, p.n_out, p.precisePower
end