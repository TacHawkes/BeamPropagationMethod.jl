"""
    mode_superposition(p, modeidxs, coeffs=ones(prod(size(modeidxs)), 1))

Calculates and returns a mode superposition of the pre-calulcated eigenmodes with the indices
given by `modeidxs`. The superposition coefficients (may be complex) can be passed using the `coeffs` parameter.
"""
function mode_superposition(m, mode_idxs, coeffs=ones(prod(size(mode_idxs)), 1))
    out = ElectricFieldProfile(
        lx = m.modes[mode_idxs[1]].lx,
        ly = m.modes[mode_idxs[1]].ly,
        field = m.modes[mode_idxs[1]].field,
        xsymmetry = lx = m.modes[mode_idxs[1]].xsymmetry,
        ysymmetry = lx = m.modes[mode_idxs[1]].ysymmetry,
    )

    for mode_idx in 1:length(mode_idxs)
        out.field += coeffs[mode_idx]*m.modes[mode_idx].field
    end

    return out
end

function find_cores(n, n_background)
    A = (single.(n) .!= single(n_background))

    nx, ny = size(A)

    B = zeros(nx,ny)
    n = 0
    for ix in 1:nx
        on = false
        for iy in 1:ny
            if !on && A[ix, iy]
                n += 1
                B[ix, iy] = n
                on = true
            elseif on && A[ix, iy]
                B[ix, iy] = n
            elseif on && !A[ix, iy]
                on = false
            end
        end
    end

    for iy in 1:ny
        for ix=2:nx
            if B[ix, iy] && B[ix-1,iy]
                B[B .== B[ix,iy]] .= B[ix-1,iy]
            end
        end
    end

   ic = indexin(unique(B), B)

   return reshape(ic .- 1, size(B))
end

function test_radial_symmetry(x,y,n,nbackground,xsymmetry,ysymmetry)
    n = n .- nbackground

    if ysymmetry != NoSymmetry
        x_c = 0
    else
        x_c = 0
        for j in eachindex(y)
            for i in eachindex(x)
                x_c += x[i] * abs2(n[i,j])
            end
        end
        x_c /= sum(abs.(n))
    end
    if xsymmetry != NoSymmetry
        y_c = 0
    else
        y_c = 0
        for j in eachindex(y)
            for i in eachindex(x)
                y_c += y[j] * abs2(n[i,j])
            end
        end
        y_c /= sum(abs.(n))
    end

    R_f(xi,yj) = √((xi - x_c)^2 + (yj - y_c)^2)
    R = [R_f(xi, yj) for (xi, yj) in Iterators.product(x,y)]
    sort_idxs = sortperm(R[:])
    nsorted = n[sort_idxs]

    monotonicity_real = sign.(diff(real(nsorted)))
    reversals_real = sum(abs.(diff(monotonicity_real[monotonicity_real .!= 0])/2))
    monotonicity_imag = sign.(diff(imag(nsorted)))
    reversals_imag = sum(abs.(diff(monotonicity_imag[monotonicity_imag .!= 0])/2))

    return reversals_real < 5 && reversals_imag < 5, x_c, y_c
end

function find_modes!(m::Model, n_modes; single_core_modes=false, sort_by_loss=false, plot_modes=true)
    if ((m.xsymmetry != NoSymmetry && !isinf(m.bending_radius_of_curvature) && sind(m.bend_direction) != 0.0)
        || (m.ysymmetry != NoSymmetry) && !isinf(m.bending_radius_of_curvature) && cosd(m.bend_direction) != 0.0)
        @error("The specified bending direction is inconsistent with the symmetry assumption")
    end

    _dx = dx(m)
    _dy = dy(m)
    _nx = nx(m)
    _ny = ny(m)
    _lx = lx(m)
    _ly = ly(m)
    _x = x(m)
    _y = y(m)

    N = _nx * _ny
    if n_modes >= N - 1
        @error("Error: The number of modes requested must be less than the pixels in the full simulation window minus one.")
    end

    k0 = 2π/m.λ

    @info "Finding modes"
    tstart = time()

    nx_source, ny_source, nz_source = size(m.n.n)
    dx_source = m.n.lx / nx_source
    dy_source = m.n.ly / ny_source
    x_source = get_grid_array(nx_source, dx_source, m.n.ysymmetry)
    y_source = get_grid_array(ny_source, dy_source, m.n.xsymmetry)
    x_source, y_source, n_source = calc_full_refractive_index(x_source, y_source, m.n.n)
    itp = Interpolations.LinearInterpolation((x_source, y_source), @view(n_source[:,:,1]); extrapolation_bc=m.nbackground)
    n = [itp(xi, yj) for (xi, yj) in Iterators.product(_x, _y)]

    anycomplex = !isreal(n)

    if single_core_modes
        core_idxs = find_cores(n, m.nbackground)
    else
        core_idxs = ones(size(n))
    end

    n_cores = maximum(core_idxs)
    m.modes = ElectricFieldProfile[]
    i_mode = 0
    for i_core in 1:n_cores
        n_core = n
        n_core[core_idxs .!= i_core] .= m.nbackground

        if !isfinite(m.bending_radius_of_curvature)
            radially_symmetric, x_c, y_c = test_radial_symmetry(
                _x,
                _y,
                n_core,
                m.nbackground,
                m.xsymmetry,
                m.ysymmetry
                )
        else
            radially_symmetric = false
        end

        n_core_bent = Matrix{Float64}(undef, length(_x), length(_y))
        for j in eachindex(_y)
            for i in eachindex(_x)
                n_core_bent[i, j] = real(n_core[i,j]) * (1 - (real(n_core[i,j])^2 *
                                    (_x[i]*cosd(m.bend_direction) + _y[j]*sind(m.bend_direction)) *
                                    m.ρe / (2*m.bending_radius_of_curvature))) *
                                    exp((_x[i]*cosd(m.bend_direction) + _y[j]*sind(m.bend_direction))/m.bending_radius_of_curvature)

            end
        end

        delta_n_2 = n_core_bent.^2 .- m.n0^2

        dz = 1e-10
        x_edge = m.lx_main * (1 + Int(m.ysymmetry != NoSymmetry))/2
        y_edge = m.ly_main * (1 + Int(m.xsymmetry != NoSymmetry))/2
        absorber = Matrix{ComplexF64}(undef, length(_x), length(_y))
        for j in eachindex(_y)
            for i in eachindex(_x)
                absorber[i, j] = exp(-dz * (max(0, max(0, max(abs(_y[j]) - y_edge, abs(_x[i]) - x_edge)))^2*m.alpha + 2π*imag(n_core[i,j])/m.λ))
            end
        end
        ax = 1.00001*dz/(_dx^2*2im*k0*m.n0)
        ay = dz/(_dy^2*2im*k0*m.n0)

        M_rhs = sparse(1:N, 1:N, absorber[1:N] + delta_n_2[1:N]*dz*k0/(2im*m.n0),N,N) +
                sparse(1:N-1,2:N,vec([repeat([fill(ax,1,_nx-1) 0],1,_ny-1) fill(ax,1,_nx-1)]),N,N) +
                sparse(2:N,1:N-1,vec([repeat([fill(ax,1,_nx-1) 0],1,_ny-1) fill(ax,1,_nx-1)]),N,N) +
                sparse(1:N-_nx,1+_nx:N,ay,N,N) +
                sparse(1+_nx:N,1:N-_nx,ay,N,N)
        M_rhs[1:N+1:N*N] = M_rhs[1:N+1:N*N] - vec(repeat([ax fill(2*ax,1,_nx-2) ax],1,_ny))
        M_rhs[1:N+1:N*N] = M_rhs[1:N+1:N*N] - vec([fill(ay,1,_nx) fill(2*ay,1,_nx*(_ny-2)) fill(ay,1,_nx)])
        absorberdiag = sparse(1:N,1:N,absorber[1:N],N,N)
        M_rhs = M_rhs*absorberdiag

        if m.xsymmetry == AntiSymmetric
            M_rhs[1:N+1:N*_nx] .= 0
            M_rhs[N*_nx:N+1:2*N*_nx] .= 0
        end
        if m.ysymmetry == AntiSymmetric
            M_rhs[1:((N+1)*_nx):end] .= 0
            M_rhs[N:((N+1)*_nx):end] .= 0
        end

        D, V = eigs(M_rhs; nev=ceil(Int, n_modes/n_cores), sigma=1.0, ncv=min(N,ceil(Int,n_modes/n_cores)*10))

        κ = @. (1 - real(D))/dz * m.λ/(4π)
        realn = @. √(m.n0^2 - 2*m.n0*imag(D/(dz*k0)))

        neff = realn[realn .> m.nbackground] + im*κ[realn .> m.nbackground]
        V = V[:,realn .> m.nbackground]
        D = D[realn .> m.nbackground]

        for i_core_mode = 1:length(D)
            i_mode += 1
            mode = ElectricFieldProfile(
                lx=_lx,
                ly=_ly,
                xsymmetry=m.xsymmetry,
                ysymmetry=m.ysymmetry
            )
            E = reshape(V[:,i_core_mode], _nx, _ny)
            E = @. E * exp(-im*angle(maximum(E)))
            mode.field = E
            if anycomplex
                mode.n_effective = neff[i_core_mode]
            else
                mode.n_effective = real(neff[i_core_mode])
            end

            if radially_symmetric
                x_full, y_full, E_full = calc_full_field(_x, _y, E)
                i_max = argmax(E_full)
                x_max = x_full[i_max[1]]
                y_max = x_full[i_max[2]]
                θ = atan(y_max - y_c, x_max - x_c)
                itp = Interpolations.LinearInterpolation((x_full, y_full), E_full; extrapolation_bc=Line())
                radialE = itp.(x_c .+ LinRange(0, max(Lx, Ly), 1000)*cos(θ), y_c .+ LinRange(0, max(Lx,Ly), 1000)*sin(θ))
                radialEpruned = radialE[abs.(radialE) .> 0.01*maximum(abs.(radialE))]
                m = sum(abs.(diff(angle.(radialEpruned*exp(im*π/2)) .> 0))) + 1

                R = √((x_max - x_c)^2 + (y_max - y_c)^2)
                azimuthalE = itp.(x_c .+ R*cos.(θ .+ LinRange(0,2π,1000)),y_c .+ R*sin.(θ .+ LinRange(0, 2π, 1000)))
                azimuthalEpruned = azimuthalE[abs.(azimuthalE) .> 0.01*maximum(abs.(azimuthalE))]
                l = sum(abs.(diff(angle.(azimuthalEpruned*exp(im*π/2)) .> 0))) / 2

                if l > 0
                    Emaxmirrored = itp(x_max, 2*y_c - y_max)
                    if real(E_full[i_max]/Emaxmirrored) < 0
                        parity = 'o'
                    else
                        parity = 'e'
                    end
                else
                    parity = ' '
                end
                mode.label = ", LP$l$m$parity"
            end

            push!(m.modes, mode)
        end
    end

    tstop = time()
    dt = round(tstop - tstart, digits=2)
    if isempty(m.modes)
        @info "Done, $dt seconds elapsed. No guided modes found."
        return
    end

    if sort_by_loss
        sorted_idxs = sortperm(imag([m.n_effective for m in m.modes]))
    else
        sorted_idxs = sortperm(real([m.n_effective for m in m.modes]); rev=true)
    end

    m.modes = m.modes[sorted_idxs[1:min(length(m.modes), n_modes)]]

    @info "Done, $dt seconds elapsed. $(length(m.modes)) guided modes found."

    figs = Figure[]
    for i_mode in eachindex(m.modes)
        m.modes[i_mode].label = "Mode $i_mode $(m.modes[i_mode].label)"
        if plot_modes
            E = m.modes[i_mode].field


            fig = Figure()

            ax = fig[1,1] = Axis(fig, aspect=DataAspect())
            heatmap!(_x, _y, abs2.(E))

            ax = fig[1,2] = Axis(fig, aspect=DataAspect())
            heatmap!(_x, _y, angle.(E), colorrange=(-π,π))

            push!(figs, fig)
        end
    end

    figs
end
