function initialize_electric_field!(m::Model,
                                callback::Function,
                                parameters::Union{Dict{Symbol,Any}, Nothing}=nothing)
    E = callback(x(m), y(m), parameters)

    power_fraction = 1/(1 + Int(m.xsymmetry != NoSymmetry)) / (1 + Int(m.ysymmetry != NoSymmetry))
    m.E.field = E ./ √(sum(abs2.(E))/power_fraction)
    m.E.lx = lx(m)
    m.E.ly = ly(m)
    m.E.xsymmetry = m.xsymmetry
    m.E.ysymmetry = m.ysymmetry

    return m
end

function initialize_refractive_index!(m::Model,
                                    callback::Function,
                                    parameters::Union{Dict{Symbol,Any}, Nothing}=nothing,
                                    nz=1)
    if nz == 1
        m.n.n = callback(x(m), y(m), m.nbackground, parameters)
    else
        dz = m.lz / (nz-1)
        z = dz*(0:nz)
        n = callback(x(m), y(m), z, m.nbackground, parameters)
        m.n.n = trim_refractive_index(n, m.nbackground)
    end

    m.n.lx = dx(m)*size(m.n.n,1)
    m.n.ly = dy(m)*size(m.n.n,2)
    m.n.xsymmetry = m.xsymmetry
    m.n.ysymmetry = m.ysymmetry

    return m
end

function trim_refractive_index(n, nbackground)
    n_background = single(nbackground)
    nx, ny = size(n,1), size(n,2)
    xmin = findfirst(any(x->x != n_background, n; dims=(2,3)))
    xmax = findlast(any(x->x != n_background, n; dims=(2,3)))
    xtrim = min(xmin - 1, nx - xmax)
    ymin = findfirst(any(x->x != n_background, n; dims=(1,3)))
    ymax = findlast(any(x->x != n_background, n; dims=(1,3)))
    ytrim = min(ymin - 1, ny - ymax)

    n_temp = n[xtrim+1:nx-xtrim,ytrim+1:ny-ytrim,:]
    if ndims(n_temp) == 3
        n = n_background * ones(eltype(n_temp), (hcat(size(n_temp)...) .+ [2 2 0])...)
    else
        n = n_background * ones(eltype(n_temp), (hcat(size(n_temp)...) .+ [2 2])...)
    end
    n[2:end-1, 2:end-1, :] .= n_temp

    return n
end

function finalize_video!(m, out_file="")
    if isempty(out_file)
        out_file = joinpath(pwd(), m.name * ".mp4")
    end

    save(out_file, m.video_stream)
end

function redlines(m::Model)
    if m.xsymmetry != NoSymmetry && m.ysymmetry != NoSymmetry
        redline_x = [0 m.lx_main m.lx_main]
        redline_y = [m.ly_main m.ly_main 0]
    elseif m.xsymmetry != NoSymmetry
        redline_x = [-m.lx_main -m.lx_main m.lx_main m.lx_main]
        redline_y = [0 m.ly_main m.ly_main 0]
    elseif m.ysymmetry != NoSymmetry
        redline_x = [0 m.lx_main m.lx_main 0]
        redline_y = [-m.ly_main -m.ly_main m.ly_main m.ly_main]
    else
        redline_x = [-m.lx_main m.lx_main m.lx_main -m.lx_main]
        redline_y = [m.ly_main m.ly_main -m.ly_main -m.ly_main]
    end

    return redline_x, redline_y
end

function update_stored_fields!(m::Model, nx, ny, E)
    if m.xsymmetry != NoSymmetry
        y0idx = 1
    else
        y0idx = round(Int, (ny-1)/2+1)
    end
    if m.ysymmetry != NoSymmetry
        x0idx = 1
    else
        x0idx = round(Int, (nx-1)/2+1)
    end

    if m.prior_data
        m.powers = [m.powers;vec(fill(NaN, m.updates))]
        push!(m.xzslice, fill(complex(NaN32), nx, m.updates))
        push!(m.yzslice, fill(complex(NaN32), ny, m.updates))
    else
        m.powers = fill(NaN, m.updates + 1)
        m.powers[1] = 1
        m.xzslice = [fill(complex(NaN32), nx, m.updates+1)]
        m.yzslice = [fill(complex(NaN32), ny, m.updates+1)]
        m.xzslice[1][:,1] = E[:,y0idx]
        m.yzslice[1][:,1] = E[x0idx,:]
    end

    if m.store_field_3D
        if m.prior_data
            push!(m.E3D, fill(complex(NaN32), nx, ny, m.updates))
        else
            m.E3D = [fill(complex(NaN32), nx, ny, m.updates+1)]
            m.E3D[1][:,:,1] = E
        end
    end
end
