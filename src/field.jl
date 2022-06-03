dx(m::Model) = m.lx_main / m.nx_main
dy(m::Model) = m.ly_main / m.ny_main
dz(m::Model) = m.lz / nz(m)

function nx(m::Model)
    target_lx = m.padfactor * m.lx_main
    nx = round(Int, target_lx/dx(m))
    if m.ysymmetry != NoSymmetry
        nx += (rem(nx, 2) != rem(m.nx_main, 2)) ? 1 : 0
    end

    return nx
end

function ny(m::Model)
    target_ly = m.padfactor * m.ly_main
    ny = round(Int, target_ly/dy(m))
    if m.xsymmetry != NoSymmetry
        ny += (rem(ny, 2) != rem(m.ny_main, 2)) ? 1 : 0
    end

    return ny
end

nz(m::Model) = max(m.updates, round(Int, m.lz/m.dz_target))

lx(m::Model) = dx(m) * nx(m)
ly(m::Model) = dy(m) * ny(m)

x(m::Model) = get_grid_array(nx(m), dx(m), m.xsymmetry)
y(m::Model) = get_grid_array(ny(m), dy(m), m.ysymmetry)

function calc_full_field(x, y, E)
    if sign(minimum(x)) == -1
        # no op
    elseif sign(minimum(x)) == 0
        x = cat(-reverse(x[2:end]; dims=1), x; dims=1)
        E = cat(-reverse(E[2:end,:]; dims=1), E; dims=1)
    elseif sign(minimum(x)) == 1
        x = cat(-reverse(x; dims=1), x; dims=1)
        E = cat(reverse(E; dims=1), E; dims=1)
    end

    if sign(minimum(y)) == -1
        # no op
    elseif sign(minimum(y)) == 0
        y = cat(-reverse(y[2:end]; dims=1), y; dims=1)
        E = cat(-reverse(E[2:end,:]; dims=2), E; dims=2)
    elseif sign(minimum(y)) == 1
        y = cat(-reverse(y; dims=1), y; dims=1)
        E = cat(reverse(E; dims=2), E; dims=2)
    end

    return x, y, E
end

function calc_full_refractive_index(x, y, n)
    if sign(minimum(x)) == -1
        # no op
    elseif sign(minimum(x)) == 0
        x = cat(-reverse(x[2:end]; dims=1), x; dims=1)
        n = cat(reverse(n[2:end,:,:]; dims=1), n; dims=1)
    elseif sign(minimum(x)) == 1
        x = cat(-reverse(x; dims=1), x; dims=1)
        n = cat(reverse(n; dims=1), n; dims=1)
    end

    if sign(minimum(y)) == -1
        # no op
    elseif sign(minimum(y)) == 0
        y = cat(-reverse(y[2:end]; dims=1), y; dims=1)
        n = cat(reverse(n[2:end,:,:]; dims=2), n; dims=2)
    elseif sign(minimum(y)) == 1
        y = cat(-reverse(y; dims=1), y; dims=1)
        n = cat(reverse(n; dims=2), n; dims=2)
    end

    return x, y, n
end
