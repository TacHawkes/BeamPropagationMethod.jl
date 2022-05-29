function ndgrid(v1::AbstractVector{T}, v2::AbstractVector{T}) where T
    m, n = length(v1), length(v2)
    v1 = reshape(v1, 1, m)
    v2 = reshape(v2, n, 1)
    (repeat(v1, n, 1), repeat(v2, 1, m))
end

function get_grid_array(nx, dx, symmetry)
    if symmetry == NoSymmetry
        x = dx*(-nx/2+1/2:nx/2-1/2)
    elseif symmetry == Symmetric
        x = dx * (1/2:nx-1/2)
    else
        x = dx*(0:nx-1)
    end

    return x
end

function get_colormap(cmap::Int)
    if cmap == 1
        return :inferno
    elseif cmap == 2
        return :jet
    elseif cmap == 3
        return :balance
    elseif cmap == 4
        return :gray
    elseif cmap == 5
        return :cividis
    end
end

get_colormap(cmap::Symbol) = cmap

function swap_pointers!(a, b)
    tmp = a
    a = b
    b = tmp
end

field_power(y) = sum(x->abs2(x), y)

sigmoid(x) = 1-1/(1+exp(x))

single(x::Float64) = convert(Float32, x)
single(x::ComplexF64) = convert(ComplexF32, x)
