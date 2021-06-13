function ndgrid(v1::AbstractVector{T}, v2::AbstractVector{T}) where T
    m, n = length(v1), length(v2)
    v1 = reshape(v1, 1, m)
    v2 = reshape(v2, n, 1)
    (repeat(v1, n, 1), repeat(v2, 1, m))
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

function field_power(x)
    y = 0.0
    @inbounds @simd for i in eachindex(x)
        y += abs2(x[i])
    end
    return y
end

sigmoid(x) = 1-1/(1+exp(x))