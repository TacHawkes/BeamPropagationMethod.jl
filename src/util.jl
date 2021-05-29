function ndgrid(v1::AbstractVector{T}, v2::AbstractVector{T}) where T
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repeat(v1, 1, n), repeat(v2, m, 1))
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