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

function swap_pointers!(a, b)
    tmp = a
    a = b
    b = tmp
end

field_power(y) = sum(x->abs2(x), y)

sigmoid(x) = 1-1/(1+exp(x))

single(x::Float64) = convert(Float32, x)
single(x::ComplexF64) = convert(ComplexF32, x)
