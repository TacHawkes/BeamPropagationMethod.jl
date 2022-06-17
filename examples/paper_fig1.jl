using Revise, BeamPropagationMethod

m = BeamPropagationMethod.Model()

# general
m.name = @__FILE__
m.use_all_cpus = false
m.use_gpu = false
m.store_field_3D = true

# visualisation
m.updates = 50

# resolution
m.lx_main = 30e-6
m.ly_main = 30e-6
m.nx_main = 400
m.ny_main = 400
m.padfactor = 1.0
m.dz_target = 0.4e-6
m.alpha = 3e14

# problem
m.λ = 800e-9
m.nbackground = 1.4533
m.n0 = 1.46
m.lz = 1e-3
m.intensity_colormap = :jet

function calcRI(x, y, n_background, n_parameters)
    n = n_background * ones(length(x), length(y), 1)
    for j in eachindex(y)
        for i in eachindex(x)
            if (x[i]^2 + y[j]^2) < 10e-6^2
                n[i,j,:] .= 1.4833
            end
        end
    end

    return n
end

function calcInitialE(x, y, Eparameters)
    w_0 = 2.5e-6
    offset = 5e-6
    E = zeros(ComplexF32, length(x), length(y))
    for j in eachindex(y)
        for i in eachindex(x)
            E[i,j] = exp(-((x[i])^2 + (y[j]-offset)^2)/w_0^2) * exp(im*(-sind(5)*x[i]/800e-9*2π))
        end
    end

    return E
end

BeamPropagationMethod.initialize_refractive_index!(m, calcRI)
BeamPropagationMethod.initialize_electric_field!(m, calcInitialE)
##
BeamPropagationMethod.fd_bpm!(m)

##
BeamPropagationMethod.finalize_video!(m)
