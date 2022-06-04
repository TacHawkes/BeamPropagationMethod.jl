using Revise, BeamPropagationMethod

m = BeamPropagationMethod.Model()

# general
m.name = @__FILE__
m.use_all_cpus = false
m.use_gpu = false

# visualisation
m.updates = 100
m.plot_field_max = 0.5

# resolution
m.lx_main = 20e-6
m.ly_main = 20e-6
m.nx_main = 200
m.ny_main = 200
m.padfactor = 1.5
m.dz_target = 1e-6
m.alpha = 3e14

# problem
m.λ = 1000e-9
m.nbackground = 1.45
m.n0 = 1.46
m.lz = 2e-3

function calcRI(x, y, n_background, n_parameters)
    n = n_background * ones(length(x), length(y), 1)
    for j in eachindex(y)
        for i in eachindex(x)
            if (x[i]^2 + y[j]^2) < 5e-6^2
                n[i,j,:] .= 1.46
            end
        end
    end

    return n
end

function calcInitialE(x, y, Eparameters)
    w_0 = 2.5e-6
    offset = 2.5e-6
    E = zeros(ComplexF32, length(x), length(y))
    for j in eachindex(y)
        for i in eachindex(x)
            E[i,j] = exp(-((x[i]-offset)^2 + y[j]^2)/w_0^2)
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
