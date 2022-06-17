using Revise, BeamPropagationMethod

m = BeamPropagationMethod.Model()

# general
m.name = @__FILE__
m.use_all_cpus = true
m.use_gpu = false

# visualisation
m.updates = 30
m.plot_field_max = 1

# resolution
m.lx_main = 20e-6
m.ly_main = 10e-6
m.nx_main = 200
m.ny_main = 100
m.padfactor = 1.5
m.dz_target = 1e-6
m.alpha = 3e14

# problem definition
m.λ = 1000e-9
m.nbackground = 1.45
m.n0 = 1.47
m.lz = 2e-4

m.xsymmetry = BeamPropagationMethod.Symmetric
m.ysymmetry = BeamPropagationMethod.NoSymmetry

function calcRI(x, y, n_background, n_parameters)
    n = n_background * ones(length(x), length(y), 1)
    for j in eachindex(y)
        for i in eachindex(x)
            if (x[i]^2 + y[j]^2) < 5e-6^2
                n[i,j,:] .= 1.47
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
BeamPropagationMethod.initialize_refractive_index!(m, calcRI)
BeamPropagationMethod.initialize_electric_field!(m, calcInitialE)
m.xsymmetry = BeamPropagationMethod.NoSymmetry
m.ysymmetry = BeamPropagationMethod.NoSymmetry
m.ly_main = 20e-6
m.ny_main = 200

##
BeamPropagationMethod.fd_bpm!(m)

##
BeamPropagationMethod.initialize_refractive_index!(m, calcRI)
BeamPropagationMethod.initialize_electric_field!(m, calcInitialE)
m.xsymmetry = BeamPropagationMethod.AntiSymmetric
m.ysymmetry = BeamPropagationMethod.Symmetric
m.lx_main = 10e-6
m.nx_main = 100
m.ly_main = 10e-6
m.ny_main = 100
##
figs = BeamPropagationMethod.find_modes!(m, 10, plot_modes=true)
##
for fig in figs
    display(fig)
end

##
m.E = BeamPropagationMethod.mode_superposition(m, 1:4)

##
BeamPropagationMethod.fd_bpm!(m)

##
m.xsymmetry = BeamPropagationMethod.NoSymmetry
m.ysymmetry = BeamPropagationMethod.NoSymmetry
m.lx_main = 20e-6
m.nx_main = 200
m.ly_main = 20e-6
m.ny_main = 200
BeamPropagationMethod.initialize_refractive_index!(m, calcRI)
m.E = BeamPropagationMethod.mode_superposition(m, 1:4)

##
BeamPropagationMethod.fd_bpm!(m)

##
BeamPropagationMethod.finalize_video!(m)
