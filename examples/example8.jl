using Revise, BeamPropagationMethod

m = BeamPropagationMethod.Model()

# general
m.name = @__FILE__
m.use_all_cpus = true
m.use_gpu = false

# visualisation
m.updates = 100

# resolution
m.lx_main = 30e-6
m.ly_main = 25e-6
m.nx_main = 300
m.ny_main = 250
m.padfactor = 1.5
m.dz_target = 5e-7
m.alpha = 3e14

# problem definition
m.λ = 1000e-9
m.nbackground = 1.45
m.n0 = 1.46
m.lz = 2e-3

function calcRI(x, y, n_background, n_parameters)
    n = n_background * ones(length(x), length(y), 1)
    corepos = [5e-6 0]
    r = 2.5e-6
    for j in eachindex(y)
        for i in eachindex(x)
            if ((x[i]+2.5e-6)^2 + (y[j])^2) < 5e-6^2
                n[i,j,:] .= 1.46
            end

            R = sqrt((x[i]-corepos[1])^2 + (y[j] - corepos[2])^2)
            if R < r
                n[i,j,:] .= n_background + (1.47 - n_background) * (1 - (R/r)^2)
            end
        end
    end

    return n
end

BeamPropagationMethod.initialize_refractive_index!(m, calcRI)

figs = BeamPropagationMethod.find_modes!(m, 10)
##
for fig in figs
    display(fig)
end
##
m.E = m.modes[9]

##
BeamPropagationMethod.fd_bpm!(m)

##
BeamPropagationMethod.finalize_video!(m)
