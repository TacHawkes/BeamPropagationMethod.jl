using Revise, BeamPropagationMethod

m = BeamPropagationMethod.Model()

# general
m.name = @__FILE__
m.use_all_cpus = true
m.use_gpu = false

# visualisation
m.updates = 30
m.plot_zoom = 1

# resolution
m.lx_main = 50e-6
m.ly_main = 50e-6
m.nx_main = 200
m.ny_main = 200
m.padfactor = 1.5
m.dz_target = 1e-6
m.alpha = 3e14

# problem definition
m.λ = 1000e-9
m.nbackground = 1.45
m.n0 = 1.46
m.lz = 2e-3
m.taper_scaling = 1
m.twist_rate = 0
m.fig_title = "Segment 1"

function calcRI(x, y, n_background, n_parameters)
    n = n_background * ones(length(x), length(y), 1)
    corepos = [2e-6 12e-6]
    r = 10e-6
    for j in eachindex(y)
        for i in eachindex(x)
            if ((x[i]+7e-6)^2 + (y[j] + 7e-6)^2) < 10e-6^2
                n[i,j,:] .= 1.46
            end

            if ((x[i]-15e-6)^2 + (y[j] + 0e-6)^2) < 1.25e-6^2
                n[i,j,:] .= 1.46
            end

            R3 = sqrt((x[i]-corepos[1])^2 + (y[j] - corepos[2])^2)
            if R3 < r
                n[i,j,:] .= n_background + (1.465 - n_background) * (1 - (R3/r)^2)
            end
        end
    end

    return n
end

function calcInitialE(x, y, Eparameters)
    w_0 = 5e-6
    E = zeros(ComplexF32, length(x), length(y))
    for j in eachindex(y)
        for i in eachindex(x)
            E[i,j] = exp(-((x[i]-1.5e-5)^2 + y[j]^2)/w_0^2) +
                    2*exp(-((x[i]+12e-6)^2 + (y[j]+7e-6)^2)/w_0^2)*exp(im*8e5*y[j])
        end
    end

    return E
end

BeamPropagationMethod.initialize_refractive_index!(m, calcRI)
BeamPropagationMethod.initialize_electric_field!(m, calcInitialE)
##
BeamPropagationMethod.fd_bpm!(m)

##
m.fig_title = "Segment 2"
m.lz = 5e-3
m.taper_scaling = 0.15
m.twist_rate = 2π/m.lz
##
BeamPropagationMethod.fd_bpm!(m)

##
m.fig_title = "Segment 3"
m.lz = 2e-3
m.taper_scaling = 1
m.twist_rate = 0

##
BeamPropagationMethod.fd_bpm!(m)

##
BeamPropagationMethod.finalize_video!(m)
