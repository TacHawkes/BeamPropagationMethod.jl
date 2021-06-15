# The parameter dictionary

Central to the propagation is a parameters dictionary which is setup at the beginning and can be reused to chain multiple propagation segments together. It is possible to change parameters in between segments.

The dictionary should be of the type `Dict{Symbol,Any}`.

## General parameters

* `:name` - A descriptive name for the model (currently unused)
* `:useGPU`- Reserved for future use of `CUDA.jl``

## Visualization parameters

* `:updates` -  The number of updates to the plot. This is usefull for live coverage of the simulation status but adds a significant plotting overhead. Set this to `1` if ou are only interested in the field at the end of the waveguide.
* `:plotEmax` - If left unset, the colorbar of the electric field display autoscales. Otherwise it will scale the colobar to the maximum of the initial field times this value.

## Resolution parameters

* `:Lx_main, :Ly_main` - The size of the main simulation region in meters for x-/y-direction respectively
* `:Nx_main, :Ny_main` - The number of grid points for the main simulation region for x-/y-direction respectively
* `dz_target` - The targeted ``\Delta z``resolution between two simulation steps in meters
* `padfactor` - This adds an additional outer area to the main simulation region which acts as absorber for field components leaving the waveguide. The `padfactor` specifies the ration between the size of the full simulation window (`Lx`, `Ly`) and the main simulation region (`Lx_main`, `Ly_main`).
* `alpha` - The absorption strength of the padding region. The absorbtion coefficient of the padding is proportional to the squared distance from the edge of the main simulation region. This means that field components which propagate further out of the padding are attenuated more strongly. The coefficient `alpha` is therefore a volumentric absortion coefficient in ``1/m^3``.

## Physical properties

* `lambda` - The wavelength of the light in meters
* `n_background` - The refractive index of the surrounding background, typically this is the cladding refractive index