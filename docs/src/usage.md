# The parameter dictionary

Central to the propagation is a parameters dictionary which is setup at the beginning and can be reused to chain multiple propagation segments together. It is possible to change parameters in between segments.

The dictionary should be of the type `Dict{Symbol,Any}`. It is recommended to look at the first example ([Simple multimode fiber](@ref)) to get a first idea how to use the dictionary.

## General parameters

* `:name` - A descriptive name for the model (currently unused)
* `:useGPU`- Reserved for future use of `CUDA.jl`

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
* `n_0` - The refractive index of the Helmholtz-equation. This should be chosed as the refractive index where the majority of the energy is propagating (i.e. close to the core refractive index).
* `L_z` - The length of the simulation segment in meters
* `shapes` - This is a 2D array which allows a simple way to define a refractive index distribution with a circular shape. Note that more sophisticated refractiv index profiles can be supplied directly using the `n` parameter or the `nFunc` parameter, where non-trivial transforming profiles (i.e. splitters/couplers, Mach-Zehnder modulators or photonic lanterns) can be simulated. The `shapes` matrix is structured as follows: The first column is the `x` coordinate and the second column the `y` coordinate. The radius of the circular core is defined in the third column. Column 4 is used to define the type of the shape. Column 5 is the refractive index of the shape and column 6 is the g parameter of gradient index (GRIN) lenses. The following shape types are supported:
    * 1. Circular step-index disk
    * 2. Anti-aliased circular step-index disk
    * 3. Parabolic graded index disk
    * 4. GRIN lens focusing in x- and y-direction
    * 5. GRIN lens focusing only in the y-direction