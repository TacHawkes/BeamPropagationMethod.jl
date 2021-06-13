# Introduction

This package is a implementation of [BPM-Matlab](https://gitlab.gbar.dtu.dk/biophotonics/BPM-Matlab) with some extensions and different routes.
It is a numerical tool for electric field propagation in optical waveguide structures using the beam propagation method. For an in-depth explanation
of the purpose and the math read the paper[^Veettikazhy2021] by the original BPM-Matlab authors.
The package is for the most part a direct reimplementation of the Matlab code, however there are some differences (see below). The most important benefit of this package is the pure Julia implementation. Therefore this package runs on any system which runs Julia 1.6. I have not performed an "scientific" benchmarks but replicating figure 1 of [^Veettikazhy2021] takes less than 30s on a i7-4850HQ with only a single update and 2 minutes with 50 updates, indicating that performance is limited by the visualization at the moment.

The syntax for using `BeamPropagationMethod.jl` is very simlar to BPM-Matlab but had to be adapted for Julia specifics. Check the examples in order to get started. It might be worthwile to check the [BPM-Matlab Readme](https://gitlab.gbar.dtu.dk/biophotonics/BPM-Matlab/blob/Release/README.md) for an explanation of the structure.

## Features

* Native Julia implementation
* Multi-Threading support
* Calculate electric field propagation for arbitrary refractive index profiles
* FD-BPM and FFT (Fresnel) propagation code for solving the Helmholtz equation (i.e. only low refractive index contrast, no vector information, polarization, ...)
* Mode solver for determining the eigenmodes of an refractive index profile
* Support for non-trivially transforming refractive mode index profiles, allowing to model waveguide building blocks (Splitters, Mach-Zehnder-Modulators, Photonic Lanterns)
* Visualization and video export using `Plots.jl`
* Progress monitor


## Installation

Install the package using the package manager:

```julia
] add BeamPropagationMethod
```

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

## Citation

As soon as I consider the package to be release ready, I will register it with the Julia registry and Zenodo. For now just have fun with it and do not forget to cite the authors of BPM-Matlab.

## Notable differences to Matlab-BPM

* No CUDA support yet (you can thank NVIDIA for the lack of MacOS CUDA drivers, so I am unable to use CUDA on my Macbook Pro with GT750M)
* Performance-optimized propagation code with more caching and code structure which can be optimized using native vector operations
* Slightly different naming for methods
* API is ready for complex waveguiding building blocks
* Less sophisticated visualization due to the lack of support in Plots.jl, I plan to improve on this point
* Mode-solver is sometimes a bit unreliable in enumerating modes

## References

[^Veettikazhy2021]: Madhu Veettikazhy, Anders Kragh Hansen, Dominik Marti, Stefan Mark Jensen, Anja Lykke Borre, Esben Ravn Andresen, Kishan Dholakia, and Peter Eskil Andersen, "BPM-Matlab: an open-source optical propagation simulation tool in MATLAB," Opt. Express 29, 11819-11832 (2021). [Link to article](https://doi.org/10.1364/OE.420493)