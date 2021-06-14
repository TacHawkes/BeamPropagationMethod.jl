var documenterSearchIndex = {"docs":
[{"location":"examples/example1/#Simple-multimode-fiber","page":"Simple multimode fiber","title":"Simple multimode fiber","text":"","category":"section"},{"location":"examples/example1/","page":"Simple multimode fiber","title":"Simple multimode fiber","text":"This example reproduces figure 1 of the BPM-Matlab [Veettikazhy2021] paper.","category":"page"},{"location":"examples/example1/#Description","page":"Simple multimode fiber","title":"Description","text":"","category":"section"},{"location":"examples/example1/","page":"Simple multimode fiber","title":"Simple multimode fiber","text":"A Gaussian beam with w_0=25mu m is launched with an offset of 5mum in the y-direction into the fiber core of a 1 mm long multimode fiber. The beam is titled by 5° towards the positive x-direction by providing a corresponding tilt phase screen for the input electric field. The wavelength is 800 nm. The fiber core has a radius of 10 mu m, a core refractive index of 14833 and a cladding refractive index of 14533. The number of grid points is chosed to result in a resolution of Delta x = Delta y = 0067 mu m and Delta z=04mu m.","category":"page"},{"location":"examples/example1/#The-code","page":"Simple multimode fiber","title":"The code","text":"","category":"section"},{"location":"examples/example1/","page":"Simple multimode fiber","title":"Simple multimode fiber","text":"Run the following code to get figure 1 of the paper. First, the module is loaded. The function calc_initial_field is used to calculate the initial gaussian beam. The parameter X and Y are a grid of points where to calculate the field. Eparameters can be used to pass additional parameters to this function.","category":"page"},{"location":"examples/example1/","page":"Simple multimode fiber","title":"Simple multimode fiber","text":"To prepare for the propagation, a Dict is setup containing all relevant parameters. Finally, the function fdbpm! is called with the dictionary as a parameter. This starts the iteration with the given parameters and a plot is generated and updated step by step. The final electric field as well as other parameters a stored within the dictionary and can be used for additional propagation steps or other purposes.","category":"page"},{"location":"examples/example1/","page":"Simple multimode fiber","title":"Simple multimode fiber","text":"using BeamPropagationMethod\n\nfunction calc_initial_field(X, Y, Eparameters)\n    w_0 = 2.5e-6\n    offset = 5e-6\n    amplitude = @. exp(-((X)^2+(Y-offset)^2)/w_0^2)\n    phase = @. -sind(5)*X/800e-9*2*π\n    E = @. amplitude*exp(im*phase)\n\n    return E\nend\n\np = Dict(    \n    :useGPU => false,\n\n    :updates => 50,    \n\n    :Lx_main => 30e-6,\n    :Ly_main => 30e-6,\n    :Nx_main => 400,\n    :Ny_main => 400,\n    :padfactor => 1.0,\n    :dz_target => 0.4e-6,\n    :alpha => 3e14,\n\n    :lambda => 800e-9,\n    :n_background => 1.4533,\n    :n_0 => 1.4533,\n    :Lz => 1e-3,\n\n    :shapes => [0 0 10e-6 2 1.4833],\n    :E => calc_initial_field,\n    \n    :Intensity_colormap => :jet\n)\n\nfdbpm!(p);","category":"page"},{"location":"examples/example1/","page":"Simple multimode fiber","title":"Simple multimode fiber","text":"[Veettikazhy2021]: Madhu Veettikazhy, Anders Kragh Hansen, Dominik Marti, Stefan Mark Jensen, Anja Lykke Borre, Esben Ravn Andresen, Kishan Dholakia, and Peter Eskil Andersen, \"BPM-Matlab: an open-source optical propagation simulation tool in MATLAB,\" Opt. Express 29, 11819-11832 (2021). Link to article","category":"page"},{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"This package is a implementation of BPM-Matlab with some extensions and different routes. It is a numerical tool for electric field propagation in optical waveguide structures using the beam propagation method. For an in-depth explanation of the purpose and the math read the paper[Veettikazhy2021] by the original BPM-Matlab authors. The package is for the most part a direct reimplementation of the Matlab code, however there are some differences (see below). The most important benefit of this package is the pure Julia implementation. Therefore this package runs on any system which runs Julia 1.6. I have not performed an \"scientific\" benchmarks but replicating figure 1 of [Veettikazhy2021] takes less than 30s on a i7-4850HQ with only a single update and 2 minutes with 50 updates, indicating that performance is limited by the visualization at the moment.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The syntax for using BeamPropagationMethod.jl is very simlar to BPM-Matlab but had to be adapted for Julia specifics. Check the examples in order to get started. It might be worthwile to check the BPM-Matlab Readme for an explanation of the structure.","category":"page"},{"location":"#Features","page":"Introduction","title":"Features","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Native Julia implementation\nMulti-Threading support\nCalculate electric field propagation for arbitrary refractive index profiles\nFD-BPM and FFT (Fresnel) propagation code for solving the Helmholtz equation (i.e. only low refractive index contrast, no vector information, polarization, ...)\nMode solver for determining the eigenmodes of an refractive index profile\nSupport for non-trivially transforming refractive mode index profiles, allowing to model waveguide building blocks (Splitters, Mach-Zehnder-Modulators, Photonic Lanterns)\nVisualization and video export using Plots.jl\nProgress monitor","category":"page"},{"location":"#Installation","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Install the package using the package manager:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"] add BeamPropagationMethod","category":"page"},{"location":"#License","page":"Introduction","title":"License","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"You should have received a copy of the GNU General Public License along with this program. If not, see [https://www.gnu.org/licenses/].","category":"page"},{"location":"#Citation","page":"Introduction","title":"Citation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"As soon as I consider the package to be release ready, I will register it with the Julia registry and Zenodo. For now just have fun with it and do not forget to cite the authors of BPM-Matlab.","category":"page"},{"location":"#Notable-differences-to-Matlab-BPM","page":"Introduction","title":"Notable differences to Matlab-BPM","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"No CUDA support yet (you can thank NVIDIA for the lack of MacOS CUDA drivers, so I am unable to use CUDA on my Macbook Pro with GT750M)\nPerformance-optimized propagation code with more caching and code structure which can be optimized using native vector operations\nSlightly different naming for methods\nAPI is ready for complex waveguiding building blocks\nLess sophisticated visualization due to the lack of support in Plots.jl, I plan to improve on this point\nMode-solver is sometimes a bit unreliable in enumerating modes","category":"page"},{"location":"#References","page":"Introduction","title":"References","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"[Veettikazhy2021]: Madhu Veettikazhy, Anders Kragh Hansen, Dominik Marti, Stefan Mark Jensen, Anja Lykke Borre, Esben Ravn Andresen, Kishan Dholakia, and Peter Eskil Andersen, \"BPM-Matlab: an open-source optical propagation simulation tool in MATLAB,\" Opt. Express 29, 11819-11832 (2021). Link to article","category":"page"}]
}
