[![Build Status Master](https://travis-ci.org/ALaDyn/piccante.png?branch=master)](https://travis-ci.org/ALaDyn/piccante "master")

``piccante`` is a massively parallel fully-relativistic electromagnetic 3D particle-in-cell code, released by the authors to the whole laser-plasma community under a GPLv3 license.  
Its strengths are related to flexibility: the user is able to define an arbitrary number of particle species with arbitrary density functions and an arbitrary number of Gaussian laser pulses.  
The code is proved to run on up to 16384 MPI processes. Some preliminary scaling tests verified a very good scalability on up to 4096 MPI tasks.  
In its present state, the code could be considered in *beta status*.  
Copyright 2014, 2015 - Andrea Sgattoni, Luca Fedeli, Stefano Sinigardi

To build it on linux, ``make``, a C++ compiler, ``boost``, ``MPI`` and ``gsl`` are required, after that a plain `make` should work (adapt the makefile according to your setup if needed).  
On MacOS X you should use [homebrew](http://brew.sh/) to install the same packages and then launch `make brew`.  
On Windows, you should use the Visual Studio solution to obtain ``boost`` as a Nuget package, while you should download [GSL for Windows](http://sourceforge.net/projects/gnuwin32/files/gsl/1.8/gsl-1.8.exe) and [Microsoft MPI (with SDK)](https://www.microsoft.com/en-us/download/details.aspx?id=49926) before building. After installing GSL you should also create its .lib files with this method:
```
1) open the VS Command Shell as an administrator
2) cd "c:\Program Files (x86)\GnuWin32\lib"
3) lib /def:libgsl.def
4) lib /def:libgslcblas.def
```
`piccante` should build cleanly after that

A paper that describes the technical features of ``piccante`` is under development now, in the meantime if you want to use the code please ask the authors for a ZENODO doi linked to the release you're using.  
A technical paper describing optimizations performed during a Preparatory PRACE can be found [here](http://www.prace-ri.eu/IMG/pdf/WP209.pdf) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.16097.svg)](http://dx.doi.org/10.5281/zenodo.16097).

### Papers containing simulations done with piccante ###

1) A. Sgattoni, S. Sinigardi, L. Fedeli, F. Pegoraro, A. Macchi, Laser-Driven Rayleigh-Taylor Instability: Plasmonics Effects and Three-Dimensional Structures, [arXiv:1404.1260](http://arxiv.org/pdf/1404.1260.pdf), [Phys. Rev. E, 91, 013106 (2015)](http://dx.doi.org/10.1103/PhysRevE.91.013106)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.10587.svg)](http://dx.doi.org/10.5281/zenodo.10587)

2) A. Sgattoni, S. Sinigardi, A. Macchi,  High Energy Gain in Three-Dimensional Simulations of Light Sail Acceleration, [arXiv:1403.2709](http://arxiv.org/pdf/1403.2709.pdf), [Appl. Phys. Lett., 105, 084105 (2014)](http://dx.doi.org/10.1063/1.4894092)
