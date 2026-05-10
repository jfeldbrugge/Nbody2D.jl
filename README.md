# Nbody2D.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfeldbrugge.github.io/Nbody2D.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfeldbrugge.github.io/Nbody2D.jl/dev/)
[![Build Status](https://github.com/jfeldbrugge/Nbody2D.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jfeldbrugge/Nbody2D.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jfeldbrugge/Nbody2D.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jfeldbrugge/Nbody2D.jl)

This package is a particle mesh code that is inspired by Johan Hidding's [nbody2d](https://zenodo.org/records/4158731) Python code. See [jhidding.github.io/nbody2d](jhidding.github.io/nbody2d) for more details. It serves as a quick environment to learn and experiment with cosmological $N$-body simulations in the two-dimensional setting. Using the package, we can sample initial conditions, evolve the initial perturbations to the current epoch and estimate the density, velocity and number of stream fields.


[![Density and number of stream fields of the Zel'dovich approximation and an $N$-body simulation.](docs/src/assets/figures/density.png)](docs/src/assets/figures/density.svg)

## Installation
The Nbody2D package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add Nbody2D
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("Nbody2D")
```

## Usage
Please have a look at the Tutorial page for details on how to use this package.

## Contributors
This code was written by:
* Job Feldbrugge ([job.feldbrugge@ed.ac.uk](mailto:job.feldbrugge@ed.ac.uk))
