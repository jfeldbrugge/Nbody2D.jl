"""
    A simple two-dimensional cosmological N-body code.
"""
module Nbody2D
    using QuadGK, FFTW, Random, SpecialFunctions
    using LinearAlgebra, DelaunayTriangulation
    using Plots

    include("scalarFields.jl")
    include("cGRF.jl")
    include("Nbody.jl")
    include("PhaseSpaceDTFE.jl")
    include("Plot.jl")

    export Box, gradient, hessian, laplacian, smooth
    export GRF, cGRF, constraintGRF, meanField, varianceField, measureConstraints
    export Cosmology, Zeldovich, LeapFrog
    export PS_DTFE, ps_density, ps_velocity, numberOfStreams
    export plotMesh, plotMesh!
end
