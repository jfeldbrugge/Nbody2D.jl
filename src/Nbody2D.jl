"""
    A simple two-dimensional cosmological N-body code.
"""
module Nbody2D
    using QuadGK, FFTW, Random, SpecialFunctions
    using LinearAlgebra, StaticArrays, DelaunayTriangulation
    using Plots

    include("Nbody.jl")
    include("PhaseSpaceDTFE.jl")
    include("Plot.jl")

    export Box, Cosmology, GRF, Zeldovich, LeapFrog
    export PS_DTFE, ps_density, ps_velocity, numberOfStreams
    export plotMesh, plotMesh!
end
