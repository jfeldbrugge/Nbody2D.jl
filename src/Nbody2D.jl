"""
    A simple two-dimensional cosmological N-body code.
"""
module Nbody2D
    using QuadGK, FFTW, Random, SpecialFunctions

    include("Nbody.jl")

    export Box, Cosmology, GRF, Zeldovich, LeapFrog
end
