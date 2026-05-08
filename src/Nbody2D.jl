"""
    A simple two-dimensional N-body code.
"""
module Nbody2D
    using QuadGK, FFTW, Random, SpecialFunctions, ProgressMeter

    include("Nbody_2D.jl")

    export Box, Cosmology, GRF, zeldovich, LeapFrog
end
