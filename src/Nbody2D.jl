"""
    A simple two-dimensional cosmological N-body code.
"""
module Nbody2D
    using QuadGK, FFTW, Random, SpecialFunctions
    using LinearAlgebra, DelaunayTriangulation
    using Plots
    using Interpolations

    include("ScalarFields.jl")
    include("cGRF.jl")
    include("Nbody.jl")
    include("PhaseSpaceDTFE.jl")
    include("Caustics.jl")
    include("Plot.jl")

    export Box, gradient, hessian, laplacian, smooth
    export GRF, cGRF, constraintGRF, meanField, varianceField, measureConstraints, generalizedMoment
    export Cosmology, Zeldovich, LeapFrog, LeapFrogSnapshots
    export PS_DTFE, ps_density, ps_velocity, numberOfStreams
    export Caustics!, eigenFields123, Eulerian
    export plotMesh, plotMesh!, plotSkeleton, plotSkeleton!
end
