```@meta
CurrentModule = Nbody2D
```

# Tutorial

In this tutorial, we demonstrate the usage of the *Nbody2D.jl*. 

## Cosmology and initial conditions
We set the cosmological background in which the $N$-body simulation operates and sample a Gaussian random field to create initial conditions for the $N$-body simulation.

```@example tutorial1
using Nbody2D, Plots

# Set the simulation box and cosmological model.
let Ni = 2^7, L = 50., H0 = 70., OmegaM = 1.0, OmegaL = 0.
    global box = Box(2, Ni, L)
    global EdS = Cosmology(H0, OmegaM, OmegaL)
end

# Sample a realization of a Gaussian random field.
let ns = 2.; Rs = 1.; α = 1.; seed = 1
    global phi = GRF(ns, Rs, α, seed, box)
end

# Plot the deformation potential.
qRange = range(0, box.L; length=box.N)
heatmap(qRange, qRange, phi, 
        aspect_ratio=:equal, 
        xlims=(0, box.L), ylims=(0, box.L), 
        title="Displacement potential")
```

## N-body simulation
Given the initial conditions, we evolve the $N$-body particles with both the Zel'dovich approximation and an $N$-body simulation.

```@example tutorial1
ai, af, δa = 0.02, 2.02, 0.02

# Zel'dovich approximation
stateZeldovich = Zeldovich(phi, af, af, box, EdS)
zeldovichPlot = plotMesh(stateZeldovich, box)

# N-body simulation
stateNbody = LeapFrog(phi, ai, af, δa, box, EdS)
nbodyPlot = plotMesh(stateNbody, box)

# Plot the particles as a mesh.
plot(zeldovichPlot, nbodyPlot, size = (1200, 600), legend = false, 
     title=["Zel'dovich" "N-body"])
```

## Phase-Space Delaunay Tessellation Density Field Estimator
We evaluate the density, velocity and number of stream fields with the Phase-Space Delaunay Tessellation Density Field Estimator (PS-DTFE).

```@example tutorial1
# Build PS-DTFE estimators.
let depth = 10
    stateInitial = Zeldovich(phi, δa, δa, box, EdS)
    global estimatorZeldovich = PS_DTFE(stateInitial, stateZeldovich, depth, box)
    global estimatorNbody = PS_DTFE(stateInitial, stateNbody, depth, box)
    nothing
end

# Evaluate density and number of stream fields.
rangeX = 0 : box.L / (8. * box.N) : box.L

densZeldovich   = [ps_density([x, y], estimatorZeldovich) for y in rangeX, x in rangeX]
numberZeldovich = [numberOfStreams([x, y], estimatorZeldovich) for y in rangeX, x in rangeX]

densNbody   = [ps_density([x, y], estimatorNbody) for y in rangeX, x in rangeX]
numberNbody = [numberOfStreams([x, y], estimatorNbody) for y in rangeX, x in rangeX]

# Plot the density and number of stream fields.
densityZeldovichPlot = heatmap(rangeX, rangeX, log10.(densZeldovich), aspect_ratio=:equal)
plotMesh!(stateZeldovich, box, :yellow)

numberZeldovichPlot = heatmap(rangeX, rangeX, log10.(numberZeldovich), aspect_ratio=:equal)
plotMesh!(stateZeldovich, box, :yellow)

densityNbodyPlot = heatmap(rangeX, rangeX, log10.(densNbody), aspect_ratio=:equal)
plotMesh!(stateNbody, box, :yellow)

numberNbodyPlot = heatmap(rangeX, rangeX, log10.(numberNbody), aspect_ratio=:equal)
plotMesh!(stateNbody, box, :yellow)

plot(densityZeldovichPlot, numberZeldovichPlot, densityNbodyPlot, numberNbodyPlot, 
     size = (1200, 1000), legend = false, 
     title=["Density Zeldovich" "Number of streams Zeldovich" "Density N-body" "Number of streams N-body"])
```