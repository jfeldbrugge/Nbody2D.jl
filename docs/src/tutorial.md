```@meta
CurrentModule = Nbody2D
```

# Tutorial

In this tutorial, we demonstrate the usage of the *Nbody2D.jl* 

## Cosmology and initial conditions
```@example tutorial1
using Nbody2D, Plots

let Ni = 2^7, L = 50., H0 = 70., OmegaM = 1.0, OmegaL = 0.
    global box = Box(2, Ni, L)
    global EdS = Cosmology(H0, OmegaM, OmegaL)
end

let ns = 2.; Rs = 1.; α = 1.; seed = 1
    global phi = GRF(ns, Rs, α, seed, box)
end

qRange = range(0, box.L; length=box.N)
heatmap(qRange, qRange, phi, aspect_ratio=:equal, xlims=(0, box.L), ylims=(0, box.L), title="Displacement potential")
```

## N-body simulation
```@example tutorial1
ai, af, δa = 0.02, 2.02, 0.02

stateZeldovich = Zeldovich(phi, af, af, box, EdS)
zeldovichPlot = plotMesh(stateZeldovich, box)

stateNbody = LeapFrog(phi, ai, af, δa, box, EdS)
nbodyPlot = plotMesh(stateNbody, box)

plot(zeldovichPlot, nbodyPlot, size = (1200, 600), legend = false, title=["Zel'dovich" "N-body"])
```

## Phase-Space Delaunay Tessellation Density Field Estimator

```@example tutorial1
let depth = 10
    stateInitial = Zeldovich(phi, δa, δa, box, EdS)
    global estimator = PS_DTFE(stateInitial, stateNbody, depth, box)
    nothing
end

rangeX = 0 : box.L / (8. * box.N) : box.L
dens   = [ps_density([x, y], estimator) for y in rangeX, x in rangeX]
number = [numberOfStreams([x, y], estimator) for y in rangeX, x in rangeX]

densityPlot = heatmap(rangeX, rangeX, log10.(dens), aspect_ratio=:equal)
plotMesh!(stateNbody, box, :yellow)

numberPlot = heatmap(rangeX, rangeX, log10.(number), aspect_ratio=:equal)
plotMesh!(stateNbody, box, :yellow)

plot(densityPlot, numberPlot, size = (1200, 500), legend = false, title=["Density" "Number of streams"])
```