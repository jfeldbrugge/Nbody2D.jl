```@meta
CurrentModule = Nbody2D
```

# Theory
We give a brief theoretical description of large-scale structure formation to provide a foundation for this package. We discuss the initial conditions, the evolution of the initial conditions into the cosmic web and the Phase-Space Delaunay Tessellation Field Estimator. [To be extended in a future version.]

## Initial conditions 
This package implements the initial conditions as a realisation of a Gaussian random field. In particular, we model the primordial displacement potential with a Gaussian random field governed by the power-law power spectrum smoothed with a Gaussian filter

```math
P_\phi(k) = \frac{α^2  4 \pi R_s^{2 + n_s}}{\Gamma\left(1 + \frac{n_s}{2}\right)}  k^{n_s - 4}  e^{-R_s^2  k^2}\,,
```

with the spectral index $n_s$, the smoothing length $R_s$, and normalization constant $\alpha$.

## Structure formation
In this code, we perform a two-dimensional cosmological $N$-body simulation, where we solve Newton's equations for gravity on an expanding cosmological background.

In a homogenous and isotropic Friedmann–Lemaître–Robertson–Walker universe, the Einstein Field equations reduce to the Friedmann equation

```math
H(t)^2 = \left(\frac{\dot{a}(t)}{a(t)}\right)^2 = H_0^2\left( \frac{\Omega_{0, r}}{a(t)^{4}} + \frac{\Omega_{0,m}}{ a(t)^{3}} + \frac{\Omega_{0,k}}{ a(t)^{2}} + \Omega_{0,\Lambda}\right)\,,
```

with the scale factor $a(t)$ describes the expansion history of our Universe as a function of time $t$ spanning the Big Bang $a=0$ to the current time $a(t_0)=1$. The Friedmann equation expression the Hubble factor $H=\frac{\dot{a}}{a}$ capturing the expansion of our Universe as a function of the current Hubble expansion rate $H_0$, the scale factor $a$, the current radiation density parameter $\Omega_{0,r}$, the current matter density parameter $\Omega_{0,m}$, the current curvature parameter $\Omega_{0,k}$ and the dark energy parameter $\Omega_{0,\Lambda}$.

In large-scale structure formation, the cosmic forms due to the gravitational collapse of small, close to Gaussian fluctuations. In an expanding universe, the gravitational potential $\phi$ is governed by the Poisson equation on an expanding background 

```math
\nabla^2 \phi = 4 \pi G \rho_u a^2 \delta\,,
```

with Newton's constant $G$, the mean matter density $\rho_u$ and the density perturbation defined by 

```math
1+\delta = \frac{\rho}{\rho_u}\,.
```

In $N$-body simulations, we follow the Lagrangian fluid formalism, where the motion of matter is governed by the Lagrangian map 

```math
\bm{x}_t(\bm{q}) = \bm{q} + \bm{s}_t(\bm{q})\,,
```

describing the motion of a mass element starting at point $\bm{q}$ as a function of time through the displacement map $\bm{s}_t$. The particles move according to the Euler equation, capturing Newton's second law for Lagrangian fluids 

```math
\frac{\partial a \bm{v}}{\partial t} = -\nabla \phi\,,
```

with the comoving velocity $\bm{v} = a \dot{\bm{x}_t}$. Working with the momentum $\bm{p} = a^2 \dot{x}=a \bm{v}$, and using the identity 

```math 
4 \pi G \rho_u = \frac{3 H_0^2 \Omega_{0,m}}{2 a^3}\,,
```

the Poisson and Euler equations assume the form 

```math 
a \nabla^2 \phi =\frac{3}{2} \Omega_{0,m} H_0^2 \delta\,,
``` 

and 

```math
\dot{p} = -\nabla \phi\,.
```

In these equations, the time dependence enters through the evolution of the scale factor. When the scale factor $a$ is a monotonically increasing function, it is natural to express the equation of motion for structure formation as 

```math
\partial_a \bm{x} = \frac{\bm{p}}{a^2 \dot{a}}\,,
``` 

and 

```math
\partial_a \bm{p} = - \frac{\nabla \phi}{\dot{a}}\,. 
```

This formalism saves us from having to solve the Friedmann equation. In this package, we solve this set of equations with a first-order Leap-Frog method, consisting of a drift and a kick operation.

## Density estimation
We use the Phase-Space Delaunay Tessellation Field Estimator (PS-DTFE). This package includes a two-dimensional implementation of the PS-DTFE method. For a detailed description, see the package [PhaseSpaceDTFE.jl](https://jfeldbrugge.github.io/PhaseSpaceDTFE.jl/stable/theory/), which implements the PS-DTFE method in 3D.
