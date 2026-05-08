```@meta
CurrentModule = Nbody2D
```

# Theory
<!-- The *Delaunay Tessellation Field Estimator* (DTFE) and its extension to phase space (*PS-DTFE*) are mathematical tools for the reconstruction of the density and velocity field of a discrete point set. We review here the derivations underlying these reconstruction methods. -->

```math
\nabla \rho = 
\begin{pmatrix}
p_{l_1}-p_{l_0}\\
\vdots\\
p_{l_d}-p_{l_0}\\
\end{pmatrix}^{-1}
\begin{pmatrix}
\rho_1-\rho_0\\
\vdots\\
\rho_d-\rho_0\\
\end{pmatrix}\,.
```
