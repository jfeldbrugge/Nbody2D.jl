# Unconstrained Gaussian random field
"""
    Generate an unconstrained Gaussian random field.
"""
function GRF(P, seed, box)
    Random.seed!(seed)

    Pk = P.(sqrt.(box.kx.^2 + box.ky.^2))
    Pk[1, 1] = 0.
    return real(ifft(sqrt.(Pk) .* fft(randn(box.N, box.N)))) * box.N^(2/2) / box.L^(2/2)
end

"Generalzed moments."
function generalizedMoment(i, ns, Rs, α, σ) 
    integrand(k) = k^(2 * i +1) * P(k, ns, Rs, α) * exp(-σ^2 * k^2)/ (2. * π)
    return sqrt(quadgk(integrand, 0, Inf, rtol=1e-10)[1])
end

# Build a constrained Gaussian Random field with the Hoffmann-Ribak algorithm
"Constrained Gaussian random field object"
struct cGRF
    Pk::Matrix{AbstractFloat}
    Hh::Vector{Matrix{ComplexF64}}
    xi_i::Array{AbstractFloat, 3}
    xi_ij::Matrix{AbstractFloat}
    xi_ij_inv::Matrix{AbstractFloat}
    bases::Array{AbstractFloat, 3}
    
    "Construct the constrained Gaussian random field object."
    function cGRF(box::Box, P, σ, p = [box.L / 2, box.L / 2], ders = [0 0; 1 0; 0 1])
        ikx, iky = im * box.kx, im * box.ky
    
        Pk = P.(sqrt.(box.kx.^2 + box.ky.^2))
        Pk[1, 1] = 0.
    
        expIKX = exp.(-1im * (box.kx * p[1] + box.ky * p[2]))
        σK = exp.(- σ^2 * (box.kx.^2 + box.ky.^2) / 2.)    
        Hh = [σK .* expIKX .* ikx.^d[1] .* iky.^d[2] for d in eachrow(ders)]
    
        xi_i = stack([real(ifft(hi .* Pk)) for hi in Hh] * box.L^-2 * box.N^2 )
        xi_ij = [real(mean(conj.(hi) .* hj .* Pk)) for hi in Hh, hj in Hh] * box.L^-2 * box.N^2
        xi_ij_inv = inv(xi_ij)
    
        bases = zeros(size(ders, 1), box.N, box.N)
        for i in axes(bases, 1), j in axes(bases, 2), k in axes(bases, 3)
            for l in axes(xi_ij_inv, 2)
                bases[i, j, k] += xi_ij_inv[i, l] * xi_i[j, k, l]
            end
        end
    
        return new(Pk, Hh, xi_i, xi_ij, xi_ij_inv, bases)
    end
end

"Mean of an Array."
mean(x::Array) = sum(x) / length(x)

"The mean field given the constraints c."
function meanField(c, cgrf::cGRF) 
    bases = cgrf.bases
    meanF = zeros(size(bases, 2), size(bases, 3))

    for j in axes(meanF, 1), k in axes(meanF, 2)
        for i in axes(bases, 1)
            meanF[j, k] += bases[i, j, k] * c[i]
        end
    end

    return meanF
end

"The variance field of the constraints."
function varianceField(cgrf::cGRF) 
    xi_i = cgrf.xi_i
    xi_ij_inv = cgrf.xi_ij_inv

    variance = zeros(size(xi_i, 1), size(xi_i, 2))
    for k in axes(variance, 1), l in axes(variance, 2)
        for i in axes(xi_ij_inv, 1), j in axes(xi_ij_inv, 2)
            variance[k, l] += xi_i[k, l, i] * xi_ij_inv[i, j] * xi_i[k, l, j]
        end
    end
    return variance
end

"Measure the constraints."
measureConstraints(f, cgrf::cGRF) = [real(mean(conj.(h) .* fft(f))) for h in cgrf.Hh]

"Generate an unconstrained Gaussian random field."
function GRF(seed, cgrf::cGRF, box::Box)
    Random.seed!(seed)
    return real(ifft(sqrt.(cgrf.Pk) .* fft(randn(box.N, box.N)))) * box.N^(2/2) / box.L^(2/2)
end

"Generate an unconstrained, constrained Gaussian random field pair."
function constraintGRF(c, seed, cgrf::cGRF, box::Box)
    Random.seed!(seed)
    f = GRF(seed, cgrf, box)
    fc = f .- meanField(measureConstraints(f, cgrf), cgrf) .+ meanField(c, cgrf)
    return (f, fc)
end

