# Box
"""
    The cosmological box, setting the volume of the box L^2 and the number of particles N^2.
"""
struct Box{T<:AbstractFloat, R<:AbstractRange{T}}
    N::Int
    L::T
    res::T
    kx::Matrix{T}
    ky::Matrix{T}
    m::T
    range::R
end

function Box(N::Integer, L::Real, m::Real = 1)
    T = float(typeof(L))
    k = repeat(T.(fftfreq(N) * N * 2. * π / L), 1, N)
    Range = range(0, T(L); length=N + 1)[begin:end - 1]
    return Box{T, typeof(Range)}(N, L, L / N, collect(k'), k, m, Range)
end

# Operations on scalar fields represented by 2D arrays
"Convolution."
convolve(a, b) = real(ifft(a .* fft(b)))

"Gradient of a scalar field"
gradient(f, box::Box) = stack([
    convolve(im * box.kx, f), 
    convolve(im * box.ky, f)])

"Hessian of a scalar field"
hessian(f, box::Box) = stack([
    convolve(-box.kx .* box.kx, f), 
    convolve(-box.kx .* box.ky, f),  
    convolve(-box.ky .* box.ky, f)])

"Laplace operator."
laplacian(f, box::Box) = convolve(-box.kx.^2 - box.ky.^2, f)

"Smooth an Array."
smooth(f, σ, box::Box) = convolve(exp.(- σ^2 * (box.kx.^2 + box.ky.^2) ./ 2.), f)

"Periodic index in 1:N, with a fast path for indices that are already in range or one period out."
@inline wrap1(i, N) = 1 <= i <= N ? i : (N < i <= 2N ? i - N : mod1(i, N))

"""
    Fourth-order finite difference gradient along dimension i (in grid units), multiplied by scale and written to out.
"""
function Grad2!(out, data, i, scale = 1)
    n1, n2 = size(data)
    c1, c2 = 1. / 12., 2. / 3.
    if i == 1
        @inbounds for k in 1:n2, j in 1:n1
            out[j, k] = scale * (c1 * data[wrap1(j - 2, n1), k] - c2 * data[wrap1(j - 1, n1), k] +
                                 c2 * data[wrap1(j + 1, n1), k] - c1 * data[wrap1(j + 2, n1), k])
        end
    else
        @inbounds for k in 1:n2
            km2, km1, kp1, kp2 = wrap1(k - 2, n2), wrap1(k - 1, n2), wrap1(k + 1, n2), wrap1(k + 2, n2)
            for j in 1:n1
                out[j, k] = scale * (c1 * data[j, km2] - c2 * data[j, km1] + c2 * data[j, kp1] - c1 * data[j, kp2])
            end
        end
    end
    return out
end

"""
    Finite diference gradient
"""
Grad2(data, i) = Grad2!(similar(data), data, i)


"""
    Two dimensional interpolation.
"""
@inline function Interp2D(data, x1::Real, x2::Real)
    i1, i2, j1, j2, m1, m2, n1, n2 = interp_weights(size(data, 1), x1, x2)
    @inbounds return data[i1, j1] * n1 * n2 + data[i2, j1] * m1 * n2 + data[i1, j2] * n1 * m2 + data[i2, j2] * m1 * m2
end

"Grid indices and bilinear weights of the point (x1, x2) in grid units on a periodic N×N grid."
@inline function interp_weights(N, x1, x2)
    i1, j1 = wrap1(floor(Int, x1) + 1, N), wrap1(floor(Int, x2) + 1, N)
    i2, j2 = wrap1(ceil(Int, x1) + 1, N), wrap1(ceil(Int, x2) + 1, N)
    m1, m2 = x1 >= 0 ? x1 - floor(x1) : mod(x1, 1.), x2 >= 0 ? x2 - floor(x2) : mod(x2, 1.)   # exact fractional part
    return i1, i2, j1, j2, m1, m2, 1. - m1, 1. - m2
end

"Two dimensional interpolation at the rows of x."
Interp2D(data, x::AbstractMatrix) = [Interp2D(data, x[i, 1], x[i, 2]) for i in axes(x, 1)]