# Box
"""
    The cosmological box, setting the volume of the box L^2 and the number of particles N^2.
"""
struct Box
    N 
    L 
    res
    kx 
    ky 
    m
    range

    function Box(N, L, m = 1) 
        k = repeat(fftfreq(N) * N * 2. * π / L, 1, N)
        Range = range(0, L; length=N + 1)[begin:end - 1]
        return new(N, L, L / N, k', k,m, Range)    
    end
end

# Operations on scalar fields represented by 2D arrays
"Convolution."
convolve(a, b) = real(ifft(conj.(a) .* fft(b)))

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

"""
    Finite diference gradient
"""
function Grad2(data, i) 
    if i == 1
        return 1. / 12. .* circshift(data, (2,0)) - 2. / 3. .* circshift(data, (1,0)) + 2. / 3. .* circshift(data, (-1,0)) -  1. / 12. .* circshift(data, (-2,0)) 
    else i == 2
        return 1. / 12. .* circshift(data, (0,2)) - 2. / 3. .* circshift(data, (0,1)) + 2. / 3. .* circshift(data, (0,-1)) -  1. / 12. .* circshift(data, (0,-2)) 
    end
end

"""
    Two dimensional interpolation.
"""
function Interp2D(data, x)
    N = size(data, 1)
    X1 = mod.(floor.(Int, x), N) .+ 1
    X2 = mod.(ceil.(Int, x),  N) .+ 1
    xm = mod.(x, 1.)
    xn = 1. .- xm

    f1 = [data[X1[i, 1], X1[i, 2]] for i in axes(X1, 1)]
    f2 = [data[X2[i, 1], X1[i, 2]] for i in axes(X1, 1)]
    f3 = [data[X1[i, 1], X2[i, 2]] for i in axes(X1, 1)]
    f4 = [data[X2[i, 1], X2[i, 2]] for i in axes(X1, 1)]

    return f1 .* xn[:,1] .* xn[:,2] + f2 .* xm[:,1] .* xn[:,2] + f3 .* xn[:,1] .* xm[:,2] + f4 .* xm[:,1] .* xm[:,2] 
end