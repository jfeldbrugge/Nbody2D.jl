# Cosmology
"""
    The cosmological background in which the N-body simulation runs with the current Hubble rate H0, the matter density OmegaM, dark energy density OmegaL and the curvature component OmegaK.
"""
struct Cosmology
    H0
    OmegaM
    OmegaL
    OmegaK
    G
    factor

    Cosmology(H0, OmegaM, OmegaL) = 
    new(H0, 
        OmegaM, 
        OmegaL, 
        1. - OmegaM - OmegaL, 
        3. / 2. * OmegaM * H0^2., 
        growing_mode_norm(H0, OmegaM, OmegaL))
end

"""
    The derivative of the scale factor.
"""
adot(a, cosmology::Cosmology) = cosmology.H0 * a * sqrt(cosmology.OmegaL + cosmology.OmegaM * a^-3. + cosmology.OmegaK * a^-2.)

"""
    The growing mode.
"""
function growing_mode(a::AbstractFloat, cosmology::Cosmology)
    if a <= 0.001
        return a 
    else
        return cosmology.factor * adot(a, cosmology) / a * (quadgk(b -> adot(b, cosmology)^-3., 0.00001, a, rtol=1e-6)[1])
    end
end

"""
    The normalization of the growing mode.
"""
function growing_mode_norm(H0, OmegaM, OmegaL)
    adot(a) = H0 * a * sqrt(OmegaL + OmegaM * a^-3. + (1. - OmegaM - OmegaL) * a^-2.)
    return 1. / (adot(1.) * quadgk(b -> adot(b)^-3., 0.00001, 1., rtol=1e-6)[1])
end

"""
    Cloud-in-Cell density estimator.
"""
function CIC(X, box)
    delta = zeros(box.N, box.N)
    pos = mod.(hcat(vec(X[:,:,1]'), vec(X[:,:,2]')), box.L) ./ box.res
    for i in 1:size(pos, 1)
        idx1, idx2 = floor(Int, pos[i, 1]), floor(Int, pos[i, 2])
        f1, f2     = pos[i, 1] - idx1, pos[i, 2] - idx2
        delta[mod1(idx1 + 1, box.N), mod1(idx2 + 1, box.N)] += (1. - f1) * (1. - f2) 
        delta[mod1(idx1 + 2, box.N), mod1(idx2 + 1, box.N)] += f1 * (1. - f2) 
        delta[mod1(idx1 + 1, box.N), mod1(idx2 + 2, box.N)] += (1. - f1) * f2
        delta[mod1(idx1 + 2, box.N), mod1(idx2 + 2, box.N)] += f1 * f2
    end
    return delta
end

# Nbody
"""
    The state of the N-body simulation in phase-space.
"""
mutable struct State 
    time 
    position
    momentum
end

"""
    Evaluate the Zel'dovich approximation.
"""
function Zeldovich(phi, a_pos, a_vel, box, cosmology::Cosmology)
    u_x, u_y = -box.N / box.L .* Grad2(phi, 1), -box.N / box.L .* Grad2(phi, 2)

    qRange = range(0. , box.L, box.N + 1)[2:end]
    q_x = repeat(qRange, 1, box.N)
    q_y = q_x'

    Dp_pos = growing_mode(a_pos, cosmology)
    Dp_vel = growing_mode(a_vel, cosmology)
    X = stack((q_x .+ Dp_pos * u_x, q_y .+ Dp_pos * u_y))
    P = stack((Dp_vel * u_x, Dp_vel * u_y))
    
    return State(a_pos, X, P)
end

"""
    A drift step.
"""
function Drift(a, δa, P, cosmology::Cosmology) 
    return δa .* P / (a^2 * adot(a, cosmology))
end

"""
    A kick step.
"""
function Kick(a, δa, X, box, m, cosmology::Cosmology)
    delta = CIC(X, box) .* m .- 1.
    phi_f = fft(delta) ./ (box.kx.^2 .+ box.ky.^2)
    phi_f[1,1] = 0.
    phi = real(ifft(phi_f)) .* cosmology.G ./ a

    u_x, u_y = box.N / box.L .* Grad2(phi, 1), box.N / box.L .* Grad2(phi, 2)
    
    acc_x = Interp2D(u_x, mod.(hcat(vec(X[:,:,1]'), vec(X[:,:,2]')), box.L) ./ box.res)
    acc_y = Interp2D(u_y, mod.(hcat(vec(X[:,:,1]'), vec(X[:,:,2]')), box.L) ./ box.res)
    acc = -stack((reshape(acc_x, (box.N ÷ 2, box.N ÷ 2))', reshape(acc_y, (box.N ÷ 2, box.N ÷ 2))'))
    return δa / adot(a, cosmology) .* acc
end

"""
    The Leap Frog integrator.
"""
function LeapFrog(phi, ai, af, δa, box, cosmology::Cosmology)
    force_box = Box(2 * box.N, box.L)
    mass = (force_box.N / box.N)^2
    state = Zeldovich(phi, ai, ai + δa / 2., box, cosmology)
    for a in ai:δa:af
        state.position += Drift(a, δa, state.momentum, cosmology)
        state.momentum -= Kick(a + δa / 2., δa, state.position, force_box, mass, cosmology)
    end
    return state
end