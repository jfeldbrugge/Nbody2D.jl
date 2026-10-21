# Cosmology
"""
    The cosmological background in which the N-body simulation runs with the current Hubble rate H0, the matter density OmegaM, dark energy density OmegaL and the curvature component OmegaK.
"""
struct Cosmology{T<:AbstractFloat}
    H0::T
    OmegaM::T
    OmegaL::T
    OmegaK::T
    G::T
    factor::T
end

Cosmology(H0, OmegaM, OmegaL) = 
    Cosmology{float(promote_type(typeof(H0), typeof(OmegaM), typeof(OmegaL)))}(
        H0, 
        OmegaM, 
        OmegaL, 
        1. - OmegaM - OmegaL, 
        3. / 2. * OmegaM * H0^2., 
        growing_mode_norm(H0, OmegaM, OmegaL))

"""
    The derivative of the scale factor.
"""
adot(a, cosmology::Cosmology) = cosmology.H0 * a * sqrt(cosmology.OmegaL + cosmology.OmegaM * a^-3. + cosmology.OmegaK * a^-2.)

"""
    The growing mode.
"""
function growing_mode(a, cosmology::Cosmology)
    if a <= 0.001
        return a 
    else
        return cosmology.factor * adot(a, cosmology) / a * (quadgk(b -> adot(b, cosmology)^-3., 0.00001, a, rtol=1e-6)[1])
    end
end

"""
    The derivative of the growing mode with respect to the scale factor.
"""
function growing_mode_derivative(a, cosmology::Cosmology)
    if a <= 0.001
        return 1.
    else
        E = cosmology.OmegaL + cosmology.OmegaM * a^-3. + cosmology.OmegaK * a^-2.
        dH = cosmology.H0 * (-3. * cosmology.OmegaM * a^-4. - 2. * cosmology.OmegaK * a^-3.) / (2. * sqrt(E))
        I = quadgk(b -> adot(b, cosmology)^-3., 0.00001, a, rtol=1e-6)[1]
        return cosmology.factor * (dH * I + adot(a, cosmology) / a * adot(a, cosmology)^-3.)
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
    Cloud-in-Cell density estimator on an N×N grid with spacing res, written to delta.
"""
function CIC!(delta, X, N::Int, L, res)
    fill!(delta, 0)
    # particles are visited row by row, which fixes the order of the floating-point sums
    @inbounds for i in axes(X, 1), j in axes(X, 2)
        p1, p2     = mod(X[i, j, 1], L) / res, mod(X[i, j, 2], L) / res
        idx1, idx2 = floor(Int, p1), floor(Int, p2)
        f1, f2     = p1 - idx1, p2 - idx2
        a1, a2, b1, b2 = wrap1(idx1 + 1, N), wrap1(idx1 + 2, N), wrap1(idx2 + 1, N), wrap1(idx2 + 2, N)
        delta[a1, b1] += (1. - f1) * (1. - f2) 
        delta[a2, b1] += f1 * (1. - f2) 
        delta[a1, b2] += (1. - f1) * f2
        delta[a2, b2] += f1 * f2
    end
    return delta
end

"""
    Cloud-in-Cell density estimator.
"""
CIC(X, box::Box) = CIC!(zeros(eltype(X), box.N, box.N), X, box.N, box.L, box.res)

# Nbody
"""
    The state of the N-body simulation in phase-space.
"""
mutable struct State{T<:Real, A<:AbstractArray{T, 3}}
    time::T
    position::A
    momentum::A
end

State(time::Real, position::AbstractArray{T, 3}, momentum::AbstractArray{T, 3}) where {T} =
    State{T, typeof(position)}(time, position, momentum)

"""
    Lagrangian grid
"""
function LagrangianGrid(box::Box)
    q_y = repeat(box.range, 1, box.N)
    q_x = q_y'
    return stack((q_x, q_y))
end

"""
    Evaluate the Zel'dovich approximation, x = q - D(a) ∇ϕ, with the positions at a_pos and the momenta p = a² da/dt dx/da at a_vel.
"""
function Zeldovich(ϕ, a_pos, a_vel, box, cosmology::Cosmology)
    Dp_pos = growing_mode(a_pos, cosmology)
    dDp_vel = a_vel^2 * adot(a_vel, cosmology) * growing_mode_derivative(a_vel, cosmology)

    s = -gradient(ϕ, box)

    X = LagrangianGrid(box) .+ Dp_pos .* s
    P = dDp_vel .* s

    return State(a_pos, X, P)
end


# """
#     Evaluate the Zel'dovich approximation.
# """
# function Zeldovich(ϕ, a_pos, a_vel, box, cosmology::Cosmology)
#     u_x, u_y = -box.N / box.L .* Grad2(ϕ, 1), -box.N / box.L .* Grad2(ϕ, 2)

#     qRange = range(0. , box.L, box.N + 1)[2:end]
#     q_x = repeat(qRange, 1, box.N)
#     q_y = q_x'

#     Dp_pos = growing_mode(a_pos, cosmology)
#     Dp_vel = growing_mode(a_vel, cosmology)
#     X = stack((q_x .+ Dp_pos * u_x, q_y .+ Dp_pos * u_y))
#     P = stack((Dp_vel * u_x, Dp_vel * u_y))
    
#     return State(a_pos, X, P)
# end

"""
    Preallocated grids and FFT plans for the particle-mesh force on the force grid `box`.
"""
struct PMWorkspace{T, B<:Box{T}, PF, PI}
    box::B
    mass::T
    delta::Matrix{T}
    δk::Matrix{Complex{T}}
    k2::Matrix{T}
    ϕ::Matrix{T}
    u_x::Matrix{T}
    u_y::Matrix{T}
    rfft_plan::PF
    irfft_plan::PI
end

function PMWorkspace(box::Box{T}, mass) where {T}
    N = box.N
    delta = zeros(T, N, N)
    rfft_plan = plan_rfft(delta)
    δk = rfft_plan * delta
    irfft_plan = plan_irfft(δk, N)
    k2 = box.kx[1:N ÷ 2 + 1, :].^2 .+ box.ky[1:N ÷ 2 + 1, :].^2
    return PMWorkspace(box, T(mass), delta, δk, k2, zeros(T, N, N), zeros(T, N, N), zeros(T, N, N), rfft_plan, irfft_plan)
end

"""
    A drift step.
"""
function Drift(a, δa, P, cosmology::Cosmology) 
    return δa .* P / (a^2 * adot(a, cosmology))
end

"In-place drift step, X += Drift(a, δa, P, cosmology)."
function Drift!(X, a, δa, P, cosmology::Cosmology)
    s = a^2 * adot(a, cosmology)
    X .+= δa .* P ./ s
    return X
end

"""
    In-place kick step, P -= Kick(a, δa, X, ws.box, ws.mass, cosmology).
"""
function Kick!(P, a, δa, X, ws::PMWorkspace, cosmology::Cosmology)
    N, L, res = ws.box.N, ws.box.L, ws.box.res

    CIC!(ws.delta, X, N, L, res)
    ws.delta .= ws.delta .* ws.mass .- 1.
    mul!(ws.δk, ws.rfft_plan, ws.delta)
    ws.δk ./= ws.k2
    ws.δk[1, 1] = 0.
    mul!(ws.ϕ, ws.irfft_plan, ws.δk)
    ws.ϕ .= ws.ϕ .* cosmology.G ./ a

    Grad2!(ws.u_x, ws.ϕ, 1, N / L)
    Grad2!(ws.u_y, ws.ϕ, 2, N / L)

    c = δa / adot(a, cosmology)
    u_x, u_y = ws.u_x, ws.u_y
    @inbounds for j in axes(X, 2), i in axes(X, 1)
        i1, i2, j1, j2, m1, m2, n1, n2 = interp_weights(N, mod(X[i, j, 1], L) / res, mod(X[i, j, 2], L) / res)
        P[i, j, 1] += c * (u_x[i1, j1] * n1 * n2 + u_x[i2, j1] * m1 * n2 + u_x[i1, j2] * n1 * m2 + u_x[i2, j2] * m1 * m2)
        P[i, j, 2] += c * (u_y[i1, j1] * n1 * n2 + u_y[i2, j1] * m1 * n2 + u_y[i1, j2] * n1 * m2 + u_y[i2, j2] * m1 * m2)
    end
    return P
end

"""
    A kick step.
"""
Kick(a, δa, X, box, m, cosmology::Cosmology) = -Kick!(zero(X), a, δa, X, PMWorkspace(box, m), cosmology)

"""
    A leap frog step moving the positions from a to a + δa and the momenta from a + δa/2 to a + 3δa/2.
"""
function LeapFrogStep!(state::State, a, δa, ws::PMWorkspace, cosmology::Cosmology)
    Drift!(state.position, a + δa / 2., δa, state.momentum, cosmology)
    Kick!(state.momentum, a + δa, δa, state.position, ws, cosmology)
    state.time = a + δa
    return state
end

LeapFrogStep(state::State, a, δa, force_box::Box, mass, cosmology::Cosmology) =
    LeapFrogStep!(state, a, δa, PMWorkspace(force_box, mass), cosmology)

"Number of leap frog steps from ai to af."
nsteps(ai, af, δa) = round(Int, (af - ai) / δa)

"""
    The Leap Frog integrator. The returned positions are at state.time = af, the momenta half a step later at af + δa/2.
"""
function LeapFrog(ϕ, ai, af, δa, box::Box, cosmology::Cosmology)
    force_box = Box(2 * box.N, box.L)
    ws = PMWorkspace(force_box, (force_box.N / box.N)^2)
    state = Zeldovich(ϕ, ai, ai + δa / 2., box, cosmology)
    for i in 0:nsteps(ai, af, δa) - 1
        LeapFrogStep!(state, ai + i * δa, δa, ws, cosmology)
    end
    return state
end

"""
    The Leap Frog integrator, returning the state after every step (at ai + δa, ai + 2δa, ..., af).
"""
function LeapFrogSnapshots(ϕ, ai, af, δa, box::Box, cosmology::Cosmology)
    n = nsteps(ai, af, δa)
    force_box = Box(2 * box.N, box.L)
    ws = PMWorkspace(force_box, (force_box.N / box.N)^2)
    state = Zeldovich(ϕ, ai, ai + δa / 2., box, cosmology)
    snapshots = Vector{typeof(state)}(undef, n)
    for i in 1:n
        LeapFrogStep!(state, ai + (i - 1) * δa, δa, ws, cosmology)
        snapshots[i] = deepcopy(state)
    end
    return snapshots
end