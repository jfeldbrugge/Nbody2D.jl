struct BVH
    data::Union{Vector{Int}, Nothing}
    box::Matrix{Float64}
    depth::Int
    left::Union{BVH, Nothing}
    right::Union{BVH, Nothing}
end

function BVH(data, box::Matrix{Float64}, points::Matrix{Float64}, simplices, depth::Int)
    if depth != 0
        dim = (depth % size(box, 1)) + 1
        @inbounds mins = vec(minimum(@view(points[@view(simplices[data, :]), dim]), dims = 2))
        @inbounds maxs = vec(maximum(@view(points[@view(simplices[data, :]), dim]), dims = 2))
        
        @inbounds div = sum(@view(box[dim,:])) / length(@view(box[dim,:]))
        L = mins .<= div
        R = maxs .>= div

        Lbox = copy(box)
        @inbounds Lbox[dim,:] = [Lbox[dim, 1], div]
        Rbox = copy(box)
        @inbounds Rbox[dim,:] = [div, Rbox[dim, 2]]

        @inbounds BVH(nothing, box, depth, BVH(data[L], Lbox, points, simplices, depth - 1), 
                                           BVH(data[R], Rbox, points, simplices, depth - 1))
    else 
        BVH(data, box, depth, nothing, nothing)
    end
end

function findCandidateSimplices(p::Vector{Float64}, BVH_tree::BVH)
    dim = (BVH_tree.depth % size(BVH_tree.box, 1)) + 1
    if BVH_tree.depth == 0
        return BVH_tree.data
    elseif p[dim] < BVH_tree.left.box[dim, 2]
        return findCandidateSimplices(p, BVH_tree.left)
    elseif p[dim] > BVH_tree.left.box[dim, 2]
        return findCandidateSimplices(p, BVH_tree.right)
    else
        return union(findCandidateSimplices(p, BVH_tree.left),
                     findCandidateSimplices(p, BVH_tree.right))
    end
end

function findBox(p::Vector{Float64}, BVH_tree::BVH)
    dim = (BVH_tree.depth % size(BVH_tree.box, 1)) + 1
    if BVH_tree.depth == 0
        println(BVH_tree.box)
    elseif p[dim] < BVH_tree.left.box[dim, 2]
        return findBox(p, BVH_tree.left)
    elseif p[dim] > BVH_tree.left.box[dim, 2]
        return findBox(p, BVH_tree.right)
    else
        findBox(p, BVH_tree.left)
        findBox(p, BVH_tree.right)
    end
end

function intersection(p::Vector{Float64}, simplex)
    @inbounds A = @view(simplex[2:end,:])' .- @view(simplex[1,:])
    det(A) == 0 && return false        # zero-area triangle
    barry = inv(A) * (p .- @view(simplex[1,:]))
    return all(barry .>= 0) & all(barry .<= 1) & (sum(barry) <= 1.)
end

function findIntersections(p::Vector{Float64}, BVH_tree::BVH, points, simplices)
    candidates = findCandidateSimplices(p, BVH_tree)
    @inbounds filter(i -> intersection(p, @view(points[@view(simplices[i,:]),:])), candidates)
end

volume(sim, points) = @inbounds abs(det(@view(points[@view(sim[2:end]),:])' .- @view(points[sim[1],:]))) /  factorial(size(points, 2))

"""
    unwrap_s(s, L)

Cosmological N-body simulations normally work with periodic boundary conditions. This function undoes the wrapping of particles around the box in the displacement field, assuming the particles did not travel more than half the size of the box.
"""
unwrap_s(s, L) = mod((s + L / 2), L) - L / 2

"""
    lattice_triangles(N)

Triangulation of the periodic N×N Lagrangian lattice. Every cell is split along the same diagonal, into
(i, j)–(i+1, j)–(i+1, j+1) and (i, j)–(i+1, j+1)–(i, j+1), so that every particle is a vertex of exactly six triangles.
Particle (i, j) has index i + (j - 1) N, matching reshape(state.position, :, 2).
"""
function lattice_triangles(N)
    idx(i, j) = mod1(i, N) + (mod1(j, N) - 1) * N
    simplices = Matrix{Int}(undef, 2N^2, 3)
    t = 0
    for j in 1:N, i in 1:N
        a, b, c, d = idx(i, j), idx(i + 1, j), idx(i, j + 1), idx(i + 1, j + 1)
        simplices[t += 1, :] .= (a, b, d)
        simplices[t += 1, :] .= (a, d, c)
    end
    return simplices
end

"""
    unwrap_triangle(x, sim, L)

Eulerian corners of the triangle with particle indices `sim` (rows of x). The second and third corners are placed at the
periodic image nearest to the first corner, and the first corner is shifted into [0, L)². This only assumes that the
edges of the triangle are shorter than L/2, so particles may have wrapped around the box any number of times.
"""
function unwrap_triangle(x, sim, L)
    corners = Matrix{Float64}(undef, 3, 2)
    @inbounds for d in 1:2
        x1 = x[sim[1], d]
        shift = fld(x1, L) * L
        corners[1, d] = x1 - shift
        for k in 2:3
            corners[k, d] = x1 + unwrap_s(x[sim[k], d] - x1, L) - shift
        end
    end
    return corners
end

"""
    PS_DTFE(state_initial, state_final, depth, box)

Phase-space Delaunay tessellation field estimator. The particles are connected by the regular triangulation of their
Lagrangian lattice (`lattice_triangles`), every triangle is unwrapped across the periodic boundaries on its own
(`unwrap_triangle`), and triangles crossing the edge of the box are stored together with their periodic copies.
Query points are taken modulo the box. A point lying exactly on an edge or vertex of the tessellation is counted in
every triangle that contains it.
"""
struct PS_DTFE
    rho::Vector{Float64}            # density at the particles
    Drho::Matrix{Float64}           # density gradient on every triangle
    Dv::Array{Float64, 3}           # velocity gradient on every triangle
    tree::BVH
    simplices::Matrix{Int}          # particle indices of every triangle (periodic copies included)
    corners::Matrix{Int}            # rows of `positions` holding the corners of every triangle
    positions::Matrix{Float64}      # Eulerian corners, three rows per triangle
    velocities::Matrix{Float64}     # velocities of the particles
    positions_initial::Matrix{Float64}
    L::Float64

    function PS_DTFE(state_initial, state_final, depth, box)
        L, m, dim = box.L, box.m, 2
        n = size(state_final.position, 1)
        size(state_final.position, 2) == n || throw(ArgumentError("PS_DTFE expects the particles on an n×n Lagrangian lattice"))
        x = reshape(state_final.position, :, 2)
        velocities = reshape(state_final.momentum, :, 2)
        triangles = lattice_triangles(n)
        tri_corners = [unwrap_triangle(x, @view(triangles[t, :]), L) for t in axes(triangles, 1)]

        rho = zeros(size(x, 1))
        @inbounds for t in axes(triangles, 1)
            vol = volume(1:3, tri_corners[t])
            for index in @view(triangles[t, :])
                rho[index] += vol
            end
        end
        rho = (1. + dim) * m ./ rho

        Drho = zeros(size(triangles, 1), dim)
        Dv   = zeros(size(triangles, 1), dim, dim)
        @inbounds for t in axes(triangles, 1)
            sim = @view(triangles[t, :])
            p = tri_corners[t]
            A = @view(p[2:end, :])' .- @view(p[1, :])
            det(A) == 0 && continue
            A_inv = inv(A)'
            Drho[t, :] = A_inv * (rho[sim[2:end]] .- rho[sim[1]])
            Dv[t, :, :] = A_inv * ((velocities[sim[2:end], :]' .- velocities[sim[1], :])')
        end

        # Keep every triangle at the image of its first corner, plus the periodic copies that overlap the box
        rows = Int[]; positions = Matrix{Float64}[]
        for t in axes(triangles, 1), sx in (-L, 0., L), sy in (-L, 0., L)
            c = tri_corners[t] .+ [sx sy]
            if maximum(@view(c[:, 1])) >= 0 && minimum(@view(c[:, 1])) <= L && maximum(@view(c[:, 2])) >= 0 && minimum(@view(c[:, 2])) <= L
                push!(rows, t); push!(positions, c)
            end
        end
        positions = reduce(vcat, positions)
        corners = permutedims(reshape(collect(1:size(positions, 1)), 3, :))

        BVH_tree = BVH(1:length(rows), [0 L; 0 L], positions, corners, depth * dim)
        return new(rho, Drho[rows, :], Dv[rows, :, :], BVH_tree, triangles[rows, :], corners, positions, velocities,
                   reshape(state_initial.position, :, 2), L)
    end
end

"Triangles containing the point p (taken modulo the box)."
streams(p, estimator::PS_DTFE) = findIntersections(mod.(p, estimator.L), estimator.tree, estimator.positions, estimator.corners)

"Density, linearly interpolated on triangle s at the point p."
function stream_density(p, s, estimator::PS_DTFE)
    @inbounds return estimator.rho[estimator.simplices[s, 1]] +
        @view(estimator.Drho[s, :])' * (p .- @view(estimator.positions[estimator.corners[s, 1], :]))
end

"Velocity, linearly interpolated on triangle s at the point p."
function stream_velocity(p, s, estimator::PS_DTFE)
    @inbounds return @view(estimator.velocities[estimator.simplices[s, 1], :]) +
        @view(estimator.Dv[s, :, :])' * (p .- @view(estimator.positions[estimator.corners[s, 1], :]))
end

function ps_density(p::Vector{Float64}, estimator::PS_DTFE)
    p = mod.(p, estimator.L)
    return sum((stream_density(p, s, estimator) for s in streams(p, estimator)), init = 0.)
end

function ps_velocity(p::Vector{Float64}, estimator::PS_DTFE)
    p = mod.(p, estimator.L)
    simplexIndices = streams(p, estimator)
    vs = zeros(length(simplexIndices), length(p))
    for (i, s) in pairs(simplexIndices)
        vs[i, :] = stream_velocity(p, s, estimator)
    end
    return vs
end

numberOfStreams(p::Vector{Float64}, estimator::PS_DTFE) = length(streams(p, estimator))

function velocitySum(p::Vector{Float64}, estimator::PS_DTFE)
    vs = ps_velocity(p, estimator)
    size(vs, 1) == 1 && return vs   # single-stream region

    # stream-density weighted sum of velocities in a multistream region
    p = mod.(p, estimator.L)
    densities = [stream_density(p, s, estimator) for s in streams(p, estimator)]
    return sum(densities .* vs, dims=1) / sum(densities)
end
