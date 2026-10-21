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
    GridIndex

Uniform G×G cell index of the triangles over the box [0, L)², in compressed form: the triangles whose bounding box
overlaps cell k are `tris[offsets[k]:offsets[k + 1] - 1]`, with cell k = i + (j - 1) G for x in cell i and y in cell j.
"""
struct GridIndex
    G::Int
    L::Float64
    offsets::Vector{Int}
    tris::Vector{Int}
end

function GridIndex(positions::Matrix{Float64}, corners::Matrix{Int}, L, G)
    T, h = size(corners, 1), L / G
    cell(x) = clamp(floor(Int, x / h) + 1, 1, G)
    lo, hi = zeros(Int, 2, T), zeros(Int, 2, T)
    @inbounds for t in 1:T, d in 1:2
        c1, c2, c3 = positions[corners[t, 1], d], positions[corners[t, 2], d], positions[corners[t, 3], d]
        lo[d, t], hi[d, t] = cell(min(c1, c2, c3)), cell(max(c1, c2, c3))
    end
    counts = zeros(Int, G^2)
    @inbounds for t in 1:T, j in lo[2, t]:hi[2, t], i in lo[1, t]:hi[1, t]
        counts[i + (j - 1) * G] += 1
    end
    offsets = vcat(1, 1 .+ cumsum(counts))
    next, tris = offsets[1:end - 1], zeros(Int, offsets[end] - 1)
    @inbounds for t in 1:T, j in lo[2, t]:hi[2, t], i in lo[1, t]:hi[1, t]
        k = i + (j - 1) * G
        tris[next[k]] = t
        next[k] += 1
    end
    return GridIndex(G, Float64(L), offsets, tris)
end

"Triangles whose bounding box overlaps the cell of the point (x, y) in [0, L)²."
@inline function cell_triangles(index::GridIndex, x, y)
    h = index.L / index.G
    k = clamp(floor(Int, x / h) + 1, 1, index.G) + (clamp(floor(Int, y / h) + 1, 1, index.G) - 1) * index.G
    return @inbounds @view index.tris[index.offsets[k]:index.offsets[k + 1] - 1]
end

"""
    PS_DTFE(state_initial, state_final, box; G = n)
    PS_DTFE(state_initial, state_final, depth, box)

Phase-space Delaunay tessellation field estimator. The particles are connected by the regular triangulation of their
Lagrangian lattice (`lattice_triangles`), every triangle is unwrapped across the periodic boundaries on its own
(`unwrap_triangle`), and triangles crossing the edge of the box are stored together with their periodic copies.

Triangles are located with a uniform G×G cell index (`GridIndex`); by default G is the size n of the particle lattice.
The second form, with a `depth`, is kept for compatibility: the depth is not used for the queries, and a bounding volume
hierarchy for debugging can be built with `BVH(estimator, depth)`.

Query points are taken modulo the box. A point lying exactly on an edge or vertex of the tessellation is counted in
every triangle that contains it. Fields can be evaluated at arbitrary points (`ps_density`, `numberOfStreams`,
`ps_velocity`, `velocitySum`) or, much faster, on regular grids (`ps_density_grid`, `numberOfStreams_grid`).
"""
struct PS_DTFE
    rho::Vector{Float64}            # density at the particles
    Drho::Matrix{Float64}           # density gradient on every triangle
    Dv::Array{Float64, 3}           # gradient of the momentum field on every triangle
    index::GridIndex                # cell index of the triangles
    simplices::Matrix{Int}          # particle indices of every triangle (periodic copies included)
    corners::Matrix{Int}            # rows of `positions` holding the corners of every triangle
    positions::Matrix{Float64}      # Eulerian corners, three rows per triangle
    origin::Matrix{Float64}         # 2 × T: first corner of every triangle
    Ainv::Matrix{Float64}           # 4 × T: inverse of the edge matrix [c2 - c1  c3 - c1] (NaN for zero-area triangles)
    original::BitVector             # true for the image of the first corner in [0, L)², false for a periodic copy
    velocities::Matrix{Float64}     # momenta p = a² dx/dt of the particles (as in State.momentum); the peculiar velocity is p / a
    positions_initial::Matrix{Float64}
    L::Float64

    function PS_DTFE(state_initial, state_final, box::Box; G::Integer = size(state_final.position, 1))
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

        # Gradients of the density and momentum on every triangle: with the edge matrix A = [c2 - c1  c3 - c1],
        # ∇f = A⁻ᵀ (f2 - f1, f3 - f1)
        Drho = zeros(size(triangles, 1), dim)
        Dv   = zeros(size(triangles, 1), dim, dim)
        @inbounds for t in axes(triangles, 1)
            p = tri_corners[t]
            a, c = p[2, 1] - p[1, 1], p[2, 2] - p[1, 2]
            b, d = p[3, 1] - p[1, 1], p[3, 2] - p[1, 2]
            Δ = a * d - b * c
            Δ == 0 && continue
            i1, i2, i3 = triangles[t, 1], triangles[t, 2], triangles[t, 3]
            f2, f3 = rho[i2] - rho[i1], rho[i3] - rho[i1]
            Drho[t, 1], Drho[t, 2] = (d * f2 - c * f3) / Δ, (a * f3 - b * f2) / Δ
            for k in 1:dim
                g2, g3 = velocities[i2, k] - velocities[i1, k], velocities[i3, k] - velocities[i1, k]
                Dv[t, 1, k], Dv[t, 2, k] = (d * g2 - c * g3) / Δ, (a * g3 - b * g2) / Δ
            end
        end

        # Keep every triangle at the image of its first corner, plus the periodic copies that overlap the box
        rows, shifts, original = Int[], Tuple{Float64, Float64}[], Bool[]
        @inbounds for t in axes(triangles, 1)
            p = tri_corners[t]
            x0, x1 = extrema((p[1, 1], p[2, 1], p[3, 1]))
            y0, y1 = extrema((p[1, 2], p[2, 2], p[3, 2]))
            for sx in (-L, 0., L), sy in (-L, 0., L)
                if x1 + sx >= 0 && x0 + sx <= L && y1 + sy >= 0 && y0 + sy <= L
                    push!(rows, t); push!(shifts, (sx, sy)); push!(original, sx == 0 && sy == 0)
                end
            end
        end
        positions = Matrix{Float64}(undef, 3length(rows), 2)
        @inbounds for (r, (t, (sx, sy))) in enumerate(zip(rows, shifts)), k in 1:3
            positions[3(r - 1) + k, 1] = tri_corners[t][k, 1] + sx
            positions[3(r - 1) + k, 2] = tri_corners[t][k, 2] + sy
        end
        corners = permutedims(reshape(collect(1:size(positions, 1)), 3, :))

        T = length(rows)
        origin, Ainv = zeros(2, T), zeros(4, T)
        @inbounds for t in 1:T
            r1, r2, r3 = corners[t, 1], corners[t, 2], corners[t, 3]
            x1, y1 = positions[r1, 1], positions[r1, 2]
            a, c = positions[r2, 1] - x1, positions[r2, 2] - y1
            b, d = positions[r3, 1] - x1, positions[r3, 2] - y1
            Δ = a * d - b * c
            origin[1, t], origin[2, t] = x1, y1
            Ainv[1, t], Ainv[2, t], Ainv[3, t], Ainv[4, t] = Δ == 0 ? (NaN, NaN, NaN, NaN) : (d / Δ, -b / Δ, -c / Δ, a / Δ)
        end

        return new(rho, Drho[rows, :], Dv[rows, :, :], GridIndex(positions, corners, L, G), triangles[rows, :], corners,
                   positions, origin, Ainv, BitVector(original), velocities, reshape(state_initial.position, :, 2), L)
    end
end

PS_DTFE(state_initial, state_final, depth::Integer, box::Box) = PS_DTFE(state_initial, state_final, box)

"""
    BVH(estimator::PS_DTFE, depth)

The bounding volume hierarchy of the estimator's triangles, as used before the cell index; for debugging, for instance
with `findBox` and `findIntersections(p, tree, estimator.positions, estimator.corners)`.
"""
BVH(estimator::PS_DTFE, depth::Integer) =
    BVH(1:size(estimator.corners, 1), [0 estimator.L; 0 estimator.L], estimator.positions, estimator.corners, 2depth)

"Offset (dx, dy) of (x, y) from the first corner of triangle t, and whether the point lies in the triangle."
@inline function locate(estimator::PS_DTFE, t, x, y)
    @inbounds begin
        dx, dy = x - estimator.origin[1, t], y - estimator.origin[2, t]
        b1 = estimator.Ainv[1, t] * dx + estimator.Ainv[2, t] * dy
        b2 = estimator.Ainv[3, t] * dx + estimator.Ainv[4, t] * dy
    end
    return dx, dy, b1 >= 0 && b2 >= 0 && b1 <= 1 && b2 <= 1 && b1 + b2 <= 1
end

"Density of triangle t at the offset (dx, dy) from its first corner."
@inline stream_density(estimator::PS_DTFE, t, dx, dy) =
    @inbounds estimator.rho[estimator.simplices[t, 1]] + estimator.Drho[t, 1] * dx + estimator.Drho[t, 2] * dy

"Triangles containing the point p (taken modulo the box)."
function streams(p, estimator::PS_DTFE)
    x, y = mod(p[1], estimator.L), mod(p[2], estimator.L)
    return [t for t in cell_triangles(estimator.index, x, y) if locate(estimator, t, x, y)[3]]
end

"Density, linearly interpolated on triangle s at the point p."
function stream_density(p, s, estimator::PS_DTFE)
    dx, dy, _ = locate(estimator, s, p[1], p[2])
    return stream_density(estimator, s, dx, dy)
end

"Momentum p = a² dx/dt, linearly interpolated on triangle s at the point p."
function stream_velocity(p, s, estimator::PS_DTFE)
    @inbounds return @view(estimator.velocities[estimator.simplices[s, 1], :]) +
        @view(estimator.Dv[s, :, :])' * (p .- @view(estimator.positions[estimator.corners[s, 1], :]))
end

"Density at the point p, summed over all streams."
function ps_density(p::AbstractVector{<:Real}, estimator::PS_DTFE)
    x, y = mod(p[1], estimator.L), mod(p[2], estimator.L)
    ρ = 0.
    for t in cell_triangles(estimator.index, x, y)
        dx, dy, inside = locate(estimator, t, x, y)
        inside && (ρ += stream_density(estimator, t, dx, dy))
    end
    return ρ
end

"""
    ps_velocity(p, estimator)

Momentum p = a² dx/dt of every stream at the point p, one row per stream; divide by the scale factor a for the
peculiar velocity.
"""
function ps_velocity(p::AbstractVector{<:Real}, estimator::PS_DTFE)
    p = mod.(p, estimator.L)
    simplexIndices = streams(p, estimator)
    vs = zeros(length(simplexIndices), length(p))
    for (i, s) in pairs(simplexIndices)
        vs[i, :] = stream_velocity(p, s, estimator)
    end
    return vs
end

"Number of streams at the point p."
function numberOfStreams(p::AbstractVector{<:Real}, estimator::PS_DTFE)
    x, y = mod(p[1], estimator.L), mod(p[2], estimator.L)
    return count(t -> locate(estimator, t, x, y)[3], cell_triangles(estimator.index, x, y))
end

"Stream-density weighted momentum p = a² dx/dt at the point p."
function velocitySum(p::AbstractVector{<:Real}, estimator::PS_DTFE)
    vs = ps_velocity(p, estimator)
    size(vs, 1) == 1 && return vs   # single-stream region

    # stream-density weighted sum of velocities in a multistream region
    p = mod.(p, estimator.L)
    densities = [stream_density(p, s, estimator) for s in streams(p, estimator)]
    return sum(densities .* vs, dims=1) / sum(densities)
end

"""
    rasterise(estimator, xs, ys)

Density and number of streams on the regular grid xs × ys (rows y, columns x), by looping once over the triangles and
over the grid points inside each of them. The grid may extend beyond the box: every triangle is placed at all its
periodic images that overlap the grid, so the result equals the point queries at mod.(p, L).
"""
function rasterise(estimator::PS_DTFE, xs::AbstractRange, ys::AbstractRange)
    L = estimator.L
    ρ, ns = zeros(length(ys), length(xs)), zeros(Int, length(ys), length(xs))
    (isempty(xs) || isempty(ys)) && return ρ, ns
    sx, sy = step(xs), step(ys)
    (sx > 0 && sy > 0) || throw(ArgumentError("the grid ranges must be increasing"))
    pixels(lo, hi, r, s) = (max(1, ceil(Int, (lo - first(r)) / s) + 1), min(length(r), floor(Int, (hi - first(r)) / s) + 1))
    @inbounds for t in axes(estimator.corners, 1)
        estimator.original[t] && !isnan(estimator.Ainv[1, t]) || continue
        r1, r2, r3 = estimator.corners[t, 1], estimator.corners[t, 2], estimator.corners[t, 3]
        x0, x1 = extrema((estimator.positions[r1, 1], estimator.positions[r2, 1], estimator.positions[r3, 1]))
        y0, y1 = extrema((estimator.positions[r1, 2], estimator.positions[r2, 2], estimator.positions[r3, 2]))
        for kx in floor(Int, (first(xs) - x1) / L):ceil(Int, (last(xs) - x0) / L), ky in floor(Int, (first(ys) - y1) / L):ceil(Int, (last(ys) - y0) / L)
            i0, i1 = pixels(x0 + kx * L, x1 + kx * L, xs, sx)
            j0, j1 = pixels(y0 + ky * L, y1 + ky * L, ys, sy)
            for i in i0:i1, j in j0:j1
                dx, dy, inside = locate(estimator, t, xs[i] - kx * L, ys[j] - ky * L)
                inside || continue
                ρ[j, i] += stream_density(estimator, t, dx, dy)
                ns[j, i] += 1
            end
        end
    end
    return ρ, ns
end

"""
    ps_density_grid(xs, ys, estimator)

Density on the regular grid xs × ys, as a matrix with rows y and columns x; equal to
`[ps_density([x, y], estimator) for y in ys, x in xs]`, but computed by rasterising the triangles.
"""
ps_density_grid(xs::AbstractRange, ys::AbstractRange, estimator::PS_DTFE) = rasterise(estimator, xs, ys)[1]

"""
    numberOfStreams_grid(xs, ys, estimator)

Number of streams on the regular grid xs × ys, as a matrix with rows y and columns x; equal to
`[numberOfStreams([x, y], estimator) for y in ys, x in xs]`, but computed by rasterising the triangles.
"""
numberOfStreams_grid(xs::AbstractRange, ys::AbstractRange, estimator::PS_DTFE) = rasterise(estimator, xs, ys)[2]
