# Eigenvalues
"Evaluate first eigenvalue and eigenvector"
function first_eigen(hess)
    H = [hess[1] hess[2];
         hess[2] hess[3]]
    F = eigen(H)
    i_max = argmax(F.values)
    (F.values[i_max], F.vectors[:, i_max])
end

function eigen123(hess)
    H = [hess[1] hess[2];
         hess[2] hess[3]]
    F = eigen(H)
    (F.values[2], F.values[1], F.vectors[:, 2], F.vectors[:, 1])
end

"Evaluate first eigenvalue and eigenvector field"
function eigenFields123(hess)
    Ni = size(hess, 1)
    λ1 = zeros(Ni, Ni)
    λ2 = zeros(Ni, Ni)
    v1 = zeros(Ni, Ni, 2)
    v2 = zeros(Ni, Ni, 2)

    for i in 1:Ni, j in 1:Ni
        eig = eigen123(hess[i,j,:])
        λ1[i, j]    = eig[1]
        λ2[i, j]    = eig[2]
        v1[i, j, :] = eig[3]
        v2[i, j, :] = eig[4]
    end
    return (λ1, λ2, v1, v2)
end

# Derivatives
"Finite difference order dx^2"
derivative(dm1, d1, Δ) = (d1 - dm1) / 2. / Δ

"Finite difference order dx^4"
derivative(dm2, dm1, d1, d2, Δ) = (4. / 6. * (d1 - dm1) -1. / 12. * ( d2 - dm2 )) / Δ

"Finite difference order dx^6"
derivative(dm3, dm2, dm1, d1, d2, d3, Δ) = (3. / 4. * (d1 - dm1) - 3. / 20. * (d2 - dm2)  + 1. / 60. * (d3 - dm3)) / Δ

"Evaluate the first order directional derivative v₁⋅∇λ₁"
function directional_derivative_1(λ1_grad, v1, i, j, i_ref, j_ref)
    return sign(dot(v1[i_ref, j_ref, :], v1[i, j, :])) * dot(v1[i, j, :], λ1_grad[i, j, :])
end

"Evaluate the second order directional derivative v₁⋅∇(v₁⋅∇λ₁)"
function directional_derivative_2(λ1_grad, v1, i, j, i_ref, j_ref, Δ)
    Ni = size(v1, 1)
    grad = [
        derivative(
            directional_derivative_1(λ1_grad, v1, mod(i - 3, 1:Ni), j, i_ref, j_ref),
            directional_derivative_1(λ1_grad, v1, mod(i - 2, 1:Ni), j, i_ref, j_ref),
            directional_derivative_1(λ1_grad, v1, mod(i - 1, 1:Ni), j, i_ref, j_ref),
            directional_derivative_1(λ1_grad, v1, mod(i + 1, 1:Ni), j, i_ref, j_ref),
            directional_derivative_1(λ1_grad, v1, mod(i + 2, 1:Ni), j, i_ref, j_ref),
            directional_derivative_1(λ1_grad, v1, mod(i + 3, 1:Ni), j, i_ref, j_ref), Δ),
        derivative(
            directional_derivative_1(λ1_grad, v1, i, mod(j - 3, 1:Ni), i_ref, j_ref),
            directional_derivative_1(λ1_grad, v1, i, mod(j - 2, 1:Ni), i_ref, j_ref),
            directional_derivative_1(λ1_grad, v1, i, mod(j - 1, 1:Ni), i_ref, j_ref),
            directional_derivative_1(λ1_grad, v1, i, mod(j + 1, 1:Ni), i_ref, j_ref),
            directional_derivative_1(λ1_grad, v1, i, mod(j + 2, 1:Ni), i_ref, j_ref),
            directional_derivative_1(λ1_grad, v1, i, mod(j + 3, 1:Ni), i_ref, j_ref), Δ)]
    
    return sign(dot(v1[i_ref, j_ref, :], v1[i, j, :])) * dot(v1[i, j, :], grad)
end

# Interpolation 
"Interpolate a vertex along a line"
VertexInterp(iso, p1, p2, v1, v2) = p1 + (iso - v1) / (v2 - v1) * (p2 - p1)

"Cubic interpolation of a regular grid"
function interpolation_cubic(q, data, Δ, g)
    xx = (q[1] - g.p[1,1]) / Δ
    yy = (q[2] - g.p[1,2]) / Δ
    
    return (
        data[1] * (1 - xx) * (1 - yy) +
        data[2] * xx       * (1 - yy) +
        data[3] * (1 - xx) * yy       +
        data[4] * xx       * yy)
end

# Caustic Skeleton
function extractCube(ID, data)
    i, j = ID[1], ID[2]
    return data[i:i+1,j:j+1][:]
end

"Caustic skeleton"
mutable struct CausticSkeleton
    A2::Vector{Matrix{Float64}}
    A3::Vector{Matrix{Float64}}
    A4::Vector{Matrix{Float64}}
    D4::Vector{Matrix{Float64}}

    function CausticSkeleton()
        return new(
            Vector{Matrix{Float64}}(undef, 0), 
            Vector{Matrix{Float64}}(undef, 0), 
            Vector{Matrix{Float64}}(undef, 0), 
            Vector{Matrix{Float64}}(undef, 0))
    end
end

"Gridcell with attributes"
struct gridcell
    p::Array
    eigen_1::Array
    inner::Array
    ID::Array

    function gridcell(i, j, λ1, λ1_grad, v1, Range)
        points = [
            Range[i]     Range[j]     
            Range[i + 1] Range[j]     
            Range[i]     Range[j + 1] 
            Range[i + 1] Range[j + 1]]
   
        eigen_1 = extractCube([i, j], λ1)

        inner = [directional_derivative_1(λ1_grad, v1, i + di, j + dj, i, j) for di in 0:1, dj in 0:1][:]

        return new(points, eigen_1, inner, [i, j])
    end
end

# Compute the caustics

# Fold caustics
"Compute the fold (A2) caustics"
function PolygoniseA2!(g::gridcell, skeleton::CausticSkeleton, iso, v)
    triindex = 0
    if g.eigen_1[v[1]] < iso triindex |= 1 end
    if g.eigen_1[v[2]] < iso triindex |= 2 end
    if g.eigen_1[v[3]] < iso triindex |= 4 end

    if     triindex == 1 || triindex == 6
        push!(skeleton.A2, hcat(
            VertexInterp(iso, g.p[v[1],:], g.p[v[3],:], g.eigen_1[v[1]], g.eigen_1[v[3]]),
            VertexInterp(iso, g.p[v[1],:], g.p[v[2],:], g.eigen_1[v[1]], g.eigen_1[v[2]]))')
    elseif triindex == 2 || triindex == 5
        push!(skeleton.A2, hcat(
            VertexInterp(iso, g.p[v[2],:], g.p[v[1],:], g.eigen_1[v[2]], g.eigen_1[v[1]]),
            VertexInterp(iso, g.p[v[2],:], g.p[v[3],:], g.eigen_1[v[2]], g.eigen_1[v[3]]))')
    elseif triindex == 3 || triindex == 4
        push!(skeleton.A2, hcat(
            VertexInterp(iso, g.p[v[1],:], g.p[v[3],:], g.eigen_1[v[1]], g.eigen_1[v[3]]),
            VertexInterp(iso, g.p[v[2],:], g.p[v[3],:], g.eigen_1[v[2]], g.eigen_1[v[3]]))')
    end
end

"Compute the fold (A2) caustics"
function PolygoniseA2!(g::gridcell, skeleton::CausticSkeleton, iso)
    PolygoniseA2!(g, skeleton, iso, [1, 2, 4])
    PolygoniseA2!(g, skeleton, iso, [1, 3, 4])
end

"Compute the cusp (A3), swallowtail (A4) and butterfly (A5) caustics"
function PolygoniseA!(g::gridcell, skeleton::CausticSkeleton, iso, λ1_grad, v1)
    # PolygoniseA!(g, skeleton, iso, [1, 2, 3], λ1_grad, v1)
    # PolygoniseA!(g, skeleton, iso, [1, 3, 4], λ1_grad, v1)
    # PolygoniseA!(g, skeleton, iso, [2, 3, 4], λ1_grad, v1)

    PolygoniseA!(g, skeleton, iso, [1, 2, 4], λ1_grad, v1)
    PolygoniseA!(g, skeleton, iso, [1, 3, 4], λ1_grad, v1)
end

"Compute the cusp (A3), swallowtail (A4) and butterfly (A5) caustics"
function PolygoniseA!(g::gridcell, skeleton::CausticSkeleton, iso, v, λ1_grad, v1)
    Δ = g.p[2,1] - g.p[1,1]

    triindex = 0
    if g.inner[v[1]] < 0 triindex |= 1 end
    if g.inner[v[2]] < 0 triindex |= 2 end
    if g.inner[v[3]] < 0 triindex |= 4 end

    function pushA4(p1, p2)
        i, j = g.ID[1], g.ID[2]

        innerD = [directional_derivative_2(λ1_grad, v1, i + di, j + dj, i, j, Δ) for di in 0:1, dj in 0:1][:]

        innerD1 = interpolation_cubic(p1, innerD, Δ, g)
        innerD2 = interpolation_cubic(p2, innerD, Δ, g)

        triindex = 0
        if innerD1 < 0 triindex |= 1 end
        if innerD2 < 0 triindex |= 2 end

        # Screen over lines
        if triindex == 1 || triindex == 2
            pp = VertexInterp(0., p1, p2, innerD1, innerD2)
            push!(skeleton.A4, pp')
        end        
    end

    function pushA3A4(v1, w1, v2, w2)
        p1 = VertexInterp(0, g.p[v1,:], g.p[w1,:], g.inner[v1], g.inner[w1])
        p2 = VertexInterp(0, g.p[v2,:], g.p[w2,:], g.inner[v2], g.inner[w2])

        eigen1 = interpolation_cubic(p1, g.eigen_1, Δ, g)
        eigen2 = interpolation_cubic(p2, g.eigen_1, Δ, g)

        triindex = 0
        if eigen1 > iso triindex |= 1 end
        if eigen2 > iso triindex |= 2 end

        # Screen over lines
        if triindex == 1
            pp2 = VertexInterp(iso, p1, p2, eigen1, eigen2)
            push!(skeleton.A3, hcat(p1, pp2)')
            pushA4(p1, pp2)
        elseif triindex == 2
            pp1 = VertexInterp(iso, p1, p2, eigen1, eigen2)
            push!(skeleton.A3, hcat(pp1, p2)')
            pushA4(pp1, p2)
        elseif triindex == 3
            push!(skeleton.A3, hcat(p1, p2)')
            pushA4(p1, p2)
        end
    end

    # Screen over triangles
    if     triindex == 1 || triindex == 6
        pushA3A4(v[1], v[3], v[1], v[2])
    elseif triindex == 2 || triindex == 5
        pushA3A4(v[2], v[1], v[2], v[3])
    elseif triindex == 3 || triindex == 4
        pushA3A4(v[1], v[3], v[2], v[3])
    end
end


"Compute the umbilic (D4) caustics"
function PolygoniseD4!(g::gridcell, skeleton::CausticSkeleton, iso, hess, c1, c2)
    i, j = g.ID[1], g.ID[2]

    val1 = extractCube(g.ID, c1) 
    val2 = extractCube(g.ID, c2)

    PolygoniseD!(g, skeleton, iso, [1, 2, 3], val1, val2, hess)
    PolygoniseD!(g, skeleton, iso, [1, 3, 4], val1, val2, hess)
    PolygoniseD!(g, skeleton, iso, [2, 3, 4], val1, val2, hess)
end

"Compute the umbilic (D4) caustics"
function PolygoniseD!(g::gridcell, skeleton::CausticSkeleton, iso, v, val1, val2, hess)
    Δ = g.p[2,1] - g.p[1,1]

    triindex = 0
    if val1[v[1]] < 0 triindex |= 1 end
    if val1[v[2]] < 0 triindex |= 2 end
    if val1[v[3]] < 0 triindex |= 4 end

    function pushD4(v1, w1, v2, w2)
        p1 = VertexInterp(0, g.p[v1,:], g.p[w1,:], val1[v1], val1[w1])
        p2 = VertexInterp(0, g.p[v2,:], g.p[w2,:], val1[v2], val1[w2])
        
        val2_1 = interpolation_cubic(p1, val2, Δ, g)
        val2_2 = interpolation_cubic(p2, val2, Δ, g)

        triindex = 0
        if val2_1 > 0 triindex |= 1 end
        if val2_2 > 0 triindex |= 2 end

        # Screen over lines
        if triindex == 1 || triindex == 2
            pp = VertexInterp(0, p1, p2, val2_1, val2_2)
            push!(skeleton.D4, pp')
        end
    end

    # Screen over triangles
    if     triindex == 1 || triindex == 6
        pushD4(v[1], v[3], v[1], v[2])
    elseif triindex == 2 || triindex == 5
        pushD4(v[2], v[1], v[2], v[3])
    elseif triindex == 3 || triindex == 4
        pushD4(v[1], v[3], v[2], v[3])
    end
end

"Compute the caustics"
function Caustics!(hess, iso, box::Box)
    skeleton = CausticSkeleton()

    (λ1, λ2, v1, v2) = eigenFields123(hess)
    
    N = box.N
    L = box.L
    Range = box.range
    λ1_grad = gradient(λ1, box)

    C1(h) = h[1] - h[3]
    C2(h) = h[2]

    c1 = [C1(hess[i,j,:]) for i in axes(hess, 1), j in axes(hess, 2)]
    c2 = [C2(hess[i,j,:]) for i in axes(hess, 1), j in axes(hess, 2)]

    for i in 1:N-1, j in 1:N-1
        g = gridcell(i, j, λ1, λ1_grad, v1, Range)
        if any(g.eigen_1 .> iso )
            PolygoniseA2!(g, skeleton, iso)
            PolygoniseA!(g,  skeleton, iso, λ1_grad, v1)
            PolygoniseD4!(g, skeleton, iso, hess, c1, c2)
        end
    end

    skeleton = transposeSkeleton(skeleton)

    return skeleton
end

function transposeSkeleton(skeleton::CausticSkeleton)
    skeleton_new = CausticSkeleton()

    skeleton_new.A2 = transposeSkeleton(skeleton.A2)
    skeleton_new.A3 = transposeSkeleton(skeleton.A3)
    skeleton_new.A4 = transposeSkeleton(skeleton.A4)
    skeleton_new.D4 = transposeSkeleton(skeleton.D4)

    return skeleton_new
end

function transposeSkeleton(simplices)
    return [[sim[:,2] sim[:,1]] for sim in simplices]
end

"Push the skeleton from Lagrangian to Eulerian space"
function Eulerian(skeleton::CausticSkeleton, state::State, box::Box)
    return Eulerian(skeleton, state.position .- LagrangianGrid(box), box::Box)
end

# Push from Lagrangian to Eulerian space
"Push the skeleton from Lagrangian to Eulerian space"
function Eulerian(skeleton::CausticSkeleton, s, box::Box)
    skeleton_Eulerian = CausticSkeleton()

    skeleton_Eulerian.A2 = Eulerian(skeleton.A2, s, box)
    skeleton_Eulerian.A3 = Eulerian(skeleton.A3, s, box)
    skeleton_Eulerian.A4 = Eulerian(skeleton.A4, s, box)
    skeleton_Eulerian.D4 = Eulerian(skeleton.D4, s, box)
    
    return skeleton_Eulerian
end

"Push the skeleton from Lagrangian to Eulerian space"
function Eulerian(simplices, s, box::Box)
    Range = box.range
    s1_itp = Interpolations.scale(interpolate(s[:,:,1], Interpolations.BSpline(Linear(Periodic()))), Range, Range)
    s2_itp = Interpolations.scale(interpolate(s[:,:,2], Interpolations.BSpline(Linear(Periodic()))), Range, Range)

    return [Eulerian(sim, s1_itp, s2_itp) for sim in simplices]
end

"Push the skeleton from Lagrangian to Eulerian space"
function Eulerian(M::Matrix, s1_itp, s2_itp)
    return stack([Eulerian(q, s1_itp, s2_itp) for q in eachrow(M)])'
end

"Push the skeleton from Lagrangian to Eulerian space"
Eulerian(q, s1_itp, s2_itp) = q .+ [
    s1_itp(q[2], q[1]), 
    s2_itp(q[2], q[1])]

function scaleSimplices!(skeleton::CausticSkeleton, L)
    scaleSimplices!(skeleton.A2, L)
    scaleSimplices!(skeleton.A3, L)
    scaleSimplices!(skeleton.A4, L)
    scaleSimplices!(skeleton.A5, L)
end

function scaleSimplices!(simplices, L)
    for i in eachindex(simplices)
        simplices[i] = simplices[i] ./ L 
    end
end