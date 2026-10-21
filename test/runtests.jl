using Nbody2D, SpecialFunctions, Plots, VisualRegressionTests
using LinearAlgebra, Test

@testset "Nbody2D.jl" begin
    println("Start test.")

    let Ni = 2^5, L = 50., H0 = 70., OmegaM = 1.0, OmegaL = 0.
        global box = Box(Ni, L)
        global EdS = Cosmology(H0, OmegaM, OmegaL)
    end

    "Power spectrum."
    function P(k, ns, Rs, α)
        return α^2 * 4. * π * Rs^(2. + ns) / gamma(1. + ns / 2.) * k^(ns - 4.) * exp(-Rs^2 * k^2)
    end

    let ns = 2., Rs = 1., α = 1., seed = 1
        global phi = GRF(k -> P(k, ns, Rs, α), seed, box)
    end

    ai, af, δa = 0.02, 2.02, 0.02
    stateZeldovich = Zeldovich(phi, af, af, box, EdS)
    stateNbody = LeapFrog(phi, ai, af, δa, box, EdS)
    @test stateNbody.time ≈ af

    # savefig(plotMesh(stateNbody, box), "testFigure.png")

    @plottest plotMesh(stateNbody, box) "baseline/testFigure.png" 

    @plottest begin
        plot()
        plotMesh!(stateNbody, box) 
    end "baseline/testFigure.png" 

    # Phase-Space Delaunay Tesselation Field Estimator
    let depth = 5
        stateInitial = Zeldovich(phi, δa, δa, box, EdS)
        global estimatorZeldovich = PS_DTFE(stateInitial, stateZeldovich, depth, box)
        global estimatorNbody = PS_DTFE(stateInitial, stateNbody, depth, box)
        nothing
    end
    
    p = [25. + 0.1, 25. + 0.3]

    @test ps_density(p, estimatorZeldovich) ≈ 0.18003857736541345
    @test ps_density(p, estimatorNbody) ≈ 0.19985896583449622
    @test ps_velocity(p, estimatorZeldovich) ≈ [115.20330374480434 -69.08250191156435]
    @test numberOfStreams(p, estimatorZeldovich) == 1

    # Constrained Gaussian Random Field
    let ns = 2., Rs = 1., α = 1., seed = 1, σ = 1., ders = [0 0; 1 0; 0 1]
        global cgrf = cGRF(box, k -> P(k, ns, Rs, α), σ, [box.L / 2, box.L / 2], ders)
        c = [1, 0, 0]
        global (f, f_c) = constraintGRF(c, seed, cgrf, box)
        global f_var = varianceField(cgrf)
    end
    mc = measureConstraints(f_c, cgrf)
    @test abs(mc[1] - 1) < 10e-4
    @test abs(mc[2]) < 10e-4
    @test abs(mc[3]) < 10e-4

    # Constraints are the derivatives ∂xⁱ∂yʲ of the field at p, including odd orders
    let b = Box(64, 20.), ders = [1 0; 0 1; 2 1; 3 0], c = [0.3, -0.2, 0.1, -0.4]
        cg = cGRF(b, k -> P(k, 2., 1., 1.), 0., [b.L / 2, b.L / 2], ders)
        _, fc = constraintGRF(c, 2, cg, b)
        g, h = gradient(fc, b), hessian(fc, b)
        i0 = b.N ÷ 2 + 1
        @test [g[i0, i0, 1], g[i0, i0, 2]] ≈ c[1:2]
        @test gradient(h[:, :, 2], b)[i0, i0, 1] ≈ c[3]    # ∂x(∂x∂y) = ∂x²∂y
        @test gradient(h[:, :, 1], b)[i0, i0, 1] ≈ c[4]    # ∂x(∂x²) = ∂x³
    end

    gra = gradient(f, box)
    hes = hessian(f, box)
    lap = laplacian(f, box)
    smo = smooth(f, 1., box)

    @test gra[box.N ÷ 2, box.N ÷ 2, 1] ≈ -0.0036668974305495428
    @test hes[box.N ÷ 2, box.N ÷ 2 , 1] ≈ 0.8255046575403882
    @test lap[box.N ÷ 2, box.N ÷ 2] ≈ 1.1604751729595333
    @test smo[box.N ÷ 2, box.N ÷ 2] ≈ -2.060513408055121

    # Spectral derivatives: x varies along the second array index
    let b = Box(64, 2π)
        q = Nbody2D.LagrangianGrid(b)
        @test gradient(sin.(q[:,:,1]), b)[:,:,1] ≈ cos.(q[:,:,1])
        @test gradient(sin.(q[:,:,2]), b)[:,:,2] ≈ cos.(q[:,:,2])
    end

    # Linear growth: on large scales the N-body displacement follows the growing mode D(a)
    for cosmology in (Cosmology(70., 1.0, 0.), Cosmology(1., 0.3, 0.7))
        b = Box(32, 50.)
        ϕ = GRF(k -> k^-2 * exp(-36 * k^2), 1, b)
        ϕ .*= 1e-3 / maximum(abs, gradient(ϕ, b))
        state = LeapFrog(ϕ, 0.02, 0.5, 0.002, b, cosmology)
        s0 = gradient(ϕ, b)
        s = Nbody2D.unwrap_s.(state.position .- Nbody2D.LagrangianGrid(b), b.L)
        @test isapprox(-dot(s, s0) / dot(s0, s0), Nbody2D.growing_mode(0.5, cosmology), rtol = 0.02)
    end

    # The velocity gradient is applied as v_b + ∂_a v_b dx_a for a rotational flow v = (y, 0)
    let b = Box(16, 1.)
        q = Nbody2D.LagrangianGrid(b)
        v = zeros(16, 16, 2); v[:,:,1] = q[:,:,2]
        est = PS_DTFE(Nbody2D.State(0., q, zeros(16, 16, 2)), Nbody2D.State(0., q, v), 3, b)
        @test vec(ps_velocity([0.5013, 0.4021], est)) ≈ [0.4021, 0.]
    end

    # Second directional derivative v₁⋅∇(v₁⋅∇λ₁) with v₁ = x̂ and v₁⋅∇λ₁ = sin(x)
    let b = Box(64, 2π)
        q = Nbody2D.LagrangianGrid(b)
        v1 = zeros(64, 64, 2); v1[:,:,1] .= 1.
        λ1_grad = zeros(64, 64, 2); λ1_grad[:,:,1] = sin.(q[:,:,1])
        @test Nbody2D.directional_derivative_2(λ1_grad, v1, 10, 7, 10, 7, b.res) ≈ cos(q[10, 7, 1]) atol = 1e-6
    end

    # Caustic skeleton: D4 points are found once, and plotting works with the defaults
    let b = Box(64, 2π)
        skeleton = Caustics!(hessian(GRF(k -> k^-3 * exp(-k^2 / 4), 3, b), b), 0., b)
        pts = [vec(m) for m in skeleton.D4]
        @test !any(norm(pts[i] - pts[j]) < b.res / 4 for i in eachindex(pts), j in eachindex(pts) if i < j)
        @test plotSkeleton(skeleton) isa Plots.Plot
    end

    # PS-DTFE: every particle is a vertex of six lattice triangles, so translations and shears give the exact density,
    # also when particles wrap around the periodic box several times
    let b = Box(24, 1.)
        q = Nbody2D.LagrangianGrid(b)
        ρ̄ = b.N^2 / b.L^2
        pts = [[0.013, 0.977], [0.5013, 0.4987], [0.731, 0.264], [0.998, 0.002], [1.37, -0.41]]   # off the lattice lines
        @test all(count(==(i), Nbody2D.lattice_triangles(b.N)) == 6 for i in 1:b.N^2)
        for s_x in (fill(2.3, 24, 24), 1.7 .* sin.(2π .* q[:, :, 2])), wrap in (false, true)
            x = copy(q); x[:, :, 1] .+= s_x
            wrap && (x = mod.(x, b.L))
            est = PS_DTFE(Nbody2D.State(0., q, zero(q)), Nbody2D.State(0., x, zero(q)), 3, b)
            @test all(ps_density(p, est) ≈ ρ̄ for p in pts)
            @test all(numberOfStreams(p, est) == 1 for p in pts)
        end
    end

    # Generalized moment: ∫ k e^{-k²} e^{-σ²k²} dk / 2π = 1 / (4π(1 + σ²))
    @test generalizedMoment(0, k -> exp(-k^2), 1.) ≈ sqrt(1 / (8π))

    println("Finish test.")
end
