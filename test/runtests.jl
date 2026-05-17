using Nbody2D, SpecialFunctions, Plots, VisualRegressionTests
using Test

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

    @test ps_density(p, estimatorZeldovich) ≈ 0.15695678684365458
    @test ps_density(p, estimatorNbody) ≈ 0.1885269531128359
    # @test ps_velocity(p, estimatorZeldovich) ≈ [-1.759013770242557 -0.09413530535089731]
    @test numberOfStreams(p, estimatorZeldovich) == 1

    # Constrained Gaussian Random Field
    let ns = 2., Rs = 1., α = 1., seed = 1, σ = 1., ders = [0 0; 1 0; 0 1]
        global cgrf = cGRF(box, k -> P(k, ns, Rs, α), σ, ders)
        c = [1, 0, 0]
        global (f, f_c) = constraintGRF(c, seed, cgrf, box)
        global f_var = varianceField(cgrf)
    end
    mc = measureConstraints(f_c, cgrf)
    @test abs(mc[1] - 1) < 10e-4
    @test abs(mc[2]) < 10e-4
    @test abs(mc[3]) < 10e-4

    gra = gradient(f, box)
    hes = hessian(f, box)
    lap = laplacian(f, box)
    smo = smooth(f, 1., box)

    @test gra[box.N ÷ 2, box.N ÷ 2, 1] ≈ 0.0036668974305495428
    @test hes[box.N ÷ 2, box.N ÷ 2 , 1] ≈ 0.8255046575403882
    @test lap[box.N ÷ 2, box.N ÷ 2] ≈ 1.1604751729595333
    @test smo[box.N ÷ 2, box.N ÷ 2] ≈ -2.060513408055121

    println("Finish test.")
end
