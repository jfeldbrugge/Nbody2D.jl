using Nbody2D, SpecialFunctions
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

    let depth = 5
        stateInitial = Zeldovich(phi, δa, δa, box, EdS)
        global estimatorZeldovich = PS_DTFE(stateInitial, stateZeldovich, depth, box)
        global estimatorNbody = PS_DTFE(stateInitial, stateNbody, depth, box)
        nothing
    end
    
    p = [25. + 0.1, 25. + 0.3]

    @test ps_density(p, estimatorZeldovich) ≈ 0.23577650370940562
    @test ps_density(p, estimatorNbody) ≈ 0.2905879496022735
    # @test ps_velocity(p, estimatorZeldovich) ≈ [-1.759013770242557 -0.09413530535089731]
    @test numberOfStreams(p, estimatorZeldovich) == 1

    println("Finish test.")
end
