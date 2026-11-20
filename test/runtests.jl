# Regression tests for the shields. Run from the repository root with
#     julia --project=. test/runtests.jl
# The reference values pin what each shield computes, so a change in them means a shield now computes something
# different. If that is intended, update the values and say why in the commit.
using Test, SHA, LazySets, LinearAlgebra, Random
include(joinpath(@__DIR__, "..", "src", "ReACTShields.jl"))
using .ReACTShields

# First 12 hex digits of the SHA-1 of the allowed bits, so a single changed cell is caught.
digest(a) = bytes2hex(sha1(reinterpret(UInt8, vec(a))))[1:12]

@testset "ReACTShields" begin
    @testset "bouncing ball through modal, granularity 0.02" begin
        sys, p = loadBouncingBallShieldedWithInput()
        msys = ModalHybridSystem(sys)
        g, δ⁻, δ⁺ = 0.02, 0.001, 0.032
        outside = [.!ReACTShields.guard_masks(Grid(sys.statespace, g), act.guard)[1] for act in sys.Act]
        # The Euclidean shields removed on 2026-11-20 gave 732669 dead cells computing every cell (exits bad) and 231284
        # online (pred_fast, which never computes cells far from a change, as if leaving were allowed). Modal has 9 cells
        # less dead in both, kept alive by their blocks, and forbids action 2 in 86 cells pred_fast never computed.
        for (exits, dead, sha) in ((true, 732660, "b81b1916d8f4"), (false, 231275, "f0c3da88e97d"))
            @testset "leaving the state space bad: $exits" begin
                run(online) = ReACTed_reachable_cell_modal(msys, p, g, δ⁻, δ⁺; online, exit_bad_for_flow=exits, exit_bad_for_actions=exits, verbose=false)[2]
                allowed = run(false)
                @test size(allowed) == (1500, 750, 1)
                @test count(==(0), allowed) == dead
                @test digest(allowed) == sha
                @test run(true) == allowed
                # Outside its guard an action does nothing, so it is allowed exactly where waiting is.
                for ai in eachindex(sys.Act)
                    @test all(i -> isodd(allowed[i, 1]) == isodd(allowed[i, 1] >> ai), findall(outside[ai]))
                end
            end
        end
        @testset "adaptive, granularity 0.04 to 0.02, leaving allowed" begin
            # At 0.04 the ball is all dead when leaving is bad, so there is nothing to refine.
            run(online, rounds=typemax(Int)) = ReACTed_reachable_cell_modal_adaptive(msys, p, 2g, 1, δ⁻, δ⁺; online, max_rounds=rounds, exit_bad_for_flow=false, exit_bad_for_actions=false, verbose=false)[2]
            allowed = run(false)
            @test count(==(0), allowed) == 523210
            @test digest(allowed) == "69bd362e1f01"
            @test run(true) == allowed
            coarse = run(true, 0)
            @test all(allowed .& coarse .== coarse) # refining only allows more
        end
    end

    @testset "modal: oil pump, granularity 0.2, δ 0.1 to 0.2" begin
        sys, p = loadOilPumpModal(; period=0.2)
        _, allowed, _, _ = ReACTed_reachable_cell_modal(sys, p, 0.2, 0.1, 0.2; online=true, verbose=false)
        @test size(allowed) == (100, 107, 10, 4)
        # 57515 and 920081310bff before cells kept only the successors their block reaches too (10 cells less dead)
        @test count(==(0), allowed) == 57505
        @test digest(allowed) == "655555fc2d5f"
        @test ReACTed_reachable_cell_modal(sys, p, 0.2, 0.1, 0.2; online=false, verbose=false)[2] == allowed
        # The pump cannot switch during the latency, so switching does nothing there and is allowed exactly when waiting is.
        for layer in (OILPUMP_OFF_LOCKED, OILPUMP_ON_LOCKED)
            @test all(b -> b == 0x00 || b == 0x07, allowed[:, :, :, layer])
        end
    end

    @testset "modal: cruise control (_t), sensor range 50, granularity 1, δ 0.25 to 1" begin
        sys, p = loadCruiseControlModal_t(; sensor_range=50.0)
        run(online) = ReACTed_reachable_cell_modal(sys, p, 1.0, 0.25, 1.0; online, act_needs_safe_flow=false, verbose=false)[2]
        allowed = run(false)
        @test size(allowed) == (50, 30, 28, 1)
        @test count(==(0), allowed) == 20362
        @test digest(allowed) == "fb1279fa54ed"
        # Online the shield only skips cells that cannot become unsafe, so it gives the same result.
        @test run(true) == allowed
    end

    @testset "modal online gives the offline result" begin
        # Random one-location linear systems on [-5, 5]²: mostly contracting dynamics with a random rotation or shear and
        # input, x >= 4 unsafe, an edge from y <= -4 resetting y to 4, and an action on x >= 0 that flips and scales x.
        rng = MersenneTwister(7)
        for _ in 1:8
            A = round.(-0.4 .* I(2) .+ 0.6 .* (rand(rng, 2, 2) .- 0.5); digits=2)
            U = Zonotope(round.(rand(rng, 2) .- 0.5; digits=2), diagm(round.(0.2 .* rand(rng, 2); digits=2)))
            edge = Edge(1, HPolyhedron([LazySets.HalfSpace([0.0, 1.0], -4.0)]), [1.0 0.0; 0.0 0.0], [0.0, 4.0])
            act = Action(1, HPolyhedron([LazySets.HalfSpace([-1.0, 0.0], 0.0)]), [round(-rand(rng) - 0.2; digits=2) 0.0; 0.0 1.0], [0.0, 0.0])
            rsys = ModalHybridSystem(Hyperrectangle([0.0, 0.0], [5.0, 5.0]), [HPolyhedron([LazySets.HalfSpace([-1.0, 0.0], -4.0)])],
                [ModalLocation(1, nothing, Matrix{Float64}(A), U, [edge], [ModalAction(act, 1)])])
            for exits in (true, false)
                run(online) = ReACTed_reachable_cell_modal(rsys, 0.5, 0.25, 0.125, 0.5; online, exit_bad_for_flow=exits, exit_bad_for_actions=exits, verbose=false)[2]
                @test run(true) == run(false)
            end
        end
    end

    @testset "modal adaptive: refining only allows more, online gives the offline result" begin
        rng = MersenneTwister(11)
        for _ in 1:4
            A = round.(-0.4 .* I(2) .+ 0.6 .* (rand(rng, 2, 2) .- 0.5); digits=2)
            U = Zonotope(round.(rand(rng, 2) .- 0.5; digits=2), diagm(round.(0.2 .* rand(rng, 2); digits=2)))
            edge = Edge(1, HPolyhedron([LazySets.HalfSpace([0.0, 1.0], -4.0)]), [1.0 0.0; 0.0 0.0], [0.0, 4.0])
            act = Action(1, HPolyhedron([LazySets.HalfSpace([-1.0, 0.0], 0.0)]), [round(-rand(rng) - 0.2; digits=2) 0.0; 0.0 1.0], [0.0, 0.0])
            rsys = ModalHybridSystem(Hyperrectangle([0.0, 0.0], [5.0, 5.0]), [HPolyhedron([LazySets.HalfSpace([-1.0, 0.0], -4.0)])],
                [ModalLocation(1, nothing, Matrix{Float64}(A), U, [edge], [ModalAction(act, 1)])])
            run(rounds; online=false) = ReACTed_reachable_cell_modal_adaptive(rsys, 0.5, 0.5, 2, 0.125, 0.5; online, max_rounds=rounds, verbose=false)[2]
            coarse, refined = run(0), run(typemax(Int))
            @test all(refined .& coarse .== coarse) # every choice allowed before refining is still allowed
            @test run(typemax(Int); online=true) == refined
        end
    end

    @testset "modal: an action outside its guard does nothing" begin
        # x' = 1, y' = 0 on [0, 10] × [0, 1]. x >= 9 is unsafe. Action 1 resets x to 0 and has the guard x >= 5.
        statespace = Hyperrectangle([5.0, 0.5], [5.0, 0.5])
        unsafe = [HPolyhedron([LazySets.HalfSpace([-1.0, 0.0], -9.0)])]
        U = Zonotope([1.0, 0.0], zeros(2, 2))
        reset = Action(1, HPolyhedron([LazySets.HalfSpace([-1.0, 0.0], -5.0)]), [0.0 0.0; 0.0 1.0], [0.0, 0.0])
        sys = ModalHybridSystem(statespace, unsafe, [ModalLocation(1, nothing, zeros(2, 2), U, Edge[], [ModalAction(reset, 1)])])
        grid, allowed, _, _ = ReACTed_reachable_cell_modal(sys, 1.0, 0.3, 0.01, 0.04; verbose=false)
        row = allowed[:, 1, 1]
        lower(i) = grid.lower[1] + (i - 1) * grid.granularity
        @test all(i -> row[i] == 0x03, findall(i -> lower(i) < 7.5, eachindex(row)))  # outside the guard acting is the same as waiting
        @test all(i -> row[i] == 0x02, findall(i -> 7.5 <= lower(i) < 8.7, eachindex(row)))  # waiting reaches x >= 9, the reset does not
        @test all(i -> row[i] == 0x00, findall(i -> lower(i) >= 8.7, eachindex(row)))
    end
end
