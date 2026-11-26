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
        @test count(==(0), allowed) == 18992
        @test digest(allowed) == "d43b7403bc19"
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

    @testset "modal: an axis narrower than a cell" begin
        # x' = y with y in [0, 0.9], one cell of g = 1. x >= 9 is unsafe, and an action resets x to 0. Over a period of 5 a
        # state moves up to 4.5, so waiting from x in [4, 5) can reach x = 9.49: unsafe. Flowed at a single y, the cell
        # missed that.
        statespace = Hyperrectangle([5.0, 0.45], [5.0, 0.45])
        unsafe = [HPolyhedron([LazySets.HalfSpace([-1.0, 0.0], -9.0)])]
        reset = Action(1, HPolyhedron([LazySets.HalfSpace([-1.0, 0.0], 0.0)]), [0.0 0.0; 0.0 1.0], [0.0, 0.0])
        sys = ModalHybridSystem(statespace, unsafe, [ModalLocation(1, nothing, [0.0 1.0; 0.0 0.0], Zonotope([0.0, 0.0], zeros(2, 2)), Edge[], [ModalAction(reset, 1)])])
        allowed = ReACTed_reachable_cell_modal(sys, 5.0, 1.0, 0.25, 0.5; act_needs_safe_flow=false, verbose=false)[2]
        @test vec(allowed[:, 1, 1]) == [0x03, 0x03, 0x03, 0x02, 0x02, 0x02, 0x02, 0x02, 0x00, 0x00]
    end

    @testset "modal: leaving through a face a guard only partly borders" begin
        # x' = -1 on [0, 10]². The guard x <= 0, y >= 5 wraps x to 10, so below y = 5 the states leave through x = 0.
        # y >= 9.5 is unsafe. In a period of 0.5 the first column crosses x = 0.
        guard = HPolyhedron([LazySets.HalfSpace([1.0, 0.0], 0.0), LazySets.HalfSpace([0.0, -1.0], -5.0)])
        sys = ModalHybridSystem(Hyperrectangle([5.0, 5.0], [5.0, 5.0]), [HPolyhedron([LazySets.HalfSpace([0.0, -1.0], -9.5)])],
            [ModalLocation(1, nothing, zeros(2, 2), Zonotope([-1.0, 0.0], zeros(2, 2)), [Edge(1, guard, [0.0 0.0; 0.0 1.0], [10.0, 0.0])], ModalAction[])])
        allowed = ReACTed_reachable_cell_modal(sys, 0.5, 1.0, 0.25, 0.5; verbose=false)[2]
        @test vec(allowed[1, :, 1]) == [0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x01, 0x01, 0x00, 0x00]
    end

    @testset "modal: clock axes at cell borders" begin
        # One axis x on [0, 8] with cells of 1 (cell c is [c-1, c), as a floor lookup sees it) and a period of 1. Each x has
        # a constant rate, so it is a clock: its successors are filtered by its exact interval. Successors are
        # (cell, layer); location 2 is layer 2.
        BE = ReACTShields.ReachabilityAnalysis.Exponentiation.BaseExp
        ge(b) = HPolyhedron([LazySets.HalfSpace([-1.0], -b)])
        rate(r) = Zonotope([r], zeros(1, 1))
        none = HPolyhedron{Float64,Vector{Float64}}[]
        space = Hyperrectangle([4.0], [4.0])
        two(r1, guard, M, j, r2) = ModalHybridSystem(space, none, [ModalLocation(1, nothing, zeros(1, 1), rate(r1), [Edge(2, guard, M, j)], ModalAction[]),
            ModalLocation(2, nothing, zeros(1, 1), rate(r2), Edge[], ModalAction[])])
        function succ(sys, c; act=0)
            ctx = ReACTShields.ModalContext(sys, 1.0, 1.0, 0.25, 0.5, BE, 5, 5, 1, 16; exit_bad_for_flow=false, exit_bad_for_actions=false)
            r = ReACTShields.box_successors(ctx, 1, CartesianIndices((c:c,)))
            return Set(Tuple(k) for k in (act == 0 ? r.flow : r.acts[act]))
        end
        # x' = 1, x >= 3 resets x to 0. All of [3, 4) is inside the guard at the start, so it jumps at once and ends at
        # exactly x = 1, in [1, 2). From [2, 3) the states jump later and end below 1.
        reset = two(1.0, ge(3.0), zeros(1, 1), [0.0], 1.0)
        @test (2, 2) in succ(reset, 4)
        @test succ(reset, 3) == Set([(1, 2), (4, 1)])
        # x' = 0, x >= 3 jumps to location 2. The border state 3 is not in [2, 3), so that cell never jumps.
        still = two(0.0, ge(3.0), ones(1, 1), [0.0], 0.0)
        @test succ(still, 3) == Set([(3, 1)])
        @test succ(still, 4) == Set([(4, 2)])
        # x' = 1, x >= 6 wraps x back by 6. A cell moves by exactly one cell, with no drift into its neighbours; the state
        # that reaches 6 at the end of the period jumps at the start of the next.
        wrap = ModalHybridSystem(space, none, [ModalLocation(1, nothing, zeros(1, 1), rate(1.0), [Edge(1, ge(6.0), ones(1, 1), [-6.0])], ModalAction[])])
        @test succ(wrap, 5) == Set([(6, 1)])
        @test succ(wrap, 6) == Set([(1, 1), (7, 1)])
        # x' = 1, an action flips x to 4 - x. [1, 2) becomes (2, 3] and ends in (3, 4]: the state from x = 1 ends at 4, in
        # [4, 5). A flip takes the closed end of the cell to the upper end, so the axis is not tracked after it.
        flip = ModalHybridSystem(space, none, [ModalLocation(1, nothing, zeros(1, 1), rate(1.0), Edge[], [ModalAction(Action(1, nothing, -ones(1, 1), [4.0]), 1)])])
        @test (4, 1) in succ(flip, 2; act=1)
        @test (5, 1) in succ(flip, 2; act=1)
    end

    @testset "modal simple: the ball as modal" begin
        # The ball has no axis that modal treats as a clock where it matters, so without them the shield is the same.
        sys = ModalHybridSystem(loadBouncingBallShieldedWithInput()[1])
        args = (sys, 0.1, 0.03, 0.001, 0.032)
        kw = (; exit_bad_for_flow=false, exit_bad_for_actions=false, verbose=false)
        @test ReACTed_reachable_cell_modal_simple(args...; kw...)[2] == ReACTed_reachable_cell_modal(args...; kw...)[2]
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
