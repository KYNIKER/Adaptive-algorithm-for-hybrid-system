# Regression tests for the shields. Run from the repository root with
#     julia --project=. test/runtests.jl
# The reference values pin what each shield computes, so a change in them means a shield now computes something
# different. If that is intended, update the values and say why in the commit.
using Test, SHA, LazySets
include(joinpath(@__DIR__, "..", "src", "ReACTShields.jl"))
using .ReACTShields

# First 12 hex digits of the SHA-1 of the allowed bits, so a single changed cell is caught.
digest(a) = bytes2hex(sha1(reinterpret(UInt8, vec(a))))[1:12]

@testset "ReACTShields" begin
    @testset "bouncing ball, granularity 0.02" begin
        sys, p = loadBouncingBallShieldedWithInput()
        g, δ⁻ = 0.02, 0.001
        for (name, run, dead, sha) in [
            ("plain", () -> ReACTed_reachable_cell_t(sys, p, g, δ⁻, 2^5 * δ⁻; verbose=false), 732669, "06cf1a99cab3"),
            ("pred_fast", () -> ReACTed_reachable_cell_online_pred_fast(sys, p, g, δ⁻, 2^5 * δ⁻; jump_hull=false, verbose=false), 231284, "d10988da9134"),
            ("adaptive_fast", () -> ReACTed_reachable_cell_adaptive_fast(sys, p, 2g, 1, δ⁻, 2^5 * δ⁻, patience=0, verbose=false), 753524, "b293994c2271"),
        ]
            @testset "$name" begin
                _, allowed, _, _ = run()
                @test size(allowed) == (1500, 750)
                @test count(==(0), allowed) == dead
                @test digest(allowed) == sha
                # Outside its guard an action does nothing, so it is allowed exactly where waiting is.
                for (ai, act) in enumerate(sys.Act)
                    outside = .!ReACTShields.guard_masks(Grid(sys.statespace, g), act.guard)[1]
                    @test all(i -> isodd(allowed[i]) == isodd(allowed[i] >> ai), findall(outside))
                end
            end
        end
    end

    @testset "modal: oil pump, granularity 0.2, δ 0.1 to 0.2" begin
        sys, p = loadOilPumpModal(; period=0.2)
        _, allowed, _, _ = ReACTed_reachable_cell_modal(sys, p, 0.2, 0.1, 0.2; online=true, verbose=false)
        @test size(allowed) == (100, 107, 10, 4)
        @test count(==(0), allowed) == 57515
        @test digest(allowed) == "920081310bff"
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
