using LazySets, LinearAlgebra, Plots, ReachabilityAnalysis, JLD2, FileIO, BenchmarkTools, ProfileView
using Plots.PlotMeasures

gr()
include(joinpath(@__DIR__, "..", "src", "ReACTShields.jl"))
using .ReACTShields
include(joinpath(@__DIR__, "..", "src", "core", "discretize.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "utilities.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "partition.jl"))
include(joinpath(@__DIR__, "..", "src", "shields", "guards.jl"))
include(joinpath(@__DIR__, "..", "src", "shields", "propagate.jl"))
include(joinpath(@__DIR__, "..", "models", "bouncingBall.jl"))
BBsys, p = loadBouncingBallShieldedWithInput()
δ⁻ = 0.001
δ⁺ = 2^2 * δ⁻
alg = ReachabilityAnalysis.Exponentiation.BaseExp
phiDict, tPhiDict, inputDict = phi_input_dict(BBsys.flowMatrix, BBsys.input, δ⁻, δ⁺, alg, 5, 5)
d = δ⁻
#dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
isInvA = isinvertible(BBsys.flowMatrix)
A = BBsys.flowMatrix
#Φ = copy(phiDict[d])
A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(BBsys.flowMatrix)
Φcache = BBsys.flowMatrix == A_abs ? phiDict[d] : nothing
P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)
jump_reachTrue = JumpReach(BBsys, p, δ⁻, δ⁺, phiDict, tPhiDict, inputDict, P2A_abs, alg, 5, 5; jump_hull=true)
jump_reachFalse = JumpReach(BBsys, p, δ⁻, δ⁺, phiDict, tPhiDict, inputDict, P2A_abs, alg, 5, 5; jump_hull=false)

Z = Zonotope([-1.0, 0.1], [0.0 0.01; 0.01 0.0])
Z1 = Zonotope([1.0, 0.0], [1.0 0.0; 0.0 0.0])
Z2 = Zonotope([0.0, 0.0], [0.0 0.0; 1.0 0.0])

C1 = Singleton([-1.0, 0.1] + [0.0, 0.01] + [0.01, 0.0])
C2 = Singleton([-1.0, 0.1] + [0.0, 0.01] - [0.01, 0.0])
C3 = Singleton([-1.0, 0.1] - [0.0, 0.01] + [0.01, 0.0])
C4 = Singleton([-1.0, 0.1] - [0.0, 0.01] - [0.01, 0.0])
C5 = Singleton([-1.0, 0.1])

Cs = []#[C1, C2, C3, C4, C5]
for i in -1:0.125:1
    for j in -1:0.125:1
        push!(Cs, Singleton([-1.0, 0.1] + ([0.0, 0.01] .* i) + ([0.01, 0.0] .* j)))

    end
end


t_impact(corner) = (-corner.element[1] - sqrt(corner.element[1]^2 - 2*-9.81*corner.element[2]))/-9.81

zxT, _ = propagate_set_jumps(jump_reachTrue, Z, nothing)
zxF, _ = propagate_set_jumps(jump_reachFalse, Z, nothing)

A = [0.0 -1.0; 1. 0.0]
#d = 0.1
eA = exp(d * A)

pis = ReachabilityAnalysis.Exponentiation.Φ₁(A, d, ReachabilityAnalysis.Exponentiation.BaseExp, false, nothing)

plt = plot(dpi=1200, thickness_scaling=1, guidefontsize=25, minorgrid=true,
    legendfont=font(6, "Times"),
    legend_position=:outertopright,
    tickfont=font(8, "Times"),
    xguidefont=font(12, "Times"),
    yguidefont=font(12, "Times"),
    bottom_margin=2mm,
    left_margin=5mm,
    right_margin=5mm,
    top_margin=2mm,
    xlabel="x", ylabel="y")

xlims!(-1.5, 1.5)
ylims!(0, 0.2)

plot!(plt, Z, c=:green, lab="original zonotope")
#@show zx
for z in zxT
    #@show reduce_order(z, 2)
    #plot!(plt, reduce_order(z, 5), c=:red, la=0)
    plot!(plt, CH(reduce_order(z, 5), minkowski_sum(minkowski_sum(symmetric_interval_hull(linear_map(ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, -0.013, alg, isInvA, Φcache), symmetric_interval_hull(LazySets.linear_map(-A, BBsys.input)))), SymmetricIntervalHull(linear_map(ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, -0.013, alg, isInvA, Φcache), SymmetricIntervalHull(linear_map(-A * -A, z))))), minkowski_sum(linear_map(exp(-0.013 * A), z), linear_map(-0.013 * I, BBsys.input)))), c=:red, la=0)

end
for z in zxF
    #@show reduce_order(z, 2)
    plot!(plt, reduce_order(z, 3), c=:blue, la=0)
end

hybridJump = []

for c in Cs
    t = t_impact(c)
    new_v = -9.81 * t + c.element[1]
    new_p = 0.5 * -9.81 * t^2 + c.element[1]*t + c.element[2]
    s1 = [-9.81 * t + c.element[1], 0.5 * -9.81 * t^2 + c.element[1]*t + c.element[2]]
    j = [-0.85 0.0; 0.0 0.] * s1
    s2 = [-9.81 * (0.1-t) + j[1], 0.5 * -9.81 * (0.1-t)^2 + j[1]*(0.1-t) + j[2]]
    plot!(plt, c)
    plot!(plt, Singleton(s1))
    #plot!(plt, Singleton(s2))
    push!(hybridJump, Singleton(s2))
    #plot!(plt, linear_map(exp((0.1-t) * BBsys.flowMatrix), linear_map([-0.85 0.0; 0.0 0.], linear_map(exp((t) * BBsys.flowMatrix), c))))
end

hul = reduce(ConvexHull!, hybridJump)

plot!(plt, hul)
lens!([1, 1.5], [0.02, 0.05], inset=(1, bbox(0.5, 0.0, 0.4, 0.4)))
#plot!(plt, z, c=:blue, lab="original zonotope")

#=
plot!(plt, LazySets.API.translate(linear_map(transpose(-A * exp(d * -A)), Z2), transpose(-A * exp(d * -A)) * Z1.center), c=:blue, lab="original zonotope")

plot!(plt, linear_map(eA, Z), c=:orange, lab="linear map of original zonotope")
plot!(plt, linear_map(A * eA, Z), c=:red, lab="linear map of original zonotope")
plot!(plt, intersection(Z, linear_map(A * eA, Z)), c=:red, lab="linear map of original zonotope")
plot!(plt, linear_map(-A * exp(-d * A) * eA, Z), c=:red, lab="linear map of original zonotope")
plot!(plt, linear_map(A, Singleton([2.0, 0.0])), c=:yellow, lab="linear map of original zonotope")
#plot!(plt, minkowski_sum(linear_map(pis, Z1), Z), c=:blue, lab="linear map of original zonotope with integral")
=#
display(plt)