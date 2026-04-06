using LazySets, Polyhedra, GLMakie, LinearAlgebra, BenchmarkTools, CSV, DataFrames, Expokit, Plots#, GLMakie #, Plots, ReachabilityAnalysis

include("Utilities.jl")
include("ReACTed.jl")
include("models/gearbox.jl")
simple = true
const fig = Plots.plot()

if simple
    G = rand(-4.0:4.0, (2, 4))
    println(G)
    Z = Zonotope([0.0, 0.0], G)

    h₁ = LazySets.HalfSpace([0.1, 1.0], 0.0) #LazySets.HalfSpace(rand([-1.0, 1.0], 2), 0.0)
    h₂ = LazySets.HalfSpace([-0.1, -1.0], 0.0001)
    h₃ = LazySets.HalfSpace([1.0, 0.1], 0.0)
    H = HPolyhedron([h₁, h₂, h₃])#, LazySets.HalfSpace([0.0, -1.0, 0.0], 0.2)])
    println(isempty(∩(H, Z)))
    println(zonotopeStripIntersection(Z, H))
    Plots.plot!(zonotopeStripIntersection(Z, H), color=:red)
    Plots.plot!(box_approximation(zonotopeStripIntersection(Z, H)), color=:white)
    Plots.plot!(Z, color=:white)
    Plots.plot!(getBoxIntersection(Z, H), color=:blue)

    Plots.plot!(h₁, color=:yellow, alpha=0.2)
    Plots.plot!(h₂, color=:yellow, alpha=0.2)
    Plots.plot!(h₃, color=:green, alpha=0.2)
    savefig(fig, joinpath("results/", "plotIntersection2D.png"))
else
    G = rand(Float64, (3, 4))

    Z = Zonotope([0.0, 0.0, 0.0], G)

    H = HPolyhedron([LazySets.HalfSpace([1., 0., 0.], 0.0)])#, LazySets.HalfSpace([0.0, -1.0, 0.0], 0.2)])
    println(isempty(∩(H, Z)))
    #plot3d(Z, alpha=0.5)
    println(zonotopeStripIntersection(Z, H))
    LazySets.plot3d(zonotopeStripIntersection(Z, H), color=:red)
end