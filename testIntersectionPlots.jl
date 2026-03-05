#example
using LazySets, ReachabilityAnalysis, LinearAlgebra, Polyhedra, Optim, Plots
using SparseArrays: sparsevec
include("Utilities.jl")


Z = Zonotope([2.0, 0.0], [1.0 2.0; 0.0 1.0])

inv = HPolyhedron([
    LazySets.HalfSpace(sparsevec([1], [1.], 2), 0.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([2], [1.], 2), 0.6), # y <= 0.6
])



fig = Plots.plot()

Plots.plot!(Z, c=:purple, lab="Original Zonotope", alpha=0.2)
Plots.plot!(inv, c=:yellow, lab="invariant", alpha=0.2)
flipped_constraint_list = map(x -> LazySets.HalfSpace(-x.a, -x.b), inv.constraints)
Plots.plot!(HPolyhedron(flipped_constraint_list), c=:green, lab="bad compliment of invariant", alpha=0.2)

intersect, rest = splitZonotope(Z, inv)

Plots.plot!(intersect, c=:black, lab="Intersect", alpha=0.5)
Plots.plot!(rest, c=:blue, lab="Rest", alpha=0.5)

println(intersect)
println(rest)

display(fig)
