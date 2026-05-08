#example
using LazySets, ReachabilityAnalysis, LinearAlgebra, Polyhedra, Optim, Plots, CDDLib
using SparseArrays: sparsevec
include("Utilities.jl")


Z = Zonotope([0.0, 0.0], [1.0 2.0 -1.0; 0.0 1.0 0.4])
input = Zonotope([3.0, 0.0], [[-0.1, 0.5]])

#=
input = HPolyhedron([
    LazySets.HalfSpace(sparsevec([1], [1.], 2), 3.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([1], [-1.], 2), -3.0),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([2], [1.], 2), 0.8), # y <= 0.6
    LazySets.HalfSpace(sparsevec([2], [-1.], 2), -0.6), # y <= 0.6
    LazySets.HalfSpace(sparsevec([1, 2], [-1., 1.], 2), -2.84), # y <= 0.6
])
=#

approx = overapproximate(Z, BoxDirections(2))
inv = HPolyhedron([
    LazySets.HalfSpace(sparsevec([1], [1.], 2), 4.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([1], [-1.], 2), 1.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([2], [1.], 2), 2.6), # y <= 0.6
    LazySets.HalfSpace(sparsevec([2], [-1.], 2), 0.6), # y <= 0.6
    LazySets.HalfSpace(sparsevec([1, 2], [-1., 1.], 2), 0.7), # y <= 0.6
])
dirs = map(x -> x.a, constraints(inv))
bounds = map(x -> x.b, constraints(inv))

inv2 = HPolyhedron([
    LazySets.HalfSpace(sparsevec([1], [1.], 2), 4.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([1], [-1.], 2), 1.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([2], [1.], 2), 2.6), # y <= 0.6
    LazySets.HalfSpace(sparsevec([2], [-1.], 2), 0.6), # y <= 0.6
    LazySets.HalfSpace(sparsevec([1, 2], -1 * [-1., 1.5], 2), -0.7), # y <= 0.6
])
dirs2 = map(x -> x.a, constraints(inv2))
bounds2 = map(x -> x.b, constraints(inv2))


offsetDistance = map(x -> ρ(x, input), dirs)
offsetDistance2 = map(x -> ρ(x, input), dirs2)

fig = Plots.plot()

Plots.plot!(input, c=:purple, lab="Input", alpha=0.6)
Plots.plot!(Z + input, c=:cyan, lab="Zonotope with input", alpha=0.2)

Plots.plot!(approx + input, c=:blue, lab="Overapproximated Z", alpha=0.2)
Plots.plot!(inv, c=:yellow, lab="invariant", alpha=0.2)
Plots.plot!(inv2, c=:orange, lab="invariant2", alpha=0.2)

#flipped_constraint_list = map(x -> LazySets.HalfSpace(-x.a, -x.b), inv.constraints)
#Plots.plot!(HPolyhedron(flipped_constraint_list), c=:green, lab="bad compliment of invariant", alpha=0.2)

#intersection = constrain(offsetDistance, approx, Z, dirs, bounds) #zonotopeStripIntersection(Z, inv)
#intersection = constrain(vcat(offsetDistance, offsetDistance2), approx, Z, vcat(dirs, dirs2), vcat(bounds, bounds2)) #zonotopeStripIntersection(Z, inv)
#intersection = revise(input, approx, Z, vcat(dirs, dirs2), vcat(bounds, bounds2)) #zonotopeStripIntersection(Z, inv)
intersection = constrain(input, approx, Z, vcat(dirs, dirs2), vcat(bounds, bounds2)) #zonotopeStripIntersection(Z, inv)

#=println(typeof(intersection))
if intersects(intersection, inv)
    println(intersection)

    intersection = zonotopeStripIntersection(intersection, inv)
    println(typeof(intersection))

end

println(intersection)=#
Plots.plot!(intersection, c=:white, lab="Intersect", alpha=0.5)
println(ρ([-0.5, 1.0], Z + input))
intersectionr = revise(intersection, minkowski_sum(Z, input), [[-0.5, 1.0]], [ρ([-0.5, 1.0], Z)])
println(isempty(intersectionr))
Plots.plot!(intersectionr, c=:black, lab="Intersect2", alpha=0.5)

A::Matrix = [0.0 1.0; 0.0 0.0]

Plots.plot!(mapPolytope(A, mapPolytope(A, approx)), c=:brown, lab="A * Intersect2", alpha=0.5)
#Plots.plot!(rest, c=:blue, lab="Rest", alpha=0.5)
xlims!(fig, (-3.0, 10.0))
ylims!(fig, (-3.0, 3.0))
#println(rest)

display(fig)
