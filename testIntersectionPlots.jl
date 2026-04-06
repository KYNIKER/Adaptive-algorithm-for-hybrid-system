#example
using LazySets, ReachabilityAnalysis, LinearAlgebra, Polyhedra, Optim, Plots
using SparseArrays: sparsevec
include("Utilities.jl")


Z = Zonotope([2.0, 0.0], [1.0 2.0 -1.0; 0.0 1.0 0.4])

inv = HPolyhedron([
    LazySets.HalfSpace(sparsevec([1], [-1.], 2), 0.6),  # x <= 0.6
    #LazySets.HalfSpace(sparsevec([1], [1.], 2), 1.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([2], [1.], 2), 0.0), # y <= 0.6
    LazySets.HalfSpace(sparsevec([2], [-1.], 2), 0.0) # y <= 0.6
])



fig = Plots.plot()

Plots.plot!(Z, c=:purple, lab="Original Zonotope", alpha=0.2)
Plots.plot!(inv, c=:yellow, lab="invariant", alpha=0.2)
#flipped_constraint_list = map(x -> LazySets.HalfSpace(-x.a, -x.b), inv.constraints)
#Plots.plot!(HPolyhedron(flipped_constraint_list), c=:green, lab="bad compliment of invariant", alpha=0.2)

intersection = zonotopeStripIntersection(Z, inv)
println(typeof(intersection))
if intersects(intersection, inv)
    println(intersection)

    intersection = zonotopeStripIntersection(intersection, inv)
    println(typeof(intersection))

end

println(intersection)
Plots.plot!(intersection, c=:black, lab="Intersect", alpha=0.5)
#Plots.plot!(rest, c=:blue, lab="Rest", alpha=0.5)

#println(rest)

display(fig)
