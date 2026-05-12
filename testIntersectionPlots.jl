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
    LazySets.HalfSpace(sparsevec([1], [-1.], 2), 1.2),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([2], [1.], 2), 2.6), # y <= 0.6
    LazySets.HalfSpace(sparsevec([2], [-1.], 2), 0.6), # y <= 0.6
    LazySets.HalfSpace(sparsevec([1, 2], [-1., 1.], 2), 0.7), # y <= 0.6
])
test = HPolytope([
    LazySets.HalfSpace(sparsevec([1], [1.], 2), 4.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([1], [-1.], 2), 1.2),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([2], [1.], 2), 2.6), # y <= 0.6
    LazySets.HalfSpace(sparsevec([2], [-1.], 2), 0.6), # y <= 0.6
    LazySets.HalfSpace(sparsevec([1, 2], [-1., 1.], 2), 0.7), # y <= 0.6
])

test3d = HPolytope([LazySets.HalfSpace(x, 1.0) for x in collect(BoxDirections(3))])
test10d = HPolytope([LazySets.HalfSpace(x, 1.0) for x in collect(BoxDirections(10))])
testHalfspace = LazySets.HalfSpace([1.0, 1.0, 0.0], 2.0)

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

vs3d, ps3d = vertexRep(test3d)
vs10d, ps10d = vertexRep(test10d)
#FromVertices(vs3d[1], ps3d[1])
#@show vs3d[1], htest
#@show vs3d
#@show ps3d


Plots.plot!(input, c=:purple, lab="Input", alpha=0.6)
Plots.plot!(Z + input, c=:cyan, lab="Zonotope with input", alpha=0.2)
Plots.plot!(approx, c=:brown, lab="approx", alpha=0.2)

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

vs, ps = vertexRep(test)
@show size(ps)
@show ps[1]
v1 = [x[1] for x in vs]
v2 = [x[2] for x in vs]
p1 = [y[1][1] for y in ps]
p2 = [y[1][2] for y in ps]

@show size(vs)

Plots.scatter!(v1, v2, c=:orange, lab="", alpha=0.5)
Plots.scatter!(p1, p2, c=:red, lab="", alpha=0.5)

A::Matrix = [1.0 0.0; 0.0 0.0]

#Plots.plot!(mapPolytope(A, approx), c=:pink, lab="A * Intersect2", alpha=0.5)
#Plots.plot!(rest, c=:blue, lab="Rest", alpha=0.5)
intPoint = sample(approx)
vsA, psA = vertexRep(approx)
@show constraints_list(A * approx)
@show vsA, psA
map!(x -> A * x, vsA)
map!(x -> map(y -> A * y, x), psA) #map(y -> A * y, x)
mIntPoint = A * intPoint
#@show psA
listHspaces::Vector{LazySets.HalfSpace} = []
for (x, y) in zip(vsA, psA)
    global listHspaces = vcat(listHspaces, halfspaceFromVertices(x, y, mIntPoint))
end
@show listHspaces
newmethod = HPolytope(listHspaces)
Plots.plot!(mapPolytope(A, approx), c=:pink, lab="A * Intersect2", alpha=1.0, lw=1.25)
Plots.plot!(newmethod, c=:black, lab="new method", alpha=0.8, lw=0.5, ls=:dash)


xlims!(fig, (-5.0, 10.0))
ylims!(fig, (-3.0, 3.0))
#println(rest)

display(fig)
