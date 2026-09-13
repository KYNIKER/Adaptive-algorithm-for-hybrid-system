#example
using LazySets, ReachabilityAnalysis, LinearAlgebra, Polyhedra, Optim, Plots
using SparseArrays: sparsevec
include("Utilities.jl")

function allunder(set, offset, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(offset, dir, bound)
        res = res && (offset[i] + ρ(dir[i], set) <= bound[i])
    end
    return res
end
function someunder(set, offset, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(offset, dir, bound)
        res = res && (offset[i] - ρ(-dir[i], set) <= bound[i])
    end
    return res
end

function someoutside(set, offset, dir, bound)
    res = false
    @inbounds @simd for i in eachindex(offset, dir, bound)
        res = res || (offset[i] - ρ(-dir[i], set) > bound[i])
    end
    return res
end
Z = Zonotope([2.0, 0.0], [1.0 2.0 -1.0; 0.0 1.0 0.4])

inv = HPolyhedron([
    LazySets.HalfSpace(sparsevec([1], [-1.], 2), 1.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([1], [1.], 2), -0.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([2], [1.], 2), 0.5), # y <= 0.6
    LazySets.HalfSpace(sparsevec([2], [-1.], 2), -0.2) # y <= 0.6
])

hs = [
    LazySets.HalfSpace(sparsevec([1], [-1.], 2), 1.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([1], [1.], 2), -0.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([2], [1.], 2), 0.5), # y <= 0.6
    LazySets.HalfSpace(sparsevec([2], [-1.], 2), -0.2) # y <= 0.6
]

IZ = Zonotope([-1.0, 0.75], [1.0 0.0; 0.0 0.25])
BZ = box_approximation(IZ)

invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(inv)

fig = Plots.plot()

Plots.plot!(Z, c=:purple, lab="Original Zonotope", alpha=0.2)
Plots.plot!(IZ, c=:blue, lab="inviariant Zonotope", alpha=0.2)

Plots.plot!(inv, c=:yellow, lab="invariant", alpha=0.6)
#flipped_constraint_list = map(x -> LazySets.HalfSpace(-x.a, -x.b), inv.constraints)
#Plots.plot!(HPolyhedron(flipped_constraint_list), c=:green, lab="bad compliment of invariant", alpha=0.2)
@time intersects(Z, inv)
@time intersects(Z, hs)
@time someoutside(Z, [0.0, 0.0, 0.0, 0.0], invarientProjVectors, invarientProjBounds)
@time intersects(Z, IZ)
@show intersects(Z, IZ)
@time intersects(Z, BZ)
@show intersects(Z, BZ)

#@show allunder(Z, [0.0, 0.0], invarientProjVectors, invarientProjBounds) #[0.0, 0.0, 0.0, 0.0]
#@show someoutside(Z, [0.0, 0.0], invarientProjVectors, invarientProjBounds)
#@show someunder(Z, [0.0, 0.0], invarientProjVectors, invarientProjBounds)

@show allunder(Z, [0.0, 0.0, 0.0, 0.0], invarientProjVectors, invarientProjBounds) #[0.0, 0.0, 0.0, 0.0]
@show someoutside(Z, [0.0, 0.0, 0.0, 0.0], invarientProjVectors, invarientProjBounds)
@show someunder(Z, [0.0, 0.0, 0.0, 0.0], invarientProjVectors, invarientProjBounds)
intersection = zonotopeStripIntersection(Z, inv)
println(intersection)
#=
if intersects(intersection, inv)
    println(intersection)

    intersection = zonotopeStripIntersection(intersection, inv)
    println(typeof(intersection))

end

println(intersection)
=#
Plots.plot!(intersection, c=:black, lab="Intersect", alpha=0.5)
#Plots.plot!(rest, c=:blue, lab="Rest", alpha=0.5)

#println(rest)

@show ρ([1.0, 0.0], Z)
@show ρ([0.0, 1.0], Z)
@show ρ([1.0, 1.0], Z)


display(fig)
