using SparseArrays: sparsevec
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityAnalysis: add_dimension
using LazySets

include("../Utilities.jl")

function loadBouncingBall()
    A = [0.0 1.0 0.0; 0.0 0.0 1.0; 0. 0. 0.]
    c = [0.0, -9.82, 0.0]
    u = Zonotope(zero(c), [zero(c)])

    guard = HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [1.], 3), 0.0), # x <= 0
        LazySets.HalfSpace(sparsevec([1], [-1.], 3), 0.0),  # x >= 0
        LazySets.HalfSpace(sparsevec([2], [1.], 3), 0.0)  # y <= 0
    ])

    jumpMatrix = [0.0 0.0 0.; 0.0 -0.75 0.; 0.0 0.0 1.0]

    edges::Vector{Edge} = [Edge(1, guard, jumpMatrix, zeros(3))]

    locations = [Location(1, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [-1.], 3), -0.0)]), A, I(3), u, nothing, edges, [])]

    H = HybridSystemV2(locations, [])

    # These values are loosely based on Colas Le Guernic Reachability Analysis of Hybrid Systems with Lineara Continuous dynamics, 2009
    X0 = Zonotope([1., 0., -1.], [[0.0, 0.0, 0.0]])
    T = 10

    return H, 1, X0, T
end

# loadBouncingBall()