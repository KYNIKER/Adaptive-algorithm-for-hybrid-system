using SparseArrays: sparsevec
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityAnalysis: add_dimension
using LazySets

include("../Utilities.jl")

function loadGirardExample(δ⁻::Float64 = 0.025)
    A1 = [-1 -4; 4 -1]
    A2 = [1 4; -4 -1]

    guard1 = HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [1.], 2), -0.5), # x == -0.5
        LazySets.HalfSpace(sparsevec([1], [-1.], 2), 0.5)
    ])
    guard2 = HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [1.], 2), -0.3), # x == -0.3
        LazySets.HalfSpace(sparsevec([1], [-1.], 2), 0.3)
    ])


    X0 = Zonotope([1., 0.], Matrix(0.1*I, 2, 2))
    μ = 0.001


    loc1 = Location(1, nothing, A1, I(2), getUFromInputUncertainty(A1, μ, δ⁻, X0), nothing, [Edge(2, guard1, I(2), zeros(2))], [])
    loc2 = Location(2, nothing, A2, I(2), getUFromInputUncertainty(A2, μ, δ⁻, X0), nothing, [Edge(1, guard2, I(2), zeros(2))], [])



    locations = [loc1, loc2]
    H = HybridSystemV2(locations, [])

    # These values are loosely based on Colas Le Guernic Reachability Analysis of Hybrid Systems with Lineara Continuous dynamics, 2009
    X0 = Zonotope([1., 0.], Matrix(0.1*I, 2, 2))
    T = 4.

    return H, 1, X0, T
end







# loadBouncingBall()