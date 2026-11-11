#using SparseArrays: sparsevec
#using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
#using ReachabilityAnalysis: add_dimension
using LazySets

include("../Utilities.jl")

function loadBouncingBallShielded()

    A = [0.0 1.0 0.0; 0.0 0.0 -9.82; 0. 0. 0.]
    S = Hyperrectangle([8.5, 0.0, 1.0], [8.5, 15.0, 0]) # p: [-1, 16], v: [-15, 15], g: [1,1]

    guard = HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [1.], 3), 0.0), # p <= 0
        LazySets.HalfSpace(sparsevec([1], [-1.], 3), 0.0),  # p >= 0
        LazySets.HalfSpace(sparsevec([2], [1.], 3), 0.0)  # v <= 0
    ])

    jumpMatrix = [0.0 0.0 0.; 0.0 -0.85 0.; 0.0 0.0 1.0]

    edges::Vector{Edge} = [Edge(1, guard, jumpMatrix, zeros(3))]

    constraint = HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [1.], 3), 0.0), # p <= 0
        #LazySets.HalfSpace(sparsevec([1], [-1.], 3), 0.0),  # p >= 0
        LazySets.HalfSpace(sparsevec([2], [1.], 3), 0.0),  # v <= 0
        LazySets.HalfSpace(sparsevec([2], [-1.], 3), 1e-5)  # v >= 1e-5
    ])

    #constraint = LazySets.API.intersection(constraint, S)

    act1 = Action(1, HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [-1.], 3), -4.0), # p >= 4
        LazySets.HalfSpace(sparsevec([2], [1.], 3), 0.0),  # v <= 0
        LazySets.HalfSpace(sparsevec([2], [-1.], 3), 4.0)  # v > -4
    ]), diagm([1.0, 0.0, 1.0]), [0.0, -4.0, 0.0])

    act2 = Action(2, HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [-1.], 3), -4.0), # p >= 4
        LazySets.HalfSpace(sparsevec([2], [-1.], 3), 0.0)  # v >= 0
    ]), diagm([1.0, 0.9, 1.0]), [0.0, -4.0, 0.0])

    Act = [act1, act2]

    H = EuclideanHybridSystem(S, Zonotope(zeros(3), zeros(3, 3)), [constraint], A, edges, Act)

    T = 0.1

    return H, T
end