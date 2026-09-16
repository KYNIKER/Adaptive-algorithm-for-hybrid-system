using SparseArrays: sparsevec
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityAnalysis: add_dimension
using LazySets

include("../Utilities.jl")

function loadBouncingBall()
    A = [0.0 1.0 0.0; 0.0 0.0 -(9.82/ℯ); 0. 0. 0.]
    c = [0.0, -9.82, 0.0]
    #u = Zonotope(c, zeros(3, 1))
    # u = overapproximate(u, Zonotope)
    u = Zonotope(zeros(3), zeros(3, 0))

    guard = HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [1.], 3), 0.0), # x <= 0
        LazySets.HalfSpace(sparsevec([1], [-1.], 3), 0.0),  # x >= 0
        LazySets.HalfSpace(sparsevec([2], [1.], 3), 0.0)  # y <= 0
    ])

    jumpMatrix = [0.0 0.0 0.; 0.0 -0.85 0.; 0.0 0.0 1.0]

    edges::Vector{Edge} = [Edge(1, guard, jumpMatrix, zeros(3))]

    locations = [Location(1, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [-1.], 3), -0.0)]), A, I(3), u, nothing, edges, [])]

    H = HybridSystemV2(locations, [])

    # These values are loosely based on Colas Le Guernic Reachability Analysis of Hybrid Systems with Lineara Continuous dynamics, 2009
    #X0 = Zonotope([1., 0., -1.], [[0.0, 0.0, 0.0]])
    X0 = Zonotope([1., 0., 1.], zeros(3, 1))

    T = 10

    return H, 1, X0, T
end

function loadBouncingBallShielded()
    A = [0.0 1.0 0.0; 0.0 0.0 -(9.82/ℯ); 0. 0. 0.]
    S = Hyperrectangle([7.5, 0.0, 1.0], [8.5, 15.0, 0])

    guard = HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [1.], 3), 0.0), # x <= 0
        #LazySets.HalfSpace(sparsevec([1], [-1.], 3), 0.0),  # x >= 0
        LazySets.HalfSpace(sparsevec([2], [1.], 3), 0.0)  # y <= 0
    ])

    jumpMatrix = [0.0 0.0 0.; 0.0 -0.85 0.; 0.0 0.0 1.0]

    edges::Vector{Edge} = [Edge(1, guard, jumpMatrix, zeros(3))]

    constraint = HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [1.], 3), 0.0), # x <= 0
        #LazySets.HalfSpace(sparsevec([1], [-1.], 3), 0.0),  # x >= 0
        LazySets.HalfSpace(sparsevec([2], [1.], 3), 0.0),  # y <= 0
        LazySets.HalfSpace(sparsevec([2], [-1.], 3), 1e-5)  # y <= 0
    ])

    #constraint = LazySets.API.intersection(constraint, S)

    act1 = Action(1, HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [-1.], 3), -4.0), # x >= 4
        LazySets.HalfSpace(sparsevec([2], [1.], 3), 0.0),  # y <= 0
        LazySets.HalfSpace(sparsevec([2], [-1.], 3), 4.0)  # y > -4
    ]), diagm([1.0, 0.0, 1.0]), [0.0, -4.0, 0.0])

    act2 = Action(2, HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [-1.], 3), -4.0), # x >= 4
        LazySets.HalfSpace(sparsevec([2], [-1.], 3), 0.0)  # y >= 0
    ]), diagm([1.0, 0.9, 1.0]), [0.0, -4.0, 0.0])

    Act = [act1, act2]

    H = EuclideanHybridSystem(S, [constraint], A, edges, Act)

    T = 0.1

    return H, T
end

function loadBouncingBallNoEdges()
    A = [0.0 -1.0; 0.0 9.81]
    c = [0.0, -9.81]
    #u = Zonotope(c, [zero(c)])
    u = Zonotope(zeros(2), [zero(c)])

    guard = HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [1.], 2), 0.0), # x <= 0
        LazySets.HalfSpace(sparsevec([1], [-1.], 2), 0.0),  # x >= 0
        LazySets.HalfSpace(sparsevec([2], [1.], 2), 0.0)  # y <= 0
    ])

    jumpMatrix = [1.0 0.0; 0.0 -0.89]

    edges::Vector{Edge} = [Edge(1, guard, jumpMatrix, zeros(2))]
    #push!(edges,)
    #locations = [Location(1, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [-1.], 3), -0.0)]), A, I(3), u, nothing, edges, [])]
    locations = [Location(1, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [-1.], 2), -0.0)]), A, I(2), u, nothing, edges, [])]

    H = HybridSystemV2(locations, [])

    # These values are loosely based on Colas Le Guernic Reachability Analysis of Hybrid Systems with Lineara Continuous dynamics, 2009
    X0 = Zonotope([10., 1.], [[0.001, 0.0]])
    T = 10

    return H, 1, X0, T
end





# loadBouncingBall()