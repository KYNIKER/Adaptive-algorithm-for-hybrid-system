using SparseArrays: sparsevec
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityAnalysis: add_dimension
using LazySets

include("../../Utilities.jl")

function loadembrake(x0, Tsample, ζ)
    # model's constants
    L = 1.e-3
    KP = 10000.0
    KI = 1000.0
    R = 0.5
    K = 0.02
    drot = 0.1
    i = 113.1167


    # state variables: [I, x, xe, xc]
    A = Matrix([-(R+K^2/drot)/L 0 KP/L KI/L;
                K/i/drot        0 0    0;
                0               0 0    0;
                0               0 0    0])
    # reset map
    Ar = sparse([1, 2, 3, 4, 4], [1, 2, 2, 2, 4], [1.0, 1.0, -1.0, -Tsample, 1.0], 4, 4)
    br = sparsevec([3, 4], [x0, Tsample * x0], 4)

    # initial condition
    X0 = Zonotope([0.0, 0., 0., 0.], [[0.0, 0.0, 0.0, 0.0]])

    # invariant
    #invariant = HPolyhedron([
    #    LazySets.HalfSpace(sparsevec([5], [1.], 5), Tsample+ζ) # T <= Tsample + ζ
    #])
  
    #edge 
    #guard = HPolyhedron([LazySets.HalfSpace(SingleEntryVector(5, 5, -1.), Tsample-ζ)]) # T >= Tsample - ζ
    #guard = HPolyhedron([LazySets.HalfSpace(sparsevec([5], [-1.], ), -(Tsample-ζ))]) # T >= Tsample - ζ
    edges::Vector{Edge} = [Edge(1, nothing, Ar, br)]

    # location
    locations = [Location(1, nothing, A, nothing, nothing, nothing, edges, [])]

    # global constraint 
    property = LazySets.HalfSpace(sparsevec([2], [1.], 4), x0)

    H = HybridSystemV2(locations, [property])

    T = 0.02
    return H, 1, X0, T
end