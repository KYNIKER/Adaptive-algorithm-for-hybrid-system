#using SparseArrays: sparsevec
#using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
#using ReachabilityAnalysis: add_dimension
using LazySets

include("../Utilities.jl")

# ego can only accelerate by +/- 2m/s^2.
# I think this can be modelled by either euclidean hybrid system or the modal system. 
function loadCruiseControl()
    ego_sensor_range = 200
    v_ego_min = -10.0
    v_ego_max = 20.0
    v_front_min = -8.0
    v_front_max = 20.0

    dims = 3
    distance = 1
    v_ego = 2
    v_front = 3
    S = Hyperrectangle([101.0, 10.0, 10.0], [(ego_sensor_range-1)/2, 10.0, 10.0])
    A = zeros(dims, dims)
    A[distance, v_ego] = -1.0
    A[distance, v_front] = 1.0


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

    backwards = Action(1, HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [-1.], 3), -4.0), # p >= 4
        LazySets.HalfSpace(sparsevec([2], [1.], 3), 0.0),  # v <= 0
        LazySets.HalfSpace(sparsevec([2], [-1.], 3), 4.0)  # v > -4
    ]), diagm([1.0, 0.0, 1.0]), [0.0, -4.0, 0.0])

    forwards = Action(2, HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [-1.], 3), -4.0), # p >= 4
        LazySets.HalfSpace(sparsevec([2], [-1.], 3), 0.0)  # v >= 0
    ]), diagm([1.0, 0.9, 1.0]), [0.0, -4.0, 0.0])

    Act = [backwards, forwards]

    H = EuclideanHybridSystem(S, Zonotope(zeros(3), zeros(3, 3)), [constraint], A, edges, Act)

    T = 1.0

    return H, T
end