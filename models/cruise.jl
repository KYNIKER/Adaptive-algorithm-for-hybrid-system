#using SparseArrays: sparsevec
#using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
#using ReachabilityAnalysis: add_dimension
using LazySets

include("../Utilities.jl")


const EGO_NEUTRAL_FRONT_NEUTRAL, EGO_NEUTRAL_FRONT_ACCELERATE, EGO_NEUTRAL_FRONT_DECELERATE,
EGO_ACCELERATE_FRONT_NEUTRAL, EGO_ACCELERATE_FRONT_ACCELERATE, EGO_ACCELERATE_FRONT_DECELERATE,
EGO_DECELERATE_FRONT_NEUTRAL, EGO_DECELERATE_FRONT_ACCELERATE, EGO_DECELERATE_FRONT_DECELERATE = 1:9

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

    e(i) = [j == i ? 1.0 : 0.0 for j in 1:dims]

    S = Hyperrectangle([8.5, 0.0, 1.0], [8.5, 15.0, 0]) # p: [-1, 16], v: [-15, 15], g: [1,1]

    maxspeedlimit_ego = HPolyhedron([
        LazySets.HalfSpace(sparsevec([v_ego], [-1.], 3), -v_ego_max) # v_ego <= 0
    ])

    minspeedlimit_ego = HPolyhedron([
        LazySets.HalfSpace(sparsevec([v_ego], [1.], 3), -v_ego_min) # v_ego >= 0
    ])

    maxspeedlimit_front = HPolyhedron([
        LazySets.HalfSpace(sparsevec([v_front], [-1.], 3), -v_front_max) # v_front <= 0
    ])

    minspeedlimit_front = HPolyhedron([
        LazySets.HalfSpace(sparsevec([v_front], [1.], 3), -v_front_min) # v_front >= 0
    ])

    speedlimit_ego = intersection(minspeedlimit_ego, maxspeedlimit_ego)

    layers = [(EGO_NEUTRAL_FRONT_NEUTRAL, 0, 0),
        (EGO_NEUTRAL_FRONT_ACCELERATE, 0, 1),
        (EGO_NEUTRAL_FRONT_DECELERATE, 0, -1),
        (EGO_ACCELERATE_FRONT_NEUTRAL, 1, 0),
        (EGO_ACCELERATE_FRONT_ACCELERATE, 1, 1),
        (EGO_ACCELERATE_FRONT_DECELERATE, 1, -1),
        (EGO_DECELERATE_FRONT_NEUTRAL, -1, 0),
        (EGO_DECELERATE_FRONT_ACCELERATE, -1, 1),
        (EGO_DECELERATE_FRONT_DECELERATE, -1, -1)]

    lookupDict = Dict()

    for (l, eg, fr) in layers
        lookupDict[(eg, fr)] = l
    end

    constraint = HPolyhedron([
        LazySets.HalfSpace(sparsevec([distance], [1.], dims), 5.0) # d <= 0
    ])

    backwards = Action(1, HPolyhedron([
        LazySets.HalfSpace(sparsevec([v_ego], [-1.], 3), -8.0) # v_ego >= -8
    ]), diagm([1.0, 0.0, 1.0]), [0.0, -4.0, 0.0])

    forwards = Action(2, HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [-1.], 3), -4.0), # p >= 4
        LazySets.HalfSpace(sparsevec([2], [-1.], 3), 0.0)  # v >= 0
    ]), diagm([1.0, 0.9, 1.0]), [0.0, -4.0, 0.0])

    zero_input = Zonotope(zeros(dims), zeros(dims, dims))
    locations = ModalLocation[]
    dia = diagm(ones(dims))
    for (i, eg, fr) in layers
        Al = copy(A)
        Al[v_ego] = eg * 2.0
        Al[v_front] = fr * 2.0

        edges = []

        eg > 0 && push!(edges, Edge(lookupDict(0, fr), maxspeedlimit_ego, dia, zeros(dims)))
        eg < 0 && push!(edges, Edge(lookupDict(0, fr), minspeedlimit_ego, dia, zeros(dims)))
        fr > 0 && push!(edges, Edge(lookupDict(eg, 0), maxspeedlimit_front, dia, zeros(dims)))
        fr < 0 && push!(edges, Edge(lookupDict(eg, 0), minspeedlimit_front, dia, zeros(dims)))
        actions = modalAction[]
        for j in -1:1:1
            j == eg && continue
            push!(actions, Action(lookupDict(j, fr), speedlimit_ego), dia, zeros(dims))
        end
    end
    #edges::Vector{Edge} = [Edge(1, guard, jumpMatrix, zeros(3))]
    #constraint = LazySets.API.intersection(constraint, S)


    Act = [backwards, forwards]

    H = EuclideanHybridSystem(S, Zonotope(zeros(3), zeros(3, 3)), [constraint], A, edges, Act)

    T = 1.0

    return H, T
end