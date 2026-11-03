using SparseArrays: sparsevec
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityAnalysis: add_dimension
using LazySets

# Automatic Synthesis of Robust and Optimal Controllers – An Industrial Case Study
# It looks like this should be modelled as switched dynamics...
# Looks like CoSHY was able to model this with a 4 dimensional state. Thus they only track time, volume of oil in the tank, state of the pump and a latency timer.
# They compute the volume of oil in the tank simply by the flow rate times the time spent in that rate. The act of turning the pump on is simply modelled subtracting the outgoing flow until next time step.
function loadOilPump()
    V_max = 25.1
    V_min = 4.9
    dims = 8
    machine_clock = 1
    machine_consumption_rate = 2
    pump_clock = 3
    pump_rate = 4
    accumulator_volume = 5
    accumulator_total_volume = 6
    inflow = 2.2
    fluctuation = 0.1
    imprecision = 0.06
    latency = 2.00

    A_1 = zeros(Float64, (dims, dims))
    A_2 = zeros(Float64, (dims, dims))

    consumption_rate_over_time = zip([0.0, 1.2, 0.0, 1.2, 2.5, 0.0, 1.7, 0.5, 0.0], collect(0.0:2.0:20.0))

    input_noise = Zonotope(zeros(Float64, dims), [sparsevec([machine_consumption_rate], [fluctuation], dims)]) # epsilon 

    constraint = HPolyhedron([
        LazySets.HalfSpace(sparsevec([accumulator_volume], [1.], dims), V_max), # v <= 25.1
        LazySets.HalfSpace(sparsevec([accumulator_volume], [-1.], dims), -V_min) # v <= 0
    ])

    A = [0.0 0.0; 1.0 0.0]
    c = [-9.81, 0.0]
    u = Zonotope(c, [zero(c)])
    #u = Zonotope(zeros(2), [zero(c)])
    statespace_radius = ntuple(i != accumulator_volume ? 10.0 : 0.0, dims)
    S = Hyperrectangle(sparsevec([accumulator_volume], [(V_max + V_min)/2], dims), statespace_radius .+ sparsevec([accumulator_volume], [(V_max - V_min)/2], dims))
    #S = Hyperrectangle([0.0, 3.5], [13., 3.5]) # p: [-1, 16], v: [-15, 15], g: [1,1]. Should be [0.0, 5.0], [14., 5.0]
    #S = Hyperrectangle([0.0, 0.25], [3., 0.25]) # TEST 

    guard = HPolyhedron([
        LazySets.HalfSpace(sparsevec([2], [1.], 2), 0.0), # p <= 0
        #LazySets.HalfSpace(sparsevec([1], [-1.], 2), 0.0),  # p >= 0
        LazySets.HalfSpace(sparsevec([1], [1.], 2), 0.0)  # v <= 0
    ])

    jumpMatrix = [-0.85 0.0; 0.0 0.]
    #jumpMatrix = [-1. 0.0; 0.0 0.0]

    edges::Vector{Edge} = [Edge(1, guard, jumpMatrix, zeros(2))]

    constraint = HPolyhedron([
        LazySets.HalfSpace(sparsevec([2], [1.], 2), 0.0), # p <= 0
        #LazySets.HalfSpace(sparsevec([1], [-1.], 2), 0.0),  # p >= 0
        LazySets.HalfSpace(sparsevec([1], [1.], 2), 1.0),  # v <= 0
        LazySets.HalfSpace(sparsevec([1], [-1.], 2), 1.0)  # v >= 1e-5
    ])

    #constraint = LazySets.API.intersection(constraint, S)

    act1 = Action(1, HPolyhedron([
        LazySets.HalfSpace(sparsevec([2], [-1.], 2), -4.0), # p >= 4
        LazySets.HalfSpace(sparsevec([1], [1.], 2), 0.0),  # v <= 0
        LazySets.HalfSpace(sparsevec([1], [-1.], 2), 4.0)  # v > -4
    ]), diagm([0.0, 1.0]), [-4.0, 0.0])

    act2 = Action(2, HPolyhedron([
        LazySets.HalfSpace(sparsevec([2], [-1.], 2), -4.0), # p >= 4
        LazySets.HalfSpace(sparsevec([1], [-1.], 2), 0.0)  # v >= 0
    ]), diagm([-0.9, 1.]), [-4.0, 0.0])

    Act = [act1, act2]

    H = EuclideanHybridSystem(S, u, [constraint], A, edges, Act)

    T = 20.0

    return H, T
end