using SparseArrays: sparsevec
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityAnalysis: add_dimension
using LazySets

# Automatic Synthesis of Robust and Optimal Controllers – An Industrial Case Study
# It looks like this should be modelled as switched dynamics...
# Looks like CoSHY was able to model this with a 4 dimensional state. Thus they only track time, volume of oil in the tank, state of the pump and a latency timer.
# They compute the volume of oil in the tank simply by the flow rate times the time spent in that rate. The act of turning the pump on is simply modelled subtracting the outgoing flow until next time step.
# I think this needs to be modelled as act -> translates set in pump_rate dimension. Then the reachability relation should be computed for each decision period in the time interval.
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

    A = zeros(Float64, (dims, dims))
    A[accumulator_total_volume, pump_rate] = log(2.2)
    A[accumulator_volume, pump_rate] = log(2.2)
    A[accumulator_volume, machine_consumption_rate] = -log(2.2)
    A[machine_clock, machine_clock] = 0.0

    A_1 = zeros(Float64, (dims, dims))
    A_2 = zeros(Float64, (dims, dims))

    consumption_rate_over_time = zip([0.0, 1.2, 0.0, 0.0, 1.2, 2.5, 0.0, 1.7, 0.5, 0.0], collect(0.0:2.0:20.0))

    input_noise = Zonotope(sparsevec([machine_clock], [1.0], dims), [sparsevec([machine_consumption_rate], [fluctuation], dims)]) # epsilon 

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

# The oil pump as a ModalHybridSystem (ReACTedShieldingModal.jl). Returns the system and the decision period.
# Continuous state: (t, v, l) = time in the consumption cycle, oil volume, time left of the pump latency.
# Layers: pump off/on x latency free/locked. Each layer has one location per piece of the consumption profile, with
# guards on t between them. The latency runs in the locked layers and ends with an edge to the free layer, the only place
# the pump can switch.
# For ReACTed_reachable_cell_continuous, where the pump can switch at any time, use latency_reset=2.0 and latency_top=2.05: the
# value just below the latency relies on decisions a fixed period apart, and l = 2 needs a cell of its own.
const OILPUMP_OFF_FREE, OILPUMP_OFF_LOCKED, OILPUMP_ON_FREE, OILPUMP_ON_LOCKED = 1, 2, 3, 4

function loadOilPumpModal(; period=0.2, latency_reset=2.0 - 1e-3, latency_top=2.0)
    time_axis, volume_axis, latency_axis = 1, 2, 3
    V_min, V_max = 4.9, 25.1
    inflow = 2.2
    fluctuation = 0.1 # l/s, only when the consumption rate is not 0
    margin = max(0.6, 3 * period) # more than the largest change in volume over a period (2.6 l/s), so no reach set leaves the grid
    # Consumption pieces (start, end, rate)
    pieces = [(0.0, 2.0, 0.0), (2.0, 4.0, 1.2), (4.0, 8.0, 0.0), (8.0, 10.0, 1.2), (10.0, 12.0, 2.5), (12.0, 14.0, 0.0), (14.0, 16.0, 1.7), (16.0, 18.0, 0.5), (18.0, 20.0, 0.0)]
    cycle = 20.0
    # By default switching sets l just below the latency of 2 s. Decisions are period apart, so the pump is free at exactly
    # the same decision, and the value stays strictly inside a cell of the grid instead of on the border between two.
    tol = 1e-6

    e(i) = [j == i ? 1.0 : 0.0 for j in 1:3]
    A = zeros(3, 3)
    keep_all = Matrix{Float64}(I, 3, 3)
    reset_latency = diagm([1.0, 1.0, 0.0])

    layers = [(OILPUMP_OFF_FREE, 0.0, true), (OILPUMP_OFF_LOCKED, 0.0, false), (OILPUMP_ON_FREE, 1.0, true), (OILPUMP_ON_LOCKED, 1.0, false)]
    location_id(layer, piece) = (layer - 1) * length(pieces) + piece

    locations = ModalLocation[]
    for (layer, pump, free) in layers, (i, (t_start, t_end, rate)) in enumerate(pieces)
        U = Zonotope([1.0, inflow * pump - rate, free ? 0.0 : -1.0], reshape(rate > 0 ? [0.0, fluctuation, 0.0] : zeros(3), 3, 1))
        # The free layers only use the cells at l = 0
        region_faces = [LazySets.HalfSpace(-e(time_axis), -t_start), LazySets.HalfSpace(e(time_axis), t_end)]
        free && push!(region_faces, LazySets.HalfSpace(e(latency_axis), tol))

        last_piece = i == length(pieces)
        edges = [Edge(location_id(layer, last_piece ? 1 : i + 1), HPolyhedron([LazySets.HalfSpace(-e(time_axis), -t_end)]), keep_all, last_piece ? -cycle .* e(time_axis) : zeros(3))]
        actions = ModalAction[]
        if free
            target_layer = pump == 0.0 ? OILPUMP_ON_LOCKED : OILPUMP_OFF_LOCKED
            push!(actions, ModalAction(Action(pump == 0.0 ? 1 : 2, nothing, reset_latency, latency_reset .* e(latency_axis)), location_id(target_layer, i)))
        else
            free_layer = pump == 0.0 ? OILPUMP_OFF_FREE : OILPUMP_ON_FREE
            push!(edges, Edge(location_id(free_layer, i), HPolyhedron([LazySets.HalfSpace(e(latency_axis), 0.0)]), reset_latency, zeros(3)))
        end
        push!(locations, ModalLocation(layer, HPolyhedron(region_faces), A, U, edges, actions))
    end

    unsafe = [HPolyhedron([LazySets.HalfSpace(e(volume_axis), V_min - tol)]), HPolyhedron([LazySets.HalfSpace(-e(volume_axis), -(V_max + tol))])]
    low_corner = [0.0, V_min - margin, 0.0]
    high_corner = [cycle, V_max + margin, latency_top]
    S = Hyperrectangle((low_corner .+ high_corner) ./ 2, (high_corner .- low_corner) ./ 2)
    return ModalHybridSystem(S, unsafe, locations), period
end