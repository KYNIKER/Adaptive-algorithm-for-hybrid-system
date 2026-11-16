#using SparseArrays: sparsevec
#using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
#using ReachabilityAnalysis: add_dimension
using LazySets



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
        LazySets.HalfSpace(sparsevec([v_ego], [1.], 3), v_ego_max) # v_ego <= 0
    ])

    minspeedlimit_ego = HPolyhedron([
        LazySets.HalfSpace(sparsevec([v_ego], [-1.], 3), v_ego_min) # v_ego >= 0
    ])

    maxspeedlimit_front = HPolyhedron([
        LazySets.HalfSpace(sparsevec([v_front], [1.], 3), v_front_max) # v_front <= 0
    ])

    minspeedlimit_front = HPolyhedron([
        LazySets.HalfSpace(sparsevec([v_front], [-1.], 3), v_front_min) # v_front >= 0
    ])

    speedlimit_ego = intersection(minspeedlimit_ego, maxspeedlimit_ego)
    speedlimit_front = intersection(minspeedlimit_front, maxspeedlimit_front)

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
        actions = ModalAction[]
        edges = []
        Al = copy(A)
        Al[v_ego] = eg * 2.0
        Al[v_front] = fr * 2.0
        if eg == 0 && fr == 0



            eg > 0 && push!(edges, Edge(lookupDict[(0, fr)], maxspeedlimit_ego, dia, zeros(dims)))
            eg < 0 && push!(edges, Edge(lookupDict[(0, fr)], minspeedlimit_ego, dia, zeros(dims)))
            fr > 0 && push!(edges, Edge(lookupDict[(eg, 0)], maxspeedlimit_front, dia, zeros(dims)))
            fr < 0 && push!(edges, Edge(lookupDict[(eg, 0)], minspeedlimit_front, dia, zeros(dims)))
            for j in -1:1:1
                j == 0 && continue
                push!(actions, ModalAction(Action(j, speedlimit_ego, dia, zeros(dims)), lookupDict[(0, j)]))
            end
        else




            eg > 0 && push!(edges, Edge(lookupDict[(0, fr)], maxspeedlimit_ego, dia, zeros(dims)))
            eg < 0 && push!(edges, Edge(lookupDict[(0, fr)], minspeedlimit_ego, dia, zeros(dims)))
            fr > 0 && push!(edges, Edge(lookupDict[(eg, 0)], maxspeedlimit_front, dia, zeros(dims)))
            fr < 0 && push!(edges, Edge(lookupDict[(eg, 0)], minspeedlimit_front, dia, zeros(dims)))

            if fr == 0
                for j in -1:1:1
                    j == fr && continue
                    push!(actions, ModalAction(Action(j, speedlimit_front, dia, zeros(dims)), lookupDict[(eg, j)]))
                end
            end
        end
        push!(locations, ModalLocation(i, nothing, Al, zero_input, edges, actions))
    end
    #edges::Vector{Edge} = [Edge(1, guard, jumpMatrix, zeros(3))]
    #constraint = LazySets.API.intersection(constraint, S)


    #Act = [backwards, forwards]

    #H = EuclideanHybridSystem(S, Zonotope(zeros(3), zeros(3, 3)), [constraint], A, edges, Act)

    T = 1.0

    return ModalHybridSystem(S, [constraint], locations), T
end

const CRUISE_DECELERATE, CRUISE_NEUTRAL, CRUISE_ACCELERATE = 1, 2, 3

# The adaptive cruise control of "SOS: Safe, Optimal and Small Strategies for Hybrid Markov Decision Processes" as a
# ModalHybridSystem (modal.jl). Returns the system and the decision period.
# Continuous state: (d, v_ego, v_front) = distance to the front car and the speeds of the two cars, with d' = v_front - v_ego.
# Layers: the ego's acceleration (decelerate, neutral, accelerate). The controller picks it every period with the action of
# the same id, and waiting keeps it. At its speed limit the ego stops accelerating: an edge sets the speed just inside the
# limit and goes to the neutral layer.
# The front car picks its acceleration in {-acceleration, 0, acceleration} at random every period, so for a shield it is
# adversarial, never an action: each layer has a location for each of the three, and the locations overlap, so the shield
# takes all three as successors. As in the paper, where the speeds change by acceleration * period from multiples of it,
# the front car only picks an acceleration that keeps its speed within its limits for the whole period, through the
# regions of the slowing and speeding locations. Their borders must be cell borders: the limits and acceleration * period
# multiples of the granularity. Both speeds are then clock axes, which the shield tracks exactly. (A time-varying input in
# [-acceleration, acceleration] on v_front is not: on a coarse grid its end interval also touches the cell below, so the
# front car could slow down faster than the ego, and every cell died.)
# Beyond the sensor range the front car is out of sight, which is not unsafe. An edge at d = sensor_range moves it back to
# the distance it can cover in a period at the largest relative speed, as if it were closer than it is, so the shield stays
# sound. Moving it to just below the range does not work: the edge would fire again at once, and a jump at the end of the
# period takes no time, so the shield would jump until max_depth and call the transition bad.
function loadCruiseControlModal(; period=1.0, acceleration=2.0, v_ego_min=-10.0, v_ego_max=20.0, v_front_min=-8.0, v_front_max=20.0, d_unsafe=5.0, sensor_range=200.0)
    n = 3
    distance, v_ego, v_front = 1, 2, 3
    e(i) = [j == i ? 1.0 : 0.0 for j in 1:n]
    keep_all = Matrix{Float64}(I, n, n)
    keep_all_but(i) = diagm([j == i ? 0.0 : 1.0 for j in 1:n])
    tol = 1e-6
    A = zeros(n, n)
    A[distance, v_ego] = -1.0
    A[distance, v_front] = 1.0

    modes = ((CRUISE_DECELERATE, -1.0), (CRUISE_NEUTRAL, 0.0), (CRUISE_ACCELERATE, 1.0))
    # The front car's choices: (sign of its acceleration, region where it may pick it). 
    change = acceleration * period
    fronts = [(-1.0, HPolyhedron([LazySets.HalfSpace(-e(v_front), -(v_front_min + change))])), (0.0, nothing),
        (1.0, HPolyhedron([LazySets.HalfSpace(e(v_front), v_front_max - change)]))]
    location_id(mode, f) = (mode - 1) * length(fronts) + f
    far_reset = sensor_range - (v_front_max - v_ego_min) * period

    locations = ModalLocation[]
    for (mode, sign) in modes, (f, (front_sign, region)) in enumerate(fronts)
        U = Zonotope([0.0, sign * acceleration, front_sign * acceleration], zeros(n, 1))
        edges = Edge[]
        sign > 0 && push!(edges, Edge(location_id(CRUISE_NEUTRAL, f), HPolyhedron([LazySets.HalfSpace(-e(v_ego), -v_ego_max)]), keep_all_but(v_ego), (v_ego_max - tol) .* e(v_ego)))
        sign < 0 && push!(edges, Edge(location_id(CRUISE_NEUTRAL, f), HPolyhedron([LazySets.HalfSpace(e(v_ego), v_ego_min)]), keep_all_but(v_ego), (v_ego_min + tol) .* e(v_ego)))
        # This one should actually be setting the speed of v_front to be in [v_ego_min, v_ego]
        push!(edges, Edge(location_id(mode, f), HPolyhedron([LazySets.HalfSpace(-e(distance), -sensor_range)]), keep_all_but(distance), far_reset .* e(distance)))

        actions = [ModalAction(Action(m, nothing, keep_all, zeros(n)), location_id(m, f)) for (m, _) in modes if m != mode]
        @show mode, sign, front_sign, actions
        push!(locations, ModalLocation(mode, region, A, U, edges, actions))
    end

    unsafe = [HPolyhedron([LazySets.HalfSpace(e(distance), d_unsafe)])]
    low_corner = [0.0, v_ego_min, v_front_min]
    high_corner = [sensor_range, v_ego_max, v_front_max]
    S = Hyperrectangle((low_corner .+ high_corner) ./ 2, (high_corner .- low_corner) ./ 2)
    return ModalHybridSystem(S, unsafe, locations), period
end

function loadCruiseControlModal_t(; period=1.0, acceleration=2.0, v_ego_min=-10.0, v_ego_max=20.0, v_front_min=-8.0, v_front_max=20.0, d_unsafe=5.0, sensor_range=200.0)
    n = 3
    distance, v_ego, v_front = 1, 2, 3
    e(i) = [j == i ? 1.0 : 0.0 for j in 1:n]
    keep_all = Matrix{Float64}(I, n, n)
    keep_all_but(i) = diagm([j == i ? 0.0 : 1.0 for j in 1:n])
    tol = 1e-2
    A = zeros(n, n)
    A[distance, v_ego] = -1.0
    A[distance, v_front] = 1.0

    modes = ((CRUISE_DECELERATE, -1.0), (CRUISE_NEUTRAL, 0.0), (CRUISE_ACCELERATE, 1.0))
    # The front car's choices: (sign of its acceleration, region where it may pick it). 
    change = acceleration * period
    fronts = [(-1.0, HPolyhedron([LazySets.HalfSpace(-e(v_front), -(v_front_min + change))])), (0.0, nothing),
        (1.0, HPolyhedron([LazySets.HalfSpace(e(v_front), v_front_max - change)]))]
    location_id(mode, f) = (mode - 1) * length(fronts) + f
    far_reset = sensor_range - (v_front_max - v_ego_min) * period

    locations = ModalLocation[]
    speed_regulation = keep_all_but(v_ego)
    speed_regulation[v_front, v_ego] = 1.0
    for (mode, sign) in modes, (f, (front_sign, region)) in enumerate(fronts)
        if sign != 0 || front_sign != 0
            U = Zonotope([0.0, sign * acceleration, front_sign * acceleration], zeros(n, 1))
            edges = Edge[]
            sign > 0 && push!(edges, Edge(location_id(CRUISE_NEUTRAL, f), HPolyhedron([LazySets.HalfSpace(-e(v_ego), -v_ego_max)]), keep_all_but(v_ego), (v_ego_max - tol) .* e(v_ego)))
            sign < 0 && push!(edges, Edge(location_id(CRUISE_NEUTRAL, f), HPolyhedron([LazySets.HalfSpace(e(v_ego), v_ego_min)]), keep_all_but(v_ego), (v_ego_min + tol) .* e(v_ego)))
            # This one should actually be setting the speed of v_front to be in [v_ego_min, v_ego]

            push!(edges, Edge(location_id(mode, 2), HPolyhedron([LazySets.HalfSpace(-e(distance), -sensor_range)]), speed_regulation, far_reset .* e(distance)))

            actions = []#[ModalAction(Action(m, nothing, keep_all, zeros(n)), location_id(m, f)) for (m, _) in modes if m != mode]
            #@show mode, sign, front_sign, actions
            push!(locations, ModalLocation(mode, region, A, U, edges, actions))
        else
            U = Zonotope([0.0, sign * acceleration, front_sign * acceleration], zeros(n, 1))
            edges = Edge[]
            #sign > 0 && push!(edges, Edge(location_id(CRUISE_NEUTRAL, f), HPolyhedron([LazySets.HalfSpace(-e(v_ego), -v_ego_max)]), keep_all_but(v_ego), (v_ego_max - tol) .* e(v_ego)))
            #sign < 0 && push!(edges, Edge(location_id(CRUISE_NEUTRAL, f), HPolyhedron([LazySets.HalfSpace(e(v_ego), v_ego_min)]), keep_all_but(v_ego), (v_ego_min + tol) .* e(v_ego)))
            # This one should actually be setting the speed of v_front to be in [v_ego_min, v_ego]
            #push!(edges, Edge(location_id(mode, f), HPolyhedron([LazySets.HalfSpace(-e(distance), -sensor_range)]), keep_all_but(distance), far_reset .* e(distance)))

            actions = [ModalAction(Action(m, nothing, keep_all, zeros(n)), location_id(m, fi)) for (m, _) in modes, (fi, _) in enumerate(fronts) if m != 2 || fi != 2]
            #@show mode, sign, front_sign, actions
            push!(locations, ModalLocation(mode, region, A, U, edges, actions))
        end
    end

    unsafe = [HPolyhedron([LazySets.HalfSpace(e(distance), d_unsafe)])]
    low_corner = [0.0, v_ego_min, v_front_min]
    high_corner = [sensor_range, v_ego_max, v_front_max]
    S = Hyperrectangle((low_corner .+ high_corner) ./ 2, (high_corner .- low_corner) ./ 2)
    return ModalHybridSystem(S, unsafe, locations), period
end