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

# The cruise control with a decision at the start of every period and no memory of the last one: a single layer. Waiting
# coasts, action CRUISE_DECELERATE brakes and action CRUISE_ACCELERATE accelerates the ego for the period. For each choice
# the front car (the adversary) brakes, coasts or accelerates, through overlapping locations as in loadCruiseControlModal.
#   - Coasting: three locations that cells start in, one per front choice, restricted to the speeds where the front car
#     may make it. They offer the two actions.
#   - Braking and accelerating: three locations each that no cell starts in (starts = false), entered by the actions. A
#     location coasting with front choice f sends the action to the location with the same f, so a cell gets every front
#     choice its speed allows. Their periods end in the layer, so the next period is a decision again.
# At its speed limit the ego coasts for the rest of the period: an edge sets the speed just inside the limit and goes to
# the coasting location with the same front choice. The sensor range is handled as in loadCruiseControlModal.
function loadCruiseControlModal_t(; period=1.0, acceleration=2.0, v_ego_min=-10.0, v_ego_max=20.0, v_front_min=-8.0, v_front_max=20.0, d_unsafe=5.0, sensor_range=200.0)
    n = 3
    distance, v_ego, v_front = 1, 2, 3
    e(i) = [j == i ? 1.0 : 0.0 for j in 1:n]
    keep_all = Matrix{Float64}(I, n, n)
    keep_all_but(i) = diagm([j == i ? 0.0 : 1.0 for j in 1:n])
    tol = 1e-6
    A = zeros(n, n)
    A[distance, v_ego] = -1.0
    A[distance, v_front] = 1.0

    # The ego's choices: (mode, sign of its acceleration); coasting is waiting
    modes = ((CRUISE_DECELERATE, -1.0), (CRUISE_NEUTRAL, 0.0), (CRUISE_ACCELERATE, 1.0))
    # The front car's choices: (sign of its acceleration, region where it may pick it)
    change = acceleration * period
    fronts = [(-1.0, HPolyhedron([LazySets.HalfSpace(-e(v_front), -(v_front_min + change))])), (0.0, nothing),
        (1.0, HPolyhedron([LazySets.HalfSpace(e(v_front), v_front_max - change)]))]
    location_id(mode, f) = (mode - 1) * length(fronts) + f
    far_reset = sensor_range - (v_front_max - v_ego_min) * period
    layer = 1

    locations = ModalLocation[]
    for (mode, sign) in modes, (f, (front_sign, region)) in enumerate(fronts)
        U = Zonotope([0.0, sign * acceleration, front_sign * acceleration], zeros(n, 1))
        edges = Edge[]
        sign > 0 && push!(edges, Edge(location_id(CRUISE_NEUTRAL, f), HPolyhedron([LazySets.HalfSpace(-e(v_ego), -v_ego_max)]), keep_all_but(v_ego), (v_ego_max - tol) .* e(v_ego)))
        sign < 0 && push!(edges, Edge(location_id(CRUISE_NEUTRAL, f), HPolyhedron([LazySets.HalfSpace(e(v_ego), v_ego_min)]), keep_all_but(v_ego), (v_ego_min + tol) .* e(v_ego)))
        push!(edges, Edge(location_id(mode, f), HPolyhedron([LazySets.HalfSpace(-e(distance), -sensor_range)]), keep_all_but(distance), far_reset .* e(distance)))
        coasting = mode == CRUISE_NEUTRAL
        actions = coasting ? [ModalAction(Action(m, nothing, keep_all, zeros(n)), location_id(m, f)) for (m, _) in modes if m != CRUISE_NEUTRAL] : ModalAction[]
        push!(locations, ModalLocation(layer, coasting ? region : nothing, A, U, edges, actions; starts=coasting))
    end

    unsafe = [HPolyhedron([LazySets.HalfSpace(e(distance), d_unsafe)])]
    low_corner = [0.0, v_ego_min, v_front_min]
    high_corner = [sensor_range, v_ego_max, v_front_max]
    S = Hyperrectangle((low_corner .+ high_corner) ./ 2, (high_corner .- low_corner) ./ 2)
    return ModalHybridSystem(S, unsafe, locations), period
end

# What the shield allows at a given ego speed, as the front car's speed (x) against the distance (y), with one panel per
# mode the ego is in now (its layer). Each cell is coloured by the modes the ego may be in next period: waiting keeps the
# current mode, and an action switches to its mode. The colours are those of plot_oil_pump, by the same bit order
# (bit 0 brake, bit 1 coast, bit 2 accelerate), with black for a dead cell.
# For loadCruiseControlModal_t, which has a single layer where waiting coasts, pass wait_mode=CRUISE_NEUTRAL.
function plot_cruise(grid, allowed, info; v_ego=0.0, wait_mode=nothing)
    g = grid.granularity
    ids = info.ctx.action_ids # action position -> mode
    j = floor(Int, (v_ego - grid.lower[2]) / g) + 1
    v_lo = grid.lower[2] + (j - 1) * g
    function next_modes(bits, current)
        bits == 0 && return 0
        m = (bits & 0x01) != 0 ? 1 << (current - 1) : 0
        for (ai, id) in enumerate(ids)
            (bits >> ai) & 0x01 != 0 && (m |= 1 << (id - 1))
        end
        return m
    end
    colors = [:black, :purple, :darkgoldenrod, :yellow, :orchid4, :orchid, :gray60, :white]
    labels = ["dead", "{brake}", "{coast}", "{brake, coast}", "{accelerate}", "{brake, accelerate}", "{coast, accelerate}", "{brake, coast, accelerate}"]
    modes = ["braking", "coasting", "accelerating"]

    # The heatmap's axes are cell indices: cell i spans [i - 0.5, i + 0.5], so a value x on axis a is at (x - lower)/g + 0.5
    at(x, a) = (x .- grid.lower[a]) ./ g .+ 0.5
    speeds = ceil(grid.lower[3]):4:floor(grid.upper[3])
    distances = ceil(grid.lower[1]):(grid.upper[1]-grid.lower[1]>100 ? 20 : 10):floor(grid.upper[1])
    xlims = (0.5, size(allowed, 3) + 0.5)
    ylims = (0.5, size(allowed, 1) + 0.5)
    panels = []
    present = Set{Int}()
    for layer in axes(allowed, 4)
        current = wait_mode === nothing ? layer : wait_mode
        m = [next_modes(allowed[i, j, k, layer], current) for i in axes(allowed, 1), k in axes(allowed, 3)]
        union!(present, m)
        name = wait_mode === nothing ? "Currently $(modes[layer])" : "At a decision"
        println("$name: ", join(("$(labels[c+1]) $(count(==(c), m))" for c in 0:7 if any(==(c), m)), ", "))
        fig = heatmap(m; c=cgrad(colors, length(colors), categorical=true), clim=(0, length(colors) - 1), colorbar=false,
            title=name, xlabel="front car speed (m/s)", ylabel="distance (m)", titlefontsize=10,
            xticks=(at(speeds, 3), string.(Int.(speeds))), yticks=(at(distances, 1), string.(Int.(distances))), xlims, ylims, legend=false)
        vline!(fig, at(speeds, 3); color=:gray, lw=0.5, alpha=0.6, label=false)
        hline!(fig, at(distances, 1); color=:gray, lw=0.5, alpha=0.6, label=false)
        push!(panels, fig)
    end
    legend = plot(; framestyle=:none, legend=:left, foreground_color_legend=nothing)
    for c in sort(collect(present))
        plot!(legend, Float64[], Float64[]; seriestype=:shape, color=colors[c+1], linecolor=:gray60, label=labels[c+1])   # Float64[] avoids a warning
    end
    k = length(panels)
    width = 500 * k + 250 # the legend keeps 250 px
    fig = plot(panels..., legend; layout=Plots.grid(1, k + 1, widths=[fill(500 / width, k); 250 / width]), size=(width, 650), left_margin=8Plots.mm, bottom_margin=6Plots.mm,
        plot_title="v_ego in [$v_lo, $(v_lo + g)) m/s", plot_titlefontsize=11)
    display(fig)
    return fig
end
