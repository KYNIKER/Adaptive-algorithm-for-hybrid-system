using SparseArrays: sparsevec
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityAnalysis: add_dimension
using LazySets


# The oil pump as a ModalHybridSystem (modal.jl). Returns the system and the decision period.
# Continuous state: (t, v, l) = time in the consumption cycle, oil volume, time left of the pump latency.
# Layers: pump off/on x latency free/locked. Each layer has one location per piece of the consumption profile, with
# guards on t between them. The latency runs in the locked layers and ends with an edge to the free layer, the only place
# the pump can switch.
# With outflow_in_A, the consumption is part of the dynamics instead of the input: a fourth state x₄ = 1 (a zero-width
# axis, so a single cell) with v' = -rate x₄ in A, so the input on v is only the pump flow and the consumption's time-varying
# fluctuation. The clocks stay in the input. allowed then has the axes (t, v, l, x₄, layer).
const OILPUMP_OFF_FREE, OILPUMP_OFF_LOCKED, OILPUMP_ON_FREE, OILPUMP_ON_LOCKED = 1, 2, 3, 4
# Consumption pieces (start, end, rate)
const OILPUMP_PIECES = [(0.0, 2.0, 0.0), (2.0, 4.0, 1.2), (4.0, 8.0, 0.0), (8.0, 10.0, 1.2), (10.0, 12.0, 2.5), (12.0, 14.0, 0.0), (14.0, 16.0, 1.7), (16.0, 18.0, 0.5), (18.0, 20.0, 0.0)]

function loadOilPumpModal(; period=0.2, latency_reset=2.0 - 1e-6, latency_top=2.0, outflow_in_A::Bool=false)
    time_axis, volume_axis, latency_axis, one_axis = 1, 2, 3, 4
    n = outflow_in_A ? 4 : 3
    V_min, V_max = 4.9, 25.1
    inflow = 2.2
    fluctuation = 0.1 # l/s, only when the consumption rate is not 0
    margin = max(0.6, 2.7 * period) # more than the largest change in volume over a period (2.6 l/s), so no reach set leaves the grid
    pieces = OILPUMP_PIECES
    cycle = 20.0
    # By default switching sets l just below the latency of 2 s. Decisions are period apart, so the pump is free at exactly
    # the same decision, and the value stays strictly inside a cell of the grid instead of on the border between two.
    tol = 1e-6

    e(i) = [j == i ? 1.0 : 0.0 for j in 1:n]
    keep_all = Matrix{Float64}(I, n, n)
    reset_latency = diagm([j == latency_axis ? 0.0 : 1.0 for j in 1:n])

    layers = [(OILPUMP_OFF_FREE, 0.0, true), (OILPUMP_OFF_LOCKED, 0.0, false), (OILPUMP_ON_FREE, 1.0, true), (OILPUMP_ON_LOCKED, 1.0, false)]
    location_id(layer, piece) = (layer - 1) * length(pieces) + piece

    locations = ModalLocation[]
    for (layer, pump, free) in layers, (i, (t_start, t_end, rate)) in enumerate(pieces)
        A = zeros(n, n)
        outflow_in_A && (A[volume_axis, one_axis] = -rate)
        U = Zonotope([1.0, inflow * pump - (outflow_in_A ? 0.0 : rate), free ? 0.0 : -1.0, zeros(n - 3)...], reshape(rate > 0 ? fluctuation .* e(volume_axis) : zeros(n), n, 1))
        # The free layers only use the cells at l = 0
        region_faces = [LazySets.HalfSpace(-e(time_axis), -t_start), LazySets.HalfSpace(e(time_axis), t_end)]
        free && push!(region_faces, LazySets.HalfSpace(e(latency_axis), tol))

        last_piece = i == length(pieces)
        edges = [Edge(location_id(layer, last_piece ? 1 : i + 1), HPolyhedron([LazySets.HalfSpace(-e(time_axis), -t_end)]), keep_all, last_piece ? -cycle .* e(time_axis) : zeros(n))]
        actions = ModalAction[]
        if free
            target_layer = pump == 0.0 ? OILPUMP_ON_LOCKED : OILPUMP_OFF_LOCKED
            push!(actions, ModalAction(Action(pump == 0.0 ? 1 : 2, nothing, reset_latency, latency_reset .* e(latency_axis)), location_id(target_layer, i)))
        else
            free_layer = pump == 0.0 ? OILPUMP_OFF_FREE : OILPUMP_ON_FREE
            push!(edges, Edge(location_id(free_layer, i), HPolyhedron([LazySets.HalfSpace(e(latency_axis), 0.0)]), reset_latency, zeros(n)))
        end
        push!(locations, ModalLocation(layer, HPolyhedron(region_faces), A, U, edges, actions))
    end

    unsafe = [HPolyhedron([LazySets.HalfSpace(e(volume_axis), V_min)]), HPolyhedron([LazySets.HalfSpace(-e(volume_axis), -(V_max))])]
    low_corner = [0.0, V_min - margin, 0.0, ones(n - 3)...]
    high_corner = [cycle, V_max + margin, latency_top, ones(n - 3)...]
    S = Hyperrectangle((low_corner .+ high_corner) ./ 2, (high_corner .- low_corner) ./ 2)
    return ModalHybridSystem(S, unsafe, locations), period
end

# What the shield allows at l = 0 over the layers (pump off/on x latency free/locked): the bits of allowed are and-ed
# over the layers where the cell is alive, so a choice is shown if it is safe in every pump and latency state from which
# the cell can still be kept safe. A cell is black only when it is dead in every layer. Where waiting is not safe but
# the cell is alive, it gets a darker shade of the colour it has when waiting is safe too. With actions that do nothing
# where they are not offered, turning the pump on when it is on is the same as waiting.
function plot_oil_pump(grid, shield, allowed, lost, info)
    g = shield.granularity
    layers = [allowed[:, :, 1, L] for L in axes(allowed, 4)]
    alive_and(bits...) = (alive=filter(!=(0), collect(bits)); isempty(alive) ? 0 : Int(reduce(&, alive)))
    # bit 0 wait, bit 1 turn on, bit 2 turn off. 8: alive, but no choice is safe in every layer where it is alive.
    m = alive_and.(layers...)
    m[(m .== 0) .& .!all(l->l==0, cat(layers...; dims=3); dims=3)[:, :, 1]] .= 8
    colors = [:black, :steelblue, :darkgoldenrod, :yellow, :orchid4, :orchid, :gray60, :white, :red]
    labels = ["dead", "{}", "{on}", "{wait, on}", "{off}", "{wait, off}",
        "must turn on or off", "anything", "no choice safe in every live layer"]
    println(join(("$l $(count(==(i - 1), m))" for (i, l) in enumerate(labels)), ", "), " of $(length(m)) cells")

    # The heatmap's axes are cell indices: cell i spans [i - 0.5, i + 0.5], so a value x on axis a is at (x - lower)/g + 0.5
    at(x, a) = (x .- shield.lower[a]) ./ g .+ 0.5
    seconds = ceil(shield.lower[1]):2:floor(shield.upper[1])
    liters = ceil(shield.lower[2]):floor(shield.upper[2])
    t_steps = [[p[1] for p in OILPUMP_PIECES]; OILPUMP_PIECES[end][2]]
    rates = [[p[3] for p in OILPUMP_PIECES]; OILPUMP_PIECES[end][3]]
    xlims = (0.5, shield.numCells[1] + 0.5)
    ylims = (0.5, shield.numCells[2] + 0.5)
    fig = heatmap(transpose(m); c=cgrad(colors, length(colors), categorical=true), clim=(0, length(colors) - 1), colorbar=false, title="", #Oil pump shield: choices allowed in every live layer (l = 0)
        xlabel="t (s)", ylabel="V (l)", xticks=(at(seconds, 1), string.(Int.(seconds))), yticks=(at(liters[1:2:end], 2), string.(Int.(liters[1:2:end]))),
        xlims, ylims, legend_position=:topleft, size=(1000, 650), left_margin=3Plots.mm, right_margin=3Plots.mm)
    vline!(fig, at(seconds, 1); color=:gray, lw=0.5, alpha=0.6, label=false)
    hline!(fig, at(liters, 2); color=:gray, lw=0.5, alpha=0.6, label=false)
    for (i, (color, label)) in enumerate(zip(colors, labels))
        any(==(i - 1), m) || continue   # only the choices that occur
        plot!(fig, Float64[], Float64[]; seriestype=:shape, label, color)   # Float64[] avoids a warning
    end
    # The consumption rate on a second axis with ticks on the right border, scaled by 5: a rate c sits 5c liters above the
    # bottom of the plot.
    plot!(twinx(fig), at(t_steps, 1), rates; seriestype=:steppost, color=:red, lw=1.5, label=false, xlims, xticks=false,
        ylims=(0.0, (ylims[2] - ylims[1]) * g / 5), ylabel="consumption (l/s)")
    display(fig)
    return fig
end

# Automatic Synthesis of Robust and Optimal Controllers – An Industrial Case Study
# It looks like this should be modelled as switched dynamics...
# Looks like CoSHY was able to model this with a 4 dimensional state. Thus they only track time, volume of oil in the tank, state of the pump and a latency timer.
# They compute the volume of oil in the tank simply by the flow rate times the time spent in that rate. The act of turning the pump on is simply modelled subtracting the outgoing flow until next time step.
# I think this needs to be modelled as act -> translates set in pump_rate dimension. Then the reachability relation should be computed for each decision period in the time interval.
