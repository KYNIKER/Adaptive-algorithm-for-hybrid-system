#=
The modal shield of modal.jl, reduced to its method. Every axis is treated alike. Every cell of every layer is computed
once, then the bad events are pushed backwards (flat.jl). Uses the system types and some helpers of modal.jl.

Cells are half-open, [l, l + g) per axis, as for a floor lookup of the state. A set reaches the cells it touches
(touching_cells), within the bounds that are known to hold for its states at the end of the period (Bounds). These
bounds are exact where they are tracked, so a cell that only touches their open end is dropped: otherwise a set that
moves by whole cells would also reach a neighbour every period. They are also what the exit check uses, as the sets
after a jump are bloated.
  - Flow and actions: the end state of x₀ in the cell is T x₀ + c, exactly, on the axes where the input is constant.
    The interval of axis i is then exact, with an open upper end when it needs some x₀ at the open end of the cell.
    An open end is moved down by OPEN_MARGIN, so all bounds are closed.
  - A guard jumps at a time that depends on the state. An axis keeps its bounds when its end value does not depend on
    that time (see jumped_bounds), and loses them otherwise.

One period of a set X in a location (reach!):
  - ReACT's adaptive time steps find the first step of δ⁻ at which the flowpipe of X touches a guard or an unsafe set.
  - From there the flowpipe is walked in pieces of δ⁻. A piece touching an unsafe set makes the transition bad. A piece
    touching a guard jumps, and the jumped set continues in the target location for its own remaining time.
  - The states that have not jumped end the period in the location.
Leaving the state space is checked on the sets at the end of the period, where a part beyond a face that lies inside a
guard has jumped. So a state that leaves and comes back within the period is not caught.
=#

const DenseZonotope = Zonotope{Float64,Vector{Float64},Matrix{Float64}}
const Halfspaces = Tuple{Vector{Vector{Float64}},Vector{Float64}}
const Affine = Tuple{Matrix{Float64},Vector{Float64}} # x = T x₀ + c, for x₀ in the start cell

const Bounds = Tuple{Vector{Float64},Vector{Float64}} # per axis, every state is in [lo, hi] at the end of the period
const PeriodEnd = Tuple{Int,DenseZonotope,Bounds} # location, set, bounds of its states

# How far an open end is moved down: more than the MODAL_TOL with which bounds are compared, in cells and absolutely.
open_margin(g::Float64) = 4MODAL_TOL * max(1.0, g)

struct SimpleLocation{D}
    dynamics::Dynamics
    layer::Int
    guards::Vector{Halfspaces}
    rate::Vector{Float64} # per axis its constant rate where it only translates, NaN for the others
    starts::BitArray{D} # the cells a period starts in this location
    discs::Dict{Matrix{Float64},CellDiscretization{D}} # per shape, the discretization of a whole period
end

struct SimpleContext{D,S}
    system::S
    grid::Grid{D}
    steps::Int
    δ⁻::Float64
    δ⁺::Float64
    alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg
    maxOrder::Int
    reduceOrder::Int
    max_depth::Int
    locations::Vector{SimpleLocation{D}}
    layers::Vector{Vector{Int}} # per layer, its locations
    unsafe::Vector{Halfspaces}
    statespace::Halfspaces
    cell::Matrix{Float64} # generators of a cell
    action_ids::Vector{Int}
    may_act::Vector{Vector{BitArray{D}}} # per location and action: some state of the cell is inside the guard
    must_act::Vector{Vector{BitArray{D}}} # every state of the cell is inside it
    covers::CoverCache{D}
    exit_bad_for_flow::Bool
    exit_bad_for_actions::Bool
end

function SimpleContext(system::ModalHybridSystem, p::Float64, g::Float64, δ⁻, δ⁺, alg, maxOrder, reduceOrder, max_depth, exit_bad_for_flow, exit_bad_for_actions)
    steps = round(Int, p / δ⁻)
    isapprox(steps * δ⁻, p; atol=MODAL_TOL) || throw(ArgumentError("The period $p must be a multiple of δ⁻ = $δ⁻"))
    isinteger(log2(δ⁺ / δ⁻)) || throw(ArgumentError("δ⁺ / δ⁻ must be a power of two"))
    grid = Grid(system.statespace, g)
    D = grid.dimension

    # On an axis of a single cell the cell is the state space: a point for a constant.
    half = [n == 1 ? (h - l) / 2 : g / 2 for (n, l, h) in zip(grid.numCells, low(system.statespace), high(system.statespace))]
    cell = Matrix{Float64}(genmat(remove_zero_generators(Zonotope(zeros(D), diagm(half)))))

    shared = Dict{Tuple{Matrix{Float64},Vector{Float64},Matrix{Float64}},Dynamics}() # locations with the same (A, U)
    locations = map(system.locations) do loc
        dyn = get!(() -> Dynamics(loc.A, loc.U, steps, δ⁻, δ⁺, alg, maxOrder, reduceOrder), shared, (loc.A, LazySets.center(loc.U), genmat(loc.U)))
        region = loc.region === nothing ? nothing : getHalfSpaceProjections(loc.region)
        starts = BitArray(loc.starts && box_overlaps_region(region, cell_low(grid, c), cell_low(grid, c) .+ g) for c in CartesianIndices(cellshape(grid)))
        guards = Halfspaces[getHalfSpaceProjections(e.guard) for e in loc.edges]
        rate = [iszero(view(loc.A, i, :)) && iszero(view(genmat(loc.U), i, :)) ? LazySets.center(loc.U)[i] : NaN for i in 1:D]
        SimpleLocation(dyn, loc.layer, guards, rate, starts, Dict{Matrix{Float64},CellDiscretization{D}}())
    end
    layers = [[li for (li, loc) in enumerate(locations) if loc.layer == layer] for layer in 1:maximum(loc.layer for loc in locations)]

    action_ids = sort!(unique(a.action.id for loc in system.locations for a in loc.actions))
    may_act = [[falses(cellshape(grid)) for _ in action_ids] for _ in locations]
    must_act = [[falses(cellshape(grid)) for _ in action_ids] for _ in locations]
    for (li, loc) in enumerate(system.locations), a in loc.actions
        all(b -> b.action.id != a.action.id || b.action.guard == a.action.guard, loc.actions) ||
            throw(ArgumentError("location $li has actions with id $(a.action.id) and different guards"))
        ai = searchsortedfirst(action_ids, a.action.id)
        may_act[li][ai], must_act[li][ai] = guard_masks(grid, a.action.guard)
    end

    statespace = getHalfSpaceProjections(hyperrectangle_to_HPolyhedron(system.statespace))
    unsafe = Halfspaces[getHalfSpaceProjections(H) for H in system.unsafe]
    return SimpleContext(system, grid, steps, δ⁻, δ⁺, alg, maxOrder, reduceOrder, max_depth, locations, layers, unsafe, statespace, cell,
        action_ids, may_act, must_act, CoverCache{D}(), exit_bad_for_flow, exit_bad_for_actions)
end

cell_low(grid::Grid, c::CartesianIndex) = grid.lower .+ (index_vector(c) .- 1) .* grid.granularity

# The locations the cell c starts a period in, in its layer.
starting(ctx::SimpleContext, c::CartesianIndex, layer::Int) = filter(li -> ctx.locations[li].starts[c], ctx.layers[layer])

can_act(ctx::SimpleContext, c::CartesianIndex, layer::Int, ai::Int) = any(li -> ctx.may_act[li][ai][c], starting(ctx, c, layer))

flowed(dyn::Dynamics, X, k::Int) = minkowski_sum(linear_map(dyn.Φ[k+1], X), dyn.V[k+1])

# The bounds of the cell Z mapped by x = T x₀ + c and flowed for a whole period. The upper end of an axis is open when
# some x₀ⱼ with Tᵢⱼ > 0 is at the open end of the cell, which an axis of one cell has not. An input that varies on an
# axis leaves it untracked: [-Inf, Inf].
function end_bounds(ctx::SimpleContext, dyn::Dynamics, (T, c)::Affine, Z)
    Φ, V = dyn.Φ[ctx.steps+1], dyn.V[ctx.steps+1]
    T′, c′ = Φ * T, Φ * c + LazySets.center(V)
    l, u = low(Z)', high(Z)'
    lo = c′ + vec(sum(min.(T′ .* l, T′ .* u); dims=2))
    hi = c′ + vec(sum(max.(T′ .* l, T′ .* u); dims=2))
    for i in axes(T′, 1)
        any(j -> T′[i, j] > 0 && ctx.grid.numCells[j] > 1, axes(T′, 2)) && (hi[i] -= open_margin(ctx.grid.granularity))
        iszero(view(genmat(V), i, :)) || ((lo[i], hi[i]) = (-Inf, Inf))
    end
    return lo, hi
end

# The bounds after a jump x := M x + j through the guard, from a location with rates `from` to one with rates `to`, for
# states with between D′ and D′ + W′ steps left. The time of the jump depends on the state, so only these axes keep
# bounds:
#   - translated at the same rate on both sides: it ends where it would have without the jump, plus j. A halfspace of the
#     guard on this axis alone also limits it, as the states were inside the guard and then moved for the time left.
#   - reset to j: it ends at j plus its rate times the time left, which a state can have exactly, so its end is closed.
# Returns nothing if no state of the bounds can take the jump.
function jumped_bounds((lo, hi)::Bounds, guard::Halfspaces, M::Matrix{Float64}, j::Vector{Float64}, from::Vector{Float64}, to::Vector{Float64}, D′::Int, W′::Int, δ⁻::Float64)
    lo′, hi′ = fill(-Inf, length(j)), fill(Inf, length(j))
    for i in eachindex(j)
        isnan(to[i]) && continue
        least, most = minmax(to[i] * D′ * δ⁻, to[i] * (D′ + W′) * δ⁻)
        if iszero(view(M, i, :))
            lo′[i], hi′[i] = j[i] + least, j[i] + most
        elseif from[i] == to[i] && M[i, i] == 1 && count(!iszero, view(M, i, :)) == 1
            l, h = lo[i], hi[i]
            for (a, c) in zip(guard...)
                (a[i] != 0 && count(!iszero, a) == 1) || continue
                a[i] > 0 ? (h = min(h, c / a[i] + most)) : (l = max(l, c / a[i] + least))
            end
            l > h + MODAL_TOL && return nothing
            lo′[i], hi′[i] = l + j[i], h + j[i]
        end
    end
    return lo′, hi′
end

# Per axis, the cells the bounds allow, as keep_lo, keep_hi of touching_cells: none for bounds below the grid.
function cell_bounds(grid::Grid{D}, (lo, hi)::Bounds) where D
    index(x, i) = floor(Int, clamp((x - grid.lower[i]) / grid.granularity + MODAL_TOL, -1, grid.numCells[i])) + 1
    keep_lo = ntuple(i -> index(lo[i], i), D)
    return keep_lo, ntuple(i -> max(index(hi[i], i), keep_lo[i]), D)
end

# ReACT's discretization of X, from the one of its shape when X flows a whole period. Only outside a jump: that one is
# overwritten by every call, and the caller still reads it after the jumps return.
function discretize(ctx::SimpleContext, li::Int, X, depth::Int, δ⁺::Float64)
    loc, dyn = ctx.locations[li], ctx.locations[li].dynamics
    if depth == 0
        G = genmat(X)
        cd = get!(() -> period_discretization(dyn, ctx.system.locations[li].U, G, ctx.δ⁻, δ⁺), loc.discs, G)
        return discretize_cell!(cd, LazySets.center(X))
    end
    return newReACTDiscretizePlus(X, ctx.δ⁻, δ⁺, dyn.A, dyn.P2A_abs, dyn.step_Φ, dyn.step_input, ctx.alg, ctx.maxOrder, ctx.reduceOrder)
end

# Appends to out what X reaches at the end of the period. X has between D and D + W steps of δ⁻ left, depending on the
# state: W > 0 after a jump. bounds hold for its states at the end of the period. Returns false if the transition is bad.
function reach!(out::Vector{PeriodEnd}, ctx::SimpleContext, li::Int, X, D::Int, W::Int, exit_bad::Bool, bounds::Bounds, depth::Int=0)
    depth > ctx.max_depth && return false # jumping again and again
    total = D + W
    total == 0 && (push!(out, (li, X, bounds)); return true)
    loc = ctx.locations[li]
    dyn = loc.dynamics
    δ⁺ = min(ctx.δ⁺, ctx.δ⁻ * 2^floor(Int, log2(total))) # no step longer than the time left
    disc = discretize(ctx, li, X, depth, δ⁺)
    k0 = first_touch(dyn, disc, Iterators.flatten((loc.guards, ctx.unsafe)), total, ctx.δ⁻, δ⁺)

    for k in min(k0, W == 0 ? total : D):(total-1)
        piece = minkowski_sum(linear_map(dyn.Φ[k+1], disc[ctx.δ⁻]), dyn.V[k+1])
        if k >= k0
            any(set -> some_inside(piece, set...), ctx.unsafe) && return false
            all_jumped = false
            for (e, guard) in zip(ctx.system.locations[li].edges, loc.guards)
                clipped = some_inside(piece, guard...) ? clip_to_polyhedron(piece, guard) : nothing
                clipped === nothing && continue
                # A state of the piece has spent between k and k + 1 steps here.
                D′, W′ = max(D - k - 1, 0), total - k - max(D - k - 1, 0)
                b = jumped_bounds(bounds, guard, e.jumpMatrix, e.jumpVector, loc.rate, ctx.locations[e.targetLoc].rate, D′, W′, ctx.δ⁻)
                b === nothing && continue
                jumped = jump_set(clipped, e.jumpMatrix, e.jumpVector, ctx.maxOrder, ctx.reduceOrder)
                reach!(out, ctx, e.targetLoc, jumped, D′, W′, exit_bad, b, depth + 1) || return false
                all_jumped |= all_inside(piece, guard...)
            end
            all_jumped && return true
        end
        W > 0 && k >= D && push!(out, (li, piece, bounds)) # the states with k + 1 steps left end in this piece
    end
    exit_bad && any(k -> leaves(ctx, loc, flowed(dyn, X, k), bounds), unique((D, total))) && return false
    W == 0 && push!(out, (li, flowed(dyn, X, total), bounds))
    return true
end

# True if part of the set lies beyond a face of the state space and outside every guard of the location: inside a guard
# it has jumped. On a tracked axis the bounds decide whether it lies beyond. The guard's halfspace that borders the face
# holds beyond it by construction, and is skipped: the clip is not exact on it.
function leaves(ctx::SimpleContext, loc::SimpleLocation, set, (lo, hi)::Bounds)
    for (a, b) in zip(ctx.statespace...)
        i = findfirst(!iszero, a)
        bound = a[i] * (a[i] > 0 ? hi[i] : lo[i])
        min(bound, ρ(a, set)) > b + MODAL_TOL || continue
        beyond = clip_to_polyhedron(set, ([-a], [-b]))
        if beyond === nothing
            isfinite(bound) && return true # the bounds are exact, the set is not
            continue
        end
        in_guard(guard) = all(((c, d),) -> (isapprox(c, -a; atol=MODAL_TOL) && isapprox(d, -b; atol=MODAL_TOL)) || ρ(c, beyond) <= d, zip(guard...))
        any(in_guard, loc.guards) || return true
    end
    return false
end

# The cells that the transitions (location, set, exit_bad, map from the cell Z) reach in one period, or nothing if one of
# them is bad.
function reached(ctx::SimpleContext, transitions, Z)
    out = PeriodEnd[]
    for (li, X, exit_bad, affine) in transitions
        reach!(out, ctx, li, X, ctx.steps, 0, exit_bad, end_bounds(ctx, ctx.locations[li].dynamics, affine, Z)) || return nothing
    end
    return unique!([CartesianIndex(c, ctx.locations[li].layer) for (li, X, b) in out for c in touching_cells(ctx.grid, X, ctx.covers, cell_bounds(ctx.grid, b)...)])
end

# The successors of cell c of a layer that starts in locations locs: of waiting, and per action nothing where it is not
# enabled or bad. An action tried outside its guard, or in a location that does not offer it, does nothing: there the
# state just flows, so the action also has the successors of waiting.
function cell_successors(ctx::SimpleContext, c::CartesianIndex, locs::Vector{Int})
    Z = Zonotope(cell_low(ctx.grid, c) .+ ctx.grid.granularity / 2, ctx.cell)
    D = length(c)
    waiting = reached(ctx, ((li, Z, ctx.exit_bad_for_flow, (Matrix{Float64}(I, D, D), zeros(D))) for li in locs), Z)
    acts = map(eachindex(ctx.action_ids)) do ai
        taking = filter(li -> ctx.may_act[li][ai][c], locs)
        isempty(taking) && return nothing
        outcomes = ((a.target, jump_set(Z, a.action.jumpMatrix, a.action.jumpVector, ctx.maxOrder, ctx.reduceOrder), ctx.exit_bad_for_actions,
                     (a.action.jumpMatrix, a.action.jumpVector))
                    for li in taking for a in ctx.system.locations[li].actions if a.action.id == ctx.action_ids[ai])
        succ = reached(ctx, outcomes, Z)
        all(li -> ctx.must_act[li][ai][c], locs) && return succ
        return succ === nothing || waiting === nothing ? nothing : unique!(vcat(succ, waiting))
    end
    return waiting, acts
end

# Shield for a ModalHybridSystem, like ReACTed_reachable_cell_modal without its optimizations and clock axes, and offline
# only. Returns grid, allowed, lost, info as every shield (see choice_bits in flat.jl); info.ctx is the context.
function ReACTed_reachable_cell_modal_simple(system::ModalHybridSystem, p, granularity, δ⁻::Float64, δ⁺::Float64,
    alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5;
    act_needs_safe_flow::Bool=true, exit_bad_for_flow::Bool=true, exit_bad_for_actions::Bool=true, max_depth::Int=16, verbose::Bool=true)
    ctx = SimpleContext(system, Float64(p), Float64(granularity), δ⁻, δ⁺, alg, maxOrder, reduceOrder, max_depth, exit_bad_for_flow, exit_bad_for_actions)
    return simple_shield(ctx, act_needs_safe_flow, verbose)
end

function simple_shield(ctx::SimpleContext{D}, act_needs_safe_flow::Bool, verbose::Bool) where D
    grid = ctx.grid
    shape = (cellshape(grid)..., length(ctx.layers))
    cells, ids = CartesianIndices(shape), LinearIndices(shape)
    n, n_acts = length(ids), length(ctx.action_ids)
    store = SuccessorStore(n, [falses(n) for _ in 1:n_acts])
    store.flow_empty_ok = !ctx.exit_bad_for_flow
    store.act_empty_ok = !ctx.exit_bad_for_actions
    s = reset!(FlatShieldState(n_acts, act_needs_safe_flow), n)
    split(k) = CartesianIndex(ntuple(i -> k[i], D)), k[D+1]
    as_ids(succ) = succ === nothing ? Int32[] : Int32[ids[k] for k in succ]

    for H in ctx.system.unsafe, c in touching_cells(grid, H), layer in eachindex(ctx.layers)
        seed_constraint!(s, Int32(ids[c, layer]))
        store.flow_none[ids[c, layer]] = true
    end
    for k in cells
        id = Int32(ids[k])
        c, layer = split(k)
        locs = starting(ctx, c, layer)
        (s.computed[id] || isempty(locs)) && continue
        waiting, acts = cell_successors(ctx, c, locs)
        store.flow_none[id] = waiting === nothing
        set_successors!(store, 0, id, as_ids(waiting))
        for ai in 1:n_acts
            store.act_enabled[ai][id] = acts[ai] !== nothing
            set_successors!(store, ai, id, as_ids(acts[ai]))
        end
        register_successors!(store, s, id)
        verbose && id % 100_000 == 0 && println("$id of $n cells : $(Dates.format(now(), "HH:MM:SS"))")
    end
    propagate_bad!(s)
    verbose && println("$n cells, dead: $(count(s.dead))")

    allowed, lost = cell_choices(shape, n_acts, s, (ai, k) -> can_act(ctx, split(cells[k])..., ai), store_acted(store))
    # Where no location of the cell can take an action, it does nothing: it is allowed exactly where waiting is.
    for k in cells
        c, layer = split(k)
        (isodd(allowed[k]) && !isempty(starting(ctx, c, layer))) || continue
        for ai in 1:n_acts
            can_act(ctx, c, layer, ai) || (allowed[k] |= one(eltype(allowed)) << ai)
        end
    end
    return grid, allowed, lost, (; ctx)
end
