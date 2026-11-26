#=
Shield synthesis for hybrid systems with several locations.
A location has its own dynamics and belongs to a layer of the grid, and can be restricted to a region of the continuous state:
  - one location, one layer, edges to itself        : a system without discrete modes (see ModalHybridSystem(::EuclideanHybridSystem))
  - one layer per location                          : discrete modes
  - several locations in a layer, different regions : dynamics that switch across the state, with edges between the regions
A cell of the grid is CartesianIndex(continuous cell..., layer), and its linear index is its id in the FlatShieldState of flat.jl.

One decision period of a cell is computed by modal_reach!:
  - If its flowpipe over [0, p] touches no guard and no unsafe set, the end point of the period is the whole answer. The check uses the adaptive time steps δ⁻..δ⁺. Nothing is assumed about how the state moves between the end points.
  - Otherwise the flowpipe is walked in pieces of δ⁻. A piece touching a guard jumps on its own and continues in the target location for its own remaining time. A piece touching an unsafe set makes the transition bad. States that have not crossed yet continue in the current location.
Safety along the flowpipe comes from the unsafe sets. Leaving the state space is only checked on the sets at the end of the period, because the flowpipe pieces are bloated and a set that starts on the border of the state space would appear to leave it. So a state that leaves and comes back within the period is not caught: what must never happen belongs in the unsafe sets. The exception is a face that a guard borders (such as the ball's floor): beyond it the set at the end of the period holds states that jumped, so instead the set at every step of δ⁻ is checked there (see leaves_through_bordered).

Clock axes: an axis whose row in A is zero and with no input generator (a clock, or anything else with a constant rate) only translates in that location. Two things keep such an axis exact, since otherwise time would drift by a cell per period:
  - Cells are half-open on it. A translated cell that merely touches a neighbouring cell at its border does not count as reaching it.
  - Its interval at the end of the period is tracked exactly through jumps that keep it rigid (same rate on both sides, jump only translates it) or reset it (see clock_after_jump and clock_kept). Any other jump, by a guard or an action, stops the tracking. The successors are filtered by that interval, so the bloat of the flowpipe pieces does not reach this axis.
=#

# Several ModalActions with the same action id in one location, and the same guard, are one choice with several outcomes, which the shield all takes as successors.
struct ModalAction
    action::Action
    target::Int
end

# A location with starts = false is only entered by actions and edges. So a layer can be a decision point that every period starts from, 
# with the outcome of each choice in locations of their own.
struct ModalLocation
    layer::Int
    region::Union{HPolyhedron,Nothing} # part of the continuous state this location covers, nothing = all of it
    A::Matrix{Float64}
    U::Zonotope{Float64,Vector{Float64},Matrix{Float64}}
    edges::Vector{Edge}
    actions::Vector{ModalAction}
    starts::Bool
end
ModalLocation(layer, region, A, U, edges, actions; starts::Bool=true) = ModalLocation(layer, region, A, U, edges, actions, starts)

struct ModalHybridSystem{S<:LazySet,C<:LazySet}
    statespace::S
    unsafe::Vector{C}
    locations::Vector{ModalLocation}
end

# A system without discrete modes, as one location. Its edges must lead back to location 1.
function ModalHybridSystem(system::EuclideanHybridSystem)
    U = Zonotope(Vector{Float64}(LazySets.center(system.input)), Matrix{Float64}(genmat(system.input)))
    location = ModalLocation(1, nothing, Matrix{Float64}(system.flowMatrix), U, system.edges, [ModalAction(act, 1) for act in system.Act])
    return ModalHybridSystem(system.statespace, system.globalConstraints, [location])
end

# wrapper for both discretized generators and cover for a set centered at the origin
mutable struct ShapeData{D}
    G::Matrix{Float64}
    disc::Union{Nothing,CellDiscretization{D}}
    covers::Union{Nothing,Tuple{ConstantZonotopeCover{D},Vector{ConstantZonotopeCover{D}},Vector{ConstantZonotopeCover{D}}}}
end

struct LocationData{D}
    dynamics::Dynamics # shared between locations with the same (A, U)
    shapes::Vector{ShapeData{D}}
    layer::Int
    starts::Bool
    region::Union{Tuple{Vector{Vector{Float64}},Vector{Float64}},Nothing}
    guards::Vector{Tuple{Vector{Vector{Float64}},Vector{Float64}}}
    fast::BitArray{D}
    clock_rate::Vector{Float64} # rate of each axis that only translates in this location, NaN for the others
    bordered::Vector{Bool} # per face of the state space: a guard of the location borders it
    # Per bordered face a ⋅ x <= b: column k + 1 of the matrix is Φ[k+1]' a, and entry k + 1 of the vector b - ρ(a, V[k+1]) + MODAL_TOL
    bordered_dirs::Vector{Matrix{Float64}}
    bordered_limits::Vector{Vector{Float64}}
end

# A set covering part of the end of the period in a location. On axes with lo[i] not NaN, the states are exactly in [lo[i], hi[i]].
struct EndSet
    location::Int
    set::Zonotope{Float64,Vector{Float64},Matrix{Float64}}
    lo::Vector{Float64}
    hi::Vector{Float64}
end

struct ModalContext{D,N,S}
    system::S
    grid::Grid{D}
    steps::Int
    δ⁻::Float64
    δ⁺::Float64
    alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg
    maxOrder::Int
    reduceOrder::Int
    locations::Vector{LocationData{D}}
    unsafe::Vector{Tuple{Vector{Vector{Float64}},Vector{Float64}}}
    statespace::Tuple{Vector{Vector{Float64}},Vector{Float64}}
    cell_generators::Matrix{Float64}
    cell_half::Vector{Float64} # half the width of a cell per axis: g / 2, or on an axis with a single cell half the state space's width (0 for a constant)
    action_ids::Vector{Int}
    may_act::Vector{Vector{BitArray{D}}} # per location, per position in action_ids: some state of the cell is inside the action's guard
    must_act::Vector{Vector{BitArray{D}}} # same, but every state of the cell is inside it
    n_layers::Int
    center::Vector{Float64}
    max_depth::Int
    covers::CoverCache{D}
    layer_starts::Vector{Vector{Int}} # per layer, the locations cells start a period in
    start_mask::Array{UInt64,N} # per (cell..., layer), bit j set if the cell starts in layer_starts[layer][j]
    sweep::Vector{Int} # cells per block and dimension
    exit_bad_for_flow::Bool
    exit_bad_for_actions::Bool
end

const MODAL_TOL = 1e-9

# The locations a cell starts a period in, as a set bit per location of its layer.
struct LocationSet
    ids::Vector{Int}
    mask::UInt64
end
Base.iterate(s::LocationSet, m::UInt64=s.mask) = m == 0 ? nothing : (s.ids[trailing_zeros(m)+1], m & (m - 1))
Base.length(s::LocationSet) = count_ones(s.mask)
Base.isempty(s::LocationSet) = s.mask == 0
Base.eltype(::Type{LocationSet}) = Int

# Exact interval at the end of the period on the clock axes of a location, for X flowing `steps` steps of δ⁻.
function clock_end(rate::Vector{Float64}, X, steps::Int, δ⁻::Float64)
    lo = fill(NaN, length(rate))
    hi = fill(NaN, length(rate))
    for i in eachindex(rate)
        if !isnan(rate[i])
            r = sum(abs, view(genmat(X), i, :))
            lo[i] = LazySets.center(X)[i] - r + rate[i] * steps * δ⁻
            hi[i] = LazySets.center(X)[i] + r + rate[i] * steps * δ⁻
        end
    end
    return lo, hi
end

# End intervals after a jump through `guard`, with between D′ and D′ + W′ steps left. For axis i:
#   the jump resets it to j[i]  : it only depends on the time left
#   the jump translates it, and the rate is the same on both sides : the time of the jump does not matter, so it just translates. A guard on this axis alone also limits it: the states that jumped were inside the guard, and then flowed for the time left.
#   otherwise                   : not tracked any more
# Returns nothing if no state of the half-open cell can have taken the jump, such as a cell whose end only touches the guard at its excluded upper border.
# The upper end of a reset, or of a guard's limit, is the image of the most time left. A state jumping in the first piece
# (at_start) can be inside the guard from the start and have exactly that time left, from anywhere in the cell, so the
# upper end is reached. It is then moved up by closed_margin, which keeps the cell above it. In a later piece the states
# with exactly that time left jump at its start, which the piece before also covers. With a negative rate the upper end
# is the least time left, which a state can have exactly, so it is moved up always.
closed_margin(g::Float64) = 2MODAL_TOL * max(1.0, g)

# Whether some state of X is inside the guard, with X half-open on the tracked clock axes (those where lo is not NaN).
function inside_at_start(X, guard, lo::Vector{Float64})
    l, h = low(X), high(X)
    for i in eachindex(lo)
        isnan(lo[i]) || (h[i] = max(l[i], h[i] - 10MODAL_TOL * max(1.0, abs(h[i]))))
    end
    # A halfspace the box lies beyond decides it exactly. Otherwise isdisjoint decides, which counts touching as meeting.
    any(((a, b),) -> sum(a[i] * (a[i] > 0 ? l[i] : h[i]) for i in eachindex(a)) > b + MODAL_TOL, zip(guard...)) && return false
    return !isdisjoint(Hyperrectangle(low=l, high=h), HPolyhedron([LazySets.HalfSpace(a, b) for (a, b) in zip(guard...)]))
end
function clock_after_jump(lo::Vector{Float64}, hi::Vector{Float64}, guard, M::Matrix{Float64}, j::Vector{Float64}, source::Vector{Float64}, target::Vector{Float64}, D′::Int, W′::Int, δ⁻::Float64, at_start::Bool, g::Float64)
    lo′ = fill(NaN, length(lo))
    hi′ = fill(NaN, length(lo))
    for i in eachindex(lo)
        isnan(target[i]) && continue
        least, most = minmax(target[i] * D′ * δ⁻, target[i] * (D′ + W′) * δ⁻)
        margin = at_start || target[i] < 0 ? closed_margin(g) : 0.0
        if all(iszero, view(M, i, :))
            lo′[i], hi′[i] = j[i] + least, j[i] + most + margin
        elseif !isnan(lo[i]) && source[i] == target[i] && clock_kept(M, i)
            l, h = lo[i], hi[i]
            for (a, b) in zip(guard...)
                if a[i] != 0.0 && count(!iszero, a) == 1
                    a[i] > 0 ? (h = min(h, b / a[i] + most + margin)) : (l = max(l, b / a[i] + least))
                end
            end
            if l > h + MODAL_TOL || (hi[i] > lo[i] + MODAL_TOL && l >= hi[i] - MODAL_TOL)
                return nothing
            end
            lo′[i], hi′[i] = l + j[i], h + j[i]
        end
    end
    return lo′, hi′
end

# Whether a jump with matrix M keeps the upper end of a cell open on axis i: it resets the axis or only translates it.
# Any other jump, such as a flip, can take the closed lower end of the cell to the upper end, so the axis is no longer
# tracked after it.
clock_kept(M::AbstractMatrix, i::Int) = all(iszero, view(M, i, :)) || (M[i, i] == 1.0 && count(!iszero, view(M, i, :)) == 1)


function leaves_statespace(ctx::ModalContext, loc::LocationData, X, D::Int, total::Int, lo::Vector{Float64}, hi::Vector{Float64})
    dyn = loc.dynamics
    for (n, k) in enumerate((D, total))
        n == 2 && k == D && break
        final = minkowski_sum(linear_map(dyn.Φ[k+1], X), dyn.V[k+1])
        for (f, (a, b)) in enumerate(zip(ctx.statespace...))
            loc.bordered[f] && continue
            axis = findfirst(!iszero, a)
            tracked = !isnan(lo[axis])
            value = tracked ? (a[axis] > 0 ? a[axis] * hi[axis] : a[axis] * lo[axis]) : ρ(a, final)
            value > b + (tracked && a[axis] > 0 ? 2closed_margin(ctx.grid.granularity) : MODAL_TOL) || continue # a tracked upper end may be closed_margin up
            beyond = clip_to_polyhedron(final, ([-a], [-b]))
            if beyond === nothing
                tracked && return true # the interval is exact, the set is not
                continue
            end
            inside_a_guard(loc, beyond, a, b) || return true
        end
    end
    return false
end

# True if piece k of the flowpipe, Φ[k+1] Z + V[k+1], reaches beyond a state-space face a ⋅ x <= b that a guard of the
# location borders: ρ(Φ[k+1]' a, Z) > b - ρ(a, V[k+1]), with both sides precomputed per face (see bordered_dirs).
function piece_beyond_bordered(loc::LocationData, k::Int, Z)
    c = LazySets.center(Z)
    G = genmat(Z)
    for (d, limit) in zip(loc.bordered_dirs, loc.bordered_limits)
        s = 0.0
        @inbounds for i in eachindex(c)
            s += d[i, k+1] * c[i]
        end
        @inbounds for j in axes(G, 2)
            t = 0.0
            for i in axes(G, 1)
                t += d[i, k+1] * G[i, j]
            end
            s += abs(t)
        end
        s > limit[k+1] && return true
    end
    return false
end

function leaves_through_bordered(ctx::ModalContext, loc::LocationData, set)
    for (f, (a, b)) in enumerate(zip(ctx.statespace...))
        loc.bordered[f] || continue
        ρ(a, set) > b + MODAL_TOL || continue
        beyond = clip_to_polyhedron(set, ([-a], [-b]))
        beyond !== nothing && !inside_a_guard(loc, beyond, a, b) && return true
    end
    return false
end

# True if `beyond`, the part of a set beyond the state-space face a ⋅ x <= b, lies inside a guard of the location. A guard's
# halfspace -a ⋅ x <= -b (it borders the face) holds there by construction, and is skipped: the clip is not exact on it.
inside_a_guard(loc::LocationData, beyond, a, b) =
    any(guard -> all(i -> (isapprox(guard[1][i], -a; atol=MODAL_TOL) && isapprox(guard[2][i], -b; atol=MODAL_TOL)) ||
        ρ(guard[1][i], beyond) <= guard[2][i], eachindex(guard[1])), loc.guards)

# Appends to `out` what X reaches, covering it at the end of the period. X starts with between D and D + W steps of δ⁻ left, depending on the state, so W > 0 after a jump. lo and hi are its exact intervals at the end of the period on the clock axes.
# Returns false if the transition is bad: it touches an unsafe set, or leaves the state space when exit_bad is true.

function modal_reach!(out::Vector{EndSet}, ctx::ModalContext, li::Int, X, D::Int, W::Int, exit_bad::Bool, lo::Vector{Float64}, hi::Vector{Float64}, depth::Int=0)
    depth > ctx.max_depth && return false # jumping again and again
    total = D + W
    if total == 0
        push!(out, EndSet(li, X, lo, hi))
        return true
    end
    loc = ctx.locations[li]
    dyn = loc.dynamics
    δ⁺ = min(ctx.δ⁺, ctx.δ⁻ * 2^floor(Int, log2(total))) # no step longer than the time left
    # A whole cell, or a cell just after an action, flowing a whole period uses the precomputed discretization of its shape.
    cd = depth == 0 && total == ctx.steps ? shape_discretization!(ctx, li, genmat(X)) : nothing
    disc = cd === nothing ? newReACTDiscretizePlus(X, ctx.δ⁻, δ⁺, dyn.A, dyn.P2A_abs, dyn.step_Φ, dyn.step_input, ctx.alg, ctx.maxOrder, ctx.reduceOrder) :
           discretize_cell!(cd, LazySets.center(X))
    k0 = first_touch(dyn, disc, Iterators.flatten((loc.guards, ctx.unsafe)), total, ctx.δ⁻, δ⁺)
    check_bordered = exit_bad && !isempty(loc.bordered_dirs)

    # Pieces from the first touch on, and the pieces that end the period when the time left varies.
    for k in (check_bordered ? 0 : min(k0, W==0 ? total : D)):(total-1)
        if check_bordered && piece_beyond_bordered(loc, k, disc[ctx.δ⁻])
            leaves_through_bordered(ctx, loc, minkowski_sum(linear_map(dyn.Φ[k+2], X), dyn.V[k+2])) && return false
        end
        (k >= k0 || (W > 0 && k >= D)) || continue
        piece = minkowski_sum(linear_map(dyn.Φ[k+1], disc[ctx.δ⁻]), dyn.V[k+1])
        if k >= k0
            any(set -> some_inside(piece, set...), ctx.unsafe) && return false
            all_jumped = false
            for (e, guard) in zip(ctx.system.locations[li].edges, loc.guards)
                some_inside(piece, guard...) || continue
                clipped = clip_to_polyhedron(piece, guard)
                clipped === nothing && continue
                # A state in this piece has spent between k and k + 1 steps here, and had between D and D + W left at the start.
                D′ = max(D - k - 1, 0)
                W′ = total - k - D′
                # After a jump (depth > 0) X is not a cell, so any state of it may be inside the guard at the start.
                at_start = k == 0 && (depth > 0 || inside_at_start(X, guard, lo))
                clocks = clock_after_jump(lo, hi, guard, e.jumpMatrix, e.jumpVector, loc.clock_rate, ctx.locations[e.targetLoc].clock_rate, D′, W′, ctx.δ⁻, at_start, ctx.grid.granularity)
                clocks === nothing && continue
                jumped = jump_set(clipped, e.jumpMatrix, e.jumpVector, ctx.maxOrder, ctx.reduceOrder)
                modal_reach!(out, ctx, e.targetLoc, jumped, D′, W′, exit_bad, clocks..., depth + 1) || return false
                all_jumped = all_jumped || all_inside(piece, guard...)
            end
            all_jumped && return true
        end
        if W > 0 && k >= D
            push!(out, EndSet(li, piece, lo, hi))
        end
    end
    exit_bad && leaves_statespace(ctx, loc, X, D, total, lo, hi) && return false
    if W == 0
        push!(out, EndSet(li, minkowski_sum(linear_map(dyn.Φ[D+1], X), dyn.V[D+1]), lo, hi))
    end
    return true
end

# True if the interior of the box [lo, hi] overlaps the region.
function box_overlaps_region(region, lo, hi)
    region === nothing && return true
    for (a, b) in zip(region...)
        low_value = sum(a[i] * (a[i] > 0 ? lo[i] : hi[i]) for i in eachindex(a))
        low_value >= b - MODAL_TOL && return false
    end
    return true
end

function ModalContext(system::ModalHybridSystem, p, granularity::Float64, δ⁻::Float64, δ⁺::Float64, alg, maxOrder::Int, reduceOrder::Int, sweep_factor, max_depth::Int;
    exit_bad_for_flow::Bool=true, exit_bad_for_actions::Bool=true)
    grid = Grid(system.statespace, granularity)
    return _ModalContext(grid, Val(ndims(grid.deadCells) + 1), system, Float64(p), δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor, max_depth, exit_bad_for_flow, exit_bad_for_actions)
end


function _ModalContext(grid::Grid{D}, ::Val{N}, system, p::Float64, δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor, max_depth, exit_bad_for_flow, exit_bad_for_actions) where {D,N}
    steps = round(Int, p / δ⁻)
    isapprox(steps * δ⁻, p; atol=MODAL_TOL) || throw(ArgumentError("The period $p must be a multiple of δ⁻ = $δ⁻"))
    isinteger(log2(δ⁺ / δ⁻)) || throw(ArgumentError("δ⁺ / δ⁻ must be a power of two"))

    faces = hyperrectangle_to_HPolyhedron(system.statespace)
    statespace = getHalfSpaceProjections(faces)
    exits = [([-Vector{Float64}(c.a)], [-c.b - MODAL_TOL]) for c in constraints_list(faces)] # each face, as the set strictly beyond it
    unsafe = [getHalfSpaceProjections(H) for H in system.unsafe]

    g = grid.granularity
    # On an axis with a single cell, the state space is at most a cell wide, the cell is the state space: a point for a constant, and as wide as the state space otherwise.
    cell_half = [n == 1 ? (h - l) / 2 : g / 2 for (n, l, h) in zip(grid.numCells, LazySets.low(system.statespace), LazySets.high(system.statespace))]
    cell_generators = Matrix{Float64}(genmat(remove_zero_generators(Zonotope(zeros(D), diagm(cell_half)))))
    action_ids = sort!(unique(a.action.id for loc in system.locations for a in loc.actions))

    sweep = sweep_factor isa Integer ? [min(sweep_factor, n) for n in grid.numCells] : [min(s, n) for (s, n) in zip(sweep_factor, grid.numCells)]
    G_big = Matrix{Float64}(genmat(remove_zero_generators(Zonotope(zeros(D), diagm(map((h, n, k) -> n == 1 ? h : g * (k / 2), cell_half, grid.numCells, sweep))))))
    cell_δ⁺ = min(δ⁺, δ⁻ * 2^floor(Int, log2(steps))) # what modal_reach! uses for a whole period
    dyn_key(l) = (l.A, LazySets.center(l.U), genmat(l.U))
    shared_dyn = Dict{Tuple{Matrix{Float64},Vector{Float64},Matrix{Float64}},Dynamics}()
    dynamics_of(l) = get!(() -> Dynamics(l.A, l.U, steps, δ⁻, δ⁺, alg, maxOrder, reduceOrder), shared_dyn, dyn_key(l))
    exits_and(guards) = vcat(guards, unsafe, exits)
    block_ranges(bidx) = ntuple(i -> (1+(bidx[i]-1)*sweep[i]):(bidx[i]*sweep[i]), D)
    block_center(bidx) = grid.lower .+ g .* (sweep ./ 2) .+ (index_vector(bidx) .- 1) .* g .* sweep
    blocks = CartesianIndices(ntuple(i -> 1:fld(grid.numCells[i], sweep[i]), D))

    locations = LocationData{D}[]
    for loc in system.locations
        dyn = dynamics_of(loc)
        li = length(locations) + 1
        guards = [getHalfSpaceProjections(e.guard) for e in loc.edges]
        targets = exit_bad_for_flow ? exits_and(guards) : vcat(guards, unsafe)
        region = loc.region === nothing ? nothing : getHalfSpaceProjections(loc.region)

        entered_at_start = loc.starts || any(a.target == length(locations) + 1 for l in system.locations for a in l.actions)
        bordered = [bordered_by(guards, a, b; atol=MODAL_TOL) for (a, b) in zip(statespace...)]
        bordered_dirs = [reduce(hcat, [transpose(dyn.Φ[k+1]) * statespace[1][f] for k in 0:steps]) for f in eachindex(bordered) if bordered[f]]
        bordered_limits = [[statespace[2][f] - ρ(statespace[1][f], dyn.V[k+1]) + MODAL_TOL for k in 0:steps] for f in eachindex(bordered) if bordered[f]]
        fast = falses(cellshape(grid))
        big = entered_at_start && !isempty(blocks) ? period_discretization(dyn, loc.U, G_big, δ⁻, cell_δ⁺) : nothing
        for bidx in blocks
            entered_at_start || break
            center = block_center(bidx)
            box_overlaps_region(region, center .- g .* (sweep ./ 2), center .+ g .* (sweep ./ 2)) || continue
            if first_touch(dyn, discretize_cell!(big, center), targets, steps, δ⁻, cell_δ⁺) == steps
                fast[block_ranges(bidx)...] .= true
            end
        end
        clock_rate = [all(iszero, view(loc.A, i, :)) && all(iszero, view(genmat(loc.U), i, :)) ? LazySets.center(loc.U)[i] : NaN for i in 1:D]
        push!(locations, LocationData{D}(dyn, ShapeData{D}[], loc.layer, loc.starts, region, guards, fast, clock_rate, bordered, bordered_dirs, bordered_limits))
    end

    # A location that does not offer an action has neither mask set for it, so the action does nothing there.
    may_act = [[falses(cellshape(grid)) for _ in action_ids] for _ in system.locations]
    must_act = [[falses(cellshape(grid)) for _ in action_ids] for _ in system.locations]
    for (li, loc) in enumerate(system.locations), a in loc.actions
        same = [b.action.guard for b in loc.actions if b.action.id == a.action.id]
        all(isequal(a.action.guard), same) || throw(ArgumentError("location $li has actions with id $(a.action.id) and different guards; the outcomes of one action share its guard"))
        ai = findfirst(==(a.action.id), action_ids)
        may_act[li][ai], must_act[li][ai] = guard_masks(grid, a.action.guard)
    end

    n_layers = maximum(loc.layer for loc in system.locations)
    layer_starts = [[li for (li, loc) in enumerate(locations) if loc.layer == layer && loc.starts] for layer in 1:n_layers]
    start_mask = start_masks(grid, locations, layer_starts)

    return ModalContext{D,N,typeof(system)}(system, grid, steps, δ⁻, δ⁺, alg, maxOrder, reduceOrder, locations, unsafe, statespace, cell_generators, cell_half, action_ids, may_act, must_act, n_layers, zeros(D), max_depth, CoverCache{D}(), layer_starts, start_mask, sweep, exit_bad_for_flow, exit_bad_for_actions)
end

# Per (cell..., layer), bit j set if the cell overlaps the region of location layer_starts[layer][j].
function start_masks(grid::Grid{D}, locations::Vector{LocationData{D}}, layer_starts::Vector{Vector{Int}}) where D
    n_layers = length(layer_starts)
    g = grid.granularity
    start_mask = zeros(UInt64, cellshape(grid)..., n_layers)
    lo, hi = zeros(D), zeros(D)
    for layer in 1:n_layers
        ids = layer_starts[layer]
        length(ids) <= 64 || throw(ArgumentError("layer $layer has $(length(ids)) locations that cells start in; at most 64 are supported"))
        if all(li -> locations[li].region === nothing, ids)
            selectdim(start_mask, D + 1, layer) .= isempty(ids) ? UInt64(0) : typemax(UInt64) >> (64 - length(ids))
            continue
        end
        for c in CartesianIndices(cellshape(grid))
            for i in 1:D
                lo[i] = grid.lower[i] + (c[i] - 1) * g
                hi[i] = lo[i] + g
            end
            m = UInt64(0)
            for (j, li) in enumerate(ids)
                box_overlaps_region(locations[li].region, lo, hi) && (m |= UInt64(1) << (j - 1))
            end
            start_mask[c, layer] = m
        end
    end

    return start_mask
end

cell_locations(ctx::ModalContext{D}, layer::Int, c::CartesianIndex{D}) where D = LocationSet(ctx.layer_starts[layer], ctx.start_mask[c, layer])

function can_act(ctx::ModalContext{D}, layer::Int, c::CartesianIndex{D}, ai::Int) where D
    for li in cell_locations(ctx, layer, c)
        ctx.may_act[li][ai][c] && return true
    end
    return false
end

function append_touching!(res::Vector{CartesianIndex{N}}, ctx::ModalContext{D,N}, out::Vector{EndSet}) where {D,N}
    grid = ctx.grid
    for e in out
        # touching_cells finds no cells for a set with zero width exactly on the lower edge of the grid on some axis, because the upper index is ceil(0) = 0. A reset to that edge (like l := 0) gives such a set, so it is moved just inside the grid.
        c = copy(LazySets.center(e.set))
        G = genmat(e.set)
        for i in 1:D
            if abs(c[i] - grid.lower[i]) <= MODAL_TOL && all(iszero, view(G, i, :))
                c[i] = grid.lower[i] + 10 * MODAL_TOL
            end
        end
        for cell in touching_cells(grid, Zonotope(c, G), ctx.covers, clock_cell_bounds(grid, e.lo, e.hi)...)
            push!(res, CartesianIndex(cell, ctx.locations[e.location].layer))
        end
    end
    return res
end

function clock_cell_bounds(grid::Grid{D}, lo::Vector{Float64}, hi::Vector{Float64}) where D
    keep_lo = ntuple(i -> isnan(lo[i]) ? 1 : floor(Int, (lo[i] - grid.lower[i]) / grid.granularity + MODAL_TOL) + 1, D)
    keep_hi = ntuple(i -> isnan(lo[i]) ? grid.numCells[i] : max(ceil(Int, (hi[i] - grid.lower[i]) / grid.granularity - MODAL_TOL), keep_lo[i]), D)
    return keep_lo, keep_hi
end

function lands_in_fast_cells(ctx::ModalContext{D}, cov::ConstantZonotopeCover{D}, loc::LocationData{D}, c::AbstractVector{Float64}) where D
    grid = ctx.grid
    for i in 1:D
        (c[i] - cov.radius[i] < grid.lower[i] || c[i] + cov.radius[i] > grid.upper[i]) && return false
    end
    cells = cells_in_cover(grid, cov, c)
    return !isempty(cells) && all(cell -> loc.fast[cell], cells)
end

function append_cover!(res::Vector{CartesianIndex{N}}, ctx::ModalContext{D,N}, cov::ConstantZonotopeCover{D}, layer::Int, lo::Vector{Float64}, hi::Vector{Float64}) where {D,N}
    for cell in cells_in_cover(ctx.grid, cov, ctx.center, clock_cell_bounds(ctx.grid, lo, hi)...)
        push!(res, CartesianIndex(cell, layer))
    end
    return res
end


function cell_shift(M::Matrix{Float64}, j::Vector{Float64}, grid::Grid{D}) where D
    for k in 1:D, i in 1:D
        M[i, k] == (i == k ? 1.0 : 0.0) || return nothing
    end
    shift = ntuple(i -> j[i] / grid.granularity, D)
    all(i -> abs(shift[i] - round(shift[i])) <= MODAL_TOL, 1:D) || return nothing
    return CartesianIndex(ntuple(i -> round(Int, shift[i]), D))
end

function sorted_unique!(v::Vector)
    sort!(v)
    j = 0
    for i in eachindex(v)
        (j == 0 || v[i] != v[j]) && (j += 1; v[j]=v[i])
    end
    return resize!(v, j)
end


function noop_choices!(allowed::Array{T,N}, ctx::ModalContext{D,N}) where {T,D,N}
    for k in CartesianIndices(allowed)
        isodd(allowed[k]) || continue
        c = CartesianIndex(ntuple(i -> k[i], D))
        locs = cell_locations(ctx, k[N], c)
        isempty(locs) && continue
        for ai in eachindex(ctx.action_ids)
            all(li -> !ctx.may_act[li][ai][c], locs) && (allowed[k] |= T(1) << ai)
        end
    end
    return allowed
end

box_cells(b::CartesianIndex{D}, size::NTuple{D,Int}, numCells) where D =
    CartesianIndices(ntuple(i -> ((b[i]-1)*size[i]+1):min(b[i]*size[i], numCells[i]), D))

function shape_data!(ctx::ModalContext{D}, li::Int, G::Matrix{Float64}) where D
    loc = ctx.locations[li]
    for sd in loc.shapes # a few shapes per location: a cell, jumped cells, the blocks, the leaves
        (sd.G === G || sd.G == G) && return sd
    end
    push!(loc.shapes, ShapeData{D}(copy(G), nothing, nothing))
    return loc.shapes[end]
end

# Discretization over a whole period of a box of shape G under the location's dynamics, from which modal_reach! gets the
# discretization of any box of that shape (see CellDiscretization).
period_discretization(dyn::Dynamics, U, G::Matrix{Float64}, δ⁻::Float64, δ⁺::Float64) =
    CellDiscretization(Zonotope(zeros(size(G, 1)), G), (; δ⁻, δ⁺, dyn.A, dyn.P2A_abs, phiDict=dyn.step_Φ, inputDict=dyn.step_input, U), Val(size(G, 1)))

function shape_discretization!(ctx::ModalContext, li::Int, G::Matrix{Float64})
    sd = shape_data!(ctx, li, G)
    if sd.disc === nothing
        sd.disc = period_discretization(ctx.locations[li].dynamics, ctx.system.locations[li].U, sd.G, ctx.δ⁻, min(ctx.δ⁺, ctx.δ⁻ * 2^floor(Int, log2(ctx.steps))))
    end
    return sd.disc
end

# The covers of a box of shape G in location li: of its end set when it is fast, and per outgoing action (as in the
# location's actions) of the jumped box at the start of the period, to find the cells it lands in, and at its end, when
# those cells are all fast in the target. Each is bloated by half a cell for the test of cell centres, and its candidate
# cells are bounded by the unbloated set, as in touching_cells: a cell that only shares a face with the set is
# not a successor.
function shape_covers!(ctx::ModalContext{D}, li::Int, G::Matrix{Float64}) where D
    sd = shape_data!(ctx, li, G)
    sd.covers === nothing || return sd.covers
    bloat = Matrix((ctx.grid.granularity / 2) * I, D, D)
    cover(H) = ConstantZonotopeCover(hcat(H, bloat), ntuple(i -> sum(abs, view(H, i, :)), D), Val(D))
    dyn = ctx.locations[li].dynamics
    starts, ends = ConstantZonotopeCover{D}[], ConstantZonotopeCover{D}[]
    for a in ctx.system.locations[li].actions
        tdyn = ctx.locations[a.target].dynamics
        G_J = Matrix{Float64}(genmat(jump_set(Zonotope(zeros(D), G), Matrix{Float64}(a.action.jumpMatrix), zeros(D), ctx.maxOrder, ctx.reduceOrder)))
        push!(starts, cover(G_J))
        push!(ends, cover(hcat(genmat(tdyn.V[ctx.steps+1]), tdyn.Φ[ctx.steps+1] * G_J)))
    end
    sd.covers = (cover(hcat(genmat(dyn.V[ctx.steps+1]), dyn.Φ[ctx.steps+1] * G)), starts, ends)
    return sd.covers
end

struct BoxSuccessors{N}
    flow::Union{Nothing,Vector{CartesianIndex{N}}}
    enabled::Vector{Bool}
    acts::Vector{Union{Nothing,Vector{CartesianIndex{N}}}}
end

# Sets ctx.center to where the centre of a box that touches nothing ends the period in location li: the box has len cells
# per axis from cell `first`.
function fast_end_center!(ctx::ModalContext{D}, li::Int, first::CartesianIndex{D}, len::NTuple{D,Int}) where D
    dyn = ctx.locations[li].dynamics
    Φ = dyn.Φ[ctx.steps+1]
    v = LazySets.center(dyn.V[ctx.steps+1])
    grid = ctx.grid
    g = grid.granularity
    for i in 1:D
        s = v[i]
        for j in 1:D
            s += Φ[i, j] * ((grid.lower[j] + len[j] * g / 2) + (first[j] - 1) * g)
        end
        ctx.center[i] = s
    end
    return ctx.center
end

# What the box of cells `cells` of a layer (a single cell, a leaf or a block) reaches in one period, per choice, flowed as
# one zonotope. It contains what each of its cells reaches. Shortcuts: a box whose cells are all fast in a location uses
# that location's cover of its shape; so does a jumped box that lands in cells that are all fast in the target, or that a
# jump moves by whole cells onto cells that are all fast there.
# flow is nothing if waiting is bad for the box. Per action position, enabled says whether some location of a cell can
# take the action there, and acts its successors, nothing if it is bad. Trying an action outside its guard, or in a
# location that does not offer it, does nothing: the state just flows. So where some state of the box does not take it,
# its successors include the waiting ones; if waiting is bad, those are flow_if_bad() (a function, nothing to give
# none), or the action is bad too.
function box_successors(ctx::ModalContext{D,N}, layer::Int, cells::CartesianIndices{D}; flow_if_bad=nothing) where {D,N}
    grid = ctx.grid
    g = grid.granularity
    len = map(length, cells.indices)
    G = all(==(1), len) ? ctx.cell_generators :
        Matrix{Float64}(genmat(remove_zero_generators(Zonotope(zeros(D), diagm(map((h, n, l) -> n == 1 ? h : g * (l / 2), ctx.cell_half, grid.numCells, len))))))
    Z = Zonotope(grid.lower .+ ((Tuple(first(cells)) .- 1) .+ len ./ 2) .* g, G)
    mask = UInt64(0)
    for c in cells
        mask |= ctx.start_mask[c, layer]
    end
    locs = LocationSet(ctx.layer_starts[layer], mask)
    out = EndSet[]
    flow = CartesianIndex{N}[]
    flow_bad = false
    for li in locs
        loc = ctx.locations[li]
        lo, hi = clock_end(loc.clock_rate, Z, ctx.steps, ctx.δ⁻)
        if all(c -> loc.fast[c], cells) # every cell touches nothing over the period, so neither does the box
            fast_end_center!(ctx, li, first(cells), len)
            append_cover!(flow, ctx, shape_covers!(ctx, li, G)[1], loc.layer, lo, hi)
        else
            empty!(out)
            modal_reach!(out, ctx, li, Z, ctx.steps, 0, ctx.exit_bad_for_flow, lo, hi) ? append_touching!(flow, ctx, out) : (flow_bad = true)
        end
    end
    n_acts = length(ctx.action_ids)
    enabled = zeros(Bool, n_acts)
    acts = Vector{Union{Nothing,Vector{CartesianIndex{N}}}}(nothing, n_acts) # a list once enabled; nothing for bad
    also_flows = zeros(Bool, n_acts)
    for li in locs, ai in 1:n_acts
        any(c -> ctx.may_act[li][ai][c], cells) || (also_flows[ai]=true; continue)
        enabled[ai] || (enabled[ai]=true; acts[ai]=CartesianIndex{N}[])
        all(c -> ctx.must_act[li][ai][c], cells) || (also_flows[ai] = true)
        acts[ai] === nothing && continue
        loc = ctx.locations[li]
        for (k, a) in enumerate(ctx.system.locations[li].actions)
            a.action.id == ctx.action_ids[ai] || continue
            jumped = jump_set(Z, a.action.jumpMatrix, a.action.jumpVector, ctx.maxOrder, ctx.reduceOrder)
            target = ctx.locations[a.target]
            lo, hi = clock_end(target.clock_rate, jumped, ctx.steps, ctx.δ⁻)
            for i in eachindex(lo)
                clock_kept(a.action.jumpMatrix, i) || (lo[i] = hi[i] = NaN)
            end
            covers = shape_covers!(ctx, li, G)
            # Fast cells may leave the state space when leaving is allowed for the flow: then only if it is for actions too.
            fast_ok = ctx.exit_bad_for_flow || !ctx.exit_bad_for_actions
            if fast_ok && lands_in_fast_cells(ctx, covers[2][k], target, LazySets.center(jumped))
                mul!(ctx.center, target.dynamics.Φ[ctx.steps+1], LazySets.center(jumped))
                ctx.center .+= LazySets.center(target.dynamics.V[ctx.steps+1])
                append_cover!(acts[ai], ctx, covers[3][k], target.layer, lo, hi)
                continue
            end
            shift = fast_ok ? cell_shift(a.action.jumpMatrix, a.action.jumpVector, grid) : nothing
            if shift !== nothing && all(c -> checkbounds(Bool, target.fast, c + shift) && target.fast[c+shift], cells)
                fast_end_center!(ctx, a.target, first(cells) + shift, len)
                append_cover!(acts[ai], ctx, shape_covers!(ctx, a.target, G)[1], target.layer, lo, hi)
            else
                empty!(out)
                modal_reach!(out, ctx, a.target, jumped, ctx.steps, 0, ctx.exit_bad_for_actions, lo, hi) ? append_touching!(acts[ai], ctx, out) : (acts[ai]=nothing; break)
            end
        end
    end
    flow = flow_bad ? nothing : sorted_unique!(flow)
    waiting = flow !== nothing ? flow : flow_if_bad === nothing ? nothing : flow_if_bad()
    for ai in 1:n_acts
        (enabled[ai] && acts[ai] !== nothing) || continue
        if also_flows[ai]
            waiting === nothing ? (acts[ai] = nothing) : append!(acts[ai], waiting)
        end
        acts[ai] === nothing || sorted_unique!(acts[ai])
    end
    return BoxSuccessors{N}(flow, collect(enabled), acts)
end

to_ids(fids, cs) = sorted_unique!(Int32[fids[v] for v in cs])

# A list of a box restricted to the list of an enclosing box (outer): outer's is a sound over-approximation of all the box
# reaches, so what the box's own over-approximation adds outside it is not reachable. Where the box's is bad and outer's
# is not, outer's stands in.
within(own, outer) = outer === nothing ? own : own === nothing ? copy(outer) : filter!(v -> insorted(v, outer), own)

# The successors of the box of cells `cells` in a layer (box_successors) as ids, per choice, within those of an
# enclosing box: its flow, acts and enabled (per action position). Returns flow and acts, nothing where a choice is bad or
# (for an action) not enabled. So they lie in the enclosing box's by construction, and are bad only where its are.
function successors_within(ctx::ModalContext, layer::Int, cells::CartesianIndices, fids, outer_flow, outer_acts, outer_enabled)
    all_cells = CartesianIndices(fids)
    waiting() = outer_flow === nothing ? nothing : [all_cells[v] for v in outer_flow]
    r = box_successors(ctx, layer, cells; flow_if_bad=waiting)
    flow = within(r.flow === nothing ? nothing : to_ids(fids, r.flow), outer_flow)
    acts = Union{Nothing,Vector{Int32}}[r.enabled[ai] ? within(r.acts[ai] === nothing ? nothing : to_ids(fids, r.acts[ai]), outer_enabled[ai] ? outer_acts[ai] : nothing) : nothing
                                        for ai in eachindex(r.enabled)]
    return flow, acts
end

# Blocks of `size` cells per axis in every layer, each flowed forward once as one box (box_successors). Block ids are the
# linear index of (block index..., layer). Per block: flow and acts (per action position) as in successors_within,
# enabled per action, and bad if some choice is bad. Online also the blocks that reach each cell id, by any choice.
struct ForwardBlocks{D,N}
    size::NTuple{D,Int}
    shape::NTuple{N,Int}
    flow::Vector{Union{Nothing,Vector{Int32}}}
    acts::Vector{Vector{Union{Nothing,Vector{Int32}}}}
    enabled::Vector{Vector{Bool}}
    bad::BitVector
    rev_off::Vector{Int}
    rev_pool::Vector{Int32}
end

function ForwardBlocks(ctx::ModalContext{D,N}, size::NTuple{D,Int}, fids, online::Bool) where {D,N}
    nc = Tuple(ctx.grid.numCells)
    shape = (ntuple(i -> cld(nc[i], size[i]), D)..., ctx.n_layers)
    nb = prod(shape)
    flow = Vector{Union{Nothing,Vector{Int32}}}(undef, nb)
    acts = Vector{Vector{Union{Nothing,Vector{Int32}}}}(undef, nb)
    enabled = Vector{Vector{Bool}}(undef, nb)
    bad = falses(nb)
    succ_off, succ_pool = ones(Int, nb + 1), Int32[]
    for (bid, k) in enumerate(CartesianIndices(shape))
        r = box_successors(ctx, k[N], box_cells(CartesianIndex(ntuple(i -> k[i], D)), size, nc))
        flow[bid] = r.flow === nothing ? nothing : to_ids(fids, r.flow)
        acts[bid] = Union{Nothing,Vector{Int32}}[r.enabled[ai] && r.acts[ai] !== nothing ? to_ids(fids, r.acts[ai]) : nothing for ai in eachindex(r.acts)]
        enabled[bid] = r.enabled
        bad[bid] = r.flow === nothing || any(ai -> r.enabled[ai] && r.acts[ai] === nothing, eachindex(r.enabled))
        if online
            reached = flow[bid] === nothing ? Int32[] : copy(flow[bid])
            foreach(a -> a === nothing || append!(reached, a), acts[bid])
            append!(succ_pool, sorted_unique!(reached))
            succ_off[bid+1] = length(succ_pool) + 1
        end
    end
    rev_off, rev_pool = online ? reverse_csr(succ_off, succ_pool, length(fids)) : (Int[], Int32[])
    return ForwardBlocks{D,N}(size, shape, flow, acts, enabled, bad, rev_off, rev_pool)
end

block_id(B::ForwardBlocks{D}, c::CartesianIndex{D}, layer::Int) where D = LinearIndices(B.shape)[CartesianIndex(ntuple(i -> cld(c[i], B.size[i]), D)..., layer)]
function block_cells(B::ForwardBlocks{D,N}, bid, numCells) where {D,N}
    k = CartesianIndices(B.shape)[bid]
    return box_cells(CartesianIndex(ntuple(i -> k[i], D)), B.size, numCells), k[N]
end
blocks_reaching(B::ForwardBlocks, v) = view(B.rev_pool, B.rev_off[v]:(B.rev_off[v+1]-1))

# The inverse of the lists succ_off/succ_pool (list b holds succ_pool[succ_off[b]:succ_off[b+1]-1]) over the ids 1:n:
# for id v, the lists holding it are rev_pool[rev_off[v]:rev_off[v+1]-1], in increasing order.
function reverse_csr(succ_off::Vector{Int}, succ_pool::Vector{Int32}, n::Int)
    rev_off = zeros(Int, n + 1)
    for v in succ_pool
        rev_off[v+1] += 1
    end
    rev_off[1] = 1
    for v in 1:n
        rev_off[v+1] += rev_off[v]
    end
    rev_pool = Vector{Int32}(undef, length(succ_pool))
    pos = rev_off[1:n]
    for b in 1:(length(succ_off)-1), i in succ_off[b]:(succ_off[b+1]-1)
        v = succ_pool[i]
        rev_pool[pos[v]] = b
        pos[v] += 1
    end
    return rev_off, rev_pool
end

# Shield for a ModalHybridSystem. Computes every cell of every layer, or with online = true only the cells whose choices
# can differ from everything allowed (see ForwardBlocks); both give the same result. Then pushes the bad events backwards.
# Leaving the state space makes waiting bad, and an action bad, unless exit_bad_for_flow or exit_bad_for_actions is false:
# then only the successors inside the grid count, and a choice whose successors all lie outside it is safe.
# Cells are tested for being fast in blocks of sweep_factor cells per dimension (an integer, or one per dimension); online
# also computes these blocks forward. It only changes the time: a fast cell has the same successors as when flowed on its
# own. On cruise control 20 finds no fast cells, 5 and 2 are the fastest, and 1 is slower again because every cell is then
# checked on its own. Every block is also computed forward as one set, and each cell keeps only the successors its block
# reaches too, per choice (see successors_within); online uses the blocks to find the cells to compute.
# Returns grid, allowed, lost, info like every shield (see choice_bits in flat.jl). grid is the continuous
# grid; allowed and lost are indexed by (continuous cell..., layer), and a cell is dead in a layer where allowed is 0.
# Trying an action outside its guard or in a location that does not offer it does nothing, so there it is allowed exactly when waiting is (see box_successors).
# A cell is enabled for an action when some location covering it can take the action there. info.ctx is the context.
function ReACTed_reachable_cell_modal(system::ModalHybridSystem, p, granularity, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5;
    online::Bool=false, act_needs_safe_flow::Bool=true, exit_bad_for_flow::Bool=true, exit_bad_for_actions::Bool=true,
    sweep_factor=5, max_depth::Int=16, verbose::Bool=true)
    timer = PhaseTimer()
    ctx = ModalContext(system, p, Float64(granularity), δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor, max_depth; exit_bad_for_flow, exit_bad_for_actions)
    return _ReACTed_reachable_cell_modal(ctx, online, act_needs_safe_flow, verbose, timer)
end

function _ReACTed_reachable_cell_modal(ctx::ModalContext{D,N}, online::Bool, act_needs_safe_flow::Bool, verbose::Bool, timer::PhaseTimer) where {D,N}
    grid = ctx.grid
    shape = (cellshape(grid)..., ctx.n_layers)
    verbose && println("$(prod(shape)) cells, shape $shape")
    n_acts = length(ctx.action_ids)
    n = prod(shape)
    n < typemax(Int32) || throw(ArgumentError("$n cells do not fit in Int32 ids"))
    cells = CartesianIndices(shape)
    ids = LinearIndices(shape)
    store = SuccessorStore(n, [falses(n) for _ in 1:n_acts])
    store.flow_empty_ok = !ctx.exit_bad_for_flow
    store.act_empty_ok = !ctx.exit_bad_for_actions
    s = reset!(FlatShieldState(n_acts, act_needs_safe_flow), n)
    buf = Int32[]
    queued = falses(shape)
    frontier = CartesianIndex{N}[]
    next_frontier = CartesianIndex{N}[]

    for H in ctx.system.unsafe, c in touching_cells(grid, H), layer in 1:ctx.n_layers
        k = CartesianIndex(c, layer)
        seed_constraint!(s, Int32(ids[k]))
        store.flow_none[ids[k]] = true
        queued[k] = true
    end
    lap!(timer, "Setup")

    blocks = ForwardBlocks(ctx, Tuple(ctx.sweep), ids, online)
    lap!(timer, "Blocks")

    # Online: compute the cells of bad blocks, and whenever a cell becomes bad, those of the blocks that reach it.
    block_queued = falses(length(blocks.bad))
    function queue_block!(out, bid)
        block_queued[bid] && return
        block_queued[bid] = true
        bcells, layer = block_cells(blocks, bid, grid.numCells)
        for c in bcells
            k = CartesianIndex(c, layer)
            queued[k] || (queued[k]=true; push!(out, k))
        end
    end
    function queue_blocks_reaching!(out, ks)
        for k in ks, bid in blocks_reaching(blocks, k)
            queue_block!(out, bid)
        end
    end
    if online
        for bid in findall(blocks.bad)
            queue_block!(next_frontier, bid)
        end
        queue_blocks_reaching!(next_frontier, findall(s.dead))
    else
        for k in CartesianIndices(shape)
            queued[k] || push!(next_frontier, k)
        end
    end
    lap!(timer, "Seeding")

    count_cells = 0
    rounds = 0
    while !isempty(next_frontier)
        frontier, next_frontier = next_frontier, empty!(frontier)
        rounds += 1
        for k in frontier
            c = CartesianIndex(ntuple(i -> k[i], D))
            layer = k[N]
            id = Int32(ids[k])
            (s.computed[id] || isempty(cell_locations(ctx, layer, c))) && continue
            bid = block_id(blocks, c, layer)
            flow, acts = successors_within(ctx, layer, CartesianIndices(ntuple(i -> c[i]:c[i], D)), ids, blocks.flow[bid], blocks.acts[bid], blocks.enabled[bid])
            store.flow_none[id] = flow === nothing
            set_successors!(store, 0, id, flow === nothing ? empty!(buf) : flow)
            for ai in 1:n_acts
                store.act_enabled[ai][id] = acts[ai] !== nothing
                set_successors!(store, ai, id, acts[ai] === nothing ? empty!(buf) : acts[ai])
            end
            register_successors!(store, s, id)
            count_cells += 1
            verbose && count_cells % 10000 == 0 && println("$count_cells cells computed : $(Dates.format(now(), "HH:MM:SS"))")
        end
        lap!(timer, "Reachability")
        propagate_bad!(s)
        online && queue_blocks_reaching!(next_frontier, unique!(s.changed))
        empty!(s.changed)
        lap!(timer, "Shield")
    end
    verbose && report(timer, "Rounds: $rounds, cells computed: $count_cells, dead: $(count(s.dead))")
    enabled(ai, k) = can_act(ctx, cells[k][N], CartesianIndex(ntuple(i -> cells[k][i], D)), ai)
    allowed, lost = cell_choices(shape, n_acts, s, enabled, store_acted(store))
    noop_choices!(allowed, ctx)
    return grid, allowed, lost, (; ctx)
end
