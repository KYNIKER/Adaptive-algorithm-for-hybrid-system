#=
Shield synthesis for hybrid systems with several locations.
A location has its own dynamics and belongs to a layer of the grid, and can be restricted to a region of the continuous state:
  - one location, one layer, edges to itself        : a system without discrete modes (see ModalHybridSystem(::EuclideanHybridSystem))
  - one layer per location                          : discrete modes
  - several locations in a layer, different regions : dynamics that switch across the state, with edges between the regions
A cell of the grid is CartesianIndex(continuous cell..., layer), and its linear index is its id in the FlatShieldState of flat.jl, as for the other shields.

One decision period of a cell is computed by modal_reach!:
  - If its flowpipe over [0, p] touches no guard and no unsafe set, the end point of the period is the whole answer. The check uses the adaptive time steps δ⁻..δ⁺. Nothing is assumed about how the state moves between the end points.
  - Otherwise the flowpipe is walked in pieces of δ⁻. A piece touching a guard jumps on its own and continues in the target location for its own remaining time. A piece touching an unsafe set makes the transition bad. States that have not crossed yet continue in the current location.
Safety along the flowpipe comes from the unsafe sets. Leaving the state space is only checked on the sets at the end of the period, because the flowpipe pieces are bloated and a set that starts on the border of the state space would appear to leave it.

Clock axes: an axis whose row in A is zero and with no input generator (a clock, or anything else with a constant rate) only translates in that location. Two things keep such an axis exact, since otherwise time would drift by a cell per period:
  - Cells are half-open on it. A translated cell that merely touches a neighbouring cell at its border does not count as reaching it.
  - Its interval at the end of the period is tracked exactly through jumps that keep it rigid (same rate on both sides, jump only translates it) or reset it (see clock_after_jump). The successors are filtered by that interval, so the bloat of the flowpipe pieces does not reach this axis.
=#

# Action taken in a location, landing in location `target`. Several ModalActions with the same action id in one location
# (and the same guard) are one choice with several outcomes, which the shield all takes as successors.
struct ModalAction
    action::Action
    target::Int
end

# A location with starts = false is only entered by actions and edges: no cell starts a period in it, but what it reaches
# at the end of the period lands in its layer. So a layer can be a decision point that every period starts from, with
# the outcome of each choice in locations of their own (see loadCruiseControlModal_t).
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

# Everything location-specific. Polyhedra are kept as the (directions, bounds) of getHalfSpaceProjections.
struct ModalLocationData{D}
    dynamics::Dynamics # shared between locations with the same (A, U)
    layer::Int
    starts::Bool
    region::Union{Tuple{Vector{Vector{Float64}},Vector{Float64}},Nothing}
    guards::Vector{Tuple{Vector{Vector{Float64}},Vector{Float64}}}
    fast::BitArray{D}
    fast_cover::ConstantZonotopeCover{D}
    clock_rate::Vector{Float64} # rate of each axis that only translates in this location, NaN for the others
    flow_lo::Vector{Float64} # bounds on how far an axis can move within a period in this location
    flow_hi::Vector{Float64}
    incoming_edges::Vector{Tuple{Int,Edge}} # (source location, edge)
    incoming_actions::Vector{Tuple{Int,Action}} # (source location, action)
    bordered::Vector{Bool} # per face of the state space: a guard of the location borders it (see left_statespace)
end

# A set covering part of the end of the period in a location. On axes with lo[i] not NaN, the states are exactly in [lo[i], hi[i]].
struct ModalEndSet
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
    locations::Vector{ModalLocationData{D}}
    unsafe::Vector{Tuple{Vector{Vector{Float64}},Vector{Float64}}}
    statespace::Tuple{Vector{Vector{Float64}},Vector{Float64}}
    cell_generators::Matrix{Float64}
    action_ids::Vector{Int}
    may_act::Vector{Vector{BitArray{D}}} # per location, per position in action_ids: some state of the cell is inside the action's guard
    must_act::Vector{Vector{BitArray{D}}} # same, but every state of the cell is inside it
    n_layers::Int
    identity_matrix::Matrix{Float64}
    center::Vector{Float64}
    max_depth::Int
    covers::CoverCache{D}
    layer_starts::Vector{Vector{Int}} # per layer, the locations cells start a period in
    start_mask::Array{UInt64,N} # per (cell..., layer), bit j set if the cell starts in layer_starts[layer][j]
end

const MODAL_TOL = 1e-9

# The locations a cell starts a period in, as a set bit per location of its layer: iterating it allocates nothing.
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
function clock_after_jump(lo::Vector{Float64}, hi::Vector{Float64}, guard, M::Matrix{Float64}, j::Vector{Float64}, source::Vector{Float64}, target::Vector{Float64}, D′::Int, W′::Int, δ⁻::Float64)
    lo′ = fill(NaN, length(lo))
    hi′ = fill(NaN, length(lo))
    for i in eachindex(lo)
        isnan(target[i]) && continue
        least, most = minmax(target[i] * D′ * δ⁻, target[i] * (D′ + W′) * δ⁻)
        if all(iszero, view(M, i, :))
            lo′[i], hi′[i] = j[i] + least, j[i] + most
        elseif !isnan(lo[i]) && source[i] == target[i] && M[i, i] == 1.0 && count(!iszero, view(M, i, :)) == 1
            l, h = lo[i], hi[i]
            for (a, b) in zip(guard...)
                if a[i] != 0.0 && count(!iszero, a) == 1
                    a[i] > 0 ? (h = min(h, b / a[i] + most)) : (l = max(l, b / a[i] + least))
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

# True if the set at the shortest or longest remaining time leaves the state space. These sets are flowed without jumps,
# so states that jumped appear in them as if they had not:
#   - a face a ⋅ x <= b that a guard of the location borders (the guard has the halfspace -a ⋅ x <= -b, such as the wrap
#     at t = 20) is skipped: beyond it are only states that jumped, and the bloat of the sets around the jump. This assumes
#     a state can only cross such a face inside that guard.
#   - beyond any other face, a part inside a guard has jumped too; any other part has left.
# Tracked clock axes use their exact interval: X may have been clipped to a guard, and that clip is not exact, so after a jump it can sit slightly outside on such an axis without any state being there.
function left_statespace(ctx::ModalContext, loc::ModalLocationData, X, D::Int, total::Int, lo::Vector{Float64}, hi::Vector{Float64})
    dyn = loc.dynamics
    for (n, k) in enumerate((D, total))
        n == 2 && k == D && break
        final = minkowski_sum(linear_map(dyn.Φ[k+1], X), dyn.V[k+1])
        for (f, (a, b)) in enumerate(zip(ctx.statespace...))
            loc.bordered[f] && continue
            axis = findfirst(!iszero, a)
            tracked = !isnan(lo[axis])
            value = tracked ? (a[axis] > 0 ? a[axis] * hi[axis] : a[axis] * lo[axis]) : ρ(a, final)
            value > b + MODAL_TOL || continue
            beyond = clip_to_polyhedron(final, ([-a], [-b]))
            if beyond === nothing
                tracked && return true # the interval is exact, the set is not
                continue
            end
            any(guard -> allunder(beyond, guard...), loc.guards) || return true
        end
    end
    return false
end

# Appends to `out` what X reaches, covering it at the end of the period. X starts with between D and D + W steps of δ⁻ left, depending on the state, so W > 0 after a jump. lo and hi are its exact intervals at the end of the period on the clock axes.
# Returns false if the transition is bad: it touches an unsafe set, or leaves the state space when exit_bad is true.
function modal_reach!(out::Vector{ModalEndSet}, ctx::ModalContext, li::Int, X, D::Int, W::Int, exit_bad::Bool, lo::Vector{Float64}, hi::Vector{Float64}, depth::Int=0)
    depth > ctx.max_depth && return false # jumping again and again: assume the worst
    total = D + W
    if total == 0
        push!(out, ModalEndSet(li, X, lo, hi))
        return true
    end
    loc = ctx.locations[li]
    dyn = loc.dynamics
    δ⁺ = min(ctx.δ⁺, ctx.δ⁻ * 2^floor(Int, log2(total))) # no step longer than the time left
    disc = newReACTDiscretizePlus(X, ctx.δ⁻, δ⁺, dyn.A, dyn.P2A_abs, dyn.phiDict, dyn.inputDict, ctx.alg, ctx.maxOrder, ctx.reduceOrder)
    k0 = first_touch(dyn, disc, Iterators.flatten((loc.guards, ctx.unsafe)), total, ctx.δ⁻, δ⁺, ctx.identity_matrix)

    # Pieces from the first touch on, and the pieces that end the period when the time left varies.
    for k in min(k0, W == 0 ? total : D):total-1
        piece = minkowski_sum(linear_map(dyn.Φ[k+1], disc[ctx.δ⁻]), dyn.V[k+1])
        if k >= k0
            any(set -> someinside(piece, set...), ctx.unsafe) && return false
            all_jumped = false
            for (e, guard) in zip(ctx.system.locations[li].edges, loc.guards)
                someinside(piece, guard...) || continue
                clipped = clip_to_polyhedron(piece, guard)
                clipped === nothing && continue
                # A state in this piece has spent between k and k + 1 steps here, and had between D and D + W left at the start.
                D′ = max(D - k - 1, 0)
                W′ = total - k - D′
                clocks = clock_after_jump(lo, hi, guard, e.jumpMatrix, e.jumpVector, loc.clock_rate, ctx.locations[e.targetLoc].clock_rate, D′, W′, ctx.δ⁻)
                clocks === nothing && continue
                jumped = jump_set(clipped, e.jumpMatrix, e.jumpVector, ctx.maxOrder, ctx.reduceOrder)
                modal_reach!(out, ctx, e.targetLoc, jumped, D′, W′, exit_bad, clocks..., depth + 1) || return false
                all_jumped = all_jumped || allunder(piece, guard...)
            end
            all_jumped && return true
        end
        if W > 0 && k >= D
            push!(out, ModalEndSet(li, piece, lo, hi))
        end
    end
    exit_bad && left_statespace(ctx, loc, X, D, total, lo, hi) && return false
    if W == 0
        push!(out, ModalEndSet(li, minkowski_sum(linear_map(dyn.Φ[D+1], X), dyn.V[D+1]), lo, hi))
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

function ModalContext(system::ModalHybridSystem, p, granularity::Float64, δ⁻::Float64, δ⁺::Float64, alg, maxOrder::Int, reduceOrder::Int, sweep_factor, max_depth::Int)
    grid = Grid(system.statespace, granularity)
    return _ModalContext(grid, Val(ndims(grid.deadCells) + 1), system, Float64(p), δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor, max_depth)
end

# Function barrier on Grid{D}. N = D + 1 is the dimension of the cell index, (continuous cell..., layer).
function _ModalContext(grid::Grid{D}, ::Val{N}, system, p::Float64, δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor, max_depth) where {D,N}
    steps = round(Int, p / δ⁻)
    isapprox(steps * δ⁻, p; atol=MODAL_TOL) || throw(ArgumentError("The period $p must be a multiple of δ⁻ = $δ⁻"))
    isinteger(log2(δ⁺ / δ⁻)) || throw(ArgumentError("δ⁺ / δ⁻ must be a power of two"))

    faces = hyperrectangle_to_HPolyhedron(system.statespace)
    statespace = getHalfSpaceProjections(faces)
    exits = [([-Vector{Float64}(c.a)], [-c.b - MODAL_TOL]) for c in constraints_list(faces)] # each face, as the set strictly beyond it
    unsafe = [getHalfSpaceProjections(H) for H in system.unsafe]

    g = grid.granularity
    cell_generators = Matrix{Float64}(genmat(remove_zero_generators(Zonotope(zeros(D), diagm(map(n -> n == 1 ? 0.0 : g / 2, grid.numCells))))))
    action_ids = sort!(unique(a.action.id for loc in system.locations for a in loc.actions))
    identity_matrix = Matrix{Float64}(I, D, D)

    shared = Dict{Tuple{Matrix{Float64},Vector{Float64},Matrix{Float64}},Dynamics}()
    sweep = sweep_factor isa Integer ? [min(sweep_factor, n) for n in grid.numCells] : [min(s, n) for (s, n) in zip(sweep_factor, grid.numCells)]
    G_big = Matrix{Float64}(genmat(remove_zero_generators(Zonotope(zeros(D), diagm(map((n, s) -> n == 1 ? 0.0 : g * (s / 2), grid.numCells, sweep))))))

    locations = ModalLocationData{D}[]
    for loc in system.locations
        dyn = get!(() -> Dynamics(loc.A, loc.U, steps, δ⁻, δ⁺, alg, maxOrder, reduceOrder), shared, (loc.A, LazySets.center(loc.U), genmat(loc.U)))
        guards = [getHalfSpaceProjections(e.guard) for e in loc.edges]
        targets = vcat(guards, unsafe, exits)
        region = loc.region === nothing ? nothing : getHalfSpaceProjections(loc.region)
        # Coarse sweep: a whole block is fast if its flowpipe touches nothing over the period.
        # Cells start a period in the location, or land in it with an action at the start of one.
        entered_at_start = loc.starts || any(a.target == length(locations) + 1 for l in system.locations for a in l.actions)
        fast = falses(cellshape(grid))
        for bidx in CartesianIndices(ntuple(i -> 1:fld(grid.numCells[i], sweep[i]), D))
            entered_at_start || break
            center = grid.lower .+ g .* (sweep ./ 2) .+ (car2vec(bidx) .- 1) .* g .* sweep
            box_overlaps_region(region, center .- g .* (sweep ./ 2), center .+ g .* (sweep ./ 2)) || continue
            disc = newReACTDiscretizePlus(Zonotope(center, G_big), δ⁻, δ⁺, dyn.A, dyn.P2A_abs, dyn.phiDict, dyn.inputDict, alg, maxOrder, reduceOrder)
            if first_touch(dyn, disc, targets, steps, δ⁻, δ⁺, identity_matrix) == steps
                fast[ntuple(i -> (1+(bidx[i]-1)*sweep[i]):(bidx[i]*sweep[i]), D)...] .= true
            end
        end
        # All fast cells reach a translate of the same zonotope, bloated by half a cell for the cover test.
        Vp = dyn.V[steps+1]
        fast_cover = ConstantZonotopeCover(remove_zero_generators(Zonotope(LazySets.center(Vp), hcat(genmat(Vp), dyn.Φ[steps+1] * cell_generators, Matrix((g / 2) * I, D, D)))), Val(D))
        clock_rate = [all(iszero, view(loc.A, i, :)) && all(iszero, view(genmat(loc.U), i, :)) ? LazySets.center(loc.U)[i] : NaN for i in 1:D]
        # How far any state of the grid can move on each axis at any time in the period, by interval arithmetic.
        flow_lo, flow_hi = fill(Inf, D), fill(-Inf, D)
        for k in 0:steps
            M = dyn.Φ[k+1] - I
            c = M * ((grid.lower .+ grid.upper) ./ 2) + LazySets.center(dyn.V[k+1])
            r = abs.(M) * ((grid.upper .- grid.lower) ./ 2) + vec(sum(abs, genmat(dyn.V[k+1]); dims=2))
            flow_lo .= min.(flow_lo, c .- r)
            flow_hi .= max.(flow_hi, c .+ r)
        end
        li = length(locations) + 1
        incoming_edges = [(src, e) for (src, l) in enumerate(system.locations) for e in l.edges if e.targetLoc == li]
        incoming_actions = [(src, a.action) for (src, l) in enumerate(system.locations) for a in l.actions if a.target == li]
        bordered = [bordered_by(guards, a, b; atol=MODAL_TOL) for (a, b) in zip(statespace...)]
        push!(locations, ModalLocationData{D}(dyn, loc.layer, loc.starts, region, guards, fast, fast_cover, clock_rate, flow_lo, flow_hi, incoming_edges, incoming_actions, bordered))
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
    start_mask = zeros(UInt64, cellshape(grid)..., n_layers)
    lo, hi = zeros(D), zeros(D)
    for layer in 1:n_layers
        ids = layer_starts[layer]
        length(ids) <= 64 || throw(ArgumentError("layer $layer has $(length(ids)) locations that cells start in; at most 64 are supported"))
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

    return ModalContext{D,N,typeof(system)}(system, grid, steps, δ⁻, δ⁺, alg, maxOrder, reduceOrder, locations, unsafe, statespace, cell_generators, action_ids, may_act, must_act, n_layers, identity_matrix, zeros(D), max_depth, CoverCache{D}(), layer_starts, start_mask)
end

# The locations a cell of the layer starts a period in, in increasing order.
cell_locations(ctx::ModalContext{D}, layer::Int, c::CartesianIndex{D}) where D = LocationSet(ctx.layer_starts[layer], ctx.start_mask[c, layer])

# True if the half-open cell overlaps the exact interval on every tracked axis. A point on a cell border belongs to the cell above it.
function cell_in_clock_intervals(grid::Grid{D}, cell::CartesianIndex{D}, lo::Vector{Float64}, hi::Vector{Float64}) where D
    for i in 1:D
        isnan(lo[i]) && continue
        first_cell = floor(Int, (lo[i] - grid.lower[i]) / grid.granularity + MODAL_TOL) + 1
        last_cell = max(ceil(Int, (hi[i] - grid.lower[i]) / grid.granularity - MODAL_TOL), first_cell)
        first_cell <= cell[i] <= last_cell || return false
    end
    return true
end

function append_touching!(res::Vector{CartesianIndex{N}}, ctx::ModalContext{D,N}, out::Vector{ModalEndSet}) where {D,N}
    grid = ctx.grid
    for e in out
        # get_touching_cell_idxs_lz finds no cells for a set with zero width exactly on the lower edge of the grid on some axis, because the upper index is ceil(0) = 0. A reset to that edge (like l := 0) gives such a set, so it is moved just inside the grid.
        c = copy(LazySets.center(e.set))
        G = genmat(e.set)
        for i in 1:D
            if abs(c[i] - grid.lower[i]) <= MODAL_TOL && all(iszero, view(G, i, :))
                c[i] = grid.lower[i] + 10 * MODAL_TOL
            end
        end
        # Only the cells in the exact intervals of the clock axes are tested against the set
        for cell in get_touching_cell_idxs_lz(grid, Zonotope(c, G), ctx.covers, clock_cell_box(grid, e.lo, e.hi)...)
            if cell_in_clock_intervals(grid, cell, e.lo, e.hi)
                push!(res, CartesianIndex(cell, ctx.locations[e.location].layer))
            end
        end
    end
    return res
end

# The box of cells that overlap the exact intervals [lo, hi] on the clock axes (all cells on the other axes), as the
# keep_lo and keep_hi of get_touching_cell_idxs_cover. The same cells as cell_in_clock_intervals accepts.
function clock_cell_box(grid::Grid{D}, lo::Vector{Float64}, hi::Vector{Float64}) where D
    keep_lo = ntuple(i -> isnan(lo[i]) ? 1 : floor(Int, (lo[i] - grid.lower[i]) / grid.granularity + MODAL_TOL) + 1, D)
    keep_hi = ntuple(i -> isnan(lo[i]) ? grid.numCells[i] : max(ceil(Int, (hi[i] - grid.lower[i]) / grid.granularity - MODAL_TOL), keep_lo[i]), D)
    return keep_lo, keep_hi
end

# Appends what a fast cell c of location li reaches at the end of the period, from the location's precomputed cover.
function append_fast!(res::Vector{CartesianIndex{N}}, ctx::ModalContext{D,N}, li::Int, c::CartesianIndex{D}, lo::Vector{Float64}, hi::Vector{Float64}) where {D,N}
    loc = ctx.locations[li]
    fast_cell_center!(ctx.center, loc.dynamics.Φ[ctx.steps+1], LazySets.center(loc.dynamics.V[ctx.steps+1]), ctx.grid, c)
    for cell in get_touching_cell_idxs_cover(ctx.grid, loc.fast_cover, ctx.center, clock_cell_box(ctx.grid, lo, hi)...)
        if cell_in_clock_intervals(ctx.grid, cell, lo, hi)
            push!(res, CartesianIndex(cell, loc.layer))
        end
    end
    return res
end

# The offset in cells of a jump that maps every cell onto another cell (identity matrix, jump vector a whole number of
# cells on every axis), or nothing.
function cell_shift(M::Matrix{Float64}, j::Vector{Float64}, grid::Grid{D}) where D
    for k in 1:D, i in 1:D
        M[i, k] == (i == k ? 1.0 : 0.0) || return nothing
    end
    shift = ntuple(i -> j[i] / grid.granularity, D)
    all(i -> abs(shift[i] - round(shift[i])) <= MODAL_TOL, 1:D) || return nothing
    return CartesianIndex(ntuple(i -> round(Int, shift[i]), D))
end

# Flow successors of a cell, or nothing if flowing from it is bad.
function flow_successors(ctx::ModalContext{D,N}, layer::Int, c::CartesianIndex{D}) where {D,N}
    res = CartesianIndex{N}[]
    out = ModalEndSet[]
    Z = Zonotope(ctx.grid.lower .+ (car2vec(c) .- 0.5) .* ctx.grid.granularity, ctx.cell_generators)
    for li in cell_locations(ctx, layer, c)
        loc = ctx.locations[li]
        lo, hi = clock_end(loc.clock_rate, Z, ctx.steps, ctx.δ⁻)
        if loc.fast[c]
            append_fast!(res, ctx, li, c, lo, hi)
        else
            empty!(out)
            modal_reach!(out, ctx, li, Z, ctx.steps, 0, true, lo, hi) || return nothing
            append_touching!(res, ctx, out)
        end
    end
    return unique!(res)
end

# Successors of taking the action at position ai of ctx.action_ids in a cell, or nothing if the action is not available, is bad, or does nothing.
# Trying an action outside its guard, or in a location that does not offer it, does nothing: the state just flows. So:
#   - a cell where it does nothing in every location covering it gets nothing, and its choice for the action is the same as waiting (see noop_choices!)
#   - a cell partly outside the guard (or covered by a location where it does nothing) reaches the action successors and the flow successors `flow`
function action_successors(ctx::ModalContext{D,N}, layer::Int, c::CartesianIndex{D}, ai::Int, exit_bad::Bool, flow) where {D,N}
    res = CartesianIndex{N}[]
    out = ModalEndSet[]
    Z = Zonotope(ctx.grid.lower .+ (car2vec(c) .- 0.5) .* ctx.grid.granularity, ctx.cell_generators)
    acts = false
    also_flows = false
    for li in cell_locations(ctx, layer, c)
        # may_act is false where the location does not offer the action
        if !ctx.may_act[li][ai][c]
            also_flows = true
            continue
        end
        ctx.must_act[li][ai][c] || (also_flows = true)
        acts = true
        # Several actions with this id in one location are one choice with several outcomes: all of them are successors
        for a in ctx.system.locations[li].actions
            a.action.id == ctx.action_ids[ai] || continue
            jumped = jump_set(Z, a.action.jumpMatrix, a.action.jumpVector, ctx.maxOrder, ctx.reduceOrder)
            lo, hi = clock_end(ctx.locations[a.target].clock_rate, jumped, ctx.steps, ctx.δ⁻)
            # A jump onto another cell that is fast in the target location: its precomputed cover, as for the flow
            shift = cell_shift(a.action.jumpMatrix, a.action.jumpVector, ctx.grid)
            if shift !== nothing && checkbounds(Bool, ctx.locations[a.target].fast, c + shift) && ctx.locations[a.target].fast[c+shift]
                append_fast!(res, ctx, a.target, c + shift, lo, hi)
                continue
            end
            empty!(out)
            modal_reach!(out, ctx, a.target, jumped, ctx.steps, 0, exit_bad, lo, hi) || return nothing
            append_touching!(res, ctx, out)
        end
    end
    acts || return nothing
    if also_flows
        flow === nothing && return nothing
        append!(res, flow)
    end
    return unique!(res)
end

# Where an action does nothing, it is allowed exactly when waiting is. An action does nothing in a cell when, in every
# location covering it, the location does not offer the action or the cell is outside its guard.
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

# The online frontier. Grid neighbours do not find the predecessors of a cell: a period can move a cell many cells away, and a jump or action can change layer and reset an axis. So it searches backward through the system with boxes.

# Shrinks the box [lo, hi] by the axis-aligned halfspaces of a polyhedron, each widened by margin. False if the box becomes empty.
function clip_box!(lo::Vector{Float64}, hi::Vector{Float64}, polyhedron, margin::Float64)
    polyhedron === nothing && return true
    for (a, b) in zip(polyhedron...)
        count(!iszero, a) == 1 || continue
        i = findfirst(!iszero, a)
        a[i] > 0 ? (hi[i] = min(hi[i], b / a[i] + margin)) : (lo[i] = max(lo[i], b / a[i] - margin))
    end
    return all(lo .<= hi)
end

# Box of states the jump can map into [lo, hi], or nothing if there are none: a translated axis moves back, a reset axis can have had any value if the reset lands in the box, and any other axis is not bounded.
function jump_preimage(lo::Vector{Float64}, hi::Vector{Float64}, M::Matrix{Float64}, j::Vector{Float64}, grid::Grid)
    lo′, hi′ = copy(grid.lower), copy(grid.upper)
    for i in eachindex(lo)
        if all(iszero, view(M, i, :))
            lo[i] - MODAL_TOL <= j[i] <= hi[i] + MODAL_TOL || return nothing
        elseif M[i, i] == 1.0 && count(!iszero, view(M, i, :)) == 1
            lo′[i], hi′[i] = lo[i] - j[i], hi[i] - j[i]
        end
    end
    return lo′, hi′
end

# What the online search has already done, so that it is not done again: queued cells are never unqueued, so a search or
# a scan contained in an earlier one adds nothing. Up to MEMO_SIZE entries are kept per location and per layer, the
# oldest replaced first; forgetting one only means a search is repeated.
const MEMO_SIZE = 32
struct FrontierMemo{D}
    searched::Vector{Vector{Tuple{Int,Vector{Float64},Vector{Float64}}}} # per location: (depth, lo, hi) searched from
    scanned::Vector{Vector{Tuple{NTuple{D,Int},NTuple{D,Int}}}} # per layer: (first, last) cell of a box of queued cells
    next::Vector{Int} # per location, then per layer: the entry to replace next once the list is full
end
FrontierMemo{D}(n_locations::Int, n_layers::Int) where D = FrontierMemo{D}([Tuple{Int,Vector{Float64},Vector{Float64}}[] for _ in 1:n_locations],
    [Tuple{NTuple{D,Int},NTuple{D,Int}}[] for _ in 1:n_layers], ones(Int, n_locations + n_layers))

function remember!(list::Vector{T}, memo::FrontierMemo, slot::Int, entry::T) where T
    if length(list) < MEMO_SIZE
        push!(list, entry)
    else
        list[memo.next[slot]] = entry
        memo.next[slot] = memo.next[slot] % MEMO_SIZE + 1
    end
end

function queue_cells!(frontier::Vector{CartesianIndex{N}}, queued::BitArray{N}, grid::Grid{D}, layer::Int, lo::Vector{Float64}, hi::Vector{Float64}, memo::FrontierMemo{D}) where {D,N}
    first = ntuple(i -> max(floor(Int, (lo[i] - grid.lower[i]) / grid.granularity) + 1, 1), D)
    last = ntuple(i -> min(ceil(Int, (hi[i] - grid.lower[i]) / grid.granularity), grid.numCells[i]), D)
    scanned = memo.scanned[layer]
    for (f, l) in scanned
        all(i -> f[i] <= first[i] && last[i] <= l[i], 1:D) && return
    end
    for c in CartesianIndices(ntuple(i -> first[i]:last[i], D))
        k = CartesianIndex(c, layer)
        if !queued[k]
            queued[k] = true
            push!(frontier, k)
        end
    end
    remember!(scanned, memo, length(memo.searched) + layer, (first, last))
end

# Queues every cell that may reach the box [lo, hi] in location li within a period: cells flowing there in li, cells taking an action into li, and, through the edges into li, cells flowing there from other locations.
# A search only adds cells to the queue, and one from a box inside an earlier box of the same location, at the same or a
# greater depth, adds none the earlier one did not (memo). So it stops there: without that, a location with an edge to
# itself (or a cycle of edges) is searched again at every depth, once per path, which grows exponentially with max_depth.
function queue_predecessors!(frontier, queued, ctx::ModalContext, li::Int, lo::Vector{Float64}, hi::Vector{Float64}, memo::FrontierMemo, depth::Int=0)
    searched = memo.searched[li]
    for (d, slo, shi) in searched
        d <= depth && all(i -> slo[i] <= lo[i] && hi[i] <= shi[i], eachindex(lo)) && return
    end
    remember!(searched, memo, li, (depth, copy(lo), copy(hi)))
    grid = ctx.grid
    g = grid.granularity
    loc = ctx.locations[li]
    from_lo, from_hi = lo .- loc.flow_hi .- g, hi .- loc.flow_lo .+ g
    clip_box!(from_lo, from_hi, loc.region, g) || return
    loc.starts && queue_cells!(frontier, queued, grid, loc.layer, from_lo, from_hi, memo)
    for (src, act) in loc.incoming_actions
        box = jump_preimage(from_lo, from_hi, act.jumpMatrix, act.jumpVector, grid)
        box === nothing && continue
        clip_box!(box..., ctx.locations[src].region, g) && queue_cells!(frontier, queued, grid, ctx.locations[src].layer, box..., memo)
    end
    depth >= ctx.max_depth && return
    for (src, e) in loc.incoming_edges
        box = jump_preimage(from_lo, from_hi, e.jumpMatrix, e.jumpVector, grid)
        box === nothing && continue
        source = ctx.locations[src]
        # A flowpipe piece can pass a guard by one step of δ⁻ before it is clipped to it.
        margin = g + maximum(max.(abs.(source.flow_lo), abs.(source.flow_hi))) / ctx.steps
        if clip_box!(box..., source.region, g) && clip_box!(box..., getHalfSpaceProjections(e.guard), margin)
            queue_predecessors!(frontier, queued, ctx, src, box..., memo, depth + 1)
        end
    end
end

# Shield for a ModalHybridSystem. Computes every cell of every layer, or with online = true only the cells that may reach an unsafe cell or a cell that has become bad, the others being safe. Then pushes the bad events backwards.
# include_border_cells also starts the online search from the border of the grid, for when leaving the state space should count as unsafe.
# With exit_bad_for_actions = false an action may leave the state space, and only its successors inside the grid count.
# Returns grid, allowed, lost, info like every shield (see choice_bits in flat.jl). grid is the continuous
# grid; allowed and lost are indexed by (continuous cell..., layer), and a cell is dead in a layer where allowed is 0.
# Trying an action outside its guard or in a location that does not offer it does nothing, so there it is allowed exactly when waiting is (see action_successors).
# A cell is enabled for an action when some location covering it can take the action there. info.ctx is the context.
function ReACTed_reachable_cell_modal(system::ModalHybridSystem, p, granularity, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5; online::Bool=false, include_border_cells::Bool=false, act_needs_safe_flow::Bool=true, exit_bad_for_actions::Bool=true, sweep_factor=20, max_depth::Int=16, verbose::Bool=true)
    timer = PhaseTimer()
    ctx = ModalContext(system, p, Float64(granularity), δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor, max_depth)
    return _ReACTed_reachable_cell_modal(ctx, online, include_border_cells, act_needs_safe_flow, exit_bad_for_actions, verbose, timer)
end

function _ReACTed_reachable_cell_modal(ctx::ModalContext{D,N}, online::Bool, include_border_cells::Bool, act_needs_safe_flow::Bool, exit_bad_for_actions::Bool, verbose::Bool, timer::PhaseTimer) where {D,N}
    grid = ctx.grid
    g = grid.granularity
    shape = (cellshape(grid)..., ctx.n_layers)
    verbose && println("$(prod(shape)) cells, shape $shape")
    n_acts = length(ctx.action_ids)
    # The shield works on ids: the linear index of (continuous cell..., layer). See flat.jl.
    n = prod(shape)
    n < typemax(Int32) || throw(ArgumentError("$n cells do not fit in Int32 ids"))
    cells = CartesianIndices(shape) # id -> cell
    ids = LinearIndices(shape) # cell -> id
    store = SuccessorStore(n, [falses(n) for _ in 1:n_acts])
    s = reset!(FlatShieldState(n_acts, act_needs_safe_flow), n)
    buf = Int32[]
    queued = falses(shape)
    frontier = CartesianIndex{N}[]
    next_frontier = CartesianIndex{N}[]

    memo = FrontierMemo{D}(length(ctx.locations), ctx.n_layers)
    function queue_predecessors_of!(c::CartesianIndex{D}, layer::Int)
        lo = grid.lower .+ (car2vec(c) .- 1) .* g
        # Where a period can end in the cell: the locations it starts in, and those that no cell starts in
        for li in Iterators.flatten((cell_locations(ctx, layer, c), (li for (li, loc) in enumerate(ctx.locations) if loc.layer == layer && !loc.starts)))
            queue_predecessors!(next_frontier, queued, ctx, li, lo, lo .+ g, memo)
        end
    end

    for H in ctx.system.unsafe, c in get_touching_cell_idxs_l(grid, H), layer in 1:ctx.n_layers
        k = CartesianIndex(c, layer)
        seed_constraint!(s, Int32(ids[k]))
        store.flow_none[ids[k]] = true
        queued[k] = true
    end
    if online
        foreach(k -> queue_predecessors_of!(CartesianIndex(ntuple(i -> cells[k][i], D)), cells[k][N]), findall(s.dead))
        if include_border_cells
            for c in border_cells(grid), layer in 1:ctx.n_layers
                queued[CartesianIndex(c, layer)] || (queued[CartesianIndex(c, layer)] = true; push!(next_frontier, CartesianIndex(c, layer)))
                queue_predecessors_of!(c, layer)
            end
        end
    else
        for k in CartesianIndices(shape)
            queued[k] || push!(next_frontier, k)
        end
    end

    count_cells = 0
    rounds = 0
    lap!(timer, "Setup")
    while !isempty(next_frontier)
        frontier, next_frontier = next_frontier, empty!(frontier)
        rounds += 1
        for k in frontier
            c = CartesianIndex(ntuple(i -> k[i], D))
            layer = k[N]
            id = Int32(ids[k])
            (s.computed[id] || isempty(cell_locations(ctx, layer, c))) && continue
            flow = flow_successors(ctx, layer, c)
            store.flow_none[id] = flow === nothing
            set_successors!(store, 0, id, flow === nothing ? empty!(buf) : map!(v -> Int32(ids[v]), resize!(buf, length(flow)), flow))
            for ai in 1:n_acts
                succ = action_successors(ctx, layer, c, ai, exit_bad_for_actions, flow)
                store.act_enabled[ai][id] = succ !== nothing
                set_successors!(store, ai, id, succ === nothing ? empty!(buf) : map!(v -> Int32(ids[v]), resize!(buf, length(succ)), succ))
            end
            register_successors!(store, s, id)
            count_cells += 1
            verbose && count_cells % 10000 == 0 && println("$count_cells cells computed : $(Dates.format(now(), "HH:MM:SS"))")
        end
        lap!(timer, "Reachability")
        propagate_bad!(s)
        if online
            foreach(k -> queue_predecessors_of!(CartesianIndex(ntuple(i -> cells[k][i], D)), cells[k][N]), unique!(s.changed))
        end
        empty!(s.changed)
        lap!(timer, "Shield")
    end
    verbose && report(timer, "Rounds: $rounds, cells computed: $count_cells, dead: $(count(s.dead))")
    function enabled(ai, k)
        c = CartesianIndex(ntuple(i -> cells[k][i], D))
        locs = cell_locations(ctx, cells[k][N], c)
        return any(li -> ctx.may_act[li][ai][c], locs)
    end
    allowed, lost = cell_choices(shape, n_acts, s, enabled, store_acted(store))
    noop_choices!(allowed, ctx)
    return grid, allowed, lost, (; ctx)
end
