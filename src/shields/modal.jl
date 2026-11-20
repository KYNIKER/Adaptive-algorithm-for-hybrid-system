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
    # Discretizations over a whole period of a cell and of each shape an incoming action gives a cell (J times it), as
    # (generators, discretization). Shared with the locations of the same dynamics.
    discs::Vector{Tuple{Matrix{Float64},CellDiscretization{D}}}
    # Per outgoing action (as in the location's actions): the cover of a cell just after the jump, to find the cells it
    # lands in, and the cover of where it is at the end of the period when those cells are all fast in the target
    act_start_cover::Vector{ConstantZonotopeCover{D}}
    act_cover::Vector{ConstantZonotopeCover{D}}
    # The same for boxes of cells (see box_successors), per box shape (its generators), made when first used: the cover of
    # a fast box's end set, and per outgoing action the covers of the jumped box at the start and the end of the period
    box_covers::Dict{Matrix{Float64},Tuple{ConstantZonotopeCover{D},Vector{ConstantZonotopeCover{D}},Vector{ConstantZonotopeCover{D}}}}
    layer::Int
    starts::Bool
    region::Union{Tuple{Vector{Vector{Float64}},Vector{Float64}},Nothing}
    guards::Vector{Tuple{Vector{Vector{Float64}},Vector{Float64}}}
    fast::BitArray{D}
    fast_cover::ConstantZonotopeCover{D}
    clock_rate::Vector{Float64} # rate of each axis that only translates in this location, NaN for the others
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
    sweep::Vector{Int} # cells per block and dimension
    exit_bad_for_flow::Bool # leaving the state space makes waiting bad
    exit_bad_for_actions::Bool # leaving the state space makes an action bad
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
function shape_discretization(loc::ModalLocationData, G)
    for (S, cd) in loc.discs
        S == G && return cd
    end
    return nothing
end

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
    # A whole cell, or a cell just after an action, flowing a whole period uses the precomputed discretization of its
    # shape. Only at depth 0: the shared one is overwritten per call, and this call still reads it after the nested calls
    # return.
    cd = depth == 0 && total == ctx.steps ? shape_discretization(loc, genmat(X)) : nothing
    disc = cd === nothing ? newReACTDiscretizePlus(X, ctx.δ⁻, δ⁺, dyn.A, dyn.P2A_abs, dyn.phiDict, dyn.inputDict, ctx.alg, ctx.maxOrder, ctx.reduceOrder) :
           discretize_cell!(cd, LazySets.center(X))
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

function ModalContext(system::ModalHybridSystem, p, granularity::Float64, δ⁻::Float64, δ⁺::Float64, alg, maxOrder::Int, reduceOrder::Int, sweep_factor, max_depth::Int;
                      exit_bad_for_flow::Bool=true, exit_bad_for_actions::Bool=true)
    grid = Grid(system.statespace, granularity)
    return _ModalContext(grid, Val(ndims(grid.deadCells) + 1), system, Float64(p), δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor, max_depth, exit_bad_for_flow, exit_bad_for_actions)
end

# Function barrier on Grid{D}. N = D + 1 is the dimension of the cell index, (continuous cell..., layer).
function _ModalContext(grid::Grid{D}, ::Val{N}, system, p::Float64, δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor, max_depth, exit_bad_for_flow, exit_bad_for_actions) where {D,N}
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

    sweep = sweep_factor isa Integer ? [min(sweep_factor, n) for n in grid.numCells] : [min(s, n) for (s, n) in zip(sweep_factor, grid.numCells)]
    G_big = Matrix{Float64}(genmat(remove_zero_generators(Zonotope(zeros(D), diagm(map((n, s) -> n == 1 ? 0.0 : g * (s / 2), grid.numCells, sweep))))))
    cell_δ⁺ = min(δ⁺, δ⁻ * 2^floor(Int, log2(steps))) # what modal_reach! uses for a whole period
    dyn_key(l) = (l.A, LazySets.center(l.U), genmat(l.U))
    shared_dyn = Dict{Tuple{Matrix{Float64},Vector{Float64},Matrix{Float64}},Dynamics}()
    dynamics_of(l) = get!(() -> Dynamics(l.A, l.U, steps, δ⁻, δ⁺, alg, maxOrder, reduceOrder), shared_dyn, dyn_key(l))
    shared_disc = Dict{Tuple{Tuple{Matrix{Float64},Vector{Float64},Matrix{Float64}},Matrix{Float64}},CellDiscretization{D}}()
    function disc_of(l, G)
        get!(shared_disc, (dyn_key(l), G)) do
            dyn = dynamics_of(l)
            reach = (; δ⁻, δ⁺=cell_δ⁺, dyn.A, dyn.P2A_abs, dyn.phiDict, dyn.inputDict, U=l.U)
            CellDiscretization(Zonotope(zeros(D), G), reach, Val(D))
        end
    end
    # The generators of a cell after the jump x := J x + j (the same for every cell)
    jumped_shape(J, G=cell_generators) = Matrix{Float64}(genmat(jump_set(Zonotope(zeros(D), G), Matrix{Float64}(J), zeros(D), maxOrder, reduceOrder)))
    exits_and(guards) = vcat(guards, unsafe, exits)
    block_ranges(bidx) = ntuple(i -> (1+(bidx[i]-1)*sweep[i]):(bidx[i]*sweep[i]), D)
    block_center(bidx) = grid.lower .+ g .* (sweep ./ 2) .+ (car2vec(bidx) .- 1) .* g .* sweep
    blocks = CartesianIndices(ntuple(i -> 1:fld(grid.numCells[i], sweep[i]), D))

    locations = ModalLocationData{D}[]
    for loc in system.locations
        dyn = dynamics_of(loc)
        li = length(locations) + 1
        incoming = [a.action.jumpMatrix for l in system.locations for a in l.actions if a.target == li]
        shapes = unique(vcat([cell_generators], [jumped_shape(J) for J in incoming]))
        discs = [(G, disc_of(loc, G)) for G in shapes]
        guards = [getHalfSpaceProjections(e.guard) for e in loc.edges]
        targets = exit_bad_for_flow ? exits_and(guards) : vcat(guards, unsafe)
        region = loc.region === nothing ? nothing : getHalfSpaceProjections(loc.region)
        # Coarse sweep: a whole block is fast if its flowpipe touches nothing over the period.
        # Cells start a period in the location, or land in it with an action at the start of one.
        entered_at_start = loc.starts || any(a.target == length(locations) + 1 for l in system.locations for a in l.actions)
        fast = falses(cellshape(grid))
        big = entered_at_start && !isempty(blocks) ? disc_of(loc, G_big) : nothing
        for bidx in blocks
            entered_at_start || break
            center = block_center(bidx)
            box_overlaps_region(region, center .- g .* (sweep ./ 2), center .+ g .* (sweep ./ 2)) || continue
            if first_touch(dyn, discretize_cell!(big, center), targets, steps, δ⁻, cell_δ⁺, identity_matrix) == steps
                fast[block_ranges(bidx)...] .= true
            end
        end
        # Per outgoing action, the covers of a jumped cell at the start and, flowing freely in the target, at the end of the period
        act_start_cover, act_cover = ConstantZonotopeCover{D}[], ConstantZonotopeCover{D}[]
        bloat = Matrix((g / 2) * I, D, D)
        cover(G) = ConstantZonotopeCover(hcat(G, bloat), ntuple(i -> sum(abs, view(G, i, :)), D), Val(D)) # bounded by the unbloated set
        for a in loc.actions
            tdyn = dynamics_of(system.locations[a.target])
            G_J = jumped_shape(a.action.jumpMatrix)
            push!(act_start_cover, cover(G_J))
            push!(act_cover, cover(hcat(genmat(tdyn.V[steps+1]), tdyn.Φ[steps+1] * G_J)))
        end
        # All fast cells reach a translate of the same zonotope, bloated by half a cell for the cover test. As in
        # get_touching_cell_idxs_lz, the candidate cells are bounded by the unbloated set, so a cell that only shares a face
        # with the end set is not a successor, the same as for a cell flowed on its own.
        Vp = dyn.V[steps+1]
        G_fast = hcat(genmat(Vp), dyn.Φ[steps+1] * cell_generators)
        fast_cover = cover(G_fast)
        clock_rate = [all(iszero, view(loc.A, i, :)) && all(iszero, view(genmat(loc.U), i, :)) ? LazySets.center(loc.U)[i] : NaN for i in 1:D]
        bordered = [bordered_by(guards, a, b; atol=MODAL_TOL) for (a, b) in zip(statespace...)]
        push!(locations, ModalLocationData{D}(dyn, discs, act_start_cover, act_cover, Dict{Matrix{Float64},Tuple{ConstantZonotopeCover{D},Vector{ConstantZonotopeCover{D}},Vector{ConstantZonotopeCover{D}}}}(), loc.layer, loc.starts, region, guards, fast, fast_cover, clock_rate, bordered))
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

    return ModalContext{D,N,typeof(system)}(system, grid, steps, δ⁻, δ⁺, alg, maxOrder, reduceOrder, locations, unsafe, statespace, cell_generators, action_ids, may_act, must_act, n_layers, identity_matrix, zeros(D), max_depth, CoverCache{D}(), layer_starts, start_mask, sweep, exit_bad_for_flow, exit_bad_for_actions)
end

# Per (cell..., layer), bit j set if the cell overlaps the region of location layer_starts[layer][j].
function start_masks(grid::Grid{D}, locations::Vector{ModalLocationData{D}}, layer_starts::Vector{Vector{Int}}) where D
    n_layers = length(layer_starts)
    g = grid.granularity
    start_mask = zeros(UInt64, cellshape(grid)..., n_layers)
    lo, hi = zeros(D), zeros(D)
    for layer in 1:n_layers
        ids = layer_starts[layer]
        length(ids) <= 64 || throw(ArgumentError("layer $layer has $(length(ids)) locations that cells start in; at most 64 are supported"))
        if all(li -> locations[li].region === nothing, ids) # every cell starts in all of them
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

# The locations a cell of the layer starts a period in, in increasing order.
cell_locations(ctx::ModalContext{D}, layer::Int, c::CartesianIndex{D}) where D = LocationSet(ctx.layer_starts[layer], ctx.start_mask[c, layer])

# True if some location the cell starts in can take the action at position ai there.
function can_act(ctx::ModalContext{D}, layer::Int, c::CartesianIndex{D}, ai::Int) where D
    for li in cell_locations(ctx, layer, c)
        ctx.may_act[li][ai][c] && return true
    end
    return false
end

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
    return append_cover!(res, ctx, loc.fast_cover, loc.layer, lo, hi)
end

# True if the cover translated to c lies in the grid and only meets cells that are fast in location loc.
function lands_in_fast_cells(ctx::ModalContext{D}, cov::ConstantZonotopeCover{D}, loc::ModalLocationData{D}, c::AbstractVector{Float64}) where D
    grid = ctx.grid
    for i in 1:D
        (c[i] - cov.radius[i] < grid.lower[i] || c[i] + cov.radius[i] > grid.upper[i]) && return false
    end
    cells = get_touching_cell_idxs_cover(grid, cov, c)
    return !isempty(cells) && all(cell -> loc.fast[cell], cells)
end

# The cells of the layer whose centre lies in the cover translated to ctx.center, within the exact clock intervals.
function append_cover!(res::Vector{CartesianIndex{N}}, ctx::ModalContext{D,N}, cov::ConstantZonotopeCover{D}, layer::Int, lo::Vector{Float64}, hi::Vector{Float64}) where {D,N}
    for cell in get_touching_cell_idxs_cover(ctx.grid, cov, ctx.center, clock_cell_box(ctx.grid, lo, hi)...)
        if cell_in_clock_intervals(ctx.grid, cell, lo, hi)
            push!(res, CartesianIndex(cell, layer))
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

# Sorts v and drops repeated entries, without hashing.
function sorted_unique!(v::Vector)
    sort!(v)
    j = 0
    for i in eachindex(v)
        (j == 0 || v[i] != v[j]) && (j += 1; v[j] = v[i])
    end
    return resize!(v, j)
end

cell_zonotope(ctx::ModalContext{D}, c::CartesianIndex{D}) where D = Zonotope(ctx.grid.lower .+ (Tuple(c) .- 0.5) .* ctx.grid.granularity, ctx.cell_generators)

# Flow successors of a cell, or nothing if flowing from it is bad.
function flow_successors(ctx::ModalContext{D,N}, layer::Int, c::CartesianIndex{D}) where {D,N}
    res = CartesianIndex{N}[]
    out = ModalEndSet[]
    Z = cell_zonotope(ctx, c)
    for li in cell_locations(ctx, layer, c)
        loc = ctx.locations[li]
        lo, hi = clock_end(loc.clock_rate, Z, ctx.steps, ctx.δ⁻)
        if loc.fast[c]
            append_fast!(res, ctx, li, c, lo, hi)
        else
            empty!(out)
            modal_reach!(out, ctx, li, Z, ctx.steps, 0, ctx.exit_bad_for_flow, lo, hi) || return nothing
            append_touching!(res, ctx, out)
        end
    end
    return sorted_unique!(res)
end

# Successors of taking the action at position ai of ctx.action_ids in a cell, or nothing if the action is not available, is bad, or does nothing.
# Trying an action outside its guard, or in a location that does not offer it, does nothing: the state just flows. So:
#   - a cell where it does nothing in every location covering it gets nothing, and its choice for the action is the same as waiting (see noop_choices!)
#   - a cell partly outside the guard (or covered by a location where it does nothing) reaches the action successors and the flow successors `flow`
function action_successors(ctx::ModalContext{D,N}, layer::Int, c::CartesianIndex{D}, ai::Int, exit_bad::Bool, flow) where {D,N}
    can_act(ctx, layer, c, ai) || return nothing
    res = CartesianIndex{N}[]
    out = ModalEndSet[]
    Z = cell_zonotope(ctx, c)
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
        for (k, a) in enumerate(ctx.system.locations[li].actions)
            a.action.id == ctx.action_ids[ai] || continue
            jumped = jump_set(Z, a.action.jumpMatrix, a.action.jumpVector, ctx.maxOrder, ctx.reduceOrder)
            target = ctx.locations[a.target]
            lo, hi = clock_end(target.clock_rate, jumped, ctx.steps, ctx.δ⁻)
            # A jumped cell that lies in cells that are all fast in the target touches nothing for the period either: the
            # precomputed cover
            # (Fast cells may leave the state space when leaving is allowed for the flow: then only if it is for actions too.)
            if (ctx.exit_bad_for_flow || !exit_bad) && lands_in_fast_cells(ctx, ctx.locations[li].act_start_cover[k], target, LazySets.center(jumped))
                mul!(ctx.center, target.dynamics.Φ[ctx.steps+1], LazySets.center(jumped))
                ctx.center .+= LazySets.center(target.dynamics.V[ctx.steps+1])
                append_cover!(res, ctx, ctx.locations[li].act_cover[k], target.layer, lo, hi)
                continue
            end
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
    return sorted_unique!(res)
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

# ---------------------------------------------------------------------------------------------------------------------
# Online: blocks reached forward
#
# Online the shield computes only the cells whose choices can differ from "everything allowed", and gives the same result
# as offline. A cell can only lose a choice if its own computation is bad (it touches an unsafe set, leaves the state
# space where that is bad, or jumps too often) or if one of its successors becomes bad. Instead of searching backwards for
# the cells that can reach a cell, every block of sweep cells per dimension is computed once forward, as one set
# (block_successors): what any of its cells reaches lies in what the block reaches. So the cells to compute are those of
# the blocks whose own computation is bad, and, whenever a cell becomes bad, those of the blocks that reach it. Only
# forward reachability is used, through guards, jumps and actions, as for the cells themselves.

# The blocks: block b covers the cells block_cells(ctx, b); the last block on an axis may be smaller.
block_grid(ctx::ModalContext{D}) where D = CartesianIndices(ntuple(i -> cld(ctx.grid.numCells[i], ctx.sweep[i]), D))
block_cells(ctx::ModalContext{D}, b::CartesianIndex{D}) where D =
    CartesianIndices(ntuple(i -> (1+(b[i]-1)*ctx.sweep[i]):min(b[i] * ctx.sweep[i], ctx.grid.numCells[i]), D))
block_of(ctx::ModalContext{D}, c::CartesianIndex{D}) where D = CartesianIndex(ntuple(i -> cld(c[i], ctx.sweep[i]), D))

# What the box of cells `cells` of a layer reaches in one period, per choice, flowed as one zonotope with the same
# shortcuts as a cell (a box whose cells are all fast, actions landing in fast cells). It contains what each of its
# cells reaches: flow is nothing if waiting is bad for the box; per action position, enabled says whether some location of
# a cell can take it there, and acts its successors (with the waiting ones where some state of the box cannot take it),
# nothing if it is bad.
# The covers of a box of shape G (its generators) in location li, made on first use: of its end set when it is fast,
# and per outgoing action of the jumped box at the start and the end of the period (see fast_cover and act_cover).
function box_covers!(ctx::ModalContext{D}, li::Int, G::Matrix{Float64}) where D
    loc = ctx.locations[li]
    get!(loc.box_covers, G) do
        bloat = Matrix((ctx.grid.granularity / 2) * I, D, D)
        cover(H) = ConstantZonotopeCover(hcat(H, bloat), ntuple(i -> sum(abs, view(H, i, :)), D), Val(D))
        dyn = loc.dynamics
        starts, ends = ConstantZonotopeCover{D}[], ConstantZonotopeCover{D}[]
        for a in ctx.system.locations[li].actions
            tdyn = ctx.locations[a.target].dynamics
            G_J = Matrix{Float64}(genmat(jump_set(Zonotope(zeros(D), G), Matrix{Float64}(a.action.jumpMatrix), zeros(D), ctx.maxOrder, ctx.reduceOrder)))
            push!(starts, cover(G_J))
            push!(ends, cover(hcat(genmat(tdyn.V[ctx.steps+1]), tdyn.Φ[ctx.steps+1] * G_J)))
        end
        (cover(hcat(genmat(dyn.V[ctx.steps+1]), dyn.Φ[ctx.steps+1] * G)), starts, ends)
    end
end

# Makes the precomputed discretization of shape G in location li (see CellDiscretization), if it has none, so that
# modal_reach! of a whole box over a whole period uses it.
function ensure_discretization!(ctx::ModalContext{D}, li::Int, G::Matrix{Float64}) where D
    loc = ctx.locations[li]
    shape_discretization(loc, G) === nothing || return
    dyn = loc.dynamics
    reach = (; ctx.δ⁻, δ⁺=min(ctx.δ⁺, ctx.δ⁻ * 2^floor(Int, log2(ctx.steps))), dyn.A, dyn.P2A_abs, dyn.phiDict, dyn.inputDict, U=ctx.system.locations[li].U)
    push!(loc.discs, (G, CellDiscretization(Zonotope(zeros(D), G), reach, Val(D))))
    return
end

struct BoxSuccessors{N}
    flow::Union{Nothing,Vector{CartesianIndex{N}}}
    enabled::Vector{Bool}
    acts::Vector{Union{Nothing,Vector{CartesianIndex{N}}}}
end

function box_successors(ctx::ModalContext{D,N}, layer::Int, cells::CartesianIndices{D}) where {D,N}
    grid = ctx.grid
    g = grid.granularity
    lo_b = grid.lower .+ (car2vec(first(cells)) .- 1) .* g
    hi_b = grid.lower .+ car2vec(last(cells)) .* g
    G = Matrix{Float64}(genmat(remove_zero_generators(Zonotope(zeros(D), diagm(map((n, r) -> n == 1 ? 0.0 : g * (length(r) / 2), grid.numCells, cells.indices))))))
    Z = Zonotope((lo_b .+ hi_b) ./ 2, G)
    mask = UInt64(0)
    for c in cells
        mask |= ctx.start_mask[c, layer]
    end
    locs = LocationSet(ctx.layer_starts[layer], mask)
    out = ModalEndSet[]
    flow = CartesianIndex{N}[]
    flow_bad = false
    for li in locs
        loc = ctx.locations[li]
        lo, hi = clock_end(loc.clock_rate, Z, ctx.steps, ctx.δ⁻)
        if all(c -> loc.fast[c], cells) # every cell touches nothing over the period, so neither does the box
            mul!(ctx.center, loc.dynamics.Φ[ctx.steps+1], LazySets.center(Z))
            ctx.center .+= LazySets.center(loc.dynamics.V[ctx.steps+1])
            append_cover!(flow, ctx, box_covers!(ctx, li, G)[1], loc.layer, lo, hi)
        else
            ensure_discretization!(ctx, li, G)
            empty!(out)
            modal_reach!(out, ctx, li, Z, ctx.steps, 0, ctx.exit_bad_for_flow, lo, hi) ? append_touching!(flow, ctx, out) : (flow_bad = true)
        end
    end
    n_acts = length(ctx.action_ids)
    enabled = falses(n_acts)
    acts = Union{Nothing,Vector{CartesianIndex{N}}}[CartesianIndex{N}[] for _ in 1:n_acts]
    also_flows = falses(n_acts)
    for li in locs, ai in 1:n_acts
        any(c -> ctx.may_act[li][ai][c], cells) || (also_flows[ai] = true; continue)
        enabled[ai] = true
        all(c -> ctx.must_act[li][ai][c], cells) || (also_flows[ai] = true)
        acts[ai] === nothing && continue
        loc = ctx.locations[li]
        for (k, a) in enumerate(ctx.system.locations[li].actions)
            a.action.id == ctx.action_ids[ai] || continue
            jumped = jump_set(Z, a.action.jumpMatrix, a.action.jumpVector, ctx.maxOrder, ctx.reduceOrder)
            target = ctx.locations[a.target]
            lo, hi = clock_end(target.clock_rate, jumped, ctx.steps, ctx.δ⁻)
            covers = box_covers!(ctx, li, G)
            if (ctx.exit_bad_for_flow || !ctx.exit_bad_for_actions) && lands_in_fast_cells(ctx, covers[2][k], target, LazySets.center(jumped))
                mul!(ctx.center, target.dynamics.Φ[ctx.steps+1], LazySets.center(jumped))
                ctx.center .+= LazySets.center(target.dynamics.V[ctx.steps+1])
                append_cover!(acts[ai], ctx, covers[3][k], target.layer, lo, hi)
            else
                ensure_discretization!(ctx, a.target, Matrix{Float64}(genmat(jumped)))
                empty!(out)
                modal_reach!(out, ctx, a.target, jumped, ctx.steps, 0, ctx.exit_bad_for_actions, lo, hi) ? append_touching!(acts[ai], ctx, out) : (acts[ai] = nothing; break)
            end
        end
    end
    flow = flow_bad ? nothing : sorted_unique!(flow)
    for ai in 1:n_acts
        (enabled[ai] && acts[ai] !== nothing) || continue
        if also_flows[ai]
            flow === nothing ? (acts[ai] = nothing) : append!(acts[ai], flow)
        end
        acts[ai] === nothing || sorted_unique!(acts[ai])
    end
    return BoxSuccessors{N}(flow, collect(enabled), acts)
end

# Every cell that a cell of the block reaches by any choice, and whether the block's own computation is bad.
function block_successors(ctx::ModalContext{D,N}, layer::Int, cells::CartesianIndices{D}) where {D,N}
    r = box_successors(ctx, layer, cells)
    bad = r.flow === nothing || any(ai -> r.enabled[ai] && r.acts[ai] === nothing, eachindex(r.enabled))
    S = r.flow === nothing ? CartesianIndex{N}[] : copy(r.flow)
    for ai in eachindex(r.acts)
        r.enabled[ai] && r.acts[ai] !== nothing && append!(S, r.acts[ai])
    end
    return sorted_unique!(S), bad
end

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
    for b in 1:length(succ_off)-1, i in succ_off[b]:(succ_off[b+1]-1)
        v = succ_pool[i]
        rev_pool[pos[v]] = b
        pos[v] += 1
    end
    return rev_off, rev_pool
end

# Shield for a ModalHybridSystem. Computes every cell of every layer, or with online = true only the cells whose choices
# can differ from everything allowed (see block_successors); both give the same result. Then pushes the bad events backwards.
# Leaving the state space makes waiting bad, and an action bad, unless exit_bad_for_flow or exit_bad_for_actions is false:
# then only the successors inside the grid count, and a choice whose successors all lie outside it is safe.
# Cells are tested for being fast in blocks of sweep_factor cells per dimension (an integer, or one per dimension); online
# also computes these blocks forward. It only changes the time: a fast cell has the same successors as when flowed on its
# own. On cruise control 20 finds no fast cells, 5 and 2 are the fastest, and 1 is slower again because every cell is then
# checked on its own. Every block is also computed forward as one set, and each cell keeps only the successors its block
# reaches too (see block_successors); online uses the blocks to find the cells to compute.
# Returns grid, allowed, lost, info like every shield (see choice_bits in flat.jl). grid is the continuous
# grid; allowed and lost are indexed by (continuous cell..., layer), and a cell is dead in a layer where allowed is 0.
# Trying an action outside its guard or in a location that does not offer it does nothing, so there it is allowed exactly when waiting is (see action_successors).
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
    # The shield works on ids: the linear index of (continuous cell..., layer). See flat.jl.
    n = prod(shape)
    n < typemax(Int32) || throw(ArgumentError("$n cells do not fit in Int32 ids"))
    cells = CartesianIndices(shape) # id -> cell
    ids = LinearIndices(shape) # cell -> id
    store = SuccessorStore(n, [falses(n) for _ in 1:n_acts])
    store.flow_empty_ok = !ctx.exit_bad_for_flow
    store.act_empty_ok = !ctx.exit_bad_for_actions
    s = reset!(FlatShieldState(n_acts, act_needs_safe_flow), n)
    buf = Int32[]
    queued = falses(shape)
    frontier = CartesianIndex{N}[]
    next_frontier = CartesianIndex{N}[]

    for H in ctx.system.unsafe, c in get_touching_cell_idxs_l(grid, H), layer in 1:ctx.n_layers
        k = CartesianIndex(c, layer)
        seed_constraint!(s, Int32(ids[k]))
        store.flow_none[ids[k]] = true
        queued[k] = true
    end
    lap!(timer, "Setup")

    # The successors of every block (see block_successors), as block -> sorted cell ids, and online also inverted as
    # cell -> blocks. A cell's successors are intersected with its block's: the block's are a sound over-approximation of
    # everything its cells reach, so what a cell's own over-approximation adds outside them is not reachable. That keeps the
    # cells' successors inside their blocks' by construction, which online relies on. A bad block's successors are
    # incomplete, so its cells keep their own (they are always computed).
    bgrid = block_grid(ctx)
    bids = LinearIndices((size(bgrid)..., ctx.n_layers))
    nb = length(bids)
    succ_off, succ_pool = ones(Int, nb + 1), Int32[]
    block_bad = falses(nb)
    for layer in 1:ctx.n_layers, b in bgrid # in increasing bid
        S, bad = block_successors(ctx, layer, block_cells(ctx, b))
        block_bad[bids[b, layer]] = bad
        append!(succ_pool, sort!([Int32(ids[v]) for v in S]))
        succ_off[bids[b, layer]+1] = length(succ_pool) + 1
    end
    # Likewise a cell's own computation can only be bad where its block's is: in a good block a bad flow or action of a cell
    # comes from its over-approximation, and the cell takes the block's successors for it instead (a superset of what it reaches).
    block_id(k::CartesianIndex{N}) = bids[block_of(ctx, CartesianIndex(ntuple(i -> k[i], D))), k[N]]
    block_range(bid) = view(succ_pool, succ_off[bid]:(succ_off[bid+1]-1))
    within_block!(buf::Vector{Int32}, bid) = block_bad[bid] ? buf : filter!(v -> insorted(v, block_range(bid)), buf)
    lap!(timer, "Blocks")

    # Online: compute the cells of bad blocks, and whenever a cell becomes bad, those of the blocks that reach it.
    rev_off, rev_pool = online ? reverse_csr(succ_off, succ_pool, n) : (ones(Int, n + 1), Int32[])
    block_queued = falses(nb)
    function queue_block!(out, bid)
        block_queued[bid] && return
        block_queued[bid] = true
        bl = CartesianIndices(bids)[bid]
        layer = bl[N]
        for c in block_cells(ctx, CartesianIndex(ntuple(i -> bl[i], D)))
            k = CartesianIndex(c, layer)
            queued[k] || (queued[k] = true; push!(out, k))
        end
    end
    function queue_blocks_reaching!(out, ks)
        for k in ks, i in rev_off[k]:(rev_off[k+1]-1)
            queue_block!(out, rev_pool[i])
        end
    end
    if online
        for bid in findall(block_bad)
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
            bid = block_id(k)
            flow = flow_successors(ctx, layer, c)
            flow === nothing && !block_bad[bid] && (flow = [cells[v] for v in block_range(bid)])
            store.flow_none[id] = flow === nothing
            set_successors!(store, 0, id, flow === nothing ? empty!(buf) : within_block!(map!(v -> Int32(ids[v]), resize!(buf, length(flow)), flow), bid))
            for ai in 1:n_acts
                succ = action_successors(ctx, layer, c, ai, ctx.exit_bad_for_actions, flow)
                if succ === nothing && !block_bad[bid] && can_act(ctx, layer, c, ai)
                    succ = [cells[v] for v in block_range(bid)] # bad here, but not for the block
                end
                store.act_enabled[ai][id] = succ !== nothing
                set_successors!(store, ai, id, succ === nothing ? empty!(buf) : within_block!(map!(v -> Int32(ids[v]), resize!(buf, length(succ)), succ), bid))
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
    verbose && report(timer, "Rounds: $rounds, cells computed: $count_cells, dead: $(count(s.dead)), block successors: $(length(succ_pool))")
    enabled(ai, k) = can_act(ctx, cells[k][N], CartesianIndex(ntuple(i -> cells[k][i], D)), ai)
    allowed, lost = cell_choices(shape, n_acts, s, enabled, store_acted(store))
    noop_choices!(allowed, ctx)
    return grid, allowed, lost, (; ctx)
end
