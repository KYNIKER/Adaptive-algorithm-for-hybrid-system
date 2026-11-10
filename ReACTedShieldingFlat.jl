using LazySets, LinearAlgebra
# Included at the end of ReACTedShielding.jl, after ReACTedPropagate.jl. Expects PartitionUtilities.jl to be included first.
#
# Building blocks shared by the fast shields (ReACTedShieldingAdaptiveFast.jl, ReACTedShieldingPredFast.jl). Both give
# every cell (pred) or leaf (adaptive) an Int32 id and keep what they know per id in flat arrays.
#   - CellReach: the reach set of each cell (or block of cells) from one zonotope translated to it
#   - reach_covers: the covers of what fast cells and actions reach, each a translate of one zonotope
#   - touching_masks: the cells touching a polyhedron, from dot products against its halfspaces
#   - SuccessorStore: the successor lists of all ids in one Int32 pool
#   - FlatShieldState: the event-driven shield over ids, with the rules of OnlinePredShieldState (ReACTedShieldingPred.jl)
#   - choice_bits and its helpers: the common results of all shields, what each cell allows and what it lost, as bits

# ---------------------------------------------------------------------------------------------------------------------
# Reach sets of cells

# One zonotope Z the size of a cell, or of a block of scale cells per dimension, with its discretization. Each
# computation translates it to one cell and back. Zero-width dimensions get no generator.
struct CellReach{N,GD,DD}
    Z::DenseZonotope
    generatorDict::GD
    discretizationDict::DD
    pc::Vector{Float64} # buffers of ReACT_discretize_combine_with_offset_vector!
    pg::Matrix{Float64}
    pig::Matrix{Float64}
    lower_offset::Vector{Float64} # centre of the first cell
    g::Float64
    scale::Vector{Int}
end

function CellReach(reach, grid::Grid{N}, scale::Vector{Int}=ones(Int, N)) where N
    (; δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict) = reach
    g = grid.granularity
    Z = remove_zero_generators(Zonotope(zeros(N), diagm(map((n, s) -> n == 1 ? 0.0 : g * (s / 2), grid.numCells, scale))))
    generatorDict = ReACT_discretize_decomposed_generators(Z, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict)
    discretizationDict, finaldims = ReACT_discretize_combine_with_offsets(Z, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict)
    return CellReach{N,typeof(generatorDict),typeof(discretizationDict)}(Z, generatorDict, discretizationDict, zeros(N),
        zeros(N, finaldims), zeros(N, finaldims), grid.lower .+ g .* (scale ./ 2), g, scale)
end

cell_offset(cr::CellReach, c::CartesianIndex) = cr.lower_offset + ((car2vec(c) .- 1) .* cr.g .* cr.scale)

# f(Z, discretization) with Z translated to cell c and discretized there. Z is translated back afterwards.
function with_cell(f, cr::CellReach, reach, c::CartesianIndex)
    of = cell_offset(cr, c)
    LazySets.API.translate!(cr.Z, of)
    ReACT_discretize_combine_with_offset_vector!(cr.pc, cr.pg, cr.pig, cr.discretizationDict, of, reach.δ⁻, reach.δ⁺, reach.A, reach.P2A_abs, reach.phiDict, reach.U, reach.inputDict, cr.generatorDict)
    res = f(cr.Z, cr.discretizationDict)
    LazySets.API.translate!(cr.Z, -of)
    return res
end

# The sets reached from cell c at the end of the period and the flag of propagate_set_jumps (1 = bad).
cell_flow(reach, cr::CellReach, c::CartesianIndex) = with_cell((Z, disc) -> propagate_set_jumps(reach.jump_reach, Z, disc), cr, reach, c)

# Whether block b of cr cannot reach a guard or leave the invariant within the period, which makes its cells fast.
is_fast_block(system, reach, cr::CellReach{N}, b::CartesianIndex{N}) where N = with_cell(cr, reach, b) do Z, disc
    propagate_set_check_hybrid_jump(Z, [0.0, reach.p], reach.δ⁻, reach.δ⁺, system, system.edges, reach.identity_matrix, reach.attempt_recorder, N, reach.P2A_abs, reach.guardProjVectors, reach.guardProjBounds, reach.invariantProjVectors, reach.invariantProjBounds, disc, reach.phiDict, reach.tPhiDict, reach.inputDict, reach.max_input, reach.max_flow, reach.guard_invariant_intersection, reach.alg, reach.maxOrder, reach.reduceOrder) == 0
end

# Centre of acting with act in cell c, in the order LazySets computes it: translate J * (cell centre) by j, map by
# max_flow, add max_input.
act_center!(out, reach, cr::CellReach, act, c::CartesianIndex) = (out .= reach.max_flow * (act.jumpMatrix * cell_offset(cr, c) + act.jumpVector) + reach.max_input.center)

# Covers of what a cell of size Z reaches, bloated by each half-cell size in halves for the cover test on that grid.
# A fast cell reaches max_flow * cell + max_input and acting reaches max_flow * (J * cell + j) + max_input, so both are
# translates of one zonotope per cell size. The action covers are those get_touching_cell_idxs_lz builds: the nonzero
# generators in the same order, with the radius of the unbloated set. Returns the fast covers per half-cell size and the
# action covers per action and half-cell size.
function reach_covers(system, reach, Z, halves::Vector{Float64}, ::Val{N}) where N
    (; max_flow, max_input) = reach
    G = hcat(genmat(max_input), genmat(linear_map(max_flow, Z)))
    fast = map(halves) do h
        Gh = hcat(G, diagm(fill(h, N)))
        ConstantZonotopeCover(Gh, ntuple(i -> sum(abs, view(Gh, i, :)), N), Val(N))
    end
    act = map(system.Act) do a
        Ga = hcat(max_flow * (a.jumpMatrix * genmat(Z)), genmat(max_input))
        Ga = Ga[:, [j for j in axes(Ga, 2) if any(!iszero, view(Ga, :, j))]]
        radius = ntuple(i -> sum(abs, view(Ga, i, :)), N)
        [ConstantZonotopeCover(hcat(Ga, Matrix(h * I, N, N)), radius, Val(N)) for h in halves]
    end
    return fast, act
end

# Appends the ids in ids whose cell centre lies in the cover translated to c (0 = no id). Visits the same cells, in the
# same order, as get_touching_cell_idxs_cover.
function cover_ids!(out::Vector{Int32}, ids::Array{Int32,N}, grid::Grid{N}, cov::ConstantZonotopeCover{N}, c::AbstractVector{Float64}) where N
    lower_bounds = ntuple(i -> max(floor(Int, abs(clamp(c[i] - cov.radius[i], grid.lower[i], grid.upper[i]) - grid.lower[i]) / grid.granularity), 1), N)
    upper_bounds = ntuple(i -> min(grid.numCells[i], ceil(Int, abs(clamp(c[i] + cov.radius[i], grid.lower[i], grid.upper[i]) - grid.lower[i]) / grid.granularity)), N)
    sweeping_dim = argmax(upper_bounds[i] - lower_bounds[i] for i in 1:N)
    next_sweep = argmax(i == sweeping_dim ? 0 : upper_bounds[i] - lower_bounds[i] for i in 1:N)
    emit(idx) = (id = ids[idx]; id != 0 && push!(out, id))
    if upper_bounds[next_sweep] == lower_bounds[next_sweep]
        foreach(emit, CartesianIndices(ntuple(i -> lower_bounds[i]:upper_bounds[i], N)))
    elseif upper_bounds[sweeping_dim] - lower_bounds[sweeping_dim] > 1
        step = CartesianIndex(ntuple(i -> i == sweeping_dim ? 1 : 0, N))
        for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], lower_bounds[i])), N))
            first = findfirst(off -> cell_center_in_cover(cov, grid, idx + step * off, c), lower_bounds[sweeping_dim]:upper_bounds[sweeping_dim])
            first === nothing && continue
            lo = lower_bounds[sweeping_dim] + first - 1
            hi = upper_bounds[sweeping_dim]
            while !cell_center_in_cover(cov, grid, idx + step * hi, c)
                hi -= 1
            end
            for off in lo:hi
                emit(idx + step * off)
            end
        end
    else
        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:upper_bounds[i], N))
            cell_center_in_cover(cov, grid, idx, c) && emit(idx)
        end
    end
    return out
end

# ---------------------------------------------------------------------------------------------------------------------
# Cells touching a polyhedron

# Halfspaces of a polyhedron as rows, for allocation-free cell tests.
struct CellHalfSpaces{N}
    P::HPolyhedron{Float64,Vector{Float64}}
    a::Matrix{Float64}
    b::Vector{Float64}
end

function CellHalfSpaces(P, ::Val{N}) where N
    cs = constraints_list(P)
    return CellHalfSpaces{N}(HPolyhedron([LazySets.HalfSpace(Vector{Float64}(h.a), Float64(h.b)) for h in cs]),
        Matrix{Float64}(reduce(vcat, (Vector{Float64}(h.a)' for h in cs))), Float64[h.b for h in cs])
end

# How the closed box of cell c meets the polyhedron: 0 misses, 1 touches, 2 inside. A box cut by exactly one halfspace
# and inside the others meets their intersection; a box cut by several is decided by an LP.
# Cell borders are sums like lower + (c - 1) * g, which rounding can move off a face they lie on, so a box within
# CELL_MEETS_TOL (relative) of a halfspace still touches it, as for the LP.
const CELL_MEETS_TOL = 1e-9
function cell_meets(H::CellHalfSpaces{N}, grid::Grid{N}, c::CartesianIndex{N}) where N
    cut = 0
    g = grid.granularity
    @inbounds for k in eachindex(H.b)
        lo = 0.0
        hi = 0.0
        for i in 1:N
            a = H.a[k, i]
            a == 0 && continue
            l = grid.lower[i] + (c[i] - 1) * g
            u = l + g
            lo += a > 0 ? a * l : a * u
            hi += a > 0 ? a * u : a * l
        end
        tol = CELL_MEETS_TOL * max(1.0, abs(H.b[k]), abs(lo), abs(hi))
        lo > H.b[k] + tol && return 0
        hi > H.b[k] + tol && (cut += 1)
    end
    cut == 0 && return 2
    cut == 1 && return 1
    box = Hyperrectangle(grid.lower .+ (car2vec(c) .- 0.5) .* g, fill(g / 2, N))
    return isdisjoint(H.P, box) ? 0 : 1
end

# Masks of the cells touching P on every level, the same cells as get_touching_cell_idxs_l per level. A finer level only
# tests the children of touching cells that are not inside P (child_cells is in ReACTedShieldingAdaptive.jl).
function touching_masks(grids::Vector{Grid{N}}, P, factor::NTuple{N,Int}) where N
    H = CellHalfSpaces(P, Val(N))
    # The candidate cells of get_touching_cell_idxs_l on a grid. Its bounds are clamped to the state space, so on a grid
    # that overhangs the state space they leave out cells there; each level's candidates lie in its parents' candidates.
    function candidates(grid)
        lo = clamp.(low(P), grid.lower, grid.upper)
        hi = clamp.(high(P), grid.lower, grid.upper)
        l = Int.(max.(floor.(abs.(lo .- grid.lower) ./ grid.granularity), 1))
        u = min.(grid.numCells, Int.(ceil.(abs.(hi .- grid.lower) ./ grid.granularity)))
        return CartesianIndices(ntuple(i -> l[i]:max(u[i], l[i]), N))
    end
    grid0 = grids[1]
    touching = falses(cellshape(grid0))
    inside = falses(cellshape(grid0))
    for c in candidates(grid0)
        k = cell_meets(H, grid0, c)
        touching[c] = k > 0
        inside[c] = k == 2
    end
    masks = [touching]
    for grid in grids[2:end]
        touching′ = falses(cellshape(grid))
        inside′ = falses(cellshape(grid))
        range = candidates(grid)
        for c in CartesianIndices(touching)
            touching[c] || continue
            for cc in child_cells(c, factor)
                cc in range || continue
                k = inside[c] ? 2 : cell_meets(H, grid, cc)
                touching′[cc] = k > 0
                inside′[cc] = k == 2
            end
        end
        push!(masks, touching′)
        touching, inside = touching′, inside′
    end
    return masks
end

# ---------------------------------------------------------------------------------------------------------------------
# Successor lists

# Successor lists of all ids in one Int32 pool, per choice ai: 0 is the flow (waiting a period), ai > 0 the action
# Act[ai], as in the bits of choice_bits. A recomputed list is appended and its old slice becomes garbage, which
# compact! removes.
mutable struct SuccessorStore
    pool::Vector{Int32}
    live::Int # pool entries still referenced
    flow_none::BitVector # the flow is bad (it hit a constraint or left the state space), so it has no successors
    flow_start::Vector{Int}
    flow_len::Vector{Int32}
    act_enabled::Vector{BitVector} # per action: the id touches its guard
    act_start::Vector{Vector{Int}}
    act_len::Vector{Vector{Int32}}
end

SuccessorStore(n::Int, act_enabled::Vector{BitVector}) = SuccessorStore(Int32[], 0, falses(n), ones(Int, n), zeros(Int32, n),
    act_enabled, [ones(Int, n) for _ in act_enabled], [zeros(Int32, n) for _ in act_enabled])

# Adds an id without successors, enabled for the actions whose mask is set at cell.
function push_id!(store::SuccessorStore, may_act::Vector{<:AbstractArray{Bool}}, cell::CartesianIndex)
    push!(store.flow_none, false)
    push!(store.flow_start, 1)
    push!(store.flow_len, 0)
    for ai in eachindex(store.act_enabled)
        push!(store.act_enabled[ai], may_act[ai][cell])
        push!(store.act_start[ai], 1)
        push!(store.act_len[ai], 0)
    end
end

successor_range(store::SuccessorStore, ai::Int, k::Integer) =
    ai == 0 ? (store.flow_start[k]:store.flow_start[k]+store.flow_len[k]-1) : (store.act_start[ai][k]:store.act_start[ai][k]+store.act_len[ai][k]-1)

# Replaces the successors of id for choice ai.
function set_successors!(store::SuccessorStore, ai::Int, id::Integer, succ::AbstractVector{Int32})
    start, len = ai == 0 ? (store.flow_start, store.flow_len) : (store.act_start[ai], store.act_len[ai])
    store.live += length(succ) - len[id]
    start[id] = length(store.pool) + 1
    len[id] = length(succ)
    append!(store.pool, succ)
end

# Rewrites the pool without garbage once garbage is most of it.
function compact!(store::SuccessorStore)
    length(store.pool) <= 2 * store.live + 2^20 && return
    pool = Vector{Int32}(undef, store.live)
    i = 1
    for (start, len) in Iterators.flatten(((store.flow_start, store.flow_len), zip(store.act_start, store.act_len)))
        for k in eachindex(start)
            copyto!(pool, i, store.pool, start[k], len[k])
            start[k] = i
            i += len[k]
        end
    end
    store.pool = pool
end

# ---------------------------------------------------------------------------------------------------------------------
# Event-driven shield over ids with the rules of OnlinePredShieldState (ReACTedShieldingPred.jl). Predecessor lists are
# linked lists in flat arrays: the predecessors of v start at edge head[v] and follow edge_next (0 ends).
mutable struct FlatShieldState
    dead::BitVector
    flow_bad::BitVector
    computed::BitVector
    act_ok::Vector{BitVector}
    n_ok::Vector{Int32}
    flow_head::Vector{Int32}
    act_head::Vector{Vector{Int32}}
    edge_from::Vector{Int32}
    edge_next::Vector{Int32}
    dead_queue::Vector{Int32}
    flow_bad_queue::Vector{Int32}
    changed::Vector{Int32}
    act_needs_safe_flow::Bool
end

FlatShieldState(n_acts::Int, act_needs_safe_flow::Bool) = FlatShieldState(BitVector(), BitVector(), BitVector(), [BitVector() for _ in 1:n_acts],
    Int32[], Int32[], [Int32[] for _ in 1:n_acts], Int32[], Int32[], Int32[], Int32[], Int32[], act_needs_safe_flow)

# Clears the state for n ids, keeping the allocated memory.
function reset!(s::FlatShieldState, n::Int)
    for x in (s.dead, s.flow_bad, s.computed, s.act_ok...)
        resize!(x, n)
        fill!(x, false)
    end
    for x in (s.n_ok, s.flow_head, s.act_head...)
        resize!(x, n)
        fill!(x, 0)
    end
    foreach(empty!, (s.edge_from, s.edge_next, s.dead_queue, s.flow_bad_queue, s.changed))
    return s
end

@inline function add_pred!(s::FlatShieldState, head::Vector{Int32}, v::Int32, k::Int32)
    push!(s.edge_from, k)
    push!(s.edge_next, head[v])
    head[v] = length(s.edge_from)
end

@inline function check_dead!(s::FlatShieldState, k::Int32)
    if !s.dead[k] && s.computed[k] && s.flow_bad[k] && s.n_ok[k] == 0
        s.dead[k] = true
        push!(s.dead_queue, k)
        push!(s.changed, k)
    end
end

@inline function mark_flow_bad!(s::FlatShieldState, k::Int32)
    if !s.flow_bad[k]
        s.flow_bad[k] = true
        push!(s.flow_bad_queue, k)
        push!(s.changed, k)
        check_dead!(s, k)
    end
end

# Makes every action that leads into v unsafe.
function invalidate_actions_into!(s::FlatShieldState, v::Int32)
    for ai in eachindex(s.act_head)
        e = s.act_head[ai][v]
        while e != 0
            k = s.edge_from[e]
            if s.act_ok[ai][k]
                s.act_ok[ai][k] = false
                s.n_ok[k] -= 1
                push!(s.changed, k)
                check_dead!(s, k)
            end
            e = s.edge_next[e]
        end
    end
end

# Constraint ids start dead without being computed.
function seed_constraint!(s::FlatShieldState, k::Int32)
    s.dead[k] = true
    s.flow_bad[k] = true
    s.computed[k] = true
end

# register_cell! for a computed id, reading its successors from the store.
function register_successors!(store::SuccessorStore, s::FlatShieldState, k::Int32)
    pool = store.pool
    if store.flow_none[k] || store.flow_len[k] == 0
        s.flow_bad[k] = true
    else
        for i in successor_range(store, 0, k)
            v = pool[i]
            add_pred!(s, s.flow_head, v, k)
            s.dead[v] && (s.flow_bad[k] = true)
        end
    end
    n = Int32(0)
    for ai in eachindex(s.act_ok)
        range = successor_range(store, ai, k)
        (store.act_enabled[ai][k] && !isempty(range)) || continue
        any(i -> s.dead[pool[i]] || (s.act_needs_safe_flow && s.flow_bad[pool[i]]), range) && continue
        for i in range
            add_pred!(s, s.act_head[ai], pool[i], k)
        end
        s.act_ok[ai][k] = true
        n += Int32(1)
    end
    s.n_ok[k] = n
    s.computed[k] = true
    if s.flow_bad[k]
        push!(s.flow_bad_queue, k)
        push!(s.changed, k)
    end
    check_dead!(s, k)
end

# propagate_bad! over the linked predecessor lists.
function propagate_bad!(s::FlatShieldState)
    while !isempty(s.dead_queue) || !isempty(s.flow_bad_queue)
        while !isempty(s.dead_queue)
            v = pop!(s.dead_queue)
            e = s.flow_head[v]
            while e != 0
                mark_flow_bad!(s, s.edge_from[e])
                e = s.edge_next[e]
            end
            invalidate_actions_into!(s, v)
        end
        while !isempty(s.flow_bad_queue)
            v = pop!(s.flow_bad_queue)
            (s.act_needs_safe_flow && !s.dead[v]) && invalidate_actions_into!(s, v)
        end
    end
end

# ---------------------------------------------------------------------------------------------------------------------
# Results. Every shield returns grid, allowed, lost, info: allowed and lost hold the choice bits of each cell (see
# choice_bits), grid.deadCells is where nothing is allowed, and info is a NamedTuple of what only that shield has.

# The smallest unsigned type with a bit for waiting and one for each of n_acts actions.
function choice_type(n_acts::Int)
    n_acts < 64 || throw(ArgumentError("at most 63 actions fit in the allowed matrix"))
    return n_acts < 8 ? UInt8 : n_acts < 16 ? UInt16 : n_acts < 32 ? UInt32 : UInt64
end

# What k (a cell or id) allows and what it lost, as bits: bit 0 (value 1) is waiting a period, bit ai (value 2^ai) is
# the action at position ai; nothing allowed means dead. s is any shield state with the fields dead, flow_bad, computed
# and act_ok. A k never computed allows waiting and every action it is enabled for (it touches the guard). An action is
# lost in an alive k that is enabled for it and acted (has successors for it), but where it is unsafe.
# enabled(ai, k) and acted(ai, k) say whether k is enabled for and has successors for action ai.
function choice_bits(::Type{T}, s, k, enabled::E, acted::A) where {T,E,A}
    s.dead[k] && return zero(T), zero(T)
    allowed = one(T)
    lost = zero(T)
    if !s.computed[k]
        for ai in eachindex(s.act_ok)
            enabled(ai, k) && (allowed |= T(1) << ai)
        end
    else
        s.flow_bad[k] && (allowed = zero(T))
        for ai in eachindex(s.act_ok)
            if s.act_ok[ai][k]
                allowed |= T(1) << ai
            elseif enabled(ai, k) && acted(ai, k)
                lost |= T(1) << ai
            end
        end
    end
    return allowed, lost
end

# enabled and acted for choice_bits from a SuccessorStore.
store_enabled(store::SuccessorStore) = (ai, k) -> store.act_enabled[ai][k]
store_acted(store::SuccessorStore) = (ai, k) -> store.act_len[ai][k] > 0

# allowed and lost of a shield whose state is indexed like a cell array of the given shape (by cell or linear index).
function cell_choices(shape, n_acts::Int, s, enabled, acted)
    T = choice_type(n_acts)
    return cell_choices!(zeros(T, shape), zeros(T, shape), s, enabled, acted)
end

# Function barrier on T
function cell_choices!(allowed::Array{T}, lost::Array{T}, s, enabled::E, acted::A) where {T,E,A}
    for k in eachindex(allowed)
        allowed[k], lost[k] = choice_bits(T, s, k, enabled, acted)
    end
    return allowed, lost
end

# allowed and lost of an adaptive shield on its finest grid: every active leaf (with the fields level, cell and active)
# writes its bits to the finest cells it covers. child_cells is in ReACTedShieldingAdaptive.jl.
function leaf_choices(shape, n_acts::Int, leaves, L::Int, factor, s, enabled, acted)
    T = choice_type(n_acts)
    return leaf_choices!(zeros(T, shape), zeros(T, shape), leaves, L, factor, s, enabled, acted)
end

# Function barrier on T
function leaf_choices!(allowed::Array{T,N}, lost::Array{T,N}, leaves, L::Int, factor::NTuple{N,Int}, s, enabled::E, acted::A) where {T,N,E,A}
    for id in Int32(1):Int32(length(leaves.level))
        leaves.active[id] || continue
        bits, lost_bits = choice_bits(T, s, id, enabled, acted)
        e = L - leaves.level[id]
        for cc in child_cells(leaves.cell[id], ntuple(i -> factor[i]^e, N))
            allowed[cc] = bits
            lost[cc] = lost_bits
        end
    end
    return allowed, lost
end

# For shields that end with successor relations keyed by cell instead of a state with act_ok (the shield_step! shields):
# the state choice_bits reads. An action is ok in an alive computed cell when it has successors and none of them is
# dead or, with act_needs_safe_flow, unable to wait. A cell is computed when it has a relation entry or a bad flow.
function relation_state(dead::BitArray{N}, flow_bad::BitArray{N}, act_rel, flow_rel, action_ids::Vector{Int}; act_needs_safe_flow::Bool=true) where N
    computed = copy(flow_bad)
    for c in keys(flow_rel)
        computed[c] = true
    end
    act_ok = [falses(size(dead)) for _ in action_ids]
    for ((id, c), succ) in act_rel
        computed[c] = true
        ok = !dead[c] && !isempty(succ) && all(v -> !dead[v] && !(act_needs_safe_flow && flow_bad[v]), succ)
        ok && (act_ok[findfirst(==(id), action_ids)][c] = true)
    end
    return (; dead, flow_bad, computed, act_ok)
end

# acted for choice_bits from an action relation keyed by (action id, cell), for linear indices into cells.
function relation_acted(act_rel, action_ids::Vector{Int}, cells::CartesianIndices)
    return function (ai, k)
        succ = get(act_rel, (action_ids[ai], cells[k]), nothing)
        return succ !== nothing && !isempty(succ)
    end
end

# The offline shield of successor relations over all cells, with the rules of FlatShieldState: loads the relations into
# a SuccessorStore with the cells' linear indices as ids, seeds the constraint cells dead and registers every other
# cell. A cell without a flow entry has a bad flow. Returns the propagated state and the store.
function relation_shield(constraint::BitArray{N}, act_rel, flow_rel, action_ids::Vector{Int}; act_needs_safe_flow::Bool=true) where N
    n = length(constraint)
    lin = LinearIndices(constraint)
    store = SuccessorStore(n, [falses(n) for _ in action_ids])
    s = reset!(FlatShieldState(length(action_ids), act_needs_safe_flow), n)
    ids(cells) = Int32[lin[v] for v in cells]
    for (k, c) in zip(Int32(1):Int32(n), CartesianIndices(constraint))
        if constraint[c]
            seed_constraint!(s, k)
            continue
        end
        flow = get(flow_rel, c, nothing)
        flow === nothing ? (store.flow_none[k] = true) : set_successors!(store, 0, k, ids(flow))
        for (ai, id) in enumerate(action_ids)
            succ = get(act_rel, (id, c), nothing)
            succ === nothing && continue
            store.act_enabled[ai][k] = true
            set_successors!(store, ai, k, ids(succ))
        end
        register_successors!(store, s, k)
    end
    propagate_bad!(s)
    return s, store
end
