using LazySets, LinearAlgebra, ReachabilityAnalysis, Dates
# Expects PartitionUtilities.jl and ReACTedShielding.jl to be included first.
#
# Adaptive (multi-resolution) online shield synthesis. The state space is partitioned into leaves of a cell tree: level
# 0 is the grid at max_granularity, and a leaf at level l can be split into the cells covering it at level l+1
# (granularity halved; zero-width axes are never split). The shield over the active leaves has the rules of
# ReACTed_reachable_cell_online_pred:
#   flow is bad   : no flow successors (hit a constraint), none at all, or a flow successor is dead
#   action is ok  : the leaf touches the action's guard and every successor is alive (and, if act_needs_safe_flow, can wait)
#   leaf is dead  : flow is bad and no action is ok
# Only leaves reached by the frontier are computed and the others are assumed safe; constraint leaves are dead and never
# computed. Every leaf's successors over-approximate its own reach set, so any mix of levels gives a valid shield and
# splitting only changes permissiveness. With all leaves on one level it computes the same shield as
# ReACTed_reachable_cell_online_pred at that granularity.
#
# How it is kept fast:
#   - Action guard and constraint masks are built top-down: level 0 tests its cells, and a finer level only tests the
#     children of touching cells that are not inside the polyhedron. Cells are tested with dot products against the
#     halfspaces, and only cells cut by two or more halfspaces fall back to the LP of isdisjoint.
#   - Fast-cell covers use the direct facet enumeration of ConstantZonotopeCover(G, radius, Val(N)), and action
#     successors use one precomputed cover per level, action and target level instead of a zonotope per leaf.
#   - Leaves, their successor lists and the shield's predecessor lists are kept in flat arrays instead of a vector or
#     Dict entry per leaf, so the live heap holds few objects and garbage collection stays cheap.
#   - The result is one matrix on the finest grid with a bit per allowed choice instead of the successor relations,
#     returned as grid, allowed, lost, info like every shield (see choice_bits in ReACTedShieldingFlat.jl).

# Per-level data: the reach of its cells, which cells are fast, the covers per target level and the masks.
struct AdaptiveFastLevel{N,CR,F}
    grid::Grid{N}
    cell_reach::CR
    fast::LazyFastCells{N,F}
    fast_covers::Vector{ConstantZonotopeCover{N}} # per target level
    act_covers::Vector{Vector{ConstantZonotopeCover{N}}} # per action, per target level
    may_act::Vector{BitArray{N}} # per action
    constraint::BitArray{N}
end

function AdaptiveFastLevel(system::EuclideanHybridSystem, reach, grid::Grid{N}, grid0::Grid{N}, L::Int, sweep_factor::Int, may_act, constraint) where N
    sweep = [min(sweep_factor, n) for n in grid.numCells]
    nblocks = ntuple(i -> fld(grid.numCells[i], sweep[i]), N)
    block_reach = CellReach(reach, grid, sweep)
    fast = LazyFastCells(cellshape(grid), Tuple(sweep), nblocks, falses(nblocks), falses(nblocks), b -> is_fast_block(system, reach, block_reach, b))
    cell_reach = CellReach(reach, grid)
    fast_covers, act_covers = reach_covers(system, reach, cell_reach.Z, [grid0.granularity / 2^m / 2 for m in 0:L], Val(N))
    return AdaptiveFastLevel(grid, cell_reach, fast, fast_covers, act_covers, may_act, constraint)
end

# ---------------------------------------------------------------------------------------------------------------------
# Leaves in flat storage, with their successor lists in a SuccessorStore indexed by leaf id. Each level maps its cells
# to active leaf ids with a dense array (0 = none), allocated when the level gets its first leaf.
struct FlatLeaves{N}
    level::Vector{Int8}
    cell::Vector{CartesianIndex{N}}
    active::BitVector
    ids::Vector{Array{Int32,N}} # per level
    shapes::Vector{NTuple{N,Int}} # per level
    n_active::Vector{Int} # per level
    succ::SuccessorStore
end

FlatLeaves{N}(shapes::Vector{NTuple{N,Int}}, n_acts::Int) where N = FlatLeaves{N}(Int8[], CartesianIndex{N}[], BitVector(),
    [zeros(Int32, ntuple(_ -> 0, N)) for _ in shapes], shapes, zeros(Int, length(shapes)), SuccessorStore(0, [BitVector() for _ in 1:n_acts]))

function deactivate_leaf!(leaves::FlatLeaves, id::Int32)
    l = leaves.level[id]
    leaves.active[id] = false
    leaves.ids[l+1][leaves.cell[id]] = 0
    leaves.n_active[l+1] -= 1
    for ai in 0:length(leaves.succ.act_enabled)
        set_successors!(leaves.succ, ai, id, Int32[])
    end
end

struct AdaptiveFastContext{N,S,R,LV}
    system::S
    reach::R
    levels::Vector{LV}
    factor::NTuple{N,Int}
    leaves::FlatLeaves{N}
    L::Int
    buffer::Vector{Int32} # successor list being computed
    covers::CoverCache{N}
end

AdaptiveFastContext(system::EuclideanHybridSystem, p, max_granularity::Float64, L::Int, δ⁻::Float64, δ⁺::Float64, alg, maxOrder::Int, reduceOrder::Int, sweep_factor::Int; jump_hull::Bool=false) =
    AdaptiveFastContext(Grid(system.statespace, max_granularity), system, p, L, δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor; jump_hull)

# Function barrier on Grid{N}
function AdaptiveFastContext(grid0::Grid{N}, system, p, L, δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor; jump_hull::Bool=false) where N
    zero_width = low(system.statespace) .== high(system.statespace)
    factor = ntuple(i -> zero_width[i] ? 1 : 2, N)
    reach = reach_setup(system, p, δ⁻, δ⁺, alg, maxOrder, reduceOrder, Val(N); jump_hull)
    grids = [level_grid(grid0, l, factor) for l in 0:L]
    may_act = [touching_masks(grids, act.guard, factor) for act in system.Act] # per action, per level
    constraint = [falses(cellshape(grid)) for grid in grids]
    for gc in system.globalConstraints
        for (l, m) in enumerate(touching_masks(grids, gc, factor))
            constraint[l] .|= m
        end
    end
    levels = [AdaptiveFastLevel(system, reach, grids[l+1], grid0, L, sweep_factor, [m[l+1] for m in may_act], constraint[l+1]) for l in 0:L]
    leaves = FlatLeaves{N}([NTuple{N,Int}(Tuple(cellshape(grid))) for grid in grids], length(system.Act))
    return AdaptiveFastContext(system, reach, levels, factor, leaves, L, Int32[], CoverCache{N}())
end

function add_leaf!(ctx::AdaptiveFastContext{N}, level::Int, cell::CartesianIndex{N}) where N
    leaves = ctx.leaves
    push!(leaves.level, Int8(level))
    push!(leaves.cell, cell)
    push!(leaves.active, true)
    push_id!(leaves.succ, ctx.levels[level+1].may_act, cell)
    id = Int32(length(leaves.level))
    isempty(leaves.ids[level+1]) && (leaves.ids[level+1] = zeros(Int32, leaves.shapes[level+1]))
    leaves.ids[level+1][cell] = id
    leaves.n_active[level+1] += 1
    return id
end

is_constraint_leaf(ctx::AdaptiveFastContext, id::Integer) = ctx.levels[ctx.leaves.level[id]+1].constraint[ctx.leaves.cell[id]]

# Appends the active leaves of every level m whose cell centre lies in covers[m+1] translated to center.
function cover_leaves!(buf::Vector{Int32}, ctx::AdaptiveFastContext, covers, center)
    for m in 0:ctx.L
        ctx.leaves.n_active[m+1] == 0 || cover_ids!(buf, ctx.leaves.ids[m+1], ctx.levels[m+1].grid, covers[m+1], center)
    end
end

# Successors of leaf id against the current active leaves, stored in the leaves' SuccessorStore.
function compute_leaf!(ctx::AdaptiveFastContext, id::Int32)
    (; leaves, reach, levels, buffer) = ctx
    succ = leaves.succ
    lev = levels[leaves.level[id]+1]
    c = leaves.cell[id]

    empty!(buffer)
    if lev.fast[c]
        fast_cell_center!(reach.center, reach.max_flow, reach.max_input.center, lev.grid, c)
        cover_leaves!(buffer, ctx, lev.fast_covers, reach.center)
        succ.flow_none[id] = false
    else
        z, f = cell_flow(reach, lev.cell_reach, c)
        succ.flow_none[id] = f == 1
        if f != 1
            for m in 0:ctx.L
                leaves.n_active[m+1] == 0 && continue
                ids = leaves.ids[m+1]
                for cc in get_touching_cell_idxs_lz(levels[m+1].grid, z, ctx.covers)
                    ids[cc] != 0 && push!(buffer, ids[cc])
                end
            end
        end
    end
    set_successors!(succ, 0, id, buffer)

    for (ai, act) in enumerate(ctx.system.Act)
        empty!(buffer)
        if succ.act_enabled[ai][id]
            act_center!(reach.center, reach, lev.cell_reach, act, c)
            cover_leaves!(buffer, ctx, lev.act_covers[ai], reach.center)
        end
        set_successors!(succ, ai, id, buffer)
    end
    return id
end

# Active leaves next to any of ids (sharing a face or corner), on every level, with the dense cell-to-leaf maps.
function neighbour_leaves(ctx::AdaptiveFastContext{N}, ids::Vector{Int32}) where N
    leaves = ctx.leaves
    levels_used = findall(>(0), leaves.n_active)
    if length(levels_used) == 1
        m = levels_used[1]
        grid = ctx.levels[m].grid
        cells = grow_indices_t(grid, [leaves.cell[id] for id in ids], offsets(grid))
        return Int32[leaves.ids[m][c] for c in cells]
    end
    res = Set{Int32}()
    for id in ids
        l = Int(leaves.level[id])
        c = leaves.cell[id]
        for m in 0:ctx.L
            leaves.n_active[m+1] == 0 && continue
            shape = leaves.shapes[m+1]
            ranges = ntuple(N) do i
                f = ctx.factor[i]
                if m >= l
                    lo, hi = (c[i] - 1) * f^(m - l) + 1, c[i] * f^(m - l)
                else
                    lo = hi = fld(c[i] - 1, f^(l - m)) + 1
                end
                max(lo-1, 1):min(hi+1, shape[i])
            end
            map_m = leaves.ids[m+1]
            for cc in CartesianIndices(ranges)
                nid = map_m[cc]
                (nid != 0 && nid != id) && push!(res, nid)
            end
        end
    end
    return collect(res)
end

function explore_leaves!(ctx::AdaptiveFastContext, s::FlatShieldState, frontier::Vector{Int32})
    n = 0
    while !isempty(frontier)
        for id in frontier
            s.computed[id] && continue
            compute_leaf!(ctx, id)
            register_successors!(ctx.leaves.succ, s, id)
            n += 1
        end
        propagate_bad!(s)
        changed = unique!(s.changed)
        frontier = isempty(changed) ? Int32[] : filter(id -> !s.computed[id], neighbour_leaves(ctx, changed))
        empty!(s.changed)
    end
    return n
end

function init_leaves_online!(ctx::AdaptiveFastContext, level::Int, include_border_cells::Bool, act_needs_safe_flow::Bool)
    leaves = ctx.leaves
    sys = ctx.system
    grid = ctx.levels[level+1].grid
    for c in CartesianIndices(cellshape(grid))
        add_leaf!(ctx, level, c)
    end
    s = reset!(FlatShieldState(length(sys.Act), act_needs_safe_flow), length(leaves.level))
    ids = leaves.ids[level+1]
    frontier = eltype(leaves.cell)[]
    for edge in sys.edges
        append!(frontier, get_touching_cell_idxs_l(grid, edge.guard))
    end
    initial_dead = eltype(leaves.cell)[]
    for gc in sys.globalConstraints, c in get_touching_cell_idxs_l(grid, gc)
        seed_constraint!(s, ids[c])
        push!(initial_dead, c)
    end
    offs = offsets(grid)
    append!(frontier, grow_indices_t(grid, initial_dead, offs))
    if include_border_cells
        border = border_cells(grid)
        append!(frontier, border)
        append!(frontier, grow_indices_t(grid, border, offs))
    end
    unique!(frontier)
    n = explore_leaves!(ctx, s, Int32[ids[c] for c in frontier])
    return s, n
end

# Splits the leaves to_split into their children, recomputes the children and the leaves whose successors were split,
# rebuilds the shield state in place and explores from what changed.
function split_leaves_online!(ctx::AdaptiveFastContext, s::FlatShieldState, to_split::AbstractVector{Int32})
    leaves = ctx.leaves
    succ = leaves.succ
    n_old = length(leaves.level)
    is_split = falses(n_old)
    for id in to_split
        @assert leaves.active[id] && s.computed[id] && leaves.level[id] < ctx.L
        is_split[id] = true
    end
    # Leaves with a successor that is split are recomputed against its children
    affected = Int32[]
    for k in Int32(1):Int32(n_old)
        (!leaves.active[k] || is_split[k] || !s.computed[k] || is_constraint_leaf(ctx, k)) && continue
        any(ai -> any(i -> is_split[succ.pool[i]], successor_range(succ, ai, k)), 0:length(succ.act_enabled)) && push!(affected, k)
    end
    for id in to_split
        l = Int(leaves.level[id])
        c = leaves.cell[id]
        deactivate_leaf!(leaves, id)
        for ch in child_cells(c, ctx.factor)
            add_leaf!(ctx, l + 1, ch)
        end
    end
    n = 0
    for id in Int32(n_old+1):Int32(length(leaves.level))
        is_constraint_leaf(ctx, id) && continue
        compute_leaf!(ctx, id)
        n += 1
    end
    for id in affected
        compute_leaf!(ctx, id)
        n += 1
    end
    # Rebuild from all computed active leaves (the previously computed ones and the new children)
    computed = vcat(s.computed, trues(length(leaves.level) - n_old))
    compact!(succ)
    reset!(s, length(leaves.level))
    for id in Int32(1):Int32(length(leaves.level))
        (leaves.active[id] && computed[id] && is_constraint_leaf(ctx, id)) && seed_constraint!(s, id)
    end
    for id in Int32(1):Int32(length(leaves.level))
        (leaves.active[id] && computed[id] && !is_constraint_leaf(ctx, id)) && register_successors!(succ, s, id)
    end
    propagate_bad!(s)
    changed = copy(unique!(s.changed))
    empty!(s.changed)
    frontier = isempty(changed) ? Int32[] : filter(id -> !s.computed[id], neighbour_leaves(ctx, changed))
    n += explore_leaves!(ctx, s, frontier)
    return s, n
end

function dead_volume(ctx::AdaptiveFastContext{N}, dead::AbstractVector{Bool}) where N
    leaves = ctx.leaves
    v = 0
    for id in eachindex(leaves.level)
        if leaves.active[id] && dead[id]
            v += prod(ntuple(i -> ctx.factor[i]^(ctx.L - leaves.level[id]), N))
        end
    end
    return v
end

# Predecessors for choice ai (0 = flow) of the active leaves k with usable(k), as one flat array: the predecessors of v
# are list[start[v]:start[v+1]-1], in increasing order.
function flat_predecessors(leaves::FlatLeaves, ai::Int, usable)
    n = length(leaves.level)
    pool = leaves.succ.pool
    start = zeros(Int, n + 1)
    for k in 1:n
        (leaves.active[k] && usable(k)) || continue
        for i in successor_range(leaves.succ, ai, k)
            start[pool[i]+1] += 1
        end
    end
    start[1] = 1
    for v in 1:n
        start[v+1] += start[v]
    end
    list = Vector{Int32}(undef, start[n+1] - 1)
    next = copy(start)
    for k in Int32(1):Int32(n)
        (leaves.active[k] && usable(k)) || continue
        for i in successor_range(leaves.succ, ai, k)
            v = pool[i]
            list[next[v]] = k
            next[v] += 1
        end
    end
    return start, list
end

# The shield over all active leaves from scratch, returning the dead leaves and the order they died in, which the
# refinement uses to pick the latest-dying leaves.
function adaptive_shield(leaves::FlatLeaves; act_needs_safe_flow::Bool=true, computed::AbstractVector{Bool})
    n = length(leaves.level)
    succ = leaves.succ
    n_acts = length(succ.act_enabled)
    flow_start, flow_list = flat_predecessors(leaves, 0, k -> true)
    act_preds = [flat_predecessors(leaves, ai, k -> succ.act_enabled[ai][k]) for ai in 1:n_acts]
    act_ok = [falses(n) for _ in 1:n_acts]
    n_ok = zeros(Int32, n)
    for k in 1:n, ai in 1:n_acts
        if leaves.active[k] && succ.act_enabled[ai][k] && succ.act_len[ai][k] > 0
            act_ok[ai][k] = true
            n_ok[k] += 1
        end
    end
    dead = falses(n)
    flow_bad = falses(n)
    order = zeros(Int32, n)
    queue_dead = Int32[]
    queue_fb = Int32[]
    function try_kill!(k)
        if !dead[k] && flow_bad[k] && n_ok[k] == 0 && computed[k]
            dead[k] = true
            push!(queue_dead, k)
            order[k] = length(queue_dead)
        end
    end
    function set_flow_bad!(k)
        if !flow_bad[k]
            flow_bad[k] = true
            push!(queue_fb, k)
            try_kill!(k)
        end
    end
    function invalidate_into!(v)
        for ai in 1:n_acts
            start, list = act_preds[ai]
            for i in start[v]:start[v+1]-1
                k = list[i]
                if act_ok[ai][k]
                    act_ok[ai][k] = false
                    n_ok[k] -= 1
                    try_kill!(k)
                end
            end
        end
    end
    for k in Int32(1):Int32(n)
        if leaves.active[k] && computed[k] && (succ.flow_none[k] || succ.flow_len[k] == 0)
            set_flow_bad!(k)
        end
    end
    hd, hf = 1, 1
    while hd <= length(queue_dead) || hf <= length(queue_fb)
        while hd <= length(queue_dead)
            v = queue_dead[hd]
            hd += 1
            for i in flow_start[v]:flow_start[v+1]-1
                set_flow_bad!(flow_list[i])
            end
            invalidate_into!(v)
        end
        while hf <= length(queue_fb)
            v = queue_fb[hf]
            hf += 1
            (act_needs_safe_flow && !dead[v]) && invalidate_into!(v)
        end
    end
    return dead, order
end

# What each finest cell allows and lost (see choice_bits).
adaptive_allowed(ctx::AdaptiveFastContext, s::FlatShieldState) = leaf_choices(cellshape(ctx.levels[end].grid), length(ctx.system.Act),
    ctx.leaves, ctx.L, ctx.factor, s, store_enabled(ctx.leaves.succ), store_acted(ctx.leaves.succ))

# Computes the online shield at max_granularity (level 0, as ReACTed_reachable_cell_online_pred, including
# include_border_cells), then refines in rounds up to granularity_steps halvings: every round splits the refine_fraction
# latest-dying dead leaves, recomputes what changed and continues the exploration. Stops when nothing is left to split,
# after `patience` rounds in a row that revive nothing, or after max_steps rounds. Returns grid, allowed, lost, info on
# the finest grid (see choice_bits), with the context as info.ctx. on_round(ctx, dead, round, leaf_computations) is
# called after every round, e.g. for experiments.
function ReACTed_reachable_cell_adaptive_fast(system::EuclideanHybridSystem, p, max_granularity, granularity_steps::Int64, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, max_steps::Int64=1000; include_border_cells::Bool=false, refine_fraction::Float64=0.4, patience::Int=1, act_needs_safe_flow::Bool=true, sweep_factor::Int=20, on_round=(ctx, dead, r, computed) -> nothing, jump_hull::Bool=false)
    recorderTimeStart = time_ns()
    ctx = AdaptiveFastContext(system, p, Float64(max_granularity), granularity_steps, δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor; jump_hull)
    setup_time = time_ns()
    reachability_time = 0.0
    shield_time = 0.0

    aux_time_recoder = time_ns()
    s, computed = init_leaves_online!(ctx, 0, include_border_cells, act_needs_safe_flow)
    reachability_time += time_ns() - aux_time_recoder
    aux_time_recoder = time_ns()
    dead, order = adaptive_shield(ctx.leaves; act_needs_safe_flow, computed=s.computed)
    shield_time += time_ns() - aux_time_recoder
    volume = dead_volume(ctx, dead)
    println("Level 0: $(computed) leaves computed, dead volume $volume")
    on_round(ctx, dead, 0, computed)

    rounds = 0
    rounds_without_gain = 0
    for r in 1:max_steps
        leaves = ctx.leaves
        candidates = [Int32(id) for id in eachindex(leaves.level) if leaves.active[id] && dead[id] && leaves.level[id] < ctx.L]
        isempty(candidates) && break
        candidates = candidates[sortperm(order[candidates]; rev=true)] # a closure over order would capture it boxed, as order is reassigned below
        chosen = candidates[1:max(1, round(Int, refine_fraction*length(candidates)))]

        aux_time_recoder = time_ns()
        s, n = split_leaves_online!(ctx, s, chosen)
        computed += n
        reachability_time += time_ns() - aux_time_recoder
        aux_time_recoder = time_ns()
        dead, order = adaptive_shield(ctx.leaves; act_needs_safe_flow, computed=s.computed)
        shield_time += time_ns() - aux_time_recoder

        rounds = r
        new_volume = dead_volume(ctx, dead)
        println("Round $r: split $(length(chosen)) leaves, active per level $(ctx.leaves.n_active), dead volume $new_volume (revived $(volume - new_volume)), $(computed) leaf computations : $(Dates.format(now(), "HH:MM:SS"))")
        on_round(ctx, dead, r, computed)
        rounds_without_gain = new_volume < volume ? 0 : rounds_without_gain + 1
        volume = new_volume
        rounds_without_gain >= patience && break
    end

    allowed, lost = adaptive_allowed(ctx, s)
    fine = ctx.levels[end].grid
    grid = Grid(length(fine.numCells), fine.granularity, copy(fine.lower), copy(fine.upper), copy(fine.numCells), allowed .== 0)
    println("Total times:
        Total: $((time_ns() - recorderTimeStart)/10^9)
        Setup: $((setup_time - recorderTimeStart) / 10^9)
        Reachability: $(reachability_time / 10^9)
        Shield: $(shield_time / 10^9)")
    println("Refinement rounds: $rounds, leaf computations: $computed")
    return grid, allowed, lost, (; ctx)
end
