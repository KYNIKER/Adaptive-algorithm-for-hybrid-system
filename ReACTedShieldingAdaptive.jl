using LazySets, LinearAlgebra, ReachabilityAnalysis, Dates
# Expects PartitionUtilities.jl and ReACTedShielding.jl to be included first.

#=
Adaptive (multi-resolution) online shield synthesis.
The state space is partitioned into leaves of a cell tree: level 0 is the grid at max_granularity, and a leaf at level l can be split into the cells covering it at level l+1 (granularity halved; zero-width axes are never split). 
The shield is computed over the active leaves with the same rules as the event-driven predecessor shield (ReACTedShieldingPred.jl):
  flow is bad   : no flow successors (hit a constraint), none at all, or a flow successor is dead
  action is ok  : the leaf touches the action's guard and every successor is alive (and, if act_needs_safe_flow, can wait)
  leaf is dead  : flow is bad and no action is ok
As in ReACTed_reachable_cell_online_pred, only leaves reached by the frontier are computed and the others are assumed safe; constraint leaves are dead and never computed. 
Every leaf's successors over-approximate its own reach set, so any mix of levels gives a valid shield; splitting only changes permissiveness. 
A leaf's successors are the active leaves its reach set touches, found with one cover test per level, so with all leaves on one level the algorithm computes exactly the same cells and shield as ReACTed_reachable_cell_online_pred at that granularity.
=#
# Fast-cell matrix for one level that classifies sweep blocks on first access (and caches them), so fine levels are only swept where cells are actually used. 
#Reads give the same values as the full sweep in ReACTed_reachable_cell_t.
mutable struct LazyFastCells{N,F} <: AbstractArray{Bool,N}
    dims::NTuple{N,Int}
    sweep::NTuple{N,Int}
    nblocks::NTuple{N,Int}
    known::BitArray{N}
    fast::BitArray{N}
    check::F # block index -> Bool
end

Base.size(A::LazyFastCells) = A.dims

# Check whether we know that the lazy cell is fast or not. If it is not in a full block for the sweep, then it will always be false, as that the hybrid jump check fails for those intitially.
function Base.getindex(A::LazyFastCells{N}, I::Vararg{Int,N}) where N
    b = CartesianIndex(ntuple(i -> fld(I[i] - 1, A.sweep[i]) + 1, N))
    for i in 1:N
        if b[i] > A.nblocks[i]
            return false # cells beyond the last full block are never fast
        end
    end
    if !A.known[b]
        A.fast[b] = A.check(b)
        A.known[b] = true
    end
    return A.fast[b]
end

# Leaves of the cell tree. Ids are never reused; a split leaf becomes inactive and its children get new ids.
mutable struct AdaptiveLeaves{N}
    level::Vector{Int8}                       # 0-based level
    cell::Vector{CartesianIndex{N}}           # cell index in its level's grid
    active::BitVector
    flow::Vector{Vector{Int32}}               # flow successors (leaf ids)
    flow_none::BitVector                      # no flow successors at all (flow hit a constraint, or constraint cell)
    act::Vector{Vector{Vector{Int32}}}        # per leaf, per action position: successors
    act_enabled::Vector{BitVector}            # per leaf: which actions' guards the leaf touches
    ids::Vector{Dict{CartesianIndex{N},Int32}} # per level: cell -> active leaf id
    n_active::Vector{Int}                     # per level
end

AdaptiveLeaves{N}(L::Int) where N = AdaptiveLeaves{N}(Int8[], CartesianIndex{N}[], BitVector(), Vector{Int32}[], BitVector(), Vector{Vector{Int32}}[], BitVector[], [Dict{CartesianIndex{N},Int32}() for _ in 0:L], zeros(Int, L + 1))

function add_leaf!(leaves::AdaptiveLeaves{N}, level::Int, cell::CartesianIndex{N}, n_acts::Int) where N
    push!(leaves.level, Int8(level))
    push!(leaves.cell, cell)
    push!(leaves.active, true)
    push!(leaves.flow, Int32[])
    push!(leaves.flow_none, false)
    push!(leaves.act, [Int32[] for _ in 1:n_acts])
    push!(leaves.act_enabled, falses(n_acts))
    id = Int32(length(leaves.level))
    leaves.ids[level+1][cell] = id
    leaves.n_active[level+1] += 1
    return id
end

function deactivate_leaf!(leaves::AdaptiveLeaves, id::Int32)
    l = leaves.level[id]
    leaves.active[id] = false
    delete!(leaves.ids[l+1], leaves.cell[id])
    leaves.n_active[l+1] -= 1
    leaves.flow[id] = Int32[]
    foreach(empty!, leaves.act[id])
end

# Grid for level l, obtained by subdividing the level-0 grid (not by Grid(statespace, g/2^l), whose ceil can break the
# nesting). Zero-width axes (factor 1) keep a single cell centred on their value.
function level_grid(grid0::Grid{N}, l::Int, factor::NTuple{N,Int}) where N
    g = grid0.granularity / 2^l
    numCells = [grid0.numCells[i] * factor[i]^l for i in 1:N]
    lower = copy(grid0.lower)
    upper = copy(grid0.upper)
    for i in 1:N
        if factor[i] == 1
            value = grid0.lower[i] + grid0.granularity / 2
            lower[i] = value - g / 2
            upper[i] = value + g / 2
        end
    end
    return Grid(N, g, lower, upper, numCells, falses(numCells...))
end

# Cells of level l+1 covering cell c of level l.
child_cells(c::CartesianIndex{N}, factor::NTuple{N,Int}) where N =
    CartesianIndices(ntuple(i -> ((c[i]-1)*factor[i]+1):(c[i]*factor[i]), N))

# Per-level reachability data (as set up in ReACTed_reachable_cell_t for a uniform grid of this granularity).
function adaptive_level(system::EuclideanHybridSystem, sh, grid0::Grid{N}, l::Int, L::Int, factor::NTuple{N,Int}, sweep_factor::Int) where N
    grid = level_grid(grid0, l, factor)
    g = grid.granularity
    (; δ⁻, δ⁺, A, P2A_abs, phiDict, tPhiDict, U, inputDict, max_input, max_flow) = sh

    # Lazy fast-cell blocks
    sweep = [min(sweep_factor, n) for n in grid.numCells]
    Z_big = remove_zero_generators(Zonotope(zeros(N), diagm(map((n, s) -> n == 1 ? 0. : g * (s/2), grid.numCells, sweep))))
    lower_offset_big = grid.lower .+ g .* (sweep ./ 2)
    generatorDict_big = ReACT_discretize_decomposed_generators(Z_big, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict)
    discretizationDict_big, finaldims_big = ReACT_discretize_combine_with_offsets(Z_big, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict_big)
    pc_big, pg_big, pig_big = zeros(N), zeros(N, finaldims_big), zeros(N, finaldims_big)
    check = function (bidx::CartesianIndex{N})
        of = lower_offset_big + ((car2vec(bidx) .- 1) .* g .* sweep)
        LazySets.API.translate!(Z_big, of)
        ReACT_discretize_combine_with_offset_vector!(pc_big, pg_big, pig_big, discretizationDict_big, of, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict_big)
        f = propagate_set_check_hybrid_jump(Z_big, [0.0, sh.p], δ⁻, δ⁺, system, system.edges, sh.identity_matrix, sh.attempt_recorder, N, P2A_abs, sh.guardProjVectors, sh.guardProjBounds, sh.invariantProjVectors, sh.invariantProjBounds, discretizationDict_big, phiDict, tPhiDict, inputDict, max_input, max_flow, sh.guard_invariant_intersection, sh.alg, sh.maxOrder, sh.reduceOrder)
        LazySets.API.translate!(Z_big, -of)
        return f == 0
    end
    nblocks = ntuple(i -> fld(grid.numCells[i], sweep[i]), N)
    fast = LazyFastCells(cellshape(grid), Tuple(sweep), nblocks, falses(nblocks), falses(nblocks), check)

    # Cell reach sets
    Z = remove_zero_generators(Zonotope(zeros(N), diagm(map(x -> x == 1 ? 0. : g/2, grid.numCells))))
    generatorDict = ReACT_discretize_decomposed_generators(Z, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict)
    discretizationDict, finaldims = ReACT_discretize_combine_with_offsets(Z, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict)

    # Fast cells of this level all reach a translate of the same zonotope; one cover per target level m (bloated by half
    # a level-m cell). For m == l this is exactly the cover used in the uniform algorithms.
    reach_gens = hcat(genmat(max_input), genmat(linear_map(max_flow, Z)))
    fast_covers = [ConstantZonotopeCover(Zonotope(max_input.center, hcat(reach_gens, diagm(fill(grid0.granularity / 2^m / 2, N)))), Val(N)) for m in 0:L]

    may_act = [falses(cellshape(grid)) for _ in system.Act]
    for (ai, act) in enumerate(system.Act), c in get_touching_cell_idxs_l(grid, act.guard)
        may_act[ai][c] = true
    end
    constraint = falses(cellshape(grid))
    for gc in system.globalConstraints, c in get_touching_cell_idxs_l(grid, gc)
        constraint[c] = true
    end
    return (; grid, g, Z, generatorDict, discretizationDict, pc=zeros(N), pg=zeros(N, finaldims), pig=zeros(N, finaldims),
        lower_offset=grid.lower .+ g / 2, fast, fast_covers, may_act, constraint)
end

struct AdaptiveContext{N,S,SH,LV}
    system::S
    shared::SH
    levels::Vector{LV}
    factor::NTuple{N,Int}
    leaves::AdaptiveLeaves{N}
    L::Int
end

function AdaptiveContext(system::EuclideanHybridSystem, p, max_granularity::Float64, L::Int, δ⁻::Float64, δ⁺::Float64, alg, maxOrder::Int, reduceOrder::Int, sweep_factor::Int; jump_hull::Bool=false)
    grid0 = Grid(system.statespace, max_granularity)
    return _AdaptiveContext(grid0, system, p, L, δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor; jump_hull)
end

# Function barrier on Grid{N}
function _AdaptiveContext(grid0::Grid{N}, system, p, L, δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor; jump_hull::Bool=false) where N
    zero_width = low(system.statespace) .== high(system.statespace)
    factor = ntuple(i -> zero_width[i] ? 1 : 2, N)
    sh = reach_setup(system, p, δ⁻, δ⁺, alg, maxOrder, reduceOrder, Val(N); jump_hull)
    levels = [adaptive_level(system, sh, grid0, l, L, factor, sweep_factor) for l in 0:L]
    return AdaptiveContext(system, sh, levels, factor, AdaptiveLeaves{N}(L), L)
end

# Active leaves touched by a set, given a function cells_at(m) returning the touched cells of level m's grid.
function touched_leaves!(out::Vector{Int32}, ctx::AdaptiveContext, cells_at)
    empty!(out)
    leaves = ctx.leaves
    for m in 0:ctx.L
        if leaves.n_active[m+1] == 0
            continue
        end
        ids = leaves.ids[m+1]
        for c in cells_at(m)
            id = get(ids, c, Int32(0))
            if id != 0
                push!(out, id)
            end
        end
    end
    return out
end

# Computes the flow and action successors of leaf id against the current set of active leaves. Not called for
# constraint leaves, which are dead and never computed.
function compute_leaf!(ctx::AdaptiveContext{N}, id::Int32) where N
    leaves, sh, sys = ctx.leaves, ctx.shared, ctx.system
    l = Int(leaves.level[id])
    c = leaves.cell[id]
    lev = ctx.levels[l+1]
    (; δ⁻, δ⁺, A, P2A_abs, phiDict, tPhiDict, U, inputDict, max_input, max_flow) = sh

    flow = leaves.flow[id]
    if lev.fast[c]
        fast_cell_center!(sh.center, max_flow, max_input.center, lev.grid, c)
        touched_leaves!(flow, ctx, m -> get_touching_cell_idxs_cover(ctx.levels[m+1].grid, lev.fast_covers[m+1], sh.center))
        leaves.flow_none[id] = false
    else
        of = lev.lower_offset + ((car2vec(c) .- 1) .* lev.g)
        LazySets.API.translate!(lev.Z, of)
        lev.pc .= 0
        lev.pg .= 0
        lev.pig .= 0
        ReACT_discretize_combine_with_offset_vector!(lev.pc, lev.pg, lev.pig, lev.discretizationDict, of, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, lev.generatorDict)
        z, f = propagate_set_jumps(sh.jump_reach, lev.Z, lev.discretizationDict)
        LazySets.API.translate!(lev.Z, -of)
        if f == 1
            leaves.flow_none[id] = true
            empty!(flow)
        else
            leaves.flow_none[id] = false
            touched_leaves!(flow, ctx, m -> get_touching_cell_idxs_lz(ctx.levels[m+1].grid, z))
        end
    end

    for (ai, act) in enumerate(sys.Act)
        if lev.may_act[ai][c]
            leaves.act_enabled[id][ai] = true
            of = lev.lower_offset + ((car2vec(c) .- 1) .* lev.g)
            LazySets.API.translate!(lev.Z, of)
            az = minkowski_sum(linear_map(max_flow, LazySets.API.translate(linear_map(act.jumpMatrix, lev.Z), act.jumpVector)), max_input)
            LazySets.API.translate!(lev.Z, -of)
            touched_leaves!(leaves.act[id][ai], ctx, m -> get_touching_cell_idxs_lz(ctx.levels[m+1].grid, az))
        else
            leaves.act_enabled[id][ai] = false
            empty!(leaves.act[id][ai])
        end
    end
    return id
end

# Event-driven shield over the active leaves (same rules as OnlinePredShieldState). Returns the dead flags and the order
# in which leaves died (0 = alive), which the refinement uses to pick the latest-dying leaves. Leaves outside the
# computed mask were never computed and are assumed safe (never seeded as bad and never dead).
function adaptive_shield(leaves::AdaptiveLeaves; act_needs_safe_flow::Bool=true, computed::Union{Nothing,AbstractVector{Bool}}=nothing)
    n = length(leaves.level)
    n_acts = isempty(leaves.act_enabled) ? 0 : length(leaves.act_enabled[1])
    pred_flow = [Int32[] for _ in 1:n]
    pred_act = [[Int32[] for _ in 1:n] for _ in 1:n_acts]
    act_ok = [falses(n) for _ in 1:n_acts]
    n_ok = zeros(Int32, n)
    for k in Int32(1):Int32(n)
        leaves.active[k] || continue
        for v in leaves.flow[k]
            push!(pred_flow[v], k)
        end
        for ai in 1:n_acts
            if leaves.act_enabled[k][ai] && !isempty(leaves.act[k][ai])
                act_ok[ai][k] = true
                n_ok[k] += 1
                for v in leaves.act[k][ai]
                    push!(pred_act[ai][v], k)
                end
            end
        end
    end
    dead = falses(n)
    flow_bad = falses(n)
    order = zeros(Int32, n)
    ndead = Ref(Int32(0))
    queue_dead = Int32[]
    queue_fb = Int32[]
    try_kill! = k -> begin
        if !dead[k] && flow_bad[k] && n_ok[k] == 0 && (computed === nothing || computed[k])
            dead[k] = true
            ndead[] += 1
            order[k] = ndead[]
            push!(queue_dead, k)
        end
    end
    set_flow_bad! = k -> begin
        if !flow_bad[k]
            flow_bad[k] = true
            push!(queue_fb, k)
            try_kill!(k)
        end
    end
    invalidate! = (ai, k) -> begin
        if act_ok[ai][k]
            act_ok[ai][k] = false
            n_ok[k] -= 1
            try_kill!(k)
        end
    end
    for k in Int32(1):Int32(n)
        if leaves.active[k] && (computed === nothing || computed[k]) && (leaves.flow_none[k] || isempty(leaves.flow[k]))
            set_flow_bad!(k)
        end
    end
    hd, hf = 1, 1
    while hd <= length(queue_dead) || hf <= length(queue_fb)
        while hd <= length(queue_dead)
            v = queue_dead[hd]
            hd += 1
            for k in pred_flow[v]
                set_flow_bad!(k)
            end
            for ai in 1:n_acts, k in pred_act[ai][v]
                invalidate!(ai, k)
            end
        end
        while hf <= length(queue_fb)
            v = queue_fb[hf]
            hf += 1
            if !act_needs_safe_flow || dead[v]
                continue
            end
            for ai in 1:n_acts, k in pred_act[ai][v]
                invalidate!(ai, k)
            end
        end
    end
    return dead, order
end

# ---------------------------------------------------------------------------------------------------------------------
# Frontier exploration. It reuses the predecessor shield (OnlinePredShieldState, register_cell!, propagate_bad! from
# ReACTedShieldingPred.jl) with leaf ids as 1-D cell indices.

const LeafIdx = CartesianIndex{1}

new_leaf_state(n::Int, n_acts::Int, act_needs_safe_flow::Bool) = OnlinePredShieldState(falses(n), n_acts, act_needs_safe_flow)

is_constraint_leaf(ctx::AdaptiveContext, id::Integer) = ctx.levels[ctx.leaves.level[id]+1].constraint[ctx.leaves.cell[id]]

# Marks a constraint leaf as dead and computed (as the predecessor algorithm does for constraint cells).
function seed_constraint_leaf!(s::OnlinePredShieldState{1}, id::Int32)
    k = LeafIdx(id)
    s.dead[k] = true
    s.flow_bad[k] = true
    s.computed[k] = true
end

# Hands a computed leaf's successors to the predecessor shield.
function register_leaf!(ctx::AdaptiveContext, s::OnlinePredShieldState{1}, id::Int32, act_buf::Vector{Union{Nothing,Vector{LeafIdx}}})
    leaves = ctx.leaves
    flow_succ = leaves.flow_none[id] ? nothing : LeafIdx.(leaves.flow[id])
    for ai in eachindex(act_buf)
        act_buf[ai] = leaves.act_enabled[id][ai] ? LeafIdx.(leaves.act[id][ai]) : nothing
    end
    register_cell!(s, LeafIdx(id), flow_succ, act_buf)
end

# Active leaves spatially adjacent to the given leaves (including diagonals). With all leaves on one level this is
# grow_indices_t on that level's grid, in the same order, so the exploration matches the predecessor algorithm.
function neighbour_leaves(ctx::AdaptiveContext{N}, ids::Vector{Int32}) where N
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
            shape = cellshape(ctx.levels[m+1].grid)
            ranges = ntuple(N) do i
                f = ctx.factor[i]
                if m >= l # Leaf or children of leaf 
                    lo, hi = (c[i] - 1) * f^(m - l) + 1, c[i] * f^(m - l)
                else # Anchestor of leaf
                    lo = hi = fld(c[i] - 1, f^(l - m)) + 1
                end
                # neighbour leaves
                max(lo-1, 1):min(hi+1, shape[i])
            end
            for cc in CartesianIndices(ranges)
                nid = get(leaves.ids[m+1], cc, Int32(0))
                (nid != 0 && nid != id) && push!(res, nid)
            end
        end
    end
    # computes a superset of the leaves we want, but this is okay because they arent active anyway.
    return collect(res)
end

# Frontier exploration rounds (the loop of ReACTed_reachable_cell_online_pred) starting from the given frontier.
# Returns the number of leaves computed.
function explore_leaves!(ctx::AdaptiveContext, s::OnlinePredShieldState{1}, frontier::Vector{Int32})
    act_buf = Vector{Union{Nothing,Vector{LeafIdx}}}(nothing, length(ctx.system.Act))
    n = 0
    while !isempty(frontier)
        for id in frontier
            s.computed[id] && continue
            compute_leaf!(ctx, id)
            register_leaf!(ctx, s, id, act_buf)
            n += 1
        end
        propagate_bad!(s)
        changed = Int32[k[1] for k in unique!(s.changed)]
        frontier = isempty(changed) ? Int32[] : filter(id -> !s.computed[id], neighbour_leaves(ctx, changed))
        empty!(s.changed)
    end
    return n
end

# Creates all leaves of one level (uncomputed), seeds constraint leaves and explores from the same initial frontier as
# the predecessor algorithm (edge guards, neighbours of constraint cells and, optionally, the grid border).
function init_leaves_online!(ctx::AdaptiveContext, level::Int, include_border_cells::Bool, act_needs_safe_flow::Bool)
    leaves = ctx.leaves
    sys = ctx.system
    grid = ctx.levels[level+1].grid
    for c in CartesianIndices(cellshape(grid))
        add_leaf!(leaves, level, c, length(sys.Act))
    end
    s = new_leaf_state(length(leaves.level), length(sys.Act), act_needs_safe_flow)
    ids = leaves.ids[level+1]
    frontier = eltype(leaves.cell)[]
    for edge in sys.edges
        append!(frontier, get_touching_cell_idxs_l(grid, edge.guard))
    end
    initial_dead = eltype(leaves.cell)[]
    for gc in sys.globalConstraints, c in get_touching_cell_idxs_l(grid, gc)
        seed_constraint_leaf!(s, ids[c])
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

# Splits leaves in online mode: children are computed (or seeded if they touch a constraint), affected predecessors are
# recomputed, the shield state is rebuilt (refinement can revive leaves, which the incremental state cannot undo) and
# exploration continues from whatever changed. Returns the new state and the number of leaves computed.
function split_leaves_online!(ctx::AdaptiveContext, s::OnlinePredShieldState{1}, to_split::AbstractVector{Int32}, act_needs_safe_flow::Bool)
    leaves = ctx.leaves
    n_acts = length(ctx.system.Act)
    n_old = length(leaves.level)
    is_split = falses(n_old)
    for id in to_split
        @assert leaves.active[id] && s.computed[id] && leaves.level[id] < ctx.L
        is_split[id] = true
    end
    affected = Int32[]
    for k in Int32(1):Int32(n_old)
        (!leaves.active[k] || is_split[k] || !s.computed[k] || is_constraint_leaf(ctx, k)) && continue
        hit = any(v -> is_split[v], leaves.flow[k]) || any(lst -> any(v -> is_split[v], lst), leaves.act[k])
        hit && push!(affected, k)
    end
    for id in to_split
        l = Int(leaves.level[id])
        c = leaves.cell[id]
        deactivate_leaf!(leaves, id)
        for ch in child_cells(c, ctx.factor)
            add_leaf!(leaves, l + 1, ch, n_acts)
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
    # Rebuild from all computed active leaves (the previously computed ones and the new children).
    computed = vcat(s.computed, trues(length(leaves.level) - n_old))
    s2 = new_leaf_state(length(leaves.level), n_acts, act_needs_safe_flow)
    act_buf = Vector{Union{Nothing,Vector{LeafIdx}}}(nothing, n_acts)
    for id in Int32(1):Int32(length(leaves.level))
        (leaves.active[id] && computed[id] && is_constraint_leaf(ctx, id)) && seed_constraint_leaf!(s2, id)
    end
    for id in Int32(1):Int32(length(leaves.level))
        (leaves.active[id] && computed[id] && !is_constraint_leaf(ctx, id)) && register_leaf!(ctx, s2, id, act_buf)
    end
    propagate_bad!(s2)
    changed = Int32[k[1] for k in unique!(s2.changed)]
    empty!(s2.changed)
    frontier = isempty(changed) ? Int32[] : filter(id -> !s2.computed[id], neighbour_leaves(ctx, changed))
    n += explore_leaves!(ctx, s2, frontier)
    return s2, n
end

# Dead volume in finest-level cells (a leaf at level l covers prod(factor .^ (L - l)) of them).
function dead_volume(ctx::AdaptiveContext{N}, dead::AbstractVector{Bool}) where N
    leaves = ctx.leaves
    v = 0
    for id in eachindex(leaves.level)
        if leaves.active[id] && dead[id]
            v += prod(ntuple(i -> ctx.factor[i]^(ctx.L - leaves.level[id]), N))
        end
    end
    return v
end

# What each finest cell allows and lost (see choice_bits in ReACTedShieldingFlat.jl). A leaf is enabled for an action
# when it touches the guard on its level.
function adaptive_allowed(ctx::AdaptiveContext{N}, s::OnlinePredShieldState{1}) where N
    leaves = ctx.leaves
    enabled = (ai, id) -> ctx.levels[leaves.level[id]+1].may_act[ai][leaves.cell[id]]
    acted = (ai, id) -> !isempty(leaves.act[id][ai])
    return leaf_choices(cellshape(ctx.levels[end].grid), length(ctx.system.Act), leaves, ctx.L, ctx.factor, s, enabled, acted)
end

# Adaptive online shield synthesis. Computes the online shield at max_granularity (level 0, exactly as
# ReACTed_reachable_cell_online_pred, including include_border_cells), then refines in rounds up to
# granularity_steps halvings: every round splits the refine_fraction latest-dying dead leaves (the tail of the death
# cascade, which calibration showed is where finer cells revive states), recomputes what changed and continues the
# exploration. Refinement stops when nothing is left to split, after `patience` consecutive rounds that revive nothing
# (no decrease in dead volume), or after max_steps rounds.
# Returns grid, allowed, lost, info like every shield (see choice_bits), on the finest grid, with the context as
# info.ctx. on_round(ctx, dead, round, leaf_computations) is called after every round, e.g. for experiments.
function ReACTed_reachable_cell_adaptive(system::EuclideanHybridSystem, p, max_granularity, granularity_steps::Int64, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, clustering=true, max_steps::Int64=1000; include_border_cells::Bool=false, refine_fraction::Float64=0.4, patience::Int=1, act_needs_safe_flow::Bool=true, sweep_factor::Int=20, on_round=(ctx, dead, r, computed) -> nothing, jump_hull::Bool=false)
    recorderTimeStart = time_ns()
    TIME_FUNCS = true
    ctx = AdaptiveContext(system, p, Float64(max_granularity), granularity_steps, δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor; jump_hull)
    fine = ctx.levels[end].grid
    @show foldl(*, fine.numCells), fine.numCells, fine.granularity, size(ctx.levels)
    setup_time = time_ns()
    reachability_time = 0.0
    shield_time = 0.0

    aux_time_recoder = time_ns()
    s, computed = init_leaves_online!(ctx, 0, include_border_cells, act_needs_safe_flow)
    reachability_time += time_ns() - aux_time_recoder
    @show reachability_time/10^9
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
        sort!(candidates; by=id -> order[id], rev=true)
        chosen = candidates[1:max(1, round(Int, refine_fraction*length(candidates)))]

        aux_time_recoder = time_ns()
        s, n = split_leaves_online!(ctx, s, chosen, act_needs_safe_flow) # includes rebuilding the shield state
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

    if TIME_FUNCS
        println("Total times:
        Total: $((time_ns() - recorderTimeStart)/10^9)
        Setup: $((setup_time - recorderTimeStart) / 10^9)
        Reachability: $(reachability_time / 10^9)
        Shield: $(shield_time / 10^9)")
        println("Refinement rounds: $rounds, leaf computations: $computed")
    end

    #@show ctx.leaves #, size(ctx.levels[1])

    allowed, lost = adaptive_allowed(ctx, s)
    fine = ctx.levels[end].grid
    grid = Grid(length(fine.numCells), fine.granularity, copy(fine.lower), copy(fine.upper), copy(fine.numCells), allowed .== 0)
    return grid, allowed, lost, (; ctx)
end
