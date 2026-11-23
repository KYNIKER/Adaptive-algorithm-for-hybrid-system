# Adaptive (multi-resolution) shield synthesis for a ModalHybridSystem.
#
# The state space is partitioned into leaves: aligned boxes of cells of the finest grid (granularity max_granularity /
# 2^L), in one layer. A leaf at level l is 2^(L-l) finest cells wide on every axis (zero-width axes are never split);
# level 0 is the coarsest, and a leaf at level L is one finest cell. Every leaf is flowed forward as one set
# (box_successors), and its successors are kept per choice
# as finest cells. So splitting a leaf never changes what the other leaves reach: each round maps the finest cells to the
# leaves that now cover them, and only the new leaves are computed.
#
# A leaf's successors are intersected with its parent's, per choice, and a choice that is bad for the leaf but not for its
# parent takes the parent's successors: the parent's are a sound over-approximation of everything the leaf reaches (the
# same rule as for the blocks of the modal shield). So splitting a leaf can only allow more, never less.
#
# Every block of level-0 leaves (about sweep_factor finest cells per axis) is also flowed forward once, and a level-0 leaf's
# successors are intersected with its block's in the same way, so every leaf's successors lie in its block's. Online, only
# the leaves of blocks that are bad, or that reach a leaf whose choices changed, are computed; every other leaf allows
# everything, which is what it allows when computed (see the online part of modal.jl). So online gives the offline result.
#
# Refinement: every round splits the leaves that lost a choice (or are dead, or touch an unsafe set) and border a leaf
# with a different set of allowed choices. Those lie on the boundary of the safe region, where a finer partition can
# separate safe states from unsafe ones; a leaf deep inside a dead or fully allowed region gains nothing from it.

# Splits of a level-l cell: f[i] finest cells per axis
level_factor(L::Int, l::Int, numCells) = ntuple(i -> numCells[i] == 1 ? 1 : 2^(L - l), length(numCells))

# Returns grid, allowed, lost, info on the finest grid like every shield (see choice_bits in flat.jl), with info.ctx the
# modal context on the finest grid, info.levels the number of active leaves per level and info.leaf_level the level of
# the leaf covering each finest cell. It runs at most max_rounds
# refinement rounds. The options are those of ReACTed_reachable_cell_modal.
function ReACTed_reachable_cell_modal_adaptive(system::ModalHybridSystem, p, max_granularity, L::Int, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5;
                                               online::Bool=false, act_needs_safe_flow::Bool=true, exit_bad_for_flow::Bool=true, exit_bad_for_actions::Bool=true, max_rounds::Int=typemax(Int), sweep_factor=5, max_depth::Int=16, verbose::Bool=true)
    timer = PhaseTimer()
    # The modal context on the finest grid. A leaf whose cells are all fast uses a cover, whatever the sweep blocks are.
    ctx = ModalContext(system, p, Float64(max_granularity) / 2^L, δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor, max_depth; exit_bad_for_flow, exit_bad_for_actions)
    lap!(timer, "Setup")
    return _ReACTed_reachable_cell_modal_adaptive(ctx, L, online, act_needs_safe_flow, max_rounds, verbose, timer)
end

function _ReACTed_reachable_cell_modal_adaptive(ctx::ModalContext{D,N}, L::Int, online::Bool, act_needs_safe_flow::Bool, max_rounds::Int, verbose::Bool, timer::PhaseTimer) where {D,N}
    grid = ctx.grid
    nc = Tuple(grid.numCells)
    n_acts = length(ctx.action_ids)
    T = choice_type(n_acts)
    fshape = (cellshape(grid)..., ctx.n_layers)
    fids = LinearIndices(fshape)
    unsafe = falses(fshape)
    for H in ctx.system.unsafe, c in touching_cells(grid, H), layer in 1:ctx.n_layers
        unsafe[c, layer] = true
    end

    # Leaves, as parallel vectors indexed by leaf id, and the active leaf covering each finest cell
    level, cell, layer_of, parent = Int[], CartesianIndex{D}[], Int[], Int32[]
    active, constraint, starts, computed = BitVector(), BitVector(), BitVector(), BitVector()
    enabled = [BitVector() for _ in 1:n_acts] # some state of the leaf can take the action
    flow_succ = Union{Nothing,Vector{Int32}}[] # finest-cell ids, nothing if waiting is bad
    act_succ = Vector{Union{Nothing,Vector{Int32}}}[] # per action position, nothing if bad (or not enabled)
    leaf_of = zeros(Int32, fshape)
    box(id) = box_cells(cell[id], level_factor(L, level[id], nc), nc)

    function add_leaf!(l, c, layer, par)
        push!(level, l); push!(cell, c); push!(layer_of, layer); push!(parent, par)
        id = Int32(length(level))
        b = box(id)
        mask = UInt64(0)
        for x in b
            leaf_of[x, layer] = id
            mask |= ctx.start_mask[x, layer]
        end
        push!(active, true)
        push!(constraint, any(x -> unsafe[x, layer], b))
        push!(starts, mask != 0)
        push!(computed, false)
        for ai in 1:n_acts
            push!(enabled[ai], any(x -> can_act(ctx, layer, x, ai), b))
        end
        push!(flow_succ, nothing)
        push!(act_succ, Union{Nothing,Vector{Int32}}[nothing for _ in 1:n_acts])
        return id
    end

    # Blocks of level-0 leaves (see ForwardBlocks), about sweep_factor finest cells per axis
    f0 = level_factor(L, 0, nc)
    blocks = ForwardBlocks(ctx, ntuple(i -> f0[i] * max(1, round(Int, ctx.sweep[i] / f0[i])), D), fids, online)
    lap!(timer, "Blocks")
    # A leaf's successors lie within its parent's, and a level-0 leaf's within its block's (see successors_within)
    function compute_leaf!(id)
        layer, par = layer_of[id], parent[id]
        outer = par != 0 ? (flow_succ[par], act_succ[par], [enabled[ai][par] for ai in 1:n_acts]) :
                (bid = block_id(blocks, first(box(id)), layer); (blocks.flow[bid], blocks.acts[bid], blocks.enabled[bid]))
        flow_succ[id], act_succ[id] = successors_within(ctx, layer, box(id), fids, outer...)
        computed[id] = true
    end

    # The shield over the active leaves. Their successors (finest cells) are mapped to the leaves covering those cells and
    # kept in store: after a split only the new leaves, and those with a successor in a split leaf, are mapped again.
    store = SuccessorStore(0, [BitVector() for _ in 1:n_acts])
    store.flow_empty_ok = !ctx.exit_bad_for_flow
    store.act_empty_ok = !ctx.exit_bad_for_actions
    s = FlatShieldState(n_acts, act_needs_safe_flow)
    function grow_store!(n)
        while length(store.flow_none) < n
            push!(store.flow_none, false); push!(store.flow_start, 1); push!(store.flow_len, 0)
            for ai in 1:n_acts
                push!(store.act_enabled[ai], false); push!(store.act_start[ai], 1); push!(store.act_len[ai], 0)
            end
        end
    end
    # preds[v]: leaves that had v among their successors when they were mapped (some may since have been remapped)
    preds = Vector{Int32}[]
    function note_preds!(id, list)
        for v in list
            push!(preds[v], id)
        end
        return list
    end
    leaves_of(list) = unique!(Int32[leaf_of[v] for v in list])
    function map_leaf!(id)
        store.flow_none[id] = flow_succ[id] === nothing
        set_successors!(store, 0, id, flow_succ[id] === nothing ? Int32[] : note_preds!(id, leaves_of(flow_succ[id])))
        for ai in 1:n_acts
            a = act_succ[id][ai]
            store.act_enabled[ai][id] = enabled[ai][id] && a !== nothing
            # a bad action keeps no list here, but its finest cells may lie in a leaf that is split later: note it too
            leaves = a === nothing ? Int32[] : note_preds!(id, leaves_of(a))
            set_successors!(store, ai, id, store.act_enabled[ai][id] ? leaves : Int32[])
        end
    end
    function unmap_leaf!(id)
        for ai in 0:n_acts
            set_successors!(store, ai, id, Int32[])
        end
    end
    # Online: the blocks reaching each finest cell, and what has been queued so far (kept across rounds: a leaf computed
    # stays computed, and its children are computed when it is split).
    block_queued = falses(length(blocks.bad))
    fine_seen = falses(length(fids))
    pending = Int32[]
    function queue_block!(bid)
        block_queued[bid] && return
        block_queued[bid] = true
        bcells, layer = block_cells(blocks, bid, nc)
        for x in bcells
            lid = leaf_of[x, layer]
            (computed[lid] || constraint[lid] || !starts[lid]) || push!(pending, lid)
        end
    end
    # Queues the blocks that reach a finest cell of leaf id
    function queue_reaching!(id)
        layer = layer_of[id]
        for x in box(id)
            v = fids[x, layer]
            fine_seen[v] && continue
            fine_seen[v] = true
            foreach(queue_block!, blocks_reaching(blocks, v))
        end
    end

    # Maps the leaves remap, then runs the shield from scratch over the active leaves. Online it then computes, until
    # nothing changes, the leaves of the blocks that reach a leaf whose choices changed (and of bad blocks).
    # Returns what each leaf allows and lost, and how many leaves it computed.
    function shield(remap)
        local n_leaves = length(level)
        grow_store!(n_leaves)
        while length(preds) < n_leaves
            push!(preds, Int32[])
        end
        foreach(map_leaf!, remap)
        compact!(store)
        reset!(s, n_leaves)
        for id in Int32(1):Int32(n_leaves)
            if active[id] && constraint[id]
                seed_constraint!(s, id)
                store.flow_none[id] = true
            end
        end
        for id in Int32(1):Int32(n_leaves)
            active[id] && !constraint[id] && computed[id] && register_successors!(store, s, id)
        end
        propagate_bad!(s)
        local n_online = 0
        if online
            foreach(queue_block!, findall(blocks.bad))
            for id in Int32(1):Int32(n_leaves)
                active[id] && constraint[id] && queue_reaching!(id)
            end
            while true
                foreach(queue_reaching!, s.changed)
                empty!(s.changed)
                isempty(pending) && break
                lap!(timer, "Shield")
                sorted_unique!(pending) # a leaf can be queued by several blocks; none is computed yet
                for id in pending
                    compute_leaf!(id)
                    n_online += 1
                end
                lap!(timer, "Reachability")
                for id in pending
                    map_leaf!(id)
                    register_successors!(store, s, id)
                end
                empty!(pending)
                propagate_bad!(s)
            end
        end
        empty!(s.changed)
        local new_bits = zeros(T, n_leaves)
        local new_lost = zeros(T, n_leaves)
        for id in 1:n_leaves
            active[id] && ((new_bits[id], new_lost[id]) = choice_bits(T, s, id, (ai, k) -> enabled[ai][k], store_acted(store)))
        end
        return new_bits, new_lost, n_online
    end
    mapped(id) = active[id] && computed[id] && !constraint[id]
    # Everything a leaf could allow: waiting and every action it can take
    everything(id) = reduce(|, (T(1) << ai for ai in 1:n_acts if enabled[ai][id]); init=one(T))

    # Level 0: every leaf
    for layer in 1:ctx.n_layers, c in CartesianIndices(ntuple(i -> cld(nc[i], level_factor(L, 0, nc)[i]), D))
        add_leaf!(0, c, layer, Int32(0))
    end
    n_computed = 0
    for id in eachindex(level)
        (online || constraint[id] || !starts[id]) && continue
        compute_leaf!(id)
        n_computed += 1
    end
    lap!(timer, "Reachability")
    bits, lost, n_new = shield(filter(mapped, eachindex(level)))
    n_computed += n_new
    lap!(timer, "Shield")
    n_dead(bits) = count(id -> active[id] && bits[id] == 0, eachindex(level))
    verbose && println("Level 0: $(length(level)) leaves, $n_computed computed, $(n_dead(bits)) dead")

    rounds = 0
    neighbourhood = CartesianIndices(ntuple(_ -> -1:1, D))
    # A leaf on the boundary of the safe region that can still be split: it allows less than everything and borders a
    # leaf that allows something else.
    function on_boundary(id, bits)
        (active[id] && level[id] < L && starts[id] && bits[id] != everything(id)) || return false
        b, layer = box(id), layer_of[id]
        for x in b, off in neighbourhood
            y = x + off
            (checkbounds(Bool, leaf_of, y, layer) && !(y in b)) || continue
            bits[leaf_of[y, layer]] != bits[id] && return true
        end
        return false
    end
    # The active leaves next to leaf id (sharing a face or corner)
    function neighbours!(out, id)
        b, layer = box(id), layer_of[id]
        for x in b, off in neighbourhood
            y = x + off
            (checkbounds(Bool, leaf_of, y, layer) && !(y in b)) && push!(out, leaf_of[y, layer])
        end
        return out
    end
    to_check = collect(Int32(1):Int32(length(level))) # leaves whose boundary status may have changed
    while rounds < max_rounds
        candidates = Int32[]
        for id in unique!(to_check)
            on_boundary(id, bits) && push!(candidates, id)
        end
        isempty(candidates) && break
        rounds += 1
        n_new = 0
        n_old = length(level)
        # Leaves with a successor in a split leaf: their successors now belong to its children
        affected = Int32[]
        for id in candidates
            append!(affected, preds[id])
        end
        for id in candidates
            mapped(id) && unmap_leaf!(id)
            active[id] = false
            l, layer = level[id], layer_of[id]
            f = ntuple(i -> nc[i] == 1 ? 1 : 2, D)
            for ch in child_cells(cell[id], f)
                # children beyond the grid (when its size is not a multiple of the leaf) do not exist
                all(i -> (ch[i] - 1) * level_factor(L, l + 1, nc)[i] < nc[i], 1:D) || continue
                cid = add_leaf!(l + 1, ch, layer, id)
                (constraint[cid] || !starts[cid]) && continue
                compute_leaf!(cid)
                n_new += 1
            end
        end
        filter!(mapped, unique!(affected)) # the split leaves are inactive now
        n_computed += n_new
        lap!(timer, "Reachability")
        old_bits = bits
        bits, lost, n_new2 = shield(vcat(affected, filter(mapped, (n_old+1):length(level))))
        n_computed += n_new2
        lap!(timer, "Shield")
        # Only leaves whose choices changed, and their neighbours, can have become (or stopped being) boundary leaves
        empty!(to_check)
        for id in Int32(1):Int32(length(level))
            active[id] || continue
            if id > n_old || bits[id] != old_bits[id]
                push!(to_check, id)
                neighbours!(to_check, id)
            end
        end
        verbose && println("Round $rounds: split $(length(candidates)) leaves, computed $n_new, remapped $(length(affected)); active per level ",
            [count(id -> active[id] && level[id] == l, eachindex(level)) for l in 0:L], ", $(n_dead(bits)) dead")
    end

    allowed, lost_f = finest_choices(ctx, leaf_of, bits, lost, unsafe)
    levels = [count(id -> active[id] && level[id] == l, eachindex(level)) for l in 0:L]
    verbose && report(timer, "Refinement rounds: $rounds, leaves computed: $n_computed, active per level $levels, dead finest cells: $(count(==(0), allowed))")
    leaf_level = Int8[level[id] for id in leaf_of] # per finest cell, the level of the leaf covering it
    return grid, allowed, lost_f, (; ctx, levels, leaf_level)
end

# What each finest cell allows and lost: its leaf's choices. A cell that no location starts in is dead if it touches an
# unsafe set and allows waiting otherwise, as in ReACTed_reachable_cell_modal (it is never computed there). Where no location of a cell can take an action, trying it
# does nothing: allowed exactly when waiting is.
function finest_choices(ctx::ModalContext{D,N}, leaf_of::Array{Int32,N}, bits::Vector{T}, lost::Vector{T}, unsafe::BitArray{N}) where {D,N,T}
    allowed = zeros(T, size(leaf_of))
    lost_f = zeros(T, size(leaf_of))
    for k in CartesianIndices(leaf_of)
        c = CartesianIndex(ntuple(i -> k[i], D))
        locs = cell_locations(ctx, k[N], c)
        if isempty(locs)
            allowed[k] = unsafe[k] ? zero(T) : one(T)
            continue
        end
        id = leaf_of[k]
        a = bits[id]
        for ai in eachindex(ctx.action_ids)
            if all(li -> !ctx.may_act[li][ai][c], locs)
                a = isodd(a) ? a | (T(1) << ai) : a & ~(T(1) << ai)
            end
        end
        allowed[k] = a
        lost_f[k] = lost[id]
    end
    return allowed, lost_f
end
