using LazySets, LinearAlgebra, ReachabilityAnalysis, Dates
# Expects PartitionUtilities.jl, ReACTedShielding.jl, ReACTedShieldingPred.jl and ReACTedShieldingAdaptive.jl to be included first.
#
# ReACTed_reachable_cell_adaptive with the same exploration, refinement and shield, but faster:
#   - Action guard and constraint masks are built top-down: level 0 tests its cells, and a finer level only tests the
#     children of touching cells that are not inside the polyhedron. Cells are tested with dot products against the
#     halfspaces, and only cells cut by two or more halfspaces fall back to the LP of isdisjoint.
#   - Fast-cell covers use the direct facet enumeration of ConstantZonotopeCover(G, radius, Val(N)), and action
#     successors use one precomputed cover per level, action and target level instead of a zonotope per leaf.
#   - Leaves, their successor lists and the shield's predecessor lists are kept in flat arrays instead of a vector or
#     Dict entry per leaf, so the live heap holds few objects and garbage collection stays cheap.
#   - The result is one matrix on the finest grid with a bit per allowed choice instead of the successor relations.

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
        lo > H.b[k] && return 0
        hi > H.b[k] && (cut += 1)
    end
    cut == 0 && return 2
    cut == 1 && return 1
    box = Hyperrectangle(grid.lower .+ (car2vec(c) .- 0.5) .* g, fill(g / 2, N))
    return isdisjoint(H.P, box) ? 0 : 1
end

# Masks of the cells touching P on every level, the same cells as get_touching_cell_idxs_l per level.
function touching_masks(grids::Vector{Grid{N}}, P, factor::NTuple{N,Int}) where N
    H = CellHalfSpaces(P, Val(N))
    grid0 = grids[1]
    lo = clamp.(low(P), grid0.lower, grid0.upper)
    hi = clamp.(high(P), grid0.lower, grid0.upper)
    l = Int.(max.(floor.(abs.(lo .- grid0.lower) ./ grid0.granularity), 1))
    u = min.(grid0.numCells, Int.(ceil.(abs.(hi .- grid0.lower) ./ grid0.granularity)))
    touching = falses(cellshape(grid0))
    inside = falses(cellshape(grid0))
    for c in CartesianIndices(ntuple(i -> l[i]:max(u[i], l[i]), N))
        k = cell_meets(H, grid0, c)
        touching[c] = k > 0
        inside[c] = k == 2
    end
    masks = [touching]
    for grid in grids[2:end]
        touching′ = falses(cellshape(grid))
        inside′ = falses(cellshape(grid))
        for c in CartesianIndices(touching)
            touching[c] || continue
            for cc in child_cells(c, factor)
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

# Per-level reachability data as in adaptive_level, with the masks given.
function adaptive_level_fast(system::EuclideanHybridSystem, sh, grid::Grid{N}, L::Int, grid0::Grid{N}, sweep_factor::Int, may_act, constraint) where N
    g = grid.granularity
    (; δ⁻, δ⁺, A, P2A_abs, phiDict, tPhiDict, U, inputDict, max_input, max_flow) = sh

    sweep = [min(sweep_factor, n) for n in grid.numCells]
    nblocks = ntuple(i -> fld(grid.numCells[i], sweep[i]), N)
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
    fast = LazyFastCells(cellshape(grid), Tuple(sweep), nblocks, falses(nblocks), falses(nblocks), check)

    Z = remove_zero_generators(Zonotope(zeros(N), diagm(map(x -> x == 1 ? 0. : g/2, grid.numCells))))
    generatorDict = ReACT_discretize_decomposed_generators(Z, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict)
    discretizationDict, finaldims = ReACT_discretize_combine_with_offsets(Z, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict)

    reach_gens = hcat(genmat(max_input), genmat(linear_map(max_flow, Z)))
    fast_covers = map(0:L) do m
        G = hcat(reach_gens, diagm(fill(grid0.granularity / 2^m / 2, N)))
        ConstantZonotopeCover(G, ntuple(j -> sum(abs, view(G, j, :)), N), Val(N))
    end
    # The reach set of acting in a cell, max_flow * (J * cell + j) + max_input, is a translate of one zonotope per level
    # and action. These are the covers get_touching_cell_idxs_lz builds for it on level m: its nonzero generators in the
    # same order, bloated by half a level-m cell, with the radius of the unbloated set.
    act_covers = map(system.Act) do act
        G = hcat(max_flow * (act.jumpMatrix * genmat(Z)), genmat(max_input))
        G = G[:, [j for j in axes(G, 2) if any(!iszero, view(G, :, j))]]
        radius = ntuple(i -> sum(abs, view(G, i, :)), N)
        [ConstantZonotopeCover(hcat(G, Matrix((grid0.granularity / 2^m / 2) * I, N, N)), radius, Val(N)) for m in 0:L]
    end
    return AdaptiveLevelFast(grid, g, Z, generatorDict, discretizationDict, zeros(N), zeros(N, finaldims), zeros(N, finaldims),
        grid.lower .+ g / 2, fast, fast_covers, may_act, constraint, act_covers, zeros(N))
end

# Fields as in the levels of adaptive_level, plus the action covers.
struct AdaptiveLevelFast{N,F,GD,DD}
    grid::Grid{N}
    g::Float64
    Z::Zonotope{Float64,Vector{Float64},Matrix{Float64}}
    generatorDict::GD
    discretizationDict::DD
    pc::Vector{Float64}
    pg::Matrix{Float64}
    pig::Matrix{Float64}
    lower_offset::Vector{Float64}
    fast::LazyFastCells{N,F}
    fast_covers::Vector{ConstantZonotopeCover{N}} # per target level
    may_act::Vector{BitArray{N}}
    constraint::BitArray{N}
    act_covers::Vector{Vector{ConstantZonotopeCover{N}}} # per action, per target level
    act_center::Vector{Float64} # buffer
end

# ---------------------------------------------------------------------------------------------------------------------
# Leaves in flat storage. The successor lists of all leaves share one Int32 pool: the flow successors of leaf k are
# pool[flow_start[k] .+ (0:flow_len[k]-1)], and likewise per action. A recomputed list is appended and its old slice
# becomes garbage, which compact! removes. Each level maps its cells to active leaf ids with a dense array (0 = none),
# allocated when the level gets its first leaf.
mutable struct FlatLeaves{N}
    level::Vector{Int8}
    cell::Vector{CartesianIndex{N}}
    active::BitVector
    flow_none::BitVector # no flow successors at all (flow hit a constraint)
    act_enabled::Vector{BitVector} # per action: the leaf touches its guard
    ids::Vector{Array{Int32,N}} # per level
    shapes::Vector{NTuple{N,Int}} # per level
    n_active::Vector{Int} # per level
    pool::Vector{Int32}
    live::Int # pool entries still referenced
    flow_start::Vector{Int}
    flow_len::Vector{Int32}
    act_start::Vector{Vector{Int}} # per action
    act_len::Vector{Vector{Int32}} # per action
end

FlatLeaves{N}(shapes::Vector{NTuple{N,Int}}, n_acts::Int) where N = FlatLeaves{N}(Int8[], CartesianIndex{N}[], BitVector(), BitVector(),
    [BitVector() for _ in 1:n_acts], [zeros(Int32, ntuple(_ -> 0, N)) for _ in shapes], shapes, zeros(Int, length(shapes)),
    Int32[], 0, Int[], Int32[], [Int[] for _ in 1:n_acts], [Int32[] for _ in 1:n_acts])

function add_leaf!(leaves::FlatLeaves{N}, level::Int, cell::CartesianIndex{N}) where N
    push!(leaves.level, Int8(level))
    push!(leaves.cell, cell)
    push!(leaves.active, true)
    push!(leaves.flow_none, false)
    push!(leaves.flow_start, 1)
    push!(leaves.flow_len, 0)
    for ai in eachindex(leaves.act_enabled)
        push!(leaves.act_enabled[ai], false)
        push!(leaves.act_start[ai], 1)
        push!(leaves.act_len[ai], 0)
    end
    id = Int32(length(leaves.level))
    if isempty(leaves.ids[level+1])
        leaves.ids[level+1] = zeros(Int32, leaves.shapes[level+1])
    end
    leaves.ids[level+1][cell] = id
    leaves.n_active[level+1] += 1
    return id
end

function deactivate_leaf!(leaves::FlatLeaves, id::Int32)
    l = leaves.level[id]
    leaves.active[id] = false
    leaves.ids[l+1][leaves.cell[id]] = 0
    leaves.n_active[l+1] -= 1
    set_successors!(leaves, leaves.flow_start, leaves.flow_len, id, Int32[])
    for ai in eachindex(leaves.act_start)
        set_successors!(leaves, leaves.act_start[ai], leaves.act_len[ai], id, Int32[])
    end
end

function set_successors!(leaves::FlatLeaves, start::Vector{Int}, len::Vector{Int32}, id::Integer, succ::Vector{Int32})
    leaves.live += length(succ) - len[id]
    start[id] = length(leaves.pool) + 1
    len[id] = length(succ)
    append!(leaves.pool, succ)
end

# Rewrites the pool without garbage once garbage is most of it.
function compact!(leaves::FlatLeaves)
    length(leaves.pool) <= 2 * leaves.live + 2^20 && return
    pool = Vector{Int32}(undef, leaves.live)
    i = 1
    for (start, len) in Iterators.flatten(((leaves.flow_start, leaves.flow_len), zip(leaves.act_start, leaves.act_len)))
        for k in eachindex(start)
            copyto!(pool, i, leaves.pool, start[k], len[k])
            start[k] = i
            i += len[k]
        end
    end
    leaves.pool = pool
end

struct AdaptiveFlatContext{N,S,SH,LV}
    system::S
    shared::SH
    levels::Vector{LV}
    factor::NTuple{N,Int}
    leaves::FlatLeaves{N}
    L::Int
    buffer::Vector{Int32} # successor list being computed
end

function AdaptiveContextFast(system::EuclideanHybridSystem, p, max_granularity::Float64, L::Int, δ⁻::Float64, δ⁺::Float64, alg, maxOrder::Int, reduceOrder::Int, sweep_factor::Int)
    return _AdaptiveContextFast(Grid(system.statespace, max_granularity), system, p, L, δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor)
end

# Function barrier on Grid{N}
function _AdaptiveContextFast(grid0::Grid{N}, system, p, L, δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor) where N
    zero_width = low(system.statespace) .== high(system.statespace)
    factor = ntuple(i -> zero_width[i] ? 1 : 2, N)
    sh = adaptive_shared(system, p, δ⁻, δ⁺, alg, maxOrder, reduceOrder, Val(N))
    grids = [level_grid(grid0, l, factor) for l in 0:L]
    may_act = [touching_masks(grids, act.guard, factor) for act in system.Act] # per action, per level
    constraint = [falses(cellshape(grid)) for grid in grids]
    for gc in system.globalConstraints
        for (l, m) in enumerate(touching_masks(grids, gc, factor))
            constraint[l] .|= m
        end
    end
    levels = [adaptive_level_fast(system, sh, grids[l+1], L, grid0, sweep_factor, [m[l+1] for m in may_act], constraint[l+1]) for l in 0:L]
    leaves = FlatLeaves{N}([NTuple{N,Int}(Tuple(cellshape(grid))) for grid in grids], length(system.Act))
    return AdaptiveFlatContext(system, sh, levels, factor, leaves, L, Int32[])
end

is_constraint_leaf(ctx::AdaptiveFlatContext, id::Integer) = ctx.levels[ctx.leaves.level[id]+1].constraint[ctx.leaves.cell[id]]

# Appends the active leaves in ids whose cell centre lies in the cover translated to c. Visits the same cells, in the
# same order, as get_touching_cell_idxs_cover.
function cover_leaves!(out::Vector{Int32}, ids::Array{Int32,N}, grid::Grid{N}, cov::ConstantZonotopeCover{N}, c::AbstractVector{Float64}) where N
    lower_bounds = ntuple(i -> max(floor(Int, abs(clamp(c[i] - cov.radius[i], grid.lower[i], grid.upper[i]) - grid.lower[i]) / grid.granularity), 1), N)
    upper_bounds = ntuple(i -> min(grid.numCells[i], ceil(Int, abs(clamp(c[i] + cov.radius[i], grid.lower[i], grid.upper[i]) - grid.lower[i]) / grid.granularity)), N)
    sweeping_dim = argmax(upper_bounds[i] - lower_bounds[i] for i in 1:N)
    next_sweep = argmax(i == sweeping_dim ? 0 : upper_bounds[i] - lower_bounds[i] for i in 1:N)
    if upper_bounds[next_sweep] == lower_bounds[next_sweep]
        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:upper_bounds[i], N))
            id = ids[idx]
            id != 0 && push!(out, id)
        end
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
                id = ids[idx+step*off]
                id != 0 && push!(out, id)
            end
        end
    else
        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:upper_bounds[i], N))
            if cell_center_in_cover(cov, grid, idx, c)
                id = ids[idx]
                id != 0 && push!(out, id)
            end
        end
    end
    return out
end

# Successors of leaf id against the current active leaves, as compute_leaf!, stored in the pool.
function compute_leaf!(ctx::AdaptiveFlatContext{N}, id::Int32) where N
    leaves, sh, sys = ctx.leaves, ctx.shared, ctx.system
    l = Int(leaves.level[id])
    c = leaves.cell[id]
    lev = ctx.levels[l+1]
    (; δ⁻, δ⁺, A, P2A_abs, phiDict, tPhiDict, U, inputDict, max_input, max_flow) = sh
    buf = ctx.buffer

    empty!(buf)
    if lev.fast[c]
        fast_cell_center!(sh.center, max_flow, max_input.center, lev.grid, c)
        for m in 0:ctx.L
            leaves.n_active[m+1] == 0 || cover_leaves!(buf, leaves.ids[m+1], ctx.levels[m+1].grid, lev.fast_covers[m+1], sh.center)
        end
        leaves.flow_none[id] = false
    else
        of = lev.lower_offset + ((car2vec(c) .- 1) .* lev.g)
        LazySets.API.translate!(lev.Z, of)
        lev.pc .= 0
        lev.pg .= 0
        lev.pig .= 0
        ReACT_discretize_combine_with_offset_vector!(lev.pc, lev.pg, lev.pig, lev.discretizationDict, of, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, lev.generatorDict)
        z, f = propagate_set_t(lev.Z, [0.0, sh.p], δ⁻, δ⁺, sys, lev.grid, lev.fast, sys.edges, sh.identity_matrix, sh.attempt_recorder, N, P2A_abs, sh.guardProjVectors, sh.guardProjBounds, sh.invariantProjVectors, sh.invariantProjBounds, lev.discretizationDict, phiDict, tPhiDict, inputDict, max_input, max_flow, sh.guard_invariant_intersection, sh.alg, sh.maxOrder, sh.reduceOrder)
        LazySets.API.translate!(lev.Z, -of)
        leaves.flow_none[id] = f == 1
        if f != 1
            for m in 0:ctx.L
                leaves.n_active[m+1] == 0 && continue
                ids = leaves.ids[m+1]
                for cc in get_touching_cell_idxs_lz(ctx.levels[m+1].grid, z)
                    ids[cc] != 0 && push!(buf, ids[cc])
                end
            end
        end
    end
    set_successors!(leaves, leaves.flow_start, leaves.flow_len, id, buf)

    for (ai, act) in enumerate(sys.Act)
        empty!(buf)
        enabled = lev.may_act[ai][c]
        leaves.act_enabled[ai][id] = enabled
        if enabled
            # Centre in the order LazySets computes it: translate J * (cell centre) by j, map by max_flow, add max_input.
            of = lev.lower_offset + ((car2vec(c) .- 1) .* lev.g)
            lev.act_center .= max_flow * (act.jumpMatrix * of + act.jumpVector) + max_input.center
            for m in 0:ctx.L
                leaves.n_active[m+1] == 0 || cover_leaves!(buf, leaves.ids[m+1], ctx.levels[m+1].grid, lev.act_covers[ai][m+1], lev.act_center)
            end
        end
        set_successors!(leaves, leaves.act_start[ai], leaves.act_len[ai], id, buf)
    end
    return id
end

# ---------------------------------------------------------------------------------------------------------------------
# Event-driven shield over the leaves with the rules of OnlinePredShieldState (ReACTedShieldingPred.jl). Predecessor
# lists are linked lists in flat arrays: the predecessors of v start at edge head[v] and follow edge_next (0 ends).
mutable struct LeafShieldState
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

LeafShieldState(n_acts::Int, act_needs_safe_flow::Bool) = LeafShieldState(BitVector(), BitVector(), BitVector(), [BitVector() for _ in 1:n_acts],
    Int32[], Int32[], [Int32[] for _ in 1:n_acts], Int32[], Int32[], Int32[], Int32[], Int32[], act_needs_safe_flow)

# Clears the state for n leaves, keeping the allocated memory.
function reset!(s::LeafShieldState, n::Int)
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

@inline function add_pred!(s::LeafShieldState, head::Vector{Int32}, v::Int32, k::Int32)
    push!(s.edge_from, k)
    push!(s.edge_next, head[v])
    head[v] = length(s.edge_from)
end

@inline function check_dead!(s::LeafShieldState, k::Int32)
    if !s.dead[k] && s.computed[k] && s.flow_bad[k] && s.n_ok[k] == 0
        s.dead[k] = true
        push!(s.dead_queue, k)
        push!(s.changed, k)
    end
end

@inline function mark_flow_bad!(s::LeafShieldState, k::Int32)
    if !s.flow_bad[k]
        s.flow_bad[k] = true
        push!(s.flow_bad_queue, k)
        push!(s.changed, k)
        check_dead!(s, k)
    end
end

@inline function invalidate_action!(s::LeafShieldState, ai::Int, k::Int32)
    if s.act_ok[ai][k]
        s.act_ok[ai][k] = false
        s.n_ok[k] -= 1
        push!(s.changed, k)
        check_dead!(s, k)
    end
end

function seed_constraint_leaf!(s::LeafShieldState, k::Int32)
    s.dead[k] = true
    s.flow_bad[k] = true
    s.computed[k] = true
end

# register_cell! for a computed leaf, reading its successors from the pool.
function register_leaf!(ctx::AdaptiveFlatContext, s::LeafShieldState, k::Int32)
    leaves = ctx.leaves
    pool = leaves.pool
    start, len = leaves.flow_start[k], leaves.flow_len[k]
    if leaves.flow_none[k] || len == 0
        s.flow_bad[k] = true
    else
        for i in start:start+len-1
            v = pool[i]
            add_pred!(s, s.flow_head, v, k)
            s.dead[v] && (s.flow_bad[k] = true)
        end
    end
    n = Int32(0)
    for ai in eachindex(s.act_ok)
        start, len = leaves.act_start[ai][k], leaves.act_len[ai][k]
        (leaves.act_enabled[ai][k] && len > 0) || continue
        ok = true
        for i in start:start+len-1
            v = pool[i]
            if s.dead[v] || (s.act_needs_safe_flow && s.flow_bad[v])
                ok = false
                break
            end
        end
        ok || continue
        for i in start:start+len-1
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
function propagate_bad!(s::LeafShieldState)
    while !isempty(s.dead_queue) || !isempty(s.flow_bad_queue)
        while !isempty(s.dead_queue)
            v = pop!(s.dead_queue)
            e = s.flow_head[v]
            while e != 0
                mark_flow_bad!(s, s.edge_from[e])
                e = s.edge_next[e]
            end
            for ai in eachindex(s.act_head)
                e = s.act_head[ai][v]
                while e != 0
                    invalidate_action!(s, ai, s.edge_from[e])
                    e = s.edge_next[e]
                end
            end
        end
        while !isempty(s.flow_bad_queue)
            v = pop!(s.flow_bad_queue)
            (!s.act_needs_safe_flow || s.dead[v]) && continue
            for ai in eachindex(s.act_head)
                e = s.act_head[ai][v]
                while e != 0
                    invalidate_action!(s, ai, s.edge_from[e])
                    e = s.edge_next[e]
                end
            end
        end
    end
end

# neighbour_leaves with the dense cell-to-leaf maps.
function neighbour_leaves(ctx::AdaptiveFlatContext{N}, ids::Vector{Int32}) where N
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

function explore_leaves!(ctx::AdaptiveFlatContext, s::LeafShieldState, frontier::Vector{Int32})
    n = 0
    while !isempty(frontier)
        for id in frontier
            s.computed[id] && continue
            compute_leaf!(ctx, id)
            register_leaf!(ctx, s, id)
            n += 1
        end
        propagate_bad!(s)
        changed = unique!(s.changed)
        frontier = isempty(changed) ? Int32[] : filter(id -> !s.computed[id], neighbour_leaves(ctx, changed))
        empty!(s.changed)
    end
    return n
end

function init_leaves_online!(ctx::AdaptiveFlatContext, level::Int, include_border_cells::Bool, act_needs_safe_flow::Bool)
    leaves = ctx.leaves
    sys = ctx.system
    grid = ctx.levels[level+1].grid
    for c in CartesianIndices(cellshape(grid))
        add_leaf!(leaves, level, c)
    end
    s = reset!(LeafShieldState(length(sys.Act), act_needs_safe_flow), length(leaves.level))
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

# split_leaves_online! on flat leaves. The shield state is rebuilt in place.
function split_leaves_online!(ctx::AdaptiveFlatContext, s::LeafShieldState, to_split::AbstractVector{Int32})
    leaves = ctx.leaves
    pool = leaves.pool
    n_old = length(leaves.level)
    is_split = falses(n_old)
    for id in to_split
        @assert leaves.active[id] && s.computed[id] && leaves.level[id] < ctx.L
        is_split[id] = true
    end
    hits(start, len, k) = any(i -> is_split[pool[i]], start[k]:start[k]+len[k]-1)
    affected = Int32[]
    for k in Int32(1):Int32(n_old)
        (!leaves.active[k] || is_split[k] || !s.computed[k] || is_constraint_leaf(ctx, k)) && continue
        if hits(leaves.flow_start, leaves.flow_len, k) || any(ai -> hits(leaves.act_start[ai], leaves.act_len[ai], k), eachindex(leaves.act_start))
            push!(affected, k)
        end
    end
    for id in to_split
        l = Int(leaves.level[id])
        c = leaves.cell[id]
        deactivate_leaf!(leaves, id)
        for ch in child_cells(c, ctx.factor)
            add_leaf!(leaves, l + 1, ch)
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
    compact!(leaves)
    reset!(s, length(leaves.level))
    for id in Int32(1):Int32(length(leaves.level))
        (leaves.active[id] && computed[id] && is_constraint_leaf(ctx, id)) && seed_constraint_leaf!(s, id)
    end
    for id in Int32(1):Int32(length(leaves.level))
        (leaves.active[id] && computed[id] && !is_constraint_leaf(ctx, id)) && register_leaf!(ctx, s, id)
    end
    propagate_bad!(s)
    changed = copy(unique!(s.changed))
    empty!(s.changed)
    frontier = isempty(changed) ? Int32[] : filter(id -> !s.computed[id], neighbour_leaves(ctx, changed))
    n += explore_leaves!(ctx, s, frontier)
    return s, n
end

function dead_volume(ctx::AdaptiveFlatContext{N}, dead::AbstractVector{Bool}) where N
    leaves = ctx.leaves
    v = 0
    for id in eachindex(leaves.level)
        if leaves.active[id] && dead[id]
            v += prod(ntuple(i -> ctx.factor[i]^(ctx.L - leaves.level[id]), N))
        end
    end
    return v
end

# Predecessors as one flat array: the predecessors of v are list[start[v]:start[v+1]-1], in increasing order.
function flat_predecessors(leaves::FlatLeaves, starts::Vector{Int}, lens::Vector{Int32}, usable)
    n = length(leaves.level)
    pool = leaves.pool
    start = zeros(Int, n + 1)
    for k in 1:n
        (leaves.active[k] && usable(k)) || continue
        for i in starts[k]:starts[k]+lens[k]-1
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
        for i in starts[k]:starts[k]+lens[k]-1
            v = pool[i]
            list[next[v]] = k
            next[v] += 1
        end
    end
    return start, list
end

# adaptive_shield on flat leaves: same rules and processing order, so the same dead leaves and the same death order,
# which the refinement uses to pick the latest-dying leaves.
function adaptive_shield(leaves::FlatLeaves; act_needs_safe_flow::Bool=true, computed::AbstractVector{Bool})
    n = length(leaves.level)
    n_acts = length(leaves.act_enabled)
    flow_start, flow_list = flat_predecessors(leaves, leaves.flow_start, leaves.flow_len, k -> true)
    act_preds = [flat_predecessors(leaves, leaves.act_start[ai], leaves.act_len[ai], k -> leaves.act_enabled[ai][k]) for ai in 1:n_acts]
    act_ok = [falses(n) for _ in 1:n_acts]
    n_ok = zeros(Int32, n)
    for k in 1:n, ai in 1:n_acts
        if leaves.active[k] && leaves.act_enabled[ai][k] && leaves.act_len[ai][k] > 0
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
    function invalidate!(ai, k)
        if act_ok[ai][k]
            act_ok[ai][k] = false
            n_ok[k] -= 1
            try_kill!(k)
        end
    end
    for k in Int32(1):Int32(n)
        if leaves.active[k] && computed[k] && (leaves.flow_none[k] || leaves.flow_len[k] == 0)
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
            for ai in 1:n_acts
                start, list = act_preds[ai]
                for i in start[v]:start[v+1]-1
                    invalidate!(ai, list[i])
                end
            end
        end
        while hf <= length(queue_fb)
            v = queue_fb[hf]
            hf += 1
            (!act_needs_safe_flow || dead[v]) && continue
            for ai in 1:n_acts
                start, list = act_preds[ai]
                for i in start[v]:start[v+1]-1
                    invalidate!(ai, list[i])
                end
            end
        end
    end
    return dead, order
end

# Allowed choices of each finest cell as bits: bit 0 (value 1) is waiting a period, bit ai (value 2^ai) is the action
# system.Act[ai]. 0 means dead. Leaves that were never computed allow waiting and every action whose guard they touch.
function adaptive_allowed(ctx::AdaptiveFlatContext{N}, s::LeafShieldState) where N
    n_acts = length(ctx.system.Act)
    n_acts < 64 || throw(ArgumentError("at most 63 actions fit in the allowed matrix"))
    T = n_acts < 8 ? UInt8 : n_acts < 16 ? UInt16 : n_acts < 32 ? UInt32 : UInt64
    return _adaptive_allowed!(zeros(T, cellshape(ctx.levels[end].grid)), ctx, s)
end

# Function barrier on the element type of allowed
function _adaptive_allowed!(allowed::Array{T,N}, ctx::AdaptiveFlatContext{N}, s::LeafShieldState) where {T,N}
    leaves = ctx.leaves
    n_acts = length(ctx.system.Act)
    for id in Int32(1):Int32(length(leaves.level))
        leaves.active[id] || continue
        bits = zero(T)
        if s.computed[id]
            if !s.dead[id]
                s.flow_bad[id] || (bits |= one(T))
                for ai in 1:n_acts
                    s.act_ok[ai][id] && (bits |= T(1) << ai)
                end
            end
        else
            lev = ctx.levels[leaves.level[id]+1]
            bits = one(T)
            for ai in 1:n_acts
                lev.may_act[ai][leaves.cell[id]] && (bits |= T(1) << ai)
            end
        end
        e = ctx.L - leaves.level[id]
        c = leaves.cell[id]
        for cc in CartesianIndices(ntuple(i -> ((c[i]-1)*ctx.factor[i]^e+1):(c[i]*ctx.factor[i]^e), N))
            allowed[cc] = bits
        end
    end
    return allowed
end

# ReACTed_reachable_cell_adaptive with the faster setup and flat storage. Returns the finest grid (deadCells set where
# nothing is allowed), the allowed matrix (see adaptive_allowed) and the context. on_round(ctx, dead, round,
# leaf_computations) is called after every round, as in ReACTed_reachable_cell_adaptive.
function ReACTed_reachable_cell_adaptive_fast(system::EuclideanHybridSystem, p, max_granularity, granularity_steps::Int64, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, max_steps::Int64=1000; include_border_cells::Bool=false, refine_fraction::Float64=0.4, patience::Int=1, act_needs_safe_flow::Bool=true, sweep_factor::Int=20, on_round=(ctx, dead, r, computed) -> nothing)
    recorderTimeStart = time_ns()
    ctx = AdaptiveContextFast(system, p, Float64(max_granularity), granularity_steps, δ⁻, δ⁺, alg, maxOrder, reduceOrder, sweep_factor)
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

    allowed = adaptive_allowed(ctx, s)
    fine = ctx.levels[end].grid
    grid = Grid(length(fine.numCells), fine.granularity, copy(fine.lower), copy(fine.upper), copy(fine.numCells), allowed .== 0)
    println("Total times:
        Total: $((time_ns() - recorderTimeStart)/10^9)
        Setup: $((setup_time - recorderTimeStart) / 10^9)
        Reachability: $(reachability_time / 10^9)
        Shield: $(shield_time / 10^9)")
    println("Refinement rounds: $rounds, leaf computations: $computed")
    return grid, allowed, ctx
end
