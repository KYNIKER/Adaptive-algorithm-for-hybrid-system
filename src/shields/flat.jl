# Building blocks shared by the shields. Each gives every cell (modal) or leaf (modal_adaptive) an Int32 id and keeps what it knows per id in flat arrays.

child_cells(c::CartesianIndex{N}, factor::NTuple{N,Int}) where N =
    CartesianIndices(ntuple(i -> ((c[i]-1)*factor[i]+1):(c[i]*factor[i]), N))


# Wrapper to group HPolyhedron and their constraints
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

# How the cell c meets the polyhedron H.P: 0 misses, 1 touches, 2 inside.
# A cell within CELL_MEETS_TOL (relative) of a halfspace still touches it, as for the LazySets.isdisjoint function.
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
    box = Hyperrectangle(grid.lower .+ (index_vector(c) .- 0.5) .* g, fill(g / 2, N))
    return isdisjoint(H.P, box) ? 0 : 1
end

# Cell ids for cells in box of P on grid
function touching_candidates(grid::Grid{N}, P) where N
    lo, hi = low(P), high(P)
    return CartesianIndices(ntuple(i -> cell_range(grid, i, lo[i], hi[i]), N))
end

# Splits the cells intersecting a guard into those touching and those inside.
guard_masks(grid::Grid, ::Nothing) = (trues(cellshape(grid)), trues(cellshape(grid)))
function guard_masks(grid::Grid{N}, P) where N
    H = CellHalfSpaces(P, Val(N))
    touching = falses(cellshape(grid))
    inside = falses(cellshape(grid))
    for c in touching_candidates(grid, P)
        k = cell_meets(H, grid, c)
        touching[c] = k > 0
        inside[c] = k == 2
    end
    return touching, inside
end


# Struct for storing successors in one contigous vector. This eliminates the need for dictionary look ups, but requires more book keeping. A recomputed list is appended, which compact! removes.
mutable struct SuccessorStore
    pool::Vector{Int32}
    live::Int # pool entries still referenced
    flow_none::BitVector # the flow is bad (it hit a constraint or left the state space), so it has no successors
    flow_start::Vector{Int}
    flow_len::Vector{Int32}
    act_enabled::Vector{BitVector} # per action: the id touches its guard
    act_start::Vector{Vector{Int}}
    act_len::Vector{Vector{Int32}}
    # Whether waiting, or an action, with no successors is safe: everything it reaches left the state space, and leaving
    # is allowed. Otherwise no successors is bad.
    flow_empty_ok::Bool
    act_empty_ok::Bool
end

SuccessorStore(n::Int, act_enabled::Vector{BitVector}) = SuccessorStore(Int32[], 0, falses(n), ones(Int, n), zeros(Int32, n),
    act_enabled, [ones(Int, n) for _ in act_enabled], [zeros(Int32, n) for _ in act_enabled], false, false)

successor_range(store::SuccessorStore, ai::Int, k::Integer) =
    ai == 0 ? (store.flow_start[k]:(store.flow_start[k]+store.flow_len[k]-1)) : (store.act_start[ai][k]:(store.act_start[ai][k]+store.act_len[ai][k]-1))

function set_successors!(store::SuccessorStore, ai::Int, id::Integer, succ::AbstractVector{Int32})
    start, len = ai == 0 ? (store.flow_start, store.flow_len) : (store.act_start[ai], store.act_len[ai])
    store.live += length(succ) - len[id]
    start[id] = length(store.pool) + 1
    len[id] = length(succ)
    append!(store.pool, succ)
end

function compact!(store::SuccessorStore)
    length(store.pool) <= 2 * store.live + 2^20 && return
    pool = Vector{Int32}(undef, store.live)
    i = 1
    for (start, len) in Iterators.flatten((((store.flow_start, store.flow_len),), zip(store.act_start, store.act_len)))
        for k in eachindex(start)
            copyto!(pool, i, store.pool, start[k], len[k])
            start[k] = i
            i += len[k]
        end
    end
    store.pool = pool
end


# Struct for storing the Shield state, that is what is the state of all cells computed. Most cell states are stored as, or vectors of, bit vectors.
# It also stores the predecessor relation for both flow and actions, as edges split into from and next. Flow_ and act_head is used to help realize this.
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

function seed_constraint!(s::FlatShieldState, k::Int32)
    s.dead[k] = true
    s.flow_bad[k] = true
    s.computed[k] = true
end

# Make the successor store and flat shield state agree. 
function register_successors!(store::SuccessorStore, s::FlatShieldState, k::Int32)
    pool = store.pool
    if store.flow_none[k] || (store.flow_len[k] == 0 && !store.flow_empty_ok)
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
        (store.act_enabled[ai][k] && (!isempty(range) || store.act_empty_ok)) || continue
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

# Lists as one pool: list b holds pool[off[b]:off[b+1]-1], nothing for an empty one.
function to_csr(lists)
    off, pool = ones(Int, length(lists) + 1), Int32[]
    for (b, list) in enumerate(lists)
        list === nothing || append!(pool, list)
        off[b+1] = length(pool) + 1
    end
    return off, pool
end

# The allowed choices of a shield that computes every id once, as the bits of choice_bits: a lean FlatShieldState that
# keeps only what they need. flow[id] holds the successors of waiting and acts[ai][id] those of action ai, nothing where
# the choice is bad (for an action, or not enabled), for the computed ids. dead holds the ids dead from the start, and
# is updated. A choice without successors is safe only where flow_empty_ok or act_empty_ok says so. An id never computed
# allows waiting. Bad events are pushed backwards along the predecessors, each edge at most once per kind.
function allowed_choices(::Type{T}, flow::Vector{Union{Nothing,Vector{Int32}}}, acts::Vector{Vector{Union{Nothing,Vector{Int32}}}},
    computed::BitVector, dead::BitVector, act_needs_safe_flow::Bool, flow_empty_ok::Bool, act_empty_ok::Bool) where T
    n, n_acts = length(flow), length(acts)
    flow_off, flow_pred = reverse_csr(to_csr(flow)..., n)
    act_off, act_pred = reverse_csr(to_csr(collect(Iterators.flatten(acts)))..., n) # list (ai - 1) n + k is action ai of k

    flow_bad = copy(dead)
    for k in findall(computed)
        succ = flow[k]
        flow_bad[k] = succ === nothing || (isempty(succ) && !flow_empty_ok) || any(v -> dead[v], succ)
    end
    ok = falses(n_acts * n)
    n_ok = zeros(Int, n)
    for ai in 1:n_acts, k in findall(computed)
        succ = acts[ai][k]
        safe = succ !== nothing && (!isempty(succ) || act_empty_ok) && !any(v -> dead[v] || (act_needs_safe_flow && flow_bad[v]), succ)
        safe && (ok[(ai-1)*n+k] = true; n_ok[k] += 1)
    end

    dead_queue, bad_queue = findall(dead), findall(flow_bad)
    die!(k) = computed[k] && !dead[k] && flow_bad[k] && n_ok[k] == 0 && (dead[k] = true; push!(dead_queue, k))
    foreach(die!, 1:n)
    function unsafe!(b) # action edge b leads to a bad id
        ok[b] || return
        ok[b] = false
        k = (b - 1) % n + 1
        n_ok[k] -= 1
        die!(k)
    end
    while !isempty(dead_queue) || !isempty(bad_queue)
        if !isempty(dead_queue)
            v = pop!(dead_queue)
            for k in view(flow_pred, flow_off[v]:flow_off[v+1]-1)
                flow_bad[k] || (flow_bad[k] = true; push!(bad_queue, k); die!(k))
            end
            foreach(unsafe!, view(act_pred, act_off[v]:act_off[v+1]-1))
        else
            v = pop!(bad_queue)
            act_needs_safe_flow && foreach(unsafe!, view(act_pred, act_off[v]:act_off[v+1]-1))
        end
    end

    allowed = zeros(T, n)
    for k in 1:n
        dead[k] && continue
        allowed[k] = computed[k] && flow_bad[k] ? zero(T) : one(T)
        for ai in 1:n_acts
            ok[(ai-1)*n+k] && (allowed[k] |= one(T) << ai)
        end
    end
    return allowed
end

# ---------------------------------------------------------------------------------------------------------------------
# Results. Every shield returns grid, allowed, lost, info: allowed and lost hold the choice bits of each cell (see
# choice_bits), grid.deadCells is where nothing is allowed, and info is a NamedTuple of what only that shield has.

function choice_type(n_acts::Int)
    n_acts < 64 || throw(ArgumentError("at most 63 actions fit in the allowed matrix"))
    return n_acts < 8 ? UInt8 : n_acts < 16 ? UInt16 : n_acts < 32 ? UInt32 : UInt64
end

# What k (a cell or id) allows and what it lost, as bits: bit 0 (value 1) is waiting a period, bit ai (value 2^ai) is
# the action at position ai; nothing allowed means dead. s is any shield state with the fields dead, flow_bad, computed
# and act_ok. A k never computed allows waiting and every action it is enabled for (it touches the guard). An action is
# lost in an alive k that is enabled for it and acted (has successors for it), but where it is unsafe.
# enabled(ai, k) and acted(ai, k) say whether k is enabled for and has successors for action ai. Where an action does
# nothing, the shields set it to the waiting bit afterwards (see noop_choices! in modal.jl).
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

store_acted(store::SuccessorStore) = (ai, k) -> store.act_len[ai][k] > 0

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
