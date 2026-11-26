# modal.jl reduced to its method: all axes alike, every cell computed once, bad events pushed back (flat.jl).
# Successors are the touching_cells of each end set, within its Bounds (modal_shared.jl).

const PeriodEnd = Tuple{Int,Zonotope{Float64,Vector{Float64},Matrix{Float64}},Bounds} # location, set, bounds

struct SimpleLocation{D}
    dynamics::Dynamics
    layer::Int
    guards::Vector{Halfspaces}
    rate::Vector{Float64} # rate of each axis that only translates, NaN otherwise
    starts::BitArray{D} # cells that start a period here
    discs::Vector{Dict{Tuple{Matrix{Float64},Float64},Union{Nothing,CellDiscretization{D}}}} # per depth, per (shape, δ⁺); nothing: seen once
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
    layers::Vector{Vector{Int}}
    unsafe::Vector{Halfspaces}
    statespace::Halfspaces
    cell::Matrix{Float64} # generators of a cell
    action_ids::Vector{Int}
    may_act::Vector{Vector{BitArray{D}}} # [location][action]: some state of the cell is in the guard
    must_act::Vector{Vector{BitArray{D}}} # [location][action]: every state of the cell is in the guard
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
    half = [n == 1 ? (h - l) / 2 : g / 2 for (n, l, h) in zip(grid.numCells, low(system.statespace), high(system.statespace))] # one cell: the whole axis
    cell = Matrix{Float64}(genmat(remove_zero_generators(Zonotope(zeros(D), diagm(half)))))

    shared = Dict{Tuple{Matrix{Float64},Vector{Float64},Matrix{Float64}},Dynamics}()
    locations = map(system.locations) do loc
        dyn = get!(() -> Dynamics(loc.A, loc.U, steps, δ⁻, δ⁺, alg, maxOrder, reduceOrder), shared, (loc.A, LazySets.center(loc.U), genmat(loc.U)))
        region = loc.region === nothing ? nothing : getHalfSpaceProjections(loc.region)
        starts = BitArray(loc.starts && box_overlaps_region(region, cell_low(grid, c), cell_low(grid, c) .+ g) for c in CartesianIndices(cellshape(grid)))
        rate = [iszero(view(loc.A, i, :)) && iszero(view(genmat(loc.U), i, :)) ? LazySets.center(loc.U)[i] : NaN for i in 1:D]
        SimpleLocation(dyn, loc.layer, Halfspaces[getHalfSpaceProjections(e.guard) for e in loc.edges], rate, starts,
            [Dict{Tuple{Matrix{Float64},Float64},Union{Nothing,CellDiscretization{D}}}() for _ in 0:max_depth])
    end
    layers = [[li for (li, loc) in enumerate(locations) if loc.layer == layer] for layer in 1:maximum(loc.layer for loc in locations)]

    action_ids = sort!(unique(a.action.id for loc in system.locations for a in loc.actions))
    may_act = [[falses(cellshape(grid)) for _ in action_ids] for _ in locations]
    must_act = deepcopy(may_act)
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
starting(ctx::SimpleContext, c::CartesianIndex, layer::Int) = filter(li -> ctx.locations[li].starts[c], ctx.layers[layer])
can_act(ctx::SimpleContext, c::CartesianIndex, layer::Int, ai::Int) = any(li -> ctx.may_act[li][ai][c], starting(ctx, c, layer))
flowed(dyn::Dynamics, X, k::Int) = minkowski_sum(linear_map(dyn.Φ[k+1], X), dyn.V[k+1])

# Built on a shape's second use; one cache per depth, as discretize_cell! overwrites what a caller still reads.
function discretize(ctx::SimpleContext, li::Int, X, depth::Int, δ⁺::Float64)
    dyn, cache, key = ctx.locations[li].dynamics, ctx.locations[li].discs[depth+1], (genmat(X), δ⁺)
    cd = get(cache, key, false)
    cd === false && (cache[key] = nothing; return newReACTDiscretizePlus(X, ctx.δ⁻, δ⁺, dyn.A, dyn.P2A_abs, dyn.step_Φ, dyn.step_input, ctx.alg, ctx.maxOrder, ctx.reduceOrder))
    cd === nothing && (cd = cache[key] = period_discretization(dyn, ctx.system.locations[li].U, key[1], ctx.δ⁻, δ⁺))
    return discretize_cell!(cd, LazySets.center(X))
end

# Appends the end sets of X, which has D to D + W steps of δ⁻ left; false if the transition is bad.
function reach!(out::Vector{PeriodEnd}, ctx::SimpleContext, li::Int, X, D::Int, W::Int, exit_bad::Bool, bounds::Bounds, depth::Int=0)
    depth > ctx.max_depth && return false # Zeno
    total = D + W
    total == 0 && (push!(out, (li, X, bounds)); return true)
    loc, dyn = ctx.locations[li], ctx.locations[li].dynamics
    δ⁺ = min(ctx.δ⁺, ctx.δ⁻ * 2^floor(Int, log2(total)))
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
                D′ = max(D - k - 1, 0) # the piece's states have spent k to k + 1 steps here
                W′ = total - k - D′
                b = jumped_bounds(bounds, guard, e.jumpMatrix, e.jumpVector, loc.rate, ctx.locations[e.targetLoc].rate, D′, W′, ctx.δ⁻, k == 0, ctx.grid.granularity)
                b === nothing && continue
                jumped = jump_set(clipped, e.jumpMatrix, e.jumpVector, ctx.maxOrder, ctx.reduceOrder)
                reach!(out, ctx, e.targetLoc, jumped, D′, W′, exit_bad, b, depth + 1) || return false
                all_jumped |= all_inside(piece, guard...)
            end
            all_jumped && return true
        end
        W > 0 && k >= D && push!(out, (li, piece, bounds)) # states with k + 1 steps left end here
    end
    exit_bad && any(k -> leaves(ctx, loc, flowed(dyn, X, k), bounds), unique((D, total))) && return false
    W == 0 && push!(out, (li, flowed(dyn, X, total), bounds))
    return true
end

# Whether part of the set is beyond a face of the state space and in no guard; the guard's own face halfspace is skipped.
function leaves(ctx::SimpleContext, loc::SimpleLocation, set, (lo, hi)::Bounds)
    for (a, b) in zip(ctx.statespace...)
        i = findfirst(!iszero, a)
        bound = a[i] * (a[i] > 0 ? hi[i] : lo[i])
        min(bound, ρ(a, set)) > b + MODAL_TOL || continue
        beyond = clip_to_polyhedron(set, ([-a], [-b]))
        beyond === nothing && (isfinite(bound) ? (return true) : continue) # exact bounds beat an inexact clip
        in_guard(guard) = all(((c, d),) -> (isapprox(c, -a; atol=MODAL_TOL) && isapprox(d, -b; atol=MODAL_TOL)) || ρ(c, beyond) <= d, zip(guard...))
        any(in_guard, loc.guards) || return true
    end
    return false
end

# Cells reached by the transitions (location, set, exit_bad, affine map of Z), or nothing if one is bad.
function reached(ctx::SimpleContext, transitions, Z)
    out = PeriodEnd[]
    for (li, X, exit_bad, affine) in transitions
        reach!(out, ctx, li, X, ctx.steps, 0, exit_bad, end_bounds(ctx, ctx.locations[li].dynamics, affine, Z)) || return nothing
    end
    return unique!([CartesianIndex(c, ctx.locations[li].layer) for (li, X, b) in out for c in touching_cells(ctx.grid, X, ctx.covers, cell_bounds(ctx.grid, b)...)])
end

# Successors of waiting and per action (nothing: disabled or bad); where an action is not taken, the state just waits.
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

# Offline shield like ReACTed_reachable_cell_modal; returns grid, allowed, nothing (no lost), (; ctx).
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
    split(k) = CartesianIndex(ntuple(i -> k[i], D)), k[D+1]
    as_ids(succ) = succ === nothing ? nothing : Int32[ids[k] for k in succ]

    computed, dead = falses(n), falses(n)
    flow = Vector{Union{Nothing,Vector{Int32}}}(nothing, n)
    acts = [Vector{Union{Nothing,Vector{Int32}}}(nothing, n) for _ in 1:n_acts]
    for H in ctx.system.unsafe, c in touching_cells(grid, H), layer in eachindex(ctx.layers)
        dead[ids[c, layer]] = true
    end
    for k in cells
        id = ids[k]
        c, layer = split(k)
        locs = starting(ctx, c, layer)
        (dead[id] || isempty(locs)) && continue
        waiting, act_succ = cell_successors(ctx, c, locs)
        flow[id] = as_ids(waiting)
        foreach(ai -> acts[ai][id] = as_ids(act_succ[ai]), 1:n_acts)
        computed[id] = true
        verbose && id % 100_000 == 0 && println("$id of $n cells : $(Dates.format(now(), "HH:MM:SS"))")
    end
    allowed = reshape(allowed_choices(choice_type(n_acts), flow, acts, computed, dead, act_needs_safe_flow, !ctx.exit_bad_for_flow, !ctx.exit_bad_for_actions), shape)
    verbose && println("$n cells, dead: $(count(dead))")
    return grid, noop_allowed!(allowed, ctx), nothing, (; ctx)
end

# An action no location can take does nothing: allowed where waiting is.
function noop_allowed!(allowed::Array{T}, ctx::SimpleContext{D}) where {T,D}
    for k in CartesianIndices(allowed)
        c, layer = CartesianIndex(ntuple(i -> k[i], D)), k[D+1]
        (isodd(allowed[k]) && !isempty(starting(ctx, c, layer))) || continue
        for ai in eachindex(ctx.action_ids)
            can_act(ctx, c, layer, ai) || (allowed[k] |= one(T) << ai)
        end
    end
    return allowed
end
