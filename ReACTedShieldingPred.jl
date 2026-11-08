using LazySets, LinearAlgebra, ReachabilityAnalysis, Dates
# Expects ReACTedShielding.jl and PartitionUtilities.jl to be included first.

# Event-driven (counter based) shield state for the online algorithm.
# Instead of re-checking every computed cell each round (shield_step!), we keep a predecessor relation that grows as
# successors are computed, and push "bad" events backwards along it. Every edge is processed at most once.
#
#   flow_bad[k]      : not acting in k is unsafe (no flow successors, hit constraint, or a flow successor is dead)
#   act_ok[ai][k]    : action Act[ai] is still a safe choice in k
#   n_ok[k]          : number of actions still safe in k
#   k is dead        : flow_bad[k] && n_ok[k] == 0
#
# If act_needs_safe_flow is true an action is also invalidated when one of its successors becomes flow_bad (same rule as
# shield_step!), i.e. after acting the system must be able to wait one period.
mutable struct OnlinePredShieldState{N}
    dead::BitArray{N}
    flow_bad::BitArray{N}
    computed::BitArray{N}
    act_ok::Vector{BitArray{N}}
    n_ok::Array{Int32,N}
    pred_flow::FlowRelation{N}
    pred_act::Vector{FlowRelation{N}} # indexed by action position, not act.id
    dead_queue::Vector{CartesianIndex{N}}
    flow_bad_queue::Vector{CartesianIndex{N}}
    changed::Vector{CartesianIndex{N}} # cells whose status changed this round, used to grow the frontier
    act_needs_safe_flow::Bool
end

function OnlinePredShieldState(deadCells::BitArray{N}, n_acts::Int, act_needs_safe_flow::Bool) where N
    dims = size(deadCells)
    return OnlinePredShieldState{N}(deadCells, falses(dims...), falses(dims...), [falses(dims...) for _ in 1:n_acts], zeros(Int32, dims...),
        FlowRelation{N}(), [FlowRelation{N}() for _ in 1:n_acts],
        CartesianIndex{N}[], CartesianIndex{N}[], CartesianIndex{N}[], act_needs_safe_flow)
end

@inline function add_pred!(d::FlowRelation{N}, v::CartesianIndex{N}, k::CartesianIndex{N}) where N
    push!(get!(() -> CartesianIndex{N}[], d, v), k)
end

@inline function check_dead!(s::OnlinePredShieldState, k)
    if !s.dead[k] && s.computed[k] && s.flow_bad[k] && s.n_ok[k] == 0
        s.dead[k] = true
        push!(s.dead_queue, k)
        push!(s.changed, k)
    end
end

@inline function mark_flow_bad!(s::OnlinePredShieldState, k)
    if !s.flow_bad[k]
        s.flow_bad[k] = true
        push!(s.flow_bad_queue, k)
        push!(s.changed, k)
        check_dead!(s, k)
    end
end

@inline function invalidate_action!(s::OnlinePredShieldState, ai::Int, k)
    if s.act_ok[ai][k]
        s.act_ok[ai][k] = false
        s.n_ok[k] -= 1
        push!(s.changed, k)
        check_dead!(s, k)
    end
end

@inline action_successor_bad(s::OnlinePredShieldState, v) = s.dead[v] || (s.act_needs_safe_flow && s.flow_bad[v])

# Registers the successors of a freshly computed cell k. flow_succ === nothing means the flow hit a constraint.
# act_succs[ai] === nothing means action ai is not enabled in k.
function register_cell!(s::OnlinePredShieldState{N}, k::CartesianIndex{N}, flow_succ, act_succs) where N
    if flow_succ === nothing || isempty(flow_succ)
        s.flow_bad[k] = true
    else
        for v in flow_succ
            add_pred!(s.pred_flow, v, k)
            if s.dead[v]
                s.flow_bad[k] = true
            end
        end
    end

    n = Int32(0)
    for (ai, succ) in enumerate(act_succs)
        if succ === nothing || isempty(succ)
            continue
        end
        ok = true
        for v in succ
            if action_successor_bad(s, v)
                ok = false
                break
            end
        end
        if ok
            for v in succ
                add_pred!(s.pred_act[ai], v, k)
            end
            s.act_ok[ai][k] = true
            n += 1
        end
    end
    s.n_ok[k] = n
    s.computed[k] = true
    # Status of k is new information for cells that already point at k.
    if s.flow_bad[k]
        push!(s.flow_bad_queue, k)
        push!(s.changed, k)
    end
    check_dead!(s, k)
end

# Pushes all pending events backwards through the known predecessor relation until nothing changes.
function propagate_bad!(s::OnlinePredShieldState)
    while !isempty(s.dead_queue) || !isempty(s.flow_bad_queue)
        while !isempty(s.dead_queue)
            v = pop!(s.dead_queue)
            preds = get(s.pred_flow, v, nothing)
            if preds !== nothing
                for k in preds
                    mark_flow_bad!(s, k)
                end
            end
            for ai in eachindex(s.pred_act)
                apreds = get(s.pred_act[ai], v, nothing)
                if apreds !== nothing
                    for k in apreds
                        invalidate_action!(s, ai, k)
                    end
                end
            end
        end
        while !isempty(s.flow_bad_queue)
            v = pop!(s.flow_bad_queue)
            # A flow_bad cell that is dead has already been handled through dead_queue.
            if !s.act_needs_safe_flow || s.dead[v]
                continue
            end
            for ai in eachindex(s.pred_act)
                apreds = get(s.pred_act[ai], v, nothing)
                if apreds !== nothing
                    for k in apreds
                        invalidate_action!(s, ai, k)
                    end
                end
            end
        end
    end
end


# Same reachability as ReACTed_reachable_cell_online, but the shield is computed with the event-driven predecessor
# relation above instead of repeated shield_step! sweeps. Returns the same tuple so it can be swapped in ShieldSandbox.jl.
function ReACTed_reachable_cell_online_pred(system::EuclideanHybridSystem, p, granularity, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, clustering=true, max_steps::Int64=100000; act_needs_safe_flow::Bool=true, include_border_cells::Bool=false)
    return _ReACTed_reachable_cell_online_pred(Grid(system.statespace, granularity), system, p, granularity, δ⁻, δ⁺, alg, maxOrder, reduceOrder, clustering, max_steps, act_needs_safe_flow, include_border_cells)
end

# Function barrier: the body is compiled for the concrete Grid{N} and system type, so N can be used in type positions.
function _ReACTed_reachable_cell_online_pred(grid::Grid{N}, system::EuclideanHybridSystem, p, granularity, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg, maxOrder::Int, reduceOrder::Int, clustering, max_steps::Int64, act_needs_safe_flow::Bool, include_border_cells::Bool) where N
    @show foldl(*, grid.numCells), grid.numCells
    sweep_factor = 20
    TIME_FUNCS = true

    recorderTimeStart = time_ns()

    phiDict, tPhiDict, inputDict = phi_input_dict(system.flowMatrix, system.input, δ⁻, δ⁺, alg, maxOrder, reduceOrder)
    d = δ⁻
    isInvA = isinvertible(system.flowMatrix)
    A = system.flowMatrix
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(system.flowMatrix)
    Φcache = system.flowMatrix == A_abs ? phiDict[d] : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)

    reachable_by_action = ActionRelation{N}()
    reachable_by_flow = FlowRelation{N}()
    degenerate_dimensions = zeros(N) # Grid centres zero-width dimensions on their value, so cells need no extra offset

    invariant = hyperrectangle_to_HPolyhedron(system.statespace)
    U = system.input
    max_input = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, p, alg, false, nothing), U)
    max_flow = exp(p .* system.flowMatrix)

    guardProjVectors, guardProjBounds = getHalfSpaceProjections(system.edges[1].guard)
    invariantProjVectors, invariantProjBounds = getHalfSpaceProjections(invariant)
    guard_invariant_intersection = intersection(system.edges[1].guard, invariant)

    # Coarse sweep marking blocks that cannot hit a guard/invariant within the period.
    fast_cells = falses(cellshape(grid))
    sweep = [min(sweep_factor, n) for n in grid.numCells] # cells per sweep block in each dimension (dimensions smaller than a block use one block)
    Z_big = remove_zero_generators(Zonotope(degenerate_dimensions, diagm(map((n, s) -> n == 1 ? 0. : granularity * (s/2), grid.numCells, sweep))))
    lower_offset_big = grid.lower .+ granularity .* (sweep ./ 2)
    generatorDict_big = ReACT_discretize_decomposed_generators(Z_big, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict)
    discretizationDict_big, finaldims_big = ReACT_discretize_combine_with_offsets(Z_big, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict_big)
    preallocated_center_big = zeros(grid.dimension)
    preallocated_genmat_big = zeros(grid.dimension, finaldims_big)
    preallocated_inter_genmat_big = zeros(grid.dimension, finaldims_big)
    ranges = ntuple(i -> 1:fld(grid.numCells[i], sweep[i]), N)

    attempt_recorder = zeros(4)
    identity_matrix = diagm(ones(grid.dimension))
    for bidx in CartesianIndices(ranges)
        of = lower_offset_big + ((car2vec(bidx) .- 1) .* granularity .* sweep)
        LazySets.API.translate!(Z_big, of)
        ReACT_discretize_combine_with_offset_vector!(preallocated_center_big, preallocated_genmat_big, preallocated_inter_genmat_big, discretizationDict_big, of, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict_big)
        f = propagate_set_check_hybrid_jump(Z_big, [0.0, p], δ⁻, δ⁺, system, system.edges, identity_matrix, attempt_recorder, grid.dimension, P2A_abs, guardProjVectors, guardProjBounds, invariantProjVectors, invariantProjBounds, discretizationDict_big, phiDict, tPhiDict, inputDict, max_input, max_flow, guard_invariant_intersection, alg, maxOrder, reduceOrder)
        if f == 0
            fast_cells[ntuple(i->(1+(bidx[i]-1)*sweep[i]):(bidx[i]*sweep[i]), N)...] .= true
        end
        LazySets.API.translate!(Z_big, -of)
    end

    Z = remove_zero_generators(Zonotope(degenerate_dimensions, diagm(map(x -> x == 1 ? 0. : granularity/2, grid.numCells))))
    count = 0
    lower_offset = grid.lower .+ (granularity/2)
    generatorDict = ReACT_discretize_decomposed_generators(Z, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict)
    discretizationDict, finaldims = ReACT_discretize_combine_with_offsets(Z, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict)
    preallocated_center = zeros(grid.dimension)
    preallocated_genmat = zeros(grid.dimension, finaldims)
    preallocated_inter_genmat = zeros(grid.dimension, finaldims)

    max_zonotope::Zonotope{Float64,Vector{Float64},Matrix{Float64}} = Zonotope(max_input.center, hcat(genmat(max_input), genmat(linear_map(max_flow, Z)), diagm(fill(granularity/2, grid.dimension))))
    #max_zonotope = minkowski_sum(max_zonotope, Zonotope(zeros(grid.dimension), diagm(fill(granularity/2, grid.dimension))))
    # Fast cells all share max_zonotope's generators, so their cover test uses its precomputed H-representation.
    fast_cover = ConstantZonotopeCover(max_zonotope, Val(N))
    fast_center = zeros(N)

    dims = grid.dimension
    offset_list = offsets(grid)
    n_acts = length(system.Act)

    # BitArray per action instead of `idx in may_act_dict[act.id]` (linear search per cell).
    may_act = [falses(cellshape(grid)) for _ in 1:n_acts]
    for (ai, act) in enumerate(system.Act)
        for gidx in get_touching_cell_idxs_l(grid, act.guard)
            may_act[ai][gidx] = true
        end
    end

    s = OnlinePredShieldState(grid.deadCells, n_acts, act_needs_safe_flow)

    frontier_list = CartesianIndex{N}[]
    for edge in system.edges
        append!(frontier_list, get_touching_cell_idxs_l(grid, edge.guard))
    end

    initial_dead = CartesianIndex{N}[]
    for global_constraint in system.globalConstraints
        for didx in get_touching_cell_idxs_l(grid, global_constraint)
            s.dead[didx] = true
            s.flow_bad[didx] = true
            s.computed[didx] = true
            push!(initial_dead, didx)
        end
    end
    append!(frontier_list, grow_indices_t(grid, initial_dead, offset_list))
    if include_border_cells
        border = border_cells(grid)
        append!(frontier_list, border)
        append!(frontier_list, grow_indices_t(grid, border, offset_list))
    end
    unique!(frontier_list)

    setup_time = time_ns()
    aux_time_recoder = time_ns()
    reachability_time = 0.0
    shield_time = 0.0
    shield_count::Int64 = 0

    act_succs = Vector{Union{Nothing,CellIdxs{N}}}(nothing, n_acts)
    shield_count_cover = zeros(Float64, cellshape(grid))

    while !isempty(frontier_list) && shield_count < max_steps
        if TIME_FUNCS
            aux_time_recoder = time_ns()
        end
        for idx in frontier_list
            if s.computed[idx]
                continue
            end

            flow_succ = nothing
            count += 1
            if !fast_cells[idx]
                of = lower_offset + ((car2vec(idx) .- 1) .* granularity)
                LazySets.API.translate!(Z, of)
                ReACT_discretize_combine_with_offset_vector!(preallocated_center, preallocated_genmat, preallocated_inter_genmat, discretizationDict, of, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict)
                z, f = propagate_set_t(Z, [0.0, p], δ⁻, δ⁺, system, grid, fast_cells, system.edges, identity_matrix, attempt_recorder, grid.dimension, P2A_abs, guardProjVectors, guardProjBounds, invariantProjVectors, invariantProjBounds, discretizationDict, phiDict, tPhiDict, inputDict, max_input, max_flow, guard_invariant_intersection, alg, maxOrder, reduceOrder)
                if f != 1
                    flow_succ = get_touching_cell_idxs_lz(grid, z)
                    reachable_by_flow[idx] = flow_succ
                end
                LazySets.API.translate!(Z, -of)
            else
                fast_cell_center!(fast_center, max_flow, max_zonotope.center, grid, idx)
                flow_succ = get_touching_cell_idxs_cover(grid, fast_cover, fast_center)
                reachable_by_flow[idx] = flow_succ
            end

            if count % 1000 == 0
                println(count)
            end

            fill!(act_succs, nothing)
            for (ai, act) in enumerate(system.Act)
                if may_act[ai][idx]
                    of = lower_offset + ((car2vec(idx) .- 1) .* granularity)
                    LazySets.API.translate!(Z, of)
                    # This should propably check if the cells are fast_cells and use a precomputed max_ai_zonotope in the case that it is and actually propagate the set otherwise.
                    succ = get_touching_cell_idxs_lz(grid, minkowski_sum(linear_map(max_flow, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector)), max_input))
                    reachable_by_action[(act.id, idx)] = succ
                    act_succs[ai] = succ
                    LazySets.API.translate!(Z, -of)
                end
            end

            register_cell!(s, idx, flow_succ, act_succs)
        end
        if TIME_FUNCS
            reachability_time += (time_ns() - aux_time_recoder)
            aux_time_recoder = time_ns()
        end

        propagate_bad!(s)

        for jdx in s.changed
            if s.dead[jdx] && haskey(reachable_by_flow, jdx) && shield_count_cover[jdx] == 0.0
                shield_count_cover[jdx] = length(reachable_by_flow[jdx]) / max(1, count_ones_at(s.dead, reachable_by_flow[jdx]))
            end
        end

        shield_count += 1
        if shield_count % 5 == 0
            println("$shield_count : $(Dates.format(now(), "HH:MM:SS"))")
        end

        frontier_list = isempty(s.changed) ? CartesianIndex{N}[] : grow_indices_t(grid, unique!(s.changed), offset_list)
        filter!(t -> !s.computed[t], frontier_list)
        empty!(s.changed)
        if TIME_FUNCS
            shield_time += (time_ns() - aux_time_recoder)
        end
    end

    if isempty(frontier_list)
        println("Fixed point found after $shield_count rounds!")
    else
        println("Stopped at max_steps = $max_steps before reaching a fixed point.")
    end

    if TIME_FUNCS
        println("Total times:
        Total: $((time_ns() - recorderTimeStart)/10^9)
        Setup: $((setup_time - recorderTimeStart) / 10^9)
        Reachability: $(reachability_time / 10^9)
        Shield: $(shield_time / 10^9)")
        println("Shield rounds: $shield_count, cells computed: $count")
    end

    return grid, reachable_by_action, reachable_by_flow, s.flow_bad, shield_count_cover
end

@inline function count_ones_at(m::BitArray, idxs)
    c = 0
    for i in idxs
        c += m[i]
    end
    return c
end
