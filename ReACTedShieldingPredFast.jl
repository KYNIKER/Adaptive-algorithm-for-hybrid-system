using LazySets, LinearAlgebra, ReachabilityAnalysis, Dates
# Expects PartitionUtilities.jl, ReACTedShielding.jl and ReACTedShieldingPred.jl to be included first.
#
# ReACTed_reachable_cell_online_pred with the same exploration and shield, but faster:
#   - Action successors use one precomputed cover per action instead of a zonotope and a cover per cell.
#   - Action guard masks come from dot products against the guard's halfspaces instead of the LP sweep.
#   - Successor lists and the shield's predecessor lists are kept in flat arrays (ReACTedShieldingFlat.jl), with the
#     cell's linear index as its id, instead of Dicts of vectors.
#   - It returns grid, allowed, lost, info like every shield (see choice_bits in ReACTedShieldingFlat.jl), with info.act_bad
#     and info.shield_count_cover as in ReACTed_reachable_cell_online_pred.

function ReACTed_reachable_cell_online_pred_fast(system::EuclideanHybridSystem, p, granularity, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, max_steps::Int64=100000; act_needs_safe_flow::Bool=true, include_border_cells::Bool=false, sweep_factor::Int=20, jump_hull::Bool=false)
    return _ReACTed_reachable_cell_online_pred_fast(Grid(system.statespace, granularity), system, p, δ⁻, δ⁺, alg, maxOrder, reduceOrder, max_steps, act_needs_safe_flow, include_border_cells, sweep_factor; jump_hull)
end

# Function barrier on Grid{N}
function _ReACTed_reachable_cell_online_pred_fast(grid::Grid{N}, system::EuclideanHybridSystem, p, δ⁻::Float64, δ⁺::Float64, alg, maxOrder::Int, reduceOrder::Int, max_steps::Int64, act_needs_safe_flow::Bool, include_border_cells::Bool, sweep_factor::Int; jump_hull::Bool=false) where N
    recorderTimeStart = time_ns()
    reach = reach_setup(system, p, δ⁻, δ⁺, alg, maxOrder, reduceOrder, Val(N); jump_hull)
    covers = CoverCache{N}()

    # Fast cells, decided per block of sweep cells per dimension
    sweep = [min(sweep_factor, n) for n in grid.numCells]
    block_reach = CellReach(reach, grid, sweep)
    fast_cells = falses(cellshape(grid))
    for b in CartesianIndices(ntuple(i -> 1:fld(grid.numCells[i], sweep[i]), N))
        is_fast_block(system, reach, block_reach, b) && (fast_cells[ntuple(i -> (1+(b[i]-1)*sweep[i]):(b[i]*sweep[i]), N)...] .= true)
    end
    cell_reach = CellReach(reach, grid)
    fast_covers, act_covers = reach_covers(system, reach, cell_reach.Z, [grid.granularity / 2], Val(N))

    n_acts = length(system.Act)
    shape = Tuple(cellshape(grid))
    n = prod(shape)
    n < typemax(Int32) || throw(ArgumentError("$n cells do not fit in Int32 ids"))
    ids = reshape(Int32.(1:n), shape) # cell -> id
    cells = CartesianIndices(shape) # id -> cell
    may_act = [touching_masks([grid], act.guard, ntuple(_ -> 1, N))[1] for act in system.Act]
    offset_list = offsets(grid)

    succ = SuccessorStore(n, [BitVector(vec(m)) for m in may_act])
    s = reset!(FlatShieldState(n_acts, act_needs_safe_flow), n)
    buf = Int32[]

    # Same initial frontier and seeding order as ReACTed_reachable_cell_online_pred
    frontier = CartesianIndex{N}[]
    for edge in system.edges
        append!(frontier, get_touching_cell_idxs_l(grid, edge.guard))
    end
    initial_dead = CartesianIndex{N}[]
    for global_constraint in system.globalConstraints
        for c in get_touching_cell_idxs_l(grid, global_constraint)
            seed_constraint!(s, ids[c])
            succ.flow_none[ids[c]] = true # never computed, so no flow entry
            push!(initial_dead, c)
        end
    end
    append!(frontier, grow_indices_t(grid, initial_dead, offset_list))
    if include_border_cells
        border = border_cells(grid)
        append!(frontier, border)
        append!(frontier, grow_indices_t(grid, border, offset_list))
    end
    unique!(frontier)

    setup_time = time_ns()
    reachability_time = 0.0
    shield_time = 0.0
    shield_count = 0
    cells_computed = 0
    shield_count_cover = zeros(Float64, shape)

    while !isempty(frontier) && shield_count < max_steps
        aux_time_recoder = time_ns()
        for c in frontier
            k = ids[c]
            s.computed[k] && continue
            cells_computed += 1
            empty!(buf)
            if fast_cells[c]
                fast_cell_center!(reach.center, reach.max_flow, reach.max_input.center, grid, c)
                cover_ids!(buf, ids, grid, fast_covers[1], reach.center)
                succ.flow_none[k] = false
            else
                z, f = cell_flow(reach, cell_reach, c)
                succ.flow_none[k] = f == 1
                if f != 1
                    for cc in get_touching_cell_idxs_lz(grid, z, covers)
                        push!(buf, ids[cc])
                    end
                end
            end
            set_successors!(succ, 0, k, buf)
            for (ai, act) in enumerate(system.Act)
                empty!(buf)
                if succ.act_enabled[ai][k]
                    act_center!(reach.center, reach, cell_reach, act, c)
                    cover_ids!(buf, ids, grid, act_covers[ai][1], reach.center)
                end
                set_successors!(succ, ai, k, buf)
            end
            register_successors!(succ, s, k)
        end
        reachability_time += time_ns() - aux_time_recoder
        aux_time_recoder = time_ns()

        propagate_bad!(s)
        changed = unique!(s.changed)
        for k in changed
            # Successors of a dead cell that are dead too: 1.0 when all of them are, 0.5 otherwise
            if s.dead[k] && !succ.flow_none[k] && shield_count_cover[k] == 0.0
                range = successor_range(succ, 0, k)
                n_dead = count(i -> s.dead[succ.pool[i]], range)
                shield_count_cover[k] = length(range) == max(1, n_dead) ? 1.0 : 0.5
            end
        end
        shield_count += 1
        frontier = isempty(changed) ? CartesianIndex{N}[] : grow_indices_t(grid, [cells[k] for k in changed], offset_list)
        filter!(c -> !s.computed[ids[c]], frontier)
        empty!(s.changed)
        shield_time += time_ns() - aux_time_recoder
    end
    println(isempty(frontier) ? "Fixed point found after $shield_count rounds!" : "Stopped at max_steps = $max_steps before reaching a fixed point.")

    # Actions of cells never computed whose jump (without the flow after it) lands on a computed dead cell
    act_bad = falses(shape)
    for (ai, act) in enumerate(system.Act)
        ActZonotope = remove_zero_generators(Zonotope(copy(act.jumpVector), act.jumpMatrix*diagm(fill(grid.granularity/2, N))))
        for c in findall(may_act[ai])
            s.computed[ids[c]] && continue
            of = act.jumpMatrix * cell_offset(cell_reach, c)
            LazySets.API.translate!(ActZonotope, of)
            any(cc -> s.computed[ids[cc]] && s.dead[ids[cc]], get_touching_cell_idxs_lz(grid, ActZonotope, covers)) && (act_bad[c] = true)
            LazySets.API.translate!(ActZonotope, -of)
        end
    end

    allowed, lost = cell_choices(shape, n_acts, s, store_enabled(succ), store_acted(succ))
    grid.deadCells .= allowed .== 0
    println("Total times:
        Total: $((time_ns() - recorderTimeStart)/10^9)
        Setup: $((setup_time - recorderTimeStart) / 10^9)
        Reachability: $(reachability_time / 10^9)
        Shield: $(shield_time / 10^9)")
    println("Shield rounds: $shield_count, cells computed: $cells_computed")
    return grid, allowed, lost, (; act_bad, shield_count_cover)
end

