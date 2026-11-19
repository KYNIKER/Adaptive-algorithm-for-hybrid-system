# Shield synthesis for an EuclideanHybridSystem. Online (the default) it computes cells from a frontier that starts at
# the guards, the constraint cells and (with include_border_cells) the border, and grows around cells whose status
# changed; cells never reached are assumed safe. With online = false it computes every cell (see ReACTed_reachable_cell_t).
# To keep it fast:
#   - Action successors use one precomputed cover per action instead of a zonotope and a cover per cell.
#   - Action guard masks come from dot products against the guard's halfspaces instead of an LP sweep.
#   - Successor lists and the shield's predecessor lists are kept in flat arrays (flat.jl), with the
#     cell's linear index as its id, instead of Dicts of vectors.
# It returns grid, allowed, lost, info like every shield (see choice_bits in flat.jl), with info.act_bad
# (actions of uncomputed cells whose jump lands on a dead cell) and info.shield_count_cover (the round each cell died in, negative when all its flow successors were dead).

function ReACTed_reachable_cell_online_pred_fast(system::EuclideanHybridSystem, p, granularity, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, max_steps::Int64=100000; act_needs_safe_flow::Bool=true, include_border_cells::Bool=false, sweep_factor::Int=20, jump_hull::Bool=false, online::Bool=true, verbose::Bool=true)
    return _ReACTed_reachable_cell_online_pred_fast(Grid(system.statespace, granularity), system, p, δ⁻, δ⁺, alg, maxOrder, reduceOrder, max_steps, act_needs_safe_flow, include_border_cells, sweep_factor; jump_hull, online, verbose)
end

# Function barrier on Grid{N}
function _ReACTed_reachable_cell_online_pred_fast(grid::Grid{N}, system::EuclideanHybridSystem, p, δ⁻::Float64, δ⁺::Float64, alg, maxOrder::Int, reduceOrder::Int, max_steps::Int64, act_needs_safe_flow::Bool, include_border_cells::Bool, sweep_factor::Int; jump_hull::Bool=false, online::Bool=true, verbose::Bool=true) where N
    timer = PhaseTimer()
    reach = reach_setup(system, p, δ⁻, δ⁺, alg, maxOrder, reduceOrder, Val(N); jump_hull)
    covers = CoverCache{N}()

    fast_cells = LazyFastCells(system, reach, grid, sweep_factor) # decided per block of sweep_factor cells per dimension
    cell_reach = CellReach(reach, grid)
    fast_covers, act_covers = reach_covers(system, reach, cell_reach.Z, [grid.granularity / 2], Val(N))
    act_covs = [covs[1] for covs in act_covers]

    n_acts = length(system.Act)
    shape = Tuple(cellshape(grid))
    n = prod(shape)
    n < typemax(Int32) || throw(ArgumentError("$n cells do not fit in Int32 ids"))
    ids = reshape(Int32.(1:n), shape) # cell -> id
    cells = CartesianIndices(shape) # id -> cell
    # An action does nothing outside its guard (see choice_bits), and only some states take it in a cell partly inside.
    masks = [guard_masks(grid, act.guard) for act in system.Act]
    may_act = [m[1] for m in masks]
    must_act = [m[2] for m in masks]
    offset_list = offsets(grid)

    succ = SuccessorStore(n, [BitVector(vec(m)) for m in may_act])
    s = reset!(FlatShieldState(n_acts, act_needs_safe_flow), n)
    buf = Int32[]
    cover!(out, cov, center) = cover_ids!(out, ids, grid, cov, center)
    touching!(out, z) = foreach(cc -> push!(out, ids[cc]), get_touching_cell_idxs_lz(grid, z, covers))

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
    online || (frontier = vec(collect(cells)))

    lap!(timer, "Setup")
    shield_count = 0
    cells_computed = 0
    shield_count_cover = zeros(Float64, shape)

    while !isempty(frontier) && shield_count < max_steps
        for c in frontier
            k = ids[c]
            s.computed[k] && continue
            cells_computed += 1
            compute_successors!(succ, buf, system, reach, cell_reach, grid, c, k, fast_cells[c], fast_covers[1], act_covs, must_act, cover!, touching!)
            register_successors!(succ, s, k)
        end
        lap!(timer, "Reachability")

        propagate_bad!(s)
        changed = unique!(s.changed)
        for k in changed
            # Successors of a dead cell that are dead too: 1.0 when all of them are, 0.5 otherwise
            if s.dead[k] && !succ.flow_none[k] && shield_count_cover[k] == 0.0
                range = successor_range(succ, 0, k)
                n_dead = count(i -> s.dead[succ.pool[i]], range)
                shield_count_cover[k] = length(range) == max(1, n_dead) ? -shield_count : shield_count
            end
        end
        shield_count += 1
        frontier = (!online || isempty(changed)) ? CartesianIndex{N}[] : grow_indices_t(grid, [cells[k] for k in changed], offset_list)
        filter!(c -> !s.computed[ids[c]], frontier)
        empty!(s.changed)
        lap!(timer, "Shield")
    end

    allowed, lost = cell_choices(shape, n_acts, s, store_enabled(succ), store_acted(succ), store_noop(succ))
    grid.deadCells .= allowed .== 0
    verbose && report(timer, isempty(frontier) ? "Fixed point after $shield_count rounds, $cells_computed cells computed" :
                             "Stopped at max_steps = $max_steps before a fixed point, $cells_computed cells computed")
    return grid, allowed, lost, (; shield_count_cover)
end

# The shield over every cell of the grid: ReACTed_reachable_cell_online_pred_fast with online = false. clustering is not used.
ReACTed_reachable_cell_t(system::EuclideanHybridSystem, p, granularity, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, clustering=true; jump_hull::Bool=false, sweep_factor::Int=25, verbose::Bool=true) =
    ReACTed_reachable_cell_online_pred_fast(system, p, granularity, δ⁻, δ⁺, alg, maxOrder, reduceOrder; jump_hull, sweep_factor, online=false, verbose)
