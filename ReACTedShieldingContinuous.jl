using LazySets, LinearAlgebra, ReachabilityAnalysis
# Expects ReACTedShieldingModal.jl to be included first.

#=
Online shield where the controller is not tied to a decision period: it may take an action at any time.
The interval T passed as the period is only how often the state is checked. In each interval a cell either does not act, as in ReACTed_reachable_cell_modal, or acts at some time in [0, T). Acting is only allowed if it is safe at every one of those times.
A controller that may act at any real time then stays safe, as long as the action was allowed at the most recent check.
T cannot simply be made small. If a period moves an axis much less than a cell, the rounding to cells can hide the movement, so a smaller T needs a smaller granularity.
=#

const NO_GUARD = (Vector{Float64}[], Float64[])

# Like modal_reach!, for X that takes action `id` at some time before the end of the period. Each piece of δ⁻ of its flowpipe takes the action on its own and continues in the action's target. Pieces touching a guard also continue through the edge, with the action still to come.
# Returns false if the action is bad at any of those times, or a location on the way does not offer it.
function anytime_action_reach!(out::Vector{ModalEndSet}, ctx::ModalContext, li::Int, id::Int, X, D::Int, W::Int, exit_bad::Bool, lo::Vector{Float64}, hi::Vector{Float64}, depth::Int=0)
    depth > ctx.max_depth && return false
    actions = ctx.system.locations[li].actions
    position = findfirst(a -> a.action.id == id, actions)
    position === nothing && return false
    a = actions[position]
    loc = ctx.locations[li]
    target = ctx.locations[a.target]
    total = D + W
    if total == 0
        jumped = apply_jump(X, a.action.jumpMatrix, a.action.jumpVector, ctx.maxOrder, ctx.reduceOrder)
        clocks = clock_after_jump(lo, hi, NO_GUARD, a.action.jumpMatrix, a.action.jumpVector, loc.clock_rate, target.clock_rate, 0, 0, ctx.δ⁻)
        return modal_reach!(out, ctx, a.target, jumped, 0, 0, exit_bad, clocks..., depth + 1)
    end
    dyn = loc.dynamics
    Ω = newReACTDiscretizePlus(X, ctx.δ⁻, ctx.δ⁺, dyn.A, dyn.P2A_abs, dyn.phiDict, dyn.inputDict, ctx.alg, ctx.maxOrder, ctx.reduceOrder)[ctx.δ⁻]
    for k in 0:total-1
        piece = minkowski_sum(linear_map(dyn.Φ[k+1], Ω), dyn.V[k+1])
        any(set -> someinside(piece, set...), ctx.unsafe) && return false
        # A state in this piece has spent between k and k + 1 steps here, and had between D and D + W left at the start.
        D′ = max(D - k - 1, 0)
        W′ = total - k - D′
        jumped = apply_jump(piece, a.action.jumpMatrix, a.action.jumpVector, ctx.maxOrder, ctx.reduceOrder)
        clocks = clock_after_jump(lo, hi, NO_GUARD, a.action.jumpMatrix, a.action.jumpVector, loc.clock_rate, target.clock_rate, D′, W′, ctx.δ⁻)
        modal_reach!(out, ctx, a.target, jumped, D′, W′, exit_bad, clocks..., depth + 1) || return false
        all_jumped = false
        for (e, guard) in zip(ctx.system.locations[li].edges, loc.guards)
            someinside(piece, guard...) || continue
            clipped = clip_to_guard(piece, guard)
            clipped === nothing && continue
            edge_clocks = clock_after_jump(lo, hi, guard, e.jumpMatrix, e.jumpVector, loc.clock_rate, ctx.locations[e.targetLoc].clock_rate, D′, W′, ctx.δ⁻)
            edge_clocks === nothing && continue
            edge_jumped = apply_jump(clipped, e.jumpMatrix, e.jumpVector, ctx.maxOrder, ctx.reduceOrder)
            anytime_action_reach!(out, ctx, e.targetLoc, id, edge_jumped, D′, W′, exit_bad, edge_clocks..., depth + 1) || return false
            all_jumped = all_jumped || allunder(piece, guard...)
        end
        all_jumped && return true
    end
    return true
end

# Successors of taking the action at position ai of ctx.action_ids at some time in the interval, or nothing if it is not available or is bad at any of those times.
function action_successors_anytime(ctx::ModalContext{D,N}, layer::Int, c::CartesianIndex{D}, ai::Int, exit_bad::Bool) where {D,N}
    res = CartesianIndex{N}[]
    out = ModalEndSet[]
    Z = Zonotope(ctx.grid.lower .+ (car2vec(c) .- 0.5) .* ctx.grid.granularity, ctx.cell_generators)
    for li in cell_locations(ctx, layer, c)
        ctx.may_act[li][ai][c] || return nothing
        empty!(out)
        lo, hi = clock_end(ctx.locations[li].clock_rate, Z, ctx.steps, ctx.δ⁻)
        anytime_action_reach!(out, ctx, li, ctx.action_ids[ai], Z, ctx.steps, 0, exit_bad, lo, hi) || return nothing
        append_touching!(res, ctx, out)
    end
    return unique!(res)
end

# The online shield with actions at any time. T is the interval between checks. Takes the same arguments and returns the same values as ReACTed_reachable_cell_modal.
ReACTed_reachable_cell_continuous(system::ModalHybridSystem, T, granularity, δ⁻::Float64, δ⁺::Float64, args...; kwargs...) =
    ReACTed_reachable_cell_modal(system, T, granularity, δ⁻, δ⁺, args...; online=true, action_reach=action_successors_anytime, anytime_predecessors=true, kwargs...)
