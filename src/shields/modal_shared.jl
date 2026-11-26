#=
What the modal shields (modal.jl, modal_adaptive.jl, modal_simple.jl) share: the system types, and the per-axis bounds of
the states at the end of a period.

Bounds: a set reaches the cells it touches (touching_cells), within per-axis bounds that hold for its states at the end of the period (end_bounds, jumped_bounds). Cells are half-open, [l, l + g) per axis, so where these bounds are exact a cell that only touches their open end is dropped: otherwise a set that moves by whole cells, such as a clock, would also reach a neighbour every period. They are also what the exit check uses, as the sets after a jump are bloated.
  - Flow and actions: the end state of x₀ in the box is T x₀ + c, exactly, on the axes where the input is constant. An open upper end is moved down by open_margin, so all bounds are closed.
  - A guard jumps at a time that depends on the state. An axis keeps its bounds when its end value does not depend on that time, and loses them otherwise: [-Inf, Inf].
=#

# Several ModalActions with the same action id in one location, and the same guard, are one choice with several outcomes, which the shield all takes as successors.
struct ModalAction
    action::Action
    target::Int
end

# A location with starts = false is only entered by actions and edges. So a layer can be a decision point that every period starts from, 
# with the outcome of each choice in locations of their own.
struct ModalLocation
    layer::Int
    region::Union{HPolyhedron,Nothing} # part of the continuous state this location covers, nothing = all of it
    A::Matrix{Float64}
    U::Zonotope{Float64,Vector{Float64},Matrix{Float64}}
    edges::Vector{Edge}
    actions::Vector{ModalAction}
    starts::Bool
end
ModalLocation(layer, region, A, U, edges, actions; starts::Bool=true) = ModalLocation(layer, region, A, U, edges, actions, starts)

struct ModalHybridSystem{S<:LazySet,C<:LazySet}
    statespace::S
    unsafe::Vector{C}
    locations::Vector{ModalLocation}
end

# A system without discrete modes, as one location. Its edges must lead back to location 1.
function ModalHybridSystem(system::EuclideanHybridSystem)
    U = Zonotope(Vector{Float64}(LazySets.center(system.input)), Matrix{Float64}(genmat(system.input)))
    location = ModalLocation(1, nothing, Matrix{Float64}(system.flowMatrix), U, system.edges, [ModalAction(act, 1) for act in system.Act])
    return ModalHybridSystem(system.statespace, system.globalConstraints, [location])
end

const MODAL_TOL = 1e-9

const Halfspaces = Tuple{Vector{Vector{Float64}},Vector{Float64}}
const Affine = Tuple{Matrix{Float64},Vector{Float64}} # x = T x₀ + c, for x₀ in the start box
const Bounds = Tuple{Vector{Float64},Vector{Float64}} # per axis, every state is in [lo, hi] at the end of the period

# How far an open end is moved down: more than the MODAL_TOL with which bounds are compared, in cells and absolutely.
open_margin(g::Float64) = 4MODAL_TOL * max(1.0, g)

# The bounds of the box Z mapped by x = T x₀ + c. The upper end of an axis is open when some x₀ⱼ with Tᵢⱼ > 0 is at the
# open end of the box, which an axis of one cell has not.
function affine_bounds(grid::Grid, T::Matrix{Float64}, c::Vector{Float64}, Z)
    l, u = low(Z)', high(Z)'
    lo = c + vec(sum(min.(T .* l, T .* u); dims=2))
    hi = c + vec(sum(max.(T .* l, T .* u); dims=2))
    for i in axes(T, 1)
        any(j -> T[i, j] > 0 && grid.numCells[j] > 1, axes(T, 2)) && (hi[i] -= open_margin(grid.granularity))
    end
    return lo, hi
end

# The bounds of the box Z mapped by x = T x₀ + c and flowed for a whole period. An input that varies on an axis leaves it
# untracked: [-Inf, Inf].
function end_bounds(ctx, dyn::Dynamics, (T, c)::Affine, Z)
    Φ, V = dyn.Φ[ctx.steps+1], dyn.V[ctx.steps+1]
    lo, hi = affine_bounds(ctx.grid, Φ * T, Φ * c + LazySets.center(V), Z)
    for i in eachindex(lo)
        iszero(view(genmat(V), i, :)) || ((lo[i], hi[i]) = (-Inf, Inf))
    end
    return lo, hi
end

# The bounds after a jump x := M x + j through the guard, from a location with rates `from` to one with rates `to`, for
# states with between D′ and D′ + W′ steps left. The time of the jump depends on the state, so only these axes keep
# bounds:
#   - translated at the same rate on both sides: it ends where it would have without the jump, plus j. A halfspace of the
#     guard on this axis alone also limits it, as the states were inside the guard and then moved for the time left.
#   - reset to j: it ends at j plus its rate times the time left.
# The most time left, D′ + W′, needs a state that jumps at the very start of the piece. In a later piece that state also
# jumps at the end of the piece before, which covers it, so with a positive rate the end it reaches is open. Only in the
# first piece, when some state is inside the guard at the start (at_start), is it closed.
# Returns nothing if no state of the bounds can take the jump.
function jumped_bounds((lo, hi)::Bounds, guard::Halfspaces, M::Matrix{Float64}, j::Vector{Float64}, from::Vector{Float64}, to::Vector{Float64}, D′::Int, W′::Int, δ⁻::Float64, at_start::Bool, g::Float64)
    lo′, hi′ = fill(-Inf, length(j)), fill(Inf, length(j))
    for i in eachindex(j)
        isnan(to[i]) && continue
        least, most = minmax(to[i] * D′ * δ⁻, to[i] * (D′ + W′) * δ⁻)
        to[i] > 0 && !at_start && (most -= open_margin(g))
        if iszero(view(M, i, :))
            lo′[i], hi′[i] = j[i] + least, j[i] + most
        elseif from[i] == to[i] && M[i, i] == 1 && count(!iszero, view(M, i, :)) == 1
            l, h = lo[i], hi[i]
            for (a, c) in zip(guard...)
                (a[i] != 0 && count(!iszero, a) == 1) || continue
                a[i] > 0 ? (h = min(h, c / a[i] + most)) : (l = max(l, c / a[i] + least))
            end
            l > h + MODAL_TOL && return nothing
            lo′[i], hi′[i] = l + j[i], h + j[i]
        end
    end
    return lo′, hi′
end

# Per axis, the cells the bounds allow, as keep_lo, keep_hi of touching_cells: none for bounds below the grid.
function cell_bounds(grid::Grid{D}, (lo, hi)::Bounds) where D
    index(x, i) = floor(Int, clamp((x - grid.lower[i]) / grid.granularity + MODAL_TOL, -1, grid.numCells[i])) + 1
    keep_lo = ntuple(i -> index(lo[i], i), D)
    return keep_lo, ntuple(i -> max(index(hi[i], i), keep_lo[i]), D)
end

# Whether the box of the bounds meets the guard. A halfspace the box lies beyond decides it exactly; otherwise isdisjoint
# does, which counts touching as meeting.
function box_meets((lo, hi)::Bounds, guard::Halfspaces)
    any(((a, b),) -> sum(a[i] * (a[i] > 0 ? lo[i] : hi[i]) for i in eachindex(a)) > b + MODAL_TOL, zip(guard...)) && return false
    return !isdisjoint(Hyperrectangle(low=lo, high=max.(lo, hi)), HPolyhedron([LazySets.HalfSpace(a, b) for (a, b) in zip(guard...)]))
end

# True if the interior of the box [lo, hi] overlaps the region.
function box_overlaps_region(region, lo, hi)
    region === nothing && return true
    for (a, b) in zip(region...)
        low_value = sum(a[i] * (a[i] > 0 ? lo[i] : hi[i]) for i in eachindex(a))
        low_value >= b - MODAL_TOL && return false
    end
    return true
end

# Discretization over a whole period of a box of shape G under the location's dynamics, from which modal_reach! gets the
# discretization of any box of that shape (see CellDiscretization).
period_discretization(dyn::Dynamics, U, G::Matrix{Float64}, δ⁻::Float64, δ⁺::Float64) =
    CellDiscretization(Zonotope(zeros(size(G, 1)), G), (; δ⁻, δ⁺, dyn.A, dyn.P2A_abs, phiDict=dyn.step_Φ, inputDict=dyn.step_input, U), Val(size(G, 1)))
