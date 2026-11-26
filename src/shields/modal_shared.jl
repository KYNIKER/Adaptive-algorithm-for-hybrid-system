# Shared by the modal shields: system types and per-axis Bounds of the states at the end of a period.
# Cells are half-open [l, l + g); bounds are closed, an open end is moved down by open_margin.

# ModalActions of one location with the same id and guard are one choice with several outcomes.
struct ModalAction
    action::Action
    target::Int
end

# A location with starts = false is only entered by actions and edges.
struct ModalLocation
    layer::Int
    region::Union{HPolyhedron,Nothing} # nothing: the whole state space
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

# A system without modes as one location; its edges must lead back to location 1.
function ModalHybridSystem(system::EuclideanHybridSystem)
    U = Zonotope(Vector{Float64}(LazySets.center(system.input)), Matrix{Float64}(genmat(system.input)))
    location = ModalLocation(1, nothing, Matrix{Float64}(system.flowMatrix), U, system.edges, [ModalAction(act, 1) for act in system.Act])
    return ModalHybridSystem(system.statespace, system.globalConstraints, [location])
end

const MODAL_TOL = 1e-9
const Halfspaces = Tuple{Vector{Vector{Float64}},Vector{Float64}}
const Affine = Tuple{Matrix{Float64},Vector{Float64}} # x = T x₀ + c
const Bounds = Tuple{Vector{Float64},Vector{Float64}} # every state is in [lo, hi]

open_margin(g::Float64) = 4MODAL_TOL * max(1.0, g) # exceeds the MODAL_TOL of comparisons

# Bounds of T Z + c; the upper end is open where some Tᵢⱼ > 0 meets the open end of a multi-cell axis.
function affine_bounds(grid::Grid, T::Matrix{Float64}, c::Vector{Float64}, Z)
    l, u = low(Z)', high(Z)'
    lo = c + vec(sum(min.(T .* l, T .* u); dims=2))
    hi = c + vec(sum(max.(T .* l, T .* u); dims=2))
    for i in axes(T, 1)
        any(j -> T[i, j] > 0 && grid.numCells[j] > 1, axes(T, 2)) && (hi[i] -= open_margin(grid.granularity))
    end
    return lo, hi
end

# Bounds of T Z + c after a whole period; axes with varying input are untracked.
function end_bounds(ctx, dyn::Dynamics, (T, c)::Affine, Z)
    Φ, V = dyn.Φ[ctx.steps+1], dyn.V[ctx.steps+1]
    lo, hi = affine_bounds(ctx.grid, Φ * T, Φ * c + LazySets.center(V), Z)
    for i in eachindex(lo)
        iszero(view(genmat(V), i, :)) || ((lo[i], hi[i]) = (-Inf, Inf))
    end
    return lo, hi
end

# Bounds after x := M x + j with D′ to D′ + W′ steps left; nothing if no state can jump. Only resets and same-rate translations keep bounds.
# The most time left needs a jump at the piece's start, which the piece before covers unless the state is in the guard at_start.
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
            for (a, c) in zip(guard...) # a guard halfspace on this axis alone bounds the state at the jump
                (a[i] != 0 && count(!iszero, a) == 1) || continue
                a[i] > 0 ? (h = min(h, c / a[i] + most)) : (l = max(l, c / a[i] + least))
            end
            l > h + MODAL_TOL && return nothing
            lo′[i], hi′[i] = l + j[i], h + j[i]
        end
    end
    return lo′, hi′
end

# keep_lo, keep_hi of touching_cells for the bounds; empty below the grid.
function cell_bounds(grid::Grid{D}, (lo, hi)::Bounds) where D
    index(x, i) = floor(Int, clamp((x - grid.lower[i]) / grid.granularity + MODAL_TOL, -1, grid.numCells[i])) + 1
    keep_lo = ntuple(i -> index(lo[i], i), D)
    return keep_lo, ntuple(i -> max(index(hi[i], i), keep_lo[i]), D)
end

# Whether the box meets the guard: exact for a separating halfspace, else isdisjoint (touching counts).
function box_meets((lo, hi)::Bounds, guard::Halfspaces)
    any(((a, b),) -> sum(a[i] * (a[i] > 0 ? lo[i] : hi[i]) for i in eachindex(a)) > b + MODAL_TOL, zip(guard...)) && return false
    return !isdisjoint(Hyperrectangle(low=lo, high=max.(lo, hi)), HPolyhedron([LazySets.HalfSpace(a, b) for (a, b) in zip(guard...)]))
end

# Whether the interior of the box [lo, hi] overlaps the region.
function box_overlaps_region(region, lo, hi)
    region === nothing && return true
    return all(((a, b),) -> sum(a[i] * (a[i] > 0 ? lo[i] : hi[i]) for i in eachindex(a)) < b - MODAL_TOL, zip(region...))
end

# Whole-period discretization of shape G, from which any box of that shape is discretized (CellDiscretization).
period_discretization(dyn::Dynamics, U, G::Matrix{Float64}, δ⁻::Float64, δ⁺::Float64) =
    CellDiscretization(Zonotope(zeros(size(G, 1)), G), (; δ⁻, δ⁺, dyn.A, dyn.P2A_abs, phiDict=dyn.step_Φ, inputDict=dyn.step_input, U), Val(size(G, 1)))
