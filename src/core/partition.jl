### Following code is adapted from "https://github.com/AstridHornBrorholt/Shielded-Learning-for-Hybrid-Systems/blob/main/Shared%20Code/Squares.jl"
# Assumed that partitioning is axis-aligned, i.e., the partitioning is done along the axes of the state space. The Grid struct represents a grid in the state space with specified granularity and bounds.
# N is the state dimension. It is a type parameter (not only the `dimension` field) so that functions taking a Grid{N}
# can use CartesianIndex{N}, BitArray{N} and NTuple{N} as concrete types. Inside such functions use N, never
# grid.dimension, in type positions (e.g. CartesianIndex{N}[] and ntuple(f, N)).
struct Grid{N}
    dimension::Int
    granularity::Float64
    lower::Vector{Float64}
    upper::Vector{Float64}
    numCells::Vector{Int}
    deadCells::BitArray{N}
end

# Number of cells per dimension as an NTuple{N,Int}. Use this when allocating per-cell arrays (falses(cellshape(grid)))
# instead of splatting grid.numCells, whose length is unknown to the compiler and makes the array dimension unknown.
cellshape(grid::Grid{N}) where N = size(grid.deadCells)::NTuple{N,Int}

function Grid(set::LazySet, granularity::Float64)
    dimension = LazySets.dim(set)

    lower_bounds = LazySets.low(set)
    upper_bounds = LazySets.high(set)

    numCells = zeros(Int, dimension)

    lower_bounds = Vector{Float64}(lower_bounds)
    upper_bounds = Vector{Float64}(upper_bounds)
    for (i, (lb, ub)) in enumerate(zip(lower_bounds, upper_bounds))
        numCells[i] = max(ceil((ub-lb)/granularity), 1)
        # A dimension narrower than a cell (zero width for a constant like gravity) becomes a single cell centred on it, so
        # the cell centre (lower + granularity/2) is the middle of the state space there and touching-cell index ranges
        # are 1:1 rather than empty. The shields give such a cell the width of the state space on that axis.
        if ub - lb < granularity
            lower_bounds[i] = (lb + ub)/2 - granularity/2
            upper_bounds[i] = (lb + ub)/2 + granularity/2
        end
    end

    deadCells = falses(numCells...)
    return Grid(dimension, granularity, lower_bounds, upper_bounds, numCells, deadCells)
end

Base.show(io::IO, grid::Grid) = println(io,
    "Grid($(grid.granularity), $(grid.lower), $(grid.upper), $(grid.numCells))")

@inline index_vector(x::CartesianIndex) = collect(Tuple(x))

@inline float_index_vector(x::CartesianIndex) = collect(Float64, Tuple(x))

function hyperrectangle_to_HPolyhedron(hyperrectangle)
    i = ones(length(hyperrectangle.center))
    D = diagm(i)
    hs = vcat(collect([LazySets.HalfSpace(a, dot(hyperrectangle.center, a) + dot(hyperrectangle.radius, a)), LazySets.HalfSpace(-a, dot(hyperrectangle.center, -a) + dot(hyperrectangle.radius, a))] for a in eachcol(D))...)
    return HPolyhedron(hs)
end

# The cells on axis i whose closed box meets [lo, hi], as first:last index. A bound on a grid line meets the cells on
# both sides of it, and CELL_BOUND_TOL (in cells) keeps a bound that rounding put just beside a grid line from losing
# the cell past it. The cells' membership tests then decide among these candidates.
const CELL_BOUND_TOL = 1e-9
@inline function cell_range(grid::Grid, i::Int, lo::Real, hi::Real)
    x_lo = (clamp(lo, grid.lower[i], grid.upper[i]) - grid.lower[i]) / grid.granularity
    x_hi = (clamp(hi, grid.lower[i], grid.upper[i]) - grid.lower[i]) / grid.granularity
    return max(ceil(Int, x_lo - CELL_BOUND_TOL), 1):min(floor(Int, x_hi + CELL_BOUND_TOL) + 1, grid.numCells[i])
end

function touching_cells(grid::Grid{N}, set::LazySet) where N
    touching_cell_idxs = CartesianIndex{N}[]

    ranges = [cell_range(grid, i, l, h) for (i, l, h) in zip(1:N, LazySets.low(set), LazySets.high(set))]
    lower_bounds, upper_bounds = first.(ranges), last.(ranges)

    sweeping_dim = argmax(upper_bounds[i]-lower_bounds[i] for i in 1:grid.dimension)
    hbox = Hyperrectangle(grid.lower .+ (grid.granularity/2), fill(grid.granularity/2, N))
    if upper_bounds[sweeping_dim] - lower_bounds[sweeping_dim] > 1

        cartesian_offset = CartesianIndex(ntuple(i -> i == sweeping_dim ? 1 : 0, N))

        sweeping_range = range(lower_bounds[sweeping_dim], upper_bounds[sweeping_dim])
        rev_sweeping_range = range(upper_bounds[sweeping_dim], lower_bounds[sweeping_dim]; step=-1)
        for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], lower_bounds[i])), N))
            lower_idx = 0
            upper_idx = 0
            for off in sweeping_range
                of = (index_vector((idx + (cartesian_offset * off))) .- 1) .* grid.granularity

                LazySets.API.translate!(hbox, of)
                if !isdisjoint(set, hbox)
                    LazySets.API.translate!(hbox, -of)

                    lower_idx = off
                    break

                end
                LazySets.API.translate!(hbox, -of)
            end
            if lower_idx == 0
                continue
            end
            for off in rev_sweeping_range

                of = (index_vector(idx + (cartesian_offset * off)) .- 1) .* grid.granularity

                LazySets.API.translate!(hbox, of)
                if !isdisjoint(set, hbox)
                    LazySets.API.translate!(hbox, -of)

                    upper_idx = off
                    break

                end
                LazySets.API.translate!(hbox, -of)
            end
            if upper_idx == 0
                continue
            end

            for elem in lower_idx:upper_idx

                push!(touching_cell_idxs, idx + (cartesian_offset * elem))
            end

        end
    else
        println("second")

        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:upper_bounds[i], N))
            of = (index_vector(idx) .- 1) .* grid.granularity
            LazySets.API.translate!(hbox, of)
            if !isdisjoint(set, hbox)
                push!(touching_cell_idxs, idx)
            end
            LazySets.API.translate!(hbox, -of)
        end

    end
    return unique(touching_cell_idxs)
end

# A zonotope whose generators never change, only its center (e.g. the reach set of a fast cell bloated by half a
# cell, max_zonotope). Its H-representation is computed once, so testing whether a point lies in a translated copy is
# a few dot products instead of an LP.
struct ConstantZonotopeCover{N}
    normals::Matrix{Float64} # N x n_facets, one facet normal per column, for the zonotope centred at the origin
    offsets::Vector{Float64} # facet i is normals[:, i]' * y <= offsets[i]
    radius::NTuple{N,Float64} # half-widths of the interval hull
    tol::Float64
end

function ConstantZonotopeCover(Z::Zonotope, ::Val{N}; tol::Float64=1e-7) where N
    G = genmat(Z)
    facets = constraints_list(Zonotope(zeros(N), G))
    normals = Matrix{Float64}(undef, N, length(facets))
    offsets = Vector{Float64}(undef, length(facets))
    for (i, h) in enumerate(facets)
        normals[:, i] .= h.a
        offsets[i] = h.b
    end
    radius = ntuple(j -> sum(abs, view(G, j, :)), N)
    return ConstantZonotopeCover{N}(normals, offsets, radius, tol)
end

# Cover for the zonotope with generator matrix G (N x m, full row rank), centred at the origin, with the given
# interval-hull radius (used only to bound the candidate cells). The facets are enumerated directly instead of through
# LazySets: every set of N-1 generators spanning a hyperplane gives a normal n (generalised cross product), and the
# pair of facets +-n' y <= sum_j |n' g_j|. A generator parallel to an earlier one only gives facets with the normal of
# a facet that comes earlier (or none), so normals come from the first generator of each direction only, while every
# offset still sums over all generators. Sets after a jump often have many parallel generators (28 with 19 parallel is
# typical for the cruise control), so this cuts the facets about tenfold.
function ConstantZonotopeCover(G::AbstractMatrix{Float64}, radius::NTuple{N,Float64}, ::Val{N}; tol::Float64=1e-7) where N
    m = size(G, 2)
    reps = Int[] # the first generator of each direction
    for j in 1:m
        nj = norm(view(G, :, j))
        nj == 0.0 && continue
        parallel = false
        for i in reps
            if abs(abs(dot(view(G, :, i), view(G, :, j))) - norm(view(G, :, i)) * nj) <= 1e-12 * norm(view(G, :, i)) * nj
                parallel = true
                break
            end
        end
        parallel || push!(reps, j)
    end
    r = length(reps)
    r < N - 1 && return ConstantZonotopeCover{N}(zeros(N, 0), Float64[], radius, tol)
    pairs = binomial(r, N - 1)
    normals = Matrix{Float64}(undef, N, 2pairs)
    offsets = Vector{Float64}(undef, 2pairs)
    count = 0
    n = zeros(N)
    minor = zeros(N - 1, N - 1)
    subset = collect(1:N-1)
    while true
        if N == 2
            n[1] = G[2, reps[subset[1]]]
            n[2] = -G[1, reps[subset[1]]]
        elseif N == 3
            # the cross product of the two generators
            i, j = reps[subset[1]], reps[subset[2]]
            n[1] = G[2, i] * G[3, j] - G[2, j] * G[3, i]
            n[2] = -(G[1, i] * G[3, j] - G[1, j] * G[3, i])
            n[3] = G[1, i] * G[2, j] - G[1, j] * G[2, i]
        else
            # n_k = (-1)^(k+1) det(G[rows != k, subset])
            for k in 1:N
                for c in 1:N-1, r in 1:N-1
                    minor[r, c] = G[r < k ? r : r + 1, reps[subset[c]]]
                end
                n[k] = (isodd(k) ? 1.0 : -1.0) * det(minor)
            end
        end
        len = norm(n)
        if len > 1e-12
            b = 0.0
            for j in 1:m
                s = 0.0
                for r in 1:N
                    s += (n[r] / len) * G[r, j]
                end
                b += abs(s)
            end
            for r in 1:N
                normals[r, count+1] = n[r] / len
                normals[r, count+2] = -(n[r] / len)
            end
            offsets[count+1] = b
            offsets[count+2] = b
            count += 2
        end
        # Next (N-1)-combination of 1:r in lexicographic order.
        i = N - 1
        while i >= 1 && subset[i] == r - (N - 1) + i
            i -= 1
        end
        if i < 1
            break
        end
        subset[i] += 1
        for k in i+1:N-1
            subset[k] = subset[k-1] + 1
        end
    end
    return ConstantZonotopeCover{N}(normals[:, 1:count], resize!(offsets, count), radius, tol)
end

# Number of facet pairs ConstantZonotopeCover(G, ...) enumerates for m generators in N dimensions.
cover_facet_pairs(m::Int, N::Int) = binomial(m, N - 1)
const MAX_COVER_FACET_PAIRS = 5000

# True if the point (centre of cell idx) lies in the cover zonotope translated to center c.
@inline function cell_center_in_cover(cov::ConstantZonotopeCover{N}, grid::Grid{N}, idx::CartesianIndex{N}, c::AbstractVector{Float64}) where N
    for i in eachindex(cov.offsets)
        s = 0.0
        for j in 1:N
            s += cov.normals[j, i] * ((grid.lower[j] + grid.granularity/2) + (idx[j] - 1) * grid.granularity - c[j])
        end
        if s > cov.offsets[i] + cov.tol * max(1.0, abs(cov.offsets[i]))
            return false
        end
    end
    return true
end

# Same cells as touching_cells_lp(grid, Z) for the translated zonotope Z with center c, using
# the precomputed H-representation instead of an LP per candidate cell. A cell touches the reach set iff its centre
# lies in the reach set bloated by half a cell, which is what cov represents.
# keep_lo and keep_hi (cell indices per axis) limit the result to the cells in that box, as filtering it afterwards would,
# without testing the cells outside it.
function cells_in_cover(grid::Grid{N}, cov::ConstantZonotopeCover{N}, c::AbstractVector{Float64},
    keep_lo::NTuple{N,Int}=ntuple(_ -> 1, N), keep_hi::NTuple{N,Int}=cellshape(grid)) where N
    touching_cell_idxs = CartesianIndex{N}[]
    full = ntuple(i -> cell_range(grid, i, c[i] - cov.radius[i], c[i] + cov.radius[i]), N)
    full_lower, full_upper = map(first, full), map(last, full)

    # The way to sweep depends on the full candidate box, so the result is the same with or without the limits.
    sweeping_dim = argmax(full_upper[i] - full_lower[i] for i in 1:N)
    next_sweep = argmax(i == sweeping_dim ? 0 : full_upper[i] - full_lower[i] for i in 1:N)
    lower_bounds = ntuple(i -> max(full_lower[i], keep_lo[i]), N)
    upper_bounds = ntuple(i -> min(full_upper[i], keep_hi[i]), N)
    any(i -> lower_bounds[i] > upper_bounds[i], 1:N) && return touching_cell_idxs

    if full_upper[next_sweep] == full_lower[next_sweep]
        append!(touching_cell_idxs, CartesianIndices(ntuple(i -> lower_bounds[i]:upper_bounds[i], N)))
    elseif full_upper[sweeping_dim] - full_lower[sweeping_dim] > 1
        # The cover is convex, so along each line in the sweeping dimension the touched cells form one interval:
        # find its first and last cell and fill in between.
        cartesian_offset = CartesianIndex(ntuple(i -> i == sweeping_dim ? 1 : 0, N))
        for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], lower_bounds[i])), N))
            lower_idx = 0
            for off in lower_bounds[sweeping_dim]:upper_bounds[sweeping_dim]
                if cell_center_in_cover(cov, grid, idx + cartesian_offset * off, c)
                    lower_idx = off
                    break
                end
            end
            if lower_idx == 0
                continue
            end
            upper_idx = 0
            for off in upper_bounds[sweeping_dim]:-1:lower_bounds[sweeping_dim]
                if cell_center_in_cover(cov, grid, idx + cartesian_offset * off, c)
                    upper_idx = off
                    break
                end
            end
            for elem in lower_idx:upper_idx
                push!(touching_cell_idxs, idx + cartesian_offset * elem)
            end
        end
    else
        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:upper_bounds[i], N))
            if cell_center_in_cover(cov, grid, idx, c)
                push!(touching_cell_idxs, idx)
            end
        end
    end
    return touching_cell_idxs
end

# Cells touching the zonotope set, without LPs: a cell touches it iff the cell centre lies in set bloated by
# half a cell, which is tested against the facets of that bloated zonotope (see ConstantZonotopeCover). Gives the same
# cells as touching_cells_lp. Falls back to it when the facet count explodes (many generators in high
# dimensions). Without LPs this is also safe to call from several threads (GLPK is not).
function touching_cells(grid::Grid{N}, set::Zonotope{Float64,Vector{Float64},Matrix{Float64}}) where N
    G = genmat(set)
    nonzero = [j for j in axes(G, 2) if any(!iszero, view(G, :, j))]
    if cover_facet_pairs(length(nonzero) + N, N) > MAX_COVER_FACET_PAIRS
        return touching_cells_lp(grid, set)
    end
    # Bound the candidate cells by the set itself (as the LP version does); test membership against the bloated set.
    radius = ntuple(i -> sum(abs, view(G, i, :)), N)
    bloated = hcat(G[:, nonzero], Matrix((grid.granularity / 2) * I, N, N))
    return cells_in_cover(grid, ConstantZonotopeCover(bloated, radius, Val(N)), set.center)
end

# Covers by their generator matrix (the set's nonzero generators and the half-cell bloat), so sets that only differ by a
# translation share one. Emptied when it reaches MAX_CACHED_COVERS.
const CoverCache{N} = Dict{Matrix{Float64},ConstantZonotopeCover{N}}
const MAX_CACHED_COVERS = 100_000

# touching_cells with the cover taken from, or added to, covers. Gives the same cells. keep_lo and keep_hi
# limit the result to a box of cells, as in cells_in_cover.
function touching_cells(grid::Grid{N}, set::Zonotope{Float64,Vector{Float64},Matrix{Float64}}, covers::CoverCache{N},
    keep_lo::NTuple{N,Int}=ntuple(_ -> 1, N), keep_hi::NTuple{N,Int}=cellshape(grid)) where N
    G = genmat(set)
    nonzero = [j for j in axes(G, 2) if any(!iszero, view(G, :, j))]
    if cover_facet_pairs(length(nonzero) + N, N) > MAX_COVER_FACET_PAIRS
        return touching_cells_lp(grid, set)
    end
    bloated = hcat(G[:, nonzero], Matrix((grid.granularity / 2) * I, N, N))
    cov = get(covers, bloated, nothing)
    if cov === nothing
        length(covers) >= MAX_CACHED_COVERS && empty!(covers)
        cov = covers[bloated] = ConstantZonotopeCover(bloated, ntuple(i -> sum(abs, view(G, i, :)), N), Val(N))
    end
    return cells_in_cover(grid, cov, set.center, keep_lo, keep_hi)
end

function touching_cells_lp(grid::Grid{N}, set::Zonotope{Float64,Vector{Float64},Matrix{Float64}}) where N
    touching_cell_idxs::Vector{CartesianIndex{N}} = CartesianIndex{N}[]

    ranges = [cell_range(grid, i, l, h) for (i, l, h) in zip(1:N, LazySets.low(set), LazySets.high(set))]
    lower_bounds::Vector{Int64}, upper_bounds::Vector{Int64} = first.(ranges), last.(ranges)

    sweeping_dim::Int64 = argmax(upper_bounds[i]-lower_bounds[i] for i in 1:grid.dimension)
    next_sweep = argmax(i == sweeping_dim ? 0 : upper_bounds[i]-lower_bounds[i] for i in 1:grid.dimension)

    #hbox::LazySets.HyperrectangleModule.Hyperrectangle{Float64,Vector{Float64},Vector{Float64}} = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, N))
    hbox = Zonotope(grid.lower .+ (grid.granularity/2), diagm(fill(grid.granularity/2, N)))
    of::Vector{Float64} = grid.lower .+ (grid.granularity/2)
    if upper_bounds[next_sweep] == lower_bounds[next_sweep]
        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:1:upper_bounds[i], N))

            push!(touching_cell_idxs, idx)

        end
    else
        if upper_bounds[sweeping_dim] - lower_bounds[sweeping_dim] > 1

            cartesian_offset::CartesianIndex{N} = CartesianIndex(ntuple(i -> i == sweeping_dim ? 1 : 0, N))

            sweeping_range::StepRange{Int64,Int64} = lower_bounds[sweeping_dim]:1:upper_bounds[sweeping_dim]
            rev_sweeping_range::StepRange{Int64,Int64} = upper_bounds[sweeping_dim]:-1:lower_bounds[sweeping_dim]
            for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], lower_bounds[i])), N))

                lower_idx::Int64 = 0
                upper_idx::Int64 = 0
                for off in sweeping_range
                    copyto!(hbox.center, (float_index_vector((idx + (cartesian_offset * off))) .- 1.0) .* grid.granularity + of)
                    if !isdisjoint(set, hbox, false)

                        lower_idx = off
                        break

                    end
                end
                if lower_idx == 0
                    continue
                end
                for off in rev_sweeping_range

                    copyto!(hbox.center, (float_index_vector(idx + (cartesian_offset * off)) .- 1.0) .* grid.granularity + of)

                    if !isdisjoint(set, hbox, false)

                        upper_idx = off
                        break

                    end
                end
                if upper_idx == 0
                    continue
                end

                for elem in lower_idx:upper_idx

                    push!(touching_cell_idxs, idx + (cartesian_offset * elem))
                end

            end

        else
            next_sweep = argmax(i == sweeping_dim ? 0 : upper_bounds[i]-lower_bounds[i] for i in 1:grid.dimension)

            for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:1:upper_bounds[i], N))
                #of::Vector{Float64} = (index_vector(idx) .- 1) .* granularity
                copyto!(hbox.center, (float_index_vector(idx) .- 1.0) .* grid.granularity + of)
                if !isdisjoint(set, hbox, false)
                    push!(touching_cell_idxs, idx)

                end
            end

        end
    end
    return touching_cell_idxs
end
