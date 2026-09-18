using LazySets, LinearAlgebra, ReachabilityAnalysis

### Following code is adapted from "https://github.com/AstridHornBrorholt/Shielded-Learning-for-Hybrid-Systems/blob/main/Shared%20Code/Squares.jl"
# Assumed that partitioning is axis-aligned, i.e., the partitioning is done along the axes of the state space. The Grid struct represents a grid in the state space with specified granularity and bounds.
struct Grid{T<:Real}
    dimension::Int
    granularity::T
    lower::Vector{T}
    upper::Vector{T}
    numCells::Vector{Int}
    deadCells::BitArray
    array
end

function Grid(convexSet::LazySet, granularity::T) where T<:Real
    dimension = LazySets.dim(convexSet)

    lower_bounds = LazySets.low(convexSet)
    upper_bounds = LazySets.high(convexSet)

    numCells = zeros(Int, dimension)

    for (i, (lb, ub)) in enumerate(zip(lower_bounds, upper_bounds))
        numCells[i] = max(ceil((ub-lb)/granularity), 1)
    end

    # Possible to check if cells actually are within the convexSet, but for now we assume they are.
    array = Array{Cell}(undef, (numCells...))  # Create an array to hold the grid cells

    deadCells = falses(numCells...)#BitArray(undef, (numCells...))  # Create an array to label the dead cells

    for id in CartesianIndices(array)
        array[id] = Cell(id, [LinearIndices(array)[id]], [])
    end

    return Grid{T}(dimension, granularity, lower_bounds, upper_bounds, numCells, deadCells, array)
end

Base.show(io::IO, grid::Grid) = println(io,
    "Grid($(grid.granularity), $(grid.lower), $(grid.upper))")

# Makes the grid iterable, returning each square in turn.
Base.length(grid::Grid) = length(grid.array)

Base.size(grid::Grid) = size(grid.array)

mutable struct Cell
    id::CartesianIndex
    sidx::Vector{Int}
    uidx::Vector{Int}
    pCells::Vector{CartesianIndex}
end

function Cell(id, sidx, uidx)
    return Cell(id, sidx, uidx, [])
end

Base.show(io::IO, cell::Cell) = println(io,
    "Cell($(cell.id), $(cell.sidx), $(cell.uidx))")



Base.IteratorSize(grid::Type{Grid}) = length(grid.array)

Base.iterate(grid::Grid) = begin
    idx = firstindex(grid.array)
    cell = grid.array[idx]
    return cell, idx
end


Base.iterate(grid::Grid, state) = begin
    indices = copy(state) + 1
    if state == lastindex(grid.array)
        return nothing
    else
        return grid.array[indices], indices
    end
end

car2vec(x::CartesianIndex) = collect(Tuple(x))

# The AAPolytope struct represents an axis-aligned polytope defined by its normals. It is used to represent the safe set in the state space.
struct AAPolytope
    normals::Vector{Vector{Float64}}
end

Base.convert(::Type{HPolytope}, poly::AAPolytope) = HPolytope(collect(Iterators.flatten((LazySets.HalfSpace(poly.normals[i], 1.0), LazySets.HalfSpace(poly.normals[i+1], -1.0)) for i in 1:2:(length(poly.normals)-1))))

#=
function initialize_grid_cells!(grid::Grid)
    for cell in grid
        grid.array[cell.id] = cell
    end
end
=#

function mark_dead_cells!(grid::Grid, unsafeSet::LazySet, unsafeSetDict::Dict{CartesianIndex,Vector{LazySet}}=Dict{CartesianIndex,Vector{LazySet}}())
    containedCells, edgeCells = get_contained_edge_cells(grid, unsafeSet)
    for containedCell in containedCells
        grid.deadCells[containedCell.id] = true
    end
    ldiag = Diagonal(grid.lower)
    for edgeCell in edgeCells
        # Two cases for such that multiple unsafeSets can be added. 
        if haskey(unsafeSetDict, edgeCell.id)
            push!(edgeCell.uidx, length(edgeCell.uidx)+1)
            offset = car2vec(edgeCell.id)
            granularityOffset = offset * grid.granularity
            lowerOffset = grid.lower
            hs = vcat(collect([LazySets.HalfSpace(a, dot(a, granularityOffset)+dot(a, lowerOffset)), LazySets.HalfSpace(-a, -dot(a, granularityOffset)-dot(a, lowerOffset)+grid.granularity)] for a in eachcol(diagm(ones(grid.dimension))))...)
            safeSet = HPolytope(hs)
            push!(unsafeSetDict[edgeCell.id], LazySets.API.intersection(unsafeSet, safeSet))
        else
            edgeCell.uidx = [1]
            offset = car2vec(edgeCell.id)
            granularityOffset = offset * grid.granularity
            lowerOffset = grid.lower
            hs = vcat(collect([LazySets.HalfSpace(a, dot(a, granularityOffset)+dot(a, lowerOffset)), LazySets.HalfSpace(-a, -dot(a, granularityOffset)-dot(a, lowerOffset)+grid.granularity)] for a in eachcol(diagm(ones(grid.dimension))))...)
            safeSet = HPolytope(hs)
            unsafeSetDict[edgeCell.id] = [LazySets.API.intersection(unsafeSet, safeSet)]
        end
    end
end

function initialize_safe_cells!(grid::Grid, safeSetDict::Dict{CartesianIndex,Vector{LazySet}}=Dict{CartesianIndex,Vector{LazySet}}())
    diag = Diagonal(fill(grid.granularity, grid.dimension))
    ldiag = Diagonal(grid.lower)
    for cell in grid
        #=offset = car2vec(cell.id)
        pos = collect(eachcol(Diagonal(offset * grid.granularity) + ldiag))
        neg = collect(eachcol(Diagonal((offset .- 1) * grid.granularity) + ldiag))

        safeSet = AAPolytope(collect(Iterators.flatten(zip(pos, neg))))=#
        offset = car2vec(cell.id)
        granularityOffset = offset * grid.granularity
        lowerOffset = grid.lower
        hs = vcat(collect([LazySets.HalfSpace(a, dot(a, granularityOffset)+dot(a, lowerOffset)), LazySets.HalfSpace(-a, -dot(a, granularityOffset)-dot(a, lowerOffset)+grid.granularity)] for a in eachcol(diagm(ones(grid.dimension))))...)
        safeSet = HPolytope(hs)
        safeSetDict[cell.id] = [safeSet]
    end
end

function car_to_HPolytope(grid, car)
    offset = car2vec(car)
    granularityOffset = offset * grid.granularity
    lowerOffset = grid.lower
    hs = vcat(collect([LazySets.HalfSpace(a, dot(a, granularityOffset)+dot(a, lowerOffset)), LazySets.HalfSpace(-a, -dot(a, granularityOffset)-dot(a, lowerOffset)+grid.granularity)] for a in eachcol(diagm(ones(grid.dimension))))...)
    return HPolytope(hs)
end

function hyperrectangle_to_HPolytope(hyperrectangle)
    i = ones(length(hyperrectangle.center))
    D = diagm(i)
    hs = vcat(collect([LazySets.HalfSpace(a, dot(hyperrectangle.center, a) + dot(hyperrectangle.radius, a)), LazySets.HalfSpace(-a, dot(hyperrectangle.center, -a) + dot(hyperrectangle.radius, a))] for a in eachcol(D))...)
    return HPolytope(hs)
end

function hyperrectangle_to_HPolyhedron(hyperrectangle)
    i = ones(length(hyperrectangle.center))
    D = diagm(i)
    hs = vcat(collect([LazySets.HalfSpace(a, dot(hyperrectangle.center, a) + dot(hyperrectangle.radius, a)), LazySets.HalfSpace(-a, dot(hyperrectangle.center, -a) + dot(hyperrectangle.radius, a))] for a in eachcol(D))...)
    return HPolyhedron(hs)
end

function initialize_zonotope_array(grid::Grid)
    zonotopeArray = Array{Zonotope}(undef, (grid.numCells...))
    generators = diagm(fill(grid.granularity / 2, grid.dimension))
    for cell in grid
        offset = car2vec(cell.id)
        center = (offset .- 0.5) * grid.granularity .+ grid.lower
        zonotopeArray[cell.id] = Zonotope(center, generators)

    end

    return zonotopeArray
end

function box(grid::Grid, state)
    indices = zeros(Int, grid.dimension)

    for i in 1:grid.dimension
        if !(grid.lower[i] <= state[i] < grid.upper[i])
            throw(ArgumentError("State is out of bounds of the grid."))
        end

        indices[i] = floor(Int, (state[i] - grid.lower[i]) / grid.granularity) + 1
    end
    return grid.array[CartesianIndex(indices...)]

end

function fastbox(grid::Grid, state)
    #t = zeros(Int, length(state))
    #axpy!(1 / grid.granularity, state - grid.lower, t)
    #difTuplet = CartesianIndex(NTuple{length(state),Int64}(t))
    #@show t
    #@show difTuplet
    @show map(x -> floor(Int64, x), (LinearAlgebra.BLAS.scal(1 / grid.granularity, state - grid.lower)))
    #t=NTuple{length(state),Integer}(LinearAlgebra.BLAS.scal(1 / grid.granularity, state - grid.lower))
    t=NTuple{length(state),Integer}(map(x -> floor(Int64, x), (LinearAlgebra.BLAS.scal(1 / grid.granularity, state - grid.lower))))
    difTuple = CartesianIndex(t)
    #@show difTuple

    try
        return grid.array[difTuple]

        #return grid.array[CartesianIndex(difTuple)]
    catch
        throw(ArgumentError("State is out of bounds of the grid."))
    end
end

function get_cell_bounds(grid::Grid, cell::Cell)
    lb = grid.lower + (car2vec(cell.id) .- 1) .* grid.granularity
    ub = lb .+ grid.granularity


    return lb, ub
end

function get_touching_cells(grid::Grid, convexSet::LazySet)
    touching_cells = []

    lower_bounds, upper_bounds = clamp.(LazySets.low(convexSet), grid.lower, grid.upper), clamp.(LazySets.high(convexSet), grid.lower, grid.upper)
    #lower_bounds = LazySets.low(convexSet)
    #lower_bounds = Int.(floor.(abs.(max.(lower_bounds, grid.lower) .- grid.lower) ./ grid.granularity) .+ 1)
    lower_bounds = Int.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity) .+ 1)

    upper_bounds = Int.(ceil.(abs.(upper_bounds .- grid.lower) ./ grid.granularity)) #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity
    #@show upper_bounds
    ranges = [lower_bounds[i]:max(upper_bounds[i], 1) for i in 1:grid.dimension]

    idxs = CartesianIndices((ranges...,))
    #@show idxs
    for idx in idxs
        cell = grid.array[idx]
        lower_bounds, upper_bounds = get_cell_bounds(grid, cell)
        cell_box = Hyperrectangle((lower_bounds + upper_bounds) / 2, (upper_bounds - lower_bounds) / 2)

        if !isempty(intersect(convexSet, cell_box))
            push!(touching_cells, cell)
        end
    end

    return touching_cells
end

function get_contained_cells(grid::Grid, convexSet::LazySet)
    contained_cells = []


    for cell in get_touching_cells(grid, convexSet)
        lower_bounds, upper_bounds = get_cell_bounds(grid, cell)
        cell_box = Hyperrectangle((lower_bounds + upper_bounds) / 2, (upper_bounds - lower_bounds) / 2)

        if issubset(cell_box, convexSet)
            push!(contained_cells, cell)

        end
    end

    return contained_cells
end

function get_contained_edge_cells(grid::Grid, convexSet::LazySet)
    contained_cells = []
    perimeter_cells = []

    for cell in get_touching_cells(grid, convexSet)
        lower_bounds, upper_bounds = get_cell_bounds(grid, cell)
        cell_box = Hyperrectangle((lower_bounds + upper_bounds) / 2, (upper_bounds - lower_bounds) / 2)

        if issubset(cell_box, convexSet)
            push!(contained_cells, cell)
        else
            push!(perimeter_cells, cell)
        end
    end

    return contained_cells, perimeter_cells
end

function get_contained_perimeter_cells(grid::Grid, convexSet::LazySet)
    contained_cells = []
    touching_cells = get_touching_cells(grid, convexSet)
    for cell in touching_cells
        lower_bounds, upper_bounds = get_cell_bounds(grid, cell)
        cell_box = Hyperrectangle((lower_bounds + upper_bounds) / 2, (upper_bounds - lower_bounds) / 2)

        if issubset(cell_box, convexSet)
            push!(contained_cells, cell)

        end
    end
    touching_cell_neighbour_ids = grow_indices(grid, collect(car2vec(c.id) for c in touching_cells))
    perimeter_cell_ids = setdiff(touching_cell_neighbour_ids, c.id for c in contained_cells)
    return collect(grid.array[CartesianIndex(x)] for x in perimeter_cell_ids)
end

function offsets(grid::Grid)
    offsets = [Tuple(id) for id in CartesianIndices(([1:1:3 for i in 1:grid.dimension]...,))]
    offsets = [t .- Tuple(fill(2, grid.dimension)) for t in offsets]
    return delete!(Set(offsets), Tuple(fill(0, grid.dimension)))
end

function grow_indices(grid::Grid, idxs)
    offset = offsets(grid)
    res = []
    for idx in idxs
        union!(res, collect(Tuple(idx .+ of) for of in offset))
    end

    intersect!(res, union(Tuple(idx) for idx in CartesianIndices(grid.array)))
    setdiff!(res, collect(Tuple(id) for id in idxs))
    return res
end

function grow_indices(grid::Grid, idxs, offset)
    res = []
    for idx in idxs
        union!(res, collect(Tuple(idx .+ of) for of in offset))
    end

    intersect!(res, union(Tuple(idx) for idx in CartesianIndices(grid.array)))
    setdiff!(res, collect(Tuple(id) for id in idxs))
    return res
end


