using LazySets, LinearAlgebra

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
        numCells[i] = ceil((ub-lb)/granularity)
    end

    # Possible to check if cells actually are within the convexSet, but for now we assume they are.
    array = Array{Cell}(undef, (numCells...))  # Create an array to hold the grid cells

    deadCells = BitArray(undef, (numCells...))  # Create an array to label the dead cells

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

struct Cell
    id::CartesianIndex
    sidx::Vector{Int}
    uidx::Vector{Int}
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

Base.convert(::Type{HPolytope}, poly::AAPolytope) = HPolytope(collect(Iterators.flatten((HalfSpace(poly.normals[i], 1.0), HalfSpace(poly.normals[i+1], -1.0)) for i in 1:2:(length(poly.normals)-1))))

function initialize_safe_cells!(grid::Grid, safeSetDict::Dict{CartesianIndex,Vector{LazySet}}=Dict{CartesianIndex,Vector{LazySet}}())
    diag = Diagonal(fill(grid.granularity, grid.dimension))
    ldiag = Diagonal(grid.lower)
    for cell in grid
        offset = car2vec(cell.id)
        pos = collect(eachcol(Diagonal(offset * grid.granularity) + ldiag))
        neg = collect(eachcol(Diagonal((offset .- 1) * grid.granularity) + ldiag))

        safeSet = AAPolytope(collect(Iterators.flatten(zip(pos, neg))))

    end
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
    difTuple = CartesianIndex(NTuple{length(state),Int64}(LinearAlgebra.BLAS.scal(1 / grid.granularity, state - grid.lower)))
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

    lower_bounds, upper_bounds = LazySets.low(convexSet), LazySets.high(convexSet)


    lower_bounds = Int.(floor.(abs.(max.(lower_bounds, grid.lower) .- grid.lower) ./ grid.granularity))
    upper_bounds = Int.(ceil.(abs.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity)) #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity

    #steps = Int.(upper_bounds .- lower_bounds)
    ranges = [lower_bounds[i]:upper_bounds[i] for i in 1:grid.dimension]

    idxs = CartesianIndices((ranges...,))

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


S = Zonotope(zeros(2), [1.0 0.0; 0.0 1.0])  # Example zonotope in 2D
granularity = 0.25  # Example granularity
grid = Grid(S, granularity)

#initialize_safe_cells!(grid)  # Initialize the safe cells in the grid

@time box(grid, [0.5, 0.5])  # Example state to find the corresponding cell
@time fastbox(grid, [0.5, 0.5])  # Example state to find the corresponding cell using fastbox
#get_cell_bounds(grid, box(grid, [-0., -1.]))  # Example to get the bounds of the corresponding cell
#@show get_touching_cells(grid, Zonotope([0., 0.], [0.5 0.0; 0.0 0.5]))
#@show get_contained_cells(grid, Zonotope([0., 0.], [0.5 0.0; 0.0 0.5]))

#AA = AAPolytope([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]])  # Example axis-aligned polytope in 2D
#@show convert(HPolytope, AA)  # Convert the axis-aligned polytope to an HPolytope




# Måske muligt i stedet for at genbruge koden fra Astrid, at bruge den samme funktion til at lave en grid som bare er en store af keys og så gemme values et andet sted. Hvis det er implementeret med et linært index
# eller som en dictionary på cartisianIndex ville man nok kunne fjerne keys som er cell'er der ikke længere har safe elementer. 
# Så ville man stadig kunne bruge bounds funktioner til at finde de tætteste celler til constraints.
# Man kunne tjekke at når en cell ikke længere er safe så fjerner man de eventuelle constraints den har medført og bytter dem ud med cell da den nok er simplere