using LazySets, LinearAlgebra

### Following code is adapted from "https://github.com/AstridHornBrorholt/Shielded-Learning-for-Hybrid-Systems/blob/main/Shared%20Code/Squares.jl"
# Assumed that partitioning is axis-aligned, i.e., the partitioning is done along the axes of the state space. The Grid struct represents a grid in the state space with specified granularity and bounds.
struct Grid{T<:Real}
    dimension::Int
    granularity::T
    lower::Vector{T}
    upper::Vector{T}
    numCells::Vector{Int}
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
    array = zeros(Int8, (numCells...))  # Create an array to hold the grid cells
    return Grid{T}(dimension, granularity, lower_bounds, upper_bounds, numCells, array)
end

Base.show(io::IO, grid::Grid) = println(io,
    "Grid($(grid.granularity), $(grid.lower), $(grid.upper))")

# Makes the grid iterable, returning each square in turn.
Base.length(grid::Grid) = length(grid.array)

Base.size(grid::Grid) = size(grid.array)

struct Cell
    grid::Grid
    indices::Vector{Int}
end

Base.show(io::IO, cell::Cell) = println(io,
    "Cell($(cell.indices))")

Base.iterate(grid::Grid) = begin
    indices = ones(Int, grid.dimension)
    cell = Cell(grid, indices)
    cell, indices
end

Base.iterate(grid::Grid, state) = begin
    indices = copy(state)
    for i in 1:grid.dimension
        indices[i] += 1
        if indices[i] <= grid.numCells[i]
            break
        else
            if i < grid.dimension
                indices[i] = 1
                # Proceed to incrementing next row
            else
                return nothing
            end
        end
    end
    return Cell(grid, indices), indices
end

function box(grid::Grid, state)
    indices = zeros(Int, grid.dimension)

    for i in 1:grid.dimension
        indices[i] += 1
        if !(grid.lower[i] <= state[i] < grid.upper[i])
            throw(ArgumentError("State is out of bounds of the grid."))
        end

        indices[i] = floor(Int, (state[i] - grid.lower[i]) / grid.granularity) + 1
    end

    Cell(grid, indices)
end

function get_cell_bounds(cell::Cell)
    grid = cell.grid
    lower_bounds = zeros(Float64, grid.dimension)
    upper_bounds = zeros(Float64, grid.dimension)

    for i in 1:grid.dimension
        lower_bounds[i] = grid.lower[i] + (cell.indices[i] - 1) * grid.granularity
        upper_bounds[i] = lower_bounds[i] + grid.granularity
    end

    return lower_bounds, upper_bounds
end

function get_touching_cells(grid::Grid, convexSet::LazySet)
    touching_cells = []

    for cell in grid
        lower_bounds, upper_bounds = get_cell_bounds(cell)
        cell_box = Hyperrectangle((lower_bounds + upper_bounds) / 2, (upper_bounds - lower_bounds) / 2)

        if !isempty(intersect(convexSet, cell_box))
            push!(touching_cells, cell)
        end
    end

    return touching_cells
end


#=
S = Zonotope(zeros(2), [1.0 0.0; 0.0 1.0])  # Example zonotope in 2D
granularity = 0.5  # Example granularity
grid = Grid(S, granularity)
box(grid, [-1., -1.])  # Example state to find the corresponding cell
get_cell_bounds(box(grid, [-1., -1.]))  # Example to get the bounds of the corresponding cell
get_touching_cells(grid, Zonotope([-0.5, -0.5], [0.5 0.0; 0.0 0.5]))  # Example to get all cells that touch the convex set
=#


# Måske muligt i stedet for at genbruge koden fra Astrid, at bruge den samme funktion til at lave en grid som bare er en store af keys og så gemme values et andet sted. Hvis det er implementeret med et linært index
# eller som en dictionary på cartisianIndex ville man nok kunne fjerne keys som er cell'er der ikke længere har safe elementer. 
# Så ville man stadig kunne bruge bounds funktioner til at finde de tætteste celler til constraints.
# Man kunne tjekke at når en cell ikke længere er safe så fjerner man de eventuelle constraints den har medført og bytter dem ud med cell da den nok er simplere