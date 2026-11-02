using LazySets, LinearAlgebra, ReachabilityAnalysis, Dates, JuMP, GLPK
include("Utilities.jl")

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
    #array
end

# Concrete types for the cell relations. Prefer these over Dict{..., AbstractArray{...}}, whose abstract value type makes
# every lookup type-unstable.
#   FlowRelation:   cell -> cells (successors of letting the system flow for one period, or predecessors when inverted)
#   ActionRelation: (action id, cell) -> cells
const CellIdxs{N} = Vector{CartesianIndex{N}}
const FlowRelation{N} = Dict{CartesianIndex{N},Vector{CartesianIndex{N}}}
const ActionRelation{N} = Dict{Tuple{Int64,CartesianIndex{N}},Vector{CartesianIndex{N}}}

# Number of cells per dimension as an NTuple{N,Int}. Use this when allocating per-cell arrays (falses(cellshape(grid)))
# instead of splatting grid.numCells, whose length is unknown to the compiler and makes the array dimension unknown.
cellshape(grid::Grid{N}) where N = size(grid.deadCells)::NTuple{N,Int}



function Grid(convexSet::LazySet, granularity::Float64)
    dimension = LazySets.dim(convexSet)

    lower_bounds = LazySets.low(convexSet)
    upper_bounds = LazySets.high(convexSet)

    numCells = zeros(Int, dimension)

    lower_bounds = Vector{Float64}(lower_bounds)
    upper_bounds = Vector{Float64}(upper_bounds)
    for (i, (lb, ub)) in enumerate(zip(lower_bounds, upper_bounds))
        numCells[i] = max(ceil((ub-lb)/granularity), 1)
        # A zero-width dimension (e.g. a constant like gravity) becomes a single cell centred on its value, so cell
        # centres (lower + granularity/2) equal the value and touching-cell index ranges are 1:1 rather than empty.
        if ub == lb
            lower_bounds[i] = lb - granularity/2
            upper_bounds[i] = ub + granularity/2
        end
    end

    # Possible to check if cells actually are within the convexSet, but for now we assume they are.
    #array = Array{Cell}(undef, (numCells...))  # Create an array to hold the grid cells

    deadCells = falses(numCells...)#BitArray(undef, (numCells...,))  # Create an array to label the dead cells falses(numCells...)#

    #for id in CartesianIndices(array)
    #    array[id] = Cell(id, [LinearIndices(array)[id]], [])
    #end

    return Grid(dimension, granularity, lower_bounds, upper_bounds, numCells, deadCells)#, array)
end

Base.show(io::IO, grid::Grid) = println(io,
    "Grid($(grid.granularity), $(grid.lower), $(grid.upper), $(grid.numCells))")

# Makes the grid iterable, returning each square in turn.
Base.length(grid::Grid) = length(grid.array)

Base.size(grid::Grid) = size(grid.array)

struct Cell
    id::CartesianIndex
    #sidx::Vector{Int}
    #uidx::Vector{Int}
    #pCells::Vector{CartesianIndex}
end

function Cell(id, sidx, uidx)
    return Cell(id)#, sidx, uidx, [])
end

Base.show(io::IO, cell::Cell) = println(io,
    "Cell($(cell.id), $(cell.sidx), $(cell.pCells))")



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

@inline car2vec(x::CartesianIndex) = collect(Tuple(x))

@inline car2vecFloat(x::CartesianIndex) = collect(Float64, Tuple(x))

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
    for idx in CartesianIndices(grid.deadCells)
        offset = car2vec(idx)
        center = (offset .- 0.5) * grid.granularity .+ grid.lower
        zonotopeArray[idx] = Zonotope(center, generators)

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

function get_cell_idx_bounds(grid::Grid, idx::CartesianIndex)
    lb = grid.lower .+ (Tuple(idx) .- 1) .* grid.granularity
    ub = lb .+ grid.granularity


    return lb, ub
end

function get_touching_cells(grid::Grid, convexSet::LazySet)
    touching_cells = []

    lower_bounds, upper_bounds = clamp.(LazySets.low(convexSet), grid.lower, grid.upper), clamp.(LazySets.high(convexSet), grid.lower, grid.upper)
    #lower_bounds = LazySets.low(convexSet)
    #lower_bounds = Int.(floor.(abs.(max.(lower_bounds, grid.lower) .- grid.lower) ./ grid.granularity) .+ 1)
    lower_bounds = Int.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity) .+ 1)

    upper_bounds = min.(grid.numCells, Int.(ceil.(abs.(upper_bounds .- grid.lower) ./ grid.granularity))) # min guards against float round-off at the top edge #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity
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

# TODO - This can be optimized significantly by utilizing the fact that Zonotopes are centrally symmetric. First halve the bounding box, then sweep from both sides.
function get_touching_cell_idxs(grid::Grid, convexSet::LazySet)
    touching_cell_idxs = []

    lower_bounds, upper_bounds = clamp.(LazySets.low(convexSet), grid.lower, grid.upper), clamp.(LazySets.high(convexSet), grid.lower, grid.upper)
    #@show lower_bounds, upper_bounds
    lower_bounds = Int.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity) .+ 1)

    upper_bounds = min.(grid.numCells, Int.(ceil.(abs.(upper_bounds .- grid.lower) ./ grid.granularity))) # min guards against float round-off at the top edge #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity
    #@show lower_bounds, upper_bounds
    #upper_bounds = Int.(floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity)
    ranges = [lower_bounds[i]:max(upper_bounds[i], lower_bounds[i]) for i in 1:grid.dimension]
    #@show ranges
    idxs = CartesianIndices((ranges...,))
    @show idxs
    hbox = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, grid.dimension))
    for idx in idxs
        #cell = grid.array[idx]
        #lower_bounds, upper_bounds = get_cell_idx_bounds(grid, idx)
        #cell_box = Hyperrectangle((lower_bounds + upper_bounds) / 2, (upper_bounds - lower_bounds) / 2)
        of = (car2vec(idx) .- 1) .* granularity
        #@show of
        LazySets.API.translate!(hbox, of)
        if !isdisjoint(convexSet, hbox)
            push!(touching_cell_idxs, idx)
        end
        LazySets.API.translate!(hbox, -of)
    end

    return touching_cell_idxs
end

function get_touching_cell_idxs_t(grid::Grid, convexSet::LazySet)
    touching_cell_idxs = []

    lower_bounds, upper_bounds = clamp.(LazySets.low(convexSet), grid.lower, grid.upper), clamp.(LazySets.high(convexSet), grid.lower, grid.upper)
    lower_bounds = Int.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity) .+ 1)

    upper_bounds = min.(grid.numCells, Int.(ceil.(abs.(upper_bounds .- grid.lower) ./ grid.granularity))) # min guards against float round-off at the top edge #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity
    #ranges = [lower_bounds[i]:max(upper_bounds[i], 1) for i in 1:grid.dimension]

    #idxs = CartesianIndices(ntuple(i -> lower_bounds[i]:max(upper_bounds[i], 1), grid.dimension))
    #@show idxs
    hbox = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, grid.dimension))
    for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:max(upper_bounds[i], 1), grid.dimension))
        #cell = grid.array[idx]
        #lower_bounds, upper_bounds = get_cell_idx_bounds(grid, idx)
        #cell_box = Hyperrectangle((lower_bounds + upper_bounds) / 2, (upper_bounds - lower_bounds) / 2)
        of = (car2vec(idx) .- 1) .* granularity
        #@show of
        LazySets.API.translate!(hbox, of)
        if !isdisjoint(convexSet, hbox)
            push!(touching_cell_idxs, idx)
        end
        LazySets.API.translate!(hbox, -of)
    end

    return touching_cell_idxs
end

# TODO - Right now uses unique to remove duplicate indices in the case that splitting_dim has an odd number of elements.   
# Also we shift the box back and forth every time, when we could just propagate along the sweeping axis... 
# Also i think the sweeping principle can be extended to the projection of the zonotope on dimensions that dont get sweeped...
# Bug: returns empty if zonotope has degenerate dimension.
function get_touching_cell_idxs_b(grid::Grid, convexSet::Zonotope)#, hbox::Hyperrectangle)
    touching_cell_idxs = []

    lower_bounds, upper_bounds = clamp.(LazySets.low(convexSet), grid.lower, grid.upper), clamp.(LazySets.high(convexSet), grid.lower, grid.upper)
    lower_bounds = Int.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity) .+ 1)

    upper_bounds = min.(grid.numCells, Int.(ceil.(abs.(upper_bounds .- grid.lower) ./ grid.granularity))) # min guards against float round-off at the top edge #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity

    splitting_dim = argmax(upper_bounds[i]-lower_bounds[i] for i in 1:grid.dimension)
    hbox = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, grid.dimension))
    if upper_bounds[splitting_dim] - lower_bounds[splitting_dim] > 1
        #println("gets used : )")
        cartesian_max = CartesianIndex(upper_bounds...)
        cartesian_min = CartesianIndex(lower_bounds...)

        upper_bounds[splitting_dim] = lower_bounds[splitting_dim] + cld(upper_bounds[splitting_dim] - lower_bounds[splitting_dim], 2) #cld(upper_bounds[splitting_dim], 2)

        sweeping_dim = argmax(upper_bounds[i]-lower_bounds[i] for i in 1:grid.dimension)
        cartesian_offset = CartesianIndex(ntuple(i -> i == sweeping_dim ? 1 : 0, grid.dimension))

        sweeping_range = range(lower_bounds[sweeping_dim], upper_bounds[sweeping_dim])
        rev_sweeping_range = range(upper_bounds[sweeping_dim], lower_bounds[sweeping_dim]; step=-1)
        #c = vec([i -> i == sweeping_dim ? Float64(lower_bounds[i] + upper_bounds[i]) / 2 : 0.0 for i in grid.dimension])
        #r = vec([i -> i == sweeping_dim ? Float64(upper_bounds[i] - lower_bounds[i]) / 2 : 0.0 for i in grid.dimension])
        #=
        c = grid.lower .+ (granularity/2)
        c[sweeping_dim] += ((lower_bounds[sweeping_dim]-1) * granularity)
        r = fill(granularity/2, grid.dimension)#grid.lower .+ (1.5 * granularity)
        r[sweeping_dim] = ((upper_bounds[sweeping_dim] - lower_bounds[sweeping_dim]) * granularity)
        sbox = Hyperrectangle(c, r)
        =#
        for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], lower_bounds[i])), grid.dimension))
            lower_idx = 0
            upper_idx = 0




            for off in sweeping_range
                of = (car2vec((idx + (cartesian_offset * off))) .- 1) .* granularity

                LazySets.API.translate!(hbox, of)
                if !isdisjoint(convexSet, hbox)
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

                of = (car2vec(idx + (cartesian_offset * off)) .- 1) .* granularity


                LazySets.API.translate!(hbox, of)
                if !isdisjoint(convexSet, hbox)
                    LazySets.API.translate!(hbox, -of)

                    upper_idx = off
                    break

                end
                LazySets.API.translate!(hbox, -of)
            end
            if upper_idx == 0
                continue
            end
            #=
            sof = (car2vec(idx)) .* granularity
            LazySets.API.translate!(sbox, sof)
            if isdisjoint(convexSet, sbox) && !((upper_idx == 0) && lower_idx == 0)
                @show (idx, sbox, sweeping_dim, lower_bounds, upper_bounds, upper_idx, lower_idx)
                LazySets.API.translate!(sbox, -sof)

                @show (idx, sbox, sweeping_dim, lower_bounds, upper_bounds, upper_idx, lower_idx)
                throw(Exception("fukcing shit"))
                continue
                #println("shit happens")
                #LazySets.API.translate!(sbox, -(car2vec(idx + cartesian_offset) .- 1) .* granularity)
            end
            
            if isdisjoint(convexSet, sbox) && (upper_idx == 0) && lower_idx == 0
                LazySets.API.translate!(sbox, -(car2vec(idx + cartesian_offset) .- 1) .* granularity)
                println("correct!")
                continue
                #@show (idx, (car2vec(idx) .- 1) .* granularity)
            elseif (upper_idx == 0) || lower_idx == 0
                println("not correct..")
                LazySets.API.translate!(sbox, -(car2vec(idx + cartesian_offset) .- 1) .* granularity)

                continue
            end
            LazySets.API.translate!(sbox, -(car2vec(idx + cartesian_offset) .- 1) .* granularity)
            =#

            for elem in lower_idx:upper_idx
                symmetric_id = cartesian_max - ((idx + (cartesian_offset * elem)) - cartesian_min)
                push!(touching_cell_idxs, idx + (cartesian_offset * elem), symmetric_id)
            end

        end
    else
        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:max(upper_bounds[i], lower_bounds[i]), grid.dimension))
            of = (car2vec(idx) .- 1) .* granularity
            LazySets.API.translate!(hbox, of)
            if !isdisjoint(convexSet, hbox)
                push!(touching_cell_idxs, idx)
            end
            LazySets.API.translate!(hbox, -of)
        end

    end
    return unique(touching_cell_idxs)
end

function get_touching_cell_idxs_l(cgrid::Grid{N}, convexSet::LazySet) where N
    touching_cell_idxs = CartesianIndex{N}[]

    lower_bounds, upper_bounds = clamp.(LazySets.low(convexSet), cgrid.lower, cgrid.upper), clamp.(LazySets.high(convexSet), cgrid.lower, cgrid.upper)
    #lower_bounds = Int.(floor.(abs.(lower_bounds .- cgrid.lower) ./ cgrid.granularity) .+ 1)
    lower_bounds = Int.(max.(floor.(abs.(lower_bounds .- cgrid.lower) ./ cgrid.granularity), 1))


    upper_bounds = min.(cgrid.numCells, Int.(ceil.(abs.(upper_bounds .- cgrid.lower) ./ cgrid.granularity))) # min guards against float round-off at the top edge #floor.(min.(upper_bounds, cgrid.upper) .- cgrid.lower) ./ cgrid.granularity

    sweeping_dim = argmax(upper_bounds[i]-lower_bounds[i] for i in 1:cgrid.dimension)
    hbox = Hyperrectangle(cgrid.lower .+ (cgrid.granularity/2), fill(cgrid.granularity/2, N))
    #@show upper_bounds, lower_bounds
    #@show hbox, sweeping_dim
    if upper_bounds[sweeping_dim] - lower_bounds[sweeping_dim] > 1
        #println("f'ing hope so")
        #upper_bounds[splitting_dim] = lower_bounds[splitting_dim] + cld(upper_bounds[splitting_dim] - lower_bounds[splitting_dim], 2) #cld(upper_bounds[splitting_dim], 2)


        cartesian_offset = CartesianIndex(ntuple(i -> i == sweeping_dim ? 1 : 0, N))

        sweeping_range = range(lower_bounds[sweeping_dim], upper_bounds[sweeping_dim])
        rev_sweeping_range = range(upper_bounds[sweeping_dim], lower_bounds[sweeping_dim]; step=-1)
        for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], lower_bounds[i])), N))
            lower_idx = 0
            upper_idx = 0
            for off in sweeping_range
                of = (car2vec((idx + (cartesian_offset * off))) .- 1) .* cgrid.granularity

                LazySets.API.translate!(hbox, of)
                if !isdisjoint(convexSet, hbox)
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

                of = (car2vec(idx + (cartesian_offset * off)) .- 1) .* cgrid.granularity


                LazySets.API.translate!(hbox, of)
                if !isdisjoint(convexSet, hbox)
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
        #idxs = CartesianIndices(ntuple(i -> lower_bounds[i]:max(upper_bounds[i], 1), grid.dimension))
        println("second")

        #for idx in CartesianIndices(ntuple(i -> lower_bounds[i] > upper_bounds[i] ? (min(lower_bounds[i], grid.numCells[i]):upper_bounds[i]) : lower_bounds[i]:min(upper_bounds[i], grid.numCells[i]), grid.dimension))
        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:upper_bounds[i], N))
            of = (car2vec(idx) .- 1) .* cgrid.granularity
            LazySets.API.translate!(hbox, of)
            if !isdisjoint(convexSet, hbox)
                push!(touching_cell_idxs, idx)
            end
            LazySets.API.translate!(hbox, -of)
        end

    end
    return unique(touching_cell_idxs)
end


function setup_model_constant_zonotope(Z::Zonotope)
    c = Z.center
    G = genmat(Z)
    N = eltype(c)
    m, n = size(G)
    lbounds = -one(N)
    ubounds = one(N)

    objective = zeros(N, n)

    model = direct_model(JuMP.optimizer_with_attributes(() -> GLPK.Optimizer(; method=GLPK.SIMPLEX, want_infeasibility_certificates=false)))
    #model = direct_model(JuMP.optimizer_with_attributes(() -> HiGHS.Optimizer()))
    set_silent(model)
    set_attribute(model, "msg_lev", GLPK.GLP_MSG_OFF)
    set_attribute(model, "presolve", GLPK.GLP_OFF)
    #@show typeof(model)
    JuMP.@variable(model, x[1:n], lower_bound=lbounds, upper_bound=ubounds)
    con_vector = JuMP.@constraint(model, con_ref, G * x .== c)
    #@show model
    #@objective(model, Min, 0.0)  # No objective function, just feasibility
    return model, con_vector
end

function change_constraints!(con_ref, rhs)
    #@show typeof(con_ref)
    #MOI.set(model, MOI.ConstraintSet(), con_ref, rhs)
    set_normalized_rhs.(con_ref, rhs)
end

function LinProg(x::Vector{Float64}, c::Vector{Float64}, model, con_ref)
    change_constraints!(con_ref, x - c)
    optimize!(model)
    return termination_status(model) == MOI.OPTIMAL
end

function get_touching_cell_idxs_test(grid::Grid{N}, convexSet::Zonotope{Float64,Vector{Float64},Matrix{Float64}}, model, con_ref) where N
    touching_cell_idxs::Vector{CartesianIndex{N}} = CartesianIndex{N}[]

    #lower_bounds = clamp.(LazySets.low(convexSet), grid.lower, grid.upper)
    #upper_bounds = clamp.(LazySets.high(convexSet), grid.lower, grid.upper)
    #lower_bounds = Int.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity) .+ 1)
    #lower_bounds = Int.(max.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity), 1))
    lower_bounds::Vector{Int64} = max.(floor.(abs.(clamp.(LazySets.low(convexSet), grid.lower, grid.upper) .- grid.lower) ./ grid.granularity), 1)


    upper_bounds::Vector{Int64} = min.(grid.numCells, Int.(ceil.(abs.(clamp.(LazySets.high(convexSet), grid.lower, grid.upper) .- grid.lower) ./ grid.granularity))) # min guards against float round-off at the top edge #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity
    #upper_bounds = Int.(ceil.(abs.(upper_bounds .- grid.lower) ./ grid.granularity)) #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity

    sweeping_dim::Int64 = argmax(upper_bounds[i]-lower_bounds[i] for i in 1:grid.dimension)
    next_sweep = argmax(i == sweeping_dim ? 0 : upper_bounds[i]-lower_bounds[i] for i in 1:grid.dimension)
    #hbox::LazySets.HyperrectangleModule.Hyperrectangle{Float64,Vector{Float64},Vector{Float64}} = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, N))
    #hbox = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, N))
    #hbox = Zonotope(grid.lower .+ (grid.granularity/2), diagm(fill(grid.granularity/2, N)))
    #convexSet = minkowski_sum(convexSet, Zonotope(zeros(grid.dimension), diagm(fill(granularity/2, N))))
    of::Vector{Float64} = grid.lower .+ (grid.granularity/2)
    tup1 = ntuple(i -> 1, grid.dimension)

    if upper_bounds[next_sweep] == lower_bounds[next_sweep]
        #throw(error("shouyld not happen!"))
        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:1:upper_bounds[i], N))

            push!(touching_cell_idxs, idx)


            #LazySets.API.translate!(hbox, -of)
        end
    else
        if upper_bounds[sweeping_dim] - lower_bounds[sweeping_dim] > 1
            #println("f'ing hope so")
            #upper_bounds[splitting_dim] = lower_bounds[splitting_dim] + cld(upper_bounds[splitting_dim] - lower_bounds[splitting_dim], 2) #cld(upper_bounds[splitting_dim], 2)


            cartesian_offset::CartesianIndex{N} = CartesianIndex(ntuple(i -> i == sweeping_dim ? 1 : 0, N))

            sweeping_range::StepRange{Int64,Int64} = lower_bounds[sweeping_dim]:1:upper_bounds[sweeping_dim]
            rev_sweeping_range::StepRange{Int64,Int64} = upper_bounds[sweeping_dim]:-1:lower_bounds[sweeping_dim]
            for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], lower_bounds[i])), N))

                #for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], lower_bounds[i])), N))
                lower_idx::Int64 = 0
                upper_idx::Int64 = 0
                for off in sweeping_range
                    #copyto!(hbox.center, (car2vecFloat((idx + (cartesian_offset * off))) .- 1.0) .* granularity + of)
                    #copyto!(hbox.center, (Tuple(idx + (cartesian_offset * off)) .- tup1) .* grid.granularity .+ of)
                    #tof = (car2vec((idx + (cartesian_offset * off))) .- 1) .* granularity
                    #tof = grid.lower .+ (car2vec((idx + (cartesian_offset * off))) .- 0.5) .* granularity
                    #LazySets.API.translate!(hbox, of)
                    #@show eltype(convexSet), eltype(hbox), promote_type(eltype(convexSet), eltype(hbox))
                    #throw(error("shshshs"))
                    #if !isdisjoint(convexSet, hbox)
                    #if !isdisjoint(convexSet, hbox, false; solver=model)
                    #if in(tof, convexSet)
                    #if LazySets.ZonotopeModule.in((Tuple(idx + (cartesian_offset * off)) .- tup1) .* grid.granularity .+ of, convexSet)
                    if LinProg((Tuple(idx + (cartesian_offset * off)) .- tup1) .* grid.granularity .+ of, convexSet.center, model, con_ref)
                        #LazySets.API.translate!(hbox, -of)
                        #@show typeof(tof), typeof(convexSet)
                        #throw(error("shshshs"))

                        lower_idx = off
                        break

                    end
                    #LazySets.API.translate!(hbox, -of)
                end
                if lower_idx == 0
                    continue
                end
                for off in rev_sweeping_range

                    #copyto!(hbox.center, (car2vecFloat(idx + (cartesian_offset * off)) .- 1.0) .* granularity + of)
                    #copyto!(hbox.center, (Tuple(idx + (cartesian_offset * off)) .- tup1) .* grid.granularity .+ of)

                    #tof = grid.lower .+ (car2vec(idx + (cartesian_offset * off)) .- 0.5) .* granularity
                    #tof = grid.lower .+ (car2vec(idx) .- 0.5) .* granularity

                    #LazySets.API.translate!(hbox, of)
                    #if !isdisjoint(convexSet, hbox, false; solver=model)
                    #if LazySets.ZonotopeModule.in((Tuple(idx + (cartesian_offset * off)) .- tup1) .* grid.granularity .+ of, convexSet)
                    if LinProg((Tuple(idx + (cartesian_offset * off)) .- tup1) .* grid.granularity .+ of, convexSet.center, model, con_ref)

                        #LazySets.API.translate!(hbox, -of)

                        upper_idx = off
                        break

                    end
                    #LazySets.API.translate!(hbox, -of)
                end
                if upper_idx == 0
                    continue
                end

                for elem in lower_idx:upper_idx

                    push!(touching_cell_idxs, idx + (cartesian_offset .* elem))
                end

            end
        else
            #idxs = CartesianIndices(ntuple(i -> lower_bounds[i]:max(upper_bounds[i], 1), N))


            #for idx in CartesianIndices(ntuple(i -> lower_bounds[i] > upper_bounds[i] ? (min(lower_bounds[i], grid.numCells[i]):upper_bounds[i]) : lower_bounds[i]:min(upper_bounds[i], grid.numCells[i]), N))
            for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:1:upper_bounds[i], N))
                #of::Vector{Float64} = (car2vec(idx) .- 1) .* granularity
                #copyto!(hbox.center, (Tuple(idx) .- tup1) .* grid.granularity .+ of)
                #tof = grid.lower .+ (car2vec(idx) .- 0.5) .* grid.granularity
                #LazySets.API.translate!(hbox, of)
                #if !isdisjoint(convexSet, hbox, false; solver=model)
                if LinProg(grid.lower .+ (car2vec(idx) .- 0.5) .* grid.granularity, convexSet.center, model, con_ref)

                    #if LazySets.ZonotopeModule.in(grid.lower .+ (car2vec(idx) .- 0.5) .* granularity, convexSet)
                    push!(touching_cell_idxs, idx)
                end
                #LazySets.API.translate!(hbox, -of)
            end

        end
    end
    return touching_cell_idxs
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
# pair of facets +-n' y <= sum_j |n' g_j|. Parallel generators give duplicate facets, which are harmless.
function ConstantZonotopeCover(G::AbstractMatrix{Float64}, radius::NTuple{N,Float64}, ::Val{N}; tol::Float64=1e-7) where N
    m = size(G, 2)
    normals = Vector{Float64}()
    offsets = Vector{Float64}()
    n = zeros(N)
    minor = zeros(N - 1, N - 1)
    subset = collect(1:N-1)
    while true
        if N == 2
            n[1] = G[2, subset[1]]
            n[2] = -G[1, subset[1]]
        else
            # n_k = (-1)^(k+1) det(G[rows != k, subset])
            for k in 1:N
                for (r, row) in enumerate(Iterators.filter(!=(k), 1:N)), (c, col) in enumerate(subset)
                    minor[r, c] = G[row, col]
                end
                n[k] = (isodd(k) ? 1.0 : -1.0) * det(minor)
            end
        end
        len = norm(n)
        if len > 1e-12
            n ./= len
            b = 0.0
            for j in 1:m
                b += abs(dot(n, view(G, :, j)))
            end
            append!(normals, n)
            push!(offsets, b)
            append!(normals, -n)
            push!(offsets, b)
        end
        # Next (N-1)-combination of 1:m in lexicographic order.
        i = N - 1
        while i >= 1 && subset[i] == m - (N - 1) + i
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
    return ConstantZonotopeCover{N}(reshape(normals, N, :), offsets, radius, tol)
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

# Centre of the reach set of fast cell idx, max_flow * (cell centre) + offset, written into out without allocating.
function fast_cell_center!(out::Vector{Float64}, max_flow::AbstractMatrix{Float64}, offset::AbstractVector{Float64}, grid::Grid{N}, idx::CartesianIndex{N}) where N
    for i in 1:N
        s = offset[i]
        for j in 1:N
            s += max_flow[i, j] * ((grid.lower[j] + grid.granularity/2) + (idx[j] - 1) * grid.granularity)
        end
        out[i] = s
    end
    return out
end

# Same cells as get_touching_cell_idxs_test(grid, Z, model, con_ref) for the translated zonotope Z with center c, using
# the precomputed H-representation instead of an LP per candidate cell. A cell touches the reach set iff its centre
# lies in the reach set bloated by half a cell, which is what cov represents.
function get_touching_cell_idxs_cover(grid::Grid{N}, cov::ConstantZonotopeCover{N}, c::AbstractVector{Float64}) where N
    touching_cell_idxs = CartesianIndex{N}[]
    lower_bounds = ntuple(i -> max(floor(Int, abs(clamp(c[i] - cov.radius[i], grid.lower[i], grid.upper[i]) - grid.lower[i]) / grid.granularity), 1), N)
    upper_bounds = ntuple(i -> min(grid.numCells[i], ceil(Int, abs(clamp(c[i] + cov.radius[i], grid.lower[i], grid.upper[i]) - grid.lower[i]) / grid.granularity)), N)

    sweeping_dim = argmax(upper_bounds[i] - lower_bounds[i] for i in 1:N)
    next_sweep = argmax(i == sweeping_dim ? 0 : upper_bounds[i] - lower_bounds[i] for i in 1:N)

    if upper_bounds[next_sweep] == lower_bounds[next_sweep]
        append!(touching_cell_idxs, CartesianIndices(ntuple(i -> lower_bounds[i]:upper_bounds[i], N)))
    elseif upper_bounds[sweeping_dim] - lower_bounds[sweeping_dim] > 1
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

# Cells touching the zonotope convexSet, without LPs: a cell touches it iff the cell centre lies in convexSet bloated by
# half a cell, which is tested against the facets of that bloated zonotope (see ConstantZonotopeCover). Gives the same
# cells as get_touching_cell_idxs_lz_lp. Falls back to it when the facet count explodes (many generators in high
# dimensions). Without LPs this is also safe to call from several threads (GLPK is not).
function get_touching_cell_idxs_lz(grid::Grid{N}, convexSet::Zonotope{Float64,Vector{Float64},Matrix{Float64}}) where N
    G = genmat(convexSet)
    nonzero = [j for j in axes(G, 2) if any(!iszero, view(G, :, j))]
    if cover_facet_pairs(length(nonzero) + N, N) > MAX_COVER_FACET_PAIRS
        return get_touching_cell_idxs_lz_lp(grid, convexSet)
    end
    # Bound the candidate cells by the set itself (as the LP version does); test membership against the bloated set.
    radius = ntuple(i -> sum(abs, view(G, i, :)), N)
    bloated = hcat(G[:, nonzero], Matrix((grid.granularity / 2) * I, N, N))
    return get_touching_cell_idxs_cover(grid, ConstantZonotopeCover(bloated, radius, Val(N)), convexSet.center)
end

function get_touching_cell_idxs_lz_lp(grid::Grid{N}, convexSet::Zonotope{Float64,Vector{Float64},Matrix{Float64}}) where N
    touching_cell_idxs::Vector{CartesianIndex{N}} = CartesianIndex{N}[]

    #lower_bounds = clamp.(LazySets.low(convexSet), grid.lower, grid.upper)
    #upper_bounds = clamp.(LazySets.high(convexSet), grid.lower, grid.upper)
    #lower_bounds = Int.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity) .+ 1)
    #lower_bounds = Int.(max.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity), 1))
    lower_bounds::Vector{Int64} = max.(floor.(abs.(clamp.(LazySets.low(convexSet), grid.lower, grid.upper) .- grid.lower) ./ grid.granularity), 1)


    upper_bounds::Vector{Int64} = min.(grid.numCells, Int.(ceil.(abs.(clamp.(LazySets.high(convexSet), grid.lower, grid.upper) .- grid.lower) ./ grid.granularity))) # min guards against float round-off at the top edge #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity
    #upper_bounds = Int.(ceil.(abs.(upper_bounds .- grid.lower) ./ grid.granularity)) #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity

    sweeping_dim::Int64 = argmax(upper_bounds[i]-lower_bounds[i] for i in 1:grid.dimension)
    next_sweep = argmax(i == sweeping_dim ? 0 : upper_bounds[i]-lower_bounds[i] for i in 1:grid.dimension)

    #hbox::LazySets.HyperrectangleModule.Hyperrectangle{Float64,Vector{Float64},Vector{Float64}} = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, N))
    #hbox = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, N))
    hbox = Zonotope(grid.lower .+ (grid.granularity/2), diagm(fill(grid.granularity/2, N)))
    #convexSet = minkowski_sum(convexSet, Zonotope(zeros(grid.dimension), diagm(fill(granularity/2, N))))
    of::Vector{Float64} = grid.lower .+ (grid.granularity/2)
    if upper_bounds[next_sweep] == lower_bounds[next_sweep]
        #throw(error("shouyld not happen!"))
        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:1:upper_bounds[i], N))

            push!(touching_cell_idxs, idx)


            #LazySets.API.translate!(hbox, -of)
        end
    else
        if upper_bounds[sweeping_dim] - lower_bounds[sweeping_dim] > 1
            #println("f'ing hope so")
            #upper_bounds[splitting_dim] = lower_bounds[splitting_dim] + cld(upper_bounds[splitting_dim] - lower_bounds[splitting_dim], 2) #cld(upper_bounds[splitting_dim], 2)


            cartesian_offset::CartesianIndex{N} = CartesianIndex(ntuple(i -> i == sweeping_dim ? 1 : 0, N))

            sweeping_range::StepRange{Int64,Int64} = lower_bounds[sweeping_dim]:1:upper_bounds[sweeping_dim]
            rev_sweeping_range::StepRange{Int64,Int64} = upper_bounds[sweeping_dim]:-1:lower_bounds[sweeping_dim]
            for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], lower_bounds[i])), N))

                #for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], lower_bounds[i])), N))
                lower_idx::Int64 = 0
                upper_idx::Int64 = 0
                for off in sweeping_range
                    copyto!(hbox.center, (car2vecFloat((idx + (cartesian_offset * off))) .- 1.0) .* grid.granularity + of)
                    #copyto!(hbox.center, (Tuple(idx + (cartesian_offset * off)) .- 1) .* granularity .+ of)
                    #tof = (car2vec((idx + (cartesian_offset * off))) .- 1) .* granularity
                    #tof = grid.lower .+ (car2vec((idx + (cartesian_offset * off))) .- 0.5) .* granularity
                    #LazySets.API.translate!(hbox, of)
                    #@show eltype(convexSet), eltype(hbox), promote_type(eltype(convexSet), eltype(hbox))
                    #throw(error("shshshs"))
                    #if !isdisjoint(convexSet, hbox)
                    if !isdisjoint(convexSet, hbox, false)
                        #if in(tof, convexSet)
                        #if LazySets.ZonotopeModule.in(grid.lower .+ (car2vec((idx + (cartesian_offset * off))) .- 0.5) .* granularity, convexSet)
                        #LazySets.API.translate!(hbox, -of)
                        #@show typeof(tof), typeof(convexSet)
                        #throw(error("shshshs"))

                        lower_idx = off
                        break

                    end
                    #LazySets.API.translate!(hbox, -of)
                end
                if lower_idx == 0
                    continue
                end
                for off in rev_sweeping_range

                    copyto!(hbox.center, (car2vecFloat(idx + (cartesian_offset * off)) .- 1.0) .* grid.granularity + of)
                    #tof = grid.lower .+ (car2vec(idx + (cartesian_offset * off)) .- 0.5) .* granularity
                    #tof = grid.lower .+ (car2vec(idx) .- 0.5) .* granularity

                    #LazySets.API.translate!(hbox, of)
                    if !isdisjoint(convexSet, hbox, false)
                        #if LazySets.ZonotopeModule.in(grid.lower .+ (car2vec(idx + (cartesian_offset * off)) .- 0.5) .* granularity, convexSet)

                        #LazySets.API.translate!(hbox, -of)

                        upper_idx = off
                        break

                    end
                    #LazySets.API.translate!(hbox, -of)
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
            #idxs = CartesianIndices(ntuple(i -> lower_bounds[i]:max(upper_bounds[i], 1), N))


            #for idx in CartesianIndices(ntuple(i -> lower_bounds[i] > upper_bounds[i] ? (min(lower_bounds[i], grid.numCells[i]):upper_bounds[i]) : lower_bounds[i]:min(upper_bounds[i], grid.numCells[i]), N))
            for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:1:upper_bounds[i], N))
                #of::Vector{Float64} = (car2vec(idx) .- 1) .* granularity
                copyto!(hbox.center, (car2vecFloat(idx) .- 1.0) .* grid.granularity + of)
                #tof = grid.lower .+ (car2vec(idx) .- 0.5) .* granularity
                #LazySets.API.translate!(hbox, of)
                if !isdisjoint(convexSet, hbox, false)
                    #if LazySets.ZonotopeModule.in(grid.lower .+ (car2vec(idx) .- 0.5) .* granularity, convexSet)
                    push!(touching_cell_idxs, idx)

                end
                #LazySets.API.translate!(hbox, -of)
            end


        end
    end
    return touching_cell_idxs
end

function get_touching_cell_idxs_lz_reuse(grid::Grid{N}, convexSet::Zonotope{Float64,Vector{Float64},Matrix{Float64}}, model) where N
    touching_cell_idxs::Vector{CartesianIndex{N}} = CartesianIndex{N}[]

    #lower_bounds = clamp.(LazySets.low(convexSet), grid.lower, grid.upper)
    #upper_bounds = clamp.(LazySets.high(convexSet), grid.lower, grid.upper)
    #lower_bounds = Int.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity) .+ 1)
    #lower_bounds = Int.(max.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity), 1))
    lower_bounds::Vector{Int64} = max.(floor.(abs.(clamp.(LazySets.low(convexSet), grid.lower, grid.upper) .- grid.lower) ./ grid.granularity), 1)


    upper_bounds::Vector{Int64} = min.(grid.numCells, Int.(ceil.(abs.(clamp.(LazySets.high(convexSet), grid.lower, grid.upper) .- grid.lower) ./ grid.granularity))) # min guards against float round-off at the top edge #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity
    #upper_bounds = Int.(ceil.(abs.(upper_bounds .- grid.lower) ./ grid.granularity)) #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity

    sweeping_dim::Int64 = argmax(upper_bounds[i]-lower_bounds[i] for i in 1:grid.dimension)
    #hbox::LazySets.HyperrectangleModule.Hyperrectangle{Float64,Vector{Float64},Vector{Float64}} = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, N))
    #hbox = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, N))
    hbox::Zonotope{Float64,Vector{Float64},Matrix{Float64}} = Zonotope(grid.lower .+ (granularity/2), diagm(fill(granularity/2, N)))
    #convexSet = minkowski_sum(convexSet, Zonotope(zeros(grid.dimension), diagm(fill(granularity/2, N))))
    if false #upper_bounds[sweeping_dim] - lower_bounds[sweeping_dim] > 1
        #println("f'ing hope so")
        #upper_bounds[splitting_dim] = lower_bounds[splitting_dim] + cld(upper_bounds[splitting_dim] - lower_bounds[splitting_dim], 2) #cld(upper_bounds[splitting_dim], 2)


        cartesian_offset::CartesianIndex{N} = CartesianIndex(ntuple(i -> i == sweeping_dim ? 1 : 0, N))

        sweeping_range::StepRange = lower_bounds[sweeping_dim]:1:upper_bounds[sweeping_dim]
        rev_sweeping_range::StepRange = upper_bounds[sweeping_dim]:-1:lower_bounds[sweeping_dim]
        for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], lower_bounds[i])), N))

            #for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], lower_bounds[i])), N))
            lower_idx::Int64 = 0
            upper_idx::Int64 = 0
            for off in sweeping_range
                #of::Vector{Float64} = (car2vec((idx + (cartesian_offset * off))) .- 1) .* granularity
                #tof = (car2vec((idx + (cartesian_offset * off))) .- 1) .* granularity
                tof::Vector{Float64} = grid.lower .+ (car2vec((idx + (cartesian_offset * off))) .- 0.5) .* granularity
                #LazySets.API.translate!(hbox, of)
                #if !isdisjoint(convexSet, hbox)
                #if !isdisjoint(convexSet, hbox, false; solver=model)
                #if in(tof, convexSet)
                if LazySets.API.in(tof, convexSet; solver=model)
                    #LazySets.API.translate!(hbox, -of)
                    #@show typeof(tof), typeof(convexSet)
                    #throw(error("shshshs"))

                    lower_idx = off
                    break

                end
                #LazySets.API.translate!(hbox, -of)
            end
            if lower_idx == 0
                continue
            end
            for off in rev_sweeping_range

                #of::Vector{Float64} = (car2vec(idx + (cartesian_offset * off)) .- 1) .* granularity
                #tof = grid.lower .+ (car2vec(idx + (cartesian_offset * off)) .- 0.5) .* granularity
                #tof = grid.lower .+ (car2vec(idx) .- 0.5) .* granularity

                #LazySets.API.translate!(hbox, of)
                #if !isdisjoint(convexSet, hbox, false; solver=model)
                if LazySets.API.in(grid.lower .+ (car2vec(idx + (cartesian_offset * off)) .- 0.5) .* granularity, convexSet; solver=model)

                    #LazySets.API.translate!(hbox, -of)

                    upper_idx = off
                    break

                end
                #LazySets.API.translate!(hbox, -of)
            end
            if upper_idx == 0
                continue
            end

            for elem in lower_idx:upper_idx

                push!(touching_cell_idxs, idx + (cartesian_offset * elem))
            end

        end
    else
        #idxs = CartesianIndices(ntuple(i -> lower_bounds[i]:max(upper_bounds[i], 1), N))


        #for idx in CartesianIndices(ntuple(i -> lower_bounds[i] > upper_bounds[i] ? (min(lower_bounds[i], grid.numCells[i]):upper_bounds[i]) : lower_bounds[i]:min(upper_bounds[i], grid.numCells[i]), N))
        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:1:upper_bounds[i], N))
            of::Vector{Float64} = (car2vec(idx) .- 1) .* granularity
            #tof = grid.lower .+ (car2vec(idx) .- 0.5) .* granularity
            LazySets.API.translate!(hbox, of)
            if !isdisjoint(convexSet, hbox, false; solver=model)
                #if LazySets.API.in(grid.lower .+ (car2vec(idx) .- 0.5) .* granularity, convexSet; solver=model)
                push!(touching_cell_idxs, idx)
            end
            LazySets.API.translate!(hbox, -of)
        end

    end
    return touching_cell_idxs
end

function get_touching_cell_idxs_bloated(grid::Grid{N}, convexSet::Zonotope{Float64,Vector{Float64},Matrix{Float64}}) where N
    touching_cell_idxs = CartesianIndex{N}[]

    #lower_bounds::Vector{Float64} = clamp.(LazySets.low(convexSet), grid.lower, grid.upper)
    #upper_bounds::Vector{Float64} = clamp.(LazySets.high(convexSet), grid.lower, grid.upper)
    #lower_bounds = Int.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity) .+ 1)
    #lower_bounds = Int.(max.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity), 1))
    lower_bounds::Vector{Int64} = Int.(max.(floor.(abs.(clamp.(LazySets.low(convexSet), grid.lower, grid.upper) .- grid.lower) ./ grid.granularity), 1))


    upper_bounds::Vector{Int64} = min.(grid.numCells, Int.(ceil.(abs.(clamp.(LazySets.high(convexSet), grid.lower, grid.upper) .- grid.lower) ./ grid.granularity))) # min guards against float round-off at the top edge #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity
    #upper_bounds = Int.(ceil.(abs.(upper_bounds .- grid.lower) ./ grid.granularity)) #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity

    sweeping_dim::Int64 = argmax(upper_bounds[i]-lower_bounds[i] for i in 1:grid.dimension)
    #hbox::LazySets.HyperrectangleModule.Hyperrectangle{Float64,Vector{Float64},Vector{Float64}} = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, N))
    #hbox = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, N))
    #zbox = Zonotope(zeros(grid.dimension), diagm(fill(granularity/2, N)))
    #convexSet = minkowski_sum(convexSet, Zonotope(zeros(grid.dimension), diagm(fill(granularity/2, N))))
    if upper_bounds[sweeping_dim] - lower_bounds[sweeping_dim] > 1
        #println("f'ing hope so")
        #upper_bounds[splitting_dim] = lower_bounds[splitting_dim] + cld(upper_bounds[splitting_dim] - lower_bounds[splitting_dim], 2) #cld(upper_bounds[splitting_dim], 2)


        cartesian_offset::CartesianIndex{N} = CartesianIndex(ntuple(i -> i == sweeping_dim ? 1 : 0, N))

        sweeping_range = range(lower_bounds[sweeping_dim], upper_bounds[sweeping_dim])
        rev_sweeping_range = range(upper_bounds[sweeping_dim], lower_bounds[sweeping_dim]; step=-1)
        for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], lower_bounds[i])), N))
            lower_idx::Int64 = 0
            upper_idx::Int64 = 0
            for off in sweeping_range
                #of = (car2vec((idx + (cartesian_offset * off))) .- 1) .* granularity
                #tof = (car2vec((idx + (cartesian_offset * off))) .- 1) .* granularity
                #tof = grid.lower .+ (car2vec((idx + (cartesian_offset * off))) .- 0.5) .* granularity
                #LazySets.API.translate!(hbox, of)
                #if !isdisjoint(convexSet, hbox)
                #if !isdisjoint(convexSet, hbox)
                #if in(tof, convexSet)
                if !isdisjoint(Singleton(grid.lower .+ (car2vec((idx + (cartesian_offset * off))) .- 0.5) .* granularity), convexSet)
                    #LazySets.API.translate!(hbox, -of)
                    #@show typeof(tof), typeof(convexSet)
                    #throw(error("shshshs"))

                    lower_idx = off
                    break

                end
                #LazySets.API.translate!(hbox, -of)
            end
            if lower_idx == 0
                continue
            end
            for off in rev_sweeping_range

                #of = (car2vec(idx + (cartesian_offset * off)) .- 1) .* granularity
                #tof = grid.lower .+ (car2vec(idx + (cartesian_offset * off)) .- 0.5) .* granularity
                #tof = grid.lower .+ (car2vec(idx) .- 0.5) .* granularity

                #LazySets.API.translate!(hbox, of)
                #if !isdisjoint(convexSet, hbox)
                if isdisjoint(Singleton(grid.lower .+ (car2vec(idx + (cartesian_offset * off)) .- 0.5) .* granularity), convexSet)

                    #LazySets.API.translate!(hbox, -of)

                    upper_idx = off
                    break

                end
                #LazySets.API.translate!(hbox, -of)
            end
            if upper_idx == 0
                continue
            end

            for elem in lower_idx:upper_idx

                push!(touching_cell_idxs, idx + (cartesian_offset * elem))
            end

        end
    else
        #idxs = CartesianIndices(ntuple(i -> lower_bounds[i]:max(upper_bounds[i], 1), N))


        #for idx in CartesianIndices(ntuple(i -> lower_bounds[i] > upper_bounds[i] ? (min(lower_bounds[i], grid.numCells[i]):upper_bounds[i]) : lower_bounds[i]:min(upper_bounds[i], grid.numCells[i]), N))
        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:upper_bounds[i], N))
            #of = (car2vec(idx) .- 1) .* granularity
            #tof = grid.lower .+ (car2vec(idx) .- 0.5) .* granularity
            #LazySets.API.translate!(hbox, of)
            #if !isdisjoint(convexSet, hbox)
            #if LazySets.ZonotopeModule.in(grid.lower .+ (car2vec(idx) .- 0.5) .* granularity, convexSet)
            if !isdisjoint(convexSet, Singleton(grid.lower .+ (car2vec(idx) .- 0.5) .* granularity))
                push!(touching_cell_idxs, idx)
            end
            #LazySets.API.translate!(hbox, -of)
        end

    end
    return unique(touching_cell_idxs)
end
#=
function get_touching_cell_idxs_static_zonotope(grid::Grid, convexSet::Zonotope{Float64,MVector{N,Float64},SMatrix{N,N,Float64,T}}) where {N,T}
    touching_cell_idxs = []

    lower_bounds, upper_bounds = clamp.(LazySets.low(convexSet), grid.lower, grid.upper), clamp.(LazySets.high(convexSet), grid.lower, grid.upper)
    #lower_bounds = Int.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity) .+ 1)
    lower_bounds = Int.(max.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity), 1))


    upper_bounds = min.(grid.numCells, Int.(ceil.(abs.(upper_bounds .- grid.lower) ./ grid.granularity))) # min guards against float round-off at the top edge #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity

    sweeping_dim = argmax(upper_bounds[i]-lower_bounds[i] for i in 1:grid.dimension)
    hbox::LazySets.HyperrectangleModule.Hyperrectangle{Float64,Vector{Float64},Vector{Float64}} = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, grid.dimension))
    #hbox = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, grid.dimension))
    if upper_bounds[sweeping_dim] - lower_bounds[sweeping_dim] > 1
        #println("f'ing hope so")
        #upper_bounds[splitting_dim] = lower_bounds[splitting_dim] + cld(upper_bounds[splitting_dim] - lower_bounds[splitting_dim], 2) #cld(upper_bounds[splitting_dim], 2)


        cartesian_offset = CartesianIndex(ntuple(i -> i == sweeping_dim ? 1 : 0, grid.dimension))

        sweeping_range = range(lower_bounds[sweeping_dim], upper_bounds[sweeping_dim])
        rev_sweeping_range = range(upper_bounds[sweeping_dim], lower_bounds[sweeping_dim]; step=-1)
        for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], lower_bounds[i])), grid.dimension))
            lower_idx = 0
            upper_idx = 0
            for off in sweeping_range
                of = (car2vec((idx + (cartesian_offset * off))) .- 1) .* granularity

                LazySets.API.translate!(hbox, of)
                if !isdisjoint(convexSet, hbox)
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

                of = (car2vec(idx + (cartesian_offset * off)) .- 1) .* granularity


                LazySets.API.translate!(hbox, of)
                if !isdisjoint(convexSet, hbox)
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
        #idxs = CartesianIndices(ntuple(i -> lower_bounds[i]:max(upper_bounds[i], 1), grid.dimension))


        #for idx in CartesianIndices(ntuple(i -> lower_bounds[i] > upper_bounds[i] ? (min(lower_bounds[i], grid.numCells[i]):upper_bounds[i]) : lower_bounds[i]:min(upper_bounds[i], grid.numCells[i]), grid.dimension))
        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:upper_bounds[i], grid.dimension))
            of = (car2vec(idx) .- 1) .* granularity
            LazySets.API.translate!(hbox, of)
            if !isdisjoint(convexSet, hbox)
                push!(touching_cell_idxs, idx)
            end
            LazySets.API.translate!(hbox, -of)
        end

    end
    return unique(touching_cell_idxs)
end
=#
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

function offsets(grid::Grid{N}) where N
    offsets = Set{CartesianIndex{N}}(CartesianIndices(ntuple(_ -> -1:1, N)))
    return delete!(offsets, zero(CartesianIndex{N}))
end

# Cells on the outer faces of the grid. Dimensions with a single cell are skipped, since every cell lies on both of
# their faces.
function border_cells(grid::Grid{N}) where N
    shape = cellshape(grid)
    res = CartesianIndex{N}[]
    for i in 1:N
        if shape[i] == 1
            continue
        end
        for side in unique((1, shape[i]))
            append!(res, CartesianIndices(ntuple(j -> j == i ? (side:side) : (1:shape[j]), N)))
        end
    end
    return unique!(res)
end

function grow_indices(grid::Grid, idxs)
    offset = offsets(grid)
    res = []
    for idx in idxs
        #@show idx, offset
        union!(res, collect(idx + CartesianIndex(of) for of in offset))
    end

    intersect!(res, union(Tuple(idx) for idx in CartesianIndices(grid.array)))
    setdiff!(res, collect(Tuple(id) for id in idxs))
    return res
end

function grow_indices(grid::Grid, idxs, offset)
    res = []
    for idx in idxs
        union!(res, collect(idx .+ of for of in offset))
    end

    intersect!(res, union(idx for idx in CartesianIndices(grid.deadCells)))
    setdiff!(res, collect(CartesianIndex(id) for id in idxs))
    #@show res
    return res
end

function grow_indices_t(grid::Grid, idxs::Vector{<:CartesianIndex{N}}, offsets::Set{<:CartesianIndex{N}}) where N
    R = CartesianIndices(grid.deadCells)
    src = Set(CartesianIndex.(idxs))          # fast membership test
    res = Set{eltype(R)}()
    for idx in src, off in offsets
        nb = idx + CartesianIndex(off)
        if nb in R && !(nb in src)
            push!(res, nb)
        end
    end
    return collect(res)                # drop sort! if order doesn't matter
end

function grow_indices_b(grid::Grid, idxs::AbstractVector{<:CartesianIndex{N}}, offsets::Set{<:CartesianIndex{N}}) where N
    mask = falses(size(grid.deadCells))
    src = CartesianIndex.(idxs)
    offs = CartesianIndex.(offsets)
    for idx in src, off in offs
        nb = idx + off
        checkbounds(Bool, mask, nb) && (mask[nb] = true)
    end
    mask[src] .= false                         # remove the originals
    return findall(mask)                       # sorted, Vector{CartesianIndex{N}}
end

# https://github.com/AstridHornBrorholt/Shielded-Learning-for-Hybrid-Systems/blob/22c9fc220ef40d55877ff1360be7286a7a506620/Shared%20Code/ShieldSynthesis.jl#L88
# Looks to terminate early..
function make_shield(grid::Grid, action_set::ActionRelation{N}, no_action_set::FlowRelation{N}, max_steps::Int, Act::Vector{Action}) where N
    if max_steps == 0
        #@show findall(p -> isempty(p), no_action_set)
        grid.deadCells[findall(p -> isempty(p), no_action_set)] .= true

        return (grid, max_steps, action_set, no_action_set, falses(cellshape(grid)))
    end

    i = max_steps
    println("$(max_steps - i) : $(Dates.format(now(), "HH:MM:SS"))")
    dims = grid.dimension
    dead = grid.deadCells
    can_act_matrix = trues(cellshape(grid))
    no_action_bad_matrix = trues(cellshape(grid))
    dead´ = nothing
    action_set´ = nothing
    no_action_set´ = nothing
    a_tombstones = 0
    n_tombstones = 0

    grid.deadCells[findall(p -> isempty(p), no_action_set)] .= true
    no_action_bad_matrix[findall(p -> !isempty(p), no_action_set)] .= false
    #=no_action_bad_matrix[findall(p -> isempty(p), no_action_set)] .= true
    for i in 290:300
        for j in 148:150
            tidx = CartesianIndex(i, j)
            @show grid.deadCells[tidx], no_action_bad_matrix[tidx]
        end
    end=#
    filter!(p -> !isempty(p.second), no_action_set)
    filter!(p -> !isempty(p.second), action_set)
    action_set::ActionRelation{N} = ActionRelation{N}(action_set)
    no_action_set::FlowRelation{N} = FlowRelation{N}(no_action_set)

    #@show length(keys(action_set))
    #@show length(keys(no_action_set))

    while i > 0
        #grid´, action_set´ = shield_step!(grid, action_set, no_action_set, Act)
        #dead´, la, ln, ld, action_set´, no_action_set´ = shield_step!(dead, can_act_matrix, no_action_bad_matrix, action_set, no_action_set, Act, dims)
        dead´, la, ln, ld = shield_step!(dead, can_act_matrix, no_action_bad_matrix, action_set, no_action_set, Act, dims)
        if i % 10 == 0
            println("$(max_steps - i) : $(Dates.format(now(), "HH:MM:SS"))")
        end
        #@show can_act_matrix == trues(grid.numCells...)
        if la == 0 && ln == 0 && ld == 0

            println("Fixed point found at $(max_steps-i) steps!")
            dead = dead´
            #action_set = action_set´
            #no_action_set = no_action_set´
            #@show shield_step!(dead, can_act_matrix, no_action_bad_matrix, action_set, no_action_set, Act, dims)
            break
        end
        a_tombstones += la
        n_tombstones += ln

        dead = dead´
        #action_set = action_set´
        #no_action_set = no_action_set´

        if a_tombstones + n_tombstones > 50
            #@show (a_tombstones, n_tombstones)
            #@show length(keys(action_set))
            #filter!(k -> !in(k.first, a_pops), action_set)
            #@show length(keys(action_set))
            #filter!(k -> !in(k.first, n_pops), no_action_set)
            action_set = ActionRelation{N}(action_set)
            no_action_set = FlowRelation{N}(no_action_set)
            a_tombstones = 0
            n_tombstones = 0

        end

        #@show action_set == action_set´
        i -= 1
    end
    #@show length(keys(action_set))
    #@show length(keys(no_action_set))

    #@show can_act_matrix == trues(grid.numCells...)
    #@show no_action_bad_matrix == falses(grid.numCells...)

    return (grid, max_steps - i, action_set, no_action_set, no_action_bad_matrix)
end

function shield_step!(deadCells::BitArray, can_act_matrix, no_action_bad_matrix, action_set::ActionRelation{N}, no_action_set::FlowRelation{N}, Act::Vector{Action}, dims::Int64) where N
    #deadCells´ = copy(deadCells)
    act_pop_keys = Tuple{Int64,CartesianIndex{N}}[]
    no_act_pop_keys = CartesianIndex{N}[]
    new_dead_cells = CartesianIndex{N}[]
    for idx in CartesianIndices(deadCells)
        if !deadCells[idx]
            #no_action_bad = any(i -> i == 0, collect(grid.array[nc].sidx[1] for nc in cell.pCells))
            #no_action_bad = any(grid.deadCells[idx] for idx in cell.pCells)
            no_action_bad = no_action_bad_matrix[idx]
            if !no_action_bad && haskey(no_action_set, idx)
                #@show collect(deadCells[idxx] for idxx in no_action_set[idx])
                #if any(deadCells[idxx] for idxx in no_action_set[idx])
                if any(deadCells[no_action_set[idx]])# for idxx in no_action_set[idx])
                    #println("has key")
                    no_action_bad = true
                    no_action_bad_matrix[idx] = true
                    push!(no_act_pop_keys, idx)
                elseif isempty(no_action_set[idx])
                    no_action_bad = true
                    no_action_bad_matrix[idx] = true
                    push!(no_act_pop_keys, idx)
                end
            end

            #can_act = can_act_matrix[idx]
            approved_action_count = 0
            if can_act_matrix[idx]
                #should_act = false
                for act in Act
                    # maybe we should check if the action is even in the action_set first, but this is probably faster than checking for every action. 
                    if haskey(action_set, (act.id, idx))
                        approved_action_count += 1
                        should_act = true
                        #if isempty(collect(deadCells[idxx] for idxx in action_set[(act.id, idx)]))
                        if isempty(collect(deadCells[action_set[(act.id, idx)]]))# for idxx in action_set[(act.id, idx)]))
                            #@show action_set[(act.id, cell.id)]
                            push!(act_pop_keys, (act.id, idx))
                            approved_action_count -= 1
                            #end

                            #elseif any(deadCells[idxx] for idxx in action_set[(act.id, idx)])
                        elseif any(deadCells[action_set[(act.id, idx)]])# for idxx in action_set[(act.id, idx)])
                            push!(act_pop_keys, (act.id, idx))
                            approved_action_count -= 1
                            #elseif any(no_action_bad_matrix[idxx] for idxx in action_set[(act.id, idx)])
                        elseif any(no_action_bad_matrix[action_set[(act.id, idx)]])# for idxx in action_set[(act.id, idx)])
                            push!(act_pop_keys, (act.id, idx))
                            approved_action_count -= 1
                        end
                    end
                end
                #if !should_act
                #can_act_matrix[idx] = false
                #push!(act_pop_keys, (0, idx))
                #end
            end
            #@show can_act
            #@show can_act, isempty(cell.pCells)
            if approved_action_count <= 0 && no_action_bad #isempty(cell.pCells) && cell.sidx[1] == 1 #no_action_bad && can_act == 0
                #grid´.array[cell.id].sidx[1] = 1
                #deadCells´[idx] = true
                #push!(no_act_pop_keys, cell.id)
                #@show grid´.deadCells[cell.id] == grid.deadCells[cell.id]
                no_action_bad_matrix[idx] = true
                can_act_matrix[idx] = false
                push!(new_dead_cells, idx)
            end
            #=else
            for act in Act
                if haskey(action_set, (act, cell.id))

                    push!(act_pop_keys, (act, cell.id))

                end
                if haskey(no_action_set, cell.id)

                    push!(no_act_pop_keys, cell.id)

                end
            end=#
        end
    end
    #grid´.deadCells[new_dead_cells...] = true
    #@show grid.deadCells == grid´.deadCells
    #@show length(act_pop_keys)
    #@show length(no_act_pop_keys)
    filter!(k -> !in(k.first, act_pop_keys), action_set)
    filter!(k -> !in(k.first, no_act_pop_keys), no_action_set)
    if !isempty(new_dead_cells)
        deadCells[new_dead_cells] .= true
    end
    #@show length(action_set)

    #return deadCells, length(act_pop_keys), length(no_act_pop_keys), length(new_dead_cells), action_set, no_action_set #, action_set
    return deadCells, length(act_pop_keys), length(no_act_pop_keys), length(new_dead_cells) #, action_set
end

# True if any of the cells in idxs is set in m. Avoids allocating m[idxs].
@inline function any_at(m::AbstractArray{Bool}, idxs)
    for i in idxs
        m[i] && return true
    end
    return false
end

# Re-evaluates the cells in `indices` against the current dead/no_action_bad state.
# Only state changes are reported, so empty return lists mean a fixed point has been reached:
#   act_idx_pop_keys : cells that just lost one of their actions (once per lost action)
#   no_act_pop_keys  : cells where not acting just became unsafe (no_action_bad_matrix went false -> true)
#   new_dead_cells   : cells that just became dead
# popped_actions holds the (action id, cell) pairs already found unsafe. Pass the same set on every call, so each lost
# action is reported only once.
function shield_step!(deadCells::BitArray{N}, indices, can_act_matrix::BitArray{N}, no_action_bad_matrix::BitArray{N}, action_set::AbstractDict{Tuple{Int64,CartesianIndex{N}}}, no_action_set::AbstractDict{CartesianIndex{N}}, Act::Vector{Action}, dims::Int64; popped_actions::Set{Tuple{Int64,CartesianIndex{N}}}) where N
    act_idx_pop_keys = CartesianIndex{N}[]
    no_act_pop_keys = CartesianIndex{N}[]
    new_dead_cells = CartesianIndex{N}[]
    for idx in indices
        if deadCells[idx]
            continue
        end

        no_action_bad = no_action_bad_matrix[idx]
        if !no_action_bad
            succ = get(no_action_set, idx, nothing)
            # No entry means the flow hit a constraint, an empty entry means it left the invariant.
            if succ === nothing || isempty(succ) || any_at(deadCells, succ)
                no_action_bad = true
                no_action_bad_matrix[idx] = true
                push!(no_act_pop_keys, idx)
            end
        end

        approved_action_count = 0
        if can_act_matrix[idx]
            for act in Act
                key = (act.id, idx)
                succ = get(action_set, key, nothing)
                if succ === nothing || key in popped_actions
                    continue
                end
                if isempty(succ) || any_at(deadCells, succ) || any_at(no_action_bad_matrix, succ)
                    push!(popped_actions, key)
                    push!(act_idx_pop_keys, idx)
                else
                    approved_action_count += 1
                end
            end
            if approved_action_count <= 0
                can_act_matrix[idx] = false
            end
        end

        if approved_action_count <= 0 && no_action_bad
            can_act_matrix[idx] = false
            push!(new_dead_cells, idx)
        end
    end
    if !isempty(new_dead_cells)
        deadCells[new_dead_cells] .= true
    end

    return deadCells, act_idx_pop_keys, no_act_pop_keys, new_dead_cells
end

# action_set: (act, idx) -> idxs. Maps action id and cell idx to the idxs of all cells that overlap the final reachset of performing the action and then propagation.
# no_action_set: (idx) -> idxs. Maps cell idx to idxs of all cells that overlap the final reachset of propagation.
# no_act_bad_matrix: bitmatrix where true denotes that not performing the action results in not being a subset of the invariant set witin the control period. 
# act_bad_matrix: bitmatrix where true denotes that performing the action and then propagating results in not being a subset of the invariant set within the control period. 
function predecessor_graph!(action_set::ActionRelation{N}, no_action_set::FlowRelation{N}, no_act_bad_matrix::BitArray, act_bad_matrix::BitArray) where N
    # New dictionaries to store the predecessor relation. Now keys map to arrays of all idxs that reach said key.
    new_action_set::ActionRelation{N} = Dict()
    new_no_action_set::FlowRelation{N} = Dict()
    no_act_bad_matrix[findall(x -> !haskey(no_action_set, x), CartesianIndices(no_act_bad_matrix))] .= true
    for (k, vx) in pairs(no_action_set)
        if isempty(vx) # If vx is empty then the final reachset didnt map to any idxs implying that it is disjoint from the invariant set.
            no_act_bad_matrix[k] = true
        else
            for v in vx # Two cases: 1) some cell has already mapped to this idx in the original successor relation, hence we append the key. 2) it hasn't happened before and we must create a new array to store the key.
                if haskey(new_no_action_set, v)
                    push!(new_no_action_set[v], k)
                else
                    new_no_action_set[v] = CartesianIndex{N}[k]
                end
            end
        end
    end

    for (k, avx) in pairs(action_set)
        if !isempty(avx)
            act_bad_matrix[k[2]] = false # act_bad_matrix is false for all cells that had an enabled action
            for v in avx
                if haskey(new_action_set, (k[1], v))
                    push!(new_action_set[(k[1], v)], k[2])
                else
                    new_action_set[(k[1], v)] = CartesianIndex{N}[k[2]]
                end
            end
        else
            act_bad_matrix[k[2]] = true
        end
    end
    return new_action_set, new_no_action_set#, act_bad_matrix, no_act_bad_matrix
end


# TODO - Does not work because the predecessor graph computes the cells reached from an action in the stored other cell. Because of this we wrongly check whether we still are able to perform an action such that we dont become unsafe. 
function shield_from_predecessor_graph!(grid, act_list, action_set::ActionRelation{N}, no_action_set::FlowRelation{N}, no_act_bad_matrix::BitArray, act_bad_matrix::BitArray) where N
    waitlist = Set{CartesianIndex{N}}()
    donelist = Set{CartesianIndex{N}}()
    dead = grid.deadCells
    for idx in findall((|).(dead, (&).(act_bad_matrix, no_act_bad_matrix)))
        dead[idx] = true
        push!(waitlist, idx)
    end
    act_idxs = collect(a.id for a in act_list)
    #may_act_matrices = Dict(collect((x, copy(act_bad_matrix)) for x in act_idxs))
    may_act_matrices = [copy(.!act_bad_matrix) for _ in act_idxs]
    #@show size(may_act_matrices)
    #@show size(act_bad_matrix)
    #return act_bad_matrix, no_act_bad_matrix
    count = 0
    while !isempty(waitlist) && count < 15000000
        count += 1
        cidx = pop!(waitlist)
        if cidx in donelist
            continue
        end

        #no_act_bad_matrix[cidx] = true
        propagate = false
        if act_bad_matrix[cidx]
            propagate = true
        elseif !(any(map(i -> i[cidx], may_act_matrices)))
            #dead[cidx] = true
            propagate = true
        end

        if propagate && no_act_bad_matrix[cidx]
            dead[cidx] = true
            push!(donelist, cidx)
            if haskey(no_action_set, cidx)
                #@show no_action_set[cidx]
                no_act_bad_matrix[no_action_set[cidx]] .= true
                push!(waitlist, no_action_set[cidx]...)
            end
            for aidx in act_idxs
                may_act_matrices[aidx][cidx] = false

                if haskey(action_set, (aidx, cidx))
                    may_act_matrices[aidx][action_set[(aidx, cidx)]] .= false
                    act_bad_matrix[action_set[(aidx, cidx)]] .= true
                    push!(waitlist, action_set[(aidx, cidx)]...)
                    # should these get added to waitlist?
                end
            end
        end

    end
    println(count)
    #no_act_bad_matrix = (|).(no_act_bad_matrix, .!reduce((x, y)->(|).(x, y), may_act_matrices))
    #return act_bad_matrix, no_act_bad_matrix
end

function split_offsets(split_factor::Int64, dim::Int64)
    return (CartesianIndex).(Iterators.product(ntuple(i -> 0:-1:(-(split_factor-1)), dim)...))
end

function split_cell_idx(split_factor::Int64, car::CartesianIndex{N}, split_offsets) where N
    return split_offsets .+ (car .* split_factor)
end

function join_cell_idx(split_factor::Int64, car::CartesianIndex{N}) where N
    return CartesianIndex(ntuple(i -> fld((car[i] - 1), split_factor) + 1, N))
end