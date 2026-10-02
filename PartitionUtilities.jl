using LazySets, LinearAlgebra, ReachabilityAnalysis
include("Utilities.jl")

### Following code is adapted from "https://github.com/AstridHornBrorholt/Shielded-Learning-for-Hybrid-Systems/blob/main/Shared%20Code/Squares.jl"
# Assumed that partitioning is axis-aligned, i.e., the partitioning is done along the axes of the state space. The Grid struct represents a grid in the state space with specified granularity and bounds.
struct Grid{T<:Real}
    dimension::Int
    granularity::T
    lower::Vector{T}
    upper::Vector{T}
    numCells::Vector{Int}
    deadCells::BitArray
    #array
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
    #array = Array{Cell}(undef, (numCells...))  # Create an array to hold the grid cells

    deadCells = falses(numCells...)#BitArray(undef, (numCells...,))  # Create an array to label the dead cells falses(numCells...)#

    #for id in CartesianIndices(array)
    #    array[id] = Cell(id, [LinearIndices(array)[id]], [])
    #end

    return Grid{T}(dimension, granularity, lower_bounds, upper_bounds, numCells, deadCells)#, array)
end

Base.show(io::IO, grid::Grid) = println(io,
    "Grid($(grid.granularity), $(grid.lower), $(grid.upper))")

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

# TODO - This can be optimized significantly by utilizing the fact that Zonotopes are centrally symmetric. First halve the bounding box, then sweep from both sides.
function get_touching_cell_idxs(grid::Grid, convexSet::LazySet)
    touching_cell_idxs = []

    lower_bounds, upper_bounds = clamp.(LazySets.low(convexSet), grid.lower, grid.upper), clamp.(LazySets.high(convexSet), grid.lower, grid.upper)
    lower_bounds = Int.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity) .+ 1)

    upper_bounds = Int.(ceil.(abs.(upper_bounds .- grid.lower) ./ grid.granularity)) #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity
    ranges = [lower_bounds[i]:max(upper_bounds[i], 1) for i in 1:grid.dimension]

    idxs = CartesianIndices((ranges...,))
    #@show idxs
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

    upper_bounds = Int.(ceil.(abs.(upper_bounds .- grid.lower) ./ grid.granularity)) #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity
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
function get_touching_cell_idxs_b(grid::Grid, convexSet::Zonotope)#, hbox::Hyperrectangle)
    touching_cell_idxs = []

    lower_bounds, upper_bounds = clamp.(LazySets.low(convexSet), grid.lower, grid.upper), clamp.(LazySets.high(convexSet), grid.lower, grid.upper)
    lower_bounds = Int.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity) .+ 1)

    upper_bounds = Int.(ceil.(abs.(upper_bounds .- grid.lower) ./ grid.granularity)) #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity

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
        for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], 1)), grid.dimension))
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
        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:max(upper_bounds[i], 1), grid.dimension))
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

function get_touching_cell_idxs_l(grid::Grid, convexSet::LazySet)
    touching_cell_idxs = []

    lower_bounds, upper_bounds = clamp.(LazySets.low(convexSet), grid.lower, grid.upper), clamp.(LazySets.high(convexSet), grid.lower, grid.upper)
    lower_bounds = Int.(floor.(abs.(lower_bounds .- grid.lower) ./ grid.granularity) .+ 1)

    upper_bounds = Int.(ceil.(abs.(upper_bounds .- grid.lower) ./ grid.granularity)) #floor.(min.(upper_bounds, grid.upper) .- grid.lower) ./ grid.granularity

    sweeping_dim = argmax(upper_bounds[i]-lower_bounds[i] for i in 1:grid.dimension)
    hbox = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, grid.dimension))
    if upper_bounds[sweeping_dim] - lower_bounds[sweeping_dim] > 1
        #println("f'ing hope so")
        #upper_bounds[splitting_dim] = lower_bounds[splitting_dim] + cld(upper_bounds[splitting_dim] - lower_bounds[splitting_dim], 2) #cld(upper_bounds[splitting_dim], 2)


        cartesian_offset = CartesianIndex(ntuple(i -> i == sweeping_dim ? 1 : 0, grid.dimension))

        sweeping_range = range(lower_bounds[sweeping_dim], upper_bounds[sweeping_dim])
        rev_sweeping_range = range(upper_bounds[sweeping_dim], lower_bounds[sweeping_dim]; step=-1)
        for idx in CartesianIndices(ntuple(i -> i == sweeping_dim ? (0:0) : (lower_bounds[i]:max(upper_bounds[i], 1)), grid.dimension))
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


        for idx in CartesianIndices(ntuple(i -> lower_bounds[i]:max(upper_bounds[i], 1), grid.dimension))
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

# https://github.com/AstridHornBrorholt/Shielded-Learning-for-Hybrid-Systems/blob/22c9fc220ef40d55877ff1360be7286a7a506620/Shared%20Code/ShieldSynthesis.jl#L88
function make_shield(grid::Grid, action_set::Dict{Tuple{Action,CartesianIndex},AbstractArray{CartesianIndex}}, no_action_set::Dict{CartesianIndex,AbstractArray{CartesianIndex}}, max_steps::Int, Act::Vector{Action})
    i = max_steps
    dims = grid.dimension
    dead = grid.deadCells
    can_act_matrix = trues(grid.numCells...)
    no_action_bad_matrix = falses(grid.numCells...)
    dead´ = nothing
    action_set´ = nothing
    no_action_set´ = nothing
    a_tombstones = 0
    n_tombstones = 0
    #filter!(p -> !isempty(p.second), no_action_set)
    filter!(p -> !isempty(p.second), action_set)
    action_set = Dict{Tuple{Action,CartesianIndex},AbstractArray{CartesianIndex}}(action_set)
    no_action_set = Dict{CartesianIndex,AbstractArray{CartesianIndex}}(no_action_set)

    @show length(keys(action_set))
    @show length(keys(no_action_set))

    while i > 0
        #grid´, action_set´ = shield_step!(grid, action_set, no_action_set, Act)
        #dead´, la, ln, ld, action_set´, no_action_set´ = shield_step!(dead, can_act_matrix, no_action_bad_matrix, action_set, no_action_set, Act, dims)
        dead´, la, ln, ld = shield_step!(dead, can_act_matrix, no_action_bad_matrix, action_set, no_action_set, Act, dims)
        #@show can_act_matrix == trues(grid.numCells...)
        if la == 0 && ln == 0 && ld == 0

            println("Fixed point found at $(max_steps-i) steps!")
            dead = dead´
            #action_set = action_set´
            #no_action_set = no_action_set´
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
            action_set = Dict{Tuple{Action,CartesianIndex},AbstractArray{CartesianIndex}}(action_set)
            no_action_set = Dict{CartesianIndex,AbstractArray{CartesianIndex}}(no_action_set)
            a_tombstones = 0
            n_tombstones = 0

        end

        #@show action_set == action_set´
        i -= 1

    end
    @show length(keys(action_set))
    @show length(keys(no_action_set))

    @show can_act_matrix == trues(grid.numCells...)
    @show no_action_bad_matrix == falses(grid.numCells...)

    return (grid, max_steps - i, action_set, no_action_set)
end

function shield_step!(deadCells::BitArray, can_act_matrix, no_action_bad_matrix, action_set::Dict{Tuple{Action,CartesianIndex},AbstractArray{CartesianIndex}}, no_action_set::Dict{CartesianIndex,AbstractArray{CartesianIndex}}, Act::Vector{Action}, dims::Int64)
    #deadCells´ = copy(deadCells)
    act_pop_keys = Tuple{Action,CartesianIndex{dims}}[]
    no_act_pop_keys = CartesianIndex{dims}[]
    new_dead_cells = CartesianIndex{dims}[]
    for idx in CartesianIndices(deadCells)
        if !deadCells[idx]
            #no_action_bad = any(i -> i == 0, collect(grid.array[nc].sidx[1] for nc in cell.pCells))
            #no_action_bad = any(grid.deadCells[idx] for idx in cell.pCells)
            no_action_bad = no_action_bad_matrix[idx]
            if !no_action_bad && haskey(no_action_set, idx)
                #@show collect(deadCells[idxx] for idxx in no_action_set[idx])
                if any(deadCells[idxx] for idxx in no_action_set[idx])
                    #println("has key")
                    no_action_bad = true
                    no_action_bad_matrix[idx] = true
                    push!(no_act_pop_keys, idx)
                end
            end

            #can_act = can_act_matrix[idx]
            approved_action_count = 0
            if can_act_matrix[idx]
                for act in Act
                    if haskey(action_set, (act, idx))
                        approved_action_count += 1
                        #=
                        if isempty(collect(deadCells[idxx] for idxx in action_set[(act, idx)]))
                            #@show action_set[(act, cell.id)]
                            push!(act_pop_keys, (act, idx))
                            approved_action_count -= 1
                        end
                        =#
                        if any(deadCells[idxx] for idxx in action_set[(act, idx)])
                            push!(act_pop_keys, (act, idx))
                            approved_action_count -= 1
                        end
                    end
                end
            end
            #@show can_act
            #@show can_act, isempty(cell.pCells)
            if approved_action_count <= 0 && no_action_bad #isempty(cell.pCells) && cell.sidx[1] == 1 #no_action_bad && can_act == 0
                #grid´.array[cell.id].sidx[1] = 1
                #deadCells´[idx] = true
                #push!(no_act_pop_keys, cell.id)
                #@show grid´.deadCells[cell.id] == grid.deadCells[cell.id]

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