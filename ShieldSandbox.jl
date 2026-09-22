using LazySets, LinearAlgebra, Plots, ReachabilityAnalysis, JLD2, FileIO
using Plots.PlotMeasures
include("PartitionUtilities.jl")
include("plotFuncs/plotFuncHelper.jl")
include("Utilities.jl")
include("models/bouncingBall.jl")
include("ReACTedShielding.jl")

fresh_grid = false
granularity = 0.25  # Example granularity
grid_name = "ball" * string(granularity)
savefile = grid_name * ".jld2"

euclideanHybridSystem, timePeriod = loadBouncingBallShieldedWithInput()

δ⁻ = 0.001

grid = nothing
reach_by_Act = nothing
if fresh_grid
    @time grid´, reach_by_Act´ = ReACTed_reachable_cell(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^0 * δ⁻)
    @show length(keys(reach_by_Act´))
    if isfile(savefile)
        #@show grid.deadCells
        rm(savefile)
        #sets = ReACTedShieldingK(euclideanHybridSystem, 2 * timePeriod, 1, granularity, 0.001, 0.004)
    end
    @save savefile grid´ reach_by_Act´
    grid = grid´
    reach_by_Act = reach_by_Act´
else
    @load savefile grid´ reach_by_Act´
    grid = grid´
    reach_by_Act = reach_by_Act´
end
@time shield, iters = make_shield(grid, reach_by_Act, 60, euclideanHybridSystem.Act)
#@show reach_by_Act
unsafe = []
jumping = []
noAct = []
@show grid.deadCells == shield.deadCells
zonotopeArray3d = initialize_zonotope_array(shield)  # Initialize the zonotope array for the grid
for cell in shield.array
    #@show cell.pCells

    if shield.deadCells[cell.id]
        push!(unsafe, zonotopeArray3d[cell.id])
    end
    #@show any(x -> haskey(reach_by_Act, (x, cell.id)), euclideanHybridSystem.Act)
    if any(x -> haskey(reach_by_Act, (x, cell.id)), euclideanHybridSystem.Act)
        #@show zonotopeArray3d[cell.id]
        push!(jumping, zonotopeArray3d[cell.id])
        #=if !isempty(cell.pCells)
            for pcell in cell.pCells
                push!(jumping, zonotopeArray3d[pcell])
            end
        end=#
    else
        push!(noAct, zonotopeArray3d[cell.id])
    end
end

@show length(unsafe)
@show length(jumping)
@show length(noAct)
dirs = [1, 2]
#grid = Grid(euclideanHybridSystem.statespace, granularity)

#unsafeDict = Dict{CartesianIndex,Vector{LazySet}}()


plt = plot(dpi=1200, thickness_scaling=1, guidefontsize=25, minorgrid=true,
    legendfont=font(12, "Times"),
    legend_position=:topright,
    tickfont=font(8, "Times"),
    xguidefont=font(12, "Times"),
    yguidefont=font(12, "Times"),
    bottom_margin=2mm,
    left_margin=5mm,
    right_margin=5mm,
    top_margin=2mm,
    xlabel="v", ylabel="p")



#mark_dead_cells!(grid, euclideanHybridSystem.globalConstraints[1], unsafeDict)

#unsafeCells, frontierCells = get_contained_edge_cells(grid, linear_map(exp(-0.01 * euclideanHybridSystem.flowMatrix), zonotopeArray3d[4, 24, 1]))
#@show unsafeCells
#@show frontierCells
#frontierZonotopes = zonotopeArray3d[collect(c.id for c in frontierCells)]
#@show frontierZonotopes
pFrontierZonotopes = get_grid_shapes(unsafe, dirs)

for z in pFrontierZonotopes
    plot!(plt, z, alpha=0.9, c=:black)
end

nFrontierZonotopes = get_grid_shapes(noAct, dirs)

for z in nFrontierZonotopes
    plot!(plt, z, c=:blue)
end

for action in euclideanHybridSystem.Act
    plot!(plt, action.guard, c=:yellow, fillstyle=://)
end


#@show jumping
eA = exp(timePeriod * euclideanHybridSystem.flowMatrix)
gFrontierZonotopes = get_grid_shapes(jumping, dirs)#[get_grid_shapes(jumping, dirs); get_grid_shapes(map(x -> linear_map(eA, x), jumping), dirs)]#get_grid_shapes(jumping, dirs)

for z in gFrontierZonotopes
    plot!(plt, z, c=:green)
end
#@show intersection(euclideanHybridSystem.globalConstraints[1], hyperrectangle_to_HPolytope(euclideanHybridSystem.statespace))
plot!(plt, intersection(euclideanHybridSystem.globalConstraints[1], hyperrectangle_to_HPolytope(euclideanHybridSystem.statespace)), c=:black, lw=5.0, fillstyle=://, lab="unsafe")
#plot!(plt, LazySets.API.project(LinearMap(exp(-timePeriod * euclideanHybridSystem.flowMatrix), intersection(euclideanHybridSystem.globalConstraints[1], hyperrectangle_to_HPolytope(euclideanHybridSystem.statespace))), dirs), c=:white)
plot!(plt, intersection(euclideanHybridSystem.edges[1].guard, hyperrectangle_to_HPolytope(euclideanHybridSystem.statespace)), c=:green)
#plot!(plt, euclideanHybridSystem.edges[1].guard, c=:green)

display(plt)  # Display the plot





#=


nextFrontierZonotopes = map(z -> linear_map(eA, z), frontierZonotopes)

nextPFrontierZonotopes = get_grid_shapes(nextFrontierZonotopes, dirs)

for z in nextPFrontierZonotopes
    plot!(plt, z, c=:blue)
end



#=
S = Zonotope([7.5, 0.0], [8.5 0.0; 0.0 15.0])  # Example zonotope in 2D
grid = Grid(S, granularity)


@time zonotopeArray = initialize_zonotope_array(grid)  # Initialize the zonotope array for the grid

#@show zonotopeArray

@show CartesianIndices(grid.array)

A = [0.0 0.0; 1. 0.0]

b = [-9.81, 0]
U = Zonotope(b, zeros(2, 2))
=#

#=
for idx in eachindex(zonotopeArray)
    if idx % 2 == 0
        plot!(plt, zonotopeArray[idx], vars=(1, 2), c=:red, alpha=0.3, lw=0.65, label="")
    else
        plot!(plt, zonotopeArray[idx], vars=(1, 2), c=:blue, alpha=0.3, lw=0.65, label="")
    end
end

eA = exp(A * 0.1)

plot!(plt, eA * Zonotope([0.0, 0.0], Diagonal(fill(grid.granularity * 3.5, grid.dimension))), vars=(1, 2), c=:black, label="")
#=
for idx in eachindex(zonotopeArray)
    if idx % 2 == 0
        plot!(plt, eA * zonotopeArray[idx], vars=(1, 2), c=:red, alpha=0.3, lw=0.65, label="")
    else
        plot!(plt, eA * zonotopeArray[idx], vars=(1, 2), c=:blue, alpha=0.3, lw=0.65, label="")
    end
end
=#

#initialize_safe_cells!(grid)  # Initialize the safe cells in the grid

#@time box(grid, [0.5, 0.5])  # Example state to find the corresponding cell
@time box(grid, [-0.195, -0.195])  # Example state to find the corresponding cell using fastbox
@time fastbox(grid, [-0.195, -0.195])  # Example state to find the corresponding cell using fastbox
#@time get_cell_bounds(grid, fastbox(grid, [0., -1.]))  # Example to get the bounds of the corresponding cell

for idx in get_touching_cells(grid, eA * Zonotope([0.0, 0.0], Diagonal(fill(grid.granularity * 3.5, grid.dimension))))
    
    plot!(plt, zonotopeArray[idx.id], vars=(1, 2), c=:blue, alpha=0.3, lw=0.65, label="")
end
for idx in get_contained_cells(grid, eA * Zonotope([0.0, 0.0], Diagonal(fill(grid.granularity * 3.5, grid.dimension))))
    plot!(plt, zonotopeArray[idx.id], vars=(1, 2), c=:white, alpha=0.3, lw=0.65, label="")
end


@time get_contained_cells(grid, Zonotope([0.0, 0.0], Diagonal(fill(grid.granularity * 3, grid.dimension))))

_, edgeCells = get_contained_edge_cells(grid, eA * Zonotope([0.0, 0.0], Diagonal(fill(grid.granularity * 3.5, grid.dimension))))

for idx in edgeCells #get_perimeter_cells(grid, eA * Zonotope([0.0, 0.0], Diagonal(fill(grid.granularity * 3.5, grid.dimension))))
    plot!(plt, zonotopeArray[idx.id], vars=(1, 2), c=:orange, alpha=0.3, lw=0.65, label="")
end
#=
=#

#AA = AAPolytope([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]])  # Example axis-aligned polytope in 2D
#@show convert(HPolytope, AA)  # Convert the axis-aligned polytope to an HPolytope


@show grow_indices(grid, [[1, 1], [1, 2]])

@show offsets(grid)

mark_dead_cells!(grid, eA * Zonotope([0.0, 0.0], Diagonal(fill(grid.granularity * 3.5, grid.dimension))), unsafeDict)

for v in values(unsafeDict)
    #@show v
    plot!(plt, v[1], c=:green)
end


inputSet = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A, 0.1, ReachabilityAnalysis.Exponentiation.BaseExp, false, nothing), U)

for v in values(unsafeDict)
    #@show v
    plot!(plt, minkowski_sum(inputSet, linear_map(eA, v[1])), c=:yellow)
end

dirs = [1, 2]
setSets = zonotopeArray3d[1, 1:2]
@show typeof(setSets)
shapes = get_grid_shapes(setSets, dirs)

for v in shapes
    #@show v
    plot!(plt, v, c=:red)
end

=#
=#
# Måske muligt i stedet for at genbruge koden fra Astrid, at bruge den samme funktion til at lave en grid som bare er en store af keys og så gemme values et andet sted. Hvis det er implementeret med et linært index
# eller som en dictionary på cartisianIndex ville man nok kunne fjerne keys som er cell'er der ikke længere har safe elementer. 
# Så ville man stadig kunne bruge bounds funktioner til at finde de tætteste celler til constraints.
# Man kunne tjekke at når en cell ikke længere er safe så fjerner man de eventuelle constraints den har medført og bytter dem ud med cell da den nok er simplere