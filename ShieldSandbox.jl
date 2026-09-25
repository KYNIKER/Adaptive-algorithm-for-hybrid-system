using LazySets, LinearAlgebra, Plots, ReachabilityAnalysis, JLD2, FileIO
using Plots.PlotMeasures
include("PartitionUtilities.jl")
include("plotFuncs/plotFuncHelper.jl")
include("Utilities.jl")
include("models/bouncingBall.jl")
include("ReACTedShielding.jl")

fresh_grid = false
make_plot = false
granularity = 0.5  # Example granularity
grid_name = "ball" * string(granularity)
savefile = grid_name * ".jld2"

euclideanHybridSystem, timePeriod = loadBouncingBallShieldedWithInput()

cscheme = palette(:default, length(euclideanHybridSystem.Act))

δ⁻ = 0.001

# TODO - Looking at the serialized reach by act dict it looks like it maps to cells that only have a idx but no other information.
#=
GC.gc(true)
grid = nothing
reach_by_Act = nothing
if fresh_grid
    if isfile(savefile)
        #@show length(keys(reach_by_Act´))

        #@show grid.deadCells
        rm(savefile)
        #sets = ReACTedShieldingK(euclideanHybridSystem, 2 * timePeriod, 1, granularity, 0.001, 0.004)
    end
    grid´, reach_by_Act´ = nothing, nothing
    @time grid´, reach_by_Act´ = ReACTed_reachable_cell(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^0 * δ⁻)
    jldsave(savefile; g=grid´, r=reach_by_Act´)
end
f = jldopen(savefile)
@show f
grid = read(f, "g")
reach_by_Act = read(f, "r")
close(f)
@show reach_by_Act
=#

@time grid, reach_by_Act, reach_by_no_Act = ReACTed_reachable_cell(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^6 * δ⁻)
#=
@show length(reach_by_Act), length(reach_by_no_Act)

@show grid.array[40, 6]
@show grid.array[40, 5]
@show grid.deadCells[40, 6]
@show grid.deadCells[40, 5]
=#
#tidx = CartesianIndex(20, 3)
#=
tidx = CartesianIndex(30, 20)
@show haskey(reach_by_no_Act, tidx)
if haskey(reach_by_no_Act, tidx)
    @show reach_by_no_Act[tidx]
end

for act in euclideanHybridSystem.Act
    @show haskey(reach_by_Act, (act, tidx))
    if haskey(reach_by_Act, (act, tidx))
        @show reach_by_Act[(act, tidx)]
    end
end
=#
@time shield, iters = make_shield(grid, reach_by_Act, reach_by_no_Act, 15, euclideanHybridSystem.Act)
#=
@show length(reach_by_Act), length(reach_by_no_Act)

@show shield.array[40, 6]
@show shield.array[40, 5]
@show shield.deadCells[40, 6]
@show shield.deadCells[40, 5]
#@show reach_by_Act
=#
if make_plot
    unsafe = []
    jumping = []
    noAct = []

    actDict = Dict()

    for act in euclideanHybridSystem.Act
        actDict[act] = []
    end

    zonotopeArray3d = initialize_zonotope_array(shield)  # Initialize the zonotope array for the grid
    invalid_cells = []
    for idx in CartesianIndices(shield.deadCells)
        #@show cell.pCells

        if shield.deadCells[idx]
            #push!(unsafe, zonotopeArray3d[cell.id])
        end
        #@show any(x -> haskey(reach_by_Act, (x, cell.id)), euclideanHybridSystem.Act)
        if any(x -> haskey(reach_by_Act, (x, idx)), euclideanHybridSystem.Act)
            for act in euclideanHybridSystem.Act
                if haskey(reach_by_Act, (act, idx))
                    push!(actDict[act], zonotopeArray3d[idx])
                end
            end
            #@show zonotopeArray3d[cell.id]
            #push!(jumping, zonotopeArray3d[idx])
            #=if !isempty(cell.pCells)
            for pcell in cell.pCells
                push!(jumping, zonotopeArray3d[pcell])
            end
        end=#
        else
            push!(noAct, zonotopeArray3d[idx])
            if !haskey(reach_by_no_Act, idx)
                push!(invalid_cells, idx)
                push!(unsafe, zonotopeArray3d[idx])
            end
        end
    end

    tidx = CartesianIndex(30, 20)

    @show length(unsafe)
    @show length(jumping)
    @show length(noAct)
    @show CartesianIndex((40, 20)) ∈ invalid_cells
    dirs = [1, 2]
    #grid = Grid(euclideanHybridSystem.statespace, granularity)

    #unsafeDict = Dict{CartesianIndex,Vector{LazySet}}()

    #rect(x, y) = Shape(x .- 1 .* granularity .+ [0, granularity, granularity, 0, 0], y .- 1 .* granularity .+ [0, 0, granularity, granularity, 0])

    plt = plot(dpi=1200, thickness_scaling=1, guidefontsize=35, minorgrid=true, ε=granularity,
        #legendfont=font(12, "Times"),
        #legend_position=:topright,
        legend=false,
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
        plot!(plt, z, alpha=0.1, c=:blue)
    end

    for (i, action) in pairs(euclideanHybridSystem.Act)
        plot!(plt, action.guard, c=cscheme[i], fillstyle=://)
    end


    #@show jumping
    eA = exp(timePeriod * euclideanHybridSystem.flowMatrix)
    #tU = linear_map((inv(euclideanHybridSystem.flowMatrix) * (eA - I)), euclideanHybridSystem.input)
    A_abs = abs.(euclideanHybridSystem.flowMatrix)
    tU = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, timePeriod, ReachabilityAnalysis.Exponentiation.BaseExp, false, nothing), euclideanHybridSystem.input)
    tz = minkowski_sum(tU, linear_map(eA, zonotopeArray3d[tidx]))
    #for id in get_touching_cell_idxs(grid, tz)
    #@show get_cell_bounds(grid, grid.array[id])
    #end
    #@show get_cell_bounds(grid, grid.array[tidx])
    #@show get_cell_bounds(grid, zonotopeArray3d[tidx])
    gFrontierZonotopes = get_grid_shapes(jumping, dirs)#[get_grid_shapes(jumping, dirs); get_grid_shapes(map(x -> linear_map(eA, x), jumping), dirs)]#get_grid_shapes(jumping, dirs)

    for (i, ac) in pairs(euclideanHybridSystem.Act)
        tempzs = get_grid_shapes(actDict[ac], dirs)
        for z in tempzs
            plot!(plt, z, c=cscheme[i])
        end
    end
    #@show intersection(euclideanHybridSystem.globalConstraints[1], hyperrectangle_to_HPolytope(euclideanHybridSystem.statespace))
    plot!(plt, intersection(euclideanHybridSystem.globalConstraints[1], hyperrectangle_to_HPolytope(euclideanHybridSystem.statespace)), c=:red, lw=1.0, lab="unsafe")
    #plot!(plt, LazySets.API.project(LinearMap(exp(-timePeriod * euclideanHybridSystem.flowMatrix), intersection(euclideanHybridSystem.globalConstraints[1], hyperrectangle_to_HPolytope(euclideanHybridSystem.statespace))), dirs), c=:white)
    plot!(plt, intersection(euclideanHybridSystem.edges[1].guard, hyperrectangle_to_HPolytope(euclideanHybridSystem.statespace)), c=:green)
    #plot!(plt, euclideanHybridSystem.edges[1].guard, c=:green)
    plot!(plt, zonotopeArray3d[tidx], c=:white)
    plot!(plt, tz, c=:white)
    display(plt)  # Display the plot

end



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