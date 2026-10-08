using LazySets, LinearAlgebra, Plots, ReachabilityAnalysis, JLD2, FileIO, BenchmarkTools, ProfileView
using Plots.PlotMeasures

gr()
include("PartitionUtilities.jl")
include("plotFuncs/plotFuncHelper.jl")
include("Utilities.jl")
include("models/bouncingBall.jl")
include("ReACTedShielding.jl")

fresh_grid = false
make_plot = true
granularity = 0.02  # Example granularity
grid_name = "ball" * string(granularity)
savefile = grid_name * ".jld2"

euclideanHybridSystem, timePeriod = loadBouncingBallShieldedWithInput()

cscheme = palette(:matter, 2^length(euclideanHybridSystem.Act)+1; rev=true)

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

#@allocations ReACTed_reachable_cell_b(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻)
@time grid, reach_by_Act, reach_by_no_Act = ReACTed_reachable_cell_b(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻)
#ProfileView.@profview _, _, _ = ReACTed_reachable_cell_b(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻)
@show length(values(reach_by_no_Act))
@show length(unique(values(reach_by_no_Act)))
@show in([], values(reach_by_no_Act))
@show length(values(reach_by_Act))
@show length(unique(values(reach_by_Act)))
@show in([], values(reach_by_Act))

#=
test_Z = Zonotope([-13.1, 0.0], diagm([1., 1.1]))
#ProfileView.@profview _ = get_touching_cell_idxs_t(grid, test_Z)
_ = get_touching_cell_idxs_b(grid, test_Z)
_ = get_touching_cell_idxs_t(grid, test_Z)
_ = get_touching_cell_idxs(grid, test_Z)
@time res1 = get_touching_cell_idxs_b(grid, test_Z)
@time res2 = get_touching_cell_idxs_t(grid, test_Z)
@time res3 = get_touching_cell_idxs(grid, test_Z)
@show res1
@show res3
@show issubset(res1, res2)
@show issubset(res2, res1)
@show issubset(res1, res3)
=#
#=
@show maximum(length, values(reach_by_Act))
@show maximum(length, values(reach_by_no_Act))
@show minimum(length, values(reach_by_Act))
@show minimum(length, values(reach_by_no_Act))
=#
#@time _, _, _ = ReACTed_reachable_cell_b(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻)
#=
_, _, _ = ReACTed_reachable_cell(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻)
@time _, _, _ = ReACTed_reachable_cell(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻)
grid, reach_by_Act, reach_by_no_Act = ReACTed_reachable_cell(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻)
=#
#ProfileView.@profview _, _, _ = ReACTed_reachable_cell_b(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^6 * δ⁻)
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
#@show (length(keys(reach_by_Act)), length(keys(reach_by_no_Act)))
@time shield, iters, act_set, no_act_set = make_shield(grid, reach_by_Act, reach_by_no_Act, 1200, euclideanHybridSystem.Act)
#@show (length(keys(act_set)), length(keys(no_act_set)))
#ProfileView.@profview shield, iters = make_shield(grid, reach_by_Act, reach_by_no_Act, 50, euclideanHybridSystem.Act)
#@time _, _ = make_shield(grid, reach_by_Act, reach_by_no_Act, 500, euclideanHybridSystem.Act)
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

    act_translation = Dict()
    for (i, act) in pairs(euclideanHybridSystem.Act)
        act_translation[act.id] = 2#^(i-1)
    end

    for act in euclideanHybridSystem.Act
        actDict[act.id] = []
    end

    #zonotopeArray3d = initialize_zonotope_array(shield)  # Initialize the zonotope array for the grid
    invalid_cells = []
    heatmap_matrix = zeros(Int64, shield.numCells...)
    #@show keys(act_set)
    for idx in CartesianIndices(shield.deadCells)
        #@show cell.pCells


        #@show any(x -> haskey(reach_by_Act, (x, cell.id)), euclideanHybridSystem.Act)
        if any(x -> haskey(act_set, (x.id, idx)), euclideanHybridSystem.Act)
            possible_acts = []
            for act in euclideanHybridSystem.Act

                if haskey(act_set, (act.id, idx))
                    #push!(actDict[act.id], zonotopeArray3d[idx])

                    #@show idx
                    push!(possible_acts, act.id)
                end

            end
            for actid in possible_acts
                heatmap_matrix[idx] += act_translation[actid]
            end
        end

        #push!(noAct, zonotopeArray3d[idx])
        if !haskey(no_act_set, idx)
            push!(invalid_cells, idx)
            #push!(unsafe, zonotopeArray3d[idx])
        else
            heatmap_matrix[idx] += 1

        end
        if shield.deadCells[idx]
            heatmap_matrix[idx] = 0
            #push!(unsafe, zonotopeArray3d[cell.id])
        end
    end

    #colors = cgrad(cscheme, length(euclideanHybridSystem.Act) + 1, categorical=true)
    #=
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
    =#
    #plot(hm)
    color_labels = ["Dead", "No action", "Act1", "Act2", "Act1 + Act2"]
    if length(color_labels) > 0
        if length(color_labels) != length(cscheme)
            throw(ArgumentError("Length of argument color_labels does not match  number of colors."))
        end
        for (color, label) in zip(cscheme, color_labels)
            # Apparently shapes are added to the legend even if the list is empty
            #plot!(plt, Float64[], Float64[], seriestype=:shape, label=label, color=color)
        end
    end
    #plot!(plt, heatmap(transpose(heatmap_matrix), c=cscheme, colorbar=nothing))

    tidx = CartesianIndex(18, 8)
    #tidx = CartesianIndex(91, 1)
    heatmap_matrix[tidx] = 4
    @show haskey(no_act_set, tidx)
    if haskey(no_act_set, tidx)
        for v in no_act_set[tidx]
            if heatmap_matrix[v] == 0
                heatmap_matrix[v] = 3
            else
                heatmap_matrix[v] = 4
            end
            #heatmap_matrix[v] = 4

        end
        @show no_act_set[tidx]
    end

    plot(heatmap(transpose(heatmap_matrix), c=cgrad([:black, :white, :purple, :yellow, :red], 5, categorical=true), clim=(0, 4), colorbar=true, xlabel="v", ylabel="p", xticks=([0.0, shield.numCells[1]/2, shield.numCells[1]], [string(shield.lower[1]), "0", string(shield.upper[1])]), yticks=([0.0, shield.numCells[2]/2, shield.numCells[2]], [string(shield.lower[2]), string((shield.lower[2]+shield.upper[2])/2), string(shield.upper[2])])); colorbar_ticks=([0.0, 1.0, 2.0, 3.0, 4.0], ["dead", "no action", "act1", "act2", "act1 + act2"]))
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