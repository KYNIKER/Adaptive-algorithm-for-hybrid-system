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
granularity::Float64 = 0.02  # Example granularity
grid_name = "ball" * string(granularity)
savefile = grid_name * ".jld2"

euclideanHybridSystem, timePeriod = loadBouncingBallShieldedWithInput()

cscheme = palette(:matter, 2^length(euclideanHybridSystem.Act)+1; rev=true)

δ⁻ = 0.001

#GC.gc(true)
# TODO - Looking at the serialized reach by act dict it looks like it maps to cells that only have a idx but no other information.
#=
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
#VSCodeServer.@profview _, _, _ = ReACTed_reachable_cell_online(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻)
@time grid, act_set, no_act_set, no_act_bad_matrix = ReACTed_reachable_cell_online(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻)
#@show grid.deadCells# == trues(grid.numCells...)
shield = grid
#@time grid, reach_by_Act, reach_by_no_Act = ReACTed_reachable_cell_t(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻)
#ProfileView.@profview _, _, _ = ReACTed_reachable_cell_t(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻)
#ProfileView.@profview _, _, _ = ReACTed_reachable_cell_b(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻)
#@show grid.deadCells
#ProfileView.@profview _, _, _, _ = make_shield(grid, reach_by_Act, reach_by_no_Act, 70, euclideanHybridSystem.Act)
#@show typeof(reach_by_Act)

#shield, iters, act_set, no_act_set = make_shield(grid, reach_by_Act, reach_by_no_Act, 2, euclideanHybridSystem.Act)

#shield = grid
#act_set = reach_by_Act
#no_act_set = reach_by_no_Act

#@show (length(keys(act_set)), length(keys(no_act_set)))
#ProfileView.@profview shield, iters = make_shield(grid, reach_by_Act, reach_by_no_Act, 50, euclideanHybridSystem.Act)
#@time _, _ = make_shield(grid, reach_by_Act, reach_by_no_Act, 500, euclideanHybridSystem.Act)

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
    heatmap_matrix = ones(Int64, shield.numCells...)
    #@show keys(act_set)
    for idx in CartesianIndices(shield.deadCells)
        #@show cell.pCells


        #@show any(x -> haskey(reach_by_Act, (x, cell.id)), euclideanHybridSystem.Act)
        if any(x -> haskey(act_set, (x.id, idx)), euclideanHybridSystem.Act)
            heatmap_matrix[idx] = 2
            #=
            possible_acts = []
            for act in euclideanHybridSystem.Act

                if haskey(act_set, (act.id, idx))
                    #push!(actDict[act.id], zonotopeArray3d[idx])

                    #@show idx
                    push!(possible_acts, act.id)
                end

            end
            for actid in possible_acts
                heatmap_matrix[idx] = 2 #+= act_translation[actid]
            end
            =#
        end

        #push!(noAct, zonotopeArray3d[idx])
        if no_act_bad_matrix[idx]
            #push!(invalid_cells, idx)
            heatmap_matrix[idx] = 2
            #push!(unsafe, zonotopeArray3d[idx])
        end
        if shield.deadCells[idx]
            heatmap_matrix[idx] = 0
            #push!(unsafe, zonotopeArray3d[cell.id])
        else
            #heatmap_matrix[idx] = 1
        end
    end
    #heatmap_matrix[findall(p -> !isempty(p.second), act_set)] .= 2
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

    #for i in 270:280
    lower_offset = grid.lower .+ (granularity/2)

    A_abs = abs.(euclideanHybridSystem.flowMatrix)
    max_input = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, timePeriod, ReachabilityAnalysis.Exponentiation.BaseExp, false, nothing), euclideanHybridSystem.input)
    degenerate_dimensions = collect(u - l == l ? l - (granularity / 2) : 0. for (l, u) in zip(grid.lower, grid.upper))
    #@show degenerate_dimensions
    aexp = exp(timePeriod * euclideanHybridSystem.flowMatrix)
    Z = remove_zero_generators(Zonotope(degenerate_dimensions, diagm(map(x -> x == 1 ? 0. : granularity/2, grid.numCells))))
    #=
    for i in 319:320
        for j in 2:2
            tidx = CartesianIndex(i, j)
            of = lower_offset + ((car2vec(tidx) .- 1) .* granularity)
            LazySets.API.translate!(Z, of)
            @show tidx, get_touching_cell_idxs_lz(grid, minkowski_sum(linear_map(aexp, Z), max_input))
            #tidx = CartesianIndex(91, 1)
            heatmap_matrix[tidx] = 4
            #@show haskey(no_act_set, tidx)
            if haskey(no_act_set, tidx)
                for v in no_act_set[tidx]
                    if heatmap_matrix[v] == 0 || heatmap_matrix[v] == 3 || heatmap_matrix[v] == 2
                        heatmap_matrix[v] = 3
                        heatmap_matrix[tidx] = 2
                    else
                        heatmap_matrix[v] = 4
                    end
                    #heatmap_matrix[v] = 4

                end

                #@show no_act_set[tidx]
            end
            for key in euclideanHybridSystem.Act
                if haskey(act_set, (key.id, tidx))
                    for v in act_set[(key.id, tidx)]
                        if heatmap_matrix[v] == 0 || heatmap_matrix[v] == 3
                            heatmap_matrix[v] = 3
                        else
                            heatmap_matrix[v] = 2
                            #@show v
                        end
                    end
                end
            end
            LazySets.API.translate!(Z, -of)

        end
    end
    =#

    fig = Plots.plot(heatmap(transpose(heatmap_matrix), c=cgrad([:black, :white, :purple, :yellow, :red], 5, categorical=true), clim=(0, 4), colorbar=true, xlabel="v", ylabel="p", xticks=([0.0, shield.numCells[1]/2, shield.numCells[1]], [string(shield.lower[1]), "0", string(shield.upper[1])]), yticks=([0.0, shield.numCells[2]/2, shield.numCells[2]], [string(shield.lower[2]), string((shield.lower[2]+shield.upper[2])/2), string(shield.upper[2])])); colorbar_ticks=([0.0, 1.0, 2.0, 3.0, 4.0], ["dead", "no action", "act1", "act2", "act1 + act2"]))
    display(fig)
    savefig(fig, "plots/shield.pdf")
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