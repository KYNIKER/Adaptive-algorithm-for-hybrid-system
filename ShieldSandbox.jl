using LazySets, LinearAlgebra, Plots, ReachabilityAnalysis, JLD2, FileIO, BenchmarkTools, ProfileView
using Plots.PlotMeasures

gr()
include("PartitionUtilities.jl")
include("plotFuncs/plotFuncHelper.jl")
include("Utilities.jl")
include("models/bouncingBall.jl")
include("ReACTedShielding.jl")
include("ReACTedShieldingPred.jl")

fresh_grid = false
make_plot = true
granularity::Float64 = 0.02  # Example granularity
grid_name = "ball" * string(granularity)
savefile = grid_name * ".jld2"

euclideanHybridSystem, timePeriod = loadBouncingBallShieldedWithInput()

cscheme = palette(:matter, 2^length(euclideanHybridSystem.Act)+1; rev=true)

δ⁻ = 0.001

include("ReACTedShieldingAdaptiveFast.jl")
include("ReACTedShieldingPredFast.jl")
#@time grid, allowed, lost, info = ReACTed_reachable_cell_t(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻)
#@time grid, allowed, lost, info = ReACTed_reachable_cell_online_pred(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻)
#@time grid, allowed, lost, info = ReACTed_reachable_cell_online_pred_fast(euclideanHybridSystem, timePeriod, granularity, δ⁻, 2^5 * δ⁻; jump_hull=false)
@time grid, allowed, lost, info = ReACTed_reachable_cell_adaptive_fast(euclideanHybridSystem, timePeriod, granularity, 1, δ⁻, 2^5 * δ⁻, patience=0)

#can_wait = (allowed .& 0x01) .!= 0
#can_act(ai) = ((allowed .>> ai) .& 0x01) .!= 0   # ai is the position in Act
# Every shield returns grid, allowed, lost, info (see choice_bits in ReACTedShieldingFlat.jl); info holds what only that shield has, here info.act_bad and info.shield_count_cover
shield = grid

#grid = grid_p
#display(plot(heatmap(transpose(shield_cover_matrix))))

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

        heatmap_matrix[idx] = allowed[idx]
        if (allowed[idx] .& 0x01) .!= 0
            heatmap_matrix[idx] = 1
        end

        if allowed[idx] == 1 && lost[idx] > 1
            heatmap_matrix[idx] = 3
        elseif allowed[idx] > 1 && heatmap_matrix[idx] != 1
            heatmap_matrix[idx] = 2
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

    #plot!(plt, heatmap(transpose(heatmap_matrix), c=cscheme, colorbar=nothing))

    #for i in 270:280
    lower_offset = grid.lower .+ (granularity/2)

    A_abs = abs.(euclideanHybridSystem.flowMatrix)
    max_input = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, timePeriod, ReachabilityAnalysis.Exponentiation.BaseExp, false, nothing), euclideanHybridSystem.input)
    degenerate_dimensions = collect(u - l == l ? l - (granularity / 2) : 0. for (l, u) in zip(grid.lower, grid.upper))
    #@show degenerate_dimensions
    aexp = exp(timePeriod * euclideanHybridSystem.flowMatrix)
    Z = remove_zero_generators(Zonotope(degenerate_dimensions, diagm(map(x -> x == 1 ? 0. : granularity/2, grid.numCells))))

    println("amount of dead cells: $(length(findall(==(0), heatmap_matrix))/foldl(*, grid.numCells))")
    println("amount of blocked cells: $(length(findall(==(2), heatmap_matrix))/foldl(*, grid.numCells))")

    y_n = floor(abs(shield.lower[2]) + abs(shield.upper[2]))
    y_nt = collect(fld(shield.numCells[2], y_n) * i for i in 1:y_n)
    y_st = collect(string(shield.lower[2] + (i * (abs(shield.lower[2]) + abs(shield.upper[2])) / y_n)) for i in 1:y_n)
    x_n = floor(abs(shield.lower[1]) + abs(shield.upper[1]))
    x_nt = collect(fld(shield.numCells[1], x_n) * i for i in 1:2:x_n)
    x_st = collect(string(shield.lower[1] + (i * (abs(shield.lower[1]) + abs(shield.upper[1])) / x_n)) for i in 1:2:x_n)
    hm = heatmap(transpose(heatmap_matrix), c=cgrad([:black, :white, :purple, :yellow, :red], 5, categorical=true), clim=(0, 4), colorbar=false, xlabel="v", ylabel="p", xticks=(x_nt, x_st), yticks=(y_nt, y_st))
    #Plots.translate!(hm.shape(), 0, 0, -15)
    fig = Plots.plot(hm; colorbar_ticks=([0.0, 1.0, 2.0, 3.0, 4.0], ["dead", "no action", "act1", "act2", "act1 + act2"]))
    color_labels = ["{}", "{hit, nohit}", "{hit}", "{nohit}"]
    for (color, label) in zip(cgrad([:black, :white, :purple, :yellow]), color_labels)
        Plots.plot!(Float64[], Float64[], seriestype=:shape,   # The Float typing has been added to avoid a weird warning
            label=label, color=color)
    end
    Plots.xlims!(0.0, shield.numCells[1])
    Plots.ylims!(1.0, shield.numCells[2])
    display(fig)
    #savefig(fig, "plots/shield.pdf")
end