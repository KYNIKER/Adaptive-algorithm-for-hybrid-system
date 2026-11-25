# Synthesizes the shield of one model and saves a plot of it. Run from the repository root:
#
#     julia --project=. scripts/plot_models.jl <model> [key=value ...]
#
# Models and their defaults (every key can be given on the command line):
#   ball      bouncing ball with input              g=0.02 dmin=0.001 dmax=0.032 online=true exit_bad=false
#   oilpump   oil pump                              g=0.1  dmin=0.2   dmax=0.4   period=0.2 online=true
#   cruise    cruise control, a layer per ego mode  g=1.0  dmin=0.25  dmax=1.0   v_ego=0.0 sensor_range=200.0 online=true
#   cruise_t  cruise control, one decision layer    g=1.0  dmin=0.25  dmax=1.0   v_ego=0.0 sensor_range=200.0 online=true
#                                                   act_needs_safe_flow=false
#   boost     dc-dc boost converter (Karamanakos)   g=0.02 dmin=0.5   dmax=4.0   R=73.0 v_min=14.5 v_max=15.5 online=true
#             (time in µs: dmin and dmax are in µs, the period is 50 µs)
# Common keys: g (granularity), dmin and dmax (δ⁻ and δ⁺), online, act_needs_safe_flow, exit_bad (leaving the state space
# is bad, for waiting and for actions; true unless given), shield (modal, or adaptive: one level of refinement from 2g
# down to g), out (the png to write; by default results/<model>_<settings>.png). For example
#
#     julia --project=. scripts/plot_models.jl cruise_t g=0.5 v_ego=4
#     julia --project=. scripts/plot_models.jl oilpump period=1.0
#
# Each run appends a row to results/plot_models_timings.csv: the date, model, settings, the time of the synthesis in
# seconds and the dead and total cells. Before the timed run, the same model and settings are synthesized once on a
# coarse grid (WARMUP_G), so the time does not include compiling the shield code.

ENV["GKSwstype"] = "100" # no plot windows, only files
using LazySets, Plots, Dates#, ProfileView
gr()
include(joinpath(@__DIR__, "..", "src", "ReACTShields.jl"))
using .ReACTShields

# shield: {"modal", "adaptive"}. The ball leaves the state space by default, as its old online shield (pred_fast) did in effect.
const DEFAULTS = Dict(
    "ball" => Dict("g" => "0.02", "dmin" => "0.001", "dmax" => "0.032", "online" => "true", "exit_bad" => "false"),
    "oilpump" => Dict("g" => "0.2", "dmin" => "0.1", "dmax" => "0.2", "period" => "0.2", "online" => "true"),
    #"cruise" => Dict("g" => "1.0", "dmin" => "0.25", "dmax" => "1.0", "v_ego" => "0.0", "sensor_range" => "200.0", "online" => "true"),
    "cruise_t" => Dict("g" => "1.0", "dmin" => "0.125", "dmax" => "1.0", "v_ego" => "0.0", "sensor_range" => "200.0", "online" => "true",
        "act_needs_safe_flow" => "false"),
    "boost" => Dict("g" => "0.02", "dmin" => "0.5", "dmax" => "4.0", "R" => "73.0", "v_min" => "14.5", "v_max" => "15.5", "online" => "true"),
)

# Granularity of the warm-up run. For the cruise control it must keep the front car's region borders (-6 and 18) on
# cell borders.
const WARMUP_G = Dict("ball" => "0.5", "oilpump" => "1.0", "cruise" => "2.0", "cruise_t" => "2.0", "boost" => "0.1")

function parse_args(args)
    if !isempty(args)
        model = args[1]
        haskey(DEFAULTS, model) || error("unknown model $model; the models are $(join(sort(collect(keys(DEFAULTS))), ", "))")
        opts = copy(DEFAULTS[model])
        haskey(opts, "act_needs_safe_flow") || (opts["act_needs_safe_flow"] = "true")
        haskey(opts, "exit_bad") || (opts["exit_bad"] = "true")
        haskey(opts, "shield") || (opts["shield"] = "modal")
        for a in args[2:end]
            k, v = split(a, "="; limit=2)
            opts[k] = v
        end
        return [model], [opts]
    else
        models = []
        opts = []
        ball_configs = [
            Dict("g" => "0.02", "dmin" => "0.001", "dmax" => "0.032", "online" => "true", "exit_bad" => "false", "shield" => "modal"),
            Dict("g" => "0.01", "dmin" => "0.001", "dmax" => "0.032", "online" => "true", "exit_bad" => "false", "shield" => "modal"),
            Dict("g" => "0.02", "dmin" => "0.001", "dmax" => "0.032", "online" => "false", "exit_bad" => "false", "shield" => "modal"),
            Dict("g" => "0.01", "dmin" => "0.001", "dmax" => "0.032", "online" => "false", "exit_bad" => "false", "shield" => "modal"),
            Dict("g" => "0.01", "dmin" => "0.001", "dmax" => "0.032", "online" => "true", "exit_bad" => "false", "shield" => "adaptive"),
            Dict("g" => "0.01", "dmin" => "0.001", "dmax" => "0.032", "online" => "false", "exit_bad" => "false", "shield" => "adaptive")
        ]
        for bc in ball_configs
            push!(models, "ball")
            haskey(bc, "act_needs_safe_flow") || (bc["act_needs_safe_flow"] = "true")
            haskey(bc, "exit_bad") || (bc["exit_bad"] = "true")
            haskey(bc, "shield") || (bc["shield"] = "modal")
            push!(opts, bc)
        end

        op_configs = [
            Dict("g" => "0.2", "dmin" => "0.1", "dmax" => "0.2", "period" => "0.2", "online" => "true", "shield"=>"modal"),
            Dict("g" => "0.2", "dmin" => "0.1", "dmax" => "0.2", "period" => "0.2", "online" => "false", "shield"=>"modal"),
            Dict("g" => "0.1", "dmin" => "0.1", "dmax" => "0.2", "period" => "0.2", "online" => "true", "shield"=>"adaptive")
        ]
        for opc in op_configs
            push!(models, "oilpump")
            haskey(opc, "act_needs_safe_flow") || (opc["act_needs_safe_flow"] = "true")
            haskey(opc, "exit_bad") || (opc["exit_bad"] = "true")
            haskey(opc, "shield") || (opc["shield"] = "modal")
            push!(opts, opc)
        end

        cc_configs = [
            Dict("g" => "1.0", "dmin" => "0.125", "dmax" => "1.0", "v_ego" => "0.0", "sensor_range" => "200.0", "online" => "true", "act_needs_safe_flow" => "false", "shield"=>"modal"),
            Dict("g" => "0.5", "dmin" => "0.125", "dmax" => "0.5", "v_ego" => "0.0", "sensor_range" => "200.0", "online" => "true", "act_needs_safe_flow" => "false", "shield"=>"modal"),
            Dict("g" => "1.0", "dmin" => "0.125", "dmax" => "1.0", "v_ego" => "0.0", "sensor_range" => "200.0", "online" => "false", "act_needs_safe_flow" => "false", "shield"=>"modal"),
            Dict("g" => "0.5", "dmin" => "0.125", "dmax" => "1.0", "v_ego" => "0.0", "sensor_range" => "200.0", "online" => "true", "act_needs_safe_flow" => "false", "shield"=>"adaptive")
        ]

        for ccc in cc_configs
            push!(models, "cruise_t")
            haskey(ccc, "act_needs_safe_flow") || (ccc["act_needs_safe_flow"] = "true")
            haskey(ccc, "exit_bad") || (ccc["exit_bad"] = "true")
            haskey(ccc, "shield") || (ccc["shield"] = "modal")
            push!(opts, ccc)
        end

        boost_configs = [
            Dict("g" => "0.02", "dmin" => "0.5", "dmax" => "4.0", "R" => "73.0", "v_min" => "14.5", "v_max" => "15.5", "online" => "true", "shield"=>"modal"),
            Dict("g" => "0.02", "dmin" => "0.5", "dmax" => "4.0", "R" => "73.0", "v_min" => "14.5", "v_max" => "15.5", "online" => "false", "shield"=>"modal"),
            Dict("g" => "0.01", "dmin" => "0.5", "dmax" => "4.0", "R" => "73.0", "v_min" => "14.5", "v_max" => "15.5", "online" => "true", "shield"=>"modal"),
            Dict("g" => "0.01", "dmin" => "0.5", "dmax" => "4.0", "R" => "73.0", "v_min" => "14.5", "v_max" => "15.5", "online" => "true", "shield"=>"adaptive")
        ]

        for boostc in boost_configs
            push!(models, "boost")
            haskey(boostc, "act_needs_safe_flow") || (boostc["act_needs_safe_flow"] = "true")
            haskey(boostc, "exit_bad") || (boostc["exit_bad"] = "true")
            haskey(boostc, "shield") || (boostc["shield"] = "modal")
            push!(opts, boostc)
        end
        return models, opts
    end
end

# The bouncing ball as a heatmap of position against velocity, coloured by the set of choices allowed in each cell.
function plot_ball(grid, allowed, sys)
    names = ["wait"; ["act $(a.id)" for a in sys.Act]]
    label(b) = b == 0 ? "dead" : "{" * join([names[i+1] for i in 0:(length(names)-1) if (b >> i) & 1 != 0], ", ") * "}"
    present = sort(unique(allowed))
    colors = [:black, :steelblue, :darkgoldenrod, :yellow, :orchid4, :white, :gray60, :white]
    idx = Dict(b => i for (i, b) in enumerate(present))
    #map!(v -> v & 0x1 == 1 ? (v > 1 ? 3 : 1) : clamp(v, 0, 2), allowed)
    m = map(b -> idx[b], allowed)
    v = grid.lower[1] .+ grid.granularity .* ((1:size(allowed, 1)) .- 0.5)
    p = grid.lower[2] .+ grid.granularity .* ((1:size(allowed, 2)) .- 0.5)
    fig = heatmap(v, p, permutedims(m); c=cgrad(colors[1:length(present)], length(present), categorical=true), clim=(0.5, length(present) + 0.5),
        colorbar=false, xlabel="v", ylabel="p", legend=:outerright, size=(1150, 650), title="Bouncing ball shield, g = $(grid.granularity)")
    for b in present
        plot!(fig, Float64[], Float64[]; seriestype=:shape, color=colors[idx[b]], linecolor=:gray60, label=label(b))
    end
    println(join(("$(label(b)) $(count(==(b), allowed))" for b in present), ", "))
    return fig
end

# The boost converter as two heatmaps of current against voltage, one per switch position, coloured by the set of
# choices allowed in each cell. The dashed box is the safe set.
function plot_boost(grid, allowed, v_min, v_max, i_max)
    names = ["keep", "turn on", "turn off"]
    label(b) = b == 0 ? "dead" : "{" * join([names[i+1] for i in 0:2 if (b >> i) & 1 != 0], ", ") * "}"
    present = sort(unique(allowed))
    palette = Dict(0 => :black, 1 => :steelblue, 2 => :darkgoldenrod, 3 => :yellow, 4 => :orchid4, 5 => :white)
    idx = Dict(b => i for (i, b) in enumerate(present))
    colors = [get(palette, b, :gray60) for b in present]
    i = grid.lower[1] .+ grid.granularity .* ((1:size(allowed, 1)) .- 0.5)
    v = grid.lower[2] .+ grid.granularity .* ((1:size(allowed, 2)) .- 0.5)
    panels = map(((BOOST_OFF, "switch off"), (BOOST_ON, "switch on"))) do (layer, title)
        m = map(b -> idx[b], allowed[:, :, layer])
        fig = heatmap(i, v, permutedims(m); c=cgrad(colors, length(present), categorical=true), clim=(0.5, length(present) + 0.5),
            colorbar=false, xlabel="inductor current (A)", ylabel="output voltage (V)", title, legend=layer == BOOST_ON ? :outerright : false)
        plot!(fig, [0, i_max, i_max, 0], [v_min, v_min, v_max, v_max]; color=:red, ls=:dash, label=false)
        for b in present
            plot!(fig, Float64[], Float64[]; seriestype=:shape, color=colors[idx[b]], linecolor=:gray60, label=label(b))
        end
        fig
    end
    println(join(("$(label(b)) $(count(==(b), allowed))" for b in present), ", "))
    return plot(panels...; layout=(1, 2), size=(1400, 600), plot_title="Boost converter shield, g = $(grid.granularity)",
        left_margin=4Plots.mm, bottom_margin=4Plots.mm)
end

# Synthesizes the shield of the model. Returns the system and what the shield returns.
function synthesize(model, opts; verbose::Bool=true)
    num(k) = parse(Float64, opts[k])
    flag(k) = parse(Bool, opts[k])
    g, δ⁻, δ⁺ = num("g"), num("dmin"), num("dmax")
    sys, p = model == "ball" ? loadBouncingBallShieldedWithInput() :
             model == "oilpump" ? loadOilPumpModal(; period=num("period")) :
             model == "cruise" ? loadCruiseControlModal(; sensor_range=num("sensor_range")) :
             model == "boost" ? loadBoostConverterKaramanakos(; R=num("R"), v_min=num("v_min"), v_max=num("v_max")) :
             loadCruiseControlModal_t(; sensor_range=num("sensor_range"))
    msys = sys isa EuclideanHybridSystem ? ModalHybridSystem(sys) : sys
    kw = (; online=flag("online"), act_needs_safe_flow=flag("act_needs_safe_flow"), exit_bad_for_flow=flag("exit_bad"), exit_bad_for_actions=flag("exit_bad"), verbose)
    shield = opts["shield"]
    res = shield == "modal" ? ReACTed_reachable_cell_modal(msys, p, g, δ⁻, δ⁺; kw...) :
          shield == "adaptive" ? ReACTed_reachable_cell_modal_adaptive(msys, p, 2g, 1, δ⁻, δ⁺; kw...) :
          error("unknown shield $shield; use modal or adaptive")
    return sys, res
end

function plot_shield(model, opts, sys, (grid, allowed, lost, info))
    model == "ball" && return plot_ball(grid, dropdims(allowed; dims=3), sys) # one layer
    model == "oilpump" && return plot_oil_pump(grid, grid, allowed, lost, info)
    model == "boost" && return plot_boost(grid, allowed, parse(Float64, opts["v_min"]), parse(Float64, opts["v_max"]), 4.0)
    return plot_cruise(grid, allowed, info; v_ego=parse(Float64, opts["v_ego"]), wait_mode=model == "cruise_t" ? CRUISE_NEUTRAL : nothing)
end

# Appends a row to results/plot_models_timings.csv, with a header if the file is new.
function record_timing(row)
    file = joinpath(@__DIR__, "..", "results", "plot_models_timings.csv")
    header = "date,model,settings,synthesis_s,dead_cells,cells"
    quote_field(x) = (x=string(x); occursin(r"[,\"]", x) ? "\"" * replace(x, "\"" => "\"\"") * "\"" : x)
    new = !isfile(file)
    open(file, "a") do io
        new && println(io, header)
        println(io, join(quote_field.(row), ","))
    end
    return abspath(file)
end

models, optss = parse_args(ARGS)
for (model, opts) in zip(models, optss)
    settings = join(("$k=$(opts[k])" for k in sort(collect(keys(opts))) if k != "out"), "_")
    out = get(opts, "out", joinpath(@__DIR__, "..", "results", "$(model)_$(settings).png"))
    mkpath(dirname(out))
    println("warm-up at g = $(WARMUP_G[model])")
    #VSCodeServer.@profview synthesize(model, merge(opts, Dict("g" => WARMUP_G[model])); verbose=false)
    synthesize(model, merge(opts, Dict("g" => WARMUP_G[model])); verbose=false)
    synth = @timed synthesize(model, opts)
    sys, res = synth.value
    savefig(plot_shield(model, opts, sys, res), out)
    allowed = res[2]
    csv = record_timing((Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), model, settings, round(synth.time, digits=3), count(==(0), allowed), length(allowed)))
    println("synthesis $(round(synth.time, digits=1)) s")
    println("saved $(abspath(out))")
    println("timing appended to $csv")
end