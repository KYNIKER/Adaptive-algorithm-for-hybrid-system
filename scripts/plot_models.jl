# Synthesizes the shield of one model and saves a plot of it. Run from the repository root:
#
#     julia --project=. scripts/plot_models.jl <model> [key=value ...]
#
# Models and their defaults (every key can be given on the command line):
#   ball      bouncing ball with input, pred_fast   g=0.02 dmin=0.001 dmax=0.032 shield=pred_fast (or t, adaptive)
#   oilpump   oil pump, modal                       g=0.1  dmin=0.2   dmax=0.4   period=0.2 online=true
#   cruise    cruise control, a layer per ego mode  g=1.0  dmin=0.25  dmax=1.0   v_ego=0.0 sensor_range=200.0 online=true
#   cruise_t  cruise control, one decision layer    g=1.0  dmin=0.25  dmax=1.0   v_ego=0.0 sensor_range=200.0 online=true
#                                                   act_needs_safe_flow=false
# Common keys: g (granularity), dmin and dmax (δ⁻ and δ⁺), online, act_needs_safe_flow, out (the png to write; by default
# results/<model>_<settings>.png). For example
#
#     julia --project=. scripts/plot_models.jl cruise_t g=0.5 v_ego=4
#     julia --project=. scripts/plot_models.jl oilpump period=1.0
#
# Each run appends a row to results/plot_models_timings.csv: the date, model, settings, the time of the synthesis in
# seconds and the dead and total cells. Before the timed run, the same model and settings are synthesized once on a
# coarse grid (WARMUP_G), so the time does not include compiling the shield code.

ENV["GKSwstype"] = "100" # no plot windows, only files
using LazySets, Plots, Dates
gr()
include(joinpath(@__DIR__, "..", "src", "ReACTShields.jl"))
using .ReACTShields

const DEFAULTS = Dict(
    "ball" => Dict("g" => "0.02", "dmin" => "0.001", "dmax" => "0.032", "shield" => "pred_fast"),
    "oilpump" => Dict("g" => "0.2", "dmin" => "0.1", "dmax" => "0.2", "period" => "0.2", "online" => "true"),
    #"cruise" => Dict("g" => "1.0", "dmin" => "0.25", "dmax" => "1.0", "v_ego" => "0.0", "sensor_range" => "200.0", "online" => "true"),
    "cruise_t" => Dict("g" => "1.0", "dmin" => "0.125", "dmax" => "1.0", "v_ego" => "0.0", "sensor_range" => "200.0", "online" => "true",
        "act_needs_safe_flow" => "false"),
)

# Granularity of the warm-up run. For the cruise control it must keep the front car's region borders (-6 and 18) on
# cell borders.
const WARMUP_G = Dict("ball" => "0.5", "oilpump" => "1.0", "cruise" => "2.0", "cruise_t" => "2.0")

function parse_args(args)
    isempty(args) && error("give a model: $(join(sort(collect(keys(DEFAULTS))), ", "))")
    model = args[1]
    haskey(DEFAULTS, model) || error("unknown model $model; the models are $(join(sort(collect(keys(DEFAULTS))), ", "))")
    opts = copy(DEFAULTS[model])
    haskey(opts, "act_needs_safe_flow") || (opts["act_needs_safe_flow"] = "true")
    for a in args[2:end]
        k, v = split(a, "="; limit=2)
        opts[k] = v
    end
    return model, opts
end

# The bouncing ball as a heatmap of position against velocity, coloured by the set of choices allowed in each cell.
function plot_ball(grid, allowed, sys)
    names = ["wait"; ["act $(a.id)" for a in sys.Act]]
    label(b) = b == 0 ? "dead" : "{" * join([names[i+1] for i in 0:(length(names)-1) if (b >> i) & 1 != 0], ", ") * "}"
    present = sort(unique(allowed))
    colors = [:black, :steelblue, :darkgoldenrod, :yellow, :orchid4, :orchid, :gray60, :white]
    idx = Dict(b => i for (i, b) in enumerate(present))
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

# Synthesizes the shield of the model. Returns the system and what the shield returns.
function synthesize(model, opts; verbose::Bool=true)
    num(k) = parse(Float64, opts[k])
    flag(k) = parse(Bool, opts[k])
    g, δ⁻, δ⁺ = num("g"), num("dmin"), num("dmax")
    if model == "ball"
        sys, p = loadBouncingBallShieldedWithInput()
        shield = opts["shield"]
        res = shield == "pred_fast" ? ReACTed_reachable_cell_online_pred_fast(sys, p, g, δ⁻, δ⁺; act_needs_safe_flow=flag("act_needs_safe_flow"), verbose) :
              shield == "t" ? ReACTed_reachable_cell_t(sys, p, g, δ⁻, δ⁺; verbose) :
              shield == "adaptive" ? ReACTed_reachable_cell_adaptive_fast(sys, p, 2g, 1, δ⁻, δ⁺; patience=0, act_needs_safe_flow=flag("act_needs_safe_flow"), verbose) :
              error("unknown shield $shield; use pred_fast, t or adaptive")
        return sys, res
    end
    sys, p = model == "oilpump" ? loadOilPumpModal(; period=num("period")) :
             model == "cruise" ? loadCruiseControlModal(; sensor_range=num("sensor_range")) :
             loadCruiseControlModal_t(; sensor_range=num("sensor_range"))
    return sys, ReACTed_reachable_cell_modal(sys, p, g, δ⁻, δ⁺; online=flag("online"), act_needs_safe_flow=flag("act_needs_safe_flow"), verbose)
end

function plot_shield(model, opts, sys, (grid, allowed, lost, info))
    model == "ball" && return plot_ball(grid, allowed, sys)
    model == "oilpump" && return plot_oil_pump(grid, grid, allowed, lost, info)
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

model, opts = parse_args(ARGS)
settings = join(("$k=$(opts[k])" for k in sort(collect(keys(opts))) if k != "out"), "_")
out = get(opts, "out", joinpath(@__DIR__, "..", "results", "$(model)_$(settings).png"))
mkpath(dirname(out))
println("warm-up at g = $(WARMUP_G[model])")
synthesize(model, merge(opts, Dict("g" => WARMUP_G[model])); verbose=false)
synth = @timed synthesize(model, opts)
sys, res = synth.value
savefig(plot_shield(model, opts, sys, res), out)
allowed = res[2]
csv = record_timing((Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), model, settings, round(synth.time, digits=3), count(==(0), allowed), length(allowed)))
println("synthesis $(round(synth.time, digits=1)) s")
println("saved $(abspath(out))")
println("timing appended to $csv")
