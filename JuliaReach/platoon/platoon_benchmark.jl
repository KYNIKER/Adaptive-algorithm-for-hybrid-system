module Platoon

using Plots, Plots.PlotMeasures, LaTeXStrings
using ReachabilityAnalysis, LazySets
using ReachabilityAnalysis.ReachabilityBase.Timing: print_timed

const model = "Platoon"
const cases = [
               # "PLAA01-BND50",
               # "PLAA01-BND50-discrete",
               # "PLAA01-BND42",
               # "PLAA01-BND42-discrete",
               "PLAD01-BND30",
               "PLAD01-BND30-discrete",
               # "PLAN01-UNB50",
               # "PLAN01-UNB50-discrete"
              ]

LazySets.deactivate_assertions()

# ==============================================================================
# Load model
# ==============================================================================

include("platoon.jl")

# ==============================================================================
# Analysis function
# ==============================================================================

function analyze_once(prob, alg, cmethod, imethod, dmin)
    # solve
    sol = solve(prob, alg=alg, clustering_method=cmethod,
                intersection_method=imethod, intersect_source_invariant=false,
                tspan = (0.0 .. 20.0))

    # verify that specification holds
    validation = (-ρ(sparsevec([1], [-1.0], 10), sol) > dmin) &&
                 (-ρ(sparsevec([4], [-1.0], 10), sol) > dmin) &&
                 (-ρ(sparsevec([7], [-1.0], 10), sol) > dmin)

    return (sol, validation)
end

function analyze(; prob, alg, cmethod, imethod, dmin, case, warmup)
    println("case: $case")

    # warm-up run
    if warmup
        analyze_once(prob, alg, cmethod, imethod, dmin)
    end

    # benchmark run
    res = @timed analyze_once(prob, alg, cmethod, imethod, dmin)
    print_timed(res)
    sol, validation = res.value
    validation || throw(ErrorException("the property is violated"))

    # report results
    io = isdefined(Main, :io) ? Main.io : stdout
    result = Int(validation)
    time = trunc(res.time, digits=3)
    print(io, "$model,$case,$result,$time\n")

    return sol
end

# ==============================================================================
# Run benchmarks
# ==============================================================================

const δ_discrete = 0.1

# ----------------------------------------
#  PLAD01-BND42
# ----------------------------------------

prob_PLAD01_BND42 = platoon(; deterministic_switching=true)
dirs = CustomDirections(collect(BoxDirections{Float64,Vector{Float64}}(10)))
imethod = TemplateHullIntersection(dirs)
dmin = -42.0

case = "PLAD01-BND42"
alg = BOX(δ=0.01)
cmethod = BoxClustering(1)
analyze(prob=prob_PLAD01_BND42, alg=alg, cmethod=cmethod, imethod=imethod,
        dmin=dmin, case=case, warmup=true)

case = "PLAD01-BND42-discrete"
alg = BOX(δ=δ_discrete, approx_model=NoBloating())
cmethod = BoxClustering(1, [3,1,1,1,1,1,1,1,1,1])
analyze(prob=prob_PLAD01_BND42, alg=alg, cmethod=cmethod, imethod=imethod,
        dmin=dmin, case=case, warmup=true)

# ----------------------------------------
#  PLAD01-BND30
# ----------------------------------------

prob_PLAD01_BND30 = platoon(; deterministic_switching=true)
cmethod = LazyClustering(1)
dirs = CustomDirections(collect(OctDirections{Float64,Vector{Float64}}(10)))
imethod = TemplateHullIntersection(dirs)
dmin = -30.0

case = "PLAD01-BND30"
alg = LGG09(δ=0.03, template=dirs, approx_model=Forward(setops=dirs))
sol_PLAD01_BND30 = analyze(prob=prob_PLAD01_BND30, alg=alg, cmethod=cmethod,
                           imethod=imethod, dmin=dmin, case=case, warmup=true)

case = "PLAD01-BND30-discrete"
alg = LGG09(δ=δ_discrete, template=dirs, approx_model=NoBloating())
analyze(prob=prob_PLAD01_BND30, alg=alg, cmethod=cmethod, imethod=imethod,
        dmin=dmin, case=case, warmup=true)

# ==============================================================================
# Plot
# ==============================================================================

const TARGET_FOLDER = isdefined(Main, :TARGET_FOLDER) ? Main.TARGET_FOLDER : @__DIR__

const fig = Plots.plot()
plot!(fig, sol_PLAD01_BND30, vars=(0, 1),
      linecolor=:blue, color=:blue, alpha=0.8, lw=1.0,
      tickfont=font(30, "Times"), guidefontsize=45,
      xlab=L"t", ylab=L"x_{1}",
      xtick=[0, 5, 10, 15, 20.], ytick=[-30, -20, -10, 0],
      xlims=(0., 20.), ylims=(-31, 7),
      bottom_margin=-6mm, left_margin=-2mm, right_margin=4mm, top_margin=0mm,
      size=(1000, 1000))
hline!(fig, [-30.0], lc=:red, ls=:dash, lw=2, lab="")

savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP25-JuliaReach-$model-PLAD01-BND30"))

sol_PLAD01_BND30 = nothing
GC.gc()

end
nothing