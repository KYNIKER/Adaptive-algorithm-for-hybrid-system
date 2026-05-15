module Spacecraft

using Plots, Plots.PlotMeasures, LaTeXStrings
using ReachabilityAnalysis, LazySets
using ReachabilityAnalysis.ReachabilityBase.Timing: print_timed

model = "Rendezvous"
cases = [
         "SRA01",
         "SRA01-discrete",
        ]
const TEST_LONG = false
if TEST_LONG
    push!(cases, "SRA04")
    # push!(cases, "SRA04-discrete")
    # push!(cases, "SRA05")
    # push!(cases, "SRA05-discrete")
    # push!(cases, "SRA06")
    # push!(cases, "SRA06-discrete")
    # push!(cases, "SRA07")
    # push!(cases, "SRA07-discrete")
    # push!(cases, "SRA08")
    # push!(cases, "SRA08-discrete")
    push!(cases, "SRU01")
    # push!(cases, "SRU01-discrete")
    push!(cases, "SRU02")
    push!(cases, "SRU02-discrete")
end

LazySets.deactivate_assertions()

# ==============================================================================
# Load model
# ==============================================================================

include("spacecraft.jl")

# ==============================================================================
# Analysis function
# ==============================================================================

const x  = 1  # x position (negative!)
const y  = 2  # y position (negative!)
const vx = 3  # x velocity
const vy = 4  # y velocity
const t  = 5  # time

# number of variables
const n = 4 + 1

# line-of-sight property in "attempt"
const cone = HPolyhedron([
    HalfSpace(sparsevec([x], [-1.0], n), 100.0),             # x >= -100
    HalfSpace(sparsevec([x, y], [tan(π/6), -1.0], n), 0.0),  # -x tan(30°) + y >= 0
    HalfSpace(sparsevec([x, y], [tan(π/6), 1.0], n), 0.0),   # -x tan(30°) - y >= 0
   ])
function line_of_sight(sol)
    all_idx = findall(L -> L == 2, location.(sol))  # attempt
    for idx in all_idx
        verif = all(set(R) ⊆ cone for R in sol[idx])
        !verif && return false
    end
    return true
end

# velocity constraint in "attempt"
const velocity = 0.055 * 60.0     # meters per minute
const cx = velocity * cos(π / 8)  # x-coordinate of the octagon's first (ENE) corner
const cy = velocity * sin(π / 8)  # y-coordinate of the octagon's first (ENE) corner
const octagon = HPolyhedron([
    HalfSpace(sparsevec([vx], [-1.0], n), cx),                # vx >= -cx
    HalfSpace(sparsevec([vx], [1.0], n), cx),                 # vx <= cx
    HalfSpace(sparsevec([vy], [-1.0], n), cx),                # vy >= -cx
    HalfSpace(sparsevec([vy], [1.0], n), cx),                 # vy <= cx
    HalfSpace(sparsevec([vx, vy], [1., 1.0], n), cy + cx),    # vx + vy <= cy + cx
    HalfSpace(sparsevec([vx, vy], [1., -1.0], n), cy + cx),   # vx - vy <= cy + cx
    HalfSpace(sparsevec([vx, vy], [-1., 1.0], n), cy + cx),   # -vx + vy <= cy + cx
    HalfSpace(sparsevec([vx, vy], [-1., -1.0], n), cy + cx)   # -vx - vy <= cy + cx
   ])
function velocity_constraint(sol)
    all_idx = findall(L -> L == 2, location.(sol))  # attempt
    for idx in all_idx
        verif = all(set(R) ⊆ octagon for R in sol[idx])
        !verif && return false
    end
    return true
end

# target set in "aborting"
const target = BallInf(zeros(2), 0.2)
function target_avoidance(sol)
   all_idx = findall(L -> L == 3, location.(sol))  # aborting
   for idx in all_idx
       # lazy:
       # verif = all(isdisjoint(set(Projection(R, [x, y])), target) for R in sol[idx])

       # concrete if R is hyperrectangular
       verif = all(isdisjoint(overapproximate(Projection(R, [x, y]), Hyperrectangle), target)
                   for R in sol[idx])
       !verif && return false
   end
   return true
end

const boxdirs = CustomDirections(collect(BoxDirections{Float64,Vector{Float64}}(5)))

function analyze_once(prob, alg, cmethod)
    # solve
    sol = solve(prob; alg=alg, clustering_method=cmethod,
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false, tspan=(0.0 .. 300.0))

    # verify that specification holds
    validation = line_of_sight(sol) && velocity_constraint(sol) && target_avoidance(sol)

    return (sol, validation)
end

function analyze(; prob, alg, cmethod, safe, case, warmup)
    println("case: $case")

    # warm-up run
    if warmup
        analyze_once(prob, alg, cmethod)
    end

    # benchmark run
    res = @timed analyze_once(prob, alg, cmethod)
    print_timed(res)
    sol, validation = res.value
    if safe
        validation || throw(ErrorException("the property is violated"))
    else
        !validation || throw(ErrorException("the property is unexpectedly verified"))
    end

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

const alg_discrete = BOX(δ=0.1, approx_model=NoBloating())

# ----------------------------------------
#  SRNA01
# ----------------------------------------

prob_SRNA01 = spacecraft(abort_time=-1.0)
cmethod = LazyClustering()

case = "SRNA01"
alg = BOX(δ=0.04)
sol_SRNA01 = analyze(prob=prob_SRNA01, alg=alg, cmethod=cmethod, safe=true, case=case, warmup=true)

case = "SRNA01-discrete"
analyze(prob=prob_SRNA01, alg=alg_discrete, cmethod=cmethod, safe=true, case=case, warmup=true)

# ----------------------------------------
#  SRA01
# ----------------------------------------

prob_SRA01 = spacecraft(abort_time=120.0)
cmethod = LazyClustering(3)

case = "SRA01"
alg_SRA01 = BOX(δ=0.04)
sol_SRA01 = analyze(prob=prob_SRA01, alg=alg_SRA01, cmethod=cmethod, safe=true, case=case, warmup=true)

case = "SRA01-discrete"
analyze(prob=prob_SRA01, alg=alg_discrete, cmethod=cmethod, safe=true, case=case, warmup=true)

# ----------------------------------------
#  SRA02
# ----------------------------------------

prob_SRA02 = spacecraft(abort_time=[120.0, 125])
cmethod = LazyClustering(3)

case = "SRA02"
alg_SRA02 = BOX(δ=0.04)
analyze(prob=prob_SRA02, alg=alg_SRA02, cmethod=cmethod, safe=true, case=case, warmup=false)

case = "SRA02-discrete"
analyze(prob=prob_SRA02, alg=alg_discrete, cmethod=cmethod, safe=true, case=case, warmup=false)

# ----------------------------------------
#  SRA03
# ----------------------------------------

prob_SRA03 = spacecraft(abort_time=[120.0, 145])
cmethod = LazyClustering(4)

case = "SRA03"
alg = BOX(δ=0.04)
analyze(prob=prob_SRA03, alg=alg, cmethod=cmethod, safe=true, case=case, warmup=false)

case = "SRA03-discrete"
analyze(prob=prob_SRA03, alg=alg_discrete, cmethod=cmethod, safe=true, case=case, warmup=false)

# ----------------------------------------
#  SRA04
# ----------------------------------------

if TEST_LONG
    prob_SRA04 = spacecraft(abort_time=240.0)

    case = "SRA04"
    alg = BOX(δ=0.01)
    cmethod = LazyClustering(40)
    analyze(prob=prob_SRA04, alg=alg, cmethod=cmethod, safe=true, case=case, warmup=false)

    #=
    case = "SRA04-discrete"
    cmethod = LazyClustering(16)
    analyze(prob=prob_SRA04, alg=alg_discrete, cmethod=cmethod, safe=true, case=case, warmup=false)
    =#
end

# ----------------------------------------
#  SRA05
# ----------------------------------------

if TEST_LONG
    #=
    prob_SRA05 = spacecraft(abort_time=[235.0, 240])
    cmethod = LazyClustering(13)

    case = "SRA05"
    alg = BOX(δ=0.04)
    analyze(prob=prob_SRA05, alg=alg, cmethod=cmethod, safe=true, case=case, warmup=false)

    case = "SRA05-discrete"
    analyze(prob=prob_SRA05, alg=alg_discrete, cmethod=cmethod, safe=true, case=case, warmup=false)
    =#
end

# ----------------------------------------
#  SRA06
# ----------------------------------------

if TEST_LONG
    #=
    prob_SRA06 = spacecraft(abort_time=[230.0, 240])
    cmethod = LazyClustering(16)

    case = "SRA06"
    alg = BOX(δ=0.04)
    analyze(prob=prob_SRA06, alg=alg, cmethod=cmethod, safe=true, case=case, warmup=false)

    case = "SRA06-discrete"
    analyze(prob=prob_SRA06, alg=alg_discrete, cmethod=cmethod, safe=true, case=case, warmup=false)
    =#
end

# ----------------------------------------
#  SRA07
# ----------------------------------------

if TEST_LONG
    #=
    prob_SRA07 = spacecraft(abort_time=[50.0, 150])
    cmethod = LazyClustering(50)

    case = "SRA07"
    alg = BOX(δ=0.04)
    analyze(prob=prob_SRA07, alg=alg, cmethod=cmethod, safe=true, case=case, warmup=false)

    case = "SRA07-discrete"
    analyze(prob=prob_SRA07, alg=alg_discrete, cmethod=cmethod, safe=true, case=case, warmup=false)
    =#
end

# ----------------------------------------
#  SRA08
# ----------------------------------------

if TEST_LONG
    #=
    prob_SRA08 = spacecraft(abort_time=[0.0, 240])

    case = "SRA08"
    alg = BOX(δ=0.04)
    cmethod = LazyClustering(800)
    analyze(prob=prob_SRA08, alg=alg, cmethod=cmethod, safe=true, case=case, warmup=false)

    case = "SRA08-discrete"
    cmethod = LazyClustering(900)
    analyze(prob=prob_SRA08, alg=alg_discrete, cmethod=cmethod, safe=true, case=case, warmup=false)
    =#
end

# ----------------------------------------
#  SRU01
# ----------------------------------------

if TEST_LONG
    prob_SRU01 = spacecraft(abort_time=260.0)
    cmethod = LazyClustering(16)

    case = "SRU01"
    # rule: use the same algorithm as for SRA01
    analyze(prob=prob_SRU01, alg=alg_SRA01, cmethod=cmethod, safe=false, case=case, warmup=false)

    #=
    case = "SRU01-discrete"
    analyze(prob=prob_SRU01, alg=alg_discrete, cmethod=cmethod, safe=false, case=case, warmup=false)
    =#
end

# ----------------------------------------
#  SRU02
# ----------------------------------------

if TEST_LONG
    prob_SRU02 = spacecraft(abort_time=[0.0, 260])
    cmethod = LazyClustering(16)

    case = "SRU02"
    # rule: use the same algorithm as for SRA02
    analyze(prob=prob_SRU02, alg=alg_SRA02, cmethod=cmethod, safe=false, case=case, warmup=false)

    case = "SRU02-discrete"
    analyze(prob=prob_SRU02, alg=alg_discrete, cmethod=cmethod, safe=false, case=case, warmup=false)
end

# ==============================================================================
# Plot
# ==============================================================================

const TARGET_FOLDER = isdefined(Main, :TARGET_FOLDER) ? Main.TARGET_FOLDER : @__DIR__

# SRNA01
fig = plot(legend=:bottomright, tickfont=font(30, "Times"), guidefontsize=45,
           xlab=L"x", ylab=L"y",
           xtick=[-800, -600, -400, -200.0, 0.0], ytick=[-400, -300, -200, -100, 0.0],
           xlims=(-1000.0, 0.0), ylims=(-450.0, 0.0),
           bottom_margin=-6mm, left_margin=-3mm, right_margin=1mm, top_margin=3mm,
           size=(1000, 1000))
for idx in findall(L -> L == 1, location.(sol_SRNA01))
    plot!(fig, sol_SRNA01[idx], vars=(1, 2), lw=0.0, color=:lightgreen, lc=:lightgreen, alpha=1.0)
end
for idx in findall(L -> L == 2, location.(sol_SRNA01))
    plot!(fig, sol_SRNA01[idx], vars=(1, 2), lw=0.0, color=:red, lc=:red, alpha=1.0)
end
for idx in findall(L -> L == 3, location.(sol_SRNA01))
    plot!(fig, sol_SRNA01[idx], vars=(1, 2), lw=0.0, color=:cyan, lc=:cyan, alpha=1.0)
end
savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP25-JuliaReach-$model-SRNA01"))

# SRA01
fig = plot(legend=:bottomright, tickfont=font(30, "Times"), guidefontsize=45,
           xlab=L"x", ylab=L"y", lw=0.0,
           xtick=[-750, -500, -250, 0, 250.0], ytick=[-400, -300, -200, -100, 0.0],
           xlims=(-1000.0, 400.0), ylims=(-450.0, 0.0),
           bottom_margin=-6mm, left_margin=-3mm, right_margin=1mm, top_margin=3mm,
           size=(1000, 1000))
for idx in findall(L -> L == 1, location.(sol_SRA01))
    plot!(fig, sol_SRA01[idx], vars=(1, 2), lw=0.0, color=:lightgreen, lc=:lightgreen, alpha=1.0)
end
for idx in findall(L -> L == 2, location.(sol_SRA01))
    plot!(fig, sol_SRA01[idx], vars=(1, 2), lw=0.0, color=:red, lc=:red, alpha=1.0)
end
for idx in findall(L -> L == 3, location.(sol_SRA01))
    plot!(fig, sol_SRA01[idx], vars=(1, 2), lw=0.0, color=:cyan, lc=:cyan, alpha=1.0)
end
savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP25-JuliaReach-$model-SRA01"))

sol_SRNA01 = nothing
sol_SRA01 = nothing
GC.gc()

end
nothing