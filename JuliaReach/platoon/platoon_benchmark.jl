using Plots, Plots.PlotMeasures, LaTeXStrings, BenchmarkTools
using ReachabilityAnalysis, LazySets
using ReachabilityAnalysis.ReachabilityBase.Timing: print_timed

LazySets.deactivate_assertions()

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

function analyze_platoon(prob, alg, cmethod, imethod, dmin, case, warmup)
    # warm-up run
    if warmup
        analyze_once(prob, alg, cmethod, imethod, dmin)
    end

    b = @benchmarkable _ = analyze_once($prob, $alg, $cmethod, $imethod, $dmin)
    y = run(b)

    # benchmark run
    res = analyze_once(prob, alg, cmethod, imethod, dmin)
    sol, validation = res

    timeList = []
    for timeVal in y.times 
        push!(timeList, timeVal * 1e-9)
    end
    df = DataFrame(name=case, avgTime=mean(timeList), medianTime=median(timeList), success=validation)

    validation || throw(ErrorException("the property is violated"))

    return sol, df
end
