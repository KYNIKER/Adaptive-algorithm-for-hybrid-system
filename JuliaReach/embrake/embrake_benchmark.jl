using Plots, Plots.PlotMeasures, LaTeXStrings, BenchmarkTools
using ReachabilityAnalysis, LazySets
using ReachabilityAnalysis.ReachabilityBase.Timing: print_timed
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector

LazySets.deactivate_assertions()



# ==============================================================================
# Analysis function
# ==============================================================================

const x0 = 0.05
const eₓ = SingleEntryVector(2, 4, 1.0)

function analyze_once(prob, alg, validate)
    # solve
    sol = solve(prob, max_jumps=1001, alg=alg)

    max_t = tend(sol)
    if validate
        # verify that specification holds
        validation = ρ(eₓ, sol) < x0
    else
        # maximum time the solution is below the x0 threshold
        validation = true
        for fp in sol, (j, Rj) in enumerate(fp)
            if ρ(eₓ, Rj) >= x0
                validation = false
                max_t = inf(tspan(fp, j))
                break
            end
        end
    end

    return (sol, validation, max_t)
end

function analyze_embrake(prob, alg, case, validate, warmup)
    # warm-up run
    if warmup
        analyze_once(prob, alg, validate)
    end

    # benchmark run
    b = @benchmarkable _ = analyze_once($prob, $alg, $validate)
    y = run(b)

    println("Proceeding to data collection...")
    res = analyze_once(prob, alg, validate)
    sol, validation, max_t = res

    timeList = []
    for timeVal in y.times 
        push!(timeList, timeVal * 1e-9)
    end

    df = DataFrame(name=case, avgTime=mean(timeList), medianTime=median(timeList), success=validation)

    # Throw an error in case of no validation
    !validate || validation || throw(ErrorException("the property is violated"))
    !validate && println("maximum time that x < x0 , case $c : $max_t")

    return sol, df
end

