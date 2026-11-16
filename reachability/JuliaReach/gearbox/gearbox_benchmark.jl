using Plots, Plots.PlotMeasures, LaTeXStrings, BenchmarkTools
using ReachabilityAnalysis, LazySets, SparseArrays
using ReachabilityAnalysis.ReachabilityBase.Timing: print_timed
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector


LazySets.deactivate_assertions()
const n = 6

function template_directions(n)
    vx, vy, px, py, θ = 1, 2, 3, 4, 0.628318530717959
    vec1 = sparsevec([px, py], [tan(θ), 1.], n)
    vec2 = sparsevec([px, py], [tan(θ), -1.], n)
    vec1perp = sparsevec([px, py], [1, -tan(θ)], n)
    vec2perp = sparsevec([px, py], [1., tan(θ)], n)
    dirs_plane_34 = [vec1, -vec1, vec2, -vec2, vec1perp, -vec1perp, vec2perp, -vec2perp]

    vec3 = sparsevec([vx, vy], [-sin(θ), -cos(θ)], n)
    vec4 = sparsevec([vx, vy], [-sin(θ), cos(θ)], n)
    vec3perp = sparsevec([vx, vy], [cos(θ), -sin(θ)], n)
    vec4perp = sparsevec([vx, vy], [cos(θ), sin(θ)], n)

    dirs_plane_12 = [vec3, -vec3, vec4, -vec4, vec3perp, -vec3perp, vec4perp, -vec4perp]
    boxdirs = BoxDirections{Float64,Vector{Float64}}(n)
    return vcat(collect(boxdirs), Vector.(dirs_plane_34), Vector.(dirs_plane_12)) |> CustomDirections
end
# ==============================================================================
# Analysis function
# ==============================================================================
const extdirs = template_directions(n)
const thull = TemplateHullIntersection(extdirs)

function analyze_once_gearbox(prob, alg)
        # solve
        sol = solve(prob, max_jumps=100, clustering_method=LazyClustering(1, convex=false),
                    intersect_source_invariant=true, intersection_source_invariant_method=thull,
                    intersection_method=thull, tspan = 0 .. 0.21, fixpoint_check=true, alg=alg)

        # verify that specification holds
        validation = true
        # property 1
        for F in array(sol.F)
            if location(F) == 1
                if tend(F) >= 0.2
                    validation = false
                    break
                end
            end
        end
        # property 2
        # would like to use n instead of 6, but for some reason n = 5 in this scope? 
        e5 = SingleEntryVector(5, 6, 1.0)
        validation &= ρ(e5, sol) < 20.0

        return (sol, validation)
    end

function analyze_gearbox(prob, alg, case, warmup)

    # warm-up run
    if warmup
        analyze_once_gearbox(prob, alg)
    end

    b = @benchmarkable _ = analyze_once_gearbox($prob, $alg)
    y = run(b)

    println("Proceeding to data collection...")
    res = analyze_once_gearbox(prob, alg)
    sol, validation = res

    timeList = []
    for timeVal in y.times 
        push!(timeList, timeVal * 1e-9)
    end
    df = DataFrame(name=case, avgTime=mean(timeList), medianTime=median(timeList), success=validation)

    validation || throw(ErrorException("the property is violated"))

    return sol, df
end
