using ReachabilityAnalysis, Plots, LazySets, BenchmarkTools, CSV, DataFrames
using ReachabilityAnalysis, LazySets


LazySets.deactivate_assertions()

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

function analyze_once_spacecraft(prob, alg, cmethod)
    # solve
    sol = solve(prob; alg=alg, clustering_method=cmethod,
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false, tspan=(0.0 .. 300.0))

    # verify that specification holds
    validation = line_of_sight(sol) && velocity_constraint(sol) && target_avoidance(sol)

    return (sol, validation)
end

function analyze_spacecraft(prob, alg, cmethod, safe, case, warmup)
    println("case: $case")

    # warm-up run
    if warmup
        analyze_once_spacecraft(prob, alg, cmethod)
    end

    b = @benchmarkable _ = analyze_once_spacecraft($prob, $alg, $cmethod)
    y = run(b)

    # benchmark run
    res = analyze_once_spacecraft(prob, alg, cmethod)
    sol, validation = res

    timeList = []
    for timeVal in y.times 
        push!(timeList, timeVal * 1e-9)
    end
    df = DataFrame(name=case, avgTime=mean(timeList), medianTime=median(timeList), success=validation)

    if safe
        validation || throw(ErrorException("the property is violated"))
    else
        !validation || throw(ErrorException("the property is unexpectedly verified"))
    end

    return sol, df
end


