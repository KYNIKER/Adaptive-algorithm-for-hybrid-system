using LazySets, LinearAlgebra, BenchmarkTools, CSV, DataFrames, Expokit, ReachabilityAnalysis

include("Utilities.jl")
include("ReACTed.jl")
include("models/gearbox.jl")
include("models/platoon.jl")
include("models/bouncingBall.jl")
include("models/powerTrain.jl")
include("models/spacecraft.jl")
include("models/embrake/solve_embrake.jl")
include("models/embrake/embrake.jl")

const RUN_FIXED = true
const RUN_ADAPTIVE = true

const MAX_ORDER = 5
const REDUCE_ORDER = 5

function RunAdaptive(name, δ⁻, δ⁺, load_func)
    clustering = true
    timeConstraintList = []
    if occursin("spacecraft", lowercase(name))
        println("Deteced SpaceCraft")
        clustering = false
    end

    if occursin("gearbox", lowercase(name))
        println("Deteced Gearbox")
        clustering = false
        timeConstraintList = [(1, 0.2)]
    end

    sys, initialState, X0, T = load_func()
    n = length(X0.center)
    LazySets.load_expokit()
    println("Running benchmark for: ", name)
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 50

    # Clear Memory
    GC.gc()

    # Warmup
    _ = ReACTed(sys, initialState, [0., T], X0, Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), [], sys.globalConstraints, δ⁻, δ⁺ , ReachabilityAnalysis.Exponentiation.BaseExp, MAX_ORDER, REDUCE_ORDER, clustering, timeConstraintList)

    # Actual Test
    b = @benchmarkable _ = ReACTed($sys, $initialState, [0., $T], $X0, $Zonotope(zeros(Float64, $n), $zeros(Float64, $n, 1)), [], $sys.globalConstraints, $δ⁻, $δ⁺ , ReachabilityAnalysis.Exponentiation.BaseExp, $MAX_ORDER, $REDUCE_ORDER, $clustering, $timeConstraintList)

    y = run(b; verbose=true)
    println("Run completed.")

    # Convert time to seconds from nanoseconds
    timeList = []
    for timeVal in y.times
        push!(timeList, timeVal / 1e9)
    end
    # Write to csv file
    df = DataFrame(δ⁺=δ⁺, δ⁻=δ⁻, avgTime=mean(timeList), medianTime=median(timeList), memory=y.memory, allocs=y.allocs)

    filename = "results/ReACT_" * name * "Results" * ".csv"
    if isfile(filename)# Check if file exists
        open(filename, "a") do File
            CSV.write(File, df, delim=";", append=true)
        end
    else
        open(filename, "w") do File
            CSV.write(File, df, delim=";", writeheader=true)
        end
    end

    println("Finished writing to file")
end

function RunBrake(name, δ⁺, δ⁻, load_func)
    LazySets.load_expokit()
    println("Running benchmark for: ", name)
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 50

    sys, initialState, X0, T = load_func()

    # Clear Memory
    GC.gc()

    # Warmup
    _ = solve_embrake(sys, initialState, X0, T, δ⁺, δ⁻, 5, 5, [], false)

    # Actual Test
    b = @benchmarkable _ = solve_embrake($sys, $initialState, $X0, $T, $δ⁺, $δ⁻, $5, $5, $[], $false)
    y = run(b; verbose=true)

    println("Run completed.")

    # Convert time to seconds from nanoseconds
    timeList = []
    for timeVal in y.times
        push!(timeList, timeVal / 1e9)
    end
    # Write to csv file
    df = DataFrame(δ⁺=δ⁺, δ⁻=δ⁻, avgTime=mean(timeList), medianTime=median(timeList), memory=y.memory, allocs=y.allocs)

    filename = "results/ReACT_" * name * "Results" * ".csv"
    if isfile(filename)# Check if file exists
        open(filename, "a") do File
            CSV.write(File, df, delim=";", append=true)
        end
    else
        open(filename, "w") do File
            CSV.write(File, df, delim=";", writeheader=true)
        end
    end

    println("Finished writing to file")
end
function RunFixed(name, δ, load_func)
    RunAdaptive(name, δ, δ, load_func)
end


names = ["brake", "gearbox", "platoon", "spaceCraft"]
minδ = [2*10^-7, 0.0008, 0.03, 0.04]
loadFunctions = [loadembrake, loadGearBox, loadPlatoon, () -> loadSpacecraft(abort_time=120.)]

# Long Versions
# names = ["gearbox-01", "gearbox-02", "platoon", "powerTrain", "spaceCraft-0", "spaceCraft-120", "spaceCraft-240"]
# minδ = [0.05, 0.008, 0.008, 0.03, 0.002, 0.04, 0.04, 0.01]
# funcs = [loadGearBox, () -> loadGearBox(2), loadPlatoon, loadPowertrain, loadSpacecraft, ()-> loadSpacecraft(abort_time=120.), ()-> loadSpacecraft(abort_time=240.)]




if RUN_FIXED
    for (name, δ, loadFunction) in zip(names, minδ, loadFunctions)
        if !occursin("embrake", lowercase(name))
            RunFixed(name * "_Fixed", δ, loadFunction)
        else
            RunBrake(name * "_Fixed", δ, δ, loadFunction)
            #RunBrake(name * "_Adaptive", 2^6*δ, δ, loadFunction)
        end
    end
end

if RUN_ADAPTIVE
    δ⁺arr = [2^1, 2^2, 2^3, 2^4, 2^6, 2^8, 2^10, 2^12]
    for d in δ⁺arr
        for (name, δ⁻, loadFunction) in zip(names, minδ, loadFunctions)
            δ⁺ = δ⁻ * d
            if !occursin("embrake", lowercase(name))
                RunAdaptive(name * "_Adaptive", δ⁻, δ⁺, loadFunction)
            else
                RunBrake(name * "_Adaptive", δ⁺, δ⁻, loadFunction)
            end
        end
    end
end
