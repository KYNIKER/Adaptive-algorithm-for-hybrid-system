using LazySets, LinearAlgebra, BenchmarkTools, CSV, DataFrames, Expokit, ReachabilityAnalysis

include("Utilities.jl")
include("ReACTed.jl")
include("models/gearbox.jl")
include("models/platoon.jl")
include("models/bouncingBall.jl")
include("models/powerTrain.jl")
include("models/spacecraft.jl")

const RUN_FIXED = true
const RUN_ADAPTIVE = true

const MAX_ORDER = 5
const REDUCE_ORDER = 5

function RunAdaptive(name, δ⁻, δ⁺, load_func, dirs)
    clustering = true
    timeConstraintList = []
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
    _ = ReACTed(sys, initialState, [0., T], X0, Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), [], sys.globalConstraints, δ⁻, δ⁺ , ReachabilityAnalysis.Exponentiation.BaseExp, MAX_ORDER, REDUCE_ORDER, timeConstraintList)

    # Actual Test
    b = @benchmarkable _ = ReACTed($sys, $initialState, [0., $T], $X0, $Zonotope(zeros(Float64, $n), $zeros(Float64, $n, 1)), [], $sys.globalConstraints, $δ⁻, $δ⁺ , ReachabilityAnalysis.Exponentiation.BaseExp, $MAX_ORDER, $REDUCE_ORDER, $timeConstraintList)

    y = run(b; verbose=true)
    println("Run completed.")

    # Convert time to seconds from nanoseconds
    timeList = []
    for timeVal in y.times
        push!(timeList, timeVal / 1e9)
    end
    # Write to csv file
    df = DataFrame(δ⁺=δ⁺, δ⁻=δ⁻, avgTime=mean(timeList), medianTime=median(timeList), memory=y.memory, allocs=y.allocs)

    filename = "results/ReACTRevised_" * name * "Results" * ".csv"
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

function RunFixed(name, δ, load_func, dirs)
    RunAdaptive(name, δ, δ, load_func, dirs)
end


names = ["gearbox", "platoon", "spaceCraft"]
minδ = [0.0008, 0.03, 0.04]
δ⁺arr = [2^1, 2^2, 2^3, 2^4, 2^6, 2^8, 2^10, 2^12]
dirs = [[3,4], [1,4,7], [1]]
loadFunctions = [loadGearBox, loadPlatoon, () -> loadSpacecraft(abort_time=120.)]

# Long Versions
# names = ["gearbox-01", "gearbox-02", "platoon", "powerTrain", "spaceCraft-0", "spaceCraft-120", "spaceCraft-240"]
# minδ = [0.05, 0.008, 0.008, 0.03, 0.002, 0.04, 0.04, 0.01]
# funcs = [loadGearBox, () -> loadGearBox(2), loadPlatoon, loadPowertrain, loadSpacecraft, ()-> loadSpacecraft(abort_time=120.), ()-> loadSpacecraft(abort_time=240.)]




if RUN_FIXED
    for (name, δ, loadFunction, dir) in zip(names, minδ, loadFunctions, dirs)
        RunFixed(name * "_Fixed", δ, loadFunction, dir)
    end
end

if RUN_ADAPTIVE
    for d in δ⁺arr
        for (name, δ⁻, loadFunction, dir) in zip(names, minδ, loadFunctions, dirs)
            δ⁺ = δ⁻ * d
            RunAdaptive(name * "_Adaptive", δ⁻, δ⁺, loadFunction, dir)
        end
    end
end

# myLoading = () -> loadSpacecraft(abort_time=120.)
# #myLoading = loadPlatoon
# myLoading = loadGearBox

# minDelta = 0.0008

# nameOfTest = "Gracie_Gearbox"

# RunFixed(nameOfTest, minDelta, myLoading)
# RunAdaptive(nameOfTest, minDelta, minDelta*2^3, myLoading)
# RunAdaptive(nameOfTest, minDelta, minDelta*2^5, myLoading)
# # RunAdaptive("NewGracieTest", minDelta, minDelta*2^7, myLoading)
# # RunAdaptive("NewGracieTest", minDelta, minDelta*2^9, myLoading)