using ReachabilityAnalysis, Plots, LazySets, BenchmarkTools, CSV, DataFrames

include("embrake/embrake.jl")
include("embrake/embrake_benchmark.jl")
include("gearbox/gearbox.jl")
include("gearbox/gearbox_benchmark.jl")
include("platoon/platoon.jl")
include("platoon/platoon_benchmark.jl")
#include("spacecraft/spacecraft_benchmark.jl")

# template directions

const extdirs = template_directions(6)
const dirs = CustomDirections(collect(OctDirections{Float64,Vector{Float64}}(10)))

cases = [
    #"BRKDC01",
    #"BRKDC01-discrete",
    #"GRBX01-MES01",
    #"GRBX01-MES01-discrete",
    "PLAD01-BND30",
    "PLAD01-BND30-discrete",
    #="SRA01",
    "SRA01-discrete"=#
]

algCheckDict = Dict(
    "BRKDC01" =>  GLGM06(δ=2e-7, max_order=1, static=true, dim=4, ngens=4, approx_model=Forward()),
    "BRKDC01-discrete" => GLGM06(δ=1e-8, max_order=1, static=true, dim=4, ngens=4, approx_model=NoBloating()),
    "GRBX01-MES01" => LGG09(δ=0.0008, template=extdirs, approx_model=Forward()),
    "GRBX01-MES01-discrete" => LGG09(δ=0.0001, template=extdirs, cache=true, approx_model=NoBloating()),
    "PLAD01-BND30" => LGG09(δ=0.03, template=dirs, approx_model=Forward(setops=dirs)),
    "PLAD01-BND30-discrete" => LGG09(δ=0.1, template=dirs, approx_model=NoBloating())
    #="SRA01" => BOX(δ=0.04),
    "SRA01-discrete" => BOX(δ=0.1, approx_model=NoBloating())=#
)

caseCheckDict = Dict(
    "BRKDC01" => "embrake",
    "BRKDC01-discrete" => "embrake",
    "GRBX01-MES01" => "gearbox",
    "GRBX01-MES01-discrete" => "gearbox",
    "PLAD01-BND30" => "platoon",
    "PLAD01-BND30-discrete" => "platoon",
    "SRA01" => "spacecraft",
    "SRA01-discrete" => "spacecraft"
)


function JuliaReachTest(name) 
    println("Running benchmark for: ", name, "...")
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 1

    GC.gc()

    alg = algCheckDict[name]

    sol, df = run(alg, name)

    # Get amonut of steps taken
    # println(length(sol))

    # Save timed data 
    filename = "results/JuliaReach/" * name * ".csv"
    if isfile(filename)
        open(filename, "a") do File
            CSV.write(File, df, delim=";", append=true)
        end
    else
        open(filename, "w") do File
            CSV.write(File, df, delim=";", writeheader=true)
        end
    end
end

function run(alg, case)
    modelname = caseCheckDict[case]
    if modelname == "embrake" 
        sol, df = runEMBrake(alg, case)
    elseif modelname == "gearbox"
        sol, df = runGearbox(alg, case)
    elseif modelname == "platoon"
        sol, df = runPlatoon(alg, case)
    elseif modelname == "spacecraft"
        sol, df = runSpacecraft(alg, case)
    end

    return sol, df
end

function runEMBrake(alg, case)
    prob_no_pv_no_jit = embrake_no_pv(ζ=0.)
    
    sol, df = analyze_embrake(prob_no_pv_no_jit, alg, case, true, true)
    
    return sol, df
end

function runGearbox(alg, case)
    X0_GRBX01h = Hyperrectangle(low=[0, 0, -0.0168, 0.0029, 0, 1],
                            high=[0, 0, -0.0166, 0.0031, 0, 1])

    prob_GRBX01 = gearbox_homog(X0=X0_GRBX01h)

    sol, df = analyze_gearbox(prob_GRBX01, alg, case, true)

    return sol, df
end

function runPlatoon(alg, case)
    prob_PLAD01_BND30 = platoon(; deterministic_switching=true)
    cmethod = LazyClustering(1)
    dirs = CustomDirections(collect(OctDirections{Float64,Vector{Float64}}(10)))
    imethod = TemplateHullIntersection(dirs)
    dmin = -30.0

    sol, df = analyze_platoon(prob_PLAD01_BND30, alg, cmethod, imethod, dmin, case, true)

    return sol, df
end

function runSpacecraft(alg, case)


    return sol, df
end

for case in cases
    JuliaReachTest(case)
end