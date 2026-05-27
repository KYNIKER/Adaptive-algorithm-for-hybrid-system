using ReachabilityAnalysis, LazySets, BenchmarkTools, CSV, DataFrames, Plots

include("embrake/embrake.jl")
include("embrake/embrake_benchmark.jl")
include("gearbox/gearbox.jl")
include("gearbox/gearbox_benchmark.jl")
include("platoon/platoon.jl")
include("platoon/platoon_benchmark.jl")
include("spacecraft/spacecraft_benchmark.jl")
include("spacecraft/spacecraft.jl")

# template directions

const extdirs = template_directions(6)
const dirs = CustomDirections(collect(OctDirections{Float64,Vector{Float64}}(10)))

cases = [
    "BRKDC01",
    "BRKDC01-discrete",
    "GRBX01-MES01",
    "GRBX01-MES01-discrete",
    "PLAD01-BND30",
    "PLAD01-BND30-discrete",
    "SRA01",
    "SRA01-discrete"
]

algCheckDict = Dict(
    "BRKDC01" =>  GLGM06(δ=2e-7, max_order=1, static=true, dim=4, ngens=4, approx_model=Forward()),
    "BRKDC01-discrete" => GLGM06(δ=1e-8, max_order=1, static=true, dim=4, ngens=4, approx_model=NoBloating()),
    "GRBX01-MES01" => LGG09(δ=0.0008, template=extdirs, approx_model=Forward()),
    "GRBX01-MES01-discrete" => LGG09(δ=0.0001, template=extdirs, cache=true, approx_model=NoBloating()),
    "PLAD01-BND30" => LGG09(δ=0.03, template=dirs, approx_model=Forward(setops=dirs)),
    "PLAD01-BND30-discrete" => LGG09(δ=0.1, template=dirs, approx_model=NoBloating()),
    "SRA01" => BOX(δ=0.04),
    "SRA01-discrete" => BOX(δ=0.1, approx_model=NoBloating())
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

    sol, df = run_algs(alg, name, true)

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

function run_algs(alg, case, doPlot)
    modelname = caseCheckDict[case]
    if modelname == "embrake" 
        sol, df = runEMBrake(alg, case, doPlot)
    elseif modelname == "gearbox"
        sol, df = runGearbox(alg, case, doPlot)
    elseif modelname == "platoon"
        sol, df = runPlatoon(alg, case, doPlot)
    elseif modelname == "spacecraft"
        sol, df = runSpacecraft(alg, case, doPlot)
    end

    return sol, df
end

function runEMBrake(alg, case, doPlot)
    prob_no_pv_no_jit = embrake_no_pv(ζ=0.)
    
    sol, df = analyze_embrake(prob_no_pv_no_jit, alg, case, true, true)
    
    if doPlot
        TARGET_FOLDER = isdefined(Main, :TARGET_FOLDER) ? Main.TARGET_FOLDER : @__DIR__

        # modify tolerance
        LazySets.set_ztol(Float64, 1e-12)

        # compute interval approximation in dimension 2 (x)
        polys = Vector{VPolygon{Float64, Vector{Float64}}}()
        for fp in sol
            for (j, Rj) in enumerate(fp)
                sfpos, sfneg = ρ(eₓ, Rj), ρ(-eₓ, Rj)
                dt = tspan(fp, j)
                ti, tf = inf(dt), sup(dt)
                p = VPolygon([[ti, sfpos], [ti, -sfneg], [tf, sfpos], [tf, -sfneg]])
                push!(polys, p)
            end
        end

        sol_no_pv_no_jit = nothing
        GC.gc()

        # ignore the first 1500 sets because the plot recipe does not like them
        fig = plot()
        plot!(fig, polys[1:500:end], color=:blue, lw=1.0, linecolor=:blue,
              tickfont=font(30, "Times"), guidefontsize=45,
              xlab=L"t", ylab=L"x",
              xtick=[0.025, 0.05, 0.075, 0.1], ytick=[0.0, 0.01, 0.02, 0.03, 0.04, 0.05],
              xlims=(0.0, 0.1), ylims=(0.0, 0.05),
              bottom_margin=-6mm, left_margin=-3mm, right_margin=12mm, top_margin=3mm,
              size=(1000, 1000))
        hline!(fig, [x0], lc=:red, ls=:dash, lw=2, lab="")

        savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP25-JuliaReach-$case.pdf"))

        # reset tolerance
        LazySets.set_tolerance(Float64)

end
    
    return sol, df
end

function runGearbox(alg, case, doPlot)
    X0_GRBX01h = Hyperrectangle(low=[0, 0, -0.0168, 0.0029, 0, 1],
                            high=[0, 0, -0.0166, 0.0031, 0, 1])

    prob_GRBX01 = gearbox_homog(X0=X0_GRBX01h)

    sol, df = analyze_gearbox(prob_GRBX01, alg, case, true)

    if doPlot
        TARGET_FOLDER = isdefined(Main, :TARGET_FOLDER) ? Main.TARGET_FOLDER : @__DIR__

        # modify tolerance
        LazySets.set_ztol(Float64, 1e-8)

        fig = plot()
        plot!(fig, sol, vars=(3, 4), ε=1e-5,
              color=:blue, alpha=0.5, lw=1.0, linecolor=:blue,
              tickfont=font(30, "Times"), guidefontsize=45,
              xlab=L"x_3", ylab=L"x_4",
              xtick=([-0.016, -0.009, -0.002],
                     [L"-0.016", L"-0.009", L"-0.002"]),
              ytick=([-0.008, -0.004, 0.0, 0.004],
                     [L"-0.008", L"-0.004", L"0", L"0.004"]),
              xlims=(-0.017, -0.0015), ylims=(-0.008, 0.004),
              bottom_margin=-5mm, left_margin=-1mm, right_margin=10mm, top_margin=3mm,
              size=(1000, 1000))
        savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP25-JuliaReach-$case.pdf"))

        # reset tolerance
        LazySets.set_tolerance(Float64)
    end

    return sol, df
end

function runPlatoon(alg, case, doPlot)
    prob_PLAD01_BND30 = platoon(; deterministic_switching=true)
    cmethod = LazyClustering(1)
    dirs = CustomDirections(collect(OctDirections{Float64,Vector{Float64}}(10)))
    imethod = TemplateHullIntersection(dirs)
    dmin = -30.0

    sol, df = analyze_platoon(prob_PLAD01_BND30, alg, cmethod, imethod, dmin, case, true)

    if doPlot
        TARGET_FOLDER = isdefined(Main, :TARGET_FOLDER) ? Main.TARGET_FOLDER : @__DIR__
        fig = plot()
        plot!(fig, sol, vars=(0, 1),
              linecolor=:blue, color=:blue, alpha=0.8, lw=1.0,
              tickfont=font(30, "Times"), guidefontsize=45,
              xlab=L"t", ylab=L"x_{1}",
              xtick=[0, 5, 10, 15, 20.], ytick=[-30, -20, -10, 0],
              xlims=(0., 20.), ylims=(-31, 7),
              bottom_margin=-6mm, left_margin=-2mm, right_margin=4mm, top_margin=0mm,
              size=(1000, 1000))
        hline!(fig, [-30.0], lc=:red, ls=:dash, lw=2, lab="")

        savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP25-JuliaReach-$case.pdf"))

    end

    return sol, df
end

function runSpacecraft(alg, case, doPlot)
    prob_SRA01 = spacecraft(abort_time=120.0)
    cmethod = LazyClustering(3)

    sol, df = analyze_spacecraft(prob_SRA01, alg, cmethod, true, case, true)

    if doPlot 
        TARGET_FOLDER = isdefined(Main, :TARGET_FOLDER) ? Main.TARGET_FOLDER : @__DIR__
        # SRA01
        fig = plot(legend=:bottomright, tickfont=font(30, "Times"), guidefontsize=45,
            xlab=L"x", ylab=L"y", lw=0.0,
            xtick=[-750, -500, -250, 0, 250.0], ytick=[-400, -300, -200, -100, 0.0],
            xlims=(-1000.0, 400.0), ylims=(-450.0, 0.0),
            bottom_margin=-6mm, left_margin=-3mm, right_margin=1mm, top_margin=3mm,
            size=(1000, 1000))
        for idx in findall(L -> L == 1, location.(sol))
            plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:lightgreen, lc=:lightgreen, alpha=1.0)
        end
        for idx in findall(L -> L == 2, location.(sol))
            plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:red, lc=:red, alpha=1.0)
        end
        for idx in findall(L -> L == 3, location.(sol))
            plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:cyan, lc=:cyan, alpha=1.0)
        end
        savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP25-JuliaReach-$case.pdf"))
    end
    return sol, df
end

for case in cases
    JuliaReachTest(case)
end