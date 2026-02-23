using Plots, LazySets, LinearAlgebra, BenchmarkTools, CSV, DataFrames, Expokit #, ReachabilityAnalysis

include("Utilities.jl")
include("ReACTed.jl")
include("models/gearbox.jl")

δ⁺ = 10^-3
δ⁻ = δ⁺ / 2^1

sys = loadGearBox()
n = length(sys.initialState.center)
res = ReACTed(sys, [0., 0.2], convert(Zonotope, sys.initialState), Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), [LazySets.HalfSpace(sparsevec([5], [1.0], 6), 20.0)], δ⁻, δ⁺)

#println(res) #
for (x, y) in res
    println(ρ(Vector(sparsevec([5], [1.0], 6)), x), " :", y)

end
#const fig = plot()
#plot!(fig, res[1], vars=(1, 2), ε=1e-5)
# plot!(fig, sol_GRBX01, vars=(3, 4), ε=1e-5,
#       color=:blue, alpha=0.5, lw=1.0, linecolor=:blue,
#       tickfont=font(30, "Times"), guidefontsize=45,
#       xlab=L"x_3", ylab=L"x_4",
#       xtick=([-0.016, -0.009, -0.002],
#              [L"-0.016", L"-0.009", L"-0.002"]),
#       ytick=([-0.008, -0.004, 0.0, 0.004],
#              [L"-0.008", L"-0.004", L"0", L"0.004"]),
#       xlims=(-0.017, -0.0015), ylims=(-0.008, 0.004),
#       bottom_margin=-5mm, left_margin=-1mm, right_margin=10mm, top_margin=3mm,
#       size=(1000, 1000))
#savefig(fig, joinpath("results/", "Sandbox.png"))

