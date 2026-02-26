using Plots, LazySets, LinearAlgebra, BenchmarkTools, CSV, DataFrames, Expokit #, ReachabilityAnalysis

include("Utilities.jl")
include("ReACTed.jl")
include("models/gearbox.jl")
include("models/bouncingBall.jl")

δ⁺ = 10^-2
δ⁻ = δ⁺ / 2^2

sys, initialState, X0, T = loadBouncingBall()
n = length(X0.center)
res = ReACTed(sys, initialState, [0., T], X0, Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), sys.globalConstraints, δ⁻, δ⁺)

#println(res) #
const fig = plot(ε=1e-5)
#plot!(LazySets.HalfSpace(sparsevec([1], [1.], 2), 0.01))
#plot!(LazySets.HalfSpace(sparsevec([1], [-1.], 2), 0.0))
#plot!(LazySets.HalfSpace(sparsevec([2], [1.], 2), 0.0))
cpallete = palette(:tab10, length(res))
i = 1
for (x, y) in res

    println(y)
    for (r, t) in x
        #box = box_approximation(r)
        v = vertices_list(r)
        dimCoords = getindex.(v, 1)
        maxcor = maximum(dimCoords)
        mincor = minimum(dimCoords)
        plot!(Shape([t[1], t[2], t[2], t[1]], [mincor, mincor, maxcor, maxcor]), c=cpallete[i], lab="", alpha=0.1)
        #plot!(r, c=cpallete[i], alpha=0.2)
    end

    h, t = x[end]
    #box = box_approximation(h)
    v = vertices_list(h)
    dimCoords = getindex.(v, 1)
    maxcor = maximum(dimCoords)
    mincor = minimum(dimCoords)
    #plot!(Shape([t[1], t[2], t[2], t[1]], [mincor, mincor, maxcor, maxcor]), c=cpallete[i], lab=i, alpha=0.2)
    h, t = x[1]
    v = vertices_list(h)
    dimCoords = getindex.(v, 1)
    maxcor = maximum(dimCoords)
    mincor = minimum(dimCoords)
    #plot!(Shape([t[1], t[2], t[2], t[1]], [mincor, mincor, maxcor, maxcor]), c=cpallete[i], lab="", alpha=0.8)
    global i += 1

end

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
savefig(fig, joinpath("results/", "Sandbox.png"))

