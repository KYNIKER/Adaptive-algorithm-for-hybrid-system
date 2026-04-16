using Plots, LazySets, LinearAlgebra, BenchmarkTools, CSV, DataFrames, Expokit, CDDLib #, ReachabilityAnalysis

include("Utilities.jl")
include("ReACTed.jl")
include("models/girardExample.jl")
include("models/gearbox.jl")
include("models/platoon.jl")
include("models/bouncingBall.jl")
include("models/simpleModel.jl")
include("models/powerTrain.jl")
include("models/spacecraft.jl")

#using Cthulhu, ProfileView

saveResult = true

δ⁺ = 10^-3 * 2
δ⁺ = 0.02
δ⁻ = δ⁺ / 2^0

#sys, initialState, X0, T = loadPlatoon()
sys, initialState, X0, T = loadBouncingBall()
#sys, initialState, X0, T = loadPowertrain(θ=3, homog = true)
# sys, initialState, X0, T = loadSpacecraft()
# sys, initialState, X0, T = loadSpacecraft(abort_time = 120.)
T = 3.5

n = length(X0.center)

#@ProfileView.profview ReACTed(sys, initialState, [0., T], X0, Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), sys.globalConstraints, δ⁻, δ⁺, ReachabilityAnalysis.Exponentiation.BaseExp, 5, 5, saveResult)
res = []

res = ReACTed(sys, initialState, [0., T], X0, Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), sys.globalConstraints, δ⁻, δ⁺, ReachabilityAnalysis.Exponentiation.BaseExp, 5, 5, saveResult)

plotProjectedFlowpipeLazy(res, 0, 1, n, joinpath("results/", "SandboxDim2.png"))

# plotProjectedFlowpipe(res, 0, 4, joinpath("results/", "SandboxDim24.png"))

# plotProjectedFlowpipe(res, 0, 7, joinpath("results/", "SandboxDim27.png"))
res = []
# plotProjectedFlowpipe(res, 0, 1, joinpath("results/", "SandboxDim1.png"))
#plotProjectedFlowpipe(res, 0, 3, joinpath("results/", "SandboxDimVx.png"))
# plotProjectedFlowpipe(res, 0, 4, joinpath("results/", "SandboxDimVy.png"))
x = 1  # x position (negative!)
y = 2  # y position (negative!)
vx = 3  # x velocity
vy = 4  # y velocity
velocity = 0.055 * 60.0     # meters per minute
cx = velocity * cos(π / 8)  # x-coordinate of the octagon's first (ENE) corner
cy = velocity * sin(π / 8)  # y-coordinate of the octagon's first (ENE) corner
n = 5


#println(loc2Constraint)

#plotProjectedFlowpipe(res, 1, 2, joinpath("results/", "SandboxDim2Constraint.png"), loc2Constraint, 0.1)

# #println(res) #
# const fig = Plots.plot(ε=1e-6)
# #plot!(LazySets.HalfSpace(sparsevec([1], [1.], 2), 0.01))
# #plot!(LazySets.HalfSpace(sparsevec([1], [-1.], 2), 0.0))
# #plot!(LazySets.HalfSpace(sparsevec([2], [1.], 2), 0.0))
# cpallete = palette(:tab10, length(res))
# i = 1

# for (x, y) in res

#     println(y)
#     for (r, t) in x
#         #box = box_approximation(r)
#         v = vertices_list(r)
#         dimCoords = getindex.(v, 1)
#         maxcor = maximum(dimCoords)
#         mincor = minimum(dimCoords)
#         Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [mincor, mincor, maxcor, maxcor]), c=cpallete[i], lab="", alpha=0.1)
#         #plot!(r, c=cpallete[i], alpha=0.2)
#     end

#     h, t = x[end]
#     #box = box_approximation(h)
#     v = vertices_list(h)
#     dimCoords = getindex.(v, 1)
#     maxcor = maximum(dimCoords)
#     mincor = minimum(dimCoords)
#     #plot!(Shape([t[1], t[2], t[2], t[1]], [mincor, mincor, maxcor, maxcor]), c=cpallete[i], lab=i, alpha=0.2)
#     h, t = x[1]
#     v = vertices_list(h)
#     dimCoords = getindex.(v, 1)
#     maxcor = maximum(dimCoords)
#     mincor = minimum(dimCoords)
#     #plot!(Shape([t[1], t[2], t[2], t[1]], [mincor, mincor, maxcor, maxcor]), c=cpallete[i], lab="", alpha=0.8)
#     global i += 1

# end
# savefig(fig, joinpath("results/", "SandboxDim1.png"))


# # const fig2 = Plots.plot(xlabel="Height", ylabel="Velocity")
# # cpallete = palette(:tab10, length(res))
# # i = 1
# # for (x, y) in res
# #     println(y)
# #     for (r, t) in x
# #         box = box_approximation(r)
# #         v = vertices_list(r)
# #         dimCoords1 = getindex.(v, 1)
# #         dimCoords2 = getindex.(v, 2)
# #         maxcor1 = maximum(dimCoords1)
# #         mincor1 = minimum(dimCoords1)
# #         maxcor2 = maximum(dimCoords2)
# #         mincor2 = minimum(dimCoords2)
# #         Plots.plot!(Shape([mincor1, maxcor1, maxcor1, mincor1], [mincor2, mincor2, maxcor2, maxcor2]), c=cpallete[i], lab="")
# #     end
# #     global i += 1
# # end
# # savefig(fig2, joinpath("results/", "SandboxPlot2Dim.png"))






# #plot!(fig, res[1], vars=(1, 2), ε=1e-5)
# # plot!(fig, sol_GRBX01, vars=(3, 4), ε=1e-5,
# #       color=:blue, alpha=0.5, lw=1.0, linecolor=:blue,
# #       tickfont=font(30, "Times"), guidefontsize=45,
# #       xlab=L"x_3", ylab=L"x_4",
# #       xtick=([-0.016, -0.009, -0.002],
# #              [L"-0.016", L"-0.009", L"-0.002"]),
# #       ytick=([-0.008, -0.004, 0.0, 0.004],
# #              [L"-0.008", L"-0.004", L"0", L"0.004"]),
# #       xlims=(-0.017, -0.0015), ylims=(-0.008, 0.004),
# #       bottom_margin=-5mm, left_margin=-1mm, right_margin=10mm, top_margin=3mm,
# #       size=(1000, 1000))
# #savefig(fig, joinpath("results/", "Sandbox.png"))

