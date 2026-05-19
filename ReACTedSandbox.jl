using Plots, LazySets, LinearAlgebra, BenchmarkTools, CSV, DataFrames, Expokit, CDDLib, Profile, PProf, JuMP, HiGHS #, ProfileView #, ReachabilityAnalysis

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
### sys, initialState, X0, T = loadPowertrain(θ=3, homog=true)
### sys, initialState, X0, T = loadSpacecraft()

dirs = [1] # Which direction to plot. Empty means no plotting

#δ⁺ = 10^-3 * 2
δ⁺ = 0.03
δ⁻ = δ⁺ / 2^0

#sys, initialState, X0, T = loadBouncingBall()
sys, initialState, X0, T = loadPlatoon()
#sys, initialState, X0, T = loadSpacecraft(abort_time=120.)
#sys, initialState, X0, T = loadGearBox()
#T = 115.5
#T = 19.2
#T = 12.0
n = length(X0.center)

#@ProfileView.profview ReACTed(sys, initialState, [0., T], X0, Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), sys.globalConstraints, δ⁻, δ⁺, ReachabilityAnalysis.Exponentiation.BaseExp, 5, 5, saveResult)
res = []
#LazySets.deactivate_assertions()
#LazySets.set_tolerance(0.00001)

res = ReACTed(sys, initialState, [0., T], X0, Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), dirs, sys.globalConstraints, δ⁻, δ⁺, ReachabilityAnalysis.Exponentiation.BaseExp, 5, 5; mustSemantics=true)
#@time ReACTed(sys, initialState, [0., T], X0, Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), dirs, sys.globalConstraints, δ⁻, δ⁺, ReachabilityAnalysis.Exponentiation.BaseExp, 5, 5)
#println(res[end])
println("Begin plotting")
plotProjectedFlowpipeLazy(res, dirs, n, joinpath("results/", "Platoon.png"); verbose=true)#; xlim=(108.0, 110.0), ylim=(2.0, 3.0)

#plotProjectedFlowpipeLazy(res, dirs, n, joinpath("results/", "SandboxDim2.png"); xlim=(0.0, T), ylim=(-30.0, 10.0))#; xlim=(108.0, 110.0), ylim=(2.0, 3.0)
