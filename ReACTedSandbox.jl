using Plots, LazySets, LinearAlgebra, BenchmarkTools, CSV, DataFrames, Expokit #, ReachabilityAnalysis

include("Utilities.jl")
include("ReACTed.jl")
include("models/gearbox.jl")

δ⁻ = 10^-3
δ⁺ = 10^-3

sys = loadGearBox()
n = length(sys.initialState.center)
res = ReACTed(sys, [0., 0.01], convert(Zonotope, sys.initialState), Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), [HalfSpace(sparsevec([5], [1.0], 6), 20.0)], δ⁻, δ⁺)
