using Plots, LazySets, LinearAlgebra, BenchmarkTools, CSV, DataFrames, Expokit #, ReachabilityAnalysis

include("Utilities.jl")
include("ReACTed.jl")
include("models/gearbox.jl")

δ⁺ = 10^-2
δ⁻ = δ⁺ / 2^3

sys = loadGearBox()
n = length(sys.initialState.center)
res = ReACTed(sys, [0., 0.05], convert(Zonotope, sys.initialState), Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), [HalfSpace(sparsevec([5], [1.0], 6), 20.0)], δ⁻, δ⁺)

println(res) #println(map(x -> ρ(Vector(sparsevec([5], [1.0], 6)), x), res))
