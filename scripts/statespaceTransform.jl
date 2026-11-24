using LazySets, LinearAlgebra, Plots, ReachabilityAnalysis, JLD2, FileIO, BenchmarkTools, ProfileView
using Plots.PlotMeasures

gr()
include(joinpath(@__DIR__, "..", "src", "ReACTShields.jl"))
using .ReACTShields
include(joinpath(@__DIR__, "..", "src", "core", "discretize.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "utilities.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "partition.jl"))
include(joinpath(@__DIR__, "..", "src", "shields", "guards.jl"))
include(joinpath(@__DIR__, "..", "src", "shields", "propagate.jl"))
include(joinpath(@__DIR__, "..", "models", "bouncingBall.jl"))


Zs = Zonotope([1.0, 0.0], [0.0 0.0; 0.0 0.0])
Hs = Zonotope([1.0, 0.0], [0.0 0.0; 4.0 0.0])
bs = Zonotope([1.0, 0.0], zeros(2, 1))
A = [0.6 1.3; -1.1 0.3]
HT = diagm([1.0, 1.0]) - (2\norm([1.0, 0.0])^2)*([1.0, 0.0] * conj(transpose([1.0, 0.0])))
u = A * [1.0, 0.0] - (norm(A * [1.0, 0.0])/norm([1.0, 0.0]))*[1.0, 0.0]
TT = diagm([1.0, 1.0]) - (2*u*transpose(u))/(transpose(u)*u)
@show HT
@show u, TT * u, TT * A * [1.0, 0.0], TT * TT
d = 0.2
fig = Plots.plot()

H1 = LazySets.Hyperplane([1.0, 0.0], 1.)
H2 = linear_map(exp(0.05*A), H1)

#plot!(fig, H1)
#plot!(fig, H2)
plot!(fig, Zs)
#plot!(fig, transpose(exp(d * A)) * Hs)
J = intersection(exp(d/20 * A) * Hs, Hs)
plot!(fig, J)
plot!(fig, intersection(Zonotope(transpose(exp(-1 * A)) * [1.0, 0.0], [[1.0 0.0; 0.0 1.0] * (transpose(exp(-1 * A))) * [2.0, 0.0]]), Hs))
I = Intersection(H1, H2)#intersection(Zonotope(transpose(exp(-d/10 * A)) * [1.0, 0.0], [[1.0 0.0; 0.0 1.0] * (transpose(exp(-d/10 * A))) * [2.0, 0.0]]), Hs)
#plot!(fig, an_element(I))
#@show I
#plot!(fig, exp(-0*A) * intersection(transpose(exp(-d * A)) * Hs, Hs))
#plot!(fig, transpose(exp(-d * A)) * bs)
#@show dot(transpose(exp(-d * A)) * [1.0, 0.0], exp(d * A) * [0.0, 4.0])
#@show dot(transpose(exp(-d * A)) * [1.0, 0.0], [0.0 -1.0; 1.0 0.0] * (transpose(exp(-d * A))) * [1.0, 0.0])
#@show ρ([1.0, 0.0], exp(-0*A) * intersection(Zonotope(transpose(exp(-d * A)) * [1.0, 0.0], [[0.0 -1.0; 1.0 0.0] * (transpose(exp(-d * A))) * [1.0, 0.0]]), Hs))
plot!(fig, A * TT * Zs, color=:blue, alpha=0.1, fa=0.0)
for i in 0.5:0.5:8
    tu = exp(i*d*A) * [1.0, 0.0] - (norm(exp(i*d*A) * [1.0, 0.0])/norm([1.0, 0.0]))*[1.0, 0.0]
    tTT = diagm([1.0, 1.0]) - (2*tu*transpose(tu))/(transpose(tu)*tu)
    @show tTT
    plot!(fig, exp(i*d*A) * Zs, color=:red, alpha=0.1)
    plot!(fig, tTT * exp(i * d * A) * Zs, color=:blue, alpha=0.1, fa=0.0)
    plot!(fig, tTT * exp(i * -d * A) * exp(i * d * A) * tTT * Zs, color=:blue, alpha=0.1, fa=0.0)
    plot!(fig, exp(i*d*A) * J, color=:yellow, alpha=0.1)
    #   @show ρ([1.0, 0.0], exp(i*d*A) * J)
end

xlims!(-4.0, 10.0)
ylims!(-2.0, 2.0)
display(fig)