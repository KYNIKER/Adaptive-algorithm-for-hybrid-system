using LazySets, ReachabilityAnalysis, Plots, LinearAlgebra
using Plots.Measures

include("Discretize.jl")

c = Float64.(rand((0:10), 2))

G = randn(2, 2)


Z = Zonotope(c, G)
Z1 = Zonotope(c, diagm(zeros(2)))
Z2 = Zonotope(zeros(2), G)

cU = randn(2)#Float64.(rand((0:5), 2))
GU = randn(2, 2)

U = Zonotope(cU, GU) #Zonotope(zeros(2), diagm(zeros(2)))#
#U = Zonotope(zeros(2), diagm(zeros(2)))

A = randn(2, 2)

d = 0.25
#dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
isInvA = false #isinvertible(A)
Φ = copy(exp(d * A))
A_abs = abs.(A)
Φcache = sum(A) == abs(sum(A)) ? Φ : nothing
P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, d, ReachabilityAnalysis.Exponentiation.BaseExp, false, Φcache)
pis = ReachabilityAnalysis.Exponentiation.Φ₁(A, d, ReachabilityAnalysis.Exponentiation.BaseExp, false, nothing)

plt = plot(dpi=1200, thickness_scaling=1, guidefontsize=25, minorgrid=true,
    legendfont=font(6, "Times"),
    legend_position=:bottomright,
    tickfont=font(8, "Times"),
    xguidefont=font(12, "Times"),
    yguidefont=font(12, "Times"),
    bottom_margin=2mm,
    left_margin=5mm,
    right_margin=5mm,
    top_margin=2mm,
    xlabel="x", ylabel="y")



#=
plot!(plt, symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z)))), c=:lightblue, lab="original bloating")
plot!(plt, symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z1)))), c=:blue, lab="0 generator bloating")
plot!(plt, symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z2)))), c=:darkblue, lab="centered bloating")
plot!(plt, minkowski_sum(symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z1)))), symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z2))))), c=:black, fillstyle=://, lab="Minkowski sum of 0 generator and centered bloating")
plot!(plt, CH(minkowski_sum(Z1, Z2), minkowski_sum(linear_map(Φ, minkowski_sum(Z1, Z2)), minkowski_sum(symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z1)))), symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z2))))))), c=:black, fillstyle=:-, lab="using zonotope decomposition")
plot!(plt, CH(Z, minkowski_sum(linear_map(Φ, Z), symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z)))))), c=:orange, lab="original method")

plot!(plt, Z, c=:green, lab="original zonotope")

plot!(plt, Z1, c=:blue, lab="0 generator zonotope")
plot!(plt, Z2, c=:yellow, lab="centered zonotope")
=#
phiDict, tPhiDict, inputDict = phi_input_dict(A, U, d, d)
phiDictp, tPhiDictp, inputDictp = phi_input_dict(A, U, d/4, d)

#@show inputDict
#@show inputDictp

isInvA = false#isinvertible(system.flowMatrix)
#Φ = copy(phiDict[d])
A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
Φcache = nothing
P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, d, ReachabilityAnalysis.Exponentiation.BaseExp, isInvA, Φcache)
P2A_absp = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, d/4, ReachabilityAnalysis.Exponentiation.BaseExp, false, nothing)

generatorDict = ReACT_discretize_decomposed_generators(Z1, d/4, d, A, P2A_absp, phiDictp, U, inputDictp)
discDict, finaldims = ReACT_discretize_combine_with_offsets(Z1, d/4, d, A, P2A_absp, phiDictp, U, inputDictp, generatorDict)
#println("new dict done")
#@show Z
oldDiscDict = newReACTDiscretizePlus(Z1, d/4, d, A, P2A_absp, phiDictp, inputDictp)
#plot!(plt, oldDiscDict[d], fa=0.1, lab="current")
#plot!(plt, oldDiscDict[d/2], fa=0.1, lab="current", c=:red)
#plot!(plt, oldDiscDict[d/2], c=:red, lab="old implementation d-")
#plot!(plt, linear_map(phiDictp[d/2], oldDiscDict[d/2]), c=:red, lab="old implementation phi * d-")
#@show Z

#=
@show LazySets.order(remove_zero_generators(discDict[d/2]))
@show LazySets.order(remove_zero_generators(oldDiscDict[d/2]))
@show LazySets.order(remove_zero_generators(discDict[d]))
@show LazySets.order(remove_zero_generators(oldDiscDict[d]))
@show discDict
=#
#plot!(plt, discDict[d/2], c=:grey, lab="implementation")
#plot!(plt, discDict[d], c=:grey, lab="implementation", fillstyle=:\)

@show finaldims
stateDict = copy(discDict)
preallocated_center = zeros(2)
preallocated_genmat = zeros(2, finaldims + 6)
preallocated_inter_genmat = zeros(2, finaldims + 6)
ReACT_discretize_state!(preallocated_center, preallocated_genmat, preallocated_inter_genmat, stateDict, c, d/4, d, A, P2A_absp, phiDictp, U, inputDictp, generatorDict)
#plot!(plt, stateDict[d/2], c=:yellow, lab="state dict implementation")
#plot!(plt, remove_zero_generators(stateDict[d]), c=:purple, lab="state dict implementation")


#plot!(plt, CH(oldDiscDict[d/2], linear_map(phiDictp[d/2], oldDiscDict[d/2])), c=:black, fillstyle=:\)
#plot!(plt, overapproximate(CH(oldDiscDict[d/2], linear_map(phiDictp[d/2], oldDiscDict[d/2])), Zonotope), c=:yellow, fillstyle=:+)
#oldDiscDictp = newReACTDiscretizePlus(Z, d, d, A, P2A_abs, phiDict, inputDict)
println("old dict done")
#plot!(plt, oldDiscDictp[d], c=:green, alpha=0.1, lab="old implementation d")

SSS = CH(Z1, minkowski_sum(minkowski_sum(linear_map(phiDictp[d/2], Z1), inputDictp[d/2]), symmetric_interval_hull(linear_map(P2A_absp, symmetric_interval_hull(linear_map(A^2, Z1))))))
SS = CH(SSS, MinkowskiSum(linear_map(phiDictp[d/2], SSS), inputDictp[d/2]))
#@show Z
#plot!(plt, SSS, lc=:white, lab="CH first step")
#plot!(plt, SS, lc=:white, lab="original method with input", fillstyle=:||)
#=
@show oldDiscDict[d]
@show discDict[d]
=#
@show oldDiscDict[d/2]
@show discDict[d/2]
@show stateDict[d/2]
@show genmat(discDict[d]) - genmat(stateDict[d])
@show LazySets.API.issubset(SS, oldDiscDict[d])
@show LazySets.API.issubset(SS, discDict[d])
@show LazySets.API.issubset(SS, stateDict[d])
#@show LazySets.API.issubset(SSS, overapproximate(CH(oldDiscDict[d/2], linear_map(phiDictp[d/2], oldDiscDict[d/2])), Zonotope))
display(plt)
#plot!(plt, symmetric_interval_hull(Z1), c=:grey)
#plot!(plt, symmetric_interval_hull(Z2), c=:grey)

