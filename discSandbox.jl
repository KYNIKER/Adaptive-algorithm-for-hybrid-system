using LazySets, ReachabilityAnalysis, Plots

c = Float64.(rand((0:10), 2))

G = rand(Float64, (2, 5))


Z = Zonotope(c, G)
Z1 = Zonotope(c, diagm(zeros(2)))
Z2 = Zonotope(zeros(2), G)

cU = Float64.(rand((0:5), 2))
GU = rand(Float64, (2, 2))

U = Zonotope(cU, GU)

A = rand(Float64, (2, 2))

d = 0.5
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




plot!(plt, symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z)))), c=:lightblue, lab="original bloating")
plot!(plt, symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z1)))), c=:blue, lab="0 generator bloating")
plot!(plt, symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z2)))), c=:darkblue, lab="centered bloating")
plot!(plt, minkowski_sum(symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z1)))), symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z2))))), c=:black, fillstyle=://, lab="Minkowski sum of 0 generator and centered bloating")
plot!(plt, CH(minkowski_sum(Z1, Z2), minkowski_sum(linear_map(Φ, minkowski_sum(Z1, Z2)), minkowski_sum(symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z1)))), symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z2))))))), c=:black, fillstyle=:-, lab="using zonotope decomposition")
plot!(plt, CH(Z, minkowski_sum(linear_map(Φ, Z), symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A^2, Z)))))), c=:orange, lab="original method")
plot!(plt, Z, c=:green, lab="original zonotope")

plot!(plt, Z1, c=:blue, lab="0 generator zonotope")
plot!(plt, Z2, c=:yellow, lab="centered zonotope")

#plot!(plt, symmetric_interval_hull(Z1), c=:grey)
#plot!(plt, symmetric_interval_hull(Z2), c=:grey)

