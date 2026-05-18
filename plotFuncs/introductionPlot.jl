include("plotFuncHelper.jl") # This also loads everything relevant
ENV["GKSwstype"] = "100"



sys, initialState, X0, T = loadBouncingBall()
T = 7
timeConstraintList = []
clustering = true
n = length(X0.center)
reduceOrder = 5
maxOrder = 5

dirs = [1]

δ⁻ = 0.01
δ⁺ = δ⁻ * 2^4


res1 = ReACTed(sys, initialState, [0., T], X0, Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), dirs, sys.globalConstraints, δ⁻, δ⁺, ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder, reduceOrder, clustering, timeConstraintList)
shapes1, _, _, maxY1, minY1 = getShapesForPlot(res1, dirs)

δ⁺ = δ⁻

res2 = ReACTed(sys, initialState, [0., T], X0, Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), dirs, sys.globalConstraints, δ⁻, δ⁺, ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder, reduceOrder, clustering, timeConstraintList)
shapes2, _, _, maxY2, minY2 = getShapesForPlot(res2, dirs)

maxVal = max(maxY1, maxY2)
minVal = min(minY1, minY2)


palette = Plots.palette(:fes10)
alp = 0.7


p = plot(dpi=1200, thickness_scaling=1, guidefontsize=25, minorgrid=true,
    legendfont=font(12, "Times"),
    legend_position=:topright,
    tickfont=font(8, "Times"),
    xguidefont=font(12, "Times"),
    yguidefont=font(12, "Times"),
    xtick=([0, 8], [L"0", L"T"]),
    ytick=([], []),
    bottom_margin=2mm,
    left_margin=5mm,
    right_margin=5mm,
    top_margin=2mm,
    ylims=(minVal, maxVal), xlims=(0, maximum(T)), xlabel=L"Time", ylabel=L"x")



for i in eachindex(shapes1)
    if i == 1
        plot!(p, shapes1[i], vars=(1, 0), c=palette[9], alpha=0.7, lw=0.05,
            label="Our approach")
    else
        plot!(p, shapes1[i], vars=(1, 0), c=palette[9], alpha=0.7, lw=0.05,
            label="")
    end
end
for i in eachindex(shapes2)
    if i == 1
        plot!(p, shapes2[i], vars=(1, 0), c=palette[6], alpha=1.0, la=0.0, lw=0.05,
            label="Fixed step")
    else
        plot!(p, shapes2[i], vars=(1, 0), c=palette[6], alpha=1.0, la=0.0, lw=0.05,
            label="")
    end
end


savefig(p, "plotResults/" * "Introduction.pdf")
savefig(p, "plotResults/" * "Introduction.png")
display(p)



