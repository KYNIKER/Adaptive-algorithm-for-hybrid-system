using Plots

include("../models/platoon.jl")
include("../ReACTed.jl")

### --- MODEL PARAMETERS --- ###
timeConstraintList = []
reduceOrder = 5
maxOrder = 5
dirs = [1]
δ⁻ = 0.03
δ⁺ = δ⁻ * 2^2
clustering = false

### --- LOAD MODEL --- ###

sys, initialState, X0, T = loadPlatoon()
n = length(X0.center)

### --- RUN --- ###

res = []
res = ReACTed(sys, initialState, [0., T], X0, Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), dirs, sys.globalConstraints, δ⁻, δ⁺, ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder, reduceOrder, clustering, timeConstraintList)

### --- PLOTTING --- ###

fig = Plots.plot(xlabel="time(s)" , ylabel="d₁", ε=1e-6)
cpallete = palette(:roma, length(res))
global i = 1
for (x, y) in res 
    sen = true
    for (d, t) in x 

        if sen
            Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [d[1], d[1], -d[2], -d[2]]), c=cpallete[i], leg=false, linealpha=0)
            sen = false
        else
            Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [d[1], d[1], -d[2], -d[2]]), c=cpallete[i], leg=false, linealpha=0)
        end
    end
    global i += 1
end

display(fig)
savefig(fig, "plots/Platoon.pdf")