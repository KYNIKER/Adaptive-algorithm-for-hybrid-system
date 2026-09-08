include("../models/spacecraft.jl")
include("../ReACTed.jl")

### --- MODEL PARAMETERS --- ###
timeConstraintList = []
reduceOrder = 5
maxOrder = 5
dirs = [1,2]
δ⁻ = 0.04
δ⁺ = δ⁻ * 2^10
clustering = false

### --- LOAD MODEL --- ###

sys, initialState, X0, T = loadSpacecraft(abort_time=120.)
n = length(X0.center)

### --- RUN --- ###

res = []
res = ReACTed(sys, initialState, [0., T], X0, Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), dirs, sys.globalConstraints, δ⁻, δ⁺, ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder, reduceOrder, clustering, timeConstraintList)

### --- PLOTTING --- ###

fig = Plots.plot(xlabel="x" , ylabel="y ", ε=1e-6)
cpallete = palette(:roma, length(res))
global i = 1
for (x, y) in res 
    for (r, t) in x 
        d1 = [r[1], -r[2]]
        d2 = [r[3], -r[4]]
        sen = true
        if sen
            Plots.plot!(Shape([d1[1], d1[2], d1[2], d1[1]], [d2[1], d2[1], d2[2], d2[2]]), c=cpallete[i], leg=false, linealpha=0)
            sen = false
        else
            Plots.plot!(Shape([d1[1], d1[2], d1[2], d1[1]], [d2[1], d2[1], d2[2], d2[2]]), c=cpallete[i], leg=false, linealpha=0)
        end
    end
    global i += 1
end

display(fig)
savefig(fig, "plots/Spacecraft2.pdf")