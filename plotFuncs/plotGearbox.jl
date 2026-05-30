include("../models/gearbox.jl")
include("../ReACTed.jl")

### --- MODEL PARAMETERS --- ###
timeConstraintList = [(1, 0.2)]
reduceOrder = 5
maxOrder = 5
dirs = [3,4]
δ⁻ = 0.0008
δ⁺ = δ⁻ * 2^4
clustering = false

### --- LOAD MODEL --- ###

sys, initialState, X0, T = loadGearBox()
n = length(X0.center)

### --- RUN --- ###

res = []
res = ReACTed(sys, initialState, [0., T], X0, Zonotope(zeros(Float64, n), zeros(Float64, n, 1)), dirs, sys.globalConstraints, δ⁻, δ⁺, ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder, reduceOrder, clustering, timeConstraintList)

### --- PLOTTING --- ###

fig = Plots.plot(xlabel="x₃" , ylabel="x₄", ε=1e-6)
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
savefig(fig, "plots/Gearbox.pdf")