using Plots

include("../models/embrake/embrake.jl")
include("../models/embrake/solve_embrake.jl")
include("../ReACTed.jl")

### --- MODEL PARAMETERS --- ###
reduceOrder = 5
maxOrder = 5
dirs = [2]
δ⁻ = 2*10^-7
δ⁺ = δ⁻ * 2^8

### --- LOAD MODEL --- ###

sys, initialState, X0, T = loadembrake()
n = length(X0.center)

### --- RUN --- ###

res = []
res = solve_embrake(sys, initialState, X0, T, δ⁺, δ⁻, maxOrder, reduceOrder, dirs, true)

### --- PLOTTING --- ###

fig = Plots.plot(xlabel="time(s)" , ylabel="x", ε=1e-6)
cpallete = palette(:roma, length(res))
global i = 1
for (x, y) in res 
    sen = true
    for (d, t) in x 

        if sen
            Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [d[1], d[1], -d[2], -d[2]]), c=cpallete[i], leg=false, linealpha=1)
            sen = false
        else
            Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [d[1], d[1], -d[2], -d[2]]), c=cpallete[i], leg=false, linealpha=1)
        end
    end
    global i += 1
end

display(fig)
savefig(fig, "plots/Brake.pdf")