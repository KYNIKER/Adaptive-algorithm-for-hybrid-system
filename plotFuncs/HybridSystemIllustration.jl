#example
using LazySets, ReachabilityAnalysis, LinearAlgebra, Polyhedra, Optim, Plots, CDDLib
using LaTeXStrings, Plots.PlotMeasures
using SparseArrays: sparsevec
include("../Utilities.jl")
#include("plotFuncHelper.jl")




#=
input = HPolyhedron([
    LazySets.HalfSpace(sparsevec([1], [1.], 2), 3.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([1], [-1.], 2), -3.0),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([2], [1.], 2), 0.8), # y <= 0.6
    LazySets.HalfSpace(sparsevec([2], [-1.], 2), -0.6), # y <= 0.6
    LazySets.HalfSpace(sparsevec([1, 2], [-1., 1.], 2), -2.84), # y <= 0.6
])
=#



b11 = LazySets.HalfSpace(sparsevec([1, 2], [-0.1, -1.], 2), 1.8)
b12 = LazySets.HalfSpace(sparsevec([1, 2], [0.2, -0.8], 2), 2.8)
t1 = LazySets.HalfSpace(sparsevec([1, 2], [0.08, 1.], 2), 4.6)
l1 = LazySets.HalfSpace(sparsevec([1, 2], [-1., -0.2], 2), 4.2)
r1 = LazySets.HalfSpace(sparsevec([1, 2], [1., 0.2], 2), 9.6)

inv1 = HPolyhedron([
    r1,  # x <= 0.6
    l1,  # x <= 0.6
    t1, # y <= 0.6
    b11, # y <= 0.6
    b12,
    LazySets.HalfSpace(sparsevec([1, 2], [-1., 1.], 2), 3.7), # y <= 0.6
])
constraint1 = HPolytope([
    LazySets.HalfSpace(sparsevec([1], [1.], 2), 8.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([1], [-1.], 2), -5.6),  # x <= 0.6
    LazySets.HalfSpace(sparsevec([2], [1.], 2), -0.6), # y <= 0.6
    LazySets.HalfSpace(sparsevec([2], [-1.], 2), 1.6), # y <= 0.6
    LazySets.HalfSpace(sparsevec([1, 2], [-1., 1.], 2), -6.7), # y <= 0.6
    LazySets.HalfSpace(sparsevec([1, 2], [1., -1.], 2), 7.7) # y <= 0.6
])


guard1 = HPolyhedron([
    r1,  # x <= 0.6
    l1,  # x <= 0.6
    t1, # y <= 0.6
    b11, # y <= 0.6
    b12,
    LazySets.HalfSpace(sparsevec([1, 2], [-1., -0.6], 2), -7.7), # y <= 0.6
])

b21 = LazySets.HalfSpace(sparsevec([1, 2], [-0.1, -1.], 2), 2.8)
b22 = LazySets.HalfSpace(sparsevec([1, 2], [0.2, -0.8], 2), 3.8)
t2 = LazySets.HalfSpace(sparsevec([1, 2], [0.08, 1.], 2), 6.6)
l2 = LazySets.HalfSpace(sparsevec([1, 2], [-1., -0.4], 2), -12.2)
r2 = LazySets.HalfSpace(sparsevec([1, 2], [1., 0.2], 2), 19.6)

inv2 = HPolyhedron([
    r2,  # x <= 0.6
    l2,  # x <= 0.6
    t2, # y <= 0.6
    b21, # y <= 0.6
    b22,
    LazySets.HalfSpace(sparsevec([1, 2], [-1., 1.], 2), 3.7), # y <= 0.6
])

guard2 = HPolyhedron([
    r2,  # x <= 0.6
    l2,  # x <= 0.6
    t2, # y <= 0.6
    b21, # y <= 0.6
    b22,
    LazySets.HalfSpace(sparsevec([1, 2], [-0.5, 1.], 2), -7.7), # y <= 0.6
])

palette = Plots.palette(:fes10)

fig = plot(dpi=1200, thickness_scaling=1, guidefontsize=25, minorgrid=false, showaxis=false,
    xtick=([], []),
    ytick=([], []),
    bottom_margin=2mm,
    left_margin=5mm,
    right_margin=5mm,
    top_margin=4mm
)

Plots.plot!(inv1, c=palette[4], lab="", alpha=0.2)
Plots.plot!(guard1, c=palette[5], lab="", alpha=0.6)

Plots.plot!(inv2, c=palette[4], lab="", alpha=0.2)
Plots.plot!(guard2, c=palette[5], lab="", alpha=0.6)

Z = Zonotope([0.0, 0.0], [2.0 -1.0; 1.0 0.4])
input1 = [1.1, 0.2]
input2 = [1.3, 0.8] * 0.6


Plots.plot!(Z, c=palette[2], lab="", alpha=1)

tstep = 0.05
A1 = [-1.0 0.1; 0.0 1.0]
A2 = [0.0 1.1; -1.0 0.0]
eAtstep = exp(A1 * tstep)
Ztransformed = copy(Z)
I1transformed = copy(input1)
intersectingSet = []

#lower = vertices_list(Z)
#@show lower
#lowest = [-1.0, -1.4]



for i in range(tstep, 0.55; step=tstep)
    global I1transformed = eAtstep * I1transformed + input1
    global Ztransformed = eAtstep * Ztransformed
    tset = minkowski_sum(Ztransformed, Singleton(I1transformed))
    if !LazySets.isdisjoint(guard1, tset)
        push!(intersectingSet, tset)
    end

    Plots.plot!(intersection(tset, inv1), c=palette[3], lab="", alpha=1, fa=0.0, ls=:dash)

end

cluster = intersection(overapproximate(foldl(CH, intersectingSet), 0.01), guard1)
Plots.plot!(cluster, c=palette[3], lab="", alpha=1, fa=0.5)

Plots.plot!(constraint1, c=palette[1], lab="", alpha=0.5)

#eAtstep = exp(A1 * tstep)
jumpSet = concretize(MinkowskiSum(LinearMap(exp([0.1 -0.3; -0.4 -0.2]), cluster), Singleton([4.5, 4.5])))
Plots.plot!(jumpSet, c=palette[2], lab="", alpha=1)

tstep = 0.05
eA2tstep = exp(A2 * tstep)
I2transformed = copy(input2)
intersectingSet = []
#@show typeof(jumpSet)
for i in range(0.0, 0.20; step=tstep)
    global I2transformed = eA2tstep * I2transformed + input2
    global jumpSet = linear_map(eA2tstep, jumpSet)


    Plots.plot!(intersection(inv2, LazySets.translate(jumpSet, I2transformed)), c=palette[3], lab="", alpha=1, fa=0.0, ls=:dash)

end

# INITIAL SET
plot!([-2, -0.5], [3, 0.1], arrow=true, color=:black, linewidth=1, annotations=(-2.3, 3.5, text(L"\mathcal{X}_0")), lab="")



# GUARD SET
plot!([9, 8.5], [-2, -0.5], arrow=true, color=:black, linewidth=1, annotations=(9, -2.5, text("Guard set", font(12, "Times"))), lab="")
plot!([9.5, 13.5], [-2, -1.1], arrow=true, color=:black, linewidth=1, lab="")



# DISCRETE SUCCESSOR
function curved_arrow(start_angle, end_angle; radius=1.45, steps=25)
    θs = range(start_angle, end_angle, length=steps)
    x = radius .* cos.(θs)
    y = radius .* sin.(θs)
    return x, y
end

# Draw the curved arrow
x, y = curved_arrow(-0.15π, 0.9π)

plot!((-1 * x) .+ 9.7, y .+ 4, arrow=true, linewidth=1, color=:black, ls=:dash, annotations=(6.9, 4.8, text(text(LaTeXString("post\$_d\$"), font(12, "Times")))), lab="")


# CONSTRAINT
plot!([0.2, 5.8], [-2, -1.4], arrow=true, color=:black, linewidth=1, annotations=(-0.2, -2.5, text("Constraint", font(12, "Times"))), lab="")

# INVARIANT

plot!([4, 4.5], [6, 3.5], arrow=true, color=:black, linewidth=1, annotations=(4, 6.5, text("Invariant set", font(12, "Times"))), lab="")
plot!([4.5, 10.5], [6, 5.6], arrow=true, color=:black, linewidth=1, lab="")


# ALGORITHM 1
plot!([16, 14.5], [6, 4.5], arrow=true, color=:black, linewidth=1, annotations=(16.5, 6.5, text("Algorithm 2", font(12, "Times"))), lab="")

# ALGORITHM 2
plot!([16, 18.5], [-1.5, 0.7], arrow=true, color=:black, linewidth=1, annotations=(16.5, -2.0, text("Algorithm 4", font(12, "Times"))), lab="")


#xlims!(fig, (-5.0, 20.0))
#ylims!(fig, (-3.0, 9.0))
#println(rest)
#default(fmt=:png)
#display(fig)
savefig(fig, "plots/" * "HybridSystemIllustration.pdf")
