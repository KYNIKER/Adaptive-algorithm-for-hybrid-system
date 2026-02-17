using ReachabilityAnalysis, SparseArrays
using ReachabilityAnalysis: add_dimension
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector


# Based on https://github.com/JuliaReach/ARCH2025_AFF_RE/blob/master/models/Gearbox/gearbox_benchmark.jl

# not used
function gearbox(; X0 = Hyperrectangle([0., 0., -0.0167, 0.003,  0., 0.],
                                       [0., 0.,  0.0001, 0.0001, 0., 0.]),
                   init = [(1, X0)])

    # variables
    vx = 1  # x velocity
    vy = 2  # y velocity
    px = 3  # x position
    py = 4  # y position
    I = 5   # accumulated impulse
    t = 6   # time

    # number of variables
    n = 5 + 1

    # constants
    Fs = 70.0  # shifting force on the sleeve [N]
    Tf = 1.0  # resisting moments on the second gear [N⋅m]
    # TODO where do the following values come from?
    ms = 3.2  # mass of the sleeve
    Rs = 0.08  # radius of the sleeve
    Jg₂ = 0.7  # inertia of the second gear
    Δp = -0.003  # horizontal (px) distance
    θ = 0.628318530717959  # included angle of the second gear [°]
    mg₂ = 18.1  # mass of the second gear
    ζ = 0.9  # coefficient of restitution

    # discrete structure (graph)
    automaton = GraphAutomaton(2)
    add_transition!(automaton, 1, 1, 1)
    add_transition!(automaton, 1, 1, 2)
    add_transition!(automaton, 1, 2, 3)

    # mode 1 ("free")
    A = zeros(n, n)
    b = zeros(n)
    A[px, vx] = 1.
    A[py, vy] = 1.
    b[vx] = Fs / ms
    b[vy] = - (Rs * Tf) / Jg₂
    b[t] = 1.
    invariant = HPolyhedron([
        HalfSpace(sparsevec([px], [1.], n), Δp),  # px <= Δp
        HalfSpace(sparsevec([px, py], [tan(θ), 1.], n), 0.),    # py <= -px * tan(θ)
        HalfSpace(sparsevec([px, py], [tan(θ), -1.], n), 0.)])  # py >= px * tan(θ)
    m_1 = @system(x' = Ax + b, x ∈ invariant)

    # mode 2 ("meshed")
    A = zeros(n, n)
    b = zeros(n)
    b[t] = 1.
    m_2 = @system(x' = Ax + b, x ∈ Universe(n))

    # common assignment matrix (requires individual modifications)
    A_template = zeros(n, n)
    for i in 3:6
        A_template[i, i] = 1.
    end
    denominator = ms * cos(θ)^2 + mg₂ * sin(θ)^2
    A_template[vx, vx] = (ms * cos(θ)^2 - mg₂ * ζ * sin(θ)^2) / denominator
    A_template[vx, vy] = (-(ζ + 1.) * mg₂ * sin(θ) * cos(θ)) / denominator
    A_template[vy, vx] = (-(ζ + 1.) * ms * sin(θ) * cos(θ)) / denominator
    A_template[vy, vy] = (mg₂ * sin(θ)^2 - ms * ζ * cos(θ)^2) / denominator
    A_template[I, vx] = ((ζ + 1.) * ms * mg₂ * sin(θ)) / denominator
    A_template[I, vy] = ((ζ + 1.) * ms * mg₂ * cos(θ)) / denominator

    # transition l1 -> l1
    # TODO what happened to the term '2nb' and the whole second constraint in the paper?
    guard = HPolyhedron([
        HalfSpace(sparsevec([px, py], [-tan(θ), -1.], n), 0.),     # py >= -px * tan(θ)
        HalfSpace(sparsevec([vx, vy], [-sin(θ), -cos(θ)], n), 0.)  # vx * sin(θ) + vy * cos(θ) >= 0
        ])
    A = copy(A_template)
    t1 = ConstrainedLinearMap(A, guard)

    # transition l1 -> l1
    # TODO same remark as with the other guard
    guard = HPolyhedron([
        HalfSpace(sparsevec([px, py], [-tan(θ), 1.], n), 0.),     # py <= px * tan(θ)
        HalfSpace(sparsevec([vx, vy], [-sin(θ), cos(θ)], n), 0.)  # vx * sin(θ) - vy * cos(θ) >= 0
        ])
    A = copy(A_template)
    A[vx, vy] *= -1.
    A[vy, vx] *= -1.
    A[I, vy] *= -1.
    t2 = ConstrainedLinearMap(A, guard)

    # transition l1 -> l2
    guard = HalfSpace(sparsevec([px], [-1.], n), -Δp)  # px >= Δp
    A = copy(A_template)
    A[vx, vx] = 0.
    A[vx, vy] = 0.
    A[vy, vx] = 0.
    A[vy, vy] = 0.
    A[I, vx] = A[I, vy] = ms
    t3 = ConstrainedLinearMap(A, guard)

    H = HybridSystem(automaton=automaton, modes=[m_1, m_2], resetmaps=[t1, t2, t3])

    return InitialValueProblem(H, init)
end

# homogenized version

function gearbox_homog(; X0 = Hyperrectangle([0., 0., -0.0167, 0.003,  0., 1.],
                                             [0., 0.,  0.0001, 0.0001, 0., 0.]),
                         init = [(1, X0)])

    # variables
    vx = 1  # x velocity
    vy = 2  # y velocity
    px = 3  # x position
    py = 4  # y position
    I = 5   # accumulated impulse
    η = 6   # accounts for the homogenized input term  TODO unused

    # number of variables
    n = 5 + 1

    # constants
    Fs = 70.0  # shifting force on the sleeve [N]
    Tf = 1.0  # resisting moments on the second gear [N⋅m]
    # TODO where do the following values come from?
    ms = 3.2  # mass of the sleeve
    Rs = 0.08  # radius of the sleeve
    Jg₂ = 0.7  # inertia of the second gear
    Δp = -0.003  # horizontal (px) distance
    θ = 0.628318530717959  # included angle of the second gear [°]
    mg₂ = 18.1  # mass of the second gear
    ζ = 0.9  # coefficient of restitution

    # discrete structure (graph)
    automaton = GraphAutomaton(2)
    add_transition!(automaton, 1, 1, 1)
    add_transition!(automaton, 1, 1, 2)
    add_transition!(automaton, 1, 2, 3)

    # mode 1 ("free")
    A = zeros(n-1, n-1)
    b = zeros(n-1)
    A[px, vx] = 1.
    A[py, vy] = 1.
    b[vx] = Fs / ms
    b[vy] = - (Rs * Tf) / Jg₂
    invariant = HPolyhedron([
        HalfSpace(sparsevec([px], [1.], n), Δp),  # px <= Δp
        HalfSpace(sparsevec([px, py], [tan(θ), 1.], n), 0.),    # py <= -px * tan(θ)
        HalfSpace(sparsevec([px, py], [tan(θ), -1.], n), 0.)])  # py >= px * tan(θ)
    Aext = add_dimension(A)
    Aext[1:n-1, n] = b
    m_1 = @system(x' = Aext * x, x ∈ invariant)

    # mode 2 ("meshed")
    A0 = zeros(n, n)
    m_2 = @system(x' = A0 * x, x ∈ Universe(n))

    # common assignment matrix (requires individual modifications)
    A_template = zeros(n, n)
    for i in 3:6
        A_template[i, i] = 1.
    end
    denominator = ms * cos(θ)^2 + mg₂ * sin(θ)^2
    A_template[vx, vx] = (ms * cos(θ)^2 - mg₂ * ζ * sin(θ)^2) / denominator
    A_template[vx, vy] = (-(ζ + 1.) * mg₂ * sin(θ) * cos(θ)) / denominator
    A_template[vy, vx] = (-(ζ + 1.) * ms * sin(θ) * cos(θ)) / denominator
    A_template[vy, vy] = (mg₂ * sin(θ)^2 - ms * ζ * cos(θ)^2) / denominator
    A_template[I, vx] = ((ζ + 1.) * ms * mg₂ * sin(θ)) / denominator
    A_template[I, vy] = ((ζ + 1.) * ms * mg₂ * cos(θ)) / denominator

    # transition l1 -> l1
    # TODO what happened to the term '2nb' and the whole second constraint in the paper?
    guard = HPolyhedron([
        HalfSpace(sparsevec([px, py], [-tan(θ), -1.], n), 0.),     # py >= -px * tan(θ)
        HalfSpace(sparsevec([vx, vy], [-sin(θ), -cos(θ)], n), 0.)  # vx * sin(θ) + vy * cos(θ) >= 0
        ])
    A = copy(A_template)
    t1 = ConstrainedLinearMap(A, guard)

    # transition l1 -> l1
    # TODO same remark as with the other guard
    guard = HPolyhedron([
        HalfSpace(sparsevec([px, py], [-tan(θ), 1.], n), 0.),     # py <= px * tan(θ)
        HalfSpace(sparsevec([vx, vy], [-sin(θ), cos(θ)], n), 0.)  # vx * sin(θ) - vy * cos(θ) >= 0
        ])
    A = copy(A_template)
    A[vx, vy] *= -1.
    A[vy, vx] *= -1.
    A[I, vy] *= -1.
    t2 = ConstrainedLinearMap(A, guard)

    # transition l1 -> l2
    guard = HalfSpace(SingleEntryVector(px, n, -1.), -Δp)  # px >= Δp
    A = copy(A_template)
    A[vx, vx] = 0.
    A[vx, vy] = 0.
    A[vy, vx] = 0.
    A[vy, vy] = 0.
    A[I, vx] = A[I, vy] = ms
    t3 = ConstrainedLinearMap(A, guard)

    H = HybridSystem(automaton=automaton, modes=[m_1, m_2], resetmaps=[t1, t2, t3])

    return InitialValueProblem(H, init)
end

# ==============================================================================
# Run the example
# ==============================================================================




using Plots, Plots.PlotMeasures, LaTeXStrings
using ReachabilityAnalysis, LazySets, SparseArrays
using ReachabilityAnalysis.ReachabilityBase.Timing: print_timed
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector

const model = "Gearbox"
# const cases = [
#                "GRBX01-MES01",
#                "GRBX01-MES01-discrete",
#                "GRBX02-MES01",
#                "GRBX02-MES01-discrete"
#               ]

LazySets.deactivate_assertions()


# ==============================================================================
# Analysis function
# ==============================================================================

const n = 6

# template directions

function template_directions(n)
    vx, vy, px, py, θ = 1, 2, 3, 4, 0.628318530717959
    vec1 = sparsevec([px, py], [tan(θ), 1.], n)
    vec2 = sparsevec([px, py], [tan(θ), -1.], n)
    vec1perp = sparsevec([px, py], [1, -tan(θ)], n)
    vec2perp = sparsevec([px, py], [1., tan(θ)], n)
    dirs_plane_34 = [vec1, -vec1, vec2, -vec2, vec1perp, -vec1perp, vec2perp, -vec2perp]

    vec3 = sparsevec([vx, vy], [-sin(θ), -cos(θ)], n)
    vec4 = sparsevec([vx, vy], [-sin(θ), cos(θ)], n)
    vec3perp = sparsevec([vx, vy], [cos(θ), -sin(θ)], n)
    vec4perp = sparsevec([vx, vy], [cos(θ), sin(θ)], n)

    dirs_plane_12 = [vec3, -vec3, vec4, -vec4, vec3perp, -vec3perp, vec4perp, -vec4perp]
    boxdirs = BoxDirections{Float64,Vector{Float64}}(n)
    return vcat(collect(boxdirs), Vector.(dirs_plane_34), Vector.(dirs_plane_12)) |> CustomDirections
end

const extdirs = template_directions(n)
const thull = TemplateHullIntersection(extdirs)

function analyze_once(prob, alg)
    # solve
    sol = solve(prob, max_jumps=100, clustering_method=LazyClustering(1, convex=false),
                intersect_source_invariant=true, intersection_source_invariant_method=thull,
                intersection_method=thull, tspan = 0 .. 0.21, fixpoint_check=true, alg=alg)

    # verify that specification holds
    validation = true
    # property 1
    for F in array(sol.F)
        if location(F) == 1
            if tend(F) >= 0.2
                validation = false
                break
            end
        end
    end
    # property 2
    e5 = SingleEntryVector(5, n, 1.0)
    validation &= ρ(e5, sol) < 20.0

    return (sol, validation)
end

function analyze(; prob, alg, case, warmup)
    println("case: $case")

    # warm-up run
    if warmup
        analyze_once(prob, alg)
    end

    # benchmark run
    res = @timed analyze_once(prob, alg)
    print_timed(res)
    sol, validation = res.value
    validation || throw(ErrorException("the property is violated"))

    # report results
    io = isdefined(Main, :io) ? Main.io : stdout
    result = Int(validation)
    time = trunc(res.time, digits=3)
    print(io, "$model,$case,$result,$time\n")

    return sol
end

# ==============================================================================
# Run benchmarks
# ==============================================================================



# ----------------------------------------
#  GRBX01
# ----------------------------------------

# initial states
# non-homogeneous system (not used)
# X0_GRBX01 = Hyperrectangle(low=[0, 0, -0.0168, 0.0029, 0, 0],
#                            high=[0, 0, -0.0166, 0.0031, 0, 0])
# homogenized system
X0_GRBX01h = Hyperrectangle(low=[0, 0, -0.0168, 0.0029, 0, 1],
                            high=[0, 0, -0.0166, 0.0031, 0, 1])

prob_GRBX01 = gearbox_homog(X0=X0_GRBX01h)

case = "GRBX01-MES01"
alg = LGG09(δ=0.0008, template=extdirs, approx_model=Forward())
alg = LGG09(δ=0.0004, template=extdirs, approx_model=Forward())



sol_GRBX01 = analyze(prob=prob_GRBX01, alg=alg, case=case, warmup=true)

# case = "GRBX01-MES01-discrete"
# alg = LGG09(δ=δ_discrete, template=extdirs, cache=true, approx_model=NoBloating())
# analyze(prob=prob_GRBX01, alg=alg, case=case, warmup=true)

# ----------------------------------------
#  GRBX02
# ----------------------------------------

# initial states
# non-homogeneous system (not used)
# X0_GRBX02 = Hyperrectangle(low=[0, 0, -0.01675, 0.00285, 0, 0],
#                            high=[0, 0, -0.01665, 0.00315, 0, 0])
# homogenized system
X0_GRBX02h = Hyperrectangle(low=[0, 0, -0.01675, 0.00285, 0, 1],
                            high=[0, 0, -0.01665, 0.00315, 0, 1])

prob_GRBX02 = gearbox_homog(X0=X0_GRBX02h)

case = "GRBX02-MES01"
alg = LGG09(δ=0.0008, template=extdirs, approx_model=Forward())
analyze(prob=prob_GRBX02, alg=alg, case=case, warmup=false)

# case = "GRBX02-MES01-discrete"
# alg = LGG09(δ=δ_discrete, template=extdirs, cache=true, approx_model=NoBloating())
# analyze(prob=prob_GRBX02, alg=alg, case=case, warmup=false)

# ==============================================================================
# Plot
# # ==============================================================================

# const TARGET_FOLDER = isdefined(Main, :TARGET_FOLDER) ? Main.TARGET_FOLDER : @__DIR__

# # modify tolerance
# LazySets.set_ztol(Float64, 1e-8)

const fig = plot()
plot!(fig, sol_GRBX01, vars=(1, 2), ε=1e-5)
# plot!(fig, sol_GRBX01, vars=(3, 4), ε=1e-5,
#       color=:blue, alpha=0.5, lw=1.0, linecolor=:blue,
#       tickfont=font(30, "Times"), guidefontsize=45,
#       xlab=L"x_3", ylab=L"x_4",
#       xtick=([-0.016, -0.009, -0.002],
#              [L"-0.016", L"-0.009", L"-0.002"]),
#       ytick=([-0.008, -0.004, 0.0, 0.004],
#              [L"-0.008", L"-0.004", L"0", L"0.004"]),
#       xlims=(-0.017, -0.0015), ylims=(-0.008, 0.004),
#       bottom_margin=-5mm, left_margin=-1mm, right_margin=10mm, top_margin=3mm,
#       size=(1000, 1000))
savefig(fig, joinpath("results/", "ARCH-COMP25-JuliaReach-$model-GRBX01.png"))

# # reset tolerance
# LazySets.set_tolerance(Float64)

# sol_GRBX01 = nothing
# GC.gc()


nothing