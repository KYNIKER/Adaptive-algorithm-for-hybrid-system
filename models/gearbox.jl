using SparseArrays: sparsevec
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityAnalysis: add_dimension
using LazySets

include("../Utilities.jl")


function loadGearBox()
    X0 = Hyperrectangle(low=[0, 0, -0.0168, 0.0029, 0, 1],
        high=[0, 0, -0.0166, 0.0031, 0, 1])
    X0 = convert(Zonotope, X0)
    X0 = Zonotope(X0.center, X0.generators)

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


    # automaton = GraphAutomaton(2)
    # add_transition!(automaton, 1, 1, 1)
    # add_transition!(automaton, 1, 1, 2)
    # add_transition!(automaton, 1, 2, 3)



    # edges = [
    #     (1, 1),
    #     (1, 1),
    #     (1, 2)
    # ]
    edgeListLoc1::Vector{Edge} = []

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
        LazySets.HalfSpace(sparsevec([px, py], [-tan(θ), -1.], n), 0.),     # py >= -px * tan(θ)
        LazySets.HalfSpace(sparsevec([vx, vy], [-sin(θ), -cos(θ)], n), 0.)  # vx * sin(θ) + vy * cos(θ) >= 0
    ])
    A = copy(A_template)
    #t1 = ConstrainedLinearMap(A, guard)

    push!(edgeListLoc1, Edge(1, guard, A, zeros(n)))

    # transition l1 -> l1
    # TODO same remark as with the other guard
    guard = HPolyhedron([
        LazySets.HalfSpace(sparsevec([px, py], [-tan(θ), 1.], n), 0.),     # py <= px * tan(θ)
        LazySets.HalfSpace(sparsevec([vx, vy], [-sin(θ), cos(θ)], n), 0.)  # vx * sin(θ) - vy * cos(θ) >= 0
    ])
    A = copy(A_template)
    A[vx, vy] *= -1.
    A[vy, vx] *= -1.
    A[I, vy] *= -1.
    #t2 = ConstrainedLinearMap(A, guard)

    push!(edgeListLoc1, Edge(1, guard, A, zeros(n)))

    # transition l1 -> l2
    guard = LazySets.HalfSpace(SingleEntryVector(px, n, -1.), -Δp)  # px >= Δp
    A = copy(A_template)
    A[vx, vx] = 0.
    A[vx, vy] = 0.
    A[vy, vx] = 0.
    A[vy, vy] = 0.
    A[I, vx] = A[I, vy] = ms
    #t3 = ConstrainedLinearMap(A, guard)

    push!(edgeListLoc1, Edge(2, guard, A, zeros(n)))


    locations = []

    # mode 1 ("free")
    A = zeros(n - 1, n - 1)
    b = zeros(n - 1)
    A[px, vx] = 1.
    A[py, vy] = 1.
    b[vx] = Fs / ms
    b[vy] = -(Rs * Tf) / Jg₂
    invariant = HPolyhedron([
        LazySets.HalfSpace(sparsevec([px], [1.], n), Δp),  # px <= Δp
        LazySets.HalfSpace(sparsevec([px, py], [tan(θ), 1.], n), 0.),    # py <= -px * tan(θ)
        LazySets.HalfSpace(sparsevec([px, py], [tan(θ), -1.], n), 0.)])  # py >= px * tan(θ)
    Aext = add_dimension(A)
    Aext[1:n-1, n] = b


    push!(locations, Location(1,                        # ID
        HPolyhedron([   # Invaariant
            LazySets.HalfSpace(sparsevec([px], [1.], n), Δp),  # px <= Δp
            LazySets.HalfSpace(sparsevec([px, py], [tan(θ), 1.], n), 0.),    # py <= -px * tan(θ)
            LazySets.HalfSpace(sparsevec([px, py], [tan(θ), -1.], n), 0.)]),  # py >= px * tan(θ)
        Aext,           # Flow matrix
        nothing, # B input
        nothing, # u
        nothing, # constant input
        edgeListLoc1)
    )

    #m_1 = @system(x' = Aext * x, x ∈ invariant)

    # mode 2 ("meshed")
    A0 = zeros(n, n)
    push!(locations, Location(2, nothing, A0, nothing, nothing, nothing, []))
    #m_2 = @system(x' = A0 * x, x ∈ Universe(n))


    H = HybridSystemV2(locations, 1, X0)

    return H
end


# loadGearBox()