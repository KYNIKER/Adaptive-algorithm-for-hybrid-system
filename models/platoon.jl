using SparseArrays: sparsevec, sparse
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityAnalysis: add_dimension
using LazySets

include("../Utilities.jl")

# Location 1
function platoon_connected(; deterministic_switching::Bool=true,
    c1=5.0)  # clock constraints
    n = 10  # 9 dimensions + time
    # x' = Ax + Bu + c
    A = Matrix{Float64}(undef, n, n)
    A[1, :] = [0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0]
    A[2, :] = [0, 0, -1.0, 0, 0, 0, 0, 0, 0, 0]
    A[3, :] = [1.6050, 4.8680, -3.5754, -0.8198, 0.4270, -0.0450, -0.1942, 0.3626, -0.0946, 0.]
    A[4, :] = [0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0,]
    A[5, :] = [0, 0, 1.0, 0, 0, -1.0, 0, 0, 0, 0]
    A[6, :] = [0.8718, 3.8140, -0.0754, 1.1936, 3.6258, -3.2396, -0.5950, 0.1294, -0.0796, 0.]
    A[7, :] = [0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0]
    A[8, :] = [0, 0, 0, 0, 0, 1.0, 0, 0, -1.0, 0]
    A[9, :] = [0.7132, 3.5730, -0.0964, 0.8472, 3.2568, -0.0876, 1.2726, 3.0720, -3.1356, 0.]
    A[10, :] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0]  # t' = 1

    if deterministic_switching
        invariant = HPolyhedron([LazySets.HalfSpace(sparsevec([n], [1.], n), c1)])  # t <= c1
    else
        invariant = nothing
    end

    # acceleration of the lead vehicle + time
    B = sparse([2], [1], [1.0], n, 1)
    u = Hyperrectangle(low=[-9.], high=[1.])



    c = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0]
    #@system(x' = A * x + B * u + c, x ∈ invariant, u ∈ U)
    return A, B, u, c, invariant
end

# Location 2
function platoon_disconnected(; deterministic_switching::Bool=true,
    c2=5.0)  # clock constraints
    n = 10  # 9 dimensions + time
    # x' = Ax + Bu + c
    A = Matrix{Float64}(undef, n, n)
    A[1, :] = [0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0]
    A[2, :] = [0, 0, -1.0, 0, 0, 0, 0, 0, 0, 0]
    A[3, :] = [1.6050, 4.8680, -3.5754, 0, 0, 0, 0, 0, 0, 0]
    A[4, :] = [0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0,]
    A[5, :] = [0, 0, 1.0, 0, 0, -1.0, 0, 0, 0, 0]
    A[6, :] = [0, 0, 0, 1.1936, 3.6258, -3.2396, 0, 0, 0, 0.]
    A[7, :] = [0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0]
    A[8, :] = [0, 0, 0, 0, 0, 1.0, 0, 0, -1.0, 0]
    A[9, :] = [0.7132, 3.5730, -0.0964, 0.8472, 3.2568, -0.0876, 1.2726, 3.0720, -3.1356, 0.]
    A[10, :] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0]  # t' = 1

    if deterministic_switching
        invariant = HPolyhedron([LazySets.HalfSpace(sparsevec([n], [1.], n), c2)])  # t <= c2
    else
        invariant = nothing
    end

    # acceleration of the lead vehicle + time
    B = sparse([2], [1], [1.0], n, 1)
    u = Hyperrectangle(low=[-9.], high=[1.])

    c = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0]
    #@system(x' = A * x + B * u + c, x ∈ invariant, u ∈ U)
    return A, B, u, c, invariant
end

function loadPlatoon(; deterministic_switching::Bool=true,
    c1=5.0,  # clock constraints
    c2=5.0,  # clock constraints
    tb=10.0,  # lower bound for loss of communication
    tc=20.0, tr=20.0)  # upper bound for loss of communication (tc) and reset time (tr)

    # three variables for each vehicle, (ei, d(et)/dt, ai) for
    # (spacing error, relative velocity, speed), and the last dimension is time
    n = 9 + 1


    edgeListLoc1 = Vector{Edge}()
    edgeListLoc2 = Vector{Edge}()

    locations = Vector{Location}()

    # common reset
    #reset = Dict(n => 0.)
    reset = Diagonal(ones(n))

    reset[n, n] = 0.  # We reset the time when taking a guard


    # transition l1 -> l2
    if deterministic_switching
        guard = HPolyhedron([#LazySets.HalfSpace(SingleEntryVector(n, n, 1.), c1), # This part is enforced by invarients
            LazySets.HalfSpace(sparsevec([n], [-1.], n), -c1)])  # t >= c1
    else
        # tb <= t <= tc
        guard = HPolyhedron([LazySets.HalfSpace(sparsevec([n], [-1.], n), -tb),
            LazySets.HalfSpace(sparsevec([n], [1.], n), tc)])
    end
    #t1 = ConstrainedResetMap(n, guard, reset)
    push!(edgeListLoc1, Edge(2, guard, reset, zeros(n)))


    # transition l2 -> l1
    if deterministic_switching
        guard = HPolyhedron([#LazySets.HalfSpace(SingleEntryVector(n, n, 1.), c2),  # This part is handled by invarients
            LazySets.HalfSpace(sparsevec([n], [-1.], n), -c2)])  # t >= c2
    else
        guard = HPolyhedron([LazySets.HalfSpace(sparsevec([n], [1.], n), tr)])  # t <= tr
    end
    #t2 = ConstrainedResetMap(n, guard, reset)

    push!(edgeListLoc2, Edge(1, guard, reset, zeros(n)))

    #resetmaps = [t1, t2]

    A1, B1, u1, constant1, invariant1 = platoon_connected(deterministic_switching=deterministic_switching, c1=c1)
    push!(locations, Location(1, invariant1, A1, B1, u1, constant1, edgeListLoc1, []))
    A2, B2, u2, constant2, invariant2 = platoon_disconnected(deterministic_switching=deterministic_switching, c2=c2)
    push!(locations, Location(2, invariant2, A2, B2, u2, constant2, edgeListLoc2, []))


    #X0 = Singleton(zeros(n))
    #X0 = convert(Zonotope, X0)
    X0 = Zonotope(zeros(n), [zeros(n)])
    # X0 = Zonotope(zeros(n))

    # Global constraints
    properties = []

    dmin = -30.  # Alternatively this is -30, -42 or -50
    push!(properties, LazySets.HalfSpace(sparsevec([1], [-1.], n), -dmin))  # >= -dmin
    push!(properties, LazySets.HalfSpace(sparsevec([4], [-1.], n), -dmin))  # >= -dmin
    push!(properties, LazySets.HalfSpace(sparsevec([7], [-1.], n), -dmin))  # >= -dmin

    H = HybridSystemV2(locations, properties)

    T = 20

    return H, 1, X0, T
    # H = HybridSystem(automaton, modes, resetmaps, [AutonomousSwitching()])

    # # initial condition is at the orgin in mode 1
    # X0 = Singleton(zeros(n))
    # initial_condition = [(1, X0)]

    # return IVP(H, initial_condition)
end

#loadPlatoon()