using SparseArrays: sparsevec 
using SparseArrays: spzeros
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityAnalysis: add_dimension
using LazySets

include("../Utilities.jl")

# This code is based on https://github.com/JuliaReach/ARCH2025_AFF_RE/tree/master/models/Powertrain
# It is adapted to our framework

const t_init = 0.2  # time to stay in the initial location
const α = 0.03  # backlash size (half of the gap width)
const τ_eng = 0.1  # engine time constant
const γ = 12.0  # gearbox ratio (dimensionless)
const u = 5.0  # requested engine torque

# moments of inertia [kg m²]
const Jₗ = 140.
const Jₘ = 0.3
const Jᵢ = 100.0  # NOTE: the paper says 0.01, but the CORA/SpaceEx models use 100

# viscous friction constants [Nm s/rad]
const bₗ = 5.6
const bₘ = 0.0
const bᵢ = 1.0

# shaft stiffness [Nm/rad]
const kᵢ = 1e5
const kₛ = 1e4

# PID parameters
const k_P = 0.5  # [Nms/rad]
const k_I = 0.5  # [Nm/rad]
const k_D = 0.5  # [Nms²/rad]


function print_dynamics(A, b, location_name)
    println("dynamics of location $location_name:")
    for i in 1:size(A, 1)-1  # ignore the last dimension (time)
        print("x_$i' = ")
        for j in 1:size(A, 2)
            if !iszero(A[i,j])
                print("$(A[i,j]) x_$j + ")
            end
        end
        println("$(b[i])\n")
    end
end


function get_dynamics(kₛ, α, u, n, θ)
    # linear dynamics
    A = spzeros(n, n)

    A[1, 7] = 1.0 / γ
    A[1, 9] = -1.0

    A[2, 1] = (-k_I * γ + k_D * kₛ / (γ * Jₘ)) / τ_eng
    A[2, 2] = (-k_D / Jₘ - 1.0) / τ_eng
    A[2, 3] = k_I * γ / τ_eng
    A[2, 4] = k_P * γ / τ_eng
    A[2, 7] = (-k_P + k_D * bₘ /Jₘ) / τ_eng
    A[2, 8] = -k_I * γ / τ_eng

    A[3, 4] = 1.0

    A[5, 6] = 1.0

    A[6, 5] = -kᵢ / Jₗ
    A[6, 6] = -bₗ / Jₗ
    A[6, 2*θ+6] = kᵢ / Jₗ

    A[7, 1] = -kₛ / (Jₘ * γ)
    A[7, 2] = 1.0 / Jₘ
    A[7, 7] = -bₘ / Jₘ

    i = 8
    while i < n-1
        A[i, i+1] = 1.0

        if i == 8
            # x9 has special dynamics
            A[i+1, 1] = kₛ / Jᵢ
            A[i+1, i] = -kᵢ / Jᵢ
        else
            A[i+1, i-2] = kᵢ / Jᵢ
            A[i+1, i] = -2. * kᵢ / Jᵢ
        end
        A[i+1, i+1] = -bᵢ / Jᵢ

        # wrap-around to x5 in the last step
        j = (i == n-2) ? 5 : i+2
        A[i+1, j] = kᵢ / Jᵢ

        i += 2
    end

    # affine vector
    b = spzeros(n)
    b[2] = k_D * (γ * u - kₛ * α / (Jₘ * γ)) / τ_eng
    b[4] = u
    b[7] = kₛ * α / (Jₘ * γ)
    b[9] = -kₛ * α / Jᵢ
    b[n] = 1.0  # time

    return A, Vector(b)
end


function get_initial_condition(n, X0_scale)
    c = Vector{Float64}(undef, n)
    g = Vector{Float64}(undef, n)
    c[1:7] = [-0.0432, -11., 0., 30., 0., 30., 360.]
    g[1:7] = [0.0056, 4.67, 0., 10., 0., 10., 120.]
    i = 8
    while i < n
        c[i] = -0.0013
        g[i] = 0.0006
        i += 1
        c[i] = 30.
        g[i] = 10.
        i += 1
    end
    c[n] = 0.0
    g[n] = 0.0
    if X0_scale < 1.0
        g = X0_scale * g
    end

    return Zonotope(c, hcat(g))
end


function _inhomogPowertrain(; θ::Int=1, X0_scale::Float64=1.0)
    θ > 0 || error("θ must be positive, but was $θ")
    (X0_scale > 0.0 && X0_scale <= 1.0) || error("scale $X0_scale ∉ (0, 1]")

    # dimension of state space (last dimension is time)
    n = 2 * θ + 7 + 1

    # hybrid automaton

    negAngleEdges = Vector{Edge}()
    deadzoneEdges = Vector{Edge}()
    posAngleEdges = Vector{Edge}()
    initAngleEdges = Vector{Edge}()

    # push!(initAngleEdges, Edge(1, HPolyhedron([LazySets.HalfSpace(sparsevec([n], [-1.], n), -t_init)]), I(n+1), zeros(n+1))) # t >= t_init
    # push!(negAngleEdges, Edge(2, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [-1.], n), α)]), I(n+1), zeros(n+1))) # x1 >= -α
    # push!(deadzoneEdges, Edge(3, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [-1.], n), -α)]), I(n+1), zeros(n+1)))  # x1 >= α
    # push!(deadzoneEdges, Edge(1, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [1.], n), -α)]), I(n+1), zeros(n+1))) # x1 <= -α
    # push!(posAngleEdges, Edge(2, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [1.], n), α)]), I(n+1), zeros(n+1))) # x1 <= α
    
    push!(initAngleEdges, Edge(1, HPolyhedron([LazySets.HalfSpace(sparsevec([n], [-1.], n), -t_init)]), I(n), zeros(n))) # t >= t_init
    push!(negAngleEdges, Edge(2, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [-1.], n), α)]), I(n), zeros(n))) # x1 >= -α
    push!(deadzoneEdges, Edge(3, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [-1.], n), -α)]), I(n), zeros(n)))  # x1 >= α
    #push!(deadzoneEdges, Edge(1, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [1.], n), -α)]), I(n), zeros(n))) # x1 <= -α
    #push!(posAngleEdges, Edge(2, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [1.], n), α)]), I(n), zeros(n))) # x1 <= α
    


    locations = Vector{Location}()

    # negAngle
    A, b = get_dynamics(kₛ, -α, u, n, θ)
    inv = HPolyhedron([LazySets.HalfSpace(sparsevec([1], [1.], n), -α)])  # x1 <= -α
    # Aext = add_dimension(A, 1)
    # Aext[1:n, n+1] .= b
    uInput = Zonotope(b, [zero(b)])
    constraints = []
    push!(locations, Location(1, inv, A, nothing, uInput, nothing, negAngleEdges, []))


    # deadzone
    A, b = get_dynamics(0., -α, u, n, θ)
    inv = HPolyhedron([LazySets.HalfSpace(sparsevec([1], [-1.], n), α),  # x1 >= -α
                     LazySets.HalfSpace(sparsevec([1], [1.], n), α)])  # x1 <= α
    # Aext = add_dimension(A, 1)
    # Aext[1:n, n+1] .= b
    #m_deadzone = @system(x' = Aext * x, x ∈ X)
    uInput = Zonotope(b, [zero(b)])
    constraints = []
    push!(locations, Location(2, inv, A, nothing, uInput, nothing, deadzoneEdges, []))


    # posAngle
    A, b = get_dynamics(kₛ, α, u, n, θ)
    inv = HPolyhedron([LazySets.HalfSpace(sparsevec([1], [-1.], n), -α)])  # x1 >= α
    # Aext = add_dimension(A, 1)
    # Aext[1:n, n+1] .= b
    #m_posAngle = @system(x' = Aext * x, x ∈ X)
    uInput = Zonotope(b, [zero(b)])
    constraints = []
    push!(locations, Location(3, inv, A, nothing, uInput, nothing, posAngleEdges, []))


    # negAngleInit
    A, b = get_dynamics(kₛ, -α, -u, n, θ)
    #X = HalfSpace(SingleEntryVector(n, n+1, 1.), t_init)  # t <= t_init
    inv = HPolyhedron([LazySets.HalfSpace(sparsevec([n], [1.], n), t_init)]) # t <= t_init
    # Aext = add_dimension(A, 1)
    # Aext[1:n, n+1] .= b
    uInput = Zonotope(b, [zero(b)])
    #m_negAngleInit = @system(x' = Aext * x, x ∈ X)
    push!(locations, Location(4, inv, A, nothing, uInput, nothing, initAngleEdges, []))


    # switching
    #switchings = [HybridSystems.AutonomousSwitching()]

    globalConstraints = []
    H = HybridSystemV2(locations, globalConstraints)

    # initial condition
    #X0_orig = get_initial_condition(n, X0_scale)
    #X0 = cartesian_product(X0_orig, Singleton([1.0]))
    X0 = get_initial_condition(n, X0_scale)

    T = 2

    # Some issues here at the end...

    return H, 4, X0, T
end


function _homogPowertrain(; θ::Int=1, X0_scale::Float64=1.0)

    θ > 0 || error("θ must be positive, but was $θ")
    (X0_scale > 0.0 && X0_scale <= 1.0) || error("scale $X0_scale ∉ (0, 1]")

    # dimension of state space (last dimension is time)
    n = 2 * θ + 7 + 1

    # hybrid automaton

    negAngleEdges = Vector{Edge}()
    deadzoneEdges = Vector{Edge}()
    posAngleEdges = Vector{Edge}()
    initAngleEdges = Vector{Edge}()

    push!(initAngleEdges, Edge(1, HPolyhedron([LazySets.HalfSpace(sparsevec([n], [-1.], n+1), -t_init)]), I(n+1), zeros(n+1))) # t >= t_init
    push!(negAngleEdges, Edge(2, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [-1.], n+1), α)]), I(n+1), zeros(n+1))) # x1 >= -α
    push!(deadzoneEdges, Edge(3, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [-1.], n+1), -α)]), I(n+1), zeros(n+1)))  # x1 >= α
    #push!(deadzoneEdges, Edge(1, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [1.], n+1), -α)]), I(n+1), zeros(n+1))) # x1 <= -α
    #push!(posAngleEdges, Edge(2, HPolyhedron([LazySets.HalfSpace(sparsevec([1], [1.], n+1), α)]), I(n+1), zeros(n+1))) # x1 <= α
    deadZoneConstraint = [LazySets.HalfSpace(sparsevec([1], [-1.], n+1), -α)]
    posAngleConstraint = [LazySets.HalfSpace(sparsevec([1], [-1.], n+1), -α)]

    deadZoneConstraint = []
    posAngleConstraint = []

    locations = Vector{Location}()

    # negAngle
    A, b = get_dynamics(kₛ, -α, u, n, θ)
    inv = HPolyhedron([LazySets.HalfSpace(sparsevec([1], [1.], n+1), -α)])  # x1 <= -α
    Aext = add_dimension(A, 1)
    Aext[1:n, n+1] .= b
    uInput = nothing
    constraints = []
    push!(locations, Location(1, inv, Aext, nothing, uInput, nothing, negAngleEdges, []))


    # deadzone
    A, b = get_dynamics(0., -α, u, n, θ)
    inv = HPolyhedron([LazySets.HalfSpace(sparsevec([1], [-1.], n+1), α),  # x1 >= -α
                     LazySets.HalfSpace(sparsevec([1], [1.], n+1), α)])  # x1 <= α
    Aext = add_dimension(A, 1)
    Aext[1:n, n+1] .= b
    #m_deadzone = @system(x' = Aext * x, x ∈ X)
    uInput = nothing
    constraints = []
    push!(locations, Location(2, inv, Aext, nothing, uInput, nothing, deadzoneEdges, deadZoneConstraint))


    # posAngle
    A, b = get_dynamics(kₛ, α, u, n, θ)
    inv = HPolyhedron([LazySets.HalfSpace(sparsevec([1], [-1.], n+1), -α)])  # x1 >= α
    Aext = add_dimension(A, 1)
    Aext[1:n, n+1] .= b
    #m_posAngle = @system(x' = Aext * x, x ∈ X)
    uInput = nothing
    constraints = []
    push!(locations, Location(3, inv, Aext, nothing, uInput, nothing, posAngleEdges, posAngleConstraint))


    # negAngleInit
    A, b = get_dynamics(kₛ, -α, -u, n, θ)
    #X = HalfSpace(SingleEntryVector(n, n+1, 1.), t_init)  # t <= t_init
    inv = HPolyhedron([LazySets.HalfSpace(sparsevec([n], [1.], n+1), t_init)]) # t <= t_init
    Aext = add_dimension(A, 1)
    Aext[1:n, n+1] .= b
    uInput = nothing
    #m_negAngleInit = @system(x' = Aext * x, x ∈ X)
    push!(locations, Location(4, inv, Aext, nothing, uInput, nothing, initAngleEdges, []))


    # switching
    #switchings = [HybridSystems.AutonomousSwitching()]

    globalConstraints = []
    H = HybridSystemV2(locations, globalConstraints)

    # initial condition
    X0_orig = get_initial_condition(n, X0_scale)
    X0 = cartesian_product(X0_orig, Singleton([1.0]))
    #X0 = get_initial_condition(n, X0_scale)

    T = 2

    # Some issues here at the end...

    return H, 4, X0, T
end


function loadPowertrain(; θ::Int=1, X0_scale::Float64=1.0, homog::Bool = true)
    if homog
        return _homogPowertrain(θ=θ, X0_scale=X0_scale)
    end
    return _inhomogPowertrain(θ=θ, X0_scale=X0_scale)
end