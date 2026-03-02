using SparseArrays: sparsevec
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityAnalysis: add_dimension
using LazySets

include("../Utilities.jl")

function loadExampleSystem()
    dampening = 0.5
    w = 1.0
    A = [0.0 1.0; -w^2 0.0]
    c = zeros(2)
    u = Zonotope(c, [zero(c)])

    guard = HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [-1.], 2), -0.6),  # x >= 0.6
    ])

    jumpMatrix = [dampening 0.0; 0.0 -1.]

    edges::Vector{Edge} = [Edge(1, guard, jumpMatrix, zeros(2))]
    #push!(edges, )

    invarient = HPolyhedron([
        LazySets.HalfSpace(sparsevec([1], [1.], 2), 0.6),  # x <= 0.6
    ])
    #invarient = guard
    #invarient = Universe(2)
    #invarient = nothing

    locations = [Location(1, invarient, A, I(2), u, nothing, edges, [])]

    H = HybridSystemV2(locations, [])

    X0 = Zonotope([-1, 0.], [[0.01, 0.0]])
    T = 5

    return H, 1, X0, T
end

function getRealExampleSystem()
    stepsPerSec = 100
    _, _, X0, T = loadExampleSystem()

    totalSteps = T * stepsPerSec

    t = range(0, T, length = totalSteps)
    x = zeros(totalSteps)
    v = zeros(totalSteps)

    x[1] = X0.center[1]
    v[1] = X0.center[2]

    

    for k in 1:totalSteps-1
        xNext = x[k] + 1/stepsPerSec * v[k]
        vNext = v[k] + 1/stepsPerSec * (x[k] * -1^2)

        if xNext >= 0.6
            xNext = 0.5 * xNext
            vNext = -1 * vNext
        end

        x[k+1] = xNext
        v[k+1] = vNext
    end
    return x, v, t
end

# loadBouncingBall()