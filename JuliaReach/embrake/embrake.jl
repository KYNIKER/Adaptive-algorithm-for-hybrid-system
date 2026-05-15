using ReachabilityAnalysis, SparseArrays

function embrake_parameters()
    L = 1.e-3
    KP = 10000.
    KI = 1000.
    R = 0.5
    K = 0.02
    drot = 0.1
    i = 113.1167
    return (L, KP, KI, R, K, drot, i)
end

function embrake_ivp(A, Tsample, ζ, x0)
    # state variables: [I, x, xe, xc]
    system = @system(x' = Ax)

    # initial condition
    I₀  = Singleton([0.0])
    x₀  = Singleton([0.0])
    xe₀ = Singleton([0.0])
    xc₀ = Singleton([0.0])
    X₀ = concretize(I₀ × x₀ × xe₀ × xc₀)

    # reset map
    Ar = sparse([1, 2, 3, 4, 4], [1, 2, 2, 2, 4], [1., 1., -1., -Tsample, 1.], 4, 4)
    br = sparsevec([3, 4], [x0, Tsample*x0], 4)
    reset_map(X) = Ar * X + br

    # hybrid system with clocked linear dynamics
    ha = HACLD1(system, reset_map, Tsample, ζ)
    return IVP(ha, X₀)
end

# model without parameter variation
function embrake_no_pv(; Tsample=1.E-4, ζ=1e-6, x0=0.05)
    L, KP, KI, R, K, drot, i = embrake_parameters()

    A = Matrix([-(R+K^2/drot)/L 0 KP/L KI/L;
                K/i/drot        0 0    0;
                0               0 0    0;
                0               0 0    0])

    return embrake_ivp(A, Tsample, ζ, x0)
end

# model with parameter variation changing only 1 coefficient
# corresponds to the Flow* settings in [SO15]
function embrake_pv_1(; Tsample=1.E-4, ζ=1e-6, Δ=3.0, x0=0.05)
    L, KP, KI, R, K, drot, i = embrake_parameters()
    p = 504. + (-Δ .. Δ)

    A = IntervalMatrix([-p       0 KP/L KI/L;
                        K/i/drot 0 0    0;
                        0        0 0    0;
                        0        0 0    0])

    return embrake_ivp(A, Tsample, ζ, x0)
end