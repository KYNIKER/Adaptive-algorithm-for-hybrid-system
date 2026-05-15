using ReachabilityAnalysis, SparseArrays
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector

function spacecraft(; abort_time::Union{Float64, Vector{Float64}}=-1.)
    # variables
    x  = 1  # x position (negative!)
    y  = 2  # y position (negative!)
    vx = 3  # x velocity
    vy = 4  # y velocity
    t  = 5  # time

    # number of variables
    n = 4 + 1

    # flag for activating the "rendezvous abort" scenario
    aborting = true
    if abort_time isa Number
        if abort_time < 0.
            aborting = false
        else
            t_abort_lower = abort_time
            t_abort_upper = abort_time
        end
    else
        if length(abort_time) != 2
            error("abort time must be a point or an interval")
        end
        t_abort_lower = abort_time[1]
        t_abort_upper = abort_time[2]
    end

    # discrete structure (graph)
    automaton = GraphAutomaton(aborting ? 3 : 2)

    # common vector of affine dynamics
    b = SingleEntryVector(t, n, 1.)

    # mode 1 ("approaching")
    A = spzeros(n, n)
    A[x, vx] = 1.
    A[y, vy] = 1.
    A[vx, x] = -0.057599765881773
    A[vx, y] = 0.000200959896519766
    A[vx, vx] = -2.89995083970656
    A[vx, vy] = 0.00877200894463775
    A[vy, x] = -0.000174031357370456
    A[vy, y] = -0.0665123984901026
    A[vy, vx] = -0.00875351105536225
    A[vy, vy] = -2.90300269286856
    invariant = HalfSpace(SingleEntryVector(x, n, 1.), -100.)  # x <= -100
    if aborting
        invariant = HPolyhedron([
            invariant,
            HalfSpace(SingleEntryVector(t, n, 1.), t_abort_upper)  # t <= t_abort_upper
           ])
    end
    m_1 = @system(x' = Ax + b, x ∈ invariant)

    # mode 2 ("rendezvous attempt")
    A = spzeros(n, n)
    A[x, vx] = 1.
    A[y, vy] = 1.
    A[vx, x] = -0.575999943070835
    A[vx, y] = 0.000262486079431672
    A[vx, vx] = -19.2299795908647
    A[vx, vy] = 0.00876275931760007
    A[vy, x] = -0.000262486080737868
    A[vy, y] = -0.575999940191886
    A[vy, vx] = -0.00876276068239993
    A[vy, vy] = -19.2299765959399
    invariant = HPolyhedron([
        HalfSpace(sparsevec([x], [-1.], n), 100.),           # x >= -100
        HalfSpace(sparsevec([x], [1.], n), 100.),            # x <= 100
        HalfSpace(sparsevec([y], [-1.], n), 100.),           # y >= -100
        HalfSpace(sparsevec([y], [1.], n), 100.),            # y <= 100
        HalfSpace(sparsevec([x, y], [-1., -1.], n), 141.1),  # x + y >= -141.1
        HalfSpace(sparsevec([x, y], [1., 1.], n), 141.1),    # x + y <= 141.1
        HalfSpace(sparsevec([x, y], [1., -1.], n), 141.1),   # -x + y >= -141.1
        HalfSpace(sparsevec([x, y], [-1., 1.], n), 141.1)    # -x + y <= 141.1
       ])
    if aborting
        addconstraint!(invariant, HalfSpace(sparsevec([t], [1.], n), t_abort_upper))  # t <= t_abort_upper
    end

    m_2 = @system(x' = Ax + b, x ∈ invariant)

    # mode 3 ("aborting")
    A = spzeros(n, n)
    A[x, vx] = 1.
    A[y, vy] = 1.
    A[vx, x] = 0.0000575894721132
    A[vx, vy] = 0.00876276
    A[vy, vx] = -0.00876276
    invariant = Universe(n)

    m_3 = @system(x' = Ax + b, x ∈ invariant)

    # modes
    modes = aborting ? [m_1, m_2, m_3] : [m_1, m_2]

    # transition 1 -> 2
    add_transition!(automaton, 1, 2, 1)
    guard = HPolyhedron([
        HalfSpace(sparsevec([x], [-1.], n), 100.),           # x >= -100
        HalfSpace(sparsevec([x], [1.], n), 100.),            # x <= 100
        HalfSpace(sparsevec([y], [-1.], n), 100.),           # y >= -100
        HalfSpace(sparsevec([y], [1.], n), 100.),            # y <= 100
        HalfSpace(sparsevec([x, y], [-1., -1.], n), 141.1),  # x + y >= -141.1
        HalfSpace(sparsevec([x, y], [1., 1.], n), 141.1),    # x + y <= 141.1
        HalfSpace(sparsevec([x, y], [1., -1.], n), 141.1),   # -x + y >= -141.1
        HalfSpace(sparsevec([x, y], [-1., 1.], n), 141.1)    # -x + y <= 141.1
       ])
    t1 = ConstrainedIdentityMap(n, guard)

    if aborting
        add_transition!(automaton, 1, 3, 2)
        guard = HalfSpace(SingleEntryVector(t, n, -1.), -t_abort_lower)  # t >= t_abort_lower
        t2 = ConstrainedIdentityMap(n, guard)

        add_transition!(automaton, 2, 3, 3)
        guard = HalfSpace(SingleEntryVector(t, n, -1.), -t_abort_lower)  # t >= t_abort_lower
        t3 = ConstrainedIdentityMap(n, guard)
    end

    # transition annotations
    resetmaps = aborting ? [t1, t2, t3] : [t1]

    H = HybridSystem(automaton, modes, resetmaps, [AutonomousSwitching()])

    # initial condition in mode 1
    X0 = Hyperrectangle([-900., -400., 0., 0., 0.],
                        [25., 25., 0., 0., 0.])
    initial_condition = [(1, X0)]

    return IVP(H, initial_condition)
end