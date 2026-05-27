using ReachabilityAnalysis, SparseArrays
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector

function loadSpacecraft(; abort_time::Union{Float64,Vector{Float64}}=-1.)
    # abort times cycle betweeen:
    # -1 (no abort)
    # 120
    # [120, 125]
    # [120, 145]
    # 240



    # variables
    x = 1  # x position (negative!)
    y = 2  # y position (negative!)
    vx = 3  # x velocity
    vy = 4  # y velocity
    t = 5  # time

    # number of variables
    n = 4 + 1

    # flag for activating the "rendezvous abort" scenario
    t_abort_lower = 0
    t_abort_upper = 0
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
    println("Aborting: $aborting")
    println("t_abort_lower = $t_abort_lower")
    println("t_abort_upper = $t_abort_upper")
    # discrete structure (graph)
    #automaton = GraphAutomaton(aborting ? 3 : 2)

    # common vector of affine dynamics
    b = SingleEntryVector(t, n, 1.)

    # mode 1 ("approaching")
    A1 = spzeros(n, n)
    A1[x, vx] = 1.
    A1[y, vy] = 1.
    A1[vx, x] = -0.057599765881773
    A1[vx, y] = 0.000200959896519766
    A1[vx, vx] = -2.89995083970656
    A1[vx, vy] = 0.00877200894463775
    A1[vy, x] = -0.000174031357370456
    A1[vy, y] = -0.0665123984901026
    A1[vy, vx] = -0.00875351105536225
    A1[vy, vy] = -2.90300269286856
    invariant1 = HPolyhedron([LazySets.HalfSpace(SingleEntryVector(x, n, 1.), -100.)])  # x <= -100
    if aborting
        invariant1 = HPolyhedron([
            LazySets.HalfSpace(SingleEntryVector(x, n, 1.), -100.),  # x <= -100
            LazySets.HalfSpace(SingleEntryVector(t, n, 1.), t_abort_upper)  # t <= t_abort_upper
        ])
    end

    #m_1 = @system(x' = Ax + b, x ∈ invariant)

    # mode 2 ("rendezvous attempt")
    A2 = spzeros(n, n)
    A2[x, vx] = 1.
    A2[y, vy] = 1.
    A2[vx, x] = -0.575999943070835
    A2[vx, y] = 0.000262486079431672
    A2[vx, vx] = -19.2299795908647
    A2[vx, vy] = 0.00876275931760007
    A2[vy, x] = -0.000262486080737868
    A2[vy, y] = -0.575999940191886
    A2[vy, vx] = -0.00876276068239993
    A2[vy, vy] = -19.2299765959399
    invariant2 = HPolyhedron([
        LazySets.HalfSpace(sparsevec([x], [-1.], n), 100.),           # x >= -100
        LazySets.HalfSpace(sparsevec([x], [1.], n), 100.),            # x <= 100
        LazySets.HalfSpace(sparsevec([y], [-1.], n), 100.),           # y >= -100
        LazySets.HalfSpace(sparsevec([y], [1.], n), 100.),            # y <= 100
        LazySets.HalfSpace(sparsevec([x, y], [-1., -1.], n), 141.1),  # x + y >= -141.1
        LazySets.HalfSpace(sparsevec([x, y], [1., 1.], n), 141.1),    # x + y <= 141.1
        LazySets.HalfSpace(sparsevec([x, y], [1., -1.], n), 141.1),   # -x + y >= -141.1
        LazySets.HalfSpace(sparsevec([x, y], [-1., 1.], n), 141.1)    # -x + y <= 141.1
    ])

    if aborting
        invariant2 = HPolyhedron([
            LazySets.HalfSpace(sparsevec([x], [-1.], n), 100.),           # x >= -100
            LazySets.HalfSpace(sparsevec([x], [1.], n), 100.),            # x <= 100
            LazySets.HalfSpace(sparsevec([y], [-1.], n), 100.),           # y >= -100
            LazySets.HalfSpace(sparsevec([y], [1.], n), 100.),            # y <= 100
            LazySets.HalfSpace(sparsevec([x, y], [-1., -1.], n), 141.1),  # x + y >= -141.1
            LazySets.HalfSpace(sparsevec([x, y], [1., 1.], n), 141.1),    # x + y <= 141.1
            LazySets.HalfSpace(sparsevec([x, y], [1., -1.], n), 141.1),   # -x + y >= -141.1
            LazySets.HalfSpace(sparsevec([x, y], [-1., 1.], n), 141.1),    # -x + y <= 141.1
            LazySets.HalfSpace(SingleEntryVector(t, n, 1.), t_abort_upper)  # t <= t_abort_upper
        ])
    end



    #m_2 = @system(x' = Ax + b, x ∈ invariant)

    # mode 3 ("aborting")
    A3 = spzeros(n, n)
    A3[x, vx] = 1.
    A3[y, vy] = 1.
    A3[vx, x] = 0.0000575894721132
    A3[vx, vy] = 0.00876276
    A3[vy, vx] = -0.00876276

    #A3 = A
    invariant3 = nothing

    #m_3 = @system(x' = Ax + b, x ∈ invariant)

    # modes
    #modes = aborting ? [m_1, m_2, m_3] : [m_1, m_2]

    locations = Vector{Location}()

    loc1edges = Vector{Edge}()
    loc2edges = Vector{Edge}()
    loc3edges = Vector{Edge}()

    # transition 1 -> 2
    guard = HPolyhedron([
        LazySets.HalfSpace(sparsevec([x], [-1.], n), 100.),           # x >= -100
        LazySets.HalfSpace(sparsevec([x], [1.], n), 100.),            # x <= 100
        LazySets.HalfSpace(sparsevec([y], [-1.], n), 100.),           # y >= -100
        LazySets.HalfSpace(sparsevec([y], [1.], n), 100.),            # y <= 100
        LazySets.HalfSpace(sparsevec([x, y], [-1., -1.], n), 141.1),  # x + y >= -141.1
        LazySets.HalfSpace(sparsevec([x, y], [1., 1.], n), 141.1),    # x + y <= 141.1
        LazySets.HalfSpace(sparsevec([x, y], [1., -1.], n), 141.1),   # -x + y >= -141.1
        LazySets.HalfSpace(sparsevec([x, y], [-1., 1.], n), 141.1)    # -x + y <= 141.1
    ])
    #t1 = ConstrainedIdentityMap(n, guard)

    push!(loc1edges, Edge(2, guard, Diagonal(ones(n)), zeros(n)))

    if aborting
        # 1 -> 3
        #add_transition!(automaton, 1, 3, 2)
        guard = LazySets.HalfSpace(SingleEntryVector(t, n, -1.), -t_abort_lower)  # t >= t_abort_lower

        # TODO put this back after testing
        # TODO put this back after testing
        push!(loc1edges, Edge(3, guard, Diagonal(ones(n)), zeros(n)))
        # TODO put this back after testing
        # TODO put this back after testing

        #t2 = ConstrainedIdentityMap(n, guard)

        # 2 -> 3

        #add_transition!(automaton, 2, 3, 3)
        guard = LazySets.HalfSpace(SingleEntryVector(t, n, -1.), -t_abort_lower)  # t >= t_abort_lower

        push!(loc2edges, Edge(3, guard, Diagonal(ones(n)), zeros(n)))
    end

    velocity = 0.055 * 60.0     # meters per minute
    cx = velocity * cos(π / 8)  # x-coordinate of the octagon's first (ENE) corner
    cy = velocity * sin(π / 8)  # y-coordinate of the octagon's first (ENE) corner

    # cx        -> 3.048802457287246
    # cx + cy   -> 4.311657784092042e

    loc2Constraint = [
        # Line of sight property
        #LazySets.HalfSpace(sparsevec([x], [-1.0], n), 100.0),             # x >= -100
        LazySets.HalfSpace(sparsevec([x, y], [tan(π / 6), -1.0], n), 0.0),  # -x tan(30°) + y >= 0
        LazySets.HalfSpace(sparsevec([x, y], [tan(π / 6), 1.0], n), 0.0),   # -x tan(30°) - y >= 0
        # Velocity / octagon property
        LazySets.HalfSpace(sparsevec([vx], [-1.0], n), cx),                # vx >= -cx
        LazySets.HalfSpace(sparsevec([vx], [1.0], n), cx),                 # vx <= cx
        LazySets.HalfSpace(sparsevec([vy], [-1.0], n), cx),                # vy >= -cx
        LazySets.HalfSpace(sparsevec([vy], [1.0], n), cx),                 # vy <= cx
        LazySets.HalfSpace(sparsevec([vx, vy], [1., 1.0], n), cy + cx),    # vx + vy <= cy + cx
        LazySets.HalfSpace(sparsevec([vx, vy], [1., -1.0], n), cy + cx),   # vx - vy <= cy + cx
        LazySets.HalfSpace(sparsevec([vx, vy], [-1., 1.0], n), cy + cx),   # -vx + vy <= cy + cx
        LazySets.HalfSpace(sparsevec([vx, vy], [-1., -1.0], n), cy + cx)   # -vx - vy <= cy + cx
    ]


    #target = BallInf(zeros(2), 0.2)
    loc3Constraint = [
        LazySets.HalfSpace(sparsevec([y], [1.], n), -0.2)   # y <= -0.2
    ]



    push!(locations, Location(1, invariant1, A1, nothing, nothing, b, loc1edges, []))
    push!(locations, Location(2, invariant2, A2, nothing, nothing, b, loc2edges, loc2Constraint))

    if aborting
        push!(locations, Location(3, invariant3, A3, nothing, nothing, b, loc3edges, loc3Constraint))
    end

    #Constraint 
    properties = []
    # verify that specification holds


    H = HybridSystemV2(locations, properties)


    # initial condition in mode 1
    X0 = Hyperrectangle([-900., -400., 0., 0., 0.],
        [25., 25., 0., 0., 0.])

    X0 = convert(Zonotope, X0)


    T = 300.

    return H, 1, X0, T
end

# loadSpacecraft(abort_time=[120., 145.])