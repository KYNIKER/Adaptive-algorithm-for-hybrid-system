# Uses ReACT_guards_t, someinside, allunder (guards.jl) and getHalfSpaceProjections.
#
# Reach set of a cell over one period of an EuclideanHybridSystem, with jumps. It replaces the jump branch of
# propagate_set_t, which flowed the whole crossing set for the longest remaining time, dropped the states that had not
# jumped by the end of the period, and stopped at a second touch. Here every piece of δ⁻ that touches a guard jumps and
# flows for its own remaining time, the states that have not jumped are kept, and further jumps are followed. This is
# modal_reach! (modal.jl) for a single location.

dense_zonotope(Z::Zonotope) = Zonotope(Vector{Float64}(LazySets.center(Z)), Matrix{Float64}(genmat(Z)))

struct JumpReach
    phiDict::Dict{Float64,Matrix{Float64}}
    tPhiDict::Dict{Float64,Matrix{Float64}}
    inputDict::Dict{Float64,Zonotope}
    A::Matrix{Float64}
    P2A_abs::Matrix{Float64}
    Φ::Vector{Matrix{Float64}} # Φ[k+1] = exp(k δ⁻ A)
    V::Vector{Zonotope{Float64,Vector{Float64},Matrix{Float64}}} # V[k+1]: the input's contribution over k δ⁻
    zero_input::Zonotope{Float64,Vector{Float64},Matrix{Float64}}
    guards::Vector{Tuple{Vector{Vector{Float64}},Vector{Float64}}} # per edge
    jumps::Vector{Tuple{Matrix{Float64},Vector{Float64}}} # per edge: jump matrix and vector
    exit_faces::Tuple{Vector{Vector{Float64}},Vector{Float64}} # state-space faces a ⋅ x <= b whose complement is not a halfspace of a guard
    steps::Int # δ⁻ steps per period
    δ⁻::Float64
    δ⁺::Float64
    alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg
    maxOrder::Int
    reduceOrder::Int
    max_depth::Int # jumps in a row before the transition is assumed bad
    jump_hull::Bool # jump the hull of the touching pieces once per guard, instead of each piece on its own
    identity_matrix::Matrix{Float64}
end

dense_projections(P) = ((v, b)=getHalfSpaceProjections(P); ([Vector{Float64}(a) for a in v], Vector{Float64}(b)))

function JumpReach(system::EuclideanHybridSystem, p, δ⁻::Float64, δ⁺::Float64, phiDict, tPhiDict, inputDict, P2A_abs, alg, maxOrder::Int, reduceOrder::Int; max_depth::Int=16, jump_hull::Bool=false)
    steps = round(Int, p / δ⁻)
    isapprox(steps * δ⁻, p; atol=1e-9) || throw(ArgumentError("The period $p must be a multiple of δ⁻ = $δ⁻"))
    A = Matrix{Float64}(system.flowMatrix)
    n = size(A, 1)
    U = system.input
    Φ = Matrix{Float64}[Matrix{Float64}(I, n, n)]
    V = [Zonotope(zeros(n), zeros(n, n))]
    for k in 1:steps
        d = k * δ⁻
        push!(Φ, exp(d .* A))
        push!(V, input_bound(A, U, d, alg)) # the input's contribution over k δ⁻, which may vary in time
    end
    guards = [dense_projections(e.guard) for e in system.edges]
    jumps = [(Matrix{Float64}(e.jumpMatrix), e.jumpVector === nothing ? zeros(n) : Vector{Float64}(e.jumpVector)) for e in system.edges]
    # A face a ⋅ x <= b is bordered by a guard with the halfspace -a ⋅ x <= -b, such as p >= 0 by the bounce guard's p <= 0.
    # This assumes a state can only cross a bordered face inside that guard, as for the ball.
    statespace = dense_projections(hyperrectangle_to_HPolyhedron(system.statespace))
    bordered(a, b) = any(((v, w),) -> any(i -> isapprox(v[i], -a; atol=1e-12) && isapprox(w[i], -b; atol=1e-12), eachindex(v)), guards)
    keep = [i for i in eachindex(statespace[1]) if !bordered(statespace[1][i], statespace[2][i])]
    exit_faces = (statespace[1][keep], statespace[2][keep])
    zero_input = Zonotope(zeros(n), zeros(size(genmat(inputDict[δ⁻]))))
    return JumpReach(phiDict, tPhiDict, inputDict, A, Matrix{Float64}(P2A_abs), Φ, V, zero_input, guards, jumps, exit_faces, steps, δ⁻, δ⁺, alg, maxOrder, reduceOrder, max_depth, jump_hull, Matrix{Float64}(I, n, n))
end

# Over-approximates the part of Z in the polyhedron with one strip intersection per halfspace that cuts it.
# Returns nothing if Z misses a halfspace.
function clip_to_polyhedron(Z, polyhedron)
    res = Z
    for (a, b) in zip(polyhedron...)
        if ρ(a, res) > b
            lo = -ρ(-a, res)
            lo > b && return nothing
            res = zonotopeStripIntersection(res, a, (lo + b) / 2, (b - lo) / 2)
        end
    end
    return res
end

function jump_set(Z, M::Matrix{Float64}, j::Vector{Float64}, maxOrder::Int, reduceOrder::Int)
    res = remove_zero_generators(Zonotope(M * LazySets.center(Z) + j, M * genmat(Z)))
    size(genmat(res), 2) == 0 && return Zonotope(LazySets.center(res), zeros(length(j), 1))
    if maxOrder > 0 && LazySets.order(res) > maxOrder
        res = reduce_order(res, reduceOrder)
    end
    return dense_zonotope(res)
end

# Appends to out sets covering what X reaches at the end of the period. X starts with between D and D + W steps of δ⁻
# left, depending on the state (W > 0 after a jump). disc is its discretization up to jr.δ⁺, or nothing to compute it here,
# only up to the largest step that fits in the time left. Returns false if the transition is bad: a set at the end crosses
# a state-space face that no guard borders, or the jumps nest deeper than max_depth.
# With jr.jump_hull, the parts of the pieces inside a guard are hulled and jump once, with the remaining-time window of
# all of them, instead of each piece jumping with its own window. Fewer flows, but a larger set over a longer window.
function reach_jumps!(out::Vector{Zonotope{Float64,Vector{Float64},Matrix{Float64}}}, jr::JumpReach, X, disc, D::Int, W::Int, depth::Int=0)
    depth > jr.max_depth && return false
    total = D + W
    if total == 0
        push!(out, X)
        return true
    end
    δ⁺ = jr.δ⁺
    if disc === nothing
        δ⁺ = min(jr.δ⁺, jr.δ⁻ * 2^floor(Int, log2(total)))
        disc = newReACTDiscretizePlus(X, jr.δ⁻, δ⁺, jr.A, jr.P2A_abs, jr.phiDict, jr.inputDict, jr.alg, jr.maxOrder, jr.reduceOrder)
    end
    # First step at which the flowpipe touches a guard
    k0 = total
    for (vectors, bounds) in jr.guards
        reachtime, _, _, flag = ReACT_guards_t(jr.δ⁻, δ⁺, [0.0, total * jr.δ⁻], jr.identity_matrix, zeros(4), jr.zero_input, copy(vectors), bounds, 2, jr.phiDict, jr.tPhiDict, disc, jr.inputDict)
        flag != 0 && (k0 = min(k0, round(Int, reachtime / jr.δ⁻)))
    end
    # jump_hull: per guard, the hull of the clipped pieces and the first and last step they came from
    hulls = jr.jump_hull ? Union{Nothing,Zonotope{Float64,Vector{Float64},Matrix{Float64}}}[nothing for _ in jr.guards] : nothing
    kfirst = jr.jump_hull ? fill(total, length(jr.guards)) : nothing
    klast = jr.jump_hull ? fill(-1, length(jr.guards)) : nothing
    all_jumped = false
    # Pieces from the first touch on, and the pieces that end the period when the time left varies
    #push!(out, minkowski_sum(linear_map(jr.Φ[k0], disc[jr.δ⁻]), jr.V[k0]))
    for k in min(k0, W==0 ? total : D):(total-1)
        piece = minkowski_sum(linear_map(jr.Φ[k+1], disc[jr.δ⁻]), jr.V[k+1])
        if k >= k0
            for (gi, (guard, (M, j))) in enumerate(zip(jr.guards, jr.jumps))
                someinside(piece, guard...) || continue
                clipped = clip_to_polyhedron(piece, guard)
                clipped === nothing && continue
                all_jumped = all_jumped || allunder(piece, guard...)
                if jr.jump_hull
                    #h = hulls[gi] === nothing ? clipped : zonotope_hull(hulls[gi], clipped)
                    h = k == k0 ? zonotope_hull(minkowski_sum(linear_map(jr.Φ[k0], disc[jr.δ⁻]), jr.V[k0]), clipped) : zonotope_hull(hulls[gi], clipped)
                    jr.maxOrder > 0 && LazySets.order(h) > jr.maxOrder && (h = reduce_order(h, jr.reduceOrder))
                    hulls[gi] = dense_zonotope(h)
                    kfirst[gi] = min(kfirst[gi], k)
                    klast[gi] = k
                    continue
                end
                # A state in this piece has spent between k and k + 1 steps, and had between D and D + W left at the start
                D′ = max(D - k - 1, 0)
                W′ = total - k - D′
                jumped = jump_set(clipped, M, j, jr.maxOrder, jr.reduceOrder)
                reach_jumps!(out, jr, jumped, nothing, D′, W′, depth + 1) || return false
            end
            all_jumped && break
        end
        W > 0 && k >= D && push!(out, piece)
    end
    if jr.jump_hull
        #push!(out, collect(clip_to_polyhedron(hulls[gi], guard) for (gi, guard) in enumerate(jr.guards) if hulls[gi] !== nothing)...)
        for (gi, (M, j)) in enumerate(jr.jumps)
            hulls[gi] === nothing && continue
            # A state in the hull has spent between kfirst and klast + 1 steps
            D′ = max(D - klast[gi] - 1, 0)
            W′ = total - kfirst[gi] - D′
            #@show D′, W′, kfirst[gi], klast[gi], total
            #push!(out, hulls[gi])
            # TODO - Determine whether the hull needs to get bloated by Phi2 
            reach_jumps!(out, jr, jump_set(zonotope_hull(minkowski_sum(linear_map(jr.Φ[k0], X), jr.V[k0]), minkowski_sum(linear_map(jr.Φ[klast[gi]+1], X), jr.V[klast[gi]+1])), M, j, jr.maxOrder, jr.reduceOrder), nothing, total-kfirst[gi]-1, 1, depth + 1) || return false
            #reach_jumps!(out, jr, jump_set(hulls[gi], M, j, jr.maxOrder, jr.reduceOrder), nothing, D′+W′-1, 1, depth + 1) || return false
            #reach_jumps!(out, jr, jump_set(hulls[gi], M, j, jr.maxOrder, jr.reduceOrder), nothing, total-kfirst[gi]-1, 1, depth + 1) || return false
        end
        #@show size(out)
    end
    all_jumped && return true
    # The sets at the shortest and longest time left are flowed without jumps, so states that jumped appear in them as if
    # they had not. Beyond a face that a guard borders are only such states (and the bloat of the sets around the jump).
    # Beyond any other face, a part inside a guard has jumped too; any other part has left the state space.
    for k in unique((D, total))
        final = minkowski_sum(linear_map(jr.Φ[k+1], X), jr.V[k+1])
        for (a, b) in zip(jr.exit_faces...)
            ρ(a, final) > b || continue
            beyond = clip_to_polyhedron(final, ([-a], [-b]))
            beyond === nothing && continue
            any(guard -> allunder(beyond, guard...), jr.guards) || return false
        end
    end
    W == 0 && push!(out, minkowski_sum(linear_map(jr.Φ[D+1], X), jr.V[D+1]))
    return true
end

# Sets covering the reach set of X at the end of the period, and the flag of propagate_set_t: 0 if the transition is safe,
# 1 if it is bad. disc is X's discretization.
function propagate_set_jumps(jr::JumpReach, X, disc)
    out = Zonotope{Float64,Vector{Float64},Matrix{Float64}}[]
    ok = reach_jumps!(out, jr, X, disc, jr.steps, 0)
    return out, ok ? 0 : 1
end

# Cells touching any of the sets, each once. With covers, the covers are taken from that cache.
function get_touching_cell_idxs_lz(grid::Grid{N}, sets::Vector{Zonotope{Float64,Vector{Float64},Matrix{Float64}}}, covers::Union{Nothing,CoverCache{N}}=nothing) where N
    touching(Z) = covers === nothing ? get_touching_cell_idxs_lz(grid, Z) : get_touching_cell_idxs_lz(grid, Z, covers)
    length(sets) == 1 && return touching(sets[1])
    cells = CartesianIndex{N}[]
    for Z in sets
        append!(cells, touching(Z))
    end
    return unique!(cells)
end

# Reachability data that a shield computes once and every cell computation reads: the discretization dictionaries, the
# whole-period flow and input used by fast cells, the projections propagate_set_check_hybrid_jump needs for the
# fast-cell sweep, and the JumpReach of propagate_set_jumps.
function reach_setup(system::EuclideanHybridSystem, p, δ⁻::Float64, δ⁺::Float64, alg, maxOrder::Int, reduceOrder::Int, ::Val{N}; jump_hull::Bool=false) where N
    phiDict, tPhiDict, inputDict = phi_input_dict(system.flowMatrix, system.input, δ⁻, δ⁺, alg, maxOrder, reduceOrder)
    A = system.flowMatrix
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    Φcache = A == A_abs ? phiDict[δ⁻] : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isinvertible(A), Φcache)
    invariant = hyperrectangle_to_HPolyhedron(system.statespace)
    U = system.input
    max_input = input_bound(A, U, p, alg)
    max_flow = exp(p .* A)
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(system.edges[1].guard)
    invariantProjVectors, invariantProjBounds = getHalfSpaceProjections(invariant)
    guard_invariant_intersection = intersection(system.edges[1].guard, invariant)
    return (; p=Float64(p), δ⁻, δ⁺, alg, maxOrder, reduceOrder, phiDict, tPhiDict, inputDict, A, P2A_abs, U, max_input, max_flow,
        guardProjVectors, guardProjBounds, invariantProjVectors, invariantProjBounds, guard_invariant_intersection,
        identity_matrix=diagm(ones(N)), attempt_recorder=zeros(4), center=zeros(N),
        jump_reach=JumpReach(system, p, δ⁻, δ⁺, phiDict, tPhiDict, inputDict, P2A_abs, alg, maxOrder, reduceOrder; jump_hull))
end
