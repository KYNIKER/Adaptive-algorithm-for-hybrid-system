using LazySets, LinearAlgebra, ReachabilityAnalysis, StaticArrays
include("Discretize.jl")
include("Utilities.jl")


function ReACTed_reachable_cell_t(system::EuclideanHybridSystem, p, granularity, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, clustering=true; jump_hull::Bool=false)
    return _ReACTed_reachable_cell_t(Grid(system.statespace, granularity), system, p, granularity, δ⁻, δ⁺, alg, maxOrder, reduceOrder, clustering; jump_hull)
end

# Function barrier: the body is compiled for the concrete Grid{N} and system type, so N can be used in type positions.
function _ReACTed_reachable_cell_t(grid::Grid{N}, system::EuclideanHybridSystem, p, granularity, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg, maxOrder::Int, reduceOrder::Int, clustering; jump_hull::Bool=false) where N
    @show foldl(*, grid.numCells), grid.numCells
    sweep_factor = 25
    TIME_FUNCS = true

    recorderTimeStart = time_ns()


    phiDict, tPhiDict, inputDict = phi_input_dict(system.flowMatrix, system.input, δ⁻, δ⁺, alg, maxOrder, reduceOrder)
    d = δ⁻
    #dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    isInvA = isinvertible(system.flowMatrix)
    A = system.flowMatrix
    #Φ = copy(phiDict[d])
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(system.flowMatrix)
    Φcache = system.flowMatrix == A_abs ? phiDict[d] : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)
    jump_reach = JumpReach(system, p, δ⁻, δ⁺, phiDict, tPhiDict, inputDict, P2A_abs, alg, maxOrder, reduceOrder; jump_hull)

    #frontierZonotopes = zonotopeArray[collect(CartesianIndex(c) for c in [gFrontierN; cFrontierN])]#zonotopeArray[collect(c.id for c in frontierN)]
    #push!(frontierZonotopes, zonotopeArray[40, 40, 1])
    #reachtimes = map(z -> propagate_set(z, [0.0, p], δ⁻, δ⁺, system, phiDict, tPhiDict, alg, maxOrder, reduceOrder), frontierZonotopes)

    #
    #   0: Didnt hit either constraint, guard or invariant
    #   1: Hit constraint
    #   2: Hit guard
    #   3: Hit invariant surface
    #
    reachable_by_action::ActionRelation{N} = ActionRelation{N}()
    reachable_by_flow::FlowRelation{N} = FlowRelation{N}()
    degenerate_dimensions = zeros(N) # Grid centres zero-width dimensions on their value, so cells need no extra offset
    invariant = hyperrectangle_to_HPolyhedron(system.statespace)
    U = system.input
    max_input = input_bound(system.flowMatrix, U, p, alg)
    max_flow = exp(p .* system.flowMatrix)
    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(system.globalConstraints[1])
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(system.edges[1].guard)
    invariantProjVectors, invariantProjBounds = getHalfSpaceProjections(invariant)

    hbox = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, N))
    guard_invariant_intersection = intersection(system.edges[1].guard, invariant)

    attempt_recorder = zeros(4)
    identity_matrix = diagm(ones(grid.dimension))

    fast_cells = falses(cellshape(grid))
    sweep = [min(sweep_factor, n) for n in grid.numCells] # cells per sweep block in each dimension (dimensions smaller than a block use one block)
    Z_big = remove_zero_generators(Zonotope(degenerate_dimensions, diagm(map((n, s) -> n == 1 ? 0. : granularity * (s/2), grid.numCells, sweep))))
    lower_offset_big = grid.lower .+ granularity .* (sweep ./ 2)
    generatorDict_big = ReACT_discretize_decomposed_generators(Z_big, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict)
    discretizationDict_big, finaldims_big = ReACT_discretize_combine_with_offsets(Z_big, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict_big)
    preallocated_center_big = zeros(grid.dimension)
    preallocated_genmat_big = zeros(grid.dimension, finaldims_big)
    preallocated_inter_genmat_big = zeros(grid.dimension, finaldims_big)
    ranges = ntuple(i -> 1:fld(grid.numCells[i], sweep[i]), N)
    #@show CartesianIndices((ranges...,))
    #solvr = LazySets.default_lp_solver(Float64)
    for bidx in CartesianIndices(ranges)

        of = lower_offset_big + ((car2vec(bidx) .- 1) .* granularity .* sweep)
        LazySets.API.translate!(Z_big, of)


        ReACT_discretize_combine_with_offset_vector!(preallocated_center_big, preallocated_genmat_big, preallocated_inter_genmat_big, discretizationDict_big, of, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict_big)
        #z, f = propagate_set(Z, [0.0, p], δ⁻, δ⁺, system, P2A_abs, phiDict, tPhiDict, inputDict, max_input, max_flow, invariant, alg, maxOrder, reduceOrder)

        f = propagate_set_check_hybrid_jump(Z_big, [0.0, p], δ⁻, δ⁺, system, system.edges, identity_matrix, attempt_recorder, grid.dimension, P2A_abs, guardProjVectors, guardProjBounds, invariantProjVectors, invariantProjBounds, discretizationDict_big, phiDict, tPhiDict, inputDict, max_input, max_flow, guard_invariant_intersection, alg, maxOrder, reduceOrder)

        if f == 0
            #@show collect(((b-1)*sweep_factor):(b*sweep_factor) for b in car2vec(bidx))
            fast_cells[ntuple(i->(1+(bidx[i]-1)*sweep[i]):(bidx[i]*sweep[i]), N)...] .= true
            #grid.deadCells[collect((1+(b-1)*sweep_factor):(b*sweep_factor) for b in car2vec(bidx))...] .= true

        end


        #grid.array[idx].sidx[1] = copy(f)
        LazySets.API.translate!(Z_big, -of)
    end
    #return grid, reachable_by_action, reachable_by_flow


    #@show degenerate_dimensions
    Z = remove_zero_generators(Zonotope(degenerate_dimensions, diagm(map(x -> x == 1 ? 0. : granularity/2, grid.numCells))))
    #dc = 0
    count = 0
    #@show grid.lower

    lower_offset = grid.lower .+ (granularity/2)
    generatorDict = ReACT_discretize_decomposed_generators(Z, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict)

    #discretizationDict = Dict()
    discretizationDict, finaldims = ReACT_discretize_combine_with_offsets(Z, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict)
    preallocated_center = zeros(grid.dimension)
    preallocated_genmat = zeros(grid.dimension, finaldims)
    preallocated_inter_genmat = zeros(grid.dimension, finaldims)
    #tarray = MArray{Tuple{grid.dimension}}(zeros(grid.dimension))
    #max_zonotope = Zonotope(MArray{Tuple{grid.dimension}}(max_input.center), SMatrix{grid.dimension,grid.dimension}(genmat(linear_map(max_flow, Z))))
    max_zonotope::Zonotope{Float64,Vector{Float64},Matrix{Float64}} = Zonotope(max_input.center, genmat(linear_map(max_flow, Z)))

    max_zonotope = minkowski_sum(max_zonotope, Zonotope(zeros(grid.dimension), diagm(fill(granularity/2, N))))
    # Fast cells all share max_zonotope's generators, so their cover test uses its precomputed H-representation.
    fast_cover = ConstantZonotopeCover(max_zonotope, Val(N))
    fast_center = zeros(N)


    for global_constraint in system.globalConstraints
        #guard_invariant_intersection = intersection(guard_invariant_intersection, global_constraint)
        for didx in get_touching_cell_idxs_l(grid, global_constraint)
            #@show didx
            grid.deadCells[didx] = true
        end
    end

    setup_time = time_ns()
    aux_time_recoder = time_ns()
    reachability_discretization_time = 0.0
    reachability_time = 0.0
    reachability_cover_time = 0.0
    act_set_time = 0.0

    for idx in CartesianIndices(grid.deadCells)
        if grid.deadCells[idx]
            #println("Skipping dead cell: ", idx)
            continue

            #else
            #@show "Processing cell: ", idx
        end

        if TIME_FUNCS
            aux_time_recoder = time_ns()
        end

        if !fast_cells[idx]
            count += 1
            of = lower_offset + ((car2vec(idx) .- 1) .* granularity)
            LazySets.API.translate!(Z, of)
            #@show Z, idx, of
            preallocated_center .= 0
            preallocated_genmat .= 0
            preallocated_inter_genmat .= 0

            ReACT_discretize_combine_with_offset_vector!(preallocated_center, preallocated_genmat, preallocated_inter_genmat, discretizationDict, of, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict)

            if TIME_FUNCS
                reachability_discretization_time += time_ns() - aux_time_recoder
            end

            z, f = propagate_set_jumps(jump_reach, Z, discretizationDict)

            if f != 1
                reachable_by_flow[idx] = get_touching_cell_idxs_lz(grid, z)
            else
                grid.deadCells[idx] = true
            end

            LazySets.API.translate!(Z, -of)
            if count % 1000 == 0
                println(count)
            end

            if TIME_FUNCS
                reachability_time += time_ns() - aux_time_recoder
            end
        else
            count += 1
            fast_cell_center!(fast_center, max_flow, max_zonotope.center, grid, idx)
            reachable_by_flow[idx] = get_touching_cell_idxs_cover(grid, fast_cover, fast_center)

            if count % 1000 == 0
                println(count)

            end
            if TIME_FUNCS
                reachability_cover_time += time_ns() - aux_time_recoder
            end
        end
    end

    if TIME_FUNCS
        aux_time_recoder = time_ns()
    end

    for act in system.Act
        #@show size(get_touching_cell_idxs_l(grid, act.guard))
        for idx in get_touching_cell_idxs_l(grid, act.guard)
            of = lower_offset + ((car2vec(idx) .- 1) .* granularity)
            LazySets.API.translate!(Z, of)

            reachable_by_action[(act.id, idx)] = get_touching_cell_idxs_lz(grid, minkowski_sum(linear_map(max_flow, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector)), max_input))

            if grid.deadCells[idx]
                grid.deadCells[idx] = false
            end
            LazySets.API.translate!(Z, -of)
        end
    end
    if TIME_FUNCS
        act_set_time += time_ns() - aux_time_recoder

        println("Total times: 
            Total: $((time_ns() - recorderTimeStart)/10^9)
            Setup: $((setup_time-recorderTimeStart) / 10^9)
            Reachability: $(reachability_time / 10^9)
                Discretization: $(reachability_discretization_time / 10^9)
            Covered reachability: $(reachability_cover_time / 10^9) 
            Total reachability: $((reachability_time+reachability_cover_time) / 10^9)
            Action relation: $(act_set_time / 10^9)")


    end

    # The offline shield over every cell, from the relations
    constraint = falses(cellshape(grid))
    for global_constraint in system.globalConstraints, c in get_touching_cell_idxs_l(grid, global_constraint)
        constraint[c] = true
    end
    s, store = relation_shield(constraint, reachable_by_action, reachable_by_flow, [a.id for a in system.Act])
    allowed, lost = cell_choices(cellshape(grid), length(system.Act), s, store_enabled(store), store_acted(store))
    grid.deadCells .= allowed .== 0
    return grid, allowed, lost, NamedTuple()
end

# TODO - Consider whether it is cheaper to check whether jumpset is part of fast_cells and does not need reachability analysis for the remaining interval. 
function propagate_jumpset(X0, interval, δ⁻::Float64, δ⁺::Float64, system, identity_matrix, attemptsRecorder, dims, Φ₂, guardProjVectors, guardProjBounds, PhiDict, TPhiDict, inputDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5)

    discretizationDict = newReACTDiscretizePlus(X0, δ⁻, δ⁺, system.flowMatrix, Φ₂, PhiDict, inputDict, alg, maxOrder, reduceOrder)
    time::Float64 = copy(minimum(interval))
    endtime::Float64 = maximum(interval)

    listOfEdges = system.edges


    zero_input = Zonotope(zeros(dims), zeros(size(genmat((inputDict[δ⁻])))))
    for edge in listOfEdges
        guard_projection_vectors = copy(guardProjVectors)
        guard_projection_bounds = copy(guardProjBounds)


        #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 
        reachtime, newSet, inputSet, flag = ReACT_guards_t(δ⁻, δ⁺, interval, identity_matrix, attemptsRecorder, zero_input, guard_projection_vectors, guard_projection_bounds, 2, PhiDict, TPhiDict, discretizationDict, inputDict)

        A_abs = abs.(system.flowMatrix)
        max_input = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, endtime - time, alg, false, nothing), system.input)

        #
        #   0: Didnt hit either constraint, guard or invariant
        #   1: Hit constraint
        #   2: Hit guard
        #   3: Hit invariant surface
        #

        if flag == 0

            A_abs = abs.(system.flowMatrix)
            max_input = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, endtime - time, alg, false, nothing), system.input)

            return (minkowski_sum(linear_map(exp((endtime - time) .* system.flowMatrix), X0), max_input), 0)


        elseif flag == 2
            return (minkowski_sum(newSet, inputSet), 0)

            intersectingSetsList = ReACT_touches_set_constant_input(newSet, inputSet, edge.guard, δ⁻, [reachtime, endtime], PhiDict, inputDict)
            tempIntersect = foldl(ConvexHull, intersectingSetsList)
            intersectedSet = convert(Zonotope, box_approximation(tempIntersect))
            jumpSet = linear_map(edge.jumpMatrix, zonotopeStripIntersection(intersectedSet, edge.guard))
            if !isnothing(edge.jumpVector)
                LazySets.translate!(jumpSet, edge.jumpVector)
            end

            return propagate_jumpset(jumpSet, [reachtime, endtime], δ⁻, δ⁺, system, identity_matrix, attemptsRecorder, dims, Φ₂, guardProjVectors, guardProjBounds, PhiDict, TPhiDict, inputDict, alg, maxOrder, reduceOrder)
        end

        return (minkowski_sum(newSet, inputSet), 1)
    end

end


function propagate_set_check_hybrid_jump(X0, interval, δ⁻::Float64, δ⁺::Float64, system, listOfEdges, identity_matrix, attemptsRecorder, dims, Φ₂, guardProjVectors, guardProjBounds, invariantProjVectors, invariantProjBounds, discretizationDict, PhiDict, TPhiDict, inputDict, max_input, max_flow, guard_intersection, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5)
    zero_input = Zonotope(zeros(dims), zeros(size(genmat((inputDict[δ⁻])))))
    for edge in listOfEdges
        guard_projection_vectors = copy(guardProjVectors)
        guard_projection_bounds = copy(guardProjBounds)


        #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 
        _, _, _, flag = ReACT_guards_t(δ⁻, δ⁺, interval, identity_matrix, attemptsRecorder, zero_input, guard_projection_vectors, guard_projection_bounds, 2, PhiDict, TPhiDict, discretizationDict, inputDict)


        #
        #   0: Didnt hit either constraint, guard or invariant
        #   1: Hit constraint
        #   2: Hit guard
        #   3: Hit invariant surface
        #
        if flag != 0
            return 1
        end
        if allunder(linear_map(max_flow, X0), max_input, invariantProjVectors, invariantProjBounds)
            continue
        else
            return 1
        end
    end
    return 0
end


# TODO - Optimize this! We know that the set will be intersecting initially, thus we can just use the fast disjoint check and actually use larger step sizes as long as we dont hit the constraint!!!
function ReACT_touches_set_constant_input(X0, U, set, δ⁻::Float64, interval, PhiDict, inputDiscretizationDict)
    # Note that in touches we always use δ⁻
    # That is, we do not adjust timestep sizes

    time::Float64 = copy(interval[1])
    endtime::Float64 = interval[2]

    intersectingSetsList = Vector{Zonotope}() # This will store all of our sets

    tempSet = minkowski_sum(X0, U)

    push!(intersectingSetsList, copy(tempSet))
    while time < endtime
        if !LazySets.API.isdisjoint(set, tempSet)
            push!(intersectingSetsList, copy(tempSet))

            U = linear_map(PhiDict[δ⁻], U) + inputDiscretizationDict[δ⁻]
            X0 = linear_map(PhiDict[δ⁻], X0)
            tempSet = minkowski_sum(X0, U)

            time += δ⁻
        else
            break
        end

    end

    return intersectingSetsList
end

function ReACT_touches_set_constant_input_b(X0, U, set, projection_vectors, projection_bounds, invariantProjVectors, invariantProjBounds, δ⁻::Float64, interval, PhiDict, inputDiscretizationDict)
    # Note that in touches we always use δ⁻
    # That is, we do not adjust timestep sizes

    time::Float64 = copy(interval[1])
    endtime::Float64 = interval[2]

    intersectingSetsList = Vector{Zonotope}() # This will store all of our sets

    tempSet = minkowski_sum(X0, U)

    push!(intersectingSetsList, copy(tempSet))
    while time < endtime

        if !allunder(tempSet, projection_vectors, projection_bounds) && someinside(tempSet, projection_vectors, projection_bounds) && someinside(tempSet, invariantProjVectors, invariantProjBounds) #touchesCheck(tempSet, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, invarientProjVectors, invarientProjBounds)

            push!(intersectingSetsList, copy(tempSet))

            U = linear_map(PhiDict[δ⁻], U) + inputDiscretizationDict[δ⁻]
            X0 = linear_map(PhiDict[δ⁻], X0)
            tempSet = minkowski_sum(X0, U)
            # i = i + 1
            time += δ⁻
        else
            break
        end

    end

    return intersectingSetsList
end


function ReACT_guards_t(δ⁻::Float64, δ⁺::Float64, interval, identity_matrix, attemptsRecorder, zero_input, guardProjVectors, guardProjBounds, STRATEGY::Integer, PhiDict, permutedphiDict, discritezationDict, inputDiscretizationDict)
    time::Float64 = copy(minimum(interval))
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(δ⁺)

    attemptsRecorder = attemptsRecorder .* 0
    origGuardProjVectors = copy(guardProjVectors)
    Φ::Matrix{Float64} = copy(identity_matrix)

    tempM = copy(identity_matrix)

    i = 1
    # The input's contribution so far. With time-varying input (inputs with generators) the contribution of a new step is
    # independent of what came before, so the two are a Minkowski sum. It stays lazy: someinside only needs its support
    # function, which is the sum of the parts' support functions. A constant input is composed exactly with the zonotope +.
    constant_input = iszero(genmat(inputDiscretizationDict[δ⁻]))
    input = copy(zero_input)

    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            # Handle if we can no longer reduce the reachset (we keep hitting something)
            if currentTimeStep >= δ⁻

                if fast_guards_check_t(discritezationDict[currentTimeStep], input, guardProjVectors, origGuardProjVectors, guardProjBounds)
                    approveFlag = true
                    if currentTimeStep + time < endtime

                        input = constant_input ? linear_map(PhiDict[currentTimeStep], input) + inputDiscretizationDict[currentTimeStep] :
                                MinkowskiSum(LinearMap(PhiDict[currentTimeStep], input), inputDiscretizationDict[currentTimeStep])

                        map!(x -> permutedphiDict[currentTimeStep] * x, guardProjVectors)

                        mul!(tempM, Φ, PhiDict[currentTimeStep])

                        copy!(Φ, tempM)

                    else
                        return (endtime, discritezationDict[δ⁻], input, 0)
                    end
                else

                    currentTimeStep /= 2
                    attempts += 1
                end

            else
                return (time, linear_map(Φ, discritezationDict[δ⁻]), input, 2)
            end

        end

        attemptsRecorder[(i%4)+1] = attempts
        i = i + 1
        time = time + currentTimeStep

        # Reset / apply strategy
        # Only do this if the current timestep is less than the initial
        if STRATEGY == 0
            # Only reduce
        elseif STRATEGY == 1
            # always try double
            if currentTimeStep < δ⁺
                currentTimeStep = currentTimeStep * 2
                #changedTimeStep = true
            end
        elseif STRATEGY == 2
            # If attemptsrecorder past 4 are successes, double timestep
            if currentTimeStep < δ⁺
                #lowest = min(4, i - 1)
                #window = @view attemptsRecorder[i-lowest:i-1]
                if all(==(1), attemptsRecorder)
                    currentTimeStep = currentTimeStep * 2
                    #changedTimeStep = true
                end

            end
        end
    end
    throw(error("ERROR!!! we should not have reached this point. "))
end

# This can be replaced with throwing an error. Currently we continue and just print
function handleHitConstraint(time, locationId)
    throw(error("ERROR!!! We have hit a constraint at loc: $(locationId) time: $time"))
    println("\nERROR!!!\n 
            ERROR!!!\n\n
            We have hit a constraint at loc: $locationId time: $time\n\n
            ERROR!!!\n
            ERROR!!!\n")
end


function guardCheck(input::Zonotope, Φ, X::Zonotope, constraint, guard, invariant) #; solver=model
    #
    #   0: Didnt hit either constraint, guard or invariant
    #   1: Hit constraint
    #   2: Hit guard
    #   3: Hit invariant surface
    #

    newR = minkowski_sum(input, linear_map(Φ, X))
    if LazySets.API.isdisjoint(newR, constraint)
        if LazySets.API.isdisjoint(newR, guard)
            return 0
        else
            return 2
        end

    else
        return 1
    end
end

function fast_guards_check(newR::Zonotope, co, constraintProjVectors, constraintProjBounds, go, guardProjVectors, guardProjBounds, io, invariantProjVectors, invariantProjBounds) #; solver=model
    if allunder(newR, co, constraintProjVectors, constraintProjBounds)
        if allunder(newR, go, guardProjVectors, guardProjBounds)
            return 0
        else
            return 2
        end

    else
        return 1
    end

    #return allunder(newR,c0, constraintProjVectors, constraintProjBounds) && someunder(newR, guardProjVectors, guardProjBounds) && someunder(newR, invarientProjVectors, invarientProjBounds)#all((input - ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))

end

function fast_guards_check_b(newR::Zonotope, input, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, invariantProjVectors, invariantProjBounds) #; solver=model
    #return all((input + ρ(x, newR)) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) && any((input + -ρ(-x, newR)) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) && all((input - ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))
    return !(someinside(newR, input, constraintProjVectors, constraintProjBounds) || someinside(newR, input, guardProjVectors, guardProjBounds)) && allunder(newR, input, invariantProjVectors, invariantProjBounds)

end
function fast_guards_check_t(newR::Zonotope, input, tguardProjVectors, guardProjVectors, guardProjBounds) #; solver=model
    #return all((input + ρ(x, newR)) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) && any((input + -ρ(-x, newR)) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) && all((input - ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))
    return !someinside(newR, input, tguardProjVectors, guardProjVectors, guardProjBounds)
    #return !someinside(minkowski_sum(newR, input), guardProjVectors, guardProjBounds)

end
function touchesCheck(newR::Zonotope, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, invarientProjVectors, invarientProjBounds) #; solver=model
    return allunder(newR, constraintProjVectors, constraintProjBounds) && someinside(newR, guardProjVectors, guardProjBounds) && someinside(newR, invarientProjVectors, invarientProjBounds)#all((input - ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))
end

function tsupfunc(d, Z, Ab, c, G)
    mul!(Z, G, d)
    Ab = sum(abs, Z)
    return dot(c, d) + Ab
end

function allunder(set, offset, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(offset, dir, bound)
        res = res && (offset[i] + ρ(dir[i], set) <= bound[i])
    end
    return res
end

function someunder(set, offset, dir, bound)
    res = !LazySets.API.isdisjoint(set, HPolyhedron(collect(LazySets.HalfSpace(a, b+c) for (a, b, c) in zip(dir, bound, offset))))
    return res
end

function someoutside(set, offset, dir, bound)
    return !allunder(set, offset, dir, bound)
end

function allunder(set, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (ρ(dir[i], set) <= bound[i])
    end
    return res
end

function allunder(set, input::LazySet, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (ρ(dir[i], set) + ρ(dir[i], input) <= bound[i])
    end
    return res
end

function someunder(set, dir, bound)
    return !LazySets.API.isdisjoint(set, HPolyhedron(collect(LazySets.HalfSpace(a, b) for (a, b) in zip(dir, bound))))
end

function someoutside(set, input::LazySet, dir, bound)
    res = false
    @inbounds @simd for i in eachindex(dir, bound)
        res = res || (-ρ(-dir[i], set)-ρ(-dir[i], input) > bound[i])
    end
    return res
end

function someinside(set, input::LazySet, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-dir[i], set)+(-ρ(-dir[i], input)) <= bound[i])
    end
    return res
end

function someinside(set, input::LazySet, tdir, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-tdir[i], set)+(-ρ(-dir[i], input)) <= bound[i])
    end
    return res
end

function someinside(set, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-dir[i], set) <= bound[i])
    end
    return res
end

function someinside_transparent(set, input::LazySet, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-dir[i], set)+(-ρ(-dir[i], input)) <= bound[i])
        @show (-ρ(-dir[i], set), -ρ(-dir[i], input), -ρ(-dir[i], set) - ρ(-dir[i], input), bound[i])
    end
    return res
end


include("ReACTedPropagate.jl")
include("ReACTedShieldingFlat.jl")
