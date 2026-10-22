using LazySets, LinearAlgebra, ReachabilityAnalysis, StaticArrays
include("Discretize.jl")
include("Utilities.jl")

export ReACTedShieldingK


function ReACTed_reachable_cell(system::EuclideanHybridSystem, p, granularity, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, clustering=true)
    grid = Grid(system.statespace, granularity)
    @show foldl(*, grid.numCells), grid.numCells
    #unsafeDict = Dict{CartesianIndex,Vector{LazySet}}()

    #zonotopeArray = initialize_zonotope_array(grid)

    #mark_dead_cells!(grid, system.globalConstraints[1], unsafeDict)

    #
    #   Take union with frontier of guards( and Act?)
    #
    #unsafeCells, cFrontierCells = get_contained_edge_cells(grid, system.globalConstraints[1])
    #cFrontierN = grow_indices(grid, map(x -> car2vec(x.id), cFrontierCells))
    #@show frontierCells
    #guardCells, gFrontierCells = get_contained_edge_cells(grid, system.edges[1].guard)
    #gFrontierCells = get_contained_perimeter_cells(grid, system.edges[1].guard)
    #@show guardCells
    #@show gFrontierCells
    #gFrontierN = grow_indices(grid, map(x -> car2vec(x.id), gFrontierCells))
    #@show gFrontierN
    #phiDicts = Dict()
    #tphiDicts = Dict()

    #=
    for x in hybridSystem.locations
        pd, tpd, id = PhiInputDict(x, δ⁻, δ⁺, alg, maxOrder, reduceOrder)
        phiDicts[x.id] = pd
        tphiDicts[x.id] = tpd
        inputDicts[x.id] = id
        timeConstraintDict[x.id] = []
        constraintDict[x.id] = vcat(x.constraints, constraint)
    end
    =#

    #phiDict = PhiDict(system.flowMatrix, δ⁻, δ⁺, alg)
    #tPhiDict = Dict(collect((k, permutedims(copy(v))) for (k, v) in pairs(phiDict)))
    phiDict, tPhiDict, inputDict = PhiInputDict(system.flowMatrix, system.input, δ⁻, δ⁺, alg, maxOrder, reduceOrder)
    d = δ⁻
    #dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    isInvA = isinvertible(system.flowMatrix)
    #Φ = copy(phiDict[d])
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(system.flowMatrix)
    Φcache = system.flowMatrix == A_abs ? phiDict[d] : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)

    #frontierZonotopes = zonotopeArray[collect(CartesianIndex(c) for c in [gFrontierN; cFrontierN])]#zonotopeArray[collect(c.id for c in frontierN)]
    #push!(frontierZonotopes, zonotopeArray[40, 40, 1])
    #reachtimes = map(z -> propagate_set(z, [0.0, p], δ⁻, δ⁺, system, phiDict, tPhiDict, alg, maxOrder, reduceOrder), frontierZonotopes)

    #
    #   0: Didnt hit either constraint, guard or invariant
    #   1: Hit constraint
    #   2: Hit guard
    #   3: Hit invariant surface
    #

    reachable_by_action = Dict()
    reachable_by_flow = Dict()
    degenerate_dimensions = collect(u - l == l ? l - (granularity / 2) : 0. for (l, u) in zip(grid.lower, grid.upper))
    #@show degenerate_dimensions
    Z = remove_zero_generators(Zonotope(degenerate_dimensions, diagm(map(x -> x == 1 ? 0. : granularity/2, grid.numCells))))
    #dc = 0
    count = 0
    #@show grid.lower
    max_input = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, p, alg, false, nothing), system.input)
    max_flow = exp(p .* system.flowMatrix)
    lower_offset = grid.lower .+ (granularity/2)
    invariant = hyperrectangle_to_HPolyhedron(system.statespace)

    #constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    #guardProjVectors, guardProjBounds = getHalfSpaceProjections(edge.guard)
    #invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(invariant)

    for idx in CartesianIndices(grid.deadCells)
        count += 1
        of = lower_offset + ((car2vec(idx) .- 1) .* granularity)
        LazySets.API.translate!(Z, of)
        #@show Z, idx, of
        for act in system.Act
            if !LazySets.API.isdisjoint(act.guard, Z)
                reachable_by_action[(act, idx)] = get_touching_cell_idxs(grid, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector))
            end
        end

        z, f = propagate_set(Z, [0.0, p], δ⁻, δ⁺, system, P2A_abs, phiDict, tPhiDict, inputDict, max_input, max_flow, invariant, alg, maxOrder, reduceOrder)

        #=if !LazySets.API.isdisjoint(z, system.edges[1].guard)
            dc += 1
        end=#

        if f != 1
            #@show LazySets.API.isdisjoint(Z, system.globalConstraints[1]), Z, idx
            #grid.deadCells[idx] = true

            #println("so deads")
            reachable_by_flow[idx] = get_touching_cell_idxs_b(grid, z)
            #else
            #grid.array[idx].pCells = get_touching_cell_idxs(grid, z)
        end
        #grid.array[idx].sidx[1] = copy(f)
        LazySets.API.translate!(Z, -of)
        if count % 1000 == 0
            println(count)
        end
    end
    #@show dc
    return grid, reachable_by_action, reachable_by_flow

end



function ReACTed_reachable_cell_b(system::EuclideanHybridSystem, p, granularity, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, clustering=true)
    grid = Grid(system.statespace, granularity)
    @show foldl(*, grid.numCells), grid.numCells
    #unsafeDict = Dict{CartesianIndex,Vector{LazySet}}()

    #zonotopeArray = initialize_zonotope_array(grid)

    #mark_dead_cells!(grid, system.globalConstraints[1], unsafeDict)

    #
    #   Take union with frontier of guards( and Act?)
    #
    #unsafeCells, cFrontierCells = get_contained_edge_cells(grid, system.globalConstraints[1])
    #cFrontierN = grow_indices(grid, map(x -> car2vec(x.id), cFrontierCells))
    #@show frontierCells
    #guardCells, gFrontierCells = get_contained_edge_cells(grid, system.edges[1].guard)
    #gFrontierCells = get_contained_perimeter_cells(grid, system.edges[1].guard)
    #@show guardCells
    #@show gFrontierCells
    #gFrontierN = grow_indices(grid, map(x -> car2vec(x.id), gFrontierCells))
    #@show gFrontierN
    #phiDicts = Dict()
    #tphiDicts = Dict()

    #=
    for x in hybridSystem.locations
        pd, tpd, id = PhiInputDict(x, δ⁻, δ⁺, alg, maxOrder, reduceOrder)
        phiDicts[x.id] = pd
        tphiDicts[x.id] = tpd
        inputDicts[x.id] = id
        timeConstraintDict[x.id] = []
        constraintDict[x.id] = vcat(x.constraints, constraint)
    end
    =#

    #phiDict = PhiDict(system.flowMatrix, δ⁻, δ⁺, alg)
    #tPhiDict = Dict(collect((k, permutedims(copy(v))) for (k, v) in pairs(phiDict)))
    phiDict, tPhiDict, inputDict = phi_input_dict(system.flowMatrix, system.input, δ⁻, δ⁺, alg, maxOrder, reduceOrder)
    d = δ⁻
    #dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    isInvA = isinvertible(system.flowMatrix)
    A = system.flowMatrix
    #Φ = copy(phiDict[d])
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(system.flowMatrix)
    Φcache = system.flowMatrix == A_abs ? phiDict[d] : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)

    #frontierZonotopes = zonotopeArray[collect(CartesianIndex(c) for c in [gFrontierN; cFrontierN])]#zonotopeArray[collect(c.id for c in frontierN)]
    #push!(frontierZonotopes, zonotopeArray[40, 40, 1])
    #reachtimes = map(z -> propagate_set(z, [0.0, p], δ⁻, δ⁺, system, phiDict, tPhiDict, alg, maxOrder, reduceOrder), frontierZonotopes)

    #
    #   0: Didnt hit either constraint, guard or invariant
    #   1: Hit constraint
    #   2: Hit guard
    #   3: Hit invariant surface
    #

    reachable_by_action = Dict{Tuple{Int64,CartesianIndex},AbstractArray{CartesianIndex}}()
    reachable_by_flow = Dict{CartesianIndex,AbstractArray{CartesianIndex}}()
    degenerate_dimensions = collect(u - l == l ? l - (granularity / 2) : 0. for (l, u) in zip(grid.lower, grid.upper))
    #@show degenerate_dimensions
    Z = remove_zero_generators(Zonotope(degenerate_dimensions, diagm(map(x -> x == 1 ? 0. : granularity/2, grid.numCells))))
    #dc = 0
    count = 0
    #@show grid.lower
    U = system.input
    max_input = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, p, alg, false, nothing), U)
    max_flow = exp(p .* system.flowMatrix)
    lower_offset = grid.lower .+ (granularity/2)
    generatorDict = ReACT_discretize_decomposed_generators(Z, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict)
    invariant = hyperrectangle_to_HPolyhedron(system.statespace)
    #discretizationDict = Dict()
    discretizationDict, finaldims = ReACT_discretize_combine_with_offsets(Z, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict)
    preallocated_center = zeros(grid.dimension)
    preallocated_genmat = zeros(grid.dimension, finaldims)
    preallocated_inter_genmat = zeros(grid.dimension, finaldims)

    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(system.globalConstraints[1])
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(system.edges[1].guard)
    invariantProjVectors, invariantProjBounds = getHalfSpaceProjections(invariant)

    hbox = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, grid.dimension))
    guard_invariant_intersection = intersection(system.edges[1].guard, invariant)
    for global_constraint in system.globalConstraints
        #guard_invariant_intersection = intersection(guard_invariant_intersection, global_constraint)
        for didx in get_touching_cell_idxs_l(grid, global_constraint)
            grid.deadCells[didx] = true
        end
    end
    for idx in CartesianIndices(grid.deadCells)
        count += 1
        of = lower_offset + ((car2vec(idx) .- 1) .* granularity)
        LazySets.API.translate!(Z, of)
        #@show Z, idx, of


        ReACT_discretize_combine_with_offset_vector!(preallocated_center, preallocated_genmat, preallocated_inter_genmat, discretizationDict, of, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict)
        #z, f = propagate_set(Z, [0.0, p], δ⁻, δ⁺, system, P2A_abs, phiDict, tPhiDict, inputDict, max_input, max_flow, invariant, alg, maxOrder, reduceOrder)

        z, f = propagate_set_b(Z, [0.0, p], δ⁻, δ⁺, system, P2A_abs, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, discretizationDict, phiDict, tPhiDict, inputDict, max_input, max_flow, guard_invariant_intersection, invariant, invariantProjVectors, invariantProjBounds, alg, maxOrder, reduceOrder)
        #=
        if idx == CartesianIndex(100, 10)
            #@show propagate_set(Z, [0.0, p], δ⁻, δ⁺, system, P2A_abs, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, discretizationDict, phiDict, tPhiDict, inputDict, max_input, max_flow, guard_invariant_intersection, invariant, invariantProjVectors, invariantProjBounds, alg, maxOrder, reduceOrder)
            #@show invariant, invariantProjVectors
            @show f
            @show z
            @show minkowski_sum(linear_map(max_flow, Z), max_input)
            @show minkowski_sum(linear_map(exp((p-δ⁻) .* system.flowMatrix), Z), linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, p-δ⁻, alg, false, nothing), U))
            @show LazySets.API.isdisjoint(z, system.globalConstraints[1])
            @show LazySets.API.isdisjoint(z, system.edges[1].guard)
            @show LazySets.API.isdisjoint(z, invariant)
            @show get_touching_cell_idxs_b(grid, z)
            @show get_touching_cell_idxs_t(grid, z)
            @show get_touching_cell_idxs(grid, z)
            reachtime, newSet, inputSet, flag = ReACT_guards_b(δ⁻, δ⁺, [0.0, p], system.statespace, system.edges[1].guard, copy(guardProjVectors), copy(guardProjBounds), system.globalConstraints[1], copy(constraintProjVectors), copy(constraintProjBounds), invariant, copy(invariantProjVectors), copy(invariantProjBounds), 2, phiDict, tPhiDict, discretizationDict, inputDict)
            @show reachtime, newSet, inputSet, flag
            @show minkowski_sum(newSet, inputSet)
        end
        =#
        #=if !LazySets.API.isdisjoint(z, system.edges[1].guard)
            dc += 1
        end=#
        #if idx == CartesianIndex(2, 9)
        #    @show z, f
        #end
        if f != 1 #&& f != 3
            #@show LazySets.API.isdisjoint(Z, system.globalConstraints[1]), Z, idx
            #grid.deadCells[idx] = true

            #println("so deads")
            #touches_list = get_touching_cell_idxs(grid, z)
            #if !isempty(touches_list)
            #    reachable_by_flow[idx] = copy(touches_list)#get_touching_cell_idxs(grid, z)
            #end

            # Maybe check if flag is 2 and handle this seperatly. 
            #if f == 3
            #@show LazySets.API.isdisjoint(z, invariant)

            #grid.deadCells[idx] = true
            #=
            elseif f == 2
                intersectingSetsList = ReACT_touches_set_constant_input_b(newSet, inputSet, edge.guard, copy(guardProjVectors), copy(guardProjBounds), invariantProjVectors, invariantProjBounds, δ⁻, [reachtime, endtime], PhiDict, inputDict)


                tempIntersect = foldl(ConvexHull, intersectingSetsList)
                intersectedSet = convert(Zonotope, box_approximation(tempIntersect))
                jumpSet = linear_map(edge.jumpMatrix, zonotopeStripIntersection(intersectedSet, guard_intersection))
                if !isnothing(edge.jumpVector)
                    LazySets.translate!(jumpSet, edge.jumpVector)
                end

                return propagate_set(jumpSet, [reachtime, endtime], δ⁻, δ⁺, system, Φ₂, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, PhiDict, TPhiDict, inputDict, max_input, max_flow, invariant, invariantProjVectors, invariantProjBounds, alg, maxOrder, reduceOrder)
                =#
            #else
            reachable_by_flow[idx] = get_touching_cell_idxs_l(grid, z)

            #end




            #reachable_by_flow[idx] = get_touching_cell_idxs_b(grid, z)
            #else
            #grid.array[idx].pCells = get_touching_cell_idxs(grid, z)
        else
            grid.deadCells[idx] = true
        end
        #grid.array[idx].sidx[1] = copy(f)
        LazySets.API.translate!(Z, -of)
        if count % 1000 == 0
            println(count)
        end
    end

    #can_act = false
    for act in system.Act
        #@show size(get_touching_cell_idxs_l(grid, act.guard))
        for idx in get_touching_cell_idxs_l(grid, act.guard)
            of = lower_offset + ((car2vec(idx) .- 1) .* granularity)
            LazySets.API.translate!(Z, of)
            #touches_list = get_touching_cell_idxs(grid, z)
            #touches_list = get_touching_cell_idxs(grid, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector))
            #if !isempty(touches_list)
            #    reachable_by_action[(act, idx)] = copy(touches_list)
            #reachable_by_flow[idx] = get_touching_cell_idxs(grid, z)
            #end
            reachable_by_action[(act.id, idx)] = get_touching_cell_idxs_l(grid, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector))

            if act.id == 1
                #@show get_touching_cell_idxs(grid, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector)), LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector)
                #reachable_by_action[(act.id, idx)] = get_touching_cell_idxs(grid, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector))
            end
            if grid.deadCells[idx]
                grid.deadCells[idx] = false
            end
            LazySets.API.translate!(Z, -of)
        end
    end

    #@show dc
    return grid, reachable_by_action, reachable_by_flow

end


function ReACTed_reachable_cell_t(system::EuclideanHybridSystem, p, granularity, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, clustering=true)
    grid = Grid(system.statespace, granularity)
    @show foldl(*, grid.numCells), grid.numCells
    sweep_factor = 20

    phiDict, tPhiDict, inputDict = phi_input_dict(system.flowMatrix, system.input, δ⁻, δ⁺, alg, maxOrder, reduceOrder)
    d = δ⁻
    #dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    isInvA = isinvertible(system.flowMatrix)
    A = system.flowMatrix
    #Φ = copy(phiDict[d])
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(system.flowMatrix)
    Φcache = system.flowMatrix == A_abs ? phiDict[d] : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)

    #frontierZonotopes = zonotopeArray[collect(CartesianIndex(c) for c in [gFrontierN; cFrontierN])]#zonotopeArray[collect(c.id for c in frontierN)]
    #push!(frontierZonotopes, zonotopeArray[40, 40, 1])
    #reachtimes = map(z -> propagate_set(z, [0.0, p], δ⁻, δ⁺, system, phiDict, tPhiDict, alg, maxOrder, reduceOrder), frontierZonotopes)

    #
    #   0: Didnt hit either constraint, guard or invariant
    #   1: Hit constraint
    #   2: Hit guard
    #   3: Hit invariant surface
    #
    reachable_by_action::Dict{Tuple{Int64,CartesianIndex{grid.dimension}},AbstractArray{CartesianIndex{grid.dimension}}} = Dict{Tuple{Int64,CartesianIndex{grid.dimension}},AbstractArray{CartesianIndex{grid.dimension}}}()
    reachable_by_flow::Dict{CartesianIndex{grid.dimension},AbstractArray{CartesianIndex{grid.dimension}}} = Dict{CartesianIndex{grid.dimension},AbstractArray{CartesianIndex{grid.dimension}}}()
    degenerate_dimensions = collect(u - l == l ? l - (granularity / 2) : 0. for (l, u) in zip(grid.lower, grid.upper))
    invariant = hyperrectangle_to_HPolyhedron(system.statespace)
    U = system.input
    max_input = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, p, alg, false, nothing), U)
    max_flow = exp(p .* system.flowMatrix)
    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(system.globalConstraints[1])
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(system.edges[1].guard)
    invariantProjVectors, invariantProjBounds = getHalfSpaceProjections(invariant)

    hbox = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, grid.dimension))
    guard_invariant_intersection = intersection(system.edges[1].guard, invariant)

    attempt_recorder = zeros(4)
    identity_matrix = diagm(ones(grid.dimension))

    fast_cells = falses(grid.numCells...)
    Z_big = remove_zero_generators(Zonotope(degenerate_dimensions, diagm(map(x -> x == 1 ? 0. : granularity * (sweep_factor/2), grid.numCells))))
    lower_offset_big = grid.lower .+ (granularity*(sweep_factor/2))
    generatorDict_big = ReACT_discretize_decomposed_generators(Z_big, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict)
    discretizationDict_big, finaldims_big = ReACT_discretize_combine_with_offsets(Z_big, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict_big)
    preallocated_center_big = zeros(grid.dimension)
    preallocated_genmat_big = zeros(grid.dimension, finaldims_big)
    preallocated_inter_genmat_big = zeros(grid.dimension, finaldims_big)
    ranges = [1:floor(Int64, grid.numCells[i]/sweep_factor) for i in 1:grid.dimension]
    #@show CartesianIndices((ranges...,))
    #solvr = LazySets.default_lp_solver(Float64)
    for bidx in CartesianIndices((ranges...,))

        of = lower_offset_big + ((car2vec(bidx) .- 1) .* (granularity*sweep_factor))
        LazySets.API.translate!(Z_big, of)



        ReACT_discretize_combine_with_offset_vector!(preallocated_center_big, preallocated_genmat_big, preallocated_inter_genmat_big, discretizationDict_big, of, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict_big)
        #z, f = propagate_set(Z, [0.0, p], δ⁻, δ⁺, system, P2A_abs, phiDict, tPhiDict, inputDict, max_input, max_flow, invariant, alg, maxOrder, reduceOrder)

        f = propagate_set_check_hybrid_jump(Z_big, [0.0, p], δ⁻, δ⁺, system, system.edges, identity_matrix, attempt_recorder, grid.dimension, P2A_abs, guardProjVectors, guardProjBounds, discretizationDict_big, phiDict, tPhiDict, inputDict, max_input, max_flow, guard_invariant_intersection, alg, maxOrder, reduceOrder)

        if f == 0
            #@show collect(((b-1)*sweep_factor):(b*sweep_factor) for b in car2vec(bidx))
            fast_cells[collect((1+(b-1)*sweep_factor):(b*sweep_factor) for b in car2vec(bidx))...] .= true
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

    max_zonotope = minkowski_sum(max_zonotope, Zonotope(zeros(grid.dimension), diagm(fill(granularity/2, grid.dimension))))
    mdl, c_ref = setup_model_constant_zonotope(max_zonotope)



    for global_constraint in system.globalConstraints
        #guard_invariant_intersection = intersection(guard_invariant_intersection, global_constraint)
        for didx in get_touching_cell_idxs_l(grid, global_constraint)
            #@show didx
            grid.deadCells[didx] = true
        end
    end
    for idx in CartesianIndices(grid.deadCells)
        if grid.deadCells[idx]
            #println("Skipping dead cell: ", idx)
            continue

            #else
            #@show "Processing cell: ", idx
        end
        if !fast_cells[idx]
            count += 1
            of = lower_offset + ((car2vec(idx) .- 1) .* granularity)
            LazySets.API.translate!(Z, of)
            #@show Z, idx, of


            ReACT_discretize_combine_with_offset_vector!(preallocated_center, preallocated_genmat, preallocated_inter_genmat, discretizationDict, of, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict)
            #z, f = propagate_set(Z, [0.0, p], δ⁻, δ⁺, system, P2A_abs, phiDict, tPhiDict, inputDict, max_input, max_flow, invariant, alg, maxOrder, reduceOrder)
            #=
            discdisc = newReACTDiscretizePlus(Z, δ⁻, δ⁺, system.flowMatrix, P2A_abs, phiDict, inputDict, alg, maxOrder, reduceOrder)
            for k in keys(discdisc)
                dumdum = Zonotope(discretizationDict[k].center, hcat(generatorDict[k], (Z.center - (phiDict[k]*Z.center))/2))
                @show discretizationDict[k].generators
                @show dumdum.generators
                @show (Z.center - (phiDict[k]*Z.center))/2, phiDict[k]
                @show low(discdisc[k]), low(discretizationDict[k]), low(dumdum)
                @show high(discdisc[k]), high(discretizationDict[k]), high(dumdum)
                @show discdisc[k].center, discretizationDict[k].center, size(genmat(discretizationDict[k])), size(generatorDict[k]), k
            end
            @show generatorDict[0]
            =#
            #throw(error("ASHJDAJSHDA"))
            z, f = propagate_set_t(Z, [0.0, p], δ⁻, δ⁺, system, system.edges, identity_matrix, attempt_recorder, grid.dimension, P2A_abs, guardProjVectors, guardProjBounds, discretizationDict, phiDict, tPhiDict, inputDict, max_input, max_flow, guard_invariant_intersection, alg, maxOrder, reduceOrder)

            if f != 1 #&& f != 3


                reachable_by_flow[idx] = get_touching_cell_idxs_lz(grid, z)


            else
                grid.deadCells[idx] = true
            end

            LazySets.API.translate!(Z, -of)
            if count % 1000 == 0
                println(count)
                #GC.gc()
            end
        else
            count += 1
            of = lower_offset + ((car2vec(idx) .- 1) .* granularity)
            of = max_flow * of
            LazySets.API.translate!(max_zonotope, of)

            # TODO - Might as well precompute the max_flow with input zonotope and just translate it back and fourth 
            #reachable_by_flow[idx] = get_touching_cell_idxs_lz(grid, minkowski_sum(max_input, linear_map(max_flow, Z)))
            #empty!(solvr)
            #reachable_by_flow[idx] = get_touching_cell_idxs_lz(grid, max_zonotope)
            reachable_by_flow[idx] = get_touching_cell_idxs_test(grid, max_zonotope, mdl, c_ref)

            LazySets.API.translate!(max_zonotope, -of)
            if count % 1000 == 0
                println(count)
                #GC.gc()
            end
        end
    end

    #can_act = false
    for act in system.Act
        #@show size(get_touching_cell_idxs_l(grid, act.guard))
        for idx in get_touching_cell_idxs_l(grid, act.guard)
            of = lower_offset + ((car2vec(idx) .- 1) .* granularity)
            LazySets.API.translate!(Z, of)
            #touches_list = get_touching_cell_idxs(grid, z)
            #touches_list = get_touching_cell_idxs(grid, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector))
            #if !isempty(touches_list)
            #    reachable_by_action[(act, idx)] = copy(touches_list)
            #reachable_by_flow[idx] = get_touching_cell_idxs(grid, z)
            #end
            reachable_by_action[(act.id, idx)] = get_touching_cell_idxs_lz(grid, minkowski_sum(linear_map(max_flow, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector)), max_input))

            #if act.id == 1
            #@show get_touching_cell_idxs(grid, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector)), LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector)
            #reachable_by_action[(act.id, idx)] = get_touching_cell_idxs(grid, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector))
            #end
            if grid.deadCells[idx]
                grid.deadCells[idx] = false
            end
            LazySets.API.translate!(Z, -of)
        end
    end

    #@show dc
    return grid, reachable_by_action, reachable_by_flow

end

function ReACTed_reachable_cell_online(system::EuclideanHybridSystem, p, granularity, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, clustering=true, max_steps::Int64=1000)
    grid = Grid(system.statespace, granularity)
    @show foldl(*, grid.numCells), grid.numCells
    #max_steps = 300
    shield_count = 0
    sweep_factor = 20

    phiDict, tPhiDict, inputDict = phi_input_dict(system.flowMatrix, system.input, δ⁻, δ⁺, alg, maxOrder, reduceOrder)
    d = δ⁻
    #dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    isInvA = isinvertible(system.flowMatrix)
    A = system.flowMatrix
    #Φ = copy(phiDict[d])
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(system.flowMatrix)
    Φcache = system.flowMatrix == A_abs ? phiDict[d] : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)

    #frontierZonotopes = zonotopeArray[collect(CartesianIndex(c) for c in [gFrontierN; cFrontierN])]#zonotopeArray[collect(c.id for c in frontierN)]
    #push!(frontierZonotopes, zonotopeArray[40, 40, 1])
    #reachtimes = map(z -> propagate_set(z, [0.0, p], δ⁻, δ⁺, system, phiDict, tPhiDict, alg, maxOrder, reduceOrder), frontierZonotopes)

    #
    #   0: Didnt hit either constraint, guard or invariant
    #   1: Hit constraint
    #   2: Hit guard
    #   3: Hit invariant surface
    #

    reachable_by_action::Dict{Tuple{Int64,CartesianIndex{grid.dimension}},AbstractArray{CartesianIndex{grid.dimension}}} = Dict{Tuple{Int64,CartesianIndex{grid.dimension}},AbstractArray{CartesianIndex{grid.dimension}}}()
    reachable_by_flow::Dict{CartesianIndex{grid.dimension},AbstractArray{CartesianIndex{grid.dimension}}} = Dict{CartesianIndex{grid.dimension},AbstractArray{CartesianIndex{grid.dimension}}}()
    degenerate_dimensions = collect(u - l == l ? l - (granularity / 2) : 0. for (l, u) in zip(grid.lower, grid.upper))

    invariant = hyperrectangle_to_HPolyhedron(system.statespace)
    U = system.input
    max_input = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, p, alg, false, nothing), U)
    max_flow = exp(p .* system.flowMatrix)
    #@show degenerate_dimensions

    guardProjVectors, guardProjBounds = getHalfSpaceProjections(system.edges[1].guard)
    guard_invariant_intersection = intersection(system.edges[1].guard, invariant)


    fast_cells = falses(grid.numCells...)
    Z_big = remove_zero_generators(Zonotope(degenerate_dimensions, diagm(map(x -> x == 1 ? 0. : granularity * (sweep_factor/2), grid.numCells))))
    lower_offset_big = grid.lower .+ (granularity*(sweep_factor/2))
    generatorDict_big = ReACT_discretize_decomposed_generators(Z_big, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict)
    discretizationDict_big, finaldims_big = ReACT_discretize_combine_with_offsets(Z_big, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict_big)
    preallocated_center_big = zeros(grid.dimension)
    preallocated_genmat_big = zeros(grid.dimension, finaldims_big)
    preallocated_inter_genmat_big = zeros(grid.dimension, finaldims_big)
    ranges = [1:floor(Int64, grid.numCells[i]/sweep_factor) for i in 1:grid.dimension]
    #@show CartesianIndices((ranges...,))
    #solvr = LazySets.default_lp_solver(Float64)

    attempt_recorder = zeros(4)
    identity_matrix = diagm(ones(grid.dimension))
    for bidx in CartesianIndices((ranges...,))

        of = lower_offset_big + ((car2vec(bidx) .- 1) .* (granularity*sweep_factor))
        LazySets.API.translate!(Z_big, of)



        ReACT_discretize_combine_with_offset_vector!(preallocated_center_big, preallocated_genmat_big, preallocated_inter_genmat_big, discretizationDict_big, of, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict_big)
        #z, f = propagate_set(Z, [0.0, p], δ⁻, δ⁺, system, P2A_abs, phiDict, tPhiDict, inputDict, max_input, max_flow, invariant, alg, maxOrder, reduceOrder)

        f = propagate_set_check_hybrid_jump(Z_big, [0.0, p], δ⁻, δ⁺, system, system.edges, identity_matrix, attempt_recorder, grid.dimension, P2A_abs, guardProjVectors, guardProjBounds, discretizationDict_big, phiDict, tPhiDict, inputDict, max_input, max_flow, guard_invariant_intersection, alg, maxOrder, reduceOrder)

        if f == 0
            #@show collect(((b-1)*sweep_factor):(b*sweep_factor) for b in car2vec(bidx))
            fast_cells[collect((1+(b-1)*sweep_factor):(b*sweep_factor) for b in car2vec(bidx))...] .= true
            #grid.deadCells[collect((1+(b-1)*sweep_factor):(b*sweep_factor) for b in car2vec(bidx))...] .= true

        end


        #grid.array[idx].sidx[1] = copy(f)
        LazySets.API.translate!(Z_big, -of)
    end

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

    max_zonotope::Zonotope{Float64,Vector{Float64},Matrix{Float64}} = Zonotope(max_input.center, genmat(linear_map(max_flow, Z)))

    max_zonotope = minkowski_sum(max_zonotope, Zonotope(zeros(grid.dimension), diagm(fill(granularity/2, grid.dimension))))
    mdl, c_ref = setup_model_constant_zonotope(max_zonotope)


    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(system.globalConstraints[1])
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(system.edges[1].guard)
    invariantProjVectors, invariantProjBounds = getHalfSpaceProjections(invariant)

    hbox = Hyperrectangle(grid.lower .+ (granularity/2), fill(granularity/2, grid.dimension))



    done_list = CartesianIndex{grid.dimension}[]
    frontier_list = CartesianIndex{grid.dimension}[]
    ln = 0
    la = 0
    ld = 0
    dims = grid.dimension
    #dead = grid.deadCells
    can_act_matrix = trues(grid.numCells...)
    no_action_bad_matrix = falses(grid.numCells...)

    offset_list = offsets(grid)

    may_act_dict::Dict{Int64,AbstractArray{CartesianIndex{grid.dimension}}} = Dict{Int64,AbstractArray{CartesianIndex{grid.dimension}}}()
    #may_act_dict::Dict{Int64,BitMatrix} = Dict{Int64,BitMatrix}()
    may_act_matrix = falses(grid.numCells...)
    for act in system.Act
        may_act_dict[act.id] = get_touching_cell_idxs_l(grid, act.guard) #falses(grid.numCells...)
        #may_act_dict[act.id][get_touching_cell_idxs_l(grid, act.guard)...] .= true
        may_act_matrix[may_act_dict[act.id]] .= true
    end

    for edge in system.edges
        for didx in get_touching_cell_idxs_l(grid, edge.guard)
            #@show didx
            push!(frontier_list, didx)
        end
    end

    for global_constraint in system.globalConstraints
        #guard_invariant_intersection = intersection(guard_invariant_intersection, global_constraint)
        for didx in get_touching_cell_idxs_l(grid, global_constraint)
            #@show didx
            grid.deadCells[didx] = true
            no_action_bad_matrix[didx] = true
            push!(done_list, didx)
        end
    end
    push!(frontier_list, grow_indices_t(grid, done_list, offset_list)...)
    #return grid, reachable_by_action, reachable_by_flow
    while shield_count < max_steps
        for idx in frontier_list
            if grid.deadCells[idx]
                #println("Skipping dead cell: ", idx)
                continue

                #else
                #@show "Processing cell: ", idx
            end
            if !fast_cells[idx]
                count += 1
                of = lower_offset + ((car2vec(idx) .- 1) .* granularity)
                LazySets.API.translate!(Z, of)
                #@show Z, idx, of


                ReACT_discretize_combine_with_offset_vector!(preallocated_center, preallocated_genmat, preallocated_inter_genmat, discretizationDict, of, δ⁻, δ⁺, A, P2A_abs, phiDict, U, inputDict, generatorDict)
                #z, f = propagate_set(Z, [0.0, p], δ⁻, δ⁺, system, P2A_abs, phiDict, tPhiDict, inputDict, max_input, max_flow, invariant, alg, maxOrder, reduceOrder)

                z, f = propagate_set_t(Z, [0.0, p], δ⁻, δ⁺, system, system.edges, identity_matrix, attempt_recorder, grid.dimension, P2A_abs, guardProjVectors, guardProjBounds, discretizationDict, phiDict, tPhiDict, inputDict, max_input, max_flow, guard_invariant_intersection, alg, maxOrder, reduceOrder)

                if f != 1 #&& f != 3


                    reachable_by_flow[idx] = get_touching_cell_idxs_lz(grid, z)


                else
                    grid.deadCells[idx] = true
                end

                LazySets.API.translate!(Z, -of)
                if count % 1000 == 0
                    println(count)
                    #GC.gc()
                end
            else
                count += 1
                of = lower_offset + ((car2vec(idx) .- 1) .* granularity)
                of = max_flow * of
                LazySets.API.translate!(max_zonotope, of)

                #reachable_by_flow[idx] = get_touching_cell_idxs_lz(grid, max_zonotope)
                reachable_by_flow[idx] = get_touching_cell_idxs_test(grid, max_zonotope, mdl, c_ref)

                LazySets.API.translate!(max_zonotope, -of)
                if count % 1000 == 0
                    println(count)
                    #GC.gc()
                end
            end

            #can_act = false
            if may_act_matrix[idx]
                for act in system.Act
                    #@show size(get_touching_cell_idxs_l(grid, act.guard))
                    if idx in may_act_dict[act.id]
                        of = lower_offset + ((car2vec(idx) .- 1) .* granularity)
                        LazySets.API.translate!(Z, of)
                        #touches_list = get_touching_cell_idxs(grid, z)
                        #touches_list = get_touching_cell_idxs(grid, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector))
                        #if !isempty(touches_list)
                        #    reachable_by_action[(act, idx)] = copy(touches_list)
                        #reachable_by_flow[idx] = get_touching_cell_idxs(grid, z)
                        #end
                        reachable_by_action[(act.id, idx)] = get_touching_cell_idxs_lz(grid, minkowski_sum(linear_map(max_flow, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector)), max_input))

                        #if act.id == 1
                        #@show get_touching_cell_idxs(grid, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector)), LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector)
                        #reachable_by_action[(act.id, idx)] = get_touching_cell_idxs(grid, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector))
                        #end
                        if grid.deadCells[idx]
                            grid.deadCells[idx] = false
                        end
                        LazySets.API.translate!(Z, -of)
                    end
                end
            end
            push!(done_list, idx)
        end
        # TODO - Maybe this should be frontier_list with done_list subset difference with deadCells
        _, la, ln, ld = shield_step!(grid.deadCells, done_list, can_act_matrix, no_action_bad_matrix, reachable_by_action, reachable_by_flow, system.Act, dims)
        if shield_count % 5 == 0
            println("$(max_steps - shield_count) : $(Dates.format(now(), "HH:MM:SS"))")
        end
        #@show can_act_matrix == trues(grid.numCells...)
        if la == 0 && ln == 0 && ld == 0

            println("Fixed point found at $(max_steps-shield_count) steps!")
            #grid.deadCells = dead´

            break
        end
        shield_count += 1
        frontier_list = []
        # TODO - This is missing the cells where no_act_bad_matrix is true
        #push!(frontier_list, grow_indices_t(grid, findall(grid.deadCells), offset_list)...)
        push!(frontier_list, grow_indices_t(grid, findall(no_action_bad_matrix), offset_list)...)
        filter!(t -> !in(t, done_list), frontier_list)
    end


    #@show frontier_list

    return grid, reachable_by_action, reachable_by_flow, no_action_bad_matrix

end

#
#   Right now this makes quite a bit of allocations.
#   TODO - Needs a fix for determining the final reachset. Right now it is overapproximated and this is worse for large timestep sizes.
#
function propagate_set(X0, interval, δ⁻::Float64, δ⁺::Float64, system, Φ₂, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, PhiDict, TPhiDict, inputDict, max_input, max_flow, invariant, invariantProjVectors, invariantProjBounds, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5)

    discretizationDict = newReACTDiscretizePlus(X0, δ⁻, δ⁺, system.flowMatrix, Φ₂, PhiDict, inputDict, alg, maxOrder, reduceOrder)
    time::Float64 = copy(minimum(interval))
    endtime::Float64 = maximum(interval)
    constraint = system.globalConstraints[1]
    listOfEdges = system.edges



    for edge in listOfEdges
        invariant_projection_vectors = copy(invariantProjVectors)
        invariant_projection_bounds = copy(invariantProjBounds)
        guard_projection_vectors = copy(guardProjVectors)
        guard_projection_bounds = copy(guardProjBounds)
        constraint_projection_vectors = copy(constraintProjVectors)
        constraint_projection_bounds = copy(constraintProjBounds)

        #Sρ = zeros(Float64, length(constraintProjVectors))
        #Gρ = zeros(Float64, length(guardProjVectors))
        #Iρ = zeros(Float64, length(invarientProjVectors))

        #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 
        reachtime, newSet, inputSet, flag = ReACT_guards_b(δ⁻, δ⁺, interval, system.statespace, edge.guard, guard_projection_vectors, guard_projection_bounds, constraint, constraint_projection_vectors, constraint_projection_bounds, invariant, invariant_projection_vectors, invariant_projection_bounds, 2, PhiDict, TPhiDict, discretizationDict, inputDict)


        #
        #   0: Didnt hit either constraint, guard or invariant
        #   1: Hit constraint
        #   2: Hit guard
        #   3: Hit invariant surface
        #
        #@show flag, time, reachtime, newSet
        if flag == 0
            #@show time
            #if time == 0.0
            #    return (minkowski_sum(linear_map(max_flow, X0), max_input), 0)
            #else
            A_abs = abs.(system.flowMatrix)
            max_input = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, endtime - time, alg, false, nothing), system.input)

            return (minkowski_sum(linear_map(exp((endtime - time) .* system.flowMatrix), X0), linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, endtime - time, alg, false, nothing), system.input)), 0)
            #end
        elseif flag == 1 # Hit constraint
            #@show LazySets.API.isdisjoint(minkowski_sum(newSet, inputSet), system.globalConstraints[1]), newSet.center, flag, reachtime

            return (nothing, 1)
        elseif flag == 2
            throw(error("Nope not happening"))
            #@show inputSet
            intersectingSetsList = ReACT_touches_set_constant_input(newSet, inputSet, edge.guard, δ⁻, [reachtime, endtime], PhiDict, inputDict)
            tempIntersect = foldl(ConvexHull, intersectingSetsList)
            intersectedSet = convert(Zonotope, box_approximation(tempIntersect))
            jumpSet = linear_map(edge.jumpMatrix, zonotopeStripIntersection(intersectedSet, edge.guard))
            if !isnothing(edge.jumpVector)
                LazySets.translate!(jumpSet, edge.jumpVector)
            end
            #res = propagate_set(jumpSet, [reachtime, endtime], δ⁻, δ⁺, system, Φ₂, PhiDict, TPhiDict, inputDict, alg, maxOrder, reduceOrder)
            #@show res
            return propagate_set(jumpSet, [reachtime, endtime], δ⁻, δ⁺, system, Φ₂, PhiDict, TPhiDict, inputDict, max_input, max_flow, invariant, alg, maxOrder, reduceOrder)
        end
        #@show LazySets.API.isdisjoint(minkowski_sum(newSet, inputSet), invariant), newSet.center, flag, reachtime
        return (minkowski_sum(newSet, inputSet), 3)
    end

end

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

        #return (minkowski_sum(linear_map(exp((endtime - time) .* system.flowMatrix), X0), max_input), 0)
        #
        #   0: Didnt hit either constraint, guard or invariant
        #   1: Hit constraint
        #   2: Hit guard
        #   3: Hit invariant surface
        #
        #@show flag, time, reachtime, newSet
        if flag == 0
            #@show time
            #if time == 0.0
            #    return (minkowski_sum(linear_map(max_flow, X0), max_input), 0)
            #else
            A_abs = abs.(system.flowMatrix)
            max_input = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, endtime - time, alg, false, nothing), system.input)

            return (minkowski_sum(linear_map(exp((endtime - time) .* system.flowMatrix), X0), max_input), 0)
            #end

        elseif flag == 2
            return (minkowski_sum(newSet, inputSet), 0)
            #return (minkowski_sum(linear_map(exp((endtime - time) .* system.flowMatrix), X0), max_input), 0)

            #@show reachtime, newSet, inputSet, time
            #throw(error("Nope not happening"))
            #@show inputSet
            intersectingSetsList = ReACT_touches_set_constant_input(newSet, inputSet, edge.guard, δ⁻, [reachtime, endtime], PhiDict, inputDict)
            tempIntersect = foldl(ConvexHull, intersectingSetsList)
            intersectedSet = convert(Zonotope, box_approximation(tempIntersect))
            jumpSet = linear_map(edge.jumpMatrix, zonotopeStripIntersection(intersectedSet, edge.guard))
            if !isnothing(edge.jumpVector)
                LazySets.translate!(jumpSet, edge.jumpVector)
            end
            #res = propagate_set(jumpSet, [reachtime, endtime], δ⁻, δ⁺, system, Φ₂, PhiDict, TPhiDict, inputDict, alg, maxOrder, reduceOrder)
            #@show res
            return propagate_jumpset(jumpSet, [reachtime, endtime], δ⁻, δ⁺, system, identity_matrix, attemptsRecorder, dims, Φ₂, guardProjVectors, guardProjBounds, PhiDict, TPhiDict, inputDict, alg, maxOrder, reduceOrder)
        end
        #@show LazySets.API.isdisjoint(minkowski_sum(newSet, inputSet), invariant), newSet.center, flag, reachtime
        return (minkowski_sum(newSet, inputSet), 1)
    end

end

function propagate_set_b(X0, interval, δ⁻::Float64, δ⁺::Float64, system, Φ₂, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, discretizationDict, PhiDict, TPhiDict, inputDict, max_input, max_flow, guard_intersection, invariant, invariantProjVectors, invariantProjBounds, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5)

    #discretizationDict = newReACTDiscretizePlus(X0, δ⁻, δ⁺, system.flowMatrix, Φ₂, PhiDict, inputDict, alg, maxOrder, reduceOrder)
    time::Float64 = copy(minimum(interval))
    endtime::Float64 = maximum(interval)
    constraint = system.globalConstraints[1]
    listOfEdges = system.edges

    for edge in listOfEdges
        #constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
        #guardProjVectors, guardProjBounds = getHalfSpaceProjections(edge.guard)
        invariant_projection_vectors = copy(invariantProjVectors)
        invariant_projection_bounds = copy(invariantProjBounds)
        guard_projection_vectors = copy(guardProjVectors)
        guard_projection_bounds = copy(guardProjBounds)
        constraint_projection_vectors = copy(constraintProjVectors)
        constraint_projection_bounds = copy(constraintProjBounds)


        #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 
        reachtime, newSet, inputSet, flag = ReACT_guards_b(δ⁻, δ⁺, interval, system.statespace, edge.guard, guard_projection_vectors, guard_projection_bounds, constraint, constraint_projection_vectors, constraint_projection_bounds, invariant, invariant_projection_vectors, invariant_projection_bounds, 2, PhiDict, TPhiDict, discretizationDict, inputDict)


        #
        #   0: Didnt hit either constraint, guard or invariant
        #   1: Hit constraint
        #   2: Hit guard
        #   3: Hit invariant surface
        #
        if flag == 0
            if time == 0.0
                return (minkowski_sum(linear_map(max_flow, X0), max_input), 0)
            else
                A_abs = abs.(system.flowMatrix)
                max_input = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, endtime - time, alg, false, nothing), system.input)

                return (minkowski_sum(linear_map(exp((endtime - time) .* system.flowMatrix), X0), linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, endtime - time, alg, false, nothing), system.input)), 0)
            end
        elseif flag == 1 # Hit constraint
            return (minkowski_sum(newSet, inputSet), 1)
        elseif flag == 2
            # Move this outside the function
            #return (minkowski_sum(newSet, inputSet), 2)

            intersectingSetsList = ReACT_touches_set_constant_input_b(newSet, inputSet, edge.guard, copy(guardProjVectors), copy(guardProjBounds), invariantProjVectors, invariantProjBounds, δ⁻, [reachtime, endtime], PhiDict, inputDict)



            tempIntersect = foldl(ConvexHull, intersectingSetsList)
            intersectedSet = convert(Zonotope, box_approximation(tempIntersect))
            #jumpSet = linear_map(edge.jumpMatrix, zonotopeStripIntersection(minkowski_sum(newSet, inputSet), guard_intersection))
            #if !LazySets.API.issubset(minkowski_sum(newSet, inputSet), invariant)
            #    @show minkowski_sum(newSet, inputSet), reachtime
            #    return (minkowski_sum(newSet, inputSet), 3)
            #end
            jumpSet = linear_map(edge.jumpMatrix, zonotopeStripIntersection(intersectedSet, guard_intersection))
            #jumpSet = linear_map(edge.jumpMatrix, getBoxIntersection(intersectedSet, guard_intersection))
            if !isnothing(edge.jumpVector)
                LazySets.translate!(jumpSet, edge.jumpVector)
            end
            #=
            if LazySets.API.isdisjoint(jumpSet, invariant)
                @show jumpSet
                @show newSet, inputSet
                @show reachtime
            end
            =#
            return propagate_set(jumpSet, [reachtime, endtime], δ⁻, δ⁺, system, Φ₂, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, PhiDict, TPhiDict, inputDict, max_input, max_flow, invariant, invariantProjVectors, invariantProjBounds, alg, maxOrder, reduceOrder)
            # Cant use propagate_set_b here because we need to recompute the discretizationDict for the new jumpSet. This is because the jumpSet may be a different shape than the original set and thus the discretization may be different. We could try to reuse the discretizationDict but this would require a lot of work and is not worth it for now.
            #return propagate_set_b(jumpSet, [reachtime, endtime], δ⁻, δ⁺, system, Φ₂, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, discretizationDict, PhiDict, TPhiDict, inputDict, max_input, max_flow, invariant, alg, maxOrder, reduceOrder)

        end
        #println("shiit")
        return (newSet, 3)
    end
end

function propagate_set_t(X0, interval, δ⁻::Float64, δ⁺::Float64, system, listOfEdges, identity_matrix, attemptsRecorder, dims, Φ₂, guardProjVectors, guardProjBounds, discretizationDict, PhiDict, TPhiDict, inputDict, max_input, max_flow, guard_intersection, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5)

    #discretizationDict = newReACTDiscretizePlus(X0, δ⁻, δ⁺, system.flowMatrix, Φ₂, PhiDict, inputDict, alg, maxOrder, reduceOrder)
    time::Float64 = copy(minimum(interval))
    endtime::Float64 = maximum(interval)
    zero_input = Zonotope(zeros(dims), zeros(size(genmat((inputDict[δ⁻])))))
    for edge in listOfEdges
        guard_projection_vectors = copy(guardProjVectors)
        guard_projection_bounds = copy(guardProjBounds)



        #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 
        reachtime, newSet, inputSet, flag = ReACT_guards_t(δ⁻, δ⁺, interval, identity_matrix, attemptsRecorder, zero_input, guard_projection_vectors, guard_projection_bounds, 2, PhiDict, TPhiDict, discretizationDict, inputDict)


        #
        #   0: Didnt hit either constraint, guard or invariant
        #   1: Hit constraint
        #   2: Hit guard
        #   3: Hit invariant surface
        #
        if flag == 0
            return (minkowski_sum(linear_map(max_flow, X0), max_input), 0)
            if time == 0.0
            else
                A_abs = abs.(system.flowMatrix)
                max_input = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, endtime - time, alg, false, nothing), inputDict[0])

                return (minkowski_sum(linear_map(exp((endtime - time) .* system.flowMatrix), X0), linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A_abs, endtime - time, alg, false, nothing), system.input)), 0)
            end
        elseif flag == 2
            # Move this outside the function
            if LazySets.API.isdisjoint(minkowski_sum(newSet, inputSet), system.statespace) && LazySets.API.isdisjoint(minkowski_sum(newSet, inputSet), guard_intersection)
                #if !LazySets.API.isdisjoint(minkowski_sum(newSet, inputSet), edge.guard) && LazySets.API.isdisjoint(minkowski_sum(newSet, inputSet), system.statespace)
                #testSet = minkowski_sum(newSet, inputSet)
                #@show reachtime, X0, low(testSet), high(testSet)
                #@show LazySets.API.isdisjoint(testSet, edge.guard)
                return (nothing, 1)
            end
            intersectingSetsList = ReACT_touches_set_constant_input(newSet, inputSet, guard_intersection, δ⁻, [reachtime, endtime], PhiDict, inputDict)



            tempIntersect = foldl(ConvexHull, intersectingSetsList)
            intersectedSet = convert(Zonotope, box_approximation(tempIntersect))
            #jumpSet = linear_map(edge.jumpMatrix, zonotopeStripIntersection(minkowski_sum(newSet, inputSet), guard_intersection))
            #if !LazySets.API.issubset(minkowski_sum(newSet, inputSet), invariant)
            #    @show minkowski_sum(newSet, inputSet), reachtime
            #    return (minkowski_sum(newSet, inputSet), 3)
            #end

            jumpSet = linear_map(edge.jumpMatrix, zonotopeStripIntersection(intersectedSet, guard_intersection))
            #jumpSet = linear_map(edge.jumpMatrix, getBoxIntersection(intersectedSet, guard_intersection))
            if !isnothing(edge.jumpVector)
                LazySets.translate!(jumpSet, edge.jumpVector)
            end
            #=
            if LazySets.API.isdisjoint(jumpSet, invariant)
                @show jumpSet
                @show newSet, inputSet
                @show reachtime
            end
            =#
            return propagate_jumpset(jumpSet, [reachtime, endtime], δ⁻, δ⁺, system, identity_matrix, attemptsRecorder, dims, Φ₂, guardProjVectors, guardProjBounds, PhiDict, TPhiDict, inputDict, alg, maxOrder, reduceOrder)
            # Cant use propagate_set_b here because we need to recompute the discretizationDict for the new jumpSet. This is because the jumpSet may be a different shape than the original set and thus the discretization may be different. We could try to reuse the discretizationDict but this would require a lot of work and is not worth it for now.
            #return propagate_set_b(jumpSet, [reachtime, endtime], δ⁻, δ⁺, system, Φ₂, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, discretizationDict, PhiDict, TPhiDict, inputDict, max_input, max_flow, invariant, alg, maxOrder, reduceOrder)

        end
    end
end

function propagate_set_check_hybrid_jump(X0, interval, δ⁻::Float64, δ⁺::Float64, system, listOfEdges, identity_matrix, attemptsRecorder, dims, Φ₂, guardProjVectors, guardProjBounds, discretizationDict, PhiDict, TPhiDict, inputDict, max_input, max_flow, guard_intersection, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5)
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
    end
    return 0
end
#
#   Given we have must semantics for guards we could probably
#
function ReACT_time_touches_set(X0, set, δ⁻::Float64, δ⁺::Float64, interval, PhiDict, Φ)
    # Note that in touches we always use δ⁻
    # That is, we do not adjust timestep sizes

    #phiDict = PhiDict
    #=
    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(guard)
    invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)
    =#
    time::Float64 = interval[1]
    endtime::Float64 = interval[2]


    #Vs = nestedInputDiscCalculate(inputDiscritezationDict, PhiDict, δ⁻, initialTime, reduce_order, max_order)
    # i = 1

    #if ismissing(Φ)
    #    Φ::Matrix{Float64} = exp(initialTime .* loc.A)
    #end

    #ϕ = similar(Φ)
    #ϕ = phiDict[δ⁻]

    # Compute current sets
    #newR = discritezationDict[δ⁻]
    #newRR = linear_map(Φ, newR)
    #V = copy(inputDiscritezationDict[δ⁻])
    #V = linear_map(Φ, V)
    #intersectingSetsList = Vector{Zonotope}() # This will store all of our sets

    tempSet = linear_map(Φ, X0)
    while time < endtime

        #if all((ρ(x, tempSet)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) && # IsSubSet
        #   all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds)) &&
        #   all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Intersects
        if !LazySets.API.isdisjoint(set, tempSet) #touchesCheck(tempSet, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, invarientProjVectors, invarientProjBounds)
            #if mapreduce(x -> intersects(newRR, x), &, guard)
            #@show time
            #push!(intersectingSetsList, copy(tempSet))

            # Main calculation. No longer changing timesteps
            #Vs = minkowski_sum(Vs, V) #Update input
            #newRR = linear_map(ϕ, newRR)
            #V = linear_map(ϕ, V)
            tempSet = linear_map(PhiDict[δ⁻], tempSet)
            # i = i + 1
            time += δ⁻
        else
            break
        end
        #=
        else # Handle hit something. We do not reduce anymore!
            #@show all((ρ(x, tempSet)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds))
            #@show all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds))
            #@show all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds))
            #if all((ρ(x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) && any(((ρ(x, tempSet)) > y) for (x, y) in zip(constraintProjVectors, constraintProjBounds))
            #    handleHitConstraint(time, loc.id)
            #end

            #continueAfter = true
            #=
            if !((all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds)) || !all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds))) == (!someunder(tempSet, guardProjVectors, guardProjBounds) || !someunder(tempSet, invarientProjVectors, invarientProjBounds)))
                println((all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds)) || !all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds))))
                @show someunder(tempSet, guardProjVectors, guardProjBounds)
                @show all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds))
                @show !someunder(tempSet, invarientProjVectors, invarientProjBounds)
            end
            =#
            if someunder(tempSet, guardProjVectors, guardProjBounds) || !someunder(tempSet, invarientProjVectors, invarientProjBounds)
                # If we stop because we are no longer intersect guards, 
                # but still intersect the invariant we try continue

                #@show all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds))
                #continueAfter = false
                return false, time, intersectingSetsList
            end

            if allunder(tempSet, constraintProjVectors, constraintProjBounds)#all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Intersects
                return true, time, intersectingSetsList
            else
                intersectionSet = LazySets.Intersection(tempSet, loc.invarient)
                if !all((ρ(x, intersectionSet)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds))  # IsSubSet
                    @show LazySets.API.high(tempSet)
                    @show LazySets.API.low(tempSet)
                    @show LazySets.API.high(Vs)
                    @show LazySets.API.low(Vs)
                    handleHitConstraint(time, loc.id)
                end
                return true, time, intersectingSetsList
            end

            #=if all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) && continueAfter

                intersectSet = zonotopeStripIntersection(tempSet, loc.invarient)
                @show isSubSet(tempSet, loc.invarient)
                tem = intersection(overapproximate(tempSet, BoxDirections(LazySets.dim(newR))), loc.invarient)
                #@show map((x, y) -> ρ(x, tem) <= y, zip(constraintProjVectors, constraintProjBounds))
                @show all((ρ(x, tempSet)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) # IsSubSet
                @show all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds))
                @show all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Intersects
                @show minimum(interval) - time
                @show LazySets.API.high(tempSet)
                @show LazySets.API.low(tempSet)
                !(all((ρ(x, tem)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds))) && handleHitConstraint(time, loc.id)
            end=#

            #return continueAfter, time, intersectingSetsList
        end
        =#
    end
    #continueAfter = false # We are at the end time horizon, therefore no continuing
    return (LazySets.API.issubset(tempSet, set), time)
end

# TODO - Optimize this! We know that the set will be intersecting initially, thus we can just use the fast disjoint check and actually use larger step sizes as long as we dont hit the constraint!!!
function ReACT_touches_set_constant_input(X0, U, set, δ⁻::Float64, interval, PhiDict, inputDiscretizationDict)
    # Note that in touches we always use δ⁻
    # That is, we do not adjust timestep sizes

    time::Float64 = copy(interval[1])
    endtime::Float64 = interval[2]

    intersectingSetsList = Vector{Zonotope}() # This will store all of our sets


    #X0 = linear_map(PhiDict[δ⁻], X0)
    #U = linear_map(PhiDict[δ⁻], U) + inputDiscretizationDict[δ⁻]
    tempSet = minkowski_sum(X0, U)

    push!(intersectingSetsList, copy(tempSet))
    while time < endtime

        #if all((ρ(x, tempSet)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) && # IsSubSet
        #   all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds)) &&
        #   all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Intersects
        if !LazySets.API.isdisjoint(set, tempSet) #touchesCheck(tempSet, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, invarientProjVectors, invarientProjBounds)
            #if mapreduce(x -> intersects(newRR, x), &, guard)
            #@show time
            push!(intersectingSetsList, copy(tempSet))

            # Main calculation. No longer changing timesteps
            #Vs = minkowski_sum(Vs, V) #Update input
            #newRR = linear_map(ϕ, newRR)
            #V = linear_map(ϕ, V)
            U = linear_map(PhiDict[δ⁻], U) + inputDiscretizationDict[δ⁻]
            X0 = linear_map(PhiDict[δ⁻], X0)
            tempSet = minkowski_sum(X0, U)
            # i = i + 1
            time += δ⁻
        else
            break
        end

    end
    #continueAfter = false # We are at the end time horizon, therefore no continuing
    return intersectingSetsList
end

function ReACT_touches_set_constant_input_b(X0, U, set, projection_vectors, projection_bounds, invariantProjVectors, invariantProjBounds, δ⁻::Float64, interval, PhiDict, inputDiscretizationDict)
    # Note that in touches we always use δ⁻
    # That is, we do not adjust timestep sizes

    time::Float64 = copy(interval[1])
    endtime::Float64 = interval[2]

    intersectingSetsList = Vector{Zonotope}() # This will store all of our sets


    #X0 = linear_map(PhiDict[δ⁻], X0)
    #U = linear_map(PhiDict[δ⁻], U) + inputDiscretizationDict[δ⁻]
    tempSet = minkowski_sum(X0, U)

    push!(intersectingSetsList, copy(tempSet))
    while time < endtime

        #if all((ρ(x, tempSet)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) && # IsSubSet
        #   all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds)) &&
        #   all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Intersects
        if !allunder(tempSet, projection_vectors, projection_bounds) && someinside(tempSet, projection_vectors, projection_bounds) && someinside(tempSet, invariantProjVectors, invariantProjBounds) #touchesCheck(tempSet, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, invarientProjVectors, invarientProjBounds)
            #@show time
            push!(intersectingSetsList, copy(tempSet))

            # Main calculation. No longer changing timesteps
            #Vs = minkowski_sum(Vs, V) #Update input
            #newRR = linear_map(ϕ, newRR)
            #V = linear_map(ϕ, V)
            U = linear_map(PhiDict[δ⁻], U) + inputDiscretizationDict[δ⁻]
            X0 = linear_map(PhiDict[δ⁻], X0)
            tempSet = minkowski_sum(X0, U)
            # i = i + 1
            time += δ⁻
        else
            break
        end

    end
    #continueAfter = false # We are at the end time horizon, therefore no continuing
    return intersectingSetsList
end

function ReACT_touches_set(X0, set, δ⁻::Float64, δ⁺::Float64, interval, PhiDict, inputDiscretizationDict)
    # Note that in touches we always use δ⁻
    # That is, we do not adjust timestep sizes

    time::Float64 = interval[1]
    endtime::Float64 = interval[2]

    intersectingSetsList = Vector{Zonotope}() # This will store all of our sets

    tempSet = X0
    while time < endtime

        #if all((ρ(x, tempSet)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) && # IsSubSet
        #   all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds)) &&
        #   all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Intersects
        if !LazySets.API.isdisjoint(set, tempSet) #touchesCheck(tempSet, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, invarientProjVectors, invarientProjBounds)
            #if mapreduce(x -> intersects(newRR, x), &, guard)
            #@show time
            push!(intersectingSetsList, copy(tempSet))

            # Main calculation. No longer changing timesteps
            #Vs = minkowski_sum(Vs, V) #Update input
            #newRR = linear_map(ϕ, newRR)
            #V = linear_map(ϕ, V)
            tempSet = minkowski_sum(linear_map(PhiDict[δ⁻], tempSet), inputDiscretizationDict[δ⁻])
            # i = i + 1
            time += δ⁻
        else
            break
        end

    end
    #continueAfter = false # We are at the end time horizon, therefore no continuing
    return (LazySets.API.issubset(tempSet, set), intersectingSetsList, time)
end


# TODO - For some reason changing this to only use the more precise method when the smallest time step is used gives a different answer than the previous version...
function ReACT_guards_b(δ⁻::Float64, δ⁺::Float64, interval, statespace, guards, guardProjVectors, guardProjBounds, constraint, constraintProjVectors, constraintProjBounds, invariant, invariantProjVectors, invariantProjBounds, STRATEGY::Integer, PhiDict, permutedphiDict, discritezationDict, inputDiscretizationDict)

    #
    #   These functions should be evaluated outside ReACT_guards and passed as parameters
    #
    #constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    #guardProjVectors, guardProjBounds = getHalfSpaceProjections(guards)
    #invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(invariant)

    dims = LazySets.API.dim(statespace)

    time::Float64 = copy(minimum(interval))
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(δ⁺)

    attemptsRecorder = zeros(4)

    Φ::Matrix{Float64} = diagm(ones(Float64, dims))

    tempM = diagm(ones(Float64, dims))

    #   Because we assume an empty input set i have trouble seeing how these can become relevant.
    #   But i am keeping them for now to reuse as much as possible

    #Sρ = zeros(Float64, length(constraintProjVectors))
    #Gρ = zeros(Float64, length(guardProjVectors))
    #Iρ = zeros(Float64, length(invarientProjVectors))


    i = 1
    input = Zonotope(zeros(dims), zeros(size(genmat((inputDiscretizationDict[δ⁻])))))
    while time < endtime

        attempts = 1
        approveFlag = false
        check = 0
        while !approveFlag
            # Handle if we can no longer reduce the reachset (we keep hitting something)
            if currentTimeStep > δ⁻
                #linear_map!(Z, Φ, discritezationDict[currentTimeStep])
                #check = guardCheck(minkowski_sum(input, linear_map(Φ, discritezationDict[currentTimeStep])), constraint, guards, invariant)
                if fast_guards_check_b(discritezationDict[currentTimeStep], input, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, invariantProjVectors, invariantProjBounds)# || (0 == guardCheck(input, Φ, discritezationDict[currentTimeStep], constraint, guards, invariant)) # check if fast or slow fails.
                    #=
                    if (0 != guardCheck(input, Φ, discritezationDict[currentTimeStep], constraint, guards, invariant))#0 == guardCheck(discritezationDict[currentTimeStep], Sρ, constraintProjVectors, constraintProjBounds, Gρ, guardProjVectors, guardProjBounds, Iρ, invarientProjVectors, invarientProjBounds)
                        println("Fast check passed but slow check failed. This should not happen. Check the fast check implementation.")
                        @show guardCheck(input, Φ, discritezationDict[currentTimeStep], constraint, guards, invariant), time, someinside_transparent(discritezationDict[currentTimeStep], input, invariantProjVectors, invariantProjBounds)
                        @show someinside_transparent(discritezationDict[currentTimeStep], input, map!(x -> permutedphiDict[δ⁻] * x, invariantProjVectors), invariantProjBounds)
                    end
                    =#
                    approveFlag = true
                    if currentTimeStep + time < endtime
                        #input = input + linear_map(Φ, inputDiscretizationDict[currentTimeStep])
                        #input = linear_map(Φ, input) + inputDiscretizationDict[currentTimeStep]
                        input = linear_map(PhiDict[currentTimeStep], input) + inputDiscretizationDict[currentTimeStep]

                        #Sρ += map(x -> ρ(x, inputDiscretizationDict[currentTimeStep]), constraintProjVectors)
                        #Gρ += map(x -> ρ(x, inputDiscretizationDict[currentTimeStep]), guardProjVectors)
                        #Iρ += map(x -> ρ(x, inputDiscretizationDict[currentTimeStep]), invarientProjVectors)


                        map!(x -> permutedphiDict[currentTimeStep] * x, constraintProjVectors)
                        map!(x -> permutedphiDict[currentTimeStep] * x, guardProjVectors)
                        map!(x -> permutedphiDict[currentTimeStep] * x, invariantProjVectors)

                        mul!(tempM, Φ, PhiDict[currentTimeStep])
                        copy!(Φ, tempM)

                    else
                        # TODO - The input gets used or referenced somewhere even with the zero flag..
                        #=
                        d = max(2^(floor(log2((endtime - time)/δ⁻))) * δ⁻, δ⁻)
                        #@show (endtime - time, d)
                        input = linear_map(PhiDict[δ⁻], input) + inputDiscretizationDict[δ⁻]
                        input = linear_map(PhiDict[d], input) + inputDiscretizationDict[d]
                        mul!(tempM, Φ, PhiDict[d])
                        copy!(Φ, tempM)
                        mul!(tempM, Φ, PhiDict[δ⁻])
                        copy!(Φ, tempM)
                        time = time + d #+ δ⁻
                        #println("from coarse")
                        =#
                        #return (time, linear_map(Φ, discritezationDict[δ⁻]), input, guardCheck(input, Φ, discritezationDict[δ⁻], constraint, guards, invariant))
                        return (endtime, discritezationDict[δ⁻], input, 0)
                    end
                else

                    currentTimeStep /= 2
                    attempts += 1
                end


            elseif currentTimeStep == δ⁻
                #elseif currentTimeStep == (δ⁻/2)
                #    currentTimeStep = δ⁻
                #linear_map!(Z, Φ, discritezationDict[currentTimeStep])
                check = guardCheck(input, Φ, discritezationDict[currentTimeStep], constraint, guards, invariant)
                #if check == 0 && !allunder(discritezationDict[δ⁻], input, invariantProjVectors, invariantProjBounds)
                #println
                #return (time, linear_map(Φ, discritezationDict[δ⁻]), input, 3)
                #end
                #if fast_guards_check_b(discritezationDict[currentTimeStep], input, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds) || (0 == guardCheck(input, Φ, discritezationDict[currentTimeStep], constraint, guards, invariant)) # check if fast or slow fails.
                if (0 == check) # check if fast or slow fails.
                    #if 0 == guardCheck(discritezationDict[currentTimeStep], Sρ, constraintProjVectors, constraintProjBounds, Gρ, guardProjVectors, guardProjBounds, Iρ, invarientProjVectors, invarientProjBounds)
                    approveFlag = true

                    #input = linear_map(Φ, input) + inputDiscretizationDict[currentTimeStep]
                    input = linear_map(PhiDict[currentTimeStep], input) + inputDiscretizationDict[currentTimeStep]

                    #Sρ += map(x -> ρ(x, inputDiscretizationDict[currentTimeStep]), constraintProjVectors)
                    #Gρ += map(x -> ρ(x, inputDiscretizationDict[currentTimeStep]), guardProjVectors)
                    #Iρ += map(x -> ρ(x, inputDiscretizationDict[currentTimeStep]), invarientProjVectors)


                    map!(x -> permutedphiDict[currentTimeStep] * x, constraintProjVectors)
                    map!(x -> permutedphiDict[currentTimeStep] * x, guardProjVectors)
                    map!(x -> permutedphiDict[currentTimeStep] * x, invariantProjVectors)

                    mul!(tempM, Φ, PhiDict[currentTimeStep])
                    copy!(Φ, tempM)


                else
                    #currentTimeStep /= 2
                    #attempts += 1

                    return (time, linear_map(Φ, discritezationDict[δ⁻]), input, check)

                end

                #else
                #println("hit something")
                #eturn (time, linear_map(Φ, discritezationDict[δ⁻]), input, check)
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
    #println("ever")
    return (endtime, linear_map(Φ, discritezationDict[δ⁻]), input, 0)
end

function ReACT_guards_t(δ⁻::Float64, δ⁺::Float64, interval, identity_matrix, attemptsRecorder, zero_input, guardProjVectors, guardProjBounds, STRATEGY::Integer, PhiDict, permutedphiDict, discritezationDict, inputDiscretizationDict)




    time::Float64 = copy(minimum(interval))
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(δ⁺)

    attemptsRecorder = attemptsRecorder .* 0

    Φ::Matrix{Float64} = copy(identity_matrix)

    tempM = copy(identity_matrix)
    lastM = copy(identity_matrix)
    i = 1
    input = copy(zero_input)
    last_input = copy(input)
    while time < endtime

        attempts = 1
        approveFlag = false
        check = 0
        while !approveFlag
            # Handle if we can no longer reduce the reachset (we keep hitting something)
            if currentTimeStep >= δ⁻

                if fast_guards_check_t(discritezationDict[currentTimeStep], input, guardProjVectors, guardProjBounds)
                    approveFlag = true
                    if currentTimeStep + time < endtime
                        #input = input + linear_map(Φ, inputDiscretizationDict[currentTimeStep])
                        #input = linear_map(Φ, input) + inputDiscretizationDict[currentTimeStep]
                        last_input = copy(input)
                        input = linear_map(PhiDict[currentTimeStep], input) + inputDiscretizationDict[currentTimeStep]

                        map!(x -> permutedphiDict[currentTimeStep] * x, guardProjVectors)

                        mul!(tempM, Φ, PhiDict[currentTimeStep])
                        copy!(lastM, Φ)
                        copy!(Φ, tempM)

                    else
                        # TODO - The input gets used or referenced somewhere even with the zero flag..
                        #=
                        d = max(2^(floor(log2((endtime - time)/δ⁻))) * δ⁻, δ⁻)
                        #@show (endtime - time, d)
                        input = linear_map(PhiDict[δ⁻], input) + inputDiscretizationDict[δ⁻]
                        input = linear_map(PhiDict[d], input) + inputDiscretizationDict[d]
                        mul!(tempM, Φ, PhiDict[d])
                        copy!(Φ, tempM)
                        mul!(tempM, Φ, PhiDict[δ⁻])
                        copy!(Φ, tempM)
                        time = time + d #+ δ⁻
                        #println("from coarse")
                        =#
                        #return (time, linear_map(Φ, discritezationDict[δ⁻]), input, guardCheck(input, Φ, discritezationDict[δ⁻], constraint, guards, invariant))
                        return (endtime, discritezationDict[δ⁻], input, 0)
                    end
                else

                    currentTimeStep /= 2
                    attempts += 1
                end

            else
                return (time, linear_map(Φ, discritezationDict[δ⁻]), input, 2)
                #return (time, linear_map(lastM, discritezationDict[δ⁻]), input, 2)

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

#=
function guardCheck(newR::Zonotope, Sρ::Vector{Float64}, constraintProjVectors::Vector{SparseArrays.SparseVector{Float64,Int}}, constraintProjBounds::Vector{Float64}, Gρ::Vector{Float64}, guardProjVectors::Vector{SparseArrays.SparseVector{Float64,Int}}, guardProjBounds::Vector{Float64}, Iρ::Vector{Float64}, invarientProjVectors::Vector{SparseArrays.SparseVector{Float64,Int}}, invarientProjBounds::Vector{Float64}) #; solver=model
    #cache = map((input, x, y) -> input + ρ(x, newR) <= y, Sρ, constraintProjVectors, constraintProjBounds)
    abssum = 0.0
    c = newR.center
    G = transpose(genmat(newR))
    tc = Vector{Float64}(undef, size(G, 1))
    #a = sum(abs, transpose(a) * G)

    res = all(input + tsupfunc(x, tc, abssum, c, G) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))

    res = res && any((input - tsupfunc(-x, tc, abssum, c, G)) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds))

    res = res && all((input - tsupfunc(-x, tc, abssum, c, G)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))

    return res
end
=#
#=
function guardCheck(newR::Zonotope, Sρ, constraintProjVectors, constraintProjBounds, Gρ, guardProjVectors, guardProjBounds, Iρ, invarientProjVectors, invarientProjBounds) #; solver=model



    if !someunder(newR, Sρ, constraintProjVectors, constraintProjBounds) && someunder(newR, Iρ, invarientProjVectors, invarientProjBounds)
        #println("didnt hit constraints")
        if someoutside(newR, Gρ, guardProjVectors, guardProjBounds)
            #println("didnt hit guard")
            return 0
        else
            return LazySets.API.isdisjoint(newR, HPolyhedron(collect(LazySets.HalfSpace(a, b+c) for (a, b, c) in zip(guardProjVectors, guardProjBounds, Gρ)))) ? 0 : 1
        end
    else
        #println("hit constraints")
        return 2
    end
end
=#

function guardCheck(input::Zonotope, Φ, X::Zonotope, constraint, guard, invariant) #; solver=model
    #
    #   0: Didnt hit either constraint, guard or invariant
    #   1: Hit constraint
    #   2: Hit guard
    #   3: Hit invariant surface
    #
    #=
    if !LazySets.API.isdisjoint(newR, guard)
        @show LazySets.API.isdisjoint(newR, constraint)
        println("GUARD TOUCHED")
    end
    if LazySets.API.isdisjoint(newR, guard)
        if LazySets.API.issubset(newR, invariant)
            if LazySets.API.isdisjoint(newR, constraint)
                return 0
            else
                return 1
            end
        else
            return 3
        end
    else
        return 2
    end
    =#
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

    #=
    if LazySets.API.issubset(newR, invariant)
        if LazySets.API.isdisjoint(newR, constraint)
            if LazySets.API.isdisjoint(newR, guard)
                return 0
            else
                return 2
            end

        else
            return 1
        end
    else
        #@show LazySets.API.isdisjoint(newR, guard)
        return 3
    end
    =#

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
function fast_guards_check_t(newR::Zonotope, input, guardProjVectors, guardProjBounds) #; solver=model
    #return all((input + ρ(x, newR)) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) && any((input + -ρ(-x, newR)) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) && all((input - ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))
    #return !someinside(newR, input, guardProjVectors, guardProjBounds)
    return !someinside(minkowski_sum(newR, input), guardProjVectors, guardProjBounds)

end
function touchesCheck(newR::Zonotope, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, invarientProjVectors, invarientProjBounds) #; solver=model
    #cache = map((input, x, y) -> input + ρ(x, newR) <= y, Sρ, constraintProjVectors, constraintProjBounds)
    #abssum = 0.0
    #c = newR.center
    #G = transpose(genmat(newR))
    #tc = Vector{Float64}(undef, size(G, 1))
    #a = sum(abs, transpose(a) * G)
    #res = all((input + ρ(x, newR)) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))#all(input + tsupfunc(x, tc, abssum, c, G) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))

    #res = res && any((input + -ρ(-x, newR)) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds))#any((input - tsupfunc(-x, tc, abssum, c, G)) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds))

    #res = res && all((input - ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))#all((input - tsupfunc(-x, tc, abssum, c, G)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))
    #(all((input + ρ(x, discritezationDict[currentTimeStep])) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) &&
    #any(sign(y) >= 0 ? (input + ρ(-x, newR)) <= y : !((input + ρ(x, newR)) < y) for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) &&
    #any((input + -ρ(-x, discritezationDict[currentTimeStep])) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) &&
    #all((input - ρ(-x, discritezationDict[currentTimeStep])) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds)))
    #return all((input + ρ(x, newR)) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) && any((input + -ρ(-x, newR)) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) && all((input - ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))
    return allunder(newR, constraintProjVectors, constraintProjBounds) && someunder(newR, guardProjVectors, guardProjBounds) && someunder(newR, invarientProjVectors, invarientProjBounds)#all((input - ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))

end

function tsupfunc(d, Z, Ab, c, G)
    #c = center(Z)
    #G = genmat(Z)
    mul!(Z, G, d)
    Ab = sum(abs, Z)
    #c = dot(center(Z), d)
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
    #res = true
    res = !LazySets.API.isdisjoint(set, HPolyhedron(collect(LazySets.HalfSpace(a, b+c) for (a, b, c) in zip(dir, bound, offset))))
    #@inbounds @simd for i in eachindex(offset, dir, bound)
    #    res = res && (offset[i] - ρ(-dir[i], set) <= bound[i])
    #end
    return res
end

function someoutside(set, offset, dir, bound)
    #=res = false
    @inbounds @simd for i in eachindex(offset, dir, bound)
        res = res || (offset[i] - ρ(-dir[i], set) > bound[i])
    end=#
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
    #=res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-dir[i], set) <= bound[i])
    end=#
    return !LazySets.API.isdisjoint(set, HPolyhedron(collect(LazySets.HalfSpace(a, b) for (a, b) in zip(dir, bound))))
end

function someoutside(set, input::LazySet, dir, bound)
    res = false
    @inbounds @simd for i in eachindex(dir, bound)
        res = res || (-ρ(-dir[i], set)-ρ(-dir[i], input) > bound[i])
    end
    return res
    #return !allunder(set, dir, bound)
end

function someinside(set, input::LazySet, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-dir[i], set)+(-ρ(-dir[i], input)) <= bound[i])
    end
    return res
    #return !allunder(set, dir, bound)
end

function someinside(set, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-dir[i], set) <= bound[i])
    end
    return res
    #return !allunder(set, dir, bound)
end

function someinside_transparent(set, input::LazySet, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-dir[i], set)+(-ρ(-dir[i], input)) <= bound[i])
        @show (-ρ(-dir[i], set), -ρ(-dir[i], input), -ρ(-dir[i], set) - ρ(-dir[i], input), bound[i])
    end
    return res
    #return !allunder(set, dir, bound)
end

function ReACTedShieldingK(system::EuclideanHybridSystem, p, k, granularity, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, clustering=true)
    grid = Grid(system.statespace, granularity)
    unsafeDict = Dict{CartesianIndex,Vector{LazySet}}()

    zonotopeArray = initialize_zonotope_array(grid)

    mark_dead_cells!(grid, system.globalConstraints[1], unsafeDict)

    #
    #   Take union with frontier of guards( and Act?)
    #
    #unsafeCells, cFrontierCells = get_contained_edge_cells(grid, system.globalConstraints[1])
    #cFrontierN = grow_indices(grid, map(x -> car2vec(x.id), cFrontierCells))
    #@show frontierCells
    #guardCells, gFrontierCells = get_contained_edge_cells(grid, system.edges[1].guard)
    #gFrontierCells = get_contained_perimeter_cells(grid, system.edges[1].guard)
    #@show guardCells
    #@show gFrontierCells
    #gFrontierN = grow_indices(grid, map(x -> car2vec(x.id), gFrontierCells))
    #@show gFrontierN
    phiDicts = Dict()
    tphiDicts = Dict()

    #=
    for x in hybridSystem.locations
        pd, tpd, id = PhiInputDict(x, δ⁻, δ⁺, alg, maxOrder, reduceOrder)
        phiDicts[x.id] = pd
        tphiDicts[x.id] = tpd
        inputDicts[x.id] = id
        timeConstraintDict[x.id] = []
        constraintDict[x.id] = vcat(x.constraints, constraint)
    end
    =#

    phiDict = PhiDict(system.flowMatrix, δ⁻, δ⁺, alg)
    tPhiDict = Dict(collect((k, permutedims(copy(v))) for (k, v) in pairs(phiDict)))
    #frontierZonotopes = zonotopeArray[collect(CartesianIndex(c) for c in [gFrontierN; cFrontierN])]#zonotopeArray[collect(c.id for c in frontierN)]
    #push!(frontierZonotopes, zonotopeArray[40, 40, 1])
    #reachtimes = map(z -> propagate_set(z, [0.0, p], δ⁻, δ⁺, system, phiDict, tPhiDict, alg, maxOrder, reduceOrder), frontierZonotopes)
    reachtimes = map(z -> propagate_set(z, [0.0, p], δ⁻, δ⁺, system, phiDict, tPhiDict, alg, maxOrder, reduceOrder), zonotopeArray)

    #@time newReACTDiscretizePlus(zonotopeArray[40, 40, 1], δ⁻, δ⁺, system.flowMatrix, phiDict, alg, maxOrder, reduceOrder)
    #@time propagate_set(zonotopeArray[40, 40, 1], [0.0, p], δ⁻, δ⁺, system, phiDict, tPhiDict, alg, maxOrder, reduceOrder)
    #@show reachtimes
    return reachtimes

    waitinglist = []
    zenoBound = 20
    transitionCount = 0

    dims = size(X0.center, 1)

    #res = auxReACTed(hybridSystem, hybridSystem.locations[loc], nothing, interval, X0, dirsVectors, constraint, δ⁻, δ⁺, flowPhiDict, alg, maxOrder, reduceOrder, missing, clustering, timeConstraintList, saveResult)

    #reachset = vcat(reachset, res)

    #return res
    push!(waitinglist, (loc, X0, nothing, interval))


    if TIMEFUNC
        println((time_ns() - startTimer) / 10^9)
    end

    while !isempty(waitinglist)
        #GC.gc()

        location, initialset, edge, interval′ = pop!(waitinglist)

        _ = auxReACTed(waitinglist, hybridSystem, hybridSystem.locations[location], edge, interval′, initialset, constraintDict[location], δ⁻, δ⁺, phiDicts[location], tphiDicts[location], inputDicts[location], dims, alg, maxOrder, reduceOrder, clustering, timeConstraintDict[location])


        if transitionCount < zenoBound
            transitionCount += 1
        else
            throw(error("Reached zeno bound"))

        end
    end

    if TIMEFUNC
        total = (time_ns() - startTimer) / 10^9
        println("Total times: 
        Total: $(total)
        TotalNotAux: $((totalDiscTime / 10^9) + (totalGuardTime / 10^9) + (totalTouchesTime / 10^9) + (totalReACTTime / 10^9))
        Disc: $(totalDiscTime / 10^9)
        Guards: $(totalGuardTime / 10^9) 
        Init Guards: $(totalInitGuardTime / 10^9) 
        Touches: $(totalTouchesTime / 10^9)
        ReACT: $(totalReACTTime / 10^9)
        Init ReACT: $(totalInitReACTTime / 10^9)")
        println("Amount of calls to auxreacted: $transitionCount")
    end

    return []
end
# AucReacted is called recursively each time we have a new starting location (after a transition)
function auxReACTedShielding(waitinglist, hybridSystem, loc::Location, currentEdge, interval, X0, constraint, δ⁻::Float64, δ⁺::Float64, PhiDict, TPhiDict, inputDiscritezationDict, dims, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, clustering=true, activeTimeConstraints=[])


    discretizationDict = newReACTDiscretizePlus(loc, X0, δ⁻, δ⁺, PhiDict, inputDiscritezationDict, alg, maxOrder, reduceOrder)
    #t2 = Base.time()
    #@show t2 - t1
    if TIMEFUNC
        global totalDiscTime += time_ns() - recorderTimeStart
    end


    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    #setOfConstraints = constraint #vcat(loc.constraints, constraint)


    # For each edge we simulate the system
    listOfEdges = loc.edges
    if !isa(currentEdge, Nothing) # This allows us to sometimes just wanna do a signular edge
        listOfEdges = [currentEdge]
    end
    #tΦ::Matrix{Float64} = diagm(ones(Float64, size(loc.A, 2)))
    for edge in listOfEdges
        # println("Handling edge at time $time from $(loc.id) -> $(edge.targetLoc)")
        #guards = edge.guard # Technically the guard is one singular HPolyhedron, but it composes the other guards


        #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 

        if TIMEFUNC
            global recorderTimeStart = time_ns()
        end

        reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, interval, edge.guard, constraint, dims, 2, PhiDict, TPhiDict, discretizationDict, inputDiscritezationDict)

        if TIMEFUNC
            global totalGuardTime += time_ns() - recorderTimeStart
        end

        #=
        if reachtime - time == 0.0
            println("zeno?")
            if !isempty(reachset)
                # println("Found an immediate transition to $(edge.targetLoc), but we are not taking it as we are scared of zeno behaviour")
            end
            if saveResult
                push!(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * " took no steps"))
            end
            #continue
        end
        =#
        #println("Everytime")
        if reachtime < endtime


            # println("Current time: $time, intersecting time start: $reachtime")
            if TIMEFUNC
                recorderTimeStart = time_ns()
            end

            tryContinueFlag, timeNotIntersected, intersectingSetsList = ReACTTouches(loc, δ⁻, [reachtime, endtime], (reachtime - time), edge.guard, constraint, PhiDict, discretizationDict[δ⁻], inputDiscritezationDict, tΦ, reduceOrder, maxOrder)

            if TIMEFUNC
                global totalTouchesTime += time_ns() - recorderTimeStart
            end

            if any((timeNotIntersected >= x) for x in activeTimeConstraints)
                handleHitConstraint(timeNotIntersected, loc.id)
            end





            #latestSet::Any = nothing
            #
            #   Here we should check whether we have reached endtime. If true we should only push the jumpSet
            #   Still need to check whether we have reached the invariant. If true we should NOT push the else branch result, only the tempReachsets[infMaxsIdx]
            #
            if !isempty(intersectingSetsList)
                if tryContinueFlag
                    #println("Continuing from previous run at time $timeNotIntersected $reachtime $(length(intersectingSetsList))")

                    #latestSet = zonotopeStripIntersection(latestSet, loc.invarient)
                    #push!(waitinglist, (edge.targetLoc, zonotopeStripIntersection(last(intersectingSetsList), loc.invarient), edge, [timeNotIntersected, endtime])) # gav tϕ foer..

                    push!(waitinglist, (edge.targetLoc, last(intersectingSetsList), edge, [timeNotIntersected, endtime])) # gav tϕ foer..
                    #branchedRun = auxReACTed(hybridSystem, loc, edge, [timeNotIntersected, endtime], latestSet, dirs, constraint, δ⁻, δ⁺, FlowPhiDict, alg, maxOrder, reduceOrder, tΦ, clustering, timeConstraintList, saveResult)
                    #if saveResult
                    #    reachset = vcat(reachset, branchedRun)
                    #end
                end



                # we first save this for later
                #if tryContinueFlag
                #    latestSet = copy(last(intersectingSetsList))
                #end
                # Afterwards we process the intersecting set
                #intersectedSet = nothing

                jumpSetsList = Vector{Zonotope}()


                for set in intersectingSetsList
                    jumpSet = set
                    #jumpSet = zonotopeStripIntersection(set, edge.guard)#set

                    if !isa(loc.invarient, Nothing)
                        jumpSet = zonotopeStripIntersection(jumpSet, loc.invarient)
                    end

                    #any(x -> ρ(x.a, jumpSet) > x.b, setOfConstraints) && handleHitConstraint(reachtime + δ⁻ * i, loc.id)
                    #=
                    if LazySets.isdisjoint(jumpSet, edge.guard)
                        break
                    end
                    =#
                    if !isa(edge.guard, Nothing)
                        jumpSet = zonotopeStripIntersection(jumpSet, edge.guard)
                    end

                    # Check invarient
                    # Push to reachset


                    # Apply jump Matrix
                    jumpSet = linear_map(edge.jumpMatrix, jumpSet)
                    if !isnothing(edge.jumpVector)
                        LazySets.translate!(jumpSet, edge.jumpVector)
                    end

                    if !isa(hybridSystem.locations[edge.targetLoc].invarient, Nothing)
                        jumpSet = zonotopeStripIntersection(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)
                        #=
                        if !isdisjoint(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)
                            jumpSet = zonotopeStripIntersection(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)
                        else


                            continue
                        end
                        =#
                    end


                    push!(jumpSetsList, jumpSet)
                end

                len = length(jumpSetsList)

                if clustering

                    if len > 1
                        tempIntersect = foldl(ConvexHull, jumpSetsList)
                        intersectedSet = convert(Zonotope, box_approximation(tempIntersect))
                        push!(waitinglist, (edge.targetLoc, intersectedSet, nothing, [reachtime, endtime]))

                    else
                        intersectedSet = jumpSetsList[1]
                        push!(waitinglist, (edge.targetLoc, intersectedSet, nothing, [reachtime, endtime]))

                    end

                    #jumpSet = intersectedSet
                    # println("Finished intersections")

                    #branchedRun = auxReACTed(hybridSystem, hybridSystem.locations[edge.targetLoc], nothing, [reachtime, endtime], jumpSet, dirs, constraint, δ⁻, δ⁺, FlowPhiDict, alg, maxOrder, reduceOrder, tΦ, clustering, timeConstraintList, saveResult)

                    #if saveResult
                    #    reachset = vcat(reachset, branchedRun)
                    #end
                else # No clustering. Do individual runningset
                    for (index, jumpSet) in enumerate(jumpSetsList)
                        # Note that we only do steps of size δ⁻ in touches
                        #timeStart = reachtime + (index - 1) * δ⁻
                        push!(waitinglist, (edge.targetLoc, jumpSet, nothing, [reachtime + (index - 1) * δ⁻, endtime])) # gav tϕ foer..
                        #branchedRun = auxReACTed(hybridSystem, hybridSystem.locations[edge.targetLoc], nothing, [timeStart, endtime], jumpSet, dirs, constraint, δ⁻, δ⁺, FlowPhiDict, alg, maxOrder, reduceOrder, tΦ, clustering, timeConstraintList, saveResult)

                        #if saveResult
                        #    reachset = vcat(reachset, branchedRun)
                        #end
                    end
                end

            end

            # If we are not encountering an invarient, try continue



        else # Reached the end time
            # println("Reached end time before intersecting guard")
            # Check if we hit activeTimeConstraints
            if any((reachtime >= x) for x in activeTimeConstraints)
                handleHitConstraint(reachtime, loc.id)
            end

        end
    end

    if isempty(loc.edges) # This means it is just a continous system from here
        # println("No edges? call ReACT")

        #tempReachset, tempInput, reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, [time, endtime], guards, setOfConstraints, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, missing, saveResult)


        #tempReachset, _, _ = ReACT(loc, δ⁻, δ⁺, interval, constraint, 2, alg, maxOrder, reduceOrder, PhiDict[loc.id])
        if TIMEFUNC
            recorderTimeStart = time_ns()
        end
        reachtime = ReACT(loc, δ⁻, δ⁺, [time, endtime], constraint, 2, TPhiDict, discretizationDict, inputDiscritezationDict)
        if TIMEFUNC
            global totalReACTTime += time_ns() - recorderTimeStart
        end
    end


end
#=
function ReACT_guards(δ⁻::Float64, δ⁺::Float64, interval, statespace, guards, constraint, invariant, STRATEGY::Integer, PhiDict, permutedphiDict, discritezationDict, inputDiscretizationDict)

    #
    #   These functions should be evaluated outside ReACT_guards and passed as parameters
    #
    #constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    #guardProjVectors, guardProjBounds = getHalfSpaceProjections(guards)
    #invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(invariant)

    dims = LazySets.API.dim(statespace)

    time::Float64 = copy(minimum(interval))
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(δ⁻)

    attemptsRecorder = zeros(4)

    Φ::Matrix{Float64} = diagm(ones(Float64, dims))

    tempM = diagm(ones(Float64, dims))

    #   Because we assume an empty input set i have trouble seeing how these can become relevant.
    #   But i am keeping them for now to reuse as much as possible
    #=
    Sρ = zeros(Float64, length(constraintProjVectors))
    Gρ = zeros(Float64, length(guardProjVectors))
    Iρ = zeros(Float64, length(invarientProjVectors))
    =#

    i = 1
    input = Zonotope(zeros(dims), zeros(size(genmat((inputDiscretizationDict[δ⁻])))))
    while time < endtime

        attempts = 1
        approveFlag = false
        check = 0
        while !approveFlag
            # Handle if we can no longer reduce the reachset (we keep hitting something)
            if currentTimeStep >= δ⁻
                #linear_map!(Z, Φ, discritezationDict[currentTimeStep])
                check = guardCheck(input, Φ, discritezationDict[currentTimeStep], constraint, guards, invariant)
                if 0 == check
                    #if 0 == guardCheck(discritezationDict[currentTimeStep], Sρ, constraintProjVectors, constraintProjBounds, Gρ, guardProjVectors, guardProjBounds, Iρ, invarientProjVectors, invarientProjBounds)
                    approveFlag = true
                    if currentTimeStep + time < endtime
                        input = input + linear_map(Φ, inputDiscretizationDict[currentTimeStep])
                        mul!(tempM, Φ, PhiDict[currentTimeStep])
                        copy!(Φ, tempM)

                    else
                        # TODO - The input gets used or referenced somewhere even with the zero flag..
                        #d = max(2^(floor(log2((endtime - time)/δ⁻))) * δ⁻, δ⁻)
                        #@show (endtime - time, d)
                        #input = linear_map(PhiDict[d], input) + inputDiscretizationDict[d]
                        #mul!(tempM, Φ, PhiDict[d])
                        #copy!(Φ, tempM)
                        #time = time + d
                        return (endtime, discritezationDict[δ⁻], input, 0)
                    end
                else

                    currentTimeStep /= 2
                    attempts += 1
                end
            else
                #@show check
                return (time, linear_map(Φ, discritezationDict[δ⁻]), input, check)
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
    #println("ever")
    return (endtime, linear_map(Φ, discritezationDict[δ⁻]), input, 0)
end
=#