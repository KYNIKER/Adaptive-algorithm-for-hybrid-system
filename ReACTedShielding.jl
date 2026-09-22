using LazySets, LinearAlgebra, ReachabilityAnalysis
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
    #@show grid.lower
    lower_offset = grid.lower .+ (granularity/2)
    for idx in CartesianIndices(grid.deadCells)
        of = lower_offset + ((car2vec(idx) .- 1) .* granularity)
        LazySets.API.translate!(Z, of)
        #@show Z, idx, of
        for act in system.Act
            if !LazySets.API.isdisjoint(act.guard, Z)
                reachable_by_action[(act, idx)] = get_touching_cell_idxs(grid, LazySets.API.translate(linear_map(act.jumpMatrix, Z), act.jumpVector))
            end
        end

        z, f = propagate_set(Z, [0.0, p], δ⁻, δ⁺, system, P2A_abs, phiDict, tPhiDict, inputDict, alg, maxOrder, reduceOrder)

        #=if !LazySets.API.isdisjoint(z, system.edges[1].guard)
            dc += 1
        end=#

        if f != 1
            #@show LazySets.API.isdisjoint(Z, system.globalConstraints[1]), Z, idx
            #grid.deadCells[idx] = true

            #println("so deads")
            reachable_by_flow[idx] = get_touching_cell_idxs(grid, z)
            #else
            #grid.array[idx].pCells = get_touching_cell_idxs(grid, z)
        end
        #grid.array[idx].sidx[1] = copy(f)
        LazySets.API.translate!(Z, -of)
    end
    #@show dc
    return grid, reachable_by_action, reachable_by_flow

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

#
#   Right now this makes quite a bit of allocations.
#
function propagate_set(X0, interval, δ⁻::Float64, δ⁺::Float64, system, Φ₂, PhiDict, TPhiDict, inputDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5)

    discretizationDict = newReACTDiscretizePlus(X0, δ⁻, δ⁺, system.flowMatrix, Φ₂, PhiDict, inputDict, alg, maxOrder, reduceOrder)
    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    constraint = system.globalConstraints[1]
    listOfEdges = system.edges
    for edge in listOfEdges
        #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 
        reachtime, newSet, inputSet, flag = ReACT_guards(δ⁻, δ⁺, interval, system.statespace, edge.guard, constraint, hyperrectangle_to_HPolyhedron(system.statespace), 2, PhiDict, TPhiDict, discretizationDict, inputDict)


        #
        #   0: Didnt hit either constraint, guard or invariant
        #   1: Hit constraint
        #   2: Hit guard
        #   3: Hit invariant surface
        #
        if flag == 0
            return (minkowski_sum(newSet, inputSet), 0)
        elseif flag == 1 # Hit constraint
            return (nothing, 1)
        elseif flag == 2
            #@show flag
            #@show inputSet
            issubset, intersectingSetsList, t = ReACT_touches_set_constant_input(newSet, inputSet, edge.guard, δ⁻, [reachtime, endtime], PhiDict, inputDict)
            tempIntersect = foldl(ConvexHull, intersectingSetsList)
            intersectedSet = convert(Zonotope, box_approximation(tempIntersect))
            jumpSet = linear_map(edge.jumpMatrix, zonotopeStripIntersection(intersectedSet, edge.guard))
            if !isnothing(edge.jumpVector)
                LazySets.translate!(jumpSet, edge.jumpVector)
            end
            #res = propagate_set(jumpSet, [reachtime, endtime], δ⁻, δ⁺, system, Φ₂, PhiDict, TPhiDict, inputDict, alg, maxOrder, reduceOrder)
            #@show res
            return propagate_set(jumpSet, [reachtime, endtime], δ⁻, δ⁺, system, Φ₂, PhiDict, TPhiDict, inputDict, alg, maxOrder, reduceOrder)
        end
        return (newSet, 3)




        if reachtime < endtime

            tryContinueFlag, timeNotIntersected, intersectingSetsList = ReACTTouches(loc, δ⁻, [reachtime, endtime], (reachtime - time), edge.guard, constraint, PhiDict, discretizationDict[δ⁻], inputDiscritezationDict, tΦ, reduceOrder, maxOrder)

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

function ReACT_touches_set_constant_input(X0, U, set, δ⁻::Float64, interval, PhiDict, inputDiscretizationDict)
    # Note that in touches we always use δ⁻
    # That is, we do not adjust timestep sizes

    time::Float64 = interval[1] + δ⁻
    endtime::Float64 = interval[2]

    intersectingSetsList = Vector{Zonotope}() # This will store all of our sets


    X0 = linear_map(PhiDict[δ⁻], X0)
    U = linear_map(PhiDict[δ⁻], U) + inputDiscretizationDict[δ⁻]
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
    return (LazySets.API.issubset(tempSet, set), intersectingSetsList, time)
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

function ReACT_guards(δ⁻::Float64, δ⁺::Float64, interval, statespace, guards, constraint, invariant, STRATEGY::Integer, PhiDict, permutedphiDict, discritezationDict, inputDiscretizationDict)

    #
    #   These functions should be evaluated outside ReACT_guards and passed as parameters
    #
    #constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    #guardProjVectors, guardProjBounds = getHalfSpaceProjections(guards)
    #invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(invariant)

    dims = LazySets.API.dim(statespace)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(δ⁺)

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
    input = Zonotope(zeros(dims), zeros(dims, dims))
    while time < endtime

        attempts = 1
        approveFlag = false
        check = 0
        while !approveFlag
            # Handle if we can no longer reduce the reachset (we keep hitting something)
            if currentTimeStep >= δ⁻
                #linear_map!(Z, Φ, discritezationDict[currentTimeStep])
                check = guardCheck(minkowski_sum(input, linear_map(Φ, discritezationDict[currentTimeStep])), constraint, guards, invariant)
                if 0 == check
                    #if 0 == guardCheck(discritezationDict[currentTimeStep], Sρ, constraintProjVectors, constraintProjBounds, Gρ, guardProjVectors, guardProjBounds, Iρ, invarientProjVectors, invarientProjBounds)
                    approveFlag = true
                    input = input + linear_map(Φ, inputDiscretizationDict[currentTimeStep])
                    mul!(tempM, Φ, PhiDict[currentTimeStep])
                    copy!(Φ, tempM)
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

    return (time, linear_map(Φ, discritezationDict[currentTimeStep]), input, 0)
end


function ReACT(loc, δ⁻::Float64, δ⁺::Float64, interval, constraint, STRATEGY::Integer, permutedphiDict, discritezationDict, inputDiscritezationDict)
    # We calculate the reachset till we reach a guard for an intersection (or till failure)
    #println("HUH")

    init1 = TIMEFUNC ? time_ns() : 0.0


    # constraintProjVectors = map(x -> x.a, constraint)
    # constraintProjBounds = ρ.(constraintProjVectors, constraint)
    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)


    #=
    permutedphiDict = Dict()
    for key in keys(phiDict)
        permutedphiDict[key] = permutedims(phiDict[key])
    end
    =#

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(δ⁻)#copy(initialTimeStep)

    attemptsRecorder = zeros(4)



    #Φ::Matrix{Float64} = diagm(ones(Float64, size(loc.A, 2)))

    #tempM = diagm(ones(Float64, size(loc.A, 2)))
    #ϕt = similar(tempM)
    #pϕt = similar(tempM)
    #V = copy(inputDiscritezationDict[initialTimeStep])

    #lastVs = copy(Vs)
    Sρ = zeros(Float64, length(constraint))
    Iρ = zeros(Float64, length(invarientProjVectors))

    #newR = discritezationDict[initialTimeStep]
    i = 1

    if TIMEFUNC
        global totalInitReACTTime += time_ns() - init1
    end

    #U = inputDiscritezationDict[0]
    while time < endtime
        # println("Time iss: $time")
        attempts = 1
        approveFlag = false

        while !approveFlag

            # Handle if we can no longer reduce the reachset (we keep hitting something)
            if currentTimeStep >= δ⁻

                #=
                if changedTimeStep
                    newR = discritezationDict[currentTimeStep]

                    V = copy(inputDiscritezationDict[currentTimeStep])
                    ϕt = phiDict[currentTimeStep]
                    pϕt = permutedphiDict[currentTimeStep]
                end
                =#



                if allunder(discritezationDict[currentTimeStep], Sρ, constraintProjVectors, constraintProjBounds) &&
                    someunder(discritezationDict[currentTimeStep], Iρ, invarientProjVectors, invarientProjBounds)



                    approveFlag = true

                    Sρ += map(x -> ρ(x, inputDiscritezationDict[currentTimeStep]), constraintProjVectors)
                    Iρ += map(x -> ρ(x, inputDiscritezationDict[currentTimeStep]), invarientProjVectors)


                    map!(x -> permutedphiDict[currentTimeStep] * x, constraintProjVectors)
                    map!(x -> permutedphiDict[currentTimeStep] * x, invarientProjVectors)


                else
                    #newR = copy(newR)
                    currentTimeStep = currentTimeStep / 2

                    attempts = attempts + 1
                end
            else# If we hit a constraint
                if !someoutside(discritezationDict[currentTimeStep], Sρ, constraintProjVectors, constraintProjBounds)

                    return time
                else
                    handleHitConstraint(time, loc.id)
                end
            end
        end

        attemptsRecorder[(i%4)+1] = attempts
        i = i + 1
        time = time + currentTimeStep
        # Reset / apply strategy
        # Only do this if the current timestep is less than the initial
        if STRATEGY == 2
            # If attemptsrecorder past 4 are successes, double timestep
            if currentTimeStep < δ⁺
                #lowest = min(4, i - 1)
                #window = @view attemptsRecorder[i-lowest:i-1]
                if all(==(1), attemptsRecorder)
                    currentTimeStep = currentTimeStep * 2
                    #changedTimeStep = true
                end

            end
        elseif STRATEGY == 1
            # always try double
            if currentTimeStep < δ⁺
                currentTimeStep = currentTimeStep * 2

            end
        elseif STRATEGY == 0
            # Only reduce

        end
    end

    return time
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

function guardCheck(newR::Zonotope, constraint, guard, invariant) #; solver=model
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

function someunder(set, dir, bound)
    #=res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-dir[i], set) <= bound[i])
    end=#
    return !LazySets.API.isdisjoint(set, HPolyhedron(collect(LazySets.HalfSpace(a, b) for (a, b) in zip(dir, bound))))
end

function someoutside(set, dir, bound)
    #=res = false
    @inbounds @simd for i in eachindex(dir, bound)
        res = res || (-ρ(-dir[i], set) > bound[i])
    end=#
    return !allunder(set, dir, bound)
end


