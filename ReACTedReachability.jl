using LazySets, LinearAlgebra, ReachabilityAnalysis
include("Discretize.jl")
include("Utilities.jl")

export ReACTedFast

const TIMEFUNC = false


totalReACTTime = 0
totalDiscTime = 0
totalGuardTime = 0
totalTouchesTime = 0
totalInitGuardTime = 0.0
totalInitReACTTime = 0.0
recorderTimeStart = 0.0
startTimer = 0.0

#=
function zonotopePrintDim(Z::Zonotope, dim)
    center = Z.center[dim]
    G = genmat(Z)
    _, genAmount = size(G)
    genContribute = sum(abs(G[dim, i]) for i in 1:genAmount)
    return "center: $center, generators: $genContribute"
end
=#

function ReACTedFast(hybridSystem::HybridSystemV2, initialLoc, interval, X0, U, constraint, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, clustering=true, timeConstraintList=[])
    loc = initialLoc
    #flowPhiDict = Dict(map(x -> x.id => PhiInputDict(x, δ⁻, δ⁺, alg, maxOrder, reduceOrder), hybridSystem.locations))
    if TIMEFUNC
        global startTimer = time_ns()
    end
    #phiDicts, inputDicts = Dict(map(x -> x.id => PhiInputDict(x, δ⁻, δ⁺, alg, maxOrder, reduceOrder), hybridSystem.locations)[1, :])
    phiDicts = Dict()
    tphiDicts = Dict()
    inputDicts = Dict()
    timeConstraintDict = Dict()
    constraintDict = Dict()
    for x in hybridSystem.locations
        pd, tpd, id = PhiInputDict(x, δ⁻, δ⁺, alg, maxOrder, reduceOrder)
        phiDicts[x.id] = pd
        tphiDicts[x.id] = tpd
        inputDicts[x.id] = id
        timeConstraintDict[x.id] = []
        constraintDict[x.id] = vcat(x.constraints, constraint)
    end


    waitinglist = []
    zenoBound = 2000
    transitionCount = 0

    dims = size(X0.center, 1)

    #res = auxReACTed(hybridSystem, hybridSystem.locations[loc], nothing, interval, X0, dirsVectors, constraint, δ⁻, δ⁺, flowPhiDict, alg, maxOrder, reduceOrder, missing, clustering, timeConstraintList, saveResult)

    #reachset = vcat(reachset, res)

    #return res
    push!(waitinglist, (loc, X0, nothing, interval))
    for (id, time) in timeConstraintList
        timeConstraintDict[id] = [time]
    end

    if TIMEFUNC
        println((time_ns() - startTimer) / 10^9)
    end


    projVectors = Vector{Vector{Float64}}()
    projBounds = Vector{Float64}()
    projRanges = Vector{UnitRange{Int64}}()

    prevLocation = -1
    while !isempty(waitinglist)
        #GC.gc()

        location, initialset, edge, interval′ = pop!(waitinglist)

        locationObject = hybridSystem.locations[location]

        # if prevLocation != location # Update constraint and buffer
        #     supportProj[1], supportBounds[1] = getHalfSpaceProjections(constraintDict[location])
        #     supportProj[2], supportBounds[2] = getHalfSpaceProjections(locationObject.invarient)
        # end

        if prevLocation != location # Update constraints
            constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraintDict[location])
            invariantProjVectors, invariantProjBounds = getHalfSpaceProjections(locationObject.invarient)
            projVectors = vcat(constraintProjVectors, invariantProjVectors)
            projBounds = vcat(constraintProjBounds, invariantProjBounds)

            amountOfConstraints = length(constraintProjBounds)
            amountOfInvariants = length(invariantProjBounds)
            projRanges = [1:amountOfConstraints, amountOfConstraints+1:(amountOfConstraints+amountOfInvariants)]
        end

        _ = auxReACTed(waitinglist, hybridSystem, locationObject, edge, interval′, initialset, constraintDict[location], δ⁻, δ⁺, phiDicts[location], tphiDicts[location], inputDicts[location], dims, projVectors, projBounds, projRanges, alg, maxOrder, reduceOrder, clustering, timeConstraintDict[location])


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
function auxReACTed(waitinglist, hybridSystem, loc::Location, currentEdge, interval, X0, constraint, δ⁻::Float64, δ⁺::Float64, PhiDict, TPhiDict, inputDiscritezationDict, dims, projVectors, projBounds, projRanges, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, clustering=true, activeTimeConstraints=[])


    if TIMEFUNC
        global recorderTimeStart = time_ns()
    end

    #PhiDict, inputDiscritezationDict = FlowPhiDict[loc.id]

    #t1 = Base.time()
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

        reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, interval, edge.guard, constraint, dims, 2, PhiDict, TPhiDict, discretizationDict, inputDiscritezationDict, projVectors, projBounds, projRanges)

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

            tryContinueFlag, timeNotIntersected, intersectingSetsList = ReACTTouches(loc, δ⁻, [reachtime, endtime], (reachtime - time), edge.guard, constraint, PhiDict, discretizationDict[δ⁻], inputDiscritezationDict, tΦ, reduceOrder, maxOrder, projVectors, projBounds, projRanges)

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

                # intersects = !isa(loc.invarient, Nothing) ? 
                #     (!isa(edge.guard, Nothing) ? HPolyhedron(vcat(loc.invarient.constraints, edge.guard.constraints)) : loc.invarient) : 
                #     (!isa(edge.guard, Nothing) ? edge.guard : nothing)
                    

                for set in intersectingSetsList
                    jumpSet = set
                    #jumpSet = zonotopeStripIntersection(set, edge.guard)#set

                    # if !isa(intersects, Nothing)
                    #     jumpSet = zonotopeStripIntersection(jumpSet, intersects)
                    # end
                    
                    if !isa(loc.invarient, Nothing)
                        jumpSet = zonotopeStripIntersection(jumpSet, loc.invarient)
                    end

                    #any(x -> ρ(x.a, jumpSet) > x.b, setOfConstraints) && handleHitConstraint(reachtime + δ⁻ * i, loc.id)

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
        reachtime = ReACT(loc, δ⁻, δ⁺, [time, endtime], constraint, 2, TPhiDict, discretizationDict, inputDiscritezationDict, projVectors, projBounds, projRanges)
        if TIMEFUNC
            global totalReACTTime += time_ns() - recorderTimeStart
        end
    end


end

function ReACTTouches(loc, δ⁻::Float64, interval, initialTime::Float64, guard, constraint, PhiDict, newR, inputDiscritezationDict, Φ, reduce_order, max_order, projVectors, projBounds, projRanges)
    # Note that in touches we always use δ⁻
    # That is, we do not adjust timestep sizes

    #phiDict = PhiDict
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(guard) 
    allProjVectors = vcat(projVectors, guardProjVectors)
    allProjBounds = vcat(projBounds, guardProjBounds)
    amountBeforeGuards = length(projBounds)
    amountTotal = length(guardProjBounds) + amountBeforeGuards

    constraintProjVectors = view(allProjVectors, projRanges[1])
    constraintProjBounds = view(allProjBounds, projRanges[1])
    invarientProjVectors = view(allProjVectors, projRanges[2])
    invarientProjBounds = view(allProjBounds, projRanges[2])
    guardProjVectors = view(allProjVectors, amountBeforeGuards+1:amountTotal)
    guardProjBounds = view(allProjBounds, amountBeforeGuards+1:amountTotal)

    # constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    # guardProjVectors, guardProjBounds = getHalfSpaceProjections(guard)
    # invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)

    time::Float64 = interval[1]
    endtime::Float64 = interval[2]


    Vs = nestedInputDiscCalculate(inputDiscritezationDict, PhiDict, δ⁻, initialTime, reduce_order, max_order)
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
    intersectingSetsList = Vector{Zonotope}() # This will store all of our sets

    tempSet = minkowski_sum(linear_map(Φ, newR), Vs)
    while time < endtime

        #if all((ρ(x, tempSet)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) && # IsSubSet
        #   all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds)) &&
        #   all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Intersects
        if touchesCheck(tempSet, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, invarientProjVectors, invarientProjBounds)
            #if mapreduce(x -> intersects(newRR, x), &, guard)
            #@show time
            push!(intersectingSetsList, copy(tempSet))

            # Main calculation. No longer changing timesteps
            #Vs = minkowski_sum(Vs, V) #Update input
            #newRR = linear_map(ϕ, newRR)
            #V = linear_map(ϕ, V)
            tempSet = minkowski_sum(linear_map(PhiDict[δ⁻], tempSet), inputDiscritezationDict[δ⁻])
            # i = i + 1
            time += δ⁻

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
    end
    #continueAfter = false # We are at the end time horizon, therefore no continuing
    return (false, time, intersectingSetsList)
end


function ReACTGuards(loc, δ⁻::Float64, δ⁺::Float64, interval, guards, constraint, dims, STRATEGY::Integer, PhiDict, permutedphiDict, discritezationDict, inputDiscritezationDict, projVectors, projBounds, projRanges)
    # We calculate the reachset till we reach a guard for an intersection (or till failure)
    init1 = TIMEFUNC ? time_ns() : 0.0


    #initialTimeStep = copy(δ⁺)
    #m = copy(δ⁻)
    #changedTimeStep = true
    #phiDict = PhiDict
    #A = copy(loc.A)
    # constraintProjVectors = map(x -> x.a, constraint)
    # constraintProjBounds = ρ.(constraintProjVectors, constraint)
    # constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    #guardProjVectors, guardProjBounds = [], []

    # guardProjVectors, guardProjBounds = getHalfSpaceProjections(guards)
    #=
    invarientProjVectors, invarientProjBounds = [], []
    if isnothing(loc.invarient)
        invarientProjVectors, invarientProjBounds = [constraintProjVectors[1]], [Inf]
    else
        invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)
    end
    =#
    # invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)

    #oldConstraintProjVectors = copy(constraintProjVectors)
    #oldGuardProjVectors = copy(guardProjVectors)
    #oldInvarientProjVectors = copy(invarientProjVectors)

    #=
    permutedphiDict = Dict()
    for key in keys(PhiDict)
        permutedphiDict[key] = permutedims(PhiDict[key])
    end
    =#
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(guards) 

    allProjVectors = vcat(projVectors, guardProjVectors)
    allProjBounds = vcat(projBounds, guardProjBounds)
    amountBeforeGuards = length(projBounds)
    amountTotal = length(guardProjBounds) + amountBeforeGuards

    constraintProjVectors = view(allProjVectors, projRanges[1])
    constraintProjBounds = view(allProjBounds, projRanges[1])
    invarientProjVectors = view(allProjVectors, projRanges[2])
    invarientProjBounds = view(allProjBounds, projRanges[2])
    guardProjVectors = view(allProjVectors, amountBeforeGuards+1:amountTotal)
    guardProjBounds = view(allProjBounds, amountBeforeGuards+1:amountTotal)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(δ⁺)

    #overapproximateIntersectingSetArray = []
    #lastNewR = []
    attemptsRecorder = zeros(4)



    Φ::Matrix{Float64} = diagm(ones(Float64, dims))

    tempM = diagm(ones(Float64, dims))
    #ϕt = similar(tempM)
    #pϕt = similar(tempM)
    #V = copy(inputDiscritezationDict[initialTimeStep])

    #lastVs = copy(Vs)
    Sρ = zeros(Float64, length(constraint))
    Gρ = zeros(Float64, length(guardProjVectors))
    Iρ = zeros(Float64, length(invarientProjVectors))

    #newR = discritezationDict[initialTimeStep]
    i = 1

    #U = inputDiscritezationDict[0]
    if TIMEFUNC
        global totalInitGuardTime += time_ns() - init1
    end
    while time < endtime
        #println("Time iss: $time")
        attempts = 1
        approveFlag = false
        while !approveFlag
            #println("Stuck?")

            # Handle if we can no longer reduce the reachset (we keep hitting something)
            if currentTimeStep >= δ⁻
                # If we hit a constraint
                #newRR = concretize(newRR)
                #if all((input + ρ(x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds)) && any((input + ρ(x, newR)) > y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))
                #    handleHitConstraint(time, loc.id)
                #end
                #return (dirVals, time, missing)


                #if changedTimeStep
                #newR = discritezationDict[currentTimeStep]

                #V = copy(inputDiscritezationDict[currentTimeStep])
                #ϕt = phiDict[currentTimeStep]
                #pϕt = permutedphiDict[currentTimeStep]

                #end



                # #println(any((input + ρ(-x, newRR)) > y for (input, x, y) in zip(Sρ, invarientProjVectors, invarientProjBounds)))
                #println(all((input + ρ(x, newRR)) <= y for (input, x, y) in zip(Sρ, invarientProjVectors, invarientProjBounds)))
                #=
                if isnothing(guards) || (all((input + ρ(x, discritezationDict[currentTimeStep])) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) &&
                                     #any(sign(y) >= 0 ? (input + ρ(-x, newR)) <= y : !((input + ρ(x, newR)) < y) for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) &&
                                     any((input + -ρ(-x, discritezationDict[currentTimeStep])) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) &&
                                     all((input - ρ(-x, discritezationDict[currentTimeStep])) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds)))
                                        =#
                if guardCheck(discritezationDict[currentTimeStep], Sρ, constraintProjVectors, constraintProjBounds, Gρ, guardProjVectors, guardProjBounds, Iρ, invarientProjVectors, invarientProjBounds)
                    #(any((input + ρ(-x, newRR)) > y for (input, x, y) in zip(-Sρ, invarientProjVectors, invarientProjBounds)) && all((input + ρ(x, newRR)) <= y for (input, x, y) in zip(Sρ, invarientProjVectors, invarientProjBounds)))
                    #all(((ρ(x, tempSet)) <= y) for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Subset
                    #(isnothing(loc.invarient) || intersects(tempSet, loc.invarient))
                    #if mapreduce(x -> intersects(newRR, x), &, guard)
                    #push!(overapproximateIntersectingSetArray, newRR)

                    #lastVs = copy(Vs)
                    #Vs = ReachabilityAnalysis.Exponentiation.Φ₁(A, time - minimum(interval), ReachabilityAnalysis.Exponentiation.BaseExp) * U

                    approveFlag = true
                    Sρ += map(x -> ρ(x, inputDiscritezationDict[currentTimeStep]), constraintProjVectors)
                    Gρ += map(x -> ρ(x, inputDiscritezationDict[currentTimeStep]), guardProjVectors)
                    Iρ += map(x -> ρ(x, inputDiscritezationDict[currentTimeStep]), invarientProjVectors)
                    #constraintProjVectors = map(x -> permutedphiDict[currentTimeStep] * x, oldConstraintProjVectors)
                    #guardProjVectors = map(x -> permutedphiDict[currentTimeStep] * x, oldGuardProjVectors)
                    #invarientProjVectors = map(x -> permutedphiDict[currentTimeStep] * x, oldInvarientProjVectors)

                    #oldConstraintProjVectors = constraintProjVectors
                    #oldGuardProjVectors = guardProjVectors
                    #oldInvarientProjVectors = invarientProjVectors
                    map!(x -> permutedphiDict[currentTimeStep] * x, constraintProjVectors)
                    map!(x -> permutedphiDict[currentTimeStep] * x, guardProjVectors)
                    map!(x -> permutedphiDict[currentTimeStep] * x, invarientProjVectors)
                    mul!(tempM, Φ, PhiDict[currentTimeStep])
                    copy!(Φ, tempM)
                else
                    #newR = copy(newR)
                    currentTimeStep /= 2
                    #@show currentTimeStep, i
                    attempts += 1
                    #changedTimeStep = true
                end
            else
                #@show time
                return (time, Φ)
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
    #intersectingSet = overapproximate(ConvexHullArray(overapproximateIntersectingSetArray), Zonotope)
    #println(lastNewR)
    return (time, Φ)
end

#=
function ReACTGuardsBrake(loc, δ⁻::Float64, δ⁺::Float64, interval, guards, constraint, dirs, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, saveResult)
    # We calculate the reachset till we reach a guard for an intersection (or till failure)
    init1 = TIMEFUNC ? time_ns() : 0.0


    initialTimeStep = copy(δ⁺)
    #m = copy(δ⁻)
    #changedTimeStep = true
    #phiDict = PhiDict
    #A = copy(loc.A)
    # constraintProjVectors = map(x -> x.a, constraint)
    # constraintProjBounds = ρ.(constraintProjVectors, constraint)
    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    guardProjVectors, guardProjBounds = [], []
    noGuards = isnothing(guards)
    if noGuards
        guardProjVectors, guardProjBounds = [constraintProjVectors[1]], [Inf]
    else
        guardProjVectors, guardProjBounds = getHalfSpaceProjections(guards)
    end
    invarientProjVectors, invarientProjBounds = [], []
    if isnothing(loc.invarient)
        invarientProjVectors, invarientProjBounds = [constraintProjVectors[1]], [Inf]
    else
        invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)
    end
    dirProjVectors = map(x -> x, dirs)

    oldDirProjVectors = copy(dirProjVectors)
    #oldConstraintProjVectors = copy(constraintProjVectors)
    #oldGuardProjVectors = copy(guardProjVectors)
    #oldInvarientProjVectors = copy(invarientProjVectors)

    permutedphiDict = Dict()
    for key in keys(PhiDict)
        permutedphiDict[key] = permutedims(PhiDict[key])
    end


    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    #overapproximateIntersectingSetArray = []
    #lastNewR = []
    attemptsRecorder = []

    dirVals = [] # This is where the plotting happens


    Φ::Matrix{Float64} = diagm(ones(Float64, size(loc.A, 2)))

    tempM = diagm(ones(Float64, size(loc.A, 2)))
    #ϕt = similar(tempM)
    #pϕt = similar(tempM)
    #V = copy(inputDiscritezationDict[initialTimeStep])

    #lastVs = copy(Vs)
    Sρ = zeros(Float64, length(constraint))
    Gρ = zeros(Float64, length(guardProjVectors))
    Iρ = zeros(Float64, length(invarientProjVectors))
    dρ = zeros(Float64, length(dirs))
    #newR = discritezationDict[initialTimeStep]
    i = 1

    #U = inputDiscritezationDict[0]
    if TIMEFUNC
        global totalInitGuardTime += time_ns() - init1
    end
    while time < endtime
        #println("Time iss: $time")
        attempts = 1
        approveFlag = false
        while !approveFlag
            #println("Stuck?")

            # Handle if we can no longer reduce the reachset (we keep hitting something)
            if currentTimeStep >= δ⁻
                # If we hit a constraint
                #newRR = concretize(newRR)
                #if all((input + ρ(x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds)) && any((input + ρ(x, newR)) > y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))
                #    handleHitConstraint(time, loc.id)
                #end
                #return (dirVals, time, missing)


                #if changedTimeStep
                #newR = discritezationDict[currentTimeStep]

                #V = copy(inputDiscritezationDict[currentTimeStep])
                #ϕt = phiDict[currentTimeStep]
                #pϕt = permutedphiDict[currentTimeStep]

                #end



                # #println(any((input + ρ(-x, newRR)) > y for (input, x, y) in zip(Sρ, invarientProjVectors, invarientProjBounds)))
                #println(all((input + ρ(x, newRR)) <= y for (input, x, y) in zip(Sρ, invarientProjVectors, invarientProjBounds)))
                #=
                if isnothing(guards) || (all((input + ρ(x, discritezationDict[currentTimeStep])) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) &&
                                     #any(sign(y) >= 0 ? (input + ρ(-x, newR)) <= y : !((input + ρ(x, newR)) < y) for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) &&
                                     any((input + -ρ(-x, discritezationDict[currentTimeStep])) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) &&
                                     all((input - ρ(-x, discritezationDict[currentTimeStep])) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds)))
                                        =#
                if isnothing(guards) || guardCheck(discritezationDict[currentTimeStep], Sρ, constraintProjVectors, constraintProjBounds, Gρ, guardProjVectors, guardProjBounds, Iρ, invarientProjVectors, invarientProjBounds)
                    #(any((input + ρ(-x, newRR)) > y for (input, x, y) in zip(-Sρ, invarientProjVectors, invarientProjBounds)) && all((input + ρ(x, newRR)) <= y for (input, x, y) in zip(Sρ, invarientProjVectors, invarientProjBounds)))
                    #all(((ρ(x, tempSet)) <= y) for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Subset
                    #(isnothing(loc.invarient) || intersects(tempSet, loc.invarient))
                    #if mapreduce(x -> intersects(newRR, x), &, guard)
                    #push!(overapproximateIntersectingSetArray, newRR)
                    if saveResult
                        push!(dirVals, (copy(dρ + map(x -> ρ(x, discritezationDict[currentTimeStep]), oldDirProjVectors)), [time, time + currentTimeStep]))
                        dρ += map(x -> ρ(x, inputDiscritezationDict[currentTimeStep]), oldDirProjVectors)
                        dirProjVectors = map(x -> permutedphiDict[currentTimeStep] * x, oldDirProjVectors)
                        oldDirProjVectors = dirProjVectors
                    end
                    #lastVs = copy(Vs)
                    #Vs = ReachabilityAnalysis.Exponentiation.Φ₁(A, time - minimum(interval), ReachabilityAnalysis.Exponentiation.BaseExp) * U

                    approveFlag = true
                    Sρ += map(x -> ρ(x, inputDiscritezationDict[currentTimeStep]), constraintProjVectors)
                    Gρ += map(x -> ρ(x, inputDiscritezationDict[currentTimeStep]), guardProjVectors)
                    Iρ += map(x -> ρ(x, inputDiscritezationDict[currentTimeStep]), invarientProjVectors)
                    #constraintProjVectors = map(x -> permutedphiDict[currentTimeStep] * x, oldConstraintProjVectors)
                    #guardProjVectors = map(x -> permutedphiDict[currentTimeStep] * x, oldGuardProjVectors)
                    #invarientProjVectors = map(x -> permutedphiDict[currentTimeStep] * x, oldInvarientProjVectors)

                    #oldConstraintProjVectors = constraintProjVectors
                    #oldGuardProjVectors = guardProjVectors
                    #oldInvarientProjVectors = invarientProjVectors
                    map!(x -> permutedphiDict[currentTimeStep] * x, constraintProjVectors)
                    map!(x -> permutedphiDict[currentTimeStep] * x, guardProjVectors)
                    map!(x -> permutedphiDict[currentTimeStep] * x, invarientProjVectors)
                    mul!(tempM, Φ, PhiDict[currentTimeStep])
                    copy!(Φ, tempM)
                else
                    #newR = copy(newR)
                    currentTimeStep = currentTimeStep / 2
                    #@show currentTimeStep, i
                    attempts = attempts + 1
                    #changedTimeStep = true
                end
            else
                #@show time
                return (dirVals, time, Φ)
            end
        end

        push!(attemptsRecorder, attempts)
        i = i + 1
        time = time + currentTimeStep
        # Reset / apply strategy
        # Only do this if the current timestep is less than the initial
        if STRATEGY == 0
            # Only reduce
        elseif STRATEGY == 1
            # always try double
            if currentTimeStep < initialTimeStep
                currentTimeStep = currentTimeStep * 2
                #changedTimeStep = true
            end
        elseif STRATEGY == 2
            # If attemptsrecorder past 4 are successes, double timestep
            if currentTimeStep < initialTimeStep
                lowest = min(4, i - 1)
                window = @view attemptsRecorder[i-lowest:i-1]
                if all(window .== 1)
                    currentTimeStep = currentTimeStep * 2
                    #changedTimeStep = true
                end

            end
        end
    end
    #intersectingSet = overapproximate(ConvexHullArray(overapproximateIntersectingSetArray), Zonotope)
    #println(lastNewR)
    return (dirVals, time, Φ)
end
=#
function ReACT(loc, δ⁻::Float64, δ⁺::Float64, interval, constraint, STRATEGY::Integer, permutedphiDict, discritezationDict, inputDiscritezationDict, projVectors, projBounds, projRanges)
    # We calculate the reachset till we reach a guard for an intersection (or till failure)
    #println("HUH")

    init1 = TIMEFUNC ? time_ns() : 0.0


    # constraintProjVectors = map(x -> x.a, constraint)
    # constraintProjBounds = ρ.(constraintProjVectors, constraint)
    constraintProjVectors = view(projVectors, projRanges[1])
    constraintProjBounds = view(projBounds, projRanges[1])
    invarientProjVectors = view(projVectors, projRanges[2])
    invarientProjBounds = view(projBounds, projRanges[2])


    # projVectors = vcat(constraintProjVectors, constraintProjBounds)
    # constraintProjVectors = 


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
function guardCheck(newR::Zonotope, Sρ, constraintProjVectors, constraintProjBounds, Gρ, guardProjVectors, guardProjBounds, Iρ, invarientProjVectors, invarientProjBounds) #; solver=model
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
    return allunder(newR, Sρ, constraintProjVectors, constraintProjBounds) && someoutside(newR, Gρ, guardProjVectors, guardProjBounds) && someunder(newR, Iρ, invarientProjVectors, invarientProjBounds)#all((input - ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))

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
    res = true
    @inbounds @simd for i in eachindex(offset, dir, bound)
        res = res && (offset[i] - ρ(-dir[i], set) <= bound[i])
    end
    return res
end

function someoutside(set, offset, dir, bound)
    res = false
    @inbounds @simd for i in eachindex(offset, dir, bound)
        res = res || (offset[i] - ρ(-dir[i], set) > bound[i])
    end
    return res
end

function allunder(set, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (ρ(dir[i], set) <= bound[i])
    end
    return res
end

function someunder(set, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-dir[i], set) <= bound[i])
    end
    return res
end

function someoutside(set, dir, bound)
    res = false
    @inbounds @simd for i in eachindex(dir, bound)
        res = res || (-ρ(-dir[i], set) > bound[i])
    end
    return res
end


#=
function ReACTTouches2(loc, δ⁻::Float64, δ⁺::Float64, interval, guard, constraint, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, Φ, accInput)
    initialTimeStep = copy(δ⁺)
    m = copy(δ⁻)
    changedTimeStep = true
    phiDict = PhiDict

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = m

    invariantIsMissing = ismissing(loc.invarient)


    V = copy(inputDiscritezationDict[initialTimeStep])
    Vs = concretize(ReachabilityAnalysis.Exponentiation.Φ₁(loc.A, time, ReachabilityAnalysis.Exponentiation.BaseExp) * inputDiscritezationDict[0])
    lastVs = copy(Vs)
    Sρ = zeros(Float64, length(constraint))
    newR = discritezationDict[initialTimeStep]
    i = 1


    if ismissing(Φ)
        Φ::Matrix{Float64} = exp((time) .* loc.A)
    end
    tempM = diagm(ones(Float64, size(loc.A, 2)))
    ϕt = similar(Φ)
    newRR = copy(newR)

    overapproximateIntersectingSetArray = []
    runningset = minkowski_sum(linear_map(Φ, discritezationDict[m]), accInput)  #minkowski_sum(linear_map(Φ * PhiDict[m], discritezationDict[m]), minkowski_sum(linear_map(PhiDict[m], inputDiscritezationDict[m]), accInput))
    preclustering = copy(runningset) #minkowski_sum(linear_map(Φ * PhiDict[m], discritezationDict[m]), minkowski_sum(linear_map(PhiDict[m], inputDiscritezationDict[m]), accInput))
    #time += m
    #println("First intersects?: ", intersects(runningset, guard), " ", isdisjoint(runningset, loc.invarient))
    #println(preclustering)
    attemptsRecorder = []


    V = copy(inputDiscritezationDict[m])
    ϕt = phiDict[m]


    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m

                if any((input + ρ(x, runningset)) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))
                    #if !reduce(&, <=(Sρ + map(x -> ρ(x, runningset), constraintProjVectors), constraintProjBounds))
                    #throw(ErrorException("Reached unsafe set."))
                    handleHitConstraint(time, loc.id)
                end
                #bigCH = foldr((x, y) -> overapproximate(CH(x, y), Zonotope), overapproximateIntersectingSetArray; init=concretize(newRR))
                #intersectingSet = overapproximate(bigCH, Zonotope)
                if ismissing(preclustering)
                    preclustering = minkowski_sum(newRR, Vs)
                end
                #println("Touches time: ", time, " ", i)
                #println(norm(intersectingSet), " ", norm(newRR))
                return (reduce_order(preclustering, 10), Sρ, time)
            end

            runningset = minkowski_sum(linear_map(ϕt, runningset), V)

            if !isnothing(loc.invarient) && isdisjoint(runningset, loc.invarient)
                #  println("Disjoint! ", i)
                return (reduce_order(preclustering, 10), Sρ, time)
            end

            if !isnothing(loc.invarient) && !isSubSet(runningset, loc.invarient)
                #println("intersection!")
                runningset = concretize(zonotopeStripIntersection(runningset, loc.invarient))
            end
            hom = map(x -> ρ(x, runningset), constraintProjVectors)
            flowpipe = map(x -> ρ(x, concretize(minkowski_sum(linear_map(Φ, discritezationDict[m]), concretize(ReachabilityAnalysis.Exponentiation.Φ₁(loc.A, time, ReachabilityAnalysis.Exponentiation.BaseExp) * inputDiscritezationDict[0])))), constraintProjVectors)

            if all(ρ(x, runningset) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) &&
               all(ρ(x, flowpipe) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) &&
               intersects(runningset, guard) &&
               (invariantIsMissing || isSubSet(runningset, loc.invarient)) #mapreduce(x -> intersects(newRR, x), &, guard)

                #if reduce(&, <=(hom, constraintProjBounds)) && reduce(&, <=(flowpipe, constraintProjBounds)) && intersects(runningset, guard) && (invariantIsMissing || isSubSet(runningset, loc.invarient)) #mapreduce(x -> intersects(newRR, x), &, guard)
                #if mapreduce(x -> intersects(newRR, x), &, guard)
                #push!(overapproximateIntersectingSetArray, concretize(minkowski_sum(newRR, Vs)))
                if ismissing(preclustering)
                    preclustering = runningset
                else
                    preclustering = overapproximate(ConvexHull!(preclustering, runningset), Zonotope)
                end

                approveFlag = true
                mul!(tempM, Φ, ϕt)
                copy!(Φ, tempM)
            else
                newR = copy(newR)
                currentTimeStep = currentTimeStep / 2
                changedTimeStep = true
                attempts = attempts + 1
            end
        end


        push!(attemptsRecorder, attempts)
        i = i + 1
        time = time + currentTimeStep
    end
    #bigCH = foldr((x, y) -> overapproximate(CH(x, y), Zonotope), overapproximateIntersectingSetArray; init=concretize(newRR))
    #intersectingSet = overapproximate(bigCH, Zonotope)
    return (preclustering, Sρ, time)
end
=#