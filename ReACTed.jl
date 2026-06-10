using LazySets, LinearAlgebra, ReachabilityAnalysis
include("Discretize.jl")
include("Utilities.jl")

const TIMEFUNCS = false
const VERBOSE = false


totalDiscTime = 0
totalGuardTime = 0
totalTouchesTime = 0

function zonotopePrintDim(Z::Zonotope, dim)
    center = Z.center[dim]
    G = genmat(Z)
    _, genAmount = size(G)
    genContribute = sum(abs(G[dim, i]) for i in 1:genAmount)
    return "center: $center, generators: $genContribute"
end

function ReACTed(hybridSystem::HybridSystemV2, initialLoc, interval, X0, U, dirs, constraint, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, clustering=true, timeConstraintList=[]) where {N}
    loc = initialLoc
    #flowPhiDict = Dict(map(x -> x.id => PhiInputDict(x, δ⁻, δ⁺, alg, maxOrder, reduceOrder), hybridSystem.locations))

    #phiDicts, inputDicts = Dict(map(x -> x.id => PhiInputDict(x, δ⁻, δ⁺, alg, maxOrder, reduceOrder), hybridSystem.locations)[1, :])
    phiDicts = Dict()
    inputDicts = Dict()
    for x in hybridSystem.locations
        pd, id = PhiInputDict(x, δ⁻, δ⁺, alg, maxOrder, reduceOrder)
        phiDicts[x.id] = pd
        inputDicts[x.id] = id
    end


    waitinglist = []
    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []
    saveResult = true
    dirsVectors = []
    zenoBound = 500
    transitionCount = 0
    dimLength = size(X0.center, 1)

    if isempty(dirs)
        saveResult = false
    else
        for dir in dirs
            template = zeros(dimLength)
            template[dir] = 1.

            push!(dirsVectors, template)
            push!(dirsVectors, -template)
        end
    end

    #res = auxReACTed(hybridSystem, hybridSystem.locations[loc], nothing, interval, X0, dirsVectors, constraint, δ⁻, δ⁺, flowPhiDict, alg, maxOrder, reduceOrder, missing, clustering, timeConstraintList, saveResult)

    #reachset = vcat(reachset, res)

    #return res
    push!(waitinglist, (loc, X0, nothing, interval, diagm(ones(dimLength))))

    if TIMEFUNCS
        recorderTimeStart = time_ns()
    end

    while !isempty(waitinglist)
        #GC.gc()

        location, initialset, edge, interval′, TΦ = pop!(waitinglist)

        #reducedPolySet = polySet #isnothing(polySet) ? polySet : reducePolytope(polySet)

        #@show LazySets.API.high(initialset)
        #@show LazySets.API.low(initialset)

        #=if !isnothing(polySet)
            @show initialset ⊆ polySet
            @show LazySets.API.high(polySet) - LazySets.API.high(initialset)
            @show LazySets.API.low(polySet) - LazySets.API.low(initialset)
        end=#
        #phiDict, inputDict = flowPhiDict[location]
        #res = auxReACTedAlternative(waitinglist, hybridSystem, hybridSystem.locations[location], edge, interval′, initialset, dirsVectors, constraint, δ⁻, δ⁺, phiDict, inputDict, alg, maxOrder, reduceOrder, TΦ, reducedPolySet, timeConstraintList, saveResult; clustering, mustSemantics)
        res = auxReACTed(waitinglist, hybridSystem, hybridSystem.locations[location], edge, interval′, initialset, dirsVectors, constraint, δ⁻, δ⁺, phiDicts[location], inputDicts[location], alg, maxOrder, reduceOrder, TΦ, clustering, timeConstraintList, saveResult)

        reachset = vcat(reachset, res)
        if transitionCount < zenoBound
            transitionCount += 1
        else
            return reachset
        end
    end
    
    if TIMEFUNCS
        println("Total times: 
        TotalAux: $((time_ns() - recorderTimeStart - totalDiscTime - totalGuardTime - totalTouchesTime) / 10^9)
        Disc: $(totalDiscTime / 10^9)
        Guards: $(totalGuardTime / 10^9) 
        Touches: $(totalTouchesTime / 10^9)")
    end

    return reachset
end


# AucReacted is called recursively each time we have a new starting location (after a transition)
function auxReACTed(waitinglist, hybridSystem, loc::Location, currentEdge, interval, X0, dirs, constraint, δ⁻::Float64, δ⁺::Float64, PhiDict, inputDiscritezationDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, Φ=missing, clustering=true, timeConstraintList=[], saveResult::Bool=true) where {N}
    activeTimeConstraints = []
    for (id, time) in timeConstraintList
        if id == loc.id
            push!(activeTimeConstraints, time)
        end
    end
    recorderTimeStart = 0 # This is used for timing functions

    dims = size(X0.center, 1)
    if TIMEFUNCS
        recorderTimeStart = time_ns()
    end

    #PhiDict, inputDiscritezationDict = FlowPhiDict[loc.id]

    #t1 = Base.time()
    discretizationDict = newReACTDiscretizePlus(loc, X0, δ⁻, δ⁺, PhiDict, inputDiscritezationDict, alg, maxOrder, reduceOrder)
    #t2 = Base.time()
    #@show t2 - t1
    if TIMEFUNCS
        global totalDiscTime += time_ns() - recorderTimeStart
    end


    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []
    setOfConstraints = vcat(loc.constraints, constraint)


    # For each edge we simulate the system
    listOfEdges = loc.edges
    if !isa(currentEdge, Nothing) # This allows us to sometimes just wanna do a signular edge
        listOfEdges = [currentEdge]
    end

    for edge in listOfEdges
        # println("Handling edge at time $time from $(loc.id) -> $(edge.targetLoc)")
        guards = edge.guard # Technically the guard is one singular HPolyhedron, but it composes the other guards


        #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 

        if TIMEFUNCS
            recorderTimeStart = time_ns()
        end
        tempReachset, reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, [time, endtime], guards, setOfConstraints, dirs, 2, PhiDict, discretizationDict, inputDiscritezationDict, saveResult)
        if TIMEFUNCS
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

            if saveResult
                push!(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(edge.targetLoc)))
            end
            # println("Current time: $time, intersecting time start: $reachtime")
            if TIMEFUNCS
                recorderTimeStart = time_ns()
            end
            tryContinueFlag, _, timeNotIntersected, intersectingSetsList = ReACTTouches(loc, δ⁻, [reachtime, endtime], (reachtime - time), guards, setOfConstraints, 2, PhiDict, discretizationDict, inputDiscritezationDict, tΦ, nothing, reduceOrder, maxOrder)
            if TIMEFUNCS
                global totalTouchesTime += time_ns() - recorderTimeStart
            end
            if any((timeNotIntersected >= x) for x in activeTimeConstraints)
                handleHitConstraint(timeNotIntersected, loc.id)
            end


            timeIntersected = timeNotIntersected - reachtime
            latestSet = nothing
            #
            #   Here we should check whether we have reached endtime. If true we should only push the jumpSet
            #   Still need to check whether we have reached the invariant. If true we should NOT push the else branch result, only the tempReachsets[infMaxsIdx]
            #
            if !isempty(intersectingSetsList)
                # we first save this for later
                if tryContinueFlag
                    latestSet = copy(last(intersectingSetsList))
                end
                # Afterwards we process the intersecting set

                intersectedSet = nothing
                plottingList = []
                jumpSetsList = Vector{Zonotope}()
                i = 0

                for set in intersectingSetsList
                    jumpSet = copy(set)

                    if !isa(loc.invarient, Nothing)
                        jumpSet = zonotopeStripIntersection(jumpSet, loc.invarient)

                    end

                    any(x -> ρ(x.a, jumpSet) > x.b, setOfConstraints) && handleHitConstraint(reachtime + δ⁻ * i, loc.id)

                    if !isa(guards, Nothing)
                        jumpSet = zonotopeStripIntersection(jumpSet, guards)
                    end

                    # Check invarient
                    # Push to reachset
                    if saveResult
                        push!(plottingList, (map(x -> ρ(x, jumpSet), dirs), [reachtime + δ⁻ * i, reachtime + δ⁻ * (i + 1)]))
                        i += 1
                    end

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


                if saveResult
                    push!(reachset, (plottingList, string(reachtime) * " - " * string(timeNotIntersected) * ": " * string(loc.id) * "->" * string(loc.id)))
                end

                len = length(jumpSetsList)
                if len != 0
                    if clustering
                        intersectedSet = nothing
                        if len == 1
                            intersectedSet = jumpSetsList[1]
                        else
                            tempIntersect = foldl(ConvexHull, jumpSetsList)
                            intersectedSet = convert(Zonotope, box_approximation(tempIntersect))
                        end

                        jumpSet = intersectedSet
                        # println("Finished intersections")
                        push!(waitinglist, (edge.targetLoc, jumpSet, nothing, [reachtime, endtime], tΦ))

                        #branchedRun = auxReACTed(hybridSystem, hybridSystem.locations[edge.targetLoc], nothing, [reachtime, endtime], jumpSet, dirs, constraint, δ⁻, δ⁺, FlowPhiDict, alg, maxOrder, reduceOrder, tΦ, clustering, timeConstraintList, saveResult)

                        #if saveResult
                        #    reachset = vcat(reachset, branchedRun)
                        #end
                    else # No clustering. Do individual runningset
                        for (index, jumpSet) in enumerate(jumpSetsList)
                            # Note that we only do steps of size δ⁻ in touches
                            #timeStart = reachtime + (index - 1) * δ⁻
                            push!(waitinglist, (edge.targetLoc, jumpSet, nothing, [reachtime + (index - 1) * δ⁻, endtime], nothing)) # gav tϕ foer..
                            #branchedRun = auxReACTed(hybridSystem, hybridSystem.locations[edge.targetLoc], nothing, [timeStart, endtime], jumpSet, dirs, constraint, δ⁻, δ⁺, FlowPhiDict, alg, maxOrder, reduceOrder, tΦ, clustering, timeConstraintList, saveResult)

                            #if saveResult
                            #    reachset = vcat(reachset, branchedRun)
                            #end
                        end
                    end
                else
                    continue
                end
            else
                tryContinueFlag = false
                # println("No intersections with guards, end branch")

            end

            # If we are not encountering an invarient, try continue
            if tryContinueFlag
                #println("Continuing from previous run at time $timeNotIntersected $reachtime $(length(intersectingSetsList))")

                latestSet = zonotopeStripIntersection(latestSet, loc.invarient)
                push!(waitinglist, (edge.targetLoc, latestSet, edge, [timeNotIntersected, endtime], tΦ)) # gav tϕ foer..
                #branchedRun = auxReACTed(hybridSystem, loc, edge, [timeNotIntersected, endtime], latestSet, dirs, constraint, δ⁻, δ⁺, FlowPhiDict, alg, maxOrder, reduceOrder, tΦ, clustering, timeConstraintList, saveResult)
                #if saveResult
                #    reachset = vcat(reachset, branchedRun)
                #end
            end


        else # Reached the end time
            # println("Reached end time before intersecting guard")
            # Check if we hit activeTimeConstraints
            if any((reachtime >= x) for x in activeTimeConstraints)
                handleHitConstraint(reachtime, loc.id)
            end
            if saveResult
                if !isa(loc.invarient, Nothing)
                    newReach = []
                    for (Z, timeInterval) in tempReachset

                        push!(newReach, (Z, timeInterval))
                    end
                    tempReachset = newReach
                    reachset = vcat(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(edge.targetLoc)))
                else
                    reachset = vcat(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(edge.targetLoc)))
                end
            end
        end
    end

    if isempty(loc.edges) # This means it is just a continous system from here
        # println("No edges? call ReACT")

        #tempReachset, tempInput, reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, [time, endtime], guards, setOfConstraints, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, missing, saveResult)


        #tempReachset, _, _ = ReACT(loc, δ⁻, δ⁺, interval, constraint, 2, alg, maxOrder, reduceOrder, PhiDict[loc.id])
        tempReachset, reachtime, _ = ReACT(loc, δ⁻, δ⁺, [time, endtime], setOfConstraints, dirs, 2, PhiDict, discretizationDict, inputDiscritezationDict, saveResult)

        if saveResult
            reachset = vcat(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(loc.id)))
        end
    end

    return reachset
end

function ReACTTouches(loc, δ⁻::Float64, interval, initialTime::Float64, guard, constraint, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, Φ, accInput, reduce_order, max_order)
    # Note that in touches we always use δ⁻
    # That is, we do not adjust timestep sizes

    phiDict = PhiDict

    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(guard)
    invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)


    Vs = nestedInputDiscCalculate(inputDiscritezationDict, PhiDict, δ⁻, initialTime, reduce_order, max_order)
    # i = 1

    if ismissing(Φ)
        Φ::Matrix{Float64} = exp(initialTime .* loc.A)
    end

    ϕ = similar(Φ)
    ϕ = phiDict[δ⁻]

    # Compute current sets
    newR = discritezationDict[δ⁻]
    newRR = linear_map(Φ, newR)
    V = copy(inputDiscritezationDict[δ⁻])
    V = linear_map(Φ, V)
    intersectingSetsList = Vector{Zonotope}() # This will store all of our sets

    while time < endtime
        tempSet = minkowski_sum(newRR, Vs)

        if all((ρ(x, tempSet)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) && # IsSubSet
           all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds)) &&
           all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Intersects
            #if mapreduce(x -> intersects(newRR, x), &, guard)
            #@show time
            push!(intersectingSetsList, copy(tempSet))

            # Main calculation. No longer changing timesteps
            Vs = minkowski_sum(Vs, V) #Update input
            newRR = linear_map(ϕ, newRR)
            V = linear_map(ϕ, V)

            # i = i + 1
            time = time + δ⁻

        else # Handle hit something. We do not reduce anymore!
            #@show all((ρ(x, tempSet)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds))
            #@show all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds))
            #@show all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds))
            #if all((ρ(x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) && any(((ρ(x, tempSet)) > y) for (x, y) in zip(constraintProjVectors, constraintProjBounds))
            #    handleHitConstraint(time, loc.id)
            #end

            continueAfter = true


            if all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds)) || !all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds))
                # If we stop because we are no longer intersect guards, 
                # but still intersect the invariant we try continue

                #@show all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds))
                continueAfter = false
            end

            if continueAfter && all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Intersects
                intersectionSet = LazySets.Intersection(tempSet, loc.invarient)
                if !all((ρ(x, intersectionSet)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds))  # IsSubSet
                    @show LazySets.API.high(tempSet)
                    @show LazySets.API.low(tempSet)
                    @show LazySets.API.high(Vs)
                    @show LazySets.API.low(Vs)
                    handleHitConstraint(time, loc.id)
                end
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

            return continueAfter, Φ, time, intersectingSetsList
        end
    end
    continueAfter = false # We are at the end time horizon, therefore no continuing
    return (continueAfter, Φ, time, intersectingSetsList)
end


function ReACTGuards(loc, δ⁻::Float64, δ⁺::Float64, interval, guards, constraint, dirs, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, saveResult)
    # We calculate the reachset till we reach a guard for an intersection (or till failure)
    initialTimeStep = copy(δ⁺)
    #m = copy(δ⁻)
    #changedTimeStep = true
    #phiDict = PhiDict
    #A = copy(loc.A)
    # constraintProjVectors = map(x -> x.a, constraint)
    # constraintProjBounds = ρ.(constraintProjVectors, constraint)
    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    guardProjVectors, guardProjBounds = [], []
    if isnothing(guards)
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
    oldConstraintProjVectors = copy(constraintProjVectors)
    oldGuardProjVectors = copy(guardProjVectors)
    oldInvarientProjVectors = copy(invarientProjVectors)

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


    Φ::Matrix{Float64} = exp(0 .* loc.A)

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
                if isnothing(guards) || (all((input + ρ(x, discritezationDict[currentTimeStep])) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) &&
                                         #any(sign(y) >= 0 ? (input + ρ(-x, newR)) <= y : !((input + ρ(x, newR)) < y) for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) &&
                                         any((input + -ρ(-x, discritezationDict[currentTimeStep])) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) &&
                                         all((input - ρ(-x, discritezationDict[currentTimeStep])) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds)))
                    #(any((input + ρ(-x, newRR)) > y for (input, x, y) in zip(-Sρ, invarientProjVectors, invarientProjBounds)) && all((input + ρ(x, newRR)) <= y for (input, x, y) in zip(Sρ, invarientProjVectors, invarientProjBounds)))
                    #all(((ρ(x, tempSet)) <= y) for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Subset
                    #(isnothing(loc.invarient) || intersects(tempSet, loc.invarient))
                    #if mapreduce(x -> intersects(newRR, x), &, guard)
                    #push!(overapproximateIntersectingSetArray, newRR)
                    if saveResult
                        push!(dirVals, (copy(dρ + map(x -> ρ(x, discritezationDict[currentTimeStep]), oldDirProjVectors)), [time, time + currentTimeStep]))
                    end
                    #lastVs = copy(Vs)
                    #Vs = ReachabilityAnalysis.Exponentiation.Φ₁(A, time - minimum(interval), ReachabilityAnalysis.Exponentiation.BaseExp) * U

                    approveFlag = true
                    dρ += map(x -> ρ(x, inputDiscritezationDict[currentTimeStep]), oldDirProjVectors)
                    Sρ += map(x -> ρ(x, inputDiscritezationDict[currentTimeStep]), constraintProjVectors)
                    Gρ += map(x -> ρ(x, inputDiscritezationDict[currentTimeStep]), guardProjVectors)
                    Iρ += map(x -> ρ(x, inputDiscritezationDict[currentTimeStep]), invarientProjVectors)
                    dirProjVectors = map(x -> permutedphiDict[currentTimeStep] * x, oldDirProjVectors)
                    constraintProjVectors = map(x -> permutedphiDict[currentTimeStep] * x, oldConstraintProjVectors)
                    guardProjVectors = map(x -> permutedphiDict[currentTimeStep] * x, oldGuardProjVectors)
                    invarientProjVectors = map(x -> permutedphiDict[currentTimeStep] * x, oldInvarientProjVectors)

                    oldDirProjVectors = dirProjVectors
                    oldConstraintProjVectors = constraintProjVectors
                    oldGuardProjVectors = guardProjVectors
                    oldInvarientProjVectors = invarientProjVectors
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

function ReACT(loc, δ⁻::Float64, δ⁺::Float64, interval, constraint, dirs, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, saveResult)
    # We calculate the reachset till we reach a guard for an intersection (or till failure)
    #println("HUH")

    initialTimeStep = copy(δ⁺)
    m = copy(δ⁻)
    changedTimeStep = true
    phiDict = PhiDict
    A = copy(loc.A)
    # constraintProjVectors = map(x -> x.a, constraint)
    # constraintProjBounds = ρ.(constraintProjVectors, constraint)
    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)
    dirProjVectors = map(x -> x, dirs)

    oldDirProjVectors = copy(dirProjVectors)
    oldConstraintProjVectors = copy(constraintProjVectors)
    oldInvarientProjVectors = copy(invarientProjVectors)

    permutedphiDict = Dict()
    for key in keys(phiDict)
        permutedphiDict[key] = permutedims(phiDict[key])
    end


    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(δ⁻)#copy(initialTimeStep)

    attemptsRecorder = []

    dirVals = [] # This is where the plotting happens


    Φ::Matrix{Float64} = exp(0 .* loc.A)

    tempM = diagm(ones(Float64, size(loc.A, 2)))
    ϕt = similar(tempM)
    pϕt = similar(tempM)
    V = copy(inputDiscritezationDict[initialTimeStep])

    #lastVs = copy(Vs)
    Sρ = zeros(Float64, length(constraint))
    Iρ = zeros(Float64, length(invarientProjVectors))
    dρ = zeros(Float64, length(dirs))
    newR = discritezationDict[initialTimeStep]
    i = 1

    U = inputDiscritezationDict[0]
    while time < endtime
        # println("Time iss: $time")
        attempts = 1
        approveFlag = false

        while !approveFlag

            # Handle if we can no longer reduce the reachset (we keep hitting something)
            if currentTimeStep < m
                # If we hit a constraint
                if any((input + ρ(x, newR)) > y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))

                    handleHitConstraint(time, loc.id)
                end
                return (dirVals, time, Φ)
            end

            if changedTimeStep
                newR = discritezationDict[currentTimeStep]

                V = copy(inputDiscritezationDict[currentTimeStep])
                ϕt = phiDict[currentTimeStep]
                pϕt = permutedphiDict[currentTimeStep]
            end


            changedTimeStep = false

            if all((input + ρ(x, newR)) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) &&
               all((input + -ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))

                if saveResult
                    push!(dirVals, (copy(dρ + map(x -> ρ(x, newR), oldDirProjVectors)), [time, time + currentTimeStep]))
                end

                approveFlag =
                    true
                dρ += map(x -> ρ(x, V), oldDirProjVectors)
                Sρ += map(x -> ρ(x, V), constraintProjVectors)
                Iρ += map(x -> ρ(x, V), invarientProjVectors)
                dirProjVectors = map(x -> pϕt * x, oldDirProjVectors)
                constraintProjVectors = map(x -> pϕt * x, oldConstraintProjVectors)
                invarientProjVectors = map(x -> pϕt * x, oldInvarientProjVectors)

                oldDirProjVectors = dirProjVectors
                oldConstraintProjVectors = constraintProjVectors
                oldInvarientProjVectors = invarientProjVectors
                mul!(tempM, Φ, ϕt)
                copy!(Φ, tempM)
            else
                #newR = copy(newR)
                currentTimeStep = currentTimeStep / 2
                changedTimeStep = true
                attempts = attempts + 1
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
                changedTimeStep = true
            end
        elseif STRATEGY == 2
            # If attemptsrecorder past 4 are successes, double timestep
            if currentTimeStep < initialTimeStep
                lowest = min(4, i - 1)
                window = @view attemptsRecorder[i-lowest:i-1]
                if all(window .== 1)
                    currentTimeStep = currentTimeStep * 2
                    changedTimeStep = true
                end
            end
        end
    end

    return (dirVals, time, Φ)
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