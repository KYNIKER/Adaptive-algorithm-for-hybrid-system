using LazySets, LinearAlgebra
include("Discretize.jl")
include("Utilities.jl")


function ReACTed(hybridSystem::HybridSystemV2, initialLoc, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, saveResult::Bool=true) where {N}
    loc = initialLoc
    flowPhiDict = Dict(map(x -> x.id => PhiDict(x.A, δ⁻, δ⁺, alg), hybridSystem.locations))




    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []

    res = auxReACTed(hybridSystem, hybridSystem.locations[loc], interval, X0, constraint, δ⁻, δ⁺, flowPhiDict, alg, maxOrder, reduceOrder, missing, saveResult)

    reachset = vcat(reachset, res)

    return res
end

function auxReACTed(hybridSystem, loc::Location, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, constraint, δ⁻::Float64, δ⁺::Float64, PhiDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, Φ=missing, saveResult::Bool = true) where {N}
    discretizationDict, inputDiscritezationDict = ReACTDiscretize(loc, X0, δ⁻, δ⁺, alg, maxOrder, reduceOrder, PhiDict[loc.id])
    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []
    setOfConstraints = vcat(loc.constraints, constraint)

    println("Smallest disc center: ", discretizationDict[δ⁻].center)


    # guards = edge.guard # Technically the guard is one singular HPolyhedron, but it composes the other guards

    #guards = constraints_list(guards)

    #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 

    tempReachset, tempInput, reachtime, tΦ, guardIntersectTattler = ReACTTattler(loc, δ⁻, δ⁺, [time, endtime], setOfConstraints, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, missing, saveResult)

    if saveResult
        listOfChosenEdges = ""
        for (count, setOfInterects) in enumerate(guardIntersectTattler)
            if !isempty(setOfInterects)
                listOfChosenEdges *= "$(loc.edges[count].targetLoc), "
            end
        end
        println("Tattling to: ", listOfChosenEdges)

        reachset = vcat(reachset, (tempReachset, string(loc.id) * "->" * " (" * listOfChosenEdges * ")"))
    end
    

    for (edgeCount, edge) in enumerate(loc.edges)

        for (intersectingSet, startTime) in guardIntersectTattler[edgeCount]

            #println(timeIntersected)
            if startTime >= 0.
                #timeIntersectedSet = concretize(overapproximateIntervalReachset(loc.A, tempReachsets[supMinsIdx], U, δ⁻, timeIntersected, alg, maxOrder, reduceOrder, flowPhiDict[loc.id]))
                #if !intersects(intersectingSet, guards)
                #throw(ErrorException("Set intersecting guard does not intersect the guard."))
                #end
                intersectedSet = intersectingSet



                if !isa(loc.invarient, Nothing) && intersects(intersectedSet, loc.invarient)
                    #println("tes")
                    tintersectedSet = getBoxIntersection(intersectedSet, loc.invarient)
                    #println("Inv intersection: ", ρ(Vector(sparsevec([5], [1.0], 6)), intersectedSet), " vs ", ρ(Vector(sparsevec([5], [1.0], 6)), tintersectedSet))
                    intersectedSet = tintersectedSet
                    #println(intersectedSet)
                end
                if isempty(intersectedSet)
                    println("Empty...")
                    continue # Go to next intersecting set
                end
                println("Inv intersectedSet: ", intersectedSet)

                # Surely this never happens right? Cause it is already the set that interdescts the guards?
                intersectedSet = getBoxIntersection(intersectedSet, edge.guard)
                if isempty(intersectedSet)
                    println("Empty..")
                    throw(error("This should not happen. The set intersecting guards is no longer intersecting guards!"))
                    return reachset
                end
                println("Guard intersectedSet: ", intersectedSet)

                #tempSet = exp(reachtime .* loc.A) * X0
                #=x, _ = tempReachset[end]
                #x = concretize(exp(reachtime .* loc.A) * X0)
                #x = getBoxIntersection(x, loc.invarient)
                if intersects(x, guards)
                    x = getBoxIntersection(x, loc.invarient)
                    if intersects(x, guards)
                        #throw(ErrorException("WHAT THE ACTUAL FUCK!!!!!"))
                    end
                    println("WHAT THE FUCK!")
                end

                if intersects(x, guards)
                    println("WHAT THE ACTUAL FUCK!")
                    throw(ErrorException("WHAT THE ACTUAL FUCK!"))
                end
                =#
                jumpSet = concretize(linear_map(edge.jumpMatrix, intersectedSet))
                #println(edge.jumpMatrix)

                if !isa(hybridSystem.locations[edge.targetLoc].invarient, Nothing) && intersects(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)
                    tjumpSet = getBoxIntersection(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)# + edge.jumpVector * intersectedSet
                    #println("Jump intersection: ", ρ(Vector(sparsevec([5], [1.0], 6)), jumpSet), " vs ", ρ(Vector(sparsevec([5], [1.0], 6)), tjumpSet))
                    jumpSet = tjumpSet
                    #println(jumpSet)
                    #println("tes2")
                end
                #push!(reachset, ([(jumpSet, [reachtime, reachtime])], string(reachtime) * ": jump(" * string(loc.id) * ")->" * string(edge.targetLoc)))
                #
                #   Here we could optimize it such that in the case where guards ⊆ timeIntersectedSet we calculate both [supMins, endtime] and [infMaxs, endtime] with guards
                #   and otherwise [supMins, endtime] with hyperplane intersection with timeIntersectedSet and [infMaxs, endtime] with guards intersection
                #   Maybe look at how input should be handled... and if we can manipulate the constraints to account for the accumulated input
                #

                #timePointInput = U  #   NEEDS FIXING
                #println(x.center)
                println("Jumpset center: ", jumpSet.center)
                #y, _ = tempReachset[1]
                #println(y.center)
                branchedRun = auxReACTed(hybridSystem, hybridSystem.locations[edge.targetLoc], [startTime, endtime], jumpSet, constraint, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder, tΦ, saveResult)
                
                if saveResult
                    reachset = vcat(reachset, branchedRun)
                end
                #reachset = vcat(reachset, nonintersectedSet) #Maybe gets the universe..

                #branchedRun = auxReACTed(loc, [timeIntersected, endtime], nonintersectedSet, intersectedInput, constraint, δ⁻, δ⁺, flowPhiDict, alg, maxOrder, reduceOrder)
                #push!(reachset, branchedRun)

            elseif timeIntersected == 0.0
                println("Nah but")
                #   When is this the case? 
                #=timePointInput = U  #   NEEDS FIXING
                tempIntectingSet = concretize(intersectingSet)
                if intersects(tempIntectingSet, constraints_list(loc.invarient))
                    tempReachsetInInv = getBoxIntersection(concretize(tempReachset), constraints_list(loc.invarient))
                    tempIntectingSet = tempReachsetInInv
                end
                #intersectedSet, nonintersectedSet = splitZonotope(intersectingSet, constraints_list(hybridSystem.locations[edge.targetLoc].invarient))
                branchedRun = auxReACTed(hybridSystem, loc, [reachtime, endtime], tempIntectingSet, U, constraint, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder)
                reachset = vcat(reachset, branchedRun)
                #=
                if ρ(timeIntersectedSet.center, guard.a) >= guard.b

                =#
                =#
            end
        end
    end

    if isempty(loc.edges) # This means it is just a continous system from here
        println("No edges? call ReACT")

        #tempReachset, _, _ = ReACT(loc, δ⁻, δ⁺, interval, constraint, 2, alg, maxOrder, reduceOrder, PhiDict[loc.id])
        tempReachset, _, _, _ = ReACT(loc, δ⁻, δ⁺, interval, constraint, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, Φ, saveResult)
        
        if saveResult
            if !isa(loc.invarient, Nothing)
                newReach = []
                for (Z, timeInterval) in tempReachset
                    if intersects(Z, loc.invarient)
                        push!(newReach, (getBoxIntersection(concretize(Z), loc.invarient), timeInterval))
                    else
                        push!(newReach, (Z, timeInterval))
                    end
                end
                reachset = vcat(reachset, (newReach, string(loc.id) * "->" * string(loc.id)))
            else
                reachset = vcat(reachset, (tempReachset, string(loc.id) * "->" * string(loc.id)))
            end
        end
    end


    return reachset
end

function ReACTTouches(loc, δ⁻::Float64, δ⁺::Float64, interval, guard, constraint, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, Φ, accInput)
    initialTimeStep = copy(δ⁺)
    m = copy(δ⁻)
    changedTimeStep = true
    phiDict = PhiDict

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)


    V = copy(inputDiscritezationDict[initialTimeStep])
    Vs = accInput
    lastVs = copy(Vs)
    Sρ = zeros(Float64, length(constraint))
    newR = discritezationDict[initialTimeStep]
    i = 1


    if ismissing(Φ)
        Φ::Matrix{Float64} = exp(time .* loc.A)
    end
    tempM = diagm(ones(Float64, size(loc.A, 2)))
    ϕt = similar(Φ)
    newRR = copy(newR)

    overapproximateIntersectingSetArray = []
    preclustering = minkowski_sum(linear_map(Φ * PhiDict[m], discritezationDict[m]), minkowski_sum(linear_map(PhiDict[m], inputDiscritezationDict[m]), accInput))
    attemptsRecorder = []

    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m

                
                if !reduce(&, <=(Sρ + map(x -> ρ(x, newRR), constraintProjVectors), constraintProjBounds))
                    #throw(ErrorException("Reached unsafe set."))
                    handleHitConstraint(time, loc.id)
                end
                bigCH = foldr((x, y) -> overapproximate(CH(x, y), Zonotope), overapproximateIntersectingSetArray; init=concretize(newRR))
                intersectingSet = overapproximate(bigCH, Zonotope)
                if ismissing(preclustering)
                    preclustering = minkowski_sum(newRR, Vs)
                end
                println("Order of input: ", LazySets.order(accInput))
                println("Touches Vs: ", accInput)
                #println(norm(intersectingSet), " ", norm(newRR))
                return (preclustering, Sρ, time)
            end

            if changedTimeStep
                newR = discritezationDict[currentTimeStep]
                V = copy(inputDiscritezationDict[currentTimeStep])
                ϕt = phiDict[currentTimeStep]
                newRR = linear_map(Φ, newR)
                V = linear_map(Φ, V)
            else
                newRR = linear_map(ϕt, newRR)
                V = linear_map(ϕt, V)
            end

            changedTimeStep = false
            hom = map(x -> ρ(x, newRR), constraintProjVectors)
            inhom = map(x -> ρ(x, V), constraintProjVectors)
            tempVs = remove_redundant_generators(minkowski_sum(Vs, V))

            if reduce(&, <=(Sρ + hom + inhom, constraintProjBounds)) && intersects(concretize(minkowski_sum(newRR, Vs)), guard) && intersects(concretize(minkowski_sum(newRR, Vs)), loc.invarient) #mapreduce(x -> intersects(newRR, x), &, guard)
                #if mapreduce(x -> intersects(newRR, x), &, guard)
                push!(overapproximateIntersectingSetArray, concretize(minkowski_sum(newRR, Vs)))
                if ismissing(preclustering)
                    preclustering = concretize(minkowski_sum(newRR, Vs))
                else
                    preclustering = overapproximate(CH(preclustering, concretize(minkowski_sum(newRR, Vs))), Zonotope)
                end
                lastVs = copy(Vs)
                Vs = reduce_order(copy(tempVs), 5)
                approveFlag = true
                Sρ += inhom
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
    bigCH = foldr((x, y) -> overapproximate(CH(x, y), Zonotope), overapproximateIntersectingSetArray; init=concretize(newRR))
    intersectingSet = overapproximate(bigCH, Zonotope)
    return (preclustering, Sρ, time)
end

function ReACTTattler(loc, δ⁻::Float64, δ⁺::Float64, interval, constraint, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, Φ, saveResult)
    # This produces a normal reachset
    # It also returns a list of times which we start intersecting guards
    # This will be of the form # [(edge, [(intersectingReachSet, startTime), ...]), ...]

    # We later use this to find intersected spots  


    amountOfEdges = length(loc.edges)


    guardIntersectTattler = [[] for e in loc.edges] # Here we store the final reachsets for each edge
    intersectBuildSpace = Vector{Union{Nothing, Zonotope}}(nothing, amountOfEdges) # Here we build the reachsets and then add them to the tattler once we stop intersecting
    startIntersectTimeTracker = Vector{Float64}(undef, amountOfEdges) # This keeps track of times we started intersecting
    currentlyIntersectingGuard = [(e.guard, false) for e in loc.edges] # This switches whether we want the set to avoid intersecting guards, or keep intersecting
    
    initialTimeStep = δ⁺
    changedTimeStep = true
    phiDict = PhiDict
    # discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    # inputDiscritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)

    #discritezationDict, inputDiscritezationDict, phiDict = ReACTDiscretize(A, X0, U, δ⁻, δ⁺, alg, maxOrder, reduceOrder, phiDict)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    attemptsRecorder = Integer[]

    #V::Zonotope{N,Vector{N},Matrix{N}} = copy(inputDiscritezationDict[initialTimeStep])
    V = copy(inputDiscritezationDict[initialTimeStep])
    Sρ = zeros(Float64, length(constraint))
    #newR::Zonotope{N,Vector{N},Matrix{N}} = discritezationDict[initialTimeStep]
    newR = discritezationDict[initialTimeStep]
    i = 1


    #Φ::Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    if ismissing(Φ)
        Φ::Matrix{Float64} = exp(0 .* loc.A)
    end
    tempM = similar(Φ)
    ϕt = similar(Φ)
    newRR = copy(newR)

    reachSets = []

    sρ = zeros(Float64, length(constraint))
    Vs = Zonotope(zeros(size(loc.A, 2)), [zeros(size(loc.A, 2))])

    while time < endtime
        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < δ⁻  # We must terminate as we are now outside an invarient
                
                # Before terminating, we look to see if the last set hits any guards
                currentStepRepresentation = concretize(minkowski_sum(newRR, Vs))

                for (count, (guard, intersectingPastStep)) in enumerate(currentlyIntersectingGuard)
                    intersectingGuardNow = intersects(currentStepRepresentation, guard)

                    if intersectingGuardNow # If we intersect now we must add it
                        println("Intersecting on last step")
                        if isnothing(intersectBuildSpace[count]) # If first step in sequence
                            intersectBuildSpace[count] = copy(currentStepRepresentation)
                        else # If not first step, we convex hull them together to 1 set
                            intersectBuildSpace[count] = overapproximate(CH(intersectBuildSpace[count], currentStepRepresentation), Zonotope)
                        end
                        if !intersectingPastStep # Remember to also update time
                            startIntersectTimeTracker[count] = time
                        end
                        intersectingPastStep = true
                    end
                    if intersectingPastStep # If we are intersecting this must be reported to the tattler
                        push!(guardIntersectTattler[count], (intersectBuildSpace[count], startIntersectTimeTracker[count]))
                    end
                end

                if isempty(reachSets)
                    println("Empty reachsets at time: ", time)
                    println("Current representation: ", currentStepRepresentation)
                    println("Homogeneous: ", newRR)
                    println("Input: ", Vs)
                end



                return reachSets, Vs, time, Φ, guardIntersectTattler
            end

            if changedTimeStep
                newR = discritezationDict[currentTimeStep]
                V = copy(inputDiscritezationDict[currentTimeStep])
                ϕt = phiDict[currentTimeStep]
                newRR = linear_map(Φ, newR)
                V = linear_map(Φ, V)
            else
                newRR = linear_map(ϕt, newRR)
                V = linear_map(ϕt, V)
            end

            changedTimeStep = false

            tempVs = remove_redundant_generators(minkowski_sum(Vs, V)) # Wny do we need this?


            notHitConstraint = false # This will be filled out later, based on the invarient
            inhom = map(x -> ρ(x, V), constraintProjVectors)

            currentStepRepresentation = concretize(minkowski_sum(newRR, Vs))

            containedInIntersection = isSubSet(currentStepRepresentation, loc.invarient)
            if !isnothing(loc.invarient) && !isSubSet(currentStepRepresentation, loc.invarient)
                #println("We are not a subset of the invarient at time $time with timestep $currentTimeStep")
                if currentTimeStep == δ⁻
                    boxIntersection = concretize(getBoxIntersection(currentStepRepresentation, loc.invarient))
                    if !isempty(boxIntersection)
                        #println("Cut strategy applied! aka we still intersect")

                        # We intersect partly, therefore we adjust the step representation to be this intersection
                        # Note that the actual computation does not use this. But we use this to represent the reachable set
                        # It it also used for jump sets

                        currentStepRepresentation = boxIntersection
                        # We must calculate whether we hit constraints with the new representation
                        projectedCurrentStepRepresentation = map(x -> ρ(x, currentStepRepresentation), constraintProjVectors)
                        notHitConstraint = reduce(&, <=(Sρ + projectedCurrentStepRepresentation, constraintProjBounds))
                    else 
                        # Even at the lowest timestep, we are completely out of the invarient. Therefore end the continous step
                        newR = copy(newR)
                        currentTimeStep = currentTimeStep / 2
                        changedTimeStep = true
                        attempts = attempts + 1
                        continue
                    end
                else
                    # We are not at lowest timestep yet
                    newR = copy(newR)
                    currentTimeStep = currentTimeStep / 2
                    changedTimeStep = true
                    attempts = attempts + 1
                    continue
                end
            else
                # We can calculate the nothitConstraint as follows, using the hom and inhom solution
                # TODO figure out if we can use currentStepRepresentation instead of hom and inhom
                hom = map(x -> ρ(x, newRR), constraintProjVectors)

                notHitConstraint = reduce(&, <=(Sρ + hom + inhom, constraintProjBounds))
            end


            # We have already checked consraint, so now we only have to check constraint and guards
            if notHitConstraint
                # We start checking through guards, to see if we need to increase of decrease time
                approveFromGuards = true # If this is true, we keep going. Otherwise we reduce down to δ⁻
                for (count, (guard, intersectingPastStep)) in enumerate(currentlyIntersectingGuard)

                    intersectingNow = intersects(currentStepRepresentation, guard)


                    if !intersectingPastStep # If we were not intersecting
                        if intersectingNow # If we are now intersecting
                            if currentTimeStep == δ⁻ # We are at the lowest timestep. And must start intersecting
                                currentlyIntersectingGuard[count] = (guard, true)
                                startIntersectTimeTracker[count] = time
                            else  # We wanna keep from intersecting till we absolutely necessary
                                approveFromGuards = false
                                break
                            end
                        end
                    else # We are previously intersecting 
                        if !intersectingNow # We are not intersecting currently. 
                            if currentTimeStep == δ⁻ # We are at the lowest timestep. And must stop intersecting
                                # We push to the tattler
                                push!(guardIntersectTattler[count], (intersectBuildSpace[count], startIntersectTimeTracker))
                                intersectBuildSpace[count] = nothing # Reset
                                currentlyIntersectingGuard[count] = (guard, false)

                            else # We reduce time till we intersect
                                approveFromGuards = false 
                                break
                            end
                        end
                    end
                end
                if approveFromGuards
                    # For all currently intersecting, we add the intersecting space
                    for (count, (_, intersectingPastStep)) in enumerate(currentlyIntersectingGuard)
                        if intersectingPastStep
                            if isnothing(intersectBuildSpace[count]) # If first step in sequence
                                intersectBuildSpace[count] = copy(currentStepRepresentation)
                            else # If not first step, we convex hull them together to 1 set
                                intersectBuildSpace[count] = overapproximate(CH(intersectBuildSpace[count], currentStepRepresentation), Zonotope)
                            end
                        end
                    end

                    if saveResult
                        push!(reachSets, (currentStepRepresentation, [time, time + currentTimeStep]))
                    end
                    approveFlag = true
                    Sρ += inhom
                    Vs = tempVs
                    mul!(tempM, Φ, ϕt)
                    copy!(Φ, tempM)
                else
                    # Reduce timestep size to fulfill guards!
                    newR = copy(newR)
                    currentTimeStep = currentTimeStep / 2
                    changedTimeStep = true
                    attempts = attempts + 1
                end
            elseif !notHitConstraint && currentTimeStep == δ⁻
                # We are hitting a constraint and cannot increase precision further
                handleHitConstraint(time, loc.id)

                newR = copy(newR)
                currentTimeStep = currentTimeStep / 2
                changedTimeStep = true
            else
                # We hit a constraint or invarient and reduce time
                newR = copy(newR)
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

    #return (newR, sρ, time)
    return reachSets, Vs, time, Φ, guardIntersectTattler
end

#function ReACT(loc, δ⁻::Float64, δ⁺::Float64, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, STRATEGY::Integer, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, PhiDict=nothing) where {N}
function ReACT(loc, δ⁻::Float64, δ⁺::Float64, interval, constraint, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, Φ, saveResult)
    initialTimeStep = δ⁺
    changedTimeStep = true
    phiDict = PhiDict
    # discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    # inputDiscritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)

    #discritezationDict, inputDiscritezationDict, phiDict = ReACTDiscretize(A, X0, U, δ⁻, δ⁺, alg, maxOrder, reduceOrder, phiDict)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    attemptsRecorder = Integer[]

    #V::Zonotope{N,Vector{N},Matrix{N}} = copy(inputDiscritezationDict[initialTimeStep])
    V = copy(inputDiscritezationDict[initialTimeStep])
    Sρ = zeros(Float64, length(constraint))
    #newR::Zonotope{N,Vector{N},Matrix{N}} = discritezationDict[initialTimeStep]
    newR = discritezationDict[initialTimeStep]
    i = 1


    #Φ::Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    if ismissing(Φ)
        Φ::Matrix{Float64} = exp(0 .* loc.A)
    end
    tempM = similar(Φ)
    ϕt = similar(Φ)
    newRR = copy(newR)

    reachSets = []

    sρ = zeros(Float64, length(constraint))
    Vs = Zonotope(zeros(size(loc.A, 2)), [zeros(size(loc.A, 2))])

    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < δ⁻ 
                return (newR, sρ, time)
            end

            if changedTimeStep
                newR = discritezationDict[currentTimeStep]
                V = copy(inputDiscritezationDict[currentTimeStep])
                ϕt = phiDict[currentTimeStep]
                newRR = linear_map(Φ, newR)
                V = linear_map(Φ, V)
            else
                newRR = linear_map(ϕt, newRR)
                V = linear_map(ϕt, V)
            end

            changedTimeStep = false
            hom = map(x -> ρ(x, newRR), constraintProjVectors)
            inhom = map(x -> ρ(x, V), constraintProjVectors)
            tempVs = remove_redundant_generators(minkowski_sum(Vs, V))


            notHitConstraint = reduce(&, <=(Sρ + hom + inhom, constraintProjBounds))

            if notHitConstraint && intersects(concretize(minkowski_sum(newRR, Vs)), loc.invarient)
                if saveResult
                    push!(reachSets, (concretize(minkowski_sum(newRR, tempVs)), [time, time + currentTimeStep]))
                end
                approveFlag = true
                Sρ += inhom
                Vs = tempVs
                mul!(tempM, Φ, ϕt)
                copy!(Φ, tempM)

            elseif !notHitConstraint && currentTimeStep == δ⁻
                # We are hitting a constraint and cannot increase precision further
                handleHitConstraint(time, loc.id)

                newR = copy(newR)
                currentTimeStep = currentTimeStep / 2
                changedTimeStep = true
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

    #return (newR, sρ, time)
    return (reachSets, Vs, time, Φ)
end

# This can be replaced with throwing an error. Currently we continue and just print
function handleHitConstraint(time, locationId)
    throw(error("ERROR!!! We have hit a constraint at loc: $(locationId) time: $time"))
    # println("\nERROR!!!\n 
    #         ERROR!!!\n\n
    #         We have hit a constraint at loc: $(loc.id) time: $time\n\n
    #         ERROR!!!\n
    #         ERROR!!!\n")
end