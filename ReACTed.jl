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

function auxReACTed(hybridSystem, loc::Location, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, constraint, δ⁻::Float64, δ⁺::Float64, PhiDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, Φ=missing, saveResult::Bool=true) where {N}
    discretizationDict, inputDiscritezationDict = ReACTDiscretizePlus(loc, X0, δ⁻, δ⁺, alg, maxOrder, reduceOrder, PhiDict[loc.id])
    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []
    setOfConstraints = vcat(loc.constraints, constraint)

    #println("Smallest disc center: ", discretizationDict[δ⁻].center)

    for edge in loc.edges
        guards = edge.guard # Technically the guard is one singular HPolyhedron, but it composes the other guards

        #guards = constraints_list(guards)

        #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 

        tempReachset, tempInput, reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, [time, endtime], guards, setOfConstraints, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, missing, saveResult)
        #println(reachtime)
        #println(reachtime - time)
        if reachtime - time == 0.0
            return tempReachset
        end
        if reachtime < endtime
            if saveResult
                push!(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(edge.targetLoc)))
            end
            intersectingSet, intersectedInput, timeNotIntersected = ReACTTouches(loc, δ⁻, δ⁺, [reachtime, endtime], guards, setOfConstraints, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, tΦ, tempInput)
            timeIntersected = timeNotIntersected - reachtime
            #
            #   Here we should check whether we have reached endtime. If true we should only push the jumpSet
            #   Still need to check whether we have reached the invariant. If true we should NOT push the else branch result, only the tempReachsets[infMaxsIdx]
            #
            #println(timeIntersected)
            if timeIntersected >= 0.
                #timeIntersectedSet = concretize(overapproximateIntervalReachset(loc.A, tempReachsets[supMinsIdx], U, δ⁻, timeIntersected, alg, maxOrder, reduceOrder, flowPhiDict[loc.id]))
                #if !intersects(intersectingSet, guards)
                #throw(ErrorException("Set intersecting guard does not intersect the guard."))
                #end
                intersectedSet = intersectingSet

                if !isa(guards, Nothing) && !isdisjoint(intersectedSet, guards)
                    intersectedSet = zonotopeStripIntersection(intersectedSet, guards)
                end
                if isempty(intersectedSet)
                    println("Empty..")
                    return reachset
                end
                #println("Guard intersectedSet: ", intersectedSet)

                if !isa(loc.invarient, Nothing) && !isdisjoint(intersectedSet, loc.invarient)
                    #println("tes")
                    tintersectedSet = zonotopeStripIntersection(intersectedSet, loc.invarient)
                    #println("Inv intersection: ", ρ(Vector(sparsevec([5], [1.0], 6)), intersectedSet), " vs ", ρ(Vector(sparsevec([5], [1.0], 6)), tintersectedSet))
                    intersectedSet = tintersectedSet
                    #println(intersectedSet)
                    #println("Inv intersectedSet: ", intersectedSet)
                end
                if isempty(intersectedSet)
                    println("Empty...")
                    return reachset
                end
                push!(reachset, ([(intersectedSet, [reachtime, timeNotIntersected])], string(reachtime) * " - " * string(timeNotIntersected) * ": " * string(loc.id) * "->" * string(loc.id)))

                if !isdisjoint(intersectedSet, guards)
                    intersectedSet = zonotopeStripIntersection(intersectedSet, guards)


                    #println("Guard intersectedSet: ", intersectedSet)


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
                    jumpSet = minkowski_sum(linear_map(edge.jumpMatrix, intersectedSet), Zonotope(edge.jumpVector, [zero(edge.jumpVector)]))
                    #println(jumpSet)

                    if !isa(hybridSystem.locations[edge.targetLoc].invarient, Nothing) && intersects(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)
                        tjumpSet = zonotopeStripIntersection(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)# + edge.jumpVector * intersectedSet
                        #println("Jump intersection: ", ρ(Vector(sparsevec([5], [1.0], 6)), jumpSet), " vs ", ρ(Vector(sparsevec([5], [1.0], 6)), tjumpSet))
                        jumpSet = tjumpSet
                        #println(LazySets.order(jumpSet))
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
                    #println("Jumpset center: ", jumpSet.center)
                    #y, _ = tempReachset[1]
                    #println(y.center)
                    branchedRun = auxReACTed(hybridSystem, hybridSystem.locations[edge.targetLoc], [reachtime, endtime], jumpSet, constraint, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder, tΦ, saveResult)

                    if saveResult
                        reachset = vcat(reachset, branchedRun)
                    end
                else
                    nonintersectedSet = minkowski_sum(linear_map(exp(timeNotIntersected .* loc.A), X0), concretize(ReachabilityAnalysis.Exponentiation.Φ₁(loc.A, timeNotIntersected, ReachabilityAnalysis.Exponentiation.BaseExp) * inputDiscritezationDict[0]))
                    if !isa(loc.invarient, Nothing) && !isdisjoint(nonintersectedSet, loc.invarient)
                        nonintersectedSet = zonotopeStripIntersection(nonintersectedSet, loc.invarient)
                    end
                    branchedRun = auxReACTed(hybridSystem, loc, [timeNotIntersected, endtime], nonintersectedSet, constraint, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder)
                    if saveResult
                        reachset = vcat(reachset, branchedRun)
                    end
                end
                #reachset = vcat(reachset, nonintersectedSet) #Maybe gets the universe..

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
        else
            println("Is here?")
            if saveResult
                if !isa(loc.invarient, Nothing)
                    newReach = []
                    for (Z, timeInterval) in tempReachset
                        if intersects(Z, loc.invarient)
                            push!(newReach, (zonotopeStripIntersection(concretize(Z), loc.invarient), timeInterval))
                        else
                            push!(newReach, (Z, timeInterval))
                        end
                    end
                    tempReachset = newReach
                    reachset = vcat(reachset, (tempReachset, string(loc.id) * "->" * string(edge.targetLoc)))
                else
                    reachset = vcat(reachset, (tempReachset, string(loc.id) * "->" * string(edge.targetLoc)))
                end
            end
        end
    end

    if isempty(loc.edges) # This means it is just a continous system from here
        println("No edges? call ReACT")

        #tempReachset, _, _ = ReACT(loc, δ⁻, δ⁺, interval, constraint, 2, alg, maxOrder, reduceOrder, PhiDict[loc.id])
        tempReachset, _, reachtime, _ = ReACT(loc, δ⁻, δ⁺, interval, constraint, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, Φ, saveResult)

        if saveResult
            if !isa(loc.invarient, Nothing)
                newReach = []
                for (Z, timeInterval) in tempReachset
                    if intersects(Z, loc.invarient)
                        push!(newReach, (zonotopeStripIntersection(concretize(Z), loc.invarient), timeInterval))
                    else
                        push!(newReach, (Z, timeInterval))
                    end
                end

                reachset = vcat(reachset, (newReach, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(loc.id)))
            else
                reachset = vcat(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(loc.id)))
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
    Vs = concretize(ReachabilityAnalysis.Exponentiation.Φ₁(loc.A, time, ReachabilityAnalysis.Exponentiation.BaseExp) * inputDiscritezationDict[0])
    lastVs = copy(Vs)
    Sρ = zeros(Float64, length(constraint))
    newR = copy(discritezationDict[initialTimeStep])
    i = 1


    if ismissing(Φ)
        Φ::Matrix{Float64} = exp(time .* loc.A)
    end
    tempM = diagm(ones(Float64, size(loc.A, 2)))
    ϕt = similar(Φ)
    newRR = copy(newR)

    overapproximateIntersectingSetArray = []
    preclustering = minkowski_sum(linear_map(Φ, discritezationDict[m]), Vs)
    attemptsRecorder = []

    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m


                #if !reduce((x, y -> x && y), <=(Sρ + map(x -> ρ(x, newRR), constraintProjVectors), constraintProjBounds))
                # Any has short-circuit, and returns false if no constraints
                if any(((input + ρ(x, newRR)) > y) for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))
                    #throw(ErrorException("Reached unsafe set."))
                    handleHitConstraint(time, loc.id)
                end
                #bigCH = foldr((x, y) -> overapproximate(CH(x, y), Zonotope), overapproximateIntersectingSetArray; init=concretize(newRR))
                #intersectingSet = overapproximate(bigCH, Zonotope)
                if ismissing(preclustering)
                    preclustering = concretize(minkowski_sum(newRR, Vs))
                end
                #println("Touches Vs: ", accInput)
                #println(norm(intersectingSet), " ", norm(newRR))
                return (reduce_order(preclustering, 5), Sρ, time)
            end

            if changedTimeStep
                newR = copy(discritezationDict[currentTimeStep])
                V = copy(inputDiscritezationDict[currentTimeStep])
                ϕt = phiDict[currentTimeStep]
                newRR = linear_map(Φ, newR)
                V = linear_map(Φ, V)
            else
                newRR = linear_map(ϕt, newRR)
                V = linear_map(ϕt, V)
            end

            changedTimeStep = false
            # hom = map(x -> ρ(x, newRR), constraintProjVectors)
            # inhom = map(x -> ρ(x, Vs), constraintProjVectors)
            #tempVs = concretize(ReachabilityAnalysis.Exponentiation.Φ₁(loc.A, time, ReachabilityAnalysis.Exponentiation.BaseExp) * inputDiscritezationDict[0])


            #if all(&, <=(hom + inhom, constraintProjBounds)) && intersects(concretize(minkowski_sum(newRR, Vs)), guard) && (isnothing(loc.invarient) || isSubSet(concretize(minkowski_sum(newRR, Vs)), loc.invarient)) #mapreduce(x -> intersects(newRR, x), &, guard)
                
            if all((ρ(x, newRR) + ρ(x, Vs)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) && intersects(concretize(minkowski_sum(newRR, Vs)), guard) && (isnothing(loc.invarient) || isSubSet(concretize(minkowski_sum(newRR, Vs)), loc.invarient)) #mapreduce(x -> intersects(newRR, x), &, guard)
                #if mapreduce(x -> intersects(newRR, x), &, guard)
                push!(overapproximateIntersectingSetArray, concretize(minkowski_sum(newRR, Vs)))
                if ismissing(preclustering)
                    preclustering = concretize(minkowski_sum(newRR, Vs))
                else
                    preclustering = overapproximate(CH(preclustering, concretize(minkowski_sum(newRR, Vs))), Zonotope)
                end
                lastVs = copy(Vs)
                Vs = concretize(ReachabilityAnalysis.Exponentiation.Φ₁(loc.A, time + currentTimeStep, ReachabilityAnalysis.Exponentiation.BaseExp) * inputDiscritezationDict[0])
                approveFlag = true
                #Sρ += inhom
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
    #bigCH = foldr((x, y) -> overapproximate(CH(x, y), Zonotope), overapproximateIntersectingSetArray; init=concretize(newRR))
    #intersectingSet = overapproximate(bigCH, Zonotope)
    return (preclustering, Sρ, time)
end

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

                if any((input + ρ(x, runningset)) <= y for (input, x, y) in zip(Sp, constraintProjVectors, constraintProjBounds))
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
                println("Disjoint! ", i)
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
                    preclustering = overapproximate(CH(preclustering, runningset), Zonotope)
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

function ReACTGuards(loc, δ⁻::Float64, δ⁺::Float64, interval, guards, constraint, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, Φ, saveResult)
    initialTimeStep = copy(δ⁺)
    m = copy(δ⁻)
    changedTimeStep = true
    phiDict = PhiDict
    A = copy(loc.A)
    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(δ⁻)#copy(initialTimeStep)

    #overapproximateIntersectingSetArray = []
    lastNewR = []
    attemptsRecorder = []


    if ismissing(Φ)
        Φ::Matrix{Float64} = exp(0 .* loc.A)
    end
    tempM = diagm(ones(Float64, size(loc.A, 2)))
    ϕt = similar(tempM)

    V = copy(inputDiscritezationDict[initialTimeStep])
    Vs = Zonotope(zeros(size(loc.A, 2)), [zeros(size(loc.A, 2))])
    lastVs = copy(Vs)
    Sρ = zeros(Float64, length(constraint))
    newR = discritezationDict[initialTimeStep]
    i = 1

    newRR = copy(newR)
    U = inputDiscritezationDict[0]

    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m
                if isempty(lastNewR) && intersects(newRR, guards)
                    newRR = concretize(newRR)
                    
                    if any((input + ρ(x, newRR)) > y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))
                    #if !reduce(&, <=(Sρ + map(x -> ρ(x, concretize(newRR)), constraintProjVectors), constraintProjBounds))
                        handleHitConstraint(time, loc.id)
                    else
                        #println("Pushing!")
                        #push!(lastNewR, (concretize(newRR), [time, time + currentTimeStep]))
                        #println(!intersects(concretize(minkowski_sum(newRR, Vs)), guards))
                        #println(guards)
                        #println((isnothing(loc.invarient) || intersects(concretize(minkowski_sum(newRR, Vs)), loc.invarient)))
                    end
                end
                #println(intersects(lastNewR[end][1], guards))
                #push!(overapproximateIntersectingSetArray, newRR)
                #intersectingSet = overapproximate(ConvexHullArray(overapproximateIntersectingSetArray), Zonotope)
                #println(norm(lastNewR), " ", time)
                # println("Guards Vs: ", concretize(minkowski_sum(newRR, Vs)))
                #println("Guards newRR: ", pop!(lastNewR))
                #newRR = minkowski_sum(newRR, Vs)
                #println("x :", ρ(sparsevec([2],[1.], 5), newRR))
                #println("t :", ρ(sparsevec([5],[1.], 5), newRR))
                #println("t :", ρ(sparsevec([5],[-1.], 5), newRR))
                #println(intersects(newRR, guards))
                return (lastNewR, Vs, time, Φ)
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

            tempSet = concretize(minkowski_sum(newRR, Vs))
            #hom = map(x -> ρ(x, tempSet), constraintProjVectors)

            if all(((ρ(x, tempSet)) <= y) for (x, y) in zip(constraintProjVectors, constraintProjBounds)) && 
                !intersects(tempSet, guards) && 
                (isnothing(loc.invarient) || intersects(tempSet, loc.invarient))
                #if mapreduce(x -> intersects(newRR, x), &, guard)
                #push!(overapproximateIntersectingSetArray, newRR)
                if saveResult
                    if isempty(lastNewR)
                        #println("first element: ", concretize(newRR).center, " ", currentTimeStep)
                    end
                    push!(lastNewR, (tempSet, [time, time + currentTimeStep]))
                end
                lastVs = copy(Vs)
                Vs = concretize(ReachabilityAnalysis.Exponentiation.Φ₁(A, time, ReachabilityAnalysis.Exponentiation.BaseExp) * U)
                approveFlag = true
                Sρ += map(x -> ρ(x, Vs), constraintProjVectors)
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
    #intersectingSet = overapproximate(ConvexHullArray(overapproximateIntersectingSetArray), Zonotope)
    #println(lastNewR)
    return (lastNewR, Vs, time, Φ)
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
                return (reachSets, Vs, time, Φ)
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
            #hom = map(x -> ρ(x, newRR), constraintProjVectors)
            inhom = map(x -> ρ(x, V), constraintProjVectors)
            tempVs = remove_redundant_generators(minkowski_sum(Vs, V))


            #notHitConstraint = reduce(&, <=(Sρ + hom + inhom, constraintProjBounds))
            # TODO should we include both Sp and V?
            notHitConstraint = all((input + input2 + ρ(x, newRR)) <= y for (input, input2, x, y) in zip(Sρ, inhom, constraintProjVectors, constraintProjBounds))

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
    #throw(error("ERROR!!! We have hit a constraint at loc: $(locationId) time: $time"))
    println("\nERROR!!!\n 
            ERROR!!!\n\n
            We have hit a constraint at loc: $locationId time: $time\n\n
            ERROR!!!\n
            ERROR!!!\n")
end