using LazySets, LinearAlgebra, ReachabilityAnalysis
include("Discretize.jl")
include("Utilities.jl")


const VERBOSE = false

function zonotopePrintDim(Z::Zonotope, dim)
    center = Z.center[dim]
    G = genmat(Z)
    _, genAmount = size(G)
    genContribute = sum(abs(G[dim, i]) for i in 1:genAmount)
    return "center: $center, generators: $genContribute"
end

function ReACTed(hybridSystem::HybridSystemV2, initialLoc, interval, X0, U, dirs, constraint, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    loc = initialLoc
    flowPhiDict = Dict(map(x -> x.id => PhiDict(x.A, δ⁻, δ⁺, alg), hybridSystem.locations))

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []
    saveResult = true
    dirsVectors = []

    if isempty(dirs)
        saveResult = false
    else
        dimLength = size(X0.center, 1)
        for dir in dirs
            template = zeros(dimLength)
            template[dir] = 1.

            push!(dirsVectors, template)
            push!(dirsVectors, -template)
        end
    end

    @show dirsVectors

    res = auxReACTed(hybridSystem, hybridSystem.locations[loc], interval, X0, dirsVectors, constraint, δ⁻, δ⁺, flowPhiDict, alg, maxOrder, reduceOrder, missing, nothing, nothing, saveResult)

    reachset = vcat(reachset, res)

    return res
end


# AucReacted is called recursively each time we have a new starting location (after a transition)
function auxReACTed(hybridSystem, loc::Location, interval, X0, dirs, constraint, δ⁻::Float64, δ⁺::Float64, PhiDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, Φ=missing, discretizationDict = nothing, inputDiscritezationDict = nothing, saveResult::Bool=true) where {N}
    locChange = false
    jumpDiscDict = Dict()
    jumpInputDict = Dict()

    if isa(discretizationDict, Nothing) && isa(inputDiscritezationDict, Nothing)
        discretizationDict, inputDiscritezationDict = ReACTDiscretizePlus(loc, X0, δ⁻, δ⁺, alg, maxOrder, reduceOrder, PhiDict[loc.id])
        locChange = true
    end
    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []
    setOfConstraints = vcat(loc.constraints, constraint)


    # For each edge we simulate the system
    for edge in loc.edges
        println("Handling edge at time $time from $(loc.id) -> $(edge.targetLoc)")
        guards = edge.guard # Technically the guard is one singular HPolyhedron, but it composes the other guards


        #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 

        tempReachset, reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, [time, endtime], guards, setOfConstraints, dirs, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, saveResult)

        if reachtime - time == 0.0
            if !isempty(reachset)
                println("Found an immediate transition to $(edge.targetLoc), but we are not taking it as we are scared of zeno behaviour")
            end
            if saveResult
                push!(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * " took no steps"))
            end
            continue
        end
        if reachtime < endtime
            if saveResult
                push!(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(edge.targetLoc)))
            end
            println("Current time: $time, intersecting time start: $reachtime")
            intersectingSet, intersectedInput, timeNotIntersected = ReACTTouches(loc, δ⁻, δ⁺, [reachtime, endtime], (reachtime - time), guards, setOfConstraints, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, tΦ, nothing)

            timeIntersected = timeNotIntersected - reachtime

            #
            #   Here we should check whether we have reached endtime. If true we should only push the jumpSet
            #   Still need to check whether we have reached the invariant. If true we should NOT push the else branch result, only the tempReachsets[infMaxsIdx]
            #


            if !isnothing(intersectingSet)
                intersectedSet = overapproximate(concretize(intersectingSet), Zonotope)

                if !isa(guards, Nothing)
                    if isdisjoint(intersectedSet, guards)
                        println("Empty intersection with guard")
                        return reachset
                    end
                    
                    intersectedSet = zonotopeStripIntersection(intersectedSet, guards)
                end


                if !isa(loc.invarient, Nothing)
                    if isdisjoint(loc.invarient, intersectedSet)
                        println("Empty intersection with invarient")
                        return reachset
                    end
                    tintersectedSet = zonotopeStripIntersection(intersectedSet, loc.invarient)
                    intersectedSet = tintersectedSet
                end
                push!(reachset, ([(map(x -> ρ(x, intersectedSet), dirs), [reachtime, timeNotIntersected])], string(reachtime) * " - " * string(timeNotIntersected) * ": " * string(loc.id) * "->" * string(loc.id)))
                if !isdisjoint(intersectedSet, guards)
                    if locChange
                        # reset_map(X) = edge.jumpMatrix * X + edge.jumpVector
                        # jumpSet = reset_map(intersectedSet) 
                        jumpSet = linear_map(edge.jumpMatrix, intersectedSet)
                        jumpSet = LazySets.translate(jumpSet, edge.jumpVector)

                        if !isa(hybridSystem.locations[edge.targetLoc].invarient, Nothing) && !isdisjoint(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)
                            jumpSet = zonotopeStripIntersection(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)
                        end 
                    
                        #
                        #   Here we could optimize it such that in the case where guards ⊆ timeIntersectedSet we calculate both [supMins, endtime] and [infMaxs, endtime] with guards
                        #   and otherwise [supMins, endtime] with hyperplane intersection with timeIntersectedSet and [infMaxs, endtime] with guards intersection
                        #   Maybe look at how input should be handled... and if we can manipulate the constraints to account for the accumulated input
                        #
                    
                        #timePointInput = U  #   NEEDS FIXING
                        #println(x.center)
                        #println("Jumpset center: ", jumpSet.center)
                        #y, _ = tempReachset[1]
                        println("Finished intersections")
                    
                        println("Going this way")
                        branchedRun = auxReACTed(hybridSystem, hybridSystem.locations[edge.targetLoc], [reachtime, endtime], jumpSet, dirs, constraint, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder, tΦ, nothing, nothing, saveResult)
                    
                        if saveResult
                            reachset = vcat(reachset, branchedRun)
                        end
                    else
                        λ = ceil((timeIntersected)/δ⁺)
                        branchNumber = 1
                        while branchNumber <= λ
                            d = δ⁻
                            jumpΦ = tΦ
                            while d <= δ⁺
                                tempSumAffine = minkowski_sum(linear_map(edge.jumpMatrix, linear_map(jumpΦ, discretizationDict[d])), edge.jumpVector)
                                tempSumInput = minkowski_sum(linear_map(edge.jumpMatrix, linear_map(jumpΦ, inputDiscritezationDict[d])), edge.jumpVector)

                                jumpDiscDict[d] = tempSumAffine
                                jumpInputDict[d] = tempSumInput
                                if !isa(hybridSystem.locations[edge.targetLoc].invarient, Nothing) && !isdisjoint(minkowski_sum(tempSumAffine, tempSumInput), hybridSystem.locations[edge.targetLoc].invarient; algorithm="sufficient")
                                    jumpDiscDict[d] = Intersection(hybridSystem.locations[edge.targetLoc].invarient, tempSumAffine)
                                    jumpInputDict[d] = Intersection(hybridSystem.locations[edge.targetLoc].invarient, tempSumInput)
                                end
                            end
                            jumpΦ = jumpΦ * PhiDict[δ⁺]
                            branchedRun = auxReACTed(hybridSystem, hybridSystem.locations[edge.targetLoc], [reachtime, endtime], nothing, dirs, constraint, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder, tΦ, jumpDiscDict, jumpInputDict, saveResult)
                        end
                    end
                else# We do not hit guards, and cannot transition
                    println("MOSHIMOSHI")
                    nonintersectedSet = MinkowskiSum(LinearMap(exp(timeNotIntersected .* loc.A), X0), ReachabilityAnalysis.Exponentiation.Φ₁(loc.A, timeNotIntersected, ReachabilityAnalysis.Exponentiation.BaseExp) * inputDiscritezationDict[0])
                    if !isa(loc.invarient, Nothing) && !isdisjoint(nonintersectedSet, loc.invarient; algorithm="sufficient")
                        nonintersectedSet = nonintersectedSet ∩ loc.invarient
                        #nonintersectedSet = zonotopeStripIntersection(nonintersectedSet, loc.invarient)
                    end
                    branchedRun = auxReACTed(hybridSystem, loc, [timeNotIntersected, endtime], nonintersectedSet, constraint, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder)
                    if saveResult
                        reachset = vcat(reachset, branchedRun)
                    end
                end
            else
                println("We have no intersecting set. Meaning we hit an invarient and have no guards fulfilled")
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

        else # Reached the end time
            println("Is here?")
            if saveResult
                if !isa(loc.invarient, Nothing)
                    newReach = []
                    for (Z, timeInterval) in tempReachset

                        push!(newReach, (Z, timeInterval))
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

        #tempReachset, tempInput, reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, [time, endtime], guards, setOfConstraints, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, missing, saveResult)


        #tempReachset, _, _ = ReACT(loc, δ⁻, δ⁺, interval, constraint, 2, alg, maxOrder, reduceOrder, PhiDict[loc.id])
        tempReachset, reachtime, _ = ReACT(loc, δ⁻, δ⁺, [time, endtime], setOfConstraints, dirs, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, saveResult)

        if saveResult
            reachset = vcat(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(loc.id)))
        end
    end




    return reachset
end

function ReACTTouches(loc, δ⁻::Float64, δ⁺::Float64, interval, initialTime::Float64, guard, constraint, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, Φ, accInput)
    initialTimeStep = copy(δ⁺)
    m = copy(δ⁻)
    changedTimeStep = true
    phiDict = PhiDict

    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(guard)
    invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    V = copy(inputDiscritezationDict[initialTimeStep])

    Vs = nestedInputDiscCalculate(inputDiscritezationDict, PhiDict, δ⁺, δ⁻, initialTime)


    #concretize(Vs)
    #lastVs = copy(Vs)
    Sρ = zeros(Float64, length(constraint))
    newR = discritezationDict[initialTimeStep]
    i = 1


    if ismissing(Φ)
        Φ::Matrix{Float64} = exp(initialTime .* loc.A)
    end

    tempM = diagm(ones(Float64, size(loc.A, 2)))
    ϕt = similar(Φ)
    #newRR = copy(newR)
    newRR = newR

    overapproximateIntersectingSetArray = []
    #preclustering = (Φ * discritezationDict[m]) ⊕ Vs
    preclustering = nothing
    attemptsRecorder = []

    while time < endtime
        #println("Touching at time $time with newRR: ", zonotopePrintDim(newRR , 5))
        println("Time $time, calculating intersect")

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m

                if any(((input + ρ(x, newRR)) > y) for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))
                    #throw(ErrorException("Reached unsafe set."))
                    handleHitConstraint(time, loc.id)
                end

                if all(((-ρ(-x, newRR ⊕ Vs)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds))
                    println("Intersecting guard, but we are going to hit an invarient")

                    if isnothing(preclustering)
                        preclustering = newRR ⊕ Vs
                    else
                        preclustering = UnionSet(preclustering, newRR ⊕ Vs)
                    end

                end

                return preclustering, Sρ, time


                #println("Touches Vs: ", accInput)
                #println(norm(intersectingSet), " ", norm(newRR))
            end

            if changedTimeStep
                newR = discritezationDict[currentTimeStep]
                V = copy(inputDiscritezationDict[currentTimeStep])
                ϕt = phiDict[currentTimeStep]
                newRR = LinearMap(Φ, newR)
                V = LinearMap(Φ, V)
            else
                newRR = LinearMap(ϕt, newRR)
                V = LinearMap(ϕt, V)
            end

            changedTimeStep = false


            tempSet = newRR ⊕ Vs

            if all((ρ(x, newRR) + ρ(x, Vs)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) && # IsSubSet
               all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds)) && 
               #!isdisjoint(tempSet, guard; algorithm="sufficient") && # intersects
               all((ρ(x, newRR) + ρ(x, Vs)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # IsSubSet
                #if mapreduce(x -> intersects(newRR, x), &, guard)
                push!(overapproximateIntersectingSetArray, newRR ⊕ Vs)
                if isnothing(preclustering)
                    preclustering = MinkowskiSum(newRR, Vs)
                else
                    preclustering = UnionSet(preclustering, newRR ⊕ Vs)
                end
                #lastVs = copy(Vs)
                Vs = Vs ⊕ V
                #Vs = LinearMap(ReachabilityAnalysis.Exponentiation.Φ₁(loc.A, time + currentTimeStep - minimum(interval), ReachabilityAnalysis.Exponentiation.BaseExp), inputDiscritezationDict[0])
                approveFlag = true
                #Sρ += inhom
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




function ReACTGuards(loc, δ⁻::Float64, δ⁺::Float64, interval, guards, constraint, dirs, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, saveResult)
    # We calculate the reachset till we reach a guard for an intersection (or till failure)
    initialTimeStep = copy(δ⁺)
    m = copy(δ⁻)
    changedTimeStep = true
    phiDict = PhiDict
    A = copy(loc.A)
    # constraintProjVectors = map(x -> x.a, constraint)
    # constraintProjBounds = ρ.(constraintProjVectors, constraint)
    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(guards)
    invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)
    dirProjVectors = map(x -> x, dirs)

    oldDirProjVectors = copy(dirProjVectors)
    oldConstraintProjVectors = copy(constraintProjVectors)
    oldGuardProjVectors = copy(guardProjVectors)
    oldInvarientProjVectors = copy(invarientProjVectors)

    permutedphiDict = Dict()
    for key in keys(phiDict)
        permutedphiDict[key] = permutedims(phiDict[key])
    end


    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(δ⁻)#copy(initialTimeStep)

    #overapproximateIntersectingSetArray = []
    #lastNewR = []
    attemptsRecorder = []

    dirVals = [] # This is where the plotting happens


    Φ::Matrix{Float64} = exp(0 .* loc.A)

    tempM = diagm(ones(Float64, size(loc.A, 2)))
    ϕt = similar(tempM)
    pϕt = similar(tempM)
    V = copy(inputDiscritezationDict[initialTimeStep])

    #lastVs = copy(Vs)
    Sρ = zeros(Float64, length(constraint))
    Gρ = zeros(Float64, length(guardProjVectors))
    Iρ = zeros(Float64, length(invarientProjVectors))
    dρ = zeros(Float64, length(dirs))
    newR = discritezationDict[initialTimeStep]
    i = 1

    U = inputDiscritezationDict[0]
    while time < endtime
        #println("Time iss: $time")
        attempts = 1
        approveFlag = false

        while !approveFlag
            #println("Stuck?")

            # Handle if we can no longer reduce the reachset (we keep hitting something)
            if currentTimeStep < m
                # If we hit a constraint
                #newRR = concretize(newRR)
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
 
            # #println(any((input + ρ(-x, newRR)) > y for (input, x, y) in zip(Sρ, invarientProjVectors, invarientProjBounds)))
            #println(all((input + ρ(x, newRR)) <= y for (input, x, y) in zip(Sρ, invarientProjVectors, invarientProjBounds)))
            if all((input + ρ(x, newR)) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) &&
                #any(sign(y) >= 0 ? (input + ρ(-x, newR)) <= y : !((input + ρ(x, newR)) < y) for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) &&
                any((input + -ρ(-x, newR)) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) &&
                all((input + -ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))
                #(any((input + ρ(-x, newRR)) > y for (input, x, y) in zip(-Sρ, invarientProjVectors, invarientProjBounds)) && all((input + ρ(x, newRR)) <= y for (input, x, y) in zip(Sρ, invarientProjVectors, invarientProjBounds)))
                #all(((ρ(x, tempSet)) <= y) for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Subset
                #(isnothing(loc.invarient) || intersects(tempSet, loc.invarient))
                #if mapreduce(x -> intersects(newRR, x), &, guard)
                #push!(overapproximateIntersectingSetArray, newRR)
                if saveResult
                    push!(dirVals, (copy(dρ + map(x -> ρ(x, newR), oldDirProjVectors)), [time, time + currentTimeStep]))
                end
                #lastVs = copy(Vs)
                #Vs = ReachabilityAnalysis.Exponentiation.Φ₁(A, time - minimum(interval), ReachabilityAnalysis.Exponentiation.BaseExp) * U

                approveFlag =
                 true
                dρ += map(x -> ρ(x, V), oldDirProjVectors)
                Sρ += map(x -> ρ(x, V), constraintProjVectors)
                Gρ += map(x -> ρ(x, V), guardProjVectors)
                Iρ += map(x -> ρ(x, V), invarientProjVectors)
                dirProjVectors = map(x -> pϕt * x, oldDirProjVectors)
                constraintProjVectors = map(x -> pϕt * x, oldConstraintProjVectors)
                guardProjVectors = map(x -> pϕt * x, oldGuardProjVectors)
                invarientProjVectors = map(x -> pϕt * x, oldInvarientProjVectors)

                oldDirProjVectors = dirProjVectors
                oldConstraintProjVectors = constraintProjVectors
                oldGuardProjVectors = guardProjVectors
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
    #intersectingSet = overapproximate(ConvexHullArray(overapproximateIntersectingSetArray), Zonotope)
    #println(lastNewR)
    return (dirVals, time, Φ)
end

function ReACT(loc, δ⁻::Float64, δ⁺::Float64, interval, constraint, dirs, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, saveResult)
    # We calculate the reachset till we reach a guard for an intersection (or till failure)
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