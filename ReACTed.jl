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

    for edge in loc.edges
        guards = edge.guard
        if !ismissing(guards)
            #guards = constraints_list(guards)

            #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 

            tempReachset, tempInput, reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, [time, endtime], guards, setOfConstraints, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, missing, saveResult)
            #println(reachtime)
            println(reachtime - time)
            if reachtime - time == 0.0
                return reachset
            end
            if reachtime < endtime
                if saveResult
                    push!(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(edge.targetLoc)))
                end
                intersectingSet, intersectedInput, timeNotIntersected = ReACTTouches(loc, δ⁻, δ⁺, [reachtime, endtime], guards, setOfConstraints, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, tΦ, tempInput)
                timeIntersected = timeNotIntersected - reachtime
                #push!(reachset, (intersectingSet, string(loc.id) * "->" * string(edge.targetLoc)))
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



                    if !isa(loc.invarient, Nothing) && intersects(intersectedSet, loc.invarient)
                        #println("tes")
                        tintersectedSet = getBoxIntersection(intersectedSet, loc.invarient)
                        #println("Inv intersection: ", ρ(Vector(sparsevec([5], [1.0], 6)), intersectedSet), " vs ", ρ(Vector(sparsevec([5], [1.0], 6)), tintersectedSet))
                        intersectedSet = tintersectedSet
                        #println(intersectedSet)
                    end
                    if isempty(intersectedSet)
                        println("Empty...")
                        return reachset
                    end
                    println("Inv intersectedSet: ", intersectedSet)
                    intersectedSet = getBoxIntersection(intersectedSet, guards)
                    if isempty(intersectedSet)
                        println("Empty..")
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
                    branchedRun = auxReACTed(hybridSystem, hybridSystem.locations[edge.targetLoc], [reachtime, endtime], jumpSet, constraint, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder, tΦ, saveResult)
                    
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
            else
                println("Is here?")
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
                        tempReachset = newReach
                        reachset = vcat(reachset, (tempReachset, string(loc.id) * "->" * string(edge.targetLoc)))
                    else
                        reachset = vcat(reachset, (tempReachset, string(loc.id) * "->" * string(edge.targetLoc)))
                    end
                end
            end
        else
            println("No guards? call ReACT")

            tempReachset, _, _ = ReACT(loc, δ⁻, δ⁺, interval, X0, U, constraint, 2, alg, maxOrder, reduceOrder, PhiDict[loc.id])
            
            if saveResult
                if intersects(concretize(tempReachset), constraints_list(loc.invarient))

                    tempReachsetInInv = getBoxIntersection(concretize(tempReachset), constraints_list(loc.invarient))
                    reachset = vcat(reachset, (tempReachsetInInv, string(loc.id) * "->" * string(edge.targetLoc)))
                else
                    reachset = vcat(reachset, (concretize(tempReachset), string(loc.id) * "->" * string(edge.targetLoc)))
                end
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

function ReACTGuards(loc, δ⁻::Float64, δ⁺::Float64, interval, guards, constraint, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, Φ, saveResult)
    initialTimeStep = copy(δ⁺)
    m = copy(δ⁻)
    changedTimeStep = true
    phiDict = PhiDict
    A = copy(loc.A)
    println(A)
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


    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m
                if isempty(lastNewR) && intersects(newRR, guards)
                    if !reduce(&, <=(Sρ + map(x -> ρ(x, concretize(newRR)), constraintProjVectors), constraintProjBounds))
                        handleHitConstraint(time, loc.id)
                    else
                        # println("Pushing!")
                        # push!(lastNewR, (concretize(newRR), [time, time + currentTimeStep]))
                    end
                end

                #push!(overapproximateIntersectingSetArray, newRR)
                #intersectingSet = overapproximate(ConvexHullArray(overapproximateIntersectingSetArray), Zonotope)
                #println(norm(lastNewR), " ", time)
                # println("Guards Vs: ", concretize(minkowski_sum(newRR, Vs)))
               #  println("Guards newRR: ", pop!(lastNewR))
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
            hom = map(x -> ρ(x, newRR), constraintProjVectors)
            inhom = map(x -> ρ(x, V), constraintProjVectors)

            tempVs = remove_redundant_generators(minkowski_sum(Vs, V))
            if reduce(&, <=(Sρ + hom + inhom, constraintProjBounds)) && !intersects(concretize(minkowski_sum(newRR, Vs)), guards) && intersects(concretize(minkowski_sum(newRR, Vs)), loc.invarient)
                #if mapreduce(x -> intersects(newRR, x), &, guard)
                #push!(overapproximateIntersectingSetArray, newRR)
                if saveResult
                    if isempty(lastNewR)
                        println("first element: ", concretize(newRR).center, " ", currentTimeStep)
                    end
                    push!(lastNewR, (concretize(minkowski_sum(newRR, Vs)), [time, time + currentTimeStep]))
                end
                lastVs = copy(Vs)
                Vs = reduce_order(copy(tempVs), 6)
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
    #intersectingSet = overapproximate(ConvexHullArray(overapproximateIntersectingSetArray), Zonotope)
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