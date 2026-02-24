using LazySets, LinearAlgebra
include("Discretize.jl")
include("Utilities.jl")


function ReACTed(hybridSystem::HybridSystemV2, initialLoc, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    loc = initialLoc
    flowPhiDict = Dict(map(x -> x.id => PhiDict(x.A, δ⁻, δ⁺, alg), hybridSystem.locations))




    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []

    res = auxReACTed(hybridSystem, hybridSystem.locations[loc], interval, X0, constraint, δ⁻, δ⁺, flowPhiDict, alg, maxOrder, reduceOrder)

    reachset = vcat(reachset, res)

    return res
end

function auxReACTed(hybridSystem, loc::Location, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, constraint, δ⁻::Float64, δ⁺::Float64, PhiDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, Φ=missing) where {N}
    discretizationDict, inputDiscritezationDict = ReACTDiscretize(loc, X0, δ⁻, δ⁺, alg, maxOrder, reduceOrder, PhiDict[loc.id])
    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []
    setOfConstraints = vcat(loc.constraints, constraint)

    for edge in loc.edges
        guards = edge.guard
        if !ismissing(guards)
            #guards = constraints_list(guards)

            #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 

            tempReachset, tempInput, reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, [time, endtime], guards, setOfConstraints, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, Φ)
            if reachtime < endtime
                push!(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(edge.targetLoc)))
                intersectingSet, intersectedInput, timeNotIntersected = ReACTTouches(loc, δ⁻, δ⁺, [reachtime, endtime], guards, setOfConstraints, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, tΦ)
                timeIntersected = timeNotIntersected - reachtime
                #push!(reachset, (intersectingSet, string(loc.id) * "->" * string(edge.targetLoc)))
                #
                #   Here we should check whether we have reached endtime. If true we should only push the jumpSet
                #   Still need to check whether we have reached the invariant. If true we should NOT push the else branch result, only the tempReachsets[infMaxsIdx]
                #

                if timeIntersected > 0.
                    #timeIntersectedSet = concretize(overapproximateIntervalReachset(loc.A, tempReachsets[supMinsIdx], U, δ⁻, timeIntersected, alg, maxOrder, reduceOrder, flowPhiDict[loc.id]))
                    if !intersects(intersectingSet, guards)
                        throw(ErrorException("Set intersecting guard does not intersect the guard."))
                    end
                    intersectedSet = getBoxIntersection(intersectingSet, guards)

                    if !isa(loc.invarient, Nothing) && intersects(intersectedSet, loc.invarient)
                        tintersectedSet = getBoxIntersection(intersectedSet, loc.invarient)
                        #println("Inv intersection: ", ρ(Vector(sparsevec([5], [1.0], 6)), intersectedSet), " vs ", ρ(Vector(sparsevec([5], [1.0], 6)), tintersectedSet))
                        intersectedSet = tintersectedSet
                        #println(intersectedSet)
                    end

                    (xx, _), _ = tempReachset
                    jumpSet = concretize(edge.jumpMatrix * intersectedSet)
                    #println(edge.jumpMatrix)

                    if !isa(hybridSystem.locations[edge.targetLoc].invarient, Nothing) && intersects(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)
                        tjumpSet = getBoxIntersection(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)# + edge.jumpVector * intersectedSet
                        #println("Jump intersection: ", ρ(Vector(sparsevec([5], [1.0], 6)), jumpSet), " vs ", ρ(Vector(sparsevec([5], [1.0], 6)), tjumpSet))
                        jumpSet = tjumpSet
                        #println(jumpSet)
                    end
                    #push!(reachset, ([(jumpSet, [reachtime, reachtime])], string(reachtime) * ": jump(" * string(loc.id) * ")->" * string(edge.targetLoc)))
                    #
                    #   Here we could optimize it such that in the case where guards ⊆ timeIntersectedSet we calculate both [supMins, endtime] and [infMaxs, endtime] with guards
                    #   and otherwise [supMins, endtime] with hyperplane intersection with timeIntersectedSet and [infMaxs, endtime] with guards intersection
                    #   Maybe look at how input should be handled... and if we can manipulate the constraints to account for the accumulated input
                    #

                    #timePointInput = U  #   NEEDS FIXING
                    branchedRun = auxReACTed(hybridSystem, hybridSystem.locations[edge.targetLoc], [reachtime, endtime], jumpSet, constraint, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder, tΦ)
                    reachset = vcat(reachset, branchedRun)
                    #reachset = vcat(reachset, nonintersectedSet) #Maybe gets the universe..

                    #branchedRun = auxReACTed(loc, [timeIntersected, endtime], nonintersectedSet, intersectedInput, constraint, δ⁻, δ⁺, flowPhiDict, alg, maxOrder, reduceOrder)
                    #push!(reachset, branchedRun)

                elseif reachtime - time < 0.0
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
                if !isa(loc.invarient, Nothing)
                    tempReachset = map((x, y) -> intersects(x, loc.invarient) ? (getBoxIntersection(concretize(x), loc.invarient), y) : (x, y), tempReachset)

                    #tempReachsetInInv = getBoxIntersection(concretize(tempReachset), loc.invarient)
                    #reachset = vcat(reachset, (tempReachset, string(loc.id) * "->" * string(edge.targetLoc)))
                else
                    #reachset = vcat(reachset, (tempReachset, string(loc.id) * "->" * string(edge.targetLoc)))
                end
            end
        else
            println("No guards? call ReACT")

            tempReachset, _, _ = ReACT(loc, δ⁻, δ⁺, interval, X0, U, constraint, 2, alg, maxOrder, reduceOrder, PhiDict[loc.id])
            if intersects(concretize(tempReachset), constraints_list(loc.invarient))

                tempReachsetInInv = getBoxIntersection(concretize(tempReachset), constraints_list(loc.invarient))
                reachset = vcat(reachset, (tempReachsetInInv, string(loc.id) * "->" * string(edge.targetLoc)))
            else
                reachset = vcat(reachset, (concretize(tempReachset), string(loc.id) * "->" * string(edge.targetLoc)))
            end

        end
    end

    return reachset
end

function ReACTTouches(loc, δ⁻::Float64, δ⁺::Float64, interval, guard, constraint, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, Φ)
    initialTimeStep = copy(δ⁺)
    m = copy(δ⁻)
    changedTimeStep = true
    phiDict = PhiDict

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    overapproximateIntersectingSetArray = []
    attemptsRecorder = []

    V = copy(inputDiscritezationDict[initialTimeStep])
    Sρ = zeros(Float64, length(constraint))
    newR = discritezationDict[initialTimeStep]
    i = 1


    if ismissing(Φ)
        Φ::Matrix{Float64} = diagm(ones(Float64, size(loc.A, 2)))
    end
    tempM = similar(Φ)
    ϕt = similar(Φ)
    newRR = copy(newR)


    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m
                if reduce(&, <=(Sρ + map(x -> ρ(x, newRR), constraintProjVectors), constraintProjBounds))
                    #throw(ErrorException("Reached unsafe set."))
                end
                bigCH = foldr((x, y) -> overapproximate(CH(x, y), Zonotope), overapproximateIntersectingSetArray; init=concretize(newRR))
                intersectingSet = overapproximate(bigCH, Zonotope)
                #println(norm(intersectingSet), " ", norm(newRR))
                return (intersectingSet, Sρ, time)
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

            if reduce(&, <=(Sρ + hom + inhom, constraintProjBounds)) && intersects(newRR, guard) #mapreduce(x -> intersects(newRR, x), &, guard)
                #if mapreduce(x -> intersects(newRR, x), &, guard)
                push!(overapproximateIntersectingSetArray, concretize(newRR))
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
    return (intersectingSet, Sρ, time)
end

function ReACTGuards(loc, δ⁻::Float64, δ⁺::Float64, interval, guards, constraint, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, Φ)
    initialTimeStep = copy(δ⁺)
    m = copy(δ⁻)
    changedTimeStep = true
    phiDict = PhiDict

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    #overapproximateIntersectingSetArray = []
    lastNewR = []
    attemptsRecorder = []


    if ismissing(Φ)
        Φ::Matrix{Float64} = diagm(ones(Float64, size(loc.A, 2)))
    end
    tempM = similar(Φ)
    ϕt = similar(Φ)

    V = copy(inputDiscritezationDict[initialTimeStep])
    Sρ = zeros(Float64, length(constraint))
    newR = discritezationDict[initialTimeStep]
    i = 1

    newRR = copy(newR)


    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m
                if ismissing(lastNewR) && intersects(newRR, guards)
                    push!(lastNewR, (concretize(newRR), [time, time + currentTimeStep]))
                elseif !reduce(&, <=(Sρ + map(x -> ρ(x, lastNewR), constraintProjVectors), constraintProjBounds))
                    throw(ErrorException("Reached unsafe set."))
                end
                #push!(overapproximateIntersectingSetArray, newRR)
                #intersectingSet = overapproximate(ConvexHullArray(overapproximateIntersectingSetArray), Zonotope)
                #println(norm(lastNewR), " ", time)
                return (lastNewR, Sρ, time, Φ)
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

            if reduce(&, <=(Sρ + hom + inhom, constraintProjBounds)) && !intersects(newRR, guards)
                #if mapreduce(x -> intersects(newRR, x), &, guard)
                #push!(overapproximateIntersectingSetArray, newRR)
                push!(lastNewR, (concretize(copy(newRR)), [time, time + currentTimeStep]))
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
    return (lastNewR, Sρ, time, Φ)
end

function ReACT(loc, δ⁻::Float64, δ⁺::Float64, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, STRATEGY::Integer, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, PhiDict=nothing) where {N}
    initialTimeStep = δ⁺
    changedTimeStep = true
    phiDict = PhiDict
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    inputDiscritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)

    discritezationDict, inputDiscritezationDict, phiDict = ReACTDiscretize(A, X0, U, δ⁻, δ⁺, alg, maxOrder, reduceOrder, phiDict)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    attemptsRecorder = Integer[]

    V::Zonotope{N,Vector{N},Matrix{N}} = copy(inputDiscritezationDict[initialTimeStep])
    Sρ = zeros(Float64, length(constraint))
    newR::Zonotope{N,Vector{N},Matrix{N}} = discritezationDict[initialTimeStep]
    i = 1


    Φ::Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    tempM = similar(Φ)
    ϕt = similar(Φ)
    newRR = copy(newR)


    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m
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

            if reduce(&, <=(Sρ + hom + inhom, constraintProjBounds))
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

    return (newR, sρ, time)
end
