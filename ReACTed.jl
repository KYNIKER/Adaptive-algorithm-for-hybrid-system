using LazySets, LinearAlgebra
include("Discretize.jl")
include("Utilities.jl")


function ReACTed(hybridSystem::HybridSystemV2, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, δ⁻, δ⁺, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    loc = hybridSystem.initialLoc
    flowPhiDict = Dict{eltype(hybridSystem.locations),Matrix{N}}(map(x -> x => PhiDict(x.A, δ⁻, δ⁺, alg), hybridSystem.locations))




    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []

    res = auxReACTed(loc, interval, X0, U, constraint, δ⁻, δ⁺, flowPhiDict[loc.id], alg, maxOrder, reduceOrder)

    #=for j in eachindex(hybridSystem.G[loc, :])
        guards = hybridSystem.G[loc, j]
        if !ismissing(guards)
            for guard in guards
                let tempReachset, tempInput, reachtime = ReACT(hybridSystem.Flow[loc], δ⁻, δ⁺, [0.0, endtime], X0, U, guard, 2, alg, maxOrder, reduceOrder, flowPhiDict[loc])
                    if reachtime < endtime
                        #=
                        timeIntersected = ReACTTouches(hybridSystem.Flow[loc], δ⁻, δ⁺, [reachtime, endtime], X0, U, guard, 2, maxOrder, reduceOrder, flowPhiDict[loc])
                        timeIntersectedSet = concretize(ReACTDiscretize(hybridSystem.Flow[loc], tempReachset, U, δ⁻, timeIntersected, alg, maxOrder, reduceOrder, flowPhiDict[loc])[timeIntersected])
                        if ρ(timeIntersectedSet.center, guard.a) >= guard.b
                            jumpSet = edge.jumpMatrix * timeIntersectedSet + edge.jumpVector * timeIntersectedSet
                            ReACT(hybridSystem.Flow[edge.targetLoc], δ⁻, δ⁺, [0.0, endtime], X0, U, guard, 2, alg, maxOrder, reduceOrder, flowPhiDict[loc])
                        nonintersectedSet, intersectedSet = intersection(timeIntersectedSet, guard)
                        newLoc = j
                        # New run of React 
                        =#
                    end
                end
            end
        end
    end=#
end

function auxReACTed(loc::Location, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, δ⁻, δ⁺, PhiDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    discretizationDict, inputDiscritezationDict = ReACTDiscretize(loc.A, X0, U, δ⁻, δ⁺, alg, maxOrder, reduceOrder, flowPhiDict[loc.id])
    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []
    for edge in loc.edges
        guards = edge.guard
        if !ismissing(guards)


            tempReachset, tempInput, reachtime = ReACTGuards(loc, δ⁻, δ⁺, [time, endtime], guard, constraint, 2, PhiDict, discretizationDict, inputDiscritezationDict)
            if reachtime < endtime
                intersectingSet, intersectedInput, timeIntersected = ReACTTouches(loc, δ⁻, δ⁺, [reachtime, endtime], guard, constraint, 2, PhiDict, discretizationDict, inputDiscritezationDict)
                timeIntersected = timeIntersected - reachtime

                #
                #   Here we should check whether we have reached endtime. If true we should only push the jumpSet
                #   Still need to check whether we have reached the invariant. If true we should NOT push the else branch result, only the tempReachsets[infMaxsIdx]
                #

                if timeIntersected >= 0
                    #timeIntersectedSet = concretize(overapproximateIntervalReachset(loc.A, tempReachsets[supMinsIdx], U, δ⁻, timeIntersected, alg, maxOrder, reduceOrder, flowPhiDict[loc.id]))
                    nonintersectedSet, intersectedSet = intersection(intersectingSet, guards)
                    jumpSet = edge.jumpMatrix * intersectedSet + edge.jumpVector * intersectedSet

                    #
                    #   Here we could optimize it such that in the case where guards ⊆ timeIntersectedSet we calculate both [supMins, endtime] and [infMaxs, endtime] with guards
                    #   and otherwise [supMins, endtime] with hyperplane intersection with timeIntersectedSet and [infMaxs, endtime] with guards intersection
                    #   Maybe look at how input should be handled... and if we can manipulate the constraints to account for the accumulated input
                    #

                    timePointInput = missing
                    branchedRun = auxReACTed(edge.targetLoc, [reachtime, endtime], jumpset, timePointInput, constraint, δ⁻, δ⁺, flowPhiDict[loc.id], alg, maxOrder, reduceOrder)
                    push!(reachset, branchedRun)
                    branchedRun = auxReACTed(loc, [timeIntersected, endtime], nonintersectedSet, intersectedInput, constraint, δ⁻, δ⁺, flowPhiDict[loc.id], alg, maxOrder, reduceOrder)

                else
                    branchedRun = auxReACTed(loc, [infMaxs, endtime], tempReachsets[infMaxsIdx], tempInputs[infMaxsIdx], constraint, δ⁻, δ⁺, flowPhiDict[loc.id], alg, maxOrder, reduceOrder)
                    push!(reachset, branchedRun)
                    #=
                    if ρ(timeIntersectedSet.center, guard.a) >= guard.b

                        =#
                end
            else
                push!(reachset, tempReachset)
            end
        end
    end
    return reachset
end

function ReACTTouches(loc, δ⁻, δ⁺, interval, guard, constraint, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict)
    initialTimeStep = δ⁺
    m = δ⁻
    changedTimeStep = true
    phiDict = PhiDict

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    overapproximateIntersectingSetArray = []
    attemptsRecorder = []

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
                if reduce(&, <=(Sρ + map(x -> ρ(x, newRR), constraintProjVectors), constraintProjBounds))
                    throw(ErrorException("Reached unsafe set."))
                end
                push!(overapproximateIntersectingSetArray, newRR)
                intersectingSet = overapproximate(ConvexHullArray(overapproximateIntersectingSetArray), Zonotope)
                return (intersectingSet, sρ, time)
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

            if reduce(&, <=(Sρ + hom + inhom, constraintProjBounds)) && mapreduce(x -> intersects(newRR, x), &, guard)
                #if mapreduce(x -> intersects(newRR, x), &, guard)
                push!(overapproximateIntersectingSetArray, newRR)
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
    intersectingSet = overapproximate(ConvexHullArray(overapproximateIntersectingSetArray), Zonotope)
    return (intersectingSet, sρ, time)
end

function ReACTGuards(loc, δ⁻, δ⁺, interval, guards, constraint, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict)
    initialTimeStep = δ⁺
    m = δ⁻
    changedTimeStep = true
    phiDict = PhiDict

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    #overapproximateIntersectingSetArray = []
    lastNewR = missing
    attemptsRecorder = []

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
                if reduce(&, <=(Sρ + map(x -> ρ(x, lastNewR), constraintProjVectors), constraintProjBounds))
                    throw(ErrorException("Reached unsafe set."))
                end
                #push!(overapproximateIntersectingSetArray, newRR)
                #intersectingSet = overapproximate(ConvexHullArray(overapproximateIntersectingSetArray), Zonotope)
                return (lastNewR, sρ, time)
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

            if reduce(&, <=(Sρ + hom + inhom, constraintProjBounds)) && !(mapreduce(x -> intersects(newRR, x), &, guards))
                #if mapreduce(x -> intersects(newRR, x), &, guard)
                #push!(overapproximateIntersectingSetArray, newRR)
                lastNewR = newRR
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
    intersectingSet = overapproximate(ConvexHullArray(overapproximateIntersectingSetArray), Zonotope)
    return (intersectingSet, sρ, time)
end

function ReACT(loc, δ⁻, δ⁺, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, STRATEGY::Integer, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, PhiDict=nothing) where {N}
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
