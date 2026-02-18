using LazySets, LinearAlgebra
include("Discretize.jl")
include("Utilities.jl")


function ReACTed(hybridSystem::HybridSystemV2, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, δ⁻, δ⁺, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    flowPhiDict = Dict{eltype(hybridSystem.locations),Matrix{N}}(map(x -> x => PhiDict(x.A, δ⁻, δ⁺, alg), hybridSystem.locations))
    loc = hybridSystem.initialLoc

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
    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []
    for edge in loc.edges
        guards = edge.guard
        if !ismissing(guards)
            mins = []
            maxs = []
            tempReachsets = []
            tempInputs = []
            for guard in constraints_list(guards)
                let tempReachset, tempInput, reachtime = ReACT(loc, δ⁻, δ⁺, [time, endtime], X0, U, guard, 2, alg, maxOrder, reduceOrder, PhiDict)
                    if reachtime < endtime
                        push!(mins, reachtime)
                        timeIntersected = ReACTTouches(loc, δ⁻, δ⁺, [reachtime, endtime], X0, U, guard, 2, maxOrder, reduceOrder, flowPhiDict[loc])
                        #=
                        if isnothing(timeIntersected) # Handle the case where the unsafe set is reached while intersecting the guard
                            return missing
                        end
                        =#
                        push!(maxs, reachtime + timeIntersected)
                    end
                    push!(mins, endtime)
                    push!(maxs, endtime)
                    push!(tempReachsets, tempReachset)
                    push!(tempInputs, tempInput)
                end
            end
            supMins, supMinsIdx = findmax(mins)
            infMaxs, infMaxsIdx = findmin(maxs)
            timeIntersected = infMaxs - supMins

            #
            #   Here we should check whether we have reached endtime. If true we should only push the jumpSet
            #   Still need to check whether we have reached the invariant. If true we should NOT push the else branch result, only the tempReachsets[infMaxsIdx]
            #

            if timeIntersected >= 0
                timeIntersectedSet = concretize(overapproximateIntervalReachset(loc.A, tempReachsets[supMinsIdx], U, δ⁻, timeIntersected, alg, maxOrder, reduceOrder, flowPhiDict[loc.id]))
                nonintersectedSet, intersectedSet = intersection(timeIntersectedSet, guard)
                jumpSet = edge.jumpMatrix * intersectedSet + edge.jumpVector * intersectedSet

                branchedRun = auxReACTed(edge.targetLoc, [supMins, endtime], jumpset, tempInputs[supMinsIdx], constraint, δ⁻, δ⁺, flowPhiDict[loc.id], alg, maxOrder, reduceOrder)
                push!(reachset, branchedRun)
            else
                branchedRun = auxReACTed(loc, [infMaxs, endtime], tempReachsets[infMaxsIdx], tempInputs[infMaxsIdx], constraint, δ⁻, δ⁺, flowPhiDict[loc.id], alg, maxOrder, reduceOrder)
                push!(reachset, branchedRun)
                #=
                if ρ(timeIntersectedSet.center, guard.a) >= guard.b

                =#
            end
        end
    end
    return reachset
end

function ReACTTouches(loc, δ⁻, δ⁺, [reachtime, endtime], X0, U, guard, 2, maxOrder, reduceOrder, flowPhiDict[loc])
    return missing
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

function ReACT(A, B, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Nothing, constraint, Digits::Integer, STRATEGY::Integer, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    XG = copy(genmat(X0))
    XDim, _ = size(XG)
    m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)   #Calculate the smallest number larger than 10^-Digits obtained by repeatedly dividing initialTimeStep by 2.
    changedTimeStep = true

    elems = (ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)
    phiDict = Dict{Float64,Matrix{Float64}}()
    sizehint!(phiDict, elems)
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    sizehint!(discritezationDict, elems)

    discritezationDict, _, phiDict = ReACTDiscretize(A, B, X0, U, m, initialTimeStep, alg, maxOrder, reduceOrder)

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    attemptsRecorder = Integer[]

    newR::Zonotope{N,Vector{N},Matrix{N}} = discritezationDict[initialTimeStep]
    i = 1


    Φ::Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    tempM = similar(Φ)
    ϕt = similar(Φ)
    newRR = copy(newR)
    initialϕ = phiDict[initialTimeStep]

    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m
                println(time)
                return false
            end

            if changedTimeStep
                newR = discritezationDict[currentTimeStep]
                ϕt = phiDict[currentTimeStep]
                newRR = linear_map(Φ, newR)
            else
                newRR = linear_map(ϕt, newRR)
            end

            changedTimeStep = false

            hom = map(x -> ρ(x, newRR), constraintProjVectors)
            if reduce(&, <=(hom, constraintProjBounds))
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

    return true
end