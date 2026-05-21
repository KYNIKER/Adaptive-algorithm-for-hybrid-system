using LazySets, LinearAlgebra, ReachabilityAnalysis, HiGHS
using LazySets.Approximations: PolygonalOverapproximation, addapproximation!
include("Discretize.jl")
include("Utilities.jl")

model = JuMP.Model(HiGHS.Optimizer)
set_string_names_on_creation(model, false)
set_attribute(model, "presolve", "off")

#set_attribute(model, "eps_abs", 1e-5)
#set_attribute(model, "eps_rel", 1e-5)
set_silent(model)
const VERBOSE = false

function zonotopePrintDim(Z::Zonotope, dim)
    center = Z.center[dim]
    G = genmat(Z)
    _, genAmount = size(G)
    genContribute = sum(abs(G[dim, i]) for i in 1:genAmount)
    return "center: $center, generators: $genContribute"
end

function ReACTed(hybridSystem::HybridSystemV2, initialLoc, interval, X0, U, dirs, constraint, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, timeConstraintList=[]; clustering=false, mustSemantics=false)
    loc = initialLoc
    flowPhiDict = Dict(map(x -> x.id => PhiDict(x.A, δ⁻, δ⁺, alg), hybridSystem.locations))
    waitinglist = []
    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []
    saveResult = true
    dirsVectors = []

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

    push!(waitinglist, (loc, X0, nothing, interval, diagm(ones(dimLength)), nothing))
    while !isempty(waitinglist)
        location, initialset, edge, interval′, TΦ, polySet = pop!(waitinglist)

        reducedPolySet = polySet #isnothing(polySet) ? polySet : reducePolytope(polySet)

        #@show LazySets.API.high(initialset)
        #@show LazySets.API.low(initialset)

        #=if !isnothing(polySet)
            @show initialset ⊆ polySet
            @show LazySets.API.high(polySet) - LazySets.API.high(initialset)
            @show LazySets.API.low(polySet) - LazySets.API.low(initialset)
        end=#
        res = auxReACTed(waitinglist, hybridSystem, hybridSystem.locations[location], edge, interval′, initialset, dirsVectors, constraint, δ⁻, δ⁺, flowPhiDict, alg, maxOrder, reduceOrder, TΦ, reducedPolySet, timeConstraintList, saveResult; clustering, mustSemantics)
        reachset = vcat(reachset, res)
    end

    return reachset
end


# AucReacted is called recursively each time we have a new starting location (after a transition)
function auxReACTed(waitlist, hybridSystem, loc::Location, currentEdge, interval, X0, dirs, constraint, δ⁻::Float64, δ⁺::Float64, PhiDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, Φ=missing, polyhedralSet=nothing, timeConstraintList=[], saveResult::Bool=true; clustering=false, mustSemantics=false)
    activeTimeConstraints = []
    #@show norm(X0)
    for (id, time) in timeConstraintList
        if id == loc.id
            push!(activeTimeConstraints, time)
        end
    end

    dims = size(X0.center, 1)
    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []
    setOfConstraints = vcat(loc.constraints, constraint)

    discretizationDict, overapproximatedDiscretizationDict, inputDiscritezationDict = ReACTDiscretizePlus(loc, X0, polyhedralSet, δ⁻, δ⁺, alg, maxOrder, reduceOrder, PhiDict[loc.id])
    #=
    discretizationDict, inputDiscritezationDict = ReACTDiscretizePlus(loc, X0, δ⁻, δ⁺, alg, maxOrder, reduceOrder, PhiDict[loc.id])

    #@show norm(discretizationDict[δ⁻])




    #t1 = Base.time()
    overapproximatedDiscretizationDict = Dict()

    intersectedSetIsNothing = isnothing(polyhedralSet)
    if !intersectedSetIsNothing
        overapproximatedDiscretizationDict, _ = ReACTDiscretizePlus(loc, polyhedralSet, δ⁻, δ⁺, alg, maxOrder, reduceOrder, PhiDict[loc.id])
    else
        for key in keys(discretizationDict)
            overapproximatedDiscretizationDict[key] = overapproximate(discretizationDict[key], BoxDirections(dims))
        end
    end
    #@show (norm(discretizationDict[δ⁻]), norm(overapproximatedDiscretizationDict[δ⁻]))
    #@show time

    #println("Finished Dicts $(loc.id) $(length(loc.edges)) in $(Base.time() - t1) s")
    =#



    # For each edge we simulate the system
    listOfEdges = loc.edges
    if !isa(currentEdge, Nothing) # This allows us to sometimes just wanna do a signular edge
        listOfEdges = [currentEdge]
    end

    for edge in listOfEdges
        Reset_Map(X) = bloatPolytope(Singleton(edge.jumpVector), edge.jumpMatrix, X)
        # println("Handling edge at time $time from $(loc.id) -> $(edge.targetLoc)")
        guards = edge.guard # Technically the guard is one singular HPolyhedron, but it composes the other guards


        #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 

        tempReachset, reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, [time, endtime], guards, setOfConstraints, dirs, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, overapproximatedDiscretizationDict, saveResult)





        if reachtime - time == 0.0
            if !isempty(reachset)
                # println("Found an immediate transition to $(edge.targetLoc), but we are not taking it as we are scared of zeno behaviour")
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
            # println("Current time: $time, intersecting time start: $reachtime")
            tryContinueFlag, _, timeNotIntersected, intersectingSetsList = ReACTTouches(loc, δ⁻, δ⁺, [reachtime, endtime], (reachtime - time), guards, setOfConstraints, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, overapproximatedDiscretizationDict, tΦ, nothing, reduceOrder, maxOrder)


            if any((timeNotIntersected >= x) for x in activeTimeConstraints)
                handleHitConstraint(timeNotIntersected, loc.id)
            end


            timeIntersected = timeNotIntersected - reachtime
            latestSet = nothing
            #
            #   Here we should check whether we have reached endtime. If true we should only push the jumpSet
            #   Still need to check whether we have reached the invariant. If true we should NOT push the else branch result, only the tempReachsets[infMaxsIdx]
            #
            invariantGuardIntersection = LazySets.intersection(loc.invarient, guards)
            if !isempty(intersectingSetsList)
                # we first save this for later
                if tryContinueFlag
                    latestSet = copy(last(intersectingSetsList))
                end
                # Afterwards we process the intersecting set

                intersectedSet = nothing
                plottingList = []
                jumpSetsList = Vector{Zonotope}()
                jumpPolySetsList = Vector{HPolytope}()
                i = 0

                for set in intersectingSetsList
                    #polySet = overapproximate(set, BoxDirections(dims))

                    if !isa(guards, Nothing)
                        if !isDisjointFast(set, guards) #!isdisjoint(guards, set)
                            set = zonotopeStripIntersection(set, guards)
                        else
                            #@show LazySets.API.high(set)
                            continue
                        end
                    end

                    # Check invarient
                    if !isa(loc.invarient, Nothing)
                        if !isDisjointFast(set, loc.invarient)
                            set = zonotopeStripIntersection(set, loc.invarient)
                        else
                            continue
                        end
                    end

                    polySet = LazySets.intersection(overapproximate(set, BoxDirections(dims)), invariantGuardIntersection; prune=false)
                    # Push to reachset
                    if saveResult
                        push!(plottingList, (map(x -> ρ(x, polySet), dirs), [reachtime + δ⁻ * i, reachtime + δ⁻ * (i + 1)]))
                        i += 1
                    end

                    # Apply jump Matrix
                    jumpSet = linear_map(edge.jumpMatrix, set)
                    LazySets.translate!(jumpSet, edge.jumpVector)
                    if !isa(hybridSystem.locations[edge.targetLoc].invarient, Nothing)
                        if !isDisjointFast(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)
                            jumpSet = zonotopeStripIntersection(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)
                        else
                            continue
                        end
                    end
                    push!(jumpSetsList, jumpSet)
                    jumpPolySet = Reset_Map(polySet)
                    push!(jumpPolySetsList, jumpPolySet)
                end


                if saveResult
                    push!(reachset, (plottingList, string(reachtime) * " - " * string(timeNotIntersected) * ": " * string(loc.id) * "->" * string(loc.id)))
                end

                len = length(jumpSetsList)
                if len == 0
                    continue
                end
                if clustering
                    intersectedSet = nothing
                    polyIntersectedSet = nothing
                    if len == 1
                        intersectedSet = jumpSetsList[1]
                        polyIntersectedSet = jumpPolySetsList[1]
                    else
                        tempIntersect = foldl(ConvexHull, jumpSetsList)
                        tempPolyIntersect = foldl(overapproximatedCH, jumpPolySetsList)
                        intersectedSet = convert(Zonotope, box_approximation(tempIntersect))
                        polyIntersectedSet = tempPolyIntersect
                    end

                    jumpSet = intersectedSet
                    # println("Finished intersections")
                    push!(waitlist, (edge.targetLoc, copy(jumpSet), nothing, [reachtime, endtime], tΦ, polyIntersectedSet))
                    #branchedRun = auxReACTed(hybridSystem, hybridSystem.locations[edge.targetLoc], nothing, [reachtime, endtime], jumpSet, dirs, constraint, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder, tΦ, clustering, timeConstraintList, saveResult)

                    #=if saveResult
                        reachset = vcat(reachset, branchedRun)
                    end=#
                else # No clustering. Do individual runningset
                    for (index, jumpSet) in enumerate(jumpSetsList)
                        # Note that we only do steps of size δ⁻ in touches
                        timeStart = reachtime + (index - 1) * δ⁻
                        push!(waitlist, (edge.targetLoc, copy(jumpSet), nothing, [timeStart, endtime], tΦ, jumpPolySetsList[index]))
                        #=branchedRun = auxReACTed(hybridSystem, hybridSystem.locations[edge.targetLoc], nothing, [timeStart, endtime], jumpSet, dirs, constraint, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder, tΦ, clustering, timeConstraintList, saveResult)

                        if saveResult
                            reachset = vcat(reachset, branchedRun)
                        end=#
                    end
                end
            else
                tryContinueFlag = false
                # println("No intersections with guards, end branch")

            end

            # If we are not encountering an invarient, try continue
            if tryContinueFlag
                # println("Continuing from previous run at time $timeNotIntersected")
                push!(waitlist, (loc.id, copy(latestSet), edge, [timeNotIntersected, endtime], tΦ, nothing))
                #=branchedRun = auxReACTed(hybridSystem, loc, edge, [timeNotIntersected, endtime], latestSet, dirs, constraint, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder, tΦ, clustering, timeConstraintList, saveResult)
                if saveResult
                    reachset = vcat(reachset, branchedRun)
                end=#
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
                    reachset = vcat(reachset, (tempReachset, string(loc.id) * "->" * string(edge.targetLoc)))
                else
                    reachset = vcat(reachset, (tempReachset, string(loc.id) * "->" * string(edge.targetLoc)))
                end
            end
        end
    end

    if isempty(loc.edges) # This means it is just a continous system from here
        # println("No edges? call ReACT")

        tempReachset, reachtime, _ = ReACT(loc, δ⁻, δ⁺, [time, endtime], setOfConstraints, dirs, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, saveResult)

        if saveResult
            reachset = vcat(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(loc.id)))
        end
    end


    return reachset
end


# Optimized
function ReACTTouches(loc, δ⁻::Float64, δ⁺::Float64, interval, initialTime::Float64, guard, constraint, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, polytopeDict, Φ, accInput, reduce_order, max_order)
    STRATEGY = 0
    initialTimeStep = copy(δ⁻)
    #initialTimeStep = copy(δ⁺)
    m = copy(δ⁻)
    #initialTimeStep = copy(δ⁺)
    changedTimeStep = true
    phiDict = PhiDict

    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(guard)
    invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)
    #dirProjVectors = map(x -> x, dirs)

    #originalConstraintProjVectors = copy(constraintProjVectors)
    #originalGuardProjVectors = copy(guardProjVectors)
    #originalInvarientProjVectors = copy(invarientProjVectors)

    #permutedphiDict = Dict()
    #for key in keys(phiDict)
    #    permutedphiDict[key] = permutedims(phiDict[key])
    #end


    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    V = copy(inputDiscritezationDict[m])
    V = linear_map(Φ, V)
    Vs = nestedInputDiscCalculate(inputDiscritezationDict, PhiDict, δ⁺, δ⁻, initialTime, reduce_order, max_order)


    #concretize(Vs)
    #lastVs = copy(Vs)
    #Sρ = map(x -> ρ(x, Vs), constraintProjVectors)
    #Gρ = map(x -> ρ(x, Vs), guardProjVectors)
    #Iρ = map(x -> ρ(x, Vs), invarientProjVectors)
    #dρ = map(x -> ρ(x, Vs), dirs)
    newR = copy(discritezationDict[initialTimeStep])
    polyNewR = polytopeDict[initialTimeStep]
    i = 1
    if ismissing(Φ)
        Φ::Matrix{Float64} = exp((initialTime) .* loc.A)
    end

    tempM = diagm(ones(Float64, size(loc.A, 2)))
    ϕt = similar(Φ)
    #pϕt = similar(tempM)
    newRR = linear_map(Φ, newR)

    attemptsRecorder = []
    intersectingSetsList = Vector{Zonotope}()

    #transposedΦ = permutedims(Φ)
    #constraintProjVectors = map(x -> transposedΦ * x, originalConstraintProjVectors)
    #guardProjVectors = map(x -> transposedΦ * x, originalGuardProjVectors)
    #invarientProjVectors = map(x -> transposedΦ * x, originalInvarientProjVectors)
    #dirProjVectors = map(x -> x, dirs)

    #oldDirProjVectors = copy(dirProjVectors)
    #oldConstraintProjVectors = copy(constraintProjVectors)
    #oldGuardProjVectors = copy(guardProjVectors)
    #oldInvarientProjVectors = copy(invarientProjVectors)


    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m

                tempSet = concretize(newRR ⊕ Vs)
                if any(((ρ(x, tempSet)) > y) for (x, y) in zip(constraintProjVectors, constraintProjBounds))
                    #throw(ErrorException("Reached unsafe set."))
                    handleHitConstraint(time, loc.id)
                end

                # Assuming we stopped because we no longer intersect guards
                continueAfter = true

                if all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds))
                    # println("Intersecting guard, but no longer intersecting invarients")
                    continueAfter = false

                    push!(intersectingSetsList, tempSet)
                end

                return continueAfter, Φ, time, intersectingSetsList


                #println("Touches Vs: ", accInput)
                #println(norm(intersectingSet), " ", norm(newRR))
            end

            if changedTimeStep
                newR = discritezationDict[currentTimeStep]
                V = copy(inputDiscritezationDict[currentTimeStep])
                ϕt = phiDict[currentTimeStep]
                #pϕt = permutedphiDict[currentTimeStep]
                newRR = linear_map(Φ, newR)
                V = linear_map(Φ, V)
            else
                newRR = linear_map(ϕt, newRR)
                V = linear_map(ϕt, V)
            end
            changedTimeStep = false

            # @show concretize(newRR)
            # #@show concretize(Vs)
            # @show concretize(V)

            tempSet = minkowski_sum(newRR, Vs)
            if all((ρ(x, tempSet)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) && # IsSubSet
               all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds)) &&
               all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Intersects

                #if all(input + ρ(x, newR) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) && # IsSubSet
                #   all(((input + -ρ(-x, newR)) <= y) for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) &&  #!isdisjoint(tempSet, guard; algorithm="sufficient") && # intersects
                #   all((input + -ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds)) # Intersects
                #if touchesCheck(newR, Sρ, constraintProjVectors, constraintProjBounds, Gρ, guardProjVectors, guardProjBounds, Iρ, invarientProjVectors, invarientProjBounds)
                #if mapreduce(x -> intersects(newRR, x), &, guard)

                push!(intersectingSetsList, copy(minkowski_sum(newRR, Vs)))
                #lastVs = copy(Vs)
                #Vs = Vs ⊕ V
                # Vs = concretize(Vs)
                # V = concretize(V) 
                V = linear_map(Φ, V)
                Vs = minkowski_sum(Vs, V)

                #Vs = LinearMap(ReachabilityAnalysis.Exponentiation.Φ₁(loc.A, time + currentTimeStep - minimum(interval), ReachabilityAnalysis.Exponentiation.BaseExp), inputDiscritezationDict[0])
                approveFlag = true
                #triedRevise = false
                #Sρ += inhom


                #Sρ += map(x -> ρ(x, V), originalConstraintProjVectors)
                #Gρ += map(x -> ρ(x, V), originalGuardProjVectors)
                #Iρ += map(x -> ρ(x, V), originalInvarientProjVectors)
                #constraintProjVectors = map(x -> pϕt * x, oldConstraintProjVectors)
                #guardProjVectors = map(x -> pϕt * x, oldGuardProjVectors)
                #invarientProjVectors = map(x -> pϕt * x, oldInvarientProjVectors)

                #oldConstraintProjVectors = constraintProjVectors
                #oldGuardProjVectors = guardProjVectors
                #oldInvarientProjVectors = invarientProjVectors


                mul!(tempM, Φ, ϕt)
                copy!(Φ, tempM)
            else
                #newR = copy(newR)
                currentTimeStep = currentTimeStep / 2
                changedTimeStep = true
                attempts = attempts + 1
                #triedRevise = false
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
    continueAfter = false # We are at the end time horizon, therefore no continuing
    return (continueAfter, Φ, time, intersectingSetsList)
end


#=
function ReACTTouches(loc, δ⁻::Float64, δ⁺::Float64, interval, initialTime::Float64, guard, constraint, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, polytopeDict, Φ, accInput, reduce_order, max_order)
    # Note that in touches we always use δ⁻
    # That is, we do not adjust timestep sizes

    phiDict = PhiDict

    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(guard)
    invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)


    Vs = nestedInputDiscCalculate(inputDiscritezationDict, PhiDict, δ⁺, δ⁻, initialTime, reduce_order, max_order)
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
        tempSet = concretize(newRR ⊕ Vs)

        if all((ρ(x, tempSet)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) && # IsSubSet
           all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds)) &&
           all((-ρ(-x, tempSet)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # Intersects
            #if mapreduce(x -> intersects(newRR, x), &, guard)

            push!(intersectingSetsList, copy(tempSet))

            # Main calculation. No longer changing timesteps
            Vs = minkowski_sum(Vs, V) #Update input
            newRR = linear_map(ϕ, newRR)
            V = linear_map(ϕ, V)

            # i = i + 1
            time = time + δ⁻

        else # Handle hit something. We do not reduce anymore!
            if any(((ρ(x, tempSet)) > y) for (x, y) in zip(constraintProjVectors, constraintProjBounds))
                handleHitConstraint(time, loc.id)
            end

            continueAfter = true
            if all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds))
                # If we stop because we are no longer intersect guards, 
                # but still intersect the invariant we try continue
                continueAfter = false
            end

            return continueAfter, Φ, time, intersectingSetsList
        end
    end
    continueAfter = false # We are at the end time horizon, therefore no continuing
    return (continueAfter, Φ, time, intersectingSetsList)
end
=#

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



function ReACTGuards(loc, δ⁻::Float64, δ⁺::Float64, interval, guards, constraint, dirs, STRATEGY::Integer, PhiDict, discritezationDict, inputDiscritezationDict, polytopeDict, saveResult)
    # We calculate the reachset till we reach a guard for an intersection (or till failure)
    initialTimeStep = copy(δ⁺)
    m = copy(δ⁻)
    changedTimeStep = true
    phiDict = PhiDict
    A = copy(loc.A)
    # constraintProjVectors = map(x -> x.a, constraint)
    # constraintProjBounds = ρ.(constraintProjVectors, constraint)
    #=constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    guardProjVectors, guardProjBounds = [], []
    if isnothing(guards)
        guardProjVectors, guardProjBounds = [constraintProjVectors[1]], [Inf]
    else
        guardProjVectors, guardProjBounds = getHalfSpaceProjections(guards)
    end

    invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)
    dirProjVectors = map(x -> x, dirs)=#
    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(guards)
    invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)
    dirProjVectors = map(x -> x, dirs)
    #=
    oldDirProjVectors = copy(dirProjVectors)
    oldConstraintProjVectors = copy(constraintProjVectors)
    oldGuardProjVectors = copy(guardProjVectors)
    oldInvarientProjVectors = copy(invarientProjVectors)
    =#
    permutedphiDict = Dict{Float64,Matrix{Float64}}()
    for key in keys(phiDict)
        permutedphiDict[key] = permutedims(phiDict[key], (2, 1))
    end


    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(δ⁻)#copy(initialTimeStep)

    #overapproximateIntersectingSetArray = []
    #lastNewR = []
    attemptsRecorder = []

    dirVals = [] # This is where the plotting happens


    Φ::Matrix{Float64} = diagm(ones(Float64, size(loc.A, 2)))

    tempM = diagm(ones(Float64, size(loc.A, 2)))
    ϕt = similar(tempM)
    pϕt = similar(tempM)
    V::Zonotope = copy(inputDiscritezationDict[initialTimeStep])

    #lastVs = copy(Vs)
    Sρ = zeros(Float64, length(constraint))
    Gρ = zeros(Float64, length(guardProjVectors))
    Iρ = zeros(Float64, length(invarientProjVectors))
    dρ = zeros(Float64, length(dirs))
    newR = copy(discritezationDict[initialTimeStep])
    polyNewR = polytopeDict[initialTimeStep]
    i = 1

    U::Zonotope = inputDiscritezationDict[0]

    #triedRevise = false
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
                if (any((input + ρ(x, newR)) > y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) && any((input + ρ(x, polyNewR; solver=model)) > y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)))
                    #@show [((input + ρ(x, newR)), y) for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)]
                    #@show [((input + ρ(x, polyNewR; solver=model)), y) for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)]

                    #@show Sρ
                    handleHitConstraint(time, loc.id)
                end
                #println("Guards i: $i")
                return (dirVals, time, Φ)
            end

            if changedTimeStep
                newR = discritezationDict[currentTimeStep]
                polyNewR = polytopeDict[currentTimeStep]
                V = copy(inputDiscritezationDict[currentTimeStep])
                ϕt = phiDict[currentTimeStep]
                pϕt = permutedphiDict[currentTimeStep]

            end

            changedTimeStep = false

            # #println(any((input + ρ(-x, newRR)) > y for (input, x, y) in zip(Sρ, invarientProjVectors, invarientProjBounds)))
            #println(all((input + ρ(x, newRR)) <= y for (input, x, y) in zip(Sρ, invarientProjVectors, invarientProjBounds)))
            #if (all(input + ρ(x, newR) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) || all((input + ρ(x, polyNewR)) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))) &&
            #   (any((input + -ρ(-x, newR)) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) || any((input + -ρ(-x, polyNewR)) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds))) &&
            #   (all((input + -ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds)) || all((input + -ρ(-x, polyNewR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds)))
            if guardCheck(newR, polyNewR, Sρ, constraintProjVectors, constraintProjBounds, Gρ, guardProjVectors, guardProjBounds, Iρ, invarientProjVectors, invarientProjBounds)
                if saveResult
                    plotDirZ = copy(dρ + map(x -> ρ(x, newR), dirProjVectors))
                    plotDirP = copy(dρ + map(x -> ρ(x, polyNewR; solver=model), dirProjVectors))
                    plotDir = [abs(x) < abs(y) ? x : y for (x, y) in zip(plotDirZ, plotDirP)]
                    push!(dirVals, (copy(plotDir), [time, time + currentTimeStep]))
                end
                #lastVs = copy(Vs)
                #Vs = ReachabilityAnalysis.Exponentiation.Φ₁(A, time - minimum(interval), ReachabilityAnalysis.Exponentiation.BaseExp) * U

                approveFlag = true

                dρ += map(x -> ρ(x, V), dirProjVectors)
                Sρ += map(x -> ρ(x, V), constraintProjVectors)
                Gρ += map(x -> ρ(x, V), guardProjVectors)
                Iρ += map(x -> ρ(x, V), invarientProjVectors)
                map!(x -> pϕt * x, dirProjVectors)
                map!(x -> pϕt * x, constraintProjVectors)
                #guardProjVectors = map(x -> pϕt * x, guardProjVectors)
                map!(x -> pϕt * x, guardProjVectors)
                #invarientProjVectors = map(x -> pϕt * x, invarientProjVectors)
                map!(x -> pϕt * x, invarientProjVectors)

                #=
                oldDirProjVectors = dirProjVectors
                oldConstraintProjVectors = constraintProjVectors
                oldGuardProjVectors = guardProjVectors
                oldInvarientProjVectors = invarientProjVectors
                =#
                mul!(tempM, Φ, ϕt)
                copy!(Φ, tempM)
            else
                #newR = copy(newR)
                currentTimeStep = currentTimeStep / 2
                changedTimeStep = true
                attempts = attempts + 1
                #triedRevise = false
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
    invarientProjVectors, invarientProjBounds = [], []
    if isnothing(loc.invarient)
        invarientProjVectors, invarientProjBounds = [constraintProjVectors[1]], [Inf]
    else
        invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient)

    end
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
    V = inputDiscritezationDict[initialTimeStep]

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
               all((input + ρ(x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))
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

function guardCheck(newR::Zonotope, polyNewR, Sρ, constraintProjVectors, constraintProjBounds, Gρ, guardProjVectors, guardProjBounds, Iρ, invarientProjVectors, invarientProjBounds) #; solver=model
    #cache = map((input, x, y) -> input + ρ(x, newR) <= y, Sρ, constraintProjVectors, constraintProjBounds)
    abssum = 0.0
    c = newR.center
    G = transpose(genmat(newR))
    tc = zeros(size(G, 1))
    #a = sum(abs, transpose(a) * G)
    res = all(input + tsupfunc(x, tc, abssum, c, G) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) || all((input + ρ(x, polyNewR)) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))
    res = res && (any((input + -tsupfunc(-x, tc, abssum, c, G)) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) || any((input + -ρ(-x, polyNewR)) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)))
    res = res && (all((input + -tsupfunc(-x, tc, abssum, c, G)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds)) || all((input + -ρ(-x, polyNewR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds)))
    return res
end

function tsupfunc(d, Z, Ab, c, G)
    #c = center(Z)
    #G = genmat(Z)
    mul!(Z, G, d)
    Ab = sum(abs, Z)
    #c = dot(center(Z), d)
    return dot(c, d) + Ab
end

function touchesCheck(newR, Sρ, constraintProjVectors, constraintProjBounds, Gρ, guardProjVectors, guardProjBounds, Iρ, invarientProjVectors, invarientProjBounds)
    return all(input + ρ(x, newR) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) && # IsSubSet
           all(((input + -ρ(-x, newR)) <= y) for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) &&
           all((input + -ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds)) # Intersects
end