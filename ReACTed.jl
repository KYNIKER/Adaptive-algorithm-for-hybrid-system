using LazySets, LinearAlgebra, ReachabilityAnalysis
using LazySets.Approximations: PolygonalOverapproximation, addapproximation!
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

function ReACTed(hybridSystem::HybridSystemV2, initialLoc, interval, X0, U, dirs, constraint, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) 
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
    if VERBOSE
        @show dirsVectors
        
    end

    push!(waitinglist, (loc, X0, interval, diagm(ones(dimLength)), nothing, nothing, nothing))
    while !isempty(waitinglist)
        location, initialset, interval′, TΦ, discDict, inputDict, nonIntersectedDiscDict = pop!(waitinglist)

        res = auxReACTed(waitinglist, hybridSystem, dimLength, hybridSystem.locations[location], interval′, initialset, dirsVectors, constraint, δ⁻, δ⁺, flowPhiDict, alg, maxOrder, reduceOrder, TΦ, discDict, inputDict, nonIntersectedDiscDict, saveResult)
        reachset = vcat(reachset, res)
    end


    return reachset
end


# AucReacted is called recursively each time we have a new starting location (after a transition)
function auxReACTed(waitlist, hybridSystem, dim, loc::Location, interval, X0, dirs, constraint, δ⁻::Float64, δ⁺::Float64, PhiDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, Φ=missing, discretizationDict=nothing, inputDiscritezationDict=nothing, polyhedralSet=nothing, saveResult::Bool=true; clustering=false)
    locChange = false
    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    reachset = []
    setOfConstraints = vcat(loc.constraints, constraint)
    jumpDiscDict = Dict()
    #jumpOverapproximatedDiscDict = Dict()
    jumpInputDict = Dict()

    if isa(discretizationDict, Nothing) && isa(inputDiscritezationDict, Nothing)
        discretizationDict, inputDiscritezationDict = ReACTDiscretizePlus(loc, X0, δ⁻, δ⁺, alg, maxOrder, reduceOrder, PhiDict[loc.id])
        locChange = true
    end

    overapproximatedDiscretizationDict = Dict()

    intersectedSetIsNothing = isnothing(polyhedralSet)
    intersectedDict = Dict()
    if intersectedSetIsNothing
        #intersectedDict = deepcopy(discretizationDict)
    else 
        #@show polyhedralSet
        intersectedDict, _ = ReACTDiscretizePlus(loc, polyhedralSet, δ⁻, δ⁺, alg, maxOrder, reduceOrder, PhiDict[loc.id])
    end


    for key in keys(discretizationDict)
        #@show isempty(discretizationDict[key])
        tempval = overapproximate(discretizationDict[key], BoxDirections(dim))
        @show norm(tempval)
        if !intersectedSetIsNothing
            @show norm(intersectedDict[key])
            tempval = intersection(tempval, intersectedDict[key])
            @show norm(tempval)
            if norm(tempval) == 0.0
                if saveResult
                    #push!(reachset, ([(map(x -> ρ(x,intersectedDict[key] ), dirs), [time, time + δ⁺])], "Guard intersection: " * string(time) * " - " * string(time) * ": " * string(loc.id) * "->" * string(loc.id)))
                    push!(reachset, ([(map(x -> ρ(x,X0 ), dirs), [time, time + δ⁺])], "Guard intersection: " * string(time) * " - " * string(time) * ": " * string(loc.id) * "->" * string(loc.id)))
                    
                end
                #return reachset
            end
            
        end
        #println(norm(tempval))

        overapproximatedDiscretizationDict[key] = tempval
    end
    #t = overapproximate(X0, BoxDirections(dim))

    println("Finished Dicts")


    # For each edge we simulate the system
    for edge in loc.edges
        reset_map(X) = MinkowskiSum(LazySets.LinearMap(copy(edge.jumpMatrix), X), Singleton(edge.jumpVector)) #edge.jumpMatrix * X + edge.jumpVector
        Reset_Map(X) = bloatPolytope(Singleton(edge.jumpVector),edge.jumpMatrix, X) 
        if VERBOSE
            println("Handling edge at time $time from $(loc.id) -> $(edge.targetLoc)")
        end
        guards = edge.guard # Technically the guard is one singular HPolyhedron, but it composes the other guards


        #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 

        tempReachset, reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, [time, endtime], guards, setOfConstraints, dirs, 2, PhiDict[loc.id], overapproximatedDiscretizationDict, discretizationDict, inputDiscritezationDict, saveResult)

        if reachtime - time == 0.0
            if !isempty(reachset)
                if VERBOSE
                    
                    println("Found an immediate transition to $(edge.targetLoc), but we are not taking it as we are scared of zeno behaviour")
                end
            end
            if saveResult
                push!(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * " took no steps"))
            end
            continue
        end
        if reachtime < endtime
            #println(string(time) * " - " * string(reachtime))
            if saveResult
                push!(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(edge.targetLoc)))
            end
            #println("Before touches")
            intersectingSet, timeNotIntersected = ReACTTouches(loc, δ⁻, δ⁺, [reachtime, endtime], (reachtime - time), guards, setOfConstraints, 2, PhiDict[loc.id], overapproximatedDiscretizationDict, discretizationDict, inputDiscritezationDict, tΦ, nothing)
            #println("After touches")
            
            if VERBOSE
                println("Current time: $time, intersecting time start: $reachtime to $timeNotIntersected")
            end
            
            timeIntersected = timeNotIntersected - reachtime

            #
            #   Here we should check whether we have reached endtime. If true we should only push the jumpSet
            #   Still need to check whether we have reached the invariant. If true we should NOT push the else branch result, only the tempReachsets[infMaxsIdx]
            #
            if clustering

            else
                @show length(intersectingSet)
                for (intersectedSet, nonIntersectedSet, startTime) in intersectingSet
                    if LazySets.isempty(intersectedSet)
                        continue
                    end
                    #return reachset
                    if saveResult
                        #return reachset
                        #push!(reachset, ([(map(x -> ρ(x, nonIntersectedSet), dirs), [startTime, startTime + δ⁻])], "Guard intersection: " * string(reachtime) * " - " * string(timeNotIntersected) * ": " * string(loc.id) * "->" * string(loc.id)))
                        #tintersectedSet = revise(intersectedSet, nonIntersectedSet, collect(BoxDirections(dim)))
                        #push!(reachset, ([(map(x -> ρ(x, tintersectedSet), dirs), [reachtime- δ⁺, reachtime + δ⁺])], "Guard intersection: " * string(reachtime) * " - " * string(timeNotIntersected) * ": " * string(loc.id) * "->" * string(loc.id)))
                        
                        #push!(reachset, ([(map(x -> ρ(x, intersectedSet), dirs), [reachtime, reachtime + δ⁺])], "Guard intersection: " * string(reachtime) * " - " * string(timeNotIntersected) * ": " * string(loc.id) * "->" * string(loc.id)))
                        
                    end
                    #continue

                    println("Before revise: $(LazySets.isempty(intersectedSet))")
                    @show (startTime, norm(nonIntersectedSet), norm(intersectedSet))
                    tintersectedSet = revise(intersectedSet, nonIntersectedSet, collect(BoxDirections(dim))) #, collect(BoxDirections(dim))
                    println("After revise: $(LazySets.isempty(tintersectedSet))")
                    if !LazySets.isempty(tintersectedSet)
                        if saveResult
                            push!(reachset, ([(map(x -> ρ(x, tintersectedSet), dirs), [reachtime, startTime])], "Guard intersection: " * string(reachtime) * " - " * string(timeNotIntersected) * ": " * string(loc.id) * "->" * string(loc.id)))
                        end
                        jumpSetLazy = reset_map(nonIntersectedSet) #nonIntersectedSet # 
                        #@show minimum(map(x -> norm(x.a) ,constraints_list(intersectedSet)))
                        #@show minimum(map(x -> abs(x.b) ,constraints_list(intersectedSet)))
                        #@show typeof(intersectedSet)
                        #jumpSetIntersected = linear_map(edge.jumpMatrix, intersectedSet)
                        jumpSetIntersected = Reset_Map(tintersectedSet) #intersection(Reset_Map(intersectedSet), hybridSystem.locations[edge.targetLoc].invarient)
                        
                        if !isnothing(hybridSystem.locations[edge.targetLoc].invarient)
                            #jumpSetLazy = Intersection( hybridSystem.locations[edge.targetLoc].invarient, jumpSetLazy) #reset_map(nonIntersectedSet) #
                            jumpSetIntersected = LazySets.intersection(jumpSetIntersected, hybridSystem.locations[edge.targetLoc].invarient)
                            
                        end
                        jumpSetIntersected = revise(jumpSetIntersected, jumpSetLazy, collect(BoxDirections(dim)))
                        @show LazySets.isempty(jumpSetIntersected)
                        push!(waitlist, (edge.targetLoc, jumpSetLazy, [startTime, endtime], missing, nothing, nothing, jumpSetIntersected))

                    else
                        if saveResult
                    
                            #push!(reachset, ([(map(x -> ρ(x, intersectedSet), dirs), [reachtime, reachtime + δ⁺])], "Guard intersection: " * string(reachtime) * " - " * string(timeNotIntersected) * ": " * string(loc.id) * "->" * string(loc.id)))
                            push!(reachset, ([(map(x -> ρ(x, nonIntersectedSet), dirs), [reachtime, reachtime + δ⁺])], "Guard intersection: " * string(reachtime) * " - " * string(timeNotIntersected) * ": " * string(loc.id) * "->" * string(loc.id)))
                            #return reachset
                            
                        end
                    end
                end
            end
        else # Reached the end time
            if VERBOSE
                println("Instant transition")
            end
            if saveResult
                if !isa(loc.invarient, Nothing)
                    newReach = []
                    for (Z, timeInterval) in tempReachset
                        #=
                        if intersects(Z, loc.invarient)
                            push!(newReach, (zonotopeStripIntersection(concretize(Z), loc.invarient), timeInterval))
                            push!(newReach, (zonotopeStripIntersection(concretize(Z), loc.invarient), timeInterval))
                        else
                            push!(newReach, (Z, timeInterval))
                        end
                        =#
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
        if VERBOSE
            
            println("No edges? call ReACT")
        end

        #tempReachset, tempInput, reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, [time, endtime], guards, setOfConstraints, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, missing, saveResult)


        #tempReachset, _, _ = ReACT(loc, δ⁻, δ⁺, interval, constraint, 2, alg, maxOrder, reduceOrder, PhiDict[loc.id])
        tempReachset, reachtime, _ = ReACT(loc, δ⁻, δ⁺, [time, endtime], setOfConstraints, dirs, 2, PhiDict[loc.id], discretizationDict, inputDiscritezationDict, saveResult)

        if saveResult
            reachset = vcat(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(loc.id)))
        end
    end




    return reachset
end

function ReACTTouches(loc, δ⁻::Float64, δ⁺::Float64, interval, initialTime::Float64, guard, constraint, STRATEGY::Integer, PhiDict, discritezationDict, lazyDiscritezationDict, inputDiscritezationDict, Φ, accInput)
    m = copy(δ⁻)
    initialTimeStep = copy(δ⁺)
    changedTimeStep = true
    phiDict = PhiDict

    #@show Φ

    # constraintProjVectors = map(x -> x.a, constraint)
    # constraintProjBounds = ρ.(constraintProjVectors, constraint)
    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint) #constraintProjBounds = map(x -> x.a, constraint) #
    guardProjVectors, guardProjBounds = getHalfSpaceProjections(guard) #map(x -> x.a, guard.constraints) #
    invarientProjVectors, invarientProjBounds = getHalfSpaceProjections(loc.invarient) #map(x -> x.a, loc.invarient.constraints) #
    @show invarientProjVectors
    @show guardProjVectors

    #oldDirProjVectors = copy(dirProjVectors)
    oldConstraintProjVectors = copy(constraintProjVectors)
    oldGuardProjVectors = copy(guardProjVectors)
    oldInvarientProjVectors = copy(invarientProjVectors)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    V = copy(inputDiscritezationDict[initialTimeStep])

    Vs = nestedInputDiscCalculate(inputDiscritezationDict, PhiDict, δ⁺, δ⁻, initialTime)
    #println(Vs)
    #lastVs = copy(Vs)
    Sρ = map(x -> ρ(x, Vs), constraintProjVectors) #zeros(Float64, length(constraint))

    Gρ = map(x -> -ρ(-x, Vs), guardProjVectors)#zeros(Float64, length(guardProjVectors))
    Iρ = map(x -> -ρ(-x, Vs),invarientProjVectors) #zeros(Float64, length(invarientProjVectors))
    #println("precomputed!")
    
    #dρ = zeros(Float64, length(dirs))
    newR = discritezationDict[initialTimeStep]
    i = 1


    if ismissing(Φ)
        Φ::Matrix{Float64} = exp((initialTime) .* loc.A)

    else
        #@show exp((initialTime) .* loc.A) - Φ
    end

    permutedphiDict = Dict()
    for key in keys(phiDict)
        permutedphiDict[key] = permutedims(phiDict[key])
    end

    tempM = diagm(ones(Float64, size(loc.A, 2)))
    ϕt = similar(Φ)
    #newRR = copy(newR)
    newRR = newR
    overapproximateIntersectingSetArray = []
    #@show Vs
    #@show newR
    #println(unique(map(x -> x.a, constraints_list(lazyDiscritezationDict[δ⁻]))))
    
    push!(overapproximateIntersectingSetArray, [constrain(Vs, linear_map(Φ, newR), LinearMap(copy(Φ), lazyDiscritezationDict[δ⁻]), vcat(invarientProjVectors), vcat(invarientProjBounds)), MinkowskiSum(Vs, LinearMap(copy(Φ), lazyDiscritezationDict[δ⁻])), time])
    #return (overapproximateIntersectingSetArray,  time)
    
    #println("first!")
    
    #return overapproximateIntersectingSetArray, Vs, time
    #preclustering = (Φ * discritezationDict[m]) ⊕ Vs
    preclustering = nothing
    attemptsRecorder = []
    triedRevise = false
    while time < endtime
        #println("Touching at time $time with newRR: ", zonotopePrintDim(newRR , 5))
        if VERBOSE

            println("Time $time, calculating intersect")
        end
        #println("Time $time, calculating intersect")

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

                if all(((-ρ(-x, newRR ⊕ Vs)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds))
                    #if !isdisjoint(newRR ⊕ Vs, guard; algorithm="sufficient")
                    if VERBOSE
                        println("Intersecting guard, but we are going to hit an invarient")
                    end
                    if isnothing(preclustering)
                        preclustering = newRR ⊕ Vs
                    else
                        preclustering = UnionSet(preclustering, newRR ⊕ Vs)
                    end
                    if isempty(overapproximateIntersectingSetArray)
                        println("Is empty_!_!")
                        push!(overapproximateIntersectingSetArray, copy([constrain(Vs, newRR, LinearMap(copy(Φ), lazyDiscritezationDict[δ⁻]), invarientProjVectors,  invarientProjBounds),MinkowskiSum(Vs, LinearMap(copy(Φ), lazyDiscritezationDict[δ⁻])), time]))
                        #push!(overapproximateIntersectingSetArray, copy([newSet, MinkowskiSum(Vs, LinearMap(Φ, lazyDiscritezationDict[currentTimeStep])), time])) #
                        
                    end
                end
                println("i: $i @$time")
                return (overapproximateIntersectingSetArray, time)

                #println("Touches Vs: ", accInput)
                #println(norm(intersectingSet), " ", norm(newRR))
            end

            if changedTimeStep
                newR = discritezationDict[currentTimeStep]
                V = inputDiscritezationDict[currentTimeStep]
                ϕt = phiDict[currentTimeStep]
                newRR = linear_map(Φ, newR)
                V = linear_map(Φ, V)
            else
                newRR = linear_map(ϕt, newRR)
                V = linear_map(ϕt, V)
            end

            changedTimeStep = false

            #if all(&, <=(hom + inhom, constraintProjBounds)) && intersects(concretize(minkowski_sum(newRR, Vs)), guard) && (isnothing(loc.invarient) || isSubSet(concretize(minkowski_sum(newRR, Vs)), loc.invarient)) #mapreduce(x -> intersects(newRR, x), &, guard)
            # println("Constraint check: ", all((ρ(x, newRR) + ρ(x, Vs)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) )    
            # println("Guard check: ", intersects(concretize(minkowski_sum(newRR, Vs)), guard))
            # println("Invarient check: ", (isnothing(loc.invarient) || isSubSet(concretize(minkowski_sum(newRR, Vs)), loc.invarient)))
            # tempSet = Minkowski_sum(newRR, Vs)
            tempSet = newRR ⊕ Vs
            #@show norm(newRR)
            #println("First")
            sen = all((ρ(x, newRR) + ρ(x, Vs)) <= y for ( x, y) in zip( constraintProjVectors, constraintProjBounds))
            #println("Second: $sen")
            #println(isempty(newRR))
            
            sen = sen && all((-ρ(-x, Vs) + (-ρ(-x, newRR)) <= y) for ( x, y) in zip( guardProjVectors, guardProjBounds))
            #println("Third: $sen")
            
            sen = sen && any((-ρ(-x, Vs) + (-ρ(-x, newRR))) <= y for ( x, y) in zip(invarientProjVectors, invarientProjBounds))
            #println("After third: $sen")
                #=
                if all((ρ(x, newRR) + ρ(x, Vs)) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds)) && # IsSubSet
               all(((-ρ(-x, tempSet)) <= y) for (x, y) in zip(guardProjVectors, guardProjBounds)) &&
               #!isdisjoint(tempSet, guard; algorithm="sufficient") && # intersects
               any((ρ(-x, newRR) + ρ(-x, Vs)) <= y for (x, y) in zip(invarientProjVectors, invarientProjBounds)) # IsSubSet
                =##if mapreduce(x -> intersects(newRR, x), &, guard)
            if sen
                #println("Before constrain")
                
                newSet = constrain(Vs, newRR, LinearMap(copy(Φ), lazyDiscritezationDict[currentTimeStep]), vcat(guardProjVectors, invarientProjVectors), vcat(guardProjBounds, invarientProjBounds))
                @show (norm(newRR), norm(newSet), norm(Vs), norm(MinkowskiSum(copy(Vs), LinearMap(copy(Φ), lazyDiscritezationDict[currentTimeStep]))))
                #println("After constrain")
                if !isempty(newSet)
                    #@show triedRevise
                    push!(overapproximateIntersectingSetArray, [newSet, MinkowskiSum(copy(Vs), LinearMap(copy(Φ), lazyDiscritezationDict[currentTimeStep])), time]) #
                    if i == 3    
                        #return (overapproximateIntersectingSetArray, time)
                    end
                    
                else
                    @show isempty(newSet)

                    return (overapproximateIntersectingSetArray, time)
                end
            
                
                if isnothing(preclustering)
                    preclustering = tempSet
                else
                    preclustering = UnionSet(preclustering, tempSet)
                end
                #lastVs = copy(Vs)
                Sρ += map(x -> ρ(x, V), constraintProjVectors)
                Gρ += map(x -> ρ(-x, V), guardProjVectors)
                Iρ += map(x -> ρ(-x, V), invarientProjVectors)
                
                #constraintProjVectors = map(x -> ϕt * x, oldConstraintProjVectors)
                #guardProjVectors = map(x -> ϕt * x, oldGuardProjVectors)
                #invarientProjVectors = map(x -> ϕt * x, oldInvarientProjVectors)

                Vs = minkowski_sum(Vs, V)
                #Vs = LinearMap(ReachabilityAnalysis.Exponentiation.Φ₁(loc.A, time + currentTimeStep - minimum(interval), ReachabilityAnalysis.Exponentiation.BaseExp), inputDiscritezationDict[0])
                approveFlag = true
                triedRevise = false
                #Sρ += inhom
                mul!(tempM, Φ, ϕt)
                copy!(Φ, tempM)
            elseif triedRevise == false
                if length(constraintProjVectors) > 0
                    changedTimeStep = true
                    #println("Before revise")
                    
                    tSet = revise(discritezationDict[currentTimeStep], lazyDiscritezationDict[currentTimeStep], map(x -> permutedims(Φ) *x, constraintProjVectors))
                    @show isempty(tSet)
                    
                    if !isempty(tSet)

                        discritezationDict[currentTimeStep] = tSet
                    elseif VERBOSE
                        println(i)
                        return (overapproximateIntersectingSetArray, time)
                    else
                        return (overapproximateIntersectingSetArray, time)
                    end
                end
                triedRevise = true
            else
                newR = copy(newR)
                currentTimeStep = currentTimeStep / 2
                changedTimeStep = true
                attempts = attempts + 1
                triedRevise = false
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

    return (overapproximateIntersectingSetArray, time)
end


function ReACTGuards(loc, δ⁻::Float64, δ⁺::Float64, interval, guards, constraint, dirs, STRATEGY::Integer, PhiDict, discritezationDict, lazyDiscritezationDict, inputDiscritezationDict, saveResult)
    # We calculate the reachset till we reach a guard for an intersection (or till failure)
    initialTimeStep = copy(δ⁺)
    m = copy(δ⁻)
    changedTimeStep = true
    phiDict = PhiDict
    A = copy(loc.A)
    # constraintProjVectors = map(x -> x.a, constraint)
    # constraintProjBounds = ρ.(constraintProjVectors, constraint)
    constraintProjVectors, constraintProjBounds = getHalfSpaceProjections(constraint)
    guardProjVectors, guardProjBounds = [], []
    if isnothing(guards)
        guardProjVectors, guardProjBounds = [constraintProjVectors[1]], [Inf]
    else
        guardProjVectors, guardProjBounds = getHalfSpaceProjections(guards)
    end

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


    Φ::Matrix{Float64} = diagm(ones(Float64, size(loc.A, 2)))

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

    triedRevise = false
    while time < endtime
        if VERBOSE

            println("Time is: $time")
        end
        attempts = 1
        approveFlag = false

        while !approveFlag
            if VERBOSE
                println("Stuck?")
            end
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

                V = inputDiscritezationDict[currentTimeStep]
                ϕt = phiDict[currentTimeStep]
                pϕt = permutedphiDict[currentTimeStep]
            end

            # constraintProjVectors = map(x -> pϕt * x, oldConstraintProjVectors)
            # guardProjVectors = map(x -> pϕt * x, oldGuardProjVectors)
            # invarientProjVectors = map(x -> pϕt * x, oldInvarientProjVectors)


            changedTimeStep = false
            #println(isempty(newR))
            #println(newR)
            if all((input + ρ(x, newR)) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) &&
               #any(sign(y) >= 0 ? (input + ρ(-x, newR)) <= y : !((input + ρ(x, newR)) < y) for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) &&
               any((input + -ρ(-x, newR)) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) &&
               all((input + ρ(x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))
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

                approveFlag = true
                triedRevise = false
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
            elseif triedRevise == false
                
                changedTimeStep = true
                tSet = revise(discritezationDict[currentTimeStep], lazyDiscritezationDict[currentTimeStep], vcat(constraintProjVectors, -1 * guardProjVectors, invarientProjVectors), vcat(constraintProjBounds, -1 * guardProjBounds, invarientProjBounds)) # - vcat(Sρ, Gρ, Iρ)
                @show isempty(tSet)
                if !isempty(tSet)
                    
                    discritezationDict[currentTimeStep] = tSet
                elseif VERBOSE
                    println(time)
                    throw(DomainError(discritezationDict[currentTimeStep], "Set is empty?"))
                else
                    throw(DomainError(discritezationDict[currentTimeStep], "Set is empty? Location: $(loc.id)"))
                end
                
                triedRevise = true
            else
                newR = copy(newR)
                currentTimeStep = currentTimeStep / 2
                changedTimeStep = true
                attempts = attempts + 1
                triedRevise = false
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
        if VERBOSE
            println("Time iss: $time")
        end
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
    #throw(error("ERROR!!! We have hit a constraint at loc: $(locationId) time: $time"))
    println("\nERROR!!!\n 
            ERROR!!!\n\n
            We have hit a constraint at loc: $locationId time: $time\n\n
            ERROR!!!\n
            ERROR!!!\n")
end