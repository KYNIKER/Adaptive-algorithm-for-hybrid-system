using LazySets, LinearAlgebra
include("../../Discretize.jl")
include("../../Utilities.jl")
include("../../ReACTed.jl")
include("embrake.jl")

function solve_embrake(δ⁺ = 2*10^-7, δ⁻ = 2*10^-7, maxOrder = 5, reduceOrder = 5, dirs = [], saveResult = true)
    time = 0

    x0 = 0.05
    Tsample = 1e-4
    ζ = 1e-6
    sys, initialState, X0, T = loadembrake(x0, Tsample, ζ)
    n = length(X0.center)
    reachset = []

    loc = sys.locations[initialState]
    alg = ReachabilityAnalysis.Exponentiation.BaseExp

    # Get the phi dict system 
    flowPhiDict = Dict(map(x -> x.id => PhiDict(x.A, δ⁻, δ⁺, alg), sys.locations))
    res = run(sys, loc, time, Tsample, ζ, T, X0, sys.globalConstraints, dirs, δ⁻, δ⁺, flowPhiDict, alg, maxOrder, reduceOrder, missing, saveResult)

    reachset = vcat(reachset, res)

    return res
end

function run(hybridSystem, loc::Location, time, Tsample, ζ, T, X0::Zonotope{N,Vector{N},Matrix{N}}, constraint, dirs, δ⁻::Float64, δ⁺::Float64, PhiDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, Φ=missing, saveResult::Bool=true) where {N}
    guardTime = Tsample - ζ
    invariantTime = Tsample + ζ
    
    reachset = []
    discretizationDict, inputDiscritezationDict = ReACTDiscretizePlus(loc, X0, δ⁻, δ⁺, alg, maxOrder, reduceOrder, PhiDict[loc.id])
    setOfConstraints = vcat(loc.constraints, constraint)
    edge = loc.edges[1]
    guards = edge.guard
    tempReachset, reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, [time, time+guardTime], nothing, setOfConstraints, dirs, 2, PhiDict[loc.id], discretizationDict, discretizationDict, inputDiscritezationDict, saveResult)

    if reachtime < T 
        if saveResult
            push!(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(edge.targetLoc)))
        end
        intersectingSet, intersectedInput, timeNotIntersected = ReACTTouches(loc, δ⁻, δ⁺, [reachtime, time+invariantTime], (reachtime - time),
                                                                            guards, setOfConstraints, 2, PhiDict[loc.id], 
                                                                            discretizationDict, discretizationDict, inputDiscritezationDict, tΦ, nothing)
        timeIntersected = timeNotIntersected - reachtime

        time = time + invariantTime # TODO : Find the correct time to add here

        if !isnothing(intersectingSet)
            intersectedSet = intersectingSet

            if !isa(guards, Nothing) && !isdisjoint(intersectedSet, guards)
                intersectedSet = zonotopeStripIntersection(intersectedSet, guards)
            end
            if isempty(intersectedSet)
                println("Empty..")
                return reachset
            end
            if !isa(loc.invarient, Nothing) && !isdisjoint(intersectedSet, loc.invarient)
                    tintersectedSet = zonotopeStripIntersection(intersectedSet, loc.invarient)
                    intersectedSet = tintersectedSet
            end
            if isempty(intersectedSet)
                println("Empty...")
                return reachset
            end
            push!(reachset, ([(intersectedSet, [reachtime, timeNotIntersected])], string(reachtime) * " - " * string(timeNotIntersected) * ": " * string(loc.id) * "->" * string(loc.id)))

            if true
                jumpSet = minkowski_sum(linear_map(edge.jumpMatrix, intersectedSet), Zonotope(edge.jumpVector, [zero(edge.jumpVector)]))

                if !isa(hybridSystem.locations[edge.targetLoc].invarient, Nothing) && intersects(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)
                    tjumpSet = zonotopeStripIntersection(jumpSet, hybridSystem.locations[edge.targetLoc].invarient)# + edge.jumpVector * intersectedSet
                    jumpSet = tjumpSet
                end

                branchedRun = run(hybridSystem, hybridSystem.locations[edge.targetLoc], time, Tsample, ζ, T, jumpSet, constraint, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder, tΦ, saveResult)

                if saveResult
                    reachset = vcat(reachset, branchedRun)
                end
            else
                nonintersectedSet = minkowski_sum(linear_map(exp(timeNotIntersected .* loc.A), X0), concretize(ReachabilityAnalysis.Exponentiation.Φ₁(loc.A, timeNotIntersected, ReachabilityAnalysis.Exponentiation.BaseExp) * inputDiscritezationDict[0]))
                if !isa(loc.invarient, Nothing) && !isdisjoint(nonintersectedSet, loc.invarient)
                    nonintersectedSet = zonotopeStripIntersection(nonintersectedSet, loc.invarient)
                end
                branchedRun = run(hybridSystem, loc, time, Tsample, ζ, T, nonintersectedSet, constraint, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder)
                if saveResult
                    reachset = vcat(reachset, branchedRun)
                end
            end
        else
            println("We have no intersecting set. Meaning we hit an invarient and have no guards fulfilled")
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

    return reachset
end

