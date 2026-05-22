using LazySets, LinearAlgebra
include("../../Discretize.jl")
include("../../Utilities.jl")
include("../../ReACTed.jl")
include("embrake.jl")

function solve_embrake(sys, initialState, X0, T, δ⁺ = 2*10^-7, δ⁻ = 2*10^-7, maxOrder = 5, reduceOrder = 5, dirs = [], saveResult = true)
    time = 0

    x0 = 0.05
    Tsample = 1e-4
    ζ = 1e-6
    n = length(X0.center)
    reachset = []
    phiDict = Dict()
    
    dirsVectors = []
    dimLength = size(X0.center, 1)
    for dir in dirs
        template = zeros(dimLength)
        template[dir] = 1.

        push!(dirsVectors, template)
        push!(dirsVectors, -template)
    end


    loc = sys.locations[initialState]
    alg = ReachabilityAnalysis.Exponentiation.BaseExp

    # Get the phi dict system 
    flowPhiDict = Dict(map(x -> x.id => PhiDict(x.A, δ⁻, δ⁺, alg), sys.locations))
    phiDict = flowPhiDict[loc.id]
    discretizationDict, inputDiscretizationDict = ReACTDiscretizePlus(loc, X0, δ⁻, δ⁺, alg, maxOrder, reduceOrder, phiDict)
    

    res = run_embrake(sys, loc, time, Tsample, ζ, x0, T, discretizationDict, inputDiscretizationDict, sys.globalConstraints, dirsVectors, δ⁻, δ⁺, phiDict, alg, maxOrder, reduceOrder, missing, saveResult)

    reachset = vcat(reachset, res)

    return res
end

function run_embrake(hybridSystem, loc::Location, time, Tsample, ζ, x0, T, discretizationDict, inputDiscretizationDict, constraint, dirs, δ⁻::Float64, δ⁺::Float64, PhiDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, Φ=missing, saveResult::Bool=true) where {N}
    guardTime = Tsample - ζ
    invariantTime = Tsample + ζ
    
    reachset = []
    setOfConstraints = vcat(loc.constraints, constraint)
    edge = loc.edges[1]

    #JumpSupport(X, dir) = ρ(transpose(edge.jumpMatrix)*dir, X) + ρ(dir, edge.jumpVector) # TODO: This is not correct
    JumpZonotope(X) = minkowski_sum(linear_map(edge.jumpMatrix, X), Singleton(edge.jumpVector))

    tempReachset, reachtime, tΦ = ReACTGuards(loc, δ⁻, δ⁺, [time, time+guardTime], nothing, setOfConstraints, dirs, 2, PhiDict, discretizationDict, inputDiscretizationDict, saveResult)
  
    #println(ρ(dirs[1], discretizationDict[δ⁺]))
    if reachtime < T 
        if saveResult
            push!(reachset, (tempReachset, string(time) * " - " * string(reachtime) * ": " * string(loc.id) * "->" * string(edge.targetLoc)))
        end
        #_, _, timeNotIntersected = ReACTTouches(loc, δ⁻, δ⁺, [reachtime, time+invariantTime], (reachtime - time),
        #                                                                    nothing, setOfConstraints, 2, PhiDict, 
        #                                                                    discretizationDict, discretizationDict, inputDiscretizationDict, tΦ, nothing)
        timeIntersected =  invariantTime - guardTime
        time = reachtime
        println("time = ", time)
        λ = ceil(timeIntersected/(δ⁺))
        i = 0
        jumpΦ = tΦ
        while i < λ
            branchedRun = []
            d = δ⁻
            jumpDiscDict = Dict()
            while d <= δ⁺
                tempSet = linear_map(jumpΦ, discretizationDict[d]) 
                #tempJump = JumpSupport(discretizationDict[d], tempPhi*dir)
                tempJump = JumpZonotope(tempSet)
                if (ρ(dirs[1], tempJump) > x0) || (ρ(dirs[1], tempSet) > x0)
                    handleHitConstraint(time, 1)
                else
                    jumpDiscDict[d] = tempJump
                end
                d = 2*d
            end
            jumpΦ = jumpΦ* PhiDict[δ⁺]
            if time < T
                branchedRun = run_embrake(hybridSystem, loc, time, Tsample, ζ, x0, T, jumpDiscDict, inputDiscretizationDict, hybridSystem.globalConstraints, dirs, δ⁻, δ⁺, PhiDict, alg, maxOrder, reduceOrder, missing, saveResult)
            end
            time = time + δ⁺
            i = i+1
            if saveResult && !isnothing(branchedRun)
                reachset = vcat(reachset, branchedRun)
            end 
        end
    end
    return reachset
end