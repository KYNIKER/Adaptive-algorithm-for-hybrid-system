export ReACTDiscretize, PhiDict
using LinearAlgebra, LazySets, ReachabilityAnalysis
using ReachabilityBase.Arrays: isinvertible

#isinvertible(x) = applicable(inv, x) && isone(inv(Matrix(x)) * x)

function ReACTDiscretize(loc, X0::Zonotope{N,Vector{N},Matrix{N}}, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, phiDict=nothing) where {N}
    XDim, _ = size(genmat(X0))
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    inputDiscritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    A = loc.A
    if isnothing(phiDict)
        phiDict = PhiDict(A, δ⁻, δ⁺, alg)
    end

    U = isnothing(loc.B) ? (isnothing(loc.u) ? Zonotope(zeros(XDim), zeros(XDim, 1)) : loc.u) : concretize(loc.B * loc.u)
    if !isnothing(loc.c)
        U = Zonotope(U.center + loc.c, genmat(U))
    end

    d = δ⁻
    dia::Matrix{Float64} = diagm(ones(XDim))
    isInvA = false #isinvertible(A)
    Φ = copy(phiDict[d])
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    Φcache = sum(A) == abs(sum(A)) ? Φ : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)
    #pis = ReachabilityAnalysis.Exponentiation.Φ₁(A, δ⁻, alg, isInvA, Φcache)

    inputDiscritezationDict[0] = U
    #if !(zeros(XDim) ∈ U) #Origin is *not* in input
    #println("Here")
    dU = overapproximate(δ⁻ * U, Zonotope)
    E_ψ = convert(Zonotope, symmetric_interval_hull(P2A_abs * symmetric_interval_hull(A * U)))
    P = minkowski_sum(dU, E_ψ) #

    lt = minkowski_sum(convert(Zonotope, phiDict[d] * X0), dU)  #minkowski_sum(convert(Zonotope, phiDict[d] * X0), dU)
    E⁺ = convert(Zonotope, symmetric_interval_hull(P2A_abs * symmetric_interval_hull(A * A * X0)))
    rt = minkowski_sum(E_ψ, E⁺)
    f = minkowski_sum(lt, rt)
    disc = overapproximate(CH(X0, f), Zonotope) # overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope) #
    #disc = Zonotope(disc.center - P̂, genmat(disc))
    #return discritezationDict, inputDiscritezationDict
    while d < δ⁺
        discritezationDict[d] = disc
        if maxOrder > 0
            if LazySets.order(P) > maxOrder
                P = reduce_order(P, reduceOrder)
            end
            if LazySets.order(disc) > maxOrder
                disc = reduce_order(disc, reduceOrder)
            end
        end
        P = concretize(ReachabilityAnalysis.Exponentiation.Φ₁(A, d, alg, isInvA, Φcache) * U)
        inputDiscritezationDict[d] = P
        #println(P)

        disc = overapproximate(CH(disc, minkowski_sum(P, linear_map(phiDict[d], disc))), Zonotope)
        d = d * 2
    end
    if maxOrder > 0
        if LazySets.order(P) > maxOrder
            P = reduce_order(P, reduceOrder)
        end
        if LazySets.order(disc) > maxOrder
            disc = reduce_order(disc, reduceOrder)
        end
    end
    discritezationDict[δ⁺] = copy(disc)
    inputDiscritezationDict[δ⁺] = P

    return discritezationDict, inputDiscritezationDict
end

function ReACTDiscretizePlus(loc, X0::Zonotope{N,Vector{N},Matrix{N}}, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, phiDict=nothing) where {N}
    #XDim, _ = size(genmat(X0))


    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    #discritezationDict = Dict()
    inputDiscritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    #inputDiscritezationDict = Dict()
    A = loc.A
    XDim = size(A, 1)
    if isnothing(phiDict)
        phiDict = PhiDict(A, δ⁻, δ⁺, alg)
    end

    #U = isnothing(loc.B) ? (isnothing(loc.u) ? Zonotope(zeros(XDim), [zeros(XDim)]) : loc.u) : concretize(loc.B * loc.u)
    U = Zonotope(zeros(XDim), zeros(XDim, 1))
    U = isnothing(loc.B) ? (isnothing(loc.u) ? Zonotope(zeros(XDim), zeros(XDim, 1)) : loc.u) : linear_map(loc.B, loc.u)
    if !isnothing(loc.c)
        U = concretize(U)
        U = Zonotope(U.center + loc.c, genmat(U))
    end

    d = δ⁻
    #dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    isInvA = false #isinvertible(A)
    Φ = copy(phiDict[d])
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    Φcache = sum(A) == abs(sum(A)) ? Φ : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)
    #pis = ReachabilityAnalysis.Exponentiation.Φ₁(A, δ⁻, alg, isInvA, Φcache)

    #X0 = Zonotope([1., 0., -1.], [[0.0, 0.0, 0.0]])

    inputDiscritezationDict[0] = U
    #if !(zeros(XDim) ∈ U) #Origin is *not* in input
    #println("Here")
    dU = overapproximate(LinearMap(δ⁻, U), Zonotope)#linear_map(dia, U)
    E_ψ = SymmetricIntervalHull(linear_map(P2A_abs, SymmetricIntervalHull(LazySets.linear_map(A, U))))
    #E_ψ = SymmetricIntervalHull(LinearMap(P2A_abs, SymmetricIntervalHull(A * U)))
    P = minkowski_sum(dU, E_ψ) #
    lt = minkowski_sum(LazySets.linear_map(phiDict[d], X0), dU)  #minkowski_sum(convert(Zonotope, phiDict[d] * X0), dU)
    E⁺ = SymmetricIntervalHull(LazySets.linear_map(P2A_abs, SymmetricIntervalHull(LazySets.linear_map(A * A, X0))))
    rt = minkowski_sum(E_ψ, E⁺)
    f = minkowski_sum(lt, rt)

    disc = overapproximate(CH(X0, f), Zonotope) # overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope) #

    # TODO maybe we can reuse this somehow?
    # if (size(genmat(disc),2)==0)
    #     disc = X0
    # end

    #disc = Zonotope(disc.center - P̂, genmat(disc))
    #P = LinearMap(ReachabilityAnalysis.Exponentiation.Φ₁(A, d, alg, isInvA, Φcache), U)
    #return discritezationDict, inputDiscritezationDict
    #inputDiscritezationDict[d] = P

    while d < δ⁺
        discritezationDict[d] = disc
        inputDiscritezationDict[d] = P
        if maxOrder > 0
            if LazySets.order(P) > maxOrder
                P = reduce_order(P, reduceOrder)
            end
            if LazySets.order(disc) > maxOrder
                disc = reduce_order(disc, reduceOrder)
            end
        end
        #println(P)
        #disc = UnionSet(disc, MinkowskiSum(P, LinearMap(phiDict[d], disc)))
        disc = overapproximate(CH(disc, minkowski_sum(P, LazySets.linear_map(phiDict[d], disc))), Zonotope)
        #P = P ⊕ LinearMap(phiDict[d], P)
        P = minkowski_sum(P, LazySets.linear_map(phiDict[d], P))
        d = d * 2
    end
    # if maxOrder > 0
    #     if LazySets.order(P) > maxOrder
    #         P = reduce_order(P, reduceOrder)
    #     end
    #     if LazySets.order(disc) > maxOrder
    #         disc = reduce_order(disc, reduceOrder)
    #     end
    # end
    discritezationDict[δ⁺] = disc
    inputDiscritezationDict[δ⁺] = concretize(P)

    return discritezationDict, inputDiscritezationDict
end

function newReACTDiscretizePlus(loc, X0::Zonotope{N,Vector{N},Matrix{N}}, δ⁻::Float64, δ⁺::Float64, phiDict, inputDiscritezationDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}


    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()

    A = loc.A
    XDim = size(A, 1)

    d = δ⁻
    #dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    isInvA = isinvertible(A)
    #Φ = copy(phiDict[d])
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    Φcache = A == A_abs ? phiDict[d] : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)


    lt = minkowski_sum(LazySets.linear_map(phiDict[d], X0), inputDiscritezationDict[d])  #minkowski_sum(convert(Zonotope, phiDict[d] * X0), dU)
    E⁺ = SymmetricIntervalHull(LazySets.linear_map(P2A_abs, SymmetricIntervalHull(LazySets.linear_map(A * A, X0))))
    f = minkowski_sum(lt, E⁺)
    #f = minkowski_sum(lt, rt)

    disc = overapproximate(CH(X0, f), Zonotope) # overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope) #

    while d < δ⁺
        discritezationDict[d] = disc

        if maxOrder > 0
            if LazySets.order(disc) > maxOrder
                disc = reduce_order(disc, reduceOrder)
            end
        end

        disc = overapproximate(CH(disc, minkowski_sum(inputDiscritezationDict[d], LazySets.linear_map(phiDict[d], disc))), Zonotope)

        d = d * 2
    end
    #if maxOrder > 0
    #    if LazySets.order(disc) > maxOrder
    #disc = reduce_order(disc, reduceOrder)
    #    end
    #end
    discritezationDict[δ⁺] = disc

    return discritezationDict
end

#=
function newReACTDiscretizePlus(X0::Zonotope{N,Vector{N},Matrix{N}}, δ⁻::Float64, δ⁺::Float64, A, P2A_abs, phiDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    d = δ⁻
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    E⁺ = SymmetricIntervalHull(LazySets.linear_map(P2A_abs, SymmetricIntervalHull(LazySets.linear_map(A * A, X0))))
    f = minkowski_sum(LazySets.linear_map(phiDict[d], X0), E⁺)

    disc = overapproximate(CH(X0, f), Zonotope)

    while d < δ⁺
        discritezationDict[d] = disc

        if maxOrder > 0
            if LazySets.order(disc) > maxOrder
                disc = reduce_order(disc, reduceOrder)
            end
        end

        disc = overapproximate(CH(disc, LazySets.linear_map(phiDict[d], disc)), Zonotope)

        d = d * 2
    end
    #if maxOrder > 0
    #    if LazySets.order(disc) > maxOrder
    #disc = reduce_order(disc, reduceOrder)
    #    end
    #end
    discritezationDict[δ⁺] = disc

    return discritezationDict
end
=#

function newReACTDiscretizePlus(X0::Zonotope{N,Vector{N},Matrix{N}}, δ⁻::Float64, δ⁺::Float64, A, P2A_absp, phiDict, inputDiscritezationDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    d = δ⁻
    discritezationDictp = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    lt = minkowski_sum(linear_map(phiDict[d], X0), inputDiscritezationDict[d])  #minkowski_sum(convert(Zonotope, phiDict[d] * X0), dU)
    E⁺ = SymmetricIntervalHull(linear_map(P2A_absp, SymmetricIntervalHull(linear_map(A * A, X0))))
    f = minkowski_sum(lt, E⁺)
    #f = minkowski_sum(lt, rt)
    #XDim = size(A, 1)
    #dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    #disc = overapproximate(CH(X0, minkowski_sum(linear_map(dia, inputDict[0]), f)), Zonotope; algorithm="mean") # overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope) #
    disc = overapproximate(CH(X0, f), Zonotope) # overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope) #
    #@show disc
    while d < δ⁺
        discritezationDictp[d] = disc
        #=
        if maxOrder > 0
            if LazySets.order(disc) > maxOrder
                disc = reduce_order(disc, reduceOrder)
            end
        end
        =#
        #@show inputDiscritezationDict[d]
        rs = minkowski_sum(linear_map(phiDict[d], disc), inputDiscritezationDict[d])
        disc = overapproximate(CH(disc, rs), Zonotope)

        d = d * 2
    end
    #if maxOrder > 0
    #    if LazySets.order(disc) > maxOrder
    #disc = reduce_order(disc, reduceOrder)
    #    end
    #end
    discritezationDictp[δ⁺] = disc

    return discritezationDictp
end

function ReACT_discretize_decomposed_generators(X0G::Zonotope{N,Vector{N},Matrix{N}}, δ⁻::Float64, δ⁺::Float64, A, P2A_abs, phiDict, U, inputDiscritezationDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    d = δ⁻
    discritezationDict = Dict{Float64,Matrix{N}}()
    UG = genmat(U) .* d
    XG = genmat(X0G)
    #@show (A * genmat(U))
    AUc = abs.(A * U.center)
    AUG = sum(abs.(A * genmat(U)), dims=2)
    #@show AUc
    #@show AUG
    #@show AUG + AUc

    U_bloat = diagm(sum(abs.(genmat(inputDiscritezationDict[0])), dims=2)[:, 1])#diagm(sum(P2A_abs * diagm((AUG+AUc)[:, 1]), dims=2)[:, 1])
    G_bloat = diagm(sum(P2A_abs * diagm(sum(abs.(A * A * XG), dims=2)[:, 1]), dims=2)[:, 1])
    @show G_bloat
    @show genmat(linear_map(phiDict[d], X0G))
    newG = hcat((XG .+ phiDict[d] * XG) .* 0.5, (XG .- phiDict[d] * XG) .* 0.5, G_bloat + U_bloat, UG)

    # Some mistake in here
    while d < δ⁺
        discritezationDict[d] = copy(newG)
        #newG´ = copy(newG)
        newG´ = copy(phiDict[d] * newG)
        newG = copy(hcat((newG .+ newG´), (newG .- newG´), genmat(inputDiscritezationDict[d])) .* 0.5)


        d = d * 2
    end

    discritezationDict[δ⁺] = newG

    return discritezationDict
end


### TODO - Figure out whether this is easier if we had passed the decomposed components instead?
function ReACT_discretize_combine_with_offsets(X0c::Zonotope{N,Vector{N},Matrix{N}}, δ⁻::Float64, δ⁺::Float64, A, P2A_abs, phiDict, U, inputDiscritezationDict, generatorDiscretizationDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    d = δ⁻
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()

    Xc = X0c.center
    Uc = inputDiscritezationDict[0].center

    eAXc = (phiDict[d] * Xc)
    newC = (Xc .+ eAXc .+ Uc) .* 0.5

    c_bloat = hcat(diagm(sum(P2A_abs * diagm(abs.(A * A * Xc)), dims=2)[:, 1]), (Xc .- (eAXc .+ Uc)) .* 0.5)

    while d < δ⁺
        discritezationDict[d] = copy(Zonotope(newC, hcat(c_bloat, generatorDiscretizationDict[d])))
        eA_c_bloat = phiDict[d] * c_bloat
        c_bloat = hcat(c_bloat .+ eA_c_bloat, c_bloat .- eA_c_bloat, (newC .- (phiDict[d] * newC))) .* 0.5
        newC = (newC .+ (phiDict[d] * newC)) .* 0.5
        #c_bloat = vcat(c_bloat, newC .- (phiDict[d] * newC) .* 0.5)

        d = d * 2
    end
    #@show generatorDiscretizationDict[d]
    discritezationDict[δ⁺] = copy(Zonotope(newC, hcat(c_bloat, generatorDiscretizationDict[δ⁺])))

    return discritezationDict
end

#=function ReACTDiscretize(A, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Nothing, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, phiDict=nothing) where {N}
    XDim, _ = size(genmat(X0))
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()

    if isnothing(phiDict)
        phiDict = PhiDict(A, δ⁻, δ⁺, alg)
    end


    d = δ⁻
    dia::Matrix{Float64} = diagm(ones(XDim))
    isInvA = isinvertible(A)
    Φ = copy(phiDict[d])
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    Φcache = sum(A) == abs(sum(A)) ? Φ : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)

    invA = inv(Matrix(A))
    û = U
    P̂ = (A \ (phiDict[d] - dia)) * û
    PZ = Zonotope(P̂, zeros(Float64, XDim, 1))
    E⁺ = convert(Zonotope, symmetric_interval_hull(P2A_abs * symmetric_interval_hull(A * A * X0)))

    f = minkowski_sum(convert(Zonotope, phiDict[d] * X0), E⁺)
    disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
    while d < δ⁺
        discritezationDict[d] = disc
        if maxOrder > 0
            if LazySets.order(disc) > maxOrder
                disc = reduce_order(disc, reduceOrder)
            end
        end
        disc = overapproximate(CH(disc, linear_map(phiDict[d], disc)), Zonotope)

        d = d * 2
    end
    if maxOrder > 0
        if LazySets.order(disc) > maxOrder
            disc = reduce_order(disc, reduceOrder)
        end
    end
    discritezationDict[δ⁺] = copy(disc)

    return discritezationDict, Nothing
end=#

function PhiDict(A, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp)
    let ϕ::Matrix{Float64} = ReachabilityAnalysis.Exponentiation._exp(A, δ⁻, alg)
        phiDict = Dict{Float64,Matrix{Float64}}()
        tempM = similar(ϕ)
        d = δ⁻
        while d < δ⁺

            phiDict[d] = copy(ϕ)

            mul!(tempM, ϕ, ϕ)
            copy!(ϕ, tempM)
            # LinearMap!(tempM, ϕ, ϕ)
            # copy!(ϕ, tempM)
            d = d * 2
        end
        phiDict[δ⁺] = copy(ϕ)
        return phiDict
    end
end

function PhiInputDict(loc, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5)
    A = loc.A
    ϕ::Matrix{Float64} = ReachabilityAnalysis.Exponentiation._exp(A, δ⁻, alg)
    phiDict = Dict{Float64,Matrix{Float64}}()
    TphiDict = Dict{Float64,Matrix{Float64}}()
    tempM = similar(ϕ)

    inputDiscritezationDict = Dict{Float64,Zonotope}()
    #inputDiscritezationDict = Dict()
    XDim = size(A, 1)

    #U = isnothing(loc.B) ? (isnothing(loc.u) ? Zonotope(zeros(XDim), [zeros(XDim)]) : loc.u) : concretize(loc.B * loc.u)
    U = Zonotope(zeros(XDim), zeros(XDim, 1))
    U = isnothing(loc.B) ? (isnothing(loc.u) ? Zonotope(zeros(XDim), zeros(XDim, 1)) : loc.u) : linear_map(loc.B, loc.u)
    if !isnothing(loc.c)
        #U = concretize(U)
        U = Zonotope(U.center + loc.c, genmat(U))
    end

    d = δ⁻
    #Φ = copy(ϕ)
    dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    isInvA = isinvertible(A)
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    Φcache = A == A_abs ? ϕ : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)
    #pis = ReachabilityAnalysis.Exponentiation.Φ₁(A, δ⁻, alg, isInvA, Φcache)

    #X0 = Zonotope([1., 0., -1.], [[0.0, 0.0, 0.0]])

    inputDiscritezationDict[0] = U
    #if !(zeros(XDim) ∈ U) #Origin is *not* in input
    #println("Here")
    dU = linear_map(dia, U)#LinearMap(δ⁻, U)#
    E_ψ = symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(LazySets.linear_map(A, U))))
    #E_ψ = SymmetricIntervalHull(LinearMap(P2A_abs, SymmetricIntervalHull(A * U)))
    P = minkowski_sum(dU, E_ψ) #


    d = δ⁻
    while d < δ⁺

        phiDict[d] = copy(ϕ)
        TphiDict[d] = copy(permutedims(ϕ))
        inputDiscritezationDict[d] = P

        if maxOrder > 0
            if LazySets.order(P) > maxOrder
                P = reduce_order(P, reduceOrder)
            end
        end

        P = minkowski_sum(P, LazySets.linear_map(phiDict[d], P))


        mul!(tempM, ϕ, ϕ)
        copy!(ϕ, tempM)
        d = d * 2
    end
    #if LazySets.order(P) > maxOrder
    #    P = reduce_order(P, reduceOrder)
    #end
    inputDiscritezationDict[δ⁺] = P
    phiDict[δ⁺] = copy(ϕ)
    TphiDict[d] = copy(permutedims(ϕ))

    return phiDict, TphiDict, inputDiscritezationDict
end

function PhiInputDict(A, U, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5)

    ϕ::Matrix{Float64} = ReachabilityAnalysis.Exponentiation._exp(A, δ⁻, alg)
    phiDict = Dict{Float64,Matrix{Float64}}()
    TphiDict = Dict{Float64,Matrix{Float64}}()
    tempM = similar(ϕ)

    inputDiscritezationDict = Dict{Float64,Zonotope}()
    #inputDiscritezationDict = Dict()
    XDim = size(A, 1)

    d = δ⁻
    #Φ = copy(ϕ)
    dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    isInvA = isinvertible(A)
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    Φcache = A == A_abs ? ϕ : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)
    pis = ReachabilityAnalysis.Exponentiation.Φ₁(A, δ⁻, alg, isInvA, Φcache)

    #X0 = Zonotope([1., 0., -1.], [[0.0, 0.0, 0.0]])

    #if !(zeros(XDim) ∈ U) #Origin is *not* in input
    #println("Here")
    dU = linear_map(dia, U)#LinearMap(δ⁻, U)#
    E_ψ = symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(LazySets.linear_map(A, U))))
    inputDiscritezationDict[0] = overapproximate(minkowski_sum(dU, E_ψ), Zonotope)
    #E_ψ = SymmetricIntervalHull(LinearMap(P2A_abs, SymmetricIntervalHull(A * U)))
    #P = overapproximate(minkowski_sum(dU, E_ψ), Zonotope) #

    P = linear_map(pis, U)
    #@show P
    d = δ⁻
    while d < δ⁺

        phiDict[d] = copy(ϕ)
        TphiDict[d] = copy(permutedims(ϕ))
        inputDiscritezationDict[d] = copy(P)

        P = plus(P, LazySets.linear_map(phiDict[d], P))
        #@show P


        mul!(tempM, ϕ, ϕ)
        copy!(ϕ, tempM)
        d = d * 2
    end
    #if LazySets.order(P) > maxOrder
    #    P = reduce_order(P, reduceOrder)
    #end
    inputDiscritezationDict[δ⁺] = P
    phiDict[δ⁺] = copy(ϕ)
    TphiDict[d] = copy(permutedims(ϕ))

    return phiDict, TphiDict, inputDiscritezationDict
end

plus(z1::Zonotope, z2::Zonotope) = Zonotope(z1.center + z2.center, z1.generators + z2.generators)