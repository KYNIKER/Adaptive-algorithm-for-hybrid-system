export ReACTDiscretize, PhiDict
using LinearAlgebra, LazySets, ReachabilityAnalysis

isinvertible(x) = applicable(inv, x) && isone(inv(Matrix(x)) * x)

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