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

    U = isnothing(loc.B) ? (isnothing(loc.u) ? nothing : loc.u) : concretize(loc.B * loc.u)


    d = δ⁻
    dia::Matrix{Float64} = diagm(ones(XDim))
    isInvA = false #isinvertible(A)
    Φ = copy(phiDict[d])
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    Φcache = sum(A) == abs(sum(A)) ? Φ : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)
    pis = ReachabilityAnalysis.Exponentiation.Φ₁(A, δ⁻, alg, isInvA, Φcache)
    if !isnothing(U)
        if !isnothing(loc.c)
            U = Zonotope(U.center + loc.c, genmat(U))
        end
        if !(zeros(XDim) ∈ U) #Origin is *not* in input
            #invA = inv(Matrix(A))
            û = copy(U.center)
            Ut = Zonotope(U.center - û, genmat(U))
            dU = overapproximate(δ⁻ * Ut, Zonotope)
            E_ψ = convert(Zonotope, symmetric_interval_hull(P2A_abs * symmetric_interval_hull(A * U)))
            P = minkowski_sum(dU, E_ψ)
            P̂ = pis * û
            #P̂ = invA * (phiDict[d] - dia) * û
            lt = minkowski_sum(convert(Zonotope, phiDict[d] * X0), dU)
            E⁺ = convert(Zonotope, symmetric_interval_hull(P2A_abs * symmetric_interval_hull(A * A * X0)))
            rt = minkowski_sum(E_ψ, E⁺)
            PZ = Zonotope(P̂, zeros(Float64, size(U.center, 1), 1))
            f = minkowski_sum(lt, rt)
            disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
            while d < δ⁺
                inputDiscritezationDict[d] = P
                P = minkowski_sum(P, linear_map(phiDict[d], P))
                discritezationDict[d] = disc
                if maxOrder > 0
                    if LazySets.order(P) > maxOrder
                        P = reduce_order(P, reduceOrder)
                    end
                    if LazySets.order(disc) > maxOrder
                        disc = reduce_order(disc, reduceOrder)
                    end
                end

                disc = overapproximate(CH(disc, linear_map(phiDict[d], disc)), Zonotope)
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
        else
            dU = overapproximate(d * U, Zonotope)
            E_ψ = convert(Zonotope, symmetric_interval_hull(P2A_abs * symmetric_interval_hull(A * U)))
            P = minkowski_sum(dU, E_ψ)
            E⁺ = convert(Zonotope, symmetric_interval_hull(P2A_abs * symmetric_interval_hull(A * A * X0)))
            lt = concretize(minkowski_sum(convert(Zonotope, phiDict[d] * X0), dU))
            rt = concretize(minkowski_sum(E_ψ, E⁺))
            f = concretize(minkowski_sum(lt, rt))
            disc = overapproximate(CH(X0, f), Zonotope)
            while d < δ⁺
                inputDiscritezationDict[d] = P
                P = minkowski_sum(P, linear_map(phiDict[d], P))
                discritezationDict[d] = copy(disc)
                if maxOrder > 0
                    if LazySets.order(P) > maxOrder
                        P = reduce_order(P, reduceOrder)
                    end
                    if LazySets.order(disc) > maxOrder
                        disc = reduce_order(disc, reduceOrder)
                    end
                end

                disc = overapproximate(CH(disc, linear_map(phiDict[d], disc)), Zonotope)
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
        end
    else
        #   Cases where no input set, but maybe a constant input
        if !isnothing(loc.c)
            #invA = inv(Matrix(A))
            û = loc.c
            println("THIS?")
            P̂ = pis * û
            PZ = Zonotope(P̂, zeros(Float64, XDim, 1))
            E⁺ = convert(Zonotope, symmetric_interval_hull(P2A_abs * symmetric_interval_hull(A * A * X0)))

            f = minkowski_sum(convert(Zonotope, phiDict[d] * X0), E⁺)
            disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
            P = Zonotope(zeros(Float64, XDim), zeros(Float64, XDim, 1))
            while d < δ⁺
                inputDiscritezationDict[d] = P
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
            inputDiscritezationDict[δ⁺] = P
        else    #   When there also is no constant input
            invA = inv(Matrix(A))
            E⁺ = convert(Zonotope, symmetric_interval_hull(P2A_abs * symmetric_interval_hull(A * A * X0)))
            rt = E⁺
            f = minkowski_sum(phiDict[d], rt)
            disc = overapproximate(CH(X0, f), Zonotope)
            P = Zonotope(zeros(Float64, XDim, 1), zeros(Float64, XDim, 1))
            while d < δ⁺
                inputDiscritezationDict[d] = P

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
            inputDiscritezationDict[δ⁺] = P
        end
    end
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
            d = d * 2
        end
        phiDict[δ⁺] = copy(ϕ)
        return phiDict
    end
end