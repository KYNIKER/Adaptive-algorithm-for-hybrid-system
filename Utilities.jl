using LazySets, ReachabilityAnalysis

export HybridSystem, HybridSystemV2, Location, Edge, overapproximateIntervalReachset, intersects

#=struct HybridSystem
    V::Vector{Any}
    E::Matrix{Int}
    G::Matrix{Vector{LazySet.HalfSpaceModule.HalfSpace}}
    I::Vector{LazySet.HalfSpaceModule.HalfSpace}
    Flow::Vector{Matrix{Float64}}
    Jump::Matrix{Tuple{Matrix{Float64},Vector{Float64}}}
end=#



struct Edge
    targetLoc::Int
    guard::HPolyhedron
    jumpMatrix::Matrix{Float64}
    jumpVector::Vector{Float64}
end

Base.show(io::Core.IO, e::Edge) = print(io, "Edge going to: ", e.targetLoc)


struct Location
    id::Int
    invarient::Union{HPolyhedron,Nothing}
    A::Matrix{Float64}
    edges::Vector{Edge}
end

Base.show(io::Core.IO, l::Location) = print(io, "Location: ", l.id, "\n invariant? ", !isnothing(l.invarient), "\n edges: ", length(l.edges))


mutable struct HybridSystemV2
    locations::Vector{Location}
    initialLoc::Int
    initialState # Fill this in later
end

Base.show(io::Core.IO, s::HybridSystemV2) = print(io, "System with ", length(s.locations), " locations.")

struct HybridSystem
    V
    E
    G
    I
    Flow
    Jump
    Init
end

function system(locations, edges, guards, invariants, flows, jumps, init)
    n = length(locations)

    E = zeros(Int8, n, n)
    for (i, j) in edges
        E[i, j] = 1
    end

    G = Matrix(missing, n, n)
    for (l, k, g) in guards
        if (G[l, k] == missing)
            G[l, k] = [g]
        else
            push!(G[l, k], g)
        end
    end

    J = Matrix(missing, n, n)
    for (l, k, j) in jumps
        if (J[l, k] == missing)
            J[l, k] = [j]
        else
            push!(J[l, k], j)
        end
    end
    return HybridSystem(locations, E, G, invariants, flows, J, init)
end

function intersection(Z, H)
    agenSum = reduce(+, abs.(genmat(Z) * H.a))
    acenSum = Z.center * H.a
    if (acenSum - agenSum <= H.b) & (H.b <= acenSum + agenSum)
        return Z
    else
        return nothing
    end
end

function intersects(Z, H)
    agenSum = reduce(+, abs.(genmat(Z) .* H.a))
    acenSum = dot(Vector(H.a), Z.center)
    return (acenSum - agenSum <= H.b) & (H.b <= acenSum + agenSum)
end

function overapproximateIntervalReachset(A, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, δ⁻, δ⁺, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, phiDict=nothing) where {N}
    XDim, _ = size(genmat(X0))
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    inputDiscritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()

    if isnothing(phiDict)
        phiDict = PhiDict(A, δ⁻, δ⁺, alg)
    end

    U = concretize(U)


    d = δ⁻
    dia::Matrix{Float64} = diagm(ones(XDim))
    isInvA = isinvertible(A)
    Φ = copy(phiDict[d])
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    Φcache = sum(A) == abs(sum(A)) ? Φ : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)

    if !(zeros(XDim) ∈ U) #Origin is *not* in input
        invA = inv(Matrix(A))
        û = copy(U.center)
        Ut = Zonotope(U.center - û, genmat(U))
        dU = overapproximate(δ⁻ * Ut, Zonotope)
        E_ψ = convert(Zonotope, symmetric_interval_hull(P2A_abs * symmetric_interval_hull(A * U)))
        P = minkowski_sum(dU, E_ψ)
        P̂ = invA * (phiDict[d] - dia) * û
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
        discritezationDict[d] = copy(disc)
        inputDiscritezationDict[d] = P
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

            disc = overapproximate(CH(disc, linear_map(ϕ, disc)), Zonotope)
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
        discritezationDict[d] = copy(disc)
        inputDiscritezationDict[d] = P
    end

    return discritezationDict, inputDiscritezationDict, phiDict
end