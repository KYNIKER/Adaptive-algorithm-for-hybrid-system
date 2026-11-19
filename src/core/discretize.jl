export PhiDict
using LinearAlgebra, LazySets, ReachabilityAnalysis
using ReachabilityBase.Arrays: isinvertible

function ReACTDiscretizePlus(loc, X0::Zonotope{N,Vector{N},Matrix{N}}, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5, phiDict=nothing) where {N}

    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    inputDiscritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    A = loc.A
    XDim = size(A, 1)
    if isnothing(phiDict)
        phiDict = PhiDict(A, δ⁻, δ⁺, alg)
    end

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

    inputDiscritezationDict[0] = U
    dU = overapproximate(LinearMap(δ⁻, U), Zonotope)#linear_map(dia, U)
    E_ψ = SymmetricIntervalHull(linear_map(P2A_abs, SymmetricIntervalHull(LazySets.linear_map(A, U))))
    P = minkowski_sum(dU, E_ψ) #
    lt = minkowski_sum(LazySets.linear_map(phiDict[d], X0), dU)  #minkowski_sum(convert(Zonotope, phiDict[d] * X0), dU)
    E⁺ = SymmetricIntervalHull(LazySets.linear_map(P2A_abs, SymmetricIntervalHull(LazySets.linear_map(A * A, X0))))
    rt = minkowski_sum(E_ψ, E⁺)
    f = minkowski_sum(lt, rt)

    disc = overapproximate(CH(X0, f), Zonotope) # overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope) #

    # TODO maybe we can reuse this somehow?

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
        disc = overapproximate(CH(disc, minkowski_sum(P, LazySets.linear_map(phiDict[d], disc))), Zonotope)
        P = minkowski_sum(P, LazySets.linear_map(phiDict[d], P))
        d = d * 2
    end
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
    discritezationDict[δ⁺] = disc

    return discritezationDict
end

# overapproximate(CH(Z1, Z2), Zonotope) for two zonotopes, written out: the same arithmetic as LazySets' default ("mean",
# Girard 2005) algorithm, so the same zonotope, without its temporaries. Always returns a new zonotope.
function zonotope_hull(Z1::Zonotope, Z2::Zonotope)
    m1, m2 = ngens(Z1), ngens(Z2)
    if m1 < m2
        Z2 = reduce_order(Z2, max(1, m1))
    elseif m1 > m2
        Z1 = reduce_order(Z1, max(1, m2))
    end
    if ngens(Z2) > ngens(Z1)
        Z1, Z2 = Z2, Z1
    end
    c1, c2 = LazySets.center(Z1), LazySets.center(Z2)
    G1, G2 = genmat(Z1), genmat(Z2)
    n, k1, k2 = length(c1), size(G1, 2), size(G2, 2)
    G = Matrix{Float64}(undef, n, k1 + k2 + 1)
    col = 0
    # Columns in LazySets' order: (G1 + G2) / 2, (c1 - c2) / 2, (G1 - G2) / 2 over the first k2 generators, then the rest of G1. Zero columns are dropped.
    for part in 1:4, j in (part == 2 ? (1:1) : part == 4 ? (k2+1:k1) : (1:k2))
        nonzero = false
        @inbounds for i in 1:n
            x = part == 1 ? (G1[i, j] + G2[i, j]) / 2 : part == 2 ? (c1[i] - c2[i]) / 2 : part == 3 ? (G1[i, j] - G2[i, j]) / 2 : G1[i, j]
            G[i, col+1] = x
            nonzero |= !iszero(x)
        end
        col += nonzero
    end
    return Zonotope((c1 + c2) / 2, col == size(G, 2) ? G : G[:, 1:col])
end

function newReACTDiscretizePlus(X0::Zonotope{N,Vector{N},Matrix{N}}, δ⁻::Float64, δ⁺::Float64, A, P2A_absp, phiDict, inputDiscritezationDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    d = copy(δ⁻)
    discritezationDictp = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    lt = minkowski_sum(linear_map(phiDict[d], X0), inputDiscritezationDict[0])  #minkowski_sum(convert(Zonotope, phiDict[d] * X0), dU)
    E⁺ = SymmetricIntervalHull(linear_map(P2A_absp, SymmetricIntervalHull(linear_map(A * A, X0))))
    f = minkowski_sum(lt, E⁺)
    #dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    disc = f isa Zonotope ? zonotope_hull(X0, f) : overapproximate(CH(X0, f), Zonotope; algorithm="mean") # overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope) #
    while d < δ⁺
        discritezationDictp[d] = disc # zonotope_hull below makes a new zonotope, so disc is not changed after this
        rs = minkowski_sum(linear_map(phiDict[d], disc), inputDiscritezationDict[d])
        disc = zonotope_hull(disc, rs)

        d = d * 2
    end
    discritezationDictp[δ⁺] = disc

    return discritezationDictp
end

function ReACT_discretize_decomposed_generators(X0G::Zonotope{N,Vector{N},Matrix{N}}, δ⁻::Float64, δ⁺::Float64, A, P2A_abs, phiDict, U, inputDiscritezationDict) where {N}
    d = copy(δ⁻)
    discritezationDict = Dict{Float64,Matrix{N}}()
    UG = genmat(U) .* d
    XG = genmat(X0G)

    U_bloat = diagm(sum(abs.(genmat(inputDiscritezationDict[-1])), dims=2)[:, 1]) #Psi #diagm(sum(P2A_abs * diagm((AUG+AUc)[:, 1]), dims=2)[:, 1])
    G_bloat = diagm(sum(P2A_abs * diagm(sum(abs.(A * A * XG), dims=2)[:, 1]), dims=2)[:, 1])
    newG = hcat((XG .+ (phiDict[d] * XG)) .* 0.5, (XG .- (phiDict[d] * XG)) .* 0.5, UG)
    discritezationDict[0] = copy(U_bloat + G_bloat)
    discritezationDict[-1] = copy(U_bloat)
    # Some mistake in here
    while d < δ⁺
        discritezationDict[d] = copy(newG)
        #newG´ = copy(newG)
        newG´ = copy(phiDict[d] * newG)
        newG = hcat((newG .+ newG´) .* 0.5, (newG .- newG´) .* 0.5, genmat(inputDiscritezationDict[d]))

        d = d * 2
    end

    discritezationDict[δ⁺] = newG

    return discritezationDict
end

### TODO - Figure out whether this is easier if we had passed the decomposed components instead?
function ReACT_discretize_combine_with_offsets(X0c::Zonotope{N,Vector{N},Matrix{N}}, δ⁻::Float64, δ⁺::Float64, A, P2A_abs, phiDict, U, inputDiscritezationDict, generatorDiscretizationDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    d = copy(δ⁻)
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()

    Xc = X0c.center
    Uc = inputDiscritezationDict[0].center

    eAXc = (phiDict[d] * Xc)
    newC = (Xc .+ eAXc .+ Uc) .* 0.5

    c_bloat = hcat(diagm(sum(P2A_abs * diagm(abs.(A * A * Xc)), dims=2)[:, 1]) + generatorDiscretizationDict[0], (Xc .- (eAXc .+ Uc)) .* 0.5)
    while d < δ⁺
        discritezationDict[d] = copy(Zonotope(newC, hcat(c_bloat, generatorDiscretizationDict[d])))
        eA_c_bloat = phiDict[d] * c_bloat
        c_bloat = hcat(c_bloat .+ eA_c_bloat, c_bloat .- eA_c_bloat, (newC .- (phiDict[d] * newC) .- inputDiscritezationDict[d].center)) .* 0.5
        newC = (newC .+ (phiDict[d] * newC) .+ inputDiscritezationDict[d].center) .* 0.5

        d = d * 2
    end
    discritezationDict[δ⁺] = copy(Zonotope(newC, hcat(c_bloat, generatorDiscretizationDict[δ⁺])))

    return discritezationDict, size(c_bloat, 2)
end

function ReACT_discretize_combine_with_offset_vector!(preallocated_center::Vector{N}, preallocated_genmat::Matrix{N}, preallocated_inter_genmat::Matrix{N}, discritezationDict, X0c::Vector{N}, δ⁻::Float64, δ⁺::Float64, A, P2A_abs, phiDict, U, inputDiscritezationDict, generatorDiscretizationDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    d = copy(δ⁻)

    Xc = X0c
    Uc = inputDiscritezationDict[-2].center # dU + Psi

    eAXc = (phiDict[d] * X0c)
    copy!(view(preallocated_center, :, :), view((X0c + eAXc + Uc) .* 0.5, :, :))
    dim = length(X0c)
    ng = dim + 1

    copy!(view(preallocated_genmat, :, 1:(ng-1)), view(diagm(sum(P2A_abs * diagm(abs.(A * A * Xc)), dims=2)[:, 1]) + generatorDiscretizationDict[0], :, :))

    copy!(view(preallocated_genmat, :, ng), (Xc - (eAXc .+ Uc)) .* 0.5)

    while d < δ⁺
        copy!(discritezationDict[d].center, preallocated_center)

        copy!(view(discritezationDict[d].generators, :, 1:(ng)), view(preallocated_genmat, :, 1:ng))

        copyto!(view(discritezationDict[d].generators, :, (ng+1):(size(discritezationDict[d].generators, 2))), generatorDiscretizationDict[d])#view(generatorDiscretizationDict[d], :, 1:(size(generatorDiscretizationDict[d], 2))))

        copyto!(view(preallocated_inter_genmat, :, 1:ng), (view(preallocated_genmat, :, 1:ng) + phiDict[d] * view(preallocated_genmat, :, 1:ng)) .* 0.5)

        copyto!(view(preallocated_genmat, :, (ng+1):(2*ng+1)), (view(preallocated_genmat, :, 1:ng) - phiDict[d] * view(preallocated_genmat, :, 1:ng)) .* 0.5)

        copyto!(view(preallocated_genmat, :, 1:ng), view(preallocated_inter_genmat, :, 1:ng))

        ng = 2 * ng + 1
        copy!(view(preallocated_genmat, :, ng), (preallocated_center .- ((phiDict[d] * preallocated_center) + inputDiscritezationDict[d].center)) .* 0.5)
        copy!(preallocated_center, (preallocated_center .+ (phiDict[d] * preallocated_center) .+ inputDiscritezationDict[d].center) .* 0.5)

        d = d * 2
    end
    copy!(discritezationDict[d].center, preallocated_center)
    copy!(view(discritezationDict[d].generators, :, 1:(ng)), view(preallocated_genmat, :, 1:ng))
    copyto!(view(discritezationDict[d].generators, :, (ng+1):(size(discritezationDict[d].generators, 2))), generatorDiscretizationDict[d])

end

function ReACT_discretize_state!(preallocated_center::Vector{N}, preallocated_genmat::Matrix{N}, preallocated_inter_genmat::Matrix{N}, discritezationStateDict, X0c::Vector{N}, δ⁻::Float64, δ⁺::Float64, A, P2A_abs, phiDict, U, inputDiscritezationDict, generatorDiscretizationDict, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    d = copy(δ⁻)

    Xc = X0c
    Uc = inputDiscritezationDict[-2].center # dU + Psi

    eAXc = (phiDict[d] * X0c)

    copy!(view(preallocated_center, :, :), view((X0c + eAXc + Uc) .* 0.5, :, :))

    dim = length(X0c)
    gnum_input = size(genmat(U), 1)
    ng = dim + 1

    copy!(view(preallocated_genmat, :, 1:(ng-1)), view(diagm(sum(P2A_abs * diagm(abs.(A * A * Xc)), dims=2)[:, 1]) + generatorDiscretizationDict[-1], :, :)) #calculate new generators from symmetric_interval_hull

    copy!(view(preallocated_genmat, :, ng), (Xc - (eAXc .+ Uc)) .* 0.5) # store generator for halfway point between state and zonotope

    # Overview:
    # Copy center and generators from preallocated matrices to Dict
    # Update the preallocated matrices like before and remember to account for the added difference of centers and the impact of input.
    if gnum_input > 0
        while d < δ⁺

            copy!(discritezationStateDict[d].center, preallocated_center)

            copy!(view(discritezationStateDict[d].generators, :, 1:(ng)), view(preallocated_genmat, :, 1:ng)) # Update generators with bloating factors and halfway point. 

            copyto!(view(preallocated_inter_genmat, :, 1:ng), (view(preallocated_genmat, :, 1:ng) + phiDict[d] * view(preallocated_genmat, :, 1:ng)) .* 0.5) # Store first set of new generators in intermediate matrix

            copyto!(view(preallocated_genmat, :, (ng+1):(2*ng+1)), (view(preallocated_genmat, :, 1:ng) - phiDict[d] * view(preallocated_genmat, :, 1:ng)) .* 0.5) # Store second set of new generators

            copyto!(view(preallocated_genmat, :, 1:ng), view(preallocated_inter_genmat, :, 1:ng)) # copy from intermediate matrix 

            ng = 2 * ng + 1 + gnum_input
            copy!(view(preallocated_genmat, :, ng-gnum_input), (preallocated_center .- ((phiDict[d] * preallocated_center) + inputDiscritezationDict[d].center)) .* 0.5)
            copy!(view(preallocated_genmat, :, (ng-gnum_input+1):ng), genmat(inputDiscritezationDict[d]))

            copy!(preallocated_center, (preallocated_center .+ (phiDict[d] * preallocated_center) .+ inputDiscritezationDict[d].center) .* 0.5)

            d = d * 2
        end
    else
        while d < δ⁺

            copy!(discritezationStateDict[d].center, preallocated_center)

            copy!(view(discritezationStateDict[d].generators, :, 1:(ng)), view(preallocated_genmat, :, 1:ng)) # Update generators with bloating factors and halfway point. 

            copyto!(view(preallocated_inter_genmat, :, 1:ng), (view(preallocated_genmat, :, 1:ng) + phiDict[d] * view(preallocated_genmat, :, 1:ng)) .* 0.5) # Store first set of new generators in intermediate matrix

            copyto!(view(preallocated_genmat, :, (ng+1):(2*ng+1)), (view(preallocated_genmat, :, 1:ng) - phiDict[d] * view(preallocated_genmat, :, 1:ng)) .* 0.5) # Store second set of new generators

            copyto!(view(preallocated_genmat, :, 1:ng), view(preallocated_inter_genmat, :, 1:ng)) # copy from intermediate matrix 

            ng = 2 * ng + 1
            copy!(view(preallocated_genmat, :, ng), (preallocated_center .- ((phiDict[d] * preallocated_center) + inputDiscritezationDict[d].center)) .* 0.5)

            copy!(preallocated_center, (preallocated_center .+ (phiDict[d] * preallocated_center) .+ inputDiscritezationDict[d].center) .* 0.5)
            d = d * 2
        end

    end
    copy!(discritezationStateDict[d].center, preallocated_center)
    copy!(view(discritezationStateDict[d].generators, :, 1:(ng)), view(preallocated_genmat, :, 1:ng))
    @show ng

end

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

function PhiInputDict(loc, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5)
    A = loc.A
    ϕ::Matrix{Float64} = ReachabilityAnalysis.Exponentiation._exp(A, δ⁻, alg)
    phiDict = Dict{Float64,Matrix{Float64}}()
    TphiDict = Dict{Float64,Matrix{Float64}}()
    tempM = similar(ϕ)

    inputDiscritezationDict = Dict{Float64,Zonotope}()
    XDim = size(A, 1)

    U = Zonotope(zeros(XDim), zeros(XDim, 1))
    U = isnothing(loc.B) ? (isnothing(loc.u) ? Zonotope(zeros(XDim), zeros(XDim, 1)) : loc.u) : linear_map(loc.B, loc.u)
    if !isnothing(loc.c)
        U = Zonotope(U.center + loc.c, genmat(U))
    end

    d = δ⁻
    #Φ = copy(ϕ)
    dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    isInvA = isinvertible(A)
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    Φcache = A == A_abs ? ϕ : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)

    inputDiscritezationDict[0] = U
    dU = linear_map(dia, U)#LinearMap(δ⁻, U)#
    E_ψ = symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(LazySets.linear_map(A, U))))
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
    inputDiscritezationDict[δ⁺] = P
    phiDict[δ⁺] = copy(ϕ)
    TphiDict[d] = copy(permutedims(ϕ))

    return phiDict, TphiDict, inputDiscritezationDict
end

# The one used in shielding
function phi_input_dict(A, U, δ⁻::Float64, δ⁺::Float64, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5)

    ϕ::Matrix{Float64} = ReachabilityAnalysis.Exponentiation._exp(A, δ⁻, alg)
    phiDict = Dict{Float64,Matrix{Float64}}()
    TphiDict = Dict{Float64,Matrix{Float64}}()
    tempM = similar(ϕ)

    inputDiscritezationDict = Dict{Float64,Zonotope}()
    XDim = size(A, 1)

    d = δ⁻
    #Φ = copy(ϕ)
    dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    isInvA = isinvertible(A)
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    Φcache = A == A_abs ? ϕ : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)
    pis = ReachabilityAnalysis.Exponentiation.Φ₁(A, δ⁻, alg, false, nothing)

    dU = linear_map(dia, U)#LinearMap(δ⁻, U)#
    E_ψ = symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(LazySets.linear_map(A, U))))
    inputDiscritezationDict[0] = minkowski_sum(dU, E_ψ)
    inputDiscritezationDict[-1] = E_ψ
    inputDiscritezationDict[-2] = dU

    # The input's contribution over δ⁻ and its doublings. u(t) ∈ U may vary in time, so only U's centre c is constant: its
    # part is exact (Φ₁ c, and c_2d = c_d + Φ_d c_d). The part of U's generators is bounded like dU ⊕ E_ψ above, and the
    # input in [d, 2d] is independent of the input in [0, d], so the doubling is a Minkowski sum, P_d ⊕ Φ_d P_d, reduced
    # to reduceOrder above maxOrder. Without generators the input is constant and P is Φ₁ U, composed with plus.
    G = genmat(U)
    constant = iszero(G)
    P = constant ? linear_map(pis, U) : Zonotope(pis * LazySets.center(U), hcat(δ⁻ .* G, Matrix(Diagonal(P2A_abs * vec(sum(abs, A * G; dims=2))))))
    d = δ⁻
    while d < δ⁺

        phiDict[d] = copy(ϕ)
        TphiDict[d] = copy(permutedims(ϕ))
        inputDiscritezationDict[d] = copy(P)

        if constant
            P = plus(P, LazySets.linear_map(phiDict[d], P))
        else
            P = Zonotope(LazySets.center(P) + phiDict[d] * LazySets.center(P), hcat(genmat(P), phiDict[d] * genmat(P)))
            if maxOrder > 0 && LazySets.order(P) > maxOrder
                P = reduce_order(P, reduceOrder)
                P = Zonotope(Vector{Float64}(LazySets.center(P)), Matrix{Float64}(genmat(P)))
            end
        end

        mul!(tempM, ϕ, ϕ)
        copy!(ϕ, tempM)
        d = d * 2
    end
    inputDiscritezationDict[δ⁺] = P
    phiDict[δ⁺] = copy(ϕ)
    TphiDict[d] = copy(permutedims(ϕ))

    return phiDict, TphiDict, inputDiscritezationDict
end

# Bound on the input's contribution over [0, t], ∫₀ᵗ e^{Ar} u(t - r) dr, when u ∈ U may vary in time. U's centre c gives
# Φ₁(A, t) c exactly. For the rest, w = u - c ∈ U₀ (U's generators G), split e^{Ar} = I + (e^{Ar} - I): ∫ w dr lies in
# t U₀ by convexity, and since e^{Ar} - I = Σ_k A^k r^{k+1}/(k+1)! A, the remainder is at most Φ₂(|A|, t) |A G| 1
# elementwise, with Φ₂(M, t) = Σ_k M^k t^{k+2}/(k+2)! (Le Guernic and Girard's E_ψ). Without generators the input is
# constant and the bound is Φ₁(A, t) U.
function input_bound(A, U, t, alg)
    P = Matrix{Float64}(ReachabilityAnalysis.Exponentiation.Φ₁(A, t, alg, false, nothing))
    G = Matrix{Float64}(genmat(U))
    c = P * LazySets.center(U)
    iszero(G) && return Zonotope(c, P * G)
    P2_abs = Matrix{Float64}(ReachabilityAnalysis.Exponentiation.Φ₂(ReachabilityAnalysis.Exponentiation.elementwise_abs(A), t, alg, false, nothing))
    return Zonotope(c, hcat(t .* G, Matrix(Diagonal(P2_abs * vec(sum(abs, A * G; dims=2))))))
end

plus(z1::Zonotope, z2::Zonotope) = Zonotope(z1.center + z2.center, z1.generators + z2.generators)