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
    isInvA = false
    Φ = copy(phiDict[d])
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    Φcache = sum(A) == abs(sum(A)) ? Φ : nothing
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)

    inputDiscritezationDict[0] = U
    dU = overapproximate(LinearMap(δ⁻, U), Zonotope)
    E_ψ = SymmetricIntervalHull(linear_map(P2A_abs, SymmetricIntervalHull(LazySets.linear_map(A, U))))
    P = minkowski_sum(dU, E_ψ)
    lt = minkowski_sum(LazySets.linear_map(phiDict[d], X0), dU)
    E⁺ = SymmetricIntervalHull(LazySets.linear_map(P2A_abs, SymmetricIntervalHull(LazySets.linear_map(A * A, X0))))
    rt = minkowski_sum(E_ψ, E⁺)
    f = minkowski_sum(lt, rt)

    disc = overapproximate(CH(X0, f), Zonotope)

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

    d = δ⁻
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, false, nothing)

    lt = minkowski_sum(LazySets.linear_map(phiDict[d], X0), inputDiscritezationDict[d])
    E⁺ = SymmetricIntervalHull(LazySets.linear_map(P2A_abs, SymmetricIntervalHull(LazySets.linear_map(A * A, X0))))
    f = minkowski_sum(lt, E⁺)

    disc = overapproximate(CH(X0, f), Zonotope)

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

# overapproximate(CH(Z1, Z2), Zonotope) for two zonotopes, written out: the same arithmetic as LazySets' default "mean".
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
    for part in 1:4, j in (part == 2 ? (1:1) : part == 4 ? ((k2+1):k1) : (1:k2))
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
    lt = minkowski_sum(linear_map(phiDict[d], X0), inputDiscritezationDict[0])
    E⁺ = SymmetricIntervalHull(linear_map(P2A_absp, SymmetricIntervalHull(linear_map(A * A, X0))))
    f = minkowski_sum(lt, E⁺)
    disc = f isa Zonotope ? zonotope_hull(X0, f) : overapproximate(CH(X0, f), Zonotope; algorithm="mean")
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

    U_bloat = diagm(sum(abs.(genmat(inputDiscritezationDict[-1])), dims=2)[:, 1])
    G_bloat = diagm(sum(P2A_abs * diagm(sum(abs.(A * A * XG), dims=2)[:, 1]), dims=2)[:, 1])
    newG = hcat((XG .+ (phiDict[d] * XG)) .* 0.5, (XG .- (phiDict[d] * XG)) .* 0.5, UG)
    discritezationDict[0] = copy(U_bloat + G_bloat)
    discritezationDict[-1] = copy(U_bloat)
    while d < δ⁺
        discritezationDict[d] = copy(newG)
        newG´ = copy(phiDict[d] * newG)
        newG = hcat((newG .+ newG´) .* 0.5, (newG .- newG´) .* 0.5, genmat(inputDiscritezationDict[d]))

        d = d * 2
    end

    discritezationDict[δ⁺] = newG

    return discritezationDict
end

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

# The discretization of a cell zonotope Z at every level δ⁻, 2δ⁻, ..., δ⁺, for any centre Xc: the same sets as newReACTDiscretizePlus on Z translated to Xc, without its order reduction. 
struct CellDiscretization{N}
    dict::Dict{Float64,Zonotope{Float64,Vector{Float64},Matrix{Float64}}}
    sets::Vector{Zonotope{Float64,Vector{Float64},Matrix{Float64}}} # per level, the same zonotopes as dict
    Kc::Vector{Matrix{Float64}} # per level, n × (2n + 1): centre = Kc * w
    KG::Vector{Matrix{Float64}} # per level, (n m) × (2n + 1): the m centre-dependent generators, column after column
    AA::Matrix{Float64}
    P2A_abs::Matrix{Float64}
    w::Vector{Float64}
    a::Vector{Float64} # |A² Xc|
end

# Merges parallel generators, given as vectors or as the maps that give them: g and α g span the same segment as
# (1 + |α|) g, so the zonotope is the same. Parallel means equal up to rounding: what is left of one after taking out its
# part along the other is at most PARALLEL_TOL of its length. Zero generators are dropped.
const PARALLEL_TOL = 1e-12
function merge_parallel(vs::AbstractVector{<:AbstractArray{Float64}})
    reps, units, scales = eltype(vs)[], Vector{Float64}[], Float64[]
    for v in vs
        x = vec(v)
        nx = norm(x)
        nx == 0 && continue
        i = findfirst(u -> norm(x - dot(u, x) * u) <= PARALLEL_TOL * nx, units)
        if i === nothing
            push!(reps, v)
            push!(units, x / nx)
            push!(scales, 1.0)
        else
            scales[i] += nx / norm(reps[i]) # |α|
        end
    end
    return [s * r for (r, s) in zip(reps, scales)]
end

function CellDiscretization(Z::Zonotope, reach, ::Val{N}) where N
    (; δ⁻, δ⁺, phiDict, inputDict) = reach
    A, P2A_abs = Matrix{Float64}(reach.A), Matrix{Float64}(reach.P2A_abs)
    shared = ReACT_discretize_decomposed_generators(Z, δ⁻, δ⁺, A, P2A_abs, phiDict, reach.U, inputDict)
    n, nw = N, 2N + 1
    Id = Matrix{Float64}(I, n, n)
    lift(M, v) = hcat(M, v, zeros(n, n)) # the map w -> M Xc + v
    # Level δ⁻: the hull of Z and Φ Z ⊕ input. The first n generators are the bloat, diag(s + b), with b the bloat of the
    # generators and the input (shared[0], diagonal); then (Xc - Φ Xc - Uc) / 2.
    Φ, Uc, b = phiDict[δ⁻], inputDict[-2].center, diag(shared[0])
    Kc = lift((Id + Φ) / 2, Uc / 2)
    # s_j is 0 for every centre when row j of P2A_abs |A²| is 0 (for instance when A² = 0): it then gets no entry
    s_used = [any(!iszero, row) for row in eachrow(P2A_abs * abs.(A * A))]
    cols = Matrix{Float64}[]
    for j in 1:n
        K = zeros(n, nw)
        K[j, n+1] = b[j]
        s_used[j] && (K[j, n+1+j] = 1.0)
        push!(cols, K)
    end
    push!(cols, lift((Id - Φ) / 2, -Uc / 2))
    dict = Dict{Float64,Zonotope{Float64,Vector{Float64},Matrix{Float64}}}()
    sets, Kcs, KGs = Zonotope{Float64,Vector{Float64},Matrix{Float64}}[], Matrix{Float64}[], Matrix{Float64}[]
    d = δ⁻
    while true
        # A zero map gives a zero generator for every centre, and parallel maps parallel generators: both stay so in the
        # generators made from them, so merging here also saves the work of the levels above.
        cols = merge_parallel(cols)
        shared_cols = merge_parallel(collect(eachcol(shared[d])))
        S = isempty(shared_cols) ? zeros(n, 0) : reduce(hcat, shared_cols)
        dict[d] = Zonotope(zeros(n), hcat(zeros(n, length(cols)), S))
        push!(sets, dict[d])
        push!(Kcs, copy(Kc))
        push!(KGs, isempty(cols) ? zeros(0, nw) : reduce(vcat, cols))
        d < δ⁺ || break
        # Next level: the hull of the set and Φ times it ⊕ the input, as in ReACT_discretize_combine_with_offsets
        Φ, u = phiDict[d], inputDict[d].center
        P, M = (Id + Φ) / 2, (Id - Φ) / 2
        cols = vcat([P * K for K in cols], [M * K for K in cols], [M * Kc - lift(zeros(n, n), u / 2)])
        Kc = P * Kc + lift(zeros(n, n), u / 2)
        d *= 2
    end
    return CellDiscretization{N}(dict, sets, Kcs, KGs, A * A, P2A_abs, zeros(nw), zeros(n))
end

function discretize_cell!(cd::CellDiscretization{N}, Xc::AbstractVector{Float64}) where N
    (; w, a) = cd
    mul!(a, cd.AA, Xc)
    a .= abs.(a)
    @inbounds for i in 1:N
        w[i] = Xc[i]
    end
    w[N+1] = 1.0
    mul!(view(w, (N+2):(2N+1)), cd.P2A_abs, a)
    @inbounds for (ℓ, Z) in enumerate(cd.sets)
        mul!(Z.center, cd.Kc[ℓ], w)
        KG = cd.KG[ℓ]
        mul!(view(vec(Z.generators), 1:size(KG, 1)), KG, w)
    end
    return cd.dict
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
    dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, false, nothing)

    inputDiscritezationDict[0] = U
    dU = linear_map(dia, U)
    E_ψ = symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(LazySets.linear_map(A, U))))
    P = minkowski_sum(dU, E_ψ)

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
    dia::Matrix{Float64} = diagm(δ⁻ * ones(XDim))
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, false, nothing)
    pis = ReachabilityAnalysis.Exponentiation.Φ₁(A, δ⁻, alg, false, nothing)

    dU = linear_map(dia, U)
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

# Without generators the input is constant and the bound is Φ₁(A, t) U. P is Φ₁(A, t) and P2_abs is Φ₂(|A|, t), which is not used without generators.
function input_bound(A, U, t, P::Matrix{Float64}, P2_abs::Matrix{Float64})
    G = Matrix{Float64}(genmat(U))
    c = P * LazySets.center(U)
    iszero(G) && return Zonotope(c, P * G)
    return Zonotope(c, hcat(t .* G, Matrix(Diagonal(P2_abs * vec(sum(abs, A * G; dims=2))))))
end

plus(z1::Zonotope, z2::Zonotope) = Zonotope(z1.center + z2.center, z1.generators + z2.generators)