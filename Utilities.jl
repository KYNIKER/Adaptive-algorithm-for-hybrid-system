using LazySets, ReachabilityAnalysis, LinearAlgebra, Polyhedra, Optim

export HybridSystem, HybridSystemV2, Location, Edge, overapproximateIntervalReachset, intersects, splitZonotope, getBoxIntersection

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
    B::Union{Nothing,Matrix{Float64}}
    u # Unsure 
    c::Union{Nothing,Vector{Float64}}
    edges::Vector{Edge}
    constraints::Vector{Union{HPolyhedron,LazySets.HalfSpace}}
end


Base.show(io::Core.IO, l::Location) = print(io, "Location: ", l.id, "\n invariant? ", !isnothing(l.invarient), "\n edges: ", length(l.edges))


mutable struct HybridSystemV2
    locations::Vector{Location}
    globalConstraints::Vector{Union{HPolyhedron,LazySets.HalfSpace}}
    #initialLoc::Int
    #initialState # Fill this in later
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

#=function intersection(Z, H)
    return Z, H
    #=agenSum = reduce(+, abs.(genmat(Z) .* H.a))
    acenSum = dot(Vector(H.a), Z.center)
    if (acenSum - agenSum <= H.b) & (H.b <= acenSum + agenSum)
    else
        return nothing
    end=#
end=#

function intersects(Z::Zonotope, H::LazySets.HalfSpace)
    agenSum = reduce(+, abs.(genmat(Z) .* H.a))
    acenSum = dot(Vector(H.a), Z.center)
    return (acenSum - agenSum <= H.b) & (H.b <= acenSum + agenSum)
end

function intersects(Z::Zonotope, H::Vector{N}) where N
    sen = true
    for h in H
        agenSum = reduce(+, abs.(genmat(Z) .* h.a))
        acenSum = dot(Vector(h.a), Z.center)
        sen = sen & (acenSum - agenSum <= h.b) & (h.b <= acenSum + agenSum)
    end
    return sen
end

function intersects(Z::Zonotope, H::Any)
    return !isempty(∩(Z, H))
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






function produceDirections(n)
    # Produce one hot encodings. 
    #TODO: Could be smarter directions?
    directions = Vector{Vector{Float64}}()
    for i in 1:n
        oneHotEncoding = zeros(n)
        oneHotEncoding[i] = 1.
        push!(directions, oneHotEncoding)
    end
    return directions
end

function project2d(Z::Zonotope, n::Vector, direction::Vector)
    c = Z.center
    G = genmat(Z)

    newC = [dot(c, n), dot(c, direction)]
    newG = [[dot(g, n), dot(g, direction)] for g in eachcol(G)]

    return Zonotope(newC, newG)
end

function segmentLineIntersection(segment1, segment2, line)
    x1, y1 = segment1
    x2, y2 = segment2

    if (x1 - line) * (x2 - line) > 0 || (x1 == x2)
        return nothing
    end
    # Use interpolation
    t = (line - x1) / (x2 - x1)
    return y1 + t * (y2 - y1)
end

function GirardGuernicAlgorithm(Z::Zonotope, line::Float64)
    # We assume it is 2 dimensional
    c = Z.center
    G = [copy(g) for g in eachcol(genmat(Z))]

    # Ensure all generators are positive
    for g in G
        if g[2] < 0 || (g[2] == 0 && g[1] < 0)
            g .*= -1
        end
    end

    # Get the lowest vertex
    P = copy(c)
    for g in G
        P .-= g
    end

    # Sort generators in trigometric order
    sort!(G, by=g -> atan(g[2], g[1]))

    m = Inf
    M = -Inf
    # Next we traverse the verticies of the zonotope, 
    # looking for intersections with the line
    # Forward scan
    for g in G
        Q = P .+ 2g
        y = segmentLineIntersection(P, Q, line)
        if y !== nothing
            m = min(m, y)
            M = max(M, y)
        end
        P = Q
    end

    # Backward scan
    for g in G
        Q = P .- 2g
        y = segmentLineIntersection(P, Q, line)
        if y !== nothing
            m = min(m, y)
            M = max(M, y)
        end
        P = Q
    end

    return m, M
end


# Zonotope intersection with hyperplane
# Gets the area where we meet the intersection.
function lineIntersection(Z::Zonotope, H::LazySets.HalfSpace)
    directions = produceDirections(length(Z.center))
    val = H.b
    n = H.a

    constraints = LazySets.HalfSpace[]

    for direction in directions
        Z2 = project2d(Z, n, direction)
        m, M = GirardGuernicAlgorithm(Z2, val)

        #Produce the constraints
        if isfinite(m)
            push!(constraints, LazySets.HalfSpace(direction, M))
            push!(constraints, LazySets.HalfSpace(-direction, -m))
        end
    end

    # enforce hyperplane equality n⋅x = γ
    push!(constraints, LazySets.HalfSpace(n, val))
    push!(constraints, LazySets.HalfSpace(-n, -val))


    HpolyRep = HPolyhedron(constraints)
    box = overapproximate(HpolyRep, Hyperrectangle)
    Zrep = convert(Zonotope, box)
    return Zrep
end


function getBoxIntersection(Z::Zonotope, H_intersection::LazySets.HalfSpace)
    #println("Ever used? ")
    S = Z ∩ H_intersection
    if !isempty(S)
        box = box_approximation(S)
        return convert(Zonotope, box)
    else
        return S
    end
end

function getBoxIntersection(Z::Zonotope, H_intersections::Vector{N}) where N
    S1 = foldr((x, y) -> ∩(x, y), H_intersections; init=Z)
    S = S1 #overapproximate(S1, Zonotope) #foldr(∩, H_intersections; init=Z)

    if !isempty(S)
        S = overapproximate(S, HPolytope, dirs=BoxDirections())
        if isbounded(S)
            box = box_approximation(S)
            return convert(Zonotope, box)
        end
        println("Warning not bounded")
        throw(ErrorException("Warning not bounded"))
    end
end

function getBoxIntersection(Z::Zonotope, H_intersections::Any)
    #S1 = foldr((x, y) -> ∩(x, y), H_intersections; init=Z)
    S = ∩(H_intersections, Z)
    if !isempty(S)
        box = box_approximation(S)#overapproximate(S, Hyperrectangle)
        return convert(Zonotope, box)
    else
        #println("EMPTY!!")
        return S
    end
end

function splitZonotope(Z::Zonotope, H_intersection::LazySets.HalfSpace)
    H_rest = LazySets.HalfSpace(-H_intersection.a, -H_intersection.b)

    # Get intersections
    Z_intersection = getBoxIntersection(Z, H_intersection)
    Z_rest = getBoxIntersection(Z, H_rest)

    return Z_intersection, Z_rest
end

function splitZonotope(Z::Zonotope, H_intersections::Vector{N}) where N
    H_rest = map(x -> LazySets.HalfSpace(-x.a, -x.b), H_intersections)

    # Get intersections
    Z_intersection = getBoxIntersection(Z, H_intersections)
    Z_rest = getBoxIntersection(Z, H_rest)

    return Z_intersection, Z_rest
end

function splitZonotope(Z::Zonotope, H_intersections::LazySets.HPolyhedronModule.HPolyhedron)
    # Get intersections
    Z_intersection = zonotopeStripIntersection(Z, H_intersections)

    println(Z_intersection)

    # Flip it!
    constraint_list = H_intersections.constraints
    flipped_constraint_list = map(x -> LazySets.HalfSpace(-x.a, -x.b), constraint_list)

    println("Flipped constraint list: ", flipped_constraint_list)

    Z_rest = zonotopeStripIntersection(Z, HPolyhedron(flipped_constraint_list))

    return Z_intersection, Z_rest
end

#   Based on Alamo et al. "Guaranteed state estimation by zonotopes" (2005)
function zonotopeStripIntersection(Z::Zonotope, H::LazySets.HyperplaneModule.Hyperplane, σ::Float64)
    G = genmat(Z)
    c = Z.center
    b = H.a
    d = H.b
    λ = (G * transpose(G) * b) / (transpose(b) * G * transpose(G) * b + σ^2)
    ĉ = c + λ * (d - transpose(b) * c)
    upper = ((I - λ * transpose(b)) * G)
    lower = (σ * λ)
    Ĝ::Matrix{eltype(G)} = hcat(upper, lower)
    return Zonotope(ĉ, Ĝ)
end

function zonotopeStripIntersection(Z::Zonotope, H::LazySets.HPolyhedronModule.HPolyhedron)
    HalfSpaces = copy(constraints_list(H))
    HSG = stack([x.a for x in HalfSpaces]; dims=1)
    res = copy(Z)
    if rank(HSG) < size(HSG, 1)
        collinear = []
        remidx = stack([false for x in HalfSpaces])
        for hs in HalfSpaces
            if any(i -> rank([i.a hs.a]) <= 1 && (i !== hs), HalfSpaces) #any(i -> abs(dot(i.a, hs.a)) == norm(i.a) * norm(hs.a) && (i !== hs), HalfSpaces)
                push!(collinear, hs)
                remidx = remidx .|| any(i -> rank([i.a hs.a]) <= 1 && (i !== hs), HalfSpaces, dims=2)
            end
        end

        HalfSpaces = deleteat!(HalfSpaces, remidx)




        for hs in HalfSpaces
            a = hs.a
            b = hs.b
            if applicable(ρ, a, H)
                y = ρ(a, H)
                x = max(ρ(a, H), ρ(a, res))
                thp = LazySets.HyperplaneModule.Hyperplane(-a, (b + x) / 2)  #   Should check the calculation of the sigma values
                res = zonotopeStripIntersection(res, thp, abs(x - b) / 2)
            else
                x = ρ(a, Z)
                thp = HyperPlane(a, (2 * b + x) / 2)
                res = zonotopeStripIntersection(res, thp, x / 2)
            end
        end
        while !isempty(collinear)
            temphs = []
            push!(temphs, pop!(collinear))
            cols = any(i -> rank([i.a temphs[1].a]) <= 1, collinear, dims=2)
            for i in eachindex(cols)
                if cols[i] == true
                    push!(temphs, collinear[i])
                end
            end

            collinear = deleteat!(collinear, cols)

            minDists = map(x -> x.b / norm(x.a), temphs)
            diff = (minimum(minDists) - maximum(minDists)) / 2
            σ = diff == 0.0 ? eps(1.0) : diff

            a = temphs[1].a
            a = a ./ norm(a)
            if applicable(ρ, a, H)
                thp = LazySets.HyperplaneModule.Hyperplane(-a, diff)  #   Should check the calculation of the sigma values
                res = zonotopeStripIntersection(res, thp, σ)
            else

                x = ρ(a, Z)
                thp = HyperPlane(a, (2 * b + x) / 2)
                res = zonotopeStripIntersection(res, thp, x / 2)
            end

        end
    else
        for hs in HalfSpaces
            a = hs.a
            b = hs.b
            if applicable(ρ, a, H)
                #println("Applicable")
                y = ρ(a, H)
                x = max(ρ(a, H), ρ(a, res))
                #println(x, b)
                thp = LazySets.HyperplaneModule.Hyperplane(-a, (b + x) / 2)  #   Should check the calculation of the sigma values
                println((x - b) / 2)
                res = zonotopeStripIntersection(res, thp, (x - b) / 2)
            else
                println("Not applicable")
                x = ρ(a, Z)
                thp = HyperPlane(a, (2 * b + x) / 2)
                res = zonotopeStripIntersection(res, thp, x / 2)
            end
        end
    end
    return res
end