using LazySets, ReachabilityAnalysis, LinearAlgebra, Polyhedra, Optim

export HybridSystem, HybridSystemV2, Location, Edge, overapproximateIntervalReachset, intersects, splitZonotope, getBoxIntersection, getHalfSpaceProjections

struct Edge
    targetLoc::Int
    guard::Union{HPolyhedron,Nothing}
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
    constraints::Vector{LazySets.HalfSpace}
    # constraints::Vector{Union{HPolyhedron,LazySets.HalfSpace}}
end


Base.show(io::Core.IO, l::Location) = print(io, "Location: ", l.id, "\n invariant? ", !isnothing(l.invarient), "\n edges: ", length(l.edges))


mutable struct HybridSystemV2
    locations::Vector{Location}
    globalConstraints::Vector{LazySets.HalfSpace}
    #globalConstraints::Vector{Union{HPolyhedron,LazySets.HalfSpace}}
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

function getHalfSpaceProjections(H::HPolyhedron)
    getHalfSpaceProjections(H.constraints) # Convert to halfspace list
end

function getHalfSpaceProjections(halfspaces::Vector{<:LazySets.HalfSpace}) # Any subtype of halfspace
    #projVectors = map(x -> x.a, halfspaces)
    projVectors = map(x -> Vector(x.a), halfspaces)
    projBounds = ρ.(projVectors, halfspaces)
    return projVectors, projBounds
end

function getHalfSpaceProjections(H::Nothing)
    return ([], []) # return empty lists
end

function sparseHPolyhedronToDense(H_sparse::HPolyhedron)
    H_dense = HPolyhedron([LazySets.HalfSpace(Vector(c.a), c.b) for c in H_sparse.constraints])
    return H_dense
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
    return isnothing(H) ? true : !isdisjoint(Z, H)
end

function isSubSet(Z::Zonotope, H::LazySets.HalfSpace)
    agenSum = reduce(+, abs.(genmat(Z) .* H.a))
    acenSum = dot(Vector(H.a), Z.center)
    return (acenSum + agenSum <= H.b)
    #return (acenSum - agenSum <= H.b) & (acenSum + agenSum <= H.b)
end

function isSubSet(Z::Zonotope, H::HPolyhedron)
    sen = true
    for h in H.constraints
        agenSum = reduce(+, abs.(genmat(Z) .* h.a))
        acenSum = dot(Vector(h.a), Z.center)
        sen = sen & (ρ(h.a, Z) <= h.b) #(acenSum + agenSum <= h.b)
        #sen = sen & (acenSum - agenSum <= h.b) & (acenSum + agenSum <= h.b)
    end
    return sen
end

function isSubSet(Z::Zonotope, H::Any)
    return ⊆(Z, H, false) # Do not return a witness, just boolean
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
        println("empty?")
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
        box = overapproximate(S, Hyperrectangle)#box_approximation(S)#
        return convert(Zonotope, box)
    else
        println("EMPTY!!")
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
    #@show isSubSet(Z, H)
    #if isSubSet(Z, H)
    #    return Z
    #end

    G = genmat(Z)
    c = Z.center
    #b = H.a
    #d = H.b
    # GG = G * transpose(G)
    # GGb = GG * b
    tb = transpose(H.a)
    GGb = G * (transpose(G) * H.a)
    #dotProduct = dot(tb, GGb) + σ^2
    # λ = GGB / (transpose(b) * GGb + σ^2)
    # λ = (G * transpose(G) * b) / (transpose(b) * G * transpose(G) * b + σ^2)
    λ = GGb / (dot(tb, GGb) + σ^2)
    ĉ = c + λ * (H.b - tb * c)
    #upper = ((I - λ * tb) * G)
    #lower = (σ * λ)
    Ĝ::Matrix{eltype(G)} = hcat(((I - λ * tb) * G), (σ * λ))
    return Zonotope(ĉ, Ĝ)
end

function zonotopeStripIntersection(Z::Zonotope, a, b, σ::Float64)
    #@show isSubSet(Z, H)
    #if isSubSet(Z, H)
    #    return Z
    #end

    G = genmat(Z)
    c = Z.center
    #b = H.a
    #d = H.b
    # GG = G * transpose(G)
    # GGb = GG * b
    tb = transpose(a)
    GGb = G * (transpose(G) * a)
    #dotProduct = dot(tb, GGb) + σ^2
    # λ = GGB / (transpose(b) * GGb + σ^2)
    # λ = (G * transpose(G) * b) / (transpose(b) * G * transpose(G) * b + σ^2)
    λ = GGb / (dot(tb, GGb) + σ^2)
    ĉ = c + λ * (b - tb * c)
    #upper = ((I - λ * tb) * G)
    #lower = (σ * λ)
    Ĝ::Matrix{eltype(G)} = hcat(((I - λ * tb) * G), (σ * λ))
    return Zonotope(ĉ, Ĝ)
end

# Takes two vectors. May be sparrse
function iscollinear(a, b; atol=1e-10)
    # If the value is zero (Or very close to zero, then they are linearly dependent)
    # Squaring rather than doing norm because norm would take longer to calc 
    return abs(dot(a, b)^2 - dot(a, a) * dot(b, b)) ≤ atol
end

function zonotopeStripIntersection(Z::Zonotope, H::LazySets.HPolyhedronModule.HPolyhedron)
    #@show isSubSet(Z, H)
    if true #!isSubSet(Z, H)


        HalfSpaces = copy(constraints_list(H))
        HSG = stack([x.a for x in HalfSpaces]; dims=1)
        res = copy(Z)
        if rank(HSG) < size(HSG, 1)
            #println("Collinear")
            #collinear = []
            collinear = Vector{LazySets.HalfSpaceModule.HalfSpace{Float64,Vector{Float64}}}()
            #remidx = stack([false for x in HalfSpaces])
            remidx = falses(length(HalfSpaces))
            for hs in HalfSpaces
                flag = false
                for (index2, i) in enumerate(HalfSpaces)
                    if i !== hs && iscollinear(hs.a, i.a)
                        flag = true
                        @inbounds remidx[index2] = true
                    end
                end
                if flag
                    push!(collinear, hs)
                end
            end

            #         if any(i -> rank([i.a hs.a]) <= 1 && (i !== hs), HalfSpaces) #any(i -> abs(dot(i.a, hs.a)) == norm(i.a) * norm(hs.a) && (i !== hs), HalfSpaces)
            #         push!(collinear, hs)
            #         remidx = remidx .|| any(i -> rank([i.a hs.a]) <= 1 && (i !== hs), HalfSpaces, dims=2)
            #     end
            # end

            HalfSpaces = deleteat!(HalfSpaces, remidx)

            #@show remidx
            # @show typeof(collinear)

            for hs in HalfSpaces

                centerOfStrip = (hs.b - ρ(-hs.a, res)) / 2
                #@show centerOfStrip
                #@show ρ(-hs.a, res)
                #@show hs.b
                #σ = abs(x + hs.b) / 2
                σ = abs(centerOfStrip - hs.b)# / 2

                σ = σ == 0.0 ? eps(1.0) : σ
                #@show (hs.b - x) / 2
                #@show σ
                #res = zonotopeStripIntersection(res, thp, σ)
                #res = zonotopeStripIntersection(res, hs.a, (hs.b - x) / 2, σ)
                res = zonotopeStripIntersection(res, hs.a, centerOfStrip, σ)
                #a = hs.a
                #b = hs.b
                #centerOffset = dot(res.center, hs.a)



                #=
                x = max(hs.b, ρ(-hs.a, res))

                #thp = LazySets.HyperplaneModule.Hyperplane(hs.a, (hs.b - x) / 2)  #   Should check the calculation of the sigma values
                σ = abs(x + hs.b) / 2
                σ = σ == 0.0 ? eps(1.0) : σ
                #@show x
                #@show (hs.b - x) / 2
                #@show σ
                #@show centerOffset
                #res = zonotopeStripIntersection(res, thp, σ)
                res = zonotopeStripIntersection(res, hs.a, (hs.b - x) / 2, σ)
                =#

            end
            while !isempty(collinear)
                #temphs = []
                temphs = Vector{LazySets.HalfSpaceModule.HalfSpace{Float64,Vector{Float64}}}()
                push!(temphs, pop!(collinear))
                cols = any(i -> rank([i.a temphs[1].a]) <= 1, collinear, dims=2)
                for i in eachindex(cols)
                    if cols[i] == true
                        push!(temphs, collinear[i])
                    end
                end

                collinear = deleteat!(collinear, cols)
                #@show temphs
                #pos = filter(x -> dot(x.a, temphs[1].a) > 0, temphs)
                #@show pos
                #minDists = map(x -> x.b / norm(x.a), temphs)
                #centerOffset = dot(res.center, temphs[1].a)
                minDists = map(x -> x.b, temphs)
                #@show minDists

                minminDists = minimum(minDists)
                maxminDists = maximum(minDists)
                maxDist = min(temphs[1].b, ρ(temphs[1].a, res))
                minDist = min(temphs[2].b, ρ(-temphs[1].a, res))


                #a = temphs[1].a
                diff = norm((minminDists * temphs[1].a + maxminDists * temphs[1].a) / 2 - minminDists * temphs[1].a)
                diff = maxDist - (maxDist + minDist) / 2 #norm((minminDists * temphs[1].a + maxminDists * temphs[1].a) / 2 - minminDists * temphs[1].a)
                #σ = maxminDists - diff
                σ = maxDist - diff
                σ = σ == 0.0 ? eps(1.0) : σ
                #@show diff
                #@show σ
                #@show temphs[1].a
                #@show maxDist
                #@show minDist
                #@show centerOffset
                #a = a ./ norm(a)
                #thp = LazySets.HyperplaneModule.Hyperplane(temphs[1].a, diff)  #   Should check the calculation of the sigma values
                #res = zonotopeStripIntersection(res, thp, σ)
                res = zonotopeStripIntersection(res, temphs[1].a, diff, σ)

                #=
                if applicable(ρ, temphs[1].a, H)
                    thp = LazySets.HyperplaneModule.Hyperplane(a, diff)  #   Should check the calculation of the sigma values
                    res = zonotopeStripIntersection(res, thp, σ)
                else
                    println("Hopefully never")
                    x = ρ(temphs[1].a, Z)
                    thp = HyperPlane(temphs[1].a, (2 * b + x) / 2)
                    res = zonotopeStripIntersection(res, thp, x / 2)
                end
                =#
            end
        else
            #println("Linear independent")
            for hs in HalfSpaces
                #a = hs.a
                #b = hs.b

                #x = max(hs.b, ρ(-hs.a, res))
                #println(sign(b) * x, " ", b)
                #thp = LazySets.HyperplaneModule.Hyperplane(hs.a, (hs.b - x) / 2)  #   Should check the calculation of the sigma values
                #println((x - b) / 2)
                centerOfStrip = (hs.b - ρ(-hs.a, res)) / 2
                #@show centerOfStrip
                #@show ρ(-hs.a, res)
                #@show hs.b
                #σ = abs(x + hs.b) / 2
                σ = abs(centerOfStrip - hs.b)# / 2

                σ = σ == 0.0 ? eps(1.0) : σ
                #@show (hs.b - x) / 2
                #@show σ
                #res = zonotopeStripIntersection(res, thp, σ)
                #res = zonotopeStripIntersection(res, hs.a, (hs.b - x) / 2, σ)
                res = zonotopeStripIntersection(res, hs.a, centerOfStrip, σ)

                #=
                if applicable(ρ, a, H)
                    #println("Applicable")
                    #y = ρ(a, H)
                    #@show y == b
                    x = max(b, ρ(-a, res))
                    #println(sign(b) * x, " ", b)
                    thp = LazySets.HyperplaneModule.Hyperplane(a, (b - x) / 2)  #   Should check the calculation of the sigma values
                    #println((x - b) / 2)
                    σ = abs(x + b) / 2
                    σ = σ == 0.0 ? eps(1.0) : σ
                    res = zonotopeStripIntersection(res, thp, σ)
                else
                    #println("Not applicable")
                    x = ρ(a, Z)
                    thp = HyperPlane(a, (2 * b + x) / 2)
                    res = zonotopeStripIntersection(res, thp, x / 2)
                end
                =#
            end
        end
        return res
    else
        println("Ever")
        return Z
    end
end

function getUFromInputUncertainty(A, μ, δ⁻, P₁)
    ANorm = norm(A, Inf)
    β = (exp(ANorm * (δ⁻)) - 1) * μ / ANorm
    #println("original area ballβ: ", area(Zonotope(zeros(dim(P₁)), ((exp(ANorm*(initialTimeStep))-1)*μ/ANorm)*I(dim(P₁)))))
    u = Zonotope(zeros(LazySets.dim(P₁)), β * I(LazySets.dim(P₁)))
    return u
end




function plotProjectedFlowpipe(flowpipe, dim1, dim2, destination, alpha=1)
    fig = Plots.plot(xlabel="dim: " * string(dim1), ylabel="dim: " * string(dim2), ε=1e-6)
    cpallete = palette(:roma, length(flowpipe))
    i = 1
    k = 0

    if dim1 != 0
        dimSize = size(genmat(flowpipe[1][1][1][1]), 1)
        projectionMatrix = zeros(Float64, dimSize, dimSize)
        projectionMatrix[dim1, dim1] = 1.0
        projectionMatrix[dim2, dim2] = 1.0


        for (x, y) in flowpipe

            println(y)
            sen = true
            for (r, t) in x
                G = genmat(r)
                c = r.center

                # projectedG = projectionMatrix * G
                # projectGDim1s = mapreduce(x -> sign(x[dim1]) * x, +, eachcol(projectedG))
                # projectGDim2s = mapreduce(x -> sign(x[dim2]) * x, +, eachcol(projectedG))
                _, genAmount = size(G)
                projectGDim1s = sum(abs(G[dim1, i]) for i in 1:genAmount)
                projectGDim2s = sum(abs(G[dim2, i]) for i in 1:genAmount)

                maxcor1s = c + projectGDim1s
                mincor1s = c - projectGDim1s
                maxcor2s = c + projectGDim2s
                mincor2s = c - projectGDim2s

                projectGDim1 = reduce(+, reduce(+, G, dims=dim1))
                projectGDim2 = reduce(+, reduce(+, G, dims=dim2))
                maxcor1 = c[dim1] + projectGDim1
                mincor1 = c[dim1] - projectGDim1
                maxcor2 = c[dim2] + projectGDim2
                mincor2 = c[dim2] - projectGDim2

                #Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [mincor, mincor, maxcor, maxcor]), c=cpallete[i], lab="", alpha=0.1)
                #Plots.plot!(Shape([mincor1s[dim1], mincor2s[dim1], maxcor2s[dim1], maxcor1s[dim1]], [mincor1s[dim2], maxcor1s[dim2], maxcor2s[dim2], mincor2s[dim2]]), c=cpallete[i], lab="") # Shape([mincor1s[dim1], mincor2s[dim1], maxcor2s[dim1], maxcor1s[dim1]], [mincor1s[dim2], mincor2s[dim2], maxcor2s[dim2], maxcor1s[dim2]])
                if sen
                    Plots.plot!(Shape([(mincor1s[dim1], mincor1s[dim2]), (mincor2s[dim1], mincor2s[dim2]), (maxcor1s[dim1], maxcor1s[dim2]), (maxcor2s[dim1], maxcor2s[dim2])]), c=cpallete[i], lab="S" * string(i), linealpha=0) # Shape([mincor1s[dim1], mincor2s[dim1], maxcor2s[dim1], maxcor1s[dim1]], [mincor1s[dim2], mincor2s[dim2], maxcor2s[dim2], maxcor1s[dim2]])
                    sen = false
                else
                    Plots.plot!(Shape([(mincor1s[dim1], mincor1s[dim2]), (mincor2s[dim1], mincor2s[dim2]), (maxcor1s[dim1], maxcor1s[dim2]), (maxcor2s[dim1], maxcor2s[dim2])]), c=cpallete[i], lab="", linealpha=0) # Shape([mincor1s[dim1], mincor2s[dim1], maxcor2s[dim1], maxcor1s[dim1]], [mincor1s[dim2], mincor2s[dim2], maxcor2s[dim2], maxcor1s[dim2]])
                end
                #plot!(r, c=cpallete[i], alpha=0.2)
            end
            #Plots.plot!(c=cpallete[i], lab=string(i))
            i += 1

        end
    else
        for (x, y) in flowpipe
            sen = true

            println(y)
            for (r, t) in x
                G = abs.(genmat(r))
                c = r.center
                #projectGDim1 = reduce(+, reduce(+, G, dims=dim1))
                _, genAmount = size(G)
                projectGDim2 = sum(abs(G[dim2, i]) for i in 1:genAmount)
                #projectGDim2 = reduce(+, reduce(+, G, dims=dim2))
                #maxcor1 = c[dim1] + projectGDim1
                #mincor1 = c[dim1] - projectGDim1
                maxcor2 = c[dim2] + projectGDim2
                mincor2 = c[dim2] - projectGDim2
                if sen
                    Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [mincor2, mincor2, maxcor2, maxcor2]), c=cpallete[i], leg=false, linealpha=0)
                    sen = false
                else
                    Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [mincor2, mincor2, maxcor2, maxcor2]), c=cpallete[i], leg=false, linealpha=0)
                end
                #Plots.plot!(Shape([mincor1, maxcor1, maxcor1, mincor1], [mincor2, mincor2, maxcor2, maxcor2]), c=cpallete[i], lab="")

                #plot!(r, c=cpallete[i], alpha=0.2)
            end
            #plot!(Shape([t[1], t[2], t[2], t[1]], [mincor, mincor, maxcor, maxcor]), c=cpallete[i], lab="", alpha=0.8)
            i += 1

        end
    end
    savefig(fig, destination)
end

function plotProjectedFlowpipeLazy(flowpipe, dims, ndim, destination, alpha=1)

    amountOfDims = length(dims)
    if amountOfDims == 1

        dim2 = dims[1]

        fig = Plots.plot(xlabel="time", ylabel="dim: " * string(dim2), ε=1e-6)
        cpallete = palette(:roma, length(flowpipe))
        i = 1
        k = 0
        for (x, y) in flowpipe
            sen = true

            println(y)
            for (d, t) in x
                #@show d

                #d = [-ρ(sparsevec([dim2], [-1.0], ndim), r), ρ(sparsevec([dim2], [1.0], ndim), r)] #r[dim2]
                if sen
                    Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [d[1], d[1], -d[2], -d[2]]), c=cpallete[i], leg=false, linealpha=1)
                    sen = false
                else
                    Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [d[1], d[1], -d[2], -d[2]]), c=cpallete[i], leg=false, linealpha=1)
                end
                #Plots.plot!(Shape([mincor1, maxcor1, maxcor1, mincor1], [mincor2, mincor2, maxcor2, maxcor2]), c=cpallete[i], lab="")

                #plot!(r, c=cpallete[i], alpha=0.2)
            end
            #plot!(Shape([t[1], t[2], t[2], t[1]], [mincor, mincor, maxcor, maxcor]), c=cpallete[i], lab="", alpha=0.8)
            i += 1
        end
    elseif amountOfDims == 2
        dim1 = dims[1]
        dim2 = dims[2]
        fig = Plots.plot(xlabel="dim: " * string(dim1), ylabel="dim: " * string(dim2), ε=1e-6)
        cpallete = palette(:roma, length(flowpipe))
        i = 1
        k = 0
        for (x, y) in flowpipe
            println(y)
            sen = true
            for (r, t) in x
                #=
                G = genmat(r)
                c = r.center

                projectedG = projectionMatrix * G
                projectGDim1s = mapreduce(x -> sign(x[dim1]) * x, +, eachcol(projectedG))
                projectGDim2s = mapreduce(x -> sign(x[dim2]) * x, +, eachcol(projectedG))

                maxcor1s = c + projectGDim1s
                mincor1s = c - projectGDim1s
                maxcor2s = c + projectGDim2s
                mincor2s = c - projectGDim2s

                projectGDim1 = reduce(+, reduce(+, G, dims=dim1))
                projectGDim2 = reduce(+, reduce(+, G, dims=dim2))
                maxcor1 = c[dim1] + projectGDim1
                mincor1 = c[dim1] - projectGDim1
                maxcor2 = c[dim2] + projectGDim2
                mincor2 = c[dim2] - projectGDim2
                =#
                #Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [mincor, mincor, maxcor, maxcor]), c=cpallete[i], lab="", alpha=0.1)
                #Plots.plot!(Shape([mincor1s[dim1], mincor2s[dim1], maxcor2s[dim1], maxcor1s[dim1]], [mincor1s[dim2], maxcor1s[dim2], maxcor2s[dim2], mincor2s[dim2]]), c=cpallete[i], lab="") # Shape([mincor1s[dim1], mincor2s[dim1], maxcor2s[dim1], maxcor1s[dim1]], [mincor1s[dim2], mincor2s[dim2], maxcor2s[dim2], maxcor1s[dim2]])
                d1 = [r[1], -r[2]]#[ρ(sparsevec([dim1], [-1.0], ndim), r), ρ(sparsevec([dim1], [1.0], ndim), r)]
                d2 = [r[3], -r[4]]#[ρ(sparsevec([dim2], [-1.0], ndim), r), ρ(sparsevec([dim2], [1.0], ndim), r)]

                if sen
                    Plots.plot!(Shape([d1[1], d1[2], d1[2], d1[1]], [d2[1], d2[1], d2[2], d2[2]]), c=cpallete[i], leg=false, linealpha=0)

                    #Plots.plot!(Shape([(d1[1], d1[1]), (d2[1], d2[1]), (d1[2], d1[2]), (d2[2], d2[2])]), c=cpallete[i], lab="S" * string(i), linealpha=0) # Shape([mincor1s[dim1], mincor2s[dim1], maxcor2s[dim1], maxcor1s[dim1]], [mincor1s[dim2], mincor2s[dim2], maxcor2s[dim2], maxcor1s[dim2]])
                    sen = false
                else
                    Plots.plot!(Shape([d1[1], d1[2], d1[2], d1[1]], [d2[1], d2[1], d2[2], d2[2]]), c=cpallete[i], leg=false, linealpha=0)

                    #Plots.plot!(Shape([(d1[1], d1[1]), (d2[1], d2[1]), (d1[2], d1[2]), (d2[2], d2[2])]), c=cpallete[i], lab="", linealpha=0) # Shape([mincor1s[dim1], mincor2s[dim1], maxcor2s[dim1], maxcor1s[dim1]], [mincor1s[dim2], mincor2s[dim2], maxcor2s[dim2], maxcor1s[dim2]])
                end
                #plot!(r, c=cpallete[i], alpha=0.2)
            end
            #Plots.plot!(c=cpallete[i], lab=string(i))
            i += 1
        end
    else
        println("Cannot plot with amount of dims $amountOfDims")
        return
    end
    default(fmt=:png)
    savefig(fig, destination)
    display(fig)
end

function nestedInputDiscCalculate(inputDict, phiDict, δ⁻, currentTime, reduceOrder=5, maxOrder=5)
    # We know that δ⁻ % currentTime == 0

    precomputedLargestStep = (length(inputDict) - 2)#log2(δ⁺ / δ⁻)#
    largestTimeStepSize = δ⁻ * 2^precomputedLargestStep

    #@show (log2(δ⁺ / δ⁻), length(inputDict))
    totalSteps = Int(round(currentTime / δ⁻))  # Steps we need to take. We round cause floats make small errors


    # Convert to bits 
    listToInclude = digits(totalSteps, base=2) # Get bit map
    #@show listToInclude
    #println(string(totalSteps, base=2))
    largestInput = copy(inputDict[largestTimeStepSize])
    ϕ = phiDict[largestTimeStepSize]
    tempM = similar(ϕ)
    outputInput = Zonotope(zeros(size(tempM, 1)), zeros(size(tempM, 1), 1))
    sen = true
    precomputed = true
    i = 0 # Iterator for precomputed
    for includeFlag in (listToInclude)
        if precomputed
            if includeFlag == 1 # If we have to add
                if sen
                    #println("first")
                    outputInput = minkowski_sum(outputInput, inputDict[2^i*δ⁻])#copy(inputDict[2^i*δ⁻])#
                    sen = false
                else
                    #println("second")

                    outputInput = minkowski_sum(linear_map(phiDict[2^i*δ⁻], outputInput), inputDict[2^i*δ⁻])
                end
            end

            #println("true")
            # Check if next step is also precomputed
            if !(i < precomputedLargestStep)
                #println("false")

                precomputed = false
            end
        else
            largestInput = minkowski_sum(largestInput, LazySets.linear_map(ϕ, largestInput))
            ϕ = ϕ * ϕ
            inputDict[2^i*δ⁻] = copy(largestInput)
            phiDict[2^i*δ⁻] = copy(ϕ)
            if includeFlag == 1 # If we have to add
                if isnothing(outputInput)
                    outputInput = largestInput
                else
                    outputInput = minkowski_sum(linear_map(ϕ, outputInput), largestInput)
                end
            end
            if maxOrder > 0
                if LazySets.order(largestInput) > maxOrder
                    largestInput = reduce_order(largestInput, reduceOrder)
                end
            end
        end
        i += 1
    end

    # if maxOrder > 0
    #     if LazySets.order(outputInput) > maxOrder
    #         outputInput = reduce_order(outputInput, reduceOrder)
    #     end
    # end
    #@show LazySets.API.high(outputInput)
    return outputInput
end

function projectReachSet(dirs, reachSet)
    projectedReachSet = [(map(x -> ρ(x, set), dirs), timings) for (set, timings) in reachSet]
    return projectedReachSet
end













Base.:+(z1::Zonotope, z2::Zonotope) = Zonotope(z1.center + z2.center, z1.generators + z2.generators)