using LazySets, ReachabilityAnalysis, LinearAlgebra, Polyhedra, Optim, JuMP, HiGHS

export HybridSystem, HybridSystemV2, Location, Edge, overapproximateIntervalReachset, intersects, splitZonotope, getBoxIntersection, getHalfSpaceProjections

Umodel = JuMP.Model(HiGHS.Optimizer)
set_attribute(Umodel, HiGHS.ComputeInfeasibilityCertificate(), false)
#Umodel = JuMP.direct_model(HiGHS.Optimizer())
set_string_names_on_creation(Umodel, false)
set_attribute(Umodel, "presolve", "off")
#set_attribute(model, "eps_abs", 1e-5)
#set_attribute(model, "eps_rel", 1e-5)
set_silent(Umodel)

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
    projVectors = map(x -> vec(x.a), halfspaces)
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

function isDisjointFast(Z::Zonotope, H)
    #res = all(map(x -> -ρ(-x.a, Z) > x.b, constraints_list(H)))
    #println(map(x -> -ρ(-x.a, Z) > x.b, constraints_list(H))[1])
    return any(map(x -> -ρ(-x.a, Z) > x.b, constraints_list(H)))
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
    if isSubSet(Z, H)
        return Z
    end

    G = genmat(Z)
    c = Z.center
    b = H.a
    d = H.b
    tb = transpose(b)
    # GG = G * transpose(G)
    # GGb = GG * b
    GGb = G * (transpose(G) * b)
    dotProduct = dot(tb, GGb) + σ^2
    # λ = GGB / (transpose(b) * GGb + σ^2)
    # λ = (G * transpose(G) * b) / (transpose(b) * G * transpose(G) * b + σ^2)
    λ = GGb / dotProduct
    ĉ = c + λ * (d - tb * c)
    upper = ((I - λ * tb) * G)
    lower = (σ * λ)
    Ĝ::Matrix{eltype(G)} = hcat(upper, lower)
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
    if isSubSet(Z, H)
        return Z
    end

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
                    found = true
                    remidx[index2] = true
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


        # @show collinear
        # @show typeof(collinear)

        for hs in HalfSpaces
            a = hs.a
            b = hs.b
            if applicable(ρ, a, H)
                y = ρ(a, H)
                x = max(ρ(a, H), ρ(-a, res))
                thp = LazySets.HyperplaneModule.Hyperplane(a, (b - x) / 2)  #   Should check the calculation of the sigma values
                σ = abs(x + b) / 2
                σ = σ == 0.0 ? eps(1.0) : σ
                res = zonotopeStripIntersection(res, thp, σ)
            else
                x = ρ(-a, Z)
                thp = HyperPlane(a, (2 * b + x) / 2)
                res = zonotopeStripIntersection(res, thp, x / 2)
            end
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

            minDists = map(x -> x.b / norm(x.a), temphs)
            minDists = map(x -> x.b, temphs)

            a = temphs[1].a
            diff = norm((minimum(minDists) * a + maximum(minDists) * a) / 2 - minimum(minDists) * a)

            σ = maximum(minDists) - diff
            σ = σ == 0.0 ? eps(1.0) : σ

            #a = a ./ norm(a)
            if applicable(ρ, a, H)
                thp = LazySets.HyperplaneModule.Hyperplane(a, diff)  #   Should check the calculation of the sigma values
                res = zonotopeStripIntersection(res, thp, σ)
            else

                x = ρ(a, Z)
                thp = HyperPlane(a, (2 * b + x) / 2)
                res = zonotopeStripIntersection(res, thp, x / 2)
            end

        end
    else
        #println("Linear independent")
        for hs in HalfSpaces
            a = hs.a
            b = hs.b
            if applicable(ρ, a, H)
                #println("Applicable")
                y = ρ(a, H)
                x = max(ρ(a, H), ρ(-a, res))
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
        end
    end
    return res
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
                    Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [mincor2, mincor2, maxcor2, maxcor2]), c=cpallete[i], linealpha=0)
                    sen = false
                else
                    Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [mincor2, mincor2, maxcor2, maxcor2]), c=cpallete[i], linealpha=0)
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

function plotProjectedFlowpipeLazy(flowpipe, dims, ndim, destination, alpha=1; xlim=nothing, ylim=nothing, verbose=false, legend=false)

    amountOfDims = length(dims)
    if amountOfDims == 1

        dim2 = amountOfDims[1]

        fig = Plots.plot(xlabel="time", ylabel="dim: " * string(dim2), ε=1e-6, legend=:outerright)
        cpallete = palette(:roma, length(flowpipe))
        i = 1
        k = 0
        for (x, y) in flowpipe
            sen = true
            if verbose
                println(y)

            end

            for (d, t) in x
                #@show d

                #d = [-ρ(sparsevec([dim2], [-1.0], ndim), r), ρ(sparsevec([dim2], [1.0], ndim), r)] #r[dim2]
                if sen
                    Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [d[1], d[1], -d[2], -d[2]]), c=cpallete[i], leg=legend, lab="S" * string(i), linealpha=0.1)
                    sen = false
                else
                    Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [d[1], d[1], -d[2], -d[2]]), c=cpallete[i], leg=legend, lab="", linealpha=0.1)
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
            if verbose
                println(y)

            end
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
                #@show r
                #Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [mincor, mincor, maxcor, maxcor]), c=cpallete[i], lab="", alpha=0.1)
                #Plots.plot!(Shape([mincor1s[dim1], mincor2s[dim1], maxcor2s[dim1], maxcor1s[dim1]], [mincor1s[dim2], maxcor1s[dim2], maxcor2s[dim2], mincor2s[dim2]]), c=cpallete[i], lab="") # Shape([mincor1s[dim1], mincor2s[dim1], maxcor2s[dim1], maxcor1s[dim1]], [mincor1s[dim2], mincor2s[dim2], maxcor2s[dim2], maxcor1s[dim2]])
                d1 = [r[1], -r[2]]#[ρ(sparsevec([dim1], [-1.0], ndim), r), ρ(sparsevec([dim1], [1.0], ndim), r)]
                d2 = [r[3], -r[4]]#[ρ(sparsevec([dim2], [-1.0], ndim), r), ρ(sparsevec([dim2], [1.0], ndim), r)]
                if sen
                    Plots.plot!(Shape([d1[1], d1[2], d1[2], d1[1]], [d2[1], d2[1], d2[2], d2[2]]), c=cpallete[i], leg=legend, lab="S" * string(i), linealpha=0)
                    sen = false
                else
                    Plots.plot!(Shape([d1[1], d1[2], d1[2], d1[1]], [d2[1], d2[1], d2[2], d2[2]]), c=cpallete[i], leg=legend, lab="", linealpha=0)
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
    if !isnothing(xlim)
        xlims!(fig, xlim)
    end
    if !isnothing(ylim)
        ylims!(fig, ylim)
    end
    default(fmt=:png)
    println("Saving..")
    savefig(fig, destination)
    println("Display..")
    display(fig)
end

function nestedInputDiscCalculate(inputDict, phiDict, δ⁺, δ⁻, currentTime, reduceOrder=5, maxOrder=5)
    # We know that δ⁻ % currentTime == 0
    n = LazySets.dim(inputDict[0])
    outputInput = Zonotope(zeros(n), zeros(n, 1))

    precomputedLargestStep = log2(δ⁺ / δ⁻)
    totalSteps = Int(round(currentTime / δ⁻))  # Steps we need to take. We round cause floats make small errors

    # Convert to bits 
    listToInclude = digits(totalSteps, base=2) # Get bit map


    largestInput = copy(inputDict[δ⁺])
    ϕ = phiDict[δ⁺]
    tempM = similar(ϕ)
    outputInput = Zonotope(zeros(size(tempM, 1)), zeros(size(tempM, 1), 1))

    precomputed = true
    i = 0
    for includeFlag in (listToInclude)
        if precomputed
            if includeFlag == 1 # If we have to add

                if isnothing(outputInput)
                    outputInput = copy(inputDict[2^i*δ⁻])

                else
                    outputInput = minkowski_sum(linear_map(phiDict[2^i*δ⁻], outputInput), inputDict[2^i*δ⁻])

                end
            end

            # Check if next step is also precomputed
            if !(i < precomputedLargestStep)
                precomputed = false
            end

        else
            # If not precomputed
            largestInput = minkowski_sum(largestInput, linear_map(ϕ, largestInput))
            ϕ = ϕ * ϕ
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
    return outputInput #reduce_order(outputInput, 5)
end

function projectReachSet(dirs, reachSet)
    projectedReachSet = [(map(x -> ρ(x, set), dirs), timings) for (set, timings) in reachSet]
    return projectedReachSet
end

function revise(approximation, lazyRepresentation, directions, bounds)
    newConstraints::Vector{LazySets.HalfSpace} = []
    #throw(ArgumentError("HECK NO WHY USED?"))
    for (idx, direction) in pairs(directions)
        n = norm(direction)
        nd = direction / n
        distance = ρ(nd, lazyRepresentation)
        if distance < bounds[idx] / n
            push!(newConstraints, LazySets.HalfSpace(nd, distance))
        end
    end
    newConstraints = vcat(approximation.constraints, newConstraints)
    return HPolytope(newConstraints)
end
function revise(approximation, lazyRepresentation, directions)
    newConstraints::Vector{LazySets.HalfSpace} = []
    for (idx, direction) in pairs(directions)
        n = norm(direction)
        nd = direction / n
        distance = ρ(nd, lazyRepresentation)
        tdistance = ρ(nd, approximation)#; solver=Umodel)
        #@show (distance, tdistance, (distance == tdistance))
        push!(newConstraints, LazySets.HalfSpace(nd, distance))

    end
    newConstraints = vcat(approximation.constraints, newConstraints)
    #@show newConstraints
    return HPolytope(newConstraints)
end

#=function revise(approximation, lazyRepresentation::Zonotope)
    newConstraints::Vector{LazySets.HalfSpace} = []
    directions = generators(lazyRepresentation)
    sumdir = reduce(+, directions)
    bounds = map(x -> ρ(sumdir - x, lazyRepresentation), directions)
    for (idx, direction) in pairs(directions)
        if ρ(direction, approximation) > bounds[idx]
            push!(newConstraints, LazySets.HalfSpace(direction, bounds[idx]))
        end
    end
    newConstraints = vcat(approximation.constraints, newConstraints)
    return HPolytope(newConstraints)
end=#

# TODO - FINISH FUNCTIONS
function revise(approximation::HPolytope, lazyRepresentation)
    newConstraints::Vector{LazySets.HalfSpace} = []
    directions = map(x -> x.a, constraints_list(approximation))
    #sumdir = reduce(+, directions)
    #bounds = map(x -> ρ(sumdir - x, lazyRepresentation), directions)
    for direction in directions
        n = norm(direction)
        nd = direction / n
        #ρ(direction, approximation)
        push!(newConstraints, LazySets.HalfSpace(nd, min(ρ(nd, approximation; solver=Umodel), ρ(nd, lazyRepresentation))))

    end
    return HPolytope(newConstraints)
end

function revise(input::LazySet, approximation, lazyRepresentation, directions, bounds)
    newConstraints::Vector{LazySets.HalfSpace} = []
    for (idx, direction) in pairs(directions)
        n = norm(direction)
        nd = direction / n
        Hdistance, Idistance = ρ(nd, lazyRepresentation; solver=Umodel), ρ(nd, input; solver=Umodel)
        if Hdistance + Idistance < bounds[idx] / n
            push!(newConstraints, LazySets.HalfSpace(nd, Hdistance))
        end
    end
    newConstraints = vcat(approximation.constraints, newConstraints)
    return HPolytope(newConstraints)
end

function revise(input::Vector{}, approximation, lazyRepresentation, directions, bounds)
    newConstraints::Vector{LazySets.HalfSpace} = []
    for (idx, direction) in pairs(directions)
        n = norm(direction)
        nd = direction / n
        distance = ρ(nd, lazyRepresentation; solver=Umodel)
        if distance + input[idx] < bounds[idx] / n
            push!(newConstraints, LazySets.HalfSpace(nd, distance))
        end
    end
    newConstraints = vcat(approximation.constraints, newConstraints)
    return HPolytope(newConstraints)
end

function constrain(approximation::HPolytope, lazyRepresentation, directions, bounds)
    newConstraints::Vector{LazySets.HalfSpace} = []
    for (idx, direction) in pairs(directions)
        n = norm(direction)
        nd = direction / n
        distance = min(ρ(nd, lazyRepresentation; solver=Umodel), bounds[idx] / n)
        #=if distance > bounds[idx]
            push!(newConstraints, LazySets.HalfSpace(direction, bounds[idx]))
        else
            println("$direction $distance")
            push!(newConstraints, LazySets.HalfSpace(direction, distance))
        end=#
        push!(newConstraints, LazySets.HalfSpace(nd, distance))

    end
    newConstraints = vcat(approximation.constraints, newConstraints)
    return HPolytope(newConstraints)
end

function constrain(input::LazySet, approximation, lazyRepresentation, directions, bounds)
    println("GETS USED???")
    newConstraints::Vector{LazySets.HalfSpace} = []
    #@show norm(approximation)
    for (idx, direction) in pairs(directions)
        n = norm(direction)
        ndir = direction / norm(direction)
        Hdistance, Idisctance = ρ(ndir, lazyRepresentation; solver=Umodel), ρ(ndir, input; solver=Umodel)
        distance = Hdistance + Idisctance
        #@show (Hdistance, Idisctance, bounds[idx] / n)
        #push!(newConstraints, LazySets.HalfSpace(direction, distance)) #bounds[idx]
        if distance > bounds[idx] / n  #(sign(bounds[idx]) == 1 ? distance > bounds[idx] : distance < bounds[idx])
            #println("$(bounds[idx])   $Hdistance $Idisctance")
            push!(newConstraints, LazySets.HalfSpace(ndir, bounds[idx] / n)) #(sign(bounds[idx]) == 1 ? bounds[idx] : distance)
        else
            push!(newConstraints, LazySets.HalfSpace(direction, distance))

        end
    end
    #println(newConstraints)
    #println(isempty(HPolytope(vcat(newConstraints, constraints(approximation)[end]))))
    #tres = HPolytope(newConstraints)
    constaintlist = map(normalize, constraints_list(approximation))
    constraintlist::Vector{LazySets.HalfSpace} = []
    dupes = Set()
    as = map(x -> x.a, constaintlist)
    if (length(as) != length(unique(as)))
        for (adx, ax) in pairs(as)
            if count(==(ax), as) > 1
                #@show (count(==(ax), as), ax)
                push!(dupes, ax)
            else
                push!(constraintlist, LazySets.HalfSpace(ax, constaintlist[adx].b))
            end
        end
        for dir in dupes
            push!(constraintlist, LazySets.HalfSpace(dir, ρ(dir, approximation; solver=Umodel)))
        end
    else
        constraintlist = constaintlist
    end
    #println(constaintlist)
    #println("Second loop")
    #@show norm(input)
    #@show LazySets.API.high(input)
    #@show LazySets.API.low(input)
    for constraint in constraintlist
        #println("Idist?")
        a = constraint.a
        n = norm(a)
        b = (constraint.b / n)
        na = a / n
        Idistance = ρ(na, input; solver=Umodel)
        #println("Idist!: $(Idistance)   $(constraint.a)")

        tconstraint = LazySets.HalfSpace(na, b + Idistance) #Idisctance
        if !isempty(HPolytope(vcat(newConstraints, [tconstraint])))
            push!(newConstraints, tconstraint)
            #println("$constraint    $Idistance $(constraint.b)")

        else
            println("NONONONONONOO $na $b $Idistance $(ρ(na, approximation)) $(ρ(na, lazyRepresentation))")
            push!(newConstraints, tconstraint)

            if isempty(approximation)

                #println(true)
            end
            #println("where $(ρ(constraint.a, approximation)) and $(ρ(-1. * constraint.a, approximation))") # $(ρ(constraint.a, tres)) and 
        end
    end
    #newConstraints = vcat(constraints(approximation), newConstraints)

    res = HPolytope(newConstraints)
    #res = remove_redundant_constraints!(res)
    #println(isempty(res))
    #println("Second loop!")

    return res
end

# Stupid but correct way
function constrain(input::LazySet, approximation, lazyRepresentation, invariantGuardIntersection)
    newConstraints::Vector{LazySets.HalfSpace} = []
    directions = map(x -> x.a, vcat(constraints_list(approximation), constraints_list(invariantGuardIntersection)))

    bounds0 = map(x -> ρ(x, input), directions) #, ρ(x, invariantGuardIntersection)
    #@show typeof(approximation)
    bounds1 = map(x -> ρ(x, approximation), directions) #, ρ(x, invariantGuardIntersection) ; solver=Umodel
    bounds2 = map(x -> ρ(x, lazyRepresentation), directions) #, ρ(x, invariantGuardIntersection)
    bounds = [min(b1, b2) for (b1, b2) in zip(bounds1, bounds2)]
    #@show maximum(bounds1)
    #@show maximum(bounds2)

    #@show maximum(bounds)

    #@show minimum(bounds1)
    #@show minimum(bounds2)
    #@show minimum(bounds)
    #@show bounds1
    #@show bounds2
    #@show bounds0
    hspaces::Vector{LazySets.HalfSpace} = [LazySets.HalfSpace(x, y + z) for (x, y, z) in zip(directions, bounds, bounds0)]

    res = HPolytope(hspaces)
    #@show norm(res)
    #@show LazySets.API.high(res)
    #@show LazySets.API.low(res)

    return res #LazySets.intersection(res, invariantGuardIntersection)
end


#=function constrain(input::Vector{}, approximation, lazyRepresentation, directions, bounds)
    newConstraints::Vector{LazySets.HalfSpace} = []
    for (idx, direction) in pairs(directions)
        Hdistance = ρ(direction, lazyRepresentation)
        distance = Hdistance + input[idx]
        push!(newConstraints, LazySets.HalfSpace(direction, bounds[idx]))
        if distance > bounds[idx]
            #println("$(bounds[idx])   $Hdistance $Idisctance")
        else
            #push!(newConstraints, LazySets.HalfSpace(direction, distance))

        end
    end
    #=
    as = map(x -> x.a, constraints_list(approximation))
    for (idx, direction) in pairs(as)
        Hdistance = ρ(direction, lazyRepresentation)
        distance = Hdistance + input[idx]
        push!(newConstraints, LazySets.HalfSpace(direction, bounds[idx]))
        if distance > bounds[idx]
            #println("$(bounds[idx])   $Hdistance $Idisctance")
        else
            #push!(newConstraints, LazySets.HalfSpace(direction, distance))

        end
    end
    =#
    #newConstraints = vcat(approximation.constraints, newConstraints)
    return HPolytope(newConstraints)
end=#

#=function constrain(input::Vector{}, approximation, directions)
    newConstraints::Vector{LazySets.HalfSpace} = []
    for (idx, direction) in pairs(directions)
        n = norm(direction)
        nd = direction / n
        Hdistance = ρ(nd, approximation)
        distance = Hdistance + input[idx]
        push!(newConstraints, LazySets.HalfSpace(direction, distance))

    end
    return HPolytope(newConstraints)
end=#


# TODO - Manage constraints that become inconsistent. ALSO FIX NORMALIZE
function bloatPolytope(input::Singleton, M, P::HPolytope)
    newP = nothing
    tempP = nothing
    newConstraints::Vector{LazySets.HalfSpace} = []

    try
        tempP = linear_map(M, P)
        #=
        hspaces = map(normalize, constraints_list(tempP))
        for hspace in hspaces
            Ma = hspace.a
            nM = norm(Ma)
            #@show nM
            b = hspace.b / nM

            push!(newConstraints, LazySets.HalfSpace(Ma / nM, b + ρ(Ma / nM, input)))
        end
        =#
        #newP = HPolytope(newConstraints)
        #@show LazySets.API.low(newP)
        newP = LazySets.translate(tempP, element(input))
    catch
        #println("Really doe?")

        #@show isinvertible(M)
        if isinvertible(M)
            inverseTransposeM = inv(transpose(M))
            hspaces = map(normalize, constraints_list(P))

            for hspace in hspaces
                Ma = inverseTransposeM * hspace.a
                nM = norm(Ma)
                b = hspace.b / nM

                push!(newConstraints, LazySets.HalfSpace(Ma / nM, b + ρ(Ma / nM, input)))
            end

            newP = HPolytope(newConstraints)
        else
            #=inputOffset = element(input)
            intPoint = LazySets.API.an_element(P)
            vs, ps = vertexRep(P)
            Mvs = map(x -> inputOffset + M * x, vs)
            Mps = map(x -> map(y -> inputOffset + M * y, x), ps)
            mIntPoint = inputOffset + M * intPoint
            listHspaces::Vector{LazySets.HalfSpace} = []
            for (x, y) in zip(Mvs, Mps)
                listHspaces = vcat(listHspaces, halfspaceFromVertices(x, y, mIntPoint))
            end
            =#
            newP = LazySets.translate(mpPol(M, P), element(input))
            #newP = bloatPolytope(input, diagm(ones(LazySets.API.dim(P))), mpPol(M, P))
            #=@show (LazySets.isempty(P), M, input)
            @show LazySets.API.high(P)
            @show LazySets.API.low(P)

            #= hspaces = map(normalize, constraints_list(P))
            push!(hspaces, LazySets.HalfSpace(SingleEntryVector(10, 10, 1.), 0.))
            push!(hspaces, LazySets.HalfSpace(SingleEntryVector(10, 10, -1.), 0.))


            newP = HPolytope(hspaces)=#

            throw(ErrorException("Tries to calculate linear map through vertices"))
            inputV = element(input)

            vertices = vertices_list(P)
            Mv = map(x -> (M * x), vertices)
            #iMv = element(input) .+ Mv
            println("Before convex_hull!..")
            #tempP = VPolytope(convex_hull!(Mv))
            tempP = VPolytope(Mv)
            tempP = minkowski_sum(tempP, input)
            println("Before tohrep...")
            newP = tohrep(tempP)
            println("HPolytope done")=#


        end
        #=oldConstraints = constraints_list(P)
        tempConstraints::Vector{LazySets.HalfSpace} = []
        ax = map(x -> x.a, oldConstraints)
        for (idx, a) in pairs(ax)
            if M * a != zero(a)
                push!(tempConstraints, LazySets.HalfSpace(a, oldConstraints[idx].b))
            else
                #push!(tempConstraints, LazySets.HalfSpace(a, 0.0))

            end
        end

        tempP = linear_map(M, HPolytope(tempConstraints))
        =#
    end
    #=
    if isinvertible(M)
        println("HECK YEA")
        inverseTransposeM = inv(transpose(M))
        hspaces = map(normalize, constraints_list(P))

        for hspace in hspaces
            Ma = inverseTransposeM * hspace.a
            b = hspace.b

            push!(newConstraints, LazySets.HalfSpace(Ma, b + ρ(Ma, input)))
        end

        newP = HPolytope(newConstraints)
    else
        println("HECK NO..")
        if applicable(x -> linear_map(M, x), P)
            tempP = nothing
            try
                tempP = linear_map(M, P)
            catch
                println("Really doe?")
                oldConstraints = constraints_list(P)
                tempConstraints::Vector{LazySets.HalfSpace} = []
                ax = map(x -> x.a, oldConstraints)
                for (idx, a) in pairs(ax)
                    if M * a != zero(a)
                        push!(tempConstraints, LazySets.HalfSpace(a, oldConstraints[idx].b))
                    else
                        #push!(tempConstraints, LazySets.HalfSpace(a, 0.0))

                    end
                end

                tempP = linear_map(M, HPolytope(tempConstraints))
            end
            hspaces = constraints_list(tempP)
            for hspace in hspaces
                a = hspace.a
                na = a / norm(a)
                b = hspace.b / norm(a)

                push!(newConstraints, LazySets.HalfSpace(na, b + ρ(na, input)))
            end
            newP = HPolytope(newConstraints)
        else
            inputV = element(input)
            vertices = vertices_list(P)
            Mv = map(x -> (M * x), vertices)
            #iMv = element(input) .+ Mv
            println("Before convex_hull!..")
            #tempP = VPolytope(convex_hull!(Mv))
            tempP = VPolytope(Mv)
            tempP = minkowski_sum(tempP, input)
            println("Before tohrep...")
            newP = tohrep(tempP)
            println("HPolytope done")
        end
    end
    #@show isempty(newP)
    =#
    return newP
end

function mpPol(M::Matrix, P::HPolytope)
    #throw(error("This should never be called! mpPol? more like.. "))
    #println("CRAZY CRAZY CRAZY")
    #@show LazySets.isempty(P), LazySets.isbounded(P)
    #@show LazySets.API.low(P), LazySets.API.high(P)

    intPoint = LazySets.API.an_element(P)

    vs, ps = vertexRep(P)
    #println("post vertexRep")
    Mvs = map(x -> M * x, vs)
    Mps = map(x -> map(y -> M * y, x), ps)
    mIntPoint = M * intPoint
    zidx = findall(==(zeros(length(vs[1]))), Mvs) # CHECK IF CORRECT
    listHspaces::Vector{LazySets.HalfSpace} = []
    #println("Before if")

    #if !isempty(zidx)
    H = LazySets.API.high(P)
    L = LazySets.API.low(P)
    C = [(x + y) / 2 for (x, y) in zip(H, L)]
    G = diagm(H - C)
    Z = Zonotope(C, G)
    #@show LazySets.API.high(Z)
    #@show LazySets.API.low(Z)
    MB = linear_map(M, Z)
    B = overapproximate(MB, BoxDirections(LazySets.dim(P)))
    #@show LazySets.isempty(B)
    #return MB
    for zdx in zidx
        #@show vs[zdx]
    end
    for xyidx in 1:length(vs)
        if !(xyidx in zidx)
            listHspaces = vcat(listHspaces, halfspaceFromVertices(Mvs[xyidx], Mps[xyidx], mIntPoint))
        else
            listHspaces = vcat(listHspaces, halfspaceFromVertices(Mvs[xyidx], Mps[xyidx], mIntPoint))

            #listHspaces = vcat(listHspaces, halfspaceFromVertices(vs[xyidx], ps[xyidx], intPoint))
            #@show halfspaceFromVertices(vs[xyidx], ps[xyidx], intPoint)
        end
    end
    #@show LazySets.API.high(MB)
    #@show LazySets.API.low(MB)
    #println("Gets pretty far")
    if isempty(listHspaces)
        tempP = B
        return tempP
    else
        tempP = HPolytope(listHspaces)
        tempP = LazySets.intersection(tempP, B; prune=false)
        return tempP
    end
    #@show ρ([0.0, 1.0, 0.0, 0.0, 0.0], tempP)
    #@show ρ([0.0, -1.0, 0.0, 0.0, 0.0], tempP)
    #for x in constraints_list(tempP)
    #println((vec(x.a) / norm(x.a)) * (x.b / norm(x.a)))
    #@show x.a
    #end

    #@show LazySets.isempty(tempP), LazySets.isbounded(tempP), LazySets.isbounded(B)
    #@show LazySets.API.high(tempP)
    #@show LazySets.API.low(tempP)
    #=
    else
    println("else") # Fails here!
    for (x, y) in zip(Mvs, Mps)
        listHspaces = vcat(listHspaces, halfspaceFromVertices(x, y, mIntPoint))
    end
    #@show listHspaces
    #println("before tempP")
    tempP = HPolytope(listHspaces)

    #@show ρ([0.0, 0.0, 0.0, 0.0, 1.0], tempP)
    return tempP
    end=#
end

function mapPolytope(M::Matrix, P::HPolytope; invertible=false)
    try
        tempP = linear_map(M, P)
        return tempP
    catch
        #println("linear_map failed")
        if invertible || isinvertible(M)
            inverseTransposeM = LinearAlgebra.inv(transpose(M))
            hspaces = constraints_list(P)
            newConstraints::Vector{LazySets.HalfSpace} = []
            #println("Made the inverse")
            for hspace in hspaces
                a = hspace.a
                b = hspace.b

                push!(newConstraints, LazySets.HalfSpace(inverseTransposeM * a, b))
            end
            #@show newConstraints
            tempP = HPolytope(map(normalize, newConstraints))
            return tempP

        else
            #@show M
            return mpPol(M, P)
            #=
            println("CRAZY CRAZY CRAZY")
            @show LazySets.isempty(P), LazySets.isbounded(P)
            @show LazySets.API.low(P), LazySets.API.high(P)
            intPoint = LazySets.API.an_element(P)
            println("post an_element")
            vs, ps = vertexRep(P)
            println("post vertexRep")
            Mvs = map(x -> M * x, vs)
            Mps = map(x -> map(y -> M * y, x), ps)
            mIntPoint = M * intPoint
            zidx = findall(==(zeros(length(vs[1]))), Mvs) # CHECK IF CORRECT
            listHspaces::Vector{LazySets.HalfSpace} = []
            if !isempty(zidx)
                H = LazySets.API.high(P)
                L = LazySets.API.low(P)
                C = [(x + y) / 2 for (x, y) in zip(H, L)]
                G = diagm(H - C)
                Z = Zonotope(C, G)
                #@show LazySets.API.high(Z)
                #@show LazySets.API.low(Z)
                MB = linear_map(M, Z)
                B = overapproximate(MB, BoxDirections(LazySets.dim(P)))
                #@show LazySets.isempty(B)
                #return MB
                for zdx in zidx
                    #@show vs[zdx]
                end
                for xyidx in 1:length(vs)
                    if !(xyidx in zidx)
                        listHspaces = vcat(listHspaces, halfspaceFromVertices(Mvs[xyidx], Mps[xyidx], mIntPoint))
                    else
                        #listHspaces = vcat(listHspaces, halfspaceFromVertices(vs[xyidx], ps[xyidx], intPoint))
                        #@show halfspaceFromVertices(vs[xyidx], ps[xyidx], intPoint)
                    end
                end
                #@show LazySets.API.high(B)
                #@show LazySets.API.low(B)
                tempP = HPolytope(listHspaces)
                #@show ρ([0.0, 1.0, 0.0, 0.0, 0.0], tempP)
                #@show ρ([0.0, -1.0, 0.0, 0.0, 0.0], tempP)
                tempP = LazySets.intersection(tempP, MB)

                #@show LazySets.API.high(tempP)
                #@show LazySets.API.low(tempP)
                @show LazySets.isempty(tempP)
                return tempP
            else
                println("Here")
                for (x, y) in zip(Mvs, Mps)
                    listHspaces = vcat(listHspaces, halfspaceFromVertices(x, y, mIntPoint))
                end
                tempP = HPolytope(listHspaces)
                @show ρ([0.0, 0.0, 0.0, 0.0, 1.0], tempP)
                return tempP
            end
            =#
        end

    end
end

function overapproximatedCH(H1::HPolytope, H2::HPolytope)
    newConstraints::Vector{LazySets.HalfSpace} = []
    directions = vcat(constraints_list(H1), constraints_list(H2))
    directions = map(x -> x.a / norm(x.a), directions)
    #directions = vcat(directions, collect(OctDirections(LazySets.dim(H1))))
    for direction in directions
        distance = max(ρ(direction, H1; solver=Umodel), ρ(direction, H2; solver=Umodel))
        push!(newConstraints, LazySets.HalfSpace(direction, distance))
    end
    #remove_redundant_constraints!(newConstraints)
    return HPolytope(newConstraints)
end

function overapproximatedCH(H1::HPolytope, H2::LazySet)
    newConstraints::Vector{LazySets.HalfSpace} = []
    directions = vcat(map(x -> x.a / norm(x.a), constraints_list(H1)))#, collect(OctDirections(LazySets.dim(H1))))
    for direction in directions
        #@show direction
        d1 = ρ(direction, H1) #; solver=Umodel
        d2 = ρ(direction, H2)
        distance = max(d1, d2)
        #@show distance
        push!(newConstraints, LazySets.HalfSpace(direction, distance))
    end
    #remove_redundant_constraints!(newConstraints)

    #println("New CH")
    return HPolytope(newConstraints)
end

function overapproximatedCH(H1::HPolyhedron, H2::LazySet)
    println("CRAZY CRAZY CRAZY")
    newConstraints::Vector{LazySets.HalfSpace} = []
    directions = vcat(map(x -> x.a, constraints_list(H1)), collect(BoxDirections(LazySets.dim(H1))))
    for direction in directions
        distance = max(ρ(direction, H1; solver=Umodel), ρ(direction, H2))
        push!(newConstraints, LazySets.HalfSpace(direction, distance))
    end
    #remove_redundant_constraints!(newConstraints)
    return HPolytope(newConstraints)
end

function vertexRep(P::HPolytope)
    hspaces = constraints_list(P)
    as = map(x -> vec(x.a), hspaces) #/ (norm(x.a)^2)
    bs = map(x -> (x.b), hspaces) #/ (norm(x.a)^2)
    #@show zip(as, bs)
    res = map(x -> (vec(x.a) / norm(x.a)) * (x.b / norm(x.a)), hspaces) #/ (norm(x.a)^2)
    #@show length(as[1])
    pes = []
    if length(as[1]) == 2
        T = [0.0 -1.0;
            1.0 0.0]
        for a in res
            push!(pes, [a + T * a / norm(T * a)])

        end
        #apes = map(y -> y + T * y / norm(T * y), res)
    else
        n = length(as[1])
        template = zeros(n)
        for a in res

            c = count(!=(0), a)
            tempv = copy(template)
            apes = []

            #F = svd(reshape(a / norm(a), 1, :); full=true)
            #V = F.V
            #@show V, a, reshape(a, 1, :)
            #push!(apes, V[:, 2:end])
            ns = nullspace(permutedims(a))
            #push!(apes, ns)
            apes = map(x -> x + a, collect(eachcol(ns)))

            #=
            if c > n - 2
                i = findfirst(!=(0), a)
                j = findfirst(!=(0), a[i:end])
                tempv[i] = a[j]
                tempv[j] = -a[i]
                push!(apes, tempv)
            else
                i = findfirst(==(0), a)
                tempv[i] = 1.0
                push!(apes, tempv)
            end



            m = 1

            while m < n - 2
                #shit
                m += 1
            end
            =#
            push!(pes, apes)

        end

    end

    #@show res, pes


    return res, pes
end

function halfspaceFromVertices(r, p, interiorPoint)
    #@show r, p
    listHspaces::Vector{LazySets.HalfSpace} = []
    n = length(r)
    k = length(p)
    #differenceVectors = map(x -> x - r, p) #stack(differenceVectors)
    D = zeros(eltype(r), k, n)
    for i in 1:k
        D[i, :] = p[i] - r
    end
    #@show D
    N = nullspace(D)
    if size(N, 2) == 1
        #@show N, D, r
        a = N[:, 1]
        a = a / norm(a)
        b = dot(a, r)

        if dot(a, interiorPoint) > b #
            a = -a
            b = -b
        end
        return [LazySets.HalfSpace(a, b)]
    else
        #@show N, size(N, 2)
        return []
        a = N[:, 1]
        a = a / norm(a)
        b = dot(a, r)
        if dot(a, interiorPoint) > b #
            a = -a
            b = -b
        end
        return [] #LazySets.HalfSpace(a, b), LazySets.HalfSpace(-a, b)

        #@show N, a, interiorPoint
    end
end

# TODO Check whether works in both directions
function reducePolytope(P::HPolytope; tolerance=0.1)
    constraints = constraints_list(P)
    listHspacesToKeep::Vector{LazySets.HalfSpace} = []
    for constraint in constraints
        if abs(ρ(constraint.a, P; solver=Umodel) - constraint.b) < tolerance
            push!(listHspacesToKeep, constraint)
        end
    end
    res = HPolytope(listHspacesToKeep)
    return res
end

function reducePolytopeFromBounding(P::HPolytope; tolerance=0.5)
    constraints = constraints_list(P)
    tconstraint::Vector{LazySets.HalfSpace} = map(x -> LazySets.HalfSpace(vec(x), ρ(x, P)), collect(BoxDirections(LazySets.dim(P))))
    res = HPolytope(tconstraint) #overapproximate(P, BoxDirections(LazySets.dim(P)))
    listHspacesToKeep::Vector{LazySets.HalfSpace} = []
    for constraint in constraints
        dir = constraint.a / norm(constraint.a)
        Pdist = ρ(dir, P; solver=Umodel)
        if abs(Pdist - ρ(dir, res)) > tolerance
            LazySets.addconstraint!(res, LazySets.HalfSpace(dir, Pdist))
        else
            #@show (ρ(constraint.a, P) - ρ(constraint.a, res))
        end
    end

    #=for (idx, dir) in pairs(mindists)
        if ρ(dir, P) == 1.0
            push!(listHspacesToKeep, constraints[idx])
        end
    end=#
    #res = HPolytope(listHspacesToKeep)
    #@show LazySets.API.high(res)
    #@show LazySets.API.low(res)
    #@show length(constraints_list(P))
    #@show length(constraints_list(res))

    return res
end

isinvertible(x::Matrix) = is_nonsingular(x) && applicable(LinearAlgebra.inv, x)
is_nonsingular(A) = !issuccess(lu(A, check=false)) ? false : true


#Base.:+(z1::Zonotope, z2::Zonotope) = Zonotope(z1.center + z2.center, z1.generators + z2.generators)