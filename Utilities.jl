using LazySets, LinearAlgebra

export HybridSystem

#=struct HybridSystem
    V::Vector{Any}
    E::Matrix{Int}
    G::Matrix{Vector{LazySet.HalfSpaceModule.HalfSpace}}
    I::Vector{LazySet.HalfSpaceModule.HalfSpace}
    Flow::Vector{Matrix{Float64}}
    Jump::Matrix{Tuple{Matrix{Float64},Vector{Float64}}}
end=#



struct Edge
    targetLoc :: Int
    guard :: HPolyhedron
    jumpMatrix :: Matrix{Float64}
    jumpVector :: Vector{Float64}
end

struct Location
    id :: Int
    invarient :: Union{HPolyhedron, Nothing}
    A :: Matrix{Float64}
    edges :: Vector{Edge}
end

struct HybridSystemV2
    locations :: Vector{Location}
    initialLoc :: Int
    initialState # Fill this in later
end


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

function project2d(Z :: Zonotope, n::Vector, direction::Vector)
    c = Z.center
    G = genmat(Z)

    newC = [dot(c, n), dot(c, direction)]
    newG =[[dot(g, n), dot(g, direction)] for g in eachcol(G)]

    return Zonotope(newC, newG)
end

function segmentLineIntersection(segment1, segment2, line)
    x1, y1 = segment1
    x2, y2 = segment2

    if (x1 - line)*(x2 - line) > 0 || (x1 == x2)
        return nothing
    end
    # Use interpolation
    t = (line - x1) / (x2 - x1)
    return y1 + t * (y2 - y1)
end

function GirardGuernicAlgorithm(Z :: Zonotope, line :: Float64)
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
function lineIntersection(Z :: Zonotope, H :: HalfSpace)
    directions = produceDirections(length(Z.center))
    val = H.b 
    n = H.a

    constraints = HalfSpace[]

    for direction in directions
        Z2 = project2d(Z, n, direction)
        m, M = GirardGuernicAlgorithm(Z2, val)

        #Produce the constraints
        if isfinite(m)
            push!(constraints, HalfSpace(direction, M))
            push!(constraints, HalfSpace(-direction, -m))
        end
    end

    # enforce hyperplane equality n⋅x = γ
    push!(constraints, HalfSpace(n, val))
    push!(constraints, HalfSpace(-n, -val))


    HpolyRep = HPolyhedron(constraints)
    box = overapproximate(HpolyRep, Hyperrectangle)
    Zrep = convert(Zonotope, box)
    return Zrep
end


function getBoxIntersection(Z :: Zonotope, H_intersection :: HalfSpace)
    S = Z ∩ H_intersection
    box = overapproximate(S, Hyperrectangle)
    return convert(Zonotope, box)
end


function splitZonotope(Z :: Zonotope, H_intersection :: HalfSpace)
    H_rest = HalfSpace(-H_intersection.a, -H_intersection.b)

    # Get intersections
    Z_intersection = getBoxIntersection(Z, H_intersection)
    Z_rest = getBoxIntersection(Z, H_rest)

    return Z_intersection, Z_rest
end
