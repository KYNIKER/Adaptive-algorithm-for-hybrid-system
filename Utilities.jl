using LazySets

export HybridSystem

#=struct HybridSystem
    V::Vector{Any}
    E::Matrix{Int}
    G::Matrix{Vector{LazySet.HalfSpaceModule.HalfSpace}}
    I::Vector{LazySet.HalfSpaceModule.HalfSpace}
    Flow::Vector{Matrix{Float64}}
    Jump::Matrix{Tuple{Matrix{Float64},Vector{Float64}}}
end=#

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