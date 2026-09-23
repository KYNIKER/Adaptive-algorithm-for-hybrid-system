using HiGHS

mutable struct LPWorkspace
    highs::Ptr{Cvoid}
    m::Int
    n::Int
    x::Vector{Float64}
end

function LPWorkspace(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64}, objSense=kHighsObjSenseMaximize)
    m, n = size(A)

    h = Highs_create()

    Highs_setBoolOptionValue(h, "output_flag", false)
    Highs_changeObjectiveSense(h, objSense)

    for j in 1:n
        Highs_addCol(h, c[j], -Inf, Inf, 0, C_NULL, C_NULL)
    end

    indices = Vector{Cint}(undef, n)
    @inbounds for j in 1:n
        indices[j] = j - 1
    end

    @inbounds for i in 1:m
        Highs_addRow(
            h,
            -Inf,
            b[i],
            n,
            indices,
            @view A[i, :]
        )
    end

    Highs_run(h)

    return LPWorkspace(
        h,
        m,
        n,
        Vector{Float64}(undef, n),
    )
end

function update_rhs!(
    lp::LPWorkspace,
    b::Vector{Float64},
)
    h = lp.highs

    @inbounds for i in 1:lp.m
        Highs_changeRowBounds(h, i - 1, -Inf, b[i])
    end

    return nothing
end

function update!(
    lp::LPWorkspace,
    A::Matrix{Float64},
    b::Vector{Float64},
    c::Vector{Float64},
)
    @assert size(A) == (lp.m, lp.n)

    h = lp.highs

    # Update A
    @inbounds for i in 1:lp.m
        for j in 1:lp.n
            Highs_changeCoeff(h, i - 1, j - 1, A[i, j])
        end
    end

    # Update b
    @inbounds for i in 1:lp.m
        Highs_changeRowBounds(h, i - 1, -Inf, b[i])
    end

    # Update objective
    @inbounds for j in 1:lp.n
        Highs_changeColCost(h, j - 1, c[j])
    end

    return nothing
end

# https://github.com/JuliaReach/LazySets.jl/blob/54f65f2da50d70fe3c4f091418c5a2bea24c30d5/src/ConcreteOperations/isdisjoint.jl#L575
function disjointness_check(lpMinimizer::LPWorkspace, G::Matrix{Float64}, c::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})
    n = size(c)
    if n <= 2
        # this implementation is slower for low-dimensional sets
        #return _isdisjoint_polyhedron(Z, P, witness; solver=solver)
    end

    h = lpMinimizer.highs

    p = size(G, 2)
    m = length(d)

    @inbounds for i in 1:size(A, 1)
        for j in 1:size(A, 2)
            Highs_changeCoeff(h, i - 1, j - 1, A[i, j])
        end
        for n in size(A, 2):lpMinimizer.n
            Highs_changeCoeff(h, i - 1, n - 1, 0)
        end
    end

    @inbounds for i in 1:n
        HiGHS_changeCoeff(h, size(A, 1)+i-1, i - 1, 1)

        for j in 1:size(G, 1)
            HiGHS_changeCoeff(h, size(A, 1)+i-1, n + j - 1, G[i, j])
        end
    end

    @inbounds for i in 1:size(b)
        Highs_changeRowBounds(h, i - 1, -Inf, Inf)
    end

    @inbounds for i in 1:size(c)
        Highs_changeRowBounds(h, size(b) + i - 1, -1, 1)
    end

end

A = rand(Float64, (5, 5))
b = rand(Float64, 5)
c = rand(Float64, 5)


lp = LPWorkspace(A, b, c)

for k in 1:100

    # Modify A, b, c in-place here
    Ap = rand(Float64, (5, 5))
    bp = rand(Float64, 5)
    cp = rand(Float64, 5)
    update!(lp, Ap, bp, cp)

    Highs_run(lp.highs)

    Highs_getSolution(
        lp.highs,
        lp.x,
        C_NULL,
        C_NULL,
        C_NULL,
    )
    @show lp.x
    # lp.x now contains the solution
end

#=
@validate_commutative function isdisjoint(Z::AbstractZonotope, P::AbstractPolyhedron,
                                          witness::Bool=false; solver=nothing)
    n = dim(Z)
    if n <= 2
        # this implementation is slower for low-dimensional sets
        return _isdisjoint_polyhedron(Z, P, witness; solver=solver)
    end

    N = promote_type(eltype(Z), eltype(P))
    c = center(Z)
    G = genmat(Z)
    C, d = tosimplehrep(P)
    p = size(G, 2)
    m = length(d)

    A = [C zeros(N, m, p);
         I(n) -G]
    b = vcat(d, c)
    obj = zeros(N, size(A, 2))

    lbounds = vcat(fill(-Inf, n), fill(-one(N), p))
    ubounds = vcat(fill(Inf, n), fill(one(N), p))
    sense = vcat(fill('<', m), fill('=', n))
    if isnothing(solver)
        solver = default_lp_solver(N)
    end

    lp = linprog(obj, A, sense, b, lbounds, ubounds, solver)

    if is_lp_optimal(lp.status)
        disjoint = false
    elseif is_lp_infeasible(lp.status)
        disjoint = true
    else
        throw(ArgumentError("unexpected LP solver status: $(lp.status)"))
    end

    if disjoint
        return _witness_result_empty(witness, true, Z, P)
    elseif witness
        w = lp.sol[1:n]
        return (false, w)
    else
        return false
    end
end
=#