using HiGHS

mutable struct LPWorkspace
    highs::Ptr{Cvoid}
    m::Int
    n::Int
    x::Vector{Float64}
end

function LPWorkspace(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64})
    m, n = size(A)

    h = Highs_create()

    Highs_setBoolOptionValue(h, "output_flag", false)
    Highs_changeObjectiveSense(h, -1)

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