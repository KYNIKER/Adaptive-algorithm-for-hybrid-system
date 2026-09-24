using LazySets
using HiGHS
using LinearAlgebra
using SparseArrays
using Plots

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


"""
    zonotope_polyhedron_intersect_direct(Z::Zonotope, P::AbstractPolyhedron)

Checks if a Zonotope and a Polyhedron intersect by passing the constraint 
matrices directly to the HiGHS C-API wrapper, completely avoiding JuMP.
"""
function zonotope_polyhedron_intersect_direct(Z::Zonotope, P::AbstractPolyhedron, highs_model)

    # 2. Extract polyhedron constraints Ax <= b
    constraints = constraints_list(P)
    num_constraints = length(constraints)

    # If the polyhedron contains no constraints, they trivially intersect
    if num_constraints == 0
        return true
    end
    # 1. Extract zonotope representations
    c = Z.center
    G = genmat(Z)
    dim_space, num_gens = size(G)

    # Pack A and b from the half-spaces into raw matrices
    A_poly = zeros(num_constraints, dim_space)
    b_poly = zeros(num_constraints)
    for (i, constraint) in enumerate(constraints)
        A_poly[i, :] = constraint.a
        b_poly[i] = constraint.b
    end

    # 3. Transform the problem to be strictly in terms of variable ξ:
    # (A_poly * G) * ξ <= b_poly - A_poly * c
    A_new = A_poly * G
    b_new = b_poly - A_poly * c



    # 5. Define Column (Variable) Parameters: ξ ∈ [-1, 1]
    num_col = Cint(num_gens)
    col_cost = zeros(Cdouble, num_gens)          # No objective function
    col_lower = fill(-1.0, num_gens)             # ξ_i >= -1
    col_upper = fill(1.0, num_gens)              # ξ_i <= 1

    # 6. Define Row (Constraint) Parameters: -Inf <= A_new * ξ <= b_new
    num_row = Cint(num_constraints)
    row_lower = fill(-Inf, num_constraints)      # No lower bound on constraints
    row_upper = Vector{Cdouble}(b_new)           # Upper bound is b_poly - A_poly * c

    # 7. Convert the dense matrix A_new to Compressed Sparse Column (CSC) format for HiGHS
    A_sparse = sparse(A_new)

    # HiGHS expects 0-indexed pointers for its C-backend matrix reading
    a_format = Cint(1)                           # 1 indicates CSC format
    a_start = Vector{Cint}(A_sparse.colptr .- 1)
    a_index = Vector{Cint}(A_sparse.rowval .- 1)
    a_value = Vector{Cdouble}(A_sparse.nzval)

    # 8. Load the entire problem structure directly into the HiGHS solver
    HiGHS.Highs_passLp(
        highs_model, num_col, num_row, Cint(length(a_value)),
        a_format, Cint(1), 0.0, col_cost, col_lower, col_upper,
        row_lower, row_upper, a_start, a_index, a_value
    )

    #@inbounds for i in 1:dim_space
    #    Highs_changeRowBounds(highs_model, i - 1, -Inf, b[i])
    #end


    # 9. Solve the linear program
    HiGHS.Highs_run(highs_model)

    # 10. Query model status (HighsModelStatus enum)
    # 7 = Optimal, 9 = Feasible
    model_status = HiGHS.Highs_getModelStatus(highs_model)

    # Destroy the C object allocation safely from memory
    #HiGHS.Highs_destroy(highs_model)

    return model_status == 7 || model_status == 9
end

# ==========================================
# Test Example
# ==========================================

# Setup a 2D Zonotope
c_z = [0.0, 0.0]
G_z = [1.0 1.2;
    0.1 1.0]
Z = Zonotope(c_z, G_z)

# 4. Initialize a raw Highs instance
highs_model = HiGHS.Highs_create()

# Mute logging output for clean terminal execution
HiGHS.Highs_setBoolOptionValue(highs_model, "output_flag", false)

# A Polyhedron that overlaps
hs_intersect = HPolyhedron([LazySets.HalfSpace([1.0, 1.0], 1.0)])

# A Polyhedron that is far away
hs_miss = HPolyhedron([LazySets.HalfSpace([-1.0, -0.0], -1.6), LazySets.HalfSpace([0.0, 1.0], 0.5)])

#lp = LPWorkspace(A, b, c_z)

println("Intersecting? ", zonotope_polyhedron_intersect_direct(Z, hs_intersect, highs_model)) # Expected: true
println("Intersecting? ", zonotope_polyhedron_intersect_direct(Z, hs_miss, highs_model))      # Expected: false


plt = plot(dpi=1200, thickness_scaling=1, guidefontsize=35, minorgrid=true,
    #legendfont=font(12, "Times"),
    #legend_position=:topright,
    legend=false)

xlims!(plt, -2.0, 10.0)
ylims!(plt, -2.0, 10.0)

bbox = Hyperrectangle([4.0, 4.0], [6.0, 6.0])

plot!(plt, Z, c=:green)
plot!(plt, intersection(hs_intersect, bbox), c=:red)
plot!(plt, intersection(hs_miss, bbox), c=:blue)

display(plt)
#=
Z = rand(Zonotope; dim=2, num_generators=20)

@time begin
    for i in 1:100
        local Z = rand(Zonotope; dim=2, num_generators=20)

        _ = zonotope_polyhedron_intersect_direct(Z, hs_intersect, highs_model)
        _ = zonotope_polyhedron_intersect_direct(Z, hs_miss, highs_model)
    end
end

@time begin
    for i in 1:100
        local Z = rand(Zonotope; dim=2, num_generators=20)

        _ = LazySets.API.isdisjoint(Z, hs_intersect)
        _ = LazySets.API.isdisjoint(Z, hs_miss)
    end
end
HiGHS.Highs_destroy(highs_model)
=#