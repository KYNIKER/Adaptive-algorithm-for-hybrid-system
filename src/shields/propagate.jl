dense_zonotope(Z::Zonotope) = Zonotope(Vector{Float64}(LazySets.center(Z)), Matrix{Float64}(genmat(Z)))

# Everything needed to flow a set under one (A, U). Φ[k+1] and V[k+1] solve the dynamics exactly for k*δ⁻.
struct Dynamics
    A::Matrix{Float64}
    # per step length of ReACT's time stepping (δ⁻ to δ⁺): the matrix exponential, its transpose and the input's contribution
    step_Φ::Dict{Float64,Matrix{Float64}}
    step_Φt::Dict{Float64,Matrix{Float64}}
    step_input::Dict{Float64,Zonotope{Float64,Vector{Float64},Matrix{Float64}}}
    P2A_abs::Matrix{Float64}
    Φ::Vector{Matrix{Float64}}
    V::Vector{Zonotope{Float64,Vector{Float64},Matrix{Float64}}}
end

function Dynamics(A::Matrix{Float64}, U, steps::Int, δ⁻::Float64, δ⁺::Float64, alg, maxOrder::Int, reduceOrder::Int)
    n = size(A, 1)
    step_Φ, step_Φt, step_input_any = phi_input_dict(A, U, δ⁻, δ⁺, alg, maxOrder, reduceOrder)
    step_input = Dict(k => dense_zonotope(Z) for (k, Z) in step_input_any)
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    P2A_abs = Matrix{Float64}(ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, false, nothing))
    Φ = Matrix{Float64}[Matrix{Float64}(I, n, n)]
    V = [Zonotope(zeros(n), zeros(n, n))]
    Φ₁ = Matrix{Float64}[zeros(n, n)] # Φ₁[k+1] is Φ₁(A, k δ⁻)
    M = Matrix{Float64}(A_abs)
    varying = !iszero(genmat(U))
    M_exp, M_Φ₁, M_Φ₂ = [copy(Φ[1])], [zeros(n, n)], [zeros(n, n)] # the same for M, only filled for a varying input
    for k in 1:steps
        d = k * δ⁻
        push!(Φ, exp(d .* A))
        # Φ₁(A, s + t) = Φ₁(A, s) + e^{As} Φ₁(A, t)
        a = k ÷ 2
        push!(Φ₁, k == 1 ? Matrix{Float64}(ReachabilityAnalysis.Exponentiation.Φ₁(A, d, alg, false, nothing)) : Φ₁[a+1] + Φ[a+1] * Φ₁[k-a+1])
        if varying
            # Φ₂(M, s + t) = Φ₂(M, s) + t Φ₁(M, s) + e^{Ms} Φ₂(M, t), and e^{Mt} = I + M Φ₁(M, t), where nothing cancels: M >= 0
            if k == 1
                push!(M_Φ₁, Matrix{Float64}(ReachabilityAnalysis.Exponentiation.Φ₁(M, d, alg, false, nothing)))
                push!(M_Φ₂, Matrix{Float64}(ReachabilityAnalysis.Exponentiation.Φ₂(M, d, alg, false, nothing)))
            else
                push!(M_Φ₁, M_Φ₁[a+1] + M_exp[a+1] * M_Φ₁[k-a+1])
                push!(M_Φ₂, M_Φ₂[a+1] + ((k - a) * δ⁻) .* M_Φ₁[a+1] + M_exp[a+1] * M_Φ₂[k-a+1])
            end
            push!(M_exp, I + M * M_Φ₁[k+1])
        end
        push!(V, input_bound(A, U, d, Φ₁[k+1], varying ? M_Φ₂[k+1] : M_Φ₂[1]))
    end
    return Dynamics(A, step_Φ, step_Φt, step_input, P2A_abs, Φ, V)
end

function first_touch(dyn::Dynamics, disc, sets, steps::Int, δ⁻::Float64, δ⁺::Float64)
    k0 = steps
    for (vectors, bounds) in sets
        time, flag = ReACT_guards(δ⁻, δ⁺, [0.0, steps * δ⁻], zeros(4), copy(vectors), bounds, dyn.step_Φt, disc, dyn.step_input)
        if flag != 0
            k0 = min(k0, round(Int, time / δ⁻))
        end
    end
    return k0
end


# True if a guard has the halfspace -a ⋅ x <= -b, so it borders the state-space face a ⋅ x <= b.
bordered_by(guards, a, b; atol) = any(((v, w),) -> any(i -> isapprox(v[i], -a; atol) && isapprox(w[i], -b; atol), eachindex(v)), guards)

# Over-approximates the part of Z in the polyhedron with one strip intersection per halfspace that cuts it.
function clip_to_polyhedron(Z, polyhedron)
    res = Z
    for (a, b) in zip(polyhedron...)
        if ρ(a, res) > b
            lo = -ρ(-a, res)
            lo > b && return nothing
            res = zonotopeStripIntersection(res, a, (lo + b) / 2, (b - lo) / 2)
        end
    end
    return res
end

function jump_set(Z, M::Matrix{Float64}, j::Vector{Float64}, maxOrder::Int, reduceOrder::Int)
    res = remove_zero_generators(Zonotope(M * LazySets.center(Z) + j, M * genmat(Z)))
    size(genmat(res), 2) == 0 && return Zonotope(LazySets.center(res), zeros(length(j), 1))
    if maxOrder > 0 && LazySets.order(res) > maxOrder
        res = reduce_order(res, reduceOrder)
    end
    return dense_zonotope(res)
end
