# What modal_reach! (modal.jl) needs to flow a set under one (A, U) and through jumps. Uses ReACT_guards (guards.jl).

dense_zonotope(Z::Zonotope) = Zonotope(Vector{Float64}(LazySets.center(Z)), Matrix{Float64}(genmat(Z)))

# Everything needed to flow a set under one (A, U). Φ[k+1] and V[k+1] solve the dynamics exactly for k*δ⁻.
struct Dynamics
    A::Matrix{Float64}
    # per step length of ReACT's time stepping (δ⁻ to δ⁺): the matrix exponential, its transpose and the input's contribution
    step_Φ::Dict{Float64,Matrix{Float64}}
    step_Φt::Dict{Float64,Matrix{Float64}}
    step_input::Dict{Float64,Zonotope}
    P2A_abs::Matrix{Float64}
    Φ::Vector{Matrix{Float64}}
    V::Vector{Zonotope{Float64,Vector{Float64},Matrix{Float64}}}
    zero_input::Zonotope{Float64,Vector{Float64},Matrix{Float64}}
end

function Dynamics(A::Matrix{Float64}, U, steps::Int, δ⁻::Float64, δ⁺::Float64, alg, maxOrder::Int, reduceOrder::Int)
    n = size(A, 1)
    step_Φ, step_Φt, step_input = phi_input_dict(A, U, δ⁻, δ⁺, alg, maxOrder, reduceOrder)
    A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
    Φcache = A == A_abs ? step_Φ[δ⁻] : nothing
    P2A_abs = Matrix{Float64}(ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isinvertible(A), Φcache))
    Φ = Matrix{Float64}[Matrix{Float64}(I, n, n)]
    V = [Zonotope(zeros(n), zeros(n, n))]
    for k in 1:steps
        d = k * δ⁻
        push!(Φ, exp(d .* A))
        push!(V, input_bound(A, U, d, alg)) # the input's contribution over k δ⁻, which may vary in time
    end
    return Dynamics(A, step_Φ, step_Φt, step_input, P2A_abs, Φ, V, Zonotope(zeros(n), zeros(size(genmat(step_input[δ⁻])))))
end

# Step at which the flowpipe of disc first touches any of the sets, or steps if it touches none within that many steps of δ⁻.
function first_touch(dyn::Dynamics, disc, sets, steps::Int, δ⁻::Float64, δ⁺::Float64, identity_matrix)
    k0 = steps
    for (vectors, bounds) in sets
        time, _, _, flag = ReACT_guards(δ⁻, δ⁺, [0.0, steps * δ⁻], identity_matrix, zeros(4), dyn.zero_input, copy(vectors), bounds, dyn.step_Φ, dyn.step_Φt, disc, dyn.step_input)
        if flag != 0
            k0 = min(k0, round(Int, time / δ⁻))
        end
    end
    return k0
end


# True if a guard has the halfspace -a ⋅ x <= -b, so it borders the state-space face a ⋅ x <= b.
bordered_by(guards, a, b; atol) = any(((v, w),) -> any(i -> isapprox(v[i], -a; atol) && isapprox(w[i], -b; atol), eachindex(v)), guards)

# Over-approximates the part of Z in the polyhedron with one strip intersection per halfspace that cuts it.
# Returns nothing if Z misses a halfspace.
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
