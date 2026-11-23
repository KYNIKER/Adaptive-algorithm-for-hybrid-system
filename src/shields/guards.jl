# ReACT's guard and constraint checks on zonotopes, and the set predicates (all_inside, some_inside, ...) built on them.

# ReACT's adaptive time stepping towards a polyhedron (dirs ⋅ x <= bounds): flows the discretization disc (per step
# length, as from newReACTDiscretizePlus) with steps from δ⁺ down to δ⁻, halving the step when a step reaches the
# polyhedron and doubling it again after four steps in a row that did not. step_Φ, step_Φt and step_input give, per step
# length, the matrix exponential, its transpose and the input's contribution. dirs is moved along with the flow (it is
# overwritten). Returns (time, set, input, flag): flag 2 with the time and set of the first δ⁻ step that reaches the
# polyhedron, or flag 0 at the end of the interval if none does.
function ReACT_guards(δ⁻::Float64, δ⁺::Float64, interval, identity_matrix, attempts_record, zero_input, dirs, bounds, step_Φ, step_Φt, disc, step_input)
    time::Float64 = copy(minimum(interval))
    end_time::Float64 = maximum(interval)
    step = copy(δ⁺)
    attempts_record = attempts_record .* 0
    orig_dirs = copy(dirs)
    Φ::Matrix{Float64} = copy(identity_matrix)
    tmp = copy(identity_matrix)
    i = 1
    # The input's contribution so far. With time-varying input (inputs with generators) the contribution of a new step is
    # independent of what came before, so the two are a Minkowski sum. It stays lazy: some_inside only needs its support
    # function, which is the sum of the parts' support functions. A constant input is composed exactly with the zonotope +.
    constant_input = iszero(genmat(step_input[δ⁻]))
    input = copy(zero_input)

    while time < end_time
        attempts = 1
        approved = false
        while !approved
            step >= δ⁻ || return (time, linear_map(Φ, disc[δ⁻]), input, 2) # reaches it even with the smallest step
            if misses_polyhedron(disc[step], input, dirs, orig_dirs, bounds)
                approved = true
                time + step < end_time || return (end_time, disc[δ⁻], input, 0)
                input = constant_input ? linear_map(step_Φ[step], input) + step_input[step] :
                        MinkowskiSum(LinearMap(step_Φ[step], input), step_input[step])
                map!(x -> step_Φt[step] * x, dirs)
                mul!(tmp, Φ, step_Φ[step])
                copy!(Φ, tmp)
            else
                step /= 2
                attempts += 1
            end
        end
        attempts_record[(i%4)+1] = attempts
        i = i + 1
        time = time + step
        # Double the step again after four steps in a row that succeeded at the first attempt
        if step < δ⁺ && all(==(1), attempts_record)
            step = step * 2
        end
    end
    error("ReACT_guards: the time stepping passed the end of the interval")
end

# True if the set plus the input misses the polyhedron, tested per halfspace: moved_dirs are the directions moved back
# along the flow so far, dirs and bounds the polyhedron's halfspaces.
misses_polyhedron(set::Zonotope, input, moved_dirs, dirs, bounds) = !some_inside(set, input, moved_dirs, dirs, bounds)

function all_inside(set, offset, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(offset, dir, bound)
        res = res && (offset[i] + ρ(dir[i], set) <= bound[i])
    end
    return res
end

function all_inside(set, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (ρ(dir[i], set) <= bound[i])
    end
    return res
end

function all_inside(set, input::LazySet, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (ρ(dir[i], set) + ρ(dir[i], input) <= bound[i])
    end
    return res
end

function some_inside(set, input::LazySet, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-dir[i], set)+(-ρ(-dir[i], input)) <= bound[i])
    end
    return res
end

function some_inside(set, input::LazySet, tdir, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-tdir[i], set)+(-ρ(-dir[i], input)) <= bound[i])
    end
    return res
end

function some_inside(set, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-dir[i], set) <= bound[i])
    end
    return res
end
