# -ρ(-d, Z) without negating d; bit-identical, as IEEE negation of products and sums is exact.
lower_ρ(d::AbstractVector, Z::Zonotope) = dot(LazySets.center(Z), d) - LazySets.abs_sum(d, genmat(Z))
lower_ρ(d::AbstractVector, X) = -ρ(-d, X)

# ReACT's adaptive time stepping over [0, end_time] for halfspaces (dirs, bounds); dirs is overwritten, moved back along the flow.
# Returns (time, 2) for the first δ⁻ step that reaches the polyhedron, or (end_time, 0).
function ReACT_guards(δ⁻::Float64, δ⁺::Float64, end_time::Float64, dirs, bounds, step_Φt, disc, step_input)
    time, step, i = 0.0, δ⁺, 1
    attempts_record = zeros(4)
    input_min = zeros(length(dirs)) # per halfspace, the least of its moved direction over the input's contribution so far
    tmp = similar(first(dirs))

    while time < end_time
        attempts = 1
        approved = false
        while !approved
            step >= δ⁻ || return (time, 2)
            if misses_polyhedron(disc[step], input_min, dirs, bounds)
                approved = true
                time + step < end_time || return (end_time, 0)
                P, Φt = step_input[step], step_Φt[step]
                for (j, d) in enumerate(dirs)
                    input_min[j] += lower_ρ(d, P)
                    copyto!(d, mul!(tmp, Φt, d))
                end
            else
                step /= 2
                attempts += 1
            end
        end
        attempts_record[(i%4)+1] = attempts
        i = i + 1
        time = time + step
        if step < δ⁺ && all(==(1), attempts_record)
            step = step * 2
        end
    end
    error("ReACT_guards: the time stepping passed the end of the interval")
end

# True if the set plus the input misses the polyhedron, tested per halfspace.
misses_polyhedron(set::Zonotope, input_min, moved_dirs, bounds) = any(i -> lower_ρ(moved_dirs[i], set) + input_min[i] > bounds[i], eachindex(bounds))

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
        res = res && (lower_ρ(dir[i], set) + lower_ρ(dir[i], input) <= bound[i])
    end
    return res
end

function some_inside(set, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (lower_ρ(dir[i], set) <= bound[i])
    end
    return res
end
