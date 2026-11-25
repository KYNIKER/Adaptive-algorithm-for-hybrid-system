# ReACT's adaptive time stepping for a set of dirs with bounds. The set and the input are propagated by support functions: dirs is moved
# back along the flow (it is overwritten). Returns (time, flag): flag 2 with the time of the first δ⁻ step that reaches the polyhedron, or flag 0 at the end.
function ReACT_guards(δ⁻::Float64, δ⁺::Float64, interval, attempts_record, dirs, bounds, step_Φt, disc, step_input)
    time::Float64 = copy(minimum(interval))
    end_time::Float64 = maximum(interval)
    step = copy(δ⁺)
    attempts_record = attempts_record .* 0
    i = 1

    # Per halfspace the minimum of its direction over the input's contribution W so far. With W(t + s) = W(t) + e^{At} W(s)
    # it grows by the minimum of the moved direction over the contribution of the new step.
    input_min = zeros(length(dirs))

    while time < end_time
        attempts = 1
        approved = false
        while !approved
            step >= δ⁻ || return (time, 2)
            if misses_polyhedron(disc[step], input_min, dirs, bounds)
                approved = true
                time + step < end_time || return (end_time, 0)
                P = step_input[step]
                for i in eachindex(dirs)
                    input_min[i] -= ρ(-dirs[i], P)
                end
                Φt = step_Φt[step]
                map!(x -> Φt * x, dirs)
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
misses_polyhedron(set::Zonotope, input_min, moved_dirs, bounds) = any(i -> -ρ(-moved_dirs[i], set) + input_min[i] > bounds[i], eachindex(bounds))

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

function some_inside(set, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-dir[i], set) <= bound[i])
    end
    return res
end
