# ReACT's adaptive time stepping for a set of dirs with bounds. Utilizes combined support function and zonotope propagation.
function ReACT_guards(δ⁻::Float64, δ⁺::Float64, interval, identity_matrix, attempts_record, zero_input, constant_input::Bool, dirs, bounds, step_Φ, step_Φt, disc, step_input)
    time::Float64 = copy(minimum(interval))
    end_time::Float64 = maximum(interval)
    step = copy(δ⁺)
    attempts_record = attempts_record .* 0
    orig_dirs = copy(dirs)
    Φ::Matrix{Float64} = copy(identity_matrix)
    tmp = copy(identity_matrix)
    i = 1

    # A constant input is composed exactly with the zonotope +.
    input = copy(zero_input)

    while time < end_time
        attempts = 1
        approved = false
        while !approved
            step >= δ⁻ || return (time, linear_map(Φ, disc[δ⁻]), input, 2)
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
        if step < δ⁺ && all(==(1), attempts_record)
            step = step * 2
        end
    end
    error("ReACT_guards: the time stepping passed the end of the interval")
end

# True if the set plus the input misses the polyhedron, tested per halfspace.
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
