# ReACT's guard and constraint checks on zonotopes, and the set predicates (allunder, someinside, ...) built on them.

function propagate_set_check_hybrid_jump(X0, interval, δ⁻::Float64, δ⁺::Float64, system, listOfEdges, identity_matrix, attemptsRecorder, dims, Φ₂, guardProjVectors, guardProjBounds, invariantProjVectors, invariantProjBounds, discretizationDict, PhiDict, TPhiDict, inputDict, max_input, max_flow, guard_intersection, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5)
    zero_input = Zonotope(zeros(dims), zeros(size(genmat((inputDict[δ⁻])))))
    for edge in listOfEdges
        guard_projection_vectors = copy(guardProjVectors)
        guard_projection_bounds = copy(guardProjBounds)


        #   Compute the reachset closest to the guard without intersecting it and not reaching the unsafe set. 
        _, _, _, flag = ReACT_guards_t(δ⁻, δ⁺, interval, identity_matrix, attemptsRecorder, zero_input, guard_projection_vectors, guard_projection_bounds, 2, PhiDict, TPhiDict, discretizationDict, inputDict)


        #
        #   0: Didnt hit either constraint, guard or invariant
        #   1: Hit constraint
        #   2: Hit guard
        #   3: Hit invariant surface
        #
        if flag != 0
            return 1
        end
        if allunder(linear_map(max_flow, X0), max_input, invariantProjVectors, invariantProjBounds)
            continue
        else
            return 1
        end
    end
    return 0
end


function ReACT_guards_t(δ⁻::Float64, δ⁺::Float64, interval, identity_matrix, attemptsRecorder, zero_input, guardProjVectors, guardProjBounds, STRATEGY::Integer, PhiDict, permutedphiDict, discritezationDict, inputDiscretizationDict)
    time::Float64 = copy(minimum(interval))
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(δ⁺)

    attemptsRecorder = attemptsRecorder .* 0
    origGuardProjVectors = copy(guardProjVectors)
    Φ::Matrix{Float64} = copy(identity_matrix)

    tempM = copy(identity_matrix)

    i = 1
    # The input's contribution so far. With time-varying input (inputs with generators) the contribution of a new step is
    # independent of what came before, so the two are a Minkowski sum. It stays lazy: someinside only needs its support
    # function, which is the sum of the parts' support functions. A constant input is composed exactly with the zonotope +.
    constant_input = iszero(genmat(inputDiscretizationDict[δ⁻]))
    input = copy(zero_input)

    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            # Handle if we can no longer reduce the reachset (we keep hitting something)
            if currentTimeStep >= δ⁻

                if fast_guards_check_t(discritezationDict[currentTimeStep], input, guardProjVectors, origGuardProjVectors, guardProjBounds)
                    approveFlag = true
                    if currentTimeStep + time < endtime

                        input = constant_input ? linear_map(PhiDict[currentTimeStep], input) + inputDiscretizationDict[currentTimeStep] :
                                MinkowskiSum(LinearMap(PhiDict[currentTimeStep], input), inputDiscretizationDict[currentTimeStep])

                        map!(x -> permutedphiDict[currentTimeStep] * x, guardProjVectors)

                        mul!(tempM, Φ, PhiDict[currentTimeStep])

                        copy!(Φ, tempM)

                    else
                        return (endtime, discritezationDict[δ⁻], input, 0)
                    end
                else

                    currentTimeStep /= 2
                    attempts += 1
                end

            else
                return (time, linear_map(Φ, discritezationDict[δ⁻]), input, 2)
            end

        end

        attemptsRecorder[(i%4)+1] = attempts
        i = i + 1
        time = time + currentTimeStep

        # Reset / apply strategy
        # Only do this if the current timestep is less than the initial
        if STRATEGY == 0
            # Only reduce
        elseif STRATEGY == 1
            # always try double
            if currentTimeStep < δ⁺
                currentTimeStep = currentTimeStep * 2
                #changedTimeStep = true
            end
        elseif STRATEGY == 2
            # If attemptsrecorder past 4 are successes, double timestep
            if currentTimeStep < δ⁺
                #lowest = min(4, i - 1)
                #window = @view attemptsRecorder[i-lowest:i-1]
                if all(==(1), attemptsRecorder)
                    currentTimeStep = currentTimeStep * 2
                    #changedTimeStep = true
                end

            end
        end
    end
    throw(error("ERROR!!! we should not have reached this point. "))
end

# This can be replaced with throwing an error. Currently we continue and just print
function handleHitConstraint(time, locationId)
    throw(error("ERROR!!! We have hit a constraint at loc: $(locationId) time: $time"))
    println("\nERROR!!!\n 
            ERROR!!!\n\n
            We have hit a constraint at loc: $locationId time: $time\n\n
            ERROR!!!\n
            ERROR!!!\n")
end


function guardCheck(input::Zonotope, Φ, X::Zonotope, constraint, guard, invariant) #; solver=model
    #
    #   0: Didnt hit either constraint, guard or invariant
    #   1: Hit constraint
    #   2: Hit guard
    #   3: Hit invariant surface
    #

    newR = minkowski_sum(input, linear_map(Φ, X))
    if LazySets.API.isdisjoint(newR, constraint)
        if LazySets.API.isdisjoint(newR, guard)
            return 0
        else
            return 2
        end

    else
        return 1
    end
end

function fast_guards_check_t(newR::Zonotope, input, tguardProjVectors, guardProjVectors, guardProjBounds) #; solver=model
    #return all((input + ρ(x, newR)) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds)) && any((input + -ρ(-x, newR)) > y for (input, x, y) in zip(Gρ, guardProjVectors, guardProjBounds)) && all((input - ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))
    return !someinside(newR, input, tguardProjVectors, guardProjVectors, guardProjBounds)
    #return !someinside(minkowski_sum(newR, input), guardProjVectors, guardProjBounds)

end
function touchesCheck(newR::Zonotope, constraintProjVectors, constraintProjBounds, guardProjVectors, guardProjBounds, invarientProjVectors, invarientProjBounds) #; solver=model
    return allunder(newR, constraintProjVectors, constraintProjBounds) && someinside(newR, guardProjVectors, guardProjBounds) && someinside(newR, invarientProjVectors, invarientProjBounds)#all((input - ρ(-x, newR)) <= y for (input, x, y) in zip(Iρ, invarientProjVectors, invarientProjBounds))
end

function allunder(set, offset, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(offset, dir, bound)
        res = res && (offset[i] + ρ(dir[i], set) <= bound[i])
    end
    return res
end

function someunder(set, offset, dir, bound)
    res = !LazySets.API.isdisjoint(set, HPolyhedron(collect(LazySets.HalfSpace(a, b+c) for (a, b, c) in zip(dir, bound, offset))))
    return res
end

function someoutside(set, offset, dir, bound)
    return !allunder(set, offset, dir, bound)
end

function allunder(set, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (ρ(dir[i], set) <= bound[i])
    end
    return res
end

function allunder(set, input::LazySet, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (ρ(dir[i], set) + ρ(dir[i], input) <= bound[i])
    end
    return res
end

function someunder(set, dir, bound)
    return !LazySets.API.isdisjoint(set, HPolyhedron(collect(LazySets.HalfSpace(a, b) for (a, b) in zip(dir, bound))))
end

function someoutside(set, input::LazySet, dir, bound)
    res = false
    @inbounds @simd for i in eachindex(dir, bound)
        res = res || (-ρ(-dir[i], set)-ρ(-dir[i], input) > bound[i])
    end
    return res
end

function someinside(set, input::LazySet, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-dir[i], set)+(-ρ(-dir[i], input)) <= bound[i])
    end
    return res
end

function someinside(set, input::LazySet, tdir, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-tdir[i], set)+(-ρ(-dir[i], input)) <= bound[i])
    end
    return res
end

function someinside(set, dir, bound)
    res = true
    @inbounds @simd for i in eachindex(dir, bound)
        res = res && (-ρ(-dir[i], set) <= bound[i])
    end
    return res
end

