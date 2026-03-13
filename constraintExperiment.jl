using LazySets, BenchmarkTools, Random

const amountOfTests = 10
const amountOfConstraints = 1
BenchmarkTools.DEFAULT_PARAMETERS.samples = 1

totalReduce = 0.
totalAny = 0.

for _ in 1:amountOfTests
    Z = rand(Zonotope)
    constraints = [rand(LazySets.HalfSpace) for _ in 1:amountOfConstraints]
    #[1:amountOfConstraints] # Produce an amount of random constraints 

    constraintProjVectors = map(x -> x.a, constraints)
    constraintProjBounds = map(x -> x.b, constraints)

    GC.gc()
    t1 = @benchmark _ = !reduce(&, <=(map(x -> ρ(x, $Z), $constraintProjVectors), $constraintProjBounds))
    # result = run(b; ve)
    global totalReduce += mean(t1.times)

    GC.gc()
    t2 = @benchmark _ = any((ρ(x, $Z) > y) for (x, y) in zip($constraintProjVectors, $constraintProjBounds))

    global totalAny += mean(t2.times)
end

    
println("Accumulated time for $amountOfConstraints constraints with $amountOfTests tests: 
        \n reduce: $totalReduce
        \n any: $totalAny")