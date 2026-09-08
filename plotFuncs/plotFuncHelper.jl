using Plots, LazySets, LinearAlgebra, BenchmarkTools, CSV, DataFrames, Expokit, CDDLib #, ReachabilityAnalysis
using LaTeXStrings, Plots.PlotMeasures

include("../Utilities.jl")
include("../ReACTed.jl")
include("../models/gearbox.jl")
include("../models/platoon.jl")
include("../models/bouncingBall.jl")
include("../models/spacecraft.jl")


function getShapesForPlot(flowpipe, dims, alpha=1)
    
    amountOfDims = length(dims)

    shapesList = Vector()
    maxX = -Inf
    maxY = -Inf
    minX = Inf
    minY = Inf

    
    if amountOfDims == 1

        dim2 = dims[1]
        
        i = 1
        k = 0

        for (x, y) in flowpipe
            sen = true

            #println(y)
            for (d, t) in x
                #@show d

                push!(shapesList, Shape([t[1], t[2], t[2], t[1]], [d[1], d[1], -d[2], -d[2]]))

                maxY = max(maxY, d[1])
                minY = min(minY, -d[2])
                #Plots.plot!(Shape([mincor1, maxcor1, maxcor1, mincor1], [mincor2, mincor2, maxcor2, maxcor2]), c=cpallete[i], lab="")

                #plot!(r, c=cpallete[i], alpha=0.2)
            end
            #plot!(Shape([t[1], t[2], t[2], t[1]], [mincor, mincor, maxcor, maxcor]), c=cpallete[i], lab="", alpha=0.8)
            i += 1
        end
    elseif amountOfDims == 2
        dim1 = dims[1]
        dim2 = dims[2]

        i = 1
        k = 0
        for (x, y) in flowpipe
            #println(y)
            for (r, t) in x
                #=
                G = genmat(r)
                c = r.center

                projectedG = projectionMatrix * G
                projectGDim1s = mapreduce(x -> sign(x[dim1]) * x, +, eachcol(projectedG))
                projectGDim2s = mapreduce(x -> sign(x[dim2]) * x, +, eachcol(projectedG))

                maxcor1s = c + projectGDim1s
                mincor1s = c - projectGDim1s
                maxcor2s = c + projectGDim2s
                mincor2s = c - projectGDim2s

                projectGDim1 = reduce(+, reduce(+, G, dims=dim1))
                projectGDim2 = reduce(+, reduce(+, G, dims=dim2))
                maxcor1 = c[dim1] + projectGDim1
                mincor1 = c[dim1] - projectGDim1
                maxcor2 = c[dim2] + projectGDim2
                mincor2 = c[dim2] - projectGDim2
                =#
                #Plots.plot!(Shape([t[1], t[2], t[2], t[1]], [mincor, mincor, maxcor, maxcor]), c=cpallete[i], lab="", alpha=0.1)
                #Plots.plot!(Shape([mincor1s[dim1], mincor2s[dim1], maxcor2s[dim1], maxcor1s[dim1]], [mincor1s[dim2], maxcor1s[dim2], maxcor2s[dim2], mincor2s[dim2]]), c=cpallete[i], lab="") # Shape([mincor1s[dim1], mincor2s[dim1], maxcor2s[dim1], maxcor1s[dim1]], [mincor1s[dim2], mincor2s[dim2], maxcor2s[dim2], maxcor1s[dim2]])
                d1 = [r[1], -r[2]]#[ρ(sparsevec([dim1], [-1.0], ndim), r), ρ(sparsevec([dim1], [1.0], ndim), r)]
                d2 = [r[3], -r[4]]#[ρ(sparsevec([dim2], [-1.0], ndim), r), ρ(sparsevec([dim2], [1.0], ndim), r)]

                maxX = max(maxX, d1[1])
                maxY = max(maxY, d2[1])
                minX = min(minX, d1[2])
                minY = min(minY, d2[2])


                push!(shapesList, Shape([d1[1], d1[2], d1[2], d1[1]], [d2[1], d2[1], d2[2], d2[2]]))

            end
            #Plots.plot!(c=cpallete[i], lab=string(i))
            i += 1
        end
    else
        println("Cannot plot with amount of dims $amountOfDims")
        return shapesList, 0, 0, 0, 0
    end

    return shapesList, maxX, minX, maxY, minY
end

