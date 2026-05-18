using Pkg
Pkg.activate(".")

# This creates a manifest, with the following packages

libraries = ["Plots", "LazySets", "LinearAlgebra", "BenchmarkTools", "CSV", "DataFrames", "Expokit", "SparseArrays", "ReachabilityAnalysis", "Polyhedra", "Optim", "CDDLib", "LaTeXStrings"]

Pkg.add(libraries)
