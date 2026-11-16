# Installs the packages of Project.toml. Run once after cloning: julia scripts/setup.jl
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()
