using Pkg

Pkg.activate(".")
println("Current project: ", Base.active_project())

# Add standard libraries and external packages
packages = [
    "LinearAlgebra",
    "DelimitedFiles",
    "ForwardDiff",
    "Plots"
]

println("Adding packages: ", packages)
Pkg.add(packages)
Pkg.status()
