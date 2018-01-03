if VERSION >= v"0.7-"
    using Test
else
    using Base.Test
end
using FLSA
using Graph

tests = ["grid_graph",
         "incidence_matrix",
         "conjugate_gradient",
         "pwl",
         "line",
         "dynamic",
         "tree"]

for t in tests
    tp = Pkg.dir("FLSA", "test", "$(t).jl")
    println("running $(tp) ...")
    Base.include(tp)
end
