using Base.Test
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
    include(tp)
end
