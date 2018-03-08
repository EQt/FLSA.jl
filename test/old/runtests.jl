if VERSION >= v"0.7-"
    using Test
    using FLSA
    using FLSA.Graph
    using LinearAlgebra     # for norm
else
    using Base.Test
    using FLSA
    using Graph
end

tests = ["grid_graph",
         "incidence_matrix",
         "conjugate_gradient",
         "line",
         "dynamic",
         "tree"]

for t in tests
    tp = Pkg.dir("FLSA", "test", "$(t).jl")
    include(tp)
end
