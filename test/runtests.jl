using Base.Test
using FLSA
using Graphs

tests = ["grid_graph", "incidence_matrix"]

for t in tests
    tp = joinpath(Pkg.dir("FLSA"),"test","$(t).jl")
    println("running $(tp) ...")
    include(tp)
end
