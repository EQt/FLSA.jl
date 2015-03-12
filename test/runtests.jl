using FLSA
using Base.Test

tests = ["grid_graph", "incidence_matrix"]

for t in tests
    tp = joinpath(Pkg.dir("FLSA"),"test","$(t).jl")
    println("running $(tp) ...")
    include(tp)
end
