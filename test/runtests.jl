tests = ["grid_graph"]


for t in tests
    tp = joinpath(Pkg.dir("FLSA"),"test","$(t).jl")
    println("running $(tp) ...")
    include(tp)
end
