"""Some tests how to compute MST in julia"""
using FLSA
using Graphs

n1, n2 = 4, 2
g = grid_graph(n1, n2)
v = vertices(g)
n = length(v)
w = rand(size(g.edges))
y = rand(n1*n2)

mst, wmst = kruskal_minimum_spantree(g, w)
# Problem: edge indices are relative to g, not gmst.
# Because Edge is immutable: add all edges by myself
gmst = graph(v, Edge{(Int64,Int64)}[]; is_directed=false)
for e in mst
    add_edge!(gmst, source(e), target(e))
end

type DPVisitor <: AbstractGraphVisitor
    y::Vector{Float64}
    lb::Vector{Float64}
    ub::Vector{Float64}
    df::Vector{PWL}
    function DPVisitor(y::Vector{Float64})
        n = length(y)
        lb = fill(-Inf, n)
        ub = fill(+Inf, n)
        df = fill(PWL([]), n)
        println("Constructing visitor")
        new(y, lb, ub, df)
    end
end

examine_neighbor!(vis::DPVisitor, u, v, color, ecolor) = println("ex n")
function discover_vertex!(vis::DPVisitor, v)
    println("discover")
    true
end
close_vertex!(vis::DPVisitor, v) = asdfprintln("close vertex: $v")

visitor = DPVisitor(y)
alg = DepthFirst()
println("Starting ...", typeof(visitor))
traverse_graph(gmst, alg, (1,1), visitor)
