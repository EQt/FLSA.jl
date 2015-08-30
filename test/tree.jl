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

## visitor = DPVisitor(y)
## alg = DepthFirst()
## println("Starting ...", typeof(visitor))
## traverse_graph(gmst, alg, (1,1), visitor)

