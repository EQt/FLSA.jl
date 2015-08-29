"""Some tests how to compute MST in julia"""
using FLSA
using Graphs

n1, n2 = 4, 2
g = grid_graph(n1, n2)
w = rand(size(g.edges))

mst, wmst = kruskal_minimum_spantree(g, w)
# Problem: edge indices are relative to g, not gmst.
# Because Edge is immutable: add all edges by myself
gmst = graph(vertices(g), Edge{(Int64,Int64)}[]; is_directed=false)
for e in mst
    add_edge!(gmst, source(e), target(e))
end

visitor = LogGraphVisitor(STDOUT)
alg = DepthFirst()
traverse_graph(gmst, alg, (1,1), visitor)
