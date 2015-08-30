"""Some tests how to compute MST in julia"""
using FLSA
using Graphs

srand(42)

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

lambda = 0.1
tm = FLSA.subtree(g, mst, (1,1))
x = FLSA.dp_tree(y, lambda, tm)
alpha = FLSA.dual_tree(y, x, tm)
abs(alpha) .<= lambda
