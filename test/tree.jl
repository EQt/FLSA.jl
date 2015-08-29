"""Some tests how to compute MST in julia"""
using FLSA
using Graphs

n1, n2 = 4, 2
g = grid_graph(n1, n2)
w = rand(size(g.edges))

mst, wmst = kruskal_minimum_spantree(g, w)
gmst = graph(vertices(g), mst; is_directed=false)

visitor = LogGraphVisitor(STDOUT)
alg = DepthFirst()
traverse_graph(gmst, alg, (1,1), visitor)
