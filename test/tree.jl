"""Some tests how to compute MST in julia"""
using FLSA
using Graphs

n1, n2 = 4, 2
g = grid_graph(n1, n2)
w = rand(size(g.edges))

mst, wmst = kruskal_minimum_spantree(g, w)
