# Some tests how to compute MST in julia

module TestTreeFunctions

import FLSA
import Graph

function test1()
    srand(42)

    n1, n2 = 4, 2
    img = FLSA.img_graph(n1, n2)
    g = img.graph
    v = Graph.vertices(g)
    n = length(v)
    w = rand(size(g.edges))
    y = rand(n1*n2)

    mst, wmst = Graph.kruskal_minimum_spantree(g, w)
    lambda = 0.1
    tm = FLSA.subtree(g, mst, 1)
    x = FLSA.dp_tree(y, lambda, tm)
    alpha = FLSA.dual_tree(y, x, tm)
    abs(alpha) .<= lambda + 1e-9
end

end

TestTreeFunctions.test1()
