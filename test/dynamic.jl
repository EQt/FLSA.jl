module TestDynamic

debug(msg) = println("DEBUG: $msg")

using Graphs
using FactCheck
using FLSA

digits = 3

facts("A random example") do
    srand(13)
    lambda = 1.0
    n1, n2 = 5, 3
    g = FLSA.igraph(FLSA.grid_graph(n1, n2))
    v = FLSA.vertices(g)
    n = length(v)
    w = rand(size(g.edges))
    y = round(10*rand(n1*n2), 1)
    root = 1

    mst, wmst = kruskal_minimum_spantree(g, w)
    t = FLSA.subtree(g, mst, root)
    x = FLSA.dp_tree_naive(y, lambda, t)

    x2 = FLSA.dp_tree(y, lambda, t)
    @fact round(x, digits) --> roughly(round(x2, digits))
end

end # module
