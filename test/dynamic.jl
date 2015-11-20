module TestDynamic

debug(msg) = println("DEBUG: $msg")

using Graphs
using FactCheck
using FLSA
import FLSA.@debug

digits = 3

facts("3 line knots") do
    sol = [2.0; 2.0; 2.5]
    y = [1.0; 2.0; 3.5]
    t = FLSA.create_tree([2,3,3])
    @debug "*"^70
    xl = FLSA.dp_line_naive(y, 1.0)
    @debug "*"^70
    x = FLSA.dp_tree_naive(y, 1.0, t)
    @fact x --> roughly(sol)
    @debug "*"^70
    x = FLSA.dp_tree(y, 1.0, t)
    @fact x --> roughly(sol)
end


facts("A random example") do
    srand(13)
    lambda = 1.0
    n1, n2 = 1, 3
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
