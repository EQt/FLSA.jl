module TestTreePWL
debug(msg) = println("DEBUG: $msg")
using Graphs
using FactCheck
using FLSA
reload("FLSA")

facts("A random 4x2 example") do
    begin
        srand(42)
        lambda = 1.0
        n1, n2 = 4, 2
        g = FLSA.igraph(FLSA.grid_graph(n1, n2))
        v = FLSA.vertices(g)
        n = length(v)
        w = rand(size(g.edges))
        y = round(10*rand(n1*n2), 1)
        root = 1

        mst, wmst = kruskal_minimum_spantree(g, w)
        t = FLSA.subtree(g, mst, root)
        tn = FLSA.PWLTree2(t, y)
        vis = FLSA.DPVisitor(y)
        x = FLSA.dp_tree_naive(y, lambda, t, vis)

        tf = FLSA.PWLTree(tn.parent, root, y, i->lambda)
        FLSA.forward_dp_treepwl(tf)
        for i=1:n
            @fact tf.nodes[i].lb --> roughly(round(vis.lb[i], 6)) "lb i = $i"
            @fact tf.nodes[i].ub --> roughly(round(vis.ub[i], 6)) "ub i = $i"
        end
    end
end

end
