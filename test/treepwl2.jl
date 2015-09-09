module TestFLSA
using Graphs
using FLSA
using FactCheck


facts("2 node line") do
    y = [3.5, 2.5]
    root = 2
    parents = [root, root]
    t = FLSA.PWLTree(parents, root, y, i->1.0)
    e1 = FLSA.create_min_event(t, 1, -1.0)
    e2 = FLSA.create_max_event(t, 1, +1.0)

    @fact e1.x --> roughly(2.5)
    @fact e1.s --> 2
    @fact e1.t --> 1
    @fact e1.slope --> roughly(1)
    @fact e1.offset --> roughly(-2.5)

    @fact e2.x --> roughly(4.5)
    @fact e2.s --> 1
    @fact e2.t --> 2
    @fact e2.slope --> roughly(-1)
    @fact e2.offset --> roughly(+4.5)
    
    n = t.nodes[2]
    n.minevs = [e1]
    n.maxevs = [e2]
    @fact FLSA.create_min_event(t, 2, -1.0).x --> roughly(2.5)
    @fact FLSA.create_max_event(t, 2, +1.0).x --> roughly(3.5)
end

facts("A 4 node, 3 level tree") do
    begin
        lambda = 1.0
        y = [3.5, 4.0, 2.5, 3.0]
        E = [(1,3), (2,4), (3,4)]
        E = [Graphs.IEdge(i, u, v) for (i, (u,v)) in enumerate(E)]
        g = Graphs.edgelist(collect(1:4), E; is_directed=false)
        tn = FLSA.subtree(g, E, 4)
        vis = FLSA.DPVisitor(y)
        x = FLSA.dp_tree_naive(y, lambda, tn, vis)

        parents = [3, 4, 4, 4]
        root = 4
        context("clip functions") do
            tf = FLSA.PWLTree(parents, root, y, i->1.0)
            @fact FLSA.create_min_event(tf, 1, -1.0).x --> 2.5
            @fact FLSA.create_max_event(tf, 1, +1.0).x --> 4.5
        end
        context("forward tree") do
            tf = FLSA.PWLTree(parents, root, y, i->1.0)
            @fact length(tf.pre_order) --> 4
            for i=1:4
                @fact tf.pre_order[i] --> less_than(5) "i = $i"
                @fact tf.pre_order[i] --> greater_than(0) "i = $i"
            end
            FLSA.forward_dp_treepwl(tf)
            for i=1:4
                if i == root; continue; end
                @fact tf.nodes[i].lb --> roughly(round(vis.lb[i], 6)) "lb i = $i"
                @fact tf.nodes[i].ub --> roughly(round(vis.ub[i], 6)) "ub i = $i"
            end
        end
    end
end

facts("A 5 nodes, 4 level tree") do
    begin
        lambda = 1.0
        y = [3.5, 4.0, 2.5, 3.0, 7]
        n = length(y)
        root = n
        E = [(1,5), (2,3), (3,4), (4,5)]
        E = [Graphs.IEdge(i, u, v) for (i, (u,v)) in enumerate(E)]
        g = Graphs.edgelist(collect(1:n), E; is_directed=false)
        tn = FLSA.subtree(g, E, root)
        vis = FLSA.DPVisitor(y)
        x = FLSA.dp_tree_naive(y, lambda, tn, vis)

        parents = [5, 3, 4, 5, 5]
        tf = FLSA.PWLTree(parents, root, y, i->lambda)
        FLSA.forward_dp_treepwl(tf)
        for i=1:4
            @fact tf.nodes[i].lb --> roughly(round(vis.lb[i], 6)) "lb i = $i"
            @fact tf.nodes[i].ub --> roughly(round(vis.ub[i], 6)) "ub i = $i"
        end
    end
end


facts("A random 4x2 example") do
    begin
        srand(42)
        lambda = 1.0
        n1, n2 = 4, 2
        g = igraph(grid_graph(n1, n2))
        v = vertices(g)
        n = length(v)
        w = rand(size(g.edges))
        y = round(10*rand(n1*n2), 1)
        root = 1
        
        mst, wmst = kruskal_minimum_spantree(g, w)
        t = FLSA.subtree(g, mst, root)
        tn = FLSA.PWLTree(t, y)
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
