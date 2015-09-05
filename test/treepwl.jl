facts("simple node (depth 1)") do
    children = Array{Int}[[], [], [], [1,2,3]]
    lb = [1, 2.5, 2, -Inf]
    ub = [3.2, 3, 4.5, Inf]
    y4 = 2.5
    context("min_knot!") do
        anode = PWLNode(children[4], y4, 4, lb, ub)
        @fact min_knot!(anode) --> 1.0
        @fact min_knot!(anode) --> 2.0
        @fact min_knot!(anode) --> 2.5
        @fact min_knot!(anode) --> 3.0
        @fact min_knot!(anode) --> 3.2
        @fact min_knot!(anode) --> 4.5
        @fact min_knot!(anode) --> Inf
    end
    context("max_knot!") do
        anode = PWLNode(children[4], y4, 4, lb, ub)
        @fact max_knot!(anode) --> 4.5
        @fact max_knot!(anode) --> 3.2
        @fact max_knot!(anode) --> 3.0
        @fact max_knot!(anode) --> 2.5
        @fact max_knot!(anode) --> 2.0
        @fact max_knot!(anode) --> 1.0
        @fact max_knot!(anode) --> -Inf
    end
end

facts("simple PWLTree") do
    parents = [4, 4, 4, 4]
    root = 4
    y = [2, 3, 2.5, 2.2]
    tree = FLSA.PWLTree(parents, root, y)
    @fact Set(tree.children[4]) --> Set(1,2,3)
    FLSA.prepare_events!(tree, 4)
    @pending map(k->k.x, tree.nodes[4].events) --> []
end


facts("A complete running example") do
    @pending begin
        srand(42)
        lambda = 0.1
        n1, n2 = 4, 2
        g = FLSA.grid_graph(n1, n2)
        v = vertices(g)
        n = length(v)
        w = rand(size(g.edges))
        y = rand(n1*n2)
        
        mst, wmst = kruskal_minimum_spantree(g, w)
        t = FLSA.subtree(g, mst, (1,1))
    end
end
