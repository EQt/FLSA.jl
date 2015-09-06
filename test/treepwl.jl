module TestFLSA
using FLSA
using FactCheck

facts("simple node (depth 1)") do
    children = Array{Int}[[], [], [], [1,2,3]]
    lb = [1, 2.5, 2, -Inf]
    ub = [3.2, 3, 4.5, Inf]
    y4 = 2.5
    context("min_knot!") do
        n = PWLNode(children[4], y4, 4, lb, ub)
        anode = PWLTree(n)
        @fact min_knot!(anode, 1) --> 1.0
        @fact min_knot!(anode, 1) --> 2.0
        @fact min_knot!(anode, 1) --> 2.5
        @fact min_knot!(anode, 1) --> 3.0
        @fact min_knot!(anode, 1) --> 3.2
        @fact min_knot!(anode, 1) --> 4.5
        @fact min_knot!(anode, 1) --> -Inf
    end
    context("max_knot!") do
        anode = PWLTree(PWLNode(children[4], y4, 4, lb, ub))
        @fact max_knot!(anode, 1) --> 4.5
        @fact max_knot!(anode, 1) --> 3.2
        @fact max_knot!(anode, 1) --> 3.0
        @fact max_knot!(anode, 1) --> 2.5
        @fact max_knot!(anode, 1) --> 2.0
        @fact max_knot!(anode, 1) --> 1.0
        @fact max_knot!(anode, 1) --> Inf
    end
end

facts("simple PWLTree (depth 1)") do
    parents = [4, 4, 4, 4]
    root = 4
    y = [2, 3, 2.5, 2.2]
    tree = FLSA.PWLTree(parents, root, y, i->1.0)
    @fact Set(tree.children[4]) --> Set(1,2,3)
    for i=1:3
        @fact tree.children[i] --> isempty "$tree.children[i]"
    end
    FLSA.prepare_events!(tree, 4)
    events = [[y[1:3] .+ 1], [y[1:3] .- 1]]
    @fact map(k->k.x, tree.nodes[4].events) --> sort!(events)
end

facts("simple PWLTree (depth 2)") do
    parents = [3, 4, 4, 4]
    root = 4
    y = [2, 3, 2.5, 2.2]
    tree = FLSA.PWLTree(parents, root, y, i->1.0)
    for i=1:2
        @fact tree.children[i] --> isempty "$tree.children[i]"
    end
    @fact Set(tree.children[3]) --> Set(1)
    @fact Set(tree.children[4]) --> Set(2,3)

    FLSA.prepare_events!(tree, 3)
    @fact map(k->k.x, tree.nodes[3].events) --> [y[1] - 1, y[1] + 1]
    FLSA.prepare_events!(tree, 4)
    events = [[y[1:2] .+ 1], [y[1:2] .- 1], y[3] + 2, y[3] - 2]
    @fact map(k->k.x, tree.nodes[4].events) --> sort!(events)
    @fact map(k->k.slope, filter(k->k.i == 1, tree.nodes[4].events)) --> [2,-2]
    @fact map(k->k.slope, filter(k->k.i != 1, tree.nodes[4].events)) --> [1,1,-1,-1]
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

end
