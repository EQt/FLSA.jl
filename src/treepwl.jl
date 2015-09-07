"""Piecewise linear function over a tree"""

type Event
    slope::Float64  # change of the slope, for lb, slope > 0, for ub, slope < 0
    x::Float64
    i::Int
end

type PWLNode
    events::Vector{Event}
    a::Int  # index of lowest events, not occured yet
    b::Int  # index of highest events, not occured yet
    slope::Int
    offset::Float64
    v::Int  # node index
    lb::Float64
    ub::Float64
    function PWLNode(children, y, v, lb, ub)
        events = [[Event(+1, lb[c], c) for c in children]
                  [Event(-1, ub[c], c) for c in children]]
        sort!(events, by=k->k.x)
        new(events, 1, length(events), 1, y, v, lb[v], ub[v])
    end
end

type PWLTree
    nodes::Vector{PWLNode}
    children::Vector{Vector{Int}}
    pre_order::Vector{Int}
    parent::Vector{Int}
    root::Int
    y::Vector{Float64}
    lam::Function

    """for convinience (testing): tree consisting of just one node"""
    function PWLTree(n::PWLNode, y=zeros(int(length(n.events)/2)), lam=i->1.0)
        new([n], Vector{Int}[[]], [1], [1], 1, [y, n.offset], lam)
    end
    

    function PWLTree(parents, root, y, lambda=i->1.0)
        n = length(parents)
        children = [Int[] for i=1:n]
        for (v,p) in enumerate(parents)
            if v != root push!(children[p], v) end
        end
        lb = [y[i] - (1+length(children[i]))*lambda(i) for i=1:n]
        ub = [y[i] + (1+length(children[i]))*lambda(i) for i=1:n]
        nodes = [PWLNode(children[i], y[i], i, lb, ub) for i in 1:n]
        pre_order = collect(1:n)
        return new(nodes, children, pre_order, parents, root, y, lambda)
    end

    PWLTree(t::ITreeSubGraph, y::Vector{Float64}, lambda=i->1.0) =
        PWLTree(t.parent, t.root, y, lambda)
end

"""Find and extract the next knot from the lower in a node; adapt v.slope"""
function min_knot!(t::PWLTree, v::Int)
    n = t.nodes[v]
    if n.a > length(n.events)
        return -Inf
    end
    e = n.events[n.a]
    n.slope += e.slope
    @assert(e.i in 1:length(t.y), "e.i = $(e.i), length(t.y) = $(length(t.y))")
    n.offset += sign(e.slope)*t.y[e.i]
    n.a += 1
    return e.x
end

"""Find and extract the next knot from the upper in a node; adapt v.slope"""
function max_knot!(t::PWLTree, v::Int)
    n = t.nodes[v]
    if n.b <= 0
        return Inf
    end
    e = n.events[n.b]
    n.slope += e.slope
    @assert(e.i in 1:length(t.y), "e.i = $(e.i), length(t.y) = $(length(t.y))")
    n.offset -= sign(e.slope)*t.y[e.i]
    n.b -= 1
    return e.x
end


"""Collect bounds from children and sort events"""
function prepare_events!(t::PWLTree, v::Int)
    node = t.nodes[v]
    for e in node.events
        if e.slope > 0
            e.x = t.nodes[e.i].lb
        else
            e.x = t.nodes[e.i].ub
        end
    end
    # append children's events
    for c in t.children[v]
        cn = t.nodes[c]
        cevents = cn.events[cn.a:cn.b]
        for e in cevents
            e.slope += sign(e.slope)*1
        end
        node.events = [node.events, cevents]
    end
    sort!(node.events, by=k->k.x)
end


"""Clip node v from below until the derivative becomes c.
Return stop position x."""
function clip_min!(t::PWLTree, v::Int, c::Float64)
    prepare_events(t, v)
    node = t.nodes[v]
    forecast() = (c + node.offset) / node.slope
    x = forecast()
    xk = min_knot!(t, v)
    while x > xk
        x = forecast()
    end
    return x
end

"""Clip node v from above until the derivative becomes c.
Return stop position x."""
function clip_max!(t::PWLTree, v::Int, c::Float64)
    prepare_events(t, v)
    node = t.nodes[v]
    forecast() = (c + node.offset) / node.slope
    x = forecast()
    xk = min_knot!(t, v)
    while x < xk
        x = forecast()
    end
    return x
end


"""Compute FLSA on a tree (fast algorithm)"""
function dp_tree_fast()
end

