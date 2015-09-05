"""Piecewise linear function over a tree"""

type Event
    islb::Bool
    x::Float64
    i::Int
end

type PWLNode
    events::Vector{Event}
    a::Int  # index of lowest events, not occured yet
    b::Int  # index of highest events, not occured yet
    slope::Int
    y::Float64
    v::Int  # node index
    lb::Float64
    ub::Float64
    function PWLNode(children, y, v, lb, ub)
        events = [[Event(true,  lb[c], c) for c in children]
                  [Event(false, ub[c], c) for c in children]]
        sort!(events, by=k->k.x)
        new(events, 1, length(events), 1, y, v, -Inf, +Inf)
    end
end

type PWLTree
    nodes::Vector{PWLNode}
    children::Vector{Vector{Int}}
    pre_order::Vector{Int}
    parent::Vector{Int}
    root::Int
    y::Vector{Float64}
    function PWLTree(parents, root, y)
        n = length(parents)
        children = fill(Int[], n)
        for (v,p) in enumerate(parents)
            if v != root push!(children[p], v) end
        end
        lb = fill(-Inf, n)
        ub = fill(+Inf, n)
        nodes = [PWLNode(children[i], y[i], i, lb, ub) for i in 1:n]
        pre_order = collect(1:n)
        return new(nodes, children, pre_order, parents, root, y)
    end
end

"""Find and extract the next knot from the lower in a node"""
function min_knot!(t::PWLTree, v::Int)
    n = t.nodes[v]
    if n.a > length(n.events)
        return -Inf
    end
    e = n.events[n.a]
    n.a += 1
    return e.x
end

"""Find and extract the next knot from the upper in a node"""
function max_knot!(t::PWLTree, v::Int)
    n = t.nodes[v]
    if n.b <= 0
        return Inf
    end
    e = n.events[n.b]
    n.b -= 1
    return e.x
end


"""Collect bounds from children and sort events"""
function prepare_events!(t::PWLTree, v::Int)
    node = t.nodes[v]
    for e in node.events
        if e.islb
            e.x = t.nodes[e.i].lb
        else
            e.x = t.nodes[e.i].ub
        end
    end
    sort!(node.events, by=k->k.x)
end


"""Clip node v from below until the derivative becomes c.
Return stop position x."""
function clip_min!(t::PWLTree, v::Int, c::Float64)
    prepare_events(t, v)
    node = t.nodes[v]
    oldx = node.y
    x = node.y - c # forecast
    df = 0.0
    xe = min_knot!(t, v)
    while x < xe
        df += node.slope * (x - oldx)
        xe = min_knot(t, v)
    end
    return x
end
