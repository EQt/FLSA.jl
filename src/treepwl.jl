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
        n = length(children)
        children = fill(Int[], n)
        for (v,p) in enumerate(parents)
            push!(children[p], v)
        end
        nodes = [PWLNode(children[i], y[i], i, -Inf, Inf) for i in 1:n]
        return new(nodes, children, [], parents, root, y)
    end
end

"""Find and extract the next knot from the lower in a node"""
function min_knot!(n::PWLNode)
    if n.a > length(n.events)
        return Inf
    end
    e = n.events[n.a]
    n.a += 1
    return e.x
end

"""Find and extract the next knot from the upper in a node"""
function max_knot!(n::PWLNode)
    if n.b <= 0
        return -Inf
    end
    e = n.events[n.b]
    n.b -= 1
    return e.x
end
