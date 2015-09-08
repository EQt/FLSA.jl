"""
Faster implementation
~~~~~~~~~~~~~~~~~~~~~


"""

"""Record what is happening, when a knot of the PWL is hit"""
type Event
    s::Int          # source node
    t::Int          # target node
    x::Float64      # position
    offset::Float64 # delta offset
    slope::Float64  # delta slope
end

"""Manage the events of a node"""
type PWLNode
    events::Vector{Event}   # should have fixed size of 2|children|
    a::Int                  # index of lowest, unprocessed event
    b::Int                  # index of highest, unprocessed event
    lb::Float64             # lower bound (computed by create_min_event)
    ub::Float64             # upper bound (computed by create_max_event)
end


type PWLTree
    nodes::Vector{PWLNode}
    children::Vector{Vector{Int}}
    pre_order::Vector{Int}
    parent::Vector{Int}
    root::Int
    y::Vector{Float64}
    lam::Function

    function PWLTree(parents::Vector{Int}, root::Int, y::Vector{Float64}, lambda=i->1.0)
        n = length(parents)
        children = [Int[] for i=1:n]
        for (v,p) in enumerate(parents)
            if v != root push!(children[p], v) end
        end
        lb = [y[i] - (1+length(children[i]))*lambda(i) for i=1:n]
        ub = [y[i] + (1+length(children[i]))*lambda(i) for i=1:n]
        nodes = [PWLNode(children[i], y, i, lb, ub) for i in 1:n]
        pre_order = zeros(Int, n)
        stack = [root]
        nr = 1
        while !isempty(stack)
            v = pop!(stack)
            pre_order[nr] = v
            nr += 1
            append!(stack, children[v])
        end
        return new(nodes, children, pre_order, parents, root, y, lambda)
    end

    PWLTree(t::ITreeSubGraph, y::Vector{Float64}, lambda=i->1.0) =
        PWLTree(t.parent, t.root, y, lambda)
end


"""Create a new event for v that corresponds to the new lower bound of v"""
function create_min_event(t, v, c::Float64)
    e = Event(v, v, 0.0, t.y[v], 1.0)
    e.offset += sum(map(i->t.lam(i), t.children[v]])
    forecast(e) = (c + e.offset) / e.slope
    e.x = forecast(e)
    xk = find_min_x(t, v)
    while x > xk
        xk = step_min_event(t, e)
        e.x = forecast(e)
    end
    t.nodes[v].lb = e.x
    return e
end

"""Find the lowest unprocessed event of node v or return None if it does not exist"""
find_min_event(t, v::Int) = find_min_event(t.nodes[v])
function find_min_event(n::PWLNode)
    if n.a > length(n.events) || n.a > n.b
        return None
    end
    return n.events[n.a]
end


"""Find the position of the lowest unprocessed event of node v"""
function find_min_x(t, v)
    e = find_min_event(t, v)
    e == None || return Inf
    return e.x
end


"""
Consume an event, by replacing the values in e.
Undefined behaviour if it does not exist.
Return x position of next event
"""
function step_min_event(t, e::Event)
    n = t.nodes[e.t]
    n.a += 1 # won't processed again
    ee = find_min_event(n)
    if ee == None
        e.t = t.parent[e.t]
        n = t.nodes[e.t]
        ee = n.events[n.a]  # = find_min_event(n)
        n.a += 1 # won't processed again
    end
    e.offset += ee.offset
    e.slope  += ee.slope
    sort_events!(t, e.s)
    return find_min_x(t, e.s)
end


function sort_events!(t, v::Int)
    n = t.nodes[v]
    sort!(n.events[n.a:n.b], by=k->k.x)
end
