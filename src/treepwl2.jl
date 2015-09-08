"""
Faster implementation
~~~~~~~~~~~~~~~~~~~~~
"""

macro debug(msg)
    :(info($msg))
end


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
    PWLNode() = new([], 1, 0, -Inf, +Inf)
end


showevents(n::PWLNode) = "$(n.events[n.a:n.b])"


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
        nodes = [PWLNode() for i in 1:n]
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

sort_events!(t, v::Int) = sort_events!(t.nodes[v])
sort_events!(n::PWLNode) = sort!(n.events[n.a:n.b], by=k->k.x)


"""Return lowest unprocessed event of node v or None if it does not exist"""
find_min_event(t, v::Int) = find_min_event(t.nodes[v])
function find_min_event(n::PWLNode)
    if n.a > length(n.events) || n.a > n.b
        throw(Inf)
    end
    return n.events[n.a]
end


"""Return lowest unprocessed event of node v or throw if it does not exist"""
find_max_event(t, v::Int) = find_max_event(t.nodes[v])
function find_max_event(n::PWLNode)
    if n.b <= 0 || n.b < n.a
        throw(-Inf)
    end
    return n.events[n.b]
end


"""Find the position of the lowest unprocessed event of node v"""
function find_min_x(t, v)
    try
        return find_min_event(t, v).x
    catch x
        return x
    end
end

function find_max_x(t, v)
    try
        return find_max_event(t, v).x
    catch x
        return x
    end
end


"""
Consume an event, by replacing the values in e.
Undefined behaviour if it does not exist.
Return x position of next event
"""
function step_min_event(t, e::Event)
    @debug "step($e)"
    n = t.nodes[e.t]
    ee = try
        find_min_event(n)
    catch
        e.t = t.parent[e.t]
        n = t.nodes[e.t]
        n.a += 1            # won't processed again
        n.events[n.a-1]     # = find_min_event(n)
    end
    n.a += 1 # won't processed again
    e.t = ee.t
    e.offset += ee.offset
    e.slope  += ee.slope
    sort_events!(t, e.s)
    return find_min_x(t, e.s)
end


function step_max_event(t, e::Event)
    n = t.nodes[e.t]
    @debug "step_max($n): BEGIN"
    ee = try
        find_max_event(n)
    catch
        e.s = t.parent[e.s]
        @debug "step_max(): no such event --> parent == $(e.s)"
        n = t.nodes[e.s]
        ee = n.events[n.b]  # = find_max_event(n)
        n.b -= 1            # won't processed again
        ee
    end
    @debug "step_max(): consuming $(ee)"
    n.b -= 1 # won't processed again
    e.s = ee.s
    e.offset -= ee.offset
    e.slope  -= ee.slope
    sort_events!(t, e.t)
    return find_min_x(t, e.t)
end


"""
Create a new event for v that corresponds to the new lower bound of v.
Requires child beeing processed
"""
function create_min_event(t, v, c::Float64)
    e = Event(v, v, 0.0, t.y[v], 1.0)
    e.offset += sum(map(i->t.lam(i), t.children[v]))
    forecast(e) = (c + e.offset) / e.slope
    e.x = forecast(e)
    xk = find_min_x(t, v)
    while e.x > xk
        xk = step_min_event(t, e)
        e.x = forecast(e)
    end
    t.nodes[v].lb = e.x
    e.offset -= t.lam(v)
    return e
end


function create_max_event(t, v, c::Float64)
    e = Event(v, v, 0.0, t.y[v], 1.0)
    e.offset -= sum(map(i->t.lam(i), t.children[v]))
    forecast(e) = (c + e.offset) / e.slope
    e.x = forecast(e)
    xk = find_max_x(t, v)
    @debug "create_max($v): events = $(showevents(t.nodes[v]))"
    @debug "create_max($v): x = $(e.x), xk = $xk, offset=$(e.offset), slope=$(e.slope)"
    while e.x < xk
        xk = step_max_event(t, e)
        e.x = forecast(e)
        @debug "create($v): events = $(showevents(t.nodes[v]))"
        @debug "create($v): x = $(e.x), xk = $xk, offset=$(e.offset), slope=$(e.slope)"
    end
    t.nodes[v].ub = e.x
    e.slope  = -e.slope
    e.offset = -e.offset - t.lam(v)
    return e
end



function forward_dp_treepwl(t)
    for i in t.pre_order[end:-1:1]
        n = t.nodes[i]
        childs = t.children[i]
        n.events = [[create_min_event(t, childs[c]) for c in childs]
                    [create_max_event(t, childs[c]) for c in childs]]
        n.a, n.b = 1, length(n.events)
        sort_events!(n)
    end
end
