"""Piecewise linear function over a tree"""

macro debug(msg)
    :(info($msg))
end


type Event
    slope::Float64  # change of the slope, for lb, slope > 0, for ub, slope < 0
    offset::Float64 # change in the offset part
    x::Float64
    i::Int
end

type PWLNode
    events::Vector{Event}
    a::Int  # index of lowest events, not occured yet
    b::Int  # index of highest events, not occured yet
    slope::Int      # current slope
    offset::Float64 # current offset
    v::Int  # node index
    lb::Float64
    ub::Float64
    function PWLNode(children, y::Vector{Float64}, v, lb, ub)
        events = [[Event(+1, +y[c], lb[c], c) for c in children]
                  [Event(-1, -y[c], ub[c], c) for c in children]]
        sort!(events, by=k->k.x)
        new(events, 1, length(events), 1, y[v], v, lb[v], ub[v])
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
    

    function PWLTree(parents, root, y::Vector{Float64}, lambda=i->1.0)
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

"""Find and extract the next knot from the lower in a node; adapt v.slope"""
function min_knot!(t::PWLTree, v::Int)
    n = t.nodes[v]
    if n.a > length(n.events)
        return Inf
    end
    @debug "min_knot!($v): n.a=$(n.a), old offset = $(n.offset)"
    e = n.events[n.a]
    if e.i == v return Inf end
    n.slope  += e.slope
    n.offset += e.offset
    n.a      += 1
    @debug "min_knot!($v): consume event e.i = $(e.i) --> new offset = $(n.offset)"
    return e.x
end

"""Find and extract the next knot from the upper in a node; adapt v.slope"""
function max_knot!(t::PWLTree, v::Int)
    n = t.nodes[v]
    if n.b <= 0
        return -Inf
    end
    e = n.events[n.b]
    if e.i == v return -Inf end
    @debug "max_knot!($v): n.b=$(n.b), old offset = $(n.offset)"
    n.slope  -= e.slope
    n.offset -= e.offset
    n.b      -= 1
    @debug "max_knot!($v): consume event e.i = $(e.i) --> new offset = $(n.offset)"
    return e.x
end


"""Collect bounds from children and sort events"""
function prepare_events!(t::PWLTree, v::Int)
    node = t.nodes[v]
    node.events = []
    for c in t.children[v]
        cn = t.nodes[c]
        node.events = [node.events, cn.events[cn.a:cn.b]]
    end
    node.a, node.b = 1, length(node.events)
    sort!(node.events, by=k->k.x)
    @debug "events($v): $([(e.i, e.x) for e in node.events])"
end


"""Clip node v from below until the derivative becomes c.
Return stop position x."""
function clip_min!(t::PWLTree, v::Int, c::Float64)
    node = t.nodes[v]
    node.slope = 1.0
    node.offset = sum([t.y[v], [t.lam(i) for i in t.children[v]]])
    @debug "clip_min!($v): node.offset = $(node.offset), c = $c, y=$(t.y[v])"
    forecast() = (c + node.offset) / node.slope
    x = forecast()
    xk = min_knot!(t, v)
    @debug "clip_min!($v): Starting at x = $x, xk = $xk"
    while x > xk
        x = forecast()
        @debug "clip_min!($v): forecast --> x = $x"
        xk = min_knot!(t, v)
        @debug "clip_min!($v): node.offset = $(node.offset), c = $c, node.slope = $(node.slope)"
        @debug "clip_min!($v): x = $x, xk = $xk"
    end
    node.a -= 1
    e = Event(1, t.y[v] - t.lam(v), x, v)
    if node.a <= 0
        node.a = 1
        node.b = max(node.b, node.a)
        unshift!(node.events, e)
    else
        node.events[node.a] = e
    end
    @debug "clip_min!($v): events --> $(node.events), [a,b] = [$(node.a), $(node.b)]"
    return x
end

"""Clip node v from above until the derivative becomes c.
Return stop position x."""
function clip_max!(t::PWLTree, v::Int, c::Float64)
    node = t.nodes[v]
    node.slope = 1.0
    node.offset = sum([t.y[v], [-t.lam(i) for i in t.children[v]]])
    @debug "clip_max!($v): node.offset = $(node.offset), c = $c, y=$(t.y[v])"
    forecast() = (c + node.offset) / node.slope
    x = forecast()
    xk = max_knot!(t, v)
    @debug "clip_max!($v): x = $x, xk = $xk"
    while x < xk
        @debug "clip_max!($v): node.offset = $(node.offset), c = $c, node.slope = $(node.slope)"
        @debug "clip_max!($v): x = $x, xk = $xk"
        x = forecast()
        xk = max_knot!(t, v)
    end
    node.b += 1
    e = Event(-1, -(t.y[v] - t.lam(v)), x, v)
    if node.b > length(node.events)
        node.a = min(node.a, node.b)
        push!(node.events, e)
    else
        node.events[node.b] = e
    end
    @debug "clip_max!($v): events --> $(node.events), [a,b] = [$(node.a), $(node.b)]"
    return x
end


"""Compute FLSA on a tree (fast algorithm)"""
function dp_treepwl(t::PWLTree)
    forward_dp_treepwl(t)
    backtrace_dp_treepwl(t)
end


function forward_dp_treepwl(t)
    for i in t.pre_order[end:-1:1]
        n = t.nodes[i]
        prepare_events!(t, i)
        n.lb = clip_min!(t, i, -t.lam(i))
        n.ub = clip_max!(t, i, +t.lam(i))
        @debug "forward_dp($i): [$(n.lb), $(n.ub)], events --> $(n.events[n.a:n.b])"
    end
end


function backtrace_dp_treepwl(t::PWLTree)
    x = zeros(y)
    x[t.root] = clip_min!(t, t.root, 0)
    for i in t.pre_order[2:end]
        x[i] = clamp(x[t.parent[i]], t.nodes[i].lb, t.nodes[i].ub)
    end
    return x
end
