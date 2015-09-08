"""Piecewise linear function over a tree"""

macro debug(msg)
#    :(info($msg))
end


type Event
    slope::Float64  # change of the slope, for lb, slope > 0, for ub, slope < 0
    offset::Float64 # change in the offset part
    x::Float64      # where it will happen
    i::Int          # node index
end

type PWLNode
    events::Vector{Event}
    a::Int          # index of lowest events, not occured yet
    b::Int          # index of highest events, not occured yet
    slope::Int      # current slope
    offset::Float64 # current offset
    v::Int          # node index
    lb::Float64     # lower bound (computed by clip_min!)
    ub::Float64     # upper bound (computed by clip_max!)
    function PWLNode(children, y::Vector{Float64}, v, lb, ub)
        events = Event[]
        resize!(events, 2*length(children))
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


"""Peek the next knot from the lower in a node v"""
function find_min_knot(t::PWLTree, v::Int)
    n = t.nodes[v]
    if n.a > length(n.events) || n.a > n.b
        return Inf
    end
    e = n.events[n.a]
    return e.x
end


"""Peek the next knot from the upper in a node v"""
function find_max_knot(t::PWLTree, v::Int)
    n = t.nodes[v]
    if n.b <= 0 || n.b < n.a
        return -Inf
    end
    e = n.events[n.b]
    return e.x
end


"""Find and extract the next knot from the lower in a node; adapt v.slope"""
function extract_min_knot!(t::PWLTree, v::Int)
    n = t.nodes[v]
    if n.a > length(n.events) || n.a > n.b
        return Inf
    end
    @debug "min_knot!($v): n.a=$(n.a), old offset = $(n.offset), slope = $(n.slope)"
    e = n.events[n.a]
    if e.i == v return Inf end
    n.slope  += e.slope
    n.offset += e.offset
    n.a      += 1
    @debug "min_knot!($v): consume event $e --> new offset = $(n.offset), slope = $(n.slope)"
    return find_min_knot(t, v)
end

"""Find and extract the next knot from the upper in a node; adapt v.slope"""
function extract_max_knot!(t::PWLTree, v::Int)
    n = t.nodes[v]
    if n.b <= 0 || n.b < n.a
        return -Inf
    end
    e = n.events[n.b]
    if e.i == v return -Inf end
    @debug "max_knot!($v): n.b=$(n.b), old offset = $(n.offset), slope = $(n.slope), e.offset = $(e.offset)"
    n.slope  -= e.slope
    n.offset -= e.offset
    n.b      -= 1
    @debug "max_knot!($v): consume event $e --> new offset = $(n.offset), slope = $(n.slope)"
    return find_max_knot(t, v)
end


"""Collect bounds from children and sort events"""
function prepare_events!(t::PWLTree, v::Int)
    node = t.nodes[v]
    for c in t.children[v]
        cn = t.nodes[c]
        node.events = [node.events, cn.events[cn.a:cn.b]]
    end
    node.a, node.b = 1, length(node.events)
    sort!(node.events, by=k->k.x)
    @debug "events($v): $(node.events)"
end

function push_event(t, v, e::Event)
    p = t.parent[v]
    if p != v
        push!(t.nodes[p].events, e)
    end
end


"""Clip node v from below until the derivative becomes c.
Return stop position x."""
function clip_min!(t::PWLTree, v::Int, c::Float64)
    node = t.nodes[v]
    node.slope = 1.0
    node.offset = sum([t.y[v], [t.lam(i) for i in t.children[v]]])
    forecast() = (c + node.offset) / node.slope
    x = forecast()
    xk = find_min_knot(t, v)
    @debug "clip_min!($v): BEGIN x = $x, xk = $xk"
    while x > xk
        @debug "clip_min!($v): forecast --> x = $x"
        xk = extract_min_knot!(t, v)
        x = forecast()
        @debug "clip_min!($v): node.offset = $(node.offset), c = $c, node.slope = $(node.slope)"
        @debug "clip_min!($v): x = $x, xk = $xk"
    end
    push_event(t, v, Event(node.slope, node.offset-t.lam(v), x, v))
    @debug "clip_min!($v): END   lb=$x, events --> $(node.events), [a,b] = [$(node.a), $(node.b)]"
    return x
end

"""Clip node v from above until the derivative becomes c.
Return stop position x."""
function clip_max!(t::PWLTree, v::Int, c::Float64)
    node = t.nodes[v]
    node.slope = 1.0
    node.offset = sum([t.y[v], [-t.lam(i) for i in t.children[v]]])
    forecast() = (c + node.offset) / node.slope
    x = forecast()
    xk = find_max_knot(t, v)
    @debug "clip_max!($v): BEGIN x = $x, xk = $xk"
    while x < xk
        @debug "clip_max!($v): node.offset = $(node.offset), c = $c, node.slope = $(node.slope)"
        @debug "clip_max!($v): x = $x, xk = $xk"
        xk = extract_max_knot!(t, v)
        x = forecast()
    end
    push_event(t, v, Event(-node.slope, -(node.offset+t.lam(v)), x, v))
    @debug "clip_max!($v): END  ub=$x, events --> $(node.events), [a,b] = [$(node.a), $(node.b)]"
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
    x = zeros(t.y)
    x[t.root] = clip_min!(t, t.root, 0.0)
    for i in t.pre_order[2:end]
        x[i] = clamp(x[t.parent[i]], t.nodes[i].lb, t.nodes[i].ub)
    end
    return x
end
