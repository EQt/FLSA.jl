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
    lb::Float64             # lower bound (computed by create_min_event)
    ub::Float64             # upper bound (computed by create_max_event)
    minevs::Vector{Event}   # should have fixed size of |children|
    maxevs::Vector{Event}   # should have fixed size of |children|
    PWLNode() = new(-Inf, +Inf, [], [])
end

forecast(e::Event, c) = (c - e.offset)/e.slope

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
sort_events!(n::PWLNode) = begin
    sort!(n.minevs, by=k->k.x, rev=false)
    sort!(n.maxevs, by=k->k.x, rev=true)
end


"""Return lowest unprocessed event of node v or None if it does not exist"""
find_min_event(t, v::Int) = find_min_event(t.nodes[v])
find_min_event(n::PWLNode) = try n.minevs[1] catch throw(+Inf) end


"""Return lowest unprocessed event of node v or throw if it does not exist"""
find_max_event(t, v::Int) = find_max_event(t.nodes[v])
find_max_event(n::PWLNode) = try n.maxevs[1] catch throw(-Inf) end


"""Find the position of the lowest unprocessed event of node v"""
find_min_x(t, v) = try find_min_event(t, v).x catch x return x end
find_max_x(t, v) = try find_max_event(t, v).x catch x return x end


"""
Consume an event, by replacing the values in e.
Undefined behaviour if it does not exist.
Return x position of next event
"""
function step_min_event(t, e::Event)
    n = t.nodes[e.t]
    @debug "step_min($e): n=$n"
    try
        ee = shift!(n.minevs)
        @debug "step_min($(e.t)): deleted $ee"
        @debug "step_min($(e.t)): n.minevs = $(n.minevs)"
        @debug "step_min($(e.t)): setting s from $(e.t) to $(ee.t)"
        e.t = ee.t
        e.offset += ee.offset
        e.slope  += ee.slope
        e.x = ee.x
        sort_events!(t, e.s)
        return find_min_x(t, e.s)
    catch
        error("This should not happen!")
    end
end


function step_max_event(t, e::Event)
    n = t.nodes[e.s]
    @debug "step_max($e): n=$n"
    try
        ee = shift!(n.maxevs)
        @debug "step_max($(e.s)): deleted $ee"
        @debug "step_max($(e.s)): n.maxevs = $(n.maxevs)"
        @debug "step_max($(e.s)): setting s from $(e.s) to $(ee.s)"
        e.s = ee.s
        e.offset += ee.offset
        e.slope  += ee.slope
        e.x = ee.x
        sort_events!(t, e.t)
        return find_max_x(t, e.t)
    catch
        error("This should not happen!")
    end
end


"""
Create a new event for v that corresponds to the new lower bound of v.
Requires child beeing processed.
Insert this event also to the corresponding child node!
"""
function create_min_event(t, v::Int, c::Float64=-t.lam(v))
    e = Event(t.parent[v], v, 0.0, -t.y[v], 1.0)
    e.offset -= sum(map(i->t.lam(i), t.children[v]))
    e.x = forecast(e, c)
    xk = find_min_x(t, v)
    @debug "create_min($v): starting e=$e, xk=$xk [y=$(t.y[v])]"
    while e.x > xk
        xk = step_min_event(t, e)
        e.x = forecast(e, c)
        @debug "create_min($v): forcast e=$e, xk=$xk"
    end
    t.nodes[v].lb = e.x
    e.offset += t.lam(t.parent[v])    # compensate status of parent
    push!(t.nodes[e.t].maxevs, e)
    @debug "create_min($v): adding e --> node[$(e.t)] = $(t.nodes[e.t])"
    return e
end

function create_max_event(t, v::Int, c::Float64=t.lam(v))
    e = Event(v, t.parent[v], 0.0, -t.y[v], 1.0)
    e.offset += sum(map(i->t.lam(i), t.children[v]))
    e.x = forecast(e, c)
    xk = find_max_x(t, v)
    @debug "create_max($v): starting e=$e, xk=$xk [y=$(t.y[v])]"
    while e.x < xk
        xk = step_max_event(t, e)
        e.x = forecast(e, c)
        @debug "create_max($v): forcast e=$e, xk=$xk"
    end
    t.nodes[v].ub = e.x
    e.slope  = -e.slope
    e.offset = -e.offset + t.lam(t.parent[v])
    push!(t.nodes[e.s].minevs, e)
    @debug "create_max($v): adding e --> node[$(e.s)] = $(t.nodes[e.s])"
    return e
end


function print_min_chain(t, v::Int)
    t = deepcopy(t)
    warn("BEGIN($v)")
    try
        while true
            n = t.nodes[v]
            e = find_min_event(t, v)
            warn("AT($(e.x)):  $(e.s) ----[ $(e.offset) / $(e.slope) ] ------>  $(e.t)")
            xk = step_min_event(t, e)
            v = e.t
        end
    catch y
        if y == Inf
            warn("END($v)")
        else
            error(y)
        end
    end
end


function print_max_chain(t, v::Int)
    t = deepcopy(t)
    warn("BEGIN($v)")
    try
        while true
            n = t.nodes[v]
            e = find_max_event(t, v)
            warn("AT($(e.x)):  $(e.s) <------[ $(e.offset) / $(e.slope) ] -----  $(e.t)")
            xk = step_max_event(t, e)
            v = e.s
        end
    catch y
        if y == -Inf
            warn("END($v)")
        else
            error(y)
        end
    end
end

function forward_dp_treepwl(t)
    for i in t.pre_order[end:-1:1]
        n = t.nodes[i]
        childs = t.children[i]
        n.minevs = [create_min_event(t, c) for c in childs]
        n.maxevs = [create_max_event(t, c) for c in childs]
        @debug "forward($i): nodes is $n"
        sort_events!(n)
        print_min_chain(t, i)
        print_max_chain(t, i)
    end
end
