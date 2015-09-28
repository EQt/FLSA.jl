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

forecast(e::Event) = - e.offset/e.slope
set_forecast!(e::Event) = (e.x = forecast(e))


"""Manage the events of a node"""
type PWLNode
    lb::Float64             # lower bound (computed by create_min_event)
    ub::Float64             # upper bound (computed by create_max_event)
    pq::PDeQue{Event}       # events
    PWLNode() = new(-Inf, +Inf, PDeQue{Event}([], e->e.x))
end


type PWLTree
    nodes::Vector{PWLNode}
    children::Vector{Vector{Int}}
    pre_order::Vector{Int}
    parent::Vector{Int}
    lbp::Vector{Int}
    ubp::Vector{Int}
    root::Int
    y::Vector{Float64}
    lam::Function

    function PWLTree(parents::Vector{Int}, root::Int, y::Vector{Float64},
                     lambda=i->1.0)
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
        return new(nodes, children, pre_order,
                   parents, copy(parents), copy(parents), root, y, lambda)
    end
    PWLTree(t::ITreeSubGraph, y::Vector{Float64}, lambda=i->1.0) =
        PWLTree(t.parent, t.root, y, lambda)
end


find_min(t, i) = try front(t.nodes[i].pq) catch Event(i, i, +Inf, 0, 0) end
find_max(t, i) = try back(t.nodes[i].pq)  catch Event(i, i, -Inf, 0, 0) end


"""Delete one front event of v and replace it by next, if not already enqueued"""
function step_min(t, v)
    pq = t.nodes[v].pq
    e = pop_front!(pq)
    @debug "step_min($v):" * " "^11 * "$e (will be deleted)"
    if e.slope > 0 # lower bound event
        t.lbp[e.t] = t.parent[v]
        @debug "bounds($v): lbp[$(e.t)] --> $(t.parent[v])"
    else
        t.ubp[e.s] = t.parent[v]
        @debug "bounds($v): ubp[$(e.s)] --> $(t.parent[v])"
    end
    u = e.t
    @debug "step_min($v): u=$u"
    pq_u = t.nodes[u].pq
    if isempty(pq_u)
        if u == v
            @debug "step_min($v): u == v"
            return
        end
        # test parent
        uold = u
        u = t.parent[u]
        pq_u = t.nodes[u].pq
        if !isempty(pq_u)
            @debug "step_min($v): u=$u (parent)"
        else
            # test bound parent
            u = t.ubp[uold]
            if u == t.parent[v] # don't want to steel parent's event
                @debug "step_min($v): finished"
                return
            end
            @debug "step_min($v): u=$u"
            pq_u = t.nodes[u].pq
            if isempty(pq_u)
                @debug "step_min($v): was last"
                return
            end
        end
        e = pop_front!(pq_u)
    else
        @debug "step_min($v): u=$u"
        e = pop_front!(pq_u)
    end
    @debug "step_min($v): moving    $e"
    push!(pq, e)
end


function step_max(t, v)
    pq = t.nodes[v].pq
    e = pop_back!(pq)
    @debug "step_max($v):" * " "^11 * "$e (will be deleted)"
    if e.slope > 0 # lower bound event
        t.lbp[e.t] = t.parent[v]
        @debug "bounds($v): lbp[$(e.t)] --> $(t.parent[v])"
    else
        t.ubp[e.s] = t.parent[v]
        @debug "bounds($v): ubp[$(e.s)] --> $(t.parent[v])"
    end
    u = e.s
    @debug "step_max($v): u=$u"
    pq_u = t.nodes[u].pq
    if isempty(pq_u)
        if u == v
            @debug "step_max($v): u == v"
            return
        end
        # test parent
        uold = u
        u = t.parent[u]
        pq_u = t.nodes[u].pq
        if !isempty(pq_u)
            @debug "step_min($v): u=$u (parent)"
        else
            # test bound parent
            u = t.lbp[uold]
            if u == t.parent[v] # don't want to steel parent's event
                @debug "step_min($v): finished"
                return
            end
            @debug "step_min($v): u=$u"
            pq_u = t.nodes[u].pq
            if isempty(pq_u)
                @debug "step_min($v): was last"
                return
            end
        end
        e = pop_back!(pq_u)
    else
        @debug "step_max($v): u=$u"
        e = pop_back!(pq_u)
    end
    @debug "step_max($v): moving    $e"
    push!(pq, e)
end


function step_min2(t, v)
    pq = t.nodes[v].pq
    e = pop_front!(pq)
    u = e.t # always the node where we go to
    pq_u = t.nodes[u].pq
    if isempty(pq_u)
        # next event moved...
        u = t.ubp[u]
        pq_u = t.nodes[u].pq
        e = pop_front(pq_u)       
    else
        e = pop_front(pq_u)
        t.lbp[u] = v
    end
    push!(pq, e)
end


function step_max_old(t, v)
    u = 6
    info("u = $u, $(t.nodes[u].pq.elements)")
    u = 5
    info("u = $u, $(t.nodes[u].pq.elements)")
    pq = t.nodes[v].pq
    e = pop_back!(pq)
    @debug "step_max($v):" * " "^11 * "$e (will be deleted)"
    u = e.s
    pq_u = t.nodes[u].pq
    if isempty(pq_u)
        t.ubp[u] = t.parent[v]  # update bound parent
        @debug "step_max($v): setting lbp of $u to $(t.ubp[u])"
        @debug "step_max($v): Going from $u to lbp ==> $(t.lbp[u])"
        u = t.lbp[u]
        if u == t.parent[v]
            @debug "step_max($v): finished"
            return
        end
        pq_u = t.nodes[u].pq
    end
    try
        #print_tree(t)
        if isempty(pq_u)
            return
        end
        e = pop_back!(pq_u)
        @debug "step_max($v): moving    $e"
        if isempty(t.nodes[e.t].pq)
            t.lbp[e.t] = v
        end
        if isempty(t.nodes[e.s].pq)
            t.ubp[e.s] = v
        end
        push!(pq, e)
    catch
        error("u = $u, $(t.nodes[u].pq.elements)")
    end
end

step_min2(t, i) = pop_front!(t.nodes[i].pq)
step_max2(t, i) = pop_back!(t.nodes[i].pq)
    

function prepare_node(t, v::Int)
    n = t.nodes[v]
    for c in t.children[v]
        nc = t.nodes[c]
        elem = nc.pq.elements
        for e in elem
            push!(n.pq, e)
        end
        resize!(elem, 0)
    end
    sort!(n.pq.elements, by=k->k.x)
end


function lower_event!(t, v::Int, c::Float64=-t.lam(v))
    p = t.parent[v]
    noffset = -t.y[v] - c - sum(map(i->t.lam(i), t.children[v]))
    e = Event(p, v, 0.0, noffset, 1.0)
    set_forecast!(e)
    ek = find_min(t, v)
    @debug "lower_event!($v): e   = $e"
    @debug "lower_event!($v): ek  = $ek"
    while ek.x < e.x
        e.offset += ek.offset
        e.slope  += ek.slope
        e.t       = ek.t
        set_forecast!(e)
        @debug "lower_event!($v): e   = $e"
        step_min(t, v)
        ek = find_min(t, v)
        @debug "lower_event!($v): ek  = $ek"
    end
    e.offset = e.offset
    @debug "lower_event!($v): final $e (pushed to $p)"
    push!(t.nodes[p].pq, e)
    t.nodes[v].lb = e.x
end


function lower_event2!(t, v::Int, c::Float64=-t.lam(v))
    p = t.parent[v]
    noffset = -t.y[v] - c - sum(map(i->t.lam(i), t.children[v]))
    e = Event(p, v, 0.0, noffset, 1.0)
    set_forecast!(e)
    ek = find_min(t, v)
    @debug "lower_event!($v): e   = $e"
    @debug "lower_event!($v): ek  = $ek"
    while ek.x < e.x
        e.offset += ek.offset
        e.slope  += ek.slope
        e.t       = ek.t
        set_forecast!(e)
        @debug "lower_event!($v): e   = $e"
        step_min2(t, v)
        ek = find_min(t, v)
        @debug "lower_event!($v): ek  = $ek"
    end
    e.offset = e.offset
    @debug "lower_event!($v): final $e (pushed to $p)"
    push!(t.nodes[p].pq, e)
    t.nodes[v].lb = e.x
end



function upper_event!(t, v::Int, c::Float64=+t.lam(v))
    p = t.parent[v]
    noffset = - t.y[v] - c + sum(map(i->t.lam(i), t.children[v]))
    e = Event(v, p, 0.0, noffset, 1.0)
    set_forecast!(e)
    ek = find_max(t, v)
    @debug "upper_event!($v): e   = $e"
    @debug "upper_event!($v): ek  = $ek"
    while ek.x > e.x
        e.offset -= ek.offset
        e.slope  -= ek.slope
        e.s       = ek.s
        set_forecast!(e)
        @debug "upper_event!($v): e   = $e"
        step_max(t, v)
        ek = find_max(t, v)
        @debug "upper_event!($v): ek  = $ek"
    end
    e.slope  = -e.slope
    e.offset = -e.offset
    @debug "upper_event!($v): final $e (pushed to $p)"
    push!(t.nodes[p].pq, e)
    t.nodes[v].ub = e.x
end


function upper_event2!(t, v::Int, c::Float64=+t.lam(v))
    p = t.parent[v]
    noffset = - t.y[v] - c + sum(map(i->t.lam(i), t.children[v]))
    e = Event(v, p, 0.0, noffset, 1.0)
    set_forecast!(e)
    ek = find_max(t, v)
    @debug "upper_event!($v): e   = $e"
    @debug "upper_event!($v): ek  = $ek"
    while ek.x > e.x
        e.offset -= ek.offset
        e.slope  -= ek.slope
        e.s       = ek.s
        set_forecast!(e)
        @debug "upper_event!($v): e   = $e"
        step_max2(t, v)
        ek = find_max(t, v)
        @debug "upper_event!($v): ek  = $ek"
    end
    e.slope  = -e.slope
    e.offset = -e.offset
    @debug "upper_event!($v): final $e (pushed to $p)"
    push!(t.nodes[p].pq, e)
    t.nodes[v].ub = e.x
end



function print_min(t, v)
    t = deepcopy(t)
    warn("BEGIN($v)")
    while true
        n = t.nodes[v]
        e = find_min(t, v)
        warn("AT($(e.x)):  $(e.s) ----[ $(e.offset) / $(e.slope) ] ------>  $(e.t)")
        if !isfinite(e.x) break end
        step_min(t, v)
    end
end


function print_max(t, v)
end


function print_tree(t)
    info("-"^70)
    for (i,n) in enumerate(t.nodes)
        @printf "\n((%d)): [% 3.1f,% 3.1f] {%d,%d}: " i n.lb n.ub t.lbp[i] t.ubp[i]
        println(join(map(string, n.pq), "\n" * (" "^26)))
    end
    info("-"^70)
end

countevents(t) = sum(map(i->length(t.nodes[i].pq.elements), 1:length(t.nodes)))

function forward_dp_treepwl(t)
    t2 = deepcopy(t)
    for i in t.pre_order[end:-1:1]
        prepare_node(t2, i)
        lower_event2!(t2, i)
        @debug "TTTT 22222"
        print_tree(t2)

        lower_event!(t, i)
        @debug "TTTT 11111"
        print_tree(t)

        countevents(t) == countevents(t2) || begin print_tree(t); print_tree(t2); error("STOP_LB($i)") end

        upper_event2!(t2, i)
        @debug "TTTT 22222"
        print_tree(t2)

        upper_event!(t, i)
        @debug "TTTT 11111"
        print_tree(t)

        countevents(t) == countevents(t2) || begin print_tree(t); print_tree(t2); error("STOP_UB($i)") end
    end
end
