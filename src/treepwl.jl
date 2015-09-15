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

function next_min_event(t::PWLTree, i::Int)
    try
        pop_front!(t.nodes[i].pq)
    catch
        i = t.ubp[i]
        pop_front!(t.nodes[i].pq)
    end
end


function next_max_event(t::PWLTree, i::Int)
    try
        pop_back!(t.nodes[i].pq)
    catch
        i = t.lbp[i]
        pop_back!(t.nodes[i].pq)
    end
end


"""Delete ek and make sure, that the next event is in the queue"""
function step_min(t, ek)
    @debug "step_min($(ek.t)): $ek"
    # first extract this event
    ekk = next_min_event(t, ek.t)
    @assert ek == ekk "ek=$ek, ekk=$ekk"
    # now, ekk is the next event
    ekk = next_min_event(t, ek.t)
    @debug "step_min($(ek.t)): ekk = $ekk (will be deleted)"
    @debug "step_min($(ek.t)): Going from $(ek.t) to $(ekk.t)"
    t.lbp[ek.t] = t.lbp[ek.s]  # update bound parent
    @debug "step_min($(ek.t)): setting lbp of $(ek.t) to $(t.lbp[ek.s])"
    push!(t.nodes[ek.s].pq, ekk)
end


function step_max(t, ek)
    @debug "step_min($(ek.s)): $ek"
    ekk = next_max_event(t, ek.s)
    @assert ek == ekk "ek=$ek, ekk=$ekk"
    ekk = next_max_event(t, ek.s)
    @debug "step_min($(ek.s)): ekk = $ekk (will be deleted)"
    @debug "step_min($(ek.s)): Going from $(ek.s) to $(ekk.s)"
    t.ubp[ek.s] = t.ubp[ek.t] # update bound parent
    @debug "step_min($(ek.s)): setting ubp of $(ek.s) to $(t.ubp[ek.t])"
    if ek.t == ekk.s
        @debug "Source and target are the same ==> even already there"
        push!(t.nodes[ek.t].pq, ekk)
        return
    end

    # "merge" ekk into ek
    ek.s = ekk.s
    ek.slope  += ekk.slope
    ek.offset += ekk.offset
    ek.x = ekk.x


    # update pdeque
    pq = t.nodes[ek.t].pq
    @debug "step_min(): push! to $(ek.s): $ek"
    push!(pq, ek)
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
        step_min(t, ek)
        ek = find_min(t, v)
        @debug "lower_event!($v): ek  = $ek"
    end
    e.offset = e.offset
    @debug "lower_event!($v): final $e"
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
        step_max(t, ek)
        ek = find_max(t, v)
        @debug "upper_event!($v): ek  = $ek"
    end
    e.slope  = -e.slope
    e.offset = -e.offset
    @debug "upper_event!($v): final $e"
    push!(t.nodes[p].pq, e)
    t.nodes[v].ub = e.x
end


function print_min(t)
end


function print_max(t)
end


function print_tree(t)
    info("-"^70)
    for (i,n) in enumerate(t.nodes)
        @printf "\n((%d)): [% 3.1f,% 3.1f] {%d,%d}: " i n.lb n.ub t.lbp[i] t.ubp[i]
        println(join(map(string, n.pq), "\n" * (" "^26)))
    end
    info("-"^70)
end


function forward_dp_treepwl(t)
    for i in t.pre_order[end:-1:1]
        lower_event!(t, i)
        print_tree(t)
        upper_event!(t, i)
        print_tree(t)
        print_min(t)
        print_max(t)
    end
end
