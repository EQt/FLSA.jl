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

forecast(e::Event, c) = (c - e.offset)/e.slope
set_forecast!(e::Event, c) = (e.x = forecast(e, c))


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
        return new(nodes, children, pre_order, parents, root, y, lambda)
    end
    PWLTree(t::ITreeSubGraph, y::Vector{Float64}, lambda=i->1.0) =
        PWLTree(t.parent, t.root, y, lambda)
end


find_min(t, i) = try front(t.nodes[i].pq) catch Event(i, i, +Inf, 0, 0) end
find_max(t, i) = try back(t.nodes[i].pq)  catch Event(i, i, -Inf, 0, 0) end

function step_min(t, ek)
    @debug "step_min($(ek.t)): $ek"
    n = t.nodes[ek.t]
    ekk = pop_front!(n.pq)
    @debug "step_min($(ek.t)): ekk = $ekk (will be deleted)"
    @debug "step_min(): Going from $(ek.t) to $(ekk.t)"
    ek.t = ekk.t
    ek.slope  += ekk.slope
    ek.offset += ekk.offset
    ek.x = ekk.x
    pq = t.nodes[ek.s].pq
    ekk = pop_front!(pq)
    if abs(ek.slope) <= 1e-6
        @debug "step_min(): trying to delete ek = $ek"
        @assert ek == ekk
        return
    else
        push!(pq, ek)
    end
end


function step_max(t, ek)
    @debug "step_max($(ek.s)): $ek"
    n = t.nodes[ek.s]
    ekk = pop_back!(n.pq)
    @debug "step_max($(ek.s)): ekk = $ekk (will be deleted)"
    @debug "step_max(): Going from $(ek.s) to $(ekk.s)"
    ek.s = ekk.s
    ek.slope  += ekk.slope
    ek.offset += ekk.offset
    ek.x = ekk.x
    pq = t.nodes[ek.t].pq
    ekk = pop_back!(pq)
    if abs(ek.slope) <= 1e-6
        @debug "step_min(): trying to delete ek = $ek"
        @assert ek == ekk
        return
    else
        push!(pq, ek)
    end
end



function lower_event!(t, v::Int, c::Float64=-t.lam(v))
    p = t.parent[v]
    e = Event(p, v, 0.0, -t.y[v], 1.0)
    e.offset -= sum(map(i->t.lam(i), t.children[v]))
    set_forecast!(e, c)
    ek = find_min(t, v)
    @debug "lower_event!($v): e  = $e"
    @debug "lower_event!($v): ek = $ek"
    while ek.x < e.x
        step_min(t, ek)
        ek = find_min(t, v)
    end
    push!(t.nodes[p].pq, e)
    t.nodes[v].lb = e.x
end


function upper_event!(t, v::Int, c::Float64=+t.lam(v))
    p = t.parent[v]
    e = Event(p, v, 0.0, -t.y[v], 1.0)
    e.offset += sum(map(i->t.lam(i), t.children[v]))
    set_forecast!(e, c)
    ek = find_max(t, v)
    @debug "upper_event!($v): e  = $e"
    @debug "upper_event!($v): ek = $ek"
    while ek.x > e.x
        step_max(t, ek)
        ek = find_max(t, v)
    end
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
        @printf "\n((%d)): [% 3.1f,% 3.1f]: " i n.lb n.ub
        println(join(map(string, n.pq), "\n" * (" "^20)))
    end
    info("-"^70)
end


function forward_dp_treepwl(t)
    for i in t.pre_order[end:-1:1]
        lower_event!(t, i)
        upper_event!(t, i)
        print_tree(t)
        print_min(t)
        print_max(t)
    end
end
